import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import kmer_selector
import tempfile, os

st.set_page_config(page_title="Diseño de sondas k-mer", layout="wide")
st.title("Diseño de panel de sondas")

def save_uploaded_file(uploaded_file) -> str:
    suffix = os.path.splitext(uploaded_file.name)[1]
    with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
        tmp.write(uploaded_file.getvalue())
        return tmp.name

def stacked_or_simple_hist(df, column, title, split_col=None, nbins=40):
    fig = go.Figure()
    if split_col is not None:
        fig.add_trace(go.Histogram(x=df.loc[df[split_col], column], name="Conservados", marker_color="#1f77b4", nbinsx=nbins))
        fig.add_trace(go.Histogram(x=df.loc[~df[split_col], column], name="Eliminados", marker_color="#d62728", nbinsx=nbins))
        fig.update_layout(barmode="stack")
    else:
        fig.add_trace(go.Histogram(x=df[column], marker_color="#1f77b4", nbinsx=nbins))
    fig.update_layout(title=title, xaxis_title=column, yaxis_title="Número de sondas", bargap=0.05)
    return fig

# ─── Paso 1: Carga de datos ──────────────────────────────────
st.header("Paso 1 · Carga de datos")
col1, col2 = st.columns(2)
with col1:
    kmer_file = st.file_uploader("Tabla de k-meros (.txt/.tsv)", type=["txt", "tsv"])
with col2:
    genome_file = st.file_uploader("Tabla de genomas (.txt/.tsv)", type=["txt", "tsv"], key="genomes")

k = st.number_input("k", min_value=1, max_value=32, value=13)

if kmer_file is not None:
    kmer_path = save_uploaded_file(kmer_file)
    if st.session_state.get("kmer_path") != kmer_path:
        with st.spinner("Leyendo tabla cruda..."):
            rows = kmer_selector.read_raw(kmer_path, k)
        st.session_state["kmer_path"] = kmer_path
        st.session_state["df_raw"] = pd.DataFrame(rows, columns=["code", "gc", "entropy", "count", "files_csv"])

if genome_file is not None:
    genome_path = save_uploaded_file(genome_file)
    st.session_state["df_genomes"] = pd.read_csv(genome_path, sep="\t", header=None,
                                                   names=["id_genoma", "archivo", "header_secuencia"])

if "df_raw" in st.session_state:
    df_raw = st.session_state["df_raw"]
    c1, c2 = st.columns(2)
    c1.metric("Candidatos cargados", len(df_raw))
    if "df_genomes" in st.session_state:
        c2.metric("Genomas cargados", len(st.session_state["df_genomes"]))

    ca, cb, cc = st.columns(3)
    ca.plotly_chart(stacked_or_simple_hist(df_raw, "gc", "GC (crudo)"), use_container_width=True)
    cb.plotly_chart(stacked_or_simple_hist(df_raw, "entropy", "Entropía (crudo)"), use_container_width=True)
    cc.plotly_chart(stacked_or_simple_hist(df_raw, "count", "Frecuencia en genomas (crudo)"), use_container_width=True)
else:
    st.info("Carga la tabla de k-meros para ver estadísticas.")

st.divider()

# ─── Paso 2: Clustering ──────────────────────────────────────
st.header("Paso 2 · Eliminación de palabras cercanas (clustering)")
if "df_raw" not in st.session_state:
    st.info("Completa el Paso 1 primero.")
else:
    max_dist = st.select_slider("Distancia mínima entre palabras (max_dist)", options=[1, 2], value=1)
    run_clustering_btn = st.button("Ejecutar clustering", type="primary")

    if run_clustering_btn:
        with st.spinner("Calculando vecinos y agrupando... (puede tardar con max_dist=2)"):
            rows = kmer_selector.run_clustering(st.session_state["kmer_path"], k, max_dist)
        df_clusters = pd.DataFrame(rows, columns=["code", "gc", "entropy", "count", "component_size", "files_csv"])
        df_clusters["n_genomas"] = df_clusters["files_csv"].apply(
        lambda s: len(s.split(",")) if s else 0
    )
        st.session_state["df_clusters"] = df_clusters

    if "df_clusters" in st.session_state:
        df_clusters = st.session_state["df_clusters"]
        n_raw = len(st.session_state["df_raw"])
        n_components = len(df_clusters)
        n_discarded = n_raw - n_components

        c1, c2, c3 = st.columns(3)
        c1.metric("Candidatos originales", n_raw)
        c2.metric("Representantes (tras clustering)", n_components)
        c3.metric("Descartados por redundancia", n_discarded, delta=f"-{n_discarded/n_raw*100:.1f}%")

st.divider()

# ─── Paso 3: Filtros + panel final ───────────────────────────
st.header("Paso 3 · Filtros y selección final")
if "df_clusters" not in st.session_state:
    st.info("Completa el Paso 2 primero.")
else:
    df_clusters = st.session_state["df_clusters"].copy()

    st.subheader("3a. Filtro de calidad (opcional)")
    use_quality_filter = st.checkbox("Activar filtro de entropía mínima")
    if use_quality_filter:
        min_entropy = st.number_input("Entropía mínima", 0.0, 2.0, 1.0)
        df_clusters["passes_quality"] = df_clusters["entropy"] >= min_entropy
    else:
        df_clusters["passes_quality"] = True

    st.subheader("3b. Representatividad en genomas (opcional)")
    use_min_presence = st.checkbox("Activar presencia mínima en N genomas")
    if use_min_presence:
        min_genomes = st.number_input("Mínimo de genomas donde debe aparecer", min_value=1, value=2)
        n_genomes_per_probe = df_clusters["files_csv"].str.split(",").apply(len)
        df_clusters["passes_representativity"] = n_genomes_per_probe >= min_genomes
    else:
        df_clusters["passes_representativity"] = True

    df_clusters["passes_all"] = df_clusters["passes_quality"] & df_clusters["passes_representativity"]
    n_qualified = df_clusters["passes_all"].sum()
    st.caption(f"{n_qualified} de {len(df_clusters)} candidatos pasan los filtros activos.")

    if use_quality_filter or use_min_presence:
        ca, cb = st.columns(2)
        ca.plotly_chart(stacked_or_simple_hist(df_clusters, "gc", "GC: antes/después de filtro", split_col="passes_all"), use_container_width=True)
        cb.plotly_chart(stacked_or_simple_hist(df_clusters, "entropy", "Entropía: antes/después de filtro", split_col="passes_all"), use_container_width=True)

    st.subheader("3c. Selección final por cuotas de GC")
    target_size = st.number_input("Tamaño objetivo del panel", min_value=1, value=12000)
    strategy = st.selectbox("Estrategia de selección dentro de cada estrato", ["entropy", "count", "random"])
    deficit = st.selectbox("Manejo de déficit", ["redistribute", "strict"])
    seed_enabled = st.checkbox("Fijar semilla (resultado reproducible)", value=True)
    seed = st.number_input("Semilla", min_value=0, value=42) if seed_enabled else None

    n_bins = st.number_input("Número de estratos de GC", min_value=1, max_value=6, value=3)
    st.caption("Con 1 estrato y % = 100, equivale a un filtro de rango GC simple.")
    bins = []
    defaults = [(0.25, 0.35, 0.25), (0.35, 0.65, 0.50), (0.65, 0.75, 0.25)]
    for i in range(n_bins):
        d = defaults[i] if i < len(defaults) else (0.0, 1.0, 0.0)
        c = st.columns(3)
        bins.append((
            c[0].number_input(f"min {i+1}", value=d[0], key=f"smin_{i}"),
            c[1].number_input(f"max {i+1}", value=d[1], key=f"smax_{i}"),
            c[2].number_input(f"% {i+1}", value=d[2], key=f"sprop_{i}"),
        ))

    run_panel_btn = st.button("Seleccionar panel final", type="primary")

    if run_panel_btn:
        qualified = df_clusters[df_clusters["passes_all"]]
        rows = list(zip(
            qualified["code"], qualified["gc"], qualified["entropy"],
            qualified["count"], qualified["component_size"],
        ))
        with st.spinner("Seleccionando panel final..."):
            result = kmer_selector.select_final_panel(
                rows=rows, target_size=target_size, bins=bins,
                strategy=strategy, deficit=deficit, seed=seed,
            )
        df_panel = pd.DataFrame(result, columns=["code", "selected"])
        df_panel = df_panel.merge(qualified, on="code")
        df_panel["sequence"] = df_panel["code"].apply(lambda c: kmer_selector.decode(c, k))
        st.session_state["df_panel"] = df_panel

    if "df_panel" in st.session_state:
        df_panel = st.session_state["df_panel"]
        n_sel = df_panel["selected"].sum()
        st.metric("Panel final", n_sel)

        st.subheader("Auditoría: GC y entropía")
        ca, cb = st.columns(2)
        ca.plotly_chart(stacked_or_simple_hist(df_panel, "gc", "GC: panel final vs. descartado por cuota", split_col="selected"), use_container_width=True)
        cb.plotly_chart(stacked_or_simple_hist(df_panel, "entropy", "Entropía: panel final vs. descartado por cuota", split_col="selected"), use_container_width=True)

        st.subheader("Auditoría: uniformidad más allá de GC/entropía")
        cc, cd = st.columns(2)
        with cc:
            st.plotly_chart(
                stacked_or_simple_hist(df_panel, "component_size", "Tamaño de cluster: seleccionados vs. descartados", split_col="selected"),
                use_container_width=True
            )
        with cd:
            long_panel = df_panel.assign(id_genoma=df_panel["files_csv"].str.split(",")).explode("id_genoma")
            long_panel = long_panel[long_panel["id_genoma"] != ""]

            coverage = long_panel.groupby("id_genoma")["selected"].agg(
                total="count", seleccionadas="sum"
            ).reset_index()
            coverage["eliminadas"] = coverage["total"] - coverage["seleccionadas"]

            fig_cov = go.Figure()
            fig_cov.add_trace(go.Bar(x=coverage["id_genoma"], y=coverage["seleccionadas"], name="Seleccionadas", marker_color="#1f77b4"))
            fig_cov.add_trace(go.Bar(x=coverage["id_genoma"], y=coverage["eliminadas"], name="Eliminadas", marker_color="#d62728"))
            fig_cov.update_layout(barmode="stack", title="Cobertura de sondas por genoma (panel final)",
                                   xaxis_title="Genoma (ID)", yaxis_title="Sondas")
            st.plotly_chart(fig_cov, use_container_width=True)

        st.subheader("Exportar resultados")
        final = df_panel[df_panel["selected"]].copy()

        # Archivo 1: mismo formato de la tabla original, con n_genomas corregido
        full_table = final[["sequence", "code", "n_genomas", "gc", "entropy", "files_csv"]]
        full_tsv = full_table.to_csv(sep="\t", index=False, header=False)

        # Archivo 2: solo la lista de secuencias, sin columnas extra
        dataset_txt = "\n".join(final["sequence"].tolist())

        col1, col2 = st.columns(2)
        with col1:
            st.download_button(
                "Descargar tabla completa (.tsv)",
                data=full_tsv,
                file_name="probes_selected.tsv",
                mime="text/tab-separated-values",
            )
        with col2:
            st.download_button(
                "Descargar lista de secuencias (dataset.txt)",
                data=dataset_txt,
                file_name="dataset.txt",
                mime="text/plain",
            )