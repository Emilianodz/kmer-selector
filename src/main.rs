mod io;
mod neighbors;
mod clustering;
mod filter;
mod export;
mod selection;

use filter::FilterConfig;
use selection::{PanelConfig, GcBin, SelectionStrategy, DeficitStrategy};
use kmer_selector::{io, clustering, filter, export, selection};
use filter::FilterConfig;
use selection::{PanelConfig, GcBin, SelectionStrategy, DeficitStrategy};

// ─── Configuración centralizada ───────────────────────────────
const K: usize = 13;
const MAX_DIST: u32 = 2;
const DEFAULT_INPUT_TSV: &str = "data/kmers.tsv";
const OUTPUT_TSV: &str = "data/kmers_selected.tsv";
const PANEL_OUTPUT_TSV: &str = "data/final_panel.tsv";

fn build_filter_config() -> FilterConfig {
    FilterConfig {
        min_gc: Some(0.25),
        max_gc: Some(0.75),
        min_entropy: None,
    }
}

fn build_panel_config() -> PanelConfig {
    PanelConfig {
        target_size: 12000,
        bins: vec![
            GcBin { min_gc: 0.25, max_gc: 0.35, proportion: 0.25 },
            GcBin { min_gc: 0.35, max_gc: 0.65, proportion: 0.50 },
            GcBin { min_gc: 0.65, max_gc: 0.75, proportion: 0.25 },
        ],
        strategy: SelectionStrategy::WeightedByEntropy,
        deficit: DeficitStrategy::Redistribute,
        seed: Some(42), // fijo mientras validas; None para producción real
    }
}
// ────────────────────────────────────────────────────────────

fn main() {
    let input_path = std::env::args().nth(1).unwrap_or_else(|| DEFAULT_INPUT_TSV.to_string());

    let candidates = match io::read_candidates(&input_path, K) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Error al leer {}: {}", input_path, e);
            std::process::exit(1);
        }
    };
    println!("Candidatos leídos: {}", candidates.len());

    let filter_config = build_filter_config();
    let filtered = filter::apply_filter(candidates, &filter_config);
    println!("Candidatos tras filtro GC/entropía: {}", filtered.len());

    let start = std::time::Instant::now();
    let components = clustering::select_representatives(&filtered, K, MAX_DIST);
    println!("Clustering: {} componentes ({:?})", components.len(), start.elapsed());

    match export::export_components(&components, K, OUTPUT_TSV) {
        Ok(()) => println!("Componentes exportados a {}", OUTPUT_TSV),
        Err(e) => { eprintln!("Error al exportar componentes: {}", e); std::process::exit(1); }
    }

    let panel_config = build_panel_config();
    let panel = selection::select_panel(&components, &panel_config);
    println!(
        "Panel final: {} sondas seleccionadas, {} eliminadas por cuota",
        panel.selected.len(), panel.eliminated.len()
    );

    match export::export_panel(&panel, &components, K, PANEL_OUTPUT_TSV) {
        Ok(()) => println!("Panel final exportado a {}", PANEL_OUTPUT_TSV),
        Err(e) => { eprintln!("Error al exportar panel: {}", e); std::process::exit(1); }
    }
}