use pyo3::prelude::*;
use pyo3::types::PyModule;
use pyo3::Bound;
use std::collections::HashSet;

pub mod io;
pub mod neighbors;
pub mod clustering;
pub mod filter;
pub mod export;
pub mod selection;

use selection::{PanelConfig, GcBin, SelectionStrategy, DeficitStrategy};

/// Paso 1: lectura cruda, sin filtro ni clustering. Solo para estadísticas
/// descriptivas iniciales sobre el universo completo de candidatos.
#[pyfunction]
fn read_raw(path: String, k: usize) -> PyResult<Vec<(u64, f32, f32, u32, String)>> {
    // (code, gc, entropy, count, files_csv)
    let candidates = io::read_candidates(&path, k)
        .map_err(pyo3::exceptions::PyIOError::new_err)?;

    Ok(candidates.iter().map(|c| {
        let files_csv = c.files.iter().map(|f| f.to_string()).collect::<Vec<_>>().join(",");
        (c.code, c.gc_content, c.entropy, c.count, files_csv)
    }).collect())
}

/// Paso 2: clustering directo sobre TODOS los candidatos crudos -- ya no
/// depende de un filtro previo. Regresa un renglón por componente
/// (representante ya elegido por mayor entropía, con tamaño del cluster).
#[pyfunction]
fn run_clustering(path: String, k: usize, max_dist: u32)
    -> PyResult<Vec<(u64, f32, f32, u32, usize, String)>>
{
    // (representative_code, gc, entropy, count, component_size, files_csv)
    let candidates = io::read_candidates(&path, k)
        .map_err(pyo3::exceptions::PyIOError::new_err)?;
    let components = clustering::select_representatives(&candidates, k, max_dist);

    Ok(components.iter().map(|c| {
        let files_csv = c.files.iter().map(|f| f.to_string()).collect::<Vec<_>>().join(",");
        (c.representative, c.representative_gc, c.representative_entropy,
         c.representative_count, c.size, files_csv)
    }).collect())
}

/// Paso 3: selección final por cuotas, sobre representantes YA filtrados
/// en Python (GC/entropía/representatividad) -- `rows` llega más pequeño
/// que el universo crudo, así que reconstruir Component aquí es barato.
#[pyfunction]
#[pyo3(signature = (rows, target_size, bins, strategy, deficit, seed=None))]
fn select_final_panel(
    rows: Vec<(u64, f32, f32, u32, usize)>, // (code, gc, entropy, count, component_size)
    target_size: usize,
    bins: Vec<(f32, f32, f32)>,
    strategy: String,
    deficit: String,
    seed: Option<u64>,
) -> PyResult<Vec<(u64, bool)>> {
    let components: Vec<clustering::Component> = rows.into_iter()
        .map(|(code, gc, entropy, count, size)| clustering::Component {
            representative: code,
            representative_count: count,
            representative_gc: gc,
            representative_entropy: entropy,
            size,
            members: vec![code],
            files: Vec::new(),
        })
        .collect();

    let strategy = match strategy.as_str() {
        "random" => SelectionStrategy::Random,
        "entropy" => SelectionStrategy::WeightedByEntropy,
        "count" => SelectionStrategy::WeightedByCount,
        other => return Err(pyo3::exceptions::PyValueError::new_err(
            format!("estrategia desconocida: {}", other)
        )),
    };
    let deficit = match deficit.as_str() {
        "redistribute" => DeficitStrategy::Redistribute,
        "strict" => DeficitStrategy::Strict,
        other => return Err(pyo3::exceptions::PyValueError::new_err(
            format!("estrategia de déficit desconocida: {}", other)
        )),
    };
    let gc_bins: Vec<GcBin> = bins.into_iter()
        .map(|(min_gc, max_gc, proportion)| GcBin { min_gc, max_gc, proportion })
        .collect();

    let config = PanelConfig { target_size, bins: gc_bins, strategy, deficit, seed };
    let panel = selection::select_panel(&components, &config);
    let selected_set: HashSet<u64> = panel.selected.into_iter().collect();

    Ok(components.iter()
        .map(|c| (c.representative, selected_set.contains(&c.representative)))
        .collect())
}

/// Decodifica un código u64 a secuencia de texto.
#[pyfunction]
fn decode(code: u64, k: usize) -> String {
    io::decode_kmer(code, k)
}

#[pymodule]
fn kmer_selector(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read_raw, m)?)?;
    m.add_function(wrap_pyfunction!(run_clustering, m)?)?;
    m.add_function(wrap_pyfunction!(select_final_panel, m)?)?;
    m.add_function(wrap_pyfunction!(decode, m)?)?;
    Ok(())
}