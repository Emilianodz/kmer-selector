use kmer_selector::{io, clustering, export};
use rand::seq::SliceRandom;
use rand::{SeedableRng, rngs::StdRng};
use std::io::Write;
use std::time::Instant;

const SAFETY_DEGREE_THRESHOLD: f64 = 0.3;

fn universe_size(k: usize) -> f64 {
    4f64.powi(k as i32)
}

fn binom(n: usize, r: usize) -> f64 {
    if r > n { return 0.0; }
    let mut result = 1.0;
    for i in 0..r {
        result *= (n - i) as f64;
        result /= (i + 1) as f64;
    }
    result
}

fn neighbors_per_word(k: usize, max_dist: u32) -> f64 {
    let mut total = 0.0;
    for d in 1..=max_dist {
        total += binom(k, d as usize) * 3f64.powi(d as i32);
    }
    total
}

fn safe_batch_size(k: usize, max_dist: u32) -> usize {
    let universe = universe_size(k);
    let neighbors = neighbors_per_word(k, max_dist);
    (SAFETY_DEGREE_THRESHOLD * universe / neighbors) as usize
}

// Formatea números grandes con separador de miles: 30292461 -> "30,292,461"
fn fmt_num(n: usize) -> String {
    let s = n.to_string();
    let bytes = s.as_bytes();
    let mut out = String::new();
    for (i, b) in bytes.iter().enumerate() {
        if i > 0 && (bytes.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(*b as char);
    }
    out
}

struct PassSummary {
    pass_num: usize,
    input: usize,
    output: usize,
}

impl PassSummary {
    fn reduction_pct(&self) -> f64 {
        (1.0 - self.output as f64 / self.input as f64) * 100.0
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!(
            "Uso: batch_cluster <archivo_entrada> [k=13] [max_dist=1] [seed=42] [batch_size=auto]"
        );
        std::process::exit(1);
    }

    let input_path = &args[1];
    let k: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(13);
    let max_dist: u32 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1);
    let seed: u64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(42);
    let batch_size: usize = args.get(5)
        .and_then(|s| s.parse().ok())
        .unwrap_or_else(|| safe_batch_size(k, max_dist));

    println!("Tamaño de lote: {} candidatos (grado promedio objetivo <= {})",
        fmt_num(batch_size), SAFETY_DEGREE_THRESHOLD);
    std::io::stdout().flush().ok();

    std::fs::create_dir_all("data").ok();

    let mut candidates = io::read_candidates(input_path, k).unwrap_or_else(|e| {
        eprintln!("Error al leer {}: {}", input_path, e);
        std::process::exit(1);
    });

    let input_total = candidates.len();
    let mut rng = StdRng::seed_from_u64(seed);
    let mut pass_num = 1;
    let mut pass_summaries: Vec<PassSummary> = Vec::new();
    let mut stop_reason = "convergencia (cabe en un solo lote seguro)";

    loop {
        let pass_start = Instant::now();
        let pass_input = candidates.len();
        println!("\n=== Pasada {} -- {} candidatos de entrada ===", pass_num, fmt_num(pass_input));
        std::io::stdout().flush().ok();

        candidates.shuffle(&mut rng);

        let n_batches = (candidates.len() + batch_size - 1) / batch_size;
        let mut representatives: Vec<io::ProbeCandidate> = Vec::new();

        // Progreso: imprime cada lote si hay pocos, o cada ~5% si hay muchos.
        let print_every = (n_batches / 20).max(1);

        for (batch_idx, batch) in candidates.chunks(batch_size).enumerate() {
            let components = clustering::select_representatives(batch, k, max_dist);
            for comp in &components {
                representatives.push(clustering::representative_as_candidate(comp));
            }

            if (batch_idx + 1) % print_every == 0 || batch_idx + 1 == n_batches {
                let elapsed = pass_start.elapsed().as_secs_f64();
                let pace = elapsed / (batch_idx + 1) as f64;
                let remaining = ((n_batches - batch_idx - 1) as f64 * pace).round() as u64;
                print!(
                    "\r  Lote {}/{} | representantes: {} | transcurrido: {:.0}s | ETA: {}s   ",
                    batch_idx + 1, n_batches, fmt_num(representatives.len()), elapsed, remaining
                );
                std::io::stdout().flush().ok();
            }
        }
        println!(); // salto de línea tras la barra de progreso en carriage return

        let intermediate_path = format!("data/intermediate_pass{}.tsv", pass_num);
        write_candidate_table(&representatives, k, &intermediate_path);

        let summary = PassSummary { pass_num, input: pass_input, output: representatives.len() };
        println!(
            "--- Resumen pasada {} ---\nEntrada:  {:>12}\nSalida:   {:>12}\nReducción: {:>10.2} %\nGuardado en: {}",
            summary.pass_num, fmt_num(summary.input), fmt_num(summary.output),
            summary.reduction_pct(), intermediate_path
        );

        let reduction_ratio = representatives.len() as f64 / candidates.len() as f64;
        candidates = representatives;
        pass_summaries.push(summary);

        if candidates.len() <= batch_size {
            stop_reason = "convergencia (cabe en un solo lote seguro)";
            break;
        }
        if reduction_ratio > 0.95 {
            stop_reason = "estancamiento (reducción < 5% en la última pasada)";
            break;
        }
        pass_num += 1;
    }

    println!("\n=== Clustering final sobre {} candidatos ===", fmt_num(candidates.len()));
    let final_components = clustering::select_representatives(&candidates, k, max_dist);
    println!("Componentes finales: {}", fmt_num(final_components.len()));

    // Archivo principal: 6 columnas, mismo formato que kmer_counter --
    // compatible directo con io::read_candidates y con la app de Python.
    let final_candidates: Vec<io::ProbeCandidate> = final_components.iter()
        .map(clustering::representative_as_candidate)
        .collect();
    write_candidate_table(&final_candidates, k, "data/probes_selected.tsv");

    // Archivo adicional (opcional): con component_size y members, para
    // auditoría manual -- NO pensado para recargarse en Python.
    if let Err(e) = export::export_components(&final_components, k, "data/probes_selected_detailed.tsv") {
        eprintln!("Aviso: no se pudo generar el archivo detallado: {}", e);
    }

    write_final_report(
        input_total, final_candidates.len(), k, max_dist, seed,
        batch_size, &pass_summaries, stop_reason,
    );
}

// Formato de 6 columnas: sequence, code, count, gc_content, entropy, files.
// Mismo layout que exporta kmer_counter -- io::read_candidates lo lee sin
// modificaciones, y la app de Streamlit lo carga igual que cualquier TSV
// de entrada normal.
fn write_candidate_table(candidates: &[io::ProbeCandidate], k: usize, path: &str) {
    let mut file = std::fs::File::create(path).expect("no se pudo crear el archivo de salida");
    for c in candidates {
        let seq = io::decode_kmer(c.code, k);
        let files_str = c.files.iter().map(|f| f.to_string()).collect::<Vec<_>>().join(",");
        writeln!(
            file, "{}\t{}\t{}\t{:.4}\t{:.4}\t{}",
            seq, c.code, c.count, c.gc_content, c.entropy, files_str
        ).unwrap();
    }
}

fn write_final_report(
    input_total: usize,
    final_total: usize,
    k: usize,
    max_dist: u32,
    seed: u64,
    batch_size: usize,
    passes: &[PassSummary],
    stop_reason: &str,
) {
    let total_reduction = (1.0 - final_total as f64 / input_total as f64) * 100.0;
    let initial_reduction = passes.first().map(|p| p.reduction_pct()).unwrap_or(0.0);
    let final_batch_reduction = passes.last().map(|p| p.reduction_pct()).unwrap_or(0.0);

    let mut report = String::new();
    report.push_str("KMER CLUSTERING REPORT\n");
    report.push_str("======================\n");
    report.push_str(&format!("Input candidates:       {:>12}\n", fmt_num(input_total)));
    report.push_str(&format!("Final candidates:       {:>12}\n", fmt_num(final_total)));
    report.push_str(&format!("Total reduction:        {:>11.3} %\n", total_reduction));
    report.push_str(&format!("k:                      {:>12}\n", k));
    report.push_str(&format!("Maximum Hamming distance:{:>11}\n", max_dist));
    report.push_str(&format!("Seed:                    {:>12}\n", seed));
    report.push_str(&format!("Batch size:              {:>12}\n", fmt_num(batch_size)));
    report.push_str(&format!("Passes:                  {:>12}\n", passes.len()));
    report.push_str(&format!("Final global clustering: {:>12}\n", "yes"));
    report.push_str(&format!("Stop reason:             {}\n", stop_reason));
    report.push_str(&format!("Initial pass reduction:  {:>10.1} %\n", initial_reduction));
    report.push_str(&format!("Final pass reduction:    {:>10.1} %\n", final_batch_reduction));
    report.push_str("\nPer-pass detail:\n");
    for p in passes {
        report.push_str(&format!(
            "  pass {:>3}: {:>12} -> {:>12}  ({:>6.2} %)\n",
            p.pass_num, fmt_num(p.input), fmt_num(p.output), p.reduction_pct()
        ));
    }

    println!("\n{}", report);
    std::fs::write("data/clustering_report.txt", &report).ok();
    println!("Reporte guardado en data/clustering_report.txt");
}