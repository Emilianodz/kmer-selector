# kmer_selector

A Rust engine for redundancy reduction and quota-based selection of candidate k-mer probes, built as the downstream companion to `kmer_counter`.

Rust library for clustering near-duplicate k-mers (Hamming-distance neighborhoods) and selecting a final probe panel under configurable GC/entropy/genome-representativity constraints. Designed as a reusable engine: the core logic lives in `lib.rs` and is exposed both as a public Rust API and as a Python extension module (via PyO3/maturin) for interactive exploration in a Streamlit application.

## Project status

- [x] TSV/TXT ingestion of `kmer_counter` output (`io.rs`)
- [x] Configurable-radius Hamming-neighbor generation (`neighbors.rs`)
- [x] Union-Find clustering with deterministic tie-breaking (`clustering.rs`)
- [x] GC/entropy quality filter (`filter.rs`)
- [x] Quota-based panel selection with GC stratification and weighted sampling (`selection.rs`)
- [x] TSV export, including per-component metadata (`export.rs`)
- [x] Parallelized neighbor search with `rayon`
- [x] Python bindings via PyO3/maturin
- [x] Interactive Streamlit application for auditing each pipeline stage
- [ ] BK-tree for arbitrary-radius interactive queries (optional, deferred)

## Pipeline

raw candidates → clustering (redundancy removal) → optional filters
(GC / entropy / genome representativity) → quota-based panel selection


Clustering runs on the full raw candidate set (not pre-filtered), so redundancy removal is not biased by quality thresholds decided afterward. Filters and quota selection are applied downstream, on the much smaller set of cluster representatives.

## Requirements

- Rust (Edition 2021) and Cargo
- Python 3.11–3.13 (PyO3 0.22 does not yet support 3.14+)
- [maturin](https://www.maturin.rs/) for building the Python extension

## Usage — Rust library / CLI

```bash
cargo test                          # Run the unit test suite
cargo run --bin kmer_selector_cli -- ./data/kmer_table.txt
```

## Usage — Python bindings + Streamlit app

```bash
python3 -m venv .venv
source .venv/bin/activate           # Linux/macOS
# .venv\Scripts\Activate.ps1        # Windows

pip install -r app/requirements.txt
maturin develop                     # builds the Rust extension into the venv

streamlit run app/streamlit_app.py
```

## Project structure

kmer_selector/
├── Cargo.toml
├── Cargo.lock
├── src/
│ ├── lib.rs # PyO3 module + public Rust API
│ ├── main.rs # standalone CLI binary
│ ├── io.rs # TSV parsing, encoding/decoding, GC/entropy
│ ├── neighbors.rs # Hamming-neighbor generation
│ ├── clustering.rs # Union-Find + representative selection
│ ├── filter.rs # GC/entropy quality filter
│ ├── selection.rs # quota-based panel selection
│ └── export.rs # TSV export
├── app/
│ ├── streamlit_app.py # interactive exploration UI
│ └── requirements.txt
└── data/ # not tracked in git — see .gitignore


## Output files

- **`probes_selected.tsv`**: full candidate table with an extra `selected` column, produced by the Streamlit app's export step. Same column layout as `kmer_counter`'s TSV output plus selection metadata.
- **`dataset.txt`**: plain list of selected probe sequences, one per line, with no additional columns.

## Author

Developed by Emiliano Díaz-López