use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
use crate::io::ProbeCandidate;
use crate::neighbors::neighbors_within;
use rayon::prelude::*;


struct UnionFind {
    parent: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind { parent: (0..n).collect() }
    }

    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }

    fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra != rb {
            self.parent[ra] = rb;
        }
    }
}

pub struct Component {
    pub representative: u64,
    pub representative_count: u32,
    pub representative_gc: f32,
    pub representative_entropy: f32,
    pub size: usize,
    pub members: Vec<u64>,
    pub files: Vec<u32>,
}

// Tamaño de bloque para procesar candidatos en fases: acota el pico de
// memoria del buffer temporal de aristas, en vez de acumular las aristas
// de TODOS los candidatos (potencialmente millones) antes de aplicar
// ninguna unión. Ajustable según la RAM disponible.
const CHUNK_SIZE: usize = 500_000;

pub fn select_representatives(
    candidates: &[ProbeCandidate],
    k: usize,
    max_dist: u32,
) -> Vec<Component> {
    let code_to_idx: HashMap<u64, usize> = candidates
        .iter()
        .enumerate()
        .map(|(i, c)| (c.code, i))
        .collect();

    let mut uf = UnionFind::new(candidates.len());

    for (chunk_idx, chunk) in candidates.chunks(CHUNK_SIZE).enumerate() {
        let chunk_start = chunk_idx * CHUNK_SIZE;

        // Fase paralela SOLO sobre este bloque: cada candidato genera sus
        // vecinos y busca matches, de forma independiente entre hilos.
        let code_to_idx_ref = &code_to_idx;
        let edges: Vec<(usize, usize)> = chunk
            .par_iter()
            .enumerate()
            .flat_map_iter(move |(local_i, candidate)| {
                let i = chunk_start + local_i;
                neighbors_within(candidate.code, k, max_dist)
                    .into_iter()
                    .filter_map(move |neighbor_code| {
                        code_to_idx_ref.get(&neighbor_code).map(|&j| (i, j))
                    })
            })
            .collect();

        // Fase secuencial: aplicar las uniones de este bloque.
        for (i, j) in edges {
            uf.union(i, j);
        }
        // `edges` sale de scope aquí -- su memoria se libera antes de
        // pasar al siguiente bloque, en vez de acumularse para todo el dataset.
    }

    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..candidates.len() {
        let root = uf.find(i);
        groups.entry(root).or_default().push(i);
    }

    groups
        .into_values()
        .map(|indices| build_component(candidates, &indices))
        .collect()
}

fn build_component(candidates: &[ProbeCandidate], indices: &[usize]) -> Component {
    let representative_idx = *indices
        .iter()
        .max_by(|&&a, &&b| {
            let ca = &candidates[a];
            let cb = &candidates[b];
            ca.entropy
                .partial_cmp(&cb.entropy)
                .unwrap()
                .then_with(|| hash_of(ca.code).cmp(&hash_of(cb.code)))
        })
        .unwrap();

    let rep = &candidates[representative_idx];

    let mut files: Vec<u32> = indices.iter()
        .flat_map(|&i| candidates[i].files.iter().copied())
        .collect();
    files.sort_unstable();
    files.dedup();

    Component {
        representative: rep.code,
        representative_count: rep.count,
        representative_gc: rep.gc_content,
        representative_entropy: rep.entropy,
        size: indices.len(),
        members: indices.iter().map(|&i| candidates[i].code).collect(),
        files,
    }
}

fn hash_of(code: u64) -> u64 {
    let mut hasher = DefaultHasher::new();
    code.hash(&mut hasher);
    hasher.finish()
}


pub fn representative_as_candidate(comp: &Component) -> ProbeCandidate {
    ProbeCandidate {
        code: comp.representative,
        count: comp.representative_count,
        gc_content: comp.representative_gc,
        entropy: comp.representative_entropy,
        files: comp.files.clone(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_candidate(code: u64, entropy: f32) -> ProbeCandidate {
        ProbeCandidate { code, count: 1, gc_content: 0.5, entropy, files: Vec::new() }
    }

    #[test]
    fn union_find_basic_union() {
        let mut uf = UnionFind::new(5);
        uf.union(0, 1);
        assert_eq!(uf.find(0), uf.find(1));
        assert_ne!(uf.find(0), uf.find(2));
    }

    #[test]
    fn transitive_chain_ends_in_same_component() {
        let a: u64 = 0b00_00_00_00_00_00_00;
        let b: u64 = 0b00_00_00_00_00_00_01;
        let c: u64 = 0b00_00_00_00_00_01_01;

        let candidates = vec![
            make_candidate(a, 1.0),
            make_candidate(b, 1.0),
            make_candidate(c, 1.0),
        ];
        let components = select_representatives(&candidates, 7, 1);
        assert_eq!(components.len(), 1);
        assert_eq!(components[0].size, 3);
    }

    #[test]
    fn isolated_candidates_form_singleton_components() {
        let a: u64 = 0;
        let far: u64 = u64::MAX;
        let candidates = vec![make_candidate(a, 1.0), make_candidate(far, 1.0)];
        let components = select_representatives(&candidates, 13, 1);
        assert_eq!(components.len(), 2);
    }

    #[test]
    fn selects_highest_entropy_as_representative() {
        let a: u64 = 0b00_00_00_00_00_00_00;
        let b: u64 = 0b00_00_00_00_00_00_01;
        let candidates = vec![
            make_candidate(a, 0.5),
            make_candidate(b, 2.0),
        ];
        let components = select_representatives(&candidates, 7, 1);
        assert_eq!(components.len(), 1);
        assert_eq!(components[0].representative, b);
    }

    #[test]
    fn representative_selection_is_order_independent() {
        let a: u64 = 0b00_00_00_00_00_00_00;
        let b: u64 = 0b00_00_00_00_00_00_01;
        let candidates_1 = vec![make_candidate(a, 1.5), make_candidate(b, 1.5)];
        let candidates_2 = vec![make_candidate(b, 1.5), make_candidate(a, 1.5)];

        let c1 = select_representatives(&candidates_1, 7, 1);
        let c2 = select_representatives(&candidates_2, 7, 1);
        assert_eq!(c1[0].representative, c2[0].representative);
    }
}