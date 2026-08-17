use std::collections::HashSet;

/// Genera todos los códigos que difieren de `code` en 1 hasta `max_dist`
pub fn neighbors_within(code: u64, k: usize, max_dist: u32) -> Vec<u64> {
    let mut neighbors = HashSet::new();
    for dist in 1..=(max_dist as usize) {
        for positions in positions_combinations(k, dist) {
            apply_substitutions(code, k, &positions, &mut neighbors);
        }
    }
    neighbors.into_iter().collect()
}

/// Todas las combinaciones de `dist` posiciones distintas dentro de 0..k
fn positions_combinations(k: usize, dist: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    let mut current = Vec::with_capacity(dist);
    combine(0, k, dist, &mut current, &mut result);
    result
}

fn combine(start: usize, k: usize, dist: usize, current: &mut Vec<usize>, result: &mut Vec<Vec<usize>>) {
    if current.len() == dist {
        result.push(current.clone());
        return;
    }
    for pos in start..k {
        current.push(pos);
        combine(pos + 1, k, dist, current, result);
        current.pop();
    }
}

/// Para un conjunto de posiciones dado, genera todos los códigos resultantes
/// de sustituir cada posición por sus 3 bases alternativas
fn apply_substitutions(code: u64, k: usize, positions: &[usize], out: &mut HashSet<u64>) {
    substitute_recursive(code, k, positions, 0, out);
}

fn substitute_recursive(code: u64, k: usize, positions: &[usize], idx: usize, out: &mut HashSet<u64>) {
    if idx == positions.len() {
        out.insert(code);
        return;
    }
    let pos = positions[idx];
    let shift = 2 * (k - 1 - pos);
    let original_base = ((code >> shift) & 0b11) as u8;

    for alt_base in (0u8..4).filter(|&b| b != original_base) {
        let mut new_code = code & !(0b11u64 << shift);
        new_code |= (alt_base as u64) << shift;
        substitute_recursive(new_code, k, positions, idx + 1, out);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_39_neighbors_at_distance_1_for_k13() {
        let code: u64 = 0; // todo A
        let neighbors = neighbors_within(code, 13, 1);
        assert_eq!(neighbors.len(), 39);
    }

    #[test]
    fn counts_741_neighbors_up_to_distance_2_for_k13() {
        let code: u64 = 0;
        let neighbors = neighbors_within(code, 13, 2);
        assert_eq!(neighbors.len(), 741); // 39 + 702
    }

    #[test]
    fn does_not_include_original_code() {
        let code: u64 = 0b1011_0100; // algún k-mero arbitrario, k=4
        let neighbors = neighbors_within(code, 4, 2);
        assert!(!neighbors.contains(&code));
    }

    #[test]
    fn no_duplicates_across_distances() {
        let code: u64 = 0;
        let neighbors = neighbors_within(code, 13, 2);
        let as_set: HashSet<u64> = neighbors.iter().copied().collect();
        assert_eq!(neighbors.len(), as_set.len());
    }

    #[test]
    fn combine_generates_expected_position_pairs() {
        let combos = positions_combinations(4, 2);
        assert_eq!(combos.len(), 6); // C(4,2) = 6
        assert!(combos.contains(&vec![0, 1]));
        assert!(combos.contains(&vec![2, 3]));
    }
}