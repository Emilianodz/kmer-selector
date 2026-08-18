use std::collections::HashSet;
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use crate::clustering::Component;

pub struct GcBin {
    pub min_gc: f32,
    pub max_gc: f32,
    pub proportion: f32,
}

#[derive(Clone, Copy)]
pub enum SelectionStrategy {
    Random,
    WeightedByEntropy,
    WeightedByCount,
}

#[derive(Clone, Copy)]
pub enum DeficitStrategy {
    Strict,
    Redistribute,
}

pub struct PanelConfig {
    pub target_size: usize,
    pub bins: Vec<GcBin>,
    pub strategy: SelectionStrategy,
    pub deficit: DeficitStrategy,
    pub seed: Option<u64>,
}

pub struct PanelResult {
    pub selected: Vec<u64>,
    pub eliminated: Vec<u64>,
}

pub fn select_panel(components: &[Component], config: &PanelConfig) -> PanelResult {
    let mut rng: StdRng = match config.seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_entropy(),
    };

    let mut bin_members: Vec<Vec<&Component>> = vec![Vec::new(); config.bins.len()];
    for c in components {
        if let Some(idx) = config.bins.iter().position(|b| {
            c.representative_gc >= b.min_gc && c.representative_gc < b.max_gc
        }) {
            bin_members[idx].push(c);
        }
    }

    let mut quotas: Vec<usize> = config.bins.iter()
        .map(|b| (config.target_size as f32 * b.proportion).round() as usize)
        .collect();

    if matches!(config.deficit, DeficitStrategy::Redistribute) {
        redistribute_deficit(&mut quotas, &bin_members);
    }

    let mut selected_ids: HashSet<u64> = HashSet::new();
    for (i, members) in bin_members.iter().enumerate() {
        let quota = quotas[i].min(members.len());
        for c in weighted_sample(members, quota, config.strategy, &mut rng) {
            selected_ids.insert(c.representative);
        }
    }

    let selected = components.iter()
        .filter(|c| selected_ids.contains(&c.representative))
        .map(|c| c.representative)
        .collect();
    let eliminated = components.iter()
        .filter(|c| !selected_ids.contains(&c.representative))
        .map(|c| c.representative)
        .collect();

    PanelResult { selected, eliminated }
}

fn redistribute_deficit(quotas: &mut [usize], bin_members: &[Vec<&Component>]) {
    let mut deficit = 0usize;
    let mut surplus: Vec<usize> = vec![0; quotas.len()];

    for i in 0..quotas.len() {
        let available = bin_members[i].len();
        if available < quotas[i] {
            deficit += quotas[i] - available;
        } else {
            surplus[i] = available - quotas[i];
        }
    }

    let total_surplus: usize = surplus.iter().sum();
    if deficit > 0 && total_surplus > 0 {
        for i in 0..quotas.len() {
            if surplus[i] > 0 {
                let share = (deficit as f32 * surplus[i] as f32 / total_surplus as f32).round() as usize;
                quotas[i] += share.min(surplus[i]);
            }
        }
    }
}

fn weighted_sample<'a>(
    members: &[&'a Component],
    quota: usize,
    strategy: SelectionStrategy,
    rng: &mut impl Rng,
) -> Vec<&'a Component> {
    let weight_of = |c: &Component| -> f32 {
        match strategy {
            SelectionStrategy::Random => 1.0,
            SelectionStrategy::WeightedByEntropy => c.representative_entropy.max(0.001),
            SelectionStrategy::WeightedByCount => c.representative_count as f32,
        }
    };

    let mut keyed: Vec<(f32, &Component)> = members.iter()
        .map(|&c| {
            let u: f32 = rng.gen_range(0.0001..1.0);
            let key = u.powf(1.0 / weight_of(c));
            (key, c)
        })
        .collect();

    keyed.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    keyed.into_iter().take(quota).map(|(_, c)| c).collect()
}


#[cfg(test)]
mod tests {
    use super::*;

    fn make_component(code: u64, gc: f32, entropy: f32, count: u32) -> Component {
        Component {
            representative: code,
            representative_count: count,
            representative_gc: gc,
            representative_entropy: entropy,
            size: 1,
            members: vec![code],
            files: Vec::new(),
        }
    }

    #[test]
    fn respects_quota_per_bin() {
        let components: Vec<Component> = (0..100)
            .map(|i| make_component(i, 0.5, 1.5, 1))
            .collect();

        let config = PanelConfig {
            target_size: 20,
            bins: vec![GcBin { min_gc: 0.0, max_gc: 1.0, proportion: 1.0 }],
            strategy: SelectionStrategy::Random,
            deficit: DeficitStrategy::Strict,
            seed: Some(42),
        };

        let result = select_panel(&components, &config);
        assert_eq!(result.selected.len(), 20);
        assert_eq!(result.eliminated.len(), 80);
    }

    #[test]
    fn redistributes_deficit_to_surplus_bin() {
        let mut components: Vec<Component> = (0..5)
            .map(|i| make_component(i, 0.2, 1.5, 1))
            .collect();
        components.extend((100..200).map(|i| make_component(i, 0.8, 1.5, 1)));

        let config = PanelConfig {
            target_size: 50,
            bins: vec![
                GcBin { min_gc: 0.0, max_gc: 0.5, proportion: 0.5 },
                GcBin { min_gc: 0.5, max_gc: 1.0, proportion: 0.5 },
            ],
            strategy: SelectionStrategy::Random,
            deficit: DeficitStrategy::Redistribute,
            seed: Some(42),
        };

        let result = select_panel(&components, &config);
        assert_eq!(result.selected.len(), 50);
    }

    #[test]
    fn strict_deficit_yields_smaller_panel() {
        let components: Vec<Component> = (0..5)
            .map(|i| make_component(i, 0.2, 1.5, 1))
            .collect();

        let config = PanelConfig {
            target_size: 20,
            bins: vec![GcBin { min_gc: 0.0, max_gc: 1.0, proportion: 1.0 }],
            strategy: SelectionStrategy::Random,
            deficit: DeficitStrategy::Strict,
            seed: Some(42),
        };

        let result = select_panel(&components, &config);
        assert_eq!(result.selected.len(), 5);
    }

    #[test]
    fn weighted_by_entropy_strongly_prefers_high_entropy() {
        let mut components = vec![make_component(0, 0.5, 0.01, 1)];
        components.extend((1..50).map(|i| make_component(i, 0.5, 2.0, 1)));

        let config = PanelConfig {
            target_size: 5,
            bins: vec![GcBin { min_gc: 0.0, max_gc: 1.0, proportion: 1.0 }],
            strategy: SelectionStrategy::WeightedByEntropy,
            deficit: DeficitStrategy::Strict,
            seed: Some(1),
        };

        let result = select_panel(&components, &config);
        assert!(!result.selected.contains(&0));
    }

    #[test]
    fn reproducible_with_same_seed() {
        let components: Vec<Component> = (0..100)
            .map(|i| make_component(i, 0.5, 1.5, 1))
            .collect();

        let config = PanelConfig {
            target_size: 30,
            bins: vec![GcBin { min_gc: 0.0, max_gc: 1.0, proportion: 1.0 }],
            strategy: SelectionStrategy::Random,
            deficit: DeficitStrategy::Strict,
            seed: Some(7),
        };

        let r1 = select_panel(&components, &config);
        let r2 = select_panel(&components, &config);
        assert_eq!(r1.selected, r2.selected);
    }
}