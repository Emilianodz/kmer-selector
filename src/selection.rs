use std::collections::HashSet;
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use crate::clustering::Component;

pub struct GcBin {
    pub min_gc: f32,
    pub max_gc: f32,
    pub proportion: f32, //target_size, ej. 0.5 = 50%
}

#[derive(Clone, Copy)]
pub enum SelectionStrategy {
    Random,
    WeightedByEntropy,
    WeightedByCount,
}

#[derive(Clone, Copy)]
pub enum DeficitStrategy {
    Strict,        // panel final puede quedar incompleto que target_size seleccionado
    Redistribute,  // el caso se utocompleta si hay faltantes con los segmentos excedentes
}

pub struct PanelConfig {
    pub target_size: usize,
    pub bins: Vec<GcBin>,
    pub strategy: SelectionStrategy,
    pub deficit: DeficitStrategy,
    pub seed: Option<u64>, // tests
}

pub struct PanelResult {
    pub selected: Vec<u64>,    // selecciones para el panel final
    pub eliminated: Vec<u64>,  // palabras descartadas en rojo
}

pub fn select_panel(components: &[Component], config: &PanelConfig) -> PanelResult {
    let mut rng: StdRng = match config.seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_entropy(),
    };

    // 1. Agrupar representantes por bin de GC
    let mut bin_members: Vec<Vec<&Component>> = vec![Vec::new(); config.bins.len()];
    for c in components {
        if let Some(idx) = config.bins.iter().position(|b| {
            c.representative_gc >= b.min_gc && c.representative_gc < b.max_gc
        }) {
            bin_members[idx].push(c);
        }
    }

    // 2. Cuota inicial por bin, según proporción configurada
    let mut quotas: Vec<usize> = config.bins.iter()
        .map(|b| (config.target_size as f32 * b.proportion).round() as usize)
        .collect();

    // 3. Redistribución de déficit, si se pidió
    if matches!(config.deficit, DeficitStrategy::Redistribute) {
        redistribute_deficit(&mut quotas, &bin_members);
    }

    // 4. Selección ponderada dentro de cada bin
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
    // Si no hay excedente suficiente en ningún otro bin, el déficit el panel queda más chico
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

    // Calculo de llave, importancia de selección
    // A-ExpJ: llave = u^(1/peso) con u ~ Uniforme(0,1). Tomar los `quota`
    // candidatos con mayor llave da muestreo ponderado sin reemplazo --
    // con Random (peso=1.0 fijo) esto se reduce exactamente a shuffle + take.
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
    use crate::clustering::Component;

    fn make_component(code: u64, gc: f32, entropy: f32, count: u32) -> Component {
        Component {
            representative: code,
            representative_count: count,
            representative_gc: gc,
            representative_entropy: entropy,
            size: 1,
            members: vec![code],
            files,
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
        // Bin bajo GC con solo 5 candidatos pero cuota de 25 (déficit=20)
        // Bin alto GC con 100 candidatos, cuota base 25, excedente=75
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
        // El bin bajo solo aporta sus 5 disponibles; el resto del target
        // (45) debe salir del bin alto vía redistribución.
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
        assert_eq!(result.selected.len(), 5); // no hay más de 5 disponibles
    }

    #[test]
    fn weighted_by_entropy_strongly_prefers_high_entropy() {
        let mut components = vec![make_component(0, 0.5, 0.01, 1)]; // entropía casi nula
        components.extend((1..50).map(|i| make_component(i, 0.5, 2.0, 1))); // entropía alta

        let config = PanelConfig {
            target_size: 5,
            bins: vec![GcBin { min_gc: 0.0, max_gc: 1.0, proportion: 1.0 }],
            strategy: SelectionStrategy::WeightedByEntropy,
            deficit: DeficitStrategy::Strict,
            seed: Some(1),
        };

        let result = select_panel(&components, &config);
        // Con peso tan desigual, es extremadamente improbable que el
        // candidato de entropía casi nula (code=0) sea elegido.
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