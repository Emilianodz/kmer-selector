use crate::io::ProbeCandidate;

pub struct FilterConfig {
    pub min_gc: Option<f32>,
    pub max_gc: Option<f32>,
    pub min_entropy: Option<f32>,
}

impl Default for FilterConfig {
    fn default() -> Self {
        FilterConfig { min_gc: None, max_gc: None, min_entropy: None }
    }
}

pub fn apply_filter(candidates: Vec<ProbeCandidate>, config: &FilterConfig) -> Vec<ProbeCandidate> {
    candidates
        .into_iter()
        .filter(|c| passes(c, config))
        .collect()
}

fn passes(c: &ProbeCandidate, config: &FilterConfig) -> bool {
    if let Some(min_gc) = config.min_gc {
        if c.gc_content < min_gc { return false; }
    }
    if let Some(max_gc) = config.max_gc {
        if c.gc_content > max_gc { return false; }
    }
    if let Some(min_entropy) = config.min_entropy {
        if c.entropy < min_entropy { return false; }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    fn make_candidate(gc: f32, entropy: f32) -> ProbeCandidate {
        ProbeCandidate { code: 0, count: 1, gc_content: gc, entropy, files: HashSet::new() }
    }

    #[test]
    fn default_config_passes_everything() {
        let candidates = vec![make_candidate(0.9, 0.1), make_candidate(0.1, 3.0)];
        let filtered = apply_filter(candidates, &FilterConfig::default());
        assert_eq!(filtered.len(), 2);
    }

    #[test]
    fn filters_by_gc_range() {
        let candidates = vec![make_candidate(0.5, 1.0), make_candidate(0.9, 1.0)];
        let config = FilterConfig { min_gc: Some(0.3), max_gc: Some(0.7), min_entropy: None };
        let filtered = apply_filter(candidates, &config);
        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].gc_content, 0.5);
    }
}