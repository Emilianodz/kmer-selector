use std::collections::HashSet;
use std::fs::File;
use std::io::{Write, BufWriter};
use std::path::Path;
use crate::clustering::Component;
use crate::io::decode_kmer;
use crate::selection::PanelResult;



pub fn export_components<P: AsRef<Path>>(
    components: &[Component],
    k: usize,
    path: P,
) -> Result<(), String> {
    let file = File::create(path).map_err(|e| e.to_string())?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "sequence\tcode\tcount\tgc_content\tentropy\tcomponent_size\tmembers")
        .map_err(|e| e.to_string())?;

    for c in components {
        let sequence = decode_kmer(c.representative, k);
        let members_str = c.members
            .iter()
            .map(|code| code.to_string())
            .collect::<Vec<_>>()
            .join(",");

        writeln!(
            writer,
            "{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{}",
            sequence,
            c.representative,
            c.representative_count,
            c.representative_gc,
            c.representative_entropy,
            c.size,
            members_str,
        ).map_err(|e| e.to_string())?;
    }

    Ok(())
}

pub fn export_panel<P: AsRef<Path>>(
    panel: &PanelResult,
    components: &[Component],
    k: usize,
    path: P,
) -> Result<(), String> {
    let file = File::create(path).map_err(|e| e.to_string())?;
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "sequence\tcode\tcount\tgc_content\tentropy\tcomponent_size\tselected"
    ).map_err(|e| e.to_string())?;

    let selected_set: HashSet<u64> = panel.selected.iter().copied().collect();

    for c in components {
        let sequence = decode_kmer(c.representative, k);
        let selected = selected_set.contains(&c.representative);

        writeln!(
            writer,
            "{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{}",
            sequence,
            c.representative,
            c.representative_count,
            c.representative_gc,
            c.representative_entropy,
            c.size,
            selected,
        ).map_err(|e| e.to_string())?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn exports_and_reads_back_correctly() {
        let components = vec![Component {
            representative: 0,
            representative_count: 5,
            representative_gc: 0.5,
            representative_entropy: 1.8,
            size: 3,
            members: vec![0, 1, 4],
            files: Vec::new(),
        }];

        let path = "test_export_output.tsv";
        export_components(&components, 4, path).unwrap();

        let content = fs::read_to_string(path).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 2); // header + 1 fila
        assert!(lines[1].starts_with("AAAA\t0\t5\t0.5000\t1.8000\t3\t0,1,4"));

        fs::remove_file(path).unwrap();
    }
}