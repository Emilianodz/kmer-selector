use std::collections::HashSet;
use std::io::{BufRead, BufReader};
use std::fs::File;
use std::path::Path;

pub const K: usize = 13;

pub struct ProbeCandidate {
    pub code: u64,
    pub count: u32,
    pub gc_content: f32,
    pub entropy: f32,
    pub files: HashSet<u32>,
}


pub fn shannon_entropy(code: u64, k: usize) -> f32 {
    let mut counts = [0u32; 4];
    let mut remaining = code;
    for _ in 0..k {
        let base = (remaining & 0b11) as usize;
        counts[base] += 1;
        remaining >>= 2;
    }
    let k_f = k as f32;
    let mut entropy = 0.0f32;
    for &count in counts.iter() {
        if count > 0 {
            let p = count as f32 / k_f;
            entropy -= p * p.log2();
        }
    }
    entropy
}

pub fn gc_content(code: u64, k: usize) -> f32 {
    let xor = code ^ (code >> 1);
    const EVEN_BITS_MASK: u64 = 0x5555_5555_5555_5555;
    let gc_bits = xor & EVEN_BITS_MASK;
    let gc_count = gc_bits.count_ones();
    gc_count as f32 / k as f32
}

pub fn read_candidates<P: AsRef<Path>>(path: P, k: usize) -> Result<Vec<ProbeCandidate>, String> {
    let file = File::open(path).map_err(|e| e.to_string())?;
    let reader = BufReader::new(file);
    let mut candidates = Vec::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| e.to_string())?;
        let candidate = parse_line(&line, k)
            .map_err(|e| format!("línea {}: {}", line_num + 1, e))?;
        candidates.push(candidate);
    }

    Ok(candidates)
}

fn parse_line(line: &str, k: usize) -> Result<ProbeCandidate, String> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() != 6 {
        return Err(format!("se esperaban 6 columnas, se encontraron {}", fields.len()));
    }

    let code = fields[1].parse::<u64>().map_err(|e| e.to_string())?;
    let count = fields[2].parse::<u32>().map_err(|e| e.to_string())?;

    let gc = parse_or_compute(fields[3], || gc_content(code, k))?;
    let entropy = parse_or_compute(fields[4], || shannon_entropy(code, k))?;

    let files = fields[5]
        .split(',')
        .map(|id| id.parse::<u32>().map_err(|e| e.to_string()))
        .collect::<Result<HashSet<u32>, String>>()?;

    Ok(ProbeCandidate { code, count, gc_content: gc, entropy, files })
}

fn parse_or_compute(field: &str, compute: impl FnOnce() -> f32) -> Result<f32, String> {
    if field == "NA" {
        Ok(compute())
    } else {
        field.parse::<f32>().map_err(|e| e.to_string())
    }
}


pub fn decode_kmer(code: u64, k: usize) -> String {
    let mut bases = Vec::with_capacity(k);

    for i in 0..k {
        let shift = (k - 1 - i) * 2;
        let bits = (code >> shift) & 0b11;
        let base = match bits {
            0b00 => b'A',
            0b01 => b'C',
            0b10 => b'G',
            0b11 => b'T',
            _ => unreachable!("una máscara de 2 bits solo puede dar 0..=3"),
        };
        bases.push(base);
    }

    String::from_utf8(bases).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shannon_entropy_max_for_uniform_distribution() {
        // k=4 con las 4 bases distintas → entropía máxima = log2(4) = 2.0
        // code = A C G T = 00 01 10 11 → 0b11100100 (orden depende de tu encode_kmer)
        let code: u64 = 0b11_10_01_00;
        let e = shannon_entropy(code, 4);
        assert!((e - 2.0).abs() < 1e-5);
    }

    #[test]
    fn shannon_entropy_zero_for_homogeneous_kmer() {
        // todo A (00) → entropía 0
        let code: u64 = 0;
        let e = shannon_entropy(code, 13);
        assert!(e.abs() < 1e-5);
    }

    #[test]
    fn gc_content_all_at() {
        // A=00, T=11 → ninguna es GC
        let code: u64 = 0b11_00_11_00; // A T A T
        let gc = gc_content(code, 4);
        assert_eq!(gc, 0.0);
    }

    #[test]
    fn gc_content_all_gc() {
        // C=01, G=10 → ambas son GC
        let code: u64 = 0b10_01_10_01; // C G C G
        let gc = gc_content(code, 4);
        assert_eq!(gc, 1.0);
    }

    #[test]
    fn decode_kmer_roundtrip_with_gc_content() {
        // ATG codificado y decodificado debe regresar "ATG"
        let code = 0b00_11_10u64; // A=00, T=11, G=10 en orden
        assert_eq!(decode_kmer(code, 3), "ATG");
    }

    #[test]
    fn parses_line_with_all_fields() {
        let line = "AAACCCGGGTTTA\t123456\t42\t0.5\t1.8\t1,2,3";
        let c = parse_line(line, K).unwrap();
        assert_eq!(c.code, 123456);
        assert_eq!(c.count, 42);
        assert_eq!(c.gc_content, 0.5);
        assert_eq!(c.entropy, 1.8);
        assert_eq!(c.files, HashSet::from([1, 2, 3]));
    }

    #[test]
    fn computes_gc_and_entropy_when_na() {
        let code = 123456u64;
        let line = format!("AAACCCGGGTTTA\t{}\t42\tNA\tNA\t1", code);
        let c = parse_line(&line, K).unwrap();
        assert_eq!(c.gc_content, gc_content(code, K));
        assert_eq!(c.entropy, shannon_entropy(code, K));
    }

    #[test]
    fn rejects_wrong_column_count() {
        let line = "AAACCCGGGTTTA\t123456\t42";
        assert!(parse_line(line, K).is_err());
    }
}