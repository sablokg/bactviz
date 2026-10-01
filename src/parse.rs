//! GenBank / FASTA readers for (possibly multi-record) bacterial genomes.

use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;

/*
Gaurav Sablok
gsablok@proton.me
*/

#[derive(Clone, Debug)]
pub struct Feature {
    pub kind: String,
    pub gene: String,
    pub locus_tag: String,
    pub product: String,
    pub strand: i8,
    /// 0-based, half-open exon ranges
    pub ranges: Vec<(usize, usize)>,
}

impl Feature {
    pub fn label(&self) -> &str {
        if !self.gene.is_empty() {
            &self.gene
        } else if !self.locus_tag.is_empty() {
            &self.locus_tag
        } else {
            &self.kind
        }
    }
    fn key(&self) -> (usize, usize, i8) {
        let s = self.ranges.iter().map(|r| r.0).min().unwrap_or(0);
        let e = self.ranges.iter().map(|r| r.1).max().unwrap_or(0);
        (s, e, self.strand)
    }
}

pub struct Record {
    pub id: String,
    pub seq: Vec<u8>,
    pub features: Vec<Feature>,
}

/// One genome = one or more replicons/contigs (chromosome, plasmids, draft contigs...).
pub struct Genome {
    pub name: String,
    pub records: Vec<Record>,
}

impl Genome {
    pub fn total_len(&self) -> usize {
        self.records.iter().map(|r| r.seq.len()).sum()
    }
    pub fn n_features(&self) -> usize {
        self.records.iter().map(|r| r.features.len()).sum()
    }
}

pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            x => x,
        })
        .collect()
}

fn sanitize(id: &str) -> String {
    id.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || "._-".contains(c) {
                c
            } else {
                '_'
            }
        })
        .collect()
}

pub fn load(path: &Path, min_record: usize) -> Result<Genome, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let stem = path
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "genome".into());
    let first = text.lines().find(|l| !l.trim().is_empty()).unwrap_or("");
    let (mut records, def) = if first.starts_with('>') {
        (parse_fasta(&text), String::new())
    } else {
        parse_genbank(&text)
    };
    records.retain(|r| !r.seq.is_empty() && r.seq.len() >= min_record);
    if records.is_empty() {
        return Err(format!(
            "{}: no sequence records found (after --min-record filter)",
            path.display()
        ));
    }
    // MUMmer identifies sequences by header, so ids must be whitespace-free and unique
    let mut seen = HashSet::new();
    for r in &mut records {
        let base = sanitize(&r.id);
        let mut id = base.clone();
        let mut k = 2;
        while !seen.insert(id.clone()) {
            id = format!("{base}_{k}");
            k += 1;
        }
        r.id = id;
    }
    let name = if def.is_empty() { stem } else { def };
    Ok(Genome { name, records })
}

fn parse_fasta(text: &str) -> Vec<Record> {
    let mut recs: Vec<Record> = Vec::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            let id = h.split_whitespace().next().unwrap_or("seq").to_string();
            recs.push(Record {
                id,
                seq: Vec::new(),
                features: vec![],
            });
        } else if let Some(r) = recs.last_mut() {
            r.seq.extend(
                line.bytes()
                    .filter(|b| b.is_ascii_alphabetic())
                    .map(|b| b.to_ascii_uppercase()),
            );
        }
    }
    recs
}

// ---------------------------------------------------------------- GenBank

struct Pending {
    kind: String,
    loc: String,
    gene: String,
    tag: String,
    product: String,
    in_quals: bool,
    last: u8, // 3 = product (may continue on next line)
}

struct Builder {
    id: String,
    def: String,
    seq: Vec<u8>,
    features: Vec<Feature>,
    section: u8, // 0 header, 1 features, 2 origin
    cur: Option<Pending>,
}

impl Builder {
    fn flush(&mut self) {
        if let Some(p) = self.cur.take() {
            if matches!(
                p.kind.as_str(),
                "CDS" | "tRNA" | "rRNA" | "tmRNA" | "ncRNA" | "misc_RNA" | "gene"
            ) {
                let (strand, ranges) = parse_loc(&p.loc);
                if !ranges.is_empty() {
                    self.features.push(Feature {
                        kind: p.kind,
                        gene: p.gene,
                        locus_tag: p.tag,
                        product: p.product.trim_matches('"').to_string(),
                        strand,
                        ranges,
                    });
                }
            }
        }
    }

    fn line(&mut self, line: &str) {
        if line.starts_with("FEATURES") {
            self.section = 1;
            return;
        }
        if line.starts_with("ORIGIN") {
            self.flush();
            self.section = 2;
            return;
        }
        match self.section {
            0 => {
                if self.def.is_empty() {
                    if let Some(d) = line.strip_prefix("DEFINITION") {
                        // "Escherichia coli str. K-12 substr. MG1655, complete genome." -> name part
                        let d = d.trim().trim_end_matches('.');
                        self.def = d.split(", ").next().unwrap_or(d).trim().to_string();
                    }
                }
            }
            1 => {
                let new_feature =
                    line.len() > 21 && line.starts_with("     ") && line.as_bytes()[5] != b' ';
                if new_feature {
                    self.flush();
                    if let (Some(k), Some(l)) = (line.get(5..21), line.get(21..)) {
                        self.cur = Some(Pending {
                            kind: k.trim().to_string(),
                            loc: l.trim().to_string(),
                            gene: String::new(),
                            tag: String::new(),
                            product: String::new(),
                            in_quals: false,
                            last: 0,
                        });
                    }
                } else if let Some(p) = self.cur.as_mut() {
                    let t = line.trim();
                    if let Some(q) = t.strip_prefix('/') {
                        p.in_quals = true;
                        p.last = 0;
                        let unq = |v: &str| v.trim_matches('"').to_string();
                        if let Some(v) = q.strip_prefix("gene=") {
                            p.gene = unq(v);
                        } else if let Some(v) = q.strip_prefix("locus_tag=") {
                            p.tag = unq(v);
                        } else if let Some(v) = q.strip_prefix("product=") {
                            p.product = unq(v);
                            p.last = 3;
                        }
                    } else if !p.in_quals {
                        p.loc.push_str(t);
                    } else if p.last == 3 {
                        p.product.push(' ');
                        p.product.push_str(t);
                    }
                }
            }
            2 => self.seq.extend(
                line.bytes()
                    .filter(|b| b.is_ascii_alphabetic())
                    .map(|b| b.to_ascii_uppercase()),
            ),
            _ => {}
        }
    }

    fn finish(mut self) -> (Record, String) {
        self.flush();
        // Bacterial annotations carry a `gene` feature next to every CDS/RNA. Keep one drawn
        // feature per locus: gene-only features (pseudogenes, unannotated genes) are kept,
        // duplicates are dropped after copying the gene symbol onto the CDS/RNA feature.
        let (genes, mut out): (Vec<Feature>, Vec<Feature>) =
            self.features.into_iter().partition(|f| f.kind == "gene");
        let symbols: HashMap<_, _> = genes.iter().map(|g| (g.key(), g.gene.clone())).collect();
        let keys: HashSet<_> = out.iter().map(|f| f.key()).collect();
        for f in &mut out {
            if f.gene.is_empty() {
                if let Some(n) = symbols.get(&f.key()) {
                    f.gene = n.clone();
                }
            }
        }
        out.extend(genes.into_iter().filter(|g| !keys.contains(&g.key())));
        (
            Record {
                id: self.id,
                seq: self.seq,
                features: out,
            },
            self.def,
        )
    }
}

/// Handles complement(...), join(...), order(...) and `<`/`>` partial markers by pulling out
/// integer pairs; the strand comes from the presence of `complement`.
fn parse_loc(loc: &str) -> (i8, Vec<(usize, usize)>) {
    let strand = if loc.contains("complement") { -1 } else { 1 };
    let mut nums = Vec::new();
    let mut cur = String::new();
    for c in loc.chars().chain(std::iter::once(' ')) {
        if c.is_ascii_digit() {
            cur.push(c);
        } else if !cur.is_empty() {
            if let Ok(n) = cur.parse::<usize>() {
                nums.push(n);
            }
            cur.clear();
        }
    }
    let ranges = nums
        .chunks_exact(2)
        .filter(|c| c[0] >= 1 && c[1] >= c[0])
        .map(|c| (c[0] - 1, c[1]))
        .collect();
    (strand, ranges)
}

fn parse_genbank(text: &str) -> (Vec<Record>, String) {
    let mut recs = Vec::new();
    let mut name = String::new();
    let mut b: Option<Builder> = None;
    let done = |x: Builder, recs: &mut Vec<Record>, name: &mut String| {
        let (r, d) = x.finish();
        if name.is_empty() {
            *name = d;
        }
        recs.push(r);
    };
    for line in text.lines() {
        if line.starts_with("//") {
            if let Some(x) = b.take() {
                done(x, &mut recs, &mut name);
            }
            continue;
        }
        if line.starts_with("LOCUS") {
            if let Some(x) = b.take() {
                done(x, &mut recs, &mut name);
            }
            b = Some(Builder {
                id: line
                    .split_whitespace()
                    .nth(1)
                    .unwrap_or("record")
                    .to_string(),
                def: String::new(),
                seq: Vec::new(),
                features: Vec::new(),
                section: 0,
                cur: None,
            });
            continue;
        }
        if let Some(x) = b.as_mut() {
            x.line(line);
        }
    }
    if let Some(x) = b.take() {
        done(x, &mut recs, &mut name);
    }
    (recs, name)
}

// ---------------------------------------------------------------- rotation

/// Re-orient and rotate a circular replicon so that gene `gene` (e.g. dnaA) starts at position 1
/// on the forward strand. Returns false if the gene is not found.
pub fn rotate_to_gene(rec: &mut Record, gene: &str) -> bool {
    let g = gene.to_lowercase();
    // Not every GenBank record carries /gene= (e.g. locus_tag-only annotations), so fall back to
    // matching the locus tag, then to the gene name appearing inside the product description.
    let pick = rec
        .features
        .iter()
        .find(|f| f.kind == "CDS" && f.gene.eq_ignore_ascii_case(gene))
        .or_else(|| {
            rec.features
                .iter()
                .find(|f| f.gene.eq_ignore_ascii_case(gene))
        })
        .or_else(|| {
            rec.features
                .iter()
                .find(|f| f.kind == "CDS" && f.locus_tag.eq_ignore_ascii_case(gene))
        })
        .or_else(|| {
            rec.features.iter().find(|f| {
                f.kind == "CDS"
                    && f.product
                        .to_lowercase()
                        .split(|c: char| !c.is_alphanumeric())
                        .any(|w| w == g)
            })
        });
    let Some(f) = pick else { return false };
    let (s, e, strand) = f.key();
    let l = rec.seq.len();
    let flip = strand < 0;
    let shift = if flip { l - e } else { s };

    let mut seq = if flip {
        revcomp(&rec.seq)
    } else {
        rec.seq.clone()
    };
    seq.rotate_left(shift);
    rec.seq = seq;

    for f in &mut rec.features {
        let mut out = Vec::new();
        for &(a, b) in &f.ranges {
            let (a, b) = if flip { (l - b, l - a) } else { (a, b) };
            let a2 = (a + l - shift) % l;
            let b2 = a2 + (b - a);
            if b2 > l {
                out.push((a2, l));
                out.push((0, b2 - l));
            } else {
                out.push((a2, b2));
            }
        }
        if flip {
            f.strand = -f.strand;
        }
        f.ranges = out;
    }
    true
}
