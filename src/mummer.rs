//! Collinear-block detection by driving MUMmer, the same pipeline pyGenomeViz uses:
//!
//!   nucmer [--mumreference|--mum|--maxmatch] -l <minmatch> -c <mincluster> -p <prefix> ref.fa qry.fa
//!   [delta-filter -1 prefix.delta]                       (optional one-to-one filtering)
//!   show-coords -T -H -l prefix.delta                    (tab-separated alignment coordinates)
//!
//! Multi-record genomes (chromosome + plasmids, draft contigs) are handled natively: every
//! alignment line carries the reference and query record ids.

use crate::parse::Genome;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/*
Gaurav Sablok
gsablok@proton.me
*/

#[derive(Clone, Debug)]
pub struct Block {
    pub r_rec: usize,
    pub r_start: usize, // 0-based, half-open
    pub r_end: usize,
    pub q_rec: usize,
    pub q_start: usize,
    pub q_end: usize,
    pub inverted: bool,
    pub identity: f64,
}

pub struct Config {
    pub dir: Option<PathBuf>,
    pub mode: String,
    pub minmatch: usize,
    pub mincluster: usize,
    pub one_to_one: bool,
    pub min_len: usize,
    pub min_ident: f64,
}

fn bin(cfg: &Config, name: &str) -> PathBuf {
    match &cfg.dir {
        Some(d) => d.join(name),
        None => PathBuf::from(name),
    }
}

fn exec(cmd: &mut Command, name: &str) -> Result<Output, String> {
    match cmd.output() {
        Ok(o) if o.status.success() => Ok(o),
        Ok(o) => Err(format!("{name} failed:\n{}", String::from_utf8_lossy(&o.stderr))),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Err(format!(
            "'{name}' not found. Install MUMmer (conda install -c bioconda mummer, or apt install mummer) \
             or point --mummer-dir at its bin directory"
        )),
        Err(e) => Err(format!("{name}: {e}")),
    }
}

pub fn write_fasta(path: &Path, g: &Genome) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    for r in &g.records {
        writeln!(w, ">{}", r.id)?;
        for chunk in r.seq.chunks(80) {
            w.write_all(chunk)?;
            w.write_all(b"\n")?;
        }
    }
    w.flush()
}

/// Align `q` against `r` (FASTA files already written) and return the collinear blocks.
pub fn run(
    cfg: &Config,
    tmp: &Path,
    idx: usize,
    r: (&Genome, &Path),
    q: (&Genome, &Path),
) -> Result<Vec<Block>, String> {
    let prefix = tmp.join(format!("pair{idx}"));
    let mut nucmer = Command::new(bin(cfg, "nucmer"));
    nucmer
        .arg(format!("--{}", cfg.mode))
        .args([
            "-l",
            &cfg.minmatch.to_string(),
            "-c",
            &cfg.mincluster.to_string(),
        ])
        .arg("-p")
        .arg(&prefix)
        .arg(r.1)
        .arg(q.1);
    exec(&mut nucmer, "nucmer")?;

    let mut delta = tmp.join(format!("pair{idx}.delta"));
    if cfg.one_to_one {
        let o = exec(
            Command::new(bin(cfg, "delta-filter")).arg("-1").arg(&delta),
            "delta-filter",
        )?;
        let filt = tmp.join(format!("pair{idx}.filt.delta"));
        fs::write(&filt, o.stdout).map_err(|e| format!("{}: {e}", filt.display()))?;
        delta = filt;
    }
    let out = exec(
        Command::new(bin(cfg, "show-coords"))
            .args(["-T", "-H", "-l"])
            .arg(&delta),
        "show-coords",
    )?;
    let text = String::from_utf8_lossy(&out.stdout);

    let rmap: HashMap<&str, usize> =
        r.0.records
            .iter()
            .enumerate()
            .map(|(i, x)| (x.id.as_str(), i))
            .collect();
    let qmap: HashMap<&str, usize> =
        q.0.records
            .iter()
            .enumerate()
            .map(|(i, x)| (x.id.as_str(), i))
            .collect();

    // columns: S1 E1 S2 E2 LEN1 LEN2 %IDY LENR LENQ [COVR COVQ] TAGR TAGQ
    let mut blocks = Vec::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 11 {
            continue;
        }
        let num = |i: usize| f[i].parse::<f64>().ok();
        let (Some(s1), Some(e1), Some(s2), Some(e2), Some(l1), Some(l2), Some(idy)) =
            (num(0), num(1), num(2), num(3), num(4), num(5), num(6))
        else {
            continue;
        };
        let (Some(&ri), Some(&qi)) = (rmap.get(f[f.len() - 2]), qmap.get(f[f.len() - 1])) else {
            continue;
        };
        if (l1.max(l2) as usize) < cfg.min_len || idy < cfg.min_ident {
            continue;
        }
        blocks.push(Block {
            r_rec: ri,
            r_start: s1.min(e1) as usize - 1,
            r_end: s1.max(e1) as usize,
            q_rec: qi,
            q_start: s2.min(e2) as usize - 1,
            q_end: s2.max(e2) as usize,
            inverted: s2 > e2,
            identity: idy,
        });
    }
    Ok(blocks)
}
