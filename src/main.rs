mod mummer;
mod parse;
mod svg;

use std::fmt::Write as _;
use std::path::PathBuf;
use std::str::FromStr;
use std::{env, fs, process};

const USAGE: &str = "\
bactviz - linear genome alignment maps for bacterial genomes (MUMmer + gene annotations)

Gaurav Sablok
gsablok@proton.me

USAGE:
    bactviz [OPTIONS] <GENOME>...

Genomes (GenBank or FASTA; multi-record files = chromosome + plasmids/contigs) are drawn top to
bottom in the order given. Each genome is aligned to the one above with MUMmer (nucmer).
Requires MUMmer on PATH (conda install -c bioconda mummer | apt install mummer).

OUTPUT:
    -o, --output <FILE>       output SVG                              [chloroviz.svg]
        --tsv <FILE>          also write all blocks as a TSV table
        --title <TEXT>        figure title
        --names <A,B,...>     override genome labels (comma separated, one per genome)
        --width <PX>          figure width                            [1800]
        --align <center|left> track alignment                         [center]
        --vmin <PCT>          lower bound of the identity colour scale [auto]

GENES:
        --no-genes            do not draw gene arrows
        --labels              label every gene that has a symbol (as many as fit). + strand genes
                              are named above the track, - strand genes below, with a leader line
        --label-genes <A,B>   always label these gene symbols (e.g. dnaA,gyrB,rpoB)
        --label-size <PX>     label font size                         [8]
                              Tip: gene names only fit when zoomed in - use a big --width
                              (SVG is vector; e.g. --width 12000) or --label-genes.

BACTERIAL OPTIONS:
        --start-gene <NAME>   re-orient and rotate the largest replicon of each genome so this gene
                              (usually dnaA) starts at position 1 on the + strand
        --min-record <BP>     ignore records/contigs shorter than this [0]

MUMMER OPTIONS:
        --minmatch <N>        nucmer -l  minimum exact match          [20]
        --mincluster <N>      nucmer -c  minimum cluster length       [65]
        --match-mode <M>      mumreference | mum | maxmatch           [mumreference]
        --one-to-one          run delta-filter -1 (one-to-one alignments)
        --min-length <BP>    drop alignments shorter than this        [1000]
        --min-identity <PCT>  drop alignments below this identity     [0]
        --mummer-dir <DIR>    directory containing nucmer/show-coords/delta-filter
        --keep-temp           keep intermediate MUMmer files

    -h, --help                show this help
";

fn die(msg: &str) -> ! {
    eprintln!("error: {msg}");
    process::exit(2);
}

fn val(it: &mut impl Iterator<Item = String>, flag: &str) -> String {
    it.next()
        .unwrap_or_else(|| die(&format!("{flag} needs a value")))
}

fn num<T: FromStr>(s: String, flag: &str) -> T {
    s.parse()
        .unwrap_or_else(|_| die(&format!("invalid value '{s}' for {flag}")))
}

fn list(s: String) -> Vec<String> {
    s.split(',')
        .map(|x| x.trim().to_string())
        .filter(|x| !x.is_empty())
        .collect()
}

fn main() {
    let mut files: Vec<PathBuf> = Vec::new();
    let mut out = String::from("chloroviz.svg");
    let mut tsv: Option<String> = None;
    let mut names: Vec<String> = Vec::new();
    let mut start_gene: Option<String> = None;
    let mut min_record = 0usize;
    let mut keep_temp = false;
    let mut cfg = mummer::Config {
        dir: None,
        mode: "mumreference".into(),
        minmatch: 20,
        mincluster: 65,
        one_to_one: false,
        min_len: 1000,
        min_ident: 0.0,
    };
    let mut o = svg::Opts {
        width: 1800.0,
        center: true,
        genes: true,
        labels: false,
        label_genes: Vec::new(),
        label_size: 8.0,
        title: String::new(),
        vmin: None,
    };

    let mut it = env::args().skip(1);
    while let Some(a) = it.next() {
        match a.as_str() {
            "-h" | "--help" => {
                print!("{USAGE}");
                return;
            }
            "-o" | "--output" => out = val(&mut it, &a),
            "--tsv" => tsv = Some(val(&mut it, &a)),
            "--title" => o.title = val(&mut it, &a),
            "--names" => names = list(val(&mut it, &a)),
            "--width" => o.width = num(val(&mut it, &a), &a),
            "--align" => match val(&mut it, &a).as_str() {
                "center" => o.center = true,
                "left" => o.center = false,
                v => die(&format!("--align must be center or left, got '{v}'")),
            },
            "--vmin" => o.vmin = Some(num(val(&mut it, &a), &a)),
            "--no-genes" => o.genes = false,
            "--labels" => o.labels = true,
            "--label-genes" => {
                o.label_genes = list(val(&mut it, &a))
                    .into_iter()
                    .map(|s| s.to_lowercase())
                    .collect()
            }
            "--label-size" => o.label_size = num(val(&mut it, &a), &a),
            "--start-gene" => start_gene = Some(val(&mut it, &a)),
            "--min-record" => min_record = num(val(&mut it, &a), &a),
            "--minmatch" => cfg.minmatch = num(val(&mut it, &a), &a),
            "--mincluster" => cfg.mincluster = num(val(&mut it, &a), &a),
            "--match-mode" => {
                cfg.mode = val(&mut it, &a);
                if !["mumreference", "mum", "maxmatch"].contains(&cfg.mode.as_str()) {
                    die("--match-mode must be mumreference, mum or maxmatch");
                }
            }
            "--one-to-one" => cfg.one_to_one = true,
            "--min-length" => cfg.min_len = num(val(&mut it, &a), &a),
            "--min-identity" => cfg.min_ident = num(val(&mut it, &a), &a),
            "--mummer-dir" => cfg.dir = Some(PathBuf::from(val(&mut it, &a))),
            "--keep-temp" => keep_temp = true,
            s if s.starts_with('-') && s.len() > 1 => {
                eprintln!("unknown option '{s}'\n\n{USAGE}");
                process::exit(2);
            }
            _ => files.push(PathBuf::from(&a)),
        }
    }
    if files.len() < 2 {
        eprint!("{USAGE}");
        die("need at least two genome files");
    }
    if !names.is_empty() && names.len() != files.len() {
        die("--names needs exactly one name per genome file");
    }

    let mut genomes: Vec<parse::Genome> = files
        .iter()
        .map(|f| parse::load(f, min_record).unwrap_or_else(|e| die(&e)))
        .collect();
    for (i, g) in genomes.iter_mut().enumerate() {
        if let Some(n) = names.get(i) {
            g.name = n.clone();
        }
        if let Some(gene) = &start_gene {
            // the replicon of interest is the chromosome = the largest record
            let big = (0..g.records.len())
                .max_by_key(|&k| g.records[k].seq.len())
                .unwrap();
            if parse::rotate_to_gene(&mut g.records[big], gene) {
                eprintln!(
                    "  {}: rotated {} to start at {gene}",
                    g.name, g.records[big].id
                );
            } else {
                eprintln!(
                    "  warning: {}: gene '{gene}' not found in {}; left unrotated",
                    g.name, g.records[big].id
                );
            }
        }
        eprintln!(
            "loaded {:<42} {:>10} bp  {:>3} record(s)  {:>5} genes",
            g.name,
            g.total_len(),
            g.records.len(),
            g.n_features()
        );
    }

    // ---- MUMmer
    let tmp = env::temp_dir().join(format!("chloroviz_{}", process::id()));
    fs::create_dir_all(&tmp).unwrap_or_else(|e| die(&format!("{}: {e}", tmp.display())));
    let fastas: Vec<PathBuf> = genomes
        .iter()
        .enumerate()
        .map(|(i, g)| {
            let p = tmp.join(format!("g{i}.fa"));
            mummer::write_fasta(&p, g).unwrap_or_else(|e| die(&format!("{}: {e}", p.display())));
            p
        })
        .collect();

    let mut links = Vec::new();
    let mut table = String::from(
        "ref_genome\tref_record\tref_start\tref_end\tquery_genome\tquery_record\tquery_start\tquery_end\tstrand\tidentity\n",
    );
    for i in 0..genomes.len() - 1 {
        let (r, q) = (&genomes[i], &genomes[i + 1]);
        let blocks = mummer::run(&cfg, &tmp, i, (r, &fastas[i]), (q, &fastas[i + 1]))
            .unwrap_or_else(|e| {
                if !keep_temp {
                    let _ = fs::remove_dir_all(&tmp);
                }
                die(&e)
            });
        eprintln!(
            "aligned {} -> {}: {} blocks ({} inverted)",
            q.name,
            r.name,
            blocks.len(),
            blocks.iter().filter(|b| b.inverted).count()
        );
        for b in &blocks {
            let _ = writeln!(
                table,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}",
                r.name,
                r.records[b.r_rec].id,
                b.r_start + 1,
                b.r_end,
                q.name,
                q.records[b.q_rec].id,
                b.q_start + 1,
                b.q_end,
                if b.inverted { '-' } else { '+' },
                b.identity
            );
        }
        links.push(blocks);
    }
    if keep_temp {
        eprintln!("kept MUMmer files in {}", tmp.display());
    } else {
        let _ = fs::remove_dir_all(&tmp);
    }

    fs::write(&out, svg::render(&genomes, &links, &o))
        .unwrap_or_else(|e| die(&format!("{out}: {e}")));
    eprintln!("wrote {out}");
    if let Some(t) = tsv {
        fs::write(&t, table).unwrap_or_else(|e| die(&format!("{t}: {e}")));
        eprintln!("wrote {t}");
    }
}
