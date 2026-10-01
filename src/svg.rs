//! SVG rendering: stacked linear genome tracks (multi-record), gene arrows for every
//! annotated gene, gene-name labels, and MUMmer collinear-block links between tracks.

use crate::mummer::Block;
use crate::parse::Genome;
use std::collections::HashSet;
use std::fmt::Write;

/*
Gaurav Sablok
gsablok@proton.me
*/

pub struct Opts {
    pub width: f64,
    pub center: bool,
    pub genes: bool,
    /// label every gene that has a gene symbol (as many as fit without overlapping)
    pub labels: bool,
    /// gene symbols that are always labelled (lower-cased)
    pub label_genes: Vec<String>,
    pub label_size: f64,
    pub title: String,
    pub vmin: Option<f64>,
}

const GREY0: [f64; 3] = [205.0, 205.0, 205.0];
const GREY1: [f64; 3] = [70.0, 70.0, 70.0];
const RED0: [f64; 3] = [250.0, 190.0, 190.0];
const RED1: [f64; 3] = [190.0, 20.0, 20.0];

fn esc(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn trunc(s: &str, n: usize) -> String {
    if s.chars().count() <= n {
        s.to_string()
    } else {
        format!(
            "{}…",
            s.chars().take(n.saturating_sub(1)).collect::<String>()
        )
    }
}

fn lerp(a: [f64; 3], b: [f64; 3], t: f64) -> String {
    let c = |i: usize| (a[i] + (b[i] - a[i]) * t).round() as u8;
    format!("#{:02x}{:02x}{:02x}", c(0), c(1), c(2))
}

fn fmt_len(n: usize) -> String {
    if n >= 1_000_000 {
        format!("{:.2} Mb", n as f64 / 1e6)
    } else if n >= 10_000 {
        format!("{:.1} kb", n as f64 / 1e3)
    } else {
        format!("{n} bp")
    }
}

fn feature_color(kind: &str) -> &'static str {
    match kind {
        "CDS" => "#f0a03c",
        "tRNA" => "#3aa17e",
        "rRNA" => "#8e6bbf",
        "tmRNA" | "ncRNA" | "misc_RNA" => "#d6609c",
        _ => "#9a9a9a", // gene-only features (pseudogenes, unannotated genes)
    }
}

/// Text with a white halo, drawn as a separate underlay so it also works in renderers that
/// ignore `paint-order` (Inkscape/cairosvg/PDF converters).
fn halo(pos: &str, extra: &str, fill: &str, txt: &str) -> String {
    format!(
        r##"<text {pos} {extra} fill="white" stroke="white" stroke-width="2.8" stroke-linejoin="round">{txt}</text><text {pos} {extra} fill="{fill}">{txt}</text>"##
    )
}

struct Label {
    x: f64,
    text: String,
    above: bool,
    prio: u8,
    forced: bool,
}

/// `links[i]` holds the blocks between genome `i` (reference, upper) and `i + 1` (query, lower).
pub fn render(gs: &[Genome], links: &[Vec<Block>], o: &Opts) -> String {
    let n = gs.len();
    let (ml, mr) = (260.0, 30.0);
    let (th, ls) = (26.0, o.label_size);
    let forced: HashSet<&str> = o.label_genes.iter().map(|s| s.as_str()).collect();
    let lab_on = o.labels || !forced.is_empty();
    let lab_h = if lab_on { 10.0 * ls * 0.62 + 10.0 } else { 0.0 };
    let top = (if o.title.is_empty() { 30.0 } else { 60.0 }) + lab_h;
    let gap_y = 70.0 + 2.0 * lab_h;
    let plot_w = o.width - ml - mr;

    // horizontal layout: records are laid side by side with a small fixed pixel gap
    let maxlen = gs.iter().map(|g| g.total_len()).max().unwrap_or(1).max(1) as f64;
    let maxgaps = gs.iter().map(|g| g.records.len() - 1).max().unwrap_or(0);
    let gap = if maxgaps == 0 {
        0.0
    } else {
        (plot_w * 0.06 / maxgaps as f64).min(6.0)
    };
    let sc = (plot_w - gap * maxgaps as f64) / maxlen;
    let mut rec_x: Vec<Vec<f64>> = Vec::new();
    let mut width: Vec<f64> = Vec::new();
    for g in gs {
        let (mut xs, mut cur) = (Vec::new(), 0.0);
        for r in &g.records {
            xs.push(cur);
            cur += r.seq.len() as f64 * sc + gap;
        }
        width.push(cur - gap);
        rec_x.push(xs);
    }
    let off = |i: usize| {
        if o.center {
            (plot_w - width[i]) / 2.0
        } else {
            0.0
        }
    };
    let x = |i: usize, rec: usize, p: usize| ml + off(i) + rec_x[i][rec] + p as f64 * sc;
    let y = |i: usize| top + i as f64 * (th + gap_y);

    let vmin = o.vmin.unwrap_or_else(|| {
        links
            .iter()
            .flatten()
            .map(|b| b.identity)
            .fold(100.0f64, f64::min)
            .floor()
            .min(99.0)
    });
    let sb_y = y(n - 1) + th + lab_h + 35.0;
    let height = sb_y + 120.0;

    let mut s = String::new();
    macro_rules! w {
        ($($t:tt)*) => { writeln!(s, $($t)*).unwrap() };
    }

    w!(r#"<?xml version="1.0" encoding="UTF-8"?>"#);
    w!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" viewBox="0 0 {} {}" font-family="Helvetica, Arial, sans-serif">"#,
        o.width,
        height,
        o.width,
        height
    );
    w!(
        r##"<defs><linearGradient id="gf"><stop offset="0" stop-color="#cdcdcd"/><stop offset="1" stop-color="#464646"/></linearGradient><linearGradient id="gi"><stop offset="0" stop-color="#fabebe"/><stop offset="1" stop-color="#be1414"/></linearGradient></defs>"##
    );
    w!(r#"<rect width="100%" height="100%" fill="white"/>"#);
    if !o.title.is_empty() {
        w!(
            r#"<text x="{}" y="32" font-size="18" font-weight="bold" text-anchor="middle">{}</text>"#,
            o.width / 2.0,
            esc(&o.title)
        );
    }

    // ---- links first, so tracks and labels paint on top
    for (i, blocks) in links.iter().enumerate() {
        let (y1, y2) = (y(i) + th, y(i + 1));
        for b in blocks {
            let t = ((b.identity - vmin) / (100.0 - vmin)).clamp(0.0, 1.0);
            let col = if b.inverted {
                lerp(RED0, RED1, t)
            } else {
                lerp(GREY0, GREY1, t)
            };
            let widen = |a: f64, b: f64| {
                if b - a < 1.4 {
                    let c = (a + b) / 2.0;
                    (c - 0.7, c + 0.7)
                } else {
                    (a, b)
                }
            };
            let (r0, r1) = widen(x(i, b.r_rec, b.r_start), x(i, b.r_rec, b.r_end));
            let (q0, q1) = widen(x(i + 1, b.q_rec, b.q_start), x(i + 1, b.q_rec, b.q_end));
            let (qa, qb) = if b.inverted { (q1, q0) } else { (q0, q1) }; // twist = inversion
            w!(
                r#"<polygon points="{:.2},{:.2} {:.2},{:.2} {:.2},{:.2} {:.2},{:.2}" fill="{}" fill-opacity="0.72" stroke="none"><title>{}:{}-{} vs {}:{}-{} ({}, {:.1}% identity)</title></polygon>"#,
                r0,
                y1,
                r1,
                y1,
                qb,
                y2,
                qa,
                y2,
                col,
                gs[i].records[b.r_rec].id,
                b.r_start + 1,
                b.r_end,
                gs[i + 1].records[b.q_rec].id,
                b.q_start + 1,
                b.q_end,
                if b.inverted { "inverted" } else { "forward" },
                b.identity
            );
        }
    }

    // soft white wash behind the label zones so names stay legible on top of the link ribbons
    if lab_on {
        for i in 0..n {
            let yt = y(i);
            if i > 0 {
                w!(
                    r##"<rect x="{}" y="{:.1}" width="{}" height="{:.1}" fill="white" fill-opacity="0.66"/>"##,
                    ml,
                    yt - lab_h,
                    plot_w,
                    lab_h
                );
            }
            if i + 1 < n {
                w!(
                    r##"<rect x="{}" y="{:.1}" width="{}" height="{:.1}" fill="white" fill-opacity="0.66"/>"##,
                    ml,
                    yt + th,
                    plot_w,
                    lab_h
                );
            }
        }
    }

    // ---- tracks
    let mut labels_out = String::new();
    for (i, g) in gs.iter().enumerate() {
        let yt = y(i);
        let mut cands: Vec<Label> = Vec::new();

        for (ri, rec) in g.records.iter().enumerate() {
            w!(
                r##"<rect x="{:.2}" y="{:.2}" width="{:.2}" height="{}" fill="#f2f2f2" stroke="#999" stroke-width="0.6"><title>{} ({})</title></rect>"##,
                x(i, ri, 0),
                yt,
                rec.seq.len() as f64 * sc,
                th,
                esc(&rec.id),
                fmt_len(rec.seq.len())
            );
            // record ids under the track when there is room (labels off) and only a few records
            if !lab_on && g.records.len() > 1 && g.records.len() <= 12 {
                let wpx = rec.seq.len() as f64 * sc;
                if wpx >= 36.0 {
                    let cx = x(i, ri, 0) + wpx / 2.0;
                    let pos = format!(
                        r#"x="{:.1}" y="{:.1}" text-anchor="middle""#,
                        cx,
                        yt + th + 10.0
                    );
                    let t = esc(&trunc(&rec.id, ((wpx / 4.6) as usize).max(4)));
                    w!("{}", halo(&pos, r#"font-size="8""#, "#444", &t));
                }
            }
            if !o.genes {
                continue;
            }
            for f in &rec.features {
                let col = feature_color(&f.kind);
                let fwd = f.strand >= 0;
                let (ya, yb) = if fwd {
                    (yt + 2.0, yt + th / 2.0 - 0.5)
                } else {
                    (yt + th / 2.0 + 0.5, yt + th - 2.0)
                };
                let ym = (ya + yb) / 2.0;
                let tip = [f.label(), &f.locus_tag, &f.product, &f.kind]
                    .iter()
                    .filter(|t| !t.is_empty())
                    .map(|t| esc(t))
                    .collect::<Vec<_>>()
                    .join(" | ");
                let (mut fmin, mut fmax) = (usize::MAX, 0usize);
                for &(a, b) in &f.ranges {
                    let b = b.min(rec.seq.len());
                    if a >= b {
                        continue;
                    }
                    fmin = fmin.min(a);
                    fmax = fmax.max(b);
                    let (x0, mut x1) = (x(i, ri, a), x(i, ri, b));
                    if x1 - x0 < 0.6 {
                        x1 = x0 + 0.6;
                    }
                    let h = ((x1 - x0) * 0.5).min(4.0);
                    let pts = if fwd {
                        format!("{x0:.1},{ya:.1} {:.1},{ya:.1} {x1:.1},{ym:.1} {:.1},{yb:.1} {x0:.1},{yb:.1}", x1 - h, x1 - h)
                    } else {
                        format!("{x1:.1},{ya:.1} {:.1},{ya:.1} {x0:.1},{ym:.1} {:.1},{yb:.1} {x1:.1},{yb:.1}", x0 + h, x0 + h)
                    };
                    w!(
                        r##"<polygon points="{}" fill="{}" stroke="#555" stroke-width="0.25"><title>{} {}:{}-{} ({})</title></polygon>"##,
                        pts,
                        col,
                        tip,
                        esc(&rec.id),
                        a + 1,
                        b,
                        if fwd { "+" } else { "-" }
                    );
                }
                // Prefer the /gene= symbol; fall back to locus_tag when the record carries no
                // gene symbols at all (common in locus_tag-only NCBI/INSDC annotations).
                let name = if !f.gene.is_empty() {
                    f.gene.as_str()
                } else {
                    f.locus_tag.as_str()
                };
                if fmin < fmax && lab_on && !name.is_empty() {
                    let lname = name.to_lowercase();
                    let is_forced = forced.contains(lname.as_str())
                        || forced.iter().any(|w| {
                            f.product
                                .to_lowercase()
                                .split(|c: char| !c.is_alphanumeric())
                                .any(|t| t == *w)
                        });
                    if o.labels || is_forced {
                        let cx = (x(i, ri, fmin) + x(i, ri, fmax)) / 2.0;
                        cands.push(Label {
                            x: cx,
                            text: trunc(name, 10),
                            above: fwd,
                            prio: if is_forced { 0 } else { 1 },
                            forced: is_forced,
                        });
                    }
                }
            }
        }

        // ---- gene-name labels: + strand above the band, - strand below it. Greedy packing:
        // forced names first, then the rest left to right, skipping anything that would overlap.
        cands.sort_by(|a, b| (a.prio, a.x).partial_cmp(&(b.prio, b.x)).unwrap());
        let spacing = ls * 1.15;
        let (mut up, mut down): (Vec<f64>, Vec<f64>) = (Vec::new(), Vec::new());
        for c in &cands {
            let placed = if c.above { &mut up } else { &mut down };
            if !c.forced && placed.iter().any(|&p| (p - c.x).abs() < spacing) {
                continue;
            }
            placed.push(c.x);
            let extra = format!(r#"font-size="{ls}" font-style="italic""#);
            let t = esc(&c.text);
            if c.above {
                let pos = format!(
                    r#"transform="translate({:.1},{:.1}) rotate(-90)""#,
                    c.x + ls * 0.35,
                    yt - 8.0
                );
                writeln!(labels_out, r##"<line x1="{x:.1}" y1="{a:.1}" x2="{x:.1}" y2="{b:.1}" stroke="#555" stroke-width="0.6"/>{}"##,
                    halo(&pos, &extra, "#111", &t), x = c.x, a = yt, b = yt - 7.0).unwrap();
            } else {
                let pos = format!(
                    r#"transform="translate({:.1},{:.1}) rotate(90)""#,
                    c.x - ls * 0.35,
                    yt + th + 8.0
                );
                writeln!(labels_out, r##"<line x1="{x:.1}" y1="{a:.1}" x2="{x:.1}" y2="{b:.1}" stroke="#555" stroke-width="0.6"/>{}"##,
                    halo(&pos, &extra, "#111", &t), x = c.x, a = yt + th, b = yt + th + 7.0).unwrap();
            }
        }

        w!(
            r#"<text x="{}" y="{:.2}" font-size="12" text-anchor="end"><title>{}</title>{}</text>"#,
            ml - 10.0,
            yt + th / 2.0 - 1.0,
            esc(&g.name),
            esc(&trunc(&g.name, 38))
        );
        w!(
            r##"<text x="{}" y="{:.2}" font-size="10" fill="#666" text-anchor="end">{} · {} record{} · {} genes</text>"##,
            ml - 10.0,
            yt + th / 2.0 + 11.0,
            fmt_len(g.total_len()),
            g.records.len(),
            if g.records.len() == 1 { "" } else { "s" },
            g.n_features()
        );
    }
    s.push_str(&labels_out); // labels last so nothing covers them

    // ---- scale bar
    let sc_y = sb_y;
    let bar = [1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6, 2e6, 5e6]
        .iter()
        .copied()
        .filter(|b| b * sc <= plot_w / 4.0)
        .last()
        .unwrap_or(1e3);
    let bw = bar * sc;
    w!(
        r##"<path d="M{ml},{a} V{sc_y} H{r} V{a}" fill="none" stroke="#333" stroke-width="1.2"/>"##,
        r = ml + bw,
        a = sc_y - 5.0
    );
    let bl = if bar >= 1e6 {
        format!("{} Mb", bar / 1e6)
    } else {
        format!("{} kb", bar / 1e3)
    };
    w!(
        r#"<text x="{:.2}" y="{}" font-size="11" text-anchor="middle">{}</text>"#,
        ml + bw / 2.0,
        sc_y + 16.0,
        bl
    );

    // ---- legend
    let ly = sc_y + 40.0;
    let items = [
        ("#f0a03c", "CDS"),
        ("#3aa17e", "tRNA"),
        ("#8e6bbf", "rRNA"),
        ("#d6609c", "other RNA"),
        ("#9a9a9a", "pseudo / other gene"),
    ];
    let mut lx = ml;
    for (c, l) in items {
        w!(
            r##"<rect x="{}" y="{}" width="12" height="12" fill="{}" stroke="#555" stroke-width="0.3"/>"##,
            lx,
            ly,
            c
        );
        w!(
            r#"<text x="{}" y="{}" font-size="11">{}</text>"#,
            lx + 17.0,
            ly + 10.0,
            l
        );
        lx += 22.0 + l.len() as f64 * 6.2;
    }
    w!(
        r##"<text x="{}" y="{}" font-size="10" fill="#555">gene arrows: + strand in upper half of the track, - strand in lower half</text>"##,
        ml,
        ly + 32.0
    );
    for (j, (g, l)) in [("gf", "forward"), ("gi", "inverted")].iter().enumerate() {
        let gx = ml + 620.0 + j as f64 * 190.0;
        w!(
            r##"<rect x="{}" y="{}" width="110" height="10" fill="url(#{})" stroke="#999" stroke-width="0.4"/>"##,
            gx,
            ly,
            g
        );
        w!(
            r#"<text x="{}" y="{}" font-size="10" text-anchor="middle">{:.0}%</text>"#,
            gx,
            ly + 22.0,
            vmin
        );
        w!(
            r#"<text x="{}" y="{}" font-size="10" text-anchor="middle">100%</text>"#,
            gx + 110.0,
            ly + 22.0
        );
        w!(
            r#"<text x="{}" y="{}" font-size="11" text-anchor="middle">{} identity</text>"#,
            gx + 55.0,
            ly - 5.0,
            l
        );
    }
    w!("</svg>");
    s
}
