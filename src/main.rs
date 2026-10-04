// seqlogo: render a WebLogo-style sequence logo (SVG) from a FASTA alignment.
//
// Usage:
//   seqlogo -i alignment.fasta -o logo.svg [options]
//
// Run `seqlogo --help` for the full option list.

use clap::{Parser, ValueEnum};
use figlet_rs::FIGfont;
use std::collections::BTreeMap;
use std::fmt::Write as FmtWrite;
use std::fs;
use std::process::exit;

/*
Gaurav Sablok
gsablok@proton.me
*/

// ---------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------

/// seqlogo - render a sequence logo (SVG) from a multiple sequence alignment
#[derive(Parser)]
#[command(
    name = "seqlogo",
    version,
    about = "render a sequence logo (SVG) from a multiple sequence alignment
      ************************************************
       Gaurav Sablok,
       Email: gsablok@proton.me
      ************************************************",
    after_help = "EXAMPLE:\n    seqlogo -i msa.fasta -o logo.svg -a protein -u bits -w 20 -H 200"
)]
struct Args {
    /// Input alignment in FASTA format (all rows same length)
    #[arg(short = 'i', long = "input", value_name = "FILE")]
    input: String,

    /// Output SVG file path
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    output: String,

    /// dna | rna | protein | auto
    #[arg(short = 'a', long = "alphabet", value_enum, default_value = "auto")]
    alphabet: Alphabet,

    /// Exclude gap ('-', '.') characters from frequency calculations
    /// (renormalizes over remaining symbols)
    #[arg(long = "ignore-gaps")]
    ignore_gaps: bool,

    /// bits | probability
    #[arg(short = 'u', long = "units", value_enum, default_value = "bits")]
    units: Units,

    /// Width per alignment column in pixels
    #[arg(short = 'w', long = "col-width", value_name = "PX", default_value_t = 24.0)]
    col_width: f64,

    /// Max stack height in pixels
    #[arg(short = 'H', long = "height", value_name = "PX", default_value_t = 220.0)]
    logo_height: f64,

    /// Title printed above the logo
    #[arg(short = 't', long = "title", value_name = "TEXT", default_value = "")]
    title: String,

    /// First alignment column to plot (1-based)
    #[arg(short = 's', long = "start", value_name = "N", default_value_t = 1)]
    start: usize,

    /// Last alignment column to plot (1-based) [default: last]
    #[arg(short = 'e', long = "end", value_name = "N")]
    end: Option<usize>,

    /// Put a position label every N columns
    #[arg(long = "tick-every", value_name = "N", default_value_t = 10)]
    tick_every: usize,

    /// Hide x-axis position numbers
    #[arg(long = "no-numbers")]
    no_numbers: bool,
}

#[derive(Clone, Copy, PartialEq, ValueEnum)]
enum Alphabet {
    Auto,
    Dna,
    Rna,
    #[value(alias = "aa", alias = "amino")]
    Protein,
}

#[derive(Clone, Copy, PartialEq, ValueEnum)]
enum Units {
    Bits,
    #[value(alias = "prob", alias = "p")]
    Probability,
}

fn print_banner() {
    if let Ok(font) = FIGfont::standard() {
        if let Some(figure) = font.convert("SeqLogo") {
            println!("{}", figure);
        }
    }
}

// ---------------------------------------------------------------------
// FASTA parsing
// ---------------------------------------------------------------------

fn parse_fasta(text: &str) -> Vec<(String, String)> {
    let mut records = Vec::new();
    let mut cur_name = String::new();
    let mut cur_seq = String::new();
    for line in text.lines() {
        let line = line.trim_end();
        if line.starts_with('>') {
            if !cur_name.is_empty() || !cur_seq.is_empty() {
                records.push((cur_name.clone(), cur_seq.clone()));
            }
            cur_name = line[1..].trim().to_string();
            cur_seq.clear();
        } else if !line.trim().is_empty() {
            cur_seq.push_str(line.trim());
        }
    }
    if !cur_name.is_empty() || !cur_seq.is_empty() {
        records.push((cur_name, cur_seq));
    }
    records
}

// ---------------------------------------------------------------------
// Alphabet + colors
// ---------------------------------------------------------------------

fn dna_alphabet() -> Vec<char> {
    vec!['A', 'C', 'G', 'T']
}
fn rna_alphabet() -> Vec<char> {
    vec!['A', 'C', 'G', 'U']
}
fn protein_alphabet() -> Vec<char> {
    "ACDEFGHIKLMNPQRSTVWY".chars().collect()
}

fn detect_alphabet(records: &[(String, String)]) -> Alphabet {
    let mut counts: BTreeMap<char, u64> = BTreeMap::new();
    for (_, seq) in records {
        for c in seq.chars() {
            let c = c.to_ascii_uppercase();
            if c == '-' || c == '.' {
                continue;
            }
            *counts.entry(c).or_insert(0) += 1;
        }
    }
    let acgtu: u64 = counts
        .iter()
        .filter(|(c, _)| matches!(c, 'A' | 'C' | 'G' | 'T' | 'U' | 'N'))
        .map(|(_, n)| *n)
        .sum();
    let total: u64 = counts.values().sum();
    if total == 0 {
        return Alphabet::Dna;
    }
    if (acgtu as f64) / (total as f64) > 0.9 {
        let has_u = counts.contains_key(&'U');
        let has_t = counts.contains_key(&'T');
        if has_u && !has_t {
            Alphabet::Rna
        } else {
            Alphabet::Dna
        }
    } else {
        Alphabet::Protein
    }
}

// Classic WebLogo-style coloring.
fn color_for(alphabet: Alphabet, symbol: char) -> &'static str {
    match alphabet {
        Alphabet::Dna | Alphabet::Rna | Alphabet::Auto => match symbol {
            'A' => "#2ca02c",       // green
            'C' => "#1f77b4",       // blue
            'G' => "#ff9900",       // orange
            'T' | 'U' => "#d62728", // red
            _ => "#7f7f7f",         // gray (N / ambiguous)
        },
        Alphabet::Protein => match symbol {
            // polar (green)
            'G' | 'S' | 'T' | 'Y' | 'C' => "#2ca02c",
            // basic (blue)
            'K' | 'R' | 'H' => "#1f77b4",
            // acidic (red)
            'D' | 'E' => "#d62728",
            // hydrophobic (black)
            'A' | 'V' | 'L' | 'I' | 'P' | 'W' | 'F' | 'M' => "#222222",
            // amide/other (purple)
            'N' | 'Q' => "#9467bd",
            _ => "#7f7f7f",
        },
    }
}

// ---------------------------------------------------------------------
// Logo math
// ---------------------------------------------------------------------

struct ColumnStack {
    // (symbol, height_in_units, color) sorted ascending (bottom to top order
    // is handled at render time)
    letters: Vec<(char, f64)>,
    total_height_units: f64,
}

fn compute_column(
    seqs_col: &[char],
    alphabet_syms: &[char],
    ignore_gaps: bool,
    units: Units,
    n_seqs: usize,
) -> ColumnStack {
    let mut counts: BTreeMap<char, u64> = BTreeMap::new();
    let mut n_counted: u64 = 0;
    for &c in seqs_col {
        let c = c.to_ascii_uppercase();
        let is_gap = c == '-' || c == '.' || c == '*';
        if is_gap && ignore_gaps {
            continue;
        }
        *counts.entry(c).or_insert(0) += 1;
        n_counted += 1;
    }
    if n_counted == 0 {
        return ColumnStack {
            letters: vec![],
            total_height_units: 0.0,
        };
    }

    let n_alpha = alphabet_syms.len() as f64;
    let max_bits = n_alpha.log2();

    // Shannon entropy over observed symbols (including gaps as their own
    // symbol unless ignored -- gaps simply don't get drawn, but they still
    // reduce the information content of the column, matching WebLogo
    // behaviour for gappy columns).
    let mut entropy = 0.0f64;
    for (_, &cnt) in counts.iter() {
        let p = cnt as f64 / n_counted as f64;
        if p > 0.0 {
            entropy -= p * p.log2();
        }
    }

    // Small-sample correction (Schneider & Stephens, 1990):
    // e(n) = (1/ln2) * (s-1)/(2n)
    let s = n_alpha;
    let correction = (1.0 / std::f64::consts::LN_2) * (s - 1.0) / (2.0 * n_seqs as f64);
    let mut r_bits = max_bits - entropy - correction;
    if r_bits < 0.0 {
        r_bits = 0.0;
    }
    if r_bits > max_bits {
        r_bits = max_bits;
    }

    let mut letters: Vec<(char, f64)> = Vec::new();
    for &sym in alphabet_syms {
        if let Some(&cnt) = counts.get(&sym) {
            let p = cnt as f64 / n_counted as f64;
            let h = match units {
                Units::Bits => p * r_bits,
                Units::Probability => p,
            };
            if h > 1e-9 {
                letters.push((sym, h));
            }
        }
    }
    // Sort ascending by height so tallest ends up on top when stacked
    // bottom-to-top... we actually want tallest on top when reading
    // top-down, so sort ascending and render from the top downward in
    // that order (see render_svg).
    letters.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let total_height_units = match units {
        Units::Bits => r_bits,
        Units::Probability => 1.0,
    };

    ColumnStack {
        letters,
        total_height_units,
    }
}

// ---------------------------------------------------------------------
// SVG rendering
// ---------------------------------------------------------------------

fn escape_xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn render_svg(
    columns: &[ColumnStack],
    alphabet_syms: &[char],
    alphabet: Alphabet,
    units: Units,
    col_width: f64,
    plot_height: f64,
    title: &str,
    first_pos: usize,
    tick_every: usize,
    show_numbers: bool,
) -> String {
    let n_alpha = alphabet_syms.len() as f64;
    let max_units = match units {
        Units::Bits => n_alpha.log2(),
        Units::Probability => 1.0,
    };

    let left_margin = 46.0;
    let right_margin = 16.0;
    let top_margin = if title.is_empty() { 20.0 } else { 44.0 };
    let bottom_margin = if show_numbers { 46.0 } else { 26.0 };

    let plot_width = col_width * columns.len() as f64;
    let width = left_margin + plot_width + right_margin;
    let height = top_margin + plot_height + bottom_margin;

    let px_per_unit = plot_height / max_units;
    let baseline_y = top_margin + plot_height;

    let mut svg = String::new();
    let _ = write!(
        svg,
        r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w:.1} {h:.1}" width="{w:.1}" height="{h:.1}" font-family="Helvetica, Arial, sans-serif">
<rect x="0" y="0" width="{w:.1}" height="{h:.1}" fill="#ffffff"/>
"##,
        w = width,
        h = height
    );

    if !title.is_empty() {
        let _ = write!(
            svg,
            r##"<text x="{cx:.1}" y="22" text-anchor="middle" font-size="16" font-weight="bold" fill="#111">{t}</text>
"##,
            cx = width / 2.0,
            t = escape_xml(title)
        );
    }

    // Y axis line + ticks
    let _ = write!(
        svg,
        r##"<line x1="{x:.1}" y1="{y1:.1}" x2="{x:.1}" y2="{y2:.1}" stroke="#333" stroke-width="1"/>
"##,
        x = left_margin,
        y1 = top_margin,
        y2 = baseline_y
    );

    let n_yticks = 4;
    for k in 0..=n_yticks {
        let frac = k as f64 / n_yticks as f64;
        let val = max_units * frac;
        let y = baseline_y - plot_height * frac;
        let label = match units {
            Units::Bits => format!("{:.1}", val),
            Units::Probability => format!("{:.2}", val),
        };
        let _ = write!(
            svg,
            r##"<line x1="{x0:.1}" y1="{y:.1}" x2="{x1:.1}" y2="{y:.1}" stroke="#ccc" stroke-width="1"/>
<text x="{lx:.1}" y="{ty:.1}" text-anchor="end" font-size="10" fill="#333">{label}</text>
"##,
            x0 = left_margin,
            x1 = left_margin + plot_width,
            y = y,
            lx = left_margin - 6.0,
            ty = y + 3.5,
            label = label
        );
    }

    let ylabel = match units {
        Units::Bits => "bits",
        Units::Probability => "probability",
    };
    let _ = write!(
        svg,
        r##"<text x="14" y="{cy:.1}" text-anchor="middle" font-size="11" fill="#333" transform="rotate(-90 14 {cy:.1})">{yl}</text>
"##,
        cy = top_margin + plot_height / 2.0,
        yl = ylabel
    );

    // Baseline (x axis)
    let _ = write!(
        svg,
        r##"<line x1="{x0:.1}" y1="{y:.1}" x2="{x1:.1}" y2="{y:.1}" stroke="#333" stroke-width="1"/>
"##,
        x0 = left_margin,
        x1 = left_margin + plot_width,
        y = baseline_y
    );

    // Columns
    for (i, col) in columns.iter().enumerate() {
        let col_x0 = left_margin + i as f64 * col_width;
        let empty_units = max_units - col.total_height_units;
        let mut y_cursor = top_margin + empty_units * px_per_unit; // top of the drawn stack

        // letters sorted ascending by height; render tallest first (topmost)
        // by iterating in reverse.
        let pad_x = col_width * 0.06;
        let draw_w = col_width - 2.0 * pad_x;
        for &(sym, h_units) in col.letters.iter().rev() {
            let h_px = h_units * px_per_unit;
            if h_px < 0.4 {
                y_cursor += h_px;
                continue;
            }
            let color = color_for(alphabet, sym);
            // Reference glyph box is 100 x 100 units; font-size 100,
            // baseline at y=88 gives a good visual fill of the box for
            // most sans-serif capitals.
            let sx = draw_w / 100.0;
            let sy = h_px / 72.0; // cap-height ~0.72 em: glyph spans y=16..88 of the 100-unit box
            let tx = col_x0 + pad_x;
            let ty = y_cursor - 16.0 * sy;
            let _ = write!(
                svg,
                r##"<g transform="translate({tx:.2},{ty:.2}) scale({sx:.4},{sy:.4})"><text x="0" y="88" font-size="100" font-weight="bold" textLength="100" lengthAdjust="spacingAndGlyphs" fill="{color}">{sym}</text></g>
"##,
                tx = tx,
                ty = ty,
                sx = sx,
                sy = sy,
                color = color,
                sym = sym
            );
            y_cursor += h_px;
        }

        // x-axis position label
        if show_numbers && (i == 0 || (i + 1) % tick_every == 0) {
            let pos = first_pos + i;
            let _ = write!(
                svg,
                r##"<text x="{cx:.1}" y="{ly:.1}" text-anchor="middle" font-size="10" fill="#333">{pos}</text>
"##,
                cx = col_x0 + col_width / 2.0,
                ly = baseline_y + 14.0,
                pos = pos
            );
            let _ = write!(
                svg,
                r##"<line x1="{cx:.1}" y1="{y0:.1}" x2="{cx:.1}" y2="{y1:.1}" stroke="#999" stroke-width="1"/>
"##,
                cx = col_x0 + col_width / 2.0,
                y0 = baseline_y,
                y1 = baseline_y + 4.0
            );
        }
    }

    svg.push_str("</svg>\n");
    svg
}

// ---------------------------------------------------------------------
// main
// ---------------------------------------------------------------------

fn main() {
    print_banner();

    let args = Args::parse();

    let text = match fs::read_to_string(&args.input) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("error: could not read '{}': {}", args.input, e);
            exit(1);
        }
    };

    let records = parse_fasta(&text);
    if records.is_empty() {
        eprintln!("error: no FASTA records found in '{}'", args.input);
        exit(1);
    }

    let aln_len = records[0].1.chars().count();
    for (name, seq) in &records {
        let l = seq.chars().count();
        if l != aln_len {
            eprintln!(
                "error: sequence '{}' has length {} but expected {} (not aligned / not equal length)",
                name, l, aln_len
            );
            exit(1);
        }
    }

    let alphabet = match args.alphabet {
        Alphabet::Auto => detect_alphabet(&records),
        other => other,
    };
    let alphabet_syms = match alphabet {
        Alphabet::Dna => dna_alphabet(),
        Alphabet::Rna => rna_alphabet(),
        Alphabet::Protein => protein_alphabet(),
        Alphabet::Auto => dna_alphabet(),
    };

    let start = args.start.max(1);
    let end = args.end.unwrap_or(aln_len).min(aln_len);
    if start > end {
        eprintln!("error: --start ({}) is greater than --end ({})", start, end);
        exit(1);
    }

    let n_seqs = records.len();
    let seq_chars: Vec<Vec<char>> = records
        .iter()
        .map(|(_, s)| s.chars().collect())
        .collect();

    let mut columns = Vec::with_capacity(end - start + 1);
    for pos in start..=end {
        let idx = pos - 1;
        let col_chars: Vec<char> = seq_chars.iter().map(|row| row[idx]).collect();
        let stack = compute_column(
            &col_chars,
            &alphabet_syms,
            args.ignore_gaps,
            args.units,
            n_seqs,
        );
        columns.push(stack);
    }

    let svg = render_svg(
        &columns,
        &alphabet_syms,
        alphabet,
        args.units,
        args.col_width,
        args.logo_height,
        &args.title,
        start,
        args.tick_every.max(1),
        !args.no_numbers,
    );

    if let Err(e) = fs::write(&args.output, svg) {
        eprintln!("error: could not write '{}': {}", args.output, e);
        exit(1);
    }

    eprintln!(
        "wrote {} ({} sequences, {} columns plotted [{}-{}], alphabet={})",
        args.output,
        n_seqs,
        columns.len(),
        start,
        end,
        match alphabet {
            Alphabet::Dna => "DNA",
            Alphabet::Rna => "RNA",
            Alphabet::Protein => "protein",
            Alphabet::Auto => "auto",
        }
    );
}
