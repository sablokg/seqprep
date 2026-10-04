# seqlogo

A small, dependency-free Rust CLI that turns a multiple sequence alignment (FASTA) into a WebLogo-style sequence logo, rendered as SVG.

## Build

```
cargo build --release
# binary at target/release/seqlogo
```

Requires only the Rust standard library — no crates.io dependencies.

## Usage

```
seqlogo -i alignment.fasta -o logo.svg [OPTIONS]
```

```
  ____                   _                             
 / ___|    ___    __ _  | |       ___     __ _    ___  
 \___ \   / _ \  / _` | | |      / _ \   / _` |  / _ \ 
  ___) | |  __/ | (_| | | |___  | (_) | | (_| | | (_) |
 |____/   \___|  \__, | |_____|  \___/   \__, |  \___/ 
                    |_|                  |___/         

render a sequence logo (SVG) from a multiple sequence alignment
      ************************************************
       Gaurav Sablok,
       Email: gsablok@proton.me
      ************************************************

Usage: seqlogo [OPTIONS] --input <FILE> --output <FILE>

Options:
  -i, --input <FILE>         Input alignment in FASTA format (all rows same length)
  -o, --output <FILE>        Output SVG file path
  -a, --alphabet <ALPHABET>  dna | rna | protein | auto [default: auto] [possible values: auto, dna, rna, protein]
      --ignore-gaps          Exclude gap ('-', '.') characters from frequency calculations (renormalizes over remaining symbols)
  -u, --units <UNITS>        bits | probability [default: bits] [possible values: bits, probability]
  -w, --col-width <PX>       Width per alignment column in pixels [default: 24]
  -H, --height <PX>          Max stack height in pixels [default: 220]
  -t, --title <TEXT>         Title printed above the logo [default: ""]
  -s, --start <N>            First alignment column to plot (1-based) [default: 1]
  -e, --end <N>              Last alignment column to plot (1-based) [default: last]
      --tick-every <N>       Put a position label every N columns [default: 10]
      --no-numbers           Hide x-axis position numbers
  -h, --help                 Print help
  -V, --version              Print version

EXAMPLE:
    seqlogo -i msa.fasta -o logo.svg -a protein -u bits -w 20 -H 200
```


## Example

```
./target/release/seqlogo -i examples/dna_aln.fasta -o dna_logo.svg -t "Demo DNA motif"
./target/release/seqlogo -i examples/protein_aln.fasta -o protein_logo.svg -a protein -w 22 -H 180
```

Sample inputs/outputs are in `examples/`.

