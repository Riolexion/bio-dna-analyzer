# bio-dna-analyzer

> Analyze DNA sequences, right from your terminal.

Biologists and students often need quick answers about a DNA sequence: how long is it, what is its GC content, and what protein does it code for? **bio-dna-analyzer** is a zero-dependency Python script that answers these questions in seconds, with no software to install.

## Features

- **Zero Dependencies:** Built entirely with Python's standard library. No `pip install` required.
- **Base Counting:** Counts A, T, G and C and reports the sequence length.
- **GC Content:** Calculates the GC percentage, a key measure in genomics.
- **Strand Tools:** Shows the reverse complement and the RNA (transcribed) sequence.
- **Translation:** Finds the first start codon (ATG) and translates to a protein using the standard genetic code.
- **Motif Search:** Finds every position of a motif such as the EcoRI site `GAATTC`.
- **FASTA Support:** Reads sequences from standard `.fasta` files.
- **Input Validation:** Clear error messages for invalid characters or missing files.

## Quickstart

1. Clone the repository:

   ```bash
   git clone https://github.com/yourusername/bio-dna-analyzer.git
   cd bio-dna-analyzer
   ```

2. Run it on the included sample file:

   ```bash
   python DNA_Analyzer.py --file sample.fasta --motif GAATTC
   ```

   On macOS/Linux you may need `python3` instead of `python`.

## Usage

| Flag | Description |
|------|-------------|
| `--seq` | Type a DNA sequence directly |
| `--file` | Path to a FASTA file |
| `--motif` | Optional motif to search for |

If you give neither `--seq` nor `--file`, the program asks you to type a sequence.

Examples:

```bash
python DNA_Analyzer.py --seq ATGGCCATTGTAATGGGCCGCTGA
python DNA_Analyzer.py --file sample.fasta
python DNA_Analyzer.py --seq ATGGAATTCGAATTC --motif GAATTC
```

## Sample output

```
==================================================
DNA SEQUENCE REPORT
==================================================
Name           : sample_gene demo sequence
Length         : 48 bp
Base counts    : A=12  T=10  G=16  C=10
GC content     : 54.17%
--------------------------------------------------
Sequence (5'-3')   : ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGGAATTCCGA
Reverse complement : TCGGAATTCCTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT
RNA (transcribed)  : AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGGAAUUCCGA
--------------------------------------------------
Protein        : MAIVMGR  (7 amino acids)
Motif GAATTC found at position(s): 40
==================================================
```

## Requirements

- Python 3.6 or newer
- No third-party packages

## How it works

1. Reads and cleans the sequence (uppercase, no spaces, only A/T/G/C allowed).
2. Counts bases and calculates GC% = (G + C) / length x 100.
3. Builds the reverse complement by reversing the sequence and swapping A-T and G-C.
4. Transcribes to RNA by replacing T with U.
5. Finds the first `ATG`, reads codons three bases at a time, and stops at a stop codon.

## Ideas for improvement

- Translate all 6 reading frames
- Calculate molecular weight of the protein
- Compare two sequences and report the percent identity
- Plot base composition with `matplotlib`

## License

Free to use and modify for learning purposes.
