import argparse
import sys

# ---------------------------------------------------------------
# Standard genetic code: built from the 64 codons in TCAG order
# ('*' marks a stop codon)
# ---------------------------------------------------------------
BASES = "TCAG"
AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {}
_index = 0
for first in BASES:
    for second in BASES:
        for third in BASES:
            CODON_TABLE[first + second + third] = AMINO_ACIDS[_index]
            _index += 1

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


def read_fasta(path):
    """Reads the first record of a FASTA file. Returns (header, sequence)."""
    header = "No header"
    sequence = []
    with open(path, "r") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if sequence:      # stop at the second record
                    break
                header = line[1:]
            else:
                sequence.append(line)
    return header, "".join(sequence)


def clean_sequence(raw):
    """Uppercases the sequence, removes spaces/newlines and checks it is valid DNA."""
    seq = "".join(raw.split()).upper()
    invalid = sorted(set(seq) - set("ACGT"))
    if not seq:
        raise ValueError("The sequence is empty.")
    if invalid:
        raise ValueError(f"Invalid DNA characters found: {', '.join(invalid)}")
    return seq


def count_bases(seq):
    """Counts A, T, G and C."""
    return {base: seq.count(base) for base in "ATGC"}


def gc_content(seq):
    """Returns the GC percentage of the sequence."""
    return (seq.count("G") + seq.count("C")) / len(seq) * 100


def reverse_complement(seq):
    """Returns the reverse complement strand."""
    return "".join(COMPLEMENT[base] for base in reversed(seq))


def transcribe(seq):
    """DNA to RNA: replaces T with U."""
    return seq.replace("T", "U")


def translate(seq):
    """Translates from the first ATG until a stop codon. Returns a protein string."""
    start = seq.find("ATG")
    if start == -1:
        return None
    protein = []
    for i in range(start, len(seq) - 2, 3):
        amino = CODON_TABLE[seq[i:i + 3]]
        if amino == "*":
            break
        protein.append(amino)
    return "".join(protein)


def find_motif(seq, motif):
    """Returns all 1-based start positions of a motif (overlaps included)."""
    positions = []
    for i in range(len(seq) - len(motif) + 1):
        if seq[i:i + len(motif)] == motif:
            positions.append(i + 1)
    return positions


def print_report(header, seq, motif=None):
    """Prints the full analysis."""
    counts = count_bases(seq)
    protein = translate(seq)

    print("=" * 50)
    print("DNA SEQUENCE REPORT")
    print("=" * 50)
    print(f"Name           : {header}")
    print(f"Length         : {len(seq)} bp")
    print(f"Base counts    : A={counts['A']}  T={counts['T']}  G={counts['G']}  C={counts['C']}")
    print(f"GC content     : {gc_content(seq):.2f}%")
    print("-" * 50)
    print(f"Sequence (5'-3')   : {seq}")
    print(f"Reverse complement : {reverse_complement(seq)}")
    print(f"RNA (transcribed)  : {transcribe(seq)}")
    print("-" * 50)
    if protein is None:
        print("Protein        : No start codon (ATG) found")
    else:
        print(f"Protein        : {protein}  ({len(protein)} amino acids)")

    if motif:
        positions = find_motif(seq, motif)
        print("-" * 50)
        if positions:
            print(f"Motif {motif} found at position(s): {', '.join(map(str, positions))}")
        else:
            print(f"Motif {motif} not found")
    print("=" * 50)


def parse_args():
    parser = argparse.ArgumentParser(
        description="bio-dna-analyzer: a simple DNA sequence analysis tool."
    )
    parser.add_argument("--seq", help="DNA sequence typed directly, e.g. ATGGCC")
    parser.add_argument("--file", help="path to a FASTA file")
    parser.add_argument("--motif", help="optional motif to search for, e.g. GAATTC")
    return parser.parse_args()


def main():
    args = parse_args()

    try:
        if args.file:
            header, raw = read_fasta(args.file)
        elif args.seq:
            header, raw = "User input", args.seq
        else:
            header = "User input"
            raw = input("Enter a DNA sequence (A, T, G, C): ")

        seq = clean_sequence(raw)
        motif = clean_sequence(args.motif) if args.motif else None
        print_report(header, seq, motif)

    except FileNotFoundError:
        print(f"Error: file '{args.file}' was not found.")
        sys.exit(1)
    except ValueError as error:
        print(f"Error: {error}")
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nStopped by user.")


if __name__ == "__main__":
    main()
