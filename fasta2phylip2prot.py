import argparse
import re
import os
from itertools import product


# 1. FASTA beolvasása
def read_fasta(filename):
    sequences = {}
    name = None

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                name = line[1:].split()[0]
                sequences[name] = ""
            elif name is not None:
                sequences[name] += line

    return sequences


# 2. Szekvenciák tisztítása
def clean_sequences(sequences):
    cleaned = {}

    for name, seq in sequences.items():
        seq = seq.upper()
        seq = re.sub(r"[^A-Z]", "", seq)
        cleaned[name] = seq

    return cleaned


# 3. Szekvenciák ellenőrzése
def validate_sequences(sequences):
    lengths = [len(seq) for seq in sequences.values()]

    if not lengths:
        raise ValueError("A FASTA fajl ures!")

    if len(set(lengths)) != 1:
        raise ValueError("A szekvenciak hossza nem azonos!")

    for name, seq in sequences.items():
        if len(seq) % 3 != 0:
            raise ValueError(
                f"A(z) {name} hossza nem oszthato 3-mal!"
            )

    print("A szekvenciak hossza megfelelo!")
    return lengths[0]


# 4. Rövid nevek létrehozása
def create_short_names(sequences):
    short_sequences = {}
    names = {}

    for i, (original_name, seq) in enumerate(
        sequences.items(), start=1
    ):
        short_name = f"seq{i}"
        short_sequences[short_name] = seq
        names[short_name] = original_name

    return short_sequences, names


# 5. Névmegfeleltetés mentése
def write_names(names, filename):
    with open(filename, "w") as file:
        for short_name, original_name in names.items():
            file.write(f"{short_name}\t{original_name}\n")


# 6. PHYLIP írása
def write_phylip(sequences, filename):
    sequence_length = len(next(iter(sequences.values())))

    with open(filename, "w") as file:
        file.write(f" {len(sequences)} {sequence_length}\n")

        for name, seq in sequences.items():
            file.write(f"{name:<10}{seq}\n")

    print("PHYLIP fajl elkeszult:", filename)


# 7. PHYLIP visszaolvasása
def read_phylip(filename):
    sequences = {}

    with open(filename, "r") as file:
        header = file.readline().split()
        number_of_sequences = int(header[0])
        sequence_length = int(header[1])

        for line in file:
            name = line[:10].strip()
            seq = line[10:].strip()

            if name:
                sequences[name] = seq

    if len(sequences) != number_of_sequences:
        raise ValueError("Hiba a PHYLIP visszaolvasasakor!")

    for seq in sequences.values():
        if len(seq) != sequence_length:
            raise ValueError("Hibas szekvenciahossz a PHYLIP fajlban!")

    print("PHYLIP visszaolvasasa sikeres!")
    return sequences


# 8. Kodontábla létrehozása
def create_codon_table():
    amino = (
        "FFLLSSSSYY**CC*W"
        "LLLLPPPPHHQQRRRR"
        "IIIMTTTTNNKKSSRR"
        "VVVVAAAADDEEGGGG"
    )

    bases = "TCAG"
    codons = [
        a + b + c
        for a, b, c in product(bases, repeat=3)
    ]

    return dict(zip(codons, amino))


# 9. DNS transzlációja
def translate_sequences(sequences):
    codon_table = create_codon_table()
    proteins = {}

    for name, dna in sequences.items():
        protein = ""

        for i in range(0, len(dna), 3):
            codon = dna[i:i+3]

            if "N" in codon:
                amino_acid = "X"
            else:
                amino_acid = codon_table.get(codon, "X")

            protein += amino_acid

        proteins[name] = protein

        if "*" in protein[:-1]:
            print(f"Figyelem: {name} belso stopkodont tartalmaz!")

    return proteins


# 10. Fehérje FASTA írása
def write_protein_fasta(proteins, names, filename):
    with open(filename, "w") as file:
        for short_name, protein in proteins.items():
            original_name = names[short_name]

            file.write(f">{original_name}\n")

            for i in range(0, len(protein), 60):
                file.write(protein[i:i+60] + "\n")

    print("Feherje FASTA elkeszult:", filename)


# 11. Főprogram
def main():
    parser = argparse.ArgumentParser(
        description="FASTA -> PHYLIP -> Protein FASTA"
    )

    parser.add_argument("input", help="Bemeneti FASTA fajl")
    parser.add_argument("output_prefix", help="Kimeneti prefix")

    args = parser.parse_args()

    try:
        sequences = read_fasta(args.input)
        print("Beolvasott genek:", len(sequences))

        sequences = clean_sequences(sequences)
        print("Tisztitas kesz!")

        validate_sequences(sequences)

        short_sequences, names = create_short_names(sequences)

        output_dir = os.path.dirname(args.output_prefix)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        phylip_file = args.output_prefix + ".phy"
        names_file = args.output_prefix + "_names.tsv"
        protein_file = args.output_prefix + "_prot.fasta"

        write_names(names, names_file)
        write_phylip(short_sequences, phylip_file)

        # PHYLIP visszaolvasása
        dna_sequences = read_phylip(phylip_file)

        # Transzláció a visszaolvasott DNS-ből
        proteins = translate_sequences(dna_sequences)

        # Fehérje FASTA írása
        write_protein_fasta(proteins, names, protein_file)

        print("\nA feldolgozas sikeresen befejezodott!")

    except (ValueError, OSError) as error:
        print("HIBA:", error)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
