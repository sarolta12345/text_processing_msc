from Bio.Seq import Seq
from Bio import SeqIO

# 1. Eredeti és rövid nevek beolvasása
names = {}

with open("output/names.tsv", "r") as file:
    for line in file:
        short_name, original_name = line.strip().split("\t")
        names[short_name] = original_name

# 2. DNS-szekvenciák beolvasása a PHYLIP visszaolvasásából
dna_sequences = {}

with open("output/genes_readback.txt", "r") as file:
    for line in file:
        short_name, dna = line.strip().split(maxsplit=1)
        dna_sequences[short_name] = dna

# 3. Saját fehérjefordítás beolvasása
my_proteins = {}

for record in SeqIO.parse("output/genes_prot.fasta", "fasta"):
    my_proteins[record.id] = str(record.seq)

# 4. Biopythonos fordítás és összehasonlítás
for short_name, dna in dna_sequences.items():

    original_name = names[short_name]

    bio_protein = str(Seq(dna).translate())
    my_protein = my_proteins[original_name]

    print("\nGén:", original_name)

    if bio_protein == my_protein:
        print("A két fordítás megegyezik!")
    else:
        print("ELTÉRÉS TALÁLHATÓ!")

        print("Saját fordítás: ", my_protein)
        print("Biopython:      ", bio_protein)

    # N bázis vizsgálata
    if "N" in dna:
        print("A DNS-szekvencia N bázist tartalmaz.")
        print("Biopythonban az ambiguis kodon általában X lesz.")

    # Stopkodon vizsgálata
    if "*" in bio_protein:
        print("Stopkodon található a fehérjeszekvenciában.")

print("\nAz ellenőrzés befejeződött.")
