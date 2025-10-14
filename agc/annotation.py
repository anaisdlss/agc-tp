import pandas as pd

# Chemin vers le fichier de sortie VSEARCH
tsv_file = "resultat.tsv"
tsv_out = "resultat.tsv"

# Colonnes BLAST6 standards
columns = [
    "query_id",
    "subject_id",
    "perc_identity",
    "align_length",
    "mismatches",
    "gap_opens",
    "q_start",
    "q_end",
    "s_start",
    "s_end",
    "evalue",
    "bit_score"
]

# Lire le TSV sans en-têtes
df = pd.read_csv(tsv_file, sep="\t", header=None, names=columns)

# Sauvegarder avec les en-têtes
df.to_csv(tsv_out, sep="\t", index=False)

print(f"Fichier annoté : {tsv_out}")
