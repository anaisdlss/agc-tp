
# AGC-TP: OTU Calculation and 16S Analysis

**Anaïs DELASSUS** – anais.delassus@etu.u-paris.fr  
Université Paris Cité

## Project Description

This project performs **OTU clustering** on amplified 16S sequences.  
The main steps are:

1.  Reading the compressed FASTA file (`.fasta.gz`)
2.  Sequence de-duplication (full-length dereplication)
3.  Identification of chimeric sequences (optional / not implemented this year)
4.  Abundance greedy clustering of non-chimeric sequences
5.  Writing the OTUs to a FASTA file
6.  Aligning the OTUs against a 16S reference bank (`mock_16S.fasta`) using **VSEARCH**, and automatic column annotation

The program is written in **Python 3.9+**, uses `nwalign3` for global alignments and `gzip` to read compressed files.

-----

## Installation

### Clone the Project

```bash
git clone https://github.com/anaisdlss/agc-tp.git
cd agc-tp
```

### Create and Activate the Conda Environment

```bash
conda env create -f environnement.yml
conda activate agc
```

-----

## Script Execution

### OTU Calculation

```bash
python3 agc/agc.py -i data/amplicon.fasta.gz -o OTU.fasta
```

You can change the minimum sequence length (400 by default) with `-s`
and change the minimum count for de-duplication (10 by default) with `-m`.

## OTU Alignment against the Reference Bank with VSEARCH

```bash
vsearch --usearch_global OTU.fasta \
--db data/mock_16S.fasta \
--id 0.8 \
--blast6out resultat.tsv
```

Where `--id` sets the minimum identity threshold and `--blast6out` generates the
results in CSV format.

## Automatic Annotation of TSV Columns

To add column headers to the VSEARCH output file `resultat.tsv`:

```bash
python3 agc/annotation.py
```

-----
