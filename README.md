# AGC-TP: OTU Clustering et Analyse 16S

## Auteurs
**Anaïs DELASSUS** – anais.delassus@etu.u-paris.fr  
Université Paris Cité

## Description du projet
Ce projet réalise un **clustering OTU** sur des séquences 16S amplifiées.  
Les étapes principales sont :

1. Lecture du fichier FASTA compressé (`.fasta.gz`)  
2. Dé-duplication des séquences (full-length dereplication)  
3. Identification des séquences chimériques (optionnelle / non implémentée cette année)  
4. Regroupement glouton des séquences non-chimériques (abundance greedy clustering)  
5. Écriture des OTUs dans un fichier FASTA  
6. Alignement des OTUs sur une banque de références 16S (`mock_16S.fasta`) avec **VSEARCH**, et annotation automatique des colonnes  

Le programme est écrit en **Python 3.9+**, utilise `nwalign3` pour les alignements globaux et `gzip` pour lire les fichiers compressés.

---

## Installation

### Cloner le projet
```bash
git clone https://github.com/anaisdlss/agc-tp.git
cd agc-tp

### Créer et activer l'environnement Conda
```bash
conda env create -f environment.yml
conda activate agc
````

---

## Exécution des scripts

### Calcul des OTUs
```bash
python3 agc/agc.py -i data/amplicon.fasta.gz -o OTU.fasta
```
Vous pouvez changer la longueur minimale des séquences (400 par défaut) avec -s
et changer le comptage minimal de la déduplication (10 par défaut) avec -m.

## Alignement OTUs sur la banque de référence avec VSEARCH
```bash
vsearch --usearch_global OTU.fasta \
--db data/mock_16S.fasta \
--id 0.8 \
--blast6out resultat.tsv
````
Avec `--id` qui fixe seuil d'identité minimal et `--blast6out` qui génère les
résultats au format csv.

## Annotation automatique des colonnes TSV
Pour ajouter les noms des colonnes au fichier VSEARCH:
```bash
python3 agc/annotation.py
```



