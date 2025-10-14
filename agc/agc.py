#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import gzip
import textwrap
from pathlib import Path
from collections import Counter
from typing import Iterator, List
import nwalign3 as nw
import numpy as np
np.int = int
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/

__author__ = "Anaïs DELASSUS"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Anaïs DELASSUS"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Anaïs DELASSUS"
__email__ = "anais.delassus@etu.u-paris.fr"
__status__ = "Developpement"


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments():  # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     usage=f"{sys.argv[0]} -h")

    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file',
                        type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen',
                        type=int, default=400,
                        help="Minimum sequence length for dereplication \
                            (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int,
                        default=10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=Path,
                        default=Path("OTU.fasta"), help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    """Read a compressed fasta and extract all fasta sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :return: A generator object that provides the Fasta sequences (str).
    """
    with gzip.open(str(amplicon_file), "rt") as fasta:
        seq_lines = []
        for line in fasta:
            line = line.strip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if seq_lines:
                    seq = "".join(s.strip() for s in seq_lines).upper()
                    if len(seq) >= minseqlen:
                        yield seq
                    seq_lines = []
                continue
            seq_lines.append(line)
        if seq_lines:
            seq = "".join(s.strip() for s in seq_lines).upper()
            if len(seq) >= minseqlen:
                yield seq


def dereplication_fulllength(amplicon_file: Path, minseqlen: int,
                             mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count] of 
    sequence with a count >= mincount and a length >= minseqlen.
    """
    counts = Counter(read_fasta(amplicon_file, minseqlen))
    filtered_sorted = sorted(
        [[seq, count] for seq, count in counts.items() if count >= mincount],
        key=lambda x: x[1], reverse=True)
    yield from filtered_sorted


def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list:  (list) A list of aligned sequences in the format 
    ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """
    if not alignment_list or len(alignment_list) < 2:
        return 0.0
    a, b = alignment_list[0], alignment_list[1]
    min_l = min(len(a), len(b))
    if min_l == 0:
        return 0.0
    identical = sum(1 for a, b in zip(a[:min_l], b[:min_l]) if a == b)
    return (identical / min_l) * 100.0


def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int,
                                mincount: int, _chunk_size: int,
                                _kmer_size: int) -> List:
    """Compute an abundance greedy clustering regarding sequence count and 
    identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """
    otu_list = []
    matrix_path = str(Path(__file__).parent / "MATCH")
    for seq, cnt in dereplication_fulllength(amplicon_file, minseqlen,
                                             mincount):
        def is_similar(otu_seq: str, query_seq=seq) -> bool:
            aln = nw.global_align(query_seq, otu_seq, gap_open=-1,
                                  gap_extend=-1, matrix=matrix_path)
            aln_pair = [aln[0], aln[1]]
            return get_identity(aln_pair) >= 97.0
        if not any(is_similar(otu_seq) for otu_seq, _ in otu_list):
            otu_list.append([seq, cnt])
    return otu_list


def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    with open(str(output_file), "w", encoding="utf-8") as out:
        for idx, (seq, count) in enumerate(OTU_list, 1):
            out.write(f">OTU_{idx} occurrence:{count}\n")
            out.write(textwrap.fill(seq, width=80) + "\n")


# ==============================================================
# Main program
# ==============================================================
def main():  # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # calculer les OTUs
    otu_list = abundance_greedy_clustering(
        amplicon_file=args.amplicon_file,
        minseqlen=args.minseqlen,
        mincount=args.mincount,
        chunk_size=100,
        kmer_size=8
    )
    print(f"Nombre d'OTUs : {len(otu_list)}")
    for i, (seq, count) in enumerate(otu_list, start=1):
        print(f"OTU_{i}: count={count}, length={len(seq)}")

    # Ecriture fichier fasta de sortie
    write_OTU(otu_list, args.output_file)
    print(f"Fichier OTU écrit : {args.output_file}")


if __name__ == '__main__':
    main()
