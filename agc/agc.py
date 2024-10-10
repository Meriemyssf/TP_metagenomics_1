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
import os
import gzip
import statistics
import textwrap
from pathlib import Path
from collections import Counter
from typing import Iterator, Dict, List
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Meriem YOUSSEF"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Meriem YOUSSEF"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Meriem YOUSSEF"
__email__ = "mariem.youssef1207@email.fr"
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


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
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
    with gzip.open(amplicon_file, "rt") as  f_in:
        seq = ""
        for line in f_in:
            if not line.startswith(">"):
                seq+=(line.strip())
            else:                
                if len(seq) >= minseqlen:
                    yield seq
                seq = ""
    if len(seq) >= minseqlen:
                yield seq    

def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count] of sequence with a count >= mincount and a length >= minseqlen.
    """
    dict_seq = dict()
    genr =  read_fasta(amplicon_file, minseqlen)
    for seq in genr:
        if seq in dict_seq:
            dict_seq[seq]+=1
        else :
            dict_seq[seq]=1
    
    for key in sorted(dict_seq, key=dict_seq.get, reverse = True):
        if dict_seq[key]>= mincount:
            yield key,dict_seq[key]
    

def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list:  (list) A list of aligned sequences in the format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """
    sequence1 = alignment_list[0]
    sequence2 = alignment_list[1]
    matches = 0
    for i in range(len(sequence1)) :
        if sequence1[i] == sequence2[i] :
            matches+=1

    id = ((matches/len(sequence1))* 100)
    return(id)

def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int, mincount: int, chunk_size: int, kmer_size: int) -> List:
    """Compute an abundance greedy clustering regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """
    OTUs_list = []
    seq_genr = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    first = True

    for seq in seq_genr:
        if first :
        # Add the first sequence as an OTU
            OTUs_list.append(seq)
            first = False
        else :
            is_OTU = True
            for otu in OTUs_list :
                # Align the sequences.
                alignment = nw.global_align(seq[0], otu[0], gap_open=-1, gap_extend=-1, matrix=str(Path(__file__).parent / "MATCH"))
                ident = get_identity(alignment)

                # Check the identity
                if ident > 97.0:
                    is_OTU=False
                    break          
            if is_OTU:
                OTUs_list.append(seq)
    return OTUs_list
    



def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    with open(output_file, "w") as f_out:
        for i in range(len(OTU_list)):
            otu, occ = OTU_list[i]
            f_out.write(f">OTU_{i+1} occurrence:{occ}\n")
            f_out.write(textwrap.fill(otu, width=80) + "\n")

#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    OTUs = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, 1, 1)
    write_OTU(OTUs, args.output_file) 

if __name__ == '__main__':
    main()
