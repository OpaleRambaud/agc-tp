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
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Opale Rambaud"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Opale Rambaud"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Opale Rambaud"
__email__ = "opale.rambaud@gmail.com"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    with gzip.open(amplicon_file, "r") as filin:
        for line in filin:
            if len(line) >= minseqlen:
                yield line.strip()
                
 def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    occ_dict = {}
    for seq in read_fasta(amplicon_file, minseqlen):
        if seq not in occ_dict:
            occ_dict[seq] = 1
        else:
            occ_dict[seq] += 1
     for seq, count in occ_dict.items():
        if count >= mincount:
            yield [seq, count]        
            
def get_chunks(sequence, chunk_size):
    chunk_num = int(len(sequence) / chunk_size)
    chunk = []
    if chunk_num < 4:
        chunk_num = int(len(sequence) / 4)
    count = 0
    for i in range(chunk_num):
        chunk.append(sequence[count:count+chunk_size])
        count += chunk_size
    return chunk

def cut_kmer(sequence, kmer_size):
    for i in range(len(sequence) - kmer_size + 1):
        yield sequence[i : i + kmer_size]

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    for kmer in cut_kmer(sequence, kmer_size):
        if kmer not in kmer_dict:
            kmer_dict[kmer] = []
        kmer_dict[kmer].append(id_seq)
    return kmer_dict            

def search_mates(kmer_dict, sequence, kmer_size):
    return [i[0] for i in Counter([ids for kmer in cut_kmer(sequence, kmer_size)
        if kmer in kmer_dict for ids in kmer_dict[kmer]]).most_common(8)]            
#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()


if __name__ == '__main__':
    main()
