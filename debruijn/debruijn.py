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

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
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
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


#==============================================================
# Main program
#==============================================================

def read_fastq(fastq):
    """Lit le fichier <fasta>.fasta passé en paramètre.

    Parametres
    ---------

    ------
    reference: str
        la séquence de référence
    """
    with open(fastq, "r") as filin:
        for line in filin:
            yield next(filin)[:-1]
            next(filin)
            next(filin)
    


def cut_kmer(sequence, kmer_size):
    """ fonction lisant une séquence est la découpant en une liste de k-mer

    Paramètres :
    ------------
    Une séquence
        issu d'un fichier fastq par la fonction read_fastq
    kmer_size
        la taille des kmers à générés
    ------------
    """
    kmer = ""
    kmerList = []
    seq_len = len(sequence)
    for i in range(seq_len-kmer_size+1):
        yield(sequence[i:i+kmer_size])

def build_kmer_dict(fastq, kmer_size):
    kmer_dict = {}
    for i in read_fastq(fastq):
        for j in cut_kmer(i, kmer_size):
            if j not in kmer_dict:
                kmer_dict[j] = 1
            else:
                kmer_dict[j] += 1
    return(kmer_dict)
            
def build_graph(kmer_dict):
    pref = []
    suff = []
    weight = []
    g = nx.DiGraph()
    for key in kmer_dict:
        pref.append(key[0:-1])
        suff.append(key[1:])
        weight.append(kmer_dict[key])
    nodes = list(zip(zip(pref, suff), weight))
    g = nx.DiGraph(nodes)
    return g
    

    



def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    kmer_size = args.kmer_size
    #fastq = read_fastq(args.fastq_file)
    kmer_dict = build_kmer_dict(args.fastq_file, kmer_size)
    graph = build_graph(kmer_dict)
    

    #cut_kmer(fastq[0], kmer_size)


if __name__ == '__main__':
    main()
