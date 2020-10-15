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
import random
from random import randint
from operator import itemgetter
import statistics
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
random.seed(9001)

__author__ = "Ferdinand Petit"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Ferdinand Petit"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Ferdinand Petit"
__email__ = "ferdinandpetit92@gmail.com"
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
    for kmer, poids in kmer_dict.items():
        g.add_edge(kmer[:-1], kmer[1:], weight=poids)
    return g
    
def get_starting_nodes(graph):
    x = list(graph.nodes())
    start_nodes = []
    for i in range(len(x)):
        if list(graph.predecessors(x[i])) == []:
            start_nodes.append(x[i])
    return start_nodes

def get_ending_nodes(graph):
    end_nodes = []
    x = list(graph.nodes())
    for i in range(len(x)):
        if list(graph.successors(x[i])) == []:
            end_nodes.append(x[i])
    return end_nodes

def get_contigs(graph, start, end):
    contigs = []
    for s in start: 
        for e in end:
            if(list(nx.all_simple_paths(graph, s, e))) != []:
                x = list(nx.all_simple_paths(graph, s, e))
                for idx in x:
                    cont = ""
                    for i, seq in enumerate(idx):
                        if i == 1:
                            cont = seq
                        else:
                            cont = cont + seq[-1]
                    tup = [cont, len(cont)]
                    contigs.append(tup)
    return contigs

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(contigs, file):
    with open(file, "w") as filout:
        for i, cont in enumerate(contigs):
            filout.write(">contig_{} len={}\n".format(i, cont[1]))
            filout.write(fill(cont[0]+"\n"))



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
    
    #plt.subplot(121)
    #nx.draw(graph, with_labels=True, font_weight='bold')
    #nx.draw(graph)
    #plt.show()
    
    #print(list(graph.nodes()))
    #for i in list(graph.nodes()):
    #    print(graph.predecessors(i))

    start = get_starting_nodes(graph)
    end = get_ending_nodes(graph)
    contigs = get_contigs(graph, start, end)
    save_contigs(contigs, "bravoloto.txt")




if __name__ == '__main__':
    main()
