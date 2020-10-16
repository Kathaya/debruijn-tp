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
#from operator import itemgetter
import statistics
#import matplotlib
#import matplotlib.pyplot as plt
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
    parser = argparse.ArgumentParser(description=__doc__, usage="{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


# ==============================================================
# Main program
# ==============================================================

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
    for i in range(len(sequence)-kmer_size+1):
        yield sequence[i:i+kmer_size]


def build_kmer_dict(fastq, kmer_size):
    """
    Fonction permetant de constuireun dictionnaire des kmer
    """
    kmer_dict = {}
    for i in read_fastq(fastq):
        for j in cut_kmer(i, kmer_size):
            if j not in kmer_dict:
                kmer_dict[j] = 1
            else:
                kmer_dict[j] += 1
    return kmer_dict


def build_graph(kmer_dict):
    """
    Construction d'un graphe orienté et pondéré à partir du dictionnaire de kmers
    """
    gdi = nx.DiGraph()
    for kmer, poids in kmer_dict.items():
        gdi.add_edge(kmer[:-1], kmer[1:], weight=poids)
    return gdi


def get_starting_nodes(graph):
    """
    Fonction renvoyant l'ensemble des noeuds d'entrée en liste
    """
    start_nodes = []
    for node in graph.nodes:
        if list(graph.predecessors(node)) == []:
            start_nodes.append(node)
    return start_nodes


def get_sink_nodes(graph):
    """
    Fonction renvoyant l'ensemble des noeuds de sortie en liste
    """
    end_nodes = []
    for node in graph.nodes:
        if list(graph.successors(node)) == []:
            end_nodes.append(node)
    return end_nodes

def get_contigs(graph, start, end):
    """
    Création des contigs à partir des noeuds d'entrée et de sortie non optimisé
    """
    contigs = []
    for debut in start:
        for fin in end:
            for path in nx.all_simple_paths(graph, debut, fin):
                cont = ""
                for i, idx in enumerate(path):
                    if i == 0:
                        cont = idx
                    else:
                        cont = cont + idx[-1]
                tup = [cont, len(cont)]
                contigs.append(tup)
    return contigs


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(contigs, file):
    """
    Fonction sauvegardant dans un fichier les contigs calculés
    """
    with open(file, "w") as filout:
        for i, cont in enumerate(contigs):
            filout.write(">contig_{} len={}\n".format(i, cont[1]))
            filout.write(fill(cont[0]+"\n"))


def std(vec_val):
    """
    Fonction calculant l'écart-type d'un vecteur
    """
    return statistics.stdev(vec_val)


def path_average_weight(graph, path):
    """
    Calcul du poids moyen d'un chemin
    """
    weigth = []
    for i in range(len(path)-1):
        weigth.append(graph.edges[path[i], path[i+1]]['weight'])
    return statistics.mean(weigth)


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """
    Supprimes tout les chemins indésirables donnés
    """
    for path in path_list:
        graph.remove_nodes_from(path[1:-1])
        if delete_entry_node:
            graph.remove_node(path[0])
        if delete_sink_node:
            graph.remove_node(path[-1])
    return graph


def select_best_path(graph, path_list, long_path_list, pm_path_list,
                     delete_entry_node=False, delete_sink_node=False):
    """
    Choisi le meilleur chemin en prennant d'abord celui avec le meilleur poids,
    Sinon deux chemins ont le même poids alors celui le plus long
    Sinon égalité alors en choisi un aléatoire.
    Puis supprime tout les path inutiles
    """
    ind_path_list = []
    # Sélection des chemin avec le meilleur poidssi variabilité
    if std(pm_path_list) != 0:
        indices = [i for i, x in enumerate(
            pm_path_list) if x == max(pm_path_list)]
        cp_path_list = list(path_list)
        pop_count = 0
        for i, path in enumerate(cp_path_list):
            if i not in indices:
                ind_path_list.append(path)
                path_list.pop(i-pop_count)
                long_path_list.pop(i-pop_count)
                pm_path_list.pop(i-pop_count)
                pop_count += 1
        if len(indices) == 1:
            return remove_paths(graph, ind_path_list, delete_entry_node, delete_sink_node)

    # Sélection des chemins avec la meilleur longeur
    if std(long_path_list) != 0:
        indices = [i for i, x in enumerate(
            long_path_list) if x == max(long_path_list)]
        cp_path_list = list(path_list)
        pop_count = 0
        for i, path in enumerate(cp_path_list):
            if i not in indices:
                ind_path_list.append(path)
                path_list.pop(i-pop_count)
                long_path_list.pop(i-pop_count)
                pm_path_list.pop(i-pop_count)
                pop_count += 1
        if len(indices) == 1:
            return remove_paths(graph, ind_path_list, delete_entry_node, delete_sink_node)

    # Sélection aléatoire du meilleur chemin
    indice = randint(0, len(path_list))
    for i, path in enumerate(path_list):
        if i not in indice:
            ind_path_list.append(path)

    return remove_paths(graph, ind_path_list, delete_entry_node, delete_sink_node)


def solve_bubble(graph, ancestor_node, descendant_node, ):
    path_list = []
    long_path_list = []
    pm_path_list = []
    for path in nx.all_simple_paths(graph, ancestor_node, descendant_node):
        path_list.append(path)
        pm_path_list.append(path_average_weight(graph, path))
        long_path_list.append(len(path))
    return select_best_path(graph, path_list, long_path_list, pm_path_list)


def simplify_bubbles(graph):
    pass


def solve_entry_tips():
    pass


def solve_out_tips():
    pass


def get_sink_nodes():
    pass


def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    kmer_size = args.kmer_size
    kmer_dict = build_kmer_dict(args.fastq_file, kmer_size)
    graph = build_graph(kmer_dict)
    # plt.subplot(121)
    #nx.draw(graph, with_labels=True, font_weight='bold')
    # nx.draw(graph)
    # plt.show()
    # print(list(graph.nodes()))
    # for i in list(graph.nodes()):
    # print(graph.predecessors(i))
    start = get_starting_nodes(graph)
    end = get_sink_nodes(graph)
    contigs = get_contigs(graph, start, end)
    save_contigs(contigs, "bravoloto.txt")


if __name__ == '__main__':
    main()
