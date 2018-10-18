#!/usr/bin/env python
import argparse, pdb

from kcollections import Kset

from Bio import SeqIO

KMER_SIZE = 19
BASE_PAIRS = ['A', 'C', 'G', 'T']


def get_kmers(seq, k=KMER_SIZE):
    for i in range(len(seq) - k + 1):
        kmer = str(seq[i:i + k])
        yield kmer


def build_de_bruijn_graph(reads, k=KMER_SIZE):
    graph = Kset(k)
    for read_file in reads:
        with open(read_file, 'rU') as read_fh:
            for read in SeqIO.parse(read_fh, 'fastq'):
                for kmer in get_kmers(read.seq):
                    graph.add(kmer.upper())
    return graph


def get_neighbors(graph, node, f_concat):
    node = node.upper()
    neighbors = []
    if node not in graph:
        return neighbors
    for neighbor in [f_concat(node, bp) for bp in BASE_PAIRS]:
        if neighbor in graph:
            neighbors += [neighbor]
    return neighbors


def get_prefix_neighbors(graph, node):
    prepend = lambda node, bp: bp + node[:-1]
    return get_neighbors(graph, node, prepend)


def get_suffix_neighbors(graph, node):
    append = lambda node, bp: node[1:] + bp
    return get_neighbors(graph, node, append)


def assemble_reads(graph):
    visited = Kset(KMER_SIZE)
    contig_num = 0
    for node in graph:
        suffix_contig = create_contig_suffix(graph, visited, node, node)
        # visited.remove(node)
        # prefix_contig = create_contig_prefix(graph, visited, node, node)
        if len(suffix_contig) > KMER_SIZE:
            print('> contig {}\n{}'.format(contig_num, suffix_contig))
            contig_num += 1


def create_contig(graph, visited, cur_node, contig, f_neighbors, f_concat):
    if cur_node in visited:
        return contig
    neighbors = f_neighbors(graph, cur_node)
    visited.add(cur_node)
    if len(neighbors) == 1:
        return create_contig(graph, visited, neighbors[0], f_concat(contig, neighbors[0]), f_neighbors, f_concat)
    return contig


def create_contig_suffix(graph, visited, cur_node, contig):
    append_suffix = lambda contig, neighbor: contig + neighbor[-1]
    return create_contig(graph, visited, cur_node, contig, get_suffix_neighbors, append_suffix)


def create_contig_prefix(graph, visited, cur_node, contig):
    prepend_prefix = lambda contig, neighbor: neighbor[0] + contig
    return create_contig(graph, visited, cur_node, contig, get_prefix_neighbors, prepend_prefix)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Assemble reads using kcollections.')
    parser.add_argument('reads', type=str, nargs='+', default='./reads.fq',
                        help='One or more paths to reads in FASTQ format')

    args = parser.parse_args()

    graph = build_de_bruijn_graph(args.reads)

    assemble_reads(graph)
