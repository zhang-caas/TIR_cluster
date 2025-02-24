from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import substitution_matrices
from collections import defaultdict

with open('cluster_tir.txt', 'r') as file:
    lines = file.readlines()

data = [line.strip().split('\t') for line in lines]
ids = [line[0] for line in data]
tir_ls = [line[1] for line in data]
tir_rs = [line[2] for line in data]

def align_sequences(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    best_alignment = alignments[0]
    identity = best_alignment[2] / best_alignment[4]
    return identity

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

with open('auto_cluster_intra4species.txt', 'w') as outfile:
    for i in range(len(data)):
        for j in range(i + 1, len(data)):
            id1, tir_l1, tir_r1 = ids[i], tir_ls[i], tir_rs[i]
            id2, tir_l2, tir_r2 = ids[j], tir_ls[j], tir_rs[j]
            if abs(len(tir_l1) - len(tir_l2)) >= 0.1 * min(len(tir_l1), len(tir_l2)):
                continue
            identity_l = align_sequences(tir_l1[:70], tir_l2[:70])
            identity_r = align_sequences(tir_r1[-70:], tir_r2[-70:])

            if identity_l >= 0.9 and identity_r >= 0.9:
                outfile.write(f"{id1},{id2},+\n")
            else:
                identity_l_rc = align_sequences(tir_l1[:70], reverse_complement(tir_r2)[:70])
                identity_r_rc = align_sequences(tir_r1[-70:], reverse_complement(tir_l2)[-70:])

                if identity_l_rc >= 0.9 and identity_r_rc >= 0.9:
                    outfile.write(f"{id1},{id2},-\n")

def find_clusters(pairs):
    clusters = []
    element_to_cluster = {}

    def find_cluster(element):
        if element in element_to_cluster:
            return element_to_cluster[element]
        return None

    for element1, element2 in pairs:
        cluster1 = find_cluster(element1)
        cluster2 = find_cluster(element2)

        if cluster1 is None and cluster2 is None:
            new_cluster = {element1, element2}
            clusters.append(new_cluster)
            element_to_cluster[element1] = new_cluster
            element_to_cluster[element2] = new_cluster
        elif cluster1 is not None and cluster2 is None:
            cluster1.add(element2)
            element_to_cluster[element2] = cluster1
        elif cluster1 is None and cluster2 is not None:
            cluster2.add(element1)
            element_to_cluster[element1] = cluster2
        elif cluster1 is not cluster2:
            cluster1.update(cluster2)
            for element in cluster2:
                element_to_cluster[element] = cluster1
            clusters.remove(cluster2)

    return clusters

def process_file(input_file, output_file):
    pairs = []
    with open(input_file, 'r') as infile:
        for line in infile:
            elements = line.strip().split(',')
            pairs.append((elements[0], elements[1]))

    clusters = find_clusters(pairs)

    with open(output_file, 'w') as outfile:
        cluster_id = 1
        for cluster in clusters:
            for element in cluster:
                outfile.write(f"{element},cluster{cluster_id}\n")
            cluster_id += 1

process_file('auto_cluster_intra4species.txt', 'auto_cluster2_intra4species.txt')