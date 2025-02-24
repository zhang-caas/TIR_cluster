from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
from collections import defaultdict

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def sequence_similarity(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    max_score = max(alignments, key=lambda x: x[2])[2]
    return max_score / min(len(seq1), len(seq2))

def process_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        lines = infile.readlines()
        data = [line.strip().split('\t') for line in lines]
        
        groups = {}
        for line in data:
            speciesid, tirtype = line[1], line[4]
            key = (speciesid, tirtype)
            if key not in groups:
                groups[key] = []
            groups[key].append(line)
        
        for key, group in groups.items():
            if len(group) > 1:
                for i in range(len(group)):
                    for j in range(i + 1, len(group)):
                        line1, line2 = group[i], group[j]
                        tir_len1, tir_len2 = int(line1[7]), int(line2[7])
                        if abs(tir_len1 - tir_len2) >= 0.1 * min(tir_len1, tir_len2):
                            continue
                        if sequence_similarity(line1[8][:70], line2[8][:70]) > 0.9 and sequence_similarity(line1[9][-70:], line2[9][-70:]) > 0.9:
                            outfile.write(f"{line1[0]},{line2[0]},+\n")
                        elif sequence_similarity(line1[8][:70], reverse_complement(line2[9])[:70]) > 0.9 and sequence_similarity(line1[9][-70:], reverse_complement(line2[8])[-70:]) > 0.9:
                            outfile.write(f"{line1[0]},{line2[0]},-\n")

process_file('auto_candidates_info.txt', 'auto_cluster.txt')


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

process_file('auto_cluster.txt', 'auto_cluster2.txt')