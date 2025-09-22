# README: Pipeline for Identifying and Analyzing TIR Transposons

## Overview
This pipeline identifies Terminal Inverted Repeat (TIR) transposons from genome assemblies, annotates them using MITE-Tracker, filters for transposase domains, groups candidates, filters active pairs based on similarity and flanking sequences, performs statistical summarization, and includes advanced clustering for TIR sequences both across and within species. It focuses on DNA transposons (Class II) while excluding retrotransposons. The process involves genome downloading, annotation, protein translation, domain searching with HMMER, grouping, filtering, final analysis, and sequence clustering.

The pipeline is designed for batch processing across multiple species (e.g., 1000+ species IDs). It assumes access to a Linux/Unix environment with sufficient computational resources (e.g., multi-core CPU for HMM searches and alignments).

**Key Features:**
- Identifies TIR elements with Target Site Duplications (TSDs) and transposase domains.
- Filters for active transposon pairs based on sequence similarity (>99% identity) and flanking sequence dissimilarity (<80% similarity).
- Groups transposons by domain, length, and direction.
- Outputs include FASTA files, domain profiles, grouped active transposons, summaries (e.g., TSD, TIR motifs, flanking sequences), and clustered TIR sequences.
- Advanced clustering for TIR sequences to identify similar families across or within species.

**Estimated Runtime:** Varies by genome size and number of species; Step 6 (pairwise comparisons) and Steps 10-11 (clustering) are computationally intensive and may need optimization (e.g., parallelization).

## Prerequisites
- **Software:**
  - Python 3.7+ with libraries: BioPython (`biopython`), Pandas (`pandas`), CSV.
  - HMMER3 (for `hmmsearch`, `hmmfetch`, `hmmpress`).
  - EMBOSS tools (for `transeq`).
  - VSEARCH (for pairwise global alignments).
  - Perl (for FASTA filtering script: `filter_fasta_by_ids.pl` – assume a standard Perl script for ID-based FASTA extraction).
  - AWK, Bash, and standard Unix tools (sort, uniq, grep, etc.).
  - Anaconda or similar for managing environments.

- **Databases/Models:**
  - Pfam-A HMM database (`Pfam-A.hmm`).
  - Custom HMM profiles for specific superfamilies (e.g., Ginger, Sola) – download from GitHub (link: xxx) or generate using HMMER.
  - List of transposase domains (e.g., `tnp_domain` file) and retrotransposon domains (`retro_domain` file).

- **Input Files:**
  - Species IDs list: `./data/species_ids.txt` (e.g., 1064species_id.txt).
  - Genome FASTA files: Downloaded to `./data/genomes/`.
  - Domain lists: `./data/domains/tnp_domain` and `./data/domains/retro_domain`.

Install dependencies:
```
conda install -c bioconda hmmer emboss vsearch biopython pandas
pip install biopython  # For clustering steps
```

## Directory Structure
Use the following structure for organization. Replace `/path/to/project/` with your actual base directory.

- `/path/to/project/data/genomes/` : Downloaded genome FASTA files (e.g., speciesid_genome.fa).
- `/path/to/project/tools/mitetracker/` : MITE-Tracker installation directory.
- `/path/to/project/results/mitetracker_results/` : MITE-Tracker output (candidates.csv).
- `/path/to/project/results/candidate_csv/` : Filtered CSV files (e.g., speciesid_withtsd.csv).
- `/path/to/project/results/candidate_fasta/` : Extracted TE FASTA files.
- `/path/to/project/results/candidate_tir/` : TIR sequences.
- `/path/to/project/results/candidate_tsd/` : TSD sequences.
- `/path/to/project/results/candidate_proteins/` : Translated protein sequences (.pep).
- `/path/to/project/results/domain_search/` : HMM search results (e.g., speciesid.tmp).
- `/path/to/project/results/extract_tnp/` : Extracted transposase-containing TEs.
- `/path/to/project/results/extract_retro/` : Extracted retrotransposons (for filtering).
- `/path/to/project/results/tnp_te_info/` : TE coordinates, lengths, TIR, TSD info.
- `/path/to/project/results/tnp_te_seq/` : Sequences of tnp-containing TEs.
- `/path/to/project/results/extract_domains/` : Domain and direction info.
- `/path/to/project/results/grouped_tes/` : Grouped TEs by domain/length.
- `/path/to/project/results/extract_grouped_fasta/` : Grouped FASTA sequences (direction-adjusted).
- `/path/to/project/results/pairwise_alignments/` : VSEARCH alignment results.
- `/path/to/project/results/flanking_filter/` : Filtered pairs based on flanking similarity.
- `/path/to/project/results/subgroups/` : Subgrouped active TEs.
- `/path/to/project/results/active_hmmsearch/` : Full Pfam HMM on active TEs.
- `/path/to/project/results/summaries/` : Final matrices, TSD/TIR extracts.
- `/path/to/project/results/tir_clustering/` : TIR clustering outputs (pairs and clusters).
- `/path/to/project/results/intra_species_clustering/` : Intra-species clustering outputs.

## Step-by-Step Instructions

### Step 1: Download Genomes
Download genome assemblies for species IDs listed in `./data/species_ids.txt`.

- Sources:
  - NCBI: https://www.ncbi.nlm.nih.gov/datasets/genome/
  - NGDC: https://ngdc.cncb.ac.cn/gwh/ncbi_assembly/

Save FASTA files as `./data/genomes/speciesid_genome.fa`.

### Step 2: Annotate TIR TEs with MITE-Tracker
Run MITE-Tracker to extract TIR TEs.

```bash
#!/bin/bash
fasta_dir="/path/to/project/data/genomes"
mitetracker_dir="/path/to/project/tools/mitetracker"
python3 ${mitetracker_dir}/MITETracker_simple.py -g ${fasta_dir}/speciesid_genome.fa -j speciesid --tsd_min_len 2 --tsd_max_len 15 --mite_min_len 2000 --mite_max_len 12000
```

Output: `./results/mitetracker_results/speciesid/candidates.csv` (columns: start,end,seq,record,mite_len,tir1_start,tir1_end,tir2_start,tir2_end,tir1_seq,tir2_seq,tsd,tsd_in,fs_left,fs_right,tir_len,candidate_id,description).

### Step 3: Organize CSV to FASTA and Translate to Proteins
Filter candidates.csv for TSD presence and no 'N' in sequences.

```bash
#!/bin/bash
candcsv_dir="/path/to/project/results/mitetracker_results"
tsdcsv_dir="/path/to/project/results/candidate_csv"
awk -F ',' '$13 == "no"' ${candcsv_dir}/speciesid/candidates.csv > ${tsdcsv_dir}/speciesid-candidates-tsd.csv
```

Batch filter for no 'N':

```bash
mapfile -t group_ids < /path/to/project/data/species_ids.txt
for group_id in "${group_ids[@]}"; do
    awk -F',' '$3 !~ /N/' ${tsdcsv_dir}/${group_id}-candidates-tsd.csv > ${tsdcsv_dir}/${group_id}_withtsd.csv
    rm ${tsdcsv_dir}/${group_id}-candidates-tsd.csv
done
```

Extract FASTA for TE, TIR, TSD; translate TE to protein.

```bash
#!/bin/bash
TE_fasta_dir="/path/to/project/results/candidate_fasta"
TE_TIR_dir="/path/to/project/results/candidate_tir"
TE_tsd_dir="/path/to/project/results/candidate_tsd"
cand_proteins="/path/to/project/results/candidate_proteins"

awk -F ',' -v OFS="\n" '{print ">"$17, $3}' ${tsdcsv_dir}/speciesid_withtsd.csv > ${TE_fasta_dir}/TE_speciesid.fa
awk 'BEGIN{Cnt=0}{if($0~/>/){Cnt=Cnt+1;tmp=">speciesid_TE"Cnt; $0=tmp; print$0}else{print$0}}' ${TE_fasta_dir}/TE_speciesid.fa > ${TE_fasta_dir}/TE_speciesid_newid.fa
mv ${TE_fasta_dir}/TE_speciesid_newid.fa ${TE_fasta_dir}/TE_speciesid.fa

# Similar for TIR1, TIR2, TSD (omitted for brevity; see original code)

transeq ${TE_fasta_dir}/TE_speciesid.fa ${cand_proteins}/TE_speciesid.pep -frame=6 -clean
```

### Step 4: HMM Search for Transposase Domains
Prepare HMM models for transposase superfamilies (see table below in query for domains like MULE, CACTA, etc.).

superfamily	Pfam	Domain name	Length of TSD (bp)	Note
MULE	PF10551	MULE	9-15	
CACTA	PF02992/PF03004/PF13963	Transposase_21/Transposase_24/Transpos_assoc	3	
hAT	PF05699	Dimer_Tnp_hAT	8	
PIF-Harbinger	PF13359	DDE_Tnp_4	3	
Tc1/Mariner	PF13358	DDE_3	2	
Merlin	PF12762	DDE_Tnp_IS1595	8	
Kolobok	PF20700	Mutator	4	
Zator	PF07592	DDE_Tnp_ISAZ013	3	
Zisupton	PF18758	KDZ	8	
Transib	PF12940	RAG1	5	
PiggyBac	PF13843	DDE_Tnp_1_7	4	
P	PF21787	TNP-like_RNaseH_N	8	
Ginger	-	ginger1, ginger2	3-5	domain profile generated by hmmer
Sola	-	sola1, sola2, sola3	4	domain profile generated by hmmer
Academ	-	Academ	6-10	domain profile generated by hmmer
Novosib	-	Novosib	8	domain profile generated by hmmer
Dada	-	DADA	6-7	domain profile generated by hmmer

```bash
hmmfetch --index /path/to/project/data/pfam/Pfam-A.hmm
hmmfetch -f /path/to/project/data/pfam/Pfam-A.hmm domain_names > /path/to/project/results/domain_search/TE-related.hmm
hmmpress /path/to/project/results/domain_search/TE-related.hmm
```

Run HMM search:

```bash
#!/bin/bash
hmmsearch --cpu 4 --domtblout /path/to/project/results/domain_search/speciesid.tmp /path/to/project/results/domain_search/TE-related.hmm /path/to/project/results/candidate_proteins/TE_speciesid.pep
```

Extract TNP-containing TEs:

```bash
#!/bin/bash
mapfile -t tnp_domains < /path/to/project/data/domains/tnp_domain
temp_record=()
while read -r line; do
    fields=($line)
    col4=${fields[3]}
    col1=${fields[0]}
    if [[ " ${tnp_domains[@]} " =~ " ${col4} " ]]; then
        temp_record+=("${col1:0:-2}")
    fi
done < /path/to/project/results/domain_search/speciesid.tmp
printf "%s\n" "${temp_record[@]}" | sort -u > /path/to/project/results/extract_tnp/speciesid_tnp_TE
```

Extract and filter retrotransposons similarly (using `retro_domain`), then subtract from TNP list:

```bash
#!/bin/bash
awk 'NR==FNR{a[$0]; next} !($0 in a)' /path/to/project/results/extract_retro/speciesid_retro_TE /path/to/project/results/extract_tnp/speciesid_tnp_TE > /path/to/project/results/extract_tnp/speciesid_tnp_filtered
```

### Step 5: Extract TNP TE Info and Group
Extract coordinates, lengths, etc.:

```python
# Python script: template_filter.py
import csv

vec = []
teid_data = []
with open('/path/to/project/results/extract_tnp/speciesid_tnp_filtered', 'r') as file:
    for line in file:
        teid_data.append(line.strip())
        if "TE" in line:
            num = line.split("TE")[1].strip()
            vec.append(int(num))

with open('/path/to/project/results/candidate_csv/speciesid_withtsd.csv', 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rows = list(reader)

with open('/path/to/project/results/tnp_te_info/speciesid_coord_len_tirlen_tsd_tsdlen', 'w', newline='') as outfile:
    writer = csv.writer(outfile, delimiter='\t')
    for idx, teid in zip(vec, teid_data):
        row = rows[idx - 1]
        selected_columns = [row[3], row[0], row[1], row[4], row[9], row[10], row[13], row[14], row[15], row[11]]
        length_col12 = len(row[11])
        output_row = [teid] + selected_columns + [length_col12]
        writer.writerow(output_row)
```

Extract domains and directions (see `extract_domain_template.py`, `uniq_combine.sh`, `domain_direction.py` in query).

Combine info:

```bash
awk -F'\t' -v OFS='\t' '{...}' /path/to/project/results/extract_domains/speciesid_domain_direction /path/to/project/results/tnp_te_info/speciesid_coord_len_tirlen_tsd_tsdlen > /path/to/project/results/tnp_te_info/speciesid_coord_len_tirlen_tsd_tsdlen_domain_direction
```

Group by domain, direction, TE/TIR lengths (see `group_template.py` in query). Filter groups with >=2 sequences.

### Step 6: Filter and Split Groups (Pairwise Comparison)
Extract grouped FASTA (direction-adjusted for '-' strands using BioPython).

```python
# retrieve_seq.py
from Bio import SeqIO
from Bio.Seq import Seq

input_file1 = "/path/to/project/results/grouped_tes/groupid"
filename1 = "groupid"
prefix = filename1.split("_group")[0]
input_file2 = f"/path/to/project/results/tnp_te_seq/{prefix}_tnp.fasta"
output_file = "/path/to/project/results/extract_grouped_fasta/groupid.fasta"

sequences_info = {}
with open(input_file1, 'r') as f:
    for line in f:
        parts = line.strip().split()
        seq_id = parts[0]
        direction = parts[9]
        sequences_info[seq_id] = direction

with open(output_file, 'w') as out_f:
    for record in SeqIO.parse(input_file2, "fasta"):
        if record.id in sequences_info:
            direction = sequences_info[record.id]
            if direction == "+":
                SeqIO.write(record, out_f, "fasta")
            else:
                record.seq = record.seq.reverse_complement()
                SeqIO.write(record, out_f, "fasta")
```

Run VSEARCH for pairwise alignments (>99% ID):

```bash
vsearch --allpairs_global /path/to/project/results/extract_grouped_fasta/groupid.fasta --minsl 0.98 --threads 20 --uc /path/to/project/results/pairwise_alignments/groupid.aln --leftjust --rightjust --idprefix 3 --idsuffix 3 --id 0.99 --iddef 1
```

Filter pairs and collect (see query for AWK and file merging).

### Step 7: Filter by Flanking Sequence Similarity
Filter pairs with flanking seq similarity <80% (using BioPython pairwise2).

```python
# filter.py
import csv
import os
from Bio.Seq import Seq
from Bio import pairwise2

def process_files(input_file, output_dir, coord_dir, tsd_dir):
    filename = os.path.basename(input_file)
    speciesid = filename.split("_group_")[0]
    with open(input_file, 'r') as file1:
        lines = file1.readlines()
    output_filename = os.path.join(output_dir, filename + "_filtered")
    with open(output_filename, 'w') as outfile:
        for line in lines:
            parts = line.strip().split()
            element1, element2 = parts[-2], parts[-1]
            row1, row2 = element1.split("TE")[-1], element2.split("TE")[-1]
            coord_file = os.path.join(coord_dir, f"{speciesid}_coord_len_tirlen_tsd_tsdlen_domain_direction")
            tsd_file = os.path.join(tsd_dir, f"{speciesid}_withtsd.csv")
            # Get directions and sequences (omitted for brevity; see query)
            # Align and filter if simil < 80

# Call with appropriate paths
```

Subgroup connected components (graph-based, see `sub_temp.py` in query).

### Step 8: Extract Active TE Sequences and Full Pfam HMM
Combine grouped FASTA by species (see `combineseq2run.sh`).

Extract and translate active TEs, run full Pfam HMM:

```bash
perl filter_fasta_by_ids.pl --id_list /path/to/project/results/subgroups/seq_id/speciesid.txt --fasta_file /path/to/project/results/extract_grouped_fasta/speciesid_diradj.fa --output_file /path/to/project/results/subgroups/seq_data/speciesid_active.fa
transeq /path/to/project/results/subgroups/seq_data/speciesid_active.fa /path/to/project/results/subgroups/seq_data/speciesid_active.pep -frame=6 -clean
hmmsearch --cpu 16 --domtblout /path/to/project/results/active_hmmsearch/speciesid_active.tmp /path/to/project/data/pfam/Pfam-A.hmm /path/to/project/results/subgroups/seq_data/speciesid_active.pep
```

### Step 9: Statistical Summary
Generate TE family-host matrix (custom script based on groups).

Extract TSD and TIR (15bp head) for type judgment (see `extract_tsd.py` and `extract_tir.py` in query).

Filter groups further by TSD length consistency to reduce false positives.

### Step 10: TIR Sequence Clustering
This step clusters Terminal Inverted Repeat (TIR) sequences based on similarity and length criteria. It processes input data from a tab-separated file, identifies similar sequence pairs with orientation (+ for same direction, - for reverse complement), and groups them into clusters.

The process is divided into two main stages:
1. Pairwise comparison and pair generation.
2. Cluster formation from pairs.

#### Features
- Groups sequences by species ID and TIR type.
- Computes sequence similarity using global alignment (via BioPython's pairwise2).
- Handles reverse complements for orientation detection.
- Builds connected components (clusters) from pairs using a union-find-like approach.

#### Input File Format
The input file (`auto_candidates_info.txt`) should be a tab-separated text file with at least the following columns (0-indexed):
- Column 0: Element ID
- Column 1: Species ID
- Column 4: TIR type
- Column 7: TIR length (integer)
- Column 8: Forward sequence
- Column 9: Reverse sequence

Example line: (Note: In our analysis, Species ID was designated using numerical identifiers.)
```
id1	1000	3	GAGCATCTCC	Pif/Harbinger	1	group_id1	59	GAGCATCTCCAGCAGATGCTTTAAAAGTGTCTCGTGCTTTAAAATCTGAGAATATAAAG	CTTTATATTCTCAGATTTTAAAGCACGAGACACTTTTAAAGCATCTGCTGGAGATGCTC
id1	1000	3	GAGCATCTCC	Pif/Harbinger	1	group_id2	25	GAGCATCTCCAACAATCCCCTTAAA	TTTAAGGGGATTGTTGGAGATGCTC
id3	1000	3	GAGCATCTCC	Pif/Harbinger	2	group_id3	20	GAGCATCTCCAACAGATTAG	CTAATCTGTTGGAGATGCTC
id4	1000	3	CACTACACCA	CACTA	1	group_id4	21	CACTACACCACAACACCAAAA	TTTTGGTGTTGTGGTGTAGTG
id5	1000	3	CACTACTACA	CACTA	1	group_id5	21	CACTACTACAAAACATGGATA	TATCCATGTTTTGTAGTAGTG
id6	1000	8	TAGGGGTGTA	Hat	1	group_id6	23	TAGGGGTGTAAACGGGGAGCGGATA	TATCCGCTCCGTTTACACCCCTA
```

#### Stage 1: Generate Pairs
This stage compares sequences within groups and outputs pairs with similarity > 0.9 on the first/last 70 bases, considering length differences < 10%.

```python
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

# Run the function
process_file('auto_candidates_info.txt', 'auto_cluster.txt')
```

#### Stage 2: Form Clusters
This stage reads the pairs and groups them into clusters.

```python
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
        elif cluster1 is None and cluster2 is None:
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

# Run the function
process_file('auto_cluster.txt', 'auto_cluster2.txt')
```

#### Output Files
- `auto_cluster.txt`: Comma-separated pairs with orientation (e.g., `element1,element2,+`).
- `auto_cluster2.txt`: Elements mapped to clusters (e.g., `element1,cluster1`).

Run both stages sequentially in `/path/to/project/results/tir_clustering/`.

### Step 11: Intra-Species TIR Sequence Clustering
This step clusters Terminal Inverted Repeat (TIR) sequences within a species based on sequence similarity and length criteria. It reads input data from a tab-separated file containing TIR sequences, identifies similar sequence pairs with orientation (+ for same direction, - for reverse complement), and groups them into clusters.

The process is divided into two main stages:
1. Pairwise comparison and pair generation for intra-species sequences.
2. Cluster formation from the identified pairs.

#### Features
- Compares TIR left (L) and right (R) sequences using global alignment (via BioPython's pairwise2).
- Checks for length differences less than 10% before alignment.
- Handles reverse complements for orientation detection.
- Builds connected components (clusters) from pairs using a union-find-like approach.

#### Input File Format
The input file (`cluster_tir.txt`) should be a tab-separated text file with the following columns (0-indexed):
- Column 0: Element ID
- Column 1: TIR left sequence (TIR_L)
- Column 2: TIR right sequence (TIR_R)

Example line:
```
cluster1     CACTAGTAAAAATCTATCGCATAGCCACT   AGTGGCTATGCGTTGGATTTTTACTAGTG
cluster2      CACTACAAGAAAACAGCGGTATTCTGACGGA TCCCTCAGTATCCCGATGTTTTCTTGTAGTG
cluster3      CACTACAAGAAAACGTAATTTTAACGAGG   CCTCGTAAAACTTGTGTTTTCTTGTAGTG
cluster4      CACTACAAGAAAACGTAACTTTAACGAGGAAGTTT     AAACTGCCCTCGTAAAACTTGTGTTTTCTTGTAGTG
cluster5      CACTACAGGAAATAAGGGTTTTGATGGCA   TGCTATTAAAACCCTTATTTCCTGTAGTG
cluster6     CACTACAAGAAAAAGAGGTTTTTATAGCGG  CCGCTAAGAATACCCTTTTTTCTTGTAGTG
```

#### Stage 1: Generate Pairs
This stage compares all pairs of sequences, aligns the first 70 bases of TIR_L and last 70 bases of TIR_R (or their reverse complements), and outputs pairs with identity >= 0.9, considering length differences < 10%.

```python
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
```

#### Stage 2: Form Clusters
This stage reads the pairs and groups them into clusters.

```python
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

# Run the function
process_file('auto_cluster_intra4species.txt', 'auto_cluster2_intra4species.txt')
```

#### Output Files
- `auto_cluster_intra4species.txt`: Comma-separated pairs with orientation (e.g., `element1,element2,+`).
- `auto_cluster2_intra4species.txt`: Elements mapped to clusters (e.g., `element1,cluster1`).

Run both stages sequentially in `/path/to/project/results/intra_species_clustering/`.

## Troubleshooting
- Slow steps: Parallelize loops over species IDs or use multi-threading for alignments in clustering.
- Memory issues: Split large FASTA or input files.
- Custom scripts: Ensure paths are updated; test on small subsets.

## Contributing
Contributions are welcome! Please open an issue or submit a pull request for improvements, bug fixes, or additional features.

## License
This project is licensed under the MIT License.

For questions, contact zhangxinyan@caas.cn
