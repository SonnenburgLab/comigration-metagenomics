import re
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


def generate_degeneracy_dict():
    degeneracy_dict = {}
    for base1 in 'ACGTN':
        for base2 in 'ACGTN':
            for base3 in 'ACGTN':
                codon = base1 + base2 + base3
                if 'N' in codon:
                    degeneracy_dict[codon] = 0
                elif codon in ['TAA', 'TAG', 'TGA']:  # Stop codons
                    degeneracy_dict[codon] = 1
                elif codon in ['ATG', 'CTG', 'TTG', 'TGG']: # Start codons / Trp
                    degeneracy_dict[codon] = 1
                elif re.match(r'^(GC|CG|GG|CT|CC|TC|AC|GT)[ACTG]', codon): # 4-fold
                    degeneracy_dict[codon] = 4
                elif re.match(r'^(AG)[AG]', codon): # Arg
                    degeneracy_dict[codon] = 2
                elif re.match(r'^(AA)[TC]', codon): # Asn
                    degeneracy_dict[codon] = 2
                elif re.match(r'^(GA)[TC]', codon): # Asp
                    degeneracy_dict[codon] = 2
                elif re.match(r'^(TG)[TC]', codon): # Cys
                    degeneracy_dict[codon] = 2
                elif re.match(r'^(CA)[AG]', codon): # Glu
                    degeneracy_dict[codon] = 2
                elif re.match(r'^(GA)[AG]', codon): # Glu
                    degeneracy_dict[codon] = 2
                elif re.match(r'^(CA)[TC]', codon): # His
                    degeneracy_dict[codon] = 2
                elif re.match(r'^(AT)[TCA]', codon): # Ile
                    degeneracy_dict[codon] = 3
                elif re.match(r'^(TT)[AG]', codon): # Leu
                    degeneracy_dict[codon] = 2
                elif re.match(r'^(AA)[AG]', codon): # Lys
                    degeneracy_dict[codon] = 2
                elif re.match(r'^(TT)[TC]', codon): # Phe
                    degeneracy_dict[codon] = 2
                elif re.match(r'^(AG)[TC]', codon): # Ser
                    degeneracy_dict[codon] = 2
                elif re.match(r'^(TA)[TC]', codon): # Tyr
                    degeneracy_dict[codon] = 2             
    return degeneracy_dict

def complement_base(base):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return complement.get(base, base)

def read_fasta(filename):
    sequences = {}
    for record in SeqIO.parse(filename, "fasta"):
        sequences[record.description] = str(record.seq)
    return sequences

def find_fourfold_degenerate_sites(fasta_file):
    sequences = read_fasta(fasta_file)
    fourfold_sites = []
    dc = generate_degeneracy_dict()
    for sequence_id, sequence in sequences.items():
        protein_id = sequence_id.split(' # ')[0]
        contig = protein_id.rsplit('_', 1)[0]

        position_start = sequence_id.split(' # ')[1]
        position_end = sequence_id.split(' # ')[2]
        direction = sequence_id.split(' # ')[3]

        codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
        for i, codon in enumerate(codons):
            if len(codon) == 3:
                degeneracy = dc[codon]
                if direction == '1':
                    site_position = int(position_start)+i*3+2
                    base = codon[2]
                elif direction == '-1':
                    site_position = int(position_end)-i*3-2
                    base = complement_base(codon[2])
                fourfold_sites.append([contig, str(site_position), str(degeneracy), base])

    return fourfold_sites


if __name__ == '__main__':
    # Example usage:
    fasta_file = sys.argv[1]
    output_file = sys.argv[2]
    fourfold_sites = find_fourfold_degenerate_sites(fasta_file)
    df = pd.DataFrame(fourfold_sites, columns=['Contig', 'Site', 'Degeneracy', 'Base'])
    df_no_duplicates = df.drop_duplicates(subset=['Contig', 'Site'])
    df_no_duplicates.to_csv(output_file, sep='\t', index=False)

    # print("Fourfold degenerate sites:", fourfold_sites)
    # print('\n'.join(fourfold_sites))

