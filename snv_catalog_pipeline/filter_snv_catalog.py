import os
import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq


def read_fasta(filename):
    sequences = {}
    for record in SeqIO.parse(filename, "fasta"):
        sequences[record.description] = str(record.seq)
    return sequences

def get_core_sites(fasta_records, core_gene_list):
    contigs = []
    coords = []

    for sequence_desc, sequence in fasta_records.items():
        gene_name = sequence_desc.split(' # ')[0]
        if gene_name in core_gene_list:
            contig = gene_name.rsplit('_', 1)[0]

            position_start = sequence_desc.split(' # ')[1]
            position_end = sequence_desc.split(' # ')[2]
            direction = sequence_desc.split(' # ')[3]

            if direction == '1':
                coord_lst = list(range(int(position_start), int(position_end)))
            elif direction == '-1':
                coord_lst = list(range(int(position_start)+1, int(position_end)+1))

            contig_lst = [contig]*len(coord_lst)

            contigs += contig_lst
            coords += coord_lst

    return (contigs, coords)



if __name__ == '__main__':
    core_gene_file = sys.argv[2]
    prodigal_file = sys.argv[1]
    output = sys.argv[3]

    with open(core_gene_file, 'r') as f:
        core_gene_list = [line.strip('\n') for line in f.readlines()]

    fasta_records = read_fasta(prodigal_file)
    (contigs, coords) = get_core_sites(fasta_records, core_gene_list)

    zipped_contigs_coords = zip(contigs, coords)

    with open(output, 'w', newline='\n') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        for contig, coord in zipped_contigs_coords:
            writer.writerow([contig, coord])

