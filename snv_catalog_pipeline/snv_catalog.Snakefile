import os
import pandas as pd

BASE_DIR = "species/"


SPECIES = os.listdir(BASE_DIR)

df = pd.read_csv('ref_genome_chunk_files/chunkax', sep='\t')

print(df)
ref_genome_dict = df.set_index('species')['rep_genome'].to_dict()

SPECIES = ref_genome_dict.keys()
print(ref_genome_dict)

genome_dict = {}
species_long = []
files_long = []
for s in SPECIES:
    v = [os.path.splitext(f)[0] for f in os.listdir(os.path.join(BASE_DIR, s, 'genomes'))]
    genome_dict[s] = v
    species_long = species_long + [s]*len(v)
    files_long = files_long + v

print(files_long)

# Define all rule dependencies
rule all:
    input:
        expand(os.path.join(BASE_DIR, "{species}", 'prodigal', "{file}.faa"), zip, species=species_long, file=files_long),
        expand(os.path.join(BASE_DIR, "{species}/output/{species}_core_genome_mask.tsv"), species=SPECIES),
        expand(os.path.join(BASE_DIR, "{species}/output/{species}_4D_sites.tsv"), species=SPECIES),
        expand(os.path.join(BASE_DIR, "{species}/output/{species}.catalog.noAuto.wtRef.tsv"), species=SPECIES)


rule prodigal:
    input:
        os.path.join(BASE_DIR, "{species}", 'genomes', "{file}.fa")
    output:
        os.path.join(BASE_DIR, "{species}", 'prodigal', "{file}.faa")
    shell:
        "prodigal -i {input} -d {output} -q > /dev/null"


# Rule to run MMseqs2 clustering
rule mmseqs:
    input:
        expand(os.path.join(BASE_DIR, "{species}", 'prodigal', "{file}.faa"), zip, species=species_long, file=files_long)
    output:
        os.path.join(BASE_DIR, "{species}/mmseqs/{species}_cluster.tsv")
    params:
        i=os.path.join(BASE_DIR, "{species}/prodigal/*"),
        o=os.path.join(BASE_DIR, "{species}/mmseqs/{species}")
    shell:
        "mmseqs easy-cluster --min-seq-id 0.9 --threads 2 {params.i} {params.o} tmp"


# Rule to determine ORF clusters present in >90% of input genomes
rule determine_clusters:
    input:
        os.path.join(BASE_DIR, "{species}/mmseqs/{species}_cluster.tsv")
    output:
        os.path.join(BASE_DIR, "{species}/mmseqs/{species}_core_cluster.tsv")
    params:
        ref_genome=lambda wildcards, input: ref_genome_dict[wildcards.species],
        species=lambda wildcards, input: str(input).split('/')[1],
        genome_dir=lambda wildcards, input: os.path.join(BASE_DIR, str(input).split('/')[1], "genomes")
    shell:
        """
        python determine_clusters.py {input} {params.genome_dir} {params.species} {params.ref_genome}
        """

rule create_core_genome_mask:
    input:
        os.path.join(BASE_DIR, "{species}/mmseqs/{species}_core_cluster.tsv")
    output:
        os.path.join(BASE_DIR, "{species}/output/{species}_core_genome_mask.tsv")
    params:
        ref_genome_faa=lambda wildcards, input: os.path.join(BASE_DIR, str(input).split('/')[1], "prodigal", ref_genome_dict[wildcards.species] + '.faa')
    shell:
        """
        python filter_snv_catalog.py {params.ref_genome_faa} {input} {output}
        """

rule determine_4d_sites:
    input: 
        expand(os.path.join(BASE_DIR, "{species}", 'prodigal', "{file}.faa"), zip, species=species_long, file=files_long)
    output:
        os.path.join(BASE_DIR, "{species}/output/{species}_4D_sites.tsv")
    params:
        ref_genome_faa=lambda wildcards, input: os.path.join(BASE_DIR, wildcards.species, "prodigal", ref_genome_dict[wildcards.species] + '.faa')
    shell:
        """
        python fourfold.py {params.ref_genome_faa} {output}
        """

rule make_snv_catalog:
    input:
        expand(os.path.join(BASE_DIR, "{species}", 'prodigal', "{file}.faa"), zip, species=species_long, file=files_long)
    output:
        os.path.join(BASE_DIR, "{species}/output/{species}.catalog.noAuto.wtRef.tsv")
    params:
        ref_genome=lambda wildcards, input: ref_genome_dict[wildcards.species],
        species=lambda wildcards, input: wildcards.species,
    shell:
        """
        bash run_all_snv_catalog.sh {params.species} {params.ref_genome}
        """


