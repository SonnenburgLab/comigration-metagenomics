from utils import snv_utils

data_batch = '240326_with_industrial'

metadata = snv_utils.load_dadi_snps_metadata(data_batch)
threshold = 10

# TODO: eventually will have different population names
mask = (metadata['Num Hadza'] > threshold) & (metadata['Num Tsimane'] > threshold) & (metadata['Num HMP'] > threshold) & (metadata['Num MetaHIT'] > threshold)
metadata_uni = metadata[mask]

print(f"Found {metadata_uni.shape[0]} species with at least {threshold} samples in each population")

species_list = metadata_uni.index.tolist()
for species in species_list:
    print(f"Converting {species} to VCF format")
    snv_utils.convert_snv_catalog_to_vcf(species, data_batch)
    print(f"VCF file written to {species}.vcf")