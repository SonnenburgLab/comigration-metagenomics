#!/bin/bash

#SPECIES="Cryptobacteroides_sp000432655"
#COMPARATOR_GENOME="Cryptobacteroides_sp000432655__TS_ADULT_100"
SPECIES=$1
COMPARATOR_GENOME=$2

SPECIES_DIR=species/$SPECIES
GENOMES_DIR=species/$SPECIES/genomes
PAIRS_IN=$SPECIES_DIR/${SPECIES}_input_pairs.tsv

OUTPUT_DIR=$SPECIES_DIR/snv_catalog/
mkdir -p $OUTPUT_DIR

while IFS= read -r pair
do
   genome1=$(echo $pair | cut -d' ' -f1)
   genome2=$(echo $pair | cut -d' ' -f2)

   GB_PATH_1=$GENOMES_DIR/$genome1.fa
   GB_PATH_2=$GENOMES_DIR/$genome2.fa


   nucmer $GB_PATH_1 $GB_PATH_2 -p $OUTPUT_DIR/$genome1-$genome2 --threads 8
   delta-filter -q -r $OUTPUT_DIR/$genome1-$genome2.delta > $OUTPUT_DIR/$genome1-$genome2.filter.delta
   show-coords $OUTPUT_DIR/$genome1-$genome2.filter.delta > $OUTPUT_DIR/$genome1-$genome2.coords
   show-snps $OUTPUT_DIR/$genome1-$genome2.filter.delta > $OUTPUT_DIR/$genome1-$genome2.snps
   show-diff $OUTPUT_DIR/$genome1-$genome2.filter.delta > $OUTPUT_DIR/$genome1-$genome2.diff
   echo "done whole genome alignment for $genome1 and $genome2"
done < $PAIRS_IN

find $OUTPUT_DIR -name "*.snps" | xargs basename -a | sed -e 's/\.snps$//' | xargs -I[] bash -c 'sed "1,5d" '$OUTPUT_DIR'/[].snps | awk "$0"' '$2 != "." && $3 != "." {printf "%s\t%s||%s||%s||%s\n", "[]", $14, $1, $2, $3}' | cut -f2 | LC_ALL=C sort -k2,2 -S20G --parallel=4 | uniq -c | awk '$1 > 1 {print $2}' > $SPECIES_DIR/$SPECIES.snps.list

python3 generate_catalog.py --shared $SPECIES_DIR/$SPECIES.snps.list --in-list <(find $OUTPUT_DIR -name "*.snps") --out $SPECIES_DIR/$SPECIES.catalog.noAuto.tsv

sed -i "1s/${COMPARATOR_GENOME}-//g" $SPECIES_DIR/$SPECIES.catalog.noAuto.tsv

echo $SPECIES | xargs -I[] bash -c 'head -n 1 species/[]/[].catalog.tsv | awk "$0"' '{for (i=1; i <= NF; ++i) {if ($i == "[]") {print "[]",i}}}' | awk '{printf "cut -f1-%s,%s- %s.catalog.tsv > %s.catalog.noAuto.tsv\n", $2-1, $2+1, $1, $1}' | xargs -I[] bash -c "[]"

# python bin_snps_by_source.py --catalog $SPECIES.catalog.noAuto.tsv --name $SPECIES --genome-info <(cut -f1,22 ./genomes_metadata.tsv) > $SPECIES.continentBin.tsv

python identify_ref_allele.py --catalog $SPECIES_DIR/$SPECIES.catalog.noAuto.tsv --name $SPECIES --coords-dir $OUTPUT_DIR --out $SPECIES_DIR/output/$SPECIES.catalog.noAuto.wtRef.tsv

rm $SPECIES_DIR/$SPECIES.snps.list
rm $SPECIES_DIR/$SPECIES.catalog.noAuto.tsv
