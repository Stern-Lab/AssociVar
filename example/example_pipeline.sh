#!/bin/bash

cd "${0%/*}"

FASTQ_DIR=./fastqs/
BLAST_PATH=/sternadi/home/volume1/shared/tools/ncbi-blast-2.2.30+/bin/
REF_GENOME=./yhr174w_reference.fasta
OUTPUT_DIR=.

# turn fastq files to fasta files and run blast on each file
for filename in $FASTQ_DIR/*.fastq; 
do sed '/^@/!d;s//>/;N' $filename > $filename.fasta;
$BLAST_PATH/makeblastdb -in $filename.fasta -dbtype nucl;
READ_COUNT=$(cat $filename.fasta | wc -l);
$BLAST_PATH/blastn -query $REF_GENOME -task blastn -db $filename.fasta -outfmt "6 sseqid qstart qend qstrand sstart send sstrand length btop" -num_alignments $(($READ_COUNT*10)) -dust no -soft_masking F -perc_identity 60 -evalue 1e-7 -out $filename.blast;
done

# parse blast files to dataframes
mkdir $OUTPUT_DIR/parsed_blasts
python ../association_tests/parse_blasts.py -i $FASTQ_DIR -o $OUTPUT_DIR/parsed_blasts

# the associations part of the script is the heaviest part of the pipeline. It is recommended to run as a few jobs in parallel. For small amounts of data or shorter genomes, it is possible to run as one job as shown here:

#mkdir $OUTPUT_DIR/association_results
#python ../association_tests/create_couples_index_file.py -s 30 -e 1284 -j 1 -o $OUTPUT_DIR/couples_index_file.csv
#python ../association_tests/association_test.py -b $OUTPUT_DIR/parsed_blasts/blasts.csv -m $OUTPUT_DIR/parsed_blasts/mutations.csv -i $OUTPUT_DIR/couples_index_file.csv -o $OUTPUT_DIR/association_results -j 0

# if we were to run this in 5 jobs, the previous commands would look like this:
mkdir $OUTPUT_DIR/association_results2
python ../association_tests/create_couples_index_file.py -s 30 -e 1284 -j 5 -o $OUTPUT_DIR/couples_index_file2.csv
# Because different HPC systems work differently, we used a loop to run the 5 association test jobs. Use your HPC system to submit these jobs
for i in {0..4}; 
do python ../association_tests/association_test.py -b $OUTPUT_DIR/parsed_blasts/blasts.csv -m $OUTPUT_DIR/parsed_blasts/mutations.csv -i $OUTPUT_DIR/couples_index_file2.csv -o $OUTPUT_DIR/association_results2 -j $i; 
done


# unite the association test results into one file and normalize
python ../association_tests/unify_association_results.py -i $OUTPUT_DIR/association_results2 -o $OUTPUT_DIR/association_results.csv
python ../association_tests/normalize_chi2.py -i $OUTPUT_DIR/association_results.csv -o $OUTPUT_DIR/association_results.normalized.csv


# no control sample, so we will use visualization to decide on a cutoff.
python ../association_tests/visualize_association_results.py -i $OUTPUT_DIR/association_results.normalized.csv -o $OUTPUT_DIR/association_results.normalized.png

# after visualization, cutoff chosen is: 100. Get the true viariants:
python ../association_tests/normalized_chi2_get_positions.py -i $OUTPUT_DIR/association_results.normalized.csv -o $OUTPUT_DIR/association_results.normalized.cutoff_100.csv -z 100
mkdir $OUTPUT_DIR/variant_association_tests
python ../association_tests/variant_association_test.py -b $OUTPUT_DIR/parsed_blasts/blasts.csv -m $OUTPUT_DIR/parsed_blasts/mutations.csv -p $OUTPUT_DIR/association_results.normalized.cutoff_100.csv -o $OUTPUT_DIR/variant_association_tests

### strains and linkage between variants
mkdir $OUTPUT_DIR/strains_and_linkage_results
python ../strains_and_linkage_between_variants/strain_analysis.py -b $OUTPUT_DIR/parsed_blasts/blasts.csv -m $OUTPUT_DIR/parsed_blasts/mutations.csv -p $OUTPUT_DIR/variant_association_tests/variants_chosen.csv -o $OUTPUT_DIR/strains_and_linkage_results/haplotypes.csv  -d 0.24 -s 0.21
python ../strains_and_linkage_between_variants/calculate_linkage.py -b $OUTPUT_DIR/parsed_blasts/blasts.csv -m $OUTPUT_DIR/parsed_blasts/mutations.csv -p $OUTPUT_DIR/variant_association_tests/variants_chosen.csv -o $OUTPUT_DIR/strains_and_linkage_results/linkage.csv 


# gzipping for upload to git
gzip $OUTPUT_DIR/association_results.normalized.csv
