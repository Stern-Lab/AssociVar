# Association Tests
This series of steps are used to run chi-square association tests between every pair of positions in the genome,
in order to separate real mutations from technical errors in MinION sequencing results. The analysis is based on
the assumption that real mutations have a tendency to either appear together or not appear together on the same genomes,
and thus are not independent of each other, while technical sequencing errors are random and so appear to be independent
of each other in the population.


### Requirements:
BLAST (v2.7.1+)

Anaconda python installation OR python (3.4+) with:
- numpy
- pandas
- tqdm
- matplotlib
- scipy


## Running chi-square tests

1. Run BLAST on fastq file(s).

   To run blast, first turn fastqs into fasta files (for example, using sed: "sed '/^@/!d;s//>/;N' FASTQ > FASTA"

   Next, run BLAST (v2.7.1+) as follows:
   
   Parameters:
     - REF_GENOME - reference genome (FASTA format)
     - FASTA - fasta file
     - MAX_NUM_ALIGNMENTS - maximum number of alignments (typically 10x the number of reads in file)
     - PC_ID_BLAST - percents identity of each alignment to the reference. suggested for minion: 60.
     - OUT_BLAST - blast output file

   Commands
     1. makeblastdb -in IN_FASTA -dbtype nucl
     2. blastn -query REF_GENOME -task blastn -db FASTA -outfmt "6 sseqid qstart qend qstrand sstart send sstrand length btop" -num_alignments MAX_NUM_ALIGNMENTS -dust no -soft_masking F -perc_identity PC_ID_BLAST -evalue 1e-7 -out OUT_BLAST


2. Run parse_blasts.py to create two output files: a dataframe of all blast outputs and a dataframe
containing every mutation per read from the blast outputs. The script gets an input directory containing blast
results file(s) (files with suffix .blast), and creates the two new csv files in the output directory.

   - usage: parse_blasts.py [-h] -i INPUT_DIR -o OUTPUT


3. Run create_couples_index_file.py to create an index file used by the association test script.
The file is a csv containing the pairs of positions to check, and also an association index for the
association test script to use when being run as more than one job in parallel (choose number of jobs to split to in NUMBER_OF_JOBS).
If planning to run as one job, choose 1 for -j, and on the next step (association_test.py) choose 0 for -j (PBS_JOB_ARRAY_ID).
Start and end positions to test associations for are recommended as 30 positons away from the start and end of the
reference sequence.

   - usage: create_couples_index_file.py [-h] -s START_POSITION -e END_POSITION
                                    -j NUMBER_OF_JOBS -o OUTPUT_FILE


4. Run association_test.py scripts, recommended as a job array. This uses the blasts dataframe (INPUT_BLAST_DF),
mutations dataframe (INPUT_MUTATION_DF) and the position couple index file (INPUT_POSITION_PAIRS_DF) created in steps 2-3. Each job runs chi square association tests for every pair of positions from the position couple index file for which the index matches the job id. It saves the contingency tables for every pair of positions, and a csv containing all of the chi2 results.
All sequenced reads are required to span a certain part of the genome, by defualt this is from the start position to the end position specified in create_couple_index_file.py, but this can also be set to a larger area using start_pos_read and end_pos_read.

   - usage: association_test.py [-h] -b INPUT_BLAST_DF -m INPUT_MUTATION_DF -i
                           INPUT_POSITION_PAIRS_DF -o OUTPUT_DIR [-s
                           START_POS_READ] [-e END_POS_READ] -j PBS_JOB_ARRAY_ID

   - This is the heaviest part of this code. For a genome of 3569 bases and results of 100,000 reads split into 100 jobs of 5000mb, this took 3 hours.

5. Unite association test results from previous step by running unify_association_tests.py. INPUT_RESULTS_DIRECTORY is the same as output_dir from previous step.

   - usage: unify_association_results.py [-h] -i INPUT_RESULTS_DIRECTORY -o
                                    OUTPUT_CSV

6. Normalize chi square results by using the modified z-test and determine for each association if it is higher than
its four neighbors (is_peak column). Gets an input path of a csv with all chi2 results and an output path to write the z-test results to.

   - usage: normalize_chi2.py [-h] -i INPUT_PATH -o OUTPUT_PATH

## Deciding on a cutoff
This can be done in two ways - with a control sample and without.

### Finding a cutoff using a control file:
Use normalized_chi2_find_cutoff.py on your control sample. This script gets a csv created by normalize_chi2.py and a confidence percentile as
to set a cutoff for. For example, for a cutoff percentile of 99.9, the cutoff is determined as the score that
identifies 0.1 percent of the positions as having significant associations. The cutoff is printed to the screen, and the
results exceeding the cutoff are written to output_csv.

   - usage: normalized_chi2_find_cutoff.py [-h] -i INPUT_CSV -o OUTPUT_CSV -c
                                     CONFIDENCE_PERCENTILE

### Choosing a cutoff without a control sample:
Use visualize_association_results.py to visualize the normalized association scores created by normalize_chi2.py and choose an appropriate cutoff. The script creates a scatter plot of all the associations per position, where different colors signify
whether the association is a local maximum or not (associations are required to be a local maximum in order to be considered associations between two real mutations).

   - usage: visualize_association_results.py [-h] -i INPUT_ASSOCIATION_RESULTS -o
                                        OUTPUT_PNG
## Getting the true variants
Use normalized_chi2_get_positions.py to get the **positions** that are identified as having real variants by AssociVar.
This script gets a csv with the pos1, pos2, modified zscores and is_peak column created by normalize_chi2.py, an
output path to write the positions to, and a modified_zscore_cutoff to use.

   - usage: normalized_chi2_get_positions.py [-h] -i INPUT_CSV -o OUTPUT_CSV -z
                                MODIFIED_ZSCORE_CUTOFF


To transform the identified **positions** into **variants** (specific nucleotide), use variant_association_test.py. The script uses the
blasts dataframe and mutations dataframe created in steps 2-3. The script also gets the positions to analyze as created by normalized_chi2_get_positions.py. The script writes the chi-square results and the final variant list to the output directory in variant_association_results.csv and chosen_variants.csv, respectively.
   - usage: variant_association_test.py [-h] -b INPUT_BLAST_DF -m INPUT_MUTATION_DF
                                   -p INPUT_CHOSEN_POSITIONS -o
                                   OUTPUT_DIR
                                   
 
 

