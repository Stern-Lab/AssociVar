### Association Tests
This series of steps are used to run chi-square association tests between every pair of positions in the genome, in order to separate real mutations from technical errors in MinION sequencing results. The analysis is based on the assumption that real mutations have a tendency to either appear together or not appear together on the same genomes, and thus are not independent of each other, while technical sequencing errors are random and so appear to be independent of each other in the population.


1. Run BLAST on all fastq files as described in [AccuNGS paper](https://github.com/SternLabTAU/AccuNGS).
To do this, either use AccuNGS script or run BLAST separately as follows:

   - BLAST (v2.7.1+) the merged fastq against the reference. 
Parameters: 
ref_genome - reference genome (FASTA format) 
in_fastq - fastq file 
max_num_alignments - maximum number of alignments (typically 10x of the length of the input file) 
pcID_blast - percents identity of each alignment to the reference. suggested for minion: 60.
out_blast - blast output file

   - Typical use case: makeblastdb -in ${in_fastq} -dbtype nucl blastn -query ${ref_genome} -task blastn -db ${in_fastq} -outfmt "6 sseqid qstart qend qstrand sstart send sstrand length btop" -num_alignments ${max_num_alignments} -dust no -soft_masking F -perc_identity ${pcID_blast} -evalue 1e-7 -out ${out_blast}";


2. Run parse_blasts.py to create two output files: a dataframe of all blast outputs and a dataframe
containing every mutation per read from the blast outputs. The script gets an input directory containing blast
results files (_.blast), and creates the two new csv files in the output directory.

   - usage: parse_blasts.py [-h] -i INPUT_DIR -o OUTPUT


3. Run create_couples_index_file.py to create an index file used by the association test script.
The file is a csv containing the pairs of positions to check, and also an association index for the
association test script to use when being run as many jobs in parallel (choose number of jobs to split to). If planning to run as one job, choose 1 for -j, and on step 4 (association_test.py) choose 0 for -j (PBS_JOB_ARRAY_ID).

   - usage: create_couples_index_file.py [-h] [-s START_POSITION] -e END_POSITION
                                    -j NUMBER_OF_JOBS -o OUTPUT_FILE


4. Run association_test.py scripts, recommended as a job array. This uses the blasts dataframe,
mutations dataframe and the position couple index file created in 2-3. Each job runs chi square association tests
for every pair of positions from the position couple index file for which the index matches the job id. It saves
the contingency tables for every pair of positions, and a csv containing all of the chi2 results. Also recieves
start_pos_read and end_pos_read, and only reads spanning these positions are used in the analysis.

   - usage: association_test.py [-h] -b INPUT_BLAST_DF -m INPUT_MUTATION_DF -i
                           INPUT_POSITION_PAIRS_DF -o OUTPUT_DIR -s
                           START_POS_READ -e END_POS_READ -j PBS_JOB_ARRAY_ID

   - This is the heaviest part of this code. For a genome of 3569 bases and results of 100,000 reads split into 100 jobs of 5000mb, this took 6 hours.
   
5. Unite association test results by running unify_association_tests.py. This creates a csv
with the chi square statistic value and the p-value for every pair of positions the test was performed for.

   - usage: unify_association_results.py [-h] -i INPUT_RESULTS_DIRECTORY -o
                                    OUTPUT_CSV

6. Normalize chi square results by using the modified z-test and determines for each association if it is higher than its two neighbors (is_peak column). Gets an input path of a csv with all chi2 results and an output path to   write the z-test results to.

   - usage: normalize_chi2.py [-h] -i INPUT_PATH -o OUTPUT_PATH

							   
If this is a control sample used to find a cutoff for the association, use normalized_chi2_find_cutoff.py. This script gets a csv with the pos1, pos2 and modified zscores created by normalize_chi2.py, 
and a confidence percentile. For example, for a cutoff percentile of 99.9, the cutoff is determined as the score that identifies 0.1 percent of the positions as having significant associations. The cutoff
printed to the screen. 

   - usage: normalized_chi2_find_cutoff.py [-h] -i INPUT_CSV -o OUTPUT_CSV -c
                                     CONFIDENCE_PERCENTILE		   

If this is not a control sample, use normalized_chi2_get_positions.py to get the positions that are identified as having real variants by AssociVar. This script gets a csv with the pos1, pos2, modified zscores and is_peak column created by normalize_chi2.py, 
an output path to write the positions to, and a modified_zscore_cutoff to use.
					 
   - usage: normalized_chi2_get_positions.py [-h] -i INPUT_CSV -o OUTPUT_CSV -z
                                MODIFIED_ZSCORE_CUTOFF

We also provide visualization tools to help identify positions with real mutations in the association test results. Helpful functions are provided in tools_to_visualize_association_results.py.

To transform the identified positions into specific variants, use variant_association_test.py. The script uses the blasts dataframe and mutations dataframe created in steps 2-3. The script also gets
the positions to analyze in the format of csv with no header where every row is a position to analyze (option p) as created by normalized_chi2_get_positions.py. The script
writes the results and the final variant list in chosen_variants.csv to the output directory.
   - usage: variant_association_test.py [-h] -b INPUT_BLAST_DF -m INPUT_MUTATION_DF
                                   -p INPUT_CHOSEN_POSITIONS -o
                                   OUTPUT_DIR
