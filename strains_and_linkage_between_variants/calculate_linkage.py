# -*- coding: utf-8 -*-

import sys
import pandas as pd
import itertools
import os
import numpy as np
import argparse


def mutation_pairings_statistics(args):
    '''
    This function calculates the conditional probabilities for every mutation in list
    to appear with any other mutation in the given list, including the porbabilities 
    of the WT variants for those positions. The function gets the blast dataframe path, 
    the mutations dataframe path, a csv path with the mutations to check and a directory to 
    save the results to.
    
    Foramt for the csv with the mutations to check linkage for: every mutation should 
    have its own row, the csv should have a header row titled "variant", and the mutations should be written 
    in the following format: "A1664.0G". For example:
    "variant
    A1664.0G
    A535.0G
    T1764.0-"
    '''
    blast_df = pd.read_csv(args.input_blast_df)
    mutations_df = pd.read_csv(args.input_mutation_df)
    recognized_mutations = pd.read_csv(args.input_chosen_mutations)['variant'].tolist()

    
    # choose only reads that were mapped only once in blast
    blast_df['read_count'] = blast_df.groupby('read')['start_ref'].transform('count')
    blast_df = blast_df[(blast_df.read_count == 1)]
  
    # conditional statistics
    matrix = pd.DataFrame(np.zeros(shape=(len(recognized_mutations)*2,len(recognized_mutations)*2)), columns = recognized_mutations + [f[1:-1].split('.')[0] for f in recognized_mutations], index = recognized_mutations + [f[1:-1].split('.')[0] for f in recognized_mutations])
    for mutation_1 in recognized_mutations:
        for mutation_2 in recognized_mutations:
            print (mutation_1, mutation_2)
            b = blast_df.copy()
            m = mutations_df.copy()
            # keep only alignments containing both positions for mutations 1 and 2.
            b = b[(b.start_ref < int(mutation_1[1:-1].split('.')[0])) & (b.end_ref > int(mutation_1[1:-1].split('.')[0])) & (b.start_ref < int(mutation_2[1:-1].split('.')[0])) & (b.end_ref > int(mutation_2[1:-1].split('.')[0]))]
            b = b[['read']]
            # look only at mutations in positions i and j
            m = m[m.position.isin([int(mutation_1[1:-1].split('.')[0]), int(mutation_2[1:-1].split('.')[0])])]
            # keep all reads that contain positions i and j (both containing mutations in i and/or j 
            # and not containing mutations there), and their appropriate mutations.
            relevant_mutations = pd.merge(b, m, on='read', how='left').fillna(0)
            
            # only keep reads with WT or known mutations. Discard reads with unknown mutations in the checked positions.
            # reads to drop:
            reads_to_drop = m[~(m.full_mutation.isin(recognized_mutations))].read.tolist()
            # dropping:
            relevant_mutations = relevant_mutations[~(relevant_mutations.read.isin(reads_to_drop))]
            
            # calculate conditional probabilities
            reads_with_mut_1 = len(relevant_mutations[relevant_mutations.full_mutation == mutation_1][['read']].drop_duplicates())
            reads_with_wt_1 = len(relevant_mutations[['read']].drop_duplicates()) - reads_with_mut_1
            reads_with_mut_1_and_mut_2 = len(pd.merge(relevant_mutations[(relevant_mutations.full_mutation == mutation_1)], relevant_mutations[(relevant_mutations.full_mutation == mutation_2)], on='read', how='inner').read.drop_duplicates())
            reads_with_mut_1_and_wt_2 = reads_with_mut_1 - reads_with_mut_1_and_mut_2
            reads_with_wt_1_and_wt_2 = len(relevant_mutations[['read']].drop_duplicates()) - len(relevant_mutations[((relevant_mutations.full_mutation == mutation_1) | (relevant_mutations.full_mutation == mutation_2))][['read']].drop_duplicates())
            reads_with_wt_1_and_mut_2 = reads_with_wt_1 - reads_with_wt_1_and_wt_2
            
            # update results csv
            matrix.at[mutation_1, mutation_2] = reads_with_mut_1_and_mut_2 / reads_with_mut_1
            matrix.at[mutation_1, mutation_2[1:-1].split('.')[0]] = reads_with_mut_1_and_wt_2 / reads_with_mut_1
            matrix.at[mutation_1[1:-1].split('.')[0], mutation_2[1:-1].split('.')[0]] = reads_with_wt_1_and_wt_2 / reads_with_wt_1
            matrix.at[mutation_1[1:-1].split('.')[0], mutation_2] = reads_with_wt_1_and_mut_2 / reads_with_wt_1
    matrix.to_csv(args.output_file)
            
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--input_blast_df", type=str, help="path to blasts df csv", required=True)
    parser.add_argument('-m', '--input_mutation_df', type=str, help='path to mutations df csv', required=True)
    parser.add_argument('-p', '--input_chosen_mutation', type=str, help='path to csv file with mutations to check linkage of. Every mutation should have its own row, a header row titled "variant", and the mutations should be written in the following format: "A1664.0G"', required=True)
    parser.add_argument("-o", "--output_file", type=str, help="a path to an output file", required=True)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)
    mutation_pairings_statistics(args)
