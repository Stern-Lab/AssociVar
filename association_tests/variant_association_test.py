#! /usr/local/python_anaconda/bin/python3.4


import pandas as pd
import numpy as np
import argparse
import scipy.stats
import itertools
from tqdm import tqdm
import os

def variant_association_test(args):
    '''
    '''    
    recognized_positions = [float(p) for p in pd.read_csv(args.input_chosen_positions, header=None)[0].tolist()] 

    blast_df = pd.read_csv(args.input_blast_df)
    # choose only reads that were mapped only once in blast
    blast_df['read_count'] = blast_df.groupby('read')['start_ref'].transform('count')
    blast_df = blast_df[(blast_df.read_count == 1)]
    blast_df = blast_df[['read', 'start_ref', 'end_ref']]
    blast_df = blast_df[(blast_df.start_ref < min(recognized_positions)) & (blast_df.end_ref > max(recognized_positions))]
    blast_df = blast_df[['read']]
    
    mutations_df = pd.read_csv(args.input_mutation_df)
    mutations_df = mutations_df[(mutations_df.ref != '-')]
    mutations_df = pd.merge(mutations_df, blast_df[['read']], how='right', on='read')
    mutations_df = mutations_df[mutations_df.position.isin(recognized_positions)]
    
    chi2_data = []
    
    combinations = [c for c in list(itertools.combinations(sorted(mutations_df[mutations_df.position.isin(recognized_positions)][['full_mutation']].drop_duplicates().full_mutation.tolist()), 2)) if c[0][1:-1] != c[1][1:-1]]
    
    for (i, j) in tqdm(combinations):
        pos_i, pos_j = float(i[1:-1]), float(j[1:-1])
        temp_matrix = pd.DataFrame(np.zeros(shape=(2,2)), columns=[j,0], index=[i,0])
        b = blast_df.copy()
        m = mutations_df.copy()
        # look only at mutations in positions i and j
        m = m[m.position.isin([pos_i, pos_j])]
        
        # drop reads containing a variation that is not the recognized mutation or the WT in the positions we are checking combinations for.
        reads_to_drop = m[(m.position.isin([pos_i, pos_j])) & ~(m.full_mutation.isin([i, j]))].read.tolist()
        b = b[~(b.read.isin(reads_to_drop))]
        # keep all reads that contain positions i and j (both containing mutations in i and/or j 
        # and not containing mutations there), and their appropriate mutations.
        relevant_mutations = pd.merge(b, m, on='read', how='left').fillna(0)
            
        # from reads containing both positions with either WT or the variant in examination, 
        # count reads with mutation i, with mutation j, with mutations in both or with no 
        # mutations and create matrix for association test.
        reads_count = len(relevant_mutations[['read']].drop_duplicates())
        reads_with_mut_i_mut_j  = len(pd.merge(relevant_mutations[(relevant_mutations.full_mutation == i)], relevant_mutations[(relevant_mutations.full_mutation == j)], on='read', how='inner').read.drop_duplicates())
        reads_with_mut_j = len(relevant_mutations[(relevant_mutations.full_mutation == j)][['read']].drop_duplicates())
        reads_with_mut_i = len(relevant_mutations[(relevant_mutations.full_mutation == i)][['read']].drop_duplicates())
        reads_with_mut_i_wt_j = reads_with_mut_i - reads_with_mut_i_mut_j
        reads_with_wt_i_mut_j = reads_with_mut_j - reads_with_mut_i_mut_j
        reads_with_wt_i_wt_j = reads_count - reads_with_mut_i_mut_j - reads_with_mut_i_wt_j - reads_with_wt_i_mut_j
        temp_matrix.at[i,j] = reads_with_mut_i_mut_j
        temp_matrix.at[i,0] = reads_with_mut_i_wt_j
        temp_matrix.at[0,j] = reads_with_wt_i_mut_j
        temp_matrix.at[0,0] = reads_with_wt_i_wt_j
        # run association test on matrix and write results to file.
        if temp_matrix.sum(axis=0).all() > 0 and temp_matrix.sum(axis=1).all() > 0:
            #temp_matrix.to_csv(args.output_dir + '/extras/' + str(i) + '_' + str(j) + '.csv')
            chi2, pvalue, dof, expected = scipy.stats.chi2_contingency(temp_matrix)
            chi2_data.append((i,j,pvalue,chi2))
            chi2_data.append((j,i,pvalue,chi2))
            
    df = pd.DataFrame(chi2_data, columns=['variant1', 'variant2', 'pvalue', 'chi2'])
    df.to_csv(args.output_dir  + '/variant_association_results.csv', index=False)
    
    a = df.groupby('variant1').chi2.mean().reset_index().sort_values('chi2')
    a['pos1'] = a.variant1.str[1:-1].astype(float)
    a['max_average_chi2'] = a.groupby('pos1').chi2.transform('max')
    a = a[a.max_average_chi2 == a.chi2]
    a[['variant1']].rename(columns={'variant1':'variant'}).to_csv(args.output_dir + '/variants_chosen.csv', index=False)           
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--input_blast_df", type=str, help="path to blasts df csv", required=True)
    parser.add_argument('-m', '--input_mutation_df', type=str, help='path to mutations df csv', required=True)
    parser.add_argument('-p', '--input_chosen_positions', type=str, help='path to csv file with positions to check associations for. Every position should have its own row, no header row.', required=False)
    parser.add_argument("-o", "--output_dir", type=str, help="a path to an output directory", required=True)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)
    variant_association_test(args)
