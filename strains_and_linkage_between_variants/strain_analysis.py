#! /usr/local/python_anaconda/bin/python3.4


import pandas as pd
import argparse
import numpy as np

def count_haplotypes(args):
    '''
    This function separates the reads into strains. A strain is defined as any combination 
    of the mutations or lack of mutations that are provided to the function.
    
    Format for the csv with the mutations to check strains for: every mutation should 
    have its own row, the csv should have a header row titled "variant", and the mutations 
    should be written in the following format: "A1664.0G". For example:
    "variant
    A1664.0G
    A535.0G
    T1764.0-"
    The output file variants_chosen.csv from association_test_variant.py can be used here.
    '''
    recognized_mutations = pd.read_csv(args.input_chosen_mutations)['variant'].tolist()
    recognized_positions = [float(p[1:-1]) for p in recognized_mutations]
    blast_df = pd.read_csv(args.input_blast_df)
    # choose only reads that were mapped only once in blast
    blast_df['read_count'] = blast_df.groupby('read')['start_ref'].transform('count')
    blast_df = blast_df[(blast_df.read_count == 1)]
    
    # choose only reads that are mapped from at least start_pos_read to end_pos_read
    blast_df = blast_df[(blast_df.start_ref < min(recognized_positions)) & (blast_df.end_ref > max(recognized_positions))]
    blast_df = blast_df[['read', 'start_ref', 'end_ref']]
    
    mutations_df = pd.read_csv(args.input_mutation_df)
    mutations_df = mutations_df[(mutations_df.ref != '-')]
    
    # drop reads containing a variation that is not the recognized mutation or the WT in the positions we are checking combinations for.
    reads_to_drop = mutations_df[(mutations_df.position.isin(recognized_positions)) & ~(mutations_df.full_mutation.isin(recognized_mutations))].read.tolist()
    blast_df = blast_df[~(blast_df.read.isin(reads_to_drop))]
    
    mutations_df = mutations_df[mutations_df.full_mutation.isin(recognized_mutations)]
    df = pd.merge(mutations_df, blast_df[['read']], how='right', on='read')
    
    # simple frequencies for every mutation, keep only mutations over mininmal mutation frequency (default 0)
    base_variants = []
    for m in df.full_mutation.drop_duplicates().tolist():
        m_freq = (len(df[df.full_mutation == m].read.drop_duplicates()) / len(df.read.drop_duplicates()))
        if (m_freq >= args.minimal_mutation_frequency) and (str(m) != 'nan'):
            base_variants.append(m)
    
    mutations_df = mutations_df[mutations_df.full_mutation.isin(base_variants)]
    df = pd.merge(mutations_df, blast_df[['read']], how='right', on='read')
    
    df = df.sort_values('position')
    df['full_mutation'] = df.full_mutation.astype(str)
    df['mutations_on_read'] = df.groupby('read')['full_mutation'].transform(', '.join)
    df = df[['mutations_on_read', 'read']].drop_duplicates()
    df_counts = df.groupby('mutations_on_read').read.count().reset_index().rename(columns={'read':'read_count'}).sort_values('read_count', ascending=False)
    df_counts['read_frequency'] = df_counts.read_count / df_counts.read_count.sum()
    df_counts['mutations_on_read'] = df_counts.mutations_on_read.str.replace('nan', 'WT')
    
    
    ####### now decide which strains are believable
    df_counts_reliable = []
    # classify strains within each base variant
    for i in base_variants:
        print(i)
        df1 = df_counts[df_counts.mutations_on_read.str.contains(i)].copy().reset_index(drop=True)
        df1['base_variant'] = i
        df1['relative_frequency'] = df1.read_frequency / df1.read_frequency.sum()
        df1 = df1.sort_values('relative_frequency', ascending=False)
        df1 = choose_believable_strains(df1, args.substitution_error_cutoff, args.deletion_error_cutoff)
        df_counts_reliable.append(df1)
    # WT as base variant
    df1 = df_counts.copy().reset_index(drop=True)
    df1['base_variant'] = 'WT'
    df1['relative_frequency'] = df1.read_frequency / df1.read_frequency.sum()
    df1 = df1.sort_values('relative_frequency', ascending=False)
    df1 = choose_believable_strains(df1, args.substitution_error_cutoff, args.deletion_error_cutoff)
    df_counts_reliable.append(df1) 
    
    df_counts_reliable = pd.concat(df_counts_reliable)
    df_counts_reliable.to_csv(args.output_file, index=False)

    # create results with only believable strains, and recalculate their relative frequency.
    df_reliable = df_counts_reliable[df_counts_reliable.believable == True][['mutations_on_read', 'read_count']].drop_duplicates().copy()
    df_reliable['frequency'] = df_reliable.read_count / df_reliable.read_count.sum()
    df_reliable = df_reliable.rename(columns={'mutations_on_read':'strain'})
    df_reliable.sort_values('frequency', ascending=False).to_csv(args.output_file + '.reliable.csv', index=False)
    return


def choose_believable_strains(df, substitution_error_cutoff, deletion_error_cutoff):
    df['believable'] = False
    df['closest_strain'] = None
    df['smallest_diff'] = None
    df.at[0, 'believable'] = True
    for i in range(1, len(df)):
        strain1 = df.at[i, 'mutations_on_read']
        df1 = df[:i].copy()
        strains1 = df1[df1.believable == True].mutations_on_read.tolist()
        closest_strain = None
        smallest_diff = None
        for s in strains1:
            diffs = list(set.symmetric_difference(set(s.split(', ')), set(strain1.split(', '))))
            if 'WT' in diffs:
                diffs.remove('WT')
            if closest_strain == None or len(smallest_diff) > len(diffs):
                closest_strain = s
                smallest_diff = diffs
            elif len(smallest_diff) == len(diffs):
                if [i[-1] for i in smallest_diff].count('-') >= [i[-1] for i in diffs].count('-'):
                    smallest_diff = diffs
        df.at[i, 'closest_strain'] = closest_strain
        df.at[i, 'smallest_diff'] = smallest_diff
        if ((substitution_error_cutoff**(len(smallest_diff) - [i[-1] for i in smallest_diff].count('-'))) * (deletion_error_cutoff**([i[-1] for i in smallest_diff].count('-')))) <= df.at[i, 'relative_frequency']:
            df.at[i, 'believable'] = True
    return df


    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--input_blast_df", type=str, help="path to blasts df csv", required=True)
    parser.add_argument('-m', '--input_mutation_df', type=str, help='path to mutations df csv', required=True)
    parser.add_argument('-p', '--input_chosen_mutations', type=str, help='path to csv file with mutations to separate into strains. Every mutation should have its own row, the header row titled "variant", and the mutations should be written in the following format: "A1664.0G". The output file variants_chosen.csv from association_test_variant.py can be used here.', required=True)
    parser.add_argument("-o", "--output_file", type=str, help="a path to an output file", required=True)
    parser.add_argument('-d', '--deletion_error_cutoff', type=float, required=True)
    parser.add_argument('-s', '--substitution_error_cutoff', type=float, required=True)
    parser.add_argument('-f', '--minimal_mutation_frequency', type=float, required=False , default=0, help='frequency cutoff for a single mutation. Only mutations that are on the list and appear at least in this frequency in the population will be included in the strain analysis.')
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)
    count_haplotypes(args)

