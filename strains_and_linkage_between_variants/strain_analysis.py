#! /usr/local/python_anaconda/bin/python3.4


import pandas as pd
import argparse

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
    mutations_for_downstream = []
    for m in df.full_mutation.drop_duplicates().tolist():
        if (len(df[df.full_mutation == m].read.drop_duplicates()) / len(df.read.drop_duplicates())) >= args.minimal_mutation_frequency:
            mutations_for_downstream.append(m)
    
    mutations_df = mutations_df[mutations_df.full_mutation.isin(mutations_for_downstream)]
    df = pd.merge(mutations_df, blast_df[['read']], how='right', on='read')
    
    df = df.sort_values('position')
    df['full_mutation'] = df.full_mutation.astype(str)
    df['mutations_on_read'] = df.groupby('read')['full_mutation'].transform(', '.join)
    df = df[['mutations_on_read', 'read']].drop_duplicates()
    df_counts = df.groupby('mutations_on_read').read.count().reset_index().rename(columns={'read':'read_count'}).sort_values('read_count', ascending=False)
    df_counts['read_frequency'] = df_counts.read_count / df_counts.read_count.sum()
    df_counts['mutations_on_read'] = df_counts.mutations_on_read.str.replace('nan', 'WT')
    
    # For every strain, calculate its percentage out of all the population containing the 
    # lowest appearing variant in the strain. This percentage can later be used as a cutoff
    # using the 90th percentile error frequency for example.
    recognized_dict = {}
    for i in recognized_mutations:
        recognized_dict[i] = df_counts[df_counts.mutations_on_read.str.contains(i)].read_frequency.sum()
    df_counts[['critical_variant_total_frequency', 'critical_variant']] = df_counts.mutations_on_read.apply(lambda x: get_critical_variant_freq(x, recognized_dict)).apply(pd.Series)
    df_counts['frequency_for_error_cutoff'] = df_counts.read_frequency / df_counts.critical_variant_total_frequency
    df_counts = choose_believable_strains(df_counts.reset_index(drop=True), args.substitution_error_cutoff, args.deletion_error_cutoff)
    df_counts.to_csv(args.output_file, index=False)
    
    # create results with only believable strains, and recalculate frequencies so they add up to 1.
    b = df_counts[df_counts.believable == True][['mutations_on_read', 'read_count']].copy()
    b['frequency'] = b.read_count / b.read_count.sum()
    b = b.rename(columns={'mutations_on_read':'strain'})
    b.to_csv(args.output_file + '.believable.csv', index=False)
    return

def get_critical_variant_freq(mutations_on_read, recognized_dict):
    mutations_on_read = mutations_on_read.split(', ')
    smallest_variant = None
    smallest_freq = 1.0
    for m in recognized_dict:
        if m in mutations_on_read:
            if recognized_dict[m] < smallest_freq:
                smallest_freq = recognized_dict[m]
                smallest_variant = m
        else:
            if (1.0 - recognized_dict[m]) < smallest_freq:
                smallest_freq = (1.0 - recognized_dict[m])
                smallest_variant = m[:-1] + m[0]
    return (smallest_freq, smallest_variant)


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
        if ((substitution_error_cutoff**(len(smallest_diff) - [i[-1] for i in smallest_diff].count('-'))) * (deletion_error_cutoff**([i[-1] for i in smallest_diff].count('-')))) <= df.at[i, 'frequency_for_error_cutoff']:
            df.at[i, 'believable'] = True
    return df


    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--input_blast_df", type=str, help="path to blasts df csv", required=True)
    parser.add_argument('-m', '--input_mutation_df', type=str, help='path to mutations df csv', required=True)
    parser.add_argument('-p', '--input_chosen_mutations', type=str, help='path to csv file with mutations to separate into strains. Every mutation should have its own row, the header row titled "variant", and the mutations should be written in the following format: "A1664.0G". The output file variants_chosen.csv from association_test_variant.py can be used here.', required=True)
    parser.add_argument("-o", "--output_file", type=str, help="a path to an output file", required=True)
    parser.add_argument('-f', '--minimal_mutation_frequency', type=float, required=False , default=0, help='frequency cutoff for a single mutation. Only mutations that are on the list and appear at least in this frequency in the population will be included in the strain analysis.')
    parser.add_argument('-d', '--deletion_error_cutoff', type=float, required=True)
    parser.add_argument('-s', '--substitution_error_cutoff', type=float, required=True)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)
    count_haplotypes(args)

