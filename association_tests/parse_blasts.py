#! /usr/local/python_anaconda/bin/python3.4

import pandas as pd
import re
import math
import os
import argparse


MUTATION_PATTERN = re.compile('[ACGTN]{2}')
INSERTION_PATTERN = re.compile('-[ACGTN]')
DELETION_PATTERN = re.compile('[ACGTN]-')
NUMBER_PATTERN = re.compile('\d+')
GENERAL_BLAST_PATTERN = re.compile('[ACGTN-]{2}|\d+')

def blast_to_df(blast_path):
    '''
    Creates pandas dataframe from blast file
    :param blast_path: input file path, needs to be blast output
    :returns: pandas dataframe
    '''
    blast_df = pd.read_csv(blast_path, sep='\t', header = None)
    blast_df.columns = ("read", "start_ref", "end_ref", "start_read", "end_read", "strand", "length", "btop")
    return blast_df

def blast_to_mutations_list(blast_path, out_csv_path=None):
    '''
    Gets a blast file and creates a file containing all single nucleotide mutations
    (point mutations, insertions, and deletions) with the matching read name.
    :param blast_path: input file path, needs to be blast output
    :param out_csv_path: output file path
    :returns: pandas dataframe
    '''
    blast = blast_to_df(blast_path)
    blast['parsed_btop'] = blast.apply(parse_btop, axis=1)
    reads = []
    positions = []
    refs = []
    bases = []
    for index, row in blast.iterrows():
        for mutations in list(row['parsed_btop'].values()):
            for i in mutations:
                reads.append(row['read'])
                positions.append(i[0])
                refs.append(i[1][0])
                bases.append(i[1][1])
    out = pd.DataFrame()
    out['read'] = reads
    out['position'] = positions
    out['ref'] = refs
    out['base'] = bases
    if out_csv_path:
        out.to_csv(out_csv_path, index=False)
    return out

def parse_btop(row):
    '''
    Gets a pandas dataframe row from a blast file and parses the btop field.
    Used in blast_to_mutations_list function.
    '''
    insertions = []
    deletions = []
    mutations = []
    location = float(row['start_ref'])
    btops = GENERAL_BLAST_PATTERN.findall(row['btop'])
    for b in btops:
        if NUMBER_PATTERN.findall(b) != []:
            location = float(math.ceil(location) + int(b))
        elif MUTATION_PATTERN.findall(b) != []:
            mutations.append((float(math.ceil(location)), b))
            location = float(math.ceil(location) + 1)
        elif INSERTION_PATTERN.findall(b) != []:
            if location.is_integer():
                location -= 1.0
            location += 0.01
            insertions.append((location, b))
        elif DELETION_PATTERN.findall(b) != []:
            deletions.append((float(math.ceil(location)), b))
            location = float(math.ceil(location) + 1) 
    return {'mutations':mutations, 'deletions':deletions, 'insertions':insertions}


def create_blasts_and_mutations_dfs(args):
    '''
    This function gets an input directory with blast files (_.blast), and created
    two new files in the output directory: a dataframe containing all of the blast 
    output files and a dataframe containing every mutation for every read in the
    blast outputs.
    '''
    blast_dfs = []
    mutations_dfs = []
    for i in [args.input_dir + '/' + f for f in os.listdir(args.input_dir) if f.endswith('.blast')]:
        print(i)
        mutations_df = blast_to_mutations_list(i)
        mutations_df['full_mutation'] = mutations_df.ref.astype(str) + mutations_df.position.astype(str) + mutations_df.base.astype(str)
        blast_df = blast_to_df(i)
        mutations_dfs.append(mutations_df)
        blast_dfs.append(blast_df)
    mutations_df = pd.concat(mutations_dfs)
    blast_df = pd.concat(blast_dfs)
    mutations_df.to_csv(args.output + '/' + 'mutations.csv', index=False)
    blast_df.to_csv(args.output + '/' + 'blasts.csv', index=False)
    return mutations_df, blast_df

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", type=str, help="input directory with blast files (_.blast)", required=True)
    parser.add_argument("-o", "--output", type=str, help="a path to an output directory", required=True)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)
    create_blasts_and_mutations_dfs(args)