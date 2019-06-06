#! /usr/local/python_anaconda/bin/python3.4


import pandas as pd
import os
import argparse
from tqdm import tqdm

    
if __name__ == "__main__":
    '''
    This script gets a directory with numbered directories, each one containing
    the chi2 results for that directory in a csv called chi2_results.csv. 
    It unifies these files into one csv.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_results_directory", type=str, help="input path to results, same path given as output path to association_results.py" , required=True)
    parser.add_argument("-o", "--output_csv", type=str, help="a path to an output_csv", required=True)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)
    dfs = []
    files = [args.input_results_directory + '/' + d  + '/chi2_results.csv' for d in os.listdir(args.input_results_directory) if d.isnumeric()]
    for f in tqdm(files):
        dfs.append(pd.read_csv(f))
    df = pd.concat(dfs)
    df = df.drop_duplicates()
    df.to_csv(args.output_csv, index=False)