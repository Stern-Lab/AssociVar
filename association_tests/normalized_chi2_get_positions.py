#! /usr/local/python_anaconda/bin/python3.4

import pandas as pd
from tqdm import tqdm
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_csv", type=str, help="path to input csv with pos1, pos2, modified_zscore and is_peak columns", required=True)
    parser.add_argument("-o", "--output_csv", type=str, help="a path to output csv to save results to", required=True)
    parser.add_argument("-z", "--modified_zscore_cutoff", type=float, help="modified z score cutoff", required=True)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)
    df = pd.read_csv(args.input_csv)
        
    df = df.sort_values('modified_zscore', ascending=False)
    df = df[df.modified_zscore > args.modified_zscore_cutoff]
    df = df[df.is_peak == True]
    df = pd.DataFrame(list(set(df.pos1.tolist() + df.pos2.tolist())), columns=['pos'])
    df.to_csv(args.output_csv, index=False, header=False)