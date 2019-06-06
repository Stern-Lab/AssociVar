#! /usr/local/python_anaconda/bin/python3.4

import pandas as pd
from tqdm import tqdm
import argparse

def is_peak(pos1, pos2, modified_zscore, df):
    a = df[((df.pos1 == pos1) & (df.pos2.isin([pos2 - 1, pos2 + 1]))) | ((df.pos2 == pos2) & (df.pos1.isin([pos1 - 1, pos1 + 1])))]
    if len(a) == 0:
        return True
    return modified_zscore > a.modified_zscore.max()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_csv", type=str, help="path to input csv with pos1, pos2 and zscore columns", required=True)
    parser.add_argument("-o", "--output_csv", type=str, help="a path to output csv to save results to", required=True)
    parser.add_argument("-z", "--modified_zscore_cutoff", type=float, help="modified z score cutoff", required=True)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)
    df = pd.read_csv(args.input_csv)
        
    df = df.sort_values('modified_zscore', ascending=False)
    df = df[df.modified_zscore > args.modified_zscore_cutoff]
    tqdm.pandas(desc="z score progress")
    df['is_peak'] = df.progress_apply(lambda x: is_peak(x.pos1, x.pos2, x.modified_zscore, df), axis=1)
    #df = df[df.is_peak == True]
    df.to_csv(args.output_csv, index=False)