#! /usr/local/python_anaconda/bin/python3.4

import pandas as pd
import argparse
import math

def is_peak(pos1, pos2, modified_zscore, df):
    a = df[((df.pos1 == pos1) & (df.pos2.isin([pos2 - 1, pos2 + 1]))) | ((df.pos2 == pos2) & (df.pos1.isin([pos1 - 1, pos1 + 1])))]
    if len(a) == 0:
        return True
    return modified_zscore > a.modified_zscore.max()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_csv", type=str, help="path to input csv with pos1, pos2 and zscore columns", required=True)
    parser.add_argument("-o", "--output_csv", type=str, help="a path to output csv to save results to", required=True)
    parser.add_argument("-c", "--confidence_percentile", type=float, help="confidence percentile, percentile of positions to be identified as not containing real mutations when using this cutoff.", required=True)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)
    df = pd.read_csv(args.input_csv)
    df = df.sort_values('modified_zscore', ascending=False)
    
    position_count = len(df.pos1.drop_duplicates())
    positions_above_cutoff = math.ceil(position_count * (100 - args.confidence_percentile) / 100)
    cutoff = 0
    for i in range(math.floor(positions_above_cutoff / 2), position_count):
        if cutoff == 0:
            df1 = df[:i].copy()
            df1['is_peak'] = df1.apply(lambda x: is_peak(x.pos1, x.pos2, x.modified_zscore, df1), axis=1)
            df1 = df1[df1.is_peak == True]
            positions1 = list(set(df1.pos1.tolist() + df1.pos2.tolist()))
            if len(positions1) >= positions_above_cutoff:
                cutoff = df1.modified_zscore.min()
    df1.to_csv(args.output_csv, index=False)
    print(cutoff)