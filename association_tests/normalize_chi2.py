#! /usr/local/python_anaconda/bin/python3.4

import pandas as pd
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_path', type=str, help='path to an input csv of chi scores, created by unify_association_results.py', required=True)
    parser.add_argument("-o", "--output_path", type=str, help="path to an output csv to save results to", required=True)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)

    df = pd.read_csv(args.input_path)
    df = df.drop_duplicates()
    df = df[(df.pos1 - df.pos2).abs() > 15] 
    
    df['pos1_median'] = df.groupby('pos1').chi2.transform('median')
    df['pos1_mad'] = df.groupby('pos1').chi2.transform('mad')
    df['modified_zscore'] = 0.6745 * (df.chi2 - df.pos1_median) / df.pos1_mad
    df = df.drop(['pos1_mad', 'pos1_median'], axis=1)
    
    # check if higher than neighbors
    df = df.sort_values(['pos1', 'pos2'])
    df['temp1'] = df.modified_zscore.diff()
    df['temp2'] = df.modified_zscore.diff(periods = -1)
    df = df.sort_values(['pos2', 'pos1'])
    df['temp3'] = df.modified_zscore.diff()
    df['temp4'] = df.modified_zscore.diff(periods = -1)
    df['is_peak'] = ((df.temp1 > 0) & (df.temp2 > 0) & (df.temp3 > 0) & (df.temp4 > 0))
    df = df.drop(['temp1', 'temp2', 'temp3', 'temp4'], axis=1)
    
    df = df.sort_values('modified_zscore', ascending=False)
    df.to_csv(args.output_path, index=False)