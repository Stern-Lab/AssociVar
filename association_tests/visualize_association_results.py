#! /usr/local/python_anaconda/bin/python3.4

import os
import pandas as pd
os.environ['QT_QPA_PLATFORM']='offscreen'
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import argparse

def association_scatter_modifed_zscore(df, out_png, proximity_limit=15, start_pos=False, end_pos=False):
    '''
    For every position, plot all the modified_zscore values, except for values of tests 
    between that position and positions [proximity_limit] bases away or closer (default 15).
    Can use start_pos and end_pos to zoom in on part of the genome.
    '''
    fig, ax = plt.subplots()
    df = df[(df.pos1 - df.pos2).abs() > proximity_limit]
    if start_pos and end_pos:
        df = df[df.pos1.isin(range(start_pos, end_pos))]
    df[df.is_peak == False].plot(x='pos1', y='modified_zscore', kind='scatter', color='#CEE9F2', ax=ax, label='local_maximum=False')
    df[df.is_peak == True].plot(x='pos1', y='modified_zscore', kind='scatter', ax=ax, label='local_maximum=True')
    fig.savefig(out_png)
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_association_results", type=str, help="input path to results, as created by normalize_chi2.py" , required=True)
    parser.add_argument("-o", "--output_png", type=str, help="a path to an output png", required=True)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)
    df = pd.read_csv(args.input_association_results)
    association_scatter_modifed_zscore(df, args.output_png, proximity_limit=15, start_pos=False, end_pos=False)
