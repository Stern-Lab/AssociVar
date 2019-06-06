#! /usr/local/python_anaconda/bin/python3.4


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

'''
This script contains functions to help visualize the association test results.
Provide pandas dataframe created by unify_association_results.py.
'''

#df = pd.read_csv([unify_association_results.py output file])

def association_scatter_chi2(df, out_png, proximity_limit=15, start_pos=False, end_pos=False):
    '''
    For every position, plot all the chi square values, except for values of tests 
    between that position and positions [proximity_limit] bases away or closer (default 15).
    Can use start_pos and end_pos to zoom in on part of the genome.
    '''
    df = df[(df.pos1 - df.pos2).abs() > proximity_limit]
    if start_pos and end_pos:
        df = df[df.pos1.isin(range(start_pos, end_pos))]
    plot = df.plot(x='pos1', y='chi2', kind='scatter')
    fig = plot.get_figure()
    fig.savefig(out_png)
    return

def association_scatter_modifed_zscore(df, out_png, proximity_limit=15, start_pos=False, end_pos=False):
    '''
    For every position, plot all the modified_zscore values, except for values of tests 
    between that position and positions [proximity_limit] bases away or closer (default 15).
    Can use start_pos and end_pos to zoom in on part of the genome.
    '''
    df = df[(df.pos1 - df.pos2).abs() > proximity_limit]
    if start_pos and end_pos:
        df = df[df.pos1.isin(range(start_pos, end_pos))]
    plot = df.plot(x='pos1', y='modified_zscore', kind='scatter')
    fig = plot.get_figure()
    fig.savefig(out_png)
    return

def matrix_heatmap(df, start_pos_a, end_pos_a, start_pos_b, end_pos_b, output_png):
    '''
    Plot a heatmap for chi square values between the range of positions chosen.
    '''
    df = df[(df.pos1 >= start_pos_a) & (df.pos1 <= end_pos_a) & (df.pos2 >= start_pos_b) & (df.pos2 <= end_pos_b)].drop_duplicates()
    df = df[(df.pos2 != df.pos1)]
    #df['chi2_log'] = np.log(df.chi2)
    pivot = df.drop_duplicates().pivot(index='pos1', columns='pos2', values='chi2')
    fig, ax = plt.subplots(figsize=(20, 20))
    sns.heatmap(pivot, ax=ax)
    fig.savefig(output_png, dpi=800)
    return

def scaled_back_full_genome_heatmap(df,output_png, scale=10):
    '''
    Plot a heatmap for the chi squares for the whole genome, but scale it down. 
    Use one value for every x by x square of positions, where x is defined by 
    scale (default 10).
    '''
    df = df.drop_duplicates()
    df['pos1_scaled'] = df.pos1 // scale * scale
    df['pos2_scaled'] = df.pos2 // scale * scale
    df['scaled_chi2_max'] = df.groupby(['pos1_scaled', 'pos2_scaled']).chi2.transform('max')
    df = df[['pos1_scaled', 'pos2_scaled', 'scaled_chi2_max']].drop_duplicates()
    df['scaled_chi2_max_log'] = np.log(df.scaled_chi2_max)
    
    df = df[~((df.pos1_scaled == df.pos2_scaled) | ((df.pos2_scaled - df.pos1_scaled <= scale) & (df.pos1_scaled - df.pos2_scaled <= scale)))]
    
    pivot = df.drop_duplicates().pivot(index='pos1_scaled', columns='pos2_scaled', values='scaled_chi2_max_log')
    fig, ax = plt.subplots(figsize=(20, 20))
    sns.heatmap(pivot, ax=ax)
    fig.savefig(output_png, dpi=800)
    return