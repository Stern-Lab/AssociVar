#! /usr/local/python_anaconda/bin/python3.4

import pandas as pd
import argparse


def main(args):
    couples = []
    # create couple for every combination of positions between start and end positions.
    for i in range(args.start_position, args.end_position + 1):
        for j in range(i, args.end_position + 1):
            couples.append((i, j))
    couples = pd.DataFrame(couples, columns=['i', 'j'])
    # split couples into groups (number of groups == number of jobs) by using the index.
    couples['association_index'] = couples.index.astype(int) % args.number_of_jobs
    couples.to_csv(args.output_file, index=False)
    return
        
      
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--start_position', type=int, help='start position for association tests', default = 1, required=False)
    parser.add_argument('-e', '--end_position', type=int, help='end position for association tests', required=True)
    parser.add_argument('-j', '--number_of_jobs', type=int, help='number of pbs jobs to split association tests into', required=True)
    parser.add_argument('-o', '--output_file', type=str, help='path to save index file to', required=True)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)
    main(args)