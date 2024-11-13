#!/usr/bin/env python

import argparse

import numpy as np
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-in_nucl', '--in_nucl',
                        type=str, required=True,
                        help="Input nucleotide results table")
    
    parser.add_argument('-in_prot', '--in_prot',
                        type=str, required=True,
                        help="Input protein results table")
    
    parser.add_argument('-in_cove', '--in_cove',
                        type=str, required=True,
                        help="Input coverage results table")
    
    parser.add_argument('-out_p', '--out_pattern',
                        type=str, required=True, 
                        help="Patter for output files")
    
    parser.add_argument('-v', '--version',
                        action='version', version="0.0.1")
        
    return parser.parse_args()


if __name__ == '__main__':

    args = parse_arguments()

    df_nucl = pd.read_csv(args.in_nucl, sep="\t")
    df_prot = pd.read_csv(args.in_prot, sep="\t")
    df_cove = pd.read_csv(args.in_cove, sep="\t")

    df = pd.merge(df_prot, df_nucl,
                  right_on='query', left_on='contig_name',
                  suffixes=("_prot", "_nucl"), how='outer'
                  ).merge(df_cove, left_on='contig_name', right_on='#rname')
    
    df.sort_values(by="believe").to_csv(args.out_pattern + "_extended_table.tsv", index=False, sep="\t")
    df_cove[['#rname', 'meandepth']].to_csv(f"{args.out_pattern}_two_col_coverage.tsv", index=False, header=False, sep="\t")
    df[~df["query"].isna()]["query"].to_csv(f"{args.out_pattern}_for_extension.txt", index=False, header=False)
