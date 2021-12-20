#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd


CHR = [str(i) for i in range(1, 23, 1)] + ['X', 'Y']


def filter_CNVs(snp_in, cells_in, cnv_in, file_out):
    with open(cells_in, 'r') as f:
        cells_raw = [i.split('\t')[0] for i in f.read().strip().split('\n')]

    cells = [f'CRC09.TC.{i.replace("-", ".")}Merged' for i in cells_raw]
    cells = [''.join(i.rsplit('.', 1)) for i in cells]
    cells.remove('CRC09.TC.XP.WGSMGMerged')

    df_cnv = pd.read_csv(cnv_in, sep='\t', index_col=[0,1,2])
    # df_cnv.drop('chrY', level='CHR', inplace=True)

    df_cnv = df_cnv[cells]
    df_filter = df_cnv[df_cnv.sum(axis=1) != 48]
    idx_filter = np.array([[i[0][3:], i[1], i[2]] for i in df_filter.index])

    idx_dict = {}
    for chr in CHR:
        idx_dict[chr] = idx_filter[idx_filter[:,0] == chr, 1:].astype(int)

    with open(file_out, 'w') as f_out:
        with open(snp_in, 'r') as f_snp:
            for line in f_snp:
                if line.count('\t') == 0 or line.count('\t') != len(cells) + 3:
                    f_out.write(line)
                else:
                    chr, pos = line.split('\t')[:2]
                    pos = int(pos)

                    if np.any((pos > idx_dict[chr][:,0]) & (pos < idx_dict[chr][:,1])):
                        print(f'Skipping: {chr} {pos}')
                    else:
                        f_out.write(line)





def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-snp', type=str, help='DataFilter candidate sites.')
    parser.add_argument('-cells', type=str, help='Relevant cells.')
    parser.add_argument('-cnv', type=str, help='CNV information.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file. Default = <SNP>.filtered.tsv.')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        merge_readCounts(snakemake.input, snakemake.output[0])
    else:
        import argparse
        args = parse_args()
        if not args.output:
            args.output = args.snp + '.filtered.tsv'
        filter_CNVs(args.snp, args.cells, args.cnv, args.output)
