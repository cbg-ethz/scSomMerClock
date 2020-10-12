#!/usr/bin/env python3

import argparse
import os
import re
import numpy as np
import pandas as pd
from pysam import VariantFile


CHROM = [str(i) for i in range(1, 23, 1)] + ['X', 'Y']
ALG_MAP = {'monovar': 0, 'sccaller': 1, 'mutect': 2}
RES_MAP = {'[1 0 0]': 0, '[0 1 0]': 1, '[0 0 1]': 2, '[1 1 0]': 3, '[1 0 1]': 4,
    '[0 1 1]': 5, '[1 1 1]': 6, '[0 0 0]': False}


def parse_args():
    parser = argparse.ArgumentParser(
        prog='QC_coverage', usage='python3 QC_coverage.py <DATA> [options]',
        description='*** Generate Lorenz curve and Gini coefficient. ***'
    )
    parser.add_argument(
        'input', type=str,  nargs='*',
        help='Absolute or relative path(s) to input VCF file'
    )
    parser.add_argument(
        '-t', '--task', type=str, choices=['summarize', 'merge'],
        default='summarize', help='Scripts task, options are: 1. summarizing the '
        'calls in a vcf file, 2. merging the previouslz generated summaries. '
        'Default = summarize'
    )
    parser.add_argument(
        '-o', '--output', type=str, default='',
        help='Path to the output directory. Default = <INPUT_DIR>.'
    )
    parser.add_argument(
        '-bn', '--bulk_normal', type=str, default='',
        help='Column name of bulk normal. Default = None.'
    )
    parser.add_argument(
        '-bt', '--bulk_tumor', nargs='*', type=str,
        help='Column name of bulk tumor. Default = None.'
    )
    parser.add_argument(
        '-q', '--quality', type=int, default=20,
        help='Minimum quality threshold. Default = 20.'
    )
    parser.add_argument(
        '-r', '--read_depth', type=int, default=10,
        help='Minimum read depth at loci. Default = 10.'
    )

    args = parser.parse_args()
    return args


def get_summary_df(args):
    vcf_in = VariantFile(args.input)

    samples = set([])
    for sample in vcf_in.header.samples:
        sample_detail = sample.split('.')
        sample_name = '.'.join(sample_detail[:-1])
        caller = sample_detail[-1]
        if caller != 'mutect':
            samples.add(sample_name)

    sc_map = {j: i for i, j in enumerate(sorted(samples))}
    for bt in args.bulk_tumor:
        sc_map[bt] = len(sc_map)

    sample_no = len(sc_map)

    # Iterate over rows
    print('Iterating calls - Start')
    if args.chr == 'all':
        germline = []
        data = []
        for chrom in CHROM:
            chr_data = vcf_in.fetch(chrom)
            chr_data = iterate_chrom(chr_data, sc_map, sample_no, chrom)

            data.extend(chr_data[0])
            germline.extend(chr_data[1])
    else:
        chr_data = vcf_in.fetch()
        data, germline = iterate_chrom(chr_data, sc_map, sample_no, args.chr)

    cols = ['CHROM', 'POS', \
        'monovar', 'sccaller', 'bulk', 'monovar_sccaller', 'monovar_bulk', \
        'sccaller_bulk', 'monovar_sccaller_bulk']
    df = pd.DataFrame(data, columns=cols)
    df.set_index(['CHROM', 'POS'], inplace=True)
    df = df.astype(int)

    out_summary = os.path.join(args.output, 'Calls.{}.DP{}_QUAL{}.tsv' \
        .format(args.chr, args.read_depth, args.quality))
    print('Writing call summary to: {}'.format(out_summary))
    df.to_csv(out_summary, sep='\t')

    if germline:
        out_germ = os.path.join(args.output, 'Germline.{}.DP{}_QUAL{}.tsv' \
            .format(args.chr, args.read_depth, args.quality))
        print('Writing germline calls to: {}'.format(out_germ))
        with open(out_germ, 'w') as f:
            f.write('\n'.join(germline))

    return df


def iterate_chrom(chr_data, sc_map, sample_size, chrom):
    data = []
    germline = []

    for i, rec in enumerate(chr_data):
        if i % 100000 == 0:
            print('Iterated records on Chr {}:\t{}'.format(chrom, i))

        # Filtered in bulk & multplie genotypes called by SCcaller
        if not 'PASS' in rec.filter:
            continue

        # 0: monovar, 1: sccaller, 2: bulk_tumor
        calls = np.zeros((sample_size, 3), dtype=int)
        # Iterate over columns (i.e. samples)
        for sample_id, sample in rec.samples.iteritems():
            # Skip all genotypes except: 0/1 | 1/1 | 0/2
            if not sample['GT'][1]:
                continue

            sample_detail = sample_id.split('.')
            if len(sample_detail) == 1:
                sample_name = sample_id
                alg = 'mutect'
            else:
                sample_name = '.'.join(sample_detail[:-1])
                alg = sample_detail[-1]

            if alg == 'mutect':
                if sample['DP'] < args.read_depth:
                    continue
            elif alg == 'sccaller':
                # Skip indels and "False" calls (only called by sccaller)
                # Skip low genotype quality calls or read depth below threshold
                if sample['SO'] != 'True' \
                        or rec.alleles[sample['GT'][1]].startswith(('+', '-')):
                    continue

                if sample['GQ'] < args.quality:
                    if min(sample['FPL'][:2]) - max(sample['FPL'][2:]) \
                            < args.quality:
                        continue
            else:
                # skip low coverage snps
                if sum(sample['AD']) < args.read_depth:
                    continue
                # skip if LL difference between WT and Hom/Het is below threshold
                if sample['GQ'] < args.quality:
                    PL_max = np.max([i for i in sample["PL"][1:] if i != None])
                    if sample['PL'][0] - PL_max < args.quality:
                        continue
                    
            try:
                sample_map_id = sc_map[sample_name]
            # Sample name not in mapping dict: bulk normal
            except KeyError:
                germline.append('{}:{}'.format(rec.chrom, rec.pos))
            else:

                calls[sample_map_id, ALG_MAP[alg]] = 1

        # only WT called
        if calls.max() == 0:
            continue

        call_data = np.zeros(7, dtype=int)
        for sample_calls in calls:
            call_data[RES_MAP[np.array2string(sample_calls)]] += 1

        rec_data = np.append([rec.chrom, rec.pos], call_data)
        data.append(rec_data)

    print('Iterating calls on Chr {} - End\n'.format(chrom))
    return data, germline


def get_summary_statistics(df, args):
    # Add sum _column
    df['sum'] = df.sum(axis=1)

    # Get singletons called either by monovar and/or SCcaller
    singletons = df[
        ((df['sum'] == 1) \
            & (df[['monovar', 'sccaller', 'monovar_sccaller']].sum(axis=1) == 1))
        | ((df['sum'] == 2) & (df['monovar'] == 1) & (df['sccaller'] == 1))
    ]
    df.drop(singletons.index, inplace=True)

    # Get overview of records called by SCCaller and Monovar but not in bulk
    sc_col_order = ['monovar', 'sccaller', 'monovar_sccaller']
    sc_callers = df[(df['bulk'] == 0) & (df['sccaller_bulk'] == 0) & \
        (df['monovar_bulk'] == 0) & (df['monovar_sccaller_bulk'] == 0)].copy()
    sc_callers['sc_unique'] = sc_callers[sc_col_order] \
        .apply(lambda x: ' '.join(str(i) for i in x ), axis=1)

    # Only called in bulk
    b = df[df['bulk'] == df['sum']]
    df.drop(b.index, inplace=True)

    # Get SNPs called only by monovar in min 2 samples
    m = df[(df['monovar'] == df['sum']) \
        | ((df[['monovar', 'sccaller']].sum(axis=1) == df['sum']) \
            & (df['sccaller'] == 1))
    ]
    df.drop(m.index, inplace=True)

    # Get SNPs called only by SCcaller in min 2 samples
    s = df[(df['sccaller'] == df['sum']) \
        | ((df[['monovar', 'sccaller']].sum(axis=1) == df['sum']) \
            & (df['monovar'] == 1))
    ]
    df.drop(s.index, inplace=True)

    # Get SNPs called by monovar and in bulk, but not by sccaller
    m_b = df[
        (((df['monovar'] > 0) & (df['bulk'] > 0)) | (df['monovar_bulk'] > 0)) \
        & (df['sccaller'] == 0) & (df['monovar_sccaller'] == 0) \
        & (df['sccaller_bulk'] == 0) & (df['monovar_sccaller_bulk'] == 0)
    ]
    df.drop(m_b.index, inplace=True)

    # Get SNPs called by sccaller and in bulk, but not by monovar
    s_b = df[
        (((df['sccaller'] > 0) & (df['bulk'] > 0)) | (df['sccaller_bulk'] != 0)) \
        & (df['monovar'] == 0) & (df['monovar_sccaller'] == 0) \
        & (df['monovar_bulk'] == 0) & (df['monovar_sccaller_bulk'] == 0)
    ]
    df.drop(s_b.index, inplace=True)

    # Get SNPs called by monovar and SCcaller, but not in bulk 
    bulk_cols = ['bulk', 'monovar_bulk', 'sccaller_bulk', 'monovar_sccaller_bulk']
    m_s = df[(df['monovar_sccaller'] >= 2) & (df[bulk_cols].sum(axis=1) == 0) \
        | ((df['monovar_sccaller'] == 1) & (df['monovar'] >=  1) \
            & (df['sccaller'] >=  1))
    ]
    df.drop(m_s.index, inplace=True)
    
    # Get SNPs called by both algorithms in SC and in bulk
    m_s_b = df[(df['monovar_sccaller_bulk'] > 0) \
        | ((df['bulk'] > 0) & ((df['monovar_sccaller'] > 0) \
            | ((df['sccaller'] > 0) & (df['monovar'] > 0))))]
    df.drop(m_s_b.index, inplace=True)

    # Get SNPs called by SC algorithms in 2 samples only
    s_shady = df[(df['monovar_sccaller'] == 1) & (df['sccaller'] > 0) \
        & (df['monovar'] == 0)]
    df.drop(s_shady.index, inplace=True)

    m_shady = df[(df['monovar_sccaller'] == 1) & (df['monovar'] > 0) \
        & (df['sccaller'] == 0)]
    df.drop(m_shady.index, inplace=True)

    if not df.empty:
        print('Unknown how to handle:\n{}'.format(df))

    sc_unique = pd.DataFrame(columns=sc_col_order + ['count'])
    for i, j in sc_callers['sc_unique'].value_counts().iteritems():
        sc_unique = sc_unique.append(
            pd.Series(np.append(np.fromstring(i, sep=' ', dtype=int), j),
                index=sc_col_order + ['count']),
            ignore_index=True
        )
    if not sc_unique.empty:
        out_QC_SConly = os.path.join(args.output, 'SConly_summary.{}.DP{}_QUAL{}.tsv' \
            .format(args.chr, args.read_depth, args.quality))
        sc_unique.astype(int).to_csv(out_QC_SConly, sep='\t', index=False)

    rel_recs = pd.concat([m_s_b, m_s, m_b, s_b])
    if not rel_recs.empty:
        rel_recs.drop('sum', axis=1, inplace=True)
        rel_SNPs = os.path.join(args.output, 'relevantSNPs.{}.DP{}_QUAL{}.tsv' \
            .format(args.chr, args.read_depth, args.quality))
        rel_recs.sort_index(ascending=False).to_csv(rel_SNPs, sep='\t')

    data = [
        ('Monovar2', m.shape[0] + m_shady.shape[0]),
        ('SCcaller2', s.shape[0] + s_shady.shape[0]),
        ('Monovar2 & SCcaller2', m_s.shape[0]),
        ('Bulk1', b.shape[0]),
        ('Monovar1 & Bulk1', m_b.shape[0]),
        ('SCcaller1 & Bulk1', s_b.shape[0]),
        ('Monovar1 & SCcaller1 & Bulk1', m_s_b.shape[0])
    ]
    out_QC = os.path.join(args.output, 'Call_summary.{}.DP{}_QUAL{}.tsv' \
        .format(args.chr, args.read_depth, args.quality)
    )
    save_summary(data, out_QC)

    return data


def save_summary(data, out_file, verbose=True):
    if verbose:
        print('Writing call-Venn-data to: {}\n\n\nSUMMARY:'.format(out_file))
    with open(out_file, 'w') as f:
        for alg, calls in data:
            f.write('{}\t{}\n'.format(calls, alg))
            if verbose:
                print('{}\t-\t{}'.format(calls, alg))


def plot_venn(data, out_dir):
    try:
        import matplotlib.pyplot as plt
        from matplotlib_venn import venn3, venn3_circles
    except ImportError:
        return

    # Hack to ensure min circle size
    # ------------------------------
    counts = np.array([i[1] for i in data])
    counts_norm = np.clip(counts / counts.sum(), 0.025, 1).round(2)
    counts_norm = counts_norm + np.arange(0.0001, 0.00071, 0.0001)

    no_formatter = {}
    for i, j in enumerate(counts_norm):
        no_formatter[j] = str(data[i][1])
        # if data[i][0] in ['Monovar', 'SCcaller', 'Monovar & SCcaller']:
        #     no_formatter[j] += r' ($\geq$2 samples)'

    def formatter(x):
        return no_formatter[x]
    # ------------------------------

    fig = plt.figure(figsize=(10, 10))
    v = venn3(
        subsets=(counts_norm),
        set_labels=('Monovar', 'SCcaller', 'Bulk'),
        set_colors=('r', 'g', 'b'),
        alpha = 0.5,
        normalize_to=counts_norm.sum(),
        subset_label_formatter=formatter
    )
    c = venn3_circles(
        subsets = (counts_norm),
        normalize_to=counts_norm.sum()
    )

    for text in v.set_labels:
        text.set_fontsize(14)
    for text in v.subset_labels:
        text.set_fontsize(16)
        
    try:
        fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)
    except AttributeError:
        pass

    out_file = os.path.join(args.output, 'SNP_counts_DP{}_QUAL{}.pdf' \
        .format(args.read_depth, args.quality))
    print('Saving call-Venn-plots to: {}'.format(out_file))
    fig.savefig(out_file, dpi=300)
    plt.close()


def load_data(in_file):
    data = []
    with open(in_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            no, alg = line.strip().split('\t')
            data.append((alg, int(no)))
    return data


def merge_summaries(args):
    counts = np.zeros(7)
    algs = ['Monovar', 'SCcaller', 'Monovar & SCcaller', 'Bulk',
        'Monovar & Bulk', 'SCcaller & Bulk', 'Monovar & SCcaller & Bulk']

    df_snp = pd.DataFrame([])
    df_sc = pd.DataFrame([],
        columns=['monovar', 'sccaller', 'monovar_sccaller', 'count'])
    df_sc.set_index(['monovar', 'sccaller', 'monovar_sccaller'], inplace=True)

    for in_file in args.input:
        with open(in_file, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                counts[i] += int(line.split('\t')[0])
    
        snp_file = in_file.replace('Call_summary', 'relevantSNPs', 1)
        if os.path.exists(snp_file):
            df_snp_new = pd.read_csv(snp_file, sep='\t', index_col=[0,1])
            df_snp = df_snp.append(df_snp_new)

        SConly_file = in_file.replace('Call_summary', 'SConly_summary', 1)
        if os.path.exists(SConly_file):
            df_sc_new = pd.read_csv(SConly_file, sep='\t', index_col=[0, 1, 2])
            for idx, row in df_sc_new.iterrows():
                try:
                    df_sc.loc[idx] += row['count']
                except KeyError:
                    df_sc = df_sc.append(row)

    data = [(algs[i], int(counts[i])) for i in range(7)]
    out_QC = os.path.join(args.output, 'Call_summary.all.DP{}_QUAL{}.tsv' \
        .format(args.read_depth, args.quality)
    )
    save_summary(data, out_QC)

    out_SNP = out_QC.replace('Call_summary', 'relevantSNPs', 1)
    df_snp.sort_index().astype(int).to_csv(out_SNP, sep='\t') 

    out_SConly = out_QC.replace('Call_summary', 'SConly_summary', 1)
    df_sc.sort_values('count', ascending=False).to_csv(out_SConly, sep='\t') 

    return data


def main(args):
    if not args.output:
        args.output = os.path.dirname(args.input[0])
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    if args.task == 'summarize':
        if len(args.input) != 1:
            raise IOError('Can only summarize 1 vcf file ({} given)' \
                .format(len(args.input)))
        args.input = args.input[0]

        try:
            args.chr = re.search('\.[0-9XY]+\.', args.input).group(0).strip('.')
        except AttributeError:
            args.chr = 'all'

        df = get_summary_df(args)
        summary = get_summary_statistics(df, args)

        if args.chr == 'all':
            plot_venn(summary, args)

        # df = pd.read_csv(args.input, sep='\t', index_col=[0,1])
        # df.astype(int, copy=False)
        # summary = load_data(args.input)
        # plot_venn(summary, args)
    else:
        summary = merge_summaries(args)
        args.depth = int(re.search('DP(\d+)_', args.input[0]).group(1))
        args.quality = int(re.search('_QUAL(\d+)\.', args.input[0]).group(1))
        plot_venn(summary, args)


if __name__ == '__main__':
    args = parse_args()
    main(args)
