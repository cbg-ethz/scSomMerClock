#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
import vcf

## SCcaller FORMAT
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=BI,Number=1,Type=Float,Description="Amplification Bias">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=FPL,Number=4,Type=Integer,Description="sequencing noise, amplification artifact, heterozygous SNV and homozygous SNV respectively">
##FORMAT=<ID=SO,Number=1,Type=String,Description="Whether it is a somatic mutation.">

## Monovar FORMAT
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">

## Mutect2 FORMAT
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles in the tumor">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=F1R2,Number=R,Type=Integer,Description="Count of reads in F1R2 pair orientation supporting each allele">
##FORMAT=<ID=F2R1,Number=R,Type=Integer,Description="Count of reads in F2R1 pair orientation supporting each allele">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">
##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phasing set (typically the position of the first variant in the set)">
##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">


def parse_args():
    parser = argparse.ArgumentParser(
        prog='QC_coverage', usage='python3 QC_coverage.py <DATA> [options]',
        description='*** Generate Lorenz curve and Gini coefficient. ***'
    )
    parser.add_argument(
        'input', type=str,
        help='Absolute or relative path(s) to input VCF file'
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
        '-bt', '--bulk_tumor', nargs='+', type=str,
        help='Column name of bulk tumor. Default = None.'
    )
    parser.add_argument(
        '-q', '--quality', type=int, default=30,
        help='Minimum quality threshold. Default = 30.'
    )
    parser.add_argument(
        '-r', '--read_depth', type=int, default=5,
        help='Minimum read depth at loci. Default = 5.'
    )
    parser.add_argument(
        '-gq', '--genotype_quality', type=int, default=0,
        help='Minimum genotype quality. Default = 0.'
    )

    args = parser.parse_args()
    return args


def get_summary_df(args):
    file_name = os.path.basename(args.input)
    vcf_reader = vcf.Reader(filename=args.input)

    samples = set([])
    for sample in vcf_reader.samples:
        sample_detail = sample.split('.')
        sample_name = '.'.join(sample_detail[:-1])
        caller = sample_detail[-1]
        if caller != 'mutect':
            samples.add(sample_name)
    # Ni8 specific
    samples.discard('P01M01E')
    samples.discard('P01P01E')

    sc_map = {j: i for i, j in enumerate(sorted(samples))}
    for bt in args.bulk_tumor:
        sc_map[bt] = len(sc_map)
    sc_size = len(sc_map)

    alg_map = {'monovar': 0, 'sccaller': 1}
    res_map = {'[1 0 0]': 0, '[0 1 0]': 1, '[0 0 1]': 2, '[1 1 0]': 3, 
        '[1 0 1]': 4, '[0 1 1]': 5, '[1 1 1]': 6}

    germline = []
    data = []
    indels = 0
    # Iterate over rows
    for i, record in enumerate(vcf_reader):
        # Skip indels (only keep snp)
        if record.var_type == 'indel':
            indels += 1
            continue

        # Only called in SC (and not in bulk) but with low quality
        try:
            if record.QUAL < args.quality and record.FILTER == None:
                continue
        except TypeError:
            # Only called in bulk contains filter
            if record.FILTER:
                continue

        # 0: monovar, 1: sccaller, 2: bulk_tumor
        calls = np.zeros((sc_size, 3), dtype=int)
        # Iterate over columns (i.e. samples)
        for sample in [i for i in record.samples if i.called]:
            # WT genotype (called by SCcaller)
            #   0 = 0/0, 1 = 0/1, 2 = 1/1
            if sample.gt_type == 0:
                continue

            sample_detail = sample.sample.split('.')
            sample_name = '.'.join(sample_detail[:-1])
            alg = sample_detail[-1]

            # Ni8 specific
            if sample_name in ['P01M01E', 'P01P01E']:
                continue

            # Skip low quality calls
            if alg == 'mutect':
                if record.FILTER:
                    continue
            else:
                if alg == 'sccaller' and sample.data != 'True':
                    continue

                if sample.data.GQ < args.genotype_quality:
                    continue
                # Skip samples with read depth below threshold
                if sum(sample.data.AD) < args.read_depth:
                    continue

            # Bulk SNV
            if alg == 'mutect':
                # Called in Normal (germline)
                if sample_name == args.bulk_normal:
                    germline.append('{}:{}'.format(record.CHROM, record.POS))
                # Called in tumor
                else:
                    calls[sample_id, 2] = 1
            # SC SNV
            else:
                sample_id = sc_map[sample_name]
                calls[sample_id, alg_map[alg]] = 1

        per_sample = np.sum(calls, axis=1)
        # WT called (by SCCaller)
        if per_sample.max() == 0:
            continue
        # All SNVs only called by 1 algorithm
        elif per_sample.max() == 1:
            call_data = np.append(calls.sum(axis=0), np.zeros(4, dtype=int))
        else:
            call_data = np.zeros(7, dtype=int)
            for sample_calls in calls:
                if sample_calls.sum() != 0:
                    call_data[res_map[np.array2string(sample_calls)]] += 1
        
        rec_data = np.append([record.CHROM, record.POS], call_data)
        data.append(rec_data)

    cols = ['CHROM', 'POS', \
        'monovar', 'sccaller', 'bulk', 'monovar_sccaller', 'monovar_bulk', \
        'sccaller_bulk', 'monovar_sccaller_bulk']
    df = pd.DataFrame(data, columns=cols)
    df.set_index(['CHROM', 'POS'], inplace=True)
    df = df.astype(int)

    out_summary = os.path.join(args.output, 'filtered_DP{}_QUAL{}_summary.{}.tsv' \
        .format(args.read_depth, args.quality, file_name))
    df.to_csv(out_summary, sep='\t')

    if germline:
        out_germ = os.path.join(args.output, 'germline_muts_DP{}_QUAL{}.{}.tsv' \
            .format(args.read_depth, args.quality, file_name))
        with open(out_germ, 'w') as f:
            f.write('\n'.join(germline))

    print('\n(Removed {} indels)\n'.format(indels))

    return df


def get_summary_statistics(df):
    # Add sum _column
    df['sum'] = df.sum(axis=1)

    # Get singleton called either by monovar or SCcaller
    singletons_single = df[
        ((df['sum'] == 1) & (df[['monovar', 'sccaller']].sum(axis=1) == 1)) \
            | ((df['sum'] == 2) & (df[['monovar', 'sccaller']].sum(axis=1) == 2))
    ]
    df.drop(singletons_single.index, inplace=True)

    # Get singleton called by monovar and SCcaller
    singletons_both = df[(df['sum'] == 1) & (df['monovar_sccaller'] == 1)]
    df.drop(singletons_both.index, inplace=True)

    # Get singelton bulk
    b = df[df['bulk'] == df['sum']]
    df.drop(b.index, inplace=True)

    # Get SNPs called only by monovar in min 2 samples
    m = df[(df['monovar'] >= 2) & (df['sccaller'] < 2) \
        & (df['monovar_sccaller'] < 2)]
    df.drop(m.index, inplace=True)

    # Get SNPs called only by SCcaller in min 2 samples
    s = df[(df['sccaller'] >= 2) & (df['monovar'] < 2) \
        & (df['monovar_sccaller'] < 2)]
    df.drop(s.index, inplace=True)

    # Get SNPs called by monovar and in bulk, but not by sccaller
    m_b = df[(df['monovar_bulk'] != 0) & (df['sccaller'] < 2) \
        & (df['monovar_sccaller'] < 2) & (df['sccaller_bulk'] == 0) \
        & (df['monovar_sccaller_bulk'] == 0)]
    df.drop(m_b.index, inplace=True)

    # Get SNPs called by sccaller and in bulk, but not by monovar
    s_b = df[(df['sccaller_bulk'] != 0) & (df['monovar'] < 2) \
        & (df['monovar_sccaller'] < 2) & (df['monovar_bulk'] == 0) \
        & (df['monovar_sccaller_bulk'] == 0)]
    df.drop(s_b.index, inplace=True)

    # Get SNPs called by monovar and SCcaller in min 2 samples, but not in bulk 
    m_s = df[(df['monovar_sccaller'] >= 2) & (df['monovar_sccaller_bulk'] == 0) ]
    df.drop(m_s.index, inplace=True)
    
    # Get SNPs called by both algorithms in SC and in bulk
    m_s_b = df[df['monovar_sccaller_bulk'] != 0]
    df.drop(m_s_b.index, inplace=True)

    # Get SNPs called by SC algorithms in 2 samples only
    s_shady = df[(df['sccaller'] == 1) & (df['monovar_sccaller'] == 1)\
        & (df['sum'] == 2)]
    df.drop(s_shady.index, inplace=True)

    m_shady = df[(df['monovar'] == 1) & (df['monovar_sccaller'] == 1)\
        & (df['sum'] == 2)]
    df.drop(m_shady.index, inplace=True)

    m_s_shady = df[(df['monovar'] == 1) & (df['sccaller'] == 1) \
        & (df['monovar_sccaller'] == 1)]
    df.drop(m_s_shady.index, inplace=True)

    if not df.empty:
        print('Unknown how to handle:\n{}'.format(df))

    data = [
        ('Monovar', m.shape[0] + m_shady.shape[0]),
        ('SCcaller', s.shape[0] + s_shady.shape[0]),
        ('Monovar & SCcaller', m_s.shape[0] + m_s_shady.shape[0]),
        ('Bulk', b.shape[0]),
        ('Monovar & Bulk', m_b.shape[0]),
        ('SCcaller & Bulk', s_b.shape[0]),
        ('Monovar & SCcaller & Bulk', m_s_b.shape[0])
    ]

    print('\n\nSUMMARY:\n')
    with open(os.path.join(args.output, 'QC_calling.tsv'), 'w') as f:
        for alg, calls in data:
            f.write('{}\t{}\n'.format(calls, alg))
            print('{}\t-\t{}'.format(calls, alg))

    return data


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
        if data[i][0] in ['Monovar', 'SCcaller', 'Monovar & SCcaller']:
            no_formatter[j] += r' ($\geq$2 samples)'

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

    out_file = os.path.join(args.output, 'SNP_counts_DP{}_QUAL{}.{}.pdf' \
        .format(args.read_depth, args.quality, os.path.basename(args.input)))
    fig.savefig(out_file, dpi=300)
    plt.close()


def main(args):
    if not args.output:
        args.output = os.path.dirname(args.input)

    df = get_summary_df(args)
    summary = get_summary_statistics(df)
    plot_venn(summary, args)


if __name__ == '__main__':
    args = parse_args()
    main(args)
