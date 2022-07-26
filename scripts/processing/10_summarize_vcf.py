#!/usr/bin/env python3

import argparse
import os
import re
import time
import numpy as np
import pandas as pd
from collections import OrderedDict
from pysam import VariantFile


CHROM = [str(i) for i in range(1, 23, 1)] + ['X', 'Y']
ALG_MAP = {'monovar': 0, 'sccaller': 1, 'mutect': 2}
MIN_READS_NAN = 5
MIN_GQ_NAN = 1


VCF_HEADER = """##fileformat=VCFv4.1
##fileDate={time.tm_year}:{time.tm_mon}:{time.tm_mday}-{time.tm_hour}:{time.tm_min}:{time.tm_sec}
##source=MolecularClockTesting_pipeline
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
{ref}
{contigs}
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}
"""


def get_summary_df(input, chr, output, read_depth, quality, bulk_tumor=[],
            prefix='', keep_sex=False, gt_sep=',', out_nexus=False):
    vcf_in = VariantFile(input)
    in_file = os.path.basename(input)

    samples = set([])
    bulk_samples = set([])
    for sample in vcf_in.header.samples:
        sample_detail = sample.split('.')
        sample_name = '.'.join(sample_detail[:-1])
        caller = sample_detail[-1]
        if caller == 'mutect':
            print(f'\tBulk samples: {sample_name}')
            bulk_samples.add(sample_name)
        else:
            samples.add(sample_name)

    sc_map = OrderedDict([(j, i) for i, j in enumerate(sorted(samples))])
    bk_map = OrderedDict()
    if bulk_tumor:
        for i, j in enumerate(sorted(bulk_tumor)):
            if j in bulk_samples:
                bk_map[j] = i
            else:
                print(f'\nWARNING: Cannot find bulk sample: {j}\n')

        for i in bulk_samples - bk_map.keys():
            print(f'\t Bulk normal sample: {i}')

    sample_maps = (sc_map, bk_map)

    # Iterate over rows
    print('Iterating calls - Start')
    if chr == 'all_chr':
        data = {'singletons': [], 
            'monovar2+': [],
            'sccaller2+': [],
            'bulk': [],
            'monovar2+_sccaller2+': [],
            'monovar1+_bulk': [],
            'sccaller1+_bulk': [],
            'monovar1+_sccaller1+_bulk': []
        }
        vcf_body = ''

        germline = []
        for chrom in CHROM:
            if not keep_sex and chrom in ['X', 'Y']:
                continue

            if prefix != '':
                chrom = prefix + chrom

            try:
                chr_data_in = vcf_in.fetch(chrom)
            except ValueError:
                raise RuntimeError(f'Missing index for file: {vcf_file}')
            chr_data, chr_vcf_body, chr_germline = iterate_chrom(
                chr_data_in, sample_maps, chrom, read_depth, quality)

            for group in chr_data:
                data[group].extend(chr_data[group])

            vcf_body += chr_vcf_body
            germline.extend(chr_germline)
    else:
        chr_data_in = vcf_in.fetch(chr)
        data, vcf_body, germline = iterate_chrom(
            chr_data_in, sample_maps, chr, read_depth, quality)

    cols = ['CHROM', 'POS', 'monovar', 'sccaller', 'bulk', 'monovar_sccaller',
        'monovar_bulk', 'sccaller_bulk', 'monovar_sccaller_bulk']
    df = pd.DataFrame([])
    for group_vals in data.values():
        df = df.append(group_vals)
    df.columns = cols
    df.set_index(['CHROM', 'POS'], inplace=True)
    df = df.astype(int)

    out_summary = os.path.join(output, 'Call_details.{}.tsv'.format(chr))
    print('Writing call summary to: {}'.format(out_summary))
    df.to_csv(out_summary, sep='\t')

    contigs = '\n'.join(['##contig=<ID={},eta=-1>'.format(i) for i in CHROM])
    samples = '\t'.join(sc_map.keys())
    ref = re.search('##reference=.*\n', str(vcf_in.header))[0].rstrip('\n')
    vcf_header = VCF_HEADER.format(time=time.localtime(), contigs=contigs,
        ref=ref, samples=samples)

    out_vcf = os.path.join(output, 'all_filtered.{}.vcf'.format(chr))
    print('Writing vcf file to: {}'.format(out_vcf))
    with open(out_vcf, 'w') as f_vcf:
        f_vcf.write(vcf_header)
        f_vcf.write(vcf_body.strip('\n'))

    if germline:
        out_germ = os.path.join(output, 'Call_germline.{}.tsv'.format(chr))
        print('Writing germline calls to: {}'.format(out_germ))
        with open(out_germ, 'w') as f:
            f.write('\n'.join(germline))

    out_QC = os.path.join(output, 'Call_summary.{}.tsv'.format(chr))
    save_summary(data, out_QC)

    return data


def iterate_chrom(chr_data, sample_maps, chrom, read_depth, quality):
    data = {'singletons': [], 
        'monovar2+': [],
        'sccaller2+': [],
        'bulk': [],
        'monovar2+_sccaller2+': [],
        'monovar1+_bulk': [],
        'sccaller1+_bulk': [],
        'monovar1+_sccaller1+_bulk': []
    }
    out_vcf = ''
    germline = []

    for idx, rec in enumerate(chr_data):
        if idx % 100000 == 0:
            print('Iterated records on Chr {}:\t{}'.format(chrom, idx))

        # Filtered in bulk and/or multplie genotypes called by SCcaller
        if not 'PASS' in rec.filter:
            continue

        sc_calls, is_bulk_snv, is_germline_snv = get_rec_summary(rec,
            sample_maps, read_depth, quality)

        if is_germline_snv:
            germline.append('{}:{}'.format(rec.chrom, rec.pos))
            continue
        # No mutation in single cells
        if sc_calls.max() <= 0:
            if is_bulk_snv:
                data['bulk'].append([rec.chrom, rec.pos, 0, 0, 1, 0, 0, 0, 0])
            continue

        # 0: monovar, 1: sccaller, 2: bulk_tumor
        monovar_only = np.sum((sc_calls[0] == 1) & (sc_calls[1] != 1))
        sccaller_only = np.sum((sc_calls[0] != 1) & (sc_calls[1] == 1))
        monovar_sccaller = np.sum((sc_calls[0] == 1) & (sc_calls[1] == 1))

        # SNV also called in bulk
        if is_bulk_snv:
            rec_data = (rec.chrom, rec.pos,
                    0, 0, 1, 0, monovar_only, sccaller_only, monovar_sccaller)
            # Detected by both algorithms and in bulk
            if (monovar_sccaller > 0) or (monovar_only >= 1 and sccaller_only >= 1):
                data['monovar1+_sccaller1+_bulk'].append(rec_data)
                rec_vcf = get_call_output(rec, sc_calls, sample_maps[0])
            # Detected by monovar and in bulk
            elif monovar_only >= 1 and sccaller_only == 0:
                data['monovar1+_bulk'].append(rec_data)
                rec_vcf = get_call_output(rec, sc_calls, sample_maps[0])
            # Detected by sccaller and in bulk
            elif sccaller_only >= 1 and monovar_only == 0:
                data['sccaller1+_bulk'].append(rec_data)
                rec_vcf = get_call_output(rec, sc_calls, sample_maps[0])
            else:
                import pdb; pdb.set_trace()
        # SNV only called in SC
        else:
            rec_data = (rec.chrom, rec.pos,
                monovar_only, sccaller_only, 0, monovar_sccaller, 0, 0, 0)
            no_snvs = sum(rec_data[2:])

            if no_snvs == 1 \
                    or (no_snvs == 2 and monovar_only == 1 and sccaller_only == 1):
                data['singletons'].append(rec_data)
                continue
            elif monovar_sccaller != 0:
                if monovar_sccaller > 1 \
                        or (monovar_only > 0 and sccaller_only > 0):
                    data['monovar2+_sccaller2+'].append(rec_data)
                    rec_vcf = get_call_output(rec, sc_calls, sample_maps[0])
                elif monovar_only > 0:
                    data['monovar2+'].append(rec_data)
                    continue
                elif sccaller_only > 0:
                    data['sccaller2+'].append(rec_data)
                    continue
                else:
                    import pdb; pdb.set_trace()
            elif sccaller_only <= 1 and monovar_only >= 2:
                data['monovar2+'].append(rec_data)
                continue
            elif monovar_only <= 1 and sccaller_only >= 2:
                data['sccaller2+'].append(rec_data)
                continue
            elif monovar_only >= 2 and sccaller_only >= 2:
                data['monovar2+_sccaller2+'].append(rec_data)
                rec_vcf = get_call_output(rec, sc_calls, sample_maps[0])
            else:
                import pdb; pdb.set_trace()
            
        out_vcf += rec_vcf

    print('Iterating calls on Chr {} - End\n'.format(chrom))
    return data, out_vcf, germline


def get_rec_summary(rec, sample_maps, depth, quality):
    dims = (2, len(sample_maps[0])) # Dim 0: Monovar, Dim 1: SCcaller
    sc_calls = np.full(dims, -1, dtype=int) # -1: nan, 0: wildtype, 1: mutation
    DP = np.zeros(dims, dtype=int)
    GQ = np.zeros(dims, dtype=int)
    snv_bulk = False
    snv_germline = False

    for sample_id, sample in rec.samples.iteritems():
        alg, sample_map_id = get_call_ids(sample_id, sample_maps)

        # Skip missing data
        if sample['GT'][0] == None:
            continue

        # Set indels and "False" calls (only called by sccaller) to NAN
        if alg == 1 and sample['SO'] != 'True':
            continue
        # Mutect calls
        elif alg == 2:
            if sample['GT'] != (0, 0):
                if sample_map_id == 'normal':
                    snv_germline = True
                else:
                    snv_bulk = True
            continue

        # Transform PL values to: wt|het|hom with min. value of 0
        if 'FPL' in sample and sample['FPL'][0] != None:
            wt_PL = min(sample['FPL'][0], sample['FPL'][1])
            PL_vals = np.array([wt_PL, sample['FPL'][2], sample['FPL'][3]], dtype=int)
        else:
            PL_vals = np.array([j for j in sample['PL'] if j != None], dtype=int)
        PL_vals -= PL_vals.min()
        GQ_s = sorted(PL_vals)[1]

        # Weird pysam/SCcaller output error (6 instead of 3 PL fields)
        try:
            sample['PL'] = [int(i) for i in PL_vals]
        except TypeError:
            try:
                sample['PL'] = [int(i) for i in PL_vals] + [None] * 3
            except TypeError:
                sample['PL'] = [int(i) for i in PL_vals] + [None] * 7
        sample['GQ'] = int(GQ_s)

        # Correct if GQ difference is between 0/1 and 1/1
        if GQ_s < quality and (np.abs(PL_vals[1] - PL_vals[2]) == GQ_s):
            wt_mut_diff = int(PL_vals[0] - GQ_s)
            if wt_mut_diff >= MIN_GQ_NAN:
                GQ_s = wt_mut_diff
                sample['GQ'] = wt_mut_diff

        DP_s = sum(sample['AD'])

        GQ[alg, sample_map_id] = GQ_s
        DP[alg, sample_map_id] = DP_s

        # Set genotype quality or read depth below minimum threshold to NAN
        if DP_s < MIN_READS_NAN or GQ_s < MIN_GQ_NAN:
            continue

        # Save genotype: 0 = wt, 1 = mut
        if sample['GT'][0] == 0 and sample['GT'][1] == 0:
            sc_calls[alg, sample_map_id] = 0
        else:
            sc_calls[alg, sample_map_id] = 1

    rel_calls = (sc_calls == 1) & (DP >= depth) & (GQ >= quality)
    # Check if at least 1 good mutation call
    if rel_calls.sum() == 0:
        sc_calls[:] = -1
    # If no bulk call, check at least two calls fullfill the DP and GQ criteria
    elif not snv_bulk:
        if not any(rel_calls.sum(axis=1) > 1) \
                and not any(rel_calls.sum(axis=0) == 2):
            sc_calls[:] = -1

    return sc_calls, snv_bulk, snv_germline


def get_call_ids(sample_id, sample_maps):
    sample_detail = sample_id.split('.')
    if len(sample_detail) == 1:
        sample_name = sample_id
        alg = 2
    else:
        sample_name = '.'.join(sample_detail[:-1])
        alg = ALG_MAP[sample_detail[-1]]

    if alg == 2:
        try:
            sample_map_id = sample_maps[1][sample_name]
        # Sample name not in mapping dict: bulk normal
        except KeyError:
            sample_map_id = 'normal'
    else:
        sample_map_id = sample_maps[0][sample_name]
    return alg, sample_map_id


def get_call_output(rec, calls, sc_map):
    best_calls = calls.argmax(axis=0)
    data_calls = np.sum(calls.max(axis=0) > -1)

    if len(rec.alts) == 1:
        alt = rec.alts[0]
    else:
        call_alls = [i['GT'][1] for i in rec.samples.values() if i['GT'][1]]
        alls, alls_count = np.unique(call_alls, return_counts=True)
        alt = rec.alts[alls[np.argmax(alls_count)] - 1]

    rec_out = f'\n{rec.chrom}\t{rec.pos}\t.\t{rec.ref}\t{alt}\t' \
        f'{min(99, rec.qual)}\tPASS\tNS={data_calls}\tGT:DP:AD:GQ:PL'
    
    for sample, sample_id in sc_map.items():
        # Prioritize SCcaller calls over monovar calls
        if calls[1, sample_id] > -1:
            alg = 'sccaller'
            gt = calls[1, sample_id]
        else:
            alg = 'monovar'
            gt = calls[0, sample_id]

        if gt == -1:
            rec_out += '\t./.:.:.:.:.'
        else:
            call = rec.samples[f'{sample}.{alg}']

            rec_out += f'\t{min(1, call["GT"][0])}/{min(1, call["GT"][1])}:' \
                f'{np.sum(call["AD"])}:{call["AD"][0]},{call["AD"][1]}:' \
                f'{call["GQ"]}:{call["PL"][0]},{call["PL"][1]},{call["PL"][2]}'

    return rec_out


def get_gt_mat_entry(gt):
    if gt[0] == 0 and gt[1] == 0:
        return '0'
    elif gt[0] == 0 and gt[1] != 0:
        return '1'
    else:
        return '2'


def save_summary(data, out_file, verbose=True):
    if verbose:
        print('Writing call-Venn-data to: {}\n\n\nSUMMARY:'.format(out_file))
    with open(out_file, 'w') as f:
        for alg, calls in data.items():
            if alg == 'singletons':
                continue
            alg_str = alg.replace('_',' & ')
            if isinstance(calls, list):
                call_no = len(calls)
            else:
                call_no = calls
            f.write('{}\t{}\n'.format(call_no, alg))
            if verbose:
                print('{}\t-\t{}'.format(call_no, alg_str))


def plot_venn(data, output):
    try:
        import matplotlib.pyplot as plt
        from matplotlib_venn import venn3, venn3_circles
    except ImportError:
        return

    # Hack to ensure min circle size
    # ------------------------------
    if isinstance(data['monovar2+'], list):
        counts = np.array([
            len(data['monovar2+']),
            len(data['sccaller2+']),
            len(data['monovar2+_sccaller2+']),
            len(data['bulk']),
            len(data['monovar1+_bulk']),
            len(data['sccaller1+_bulk']),
            len(data['monovar1+_sccaller1+_bulk']),
        ])
    else:
        counts = np.array([data['monovar2+'], data['sccaller2+'],
            data['monovar2+_sccaller2+'], data['bulk'], data['monovar1+_bulk'],
            data['sccaller1+_bulk'], data['monovar1+_sccaller1+_bulk']])
    counts_norm = np.clip(counts / counts.sum(), 0.025, 1).round(2)
    counts_norm = counts_norm + np.arange(0.0001, 0.00071, 0.0001)

    venn_labels = [i.replace('_', ' & ') for i in data.keys() if i != 'singletons']
    no_formatter = {}
    for i, j in enumerate(counts_norm):
        no_formatter[j] = counts[i]

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

    out_file = os.path.join(output, 'SNP_counts.pdf')
    print('Saving call-Venn-plots to: {}'.format(out_file))
    fig.savefig(out_file, dpi=300)
    plt.close()


def parse_args():
    parser = argparse.ArgumentParser(
        description='Merge and filter vcf files from different SNVcallers.')
    parser.add_argument('input', type=str,
        help='Absolute or relative path(s) to input VCF file')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Path to the output directory. Default = <INPUT_DIR>.')
    parser.add_argument('-bt', '--bulk_tumor', nargs='*', type=str,
        help='Column name of bulk tumor. Default = None.')
    parser.add_argument('-ks', '--keep_sex', action='store_true',
        help='If flag is set, keep the sex chromosomes in the all file.')
    parser.add_argument('-q', '--quality', type=int, default=13,
        help='Minimum quality threshold. Default = 13.')
    parser.add_argument('-r', '--read_depth', type=int, default=10,
        help='Minimum read depth at loci. Default = 10.')
    parser.add_argument('-nq', '--nan_quality', type=int, default=1,
        help='Min. locus GQ to be not reported as GQ. Default = 1.')
    parser.add_argument('-nr', '--nan_read_depth', type=int, default=5,
        help='Min. locus read depth to be not reported as missing. Default = 5.')
    parser.add_argument('-p', '--prefix', type=str, default='',
        help='Prefix for chromosome, e.g. "chr". Default = "".')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        get_summary_df(
            input=snakemake.input,
            chr=snakemake.wildcards.chr,
            output=os.path.dirname(snakemake.output[0]),
            read_depth=snakemake.params.filter_DP,
            quality=snakemake.params.filter_QUAL,
            bulk_tumor=snakemake.params.bulk_tumor,
            prefix=snakemake.params.pref,
        )
    else:
        args = parse_args()

        if not args.output:
            args.output = os.path.dirname(args.input)
        if not os.path.exists(args.output):
            os.makedirs(args.output)

        MIN_READS_NAN = args.nan_read_depth
        MIN_GQ_NAN = args.nan_quality

        try:
            args.chr = re.search('\.[0-9XY]+\.', args.input).group(0).strip('.')
        except AttributeError:
            args.chr = 'all_chr'

        data = get_summary_df(
            args.input,
            args.chr,
            args.output,
            args.read_depth,
            args.quality,
            bulk_tumor=args.bulk_tumor,
            prefix=args.prefix,
            keep_sex=args.keep_sex,
        )
        if args.chr == 'all_chr':
            plot_venn(data, args.output)

