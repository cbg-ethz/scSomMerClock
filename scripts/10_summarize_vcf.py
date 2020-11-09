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
SC_COLS = [0, 1, 3]
BULK_COLS = [2, 4, 5, 6]


VCF_HEADER = """##fileformat=VCFv4.1
##fileDate={time.tm_year}:{time.tm_mon}:{time.tm_mday}-{time.tm_hour}:{time.tm_min}:{time.tm_sec}
##source=MolecularClockTesting_pipeline
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
{ref}
{contigs}
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}
"""

NEXUS_TEMPLATE = """#NEXUS
Begin TAXA;
  Dimensions ntax={sample_no};
  TaxLabels {sample_labels}
End;

Begin data;
  Dimensions nchar={rec_no};
  Format datatype=dna missing=? gap=-;
  Matrix
{matrix}
  ;
End;

BEGIN TREES;
END;
"""


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
    parser.add_argument(
        '-s', '--gt_sep', type=str, default=',',
        help='Separator for genotype matrix. Default = ",".'
    )

    args = parser.parse_args()
    return args


def get_summary_df(args):
    vcf_in = VariantFile(args.input)
    in_file = os.path.basename(args.input)

    samples = set([])
    for sample in vcf_in.header.samples:
        sample_detail = sample.split('.')
        sample_name = '.'.join(sample_detail[:-1])
        caller = sample_detail[-1]
        if caller != 'mutect':
            samples.add(sample_name)

    sc_map = OrderedDict([(j, i) for i, j in enumerate(sorted(samples))])
    bk_map = OrderedDict([(j, i) for i, j in enumerate(sorted(args.bulk_tumor))])
    sample_maps = (sc_map, bk_map)

    # Iterate over rows
    print('Iterating calls - Start')
    if args.chr == 'all_chr':
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
        gt_mat = ''
        germline = []
        for chrom in CHROM:
            chr_data_in = vcf_in.fetch(chrom)
            chr_data, chr_vcf_body, chr_gt_mat, chr_all_mat, chr_germline = \
                iterate_chrom(chr_data_in, sample_maps, chrom, args.gt_sep)

            for group in chr_data:
                data[group].extend(chr_data[group])
            try:
                all_mat = np.concatenate([all_mat, chr_all_mat], axis=1)
            except NameError:
                all_mat = chr_all_mat
            except ValueError:
                pass
            vcf_body += chr_vcf_body
            gt_mat += chr_gt_mat
            germline.extend(chr_germline)
    else:
        chr_data = vcf_in.fetch(args.chr)
        data, vcf_body, gt_mat, all_mat, germline = \
            iterate_chrom(chr_data, sample_maps, args.chr, args.gt_sep)

    cols = ['CHROM', 'POS', \
        'monovar', 'sccaller', 'bulk', 'monovar_sccaller', 'monovar_bulk', \
        'sccaller_bulk', 'monovar_sccaller_bulk']
    df = pd.DataFrame([])
    for group_vals in data.values():
        df = df.append(group_vals)
    df.columns = cols
    df.set_index(['CHROM', 'POS'], inplace=True)
    df = df.astype(int)

    out_summary = os.path.join(args.output, 'Call_details.{}.tsv' \
        .format(args.chr))
    print('Writing call summary to: {}'.format(out_summary))
    df.to_csv(out_summary, sep='\t')


    contigs = '\n'.join(['##contig=<ID={},eta=-1>'.format(i) for i in CHROM])
    samples = '\t'.join(sc_map.keys())
    ref = re.search('##reference=.*\n', str(vcf_in.header))[0].rstrip('\n')
    vcf_header = VCF_HEADER.format(time=time.localtime(), 
        contigs=contigs, ref=ref, samples=samples)

    out_vcf = os.path.join(args.output, 'all.{}.filtered.vcf'.format(args.chr))
    print('Writing vcf file to: {}'.format(out_vcf))
    with open(out_vcf, 'w') as f_vcf:
        f_vcf.write(vcf_header)
        f_vcf.write(vcf_body.strip('\n'))

    out_gt = os.path.join(args.output, 'Genotype_matrix.{}.csv'.format(args.chr))
    print('Writing genotype matrix to: {}'.format(out_gt))
    with open(out_gt, 'w') as f_gt:
        f_gt.write('chrom:pos{}{}' \
            .format(args.gt_sep, args.gt_sep.join(sample_maps[0].keys())))
        f_gt.write(gt_mat.rstrip('\n'))


    out_nexus = os.path.join(args.output, 'Genotype_matrix.{}.nex'.format(args.chr))
    print('Writing NEXUS file to: {}'.format(out_nexus))
    nex_labels = ['REF'] + list(sc_map.keys())
    nex_matrix = ''
    if isinstance(all_mat, bool):
        rec_no = 0
    else:
        for i, all_row in enumerate(all_mat):
            nex_matrix += '\t{}\t{}\n'.format(nex_labels[i], ''.join(all_row))
        rec_no = all_mat.shape[1]
    with open(out_nexus, 'w') as f_nex:
        f_nex.write(NEXUS_TEMPLATE.format(sample_no=len(sc_map) + 1,
            sample_labels=' '.join(nex_labels), rec_no=rec_no,
            matrix=nex_matrix.strip('\n')))

    if germline:
        out_germ = os.path.join(args.output, 'Call_germline.{}.tsv'.format(args.chr))
        print('Writing germline calls to: {}'.format(out_germ))
        with open(out_germ, 'w') as f:
            f.write('\n'.join(germline))

    out_QC = os.path.join(args.output, 'Call_summary.{}.tsv'.format(args.chr))
    save_summary(data, out_QC)

    return data


def iterate_chrom(chr_data, sample_maps, chrom, sep=','):
    out_vcf = ''
    gt_mat = ''
    all_mat = []
    data = {'singletons': [], 
        'monovar2+': [],
        'sccaller2+': [],
        'bulk': [],
        'monovar2+_sccaller2+': [],
        'monovar1+_bulk': [],
        'sccaller1+_bulk': [],
        'monovar1+_sccaller1+_bulk': []
    }
    germline = []
    filtered = 0
    only_wt = 0

    for idx, rec in enumerate(chr_data):
        if idx % 100000 == 0:
            print('Iterated records on Chr {}:\t{}'.format(chrom, idx))

        # Filtered in bulk & multplie genotypes called by SCcaller
        if not 'PASS' in rec.filter:
            filtered += 1
            continue

        sc_calls, is_bulk_snv, is_germline_snv = get_call_summary(rec.samples,
            sample_maps, args.read_depth, args.quality)
        if is_germline_snv:
            germline.append('{}:{}'.format(rec.chrom, rec.pos))
            continue
        if sc_calls.max() <= 0:
            continue

        # 0: monovar, 1: sccaller, 2: bulk_tumor
        monovar_only = np.sum((sc_calls[:,0] == 1) & (sc_calls[:,1] != 1))
        sccaller_only = np.sum((sc_calls[:,0] != 1) & (sc_calls[:,1] == 1))
        monovar_sccaller = np.sum((sc_calls[:,0] == 1) & (sc_calls[:,1] == 1))
        # SNV also called in bulk
        if is_bulk_snv:
            rec_data = [rec.chrom, rec.pos,
                    0, 0, 1, 0, monovar_only, sccaller_only, monovar_sccaller]
            # Only detected in bulk:
            if sum(rec_data[2:]) == 1:
                data['bulk'].append(rec_data)
                continue
            # Detected by both algorithms and in bulk
            if (monovar_sccaller > 0) or (monovar_only > 1 and sccaller_only > 1):
                data['monovar1+_sccaller1+_bulk'].append(rec_data)
                rec_vcf, gt_row = get_call_output(rec, sc_calls, sample_maps[0])
            # Detected by monovar and in bulk
            elif monovar_only > 1 and sccaller_only == 0:
                data['monovar1+_bulk'].append(rec_data)
                rec_vcf, gt_row = get_call_output(rec, sc_calls, sample_maps[0])
            # Detected by sccaller and in bulk
            elif sccaller_only > 1 and sccaller_only == 0:
                data['sccaller1+_bulk'].append(rec_data)
                rec_vcf, gt_row = get_call_output(rec, sc_calls, sample_maps[0])
        # SNV only called in SC
        else:
            rec_data = [rec.chrom, rec.pos,
                monovar_only, sccaller_only, 0, monovar_sccaller, 0, 0, 0]
            no_snvs = sum(rec_data[2:])

            if no_snvs == 0 :
                only_wt += 1
                continue
            elif no_snvs == 1 \
                    or (no_snvs == 2 and monovar_only == 1 and sccaller_only == 1):
                data['singletons'].append(rec_data)
                continue
            elif monovar_sccaller != 0:
                if monovar_sccaller > 1 \
                        or (monovar_only > 0 and sccaller_only > 0):
                    data['monovar2+_sccaller2+'].append(rec_data)
                    rec_vcf, gt_row = get_call_output(rec, sc_calls, sample_maps[0])
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
                rec_vcf, gt_row = get_call_output(rec, sc_calls, sample_maps[0])
            else:
                import pdb; pdb.set_trace()
             
        out_vcf += rec_vcf
        gt_mat += '\n{}:{}{}{}' \
            .format(rec.chrom, rec.pos, sep, sep.join(gt_row))
        all_row = np.array([rec.ref] + [rec.alts[0]] * len(sample_maps[0]))
        all_row[np.argwhere(gt_row == '0') + 1] = rec.ref
        all_row[np.argwhere(gt_row == '3') + 1] = '?'
        all_mat.append(all_row)
        
    if len(all_mat) > 0:
        all_mat_t = np.stack(all_mat).T
    else:
        all_mat_t = False
    print('Iterating calls on Chr {} - End\n'.format(chrom))
    return data, out_vcf, gt_mat, all_mat_t, germline
    


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


def get_call_summary(samples, sample_maps, depth, quality):
    sc_calls = np.full((len(sample_maps[0]), 2), -1, dtype=int)
    snv_bulk = False
    snv_germline = False

    # 0: monovar, 1: sccaller, 2: bulk_tumor
    for sample_id, sample in samples.iteritems():
        alg, sample_map_id = get_call_ids(sample_id, sample_maps)

        # Skip all genotypes except: 0/1 | 1/1 | 0/2
        if not sample['GT'][1]:
            if sample['GT'][0] == 0 and alg < 2:
                sc_calls[sample_map_id, alg] = 0
            continue

        if alg == 2:
            if sample['DP'] < depth:
                continue
            if sample_map_id == 'normal':
                snv_germline = True
            else:
                snv_bulk = True
        else:
            # Skip indels and "False" calls (only called by sccaller)
            if alg == 1 and sample['SO'] != 'True':
                continue

            # Skip low genotype quality calls or read depth below threshold
            if sum(sample['AD']) < depth:
                continue

            if sample['GQ'] < quality:
                if 'FPL' in sample and sample['FPL'][0] != None:
                    if min(sample['FPL'][:2]) - max(sample['FPL'][2:]) \
                            < quality:
                        continue
                else:
                    PL_max = np.max([j for j in sample['PL'][1:] if j != None])
                    if sample['PL'][0] - PL_max < quality:
                        continue           
                
            sc_calls[sample_map_id, alg] = 1

    return sc_calls, snv_bulk, snv_germline


def get_call_output(rec, calls, sc_map):
    best_calls = calls.argmax(axis=1)
    data_calls = np.sum(calls.max(axis=1) > -1)

    rec_out = f'\n{rec.chrom}\t{rec.pos}\t.\t{rec.ref}\t{rec.alts[0]}\t' \
        f'{min(99, rec.qual)}\tPASS\tNS={data_calls}\tGT:AD:GQ:PL'
    
    gt_mat_row = np.zeros(len(sc_map), dtype=str)
    for sample, sample_id in sc_map.items():
        alg = best_calls[sample_id]
        gt = calls[sample_id, alg]
        if gt == -1:
            gt_mat_row[sample_id] = '3'
            rec_out += '\t./.:.:.:.'
        else:
            if alg == 0:
                call = rec.samples['{}.monovar'.format(sample)]
            else:
                call = rec.samples['{}.sccaller'.format(sample)]
            
            if 'FPL' in call and call['FPL'][0] != None:
                rec_out += f'\t{call["GT"][0]}/{call["GT"][1]}:' \
                    f'{call["AD"][0]},{call["AD"][1]}:{call["GQ"]}:' \
                    f'{call["FPL"][0]},{call["FPL"][2]},{call["FPL"][2]}'
                gt_mat_row[sample_id] = get_gt_mat_entry(call["GT"])
            else:
                PL = [i for i in call['PL'] if i != None]
                if len(PL) == 0:
                    rec_out += '\t./.:.:.:.'
                    gt_mat_row[sample_id] = '3'
                elif len(PL) == 3:
                    rec_out += f'\t{call["GT"][0]}/{call["GT"][1]}:' \
                        f'{call["AD"][0]},{call["AD"][1]}:{call["GQ"]}:' \
                        f'{PL[0]},{PL[1]},{PL[2]}'
                    gt_mat_row[sample_id] = get_gt_mat_entry(call["GT"])
                else:
                    import pdb; pdb.set_trace()

    if "None" in rec_out: import pdb; pdb.set_trace()
    return rec_out, gt_mat_row


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


def plot_venn(data, out_dir):
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

    out_file = os.path.join(args.output, 'SNP_counts.pdf')
    print('Saving call-Venn-plots to: {}'.format(out_file))
    fig.savefig(out_file, dpi=300)
    plt.close()


def merge_summaries(args):
    counts = {}
    gt_mat = ''
    nex_mat = {}
    vcf_out = ''
    vcf_map = {os.path.basename(i).split('.')[1]: i for i in args.input}
    sorted_chr = sorted(vcf_map.keys(),
        key = lambda x: int(x) if x not in ['X', 'Y'] else 23)

    for i, chr_no in enumerate(sorted_chr):
        vcf_file = vcf_map[chr_no]
        base_dir = os.path.dirname(vcf_file)

        sum_file = os.path.join(base_dir, 'Call_summary.{}.tsv'.format(chr_no))
        with open(sum_file, 'r') as f_cnt:
            lines = f_cnt.readlines()
            for line in lines:
                alg_counts, alg = line.strip().split('\t')
                try:
                    counts[alg] += int(alg_counts)
                except KeyError:
                    counts[alg] = int(alg_counts)

        gt_file = os.path.join(base_dir, 'Genotype_matrix.{}.csv'.format(chr_no))
        with open(gt_file, 'r') as f_gt:
            if i == 0:
                gt_mat += f_gt.read()
            else:
                header = f_gt.readline()
                gt_mat += '\n' + f_gt.read()     

        nex_file = os.path.join(base_dir, 'Genotype_matrix.{}.nex'.format(chr_no))
        with open(nex_file, 'r') as f_nex:
            nex_str = f_nex.read()
            start = nex_str.find('Matrix\n')
            end = nex_str.find(';', start)
            chr_mat = nex_str[start+6:end].strip()
            if chr_mat:
                for taxa in chr_mat.split('\n'):
                    taxa_info = taxa.strip().split('\t')
                    try:
                        nex_mat[taxa_info[0]] += taxa_info[1]
                    except KeyError:
                        nex_mat[taxa_info[0]] = taxa_info[1]

        vcf = VariantFile(vcf_file)
        if i == 0:
            contigs = ''
            for contig, contig_obj in vcf.header.contigs.items():
                if contig in sorted_chr:
                    contigs += str(contig_obj.header_record)
            samples = '\t'.join([i for i in vcf.header.samples])
            ref = re.search('##reference=.*\n', str(vcf.header))[0].rstrip('\n')
            vcf_out += VCF_HEADER.format(time=time.localtime(), 
                contigs=contigs.rstrip('\n'), ref=ref, samples=samples)

        for rec in vcf.fetch():
            vcf_out += str(rec)

    nex_mat_str = ''
    for sample_row in nex_mat.items():
        nex_mat_str += '\t{}\t{}\n'.format(*sample_row)

    nex_out_file = os.path.join(args.output, 'Genotype_matrix.all.nex')
    with open(nex_out_file, 'w') as f_nex:
        f_nex.write(NEXUS_TEMPLATE.format(sample_no=len(nex_mat),
            sample_labels=' '.join(nex_mat.keys()),
            rec_no=len(nex_mat[taxa_info[0]]), matrix=nex_mat_str.strip('\n')))

    vcf_out_file = os.path.join(args.output, 'all.filtered.vcf')
    with open(vcf_out_file, 'w') as f_vcf:
        f_vcf.write(vcf_out.strip('\n'))

    gt_out_file = os.path.join(args.output, 'Genotype_matrix.all.csv')
    with open(gt_out_file, 'w') as f_gt:
        f_gt.write(gt_mat.strip('\n'))

    out_QC = os.path.join(args.output, 'Call_summary.all.tsv' )
    save_summary(counts, out_QC)

    return counts


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
            args.chr = 'all_chr'

        data = get_summary_df(args)
        if args.chr == 'all_chr':
            plot_venn(data, args)

    else:
        summary = merge_summaries(args)
        plot_venn(summary, args)


if __name__ == '__main__':
    args = parse_args()
    main(args)


