#!/usr/bin/env python3

import argparse
import gzip
import numpy as np

import re


READS = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
READS_INV = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
GT_MAPPING = {0: 'AA', 1: 'AC', 2: 'AG', 3: 'AT', 4: 'CC', 5: 'CG', 6: 'CT',
    7: 'GG', 8: 'GT', 9: 'TT'}


def postprocess_vcf(vcf_file, out_file, minDP=1, minGQ=0, s_minDP=5,
            s_minAlt=1, s_filter=False, stats_file='', ADO_file=''):
    if vcf_file.endswith('gz'):
        file_stream = gzip.open(vcf_file, 'rb')
    else:
        file_stream = open(vcf_file, 'r')

    header = ''
    body = ''

    rows_skipped = 0
    mut_no = np.zeros(4) # muts cells, muts outg, FP cells/outg, ISA violations
    stats = np.zeros(5) # TP, FP, TN, FN, MS
    monovar = False
    with file_stream as f_in:
        for line in f_in:
            try:
                line = line.decode()
            except AttributeError:
                pass
            # Skip VCF header lines
            if line.startswith('#'):
                # Safe column headers
                if line.startswith('##source') and 'MonoVar' in line:
                    monovar = True
                elif line.startswith('#CHROM'):
                    header += '##FILTER=<ID=singleton,Description="SNP ' \
                        'is singleton">\n##FILTER=<ID=wildtype,Description="' \
                        'Only 0|0 called">\n##droppedRows=0\n'
                    sample_no = len(line.strip().split('\t')[9:])
                    if re.search('_bulk(\d+)x_', vcf_file):
                        sample_no -= 1
                        # Remove doubled outgrp from sample list
                        line = '\t'.join(
                            line.split('\t')[:-2] + [line.split('\t')[-1]])
                    ADOs = np.zeros((2, sample_no), dtype=int)
                    header += '##FORMAT=<ID=GQ,Number=1,Type=Integer,' \
                        'Description="Genotype Quality">\n'
                    format_short = 'GT:DP:RC:GQ:TG'
                    header += '\t'.join([i.strip() for i in line.split('\t')]) \
                        + '\n'

                header += line
                continue
            elif line.strip() == '':
                break

            # VCF records
            line_cols = line.strip().split('\t')

            ref = line_cols[3]
            alts = line_cols[4].split(',')

            bases = {ref: 0}
            for i, j in enumerate(alts, 1):
                bases[j] = i

            FORMAT_col = line_cols[8].split(':')
            if monovar:
                if line_cols[6] == '.':
                    line_cols[6] = 'PASS'
            else:
                PLN_col = FORMAT_col.index('PLN')
                true_gt = np.zeros((sample_no, 2), dtype=int)

            new_line = np.zeros(sample_no, dtype=object)
            genotypes = np.full((sample_no, 2), np.nan, dtype=float)
            reads = np.zeros((sample_no, 4), dtype=int)

            for s_i, s_rec_raw in enumerate(line_cols[9:9+sample_no]):
                s_rec = s_rec_raw.split(':')
                gt = sorted([int(i) for i in s_rec[0].replace('.', '-1').split('|')])
                dp = int(s_rec[1])
                rc = s_rec[2]

                try:
                    tgt = sorted([bases[i] for i in s_rec[-1].split('|')])
                # True alt in vcf
                except KeyError:
                    new_base = [i for i in s_rec[-1].split('|') \
                        if i not in bases][0]
                    bases[new_base] = len(bases)
                    tgt = sorted([bases[i] for i in s_rec[-1].split('|')])
                assert len(tgt) == 2, f'Wrong true GT field: {tgt}'
                true_gt[s_i] = tgt

                if len(FORMAT_col) == len(s_rec):
                    pln = np.array([-float(i) \
                        for i in s_rec[PLN_col].split(',')])
                    gq = min(99, sorted(pln - pln.min())[1])
                else:
                    gq = -1

                if gt == [-1, -1] and dp > 0:
                    gt_raw = GT_MAPPING[
                        np.array(s_rec[3].split(','), dtype=float).argmax()]
                    gt = sorted([bases[gt_raw[0]], bases[gt_raw[1]]])

                if (tgt[1] > 0) and (gt[1] > 0): #TP
                    stats[0] += 1
                elif (tgt[1] == 0) and (gt[1] > 0): #FP
                    stats[1] += 1
                elif  (tgt[1] == 0) and (gt[1] == 0): # TN
                    stats[2] += 1
                elif (tgt[1] > 0) and (gt[1] == 0): # FN
                    stats[3] += 1
                elif gt == [-1, -1]:
                    stats[4] += 1
                else:
                    raise RuntimeError('Unknown case for stats update')

                if (dp < minDP) or (dp == 0) or (gq < minGQ):
                    new_line[s_i] = f'.|.:{dp}:{rc}:{gq:.0f}:{tgt[0]}|{tgt[1]}'
                else:
                    if monovar:
                        new_line[s_i] = s_rec_raw
                    else:
                        if gt == [-1, -1]:
                            gt_out = '.|.'
                        else:
                            gt_out = f'{gt[0]}|{gt[1]}'
                        new_line[s_i] = f'{gt_out}:{dp}:{rc}:{gq:.0f}:' \
                            f'{tgt[0]}|{tgt[1]}'
                    genotypes[s_i] = gt

                reads[s_i] = rc.split(',')

            bases_rev = {j: i for i,j in bases.items()}
            alt_bases = np.unique(true_gt[:,1])
            if alt_bases.size == 1:
                pass
            elif alt_bases.size == 2:
                alt_base = alt_bases[1]
                het = np.all(true_gt == np.array([0, alt_base]), axis=1)
                # ADO on mut allele: wt|mut -> wt|wt
                ADOs[0] += het & (reads[:,READS[bases_rev[alt_base]]] == 0)
                # ADO on wt allele: wt|mut -> mut|mut
                ADOs[0] += het & (reads[:,READS[ref]] == 0)
                # Total number of true het muts
                ADOs[1] += het
            else:
                print(f'Unknown ADO calc for true genotypes: {reads}')

            called_gt = genotypes[~np.isnan(genotypes.sum(axis=1))]
            # Count true mutations
            unique_genos = np.unique(true_gt[:-1], axis=0) # exclude outgroup
            # All wt or all mut
            if unique_genos.shape[0] == 1:
                all_wt = (unique_genos[0] == 0).all()
                outg_wt = all(true_gt[-1] == 0)
                # Some wildtype and some mutations
                if not all_wt:
                    mut_no[0] += 1
                # False positive: all true ones are wildtype
                elif all_wt and outg_wt:
                    mut_no[2] += 1
                # Outgroup mutation: wildtype in all cells except outgroup
                elif all_wt and not outg_wt:
                    mut_no[1] += 1
                else:
                    import pdb; pdb.set_trace()

            # Some wt and some muts
            elif unique_genos.shape[0] == 2:
                # ISA violation: different mutations at locus
                if not any((unique_genos == 0).sum(axis=0) == 2):
                    mut_no[3] += 1
                mut_no[0] += 1
            # ISA violation: different mutations at locus
            else:
                import pdb; pdb.set_trace()

            # No signal above filtering threshold for position: skip
            if called_gt.size == 0:
                rows_skipped += 1
                continue

            filter_str = line_cols[6]

            # Only wildtype called
            if called_gt.sum() == 0:
                filter_str = 'wildtype'
            # >1 diff. genotypes detected, and check for singleton
            elif s_minDP > 0 or s_minAlt > 0:
                diff_gt, diff_count = np.unique(called_gt,
                    return_counts=True, axis=0)

                # Check if any non-wildtype is called more than once
                is_singleton = True
                for i, gt_i in enumerate(diff_gt):
                    if gt_i.sum() != 0 and diff_count[i] > 1:
                        is_singleton = False
                        break

                if is_singleton:
                    cell_no = np.argwhere((genotypes[:,0] == gt_i[0]) \
                            & (genotypes[:,1] == gt_i[1])).flatten()[0]
                    cell_rec = line_cols[9 + cell_no].split(':')
                    cell_reads = np.array(cell_rec[2].split(','), dtype=int)
                    alt_reads = cell_reads[READS[alts[int(gt_i[1] - 1)]]]

                    if int(cell_rec[1]) < s_minDP or alt_reads < s_minAlt:
                        new_line[cell_no] = f'.|.:{dp}:{rc}:{gq:.0f}:' \
                            f'{tgt[0]}|{tgt[1]}'
                        filter_str = 'wildtype'
                    elif s_filter:
                        filter_str = 'singleton'

                # If two different singleton het muts are called, filter these
                two_diff = False
                diff_bp, bp_count = np.unique(called_gt[:,1], return_counts=True)
                if (bp_count == 1).sum() == 2:
                    rel_cells = np.argwhere(
                        np.ma.masked_invalid(genotypes).sum(axis=1) > 0)
                    for cell_no in rel_cells:
                        cell_rec = line_cols[9 + cell_no[0]].split(':')
                        if cell_rec[0].startswith('.'):
                            cell_rec[1] = -1
                            alt_reads = -1
                        else:
                            cell_reads = np.array(cell_rec[2].split(','), dtype=int)
                            alt_id = READS[alts[int(cell_rec[0][-1]) - 1]]
                            alt_reads = cell_reads[alt_id]

                        if int(cell_rec[1]) < s_minDP or alt_reads < s_minAlt:
                            new_line[cell_no[0]] = '.|.' \
                                + new_line[cell_no[0]][3:]

            body += '\t'.join(line_cols[:6]) + '\t{}\t{}\t{}\t' \
                    .format(filter_str, line_cols[7], format_short) \
                + '\t'.join(new_line) + '\n'

    if rows_skipped:
        header = re.sub('(?<=droppedRows\=)0(?=\n)', str(rows_skipped), header)

    with gzip.open(out_file, 'wb') as f_out:
        f_out.write(f'{header}{body}'.encode())

    run_no = vcf_file.split('.')[-1]

    if ADO_file:
        with open(ADO_file, 'a+') as f:
            ADO_rate = (ADOs[0] / ADOs[1]).round(4)
            f.write(f'{run_no}\t' + '\t'.join(ADO_rate.astype(str)) +'\n')

    if stats_file:
        out_line = f'{run_no}\t{mut_no[0]}\t{mut_no[1]}\t{mut_no[2]}\t' \
            f'{np.sum(stats)}\t{stats[0]}\t{stats[1]}\t{stats[2]}\t{stats[3]}\t' \
            f'{stats[4]}\n'
        with open(stats_file, 'a+') as f:
            f.write(out_line)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str,
        help='Absolute or relative path(s) to input file(s)')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Path to the output file. Default = <INPUT_DIR>.')
    parser.add_argument('-dp', '--minDP', type=int, default=1,
        help='Min. reads to include a locus (else missing). Default = 1.')
    parser.add_argument('-gq', '--minGQ', type=int, default=0,
        help='Min. Genotype Quality to include a locus (else missing). '
            'Default = 0.')
    parser.add_argument('-sdp', '--singeltonMinDP', type=int, default=1,
        help='Min. reads to include a locus (else missing). Default = 1.')
    parser.add_argument('-sgq', '--singeltonMinAlt', type=int, default=1,
        help='Min. alternative reads to include a locus (else missing). '
            'Default = 0.')
    parser.add_argument('-fs', '--filter_singletons', action='store_true',
        help='If set, singleton SNPs are filtered out.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():

        if not snakemake.params['s_dp']:
            s_dp = snakemake.params['dp']
        else:
            s_dp = snakemake.params['s_dp']
        if not snakemake.params.get('s_alt', 0):
            s_alt = 0
        else:
            s_alt = snakemake.params['s_alt']

        postprocess_vcf(snakemake.input[0], snakemake.output[0],
            minDP=snakemake.params['dp'],
            minGQ=snakemake.params['gq'],
            s_minDP=s_dp,
            s_minAlt=s_alt,
            s_filter=snakemake.params['singletons'],
            stats_file=snakemake.input[1],
            ADO_file=snakemake.input[2])
    else:
        args = parse_args()
        if not args.output:
            args.output = args.input + '.final'
        postprocess_vcf(args.input, args.output,
            minDP=args.minDP,
            minGQ=args.minGQ,
            s_minDP=args.singeltonMinDP,
            s_minAlt=args.singeltonMinAlt,
            s_filter=args.filter_singletons)