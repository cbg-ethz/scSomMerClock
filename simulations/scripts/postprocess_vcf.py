#!/usr/bin/env python3

import re
import argparse
import numpy as np


def postprocess_vcf(vcf_file, out_file, minDP=1, minGQ=0, s_filter=False):
    if vcf_file.endswith('gz'):
        file_stream = gzip.open(vcf_file, 'rb')
    else:
        file_stream = open(vcf_file, 'r')

    header = ''
    body = ''

    rows_skipped = 0
    mut_no = np.zeros(2)
    monovar = False
    with file_stream as f_in:
        for line in f_in:
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
                    if monovar:
                        format_short = 'GT:AD:DP:GQ:PL'
                        missing = '.|.:.,.:.:.:.'
                    else:
                        header += '##FORMAT=<ID=GQ,Number=1,Type=Integer,' \
                            'Description="Genotype Quality">\n'
                        format_short = 'GT:DP:RC:G10:PL:GQ:TG'
                        missing = '.|.:.:.,.,.,.:.:.:.:'

                header += line
                continue
            elif line.strip() == '':
                break
            
            # VCF records
            line_cols = line.strip().split('\t')

            ref = line_cols[3]
            alts = line_cols[4].split(',')
            FORMAT_col = line_cols[8].split(':')
            if monovar:
                if line_cols[6] == '.':
                    line_cols[6] = 'PASS'
            else:
                PLN_col = FORMAT_col.index('PLN')
                true_gt = np.chararray(sample_no, itemsize=3)

            new_line = np.zeros(sample_no, dtype=object)
            genotypes = np.chararray(sample_no, itemsize=3)

            for s_i, s_rec_raw in enumerate(line_cols[9:]):
                s_rec = s_rec_raw.split(':')

                if monovar:
                    gt = s_rec[0]
                    gq = int(s_rec[3])
                    dp = s_rec[2]
                    tg = ''
                else:
                    gt = s_rec[0]
                    dp = s_rec[1]
                    rc = s_rec[2]
                    g10 = s_rec[3]
                    tg = s_rec[-1]
                    if len(FORMAT_col) == len(s_rec):
                        pln = np.array([-float(i) \
                            for i in s_rec[PLN_col].split(',')])
                        gq = min(99, sorted(pln - pln.min())[1])
                        pl = ','.join([str(max(int(i), -999)) \
                            for i in s_rec[PLN_col - 1].split(',')])
                    else:
                        gq = -1
                        pl = '.'
                    true_gt[s_i] = tg
                
                if int(dp) < minDP or gq < minGQ:
                    new_line[s_i] = missing + tg
                    genotypes[s_i] = '.|.'
                else:
                    if monovar:
                        new_line[s_i] = s_rec_raw
                    else:
                        new_line[s_i] = '{}:{}:{}:{}:{}:{:.0f}:{}' \
                            .format(gt, dp, rc, g10, pl, gq, tg)
                    genotypes[s_i] = gt

            called_gt = genotypes[genotypes != b'.|.']
            if not monovar:
                if np.unique(true_gt).size > 1:
                    # import pdb; pdb.set_trace()
                    mut_no[0] += 1
                else:
                    all_gt = true_gt[0].decode()
                    if all_gt[0] == ref and all_gt[-1] == ref:
                        mut_no[1] += 1
                    else:
                        import pdb; pdb.set_trace()

            # No signal above filtering threshold for position: skip
            if called_gt.size == 0:
                rows_skipped += 1
                continue

            filter_str = line_cols[6]
            # Only wildtype called
            if np.all(called_gt == b'0|0'):
                filter_str = 'wildtype'
            elif s_filter:
                # One or more different genotypes detected
                diff_gt, diff_count = np.unique(called_gt, return_counts=True)
                # Check if any non-wildtype is called more than once
                is_singleton = True
                for i in np.argwhere(diff_gt != b'0|0'):
                    if diff_count[i] > 1:
                        is_singleton = False
                        break

                if is_singleton:
                    filter_str = 'singleton'

            body += '\t'.join(line_cols[:6]) + '\t{}\t{}\t{}\t' \
                    .format(filter_str, line_cols[7], format_short) \
                + '\t'.join(new_line) + '\n'

    if rows_skipped:
        header = re.sub('(?<=droppedRows\=)0(?=\n)', str(rows_skipped), header)

    with open(out_file, 'w') as f_out:
        f_out.write('{}{}'.format(header, body))


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
            'Default = 1.')
    parser.add_argument('-fs', '--filter_singletons', action='store_true',
        help='If set, singleton SNPs are filtered out.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        postprocess_vcf(snakemake.input[0], snakemake.output[0],
            minDP=snakemake.params['dp'], minGQ=snakemake.params['gq'],
            s_filter=snakemake.params['singletons'])
    else:
        args = parse_args()
        if not args.output:
            args.output = args.input + '.final'
        postprocess_vcf(args.input, args.output, minDP=args.minDP,
            minGQ=args.minGQ, s_filter=args.filter_singletons)