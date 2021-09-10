#!/usr/bin/env python3

import argparse
import gzip
import numpy as np
import re


def postprocess_vcf(vcf_file, out_file, minDP=1, minGQ=0, s_minDP=5,
            s_minGQ=1, s_filter=False, stats_file=''):
    if vcf_file.endswith('gz'):
        file_stream = gzip.open(vcf_file, 'rb')
    else:
        file_stream = open(vcf_file, 'r')

    header = ''
    body = ''

    rows_skipped = 0
    mut_no = np.zeros(2)
    stats = [0, 0, 0, 0, 0] # TP, FP, TN, FN, MS
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
                    if monovar:
                        format_short = 'GT:AD:DP:GQ:PL'
                        missing = '.|.:.,.:.:.:.'
                    else:
                        header += '##FORMAT=<ID=GQ,Number=1,Type=Integer,' \
                            'Description="Genotype Quality">\n'
                        format_short = 'GT:DP:RC:GQ:TG'
                        missing = '.|.:.:.,.,.,.:.:.|.'

                header += line
                continue
            elif line.strip() == '':
                break

            # VCF records
            line_cols = line.strip().split('\t')

            bases = {line_cols[3]: 0}
            ref = line_cols[3]
            alts = line_cols[4].split(',')
            for i, j in enumerate(alts, 1):
                bases[j] = i

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
                    gt = '|'.join(sorted([i for i in s_rec[0].split('|')]))
                    dp = s_rec[1]
                    rc = s_rec[2]
                    # g10 = s_rec[3]
                    try:
                        tgt = '|'.join(sorted([str(bases[i]) \
                            for i in s_rec[-1].split('|')]))
                    except KeyError:
                        tgt = s_rec[-1].replace(ref, '0')
                        if tgt[0] != '0':
                            bases[tgt[0]] = len(bases)
                            tgt = tgt.replace(tgt[0], str(len(bases)))
                        if tgt[-1] != '0':
                            tgt = tgt.replace(tgt[-1], str(len(bases)))

                    assert len(tgt) == 3, f'Wrong true GT field: {tgt}'
                    if len(FORMAT_col) == len(s_rec):
                        pln = np.array([-float(i) \
                            for i in s_rec[PLN_col].split(',')])
                        gq = min(99, sorted(pln - pln.min())[1])
                    else:
                        gq = -1
                    true_gt[s_i] = tgt

                    if gt == '.|.':
                        stats[4] += 1
                    elif (int(tgt[-1]) > 0) and (int(gt[-1]) > 0): #TP
                        stats[0] += 1
                    elif (int(tgt[-1]) == 0) and (int(gt[-1]) > 0): #FP
                        stats[1] += 1
                    elif  (int(tgt[-1]) == 0) and (int(gt[-1]) == 0): # TN
                        stats[2] += 1
                    elif (int(tgt[-1]) > 0) and (int(gt[-1]) == 0): # FN
                        stats[3] += 1
                    else:
                        raise RuntimeError('Unknown case for stats update')

                if int(dp) < minDP or gq < minGQ:
                    new_line[s_i] = missing + tgt
                    genotypes[s_i] = '.|.'
                else:
                    if monovar:
                        new_line[s_i] = s_rec_raw
                    else:
                        new_line[s_i] = f'{gt}:{dp}:{rc}:{gq:.0f}:{tgt}'
                    genotypes[s_i] = gt

            called_gt = genotypes[genotypes != b'.|.']
            # Count true mutations
            if not monovar:
                if np.unique(true_gt).size > 1:
                    # import pdb; pdb.set_trace()
                    mut_no[0] += 1
                else:
                    all_gt = true_gt[0].decode()
                    if all_gt[0] == '0' and all_gt[-1] == '0':
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
            else:
                # One or more different genotypes detected
                diff_gt, diff_count = np.unique(called_gt, return_counts=True)
                # Check if any non-wildtype is called more than once
                is_singleton = True
                for i in np.argwhere(diff_gt != b'0|0'):
                    if diff_count[i] > 1:
                        is_singleton = False
                        break

                if is_singleton:
                    cell_no = np.argwhere(
                        (genotypes != b'.|.') & (genotypes != b'0|0')).flatten()[0]
                    cell_rec = line_cols[9 + cell_no].split(':')
                    s_pln = np.array([-float(i) \
                        for i in cell_rec[PLN_col].split(',')])
                    s_gq = min(99, sorted(pln - pln.min())[1])

                    if int(cell_rec[1]) < s_minDP or s_gq < s_minGQ:
                        new_line[cell_no] = missing + tgt
                        filter_str = 'wildtype'
                    elif s_filter:
                        filter_str = 'singleton'

            body += '\t'.join(line_cols[:6]) + '\t{}\t{}\t{}\t' \
                    .format(filter_str, line_cols[7], format_short) \
                + '\t'.join(new_line) + '\n'

    if rows_skipped:
        header = re.sub('(?<=droppedRows\=)0(?=\n)', str(rows_skipped), header)

    with gzip.open(out_file, 'wb') as f_out:
        f_out.write(f'{header}{body}'.encode())

    if stats_file:
        run_no = vcf_file.split('.')[-1]

        with open(stats_file, 'a+') as f:
            f.write(f'{run_no}\t{mut_no.sum()}\t{np.sum(stats)}\t{stats[0]}\t' \
                f'{stats[1]}\t{stats[2]}\t{stats[3]}\t{stats[4]}\n')


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
    parser.add_argument('-sgq', '--singeltonMinGQ', type=int, default=0,
        help='Min. Genotype Quality to include a locus (else missing). '
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
        if not snakemake.params['s_gq']:
            s_gq = snakemake.params['gq']
        else:
            s_gq = snakemake.params['s_gq']

        postprocess_vcf(snakemake.input[0], snakemake.output[0],
            minDP=snakemake.params['dp'],
            minGQ=snakemake.params['gq'],
            s_minDP=s_dp,
            s_minGQ=s_gq,
            s_filter=snakemake.params['singletons'],
            stats_file=snakemake.input[1])
    else:
        args = parse_args()
        if not args.output:
            args.output = args.input + '.final'
        postprocess_vcf(args.input, args.output,
            minDP=args.minDP,
            minGQ=args.minGQ,
            s_minDP=args.singeltonMinDP,
            s_minGQ=args.singeltonMinGQ,
            s_filter=args.filter_singletons)