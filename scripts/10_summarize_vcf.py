#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
import vcf

##                     < FOORMAT
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


def main(args):
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

        # Only called in SC with low quality
        try:
            if record.QUAL < args.quality and record.FILTER == None:
                continue
        except TypeError:
            # Only called by Mutect2 but didnt pass filtering
            if 'PASS' not in record.FILTER:
                continue

        # if qual_record >= args.quality:
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
                if 'PASS' not in record.FILTER:
                    continue
            else:
                if sample.data.GQ < args.genotype_quality:
                    continue
                # Skip samples with read depth below threshold
                if sum(sample.data.AD) < args.read_depth:
                    continue

            # Bulk SNV
            if alg == 'mutect':
                # Called in Normal (germline)
                if sample_name == args.bulk_normal:
                    germline.append(f'{record.CHROM}:{record.POS}')
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
    df.to_csv(os.path.join(args.output, f'filtered_summary.{file_name}.tsv'),
        sep='\t')

    if germline:
        germ_file = os.path.join(args.output, f'germline_muts.{file_name}.tsv')
        with open(germ_file, 'w') as f:
            f.write('\n'.join(germline))


if __name__ == '__main__':
    args = parse_args()
    main(args)
