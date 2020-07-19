#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
import vcf

## SCcaller FOORMAT
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
        '-bt', '--bulk_tumor', nargs='+', type=str, default='',
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
        '-gq', '--genotype_quality', type=int, default=20,
        help='Minimum read depth at loci. Default = 40.'
    )

    args = parser.parse_args()
    return args


def main(args):
    vcf_reader = vcf.Reader(filename=args.input)

    samples = {'sc': set([])}
    for sample in vcf_reader.samples:
        sample_detail = sample.split('.')
        sample_name = '.'.join(sample_detail[:-1])
        caller = sample_detail[-1]
        if caller != 'mutect':
            samples['sc'].add(sample_name)

    alg_map = {'monovar': 0, 'sccaller': 1, args.bulk_normal: 2}
    for bt in args.bulk_normal:
        alg_map[bt] = len(alg_map)

    
    sc_map = {j: i for i, j in enumerate(sorted(samples['sc']))}
    sc_size = len(sc_map)

    data = []
    # Iterate over rows
    for record in vcf_reader:
        if record.QUAL == None: import pdb; pdb.set_trace()
        # Skip rows with quality below a threshold
        if record.QUAL >= args.quality:
            rec_data = [record.CHROM, record.POS, record.var_type]
            calls = np.zeros((sc_size, 4))
            # Iterate over columns (i.e. samples)
            for sample in record.samples:
                # Skip samples where record is not called
                QG = sample.data.GQ
                if QG == None and sample.called: import pdb; pdb.set_trace()

                if sample.called and QG > args.genotype_quality:
                    # Skip samples with read depth below threshold
                    try: 
                        depth = sample.data.DP
                        if depth == None:
                            raise AttributeError
                    except AttributeError:
                        depth = sum(sample.data.AD)
                    if depth < args.read_depth:
                            continue

                    sample_detail = sample.sample.split('.')
                    sample_id = sc_map['.'.join(sample_detail[:-1])]
                    alg = sample_detail[-1]
                    if alg != 'mutect':
                        calls[sample_id, alg_map[alg]] = 1
                    else:
                        import pdb; pdb.set_trace()
                        calls[sample_id, alg_map[sample_id]] = 1

            per_sample = np.sum(calls, axis=1)
            # Nothing called (why reported?)
            if per_sample.max() == 0:
                continue
            # Algorithm only called in 1 
            elif per_sample.max() == 1:
                call_data = calls.sum(axis=0)
            else:
                import pdb; pdb.set_trace()

            rec_data.extend(call_data)
            data.append(rec_data)

    cols = ['CHROM', 'POS', 'type', 'monovar', 'sccaller', 'mutect']
    summary = pd.DataFrame(columns=cols)


if __name__ == '__main__':
    args = parse_args()
    main(args)