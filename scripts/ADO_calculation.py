#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
from pysam import VariantFile, FastaFile, AlignmentFile


def parse_args():
    parser = argparse.ArgumentParser(
        description='*** Calculate ADO rate per single cell. ***')
    parser.add_argument('input', type=str,  nargs='*',
        help='Path(s) to input bam files.')
    parser.add_argument('-b', '--bulk', type=str, required=True,
        default='Path to vcf file containing mutations called in tumor bulk.')
    parser.add_argument('-d', '--dbsnp', type=str, required=True,
        help='Path to zipped and indexed DBSNP file.')
    parser.add_argument('-o', '--outfile', type=str, default='',
        help='Ooutput file name. Default = <BULK DIR>/ADO_rates.tsv.')
    args = parser.parse_args()
    return args


def get_high_confidence_heterozygous(args):
    bulk = VariantFile(args.bulk)
    dbsnp = VariantFile(args.dbsnp)

    hc_het = {}
    for bulk_rec in bulk.fetch():
        # Skip deletion and insertions
        if len(bulk_rec.ref) > 1 or len(bulk_rec.alts) > 1 \
                or len(bulk_rec.alts[0]) > 1:
            continue

        is_hc = False
        for hc_rec in dbsnp.fetch(bulk_rec.chrom, bulk_rec.pos-1, bulk_rec.pos):
            if bulk_rec.alts[0] in hc_rec.alts:
                if hc_rec.info['VC'] == 'SNV':
                    is_hc = True
                elif hc_rec.info['VC'] == 'DIV':
                    continue
                else:
                    import pdb; pdb.set_trace()
            else:
                continue

        if is_hc:
            hc_het[(bulk_rec.chrom, bulk_rec.pos)] = bulk_rec.alts[0]
    
    return hc_het


def main(args):
    hc_hets = get_high_confidence_heterozygous(args)

    df = pd.DataFrame([], index=hc_hets.keys(), dtype=int)
    df.index.set_names(['chrom', 'pos'], inplace=True)
    # Iterate over input files
    for bam_file_path in sorted(args.input):
        cell = os.path.basename(bam_file_path).split('.')[0]

        bam_file = AlignmentFile(bam_file_path, 'rb')
        # Iterate over high-confidence heterozygous ones
        for (chrom, pos), alt in hc_hets.items():
            pileup_data = bam_file.pileup(chrom, pos, pos + 1, truncate=True)
            for i, pileup in enumerate(pileup_data):
                alt_no = np.sum([i == 'A' for i in pileup.get_query_sequences()])
                total_no = pileup.nsegments
                if i > 0: import pdb; pdb.set_trace()
            df.loc[(chrom, pos), f'{cell}_alt'] = alt_no
            df.loc[(chrom, pos), f'{cell}_n'] = total_no
    
    if args.outfile == '':
        args.outfile = os.path.join(os.path.dirname(args.bulk), 'ADO_rates.tsv')
    df.astype(int).to_csv(args.outfile, sep='\t')


if __name__ == '__main__':
    args = parse_args()
    main(args)