#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
from pysam import VariantFile, FastaFile, AlignmentFile

CHROM = [str(i) for i in range(1, 23, 1)] + ['X', 'Y']


def parse_args():
    parser = argparse.ArgumentParser(
        description='*** Calculate ADO rate per single cell. ***')
    parser.add_argument('input', type=str,  nargs='*',
        help='Path(s) to input bam files.')
    parser.add_argument('-b', '--bulk', type=str, required=True,
        default='Path to vcf file containing mutations called in tumor bulk.')
    parser.add_argument('-d', '--dbsnp', type=str, required=True,
        help='Path to zipped and indexed DBSNP file.')
    parser.add_argument('-r', '--ref', type=str, default='',
        help='Path to fasta reference file.')
    parser.add_argument('-o', '--outfile', type=str, default='',
        help='Output file name. Default = <BULK DIR>/ADO_rates.tsv.')
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


def get_high_conf_het_in_samples(args):
    hc_hets = get_high_confidence_heterozygous(args)

    pileup_arg = {'stepper': 'samtools', 'adjust_capq_threshold': 50,
        'min_mapping_quality': 20, 'min_base_quality': 13, 'truncate': True}

    if args.ref:
        fasta_file = FastaFile(args.ref)
        pileup_arg['fastafile'] = fasta_file

    df = pd.DataFrame([], index=hc_hets.keys(), dtype=int)
    df.index.set_names(['chrom', 'pos'], inplace=True)
    # Iterate over input files
    cell_names = set([])
    for bam_file_path in sorted(args.input):
        sample_name = os.path.basename(bam_file_path).split('.')
        cell = sample_name[0]
        if sample_name[2] in CHROM:
            s_chrom = [sample_name[2]]
        else:
            s_chrom = df.index.get_level_values('chrom').unique()
        cell_names.add(cell)

        bam_file = AlignmentFile(bam_file_path, 'rb')
        # Iterate over high-confidence heterozygous ones
        
        for chrom, pos in df.loc[s_chrom].index:
            alt = hc_hets[(chrom, pos)]
            pileup_arg['contig'] = chrom
            pileup_arg['start'] = pos
            pileup_arg['stop'] = pos + 1
            pileup_data = bam_file.pileup(**pileup_arg)
            for i, pileup in enumerate(pileup_data):
                seq = ''.join(pileup.get_query_sequences()).upper()
                alt_no = seq.count(alt)
                total_no = pileup.get_num_aligned()
                if i > 0: import pdb; pdb.set_trace()
            try:
                df.loc[(chrom, pos), cell + '_alt'] = alt_no
                df.loc[(chrom, pos), cell + '_n'] = total_no
            except NameError:
                pass
  
    out_dir = os.path.dirname(args.outfile)
    if len(cell_names) == 1:
        out_cell_name = cell_names.pop()
        out_file = os.path.join(out_dir,
            'ADO_overview_{}.tsv'.format(out_cell_name))
    else:
        out_file = os.path.join(out_dir, 'ADO_overview.tsv')
    df.to_csv(out_file, sep='\t')
    return df


def calc_ADO_rates(df, args):
    total_cols = [i for i in df.columns if i.endswith('_n')]
    n_cols = [i for i in df.columns if i.endswith('_alt')]
    import pdb; pdb.set_trace()


def main(args):
    if args.outfile == '':
        args.outfile = os.path.join(os.path.dirname(args.bulk), 'ADO_rates.tsv')

    df = get_high_conf_het_in_samples(args)
    # import pdb; pdb.set_trace()
    # int_file = '/home/hw/Desktop/molClock_project/scDNA_data/Ni9/QC/ADO_overview.tsv'
    # df = pd.read_csv(int_file, sep='\t', index_col=[0, 1])
    calc_ADO_rates(df, args)


if __name__ == '__main__':
    args = parse_args()
    main(args)