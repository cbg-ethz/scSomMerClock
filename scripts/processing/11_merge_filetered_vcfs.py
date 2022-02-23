#!/usr/bin/env python3

import argparse
import os
import re
import time
from pysam import VariantFile


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

begin data;
    dimensions ntax={sample_no} nchar={rec_no};
    format datatype=dna missing=? gap=-;
    matrix
{matrix}
    ;
end;
"""


def merge_summaries(in_files, out_dir, keep_sex=False, nexus=False):
    counts = {}
    gt_mat = ''
    nex_mat = {}
    vcf_out = ''
    vcf_map = {os.path.basename(i).split('.')[1]: i for i in in_files}
    sorted_chr = sorted(vcf_map.keys(),
        key = lambda x: int(x) if x not in ['X', 'Y'] else 23)

    for i, chrom in enumerate(sorted_chr):
        if not keep_sex and chrom in ['X', 'Y']:
            continue

        vcf_file = vcf_map[chrom]
        base_dir = os.path.dirname(vcf_file)

        sum_file = os.path.join(base_dir, f'Call_summary.{chrom}.tsv')
        with open(sum_file, 'r') as f_cnt:
            lines = f_cnt.readlines()
            for line in lines:
                alg_counts, alg = line.strip().split('\t')
                try:
                    counts[alg] += int(alg_counts)
                except KeyError:
                    counts[alg] = int(alg_counts)

        gt_file = os.path.join(base_dir, f'Genotype_matrix.{chrom}.csv')
        with open(gt_file, 'r') as f_gt:
            if i == 0:
                gt_mat += f_gt.read()
            else:
                header = f_gt.readline()
                gt_mat += '\n' + f_gt.read()     

        if nexus:
            nex_file = os.path.join(base_dir, f'Genotype_matrix.{chrom}.nex')
            if os.path.exists(nex_file):
                with open(nex_file, 'r') as f_nex:
                    nex_str = f_nex.read()
                    start = nex_str.find('matrix\n')
                    end = nex_str.find('    ;', start)
                    chr_mat = nex_str[start + 6:end].strip()
                    if chr_mat:
                        for taxa in chr_mat.split('\n'):
                            taxa_info = taxa.strip().split('    ')
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

    vcf_out_file = os.path.join(out_dir, 'all_filtered.vcf')
    with open(vcf_out_file, 'w') as f_vcf:
        f_vcf.write(vcf_out.strip('\n'))

    gt_out_file = os.path.join(out_dir, 'Genotype_matrix.all.csv')
    with open(gt_out_file, 'w') as f_gt:
        f_gt.write(gt_mat.strip('\n'))

    out_QC = os.path.join(out_dir, 'Call_summary.all.tsv' )
    with open(out_QC, 'w') as f:
        for alg, calls in counts.items():
            if alg == 'singletons':
                continue
            if isinstance(calls, list):
                call_no = len(calls)
            else:
                call_no = calls
            f.write(f'{call_no}\t{alg}\n')

    if nexus:
        nex_mat_str = ''
        for sample_row in nex_mat.items():
            nex_mat_str += '{}    {}\n'.format(*sample_row)

        nex_out_file = os.path.join(out_dir, 'Genotype_matrix.all.nex')
        with open(nex_out_file, 'w') as f_nex:
            f_nex.write(NEXUS_TEMPLATE.format(sample_no=len(nex_mat),
                sample_labels=' '.join(nex_mat.keys()),
                rec_no=len(nex_mat[taxa_info[0]]), matrix=nex_mat_str.strip('\n')))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, nargs='*', help='Input files')
    parser.add_argument('-o', '--outdir', type=str, default='', 
        help='Output file. Default = <INPUT_DIR>.')
    parser.add_argument('-on', '--output_nexus', action='store_true',
        help='Write data additionally as nexus file. Default = False.')
    parser.add_argument('-s', '--keep_sex', action='store_true',
        help='If set, sex chromosomes are kept.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        outdir = os.path.dirname(snakemake.output[0])
        merge_summaries(snakemake.input, outdir)
    else:
        args = parse_args()
        if not args.outdir:
            args.outdir = os.path.dirname(args.input[0])
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)

        merge_summaries(args.input, args.outdir, keep_sex=args.keep_sex,
            nexus=args.output_nexus)