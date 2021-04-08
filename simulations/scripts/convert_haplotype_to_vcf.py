#!/usr/bin/env python3

import os
import argparse


VCF_TEMPLATE = """##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="True genotype">
##contig=<ID=1,length={contig_length}>
##source=germline variants from CellCoal simulation
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT{cells}
{rows}
"""


def save_fasta(ref, out_dir, run_no):
    lines = ((len(ref) - 1) // 60) + 1
    lines_str = [ref[i * 60: (i + 1) * 60] for i in range(lines)]
    fa_str = '\n'.join(lines_str)

    out_file = os.path.join(out_dir, 'reference.{}.fa'.format(run_no))
    with open(out_file, 'w') as f:
        f.write('>1 - replicate.{}\n{}'.format(run_no, fa_str))


def save_germline_vcf(data, ref, out_dir, run_no):
    gl = []
    for al, al_bases in data.items():
        gl.extend([(i, j) for i,j in enumerate(al_bases) if j != ref[i]])

    out_str = ''
    for i, (pos, alt) in enumerate(gl):
        out_str += '1\t{pos}\tgermline{id}\t{ref}\t{alt}\t.\tPASS\t.\t.\n' \
            .format(pos=pos, id=i, ref=ref[pos], alt=alt)

    out_file = os.path.join(out_dir, 'germlines.{}.vcf'.format(run_no))
    with open(out_file, 'w') as f:
        f.write(VCF_TEMPLATE \
            .format(contig_length=len(ref), rows=out_str.rstrip(), cells=''))


def save_true_vcf(data, ref, out_file):
    out_str = ''
    cells = data.keys()
    for i, ref_al in enumerate(ref):
        is_SNV = False
        alts = []
        recs = []
        for cell in cells:
            al1 = data[cell]['m'][i]
            al2 = data[cell]['p'][i]
            # Homozygous
            if al1 != ref_al and al2 != ref_al:
                if not al1 in alts:
                    alts.append(al1)
                recs.append('{i}|{i}'.format(i=alts.index(al1)))
                is_SNV = True
            # Heterozygous maternal allele
            elif al1 != ref_al:       
                if not al1 in alts:
                    alts.append(al1)
                recs.append('0|{}'.format(alts.index(al1)+1))
                is_SNV = True
            # Heterozygous paternal allele
            elif al2 != ref_al:
                if not al2 in alts:
                    alts.append(al2)
                recs.append('0|{}'.format(alts.index(al2)+1))
                is_SNV = True
            # Wildtype
            else:
                recs.append('0|0')

        if is_SNV:
            alt_str = ','.join(alts)
            rec_str = '\t'.join(recs)
            out_str += '1\t{pos}\tsnv{id:06d}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{recs}\n' \
                .format(pos=i+1, id=i+1, ref=ref_al, alt=alt_str, recs=rec_str)

    with open(out_file, 'w') as f:
        f.write(VCF_TEMPLATE\
            .format(contig_length=len(ref), rows=out_str.rstrip(),
                cells='\t{}'.format('\t'.join(cells))))


def haplotypes_to_vcf(true_hap_file, out_file):
    run_no = os.path.basename(true_hap_file).split('.')[-1]

    ref = []
    with open(true_hap_file, 'r') as f:
        haps = {}
        for cell_hap in f.read().strip().split('\n')[1:]:
            cell_name, cell_bases = [i for i in cell_hap.split(' ') if i]
            base_array = cell_bases
            if cell_name.startswith('outgcell'):
                ref = base_array
                continue

            try:
                haps[cell_name[:-1]][cell_name[-1]] = base_array
            except KeyError:
                haps[cell_name[:-1]] = {cell_name[-1]: base_array}

    if len(ref) == 0:
        raise IOError('true_hap need to contain outgcell! (OUTPUT arg: -5)')

    out_dir = os.path.dirname(out_file)
    if not os.path.exists(out_dir):
        os.path.mkdir(out_dir)

    save_fasta(ref, out_dir, run_no)
    outgr = haps.pop('outgroot')
    save_germline_vcf(outgr, ref, out_dir, run_no)
    save_true_vcf(haps, ref, out_file)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='Input file')
    parser.add_argument('-o', '--output', type=str, help='Output file.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        haplotypes_to_vcf(snakemake.input[0], snakemake.output[0])
    else:
        args = parse_args()
        if args.output == '':
            run_id = os.path.basename(args.input).split('.')[-1]
            args.output = os.path.join(os.path.dirname(args.input),
                'true_vcf.{}'.format(run_id))
        haplotypes_to_vcf(args.input, args.output)