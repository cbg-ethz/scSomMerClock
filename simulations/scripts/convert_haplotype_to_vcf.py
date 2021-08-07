#!/usr/bin/env python3

import os
import argparse


VCF_TEMPLATE = """##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Called genotype">
##FORMAT=<ID=TG,Number=1,Type=String,Description="True genotype">
##contig=<ID=1,length={contig_length}>
##source=Haplotypes from CellCoal simulation
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


def get_gt(ref_als, als, alts):
    al_m, al_p = als
    ref_al_m, ref_al_p = ref_als
    # Missing
    if al_m == 'N' and al_p == 'N':
        return '.|.'
    # Homozygous
    elif al_m == al_p and al_m != ref_al_m and al_p != ref_al_p:
        if al_m != ref_al_m and al_m not in alts:
            alts.append(al_m)
        return '{i}|{i}'.format(i=alts.index(al_m))
    # Heterozygous maternal allele
    elif al_m != ref_al_m and al_p == ref_al_p:
        if not al_m in alts:
            alts.append(al_m)
        return '0|{i}'.format(i=alts.index(al_m))
    # Heterozygous paternal allele
    elif al_p != ref_al_p and al_m == ref_al_m:
        if not al_p in alts:
            alts.append(al_p)
        return '0|{i}'.format(i=alts.index(al_p))
    # Both alleles not ref
    elif al_m != ref_al_m and al_p != ref_al_p:
        if not al_m in alts:
            alts.append(al_m)
        if not al_p in alts:
            alts.append(al_p)
        return '|'.join(sorted([str(alts.index(al_m)), str(alts.index(al_p))]))
    # Wildtype
    else:
        return '|'.join([str(alts.index(al_m)), str(alts.index(al_p))])


def save_true_vcf(haps, out_file):
    SNVs = set([])
    ref = haps.pop('ingrroot')
    cells = []
    for cell, cell_dict in haps.items():
        cells.append(cell)
        for al, bp_dict in cell_dict.items():
            for pos, base in enumerate(bp_dict['true']):
                if base != ref[al][pos] or base != bp_dict['full'][pos]:
                    SNVs.add(pos)

    out_str = ''
    stats = [0, 0, 0, 0, 0] # TP, FP, TN, FN, MS
    hom_errors = 0
    for i in sorted(list(SNVs)):
        ref_als = ''.join(sorted([ref["m"][i], ref["p"][i]]))

        if ref_als[0] == ref_als[1]:
            alts = [ref_als[0]]
        else:
            alts = [ref_als[0], ref_als[1]]
        gts = []
        true_gts = []
        for cell, cell_dict in haps.items():
            als = ''.join(
                sorted([cell_dict["m"]["full"][i], cell_dict["p"]["full"][i]]))
            true_als = ''.join(
                sorted([cell_dict["m"]["true"][i], cell_dict["p"]["true"][i]]))

            # MS
            if als[0] == 'N' or als[1] == 'N':
                stats[4] += 1
            # TP or TN
            elif als[0] == true_als[0] and als[1] == true_als[1]:
                # Same as ref: TN
                if al[0] == ref_als[0] and als[1] == ref_als[1]:
                    stats[2] += 1
                # Different from ref: TP
                else:
                    stats[0] += 1
            # FP or FN
            else:
                # True hap is same as ref but mutation called: FP
                if true_als[0] == ref_als[0] and true_als[1] == ref_als[1]:
                    stats[1] += 1
                # Ref called although the true hap should be different: FN
                elif als[0] == ref_als[0] and als[1] == ref_als[1]:
                    stats[3] += 1
                # True hap is different from ref but also different mut called: TP
                elif true_als[0] == ref_als[0] or true_als[1] == ref_als[1]:
                    stats[0] += 1
                    hom_errors += 1
                # True is different from ref: FN
                else:
                    print(ref_als, true_als, als)
                    stats[3] += 1
            gts.append(get_gt(ref_als, als, alts))
            true_gts.append(get_gt(ref_als, true_als, alts))

        alt_str = ','.join(alts[1:])
        format_str = '\t'.join([f'{i[0]}:{i[1]}' for i in zip(gts, true_gts)])
        out_str += '1\t{p}\tsnv{p:06d}\t{r}\t{a}\t.\tPASS\t.\tGT:TG\t{f}\n' \
                .format(p=i+1, r=ref_als[0], a=alt_str, f=format_str)

    bp_total = len(ref['m']) * len(cells)
    stats[2] = bp_total - sum(stats)
    stats = [i / bp_total for i in stats]

    for key, val in zip(['TP', 'FP', 'TN', 'FN', 'MS'], stats):
        print(f'{key}:\t{val:.4f}')
    print(f'FP:\t{hom_errors / bp_total:.4f}\t(het -> hom)')
    import pdb; pdb.set_trace()
    with open(out_file, 'w') as f:
        f.write(VCF_TEMPLATE\
            .format(contig_length=len(ref['m']), rows=out_str.rstrip(),
                cells='\t{}'.format('\t'.join(cells))))


def haplotypes_to_vcf(full_hap_file, true_hap_file, out_file):
    run_no = os.path.basename(true_hap_file).split('.')[-1]

    haps = {}
    with open(full_hap_file, 'r') as f:
        full_haps_raw = f.read().strip().split('\n')
        for cell_hap in full_haps_raw[1:]:
            cell_name, cell_bases = [i for i in cell_hap.split('  ') if i]
            if cell_name.startswith('cell') or cell_name.startswith('tumcell'):
                al = cell_name[-1]
                try:
                    haps[cell_name[:-1]][al] = {'full': cell_bases}
                except KeyError:
                    haps[cell_name[:-1]] = {al: {'full': cell_bases}}

    with open(true_hap_file, 'r') as f:
        true_haps_raw = f.read().strip().split('\n')
        for cell_hap in true_haps_raw[1:]:
            cell_name, cell_bases = [i for i in cell_hap.split('  ') if i]
            if cell_name.startswith('cell') or cell_name.startswith('tumcell'):
                al = cell_name[-1]
                haps[cell_name[:-1]][al]['true'] = cell_bases
            elif cell_name.startswith('ingrroot'):
                al = cell_name[-1]
                try:
                    haps[cell_name[:-1]][al] = cell_bases
                except KeyError:
                    haps[cell_name[:-1]] = {al: cell_bases}

    out_dir = os.path.dirname(out_file)
    if not os.path.exists(out_dir):
        os.path.mkdir(out_dir)

    # save_fasta(ref, out_dir, run_no)
    # save_germline_vcf(outgr, ref, out_dir, run_no)
    save_true_vcf(haps, out_file)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--full', type=str, help='Full haplotype file')
    parser.add_argument('-t', '--true', type=str, help='True haplotype file')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        haplotypes_to_vcf(snakemake.input[0], snakemake.output[0])
    else:
        args = parse_args()
        if args.output:
            out_file = args.output
        else:
            run_id = os.path.basename(args.true).split('.')[-1]
            out_file = os.path.join(os.path.dirname(args.true),
                'true_vcf.{}'.format(run_id))
            print(f'Save output to: {out_file}')
        haplotypes_to_vcf(args.full, args.true, out_file)