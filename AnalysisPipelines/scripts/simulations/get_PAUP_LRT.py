#!/usr/bin/env python3

import os
import re
import argparse
import tempfile
import subprocess
from utils import get_sample_dict_from_vcf, change_tree_root


NEXUS_TEMPLATE = """#NEXUS

begin data;
    dimensions ntax={sample_no} nchar={rec_no};
    format datatype=dna missing=? gap=-;
    matrix
{matrix}
    ;
end;

begin trees;
    {tree_str}
end;

begin PAUP;
    Set autoclose=yes warnreset=no warntree=no warntsave=No;
    Set criterion=like;
    Outgroup {outg};
    exclude constant;
    LSet nst=1 rates=gamma shape=est ncat=4 condvar=no;

    ClockChecker;
    quit;
end;
    
"""


def vcf_to_nex(vcf_file, tree, out_file, paup_exe=None, exclude='', include='',
            outg='healthycell'):

    samples, sample_names = get_sample_dict_from_vcf(vcf_file, include=include,
        exclude=exclude)

    mat_str = ''
    for sample_idx, genotypes in samples.items():
        mat_str += '{}    {}\n'.format(sample_names[sample_idx],  genotypes)
    sample_names = [j for i, j in enumerate(sample_names) if i in samples]

    tree_str_rooted, tree_str_unrooted = change_tree_root(tree, sample_names)

    tree_str = f'tree simulatedTree = [&R] {tree_str_rooted}'

    nex_str = NEXUS_TEMPLATE.format(sample_no=len(sample_names),
        rec_no=len(genotypes), matrix=mat_str.strip('\n'), tree_str=tree_str,
        outg=outg)

    nex_file = tempfile.NamedTemporaryFile(delete=False)
    nex_file.write(str.encode(nex_str))
    nex_file.close()

    shell_cmd = ' '.join([paup_exe, '-n', nex_file.name, '>', out_file])
    paup = subprocess.Popen(shell_cmd, shell=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    stdout, stderr = paup.communicate()
    paup.wait()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='Input file')
    parser.add_argument('tree', type=str, help='Newick tree file')
    parser.add_argument('-o', '--output', type=str, help='Output file.')
    parser.add_argument('-e', '--exe', type=str, help='Path to PAUP exe.')
    parser.add_argument('-ex', '--exclude', type=str, default='',
        help='Regex pattern for samples to exclude from LRT test,')
    parser.add_argument('-in', '--include', type=str, default='',
        help='Regex pattern for samples to include from LRT test,')
    parser.add_argument('--outg', type=str, default='healthycell',
        help='Outgroup for PAUP block')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        if snakemake.params.exclude == None:
            snakemake.params.exclude = ''
        if snakemake.params.include == None:
            snakemake.params.include = ''

        vcf_to_nex(snakemake.input.vcf, snakemake.input.tree, snakemake.output[0],
            paup_exe=snakemake.params.paup_exe, exclude=snakemake.params.exclude,
            include=snakemake.params.include)
    else:
        args = parse_args()
        vcf_to_nex(args.input, args.tree, args.output, paup_exe=args.exe,
            exclude=args.exclude, include=args.include, outg=args.outg)