#!/usr/bin/env python3

import os
import re
import argparse
import tempfile
import subprocess
from utils import get_sample_dict_from_vcf, change_newick_tree_root


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

begin mrbayes;
    set autoclose=yes nowarnings=yes;
    lset nst=1 rates=gamma;
    lset ngammacat=2;
    lset rates=propinv;
    prset pinvarpr=fixed(0);
    outgroup healthycell;
    prset brlenspr={brlen_prior};
    {fixed_tree}
    {alg} ngen={ngen};
    sump outputname={mrbayes_out};
end;

begin PAUP;
    Set autoclose=yes warnreset=no warntree=no warntsave=No;
    Set criterion=like;
    {outg_cmd}
    {paup_const}
    LSet nst=1 rates=gamma shape=est ncat=4 condvar={paup_corr};

    {paup_tree}
    quit;
end;
    
"""


def get_sample_dict_from_FG(vcf_file):
    ref = {}
    with open(vcf_file, 'r') as f_in:
        for line in f_in:
            if line.startswith('#') or line.strip() == '':
                continue
            line_cols = line.strip().split('\t')
            pos = int(line.strip().split('\t')[1])
            ref_base = line.strip().split('\t')[3]
            ref[pos] = ref_base

    run_no = re.search('\.(\d+)\.', os.path.basename(vcf_file)).group(1)
    base_dir = os.path.sep.join(vcf_file.split(os.path.sep)[:-3])
    FG_file = os.path.join(base_dir, 'full_genotypes_dir',
        'full_gen.{}'.format(run_no))

    with open(FG_file, 'r') as f_in:
        samples = {}
        sample_names = []
        for i, line in enumerate(f_in):
            if i == 0:
                continue
            sample_name, BPs = line.split('  ')
            GT = ''
            for j, BP in enumerate(BPs.split(' '), 1):
                # Reference or no mutation detected
                if BP[0] == BP[1] or j not in ref:
                    GT += BP[0]
                elif BP[0] == ref[j]:
                    GT += BP[1]
                else:
                    GT += BP[0]


            if sample_name.startswith('cell'):
                sample_names.append(sample_name.replace('cell', 'tumcell'))
            else:
                sample_names.append('healthycell')

            samples[i - 1] = GT
    return samples, sample_names


def vcf_to_nex(vcf_file, out_files, ngen, ss_flag=False, tree=None,
            paup_exe=None, learn_tree=False, full_GT=False, GT=False, exclude='',
            include='', outg='healthycell'):

    if full_GT:
        samples, sample_names = get_sample_dict_from_FG(vcf_file)
        paup_corr = ('no', '')
    else:
        samples, sample_names = get_sample_dict_from_vcf(vcf_file,
            GT=GT, include=include, exclude=exclude)
        paup_corr = ('no', 'exclude constant;')

    mat_str = ''
    for sample_idx, genotypes in samples.items():
        mat_str += '{}    {}\n'.format(sample_names[sample_idx],  genotypes)
    sample_names = [j for i, j in enumerate(sample_names) if i in samples]
    rec_no = len(genotypes)

    if ss_flag:
        alg = 'ss'
    else:
        alg = 'mcmc'

    if tree == None:
        tree_rooted = ''
        tree_unrooted = ''
        fixed_tree = ''
    else:
        if 'scite' in tree:
            tree_str_rooted, tree_str_unrooted = change_newick_tree_root(
                tree, paup_exe, root=False, sample_names=sample_names,
                outg=outg)
            tree_name = 'treeSCITE'
        elif 'cellphy' in tree:
            tree_str_unrooted, tree_str_rooted = change_newick_tree_root(
                tree, paup_exe, root=True, outg=outg)
            tree_name = 'treeCellPhy'
        else:
            tree_str_rooted, tree_str_unrooted = change_newick_tree_root(
                tree, paup_exe, root=False, outg=outg)
            tree_name = 'treeCellCoal'

        tree_rooted = 'tree {} = [&R] {}'.format(tree_name, tree_str_rooted)
        tree_unrooted = 'tree {} = [&U] {}'.format(tree_name, tree_str_unrooted)
        fixed_tree = 'prset topologypr=fixed({});'.format(tree_name) 

    if learn_tree:
        paup_tree_raw = 'Lset clock=no;\n' \
            '    Hsearch;\n' \
            '    RootTrees rootMethod=outgroup outroot=monophyl;\n'
    else:
        paup_tree_raw = ''

    if outg in tree_str_rooted:
        outg_cmd = 'Outgroup ;'.format(outg)
    else:
        outg_cmd = ''


    for out_file in out_files:
        nxs_dir, nxs_file = os.path.split(out_file)

        model = nxs_file.split('.')[-1]
        tree_alg = tree_name[4:].lower()

        mrbayes_out = os.path.join('mrbayes_dir_{}'.format(tree_alg), 
            nxs_file.replace('nxs', 'mrbayes'))
        paup_dir = os.path.join('..', 'paup_dir_{}'.format(tree_alg))
        paup_out = os.path.join(paup_dir, nxs_file.replace('nxs', 'paup'))
        os.makedirs(os.path.join(nxs_dir, paup_dir), exist_ok=True)

        if model == 'clock':
            brlen_prior = 'clock:uniform'
            clock_str='yes'
            tree_str = tree_rooted
        else:
            brlen_prior = 'unconstrained:exp(1.0)'
            clock_str='no'
            tree_str = tree_unrooted

        if learn_tree:
            paup_tree = paup_tree_raw + \
                'log file={out_file}.PAUP.score start=yes replace=yes;\n' \
                '    clockChecker tree=all lrt=yes;\n' \
                '    log stop;'.format(paup_out)
        else:
            paup_tree = paup_tree_raw + \
                'LSet clock={clock_str};\n' \
                '    lscores 1/ clock={clock_str} scorefile={out_file}.PAUP.score append;'\
                .format(clock_str=clock_str, out_file=paup_out)

        nex_str = NEXUS_TEMPLATE.format(sample_no=len(sample_names),
            rec_no=rec_no, sample_labels=' '.join(sample_names),
            matrix=mat_str.strip('\n'), tree_str=tree_str, fixed_tree=fixed_tree,
            mrbayes_out=mrbayes_out, alg=alg, ngen=ngen, brlen_prior=brlen_prior, 
            paup_tree=paup_tree, outg_cmd=outg_cmd, paup_corr=paup_corr[0],
            paup_const=paup_corr[1])

        with open(out_file, 'w') as f_out:
            f_out.write(nex_str)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='Input file')
    parser.add_argument('-t', '--tree', type=str, default='',
        help='Newick tree file')
    parser.add_argument('-o', '--output', type=str, help='Output directory.')
    parser.add_argument('-n', '--ngen', type=int, default=1e6,
        help='Number of MCMC/SS steps in NEXUS MrBayes block. Default = 1e6.')
    parser.add_argument('-ss', '--stepping_stone', action='store_true',
        help='Use stepping stone sampling instead of MCMC.')
    parser.add_argument('-lt', '--learn_tree', action='store_true',
        help='Learn the tree under no constraints model.')
    parser.add_argument('-fg', '--full_GT', action='store_true',
        help='Use the full genotype instead of only SNPs for inference.')
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
        try:
            ngen = snakemake.wildcards.ngen
        except AttributeError:
            ngen = 1e6
        if snakemake.params.exclude == None:
            snakemake.params.exclude = ''
        if snakemake.params.include == None:
            snakemake.params.include = ''

        vcf_to_nex(snakemake.input.vcf, snakemake.output, ngen,
            ss_flag=snakemake.params.ss, tree=snakemake.input.tree,
            paup_exe=snakemake.params.paup_exe, exclude=snakemake.params.exclude,
            include=snakemake.params.include)
    else:
        args = parse_args()
        out_files = [os.path.join(args.output, 'nxs.{}.{}'.format(args.ngen, i))
            for i in ['clock', 'noClock']]
        vcf_to_nex(args.input, out_files, args.ngen,
            ss_flag=args.stepping_stone, tree=args.tree,
            learn_tree=args.learn_tree, full_GT=args.full_GT, paup_exe=args.exe,
            exclude=args.exclude, include=args.include, outg=args.outg)