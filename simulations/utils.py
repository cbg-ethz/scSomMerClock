#!/usr/bin/env python3

import os
import re
import gzip
import math
import argparse
import subprocess
from pprint import pprint as pp
from statistics import mean, stdev

TICK_FONTSIZE = 12
LABEL_FONTSIZE = 16


BASES = ['A', 'C', 'G', 'T']
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
    lset nst=6 rates=invgamma;
    [lset ngammacat=2;]
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
    Outgroup healthycell;
    {paup_const}
    LSet nst=1 base=equal rates=gamma shape=est condvar={paup_corr};

    {paup_tree}
    quit;
end;
    
"""


VCF_TEMPLATE = """##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="True genotype">
##contig=<ID=1,length={contig_length}>
##source=germline variants from CellCoal simulation
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT{cells}
{rows}
"""


def get_out_dir(config):
    if config['cellcoal']['model'].get('branch_rate_var', 0):
        model = 'clock{}'.format(config['cellcoal']['model']['branch_rate_var'])
    else:
        model = 'clock0'

    if config['cellcoal']['scWGA'].get('errors', False):
        sim_scWGA = 'WGA{}-{}-{}'.format(
            config['cellcoal']['scWGA']['ADO_rate'],
            config['cellcoal']['scWGA']['doublet_rate'][0],
            config['cellcoal']['scWGA']['ampl_error'][0])
    else:
        sim_scWGA = 'WGA0-0-0'

    sim_NGS = 'NGS{}-'.format(config['cellcoal']['NGS']['seq_cov'])
    if config['cellcoal']['NGS'].get('errors', False):
        sim_NGS += '{}-{}'.format(config['cellcoal']['NGS']['seq_overdis'],
            config['cellcoal']['NGS']['seq_error'])
    else:
        sim_NGS += '0-0'

    if not config.get('static', {}).get('out_dir', False):
        out_dir = os.path.dirname(os.path.relpath(__file__))
    else:
        out_dir = config['static']['out_dir']

    sim_dir = os.path.join(out_dir, 'res_{}_{}_{}' \
        .format(model, sim_scWGA, sim_NGS))

    filter_dir =  'minDP{}-minGQ{}'.format(
        config.get('SNP_filter', {}).get('depth', 1),
        config.get('SNP_filter', {}).get('quality', 0))
    if config.get('SNP_filter', {}).get('singletons', False):
        filter_dir += '-noSingletons'

    return sim_dir, filter_dir


def get_cellcoal_config(config, template_file, out_dir):
    with open(template_file, 'r') as f:
        templ = f.read()

    templ = re.sub('{no_rep}', str(config['cellcoal'].get('no_rep', 10)), templ) 
    templ = re.sub('{no_cells}',
        str(config['cellcoal']['model'].get('no_cells', 30)), templ)
    templ = re.sub('{no_sites}',
        str(config['cellcoal']['model'].get('no_sites', 10000)), templ)
    templ = re.sub('{pop_size}',
        str(config['cellcoal']['model'].get('pop_size', 10000)), templ)
    templ = re.sub('{out_dir}', out_dir, templ)
    templ = re.sub('{seq_cov}',
        str(config['cellcoal'].get('NGS', {}).get('seq_cov', 20)), templ)
    templ = re.sub('{outgroup_branch_length}',
        str(config['cellcoal']['model'].get('outgroup_branch_length', 1)), templ)

    if config['cellcoal']['model'].get('no_muts', None):
        templ = re.sub('{no_muts}',
            'j{}'.format(config['cellcoal']['model']['no_muts']), templ)
    else:
        templ = re.sub('{no_muts}', '', templ)

    if config['cellcoal']['model'].get('branch_rate_var', None):
        templ = re.sub('{branch_rate_var}',
            'i{}'.format(config['cellcoal']['model']['branch_rate_var']), templ)
    else:
        templ = re.sub('{branch_rate_var}', '', templ)

    if config['cellcoal']['model'].get('germline_rate', 0.0) > 0:
        templ = re.sub('{germline_rate}',
            'c{}'.format(config['cellcoal']['model']['germline_rate']), templ)
    else:
        templ = re.sub('{germline_rate}', '', templ)

    if config['cellcoal']['model'].get('binary_alphabet', False):
        templ = re.sub('{alphabet}', '0', templ)
    else:
        templ = re.sub('{alphabet}', '1', templ)

    if config['cellcoal']['scWGA'].get('errors', False):
        templ = re.sub('{ADO_rate}',
            'D{}'.format(config['cellcoal']['scWGA']['ADO_rate']), templ)
        templ = re.sub('{ampl_error}',
            'A{} {} {}'.format(*config['cellcoal']['scWGA']['ampl_error']), templ)
        templ = re.sub('{doublet_rate}',
            'B{} {}'.format(*config['cellcoal']['scWGA']['doublet_rate']), templ)
    else:
        templ = re.sub('{ADO_rate}', '', templ)
        templ = re.sub('{ampl_error}', '', templ)
        templ = re.sub('{doublet_rate}', '', templ)

    if config['cellcoal']['NGS'].get('errors', False):
        templ = re.sub('{seq_overdis}',
            'V{}'.format(config['cellcoal']['NGS']['seq_overdis']), templ)
        templ = re.sub('{seq_error}',
            'E{}'.format(config['cellcoal']['NGS']['seq_error']), templ)
    else:
        templ = re.sub('{seq_overdis}', '', templ)
        templ = re.sub('{seq_error}', '', templ)

    templ = re.sub('{out_tree}', '6', templ)

    if config['cellcoal'].get('output', {}).get('full_GT', False) \
            or config.get('paup', {}).get('full_GT', False) \
            or config.get('sieve', {}).get('run', False):
        templ = re.sub('{out_full_GT}', '3', templ)
    else:
        templ = re.sub('{out_full_GT}', '', templ)

    if config['cellcoal'].get('output', {}).get('true_haplotype', False):
        templ = re.sub('{out_true_hap}', '9', templ)
    else:
        templ = re.sub('{out_true_hap}', '', templ)

    return templ


def load_cellcoal_config(configfile):
    cc_params = {}
    with open(configfile, 'r') as f:
        for line in f.read().strip().split('\n'):
            pair = line.strip().split('] ')
            if len(pair) < 2 or pair[1].startswith('['):
                continue
            cc_params[pair[0].lstrip('[')] = pair[1][1:]
    return cc_params


def vcf_to_pileup(vcf_file, out_pileup, out_samples='', out_sample_types=''):
    if vcf_file.endswith('gz'):
        file_stream = gzip.open(vcf_file, 'rb')
    else:
        file_stream = open(vcf_file, 'r')

    pileup = ''
    with file_stream as f_in:
        for line in f_in:
            # Skip VCF header lines
            if line.startswith('#'):
                # Safe cell/sample names
                if line.startswith('#CHROM'):
                    sample_names = line.strip().split('\t')[9:]
                continue
            # VCF records
            cols = line.strip().split('\t')
            ref = cols[3]

            new_pileup = ''
            for s_rec in cols[9:]:
                s_rec_format = s_rec.split(':')
                try:
                    DP = int(s_rec_format[1])
                except ValueError:
                    DP = 0
                    import pdb; pdb.set_trace()
                else:
                    if DP == 0:
                        new_pileup += '0\t*\t*\t'
                        continue
                s_bases = ''
                for base_i, base_no in enumerate(s_rec_format[2].split(',')):
                    if BASES[base_i] == ref:
                        s_bases += '.' * int(base_no)
                    else:
                        s_bases += BASES[base_i] * int(base_no)

                new_pileup += '{}\t{}\t{}\t'.format(DP, s_bases, '~' * DP)
            pileup += '{}\t{}\t{}\t{}\n' \
                .format(cols[0], cols[1], ref, new_pileup.rstrip())
    
    with open(out_pileup, 'w') as f_pileup:
        f_pileup.write(pileup.rstrip())

    if not out_samples:
        out_samples = '{}.SampleNames.txt'.format(vcf_file)
    with open(out_samples, 'w') as f_names:
        f_names.write('\n'.join(sample_names))

    if not out_sample_types:
        out_samples = '{}.SampleTypes.txt'.format(vcf_file)
    with open(out_samples, 'w') as f_names:
        f_names.write('\n'.join(
            ['{}\tCT'.format(i) if not i.startswith('healthy') \
                else '{}\tCN'.format(i) for i in sample_names])
        )


def postprocess_vcf(vcf_file, out_file, minDP=1, minGQ=0, s_filter=False):
    import numpy as np

    if vcf_file.endswith('gz'):
        file_stream = gzip.open(vcf_file, 'rb')
    else:
        file_stream = open(vcf_file, 'r')

    header = ''
    body = ''

    monovar = False
    with file_stream as f_in:
        for line in f_in:
            # Skip VCF header lines
            if line.startswith('#'):
                # Safe column headers
                if line.startswith('##source') and 'MonoVar' in line:
                    monovar = True
                elif line.startswith('#CHROM'):
                    header += '##FILTER=<ID=singleton,Description="SNP ' \
                        'is singleton">\n'
                    sample_no = len(line.strip().split('\t')[9:])
                    if monovar:
                        format_short = 'GT:AD:DP:GQ:PL'
                        missing = '.|.:.,.:.:.:.'
                    else:
                        header += '##FORMAT=<ID=GQ,Number=1,Type=Integer,' \
                            'Description="Genotype Quality">\n'
                        format_short = 'GT:DP:RC:G10:PL:GQ:TG'
                        missing = '.|.:.:.,.,.,.:.:.:.:'

                header += line
                continue
            elif line.strip() == '':
                break
            
            # VCF records
            line_cols = line.strip().split('\t')

            ref = line_cols[3]
            alts = line_cols[4].split(',')
            FORMAT_col = line_cols[8].split(':')
            if monovar:
                if line_cols[6] == '.':
                    line_cols[6] = 'PASS'
            else:
                PLN_col = FORMAT_col.index('PLN')

            new_line = np.zeros(sample_no, dtype=object)
            genotypes = np.chararray(sample_no, itemsize=3)

            for s_i, s_rec_raw in enumerate(line_cols[9:]):
                s_rec = s_rec_raw.split(':')

                if monovar:
                    gq = int(s_rec[3])
                    dp = s_rec[2]
                    tg = ''
                else:
                    gt = s_rec[0]
                    dp = s_rec[1]
                    rc = s_rec[2]
                    g10 = s_rec[3]
                    tg = s_rec[-1]
                    if len(FORMAT_col) == len(s_rec):
                        pln = np.array([-float(i) \
                            for i in s_rec[PLN_col].split(',')])
                        gq = min(99, sorted(pln - pln.min())[1])
                        pl = ','.join([str(max(int(i), -999)) \
                            for i in s_rec[PLN_col - 1].split(',')])
                    else:
                        gq = -1
                        pl = '.'
                
                if int(dp) < minDP or gq < minGQ:
                    new_line[s_i] = missing + tg
                else:
                    if monovar:
                        new_line[s_i] = s_rec_raw
                    else:
                        new_line[s_i] = '{}:{}:{}:{}:{}:{:.0f}:{}' \
                            .format(gt, dp, rc, g10, pl, gq, tg)

            filter_str = line_cols[6]
            if s_filter:
                diff_gt, diff_count = np.unique(genotypes[:-1], return_counts=True)
                # Only one genotype detected
                if len(diff_gt) == 1:
                    if diff_gt[0] == b'0|0':
                        filter_str = 'singleton'
                # Two different genotypes detected
                elif len(diff_gt) == 2:
                    # But one out of the two is just detected in one cell
                    if min(diff_count) == 1:
                        filter_str = 'singleton'
                else:
                    if not any([i > 1 for i in sorted(diff_count)[:-1]]):
                        filter_str = 'singleton'

            body += '\t'.join(line_cols[:6]) + '\t{}\t{}\t{}\t' \
                    .format(filter_str, line_cols[7], format_short) \
                + '\t'.join(new_line) + '\n'

    with open(out_file, 'w') as f_out:
        f_out.write('{}{}'.format(header, body))


def change_newick_tree_root(in_file, out_file, paup_exe, root=True,
            outg='healthycell'):

    in_file_name = os.path.basename(in_file)
    out_file_name = os.path.basename(out_file)
    run_no = re.search('\d{4}', in_file_name).group()
    res_dir = os.path.dirname(in_file)

    if 'scite' in in_file_name:
        with open(in_file, 'r') as f_scite:
            tree = f_scite.read().strip()
        sem_count = tree.count(';')
        if sem_count == 0:
            with open(in_file, 'a') as f_scite_new:
                f_scite_new.write(';')
        elif sem_count > 1:
            with open(in_file, 'a') as f_scite_new:
                out_str = tree.replace(';', '').strip() + ';'
                f_scite_new.write(out_str)
    if root:
        root_cmd = 'DerootTrees;\nRootTrees rootMethod=outgroup outroot=monophyl'
        root = 'yes'
        paup_file = os.path.join(res_dir, 'paup_rooting.{}'.format(run_no))
    else:
        root_cmd = 'DerootTrees;\noutgroup {}'.format(outg)
        root = 'no'
        paup_file = os.path.join(res_dir, 'paup_unrooting.{}'.format(run_no))
    paup_cmd = 'getTrees file={i};\n' \
        'outgroup {g};\n' \
        '{c};\n' \
        'saveTrees format=Newick root={r} file={o};\n' \
        'quit;'.format(i=in_file_name, g=outg, c=root_cmd, r=root, o=out_file_name)

   
    with open(paup_file, 'w') as f_paup:
        f_paup.write(paup_cmd)

    shell_cmd = ' '.join([paup_exe, '-n', paup_file, '>', '/dev/null'])
    paup = subprocess.Popen(shell_cmd, shell=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    stdout, stderr = paup.communicate()
    paup.wait()

    if stderr:
        for i in str(stdout).split('\\n'):
            print(i)
        raise RuntimeError('PAUP error im command: {}'.format(shell_cmd))
    
    os.remove(paup_file)



def get_tree(tree_file, sample_names=[]):
    with open(tree_file, 'r') as f_tree:
        tree = f_tree.read().strip()

    if 'scite_dir' in tree_file:
        for s_i, s_name in enumerate(sample_names):
            pat = '(?<=[\(\),]){}(?=[,\)\)])'.format(s_i + 1)
            tree = re.sub(pat, s_name, tree)
    elif 'trees_dir' in tree_file:
        tree = tree.replace('cell', 'tumcell') \
            .replace('outgtumcell', 'healthycell')
    
    return tree


def vcf_to_nex(vcf_file, out_files, ngen, ss_flag=False, trees=None,
            learn_tree=False, full_GT=False):
    if full_GT:
        samples, sample_names = get_sample_dict_from_FG(vcf_file)
        paup_corr=('yes', '')
    else:
        samples, sample_names = get_sample_dict_from_vcf(vcf_file)
        paup_corr=('no', 'exclude constant;')

    mat_str = ''
    for sample_idx, genotypes in samples.items():
        mat_str += '{}    {}\n'.format(sample_names[sample_idx],  genotypes)
    rec_no = len(samples[0])

    if ss_flag:
        alg = 'ss'
    else:
        alg = 'mcmc'

    if trees == None:
        tree_rooted = ''
        tree_unrooted = ''
        fixed_tree = ''
    else:
        if 'scite_dir' in trees[0]:
            tree_name = 'treeSCITE'
        elif 'trees_dir' in trees[0]:
            tree_name = 'treeCellCoal'
        elif 'cellphy_dir' in trees[0]:
            tree_name = 'treeCellPhy'

        tree_rooted = 'tree {} = [&R] {}' \
            .format(tree_name, get_tree(trees[0], sample_names))
        tree_unrooted = 'tree {} = [&U] {}' \
            .format(tree_name, get_tree(trees[1], sample_names))

        fixed_tree = 'prset topologypr=fixed({});'.format(tree_name) 

    if learn_tree:
        paup_tree_raw = 'Lset clock=no;\n' \
            '    Hsearch;\n' \
            '    RootTrees rootMethod=outgroup outroot=monophyl;\n'
    else:
        paup_tree_raw = ''

    for out_file in out_files:
        nxs_file = os.path.basename(out_file)

        model = nxs_file.split('.')[-1]
        tree_alg = tree_name[4:].lower()

        mrbayes_out = os.path.join('mrbayes_dir_{}'.format(tree_alg), 
            nxs_file.replace('nxs', 'mrbayes'))
        paup_out = os.path.join('..', 'paup_dir_{}'.format(tree_alg), 
            nxs_file.replace('nxs', 'paup'))

        if model == 'clock':
            brlen_prior = 'clock:uniform'
            clock_str='yes'
            tree_str = tree_rooted
        else:
            brlen_prior = 'unconstrained:exp(10.0)'
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
            paup_tree=paup_tree, paup_corr=paup_corr[0], paup_const=paup_corr[1])

        with open(out_file, 'w') as f_out:
            f_out.write(nex_str)
   

def get_sample_dict_from_vcf(vcf_file, GT=False):
    if vcf_file.endswith('gz'):
        file_stream = gzip.open(vcf_file, 'rb')
    else:
        file_stream = open(vcf_file, 'r')

    full_gt_file = vcf_file.replace('vcf_dir', 'full_genotypes_dir') \
        .replace('vcf', 'full_gen')
    if os.path.exists(full_gt_file):
        get_bg = True
        mut_i = []
    else:
        get_bg = False

    if GT:
        missing = '3'
    else:
        missing = '?' 

    with file_stream as f_in:
        for line in f_in:
            # Skip VCF header lines
            if line.startswith('#'):
                # Safe column headers
                if line.startswith('#CHROM'):
                    sample_names = line.strip().split('\t')[9:]
                    samples = {int(i): '' for i in range(len(sample_names))}
                continue
            elif line.strip() == '':
                continue
            # VCF records
            line_cols = line.strip().split('\t')
            # Check if filter passed
            if not 'PASS' in line_cols[6]:
                continue

            if get_bg:
                mut_i.append(int(line_cols[1]))

            ref = line_cols[3]
            alts = line_cols[4].split(',')
            FORMAT_col = line_cols[8].split(':')

            try:
                GQ_col = FORMAT_col.index('PLN')
                monovar = False
            except ValueError:
                GQ_col = FORMAT_col.index('PL')
                monovar = True

            for s_i, s_rec in enumerate(line_cols[9:]):
                try:
                    gt = s_rec[:s_rec.index(':')]
                # Missing in Monovar output format
                except ValueError:
                    samples[s_i] += missing
                    continue

                s_rec_ref, s_rec_alt = re.split('[/\|]', gt)[:2]

                # Missing
                if s_rec_ref == '.':
                    samples[s_i] += missing
                # Wildtype
                elif s_rec_ref == '0' and s_rec_alt == '0':
                    if GT:
                        samples[s_i] += '0'
                    else:
                        samples[s_i] += ref
                else:
                    if GT:
                        samples[s_i] += '1'
                    else:
                        samples[s_i] += alts[max(int(s_rec_ref), int(s_rec_alt)) - 1]

    if get_bg:
        true_GT = tail(open(full_gt_file, 'r'), 1).strip().split(' ')[2:]
        cnts = {}
        for i, j in enumerate(true_GT):
            if len(mut_i) > 0 and i == mut_i[0] -1:
                mut_i.pop(0)
            else:
                try:
                    cnts[j] += 1
                except KeyError:
                    cnts[j] = 1

        bg_out_dir = full_gt_file.replace('full_gen.', 'background.')
        with open(bg_out_dir, 'w') as bg_out:
            bg_out.write(' ' \
                .join([str(cnts[i]) for i in ['AA', 'CC', 'GG', 'TT']]))

    return samples, sample_names


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

    FG_file = vcf_file.replace('vcf_dir', 'full_genotypes_dir') \
        .replace('vcf', 'full_gen')

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


def tail(f, lines=1, _buffer=4098):
    """From https://stackoverflow.com/questions/136168"""
    lines_found = []
    block_counter = -1

    while len(lines_found) < lines:
        try:
            f.seek(block_counter * _buffer, os.SEEK_END)
        except IOError:
            f.seek(0)
            lines_found = f.readlines()
            break

        lines_found = f.readlines()
        block_counter -= 1

    return '\n'.join(lines_found[-lines:])


def get_Bayes_factor(in_files, out_file, ss=False):
    scores = {}
    runtimes = {}
    for in_file in in_files:
        _, run, steps, model, _ = os.path.basename(in_file).split('.')
        steps = int(steps)
        with open(in_file, 'r') as f_score:
            score_raw = f_score.read().strip().split('\n')
        score = float(score_raw[-1].split('\t')[2])
        score_diff = abs(float(score_raw[-2].split('\t')[2]) - \
            float(score_raw[-3].split('\t')[2]))
        
        try:
            log_tail = tail(open(in_file.replace('.lstat', '.log'), 'r'), 1000)
        except FileNotFoundError:
            runtime = (0, -1)
        else:
            runtime_start = log_tail.index('Analysis completed in')
            runtime_end = log_tail[runtime_start:].index('seconds') + runtime_start
            runtime = re.findall('\d+', log_tail[runtime_start: runtime_end])

        if ss:
            if not log_tail:
                raise IOError('Log file with SS score not found: {}' \
                    .format(in_file.replace('.lstat', '.log')))
            ss_start = log_tail.index('Marginal likelihood (ln)')
            ss_end = log_tail[ss_start:].index('More statistics on stepping-stone')
            ss_raw = [i.strip() for i in \
                    log_tail[ss_start:ss_start + ss_end].strip().split('\n') \
                if i.strip() and not i.strip().startswith('-')]

            ss_mean = float(ss_raw[-1].split(' ')[-1])
            ss_diff = abs(float(ss_raw[-2].split(' ')[-1]) \
                - float(ss_raw[-3].split(' ')[-1]))
        else:
            ss_mean = None
            ss_diff = None

        # HH:MM:SS
        if len(runtime) == 3:
            runtime = 60 * 60 * int(runtime[0]) + 60 * int(runtime[1]) \
                + int(runtime[2])
        # MM:SS
        elif len(runtime) == 2:
            runtime = 60 * int(runtime[0]) + int(runtime[1])
        # SS
        else:
            runtime = int(runtime[0])
            

        if not steps in scores:
            scores[steps] = {run: {'clock': [-1] * 4, 'noClock': [-1] * 4}}
            runtimes[steps] = []
        if not run in scores[steps]:
            scores[steps][run] = {'clock': [-1] * 4, 'noClock': [-1] * 4}
        scores[steps][run][model] = (score, score_diff, ss_mean, ss_diff)
        runtimes[steps].append(runtime)

    out_str = ''
    summary_str = ''
    for step, step_info in sorted(scores.items()):
        summary_df = {'harmonic': {'h0': [], 'h0_diff': [], 'h1': [],
                'h1_diff': [], 'logB_01': [], 'evidence': []},
            'ss': {'h0': [], 'h0_diff': [], 'h1': [], 'h1_diff': [],
                'logB_01': [], 'evidence': []}}
        
        for run, run_info in step_info.items():
            h0 = run_info['clock']
            h1 = run_info['noClock']
            
            if h0[0] == -1 or h1[0] == -1:
                continue

            summary_df['harmonic']['h0'].append(h0[0])
            summary_df['harmonic']['h0_diff'].append(h0[1])
            summary_df['harmonic']['h1'].append(h1[0])
            summary_df['harmonic']['h1_diff'].append(h1[1])

            new_out_str, logB_01, evidence_flag = \
                get_marginal_ll_str(h0[:2], h1[:2])

            summary_df['harmonic']['logB_01'].append(logB_01)
            summary_df['harmonic']['evidence'].append(evidence_flag)

            out_str += f'{step}\t{run}\t' + new_out_str

            if ss:
                summary_df['ss']['h0'].append(h0[0])
                summary_df['ss']['h0_diff'].append(h0[1])
                summary_df['ss']['h1'].append(h1[0])
                summary_df['ss']['h1_diff'].append(h1[1])

                new_out_str_ss, logB_01_ss, evidence_flag_ss = \
                    get_marginal_ll_str(h0[2:], h1[2:])

                summary_df['ss']['logB_01'].append(logB_01_ss)
                summary_df['ss']['evidence'].append(evidence_flag_ss)

                out_str += '\t' + new_out_str_ss

            out_str += '\n'

        if len(step_info) < 2:
            continue

        summary_str += f"{step:.0E}\t{max(runtimes[step]):.2f}\t" \
            f"{mean(summary_df['harmonic']['h0']):.1f}\t" \
            f"{stdev(summary_df['harmonic']['h0']):.1f}\t" \
            f"{mean(summary_df['harmonic']['h0_diff']):.0f}\t" \
            f"{stdev(summary_df['harmonic']['h0_diff']):.0f}\t" \
            f"{mean(summary_df['harmonic']['h1']):.1f}\t" \
            f"{stdev(summary_df['harmonic']['h1']):.1f}\t" \
            f"{mean(summary_df['harmonic']['h1_diff']):.0f}\t" \
            f"{stdev(summary_df['harmonic']['h1_diff']):.0f}\t" \
            f"{mean(summary_df['harmonic']['logB_01']):.1f}\t" \
            f"{stdev(summary_df['harmonic']['logB_01']):.1f}\t" \
            f"{sum(summary_df['harmonic']['evidence'])}/" \
            f"{len(summary_df['harmonic']['evidence'])}"

        if ss:
            summary_str += f"\t{mean(summary_df['ss']['h0']):.1f}\t" \
            f"{stdev(summary_df['ss']['h0']):.1f}\t" \
            f"{mean(summary_df['ss']['h0_diff']):.0f}\t" \
            f"{stdev(summary_df['ss']['h0_diff']):.0f}\t" \
            f"{mean(summary_df['ss']['h1']):.1f}\t" \
            f"{stdev(summary_df['ss']['h1']):.1f}\t" \
            f"{mean(summary_df['ss']['h1_diff']):.0f}\t" \
            f"{stdev(summary_df['ss']['h1_diff']):.0f}\t" \
            f"{mean(summary_df['ss']['logB_01']):.1f}\t" \
            f"{stdev(summary_df['ss']['logB_01']):.1f}\t" \
            f"{sum(summary_df['ss']['evidence'])}/" \
            f"{len(summary_df['ss']['evidence'])}"
        summary_str += '\n'

    with open(out_file, 'w') as f_out:
        if ss:
            f_out.write('steps\trun\t' \
                'H0:clock (harmonic mean)\tH0 LL diff (harmonic mean)\t'
                'H1:noClock (harmonic mean)\tH1 LL diff (harmonic mean)\t' \
                '2log_e(B_01) (harmonic mean)\tB_01 (harmonic mean)\t' \
                'Evidence (harmonic mean)\t' \
                'H0:clock (ss)\tH0 LL diff (ss)\tH1:noClock (ss)\t'
                'H1 LL diff (ss)\t2log_e(B_01) (ss)\tB_01 (ss)\tEvidence (ss)\n')
        else:
            f_out.write('steps\trun\t' \
                'H0:clock (harmonic mean)\tH0 LL diff (harmonic mean)\t'
                'H1:noClock (harmonic mean)\tH1 LL diff (harmonic mean)\t' \
                '2log_e(B_01) (harmonic mean)\tB_01 (harmonic mean)\t' \
                'Evidence (harmonic mean)\n')
        f_out.write(out_str.strip())

    if summary_str:
        with open(out_file.replace('.tsv', '.short.tsv'), 'w') as f_out:
            if ss:
                f_out.write('steps\tMax. runtime [secs]\t'
                    'Avg. H0:clock (harmonic mean)\tStd. H0:clock (harmonic mean)\t'
                    'Avg. H0 LL diff (harmonic mean)\tStd. H0 LL diff (harmonic mean)\t'
                    'Avg. H1:noClock (harmonic mean)\tStd. H1:noClock (harmonic mean)\t'
                    'Avg. H1 LL diff (harmonic mean)\tStd. H1 LL diff (harmonic mean)\t'
                    'Avg. 2log_e(B_01) (harmonic mean)\t'
                    'Std. 2log_e(B_01) (harmonic mean)\tEvidence (harmonic mean)\t'
                    'Avg. H0:clock (ss)\tStd. H0:clock (ss)\tAvg. H0 LL diff (ss)\t'
                    'Std. H0 LL diff (ss)\tAvg. H1:noClock (ss)\tStd. H1:noClock (ss)\t'
                    'Avg. H1 LL diff (ss)\tStd. H1 LL diff (ss)\t'
                    'Avg. 2log_e(B_01) (ss)\tStd. 2log_e(B_01) (ss)\t'
                    'Evidence (ss)\n')
            else:
                f_out.write('steps\tMax. runtime [secs]\t'
                    'Avg. H0:clock (harmonic mean)\tStd. H0:clock (harmonic mean)\t'
                    'Avg. H0 LL diff (harmonic mean)\tStd. H0 LL diff (harmonic mean)\t'
                    'Avg. H1:noClock (harmonic mean)\tStd. H1:noClock (harmonic mean)\t'
                    'Avg. H1 LL diff (harmonic mean)\tStd. H1 LL diff (harmonic mean)\t'
                    'Avg. 2log_e(B_01)\tStd. 2log_e(B_01)\tEvidence\n')
            f_out.write(summary_str.strip())


def get_marginal_ll_str(h0, h1):
    diff = h1[0] - h0[0]
    if diff > 50:
        B01 = math.inf
    elif diff < -50:
        B01 = math.inf
    else:
        B01 = math.exp(diff)

    logB_01 = 2 * diff

    if logB_01 < 2:
        evidence = 'None'
        evidence_flag = 0
    elif logB_01 < 6:
        evidence = 'Positive'
        evidence_flag = 1
    elif logB_01 < 10:
        evidence = 'Strong'
        evidence_flag = 1
    else:
        evidence = 'Very Strong'
        evidence_flag = 1

    out_str = f'{h0[0]:.1f}\t{h0[1]:.0f}\t{h1[0]:.1f}\t{h1[1]:.0f}\t' \
        f'{logB_01:.1f}\t{B01:.0f}\t{evidence}'
    return out_str, logB_01, evidence_flag



def get_LRT(masterfile, out_file, cell_no, alpha=0.05):
    from scipy.stats.distributions import chi2

    with open(masterfile, 'r') as f_master:
        in_files = f_master.read().strip().split('\n')

    scores = {}
    for in_file in in_files:
        _, run, _, model, _, _ = os.path.basename(in_file).split('.')
        with open(in_file, 'r') as f_score:
            score_raw = f_score.read().strip().split('\n')
        try:
            score = -float(score_raw[-1].split('\t')[1])
        except (ValueError, IndexError):
            print(f'Cannot read: {in_file}')
            continue

        if not run in scores:
            scores[run] ={'clock': -1, 'noClock': -1}
        scores[run][model] = score

    out_str = ''
    avg = [[], [], [], 0, 0]
    for run, run_info in sorted(scores.items()):
        h0 = run_info['clock']
        h1 = run_info['noClock']
        if h0 == -1 or h1 == -1:
            continue

        LR = -2 * (h0 - h1)
        p_val = chi2.sf(LR, cell_no - 2)
        if p_val < alpha:
            hyp = 'H1'
            avg[4] += 1
        else:
            hyp = 'H0'
            avg[3] += 1

        avg[0].append(h0)
        avg[1].append(h1)
        avg[2].append(p_val)

        out_str += f'{run}\t{h0:.1f}\t{h1:.1f}\t{p_val:.2E}\t{hyp}\n'

    with open(out_file, 'w') as f_out:
        f_out.write('run\tH0:clock\tH1:noClock\tp-value\thypothesis\n')
        f_out.write(out_str)
        f_out.write(f'\nAvg.\t{mean(avg[0])}\t{mean(avg[1])}\t{mean(avg[2])}\t'
            f'H0:{avg[3]};H1:{avg[4]}\n')


def run_scite_subprocess(exe, steps, vcf_file, fd=0.001, ad=0.2, verbose=False):
    import numpy as np

    # RUN SiCloneFit
    if not os.path.exists(exe):
        raise RuntimeError(
            'SCITE not compiled: run "g++ *.cpp -o SCITE -std=c++11" inside '
            '{}'.format(os.path.dirname(exe))
        )

    run_no = re.search('\d\d\d\d', os.path.basename(vcf_file)).group()
    out_dir = os.path.sep.join(vcf_file.split(os.path.sep)[:-2] + ['scite_dir'])
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    data_raw, cells = get_sample_dict_from_vcf(vcf_file, GT=True)
    data_list = []
    for cell_data in data_raw.values():
        data_list.append([int(i) for i in cell_data])
    data = np.array(data_list).T

    data_file = os.path.join(out_dir, 'SCITE.{}.csv'.format(run_no))
    np.savetxt(data_file, data.astype(int), delimiter=' ', newline='\n', fmt='%d')
    no_muts, no_cells = data.shape

    mut_file = os.path.join(out_dir, 'muts.txt')
    if not os.path.exists(mut_file):
        with open(mut_file, 'w') as f_mut:
            f_mut.write('\n'.join(['m{}'.format(i) for i in range(no_muts)]))

    out_files = os.path.join(out_dir, 'scite_tree.{}'.format(run_no))
    cmmd = ' '.join(
        [exe, '-i', data_file, '-transpose', '-r 1', '-n', str(no_muts),
        '-m', str(no_cells), '-l', str(steps), '-fd', str(fd), '-ad', str(ad),
        '-e 0.1', '-a',  '-o', out_files, '-names', mut_file,
        '-max_treelist_size 1']
    )

    if verbose:
        print('output directory:\n{}'.format(out_dir))
        print('\nShell command:\n{}\n'.format(cmmd))

    SCITE = subprocess.Popen(
        cmmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = SCITE.communicate()
    SCITE.wait()

    if stderr:
        for i in str(stdout).split('\\n'):
            print(i)
        raise RuntimeError('SCITE Error')
    os.remove(data_file)


def get_sieve_tree(in_file, out_file, cells):
    with open(in_file, 'r') as f_in:
        tree = f_in.read().strip()

    if 'scite_dir' in in_file:
        for i in range(1, cells + 1, 1):
            tree = re.sub('(?<=[\(\),]){}(?=[,\)\)])'.format(i),
                'tumcell{:0>4}'.format(i), tree)
        tree = re.sub('(?<=[\(\),]){}(?=[,\)\)])'.format(cells + 1), 
            'healthycell', tree)
        tree += ';'
    elif 'trees_dir' in in_file:
        tree = tree.replace('cell', 'tumcell') \
            .replace('outgtumcell', 'healthycell')

    tree_cut = re.sub(',healthycell:\d+.\d+', '', tree[1:-2]) + ';'
    with open(out_file, 'w') as f_out:
        f_out.write(tree_cut)


def get_sieve_xml(template_file, tree_file, samples_file, model, steps,
            out_file):
    with open(template_file, 'r') as f:
        templ = f.read()

    run = re.search('sieve_tree.(\d+)', tree_file).group(1)

    # bg_file = os.path.join(os.path.sep.join(tree_file.split(os.path.sep)[:-2]),
    #     'full_genotypes_dir', 'background.{}'.format(run))

    # with open(bg_file, 'r') as f_bg:
    #     bg_bases = f_bg.read().strip()

    # background_str = '<data id="alignment" spec="FilteredAlignment" filter="-" ' \
    #     'data="@original-alignment" constantSiteWeights="{}"/>'.format(bg_bases)
    background_str = ''
    templ = re.sub('{background}', background_str, templ)

    with open(tree_file, 'r') as f_tree:
        newick = f_tree.read().strip()

        sample_no = 0
        with open(samples_file, 'r') as f_smpl:
            for s_i, s_raw in enumerate(f_smpl.read().strip().split('\n')):
                s_name, _ = s_raw.split('\t')
                newick = re.sub('(?<=[\(\),]){}(?=[,\)\)])'.format(s_i + 1), s_name, newick)
        # Replace names
        newick += ';'

    tree_node = '<tree id="tree" name="stateNode" ' \
        'nodetype="beast.evolution.tree.ScsNode" ' \
        'spec="beast.evolution.tree.ScsTree" ' \
        'treeFileName="{}"/>'.format(tree_file)

    # TODO <NB> Remove curly brackets in XML
    templ = re.sub('{branch_no}', str(2 * (s_i + 1) - 1), templ)

    # '<stateNode spec="beast.util.TreeParser" id="randomTree" ' \
    #     'IsLabelledNewick="true" adjustTipHeights="false" taxa="@alignment" ' \
    #     'newick="{}"/>'.format(newick)

    templ = re.sub('{tree}', tree_node, templ)


    templ = re.sub('{steps}', str(steps), templ)

    if model == 'clock':
        prior_node = ''
        model_node = """<branchRateModel id='strictClock' spec='beast.evolution.branchratemodel.StrictClockModel'>
                                                                        
                        <parameter estimate='false' id='clockRate' name='clock.rate' spec='parameter.RealParameter'>1.0</parameter>
                                                                    
                    </branchRateModel>"""

        op_node = ''
        log_node = ''
        br_rate_node = '<log id="TreeWithMetaDataLogger" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@tree" />'
    else:
        prior_node = """<prior id="ucldStdevPrior" name="distribution" x="@ucldStdev">
                    <Gamma id="Gamma.0" name="distr">
                        <parameter estimate="false" id="RealParameter.6" name="alpha" spec="parameter.RealParameter">0.5396</parameter>
                        <parameter estimate="false" id="RealParameter.7" name="beta" spec="parameter.RealParameter">0.3819</parameter>
                    </Gamma>
                </prior>"""
        model_node = """<branchRateModel id="RelaxedClock" spec="beast.evolution.branchratemodel.UCRelaxedClockModel" rateCategories="@rateCategories" tree="@tree">
                        <LogNormal id="LogNormalDistributionModel" S="@ucldStdev" meanInRealSpace="true" name="distr">
                            <parameter id="RealParameter.5" spec="parameter.RealParameter" estimate="false" lower="0.0" name="M" upper="1.0">1.0</parameter>
                        </LogNormal>
                        <parameter id="ucldMean" spec="parameter.RealParameter" estimate="false" name="clock.rate">1.0</parameter>
                    </branchRateModel>"""
        op_node = """<operator id="ucldStdevScaler" spec="ScaleOperator" parameter="@ucldStdev" scaleFactor="0.5" weight="3.0"/>

        <operator id="CategoriesRandomWalk" spec="IntRandomWalkOperator" parameter="@rateCategories" weight="5.0" windowSize="1"/>

        <operator id="CategoriesSwapOperator" spec="SwapOperator" intparameter="@rateCategories" weight="5.0"/>

        <operator id="CategoriesUniform" spec="UniformOperator" parameter="@rateCategories" weight="5.0"/>
"""

        log_node = """<log idref="ucldStdev"/>
            <log id="rate" spec="beast.evolution.branchratemodel.RateStatistic" branchratemodel="@RelaxedClock" tree="@tree"/>
"""
        br_rate_node = '<log id="TreeWithMetaDataLogger" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@tree" branchratemodel="@RelaxedClock"/>'

    templ = re.sub('{priors}', prior_node, templ)
    templ = re.sub('{model}', model_node, templ)
    templ = re.sub('{relaxed_clock_op}', op_node, templ)
    templ = re.sub('{relaxed_clock_log}', log_node, templ)
    templ = re.sub('{br_rate_log}', br_rate_node, templ)

    with open(out_file, 'w') as f_out:
        f_out.write(templ)


def generate_mrbayes_plots(in_file, out_file, regress=False):
    import numpy as np
    import pandas as pd
    from scipy.stats import linregress
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from matplotlib.ticker import FormatStrFormatter

    def get_ratio(ev_str):
        success, total = ev_str.split('/')
        return float(success) / float(total)

    df = pd.read_csv(in_file, sep='\t')
    df['ratio'] = df['Evidence'].apply(lambda x: get_ratio(x))
    df['runtime'] = df['Avg. runtime [secs]'].apply(lambda x: x / 60 / 60)

    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(3, 1)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[2, 0])

    ax0.errorbar(df['steps'], df['Avg. H_0:clock'], yerr= df['Std. H_0:clock'],
        label=r'$H_0$', color='#1F78B4', capsize=4, ls='--', marker='x')
    ax0.errorbar(df['steps'], df['Avg. H_1:noClock'], yerr= df['Std. H_1:noClock'],
        label=r'$H_1$', color='#33A02C', capsize=4, ls='--', marker='x')
    ax0.set_ylabel('Log-likelihood', fontsize=LABEL_FONTSIZE)
    ax0.set_xlabel('MCMC steps', fontsize=LABEL_FONTSIZE)
    ax0.xaxis.set_major_formatter(FormatStrFormatter('%.0E'))
    ax0.legend(fontsize=TICK_FONTSIZE)

    ax1.errorbar(df['steps'], df['Avg. 2log_e(B_01)'],
        yerr= df['Std. 2log_e(B_01)'], color='#FF7F00', capsize=4, ls='--',
        marker='x')
    ax1.set_ylabel(r'$2 log_e(B_{01})$', fontsize=LABEL_FONTSIZE)
    ax1.set_xlabel('MCMC steps', fontsize=LABEL_FONTSIZE)
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0E'))

    ax2.plot(df['steps'], df['ratio'], color='#E31A1C', ls='', marker='x',
        label='real')

    if regress:
        reg_line = linregress(df['steps'], df['ratio'])
        reg_x_alpha = (0.05 - reg_line.intercept) / reg_line.slope
        reg_x = np.linspace(df['steps'].min(), math.ceil(reg_x_alpha / 1e7) * 1e7, 20)
        reg_y = reg_line.intercept + reg_line.slope * reg_x

        ax2.plot(reg_x, reg_y, color='#6A3D9A', ls='--', marker='x',
            label=f'{reg_line.intercept:.2f} + {reg_line.slope:.2E} * x    '
                f'($R^2$={reg_line.rvalue:.2f})')
        ax2.axhline(0.05, ls=':')
        ax2.axvline(reg_x_alpha, ls=':')
        ax2.text(reg_x_alpha, 0.05, f'({reg_x_alpha:.2E}, 0.05)', va='top',
            ha='right', fontsize=TICK_FONTSIZE)
        ax2.legend(fontsize=TICK_FONTSIZE)

    ax2.set_ylabel(r'$H_0$ rejected [%]', fontsize=LABEL_FONTSIZE)
    ax2.set_xlabel('MCMC steps', fontsize=LABEL_FONTSIZE)
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0E'))
    ax2.set_ylim(-.05, 1.05)

    reg_line2 = linregress(df['steps'], df['runtime'])
    def mcmc_to_runtime(x):
        return reg_line2.intercept + reg_line2.slope * x
    
    secax = ax2.secondary_xaxis('top', functions=(mcmc_to_runtime, mcmc_to_runtime))
    secax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    secax.set_xlabel('Runtime [h]')   

    fig.suptitle('Simulations: no WGS, no NGS', fontsize=LABEL_FONTSIZE * 1.5)

    fig.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.92,
        hspace=0.5)
    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()
    plt.close()


def generate_pval_plot(in_file, out_file):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    df = pd.read_csv(in_file, sep='\t')    
    fig, ax = plt.subplots(figsize=(16, 12))

    ax.hist(df['p-value'], bins=100, range=(0, 1))
    ax.set_ylabel(f'counts (n={df.shape[0]})', fontsize=LABEL_FONTSIZE)
    ax.set_xlabel('p-values', fontsize=LABEL_FONTSIZE)
    
    fig.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.92,
        hspace=0.5)
    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()
    plt.close()


def convert_steps(in_dir, steps):
    runs = set([])
    for in_file in os.listdir(in_dir):
        _, run, steps, model, _ = in_file.split('.')
        if run in runs:
            break
        else:
            runs.add(run)

        for step in steps:
            pass


def parse_args():
    parser = argparse.ArgumentParser(prog='utils.py',
        usage='python3 utils.py <DATA> [options]',
        description='*** This script does one of the following: 1. VCF to NEXUS'
            '\n\t2. VCF to mpileup\n\t3. true_hap to ref & germline\n\t4. summarize '
            'mrbayes .lstat files\n\t5. plot summarized mrbayes output\n***')
    parser.add_argument('input', nargs='+', type=str,
        help='Absolute or relative path(s) to input file(s)')
    parser.add_argument('-f', '--format', type=str,  default='nxs',
        choices=['nxs', 'mpileup', 'ref', 'bayes', 'LRT', 'plot', 'steps',
            'pval', 'scite', 'post'],
        help='Output format to convert to. Default = "nxs".')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Path to the output directory/file. Default = <INPUT_DIR>.')
    parser.add_argument('-n', '--ngen', type=int, default=1e6,
        help='Number of MCMC steps in NEXUS MrBayes block. Default = 1e6.')
    parser.add_argument('-nc', '--no_cells', type=int, default=-1,
        help='Number of cells. Required for LRT.')
    parser.add_argument('-ss', '--stepping_stone', action='store_true',
        help='Use stepping stone sampling instead of MCMC.')
    parser.add_argument('-t', '--use_tree', action='store_true',
        help='Add the true cellcoal tree to the nexus file.')
    parser.add_argument('-lt', '--learn_tree', action='store_true',
        help='Learn the tree under no constraints model.')
    parser.add_argument('-fg', '--full_GT', action='store_true',
        help='Use the full genotype instead of only SNPs for inference.')
    parser.add_argument('-dp', '--minDP', type=int, default=1,
        help='Min. reads to include a locus (else missing). Default = 1.')
    parser.add_argument('-gq', '--minGQ', type=int, default=0,
        help='Min. Genotype Quality to include a locus (else missing). '
            'Default = 1.')
    parser.add_argument('-fs', '--filter_singletons', action='store_true',
        help='If set, singleton SNPs are filtered out.')
    parser.add_argument('-s', '--steps', nargs='+', type=int,
        help='Adjust the number of mcmc/ss steps for all runs given nxs dir.')
    parser.add_argument('-e', '--exe', type=str, default='',
        help='Path to exe to run in subprocess. Default = None.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    if not args.output:
        args.output = os.path.dirname(args.input[0])

    if args.format == 'nxs':
        out_files = [os.path.join(args.output, 'nxs.{}.{}'.format(args.ngen, i))
            for i in ['clock', 'noClock']]
        vcf_to_nex(args.input[0], out_files, args.ngen,
            ss_flag=args.stepping_stone, tree=args.use_tree,
            learn_tree=args.learn_tree, full_GT=args.full_GT)
          
    elif args.format == 'mpileup':
        out_file = os.path.join(args.output, '{}.mpileup'.format(args.input[0]))
        vcf_to_pileup(args.input[0], out_file)
    elif args.format == 'plot':
        if args.output == os.path.dirname(args.input[0]):
            out_file = ''
        else:
            out_file = args.output
        generate_mrbayes_plots(args.input[0], out_file)
    elif args.format == 'pval':
        if args.output == os.path.dirname(args.input[0]):
            out_file = ''
        else:
            out_file = args.output
        generate_pval_plot(args.input[0], out_file)
    elif args.format == 'ref':
        run_id = os.path.basename(args.input[0]).split('.')[-1]
        out_file = os.path.join(args.output, 'true_vcf.{}'.format(run_id))
        haplotypes_to_vcf(args.input[0], out_file)
    elif args.format == 'steps':
        convert_steps(args.input, args.steps)
    elif args.format == 'LRT':
        if args.no_cells < 1:
            raise IOError('For LRT, the number of cells must be specified!')
        get_LRT(args.input[0], args.output, args.no_cells + 1)
    elif args.format == 'bayes':
        if args.output == os.path.dirname(args.input[0]):
            out_file = os.path.join(args.output, 'clock_test_summary.tsv')
        else:
            out_file = args.output
        get_Bayes_factor(args.input, out_file, args.stepping_stone)
    elif args.format == 'scite':
        run_scite_subprocess(args.exe, args.steps[0], args.input[0])
    elif args.format == 'post':
        postprocess_vcf(args.input[0], args.output, minDP=args.minDP,
            minGQ=args.minGQ, s_filter=args.filter_singletons)
    else:
        raise IOError('Unknown format type: {}'.format(args.format))