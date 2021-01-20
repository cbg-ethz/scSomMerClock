#!/usr/bin/env python3

import os
import re
import gzip
import math
import argparse
from statistics import mean, stdev


BASES = ['A', 'C', 'G', 'T']
NEXUS_TEMPLATE = """#NEXUS

begin data;
    dimensions ntax={sample_no} nchar={rec_no};
    format datatype=dna missing=? gap=-;
    matrix
{matrix}
    ;
end;

begin mrbayes;
    set autoclose=yes nowarnings=yes;
    lset nst=6 rates=invgamma;
    prset brlenspr={brlen_prior};
    {alg} ngen={ngen};
    sump outputname={out_file};
end;
"""


def get_out_dir(config):
    cc_brv = config['cellcoal']['model'].get('brach_rate_var', None)
    if cc_brv:
        model = 'noClock{}'.format(cc_brv)
    else:
        model = 'clock'

    if not config['cellcoal']['scWGA'].get('simulate', False):
        sim_scWGA = 'noScWGA'
    else:
        sim_scWGA = 'scWGA'

    if not config['cellcoal']['NGS'].get('simulate', False):
        sim_NGS = 'noNGS'
    else:
        sim_NGS = 'NGS'

    mb_ngen = '-'.join(['{:.0E}'.format(i) for i in config['mrbayes']['ngen']])
    if config['mrbayes'].get('ss', False):
        mb_sampling = 'ss'
    else:
        mb_sampling = 'mcmc'

    if not config.get('static', {}).get('out_dir', False):
        out_dir = os.path.dirname(os.path.abspath(__file__))
    else:
        out_dir = config['static']['out_dir']

    return os.path.join(out_dir, 'results_{}_{}_{}_{}{}' \
        .format(model, sim_scWGA, sim_NGS, mb_sampling, mb_ngen))


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

    if config['cellcoal']['model'].get('no_muts', None):
        templ = re.sub('{no_muts}',
            'j{}'.format(config['cellcoal']['model']['no_muts']), templ)
    else:
        templ = re.sub('{no_muts}', '', templ)

    if config['cellcoal']['model'].get('branch_rate_var', None):
        templ = re.sub('{no_muts}',
            'i{}'.format(config['cellcoal']['model']['branch_rate_var']), templ)
    else:
        templ = re.sub('{branch_rate_var}', '', templ)

    if config['cellcoal']['scWGA'].get('simulate', False):
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

    if config['cellcoal']['NGS'].get('simulate', False):
        templ = re.sub('{seq_cov}',
            'C{}'.format(config['cellcoal']['NGS']['seq_cov']), templ)
        templ = re.sub('{seq_overdis}',
            'V{}'.format(config['cellcoal']['NGS']['seq_overdis']), templ)
        templ = re.sub('{seq_error}',
            'E{}'.format(config['cellcoal']['NGS']['seq_error']), templ)
    else:
        templ = re.sub('{seq_cov}', '', templ)
        templ = re.sub('{seq_overdis}', '', templ)
        templ = re.sub('{seq_error}', '', templ)

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


def vcf_to_pileup(vcf_file, out_pileup, out_samples=''):
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


def vcf_to_nex(vcf_file, out_files, ngen, ss_flag, minDP=10):
    if vcf_file.endswith('gz'):
        file_stream = gzip.open(vcf_file, 'rb')
    else:
        file_stream = open(vcf_file, 'r')

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

            ref = line_cols[3]
            alts = line_cols[4].split(',')
            DP_col = line_cols[8].split(':').index('DP')

            for s_i, s_rec in enumerate(line_cols[9:]):
                try:
                    gt = s_rec[:s_rec.index(':')]
                # Missing in Monovar output format
                except ValueError:
                    samples[s_i] += '?'
                    continue

                if int(s_rec.split(':')[DP_col]) < minDP:
                    samples[s_i] += '?'
                    continue

                if len(gt) < 3:
                    true_all1, true_all2 = s_rec[-3:].split('|')
                    if true_all1 == ref:
                        s_rec_ref = '0'
                    else:
                        s_rec_ref = str(alts.index(true_all1) + 1)

                    if true_all2 == ref:
                        s_rec_alt = '0'
                    else:
                        s_rec_alt = str(alts.index(true_all2) + 1)
                    
                else:
                    s_rec_ref, s_rec_alt = re.split('[/\|]', gt)[:2]

                # Missing
                if s_rec_ref == '.':
                    samples[s_i] += '?'
                # Wildtype
                elif s_rec_ref == '0' and s_rec_alt == '0':
                    samples[s_i] += ref
                # Heterozygous/Homozygous
                else:
                    samples[s_i] += alts[max(int(s_rec_ref), int(s_rec_alt)) - 1]
    
    mat_str = ''
    for sample_idx, genotypes in samples.items():
        mat_str += '{}    {}\n'.format(sample_names[sample_idx],  genotypes)
    rec_no = len(samples[0])

    if ss_flag:
        alg = 'ss'
    else:
        alg = 'mcmc'

    for out_file in out_files:
        model = out_file.split('.')[-1]
        if model == 'clock':
            brlen_prior = 'clock:uniform'
        else:
            brlen_prior = 'unconstrained:exp(10.0)'

        nex_str = NEXUS_TEMPLATE.format(sample_no=len(sample_names),
            rec_no=rec_no, sample_labels=' '.join(sample_names),
            matrix=mat_str.strip('\n'), out_file=os.path.basename(out_file),
            alg=alg, ngen=ngen, brlen_prior=brlen_prior)

        with open(out_file, 'w') as f_out:
            f_out.write(nex_str)
    

def format_record(format_str, s_rec):
    s_rec_vals = s_rec.split(':')
    return {j: s_rec_vals[i] for i,j in enumerate(format_str.split(':'))}
    

def parse_with_pysam(vcf_in):
    from pysam import VariantFile

    vcf = VariantFile(vcf_in)
    samples = {i: '' for i in vcf.header.samples}
    
    for rec_idx, rec in enumerate(vcf.fetch()):
        for sample_id, sample in rec.samples.iteritems():
            s_ref, s_alt = sample['GT']
            if s_ref == None:
                samples[sample_id] += '?'
            elif s_ref == 0 and s_alt == 0:
                samples[sample_id] += rec.ref
            else:
                samples[sample_id] += rec.alts[max(s_ref, s_alt) - 1]

    mat_str = ''
    sample_names = []
    for i, (sample_name, genotypes) in enumerate(samples.items()):
        sample_names.append(sample_name)
        if i == 0:
            rec_no = len(genotypes)
        mat_str += '{}    {}\n'.format(sample_name,  genotypes)

    return mat_str, rec_no, sample_names


def snv_gen_to_new(in_file):
    import numpy as np
    
    with open(in_file, 'r') as f_in:
        data_raw = f_in.read().strip().split('\n')

    sample_no = int(data_raw[0].split(' ')[0])
    rec_no = int(data_raw[0].split(' ')[1])

    df = np.empty((sample_no, rec_no), dtype='S2')
    cells = []
    for i, line in enumerate(data_raw[2:]):
        cell_name, cell_snvs = line.split('  ')
        cells.append(cell_name)
        df[i] = cell_snvs.split(' ')
    import pdb; pdb.set_trace()


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


def get_Bayes_factor(in_files, out_file):
    scores = {}
    runtimes = {}
    for in_file in in_files:
        _, run, steps, model, _ = os.path.basename(in_file).split('.')
        steps = int(steps)
        with open(in_file, 'r') as f_score:
            score_raw = f_score.read().strip().split('\n')
        score = float(score_raw[-1].split('\t')[2])
        
        log_tail = tail(open(in_file.replace('.lstat', '.log'), 'r'), 1000)
        runtime_str = re.search(
            'Analysis completed in (\d+) hours (\d+) mins (\d+) seconds',
            log_tail)
        runtime = int(runtime_str[3]) \
            + 60 * (int(runtime_str[2]) + 60 * int(runtime_str[1]))

        if not steps in scores:
            scores[steps] = {run: {'clock': -1, 'noClock': -1}}
            runtimes[steps] = []
        if not run in scores[steps]:
            scores[steps][run] = {'clock': -1, 'noClock': -1}
        scores[steps][run][model] = score
        runtimes[steps].append(runtime)

    out_str = ''
    summary_str = ''
    for step, step_info in sorted(scores.items()):
        h0_all = []
        h1_all = []
        logB_01_all = []
        evidence_all = []
        for run, run_info in step_info.items():
            h0 = run_info['clock']
            h0_all.append(h0)
            h1 = run_info['noClock']
            h1_all.append(h1)
            diff = min(max(-99, h1 - h0), 99)

            logB_01 = 2 * diff
            logB_01_all.append(logB_01)
            if logB_01 < 2:
                evidence = 'None'
                evidence_all.append(0)
            elif logB_01 < 6:
                evidence = 'Positive'
                evidence_all.append(1)
            elif logB_01 < 10:
                evidence = 'Strong'
                evidence_all.append(1)
            else:
                evidence = 'Very Strong'
                evidence_all.append(1)

            out_str += f'{step}\t{run}\t{h0:.2f}\t{h1:.2f}\t{logB_01:.2f}\t' \
                f'{math.exp(diff):.0f}\t{evidence}\n'

        if len(step_info) < 2:
            continue

        summary_str += f'{step:.0E}\t{mean(runtimes[step]):.2f}\t' \
            f'{stdev(runtimes[step])}\t{mean(h0_all):.2f}\t{stdev(h0_all):.2f}\t' \
            f'{mean(h1_all):.2f}\t{stdev(h1_all):.2f}\t{mean(logB_01_all):.2f}\t' \
            f'{stdev(logB_01_all):.2f}\t{sum(evidence_all)}/{len(evidence_all)}\n'

    with open(out_file, 'w') as f_out:
        f_out.write('steps\trun\tH_0:clock\tH_1:noClock\t2log_e(B_01)\tB_01\t'
            'Evidence\n')
        f_out.write(out_str.strip())

    if summary_str:
        with open(out_file.replace('.tsv', '.short.tsv'), 'w') as f_out:
            f_out.write('steps\tAvg. runtime [secs]\tStd. runtime [secs]\t'
                'Avg. H_0:clock\tStd. H_0:clock\tAvg. H_1:noClock\t'
                'Std. H_1:noClock\tAvg. 2log_e(B_01)\tStd. 2log_e(B_01)\tEvidence\n')
            f_out.write(summary_str.strip())



def generate_mrbayes_plots(in_file, out_file):
    import numpy as np
    import pandas as pd
    from scipy.stats import linregress
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from matplotlib.ticker import FormatStrFormatter

    TICK_FONTSIZE = 12
    LABEL_FONTSIZE = 16

    def get_ratio(ev_str):
        success, total = ev_str.split('/')
        return float(success) / float(total)

    df = pd.read_csv(in_file, sep='\t')
    df['ratio'] = df['Evidence'].apply(lambda x: get_ratio(x))
    df['runtime'] = df['Avg. runtime [secs]'].apply(lambda x: x / 60 / 60)

    fig = plt.figure(figsize=(10, 12))
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
        yerr= df['Std. 2log_e(B_01)'], color='#E31A1C', capsize=4, ls='--',
        marker='x')
    ax1.set_ylabel(r'$2 log_e(B_{01})$', fontsize=LABEL_FONTSIZE)
    ax1.set_xlabel('MCMC steps', fontsize=LABEL_FONTSIZE)
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0E'))

    reg_line = linregress(df['steps'], df['ratio'])
    reg_x_alpha = (0.05 - reg_line.intercept) / reg_line.slope
    reg_x = np.linspace(df['steps'].min(), math.ceil(reg_x_alpha / 1e7) * 1e7, 20)
    reg_y = reg_line.intercept + reg_line.slope * reg_x

    ax2.plot(reg_x, reg_y, color='#6A3D9A', ls='--', marker='x',
        label=f'{reg_line.intercept:.2f} + {reg_line.slope:.2E} * x    ($R^2$={reg_line.rvalue:.2f})')
    ax2.plot(df['steps'], df['ratio'], color='#FF7F00', ls='', marker='x',
        label='real')
    ax2.axhline(0.05, ls=':')
    ax2.axvline(reg_x_alpha, ls=':')
    ax2.text(reg_x_alpha, 0.05, f'({reg_x_alpha:.2E}, 0.05)', va='top',
        ha='right', fontsize=TICK_FONTSIZE)
    ax2.set_ylabel(r'$H_0$ rejected [%]', fontsize=LABEL_FONTSIZE)
    ax2.set_xlabel('MCMC steps', fontsize=LABEL_FONTSIZE)
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0E'))
    ax2.set_ylim(-.05, 1.05)
    ax2.legend(fontsize=TICK_FONTSIZE)

    reg_line2 = linregress(df['steps'], df['runtime'])
    def mcmc_to_runtime(x):
        return reg_line2.intercept + reg_line2.slope * x
    
    secax = ax2.secondary_xaxis('top', functions=(mcmc_to_runtime, mcmc_to_runtime))
    secax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    secax.set_xlabel('Runtime [h]')   

    fig.suptitle('Simulations: no WGS, no NGS', fontsize=LABEL_FONTSIZE * 1.5)

    plt.show()

# REF: C
# ALT: A,G,T

# GT: 0|0
# DP:18
# RC: 0,18,0,0
# G10: -44.58,-1.08,-44.46,-44.46,-0.08,-1.08,-1.08,-44.58,-44.46,-44.58
# G10N: -44.51,-1.00,-44.38,-44.38,0.00,-1.00,-1.00,-44.51,-44.38,-44.51
# GL: -0.08,-1.08,-44.58,-1.08,-44.46,-44.58,-1.08,-44.46,-44.46,-44.58
# GLN: 0.00,-1.00,-44.51,-1.00,-44.38,-44.51,-1.00,-44.38,-44.38,-44.51
# PL: -1,-11,-446,-11,-445,-446,-11,-445,-445,-446
# PLN: 0,-10,-445,-10,-444,-445,-10,-444,-444,-445
# ML: C/C
# NG: C|C
# DG: ./.
# TG: C|C


def calc_GT_likelihood(reads, eps=None, delta=None, gamma=None):
    # Line 4647 in cellcoal
    import numpy as np
    # CellCoal order: AA AC AG AT CC CG CT GG GT TT
    ll_GT = {}
    for ia1, A1 in enumerate(['A', 'C', 'G', 'T']):
        for ia2, A2 in [(0, 'A'), (1, 'C'), (2, 'G'), (3, 'T')][ia1:]:
            ll = 0
            ll2 = 0
            for ib, read in enumerate(reads):
                if gamma:
                    p_bA1 = four_temp_ampl(ib == ia1, eps, gamma) 
                    p_bA2 = four_temp_ampl(ib == ia2, eps, gamma)
                else:
                    p_bA1 = GATK(ib == ia1, eps) 
                    p_bA2 = GATK(ib == ia2, eps)

                if ia1 == ia2:
                    ll += read * math.log10(p_bA1)
                    ll2 += read * math.log10(p_bA1)
                    continue

                if delta:
                    ll += read * math.log10(
                        (1 - delta) * (0.5 * p_bA1 + 0.5 * p_bA2) +
                        delta/2 * p_bA1 + delta/2 * p_bA2)

                    t1 = read * np.log10(0.5 * p_bA1 + 0.5 * p_bA2)
                    t2 = read * np.log10(p_bA1)
                    t3 = read * np.log10(p_bA2)
                    t_max = np.max([t1, t2, t3])
                    num1 = (1 - delta) * math.pow(10, t1 - t_max)
                    num2 = delta/2 * math.pow(10, t2 - t_max)
                    num3 = delta/2 * math.pow(10, t3 - t_max)
                    gl = num1 + num2 + num3
                    ll2 += (t_max + math.log10(gl))
                else:
                    ll += read * math.log10(0.5 * p_bA1 + 0.5 * p_bA2)
                # print(ll, ll2)
                # import pdb; pdb.set_trace()
            ll_GT[A1 + A2] = round(ll, 2)
            print(f'{A1}{A2}:\t{ll:02.4f}\t{ll2:.4f}')
    print(ll_GT)
    import pdb; pdb.set_trace()


def GATK(is_same, eps):
    if is_same:
        return 1 - eps
    else:
        return eps / 3


# Line 4675
def four_temp_ampl(is_same, eps, gamma):
    if is_same:
        return (1 - gamma) * (1 - eps) + gamma * eps / 3
    else:
        return (1 - gamma) * eps / 3 + (gamma / 3) * (1 - eps / 3)


# def two_temp_ampl(is_same, eps, gamma):
#     3 * ((1 - gamma) * eps + gamma * eps)

#     if is_same:
#         return (1 - gamma) * (1 - eps) + gamma * eps / 3
#     else:
#         return (1 - gamma) * eps / 3 + gamma * (1 - eps / 3) / 3



def parse_args():
    parser = argparse.ArgumentParser(prog='utils.py',
        usage='python3 utils.py <DATA> [options]',
        description='*** Convert VCF to NEXUS or mpileup file ***')
    parser.add_argument('input', nargs='+', type=str, 
        help='Absolute or relative path(s) to input VCF file(s)')
    parser.add_argument('-f', '--format', type=str, 
        choices=['nxs', 'mpileup', 'bayes', 'plot'], default='nxs',
        help='Output format to convert to. Default = "nxs".')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Path to the output directory. Default = <INPUT_DIR>.')
    parser.add_argument('-n', '--ngen', type=int, default=1e6,
        help='Number of MCMC steps in NEXUS MrBayes block. Default = 1e6.')
    parser.add_argument('-ss', '--stepping_stone', action='store_true',
        help='Use stepping stone sampling instead of MCMC.')
    parser.add_argument('-dp', '--minDP', type=int, default=10,
        help='Minimum reads to include a locus (else: missing). Default = 10.')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    # calc_GT_likelihood([0, 18, 0, 0], eps=1e-02, delta=0.2, gamma=1e-02)

    args = parse_args()
    if not args.output:
        args.output = os.path.dirname(args.input[0])

    if args.format == 'nxs':
        out_files = [os.path.join(args.output, 'nxs.{}.{}'.format(args.ngen, i))
            for i in ['clock', 'noClock']]
        vcf_to_nex(args.input[0], out_files, args.ngen, args.stepping_stone,
            args.minDP)
    elif args.format == 'mpileup':
        out_file = os.path.join(args.output, '{}.mpileup'.format(args.input[0]))
        vcf_to_pileup(args.input[0], out_file)
    elif args.format == 'plot':
        out_file = os.path.join(args.output, '{}.summary.pdf'.format(args.input[0]))
        generate_mrbayes_plots(args.input[0], out_file)
    else:
        out_file = os.path.join(args.output, 'clock_test_summary.tsv')
        get_Bayes_factor(args.input, out_file)