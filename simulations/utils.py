#!/usr/bin/env python3

import os
import re
import gzip
import math
import argparse
from pprint import pprint as pp
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


VCF_TEMPLATE = """##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="True genotype">
##contig=<ID=1,length={contig_length}>
##source=germline variants from CellCoal simulation
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT{cells}
{rows}
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
        out_dir = os.path.dirname(os.path.relpath(__file__))
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
    templ = re.sub('{seq_cov}',
        str(config['cellcoal'].get('NGS', {}).get('seq_cov', 20)), templ)

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
            FORMAT_col = line_cols[8].split(':')
            DP_col = FORMAT_col.index('DP')


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
                    
                # s_recs = s_rec.split(':')
                # ml_gt = ''.join(sorted(s_recs[-4].split('/')))
                # tr_gt = ''.join(sorted(s_recs[-1].split('|')))
                # if ml_gt != tr_gt:
                #     rec_dict = {j: s_recs[i] for i,j in enumerate(FORMAT_col)}
                #     print('read counts: {}'.format(rec_dict['RC']))
                #     print('GT ll norm: {}'.format(rec_dict['G10N']))
                #     print(f'predicted: {ml_gt}\t true: {tr_gt}')
                #     read_counts = [int(i) for i in rec_dict['RC'].split(',')]
                #     calc_G10N_likelihood(read_counts, eps=1e-02, delta=0.2, gamma=1e-02)
                #     import pdb; pdb.set_trace()


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
            f'{stdev(runtimes[step]):.2f}\t{mean(h0_all):.2f}\t{stdev(h0_all):.2f}\t' \
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



def generate_mrbayes_plots(in_file, out_file, regress=False):
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


def calc_G10N_likelihood(reads, eps=None, delta=None, gamma=None):
    # Line 4647 in cellcoal
    import numpy as np
    # CellCoal order: AA AC AG AT CC CG CT GG GT TT
    ll_GT = {}
    for ia1, A1 in enumerate(['A', 'C', 'G', 'T']):
        for ia2, A2 in [(0, 'A'), (1, 'C'), (2, 'G'), (3, 'T')][ia1:]:
            ll = 0

            g0 = 0
            g1 = 0
            g2 = 0
            for ib, read in enumerate(reads):
                if gamma:
                    p_bA1 = four_temp_ampl(ib == ia1, eps, gamma) 
                    p_bA2 = four_temp_ampl(ib == ia2, eps, gamma)
                else:
                    p_bA1 = GATK(ib == ia1, eps) 
                    p_bA2 = GATK(ib == ia2, eps)
                if delta:
                    # Lines 4688f
                    g0 += read * np.log10(0.5 * p_bA1 + 0.5 * p_bA2)
                    g1 += read * np.log10(p_bA1)
                    g2 += read * np.log10(p_bA2)
                else:
                    ll += read * np.log10(0.5 * p_bA1 + 0.5 * p_bA2)

            if delta:
                g0 = np.log10(1 - delta) + g0
                g1 = np.log10(delta/2) + g1
                g2 = np.log10(delta/2) + g2

                max_t = np.max([g0, g1, g2])

                # t0 = (1 - delta) * g0
                # t1 = (delta/2) * g1
                # t2 = (delta/2) * g2

                t0 = g0
                t1 = g1
                t2 = g2

                ll = max_t + np.log10(np.power(10, t0 - max_t) \
                    + np.power(10, t1 - max_t) + np.power(10, t2 - max_t))
                
                ll_GT[A1 + A2] = round(ll, 2)
    import pdb; pdb.set_trace()
    x = np.array(list(ll_GT.values()))
    ll_norm = _normalize_log(np.array(list(ll_GT.values())))
    

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
    parser.add_argument('-f', '--format', type=str, 
        choices=['nxs', 'mpileup', 'ref', 'bayes', 'plot', 'steps'], default='nxs',
        help='Output format to convert to. Default = "nxs".')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Path to the output directory/file. Default = <INPUT_DIR>.')
    parser.add_argument('-n', '--ngen', type=int, default=1e6,
        help='Number of MCMC steps in NEXUS MrBayes block. Default = 1e6.')
    parser.add_argument('-ss', '--stepping_stone', action='store_true',
        help='Use stepping stone sampling instead of MCMC.')
    parser.add_argument('-dp', '--minDP', type=int, default=10,
        help='Minimum reads to include a locus (else: missing). Default = 10.')
    parser.add_argument('-s', '--steps', nargs='+', type=int,
        help='Adjust the number of mcmc/ss steps for all runs given nxs dir.')

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
        if args.output == os.path.dirname(args.input[0]):
            out_file = ''
        else:
            out_file = args.output
        generate_mrbayes_plots(args.input[0], out_file)
    elif args.format == 'ref':
        run_id = os.path.basename(args.input[0]).split('.')[-1]
        out_file = os.path.join(args.output, 'true_vcf.{}'.format(run_id))
        haplotypes_to_vcf(args.input[0], out_file)
    elif args.format == 'steps':
        convert_steps(args.input, args.steps)
    else:
        out_file = os.path.join(args.output, 'clock_test_summary.tsv')
        get_Bayes_factor(args.input, out_file)