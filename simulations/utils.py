#!/usr/bin/env python3

import os
import re
import gzip
import math
import argparse


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
        config['static']['out_dir'] == os.path.dirname(os.path.abspath(__file__))

    return os.path.join(config['static']['out_dir'], 'results_{}_{}_{}_{}{}' \
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


def vcf_to_nex(vcf_file, out_files, ngen, ss_flag):
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
            ref = line_cols[3]
            alts = line_cols[4].split(',')
            for s_i, s_rec in enumerate(line_cols[9:]):
                try:
                    gt = s_rec[:s_rec.index(':')]
                # Missing in Monovar output format
                except ValueError:
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
            matrix=mat_str.strip('\n'), out_file=out_file, alg=alg, ngen=ngen,
            brlen_prior=brlen_prior)

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


def get_Bayes_factor(in_files, out_file):
    scores = {}
    for in_file in in_files:
        _, run, model, _ = os.path.basename(in_file).split('.')
        with open(in_file, 'r') as f_score:
            score_raw = f_score.read().strip().split('\n')
        score = float(score_raw[-1].split('\t')[2])
        
        if not run in scores:
            scores[run] = {'clock': -1, 'noClock': -1}
        scores[run][model] = score

    out_str = ''
    for run, run_info in scores.items():
        h0 = run_info['clock']
        h1 = run_info['noClock']
        diff = min(max(-99, h1 - h0), 99)

        logB_01 = 2 * diff
        if logB_01 < 2:
            evidence = 'None'
        elif logB_01 < 6:
            evidence = 'Positive'
        elif logB_01 < 10:
            evidence = 'Strong'
        else:
            evidence = 'Very Strong'

        out_str += f'{run}\t{h0}\t{h1}\t{logB_01}\t{math.exp(diff)}\t{evidence}\n'

    with open(out_file, 'w') as f_out:
        f_out.write('run\tH_0:clock\tH_1:noClock\t2log_e(B_01)\tB_01\tEvidence\n')
        f_out.write(out_str.strip())
    


def parse_args():
    parser = argparse.ArgumentParser(prog='utils.py',
        usage='python3 utils.py <DATA> [options]',
        description='*** Convert VCF to NEXUS or mpileup file ***')
    parser.add_argument('input', type=str,
        help='Absolute or relative path(s) to input VCF file')
    parser.add_argument('-f', '--format', type=str, choices=['nxs', 'mpileup'],
        default='nxs', help='Output format to convert to. Default = "nxs".')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Path to the output directory. Default = <INPUT_DIR>.')
    parser.add_argument('-n', '--ngen', type=int, default=1e6,
        help='Number of MCMC steps in NEXUS MrBayes block. Default = 1e6.')
    parser.add_argument('-ss', '--stepping_stone', action='store_true',
        help='Use stepping stone sampling instead of MCMC.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    if not args.output:
        args.output = os.path.dirname(args.input)

    if args.format == 'nxs':
        out_files = [os.path.join(args.output, 'nxs.{}.{}'.format(args.ngen, i))
            for i in ['clock', 'noClock']]
        vcf_to_nex(args.input, out_files, args.ngen, args.stepping_stone)
    else:
        out_file = '{}.mpileup'.format(args.output)
        vcf_to_pileup(args.input, out_file)
    
