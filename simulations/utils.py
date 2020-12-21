#!/usr/bin/env python3

import os
import re
import math
import argparse


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
    ss ngen={ngen};
    sump outputname={out_file};
end;
"""


def load_config(configfile):
    cc_params = {}
    with open(configfile, 'r') as f:
        for line in f.read().strip().split('\n'):
            pair = line.strip().split('] ')
            if len(pair) < 2 or pair[1].startswith('['):
                continue
            cc_params[pair[0].lstrip('[')] = pair[1][1:]
    return cc_params


def vcf_to_nex(vcf_file, out_files, ngen):
    header = ''
    if vcf_file.endswith('gz'):
        import gzip
        file_stream = gzip.open(vcf_file, 'rb')
        byte = True
    else:
        file_stream = open(vcf_file, 'r')
        byte = False

    with file_stream as f_in:
        for line in f_in:
            if byte:
                line = line.decode("utf-8")
            # VCF header
            if line.startswith('##'):
                header += line
            # Column headers
            elif line.startswith('#'):
                sample_names = line.strip().split('\t')[9:]
                samples = {int(i): '' for i in range(len(sample_names))}
            # VCF records
            else:
                line_cols = line.strip().split('\t')
                ref = line_cols[3]
                alts = line_cols[4].split(',')
                for s_i, s_rec in enumerate(line_cols[9:]):
                    gt = s_rec[:s_rec.index(':')]
                    if len(gt) < 3:
                        s_rec_ref = '0'
                        s_rec_alt = re.split('[/\|]', gt)
                    else:
                        s_rec_ref, s_rec_alt = re.split('[/\|]', gt)[:2]

                    if s_rec_ref == '.':
                        samples[s_i] += '?'
                    elif s_rec_ref == '0' and s_rec_alt == '0':
                        samples[s_i] += ref
                    else:
                        samples[s_i] += alts[max(int(s_rec_ref), int(s_rec_alt)) - 1]
    
    mat_str = ''
    for sample_idx, genotypes in samples.items():
        mat_str += '{}    {}\n'.format(sample_names[sample_idx],  genotypes)
    rec_no = len(samples[0])

    for out_file in out_files:
        model = out_file.split('.')[-1]
        if model == 'clock':
            brlen_prior = 'clock:uniform'
        else:
            brlen_prior = 'unconstrained:exp(10.0)'

        nex_str = NEXUS_TEMPLATE.format(sample_no=len(sample_names),
            rec_no=rec_no, sample_labels=' '.join(sample_names),
            matrix=mat_str.strip('\n'), out_file=out_file, ngen=ngen,
            brlen_prior=brlen_prior)

        with open(out_file, 'w') as f_out:
            f_out.write(nex_str)
    

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
        diff = max(-99, h1 - h0)

        logB_12 = 2 * diff
        if logB_12 < 2:
            evidence = 'None'
        elif logB_12 < 6:
            evidence = 'Positive'
        elif logB_12 < 10:
            evidence = 'Strong'
        else:
            evidence = 'Very Strong'

        out_str += f'{run}\t{h1}\t{h2}\t{logB_12}\t{math.exp(diff)}\t{evidence}\n'

    with open(out_file, 'w') as f_out:
        f_out.write('run\tH_0:clock\tH_1:noClock\t2log_e(B_01)\tB_01\tEvidence\n')
        f_out.write(out_str.strip())
    


def parse_args():
    parser = argparse.ArgumentParser(prog='utils.py',
        usage='python3 utils.py <DATA> [options]',
        description='*** Convert VCF to NEXUS file ***')
    parser.add_argument('input', type=str,
        help='Absolute or relative path(s) to input VCF file')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Path to the output directory. Default = <INPUT_DIR>.')
    parser.add_argument('-n', '--ngen', type=int, default=1e6,
        help='Number of MCMC steps in NEXUS MrBayes block. Default = 1e6.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    vcf_to_nex(args.input, args.output, args.ngen)
    
