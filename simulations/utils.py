#!/usr/bin/env python3

import numpy as np

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
    ss ngen={ngen};
    sump outputname={out_file}.noClock;
    prset brlenspr=clock:uniform;
    ss ngen={ngen};
    sump outputname={out_file}.clock;
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


def vcf_to_nex(vcf_file, out_file, ngen):
    header = ''
    with open(vcf_file, 'r') as f_in:
        for line in f_in:
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
                        s_rec_alt = gt.strip('|')
                    else:
                        s_rec_ref, s_rec_alt = gt.split('|')

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

    nex_str = NEXUS_TEMPLATE.format(sample_no=len(sample_names), rec_no=rec_no,
        sample_labels=' '.join(sample_names), matrix=mat_str.strip('\n'),
        out_file=out_file, ngen=ngen)

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