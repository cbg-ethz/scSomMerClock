#!/usr/bin/env python3

import os
import re
import gzip
import math
import argparse
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
    outgroup healthycell;
    prset brlenspr={brlen_prior};
    {fixed_tree}
    {alg} ngen={ngen};
    sump outputname={out_file};
end;

begin PAUP;
    Set autoclose=yes warnreset=no warntree=no warntsave=No;
    Set criterion=like;
    Outgroup healthycell;
    LSet nst=1 base=equal rates=gamma shape=est;

    {paup_tree}
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
    cc_brv = config['cellcoal']['model'].get('branch_rate_var', None)
    if cc_brv:
        model = 'noClock{}'.format(cc_brv)
    else:
        model = 'clock'

    if not config['cellcoal']['scWGA'].get('errors', False):
        sim_scWGA = 'noScWGA'
    else:
        sim_scWGA = 'scWGA'

    if not config['cellcoal']['NGS'].get('errors', False):
        sim_NGS = 'noNGS'
    else:
        sim_NGS = 'NGS'

    if config.get('mrbayes', {}).get('run', False):
        mb_ngen = '-'.join(['{:.0E}'.format(i) for i in config['mrbayes']['ngen']])
        if config['mrbayes'].get('ss', False):
            sampling = 'ss'
        else:
            sampling = 'mcmc'
    elif config.get('paup', {}).get('run', False):
        mb_ngen = 1
        sampling = 'ML'
    else:
        mb_ngen = -1
        sampling = ''

    if not config.get('static', {}).get('out_dir', False):
        out_dir = os.path.dirname(os.path.relpath(__file__))
    else:
        out_dir = config['static']['out_dir']

    if config.get('paup', {}).get('run', False):
        if config['paup'].get('full_GT', False):
            data_type = '_fullGT'
        else:
            data_type = '_onlySNP'
    else:
        data_type = ''

    filters = ''
    if config['cellcoal'].get('SNP_filter', {}).get('depth', False):
        filters += '_minDP{}'.format(config['cellcoal']['SNP_filter']['depth'])
    if config['cellcoal'].get('SNP_filter', {}).get('quality', False): 
        filters += '_minGQ{}'.format(config['cellcoal']['SNP_filter']['quality'])

    return os.path.join(out_dir, 'results_{}_{}_{}{}{}_{}{}' \
        .format(model, sim_scWGA, sim_NGS, data_type, filters, sampling, mb_ngen))


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

    if config.get('mrbayes', {}).get('use_tree', False) \
            or (config.get('paup', {}).get('run', False) \
                and not config.get('paup', {}).get('learn_tree', False)) \
            or config['cellcoal'].get('output', {}).get('tree', False):
        templ = re.sub('{out_tree}', '6', templ)
    else:
        templ = re.sub('{out_tree}', '', templ)

    if config['cellcoal'].get('output', {}).get('full_GT', False) \
            or config.get('paup', {}).get('full_GT', False) :
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


def vcf_to_nex(vcf_file, out_files, ngen, ss_flag=False, tree=False,
            learn_tree=False, full_GT=False, minDP=1, minGQ=1):
    if full_GT:
        samples, sample_names = get_sample_dict_from_FG(vcf_file)
    else:
        samples, sample_names = get_sample_dict_from_vcf(vcf_file, minDP, minGQ)

    mat_str = ''
    for sample_idx, genotypes in samples.items():
        mat_str += '{}    {}\n'.format(sample_names[sample_idx],  genotypes)
    rec_no = len(samples[0])

    if ss_flag:
        alg = 'ss'
    else:
        alg = 'mcmc'

    if tree:
        run_nr = os.path.basename(vcf_file).split('.')[-1]
        tree_file = os.sep.join(vcf_file.split(os.sep)[:-2] \
            + ['trees_dir', 'trees.{}'.format(run_nr)])
        with open(tree_file, 'r') as f_tree:
            tree_nexus = f_tree.read().strip()
        tree_nexus = tree_nexus.replace('cell', 'tumcell')
        tree_nexus = tree_nexus.replace('outgtumcell', 'healthycell')

        tree_name = 'treeCellCoal'
        fixed_tree = 'prset topologypr=fixed({});'.format(tree_name)
    else:
        fixed_tree = ''

    if learn_tree:
        paup_tree_raw = 'Lset clock=no;\n' \
            '    log file={out_file}.PAUP.score start=yes replace=yes;\n' \
            '    Hsearch;\n' \
            '    RootTrees rootMethod=outgroup outroot=monophyl;\n' \
            '    clockChecker tree=all lrt=yes;\n' \
            '    log stop;' \
            '    quit;'
    else:
        paup_tree_raw = 'LSet clock={clock_str};\n' \
            '    lscores 1/ clock={clock_str} scorefile={out_file}.PAUP.score append;\n'\
            '    quit;'

    for out_file in out_files:
        model = out_file.split('.')[-1]
        if model == 'clock':
            brlen_prior = 'clock:uniform'
            clock_str='yes'
            if tree:
                tree_str = 'tree {} = [&R] {}'.format(tree_name, tree_nexus)
            else:
                tree_str = ''
        else:
            brlen_prior = 'unconstrained:exp(10.0)'
            clock_str='no'
            if tree:
                tree_nexus_unrooted = re.sub('\):\d+.\d+(?=,healthycell)', '',
                    tree_nexus[1:])
                tree_str = 'tree {} = [&U] {}'\
                    .format(tree_name, tree_nexus_unrooted)
            else:
                tree_str = ''

        if learn_tree:
            paup_tree = paup_tree_raw.format(out_file=os.path.basename(out_file))
        else:
            paup_tree = paup_tree_raw \
                .format(clock_str=clock_str, out_file=os.path.basename(out_file))

        nex_str = NEXUS_TEMPLATE.format(sample_no=len(sample_names),
            rec_no=rec_no, sample_labels=' '.join(sample_names),
            matrix=mat_str.strip('\n'), tree_str=tree_str, fixed_tree=fixed_tree,
            out_file=os.path.basename(out_file), alg=alg, ngen=ngen,
            brlen_prior=brlen_prior, paup_tree=paup_tree)

        with open(out_file, 'w') as f_out:
            f_out.write(nex_str)
   

def get_sample_dict_from_vcf(vcf_file, minDP=1, minGQ=1):
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
            GQ_col = FORMAT_col.index('PLN')

            for s_i, s_rec in enumerate(line_cols[9:]):
                try:
                    gt = s_rec[:s_rec.index(':')]
                # Missing in Monovar output format
                except ValueError:
                    samples[s_i] += '?'
                    continue
                s_rec_details = s_rec.split(':')
                

                if len(FORMAT_col) == len(s_rec_details):
                    s_rec_GQ_col = GQ_col
                else:
                    s_rec_GQ_col = GQ_col - len(FORMAT_col) + len(s_rec_details)
                s_rec_GQ_all = [-float(i) for i in s_rec_details[s_rec_GQ_col] \
                    .split(',')]
                s_rec_GQ = sorted(s_rec_GQ_all)[1]

                if s_rec_GQ >= minGQ:
                    GQ_skip_flag = False
                else:
                    GQ_skip_flag = True
                    GT_max = s_rec_GQ_all.index(0)
                    # 00,01,11,02,12,22,03,13,23,33
                    if GT_max != 0 \
                            and s_rec_GQ_all[0] - max(s_rec_GQ_all[1:]) >= minGQ:
                        GQ_skip_flag = False

                if int(s_rec_details[DP_col]) < minDP or GQ_skip_flag:
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
                else:
                    samples[s_i] += alts[max(int(s_rec_ref), int(s_rec_alt)) - 1]
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

        if len(runtime) == 3:
            runtime = 60 * 60 * int(runtime[0]) + 60 * int(runtime[1]) \
                + int(runtime[2])
        else:
            runtime = 60 * int(runtime[0]) + int(runtime[1])

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



def get_LRT(in_files, out_file, cell_no, alpha=0.05):
    from scipy.stats.distributions import chi2

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
    parser.add_argument('-f', '--format', type=str,  default='nxs',
        choices=['nxs', 'mpileup', 'ref', 'bayes', 'LRT', 'plot', 'steps', 'pval'],
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
    parser.add_argument('-gq', '--minGQ', type=int, default=1,
        help='Min. Genotype Quality to include a locus (else missing). '
            'Default = 1.')
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
        vcf_to_nex(args.input[0], out_files, args.ngen,
            ss_flag=args.stepping_stone, tree=args.use_tree,
            learn_tree=args.learn_tree, full_GT=args.full_GT, minDP=args.minDP,
            minGQ=args.minGQ)
          
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
        get_LRT(args.input, args.output, args.no_cells + 1)
    else:
        if args.output == os.path.dirname(args.input[0]):
            out_file = os.path.join(args.output, 'clock_test_summary.tsv')
        else:
            out_file = args.output
        get_Bayes_factor(args.input, out_file, args.stepping_stone)