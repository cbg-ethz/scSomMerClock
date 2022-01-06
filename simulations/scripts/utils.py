#!/usr/bin/env python3

import os
import re
import gzip
import tempfile
import subprocess


def get_out_dir(config, bulk=False):
    if config['cellcoal']['model'].get('branch_rate_var', 0):
        model = 'clock{}'.format(config['cellcoal']['model']['branch_rate_var'])
    else:
        model = 'clock0'

    scWGA = config['cellcoal']['scWGA']
    if scWGA.get('errors', False):
        if isinstance(scWGA["ADO_rate"], float):
            sim_scWGA  = f'WGA{scWGA["ADO_rate"]}'
        else:
            if len(scWGA['ADO_rate']) == 1:
                sim_scWGA  = f'WGA{scWGA["ADO_rate"][0]}'
            else:
                max_var = scWGA["ADO_rate"][0] * (1 - scWGA["ADO_rate"][0])
                max_var = int(max_var * 100) / 100 - 0.001
                scWGA["ADO_rate"][1] = min(scWGA["ADO_rate"][1], max_var)
                sim_scWGA = f'WGA{scWGA["ADO_rate"][0]},{scWGA["ADO_rate"][1]}'
        sim_scWGA += f'-{scWGA["doublet_rate"][0]}-{scWGA["ampl_error"][0]}'
    else:
        sim_scWGA = 'WGA0-0-0'

    NGS = config['cellcoal']['NGS']
    sim_NGS = f'NGS{NGS["seq_cov"]}-'
    if NGS.get('errors', False):
        sim_NGS += f'{NGS["seq_overdis"]}-{NGS["seq_error"]}'
    else:
        sim_NGS += '0-0'

    if not config.get('static', {}).get('out_dir', False):
        out_dir = os.path.dirname(os.path.dirname(os.path.relpath(__file__)))
    else:
        out_dir = config['static']['out_dir']

    if bulk:
        depth = int(config['cellcoal']['model']['no_cells'] \
            * config['cellcoal']['NGS']['seq_cov'])
        sim_dir = os.path.join(out_dir, f'res_{model}_bulk{depth}x')
    else:
        sim_dir = os.path.join(out_dir, f'res_{model}_{sim_scWGA}_{sim_NGS}')

    filter_dir =  'minDP{}-minGQ{}'.format(
        config.get('SNP_filter', {}).get('depth', 1),
        config.get('SNP_filter', {}).get('quality', 0))
    if config.get('SNP_filter', {}).get('singletons', False):
        filter_dir += '-noSingletons'

    return sim_dir, filter_dir


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


def get_sample_dict_from_vcf(vcf_file, GT=False, include='', exclude=''):
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

    exclude_i = []
    with file_stream as f_in:
        for line in f_in:
            try:
                line = line.decode()
            except AttributeError:
                pass
            # Skip VCF header lines
            if line.startswith('#'):
                # Safe column headers
                if line.startswith('#CHROM'):
                    sample_names = line.strip().split('\t')[9:]
                    samples = {int(i): '' for i in range(len(sample_names))}
                    if exclude != '':
                        print('Exluded:')
                        for s_i, sample in enumerate(sample_names):
                            if re.fullmatch(exclude, sample):
                                exclude_i.append(s_i)
                                print('\t{}'.format(sample))
                        if len(exclude_i) == 0:
                            print('\nWARNING: no samples with pattern {} in vcf!\n' \
                                .format(exclude))
                    if include != '':
                        print('Include:')
                        for s_i, sample in enumerate(sample_names):
                            if not re.fullmatch(include, sample):
                                exclude_i.append(s_i)
                            else:
                                if not re.fullmatch(exclude, sample):
                                    print('\t{}'.format(sample))
                        if len(exclude_i) == 0:
                            print('\nWARNING: no samples with pattern {} in vcf!\n' \
                                .format(include))

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

    # Remove excluded samples
    for i in sorted(exclude_i, reverse=True):
        samples.pop(i)
    # Sanity check: remove sample without name
    try:
        samples.pop(sample_names.index(''))
    except (ValueError, KeyError):
        pass


    return samples, [i.replace('-', '_') for i in sample_names]


def change_newick_tree_root(in_file, paup_exe, root=True, outg='healthycell',
        sample_names=[], br_length=False):
    paup_file = tempfile.NamedTemporaryFile(delete=False)
    out_file = tempfile.NamedTemporaryFile(delete=False)
    temp_tree_file = tempfile.NamedTemporaryFile(delete=False)

    with open(in_file, 'r') as f_tree:
        tree = f_tree.read().strip()

    if 'scite' in in_file:
        # Add missing semicolon
        sem_count = tree.count(';')
        if sem_count == 0:
            tree += ';'
        elif sem_count > 1:
            tree = tree.split(';')[0].strip() + ';'

        nodes = [int(i) for i in re.findall('(?<=[\(\),])\d+(?=[,\)\(;])', tree)]
        cells = len(sample_names)
        for i, s_i in enumerate(sorted(nodes)):
            if i < cells:
                pat = '(?<=[\(\),]){}(?=[,\)\(;)])'.format(s_i)
                repl = sample_names[i]
                if br_length:
                    repl = '{}:0.1'.format(sample_names[i])
                else:
                    repl = str(sample_names[i])
            # elif i == cells:
            #     import pdb; pdb.set_trace()
            #     pat = '(?<=[\(\),]){},?(?=[\)\(;)])'.format(s_i)
            #     if br_length:
            #         repl = ':0.1'
            #     else:
            #         repl = ''
            else:
                pat = '(?<=[\(\),]){}(?=[,\)\(;)])'.format(s_i)
                if br_length:
                    repl = f'{i+1}:0.1'
                else:
                    repl = ''
            tree = re.sub(pat, repl, tree)
    else:
        if not br_length:
            tree = re.sub(':0.\d+', '', tree)

        if 'cellphy' in in_file:
            tree =re.sub('\[\d+\]', '', tree)
        else:
            tree = tree.replace('cell', 'tumcell') \
                .replace('outgtumcell', 'healthycell')
    
    temp_tree_file.write(str.encode(tree))
    temp_tree_file.close()

    if outg == '' or not outg in tree:
        outg_cmd = ''
    else:
        outg_cmd = 'outgroup {};\n'.format(outg)

    if root:
        root_cmd = 'DerootTrees;\nRootTrees rootMethod=outgroup outroot=monophyl'
        root = 'yes'
    else:
        root_cmd = 'DerootTrees;\n{}'.format(outg_cmd)
        root = 'no'
    if br_length:
        save_brLens = 'user'
    else:
        save_brLens = 'no'
        
    paup_cmd = 'getTrees file={i} allBlocks=yes;\n' \
        '{g}' \
        '{c};\n' \
        'saveTrees format=Newick root={r} brLens={b} file={o} supportValues=nodeLabels;\n' \
        'quit;'.format(i=temp_tree_file.name, g=outg_cmd, c=root_cmd, r=root,
            b=save_brLens, o=out_file.name)

    paup_file.write(str.encode(paup_cmd))
    paup_file.close()

    shell_cmd = ' '.join([paup_exe, '-n', paup_file.name ])#, '>', '/dev/null'])
    paup = subprocess.Popen(shell_cmd, shell=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    stdout, stderr = paup.communicate()
    paup.wait()

    if stderr != b'' or not 'tree saved to file ' in str(stdout):
        raise RuntimeError(f'{stdout}\n{stderr}')
    
    with open(out_file.name, 'r') as f_tree:
        tree_new = f_tree.read().strip()
    out_file.close()

    assert tree != '', 'Could not read tree from {}'.format(in_file)
    assert tree_new != '', 'Failed to root/unroot tree with PAUP'

    return tree, tree_new


if __name__ == '__main__':
    print('Nothing to do...')