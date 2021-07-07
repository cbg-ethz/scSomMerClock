#!/usr/bin/env python3

import os
import re
import tempfile
import subprocess
from io import StringIO

import numpy as np
np.seterr(all='raise')
np.set_printoptions(suppress=True)
import pandas as pd
from scipy.stats.distributions import chi2
import scipy.special
from scipy.optimize import minimize
from scipy.stats import poisson

from Bio import Phylo


LAMBDA_MIN = 1e-6
LAMBDA_MAX = np.inf


def change_newick_tree_root(in_file, paup_exe, root=True, tree_format='scite'):
    paup_file = tempfile.NamedTemporaryFile(delete=False)
    out_file = tempfile.NamedTemporaryFile(delete=False)
    temp_tree_file = tempfile.NamedTemporaryFile(delete=False)

    with open(in_file, 'r') as f_tree:
        tree = f_tree.read().strip()

    if tree_format == 'scite':
        # Add missing semicolon
        sem_count = tree.count(';')
        if sem_count == 0:
            tree += ';'
        elif sem_count > 1:
            tree = tree.split(';')[0].strip() + ';'

        nodes = [int(i) for i in re.findall('(?<=[\(\),])\d+(?=[,\)\(;])', tree)]
        sample_names = [f'cell{i:0>4d}' for i in range(max(nodes) // 2 + 1)]
        cells = len(sample_names)

        for i, s_i in enumerate(sorted(nodes)):
            if i < cells:
                pat = '(?<=[\(\),]){}(?=[,\)\(;)])'.format(s_i)
                repl = sample_names[i]
                repl = '{}:0.1'.format(sample_names[i])
            else:
                pat = '(?<=[\(\),]){}(?=[,\)\(;)])'.format(s_i)
                repl = ':0.1'

            tree = re.sub(pat, repl, tree)

    return tree, None

    # temp_tree_file.write(str.encode(tree))
    # temp_tree_file.close()


    # if outg == '' or not outg in tree:
    #     outg_cmd = ''
    # else:
    #     outg_cmd = 'outgroup {};\n'.format(outg)

    # if root:
    #     root_cmd = 'DerootTrees;\nRootTrees rootMethod=outgroup outroot=monophyl'
    #     root = 'yes'
    # else:
    #     root_cmd = 'DerootTrees;\n{}'.format(outg_cmd)
    #     root = 'no'
    # if br_length:
    #     save_brLens = 'user'
    # else:
    #     save_brLens = 'no'

    # paup_cmd = 'getTrees file={i};\n' \
    #     '{g}' \
    #     '{c};\n' \
    #     'saveTrees format=Newick root={r} brLens={b} file={o};\n' \
    #     'quit;'.format(i=temp_tree_file.name, g=outg_cmd, c=root_cmd, r=root,
    #         b=save_brLens, o=out_file.name)

    # paup_file.write(str.encode(paup_cmd))
    # paup_file.close()

    # shell_cmd = ' '.join([paup_exe, '-n', paup_file.name ])#, '>', '/dev/null'])
    # paup = subprocess.Popen(shell_cmd, shell=True, stdout=subprocess.PIPE,
    #     stderr=subprocess.PIPE)
    # stdout, stderr = paup.communicate()
    # paup.wait()

    # assert stderr == b'', str(stdout) + '\n' + str(stderr)

    # with open(out_file.name, 'r') as f_tree:
    #     tree_new = f_tree.read().strip()
    # out_file.close()

    # assert tree != '', 'Could not read tree from {}'.format(in_file)
    # assert tree_new != '', 'Failed to root/unroot tree with PAUP'

    # return tree, tree_new


def get_mut_df(vcf_file):
    if vcf_file.endswith('gz'):
        file_stream = gzip.open(vcf_file, 'rb')
    else:
        file_stream = open(vcf_file, 'r')

    muts_raw = []
    with file_stream as f_in:
        for line in f_in:
            # Skip VCF header lines
            if line.startswith('#'):
                # Safe column headers
                if line.startswith('#CHROM'):
                    sample_names = line.strip().split('\t')[9:]
                    sample_no = len(sample_names)
                    samples = [0 for i in range(sample_no)]
                continue

            elif line.strip() == '':
                continue

            # VCF records
            line_cols = line.strip().split('\t')
            line_muts = np.zeros(sample_no)
            for s_i, s_rec in enumerate(line_cols[9:]):
                try:
                    gt = s_rec[:s_rec.index(':')]
                # Missing in Monovar output format
                except ValueError:
                    line_muts[s_i] = np.nan
                    continue

                s_rec_ref = gt[0]
                s_rec_alt = gt[-1]
                if s_rec_ref == '.' or s_rec_alt == '.':
                    continue

                if s_rec_alt != '0' or s_rec_ref != '0':
                    line_muts[s_i] = 1
            muts_raw.append(line_muts)

    muts = pd.DataFrame(muts_raw, columns=sample_names)
    include = [i for i in sample_names if i]
    
    return muts[include]


def get_tree_dict(tree_file, muts, paup_exe,  exclude='', include='',
        tree_format='scite'):
    if tree_format == 'cellphy':
        _, tree_str = change_newick_tree_root(tree_file, paup_exe, root=True,
            tree_format='cellphy')
    elif tree_format == 'scite':
        tree_str, _ = change_newick_tree_root(tree_file, paup_exe, root=False,
            tree_format='scite')
    else:
        raise RuntimeError(f'Unknown tree format: {tree_format}')


    # With BioPython package
    tree = Phylo.read(StringIO(tree_str), 'newick')

    # Get max age if branch lengths display age instead of mutation number
    max_age = muts.sum(axis=0).max()

    # Add number of mutations to terminal nodes
    for leaf_node in tree.get_terminals():
        cell_name = muts.columns[int(leaf_node.name[4:])]
        leaf_node.name = cell_name
        leaf_node.muts = set(muts.loc[muts[cell_name] == 1].index)

    # Exclude\include cells based on regex
    exclude_cnt = 0
    if exclude != '':
        print('Exluded:')
        for leaf_node in tree.get_terminals():
            if re.fullmatch(exclude, leaf_node.name):
                exclude_cnt += 1
                parent = tree.prune(leaf_node)
                print('\t{}'.format(leaf_node.name))
        if exclude_cnt == 0:
            print('\nWARNING: no samples with pattern {} in vcf!\n' \
                .format(exclude))
    if include != '':
        print('Include:')
        for leaf_node in tree.get_terminals():
            if not re.fullmatch(include, leaf_node.name):
                exclude_cnt += 1
                tree.prune(leaf_node)
            else:
                if not re.fullmatch(exclude, leaf_node.name):
                    print('\t{}'.format(leaf_node.name))
        if exclude_cnt == 0:
            print('\nWARNING: no samples with pattern {} in vcf!\n' \
                .format(include))

    # Add number of mutations to internal nodes
    for int_node in tree.get_nonterminals(order='postorder'):
        int_node.name = '+'.join([i.name for i in int_node.get_terminals()])

        try:
            child0 = int_node.clades[0]
            child1 = int_node.clades[1]
        except:
            import pdb; pdb.set_trace()

        int_node.muts = child0.muts.intersection(child1.muts)
        child0.age = max_age - len(child0.muts)
        child1.age = max_age - len(child1.muts)
        child0.mut_no = len(child0.muts.difference(int_node.muts))
        child1.mut_no = len(child1.muts.difference(int_node.muts))

    tree.root.mut_no = len(tree.root.muts)
    tree.root.age = max_age
    tree.root.lambd = 'root'

    # write_tree(tree, 'example.newick')
    # import pdb; pdb.set_trace()

    return tree


def write_tree(tree, out_file='example.newick'):
    for i in tree.find_clades():
        i.branch_length = i.mut_no
    tree.ladderize(reverse=False)
    Phylo.write(tree, out_file, 'newick')


def show_tree(tree, dendro=False, br_mut=True):
    tree.ladderize(reverse=False)

    if dendro:
        max_depth = max(tree.depths().values())
        for i in tree.find_clades(terminal=False, order='postorder'):
            int_len = i.name.count('+')
            for child in i.clades:
                child.branch_length = int_len - child.name.count('+')

    if br_mut:
        br_labels = lambda c: c.mut_no
    else:
        br_labels = lambda c: c.branch_length

    Phylo.draw(tree, branch_labels=br_labels,
        label_func=lambda c: c.name if c.name.count('+') == 0 else '',
        subplots_adjust=({'left': 0.01, 'bottom': 0.01, 'right': 0.99, 'top': 0.99})
    )


def get_sign(x):
    if x < 0:
        return '-'
    else:
        return '+'


def get_model_data(tree, min_dist=0):
    # Lambdas are assigned branchwise from the tips to the root, following the 
    #   classical strict molecular clock

    internal_nodes = [i for i in tree.find_clades(terminal=False)]

    br_indi_no = tree.count_terminals() - 1 
    br_no = br_indi_no + len(internal_nodes)
    
    X = np.zeros((br_no, br_indi_no), dtype=int)
    Y = np.zeros(br_no, dtype=float)
    constr = np.zeros((br_indi_no, br_no), dtype=int)
    init = np.zeros(br_no)
    init_p = np.zeros(br_no)

    # Get nodes in reverse depth first order
    br_dist = []
    for int_node in internal_nodes:
        dist = tree.distance(int_node)
        no_desc = int_node.name.count('+')
        br_id = dist + (1 - no_desc / 1000)
        br_dist.append((br_id, int_node))
    sorted_br = sorted(br_dist, key=lambda x: x[0], reverse=True)

    def get_traversal(start_node, tree, X, br_cells):
        X_traverse = np.zeros(br_indi_no, dtype=int)
        traverse = start_node
        while True:
            # Exit case: terminal node
            try:
                child0 = X[br_cells[traverse.clades[0]]]
            except IndexError:
                break
            else:
                child1 = X[br_cells[traverse.clades[1]]]

            if sum(child1 == -1) > sum(child0 == -1):
                X_traverse += child0
                traverse = traverse.clades[0]
            else:
                X_traverse += child1
                traverse = traverse.clades[1]
        return X_traverse

    # Iterate over internal nodes
    node_idx = 0
    br_cells = {}
    lambdas = {}
    for l_idx, (mut_no, int_node) in enumerate(sorted_br):
        is_terminal = [i.is_terminal() for i in int_node.clades]
        # Both terminal nodes
        if all(is_terminal):
            # Assign new lambda to terminals
            for i, terminal in enumerate(int_node.clades):
                Y[node_idx] = max(min_dist, terminal.mut_no)
                X[node_idx, l_idx] = 1
                br_cells[terminal] = node_idx
                terminal.lambd = f'+{l_idx}'
                if i == 0:
                    lambdas[l_idx] = node_idx
                    init[node_idx] = max(min_dist, Y[node_idx], 2 * LAMBDA_MIN)
                    node_idx += 1
            
        # One internal, one terminal node
        elif any(is_terminal):
            # Assign new lambda to terminal
            terminal = int_node.clades[is_terminal.index(True)]
            Y[node_idx] = max(min_dist, terminal.mut_no)
            X[node_idx, l_idx] = 1
            br_cells[terminal] = node_idx
            terminal.lambd = f'+{l_idx}'
            lambdas[l_idx] = node_idx
            init[node_idx] = max(min_dist, Y[node_idx], 2 * LAMBDA_MIN)
            node_idx += 1
            # Assign lambda sum to internal
            internal = int_node.clades[is_terminal.index(False)]
            Y[node_idx] = max(min_dist, internal.mut_no)
            X[node_idx, l_idx] = 1
            X[node_idx] -= get_traversal(internal, tree, X, br_cells)
            br_cells[internal] = node_idx
            internal.lambd = ' '.join(sorted([f'{get_sign(j)}{i}' \
                for i, j in enumerate(X[node_idx]) if j != 0]))
            
        # Both internal nodes
        else:
            shorter, longer = sorted(int_node.clades, key=lambda x: x.branch_length)
            # Assign new lambda to shorter branch
            Y[node_idx] = max(min_dist, shorter.mut_no)
            X[node_idx, l_idx] = 1
            br_cells[shorter] = node_idx
            shorter.lambd = f'+{l_idx}'
            lambdas[l_idx] = node_idx
            init[node_idx] = max(min_dist, Y[node_idx], 2 * LAMBDA_MIN)
            node_idx += 1
            # Assign lambda sum to longer
            Y[node_idx] = max(min_dist, longer.mut_no)
            X[node_idx, l_idx] = 1
            # Add (shortest) lambda traversal over shorter branch
            X[node_idx] += get_traversal(shorter, tree, X, br_cells)
            # Subtract (shortest) lambda traversal over longer branch
            X[node_idx] -= get_traversal(longer, tree, X, br_cells)
            br_cells[longer] = node_idx
            longer.lambd = ' '.join(sorted([f'{get_sign(j)}{i}' \
                for i, j in enumerate(X[node_idx]) if j != 0]))

        constr_idx = node_idx // 2
        constr[constr_idx, node_idx] = 1
        for i, j in enumerate(X[node_idx]):
            if j:
                constr[node_idx // 2, i*2] = -j
        init_val = np.matmul(X[node_idx], init[::2])
        # Make sure init is within bounds
        if init_val >= max(2 * LAMBDA_MIN, min_dist):
            init[node_idx,] = init_val
        else:
            init[node_idx-1,] += -init_val + (1 + min_dist)
            init[node_idx,] = np.matmul(X[node_idx], init[::2])
        node_idx += 1

    assert np.allclose(constr @ init, 0, atol=LAMBDA_MIN), \
        'Constraints not fulfilled for x_0'
    assert (init >= max(min_dist, 2 * LAMBDA_MIN)).all(), \
        f'Init value smaller than min. distance: {init.min()}'

    return Y, constr, init


def get_LRT_poisson(Y, constr, init, short=True):

    if short:
        ll_H1 = np.nansum(Y * np.log(np.where(Y > 0, Y, np.nan)) \
            - np.where(Y > 0, Y, np.nan))
        scale_fac = (ll_H1 // 1000) * 1000
        def fun_opt(l, Y):
            return -np.nansum(Y * np.log(np.where(l>LAMBDA_MIN, l, np.nan)) - l) \
                / scale_fac
    else:
        ll_H1 = np.sum(poisson.logpmf(Y, Y))
        scale_fac = (-ll_H1 // 10) * 10
        def fun_opt(l, Y):
            l[:] = np.clip(l, LAMBDA_MIN, None)
            return -np.sum(poisson.logpmf(Y, l)) / scale_fac


    def fun_jac(l, Y):
        return -(Y / np.clip(l, LAMBDA_MIN, None) - 1) / scale_fac
                

    def fun_hess(l, Y):
        return np.identity(Y.size) * (Y / np.clip(l, LAMBDA_MIN, None)**2) / scale_fac


    const = [{'type': 'eq', 'fun': lambda x: np.matmul(constr, x)}]
    init = np.clip(init, LAMBDA_MIN, None)
    bounds = np.full((Y.size, 2), (LAMBDA_MIN, Y.sum()))

    opt = minimize(fun_opt, init, args=(Y,), constraints=const,
        bounds=bounds, method='trust-constr', jac=fun_jac, hess=fun_hess,
        options={'disp': False, 'maxiter': 10000, 'barrier_tol': 1e-5})

    # opt = minimize(fun_opt, init, args=(Y,), constraints=const,
    #     bounds=bounds, method='SLSQP', jac=fun_jac,
    #     options={'disp': False, 'maxiter': 10000,})

    if not opt.success:
        print('\nFAILED POISSON OPTIMIZATION with trust-constr\n')
        opt = minimize(fun_opt, init, args=(Y,), constraints=const,
            bounds=bounds, method='SLSQP', jac=fun_jac,
            options={'disp': False, 'maxiter': 10000,})
        if not opt.success:
            print('\nFAILED POISSON OPTIMIZATION with SLSQP\n')
            return np.nan, np.nan, np.nan, np.nan, np.nan,

    if short:
        ll_H0 = np.nansum(Y * np.log(opt.x) - opt.x)
    else:
        ll_H0 = np.sum(poisson.logpmf(Y, opt.x))

    dof = Y.size - constr.shape[0]
    LR = -2 * (ll_H0 - ll_H1)
    # LR_test = -2 * (np.nansum(Y * np.log(opt2.x) - opt2.x) - ll_H1)

    on_bound = np.sum(opt.x <= np.sqrt(LAMBDA_MIN))
    if on_bound > 0:
        dof_diff = np.arange(on_bound + 1)
        weights = scipy.special.binom(on_bound, dof_diff)
        p_vals = np.clip(chi2.sf(LR, dof - dof_diff), 1e-100, 1)
        p_val = np.average(p_vals, weights=weights)
    else:
        p_val = chi2.sf(LR, dof)

    return ll_H0, ll_H1, LR, dof - on_bound / 2, p_val


def get_LRT_poisson_relaxed(Y, constr, init, relax_idx):
    scale_fac = (np.nansum(Y * np.log(np.where(Y > 0, Y, np.nan)) \
        - np.where(Y > 0, Y, np.nan)) // 1000) * 1000

    def fun_opt(l, Y):
        return -np.nansum(Y * np.log(np.where(l>LAMBDA_MIN, l, np.nan)) - l) \
            / scale_fac

    def fun_jac(l, Y):
        return -(Y / np.clip(l, LAMBDA_MIN, None) - 1) / scale_fac

    def fun_hess(l, Y):
        return np.identity(Y.size) * (Y / np.clip(l, LAMBDA_MIN, None)**2) / scale_fac

    constr_H1 = np.delete(constr, relax_idx, axis=0)
    const_H1 = [{'type': 'eq', 'fun': lambda x: np.matmul(constr_H1, x)}]
    init = np.clip(init, LAMBDA_MIN, None)
    bounds = np.full((Y.size, 2), (LAMBDA_MIN, Y.sum()))

    opt_H1 = minimize(fun_opt, init, args=(Y,), constraints=const_H1,
        bounds=bounds, method='trust-constr', jac=fun_jac, hess=fun_hess,
        options={'disp': False, 'maxiter': 10000, 'barrier_tol': 1e-5})

    if not opt_H1.success:
        print('\nFAILED OPTIMIZATION with trust-constr\n')
        opt_H1 = minimize(fun_opt, init, args=(Y,), constraints=const_H1,
            bounds=bounds, method='SLSQP', jac=fun_jac,
            options={'disp': False, 'maxiter': 10000,})
        if not opt_H1.success:
            print('\nFAILED POISSON OPTIMIZATION with SLSQP\n')
            return np.nan, np.nan, np.nan, np.nan, np.nan,

    return np.nansum(Y * np.log(opt_H1.x) - opt_H1.x)


def test_data(vcf_file, tree_file, out_file, paup_exe, exclude='', include='',
            alpha=0.05, tree_format='scite'):

    muts = get_mut_df(vcf_file)
    tree = get_tree_dict(tree_file, muts, paup_exe,  exclude, include, tree_format)
    
    Y, constr, init = get_model_data(tree, 0)
    write_tree(tree, out_file.replace('tsv', 'newick'))


    cols = ['alt. model', 'H0', 'H1', '-2logLR', 'dof', 'p-value', 'q-value',
        'hypothesis']
    header_str = '\t'.join(cols)
    model_str = ''

    # Strict clock
    ll_H0, ll_H1, LR, dof, p_val = get_LRT_poisson(Y, constr, init)
    hyp = f'H{int(p_val < alpha)}'
    model_str += f'constraint free\t{ll_H0:0>5.2f}\t{ll_H1:0>5.2f}\t' \
        f'{LR:0>5.2f}\t{dof}\t{p_val:.2E}\t{p_val:.2E}\t{hyp}'

    n = constr.shape[0]
    for i in np.arange(n):
        print(f'Testing the relaxed clock\t({i + 1:0>2} / {n})')
        ll_H1 = get_LRT_poisson_relaxed(Y, constr, init, i)
        dof = 1
        LR = -2 * (ll_H0 - ll_H1)
        p_val = chi2.sf(LR, dof)
        # Multiple testing correction with Bonferroni
        q_val = min(p_val * n, 1)
        hyp = f'H{int(q_val < alpha)}'
        model_str += f'\nrelaxed clock ({i:0>2})\t{ll_H0:0>5.2f}\t{ll_H1:0>5.2f}\t' \
        f'{LR:0>5.2f}\t{dof}\t{p_val:.2E}\t{q_val:.2E}\t{hyp}'

    with open(out_file, 'w') as f_out:
        f_out.write(f'{header_str}\n{model_str}')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf', type=str, help='SNP file in vcf format')
    parser.add_argument('tree', type=str, help='Tree file in newick format')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-e', '--exe', type=str, help='Path to PAUP exe.')
    parser.add_argument('-a', '--alpha', type=float, default=0.05,
        help='Significance threshold. Default = 0.05.')
    parser.add_argument('-tf', '--tree_format', type=str, default='scite',
        choices=['scite', 'cellphy'], help='Tree format. Default = scite')
    parser.add_argument('-excl', '--exclude', type=str, default='',
        help='Regex pattern for samples to exclude from LRT test,')
    parser.add_argument('-incl', '--include', type=str, default='',
        help='Regex pattern for samples to include from LRT test,')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        if snakemake.params.exclude == None:
            snakemake.params.exclude = ''
        if snakemake.params.include == None:
            snakemake.params.include = ''
        test_data(snakemake.input.vcf, snakemake.input.tree, 
            snakemake.output[0], snakemake.params.paup_exe,
            snakemake.params.exclude, snakemake.params.include,
            tree_format='scite')
    else:
        import argparse
        args = parse_args()
        if not args.output:
            args.output = os.path.join(os.path.dirname(args.vcf),
                'poisson_tree_LRT.tsv')
        test_data(args.vcf, args.tree, args.output, args.exe, args.exclude,
            args.include, args.alpha, args.tree_format)