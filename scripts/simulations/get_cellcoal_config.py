#!/usr/bin/env python3

import os
import re
import yaml


def get_cellcoal_config(config, template_file, out_dir, user_tree=''):
    with open(template_file, 'r') as f:
        templ = f.read()

    model = config['cellcoal']['model']

    templ = re.sub('{no_sites}', str(model.get('no_sites', 10000)), templ)
    templ = re.sub('{pop_size}', str(model.get('pop_size', 10000)), templ)
    templ = re.sub('{out_dir}', out_dir, templ)

    cell_no = model.get('no_cells', 30)
    tree_flag = '6'
    if user_tree:
        rep_no = 1
        cell_depth = config['sc_from_bulk'].get('NGS', {}).get('seq_cov', 20)
        cell_no += 1
        templ = re.sub('{outgroup_branch_length}', '0', templ)
        templ = re.sub('{root_branch_length}', '0', templ)
    else:
        rep_no = config['cellcoal'].get('no_rep', 10)
        cell_depth = config['cellcoal'].get('NGS', {}).get('seq_cov', 20)
        templ = re.sub('{outgroup_branch_length}',
            str(model.get('outgroup_branch_length', 1)), templ)
        templ = re.sub('{root_branch_length}',
            str(model.get('root_branch_length', 0.5)), templ)

    templ = re.sub('{no_rep}', str(rep_no), templ)
    templ = re.sub('{seq_cov}', str(cell_depth), templ)
    templ = re.sub('{no_cells}', str(cell_no), templ)
    templ = re.sub('{out_tree}', tree_flag, templ)

    if model.get('no_muts', None):
        templ = re.sub('{no_muts}', f'j{model["no_muts"]}', templ)
    else:
        templ = re.sub('{no_muts}', '', templ)

    if model.get('mut_rate', None):
        templ = re.sub('{mut_rate}', f'u{model["mut_rate"]}', templ)
    else:
        templ = re.sub('{mut_rate}', 'u1e-6', templ)

    if model.get('branch_rate_var', None):
        templ = re.sub('{branch_rate_var}', f'i{model["branch_rate_var"]}', templ)
    else:
        templ = re.sub('{branch_rate_var}', '', templ)

    if model.get('branch_rate_switch', None):
        switches = model['branch_rate_switch']
        if not isinstance(switches, list):
            switches = [switches]
        switch_str = f'N{len(switches)} {" ".join([f"{i:.1f}" for i in switches])}'
        templ = re.sub('{rate_switches}', switch_str, templ)
    else:
        templ = re.sub('{rate_switches}', '', templ)

    if model.get('germline_rate', 0.0) > 0:
        templ = re.sub('{germline_rate}', f'c{model["germline_rate"]}', templ)
    else:
        templ = re.sub('{germline_rate}', '', templ)

    if model.get('binary_alphabet', False):
        templ = re.sub('{alphabet}', '0', templ)
    else:
        templ = re.sub('{alphabet}', '1', templ)

    scWGA = config['cellcoal']['scWGA']
    if scWGA.get('errors', False) or \
            (user_tree and config.get('sc_from_bulk', {}).get('errors', True)):
        if isinstance(scWGA['ADO_rate'], float):
            templ = re.sub('{ADO_rate}', f'D{scWGA["ADO_rate"]}', templ)
            templ = re.sub('{ADO_rate_var}', '', templ)
        else:
            if len(scWGA["ADO_rate"]) == 1:
                templ = re.sub('{ADO_rate}', f'D{scWGA["ADO_rate"][0]}', templ)
                templ = re.sub('{ADO_rate_var}', '', templ)
            else:
                templ = re.sub('{ADO_rate}', '', templ)
                templ = re.sub('{ADO_rate_var}',
                    f'P{scWGA["ADO_rate"][0]} {scWGA["ADO_rate"][1]}', templ)
        templ = re.sub('{ampl_error}', f'A{scWGA["ampl_error"][0]} ' \
            f'{scWGA["ampl_error"][1]} {scWGA["ampl_error"][2]}', templ)
        templ = re.sub('{doublet_rate}',
            f'B{scWGA["doublet_rate"][0]} {scWGA["doublet_rate"][1]}', templ)
    else:
        templ = re.sub('{ADO_rate}', '', templ)
        templ = re.sub('{ADO_rate_var}', '', templ)
        templ = re.sub('{ampl_error}', '', templ)
        templ = re.sub('{doublet_rate}', '', templ)

    NGS = config['cellcoal']['NGS']
    if NGS.get('errors', False) or \
            (user_tree and config.get('sc_from_bulk', {}).get('errors', True)):
        if not user_tree and NGS["seq_overdis"] > 0:
            templ = re.sub('{seq_overdis}', f'V{NGS["seq_overdis"]}', templ)
        elif user_tree and config['sc_from_bulk']['NGS']['seq_overdis'] > 0:
            templ = re.sub('{seq_overdis}',
                f'V{config["sc_from_bulk"]["NGS"]["seq_overdis"]}', templ)
        else:
            templ = re.sub('{seq_overdis}', '', templ)
        templ = re.sub('{seq_error}', f'E{NGS["seq_error"]}', templ)
    else:
        templ = re.sub('{seq_overdis}', '', templ)
        templ = re.sub('{seq_error}', '', templ)

    if config['cellcoal'].get('output', {}).get('observed_haplotype', False):
        templ = re.sub('{out_full_hap}', '4', templ)
    else:
        templ = re.sub('{out_full_hap}', '', templ)

    if config['cellcoal'].get('output', {}).get('true_haplotype', False):
        templ = re.sub('{out_true_hap}', 'Y', templ)
    else:
        templ = re.sub('{out_true_hap}', '', templ)

    if user_tree:
        templ = re.sub('{tree_file}', f'T{user_tree}', templ)
    else:
        templ = re.sub('{tree_file}', '', templ)

    return templ


if __name__ == '__main__':
    if not 'snakemake' in globals():
        raise IOError('Script only works with snakemake object')


    if snakemake.input.get('tree', False):
        tree = snakemake.input.tree
        out_dir = f'{snakemake.params.out_dir}.{snakemake.wildcards.run}'
    else:
        tree = ''
        out_dir = snakemake.params.out_dir

    cc_config = get_cellcoal_config(snakemake.config,
        snakemake.params.template, out_dir, tree)

    if snakemake.output.get('snakemake', False):
        with open(snakemake.output.snakemake, 'w') as f_yaml:
            yaml.dump(snakemake.config, f_yaml)

    with open(snakemake.output.cellcoal, 'w') as f:
        f.write(cc_config)