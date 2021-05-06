#!/usr/bin/env python3

import os
import re
import yaml


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

    if config['cellcoal']['model'].get('mut_rate', None):
        templ = re.sub('{mut_rate}',
            'u{}'.format(config['cellcoal']['model']['mut_rate']), templ)
    else:
        templ = re.sub('{mut_rate}', 'u1e-6', templ)

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


if __name__ == '__main__':
    if not 'snakemake' in globals():
        raise IOError('Script only works with snakemake object')
    
    cc_config = get_cellcoal_config(snakemake.config, snakemake.params.template,
        snakemake.params.out_dir)
    with open(snakemake.output.cellcoal, 'w') as f:
        f.write(cc_config)
    with open(snakemake.output.snakemake, 'w') as f_yaml:
        yaml.dump(snakemake.config, f_yaml)