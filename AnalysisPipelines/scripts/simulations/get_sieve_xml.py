#!/usr/bin/env python3

import os
import re


def get_sieve_xml(template_file, tree_file, samples_file, model, steps, base_no,
            out_file):
    with open(template_file, 'r') as f_temp:
        templ = f_temp.read()

    with open(samples_file, 'r') as f_samp:
        samples = f_samp.read().split('\n')

    tree_node = '<tree id="tree" name="stateNode" ' \
        'nodetype="beast.evolution.tree.ScsNode" ' \
        'spec="beast.evolution.tree.ScsTree" ' \
        'treeFileName="{}"/>'.format(tree_file)

    # TODO <NB> Remove curly brackets in XML
    templ = re.sub('{branch_no}', str(2 * (len(samples) -1) - 1), templ)
    templ = re.sub('{tree}', tree_node, templ)
    templ = re.sub('{bgSites}', str(base_no), templ)
    templ = re.sub('{steps}', str(steps), templ)

    if model == 'clock':
        prior_node = ''
        model_node = 'id="strictClock" ' \
            'spec="beast.evolution.branchratemodel.StrictClockModel"'
        param_node = '<parameter estimate="false" id="clockRate" ' \
            'name="clock.rate" spec="parameter.RealParameter">1.0</parameter>'
        op_node = ''
        log_node = ''
        br_rate_node = '<log id="TreeWithMetaDataLogger" '\
            'spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@tree" />'
    else:
        prior_node = """<prior id="ucldStdevPrior" name="distribution" x="@ucldStdev">
                    <Gamma id="ucldStdevDist" name="distr">
                        <parameter estimate="false" id="ucldStdevDistRealParameter.1" name="alpha" spec="parameter.RealParameter">0.5396</parameter>
                        <parameter estimate="false" id="ucldStdevDistRealParameter.2" name="beta" spec="parameter.RealParameter">0.3819</parameter>
                    </Gamma>
                </prior>"""
        model_node = 'id="RelaxedClock" spec="beast.evolution.branchratemodel.UCRelaxedClockModel" rateCategories="@rateCategories" tree="@tree"'
        param_node = """<LogNormal id="LogNormalDistributionModel" S="@ucldStdev" meanInRealSpace="true" name="distr">
                            <parameter id="LogNormalDistributionModelRealParameter" spec="parameter.RealParameter" estimate="false" lower="0.0" name="M" upper="1.0">1.0</parameter>
                        </LogNormal>
                        <parameter id="ucldMean" spec="parameter.RealParameter" estimate="false" name="clock.rate">1.0</parameter>"""
        op_node = """<operator id="ucldStdevScaler" spec="ScaleOperator" parameter="@ucldStdev" scaleFactor="0.5" weight="3.0"/>

        <operator id="CategoriesRandomWalk" spec="IntRandomWalkOperator" parameter="@rateCategories" weight="10.0" windowSize="1"/>

        <operator id="CategoriesSwapOperator" spec="SwapOperator" intparameter="@rateCategories" weight="10.0"/>

        <operator id="CategoriesUniform" spec="UniformOperator" parameter="@rateCategories" weight="10.0"/>
"""

        log_node = """<log idref="ucldStdev"/>
            <log id="rate" spec="beast.evolution.branchratemodel.RateStatistic" branchratemodel="@RelaxedClock" tree="@tree"/>
"""
        br_rate_node = '<log id="TreeWithMetaDataLogger" ' \
            'spec="beast.evolution.tree.TreeWithMetaDataLogger" ' \
            'substitutions="true" branchratemodel="@RelaxedClock" tree="@tree"/>'

    templ = re.sub('{priors}', prior_node, templ)
    templ = re.sub('{model}', model_node, templ)
    templ = re.sub('{model_params}', param_node, templ)
    templ = re.sub('{relaxed_clock_op}', op_node, templ)
    templ = re.sub('{relaxed_clock_log}', log_node, templ)
    templ = re.sub('{br_rate_log}', br_rate_node, templ)

    with open(out_file, 'w') as f_out:
        f_out.write(templ)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, nargs='+', help='Input xml files')
    parser.add_argument('-t', '--tree', type=str, help='Input tree file.')
    parser.add_argument('-s', '--samples', type=str, help='Output file.')
    parser.add_argument('-m', '--model', type=str, default='clock',
        choices=['clock', 'noClock'], help='Model to chose.')
    parser.add_argument('-n', '--steps', type=int, default=1E6,
        help='MCMC step number.')
    parser.add_argument('-b', '--background_sites', type=int, default=-1,
        help='Number of total sequences/simulated bases.')
    parser.add_argument('-o', '--output', type=str, help='Output file.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        bg = False
        with open(snakemake.input.snps, 'r') as f:
            for line in f:
                if bg:
                    base_no = int(line)
                    break
                if line == '=numBackgroundSites=\n':
                    bg = True

        get_sieve_xml(snakemake.params.xml, snakemake.input.tree,
            snakemake.input.samples, snakemake.wildcards.model,
            snakemake.params.steps, base_no, snakemake.output[0])
    else:
        import argparse
        args = parse_args()
        get_sieve_xml(args.input, args.tree, args.samples, args.model,
            args.steps, args.background_sites, args.output)