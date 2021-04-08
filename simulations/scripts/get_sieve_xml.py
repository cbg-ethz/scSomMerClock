#!/usr/bin/env python3

import os
import re


def get_sieve_xml(template_file, tree_file, samples_file, model, steps,
            out_file):
    with open(template_file, 'r') as f:
        templ = f.read()

    run = re.search('\.(\d+)', tree_file).group(1)

    # bg_file = os.path.join(os.path.sep.join(tree_file.split(os.path.sep)[:-2]),
    #     'full_genotypes_dir', 'background.{}'.format(run))

    # with open(bg_file, 'r') as f_bg:
    #     bg_bases = f_bg.read().strip()

    # background_str = '<data id="alignment" spec="FilteredAlignment" filter="-" ' \
    #     'data="@original-alignment" constantSiteWeights="{}"/>'.format(bg_bases)
    background_str = ''
    templ = re.sub('{background}', background_str, templ)

    with open(tree_file, 'r') as f_tree:
        newick = f_tree.read().strip()

        sample_no = 0
        with open(samples_file, 'r') as f_smpl:
            for s_i, s_raw in enumerate(f_smpl.read().strip().split('\n')):
                s_name, _ = s_raw.split('\t')
                newick = re.sub('(?<=[\(\),]){}(?=[,\)\)])'.format(s_i + 1), s_name, newick)
        # Replace names
        newick += ';'

    tree_node = '<tree id="tree" name="stateNode" ' \
        'nodetype="beast.evolution.tree.ScsNode" ' \
        'spec="beast.evolution.tree.ScsTree" ' \
        'treeFileName="{}"/>'.format(tree_file)

    # TODO <NB> Remove curly brackets in XML
    templ = re.sub('{branch_no}', str(2 * (s_i + 1) - 1), templ)

    # '<stateNode spec="beast.util.TreeParser" id="randomTree" ' \
    #     'IsLabelledNewick="true" adjustTipHeights="false" taxa="@alignment" ' \
    #     'newick="{}"/>'.format(newick)

    templ = re.sub('{tree}', tree_node, templ)


    templ = re.sub('{steps}', str(steps), templ)

    if model == 'clock':
        prior_node = ''
        model_node = """<branchRateModel id='strictClock' spec='beast.evolution.branchratemodel.StrictClockModel'>
                                                                        
                        <parameter estimate='false' id='clockRate' name='clock.rate' spec='parameter.RealParameter'>1.0</parameter>
                                                                    
                    </branchRateModel>"""

        op_node = ''
        log_node = ''
        br_rate_node = '<log id="TreeWithMetaDataLogger" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@tree" />'
    else:
        prior_node = """<prior id="ucldStdevPrior" name="distribution" x="@ucldStdev">
                    <Gamma id="Gamma.0" name="distr">
                        <parameter estimate="false" id="RealParameter.6" name="alpha" spec="parameter.RealParameter">0.5396</parameter>
                        <parameter estimate="false" id="RealParameter.7" name="beta" spec="parameter.RealParameter">0.3819</parameter>
                    </Gamma>
                </prior>"""
        model_node = """<branchRateModel id="RelaxedClock" spec="beast.evolution.branchratemodel.UCRelaxedClockModel" rateCategories="@rateCategories" tree="@tree">
                        <LogNormal id="LogNormalDistributionModel" S="@ucldStdev" meanInRealSpace="true" name="distr">
                            <parameter id="RealParameter.5" spec="parameter.RealParameter" estimate="false" lower="0.0" name="M" upper="1.0">1.0</parameter>
                        </LogNormal>
                        <parameter id="ucldMean" spec="parameter.RealParameter" estimate="false" name="clock.rate">1.0</parameter>
                    </branchRateModel>"""
        op_node = """<operator id="ucldStdevScaler" spec="ScaleOperator" parameter="@ucldStdev" scaleFactor="0.5" weight="3.0"/>

        <operator id="CategoriesRandomWalk" spec="IntRandomWalkOperator" parameter="@rateCategories" weight="5.0" windowSize="1"/>

        <operator id="CategoriesSwapOperator" spec="SwapOperator" intparameter="@rateCategories" weight="5.0"/>

        <operator id="CategoriesUniform" spec="UniformOperator" parameter="@rateCategories" weight="5.0"/>
"""

        log_node = """<log idref="ucldStdev"/>
            <log id="rate" spec="beast.evolution.branchratemodel.RateStatistic" branchratemodel="@RelaxedClock" tree="@tree"/>
"""
        br_rate_node = '<log id="TreeWithMetaDataLogger" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@tree" branchratemodel="@RelaxedClock"/>'

    templ = re.sub('{priors}', prior_node, templ)
    templ = re.sub('{model}', model_node, templ)
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
    parser.add_argument('-o', '--output', type=str, help='Output file.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        get_sieve_xml(snakemake.input.xml, snakemake.input.tree,
            snakemake.input.samples, snakemake.wildcards.model,
            snakemake.params.steps, snakemake.output[0])
    else:
        import argparse
        args = parse_args()
        get_sieve_xml(args.input, args.tree, args.samples, args.model,
            args.steps, args.output)