#!/usr/bin/env python3

import os
import sys
from itertools import product

BASE_DIR = workflow.basedir
DATA_DIR = config['specific']['data_path']
RES_PATH = config['static']['resources_path']
workdir: DATA_DIR

if not os.path.exists('logs'):
    os.mkdir('logs')


cell_map = {}
with open(config['specific']['cellnames'], 'r') as f:
    lines = f.read().strip().split('\n')
    for line in lines:
        row = line.split('\t')
        if len(row) == 1:
            raise IOError('cellnames file contains only 1 columns.')
        elif len(row) > 3:
            raise NotImplementedError(
                'Only implemented for 1 or 2 samples per cell.')
        else:
            cell_map[row[-1]] = row[:-1]

# Get samples to exclude for Monovar SNV calling
if config['specific'].get('bulk_normal', False):
    bulk_samples = set([config['specific']['bulk_normal']])
else:
    bulk_samples = set([])

if config['specific'].get('bulk_samples', False):
    bulk = config['specific']['bulk_samples']
    if isinstance(bulk, str) :
        bulk_samples.add(bulk)
    elif isinstance(bulk, list):
        bulk_samples.union(bulk)


def get_corr_samples(wildcards):
    return [os.path.join('Processing', f'{i}.sorted.bam') \
        for i in cell_map[wildcards.cell]]


def get_final_vcfs(wildcards):
    chrom = [i for i in range(1, 23, 1)] + ['X', 'Y']
    final_files = ['QC_sequencing.tsv']
    if config.get('SCcaller', {}).get('run', False):
        sccaller = [os.path.join('Calls', f'{i[0]}.real.{i[1]}.sccallerlab.vcf') \
            for i in product(cell_map, chrom) ]
        final_files.extend(sccaller)
    if config.get('monovar', {}).get('run', False):
        monovar = os.path.join('Calls', 'all.monovar.vcf.gz')
        final_files.append(monovar)
    if bulk_samples:
        mutect = [os.path.join('Calls', f'{i}.mutect.vcf') for i in chrom]
        final_files.extend(mutect)
    return final_files


rule all:
    input:
        get_final_vcfs


rule adapter_cutting:
    input:
        os.path.join('Raw_Data', '{sample}_1.fastq.gz')
    output:
        os.path.join('Processing', '{sample}.trimmed_1.fastq.gz')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('cutadapt', ['cutadapt'])]),
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        pair_end = ' ' if config['specific'].get('pair_end', True) else '-se',
        WGA_lib = config['specific']['WGA_library']
    shell:
        '{params.base_dir}/scripts/1_fastqc.sh {params.modules} '
        '-s {wildcards.sample} -r {params.ref_genome} -l {params.WGA_lib} '
        '{params.pair_end}'


rule allignment1:
    input:
        os.path.join('Processing', '{sample}.trimmed_1.fastq.gz')
    output:
        os.path.join('Processing', '{sample}.sam')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('bwa', ['bwa'])]),
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        pair_end = ' ' if config['specific'].get('pair_end', True) else '-se',
        WGA_lib = config['specific']['WGA_library']
    shell:
        '{params.base_dir}/scripts/2.1_bwa.sh {params.modules} '
        '-s {wildcards.sample} -r {params.ref_genome} -l {params.WGA_lib} '
        '{params.pair_end}'


rule allignment2:
    input:
        os.path.join('Processing', '{sample}.sam')
    output:
        os.path.join('Processing', '{sample}.sorted.bam')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('picard', ['picard'])])
    shell:
        '{params.base_dir}/scripts/2.2_bwa.sh {params.modules} '
        '-s {wildcards.sample}'


rule remove_duplicates:
    input:
        get_corr_samples
    output:
        os.path.join('Processing', '{cell}.dedup.bam')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('picard', ['picard'])])
    shell:
        '{params.base_dir}/scripts/3_md_merge_rename.sh {input} '
        '{params.modules} -s {wildcards.cell}'


rule sanity_check_bam:
    input:
        expand(os.path.join('Processing', '{cell}.dedup.bam'),
            cell=cell_map.keys())
    output:
        'bad_bams.fofn'
    params:
        samtools = config['modules'].get('samtools', ['samtools'])
    shell:  
        'module load {params.samtools};'
        'samtools quickcheck -v {input} > {output[0]};'
        'xargs rm < bad_bams.fofn'


rule base_recal1:
    input:
        os.path.join('Processing', '{cell}.dedup.bam')
    output:
        os.path.join('Processing', '{cell}.recal.table')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('gatk', ['gatk'])]),
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        dbsnp = os.path.join(RES_PATH, config['static']['dbsnp']),
        indels1 = os.path.join(RES_PATH, config['static']['indel_db1'])
    shell:
        '{params.base_dir}/scripts/4.1_base_recal.sh {params.modules} '
        '-s {wildcards.cell} -r {params.ref_genome} -d {params.dbsnp} '
        '-i {params.indels1}'


rule base_recal2:
    input:
        os.path.join('Processing', '{cell}.recal.table')
    output:
        os.path.join('Processing', '{cell}.recal.bam')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('gatk', ['gatk'])]),
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref'])
    shell:
        '{params.base_dir}/scripts/4.2_base_recal.sh {params.modules} '
        '-s {wildcards.cell} -r {params.ref_genome}'


rule indel_reallignment0:
    input:
        bams = expand(os.path.join('Processing', '{cell}.recal.bam'),
            cell=cell_map.keys())
    output:
        map_file = os.path.join('Reallignment', '{chr}.map')
    run:
        with open(output.map_file, 'w') as f:
            for bam_full in input.bams:
                bam = os.path.basename(bam_full)
                cell_id = bam.split('.')[0]
                f.write(f'{bam}\tProcessing/{cell_id}.real.{{chr}}.bam\n')


rule indel_reallignment1:
    input:
        bams = expand(os.path.join('Processing', '{cell}.recal.bam'),
            cell=cell_map.keys())
    output:
        os.path.join('Reallignment', '{chr}.intervals')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('gatk', ['gatk'])]),
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        indels1 = os.path.join(RES_PATH, config['static']['indel_db1']),
        indels2 = os.path.join(RES_PATH, config['static']['indel_db2'])
    shell:
        '{params.base_dir}/scripts/5.1_indel_realign.sh {input.bams} '
        '{params.modules} -c {wildcards.chr} -r {params.ref_genome} '
        '-i1 {params.indels1} -i2 {params.indels2}'


rule indel_reallignment2:
    input:
        bams = expand(os.path.join('Processing', '{cell}.recal.bam'),
            cell=cell_map.keys()),
        intervals = os.path.join('Reallignment', '{chr}.intervals'),
        map_file = os.path.join('Reallignment', '{chr}.map')
    output:
        expand(os.path.join('Processing', '{cell}.real.{{chr}}.bam'),
            cell=cell_map.keys())
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('gatk', ['gatk'])]),
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        indels1 = os.path.join(RES_PATH, config['static']['indel_db1']),
        indels2 = os.path.join(RES_PATH, config['static']['indel_db2'])
    shell:
        '{params.base_dir}/scripts/5.2_indel_realign.sh {input.bams} '
        '{params.modules} -c {wildcards.chr} -r {params.ref_genome} '
        '-i1 {params.indels1} -i2 {params.indels2}'


rule SCcaller:
    input:
        os.path.join('Processing', '{cell}.real.{chr}.bam')
    output:
        os.path.join('Calls', '{cell}.real.{chr}.sccallerlab.vcf')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('SCcaller', ['pysam', 'numpy'])]),
        bulk = config['specific']['bulk_normal'],
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        dbsnp = os.path.join(RES_PATH, config['static']['dbsnp']),
        sccaller = config['SCcaller']['exe']
    shell:
        '{params.base_dir}/scripts/6_sccallerlab.sh {params.modules} '
        '-s {wildcards.cell} -c {wildcards.chr} -b {params.bulk} '
        '-r {params.ref_genome} -d {params.dbsnp} -e {params.sccaller}'


rule monovar0:
    input:
        expand(os.path.join('Processing', '{cell}.real.{{chr}}.bam'),
            cell=[i for i in cell_map.keys() if i not in bulk_samples])
    output:
        os.path.join('Processing', '{chr}.bamspath.txt')
    run:
        with open(output[0], 'w') as f:
            for bam_file in input:
                f.write(f'{bam_file}\n')
        

rule monovar1:
    input:
        os.path.join('Processing', '{chr}.bamspath.txt')
    output:
        os.path.join('Calls', '{chr}.monovar.vcf.gz')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('monovar', ['monovar'])]),
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),

    shell:
        '{params.base_dir}/scripts/7-1_monovar.sh {params.modules} '
        '-c {wildcards.chr} -r {params.ref_genome}'


rule monovar2:
    input:
        expand(os.path.join('Calls', '{chr}.monovar.vcf.gz'),
            chr=[i for i in range(1, 23, 1)] + ['X', 'Y'])
    output:
        os.path.join('Calls', 'all.monovar.vcf.gz')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('bcftools', ['bcftools'])]),
    shell:
        '{params.base_dir}/scripts/7.2_monovar.sh {input} {params.modules} '
        ' -o {output[0]}'


rule mutect:
    input: 
        expand(os.path.join('Processing', '{cell}.real.{{chr}}.bam'), 
            cell=bulk_samples)
    output:
        os.path.join('Calls', '{chr}.mutect.vcf')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('gatk', ['gatk'])]),
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        germ_res = os.path.join(RES_PATH, config['static']['germline']),
        pon = os.path.join(RES_PATH, config['specific']['PON']),
        bulk = ' '.join([f'-b {i}' for i in bulk_samples]),
        normal = f'-n {cell_map[config["specific"]["bulk_normal"]][0]}'
    shell:
        '{params.base_dir}/scripts/8_mutect.sh {params.modules} '
        '-c {wildcards.chr} -r {params.ref_genome} -g {params.germ_res} '
        '-p {params.pon} {params.bulk} {params.normal}'

# ------------------------------------------------------------------------------
# ------------------------------ SEQUENCING QC ---------------------------------
# ------------------------------------------------------------------------------


rule create_bed:
    input:
        os.path.join('Processing', '{cell}.dedup.bam')
    output:
        os.path.join('Processing', '{cell}.genome.tsv')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('bedtools', ['bedtools'])]),
        seq = config['specific']['SEQ'],
        target = os.path.join(RES_PATH, 
            config['specific'].get('WES_target', '-1')),
        genome = os.path.join(RES_PATH, 
            config['specific'].get('WES_target_genome', '-1'))
    shell:
        '{params.base_dir}/scripts/QC_cov.sh {params.modules} '
        '-i {input} -o {output} --seq {params.seq} -e {params.target} '
        '-g {params.genome}'


rule QC_sequencing:
    input:
        expand(os.path.join('Processing', '{cell}.genome.tsv'),
            cell=cell_map.keys())
    output:
        'QC_sequencing.tsv'
    params:
        base_dir = BASE_DIR,
        modules = ' '.join( \
                config['modules'].get('QC_seq', ['pandas', 'matplotlib']))
    shell:
        'module load {params.modules} && '
        'python {params.base_dir}/scripts/QC_coverage.py {input} -o "./"'
