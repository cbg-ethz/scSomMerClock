#!/usr/bin/env python3

import os
import sys
from itertools import product

BASE_DIR = workflow.basedir
DATA_DIR = config['static_data']['data_path']
NAME = os.path.basename(DATA_DIR)
workdir: DATA_DIR

if not os.path.exists('logs'):
    os.mkdir('logs')


cell_map = {}
with open(config['static_data']['cellnames'], 'r') as f:
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
if config['static_data'].get('bulk_normal', False):
    cells_exclude = [config['static_data']['bulk_normal']]
else:
    cells_exclude = []

if config['static_data'].get('exclude_samples', False):
    cells_ex = config['static_data']['exclude_samples']
    if isinstance(cells_ex, str):
        cells_exclude.append(cells_ex)
    elif isinstance(cells_ex, list):
        cells_exclude.extend(cells_ex)


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
        monovar = [os.path.join('Calls', f'{NAME}.real.{i}.monovar.vcf') \
            for i in chrom]
        final_files.extend(monovar)
    return final_files

rule all:
    input:
        get_final_vcfs


rule adapter_cutting:
    input:
        os.path.join('Raw_Data', '{sample}_1.fastq.gz'),
        os.path.join('Raw_Data', '{sample}_2.fastq.gz')
    output:
        os.path.join('Processing', '{sample}.trimmed_1.fastq.gz'),
        os.path.join('Processing', '{sample}.trimmed_2.fastq.gz')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('cutadapt', ['cutadapt'])]),
        ref_genome = os.path.join(config['static_data']['resources_path'],
            config['static_data']['WGA_ref']),
        WGA_lib = config['static_data']['WGA_library']
    shell:
        '{params.base_dir}/scripts/1_fastqc.sh {params.modules} '
        '-s {wildcards.sample} -r {params.ref_genome} -l {params.WGA_lib}'
    

rule allignment1:
    input:
        os.path.join('Processing', '{sample}.trimmed_1.fastq.gz'),
        os.path.join('Processing', '{sample}.trimmed_2.fastq.gz')
    output:
        os.path.join('Processing', '{sample}.sam')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('bwa', ['bwa'])]),
        ref_genome = os.path.join(config['static_data']['resources_path'],
            config['static_data']['WGA_ref']),
        WGA_lib = config['static_data']['WGA_library']
    shell:
        '{params.base_dir}/scripts/2.1_bwa.sh {params.modules} '
        '-s {wildcards.sample} -r {params.ref_genome} -l {params.WGA_lib}'


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


rule base_recal1:
    input:
        os.path.join('Processing', '{cell}.dedup.bam')
    output:
        os.path.join('Processing', '{cell}.recal.table')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('gatk', ['gatk'])]),
        ref_genome = os.path.join(config['static_data']['resources_path'],
            config['static_data']['WGA_ref']),
        dbsnp = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['dbsnp']),
        indels1 = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['indel_db1'])
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
        ref_genome = os.path.join(config['static_data']['resources_path'],
            config['static_data']['WGA_ref'])
    shell:
        '{params.base_dir}/scripts/4.2_base_recal.sh {params.modules} '
        '-s {wildcards.cell} -r {params.ref_genome}'


rule indel_reallignment0:
    input:
        bams = expand(os.path.join('Processing', '{cell}.recal.bam'),
            cell=cell_map.keys())
    output:
        map_file = os.path.join('Reallignment', f'{NAME}.{{chr}}.map')
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
        os.path.join('Reallignment', f'{NAME}.{{chr}}.intervals')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('gatk', ['gatk'])]),
        name = NAME, 
        ref_genome = os.path.join(config['static_data']['resources_path'],
            config['static_data']['WGA_ref']),
        indels1 = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['indel_db1']),
        indels2 = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['indel_db2'])
    shell:
        '{params.base_dir}/scripts/5.1_indel_realign.sh {input.bams} '
        '{params.modules} -n {params.name} -c {wildcards.chr} '
        '-r {params.ref_genome} -i1 {params.indels1} -i2 {params.indels2}'


rule indel_reallignment2:
    input:
        bams = expand(os.path.join('Processing', '{cell}.recal.bam'),
            cell=cell_map.keys()),
        intervals = os.path.join('Reallignment', f'{NAME}.{{chr}}.intervals'),
        map_file = os.path.join('Reallignment', f'{NAME}.{{chr}}.map')
    output:
        expand(os.path.join('Processing', '{cell}.real.{{chr}}.bam'),
            cell=cell_map.keys())
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('gatk', ['gatk'])]),
        name = NAME, 
        ref_genome = os.path.join(config['static_data']['resources_path'],
            config['static_data']['WGA_ref']),
        indels1 = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['indel_db1']),
        indels2 = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['indel_db2'])
    shell:
        '{params.base_dir}/scripts/5.2_indel_realign.sh {input.bams} '
        '{params.modules} -n {params.name} -c {wildcards.chr} '
        '-r {params.ref_genome} -i1 {params.indels1} -i2 {params.indels2}'


rule SCcaller:
    input:
        os.path.join('Processing', '{cell}.real.{chr}.bam')
    output:
        os.path.join('Calls', '{cell}.real.{chr}.sccallerlab.vcf')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('SCcaller', ['pysam', 'numpy'])]),
        bulk = os.path.join(config['static_data']['bulk_normal']),
        ref_genome = os.path.join(config['static_data']['resources_path'],
            config['static_data']['WGA_ref']),
        dbsnp = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['dbsnp']),
        sccaller = config['SCcaller']['exe']
    shell:
        '{params.base_dir}/scripts/6_sccallerlab.sh {params.modules} '
        '-s {wildcards.cell} -c {wildcards.chr} -b {params.bulk} '
        '-r {params.ref_genome} -d {params.dbsnp} -e {params.sccaller}'


rule monovar0:
    input:
        expand(os.path.join('Processing', '{cell}.real.{{chr}}.bam'),
            cell=[i for i in cell_map.keys() if i not in cells_exclude])
    output:
        os.path.join('Processing', '{chr}.bamspath.txt')
    run:
        with open(output[0], 'w') as f:
            for bam_file in input:
                f.write(f'{bam_file}\n')
        

rule monovar:
    input:
        os.path.join('Processing', '{chr}.bamspath.txt')
    output:
        os.path.join('Calls', '{chr}.monovar.vcf')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('monovar', ['monovar'])]),
        name = NAME,
        ref_genome = os.path.join(config['static_data']['resources_path'],
            config['static_data']['WGA_ref']),
    shell:
        '{params.base_dir}/scripts/7_monovar.sh {params.modules} '
        '-c {wildcards.chr} -r {params.ref_genome}'


# ------------------------------------------------------------------------------
# ------------------------------ SEQUENCING QC ---------------------------------
# ------------------------------------------------------------------------------


rule create_bed:
    input:
        os.path.join('Processing', '{cell}.dedup.bam')
    output:
        os.path.join('Processing', '{cell}.genome.bed')
    params:
        base_dir = BASE_DIR,
        seq = config['static_data']['SEQ'],
        target = os.path.join(config['static_data']['resources_path'],
            config.get('WES', {}).get('target_path', ''))
    shell:
        '{params.base_dir}/scripts/QC_cov.sh {input} {output} {params.seq} '
        '{params.target}'


rule QC_sequencing:
    input:
        expand(os.path.join('Processing', '{cell}.genome.bed'),
            cell=cell_map.keys())
    output:
        'QC_sequencing.tsv'
    params:
        base_dir = BASE_DIR,
    shell:
        'module load python/3.7.7 numpy/1.18.1-python-3.7.7 '
        'matplotlib/3.1.3-python-3.7.7 pandas/1.0.1-python-3.7.7 && '
        'python3 {params.base_dir}/scripts/QC_coverage.py {input} -o "./"'
