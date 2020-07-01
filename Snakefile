#!/usr/bin/env python3

import os
import sys


BASE_DIR = workflow.basedir
DATA_DIR = config['static_data']['data_path']
NAME = os.path.basename(DATA_DIR), 
workdir: DATA_DIR


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
                    

def get_corr_samples(wildcards):
    return [os.path.join('Processing', f'{i}.sorted.bam') \
        for i in cell_map[wildcards.cell]]


rule all:
    input:
        expand(
            os.path.join('Processing', '{cell}.real.{chr}.sccallerlab.vcf'),
            cell=cell_map.keys(), 
            chr=[i for i in range(1, 23, 1)] + ['X', 'Y']
        ),
        'QC_sequencing.tsv'


rule adapter_cutting:
    input:
        os.path.join('Raw_Data', '{sample}_1.fastq.gz'),
        os.path.join('Raw_Data', '{sample}_2.fastq.gz')
    output:
        os.path.join('Processing', '{sample}.trimmed_1.fastq.gz'),
        os.path.join('Processing', '{sample}.trimmed_2.fastq.gz')
    params:
        base_dir = BASE_DIR,
        ref_genome = os.path.join(config['static_data']['resources_path'],
            config['static_data']['WGA_ref']),
        WGA_lib = config['static_data']['WGA_library']
    shell:
        '{params.base_dir}/scripts/1_fastqc.sh {wildcards.sample} '
        '{params.ref_genome} {params.WGA_lib}'
    

rule allignment:
    input:
        os.path.join('Processing', '{sample}.trimmed_1.fastq.gz'),
        os.path.join('Processing', '{sample}.trimmed_2.fastq.gz')
    output:
        os.path.join('Processing', '{sample}.sorted.bam')
    params:
        base_dir = BASE_DIR,
        ref_genome = os.path.join(config['static_data']['resources_path'],
            config['static_data']['WGA_ref']),
        WGA_lib = config['static_data']['WGA_library']
    shell:
        '{params.base_dir}/scripts/2_bwa.sh {wildcards.sample} '
        '{params.ref_genome} {params.WGA_lib}'


rule remove_duplicates:
    input:
        get_corr_samples
    output:
        os.path.join('Processing', '{cell}.dedup.bam')
    params:
        base_dir = BASE_DIR
    shell:
        '{params.base_dir}/scripts/3_md_merge_rename.sh {input} {wildcards.cell}'


rule base_recal1:
    input:
        os.path.join('Processing', '{cell}.dedup.bam')
    output:
        os.path.join('Processing', '{cell}.recal.table')
    params:
        base_dir = BASE_DIR,
        ref_genome = os.path.join(config['static_data']['resources_path'],
            config['static_data']['WGA_ref']),
        dbsnp = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['dbsnp']),
        indels1 = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['indel_db1'])
    shell:
        '{params.base_dir}/scripts/4.1_base_recal.sh {wildcards.cell} '
        '{params.ref_genome} {params.dbsnp} {params.indels1}'


rule base_recal2:
    input:
        os.path.join('Processing', '{cell}.recal.table')
    output:
        os.path.join('Processing', '{cell}.recal.bam')
    params:
        base_dir = BASE_DIR,
        ref_genome = os.path.join(config['static_data']['resources_path'],
            config['static_data']['WGA_ref'])
    shell:
        '{params.base_dir}/scripts/4.2_base_recal.sh {wildcards.cell} '
        '{params.ref_genome}'


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
        name = NAME, 
        ref_genome = os.path.join(config['static_data']['resources_path'],
            config['static_data']['WGA_ref']),
        indels1 = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['indel_db1']),
        indels2 = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['indel_db2'])
    shell:
        '{params.base_dir}/scripts/5.1_indel_realign.sh {input.bams} '
        '-n {params.name} -c {wildcards.chr} -r {params.ref_genome} '
        '-i1 {params.indels1} -i2 {params.indels2}'


rule indel_reallignment2:
    input:
        bams = expand(os.path.join('Processing', '{cells}.recal.bam'),
            cell=cell_map.keys()),
        intervals = os.path.join('Reallignment', f'{NAME}.{{chr}}.intervals'),
        map_file = os.path.join('Reallignment', f'{NAME}.real.{{chr}}.map')
    output:
        expand(os.path.join('Processing', '{cell}.real.{{chr}}.bam'),
            cell=cell_map.keys())
    params:
        base_dir = BASE_DIR,
        name = NAME, 
        ref_genome = os.path.join(config['static_data']['resources_path'],
            config['static_data']['WGA_ref']),
        indels1 = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['indel_db1']),
        indels2 = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['indel_db2'])
    shell:
        '{params.base_dir}/scripts/5.2_indel_realign.sh {input.bams} '
        '-n {params.name} -c {wildcards.chr} -r {params.ref_genome} '
        '-i1 {params.indels1} -i2 {params.indels2}'


rule SCCaller:
    input:
        os.path.join('Processing', '{cell}.real.{chr}.bam')
    output:
        os.path.join('Processing', '{cell}.real.{chr}.sccallerlab.vcf')
    params:
        base_dir = BASE_DIR,
        bulk = os.path.join(config['static_data']['bulk_name']),
        ref_genome = os.path.join(config['static_data']['resources_path'],
            config['static_data']['WGA_ref']),
        dbsnp = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['dbsnp']),
        sccaller = config['SCCaller']['exe']
    shell:
        '{params.base_dir}/scripts/6_sccallerlab.sh {wildcards.cell} '
        '{wildcards.chr} {params.bulk} {params.ref_genome} {params.dbsnp} '
        '{params.sccaller}'

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
