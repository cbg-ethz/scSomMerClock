#!/usr/bin/env python3

import sys
import os

BASE_DIR = config['static_data']['data_path']

cell_map = {}
sample_list = []
with open(config['static_data']['cellnames'], 'r') as f:
    lines = f.read().strip().split('\n')
    for line in lines:
        row = line.split('\t')
        if len(row) == 1:
            raise IOError('cellnames file contains only 1 columns')
        elif len(row) > 3:
            raise NotImplementedError(
                'Pipeline only implemented for 1 or 2 samples per cell')
        else:
            cell_map[row[-1]] = row[:-1]
            sample_list.extend(row[:-1])
                    

rule all:
    input:
        expand(
            os.path.join('Processing', '{cell_real}.real.{chr}.sccallerlab.vcf'),
            cell_real=cell_map.keys(), 
            chr=[i for i in range(1,23,1)] + ['X', 'Y']
        )


rule adapter_cutting:
    output:
        os.path.join('slurm_files', '{sample}_CutadaptWGA.txt')
    params:
        ref_genome = os.path.join(config['static_data']['resources_path'],
            config['static_data']['WGA_ref']),
        WGA_lig = os.path.join(config['static_data']['resources_path'],
            config['static_data']['WGA_library'])
    shell:
        'scripts/1_fastqc.sh {sample} {params.ref_genome} {params.WGA_lig}'
    

rule allignment:
    input:
        os.path.join('Processing', '{sample}.trimmed_1.fastq.gz'),
        os.path.join('Processing', '{sample}.trimmed_2.fastq.gz')
    output:
        os.path.join('Processing', '{sample}.sorted.bam')
    params:
        ref_genome = os.path.join(config['static_data']['resources_path'],
            config['static_data']['WGA_ref'])
    shell:
        'scripts/2_bwa.sh {sample} {params.ref_genome}'


rule remove_duplicates:
    input:
        expand(os.path.join('Processing', '{sample}.sorted.bam'),
            sample=cell_map['{cell}'])
    output:
        os.path.join('Processing', '{cell}.dedup.bam')
    params:
        ref_genome = os.path.join(config['static_data']['resources_path'],
            config['static_data']['WGA_ref'])
        corr_samples = ' '.join(cell_map['{cell}'])
    shell:
        'scripts/3_md_merge_rename.sh {params.corr_samples} {cell}'


rule base_recal:
    input:
        os.path.join('Processing', '{cell}.dedup.bam')
    output:
        os.path.join('Processing', '{cell}.recal.bam')
    params:
        ref_genome = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['WGA_ref'])
        dbsnp = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['dbsnp'])
        indels1 = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['indel_db1'])
    shell:
        'scripts/4_base_recal.sh {cell} {params.ref_genome} {params.dbsnp} '
        '{params.indel1}'

rule indel_reallignment:
    input:
        expand(os.path.join('Processing', '{cell}.recal.bam'),
            cell=cell_map.keys())
    output:
        os.path.join('Processing', '{cell_real}.real.{chr}.bam')
    params:
        ref_genome = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['WGA_ref'])
        indels1 = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['indel_db1']),
        indels2 = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['indel_db2'])
    shell:
        'scripts/5_indel_realign.sh {input} -r {params.ref_genome} '
        '-i1 {params.indels1} -i2 {params.indel2}'

# rule QC_alligned:
#     input:
#     output:
#     shell:


rule SCCaller:
    input:
        os.path.join('Processing', '{cell_real}.real.{chr}.bam')
    output:
        os.path.join('Processing', '{cell_real}.real.{chr}.sccallerlab.vcf')
    params:
        bulk = os.path.join(config['static_data']['bulk_name']),
        ref_genome = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['WGA_ref'])
        dbsnp = os.path.join(config['static_data']['resources_path'],
            config['base_recal']['dbsnp'])
        sccaller = config['SCCaller']['exe']
    shell:
        'scripts/6_sccallerlab.sh {cell_real} {chr} {params.bulk} '
        '{params.ref_genome} {params.dbsnp} {params.sccaller}'