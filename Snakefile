#!/usr/bin/env python3

import os


workdir: config['specific']['data_path']
RES_PATH = config['static']['resources_path']
SCRIPT_DIR = os.path.join(workflow.basedir, 'scripts', 'processing')
CHROM = [i for i in range(1, 23, 1)] + ['X', 'Y']


if not os.path.exists('logs'):
    os.mkdir('logs')


cell_map = {}
with open(config['specific']['cellnames'], 'r') as f:
    lines = f.read().strip().split('\n')
    for line in lines:
        row = line.split('\t')
        if len(row) == 1:
            raise IOError('cellnames file contains only 1 columns.')
        else:
            if row[0].startswith('SRR'):
                if len(row) == 2:
                    cell_map[row[1]] = [row[0]]
                elif row[1].startswith('SRR') and len(row) == 3:
                    cell_map[row[-1]] = row[:-1]
                else:
                    raise IOError('Cannot handle row in cellnames file: {}' \
                        .format(row))
            else:
                cell_map[row[0]] = row[1:]

# Get samples to exclude for Monovar SNV calling
bulk_samples = {'normal': None, 'tumor': set([]), 'all': set([])}
if config['specific'].get('bulk_normal', False):
    bulk_samples['normal'] = config['specific']['bulk_normal']
    bulk_samples['all'].add(config['specific']['bulk_normal'])

if config['specific'].get('bulk_samples', False):
    bulk = config['specific']['bulk_samples']
    if isinstance(bulk, str) :
        bulk_samples['tumor'].add(bulk)
        bulk_samples['all'].add(bulk)
    elif isinstance(bulk, list):
        bulk_samples['tumor'] = bulk_samples['tumor'].union(bulk)
        bulk_samples['all'] = bulk_samples['all'].union(bulk)

ss_samples = list(set(cell_map.keys()).difference(bulk_samples['all']))
ss_samples.sort()


rule all:
    input:
        os.path.join('Calls', 'all_filtered.vcf')
            

rule adapter_cutting:
    input:
        os.path.join('Raw_Data', '{sample}_1.fastq.gz')
    output:
        temp(os.path.join('Processing', '{sample}.trimmed_1.fastq.gz'))
    envmodules:
        'pigz',
        'cutadapt/1.18-python-3.7.0',
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16384
    params:
        pair_end = ' ' if config['specific'].get('pair_end', True) else '-se',
        WGA_lib = config['specific']['WGA_library'],
    shell:
        f'{SCRIPT_DIR}/01_fastqc.sh '
        '-s {wildcards.sample} -l {params.WGA_lib} {params.pair_end}'


rule alignment1:
    input:
        os.path.join('Processing', '{sample}.trimmed_1.fastq.gz')
    output:
        temp(os.path.join('Processing', '{sample}.sam'))
    envmodules:
        'gcc/6.4.0',
        'bwa/0.7.17',
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16384,
    threads: 8
    params:
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        pair_end = ' ' if config['specific'].get('pair_end', True) else '-se',
        WGA_lib = config['specific']['WGA_library']
    shell:
        '{SCRIPT_DIR}/02.1_bwa.sh -s {wildcards.sample} '
        '-r {params.ref_genome} -l {params.WGA_lib} {params.pair_end}'


rule alignment2:
    input:
        os.path.join('Processing', '{sample}.sam')
    output:
        temp(os.path.join('Processing', '{sample}.sorted.bam'))
    envmodules:
        'picard/2.18.14',
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16384,
    shell:
        '{SCRIPT_DIR}/02.2_bwa.sh -s {wildcards.sample}'


def get_dedup_samples(wildcards):
    return [os.path.join('Processing', f'{i}.sorted.bam') \
        for i in cell_map[wildcards.cell]]


rule remove_duplicates:
    input:
        get_dedup_samples
    output:
        os.path.join('Processing', '{cell}.dedup.bam')
    envmodules:
        'picard/2.18.14',
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32768,
    shell:
        '{SCRIPT_DIR}/03_md_merge_rename.sh {input} '
        '-s {wildcards.cell}'


rule sanity_check_bam:
    input:
        expand(os.path.join('Processing', '{cell}.dedup.bam'),
            cell=cell_map.keys())
    output:
        'bad_bams.fofn'
    envmodules:
        'samtools'
    shell:  
        'samtools quickcheck -v {input} > {output}; xargs rm < bad_bams.fofn'


rule indel_realignment0:
    input:
        bams = expand(os.path.join('Processing', '{cell}.dedup.bam'),
            cell=cell_map.keys())
    output:
        map_file = os.path.join('Realignment', '{chr}.map')
    run:
        with open(output.map_file, 'w') as f:
            for bam_full in input.bams:
                bam = os.path.basename(bam_full)
                cell_id = bam.split('.')[0]
                chrom = wildcards.chr
                f.write(f'{bam}\tProcessing/{cell_id}.real.{chrom}.bam\n')


rule indel_realignment1:
    input:
        expand(os.path.join('Processing', '{cell}.dedup.bam'),
            cell=cell_map.keys())
    output:
        os.path.join('Realignment', '{chr}.intervals')
    envmodules:
        'gatk/3.7-0-gcfedb67',
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16384
    threads: 4
    params:
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        pref = 'chr' if config['static']['WGA_ref'].startswith('hg19') else '',
        indels1 = os.path.join(RES_PATH, config['static']['indel_db1']),
        indels2 = '' if not config['static'].get('indel_db2', False) else
            '-i2 {}'.format(os.path.join(RES_PATH, config['static']['indel_db2']))
    shell:
        '{SCRIPT_DIR}/05.1_indel_realign.sh {input} '
        '-c {params.pref}{wildcards.chr} -o {output} '
        '-r {params.ref_genome} -i1 {params.indels1} {params.indels2}'


rule indel_realignment2:
    input:
        bams = expand(os.path.join('Processing', '{cell}.dedup.bam'),
            cell=cell_map.keys()),
        target = os.path.join('Realignment', '{chr}.intervals'),
        maps = os.path.join('Realignment', '{chr}.map')
    output:
        temp(expand(os.path.join('Processing', '{cell}.real.{{chr}}.bam'),
            cell=cell_map.keys()))
    envmodules:
        'gatk/3.7-0-gcfedb67',
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32768
    params:
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        pref = 'chr' if config['static']['WGA_ref'].startswith('hg19') else '',
        indels1 = os.path.join(RES_PATH, config['static']['indel_db1']),
        indels2 = '' if not config['static'].get('indel_db2', False) else
            '-i2 {}'.format(os.path.join(RES_PATH, config['static']['indel_db2']))
    shell:
        '{SCRIPT_DIR}/05.2_indel_realign.sh {input.bams} '
        ' -c {params.pref}{wildcards.chr} -r {params.ref_genome} '
        '-t {input.target} -ma {input.maps} -i1 {params.indels1} '
        '{params.indels2}'


rule base_recal:
    input:
        expand(os.path.join('Processing', '{{cell}}.real.{chr}.bam'), chr=CHROM)
    output:
        temp(os.path.join('Processing', '{cell}.real.bam'))
    envmodules:
        'samtools',
    shell:
        'samtools merge {output} {input}'


rule base_recal1:
    input:
        os.path.join('Processing', '{cell}.real.bam')
    output:
        os.path.join('Processing', '{cell}.recal.table')
    envmodules:
        'gatk/4.0.10.0'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32768
    params:
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        dbsnp = os.path.join(RES_PATH, config['static']['dbsnp']),
        indels1 = os.path.join(RES_PATH, config['static']['indel_db1']),
    shell:
        '{SCRIPT_DIR}/04.1_base_recal.sh -i {input} -o {output} '
        '-r {params.ref_genome} -d {params.dbsnp} -id {params.indels1}'


rule base_recal2:
    input:
        bam = os.path.join('Processing', '{cell}.real.bam'),
        table = os.path.join('Processing', '{cell}.recal.table')
    output:
        os.path.join('Processing', '{cell}.recal.bam')
    envmodules:
        'gatk/4.0.10.0'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32768
    params:
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
    shell:
        '{SCRIPT_DIR}/04.2_base_recal.sh -i {input.bam} '
        '-t {input.table} -o {output} -r {params.ref_genome}'


rule base_recal3:
    input:
        os.path.join('Processing', '{cell}.recal.bam')
    output:
        expand(os.path.join('Processing', '{{cell}}.recal.{chr}.bam'),
            chr=CHROM)
    envmodules:
        'samtools'
    params:
        chrom = CHROM,
    shell:
        'for chr in {params.chrom}; do '
        'samtools view {input} ${{chr}} -b > Processing/{wildcards.cell}.recal.${{chr}}.bam; '
        'samtools index Processing/{wildcards.cell}.recal.${{chr}}.bam; '
        'done'
        
    
# ------------------------------------------------------------------------------
# ----------------------------- MUTATION CALLING -------------------------------
# ------------------------------------------------------------------------------

rule SCcaller1:
    input:
        os.path.join('Processing', '{cell}.recal.{chr}.bam')
    output:
        temp(os.path.join('Calls', '{cell}.{chr}.sccaller.vcf'))
    envmodules:
        'pysam/0.15.4-python-2.7.17',
        'numpy/1.16.6-python-2.7.17',
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32768
    params:
        bulk = config['specific'].get('bulk_normal', ''),
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        dbsnp = os.path.join(RES_PATH, config['static']['dbsnp']),
        min_depth = config['filters'].get('depth', 10),
        sccaller = config['SCcaller']['exe']
    shell:
        '{SCRIPT_DIR}/06.1_sccallerlab.sh -s {wildcards.cell} '
        '-c {wildcards.chr} -b {params.bulk} -r {params.ref_genome} '
        '-d {params.dbsnp} -e {params.sccaller} -md {params.min_depth}'


rule SCcaller2:
    input:
        expand(os.path.join('Calls', '{{cell}}.{chr}.sccaller.vcf'),
            chr=CHROM)
    output:
        os.path.join('Calls', '{cell}.sccaller.vcf.gz')
    envmodules:
        'bcftools',
    shell:
        '{SCRIPT_DIR}/06.2_sccallerlab.sh {input} -o {output[0]}'


rule SCcaller3:
    input:
        expand(os.path.join('Calls', '{cell}.sccaller.vcf.gz'), cell=ss_samples)
    output:
        os.path.join('Calls', 'all.sccaller.vcf.gz')
    envmodules:
        'bcftools',
    shell:
        '{SCRIPT_DIR}/06.3_sccallerlab.sh {input} -o {output[0]}'


rule monovar0:
    input:
        expand(os.path.join('Processing', '{cell}.recal.{{chr}}.bam'),
            cell=ss_samples)
    output:
        temp(os.path.join('Processing', '{chr}.bamspath.txt'))
    run:
        with open(output[0], 'w') as f:
            for bam_file in input:
                f.write(f'{bam_file}\n')
        

rule monovar1:
    input:
        os.path.join('Processing', '{chr}.bamspath.txt')
    output:
        temp(os.path.join('Calls', '{chr}.monovar.vcf'))
    envmodules:
        'samtools',
        'pysam/0.15.4-python-2.7.17',
        'numpy/1.16.6-python-2.7.17',
        'scipy/1.1.0-python-2.7.15',
    threads: 3
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32768
    params:
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        pref = 'chr' if config['static']['WGA_ref'].startswith('hg19') else '',
        monovar = config['monovar']['exe'],
    shell:
        '{SCRIPT_DIR}/07.1_monovar.sh -c {wildcards.chr} '
        '-r {params.ref_genome} -p {params.pref} -e {params.monovar}'


rule monovar2:
    input:
        expand(os.path.join('Calls', '{chr}.monovar.vcf'),
            chr=CHROM)
    output:
        os.path.join('Calls', 'all.monovar.vcf.gz')
    envmodules:
        'bcftools',
    params:
        min_depth = config['filters'].get('depth', 10),
    shell:
        '{SCRIPT_DIR}/07.2_monovar.sh {input} -o {output}'
        '-md {params.min_depth}'


rule mutect1:
    input: 
        expand(os.path.join('Processing', '{cell}.recal.{{chr}}.bam'), 
            cell=bulk_samples['all'])
    output:
        os.path.join('Calls', '{chr}.mutect.vcf'),
        os.path.join('Calls', '{chr}.f1r2.mutect.tar.gz')
    envmodules:
        'picard/2.18.14',
        'gatk/4.1.1.0',
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32768
    params:
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        germ_res = os.path.join(RES_PATH, config['static']['germline']),
        pon = os.path.join(RES_PATH, config['specific']['PON']),
        normal = f'-n {config["specific"]["bulk_normal"]}'
            if config['specific'].get('bulk_normal', False) else ''
    shell:
        '{SCRIPT_DIR}/08.1_mutect.sh {input} -c {wildcards.chr} '
        '-r {params.ref_genome} -g {params.germ_res} -p {params.pon} '
        '{params.normal}'


rule mutect2:
    input: 
        expand(os.path.join('Calls', '{chr}.mutect.vcf'), chr=CHROM)
    output:
        os.path.join('Calls', 'all.mutect.vcf.gz')
    envmodules:
        'bcftools',
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32768
    shell:
        '{SCRIPT_DIR}/08.2_mutect.sh {input} -o {output}'


if config.get('mutect', {}).get('filter', '') == 'simple':
    rule mutect4_simple:
        input: 
            os.path.join('Calls', 'all.mutect.vcf.gz'),
        output:
            os.path.join('Calls', 'all.mutect.filtered.vcf.gz')
        envmodules:
            'gatk/4.1.1.0',
        threads: 2
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 16384
        params:
            ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        shell:
            '{SCRIPT_DIR}/08.4_mutect_simple.sh -i {input} '
            '-r {params.ref_genome} -o {output}'

else:
    rule mutect3:
        input:
            tumor = expand(os.path.join('Processing', '{cell}.recal.bam'),
                cell=bulk_samples['tumor']),
            normal = os.path.join('Processing', 
                '{}.recal.bam'.format(bulk_samples['normal'])),
            f1r2 = expand(os.path.join('Calls', '{chr}.f1r2.mutect.tar.gz'),
                chr=CHROM)
        output:
            expand(os.path.join('Calls', '{cell}.contamination.table'), 
                cell=bulk_samples['tumor']),
            os.path.join('Calls', 'read-orientation-model.tar.gz')
        envmodules:
            'gatk/4.1.1.0',
        threads: 2
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 16384
        params:
            gnomAD = os.path.join(RES_PATH, config['static']['gnomAD']),
        shell:
            '{SCRIPT_DIR}/08.3_mutect.sh {input.tumor} '
            '-n {input.normal} -gAD {params.gnomAD}'


    rule mutect4:
        input: 
            vcf = os.path.join('Calls', 'all.mutect.vcf.gz'),
            cont_tables = expand(os.path.join('Calls', '{cell}.contamination.table'), 
                cell=bulk_samples['tumor']),
            rom = os.path.join('Calls', 'read-orientation-model.tar.gz')
        output:
            os.path.join('Calls', 'all.mutect.filtered.vcf.gz')
        envmodules:
            'gatk/4.1.1.0',
        threads: 2
        params:
            ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        shell:
            '{SCRIPT_DIR}/08.4_mutect.sh {input.cont_tables} '
            '-i {input.vcf} -rom {input.rom} -r {params.ref_genome} -o {output}'


def get_final_vcfs(wildcards):
    final_files = []
    if config.get('SCcaller', {}).get('run', False):
        final_files.append(os.path.join('Calls', 'all.sccaller.vcf.gz'))
    if config.get('monovar', {}).get('run', False):
        final_files.append(os.path.join('Calls', 'all.monovar.vcf.gz'))
    if bulk_samples['all']:
        final_files.append(os.path.join('Calls', 'all.mutect.filtered.vcf.gz'))
    return final_files


rule merge_calls:
    input:
        get_final_vcfs
    output:
        all_vcf = os.path.join('Calls', 'all.vcf.gz'),
        chr_vcf = expand(os.path.join('Calls', 'all.{chr}.vcf.gz'), chr=CHROM)
    envmodules:
        'bcftools',
    params:
        out_dir = 'Calls',
    shell:
        '{SCRIPT_DIR}/09_merge_vcfs.sh {input} -o {params.out_dir}'


# ------------------------------------------------------------------------------
# ----------------------------- FINAL FILTERING --------------------------------
# ------------------------------------------------------------------------------


rule filter_calls_chr:
    input:
        os.path.join('Calls', 'all.{chr}.vcf.gz')
    output:
        os.path.join('Calls', 'all_filtered.{chr}.vcf.gz')
    envmodules:
        'cesga/2018',
        'gcccore/6.4.0',
        'pysam/0.16.0.1-python-3.8.1', 
        'pandas/1.0.0-python-3.8.1',
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8192
    params:
        bulk_tumor = bulk_samples['tumor'],
        filter_DP = config.get('filters', {}).get('depth', 10),
        filter_QUAL = config.get('filters', {}).get('qual', 20),
        pref = '-p chr' if config['static']['WGA_ref'].startswith('hg19') else '',
    shell:
        '{SCRIPT_DIR}/10_summarize_vcf.py'


rule merge_filtered_calls:
    input:
        expand(os.path.join('Calls', 'all_filtered.{chr}.vcf.gz'), chr=CHROM)
    output:
        os.path.join('Calls',  'all_filtered.vcf')
    envmodules:  
        'pysam/0.16.0.1-python-3.7.7',
        'pandas/1.0.1-python-3.7.7',
    script:
        f'{SCRIPT_DIR}/11_merge_filtered_vcfs.py'