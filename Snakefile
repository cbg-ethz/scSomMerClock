#!/usr/bin/env python3

import os


BASE_DIR = workflow.basedir
DATA_DIR = config['specific']['data_path']
RES_PATH = config['static']['resources_path']
SCRIPT_DIR = os.path.join(BASE_DIR, 'scripts')
MODELS = ['clock', 'noClock']
workdir: DATA_DIR

if not os.path.exists('logs'):
    os.mkdir('logs')

CHROM = [i for i in range(1, 23, 1)] + ['X', 'Y']

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

def get_corr_samples(wildcards):
    return [os.path.join('Processing', f'{i}.sorted.bam') \
        for i in cell_map[wildcards.cell]]


def get_all_files(wildcards):
    files = [os.path.join('Calls', 'all.vcf.gz'),
        'PAUP_LRT_test.tsv', 'poisson_dispersion_test.all.tsv']
    
    if config['static'].get('QC_seq', False):
        files.append(os.path.join('QC', 'QC_sequencing.tsv'))

    return files


rule all:
    input:
        get_all_files
            

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
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16384
    threads: 8
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('bwa', ['bwa'])]),
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        pair_end = ' ' if config['specific'].get('pair_end', True) else '-se',
        WGA_lib = config['specific']['WGA_library']
    shell:
        '{params.base_dir}/scripts/02.1_bwa.sh {params.modules} '
        '-s {wildcards.sample} -r {params.ref_genome} -l {params.WGA_lib} '
        '{params.pair_end}'


rule alignment2:
    input:
        os.path.join('Processing', '{sample}.sam')
    output:
        temp(os.path.join('Processing', '{sample}.sorted.bam'))
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16384
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('picard', ['picard'])])
    shell:
        '{params.base_dir}/scripts/02.2_bwa.sh {params.modules} '
        '-s {wildcards.sample}'


rule remove_duplicates:
    input:
        get_corr_samples
    output:
        os.path.join('Processing', '{cell}.dedup.bam')
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32768
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('picard', ['picard'])])
    shell:
        '{params.base_dir}/scripts/03_md_merge_rename.sh {input} '
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
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16384
    threads: 4
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('gatk3', ['gatk/3'])]),
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        pref = 'chr' if config['static']['WGA_ref'].startswith('hg19') else '',
        indels1 = os.path.join(RES_PATH, config['static']['indel_db1']),
        indels2 = '' if not config['static'].get('indel_db2', False) else
            '-i2 {}'.format(os.path.join(RES_PATH, config['static']['indel_db2']))
    shell:
        '{params.base_dir}/scripts/05.1_indel_realign.sh {input} '
        '{params.modules} -c {params.pref}{wildcards.chr} -o {output} '
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
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32768
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('gatk3', ['gatk/3'])]),
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        pref = 'chr' if config['static']['WGA_ref'].startswith('hg19') else '',
        indels1 = os.path.join(RES_PATH, config['static']['indel_db1']),
        indels2 = '' if not config['static'].get('indel_db2', False) else
            '-i2 {}'.format(os.path.join(RES_PATH, config['static']['indel_db2']))
    shell:
        '{params.base_dir}/scripts/05.2_indel_realign.sh {input.bams} '
        '{params.modules} -c {params.pref}{wildcards.chr} -r {params.ref_genome} '
        '-t {input.target} -ma {input.maps} -i1 {params.indels1} '
        '{params.indels2}'


rule base_recal:
    input:
        expand(os.path.join('Processing', '{{cell}}.real.{chr}.bam'),
            chr=CHROM)
    output:
        temp(os.path.join('Processing', '{cell}.real.bam'))
    params:
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('samtools', ['samtools'])]),
    shell:
        'module load samtools && samtools merge {output} {input}'


rule base_recal1:
    input:
        os.path.join('Processing', '{cell}.real.bam')
    output:
        os.path.join('Processing', '{cell}.recal.table')
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32768
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('gatk4', ['gatk/4'])]),
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        dbsnp = os.path.join(RES_PATH, config['static']['dbsnp']),
        indels1 = os.path.join(RES_PATH, config['static']['indel_db1'])
    shell:
        '{params.base_dir}/scripts/04.1_base_recal.sh {params.modules} '
        '-i {input} -o {output} -r {params.ref_genome} -d {params.dbsnp} '
        '-id {params.indels1}'


rule base_recal2:
    input:
        bam = os.path.join('Processing', '{cell}.real.bam'),
        table = os.path.join('Processing', '{cell}.recal.table')
    output:
        os.path.join('Processing', '{cell}.recal.bam')
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32768
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('gatk4', ['gatk/4'])]),
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref'])
    shell:
        '{params.base_dir}/scripts/04.2_base_recal.sh {params.modules} '
        '-i {input.bam} -t {input.table} -o {output} -r {params.ref_genome}'


rule base_recal3:
    input:
        os.path.join('Processing', '{cell}.recal.bam')
    output:
        expand(os.path.join('Processing', '{{cell}}.recal.{chr}.bam'),
            chr=CHROM)
    params:
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('samtools', ['samtools'])]),
        chrom=CHROM
    shell:
        'module purge && module load samtools && '
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
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32768
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('SCcaller', ['pysam', 'numpy'])]),
        bulk = config['specific'].get('bulk_normal', ''),
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        dbsnp = os.path.join(RES_PATH, config['static']['dbsnp']),
        min_depth = config['filters'].get('depth', 10),
        sccaller = config['SCcaller']['exe']
    shell:
        '{params.base_dir}/scripts/06.1_sccallerlab.sh {params.modules} '
        '-s {wildcards.cell} -c {wildcards.chr} -b {params.bulk} '
        '-r {params.ref_genome} -d {params.dbsnp} -e {params.sccaller} '
        '-md {params.min_depth}'


rule SCcaller2:
    input:
        expand(os.path.join('Calls', '{{cell}}.{chr}.sccaller.vcf'),
            chr=CHROM)
    output:
        os.path.join('Calls', '{cell}.sccaller.vcf.gz')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('bcftools', ['bcftools'])])
    shell:
        '{params.base_dir}/scripts/06.2_sccallerlab.sh {input} {params.modules} '
        '-o {output[0]}'


rule SCcaller3:
    input:
        expand(os.path.join('Calls', '{cell}.sccaller.vcf.gz'), cell=ss_samples)
    output:
        os.path.join('Calls', 'all.sccaller.vcf.gz')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('bcftools', ['bcftools'])])
    shell:
        '{params.base_dir}/scripts/06.3_sccallerlab.sh {input} {params.modules} '
        '-o {output[0]}'


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
    threads: 3
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32768
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('monovar', ['monovar', 'samtools'])]),
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        pref = 'chr' if config['static']['WGA_ref'].startswith('hg19') else '',
        monovar = config['monovar']['exe']
    shell:
        '{params.base_dir}/scripts/07.1_monovar.sh {params.modules} '
        '-c {wildcards.chr} -r {params.ref_genome} -p {params.pref} '
        '-e {params.monovar}'


rule monovar2:
    input:
        expand(os.path.join('Calls', '{chr}.monovar.vcf'),
            chr=CHROM)
    output:
        os.path.join('Calls', 'all.monovar.vcf.gz')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('bcftools', ['bcftools'])]),
        min_depth = config['filters'].get('depth', 10),
    shell:
        '{params.base_dir}/scripts/07.2_monovar.sh {input} {params.modules} '
        '-o {output[0]} -md {params.min_depth}'


rule mutect1:
    input: 
        expand(os.path.join('Processing', '{cell}.recal.{{chr}}.bam'), 
            cell=bulk_samples['all'])
    output:
        os.path.join('Calls', '{chr}.mutect.vcf'),
        os.path.join('Calls', '{chr}.f1r2.mutect.tar.gz')
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32768
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('mutect1', ['gatk/4.1', 'picard'])]),
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        germ_res = os.path.join(RES_PATH, config['static']['germline']),
        pon = os.path.join(RES_PATH, config['specific']['PON']),
        normal = f'-n {config["specific"]["bulk_normal"]}'
            if config['specific'].get('bulk_normal', False) else ''
    shell:
        '{params.base_dir}/scripts/08.1_mutect.sh {input} {params.modules} '
        '-c {wildcards.chr} -r {params.ref_genome} -g {params.germ_res} '
        '-p {params.pon} {params.normal}'


rule mutect2:
    input: 
        expand(os.path.join('Calls', '{chr}.mutect.vcf'), chr=CHROM)
    output:
        os.path.join('Calls', 'all.mutect.vcf.gz')
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32768
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('bcftools', ['bcftools'])])
    shell:
        '{params.base_dir}/scripts/08.2_mutect.sh {input} {params.modules} '
        '-o {output[0]}'


if config.get('mutect', {}).get('filter', '') == 'simple':
    rule mutect4_simple:
        input: 
            os.path.join('Calls', 'all.mutect.vcf.gz'),
        output:
            os.path.join('Calls', 'all.mutect.filtered.vcf.gz')
        threads: 2
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 16384
        params:
            base_dir = BASE_DIR,
            modules = ' '.join([f'-m {i}' for i in \
                config['modules'].get('gatk41', ['gatk/4.1'])]),
            ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref']),
        shell:
            '{params.base_dir}/scripts/08.4_mutect_simple.sh {params.modules} '
            '-i {input} -r {params.ref_genome} -o {output[0]}'

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
        threads: 2
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 16384
        params:
            base_dir = BASE_DIR,
            modules = ' '.join([f'-m {i}' for i in \
                config['modules'].get('gatk41', ['gatk/4.1'])]),
            gnomAD = os.path.join(RES_PATH, config['static']['gnomAD'])
        shell:
            '{params.base_dir}/scripts/08.3_mutect.sh {input.tumor} '
            '{params.modules} -n {input.normal} -gAD {params.gnomAD}'


    rule mutect4:
        input: 
            vcf = os.path.join('Calls', 'all.mutect.vcf.gz'),
            cont_tables = expand(os.path.join('Calls', '{cell}.contamination.table'), 
                cell=bulk_samples['tumor']),
            rom = os.path.join('Calls', 'read-orientation-model.tar.gz')
        output:
            os.path.join('Calls', 'all.mutect.filtered.vcf.gz')
        threads: 2
        params:
            base_dir = BASE_DIR,
            modules = ' '.join([f'-m {i}' for i in \
                config['modules'].get('gatk41', ['gatk/4.1'])]),
            ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref'])
        shell:
            '{params.base_dir}/scripts/08.4_mutect.sh {input.cont_tables} '
            '{params.modules} -i {input.vcf} -rom {input.rom} '
            '-r {params.ref_genome} -o {output[0]}'


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
    params:
        base_dir = BASE_DIR,
        modules = ' '.join([f'-m {i}' for i in \
            config['modules'].get('bcftools', ['bcftools'])]),
        out_dir = 'Calls',
    shell:
        '{params.base_dir}/scripts/09_merge_vcfs.sh {input} {params.modules} '
        '-o {params.out_dir}'


# ------------------------------------------------------------------------------
# -------------------------------- CALLING QC ----------------------------------
# ------------------------------------------------------------------------------


rule QC_calling_chr:
    input:
        os.path.join('Calls', 'all.{chr}.vcf.gz')
    output:
        temp(os.path.join('QC', 'all.{chr}.filtered.vcf')),

    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8192
    params:
        base_dir = BASE_DIR,
        modules = ' '.join(config['modules'] \
            .get('QC_calling', ['pysam', 'pandas'])),
        bulk_normal = config['specific'].get('bulk_normal', ''),
        bulk_tumor = ' '.join([' '.join(i) for i in bulk_samples['tumor']]),
        filter_DP = config.get('filters', {}).get('depth', 10),
        filter_QUAL = config.get('filters', {}).get('QUAL', 20),
        pref = '-p chr' if config['static']['WGA_ref'].startswith('hg19') else '',
    shell:
        'module load {params.modules} && '
        'python {params.base_dir}/scripts/10_summarize_vcf.py {input} -o QC '
        '-bn {params.bulk_normal} -bt {params.bulk_tumor} '
        '-q {params.filter_QUAL} -r {params.filter_DP} {params.pref}'


rule QC_calling_all:
    input:
        expand(os.path.join('QC', 'all.{chr}.filtered.vcf'), chr=CHROM)
    output:
        os.path.join('QC',  'all.filtered.vcf')
    envmodules:  
        'pysam/0.16.0.1-python-3.7.7',
        'pandas/1.0.1-python-3.7.7',
    script:
        f'{SCRIPT_DIR}/merge_summary_vcfs.py'


# ------------------------------------------------------------------------------
# ------------------------------ Clock testing ---------------------------------
# ------------------------------------------------------------------------------

# -------------------------- POISSON DISPERSION --------------------------------

rule run_poisson_test:
    input:
        os.path.join('QC',  'all.filtered.vcf')
    output:
        'poisson_dispersion_test.all.tsv'
    envmodules:  
        'scipy/1.4.1-python-3.7.7',
    resources:
        mem_mb = 2048
    params:
        exclude = config['specific'].get('sc_normal', '')
    script:
        f'simulations/{SCRIPT_DIR}/get_poisson_LRT.py'


# ------------------------------- PAUP LRT -------------------------------------


rule run_scite:
    input:
        os.path.join('QC', 'all.filtered.vcf')
    output:
        os.path.join('QC', 'scite_tree_ml0.newick')
    resources:
        mem_mb = 4096
    envmodules:
        'numpy/1.18.1-python-3.7.7'
    params:
        steps = config.get('scite', {}).get('steps', 1E6),
        exe = config.get('scite', {}).get('exe', './infSCITE'),
    script:
        f'simulations/{SCRIPT_DIR}/run_scite.py'


rule vcf_to_nex:
    input:
        vcf = os.path.join('QC',  'all.filtered.vcf'),
        tree = os.path.join('QC', 'scite_tree_ml0.newick'),
    output:
        expand(os.path.join('QC', 'all.nxs.{model}'), model=MODELS)
    resources:
        mem_mb = 1024
    params:
        ss = config.get('mrbayes', {}).get('ss', False),
        paup_exe = config.get('paup', {}).get('exe', 'paup'),
    script:
        f'simulations/{SCRIPT_DIR}/convert_vcf_to_nexus.py'


rule run_PAUP:
    input:
        os.path.join('QC', 'all.nxs.{model}')
    output:
        os.path.join('QC', 'all.paup.{model}.PAUP.score')
    params:
        paup_exe = config.get('paup', {}).get('exe', 'paup'),
    shell:
        '{params.paup_exe} -n {input} > QC/paup.{wildcards.model}.log'


rule merge_paup_results:
    input:
        expand(os.path.join('QC', 'all.paup.{model}.PAUP.score'), model=MODELS)
    output:
        'PAUP_LRT_test.tsv'
    conda: 'envs/monovar.yaml'
    envmodules:
        'scipy/1.4.1-python-3.7.7',
    resources:
        mem_mb = 2048
    params:
        no_cells = config['cellcoal']['model']['no_cells'] + 1,
    script:
        f'{SCRIPT_DIR}/get_PAUP_LRT.py'


# ------------------------------------------------------------------------------
# ----------------------------- ADO Calculation --------------------------------
# ------------------------------------------------------------------------------


rule ADO_calculation:
    input:
        bulk = os.path.join('Calls', 'all.mutect.filtered.vcf.gz'),
        ss = expand(os.path.join('Processing', '{cell}.recal.{chr}.bam'),
            chr=CHROM, cell=ss_samples)
    output:
        os.path.join('QC',  'ADO_rates.txt')
    params:
        base_dir = BASE_DIR,
        modules = ' '.join(config['modules'] \
                .get('QC_calling', ['pysam', 'pandas'])),
        dbsnp = os.path.join(RES_PATH, config['static']['dbsnp']),
        ref_genome = os.path.join(RES_PATH, config['static']['WGA_ref'])
    shell:
        'module load {params.modules} && '
        'python {params.base_dir}/scripts/ADO_calculation.py {input} '
        '--bulk {input.bulk} --dbsnp {params.dbsnp} -r {params.ref_genome} '
        '-o {output}'


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
        expand(os.path.join('Processing', '{cell}.genome.tsv'), cell=ss_samples)
    output:
        'QC_sequencing.tsv'
    params:
        base_dir = BASE_DIR,
        modules = ' '.join( \
                config['modules'].get('QC_seq', ['pandas', 'matplotlib']))
    shell:
        'module load {params.modules} && '
        'python {params.base_dir}/scripts/QC_coverage.py {input} -o "./"'
