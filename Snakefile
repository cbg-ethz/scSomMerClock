#!/usr/bin/env python3

import sys
import os

BASE_DIR = config['static_data']['data_path']
samples = os.listdir(os.path.join(BASE_DIR, 'Raw_Data'))

if config['static_data'].get('cellnames', False):
    sample_mapping = pd.read_csv(config['static_data']['cellnames'], sep='\t')
    samples = sample_mapping.loc[samples]
else:
    samples = samples


rule all:
    input:
        expand(
            os.path.join(BASE_DIR, 'Processed', '{sample}.calls.tsv'),
            sample=samples
        )


rule adapter_cutting:
    params:
        ref_genome = os.path.join(
            config['static_data']['resources_path'], config['static_data']['WGA_ref']),
        WGA_lig = os.path.join(
            config['static_data']['resources_path'], config['static_data']['WGA_library]),
    output:
    run:
  

rule allignment:
  

rule remove_duplicates:
  

rule base_recal:
  

rule indel_reallignment:


rule QC_alligned:
  

rule SCCaller:

