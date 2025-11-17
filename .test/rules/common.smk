# File: common.smk
# Author: Tomas Bencomo
# Email: tjbencomo@gmail.com
# Description:
# Initialization file that loads config information and
# stores functions to aide with wildcard creation

import sys
import os
from pathlib import Path
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.9.1")

configfile: "config/config.yaml"

# Patient/Sample Info
# validate(config, schema = "../schemas/config.schema.yaml")
patients = pd.read_csv(config['patients'])['patient']
units = pd.read_csv(config['units'], dtype=str).set_index(["patient", "sample", "readgroup"], drop=False)
units = units.sort_index()
seqtype= config['sequencing_type']

# Check that there are no duplicate patients
# Duplicate patients cause patient MAFs to be repeated as input to concat_mafs and duplicate SNVs
if patients.shape[0] != patients.unique().shape[0]:
    raise ValueError(f"Duplicate patients present in {config['patients']}. Remove duplicates!")

# # Logs
# slurm_logdir = config['slurm_log_dir']
# logpath = Path(slurm_logdir)
# logpath.mkdir(parents=True, exist_ok=True) 

# Reference Files
ref_fasta = config['resources']['reference_fasta']
file_suffixes= ['amb', 'ann', 'bwt', 'pac', 'sa']
# References
# ref_dir = config['ref_dir']
# ref_fasta = os.path.join(ref_dir, config['ref_fasta'])
# ref_dict = ref_fasta.replace('.fasta', '.dict') #used for BedToIntervals
# known_sites = config['known_sites'].replace(' ', '').split(',')
# known_sites = [os.path.join(ref_dir, s) for s in known_sites]
# germline_resource = config['germline_resource']
# contamination_resource = config['contamination_resource']


# Mutect2
num_workers = config['params']['mutect2']['num_workers']
# mutect_flags = config['mutect2_flags']
# filter_flags = config['filter_flags']
tumor_only = config['tumor_only']
# stringent_filtering = config['stringent_filtering']

if tumor_only:
    sample_types = ['tumor']
else:
    sample_types = ['normal', 'tumor']


# TMP Directory
tmp_dir = config['tmp_dir']
if tmp_dir == 'None':
    tmp_dir = 'null'

# Genomic Regions
regions_bed = config['resources']['genomic_regions']
regions_gatk = os.path.basename(regions_bed).replace('.bed', '.interval_list')
regions_gatk = os.path.join(config['output_folder'], 'interval-files', regions_gatk)


wildcard_constraints:
    patient="|".join(patients),
    sample_type="|".join(sample_types)


def get_deepsomatic_input(wildcards):
    files = {}
    files['cram'] = os.path.join(
        config['output_folder'],
        "bams",
        f"{wildcards.patient}.tumor.cram")
    files['crai'] = os.path.join(
        config['output_folder'],
        "bams",
        f"{wildcards.patient}.tumor.cram.crai")
    return files
    
def get_dedup_input(wildcards):
    rgs = units.loc[(wildcards.patient, wildcards.sample_type), 'readgroup'].unique().tolist()
    bams = [f"{config['output_folder']}bams/{wildcards.patient}.{wildcards.sample_type}.{rg}.merged.bam" for rg in rgs]
    return bams

def get_platform(wildcards):
    return units.loc[(wildcards.patient, wildcards.sample_type, wildcards.readgroup), 'platform']

def get_mutect2_input(wildcards):
    files = {}
    files['tumor'] = f"{config['output_folder']}bams/{wildcards.patient}.tumor.cram"
    if not tumor_only:
        files['normal'] = f"{config['output_folder']}bams/{wildcards.patient}.normal.cram"
    return files
    # if use_pon:
    #     return {
    #         'normal' : f"bams/{wildcards.patient}.normal.bam",
    #         'tumor' : f"bams/{wildcards.patient}.tumor.bam",
    #         'pon' : pon_vcf
    #     }
    # else:
    #     return {
    #             'normal' : f"bams/{wildcards.patient}.normal.bam",
    #             'tumor' : f"bams/{wildcards.patient}.tumor.bam",
    #     }



def get_contamination_input(wildcards):
    out = {}
    out['tumor'] = f'{config['output_folder']}qc/{wildcards.patient}_tumor_pileupsummaries.table'
    if not tumor_only:
        out['normal'] = f'{config['output_folder']}/{wildcards.patient}_normal_pileupsummaries.table'
    return out

def get_coverage_input(wildcards):
    files = {}
    files['bam'] = f"{config['output_folder']}bams/{wildcards.patient}.{wildcards.sample_type}.bam"
    if seqtype == "WES":
        files['regions'] = regions_bed
    return files

def isWGS(wildcards):
    seqtype = units.loc[(wildcards.patient, wildcards.sample_type), 'seqtype'][0]
    return seqtype == "WGS"


def get_intervals():
    ints = []
    for i in range(num_workers):
        num_zeros = 4 - len(str(i))
        interval = '0' * num_zeros + str(i)
        ints.append(interval)
    return ints
   
def get_interval_files():
    ints = get_intervals()
    files = [f'{i}-scattered.interval_list' for i in ints]
    files = [os.path.join(config['output_folder'],"interval-files", f) for f in files]
    return files

def get_orientationbias_input(wildcards):
    intervals = get_intervals()
    files = [f"{config['output_folder']}vcfs/{wildcards.patient}.{i}.f1r2.tar.gz" for i in intervals]
    return files

def get_mergevcfs_input(wildcards):
    intervals = get_intervals()
    files = [f"{config['output_folder']}vcfs/{wildcards.patient}.{i}.unfiltered.vcf" for i in intervals]
    return files

def get_mergestats_input(wildcards):
    intervals = get_intervals()
    files = [f"{config['output_folder']}vcfs/{wildcards.patient}.{i}.unfiltered.vcf.stats" for i in intervals]
    return files

def get_gather_bam_input(wildcards):
    intervals = get_intervals()
    files = [f"{config['output_folder']}bams/{wildcards.patient}.{wildcards.sample_type}.{i}.bam" for i in intervals]
    return files

def get_bqsr_output(wildcards):
    intervals = get_intervals()
    recal = [f"{config['output_folder']}qc/{wildcards.patient}.{wildcards.sample_type}.{i}.recal_data.table" for i in intervals]
    return recal


def get_gather_bqsr_reports(wildcards):
    intervals = get_intervals()
    recal=[config['output_folder']
            +"-I qc/{patient}.{sample_type}.{interval}.recal_data.table" for interval in intervals]
    return recal

interval_files = get_interval_files()
