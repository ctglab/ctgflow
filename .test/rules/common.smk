import os
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.9.1")

configfile: "config/config.yaml"

# Patient/Sample Info
validate(config, schema = "../../workflow/schemas/config.schema.yaml")
units = pd.read_csv(config['units'], dtype=str).set_index(["patient", "sample", "readgroup"], drop=False)
units = units.sort_index()
patients = units.index.get_level_values('patient').unique()
seqtype= config['sequencing_type']

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
    out['tumor'] = f"{config['output_folder']}qc/{wildcards.patient}_tumor_pileupsummaries.table"
    if not tumor_only:
        out['normal'] = f"{config['output_folder']}/{wildcards.patient}_normal_pileupsummaries.table"
    return out

def get_intervals():
    ints = []
    for i in range(num_workers):
        num_zeros = 4 - len(str(i))
        interval = '0' * num_zeros + str(i)
        ints.append(interval)
    return ints
   
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

def get_mosdepth_input(wildcards):
    if config['sequencing_type'] == "WGS":
        # no need for intervals, just bam
        return {
            'bam': f"{config['output_folder']}bams/{wildcards.patient}.{wildcards.sample_type}.bam"
        }
    else:
        return {
            'bam': f"{config['output_folder']}bams/{wildcards.patient}.{wildcards.sample_type}.bam",
            'regions': regions_bed
        }

def get_multiqc_inputs(wildcards):
    """Generate input files for multiqc based on execution mode."""
    ps_index = units.index.droplevel('readgroup').unique()
    inputs = {
        "fastqc": expand(
            os.path.join(
                config["output_folder"],
                "qc",
                "fastqc",
                "{patient}.{sample_type}.{readgroup}_fastqc.html"
            ),
            zip,
            patient=units.index.get_level_values('patient'),
            sample_type=units.index.get_level_values('sample'),
            readgroup=units.index.get_level_values('readgroup'),
        ),
        "samtools": expand(
            os.path.join(
                config["output_folder"],
                "qc",
                "samtools",
                "{patient}.{sample_type}.stats"
            ),
            zip,
            patient=ps_index.get_level_values('patient'),
            sample_type=ps_index.get_level_values('sample'),
        ),
        "mosdepth": expand(
            os.path.join(
                config["output_folder"],
                "qc",
                "mosdepth",
                "{patient}.{sample_type}.summary.txt"
            ),
            zip,
            patient=ps_index.get_level_values('patient'),
            sample_type=ps_index.get_level_values('sample'),
        ),
        "markdup": expand(
            os.path.join(
                config["output_folder"], "qc", "{patient}.{sample_type}.markdup_metrics"
            ),
            zip,
            patient=ps_index.get_level_values('patient'),
            sample_type=ps_index.get_level_values('sample'),
        ),
    }
    return inputs
