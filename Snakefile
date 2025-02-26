singularity: "docker://continuumio/miniconda3"
configfile: "config/config_main.yaml"
include: "rules/common.smk"

rule targets:
    input:
        vcf=expand(
            config['output_folder']
            + "vcfs/{patient}.vep.vcf.gz",
            patient=patients
            ),
        idx=expand(
            config['output_folder']
            + "vcfs/{patient}.vep.vcf.gz.tbi",
            patient=patients
            )

include: "rules/preprocessing.smk"
include: "rules/somatic.smk"
