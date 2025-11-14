from pathlib import Path
import os


rule combine_fqs:
    input:
        r1=lookup(
            query=lambda w: f"patient == '{w.patient}' and sample == '{w.sample_type}' and readgroup == '{w.readgroup}'",
            cols=["fq1"],
            within=units
        ),
        r2=lookup(
            query=lambda w: f"patient == '{w.patient}' and sample == '{w.sample_type}' and readgroup == '{w.readgroup}'",
            cols=["fq1"],
            within=units
        ),
    output:
        temp(
            os.path.join(
                config["output_folder"],
                "bams",
                "{patient}.{sample_type}.{readgroup}.unaligned.bam",
            )
        ),
    params:
        pl=get_platform,
        tmp=config['tmp_dir'],
    conda:
        "../envs/gatk4.yml"
    container:
        config["containers"]["ctgflow_core"]
    log:
        os.path.join(
            config["log_folder"],
            "combine_fqs",
            "{patient}.{sample_type}.{readgroup}.log",
        ),
    shell:
        """
        gatk FastqToSam -F1 {input.r1} -F2 {input.r2} \
            -O {output} \
            -SM {wildcards.patient}.{wildcards.sample_type} \
            -RG {wildcards.patient}.{wildcards.sample_type}.{wildcards.readgroup} \
            --TMP_DIR {params.tmp} \
            -PL {params.pl}
        """



rule bwa_index:
    input:
        os.path.join(config["resources"]["reference_fasta"]),
    output:
        [
                f"{config['resources']['reference_fasta']}.{suf}"
                for suf in file_suffixes
            ],
    conda:
        "../envs/gatk4.yml"
    container:
        config["containers"]["ctgflow_core"]
    log:
        os.path.join(config["log_folder"], "bwa_index.log"),
    shell:
        """
        bwa index {input}
        """


rule bwa:
    input:
        [
                f"{config['resources']['reference_fasta']}.{suf}"
                for suf in file_suffixes
            ],
        bam=os.path.join(
            config["output_folder"],
            "bams",
            "{patient}.{sample_type}.{readgroup}.unaligned.bam",
        ),
        ref_fasta=config["resources"]["reference_fasta"],
    output:
        temp(
            os.path.join(
                config["output_folder"],
                "bams",
                "{patient}.{sample_type}.{readgroup}.aligned.bam",
            )
        ),
    threads: 8
    conda:
        "../envs/gatk4.yml"
    container:
        config["containers"]["ctgflow_core"]
    log:
        os.path.join(
            config["log_folder"], "bwa", "{patient}.{sample_type}.{readgroup}.log"
        ),
    shell:
        """
        gatk SamToFastq -I {input.bam} -F /dev/stdout -INTER true -NON_PF true \
        | \
        bwa mem -p -v 3 -t {threads} -T 0 \
            {input.ref_fasta} /dev/stdin - 2> >(tee {log} >&2) \
        | \
        samtools view -Shb -o {output}
        """