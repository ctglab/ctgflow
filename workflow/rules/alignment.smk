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
        "envs/gatk4.yml"
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


if config["viral_integrated"]:

    rule integrate_genome:
        input:
            ref=config["resources"]["reference_fasta"],
            viral=config["resources"]["viral_genome"],
        output:
            fasta=os.path.join(
                config["output_folder"], "reference", "viral_integrated.fasta"
            ),
            fasta_dict=os.path.join(
                config["output_folder"], "reference", "viral_integrated.dict"
            )
        conda:
            "envs/gatk4.yml"
        container:
            config["containers"]["ctgflow_core"]
        params:
            virus_name=config["virus_name"],
        log:
            os.path.join(config["log_folder"], "genome_integration.log"),
        shell:
            """
            cat {input.ref} > {output.fasta}
            echo ">{params.virus_name}" >> {output.fasta}
            sed -n '2,$p' {input.viral} >> {output.fasta}
            gatk CreateSequenceDictionary -R {output.fasta}
            """


rule bwa_index:
    input:
        branch(
            config["viral_integrated"],
            then=os.path.join(
                config["output_folder"], "reference", "viral_integrated.fasta"
            ),
            otherwise=os.path.join(config["resources"]["reference_fasta"]),
        ),
    output:
        branch(
            config["viral_integrated"],
            then=[
                os.path.join(
                    config["output_folder"],
                    "reference",
                            f"viral_integrated.fasta.{suf}",
                        )
                        for suf in file_suffixes
                    ],
                    otherwise=[
                f"{config['resources']['reference_fasta']}.{suf}"
                for suf in file_suffixes
            ],
        ),
    conda:
        "envs/gatk4.yml"
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
        branch(
            config["viral_integrated"],
            then=[
                os.path.join(
                    config["output_folder"],
                    "reference",
                    f"viral_integrated.fasta.{suf}")
                    for suf in file_suffixes
                    ],
                otherwise=[
                f"{config['resources']['reference_fasta']}.{suf}"
                for suf in file_suffixes
            ],
        ),
        bam=os.path.join(
            config["output_folder"],
            "bams",
            "{patient}.{sample_type}.{readgroup}.unaligned.bam",
        ),
        ref_fasta=branch(
            config["viral_integrated"],
            then=os.path.join(
                config["output_folder"], "reference", "viral_integrated.fasta"
            ),
            otherwise=config["resources"]["reference_fasta"],
        ),
    output:
        temp(
            os.path.join(
                config["output_folder"],
                "bams",
                "{patient}.{sample_type}.{readgroup}.unfiltered.bam",
            )
        ),
    threads: 8
    conda:
        "envs/gatk4.yml"
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


if config["viral_integrated"]:

    rule sort_unaligned:
        input:
            bam=os.path.join(
                config["output_folder"],
                "bams",
                "{patient}.{sample_type}.{readgroup}.unfiltered.bam",
            ),
        output:
            bam=temp(
                os.path.join(
                    config["output_folder"],
                    "bams",
                    "{patient}.{sample_type}.{readgroup}.sorted.bam",
                )
            ),
            bai=temp(
                os.path.join(
                    config["output_folder"],
                    "bams",
                    "{patient}.{sample_type}.{readgroup}.sorted.bam.bai",
                )
            ),
        conda:
            "envs/gatk4.yml"
        container:
            config["containers"]["ctgflow_core"]
        log:
            os.path.join(
                config["log_folder"],
                "sort_unfiltered",
                "{patient}.{sample_type}.{readgroup}.log"
            ),
        threads: 4
        shell:
            """
            samtools sort -@ {threads} -o {output.bam} {input.bam}
            samtools index {output.bam}
            """

    rule extract_viral:
        input:
            bam=os.path.join(
                config["output_folder"],
                "bams",
                "{patient}.{sample_type}.{readgroup}.sorted.bam",
            ),
            bai=os.path.join(
                config["output_folder"],
                "bams",
                "{patient}.{sample_type}.{readgroup}.sorted.bam.bai",
            ),
            ref_fasta=os.path.join(
                config["output_folder"], "reference", "viral_integrated.fasta"
            ),
        output:
            viral_cram=temp(
                os.path.join(
                    config["output_folder"],
                    "bams",
                    "{patient}.{sample_type}.{readgroup}.viral.cram",
                )
            ),
            non_viral_bam=temp(
                os.path.join(
                    config["output_folder"],
                    "bams",
                    "{patient}.{sample_type}.{readgroup}.non_viral.bam",
                )
            ),
        params:
            virus_name=config["virus_name"],
        conda:
            "envs/gatk4.yml"
        container:
            config["containers"]["ctgflow_core"]
        log:
            os.path.join(
                config["log_folder"],
                "extract_viral",
                "{patient}.{sample_type}.{readgroup}.log"
            ),
        shell:
            """
            samtools view -b {input.bam} -C -T {input.ref_fasta} -o {output.viral_cram} {params.virus_name} ;
            samtools view -h {input.bam} | grep -v '{params.virus_name}' | samtools view -b - > {output.non_viral_bam}
            samtools index {output.non_viral_bam}
            """
