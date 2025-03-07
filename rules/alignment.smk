from pathlib import Path

rule combine_fqs:
    input:
        unpack(get_fastq)
    output:
        temp(
            config['output_folder']
            + "bams/{patient}.{sample_type}.{readgroup}.unaligned.bam")
    params:
        pl=get_platform,
        tmp=tmp_dir
    container: config['containers']['ctgflow_core']
    resources:
    log:
        config['log_folder']
        + "combine_fqs"
        + "/"
        + "{patient}.{sample_type}.{readgroup}.log"
    shell:
        """
        gatk FastqToSam -F1 {input.r1} -F2 {input.r2} \
            -O {output} \
            -SM {wildcards.patient}.{wildcards.sample_type} \
            -RG {wildcards.patient}.{wildcards.sample_type}.{wildcards.readgroup} \
            --TMP_DIR {params.tmp} \
            -PL {params.pl}
        """

if config['viral_integrated']:
    rule integrate_genome:
        input:
            ref=config['resources']['reference_fasta'],
            viral=config['resources']['viral_genome'],
        output:
            fasta=config['output_folder']
                + "reference"
                + "/"
                + "viral_integrated.fasta",
        container: config['containers']['ctgflow_core']
        resources:
        params: 
            virus_name=config['virus_name'],
        log:
            config['log_folder']
            + "genome_integration.log"
        shell:
            """
            cat {input.ref} > {output.fasta}
            echo ">{params.virus_name}" >> {output.fasta}
            sed -n '2,$p' {input.viral} >> {output.fasta}
            """
    
rule bwa_index:
    input:
        branch(
            config['viral_integrated'],
            then=config['output_folder']
                + "reference"
                + "/"
                + "viral_integrated.fasta",
            otherwise=config['resources']['reference_fasta']
        ),
    output:
        branch(
            config['viral_integrated'],
            then=[f"{config['output_folder']}reference/viral_integrated.fasta.{suf}" for suf in file_suffixes],
            otherwise=[f"{config['resources']['reference_fasta']}.{suf}" for suf in file_suffixes]
        ),
    container: config['containers']['ctgflow_core']
    resources:
    log:
        config['log_folder']
        + "bwa_index.log"
    shell:
        """
        bwa index {input}
        """

rule bwa:
    input:
        branch(
            config['viral_integrated'],
            then=[
                f"{config['output_folder']}reference/viral_integrated.fasta.{suf}" for suf in file_suffixes
            ],
            otherwise=[
                f"{config['resources']['reference_fasta']}.{suf}" for suf in file_suffixes
            ]
        ),
        bam=config['output_folder']
            + "bams/{patient}.{sample_type}.{readgroup}.unaligned.bam",
        ref_fasta=config['output_folder']
            + "reference"
            + "/"
            + "viral_integrated.fasta",
    output:
        temp(
            config['output_folder']
            + "bams/{patient}.{sample_type}.{readgroup}.unfiltered.bam")
    threads: 8
    container: config['containers']['ctgflow_core']
    resources:
    log:
        config['log_folder']
        + "bwa"
        + "/"
        + "{patient}.{sample_type}.{readgroup}.log"
    shell:
        """
        gatk SamToFastq -I {input.bam} -F /dev/stdout -INTER true -NON_PF true \
        | \
        bwa mem -p -v 3 -t {threads} -T 0 \
            {input.ref_fasta} /dev/stdin - 2> >(tee {log} >&2) \
        | \
        samtools view -Shb -o {output}
        """

if config['viral_integrated']:

    rule sort_unaligned:
        input:
            bam=config['output_folder']
                + "bams/{patient}.{sample_type}.{readgroup}.unfiltered.bam"
        output:
            temp(
                config['output_folder']
                + "bams/{patient}.{sample_type}.{readgroup}.sorted.bam"),
            temp(
                config['output_folder']
                + "bams/{patient}.{sample_type}.{readgroup}.sorted.bam.bai"),
        container: config['containers']['ctgflow_core']
        threads: 4
        resources:
        shell:
            """
            samtools sort --write-index -@ {threads} -o {output[0]} {input.bam}
            """

    rule extract_viral:
        input:
            bam=config['output_folder']
                + "bams/{patient}.{sample_type}.{readgroup}.sorted.bam",
            bai=config['output_folder']
                + "bams/{patient}.{sample_type}.{readgroup}.sorted.bam.bai",
            ref_fasta=config['output_folder']
                + "reference"
                + "/"
                + "viral_integrated.fasta"
        output:
            viral_cram=temp(
                config['output_folder']
                + "bams/{patient}.{sample_type}.{readgroup}.viral.cram"),
            non_viral_bam=temp(
                config['output_folder']
                + "bams/{patient}.{sample_type}.{readgroup}.aligned.bam"),
        params:
            virus_name=config['virus_name'],
        container: config['containers']['ctgflow_core']
        resources:
        log:
            config['log_folder']
            + "extract_viral"
            + "/"
            + "{patient}.{sample_type}.{readgroup}.log"
        shell:
            """
            samtools view -b {input.bam} -C -T {input.ref_fasta} -o {output.viral_cram} {params.virus_name} ;
            samtools view -h {input.bam} | awk '($3!={params.virus_name})' | samtools view -b - > {output.non_viral_bam}
            """
