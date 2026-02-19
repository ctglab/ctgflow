rule fastqc:
    input:
        os.path.join(
                config["output_folder"],
                "bams",
                "{patient}.{sample_type}.{readgroup}.unaligned.bam",
            ),
    output:
        html=os.path.join(
            config["output_folder"], "qc", "fastqc", "{patient}.{sample_type}.{readgroup}_fastqc.html"),
        fastqc_zip=temp(os.path.join(
            config["output_folder"], "qc", "fastqc", "{patient}.{sample_type}.{readgroup}_fastqc.zip")),
    conda:
        "../envs/qc.yml",
    container:
        config["containers"]["ctgflow_core"],
    log:
        os.path.join(
            config["log_folder"], "stats", "{patient}.{sample_type}.{readgroup}.fastqc.log"),
    shell:
        """
        fastqc {input} 
        """


rule samtools_stats:
    input:
        bam=os.path.join(
                config["output_folder"], "bams", "{patient}.{sample_type}.bam"
            ),
    output:
        stats = os.path.join(
            config["output_folder"], "qc", "samtools", "{patient}.{sample_type}.stats"
        ),
    params:
        threads = config["params"]["samtools"]["threads"],
    conda:
        "../envs/gatk4.yml"
    container:
        config["containers"]["ctgflow_core"]
    log:
        os.path.join(
            config["log_folder"], "stats", "{patient}.{sample_type}.samtools.log"),
    shell:
        """
        samtools stats -@ {params.threads} {input.bam} > {output.stats}
        """
    
rule mosdepth:
    input:
        unpack(get_mosdepth_input),
    output:
        [
            os.path.join(
            config["output_folder"], "qc", "mosdepth", f"{{patient}}.{{sample_type}}.{suf}"
            ) for suf in ["mosdepth.global.dist.txt", "mosdepth.regions.dist.txt", "summary.txt"]
        ],
    params:
        regions = lambda wc, input: f"--by {input.regions}" if "regions" in input else "",
        threads = config["params"]["samtools"]["threads"],
        extra = config["params"]["mosdepth"]["extra"],
        prefix = os.path.join(
            config["output_folder"], "qc", "mosdepth", "{patient}.{sample_type}"
        ),
    conda:
        "../envs/gatk4.yml"
    container: "docker://brentp/mosdepth:v0.3.3"
    log:
        os.path.join(
            config["log_folder"], "stats", "{patient}.{sample_type}.mosdepth.log"),
    shell:
        """
        mosdepth -t {params.threads} {params.regions} {params.extra} {params.prefix} {input.bam}
        """

rule multiqc:
    input:
        unpack(get_multiqc_inputs),
    output:
        os.path.join(config["output_folder"], "qc", "multiqc", "multiqc_report.html"),
    conda:
        "../envs/qc.yml"
    container:
        config["containers"]["ctgflow_core"],
    log:
        os.path.join(config["log_folder"], "stats", "multiqc.log"),
    shell:
        """
        multiqc --filename {output} {input.fastqc} {input.samtools} {input.mosdepth} {input.markdup}
        """

