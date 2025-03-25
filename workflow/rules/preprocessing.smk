rule merge_bams:
    input:
        unaligned=os.path.join(
            config["output_folder"],
            "bams",
            "{patient}.{sample_type}.{readgroup}.unaligned.bam"
        ),
        aligned=branch(
            config["viral_integrated"],
            then=os.path.join(
                config["output_folder"],
                "bams",
                "{patient}.{sample_type}.{readgroup}.non_viral.bam",
            ),
            otherwise=os.path.join(
                config["output_folder"],
                "bams",
                "{patient}.{sample_type}.{readgroup}.aligned.bam",
            ),
        ),
        ref=config["resources"]["reference_fasta"],
    output:
        temp(
            os.path.join(
                config["output_folder"],
                "bams",
                "{patient}.{sample_type}.{readgroup}.merged.bam",
            ),
        ),
    params:
        tmp=config["tmp_dir"],
    conda:
        "../envs/gatk4.yml"
    container:
        config["containers"]["ctgflow_core"]
    log:
        os.path.join(
            config["log_folder"],
            "merge_bams",
            "{patient}.{sample_type}.{readgroup}.log",
        ),
    shell:
        """
        gatk MergeBamAlignment \
            -R {input.ref} \
            -UNMAPPED {input.unaligned} \
            -ALIGNED {input.aligned} \
            -O {output} \
            --SORT_ORDER queryname \
            --TMP_DIR {params.tmp}
        """


rule markdups_sort:
    input:
        lambda wc: [
            os.path.join(
                config["output_folder"],
                "bams",
                f"{wc.patient}.{wc.sample_type}.{rg}.merged.bam",
            )
            for rg in units.loc[(wc.patient, wc.sample_type), "readgroup"].unique()
        ],
    output:
        bam=temp(
            os.path.join(
                config["output_folder"], "bams", "{patient}.{sample_type}.sorted.bam"
            )
        ),
        metrics=os.path.join(
            config["output_folder"], "qc", "{patient}.{sample_type}.markdup_metrics"
        ),
    params:
        inbams=lambda wildcards, input: " -I  ".join(input),
        tmp=config["tmp_dir"],
    threads: 4
    conda:
        "../envs/gatk4.yml"
    container:
        config["containers"]["ctgflow_core"]
    log:
        os.path.join(
            config["log_folder"], "markdups_sort", "{patient}.{sample_type}.log"
        ),
    shell:
        """
        gatk --java-options "-XX:ParallelGCThreads={threads}" MarkDuplicates \
            -I {params.inbams} \
            -O {output.bam} \
            -M {output.metrics} \
            --TMP_DIR {params.tmp} 
        """


rule sort_mrkdups:
    input:
        bam=os.path.join(
            config["output_folder"], "bams", "{patient}.{sample_type}.sorted.bam"
        ),
        ref=config["resources"]["reference_fasta"],
    output:
        bam=temp(
            os.path.join(
                config["output_folder"],
                "bams",
                "{patient}.{sample_type}.sorted.markdup.bam",
            )
        ),
        bai=temp(
            os.path.join(
                config["output_folder"],
                "bams",
                "{patient}.{sample_type}.sorted.markdup.bam.bai",
            )
        ),
    conda:
        "../envs/gatk4.yml"
    container:
        config["containers"]["ctgflow_core"]
    log:
        os.path.join(
            config["log_folder"], "sort_mrkdups", "{patient}.{sample_type}.log"
        ),
    shell:
        """
        samtools sort -O bam -o {output.bam} {input.bam} ;
        samtools index {output.bam}
        """


# @TODO: scatter on multiple intervals like Mutect2 follwing
# https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/PairedEndSingleSampleWf-fc-hg38.wdl
rule bqsr:
    input:
        bam=os.path.join(
            config["output_folder"],
            "bams",
            "{patient}.{sample_type}.sorted.markdup.bam",
        ),
        bai=os.path.join(
            config["output_folder"],
            "bams",
            "{patient}.{sample_type}.sorted.markdup.bam.bai",
        ),
        interval=os.path.join(
            config["output_folder"],
            "interval-files",
            "{interval}-scattered.interval_list",
        ),
        ref=branch(
            config["viral_integrated"],
            then=os.path.join(
                config["output_folder"], "reference", "viral_integrated.fasta"
            ),
            otherwise=config["resources"]["reference_fasta"],
        ),
    output:
        recal=temp(
            os.path.join(
                config["output_folder"],
                "qc",
                "{patient}.{sample_type}.{interval}.recal_data.table",
            )
        ),
    params:
        ks=lambda wildcards, input: " ".join(
            [
                f"--known-sites {config['resources'][ks]}"
                for ks in ["dbsnps", "known-indels", "1000g"]
            ]
        ),
        tmp=config["tmp_dir"],
        ram=config["params"]["gatk"]["RAM"],
    conda:
        "../envs/gatk4.yml"
    container:
        config["containers"]["ctgflow_core"]
    log:
        os.path.join(
            config["log_folder"], "bqsr", "{patient}.{sample_type}.{interval}.log"
        ),
    shell:
        """
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 \
        -Xms{params.ram}m" \
        BaseRecalibrator \
        -I {input.bam} -R {input.ref} \
        --use-original-qualities \
        -O {output.recal} \
        -L {input.interval} \
        {params.ks} --tmp-dir {params.tmp}
        """


rule GatherBQSRReports:
    input:
        lambda wc: [
            os.path.join(
                config["output_folder"],
                "qc",
                f"{wc.patient}.{wc.sample_type}.{interval}.recal_data.table",
            )
            for interval in get_intervals()
        ],
    output:
        os.path.join(
            config["output_folder"], "qc", "{patient}.{sample_type}.recal_data.table"
        ),
    conda:
        "../envs/gatk4.yml"
    container:
        config["containers"]["ctgflow_core"]
    log:
        os.path.join(
            config["log_folder"], "gather_bqsr_reports", "{patient}.{sample_type}.log"
        ),
    shell:
        """
        gatk --java-options "-Xms3000m" \
            GatherBQSRReports \
            {input} \
            -O {output}
        """


rule apply_bqsr:
    input:
        bam=os.path.join(
            config["output_folder"],
            "bams",
            "{patient}.{sample_type}.sorted.markdup.bam",
        ),
        bai=os.path.join(
            config["output_folder"],
            "bams",
            "{patient}.{sample_type}.sorted.markdup.bam.bai",
        ),
        recal=os.path.join(
            config["output_folder"],
            "qc",
            "{patient}.{sample_type}.{interval}.recal_data.table",
        ),
        interval=os.path.join(
            config["output_folder"],
            "interval-files",
            "{interval}-scattered.interval_list",
        ),
        ref=config["resources"]["reference_fasta"],
    output:
        bam=temp(
            os.path.join(
                config["output_folder"],
                "bams",
                "{patient}.{sample_type}.{interval}.bam",
            )
        ),
        bai=temp(
            os.path.join(
                config["output_folder"],
                "bams",
                "{patient}.{sample_type}.{interval}.bai",
            )
        ),
    params:
        tmp=config["tmp_dir"],
    conda:
        "../envs/gatk4.yml"
    container:
        config["containers"]["ctgflow_core"]
    log:
        os.path.join(
            config["log_folder"],
            "apply_bqsr",
            "{patient}.{sample_type}.{interval}.log",
        ),
    shell:
        """
        gatk --java-options "-Xms3000m -Xmx10G" \
            ApplyBQSR \
            -I {input.bam} \
            -R {input.ref} \
            -O {output.bam} \
            -L {input.interval} \
            --use-original-qualities \
            -bqsr {input.recal} \
            --add-output-sam-program-record \
            --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
            --tmp-dir {params.tmp}
        """


rule GatherSortedBam:
    input:
        lambda wc: [
            os.path.join(
                config["output_folder"],
                "bams",
                f"{wc.patient}.{wc.sample_type}.{interval}.bam",
            )
            for interval in get_intervals()
        ],
    output:
        bam=temp(
            os.path.join(
                config["output_folder"], "bams", "{patient}.{sample_type}.bam"
            )
        ),
    params:
        tmp=config["tmp_dir"],
        bams=lambda wildcards, input: " ".join([f"-I {f}" for f in input]),
    conda:
        "../envs/gatk4.yml"
    container:
        config["containers"]["ctgflow_core"]
    log:
        os.path.join(config["log_folder"], "gather_bam", "{patient}.{sample_type}.log"),
    shell:
        """
        gatk --java-options "-Xms2000m -Xmx2500m" \
            GatherBamFiles \
            {params.bams} \
            -O {output.bam}
        """


rule sortGather:
    input:
        bam=os.path.join(config["output_folder"], "bams", "{patient}.{sample_type}.bam"),
        ref=config["resources"]["reference_fasta"],
    output:
        cram=os.path.join(
            config["output_folder"], "bams", "{patient}.{sample_type}.cram"
        ),
        crai=os.path.join(
            config["output_folder"], "bams", "{patient}.{sample_type}.cram.crai"
        ),
    conda:
        "../envs/gatk4.yml"
    container:
        config["containers"]["ctgflow_core"]
    log:
        os.path.join(config["log_folder"], "sortGather", "{patient}.{sample_type}.log"),
    shell:
        """
        samtools sort -O cram \
        -T {wildcards.patient}.{wildcards.sample_type} \
        -o {output.cram} {input.bam};
        samtools index {output.cram}
        """

rule split_intervals:
    input:
        ref=config["resources"]["reference_fasta"],
        intervals=regions_gatk,
    output:
        [
            os.path.join(
                config["output_folder"],
                "interval-files",
                f"{i:04d}-scattered.interval_list",
            )
            for i in range(num_workers)
        ],
    params:
        N=num_workers,
        d=lambda w, input: os.path.join(
            Path(input.intervals).parents[0]
        ),
    conda:
        "../envs/gatk4.yml"
    container:
        config["containers"]["ctgflow_core"]
    log:
        os.path.join(config["log_folder"], "split_intervals", "split_intervals.log"),
    shell:
        """
        gatk SplitIntervals -R {input.ref} -L {input.intervals} \
            --scatter-count {params.N} -O {params.d} \
            --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION
        """


rule make_gatk_regions:
    input:
        bed=regions_bed,
        d=lambda x: config["resources"]["reference_fasta"].replace(".fasta", ".dict"),
    output:
        intlist=regions_gatk,
    conda:
        "../envs/gatk4.yml"
    container:
        config["containers"]["ctgflow_core"]
    log:
        os.path.join(config["log_folder"], "make_gatk_regions", "make_gatk_regions.log"),
    shell:
        """
        gatk BedToIntervalList \
            -I {input.bed} \
            -SD {input.d} \
            -O {output}
        """
