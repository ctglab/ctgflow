
rule merge_bams:
    input:
        unaligned=config['output_folder']
            +"bams/{patient}.{sample_type}.{readgroup}.unaligned.bam",
        aligned=config['output_folder']
            +"bams/{patient}.{sample_type}.{readgroup}.aligned.bam",
        ref=ref_fasta
    output:
        temp(
            config['output_folder']
            + "bams/{patient}.{sample_type}.{readgroup}.merged.bam")
    params:
        tmp=tmp_dir
    container: config['containers']['ctgflow_core']
    resources:
    log:
        config['log_folder']
        + "merge_bams"
        + "/"
        + "{patient}.{sample_type}.{readgroup}.log"
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
        get_dedup_input
    output:
        bam=temp(
            config['output_folder']
            + "bams/{patient}.{sample_type}.sorted.bam"),
        metrics=config['output_folder']
            + "qc/{patient}.{sample_type}.markdup_metrics",
    params:
        inbams=lambda wildcards, input: " -I  ".join(input),
        tmp=tmp_dir
    threads: 4
    container: config['containers']['ctgflow_core']
    resources:
    log:
        config['log_folder']
        + "markdups_sort"
        + "/"
        + "{patient}.{sample_type}.log"
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
        bam=config['output_folder']
            + "bams/{patient}.{sample_type}.sorted.bam",
        ref=ref_fasta
    output:
        bam=temp(
            config['output_folder']
            + "bams/{patient}.{sample_type}.sorted.markdup.bam"),
        bai=config['output_folder']
            + "bams/{patient}.{sample_type}.sorted.markdup.bam.bai"
    container: config['containers']['ctgflow_core']
    resources:
    log:
        config['log_folder']
        + "sort_mrkdups"
        + "/"
        + "{patient}.{sample_type}.log"
    shell:
        """
        samtools sort -@ 10 -m 2G -O bam -o {output.bam} {input.bam} \
        && samtools index {output.bam}
        """


#@TODO: scatter on multiple intervals like Mutect2 follwing 
# https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/PairedEndSingleSampleWf-fc-hg38.wdl
rule bqsr:
    input:
        bam=config['output_folder']
            +"bams/{patient}.{sample_type}.sorted.markdup.bam",
        interval=config["output_folder"]
            +"interval-files/{interval}-scattered.interval_list",
        ref=ref_fasta
    output:
        recal=temp(
            config['output_folder']
            +"qc/{patient}.{sample_type}.{interval}.recal_data.table"
        )
    params:
        ks=lambda wildcards, input: " ".join(
        [ f"--known-sites {config['resources'][ks]}" for ks in ['dbsnps', 'known-indels', '1000g']]),
        tmp=tmp_dir
    container: config['containers']['ctgflow_core']
    resources:
    log:
        config['log_folder']
        + "bqsr"
        + "/"
        + "{patient}.{sample_type}.{interval}.log"
    shell:
        """
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms4000m" \
        BaseRecalibrator \
        -I {input.bam} -R {input.ref} \
        --use-original-qualities \
        -O {output.recal} \
        -L {input.interval} \
        {params.ks} --tmp-dir {params.tmp}
        """

rule GatherBQSRReports:
    input:
        get_gather_bqsr_reports
    output:
        config['output_folder']
            +"qc/{patient}.{sample_type}.recal_data.table"
    container: config['containers']['ctgflow_core']
    resources:
    log:
        config['log_folder']
        + "gather_bqsr_reports"
        + "/"
        + "{patient}.{sample_type}.log"
    shell:
        """
        gatk --java-options "-Xms3000m" \
            GatherBQSRReports \
            {input} \
            -O {output}
        """

rule apply_bqsr:
    input:
        bam=config['output_folder']
            +"bams/{patient}.{sample_type}.sorted.markdup.bam",
        recal=config['output_folder']
            +"qc/{patient}.{sample_type}.{interval}.recal_data.table",
        interval=config['output_folder'] +
        "interval-files/{interval}-scattered.interval_list",
        ref=ref_fasta
    output:
        bam=temp(
            config['output_folder']
            +"bams/{patient}.{sample_type}.{interval}.bam"),
        bai=config['output_folder']
            +"bams/{patient}.{sample_type}.{interval}.bai",
    params:
        tmp=tmp_dir
    container: config['containers']['ctgflow_core']
    resources:
    log:
        config['log_folder']
        + "apply_bqsr"
        + "/"
        + "{patient}.{sample_type}.{interval}.log"
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
        get_gather_bam_input
    output:
        bam=temp(config['output_folder']
            +"bams/{patient}.{sample_type}.bam"),
    params:
        tmp=tmp_dir,
        bams=lambda wildcards, input: " ".join([f"-I {f}" for f in input])
    container: config['containers']['ctgflow_core']
    resources:
    log:
        config['log_folder']
        + "gather_bam"
        + "/"
        + "{patient}.{sample_type}.log"
    shell:
        """
        gatk --java-options "-Xms2000m -Xmx2500m" \
            GatherBamFiles \
            {params.bams} \
            -O {output.bam}
        """

rule sortGather:
    input:
        bam=config['output_folder']
            +"bams/{patient}.{sample_type}.bam",
        ref=ref_fasta,
    output:
        cram=config['output_folder']
            +"bams"
            + "/"
            + "{patient}.{sample_type}.cram",
        crai=config['output_folder']
            +"bams/{patient}.{sample_type}.cram.crai"
    params:
    container: config['containers']['ctgflow_core']
    resources:
    log:
        config['log_folder']
        + "sortGather"
        + "/"
        + "{patient}.{sample_type}.log"
    shell:
        """
        samtools sort -@ 10 -m 2G -O cram \
        -T {wildcards.patient}.{wildcards.sample_type} \
        -o {output.cram} {input.bam};
        samtools index {output.cram}
        """


rule coverage:
    input:
        bam=config['output_folder']
            + "bams"
            + "/"
            + "{patient}.{sample_type}.cram",
        regions=config['resources']['genomic_regions'],
    output:
        config['output_folder']
            + "qc/{patient}.{sample_type}.mosdepth.region.dist.txt",
        config['output_folder']
            + "qc/{patient}.{sample_type}.regions.bed.gz",
        config['output_folder']
            + "qc/{patient}.{sample_type}.mosdepth.global.dist.txt",
        config['output_folder']
            + "qc/{patient}.{sample_type}.mosdepth.summary.txt"
    threads: 4
    params:
        by=lambda wildcards, input: '500' if seqtype == 'WGS' else input.regions
    container: config['containers']['ctgflow_core']
    resources:
    log:
        config['log_folder']
        + "coverage"
        + "/"
        + "{patient}.{sample_type}.log"
    shell:
        """
        mosdepth --by {params.by} -t {threads} \
        qc/{wildcards.patient}.{wildcards.sample_type} \
        {input.bam}
        """

rule stats:
    input:
        config['output_folder']
            + "bams/{patient}.{sample_type}.bam"
    output:
        config['output_folder']
            + "qc/{patient}.{sample_type}.flagstat"
    container: config['containers']['ctgflow_core']
    resources:
    log:
        config['log_folder']
        + "samtools_stats"
        + "/"
        + "{patient}.{sample_type}.log"
    shell:
        """
        samtools flagstat {input} > {output}
        """

rule fastqc:
    input:
        config['output_folder']
            + "bams/{patient}.{sample_type}.bam"
    output:
        html=config['output_folder']
            + "qc/fastqc/{patient}.{sample_type}_fastqc.html",
        zipdata=config['output_folder']
            + "qc/fastqc/{patient}.{sample_type}_fastqc.zip"
    container: config['containers']['ctgflow_core']
    resources:
    log:
        config['log_folder']
        + "fastqc"
        + "/"
        + "{patient}.{sample_type}.log"
    shell:
        """
        tmpdir=qc/fastqc/{wildcards.patient}-{wildcards.sample_type}.tmp 
        mkdir $tmpdir 
        fastqc --outdir $tmpdir {input} 
        mv $tmpdir/{wildcards.patient}.{wildcards.sample_type}_fastqc.html {output.html} 
        mv $tmpdir/{wildcards.patient}.{wildcards.sample_type}_fastqc.zip {output.zipdata} 
        rm -r $tmpdir
        """

rule multiqc:
    input:
        expand(
            config['output_folder']
            + "qc/fastqc/{patient}.{sample_type}_fastqc.zip", patient=patients, sample_type=sample_types),
        expand(
            config['output_folder']
            + "qc/{patient}.{sample_type}.mosdepth.region.dist.txt", patient=patients, sample_type=sample_types),
        expand(
            config['output_folder']
            + "qc/{patient}.{sample_type}.flagstat", patient=patients, sample_type=sample_types)
    output:
        config['output_folder']
            + "qc/multiqc_report.html"
    log:
        config['output_folder']
            + "logs/multiqc.log"
    container: config['containers']['ctgflow_core']
    resources:
    log:
        config['log_folder']
        + "multiqc"
        + "/"
        + "{patient}.{sample_type}.log"
    wrapper:
        "0.50.4/bio/multiqc"

rule seq_depths:
    input:
        expand(
            config['output_folder']
            + "qc/{patient}.{sample_type}.mosdepth.summary.txt", patient=patients, sample_type=sample_types)
    output:
        config['output_folder']
            + "qc/depths.csv"
    container: config['containers']['ctgflow_core']
    resources:
    log:
        config['log_folder']
        + "seq_depths"
        + "/"
        + "depths_aggregate.log"
    script:
        "../scripts/gather_depths.py"

rule plot_depths:
    input:
        config['output_folder']
            + "qc/depths.csv"
    output:
        config['output_folder']
            + "qc/depths.svg"
    container: config['containers']['ctgflow_core']
    resources:
    log:
        config['log_folder']
        + "plot_depths"
        + "/"
        + "plot_depths.log"
    script:
        "../scripts/plot_depth.R"

rule split_intervals:
    input:
        ref=ref_fasta,
        intervals=regions_gatk
    output:
        interval_files
    params:
        N=num_workers,
        d=lambda w, input: os.path.join(
            Path(input.intervals).parents[0], "split_intervals")
    container: config['containers']['ctgflow_core']
    resources:
    log:
        config['log_folder']
        + "split_intervals"
        + "/"
        + "split_intervals.log"
    shell:
        """
        gatk SplitIntervals -R {input.ref} -L {input.intervals} \
            --scatter-count {params.N} -O {params.d} \
            --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION
        """
rule make_gatk_regions:
    input:
        bed=regions_bed,
        d=lambda x: config['resources']['reference_fasta'].replace('.fasta', '.dict')
    output:
        intlist=regions_gatk
    container: config['containers']['ctgflow_core']
    resources:
    log:
        config['log_folder']
        + "make_gatk_regions"
        + "/"
        + "make_gatk_regions.log"
    shell:
        """
        gatk BedToIntervalList \
            -I {input.bed} \
            -SD {input.d} \
            -O {output}
        """