#TODO: add the deepsomatic rule 
# rule deepsomatic:
#     input:
#         unpack(get_deepsomatic_input)
#     output:
#         vcf=
#     params:
#     container: config['containers']['deepsomatic']
#     resources:
#     log:
#     shell:

rule mutect2:
    input:
        unpack(get_mutect2_input),
        ref=config['resources']['reference_fasta'],
        gnomad=config['resources']['gnomad'],
        interval=config['output_folder'] +
        "interval-files/{interval}-scattered.interval_list"
    output:
        vcf=config['output_folder'] +
        temp("vcfs/{patient}.{interval}.unfiltered.vcf"),
        stats=config['output_folder'] +
        temp("vcfs/{patient}.{interval}.unfiltered.vcf.stats"),
        f1r2tar=config['output_folder'] +
        temp("vcfs/{patient}.{interval}.f1r2.tar.gz")
    params:
        tumor="{patient}.tumor",
        normalname= ' ' if tumor_only else '-normal ' + "{patient}.normal",
        normal_input=lambda wildcards, input: ' ' if tumor_only else "-I " + input.normal,
        pon="--panel-of-normals " + config['resources']['PoN'],
        extra=config['params']['mutect2']['args'],
    container: config['containers']['ctgflow_core']
    resources:
    log:
    shell:
        """
        gatk Mutect2 -R {input.ref} -I {input.tumor} \
            -tumor {params.tumor} -O {output.vcf} \
            --germline-resource {input.gnomad} \
            --f1r2-tar-gz {output.f1r2tar} \
            -L {input.interval} \
            {params.normal_input} {params.normalname} \
            {params.pon} {params.extra}
        """

rule orientation_bias:
    input:
        get_orientationbias_input
    output:
        temp(
            config['output_folder']
            + "vcfs/{patient}.read_orientation_model.tar.gz")
    params:
        i=lambda wildcards, input: ['-I ' + d for d in input]
    container: config['containers']['ctgflow_core']
    resources:
    log:
    shell:
        """
        gatk LearnReadOrientationModel {params.i} -O {output}
        """

rule pileup_summaries:
    input:
        bam=config['output_folder']
        + "bams/{patient}.{sample_type}.cram",
        ref=config['resources']['reference_fasta'],
        germ_res=config['resources']['contamination']
    output:
        config['output_folder'] +
        "qc/{patient}_{sample_type}_pileupsummaries.table"
    container: config['containers']['ctgflow_core']
    resources:
    log:
    shell:
        """
        gatk GetPileupSummaries -I {input.bam} \
        -V {input.germ_res} \
        -L {input.germ_res} -O {output} \
        -R {input.ref}
        """

rule calculate_contamination:
    input:
        unpack(get_contamination_input)
    output:
        config['output_folder'] +
        "qc/{patient}_contamination.table"
    params:
        matched=lambda wildcards, input:'' if tumor_only else '-matched ' + input.normal
    container: config['containers']['ctgflow_core']
    resources:
    log:
    shell:
        """
        gatk CalculateContamination -I {input.tumor}  \
            {params.matched} \
            -O {output}
        """

rule merge_vcfs:
    input:
        get_mergevcfs_input
    output:
        vcf=config['output_folder']
        + "vcfs/{patient}.unfiltered.vcf",
    params:
        i=lambda wildcards, input: ['-I ' + vcf for vcf in input]
    container: config['containers']['ctgflow_core']
    resources:
    log:
    shell:
        """
        gatk MergeVcfs {params.i} -O {output.vcf}
        """

rule merge_stats:
    input:
        get_mergestats_input
    output:
        stats=config['output_folder']
        + "vcfs/{patient}.unfiltered.vcf.stats",
    params:
        i=lambda wildcards, input: ['-stats ' + s for s in input]
    container: config['containers']['ctgflow_core']
    resources:
    log:
    shell:
        """
        gatk MergeMutectStats {params.i} -O {output.stats} 
        """

rule filter_calls:
    input:
        vcf=config['output_folder']
        + "vcfs/{patient}.unfiltered.vcf",
        ref=ref_fasta,
        contamination=config['output_folder']
        + "qc/{patient}_contamination.table",
        stats=config['output_folder']
        + "vcfs/{patient}.unfiltered.vcf.stats",
        f1r2model=config['output_folder']
        + "vcfs/{patient}.read_orientation_model.tar.gz"
    output:
        vcf=temp(
            config['output_folder']
            + "vcfs/filtered/{patient}.filtered.vcf"),
        idx=temp(
            config['output_folder']
            + "vcfs/filtered/{patient}.filtered.vcf.idx"),
        inter_stats=config['output_folder']
        + "vcfs/filtered/{patient}.filtered.vcf.filteringStats.tsv",
    params:
        extra=config['params']['mutect2']['filtering']
    container: config['containers']['ctgflow_core']
    resources:
    log:
    shell:
        """
        gatk FilterMutectCalls -V {input.vcf} -R {input.ref} \
            --contamination-table {input.contamination} \
            --stats {input.stats} \
            -ob-priors {input.f1r2model} \
            -O {output.vcf} \
            {params.extra}
        """

rule compress_calls:
    input:
        vcf=config['output_folder']
        + "vcfs/filtered/{patient}.filtered.vcf",
    output:
        vcf=config['output_folder']
        + "vcfs/filtered/{patient}.withFilters.vcf.gz"
    params:
    container:
        config["containers"]["ctgflow_core"],
    resources:
    log:
    shell:
        """
        bgzip -c {input.vcf} > {output.vcf} 
        tabix -p vcf {output.vcf}
        """


rule select_calls:
    input:
        vcf=config['output_folder']
        + "vcfs/filtered/{patient}.withFilters.vcf.gz",
    output:
        vcf=temp(
        config['output_folder']
        + "vcfs/filtered/{patient}.filtered.vcf.gz"),
        idx=temp(
        config['output_folder']
        + "vcfs/filtered/{patient}.filtered.vcf.gz.tbi"),
    params:
    container:
        config['containers']['ctgflow_core']
    resources:
        # use default resources here
    log:
    shell:
        """
        bcftools view -f .,PASS {input.vcf} \
        -Oz -o {output.vcf} ;
        tabix -p vcf {output.vcf}
        """

rule vep:
    input:
        vcf=config['output_folder']
        + "vcfs/filtered/{patient}.filtered.vcf.gz",
        idx=config['output_folder']
        + "vcfs/filtered/{patient}.filtered.vcf.gz.tbi",
        fasta=config['resources']['reference_fasta'],
        vep_cache=config['resources']['vep_cache'],
    output:
        vcf=temp(
            config['output_folder']
            + "vcfs/{patient}.vep.vcf"),
    params:
        extra=config['params']['vep']['extra'],
        assembly=config['params']['vep']['assembly'],
    container:
        config['containers']['ctgflow_core']
    resources:
    log:
    shell:
        """
        vep --input_file {input.vcf} --output_file {output.vcf} \
            --fasta {input.fasta} --cache --dir {input.vep_cache} \
            --assembly {params.assembly} {params.extra}
        """

rule compress_vcf:
    input:
        vcf=config['output_folder']
        + "vcfs/{patient}.vep.vcf",
    output:
        vcf=config['output_folder']
        + "vcfs/{patient}.vep.vcf.gz",
        idx=config['output_folder']
        + "vcfs/{patient}.vep.vcf.gz.tbi",
    params:
    container:
        config['containers']['ctgflow_core']
    resources:
    log:
    shell:
        """
        bgzip -c {input.vcf} > {output.vcf} ;
        tabix -p vcf {output.vcf}
        """
