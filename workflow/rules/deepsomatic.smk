#@TODO:
# This rule should be adapted to work in either tumor-only or tumor-normal mode
rule deepsomatic:
    input:
        unpack(get_deepsomatic_input),
        ref=config['resources']['reference_fasta'],
        regions=regions_bed,
    output:
        vcf=os.path.join(
            config['output_folder'],
            "vcfs",
            "{patient}.deepsomatic.vcf.gz"
            ),
    params:
        tumor="{patient}.tumor",
        logging_dir=os.path.join(
            config['log_folder'],
            "deepsomatic",
            "{patient}_log"
        ),
        intermediate_results_dir=os.path.join(
            config['output_folder'],
            "deepsomatic_intermediate",
            "{patient}"
        ),
        threads=config['params']['deepsomatic']['threads'],
        extra=config['params']['mutect2']['args'],
        model_type=config['params']['deepsomatic']['model_type'],
    container: config['containers']['deepsomatic']
    log:
        os.path.join(
            config['log_folder'],
            "deepsomatic",
            "{patient}.overall.log"
        )
    shell:
        """
        run_deepsomatic \
            --model_type={params.model_type} \
            --ref={input.ref} \
            --reads_tumor={input.cram} \
            --output_vcf={output.vcf} \
            --sample_name_tumor={params.tumor} \
            --num_shards={params.threads} \
            --logging_dir={params.logging_dir} \
            --intermediate_results_dir={params.intermediate_results_dir} \
            --use_default_pon_filtering=true \
            --regions={input.regions} 
        """
