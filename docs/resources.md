This pipeline uses a reference genome and some additional required references that must be downloaded are accessible for reading and writing by the user who's running the pipeline. The reference genome in use is the GRCh38. The suggested version to use for variant calling is the version [without alternative contigs](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) and, more in detail, the GIAB v3 version. This could be downloaded using:

```bash
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz && gzip -d GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
```

The PATH to the obtained genome must be placed in `config/config_main.yaml` under the `resources/reference_fasta` field. 

### Known alterations

Multiple files reporting known source are required for alignment preprocessing and variant calling. These files could be huge in size, so be sure to store them in a folder with enough space and with the right reading permission. Use the following table to fill the `config/config_main.yaml` file with the right paths to the files:

| File | config key | VCF | Index| 
| --- | --- | --- | --- |
| Known SNPs | dbsnps | ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz | ftp://ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi |
| Known indels 1000G | known_indels | https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz | https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi |
| Known snps 1000G | 1000g | https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz | https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi |
| gnomAD v3 | gnomad | https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz | https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi |
| small exac | contamination | https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz | https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi |
| Panel of Normals | PoN | https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz | https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi |

### Cache for vep

The pipeline uses VEP in version 105 for the annotation of detected variants. It's better to have a local copy of the cache files required for VEP. The cache files can be downloaded using the following command:

```bash
wget https://ftp.ensembl.org/pub/release-105/variation/indexed_vep_cache/homo_sapiens_vep_105_GRCh38.tar.gz
tar -xvzf homo_sapiens_vep_105_GRCh38.tar.gz
```

Write then the folder containing the cache under `resources/vep_cache` in `config/config_main.yaml`. 

