'''
PREAMBLE
    Maybe useful:
        `ls reads | cut -d '_' -f 1 | uniq >> samples.in`

'''
import itertools
import os

configfile: "test.yaml"


# Rule to target full run
rule all:
    input:
        "filtered_vcf/cohort.vcf"

# Trim adapters off single-end reads
rule fastq_trim:
    input:
        lambda wildcards: expand("reads/{end}.fastq.gz", end= config["reads"][wildcards.read]["ends"])
    params:
        ref= "bbmap/resources/adapters.fa",
        paired= lambda wildcards: config['reads'][wildcards.read]['paired']
    output:
        pe= {'trimmed_reads/{read}_1.fastq',
            'trimmed_reads/{read}_2.fastq'},
        se= 'trimmed_reads/{read}.fastq'
    run:
        if {params.paired}:
            shell("bbduk.sh in={input} "
                  "out={output.pe} ref= {params.ref} "
                  "ktrim=r k=25")
        else:
            shell("bbduk.sh in={input} "
                  "out={output.se} ref= {params.ref} "
                  "ktrim=r k=25")

# Merge paired ends into one read
rule align_ends:
    input:
        lambda wildcards: expand("trimmed_reads/{read}.fastq", read= config["reads"][wildcards.read]["ends"])
    params:
        "Ptrichocarpa_444_v3.0.fa"
    output:
        "mapped_reads/{read}.bam"
    shell:
        "bwa mem -t 12 -R {params} {input} | "
        "samtools sort -@12 -o {output} –"

# Merge reads for one sample based on config file
rule merge_sort_bam:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "merged_bam/{sample}.bam"
    shell:
        "samtools merge {output} {input} | "
        "samtools sort -@12 -o {output}"

rule read_groups:
        input:
            "merged_bam/{sample}.bam"
        output:
            "labeled_bam/{sample}.bam"
        shell:
            "picard AddOrReplaceReadGroups I={input} "
            "O={output} RGID={wildcards.sample} RGLB=NA "
            "RGPL=illumina RGPU=NA RGSM={wildcards.sample}"

# Marks alleles that may have been falsely amplified during PCR
rule mark_duplicates:
    input:
        "labeled_bam/{sample}.bam"
    output:
        "duplicate_marked/{sample}.bam"
    shell:
        "sambamba markdup {input} {output}"

## Scott hasn't seen much benefit to recalibrating and BQSR
rule recalibrate:
    input:
        ref= "Ptrichocarpa_444_v3.0.fa",
        reads= "duplicate_marked/{sample}.bam"
    output:
        "recal/{sample}_report.table"
    shell:
        "gatk BaseRecalibrator -R {input.ref} "
        "-I {input.reads} -o {output}"

# Parameterizes to detect false positives
rule BQSR:
    input:
        ref= "Ptrichocarpa_444_v3.0.fa",
        reads= "duplicate_marked/{sample}.bam",
        recal= "recal/{sample}_report.table"
    output:
        "recal/{sample}.bam"
    shell:
        "gatk ApplyBQSR -R {input.ref} "
        "-I {input.reads} --bqsr-recal-file "
        "{input.recal} -O {output}"

# samtools is faster than HaplotypeCaller
rule haplotype_caller:
    input:
        ref= "Ptrichocarpa_444_v3.0.fa",
        reads= "recal/{sample}.bam"
    output:
        "GVCF/{sample}.g.vcf"
    shell:
        "gatk HaplotypeCaller -R "
        "{input.ref} -T "
        "HaplotypeCaller -I {input.reads} "
        "--emitRefConfidence GVCF "
        "-O {output}"

# TODO will it work passing individually?
rule consolidate_gvcf:
    input:
        expand("GVCF/{sample}.g.vcf", sample= config['samples'])
    params:
        "Ptrichocarpa_444_v3.0.bed"
    output:
        directory("cohort_database")
    shell:
        "gatk GenomicsDBImport "
        "--genomicsdb-workspace-path {output} "
        "--batch-size 50 -L {params} "
        "-V {input} "
        "--reader-threads 5"

rule genotype_gvcf:
    input:
        ref= "Ptrichocarpa_444_v3.0.fa",
        db= "cohort_database"
    output:
        "unfiltered_vcf/cohort.vcf"
    shell:
        "gatk GenotypeGVCFs -R {input.ref} "
        "-V {input.db} "
        "-newQual -O {output}"

#TODO: insert quality filters flag
rule filter_vcf:
    input:
        "unfiltered_vcf/cohort.vcf"
    output:
        "filtered_vcf/cohort.vcf"
    shell:
        "picard FilterVcf -I {input} "
        "-O {output} "
        "-[insert quality filters here]"
