'''
PREAMBLE
    Maybe useful:
        `ls reads | cut -d '_' -f 1 | uniq >> samples.in`

'''
import itertools
import os

xx = sorted(glob_wildcards('test_reads_100k/{sample}.fastq.gz').sample)

# print(xx)

fq1 = [x for x in xx if '_1' in str(x)]
fq2 = [x for x in xx if '_2' in str(x)]

paired = [x.replace('_1', '') for x in fq1]

# print(paired)

fq3 = fq1 + fq2

single = [x for x in xx if x not in fq3]
# print(single)

configfile: "test.yaml"


# Rule to target full run
rule all:
    input:
        "filtered_vcf/cohort.vcf"

# Trim adapters off single-end reads
rule fastq_trim_pe:
    input:
        fq1 = "test_reads_100k/{sample}_1.fastq.gz",
        fq2 = "test_reads_100k/{sample}_2.fastq.gz"
    params:
        ref= "adapters"

    output:
        fq1o = 'trimmed_reads/{sample}_1.fastq',
        fq2o = 'trimmed_reads/{sample}_2.fastq'

    shell:
        "bbduk.sh in1={input.fq1} in2={input.fq2} "
        "out1={output.fq1o} out2={output.fq2o} ref={params.ref} "
        "ktrim=r k=25"

rule fastq_trim_se:
    input:
        fq = "test_reads_100k_se/{sample}.fastq.gz"
    params:
        ref = "adapters"
    output:
        fqo = "trimmed_reads/{sample}.fastq"
    shell:
        "bbduk.sh in={input.fq} "
        "out={output.fqo} ref={params.ref} "
        "ktrim=r k=25"

    # run:
    #     # out1_index = [x for x in output.pe if '_1.fastq' in str(x)]
    #     # out2_index = [x for x in output.pe if '_2.fastq' in str(x)]
    #     input1=str(input.fq1[0])
    #     input2=str(input.fq2[0])
    #     out1=str(output.fq1o[0])
    #     out2=str(output.fq2o[0])
    #     # out1=str(out1_index[0])
    #     # out2=str(out2_index[0])
    #     # print(input1)
    #     # print(input2)
    #     print(out1)
    #     print(out2)
    #     shell("bbduk.sh in1={input.fq1) in2={input2} "
    #               "out1={out1} out2={out2} ref={params.ref} "
    #               "ktrim=r k=25")


# rule fastq_trim_se:
#     input:
#         fq1 = expand(trimmed_reads/)\
# else:
#     shell("bbduk.sh in={input} "
#           "out={output.se} ref= {params.ref} "
#           "ktrim=r k=25")

# Merge paired ends into one read
rule align_ends:
    input:
        lambda wildcards: expand("trimmed_reads/{read}.fastq", read= config["reads"][wildcards.read]["ends"])
    params:
        "Ptrichocarpa_444_v3.0.fa"
    output:
        "mapped_reads/{read}.bam"
    shell:
        "bwa mem -t 1 {params} {input} | "
        "samtools sort -@1 -o {output}"

# Merge reads for one sample based on config file
rule merge_sort_bam:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "merged_bam/{sample}.bam"
    shell:
        "samtools merge {output} {input} "

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
        "-I {input.reads} -O {output} --known-sites Chr13_subset_882.vcf"

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
        "{input.ref} "
        "-I {input.reads} "
        "--emit-ref-confidence GVCF "
        "-O {output}"

#https://bitbucket.org/snakemake/snakemake/issues/895/combine-multiple-files-for-input-but
rule consolidate_gvcf:
    input:
        expand("GVCF/{sample}.g.vcf", sample= config['samples'])
    params:
        # ref = "Ptrichocarpa_444_v3.0.bed"
        lambda wildcards, input: " -V ".join(input)
    output:
        directory("cohort_database")
    shell:
        "gatk GenomicsDBImport "
        "--genomicsdb-workspace-path {output} "
        "--batch-size 50 -L Ptrichocarpa_444_v3.0.bed "
        "-V {params} "
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
