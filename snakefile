'''
PREAMBLE
    1. Store all SR-ID read files in format 'reads/{SR-ID}_{end}.{suffix}'
    2. run `ls reads | cut -d '_' -f 1 | uniq >> samples.in`

'''
configfile: "test.yaml"

READS= ["SRR1801135", "SRR1801136"]

# Rule to target full run
rule all:
    input:
        "filtered_vcf/cohort.vcf"
'''
# Trim adapters off single-end reads
rule se_fastq_trim:
    input:
        'reads/{sample}.fastq.gz'
    params:
        "bbmap/resources/adapters.fa"
    output:
        'trimmed_reads/{sample}.fastq'
    shell:
        "bbduk.sh in={input} "
        "out={output} ref= {params} "
        "ktrim=r k=25"
'''

# Trim adapters off paired-end reads
rule pe_fastq_trim:
    input:
        'reads/{read}_1.fastq.gz',
        'reads/{read}_2.fastq.gz'
    params:
        "bbmap/resources/adapters.fa"
    output:
        'trimmed_reads/{read}_1.fastq',
        'trimmed_reads/{read}_2.fastq'
    shell:
        "bbduk.sh in={input} "
        "out={output} ref= {params} "
        "ktrim=r k=25"

# Merge paired ends into one read
rule align_ends:
    input:
        'trimmed_reads/{read}_1.fastq',
        'trimmed_reads/{read}_2.fastq'
        #lambda wildcards: expand("trimmed/{read}.fastq", read= config["reads"][wildcards.read])
        #lambda wildcards: config["reads"][wildcards.read]
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
        reads= expand("recal/{sample}.bam", sample=SAMPLES)
    output:
        "GVCF/{sample}.g.vcf"
    shell:
        "gatk HaplotypeCaller -R "
        "{input.ref} -T "
        "HaplotypeCaller -I {input.reads} "
        "--emitRefConfidence GVCF "
        "-O {output}"

rule create_map:
    input:
        complete= "map.done",
        dir= directory("GVCF")
    output:
        "test.sample_map"
    run:
        for sample in input.dir:
            pieces= sample.split(".")
            with open(output, "w") as out:
                out.write("{} \t {}".format(pieces[0]))

# TODO: Add step to create map file for gvcf consolidation
# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.1/org_broadinstitute_hellbender_tools_genomicsdb_GenomicsDBImport.php
# TODO: Consider adding step for automated .bed file creation
# TODO: how is this dependent on previous step?
# This step can be a pain
rule consolidate_gvcf:
    input:
        "test.sample_map"
    params:
        "Ptrichocarpa_444_v3.0.bed"
    output:
        touch("map.done"),
        directory("cohort_database")
    shell:
        "gatk GenomicsDBImport "
        "--genomicsdb-workspace-path {output} "
        "--batch-size 50 -L {params} "
        "--sample-name-map {input} "
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
