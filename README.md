# Variant Calling Pipeline

## Running the Pipeline
1. The config file (currently using `test.yaml`) must be updated to run with the samples you want until we create a script to automate this. The format should be clear from the current file but ask if it is not.

2. To run the entire pipeline simply execute `snakemake` and it will process for any samples specified in the config file. If you would like to only partly execute the pipeline, specify the output of the last rule you would like executed, i.e. `snakemake duplicate_marked/{sample}.bam` with {sample} replaced with the desired samples to process.

3. Update README as we add features
