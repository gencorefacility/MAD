# Minor Allele Detection Pipeline

Identifying minor-allele single nucleotide variants (SNVs) in SARS-CoV2 is essential for understanding the evolution of the virus, and can have implications in public health, vaccine development, and treatment. One important question in the emergence of new mutations is whether we can effectively identify new variants that arise in an individual before they become fixed in the population. To that end, we developed the MAD (Minor Allele Detection) Pipeline using Nextflow to evaluate the performance of various different variant callers (iVar, VarScan, HaplotypeCaller, Mutect2, lofreq, freebayes, and our in-house pipeline timo) in identifying minority variants and their allele frequencies.

The MAD pipeline is a fully containerized, out of the box solution for evaluating variant callers for any organism, requiring the user to provide a set of SNPs to build a mutation model (in VCF format), an aligned bam file (for modeling), paired-end reads (in fastq format), and a reference genome. Given this organism and sequencer specific input, MAD will build the appropriate models, simulate SNVs and fastq reads, downsample the reads, perform alignment, deduplication, variant calling, and produce plots summarizing the results in an HTML report. 

## Requirements:
Nextflow v20.07.1

## Docker: 
All required software is packeged in a docker container available at: https://hub.docker.com/r/gencorefacility/mad

## Config:
Download and edit the following parameters in the config file (nextflow.config)

`params.ref`: path to reference genome fasta. BWA index, fasta index, and picard reference dictionary must exist in the same dir.

`params.outdir`: output path

`params.mut_model_vcf`: path to VCF file to build mutation model

`params.error_model_fq_read1`: Read 1 of paired end fastq reads for error model

`params.error_model_fq_read2`: Read 2 of paired end fastq reads for error model

`params.readsim_model_bam`: Bam file for modeling

`params.mut_rate`: mutation rate (default: 0.45%)

`params.readsim_cov`: Simulation coverage (default: 10000)

`params.readsim_downsample_fracs`: Simulation downsampling fractions [random seed, fraction]

`params.readsim_allele_fracs`: Simulation allele frequencies


`process`: configuration for scheduler. Replace 'slurm' with execution scheduler. Resources required will vary greatly depending on input, however reasonable defaults have been provided by default. 

## Running the Pipeline:
Call the pipeline directly from github using nextflow, providing the path to the config, and specifying the -with-docker or -with-singularity parameter along with the docker repo.

`nextflow run gencorefacility/mad -c <path_to_config> -with-docker gencorefacility/mad:1`

## Output:
An HTML report summarizing the results will be stored in `params.outdir/out/reports`


