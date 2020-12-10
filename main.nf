/*  Minor Allele Simulation + Detection Pipeline
 *  Usage: nextflow run /path/to/main.nf
 *
 *  Author: Mohammed Khalfan < mkhalfan@nyu.edu >
 *  NYU Center for Genetics and System Biology 2020
 */

// Setting some defaults here,
// can be overridden in config or via command line
params.out = "${params.outdir}/out"
params.tmpdir = "${params.outdir}/gatk_temp"
params.snpEff_config = "static/snpEff.config"

println "ref: $params.ref"
println "outdir: $params.out"

// Stage some files we will need
ref = file(params.ref)
error_model_fq_read1 = file(params.error_model_fq_read1)
error_model_fq_read2 = file(params.error_model_fq_read2)
mut_model_vcf = file(params.mut_model_vcf)
readsim_model_bam = file(params.readsim_model_bam)

// Prepare the fastq read pairs for input.
// Use the size parameter to not auto-group, and instead
// use the mapping through getBaseName() and subtract
// two regexs to get the ID.
// This enables support for CGSB sequence data file naming format
Channel
    .fromFilePairs( params.reads, size: -1)
    { file -> file.getBaseName() - ~/${params.read_pair_regex}/ - ~/.f*q/ }
    .set { read_pairs_ch }

Channel
    .fromFilePairs( params.bams, size: 1) 
    { file -> file.getBaseName() - ~/.bam/ }
    .set { bams_in_ch }

process genMutModel{
    output:
    file('MutModel.p') into mut_model_ch

    script:
    """
    python2 /apps/neat-genreads/2.0/utilities/genMutModel.py \
        -r $ref \
        -m $mut_model_vcf \
        -o MutModel.p \
	--no-whitelist
    """
}

process seqErrorModel{
    output:
    file('SeqErrorModel.p') into seq_err_model_ch
    file('SeqErrorModel.p') into seq_err_model_ch_2

    script:
    """
    python2 /apps/neat-genreads/2.0/utilities/genSeqErrorModel.py \
        -i $error_model_fq_read1 \
        -o SeqErrorModel.p \
        -i2 $error_model_fq_read2    
    """
}

process gcModel{
    output:
    file('gc_model.p') into gc_model_ch


    script:
    """
    bedtools genomecov -d -ibam $readsim_model_bam > ${readsim_model_bam}.genomecov
    python2 /apps/neat-genreads/2.0/utilities/computeGC.py -r $ref -i ${readsim_model_bam}.genomecov -o gc_model.p
    """
}

process fraglenModel{
    output:
    file('fraglen.p') into fraglen_model_ch
    file('fraglen.p') into fraglen_model_ch_2

    script:
    """
    samtools view $readsim_model_bam | python2 /apps/neat-genreads/2.0/utilities/computeFraglen.py    
    """
}

process simulate_snvs{
    publishDir "${params.out}/readsim_1", mode:'copy'

    input:
    file(mut_model) from mut_model_ch
    file(seq_err_model) from seq_err_model_ch
    file(fraglen_model) from fraglen_model_ch

    output:
    set val(pair_id),
	file("${pair_id}_golden.vcf") into readsim_1_out_ch
    file('*') into readsim_out

    when:
    params.do_sim_reads

    script:
    pair_id = params.fcid
    """
    python2 /apps/neat-genreads/2.0/genReads.py \
	-r $ref \
	-p 100 \
	-R 151 \
	-o $pair_id \
	-e $seq_err_model \
	-m $mut_model \
	-M $params.mut_rate \
	--vcf \
	--pe-model $fraglen_model \
	--no-fastq \
	-c $params.readsim_cov
    """
}

process simulate_reads{
    publishDir "${params.out}/readsim_2", mode:'copy'

    input:
    set val(pair_id), 
    file(vcf) from readsim_1_out_ch
    file(seq_err_model) from seq_err_model_ch_2
    file(fraglen_model) from fraglen_model_ch_2
    file(gc_model) from gc_model_ch
    each af from params.readsim_allele_fracs

    output:
    set val(pair_id),
        file("${pair_id}_read1.fq"),
	file("${pair_id}_read2.fq"),
	file("${pair_id}.vcf") \
	into readsim_out_ch
    file("${pair_id}.vcf") into analyze_af_report_vcf

    script:
    pair_id = pair_id + "_AF_${af}"
    """
    # This script will edit the GT of all snps to $af
    # for consumption by genReads.py in the next step
    prepare_neat_vcf.py $vcf $ref $af > ${pair_id}.vcf
    
    # Simulate reads inserting snps from the output
    # of the above step directly into the reads
    python2 /apps/neat-genreads/2.0/genReads.py \
	-r $ref \
	-p 100 \
	-R 151 \
	-o ${pair_id} \
	-e $seq_err_model \
	-v ${pair_id}.vcf \
	--vcf \
	--pe-model $fraglen_model \
	--gc-model $gc_model \
	-c $params.readsim_cov \
	-M 0


    # Simulate reads inserting snps from the output
    # of the above step directly into the reads 
    # using ReSeq
    #reseq illuminaPE \
    #    -r $ref \
    #    -b $readsim_model_bam \
    #    -V ${pair_id}.vcf \
    #    -1 ${pair_id}_read1.fq \
    #    -2 ${pair_id}_read2.fq \
    #    -c $params.readsim_cov \
    #	--noBias
    """
}

process downsample_readsim_fq{
    publishDir "${params.out}/downsampled_fastqs", mode:'copy'

    input:
    set pair_id, 
	file(read1), 
	file(read2),
	file(vcf) \
	from readsim_out_ch
    each seed_frac_pair from params.readsim_downsample_fracs

    output:
    set val(pair_id),
        file("${pair_id}_read[12].fq") \
        into readsim_downsampled_ch
    file("${pair_id}_golden.vcf") into downsample_bzip_tabix_vcf_ch
    file("${pair_id}_golden.vcf") into golden_vcf_comp_ch
    
    script:
    seed = seed_frac_pair[0]
    frac = seed_frac_pair[1]
    pair_id = pair_id + "_frac_" + frac
    downsampled_dp = params.readsim_cov * frac


    if (frac < 1.0)
    """
    seqtk sample -s${seed} ${read1} $frac > ${pair_id}_read1.fq
    seqtk sample -s${seed} ${read2} $frac > ${pair_id}_read2.fq
    modify_neat_dp.py $vcf $downsampled_dp > ${pair_id}_golden.vcf
    """

    else
    """
    cp ${read1} ${pair_id}_read1.fq
    cp ${read2} ${pair_id}_read2.fq
    modify_neat_dp.py $vcf $downsampled_dp > ${pair_id}_golden.vcf
    """
}

process trim {
    publishDir "${params.out}/trimmed", mode:'copy'

    input:
    set pair_id,
        file(reads) from read_pairs_ch
	.mix(readsim_downsampled_ch)

    output:
    set val(pair_id),
	file("${pair_id}_trimmed_1.fq.gz"),
	file("${pair_id}_trimmed_2.fq.gz") \
	into trimmed_ch

    script:
    """
    java -jar \$TRIMMOMATIC_JAR \
	PE \
	-phred33 \
	-threads ${task.cpus} \
	${reads[0]} \
	${reads[1]} \
	${pair_id}_trimmed_1.fq.gz \
	${pair_id}.unpair_trimmed_1.fq.gz \
	${pair_id}_trimmed_2.fq.gz \
	${pair_id}.unpair_trimmed_2.fq.gz \
	#ILLUMINACLIP:${params.adapters}:2:30:10:8:true \
	LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:20
    """
}

process align {
    publishDir "${params.out}/aligned_reads", mode:'copy'
	
    input:
    set pair_id, 
	file(read_1),
	file(read_2) from trimmed_ch
     
    output:
    set val(pair_id), 
	file("${pair_id}_aligned_reads.bam") \
	into aligned_reads_ch
	
    script:
    readGroup = "@RG\\tID:${pair_id}\\tLB:${pair_id}\\tPL:${params.pl}\\tPM:${params.pm}\\tSM:${pair_id}"
    """
    bwa mem \
	-K 100000000 \
	-v 3 -t ${task.cpus} \
	-Y \
	-R \"${readGroup}\" \
	$ref \
	$read_1 \
	$read_2 \
	> ${pair_id}_aligned_reads.sam

    java -jar \$PICARD_JAR SortSam \
	I=${pair_id}_aligned_reads.sam \
	O=${pair_id}_aligned_reads.bam \
	SORT_ORDER=coordinate \
	CREATE_INDEX=true
    """
}

process markDuplicatesSpark  {
    publishDir "${params.out}/sorted", mode:'copy'

    input:
    set val(sample_id), 
	file(bam) from aligned_reads_ch
	.mix(bams_in_ch)

    output:
    set val(sample_id),
	file("${sample_id}_sorted_dedup.bam") \
	into sorted_dedup_bam_ch, sorted_dedup_ch_for_metrics, \
	downsample_bam_ch, pilon_ch, bcftools_ch, mutect2_ch, \
	tims_pipeline_ch, varscan_ch, ivar_ch, cliquesnv_ch
    set val(sample_id),
        file("${sample_id}_sorted_dedup.bam"),
        file("${sample_id}_sorted_dedup.bam.bai") \
	into bw_ch, freebayes_ch, lofreq_ch
    set val(sample_id),
	file("${sample_id}_dedup_metrics.txt") into dedup_qc_ch
    val(sample_id) into pair_id_ch

    script:
    """
    #mkdir -p $params.tmpdir/$workflow.runName/$sample_id
    #gatk --java-options "-Djava.io.tmpdir=${params.tmpdir}/${workflow.runName}/${sample_id}"
    gatk \
	MarkDuplicatesSpark \
	-I ${bam} \
	-M ${sample_id}_dedup_metrics.txt \
	-O ${sample_id}_sorted_dedup.bam
    #rm -r $params.tmpdir/$workflow.runName/$sample_id
    """ 
}

process getMetrics {
    publishDir "${params.out}/metrics", mode:'copy'

    input:
    set val(sample_id),
	file(sorted_dedup_reads) from sorted_dedup_ch_for_metrics

    output:
    set val(sample_id), 
            file("${sample_id}_alignment_metrics.txt"),
            file("${sample_id}_insert_metrics.txt"),
            file("${sample_id}_insert_size_histogram.pdf"),
            file("${sample_id}_depth_out.txt") \
            into metrics_output

    script:
    """
    java -jar \$PICARD_JAR \
        CollectAlignmentSummaryMetrics \
	R=${params.ref} \
        I=${sorted_dedup_reads} \
	O=${sample_id}_alignment_metrics.txt
    java -jar \$PICARD_JAR \
        CollectInsertSizeMetrics \
        INPUT=${sorted_dedup_reads} \
	OUTPUT=${sample_id}_insert_metrics.txt \
        HISTOGRAM_FILE=${sample_id}_insert_size_histogram.pdf 
    samtools depth -a ${sorted_dedup_reads} > ${sample_id}_depth_out.txt
    """
}

process timo{
    publishDir "${params.out}/tims-pipeline", mode:'copy'

    input:
    set val(sample_id),
        file(preprocessed_bam) from tims_pipeline_ch

    output:
    file("${sample_id}_tims.vcf") into tims_bzip_tabix_vcf_ch
    set val(sample_id),
        file("${sample_id}_tims.vcf") \
        into tims_vcf_ch
    file '*' into tims_out_ch

    script:
    """
    samtools view -bq 20 -F 1284 -o filtered.bam $preprocessed_bam
    samtools index filtered.bam
    readreport_v4_2.py \
	--infile filtered.bam \
	--ref $ref \
	-c 0.001 \
	-C 1
    ## parse_tims_output.py will look for all files created by 
    ## readreport_v4_2.py in the working directory and convert
    ## them into a single VCF file
    parse_tims_output.py $ref ${sample_id}
    mv ${sample_id}.vcf ${sample_id}_tims.vcf
    """

}

process varscan {
    publishDir "${params.out}/varscan", mode:'copy'

    input:
    set val(sample_id),
        file(preprocessed_bam) from varscan_ch

    output:
    file("${sample_id}_varscan.vcf") into varscan_bzip_tabix_vcf_ch
    set val(sample_id),
        file("${sample_id}_varscan.vcf") \
        into varscan_vcf_ch

    script:
    """
    samtools mpileup -f $ref $preprocessed_bam |\
	java -jar \$VARSCAN_JAR pileup2snp \
	--min-coverage 1 \
	--min-var-freq 0.001 > ${sample_id}_varscan.snp
    #vs_to_vcf.py ${sample_id}_varscan.snp true
    #mv ${sample_id}_varscan.vcf ${sample_id}_varscan_unfiltered.vcf
    vs_to_vcf.py ${sample_id}_varscan.snp
    """
}

process ivar{
    publishDir "${params.out}/ivar", mode:'copy'

    input:
    set val(sample_id),
        file(preprocessed_bam) from ivar_ch

    output:
    file("${sample_id}_ivar.vcf") into ivar_bzip_tabix_vcf_ch
    set val(sample_id), 
	file("${sample_id}_ivar.vcf") \
	into ivar_vcf_ch

    script:
    """
    samtools mpileup -aa -A -d 0 -B -Q 0 ${preprocessed_bam} | ivar variants -p ${sample_id}_ivar -t 0.001 -r $ref
    ivar_to_vcf.py ${sample_id}_ivar.tsv
    """
}

process lofreq{
    publishDir "${params.out}/lofreq", mode:'copy'

    input:
    set val(sample_id),
        file(preprocessed_bam),
	file(preprocessed_bam_index) from lofreq_ch

    output:
    file("${sample_id}_lofreq.vcf") into lofreq_bzip_tabix_vcf_ch
    set val(sample_id),
        file("${sample_id}_lofreq.vcf") \
        into lofreq_vcf_ch
    file '*' into lofreq_out_ch

    script:
    """
    lofreq call-parallel \
	--pp-threads ${task.cpus} \
	-f $ref \
	-o ${sample_id}_lofreq.vcf \
	$preprocessed_bam
    """
}

process cliquesnv{
    publishDir "${params.out}/cliquesnv", mode:'copy'

    input:
    set val(sample_id),
        file(preprocessed_bam) from cliquesnv_ch

    output:
    file("${sample_id}_cliquesnv.vcf") into cliquesnv_bzip_tabix_vcf_ch
    set val(sample_id),
        file("${sample_id}_cliquesnv.vcf") \
        into cliquesnv_vcf_ch
    file '*' into cliquesnv_out_ch

    when:
    false

    script:
    """
    java  -Xmx58G -jar \$CLIQUESNV_JAR \
	-m snv-illumina-vc \
	-in $preprocessed_bam \
	-outDir vcf_out/ \
	-t 1 \
	-tf 0.005
    mv vcf_out/${preprocessed_bam.baseName}.vcf ${sample_id}_cliquesnv.vcf
    
    # CliqueSNV names contig as 'ref' for some reason, rename it to 'SARS-CoV2'
    # here for downstream analysis (ex: bcftools isec). Todo: Make this param
    sed -i 's/ref/SARS-CoV2/g' ${sample_id}_cliquesnv.vcf
    """
}

process freebayes{
    publishDir "${params.out}/freebayes", mode:'copy'

    input:
    set val(sample_id),
        file(preprocessed_bam),
        file(preprocessed_bam_index) from freebayes_ch
    each freebayes_config from params.freebayes_configs

    output:
    file("${sample_id}_freebayes_${name}.vcf") into freebayes_bzip_tabix_vcf_ch
    set val(sample_id),
        file("${sample_id}_freebayes_${name}.vcf") \
        into freebayes_vcf_ch
    file '*' into freebayes_out_ch

    script:
    name = freebayes_config[0]
    fb_params = freebayes_config[1]
    """
    freebayes $fb_params -f $ref $preprocessed_bam > ${sample_id}_freebayes_${name}.vcf
    """
}


process mutect2{
    publishDir "${params.out}/mutect2", mode:'copy'

    input:
    set val(sample_id),
        file(preprocessed_bam) from mutect2_ch

    output:
    set val(sample_id), 
	file("${sample_id}_mutect2_filtered_pf.vcf") \
	into mutect2_vcf_ch
    file("${sample_id}_mutect2_filtered_pf.vcf") into mutect2_bzip_tabix_vcf_ch
    file '*' into mutect2_out_ch

    script:
    """
    gatk Mutect2 -R $ref -I $preprocessed_bam -O ${sample_id}_mutect2.vcf
    gatk FilterMutectCalls -R $ref -V ${sample_id}_mutect2.vcf -O ${sample_id}_mutect2_filtered.vcf
    gatk SelectVariants -R $ref -V ${sample_id}_mutect2_filtered.vcf --exclude-filtered -O ${sample_id}_mutect2_filtered_pf.vcf
    """
}

process pilon{
    publishDir "${params.out}/pilon", mode:'copy'

    input:
    set val(sample_id),
	file(preprocessed_bam) from pilon_ch

    output:
    file("${sample_id}_pilon.vcf") into pilon_bzip_tabix_vcf_ch
    file '*' into pilon_out_ch

    script:
    """
    java -Xmx16G -jar \$PILON_JAR \
	--genome $ref \
	--bam $preprocessed_bam \
	--fix bases \
	--changes \
	--vcf \
	--threads ${task.cpus} \
	--mindepth 10 \
	--output ${sample_id}_pilon_g
    
    gatk SelectVariants \
	-V ${sample_id}_pilon_g.vcf \
	-O ${sample_id}_pilon.vcf \
	--exclude-non-variants \
	--exclude-filtered
    """
}

process bcftools{
    publishDir "${params.out}/bcftools", mode:'copy'

    input:
    set val(sample_id),
        file(preprocessed_bam) from bcftools_ch

    output:
    file("${sample_id}_bcftools.vcf") into bcftools_bzip_tabix_vcf_ch

    script:
    """
    bcftools mpileup \
	--redo-BAQ \
	--adjust-MQ 50 \
	--gap-frac 0.05 \
	--max-depth 10000 \
	--max-idepth 200000 \
	--fasta-ref $ref \
	$preprocessed_bam | bcftools call \
	--ploidy 1 \
	--keep-alts \
	--multiallelic-caller \
	--variants-only \
	--output ${sample_id}_bcftools.vcf
    """
}

process haplotypeCaller {
    input:
    set val(sample_id), 
	file(preprocessed_bam) from sorted_dedup_bam_ch

    output:
    set val(sample_id),
        file("${sample_id}_raw_snps.vcf") \
        into raw_snps_ch, raw_snps_qc_ch

    script:
    """
    gatk HaplotypeCaller \
	-R $ref \
	-I $preprocessed_bam \
	-O ${sample_id}_raw_variants.vcf \
	-ploidy 100
    gatk SelectVariants \
        -R $ref \
        -V ${sample_id}_raw_variants.vcf \
        -select-type SNP \
        -O ${sample_id}_raw_snps.vcf
    """
}

process filterSnps {
    publishDir "${params.out}/filtered_snps", mode:'copy'
    
    input:
    set val(sample_id), 
	file(raw_snps) from raw_snps_ch

    output:
    set val(sample_id),
        file("${sample_id}_filtered_snps.vcf") \
        into filtered_snps_qc_ch, hc_vcf_ch

    script:
    """
    gatk VariantFiltration \
	-R $ref \
	-V $raw_snps \
	-O ${sample_id}_filtered_snps.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

    # This script generates the _consensus_snps.vcf
    # and _eaf.vcf (empirical AF) files
    filter_variants.py ${sample_id}
    """
}

process make_bw{
    publishDir "${params.out}/bigwig", mode:'copy'

    input:
    set val(id), 
	file(bam),
	file(bam_index) \
	from bw_ch

    output:
    file("${id}_coverage.bam.bw") into jbrowse_bw_ch 

    when:
    false

    script:
    """
    bamCoverage \
        -p max  \
        --bam $bam \
	--binSize 1 \
	--ignoreDuplicates \
	--minMappingQuality 20 \
        -o ${id}_coverage.bam.bw
    """
}

process downsample_bam{
    input:
    set val(sample_id), file(bam) from downsample_bam_ch

    output:
    set file("${sample_id}_downsampled.bam"),
        file("${sample_id}_downsampled.bam.bai") into jbrowse_bam_ch

    when:
    false

    script:
    """
    java -jar \$SORTSAMREFNAME_JAR \
        --samoutputformat BAM \
        $bam |\
        java -jar \$BIOSTAR_JAR \
        -n 75 \
        --samoutputformat BAM |\
        samtools sort -o ${sample_id}_downsampled.bam
    samtools index ${sample_id}_downsampled.bam
    """
}

process bzip_tabix_vcf{
    input:
    file(vcf) from pilon_bzip_tabix_vcf_ch
	.mix(mutect2_bzip_tabix_vcf_ch)
	.mix(freebayes_bzip_tabix_vcf_ch)
	.mix(downsample_bzip_tabix_vcf_ch)
	.mix(tims_bzip_tabix_vcf_ch)
	.mix(varscan_bzip_tabix_vcf_ch)
	.mix(ivar_bzip_tabix_vcf_ch)
        .mix(lofreq_bzip_tabix_vcf_ch)
	.mix(cliquesnv_bzip_tabix_vcf_ch)

    output:
    file("*.vcf.gz*") into jbrowse_vcf_ch

    when:
    false

    script:
    """
    bgzip -c ${vcf} > ${vcf}.gz
    tabix -p vcf ${vcf}.gz
    """
}

/* Collect all the pair_ids and send them to pair_id_list_ch for use in jbrowse.
 */
pair_id_ch.collectFile(storeDir: "${params.out}/reports") { item ->
       [ "sample_ids.txt", item + '\n' ]
}.tap{pair_id_list_ch}

process jbrowse{
    cache false

    publishDir "${params.out}/trackList", mode:'copy'

    input:
    file '*' from jbrowse_bw_ch.collect()
    file '*' from jbrowse_bam_ch.collect()
    file '*' from jbrowse_vcf_ch.collect()
    file(pair_id_list) from pair_id_list_ch

    when:
    false

    script:
    """
    cgsb_upload2jbrowse -p MAD -d $params.fcid -s $pair_id_list -f . $ref
    """
}

process compare_afs{
    cache false

    publishDir "${params.out}/compare_afs", mode:'copy'

    input:
    file '*' from golden_vcf_comp_ch.collect()
    set val(pair_id), 
	file(workflow_vcf) \
	from ivar_vcf_ch
	.mix(freebayes_vcf_ch)
	.mix(mutect2_vcf_ch)
	.mix(tims_vcf_ch)
	//.mix(tims_nb_vcf_ch)
	.mix(hc_vcf_ch)
	.mix(varscan_vcf_ch)
        .mix(lofreq_vcf_ch)
        .mix(cliquesnv_vcf_ch)

    output:
    file("${pair_id}_af_report.csv") into make_af_csv_output

    script:
    // This is important. The pair id is used to set
    // the golden vcf for this comparison. 
    // the pair id contains all the variables
    // e.g. the AF and DP info
    golden_vcf=pair_id + "_golden.vcf"
    """
    # need to bzip and index the vcf for bcftools
    bgzip -c ${golden_vcf} > ${golden_vcf}.gz
    tabix -p vcf ${golden_vcf}.gz
    bgzip -c ${workflow_vcf} > ${workflow_vcf}.gz
    tabix -p vcf ${workflow_vcf}.gz

    bcftools isec -c none ${golden_vcf}.gz ${workflow_vcf}.gz -p isec

    compare_afs.py --golden_vcf ${golden_vcf}.gz --workflow_vcf ${workflow_vcf}.gz --out ${pair_id}_af_report
    """
}

process qc {
    input:
    set val(sample_id),
	file("${sample_id}_alignment_metrics.txt"),
	file("${sample_id}_insert_metrics.txt"),
	file("${sample_id}_insert_size_histogram.pdf"),
	file("${sample_id}_depth_out.txt"),
	file("${sample_id}_dedup_metrics.txt"),
	file("${sample_id}_raw_snps.vcf"),
        file("${sample_id}_filtered_snps.vcf") \
	from metrics_output
	.join(dedup_qc_ch)
	.join(raw_snps_qc_ch)
	.join(filtered_snps_qc_ch)

    output:
    file("${sample_id}_report.csv") into parse_metrics_output

    when:
    false

    script:
    """
    parse_metrics.sh ${sample_id} > ${sample_id}_report.csv 
    """
}

/* Process qc above creates a report for each sample.
 * Below we compile these into a single report.
 * Same for af_csv files
 */
//parse_metrics_output.collectFile(name: "${workflow.runName}_report.csv", keepHeader: true, storeDir: "${params.out}/reports")
make_af_csv_output.collectFile(name: "${workflow.runName}_af_data.csv", keepHeader: true, storeDir: "${params.out}/reports").tap{af_report_in}

process analyze_af_report {
    publishDir "${params.out}/reports", mode:'copy'	

    input:
    file(af_report) from af_report_in
    file(golden_vcf) from analyze_af_report_vcf.take(1)

    output:
    file("*") into analyze_af_repot_out

    when:false

    script:
    """
    R -e "af_report = '${af_report}'; golden_vcf = '${golden_vcf}'; rmarkdown::render('analyze_af_report.Rmd')"
    """
}
