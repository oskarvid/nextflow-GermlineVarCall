#! /usr/bin/env nextflow

params.gatk4 = "/home/oskar/wdl_pipeline_bak/tools/GATK-4.jar"
params.gatk3 = "/home/oskar/wdl_pipeline_bak/tools/GenomeAnalysisTK.jar"
params.read1 = "/home/oskar/01-workspace/00-temp/nextflow-testing/sample.markdup.sortsam.bwa.bam"
params.fasta_ref = "/home/oskar/01-workspace/04-pipelines/GATK-Ghislain/ref_filer/human_g1k_v37_decoy.fasta"
params.dbsnp_vcf = "/home/oskar/01-workspace/04-pipelines/GATK-Ghislain/ref_filer/dbsnp.vcf"
params.dbsnp_vcf_index = "/home/oskar/01-workspace/04-pipelines/GATK-Ghislain/ref_filer/dbsnp.vcf.idx"
params.mills_vcf = "/home/oskar/01-workspace/04-pipelines/GATK-Ghislain/ref_filer/Mills_and_1000G_gold_standard.indels.b37.vcf"
params.mills_vcf_index = "/home/oskar/01-workspace/04-pipelines/GATK-Ghislain/ref_filer/Mills_and_1000G_gold_standard.indels.b37.vcf.idx"
params.hapmap_vcf = "/home/oskar/01-workspace/04-pipelines/GATK-Ghislain/ref_filer/hapmap.vcf"
params.hapmap_vcf_index = "/home/oskar/01-workspace/04-pipelines/GATK-Ghislain/ref_filer/hapmap.vcf.idx"
params.v1000g_vcf = "/home/oskar/01-workspace/04-pipelines/GATK-Ghislain/ref_filer/1000g.vcf"
params.v1000g_vcf_index = "/home/oskar/01-workspace/04-pipelines/GATK-Ghislain/ref_filer/1000g.vcf.idx"
params.omni_vcf = "/home/oskar/01-workspace/04-pipelines/GATK-Ghislain/ref_filer/omni.vcf"
params.omni_vcf_index = "/home/oskar/01-workspace/04-pipelines/GATK-Ghislain/ref_filer/omni.vcf.idx"
params.contigs = "/home/oskar/wdl_pipeline_bak/intervals/groups.list"

//contigs = file(params.contigs)
gatk4 = file(params.gatk4)
gatk3 = file(params.gatk3)
read1 = file(params.read1)
fasta_ref = file(params.fasta_ref)
fasta_ref_fai = file( params.fasta_ref+'.fai' )
fasta_ref_sa = file( params.fasta_ref+'.sa' )
fasta_ref_bwt = file( params.fasta_ref+'.bwt' )
fasta_ref_ann = file( params.fasta_ref+'.ann' )
fasta_ref_amb = file( params.fasta_ref+'.amb' )
fasta_ref_pac = file( params.fasta_ref+'.pac' )
fasta_ref_dict = file( params.fasta_ref.replace(".fasta",".dict") )
dbsnp_vcf = file(params.dbsnp_vcf)
v1000g_vcf = file(params.v1000g_vcf)
mills_vcf = file(params.mills_vcf)
dbsnp_vcf_index = file(params.dbsnp_vcf_index)
v1000g_vcf_index = file(params.v1000g_vcf_index)
mills_vcf_index = file(params.mills_vcf_index)
omni_vcf = file(params.omni_vcf)
omni_vcf_index = file(params.omni_vcf_index)
hapmap_vcf = file(params.hapmap_vcf)
hapmap_vcf_index = file(params.hapmap_vcf_index)

Channel
    .fromPath( params.read1 )
    .ifEmpty { error "Cannot find any reads matching: ${params.read1}" }
	.set { read_id }

Channel.fromPath(params.contigs)
        .splitText()
        .map { file(it) }
		.set { intervals }

process BaseRecalibrator {
	echo true

	input:
	file gatk4
	file fasta_ref
	file fasta_ref_fai
	file fasta_ref_dict
	file dbsnp_vcf
	file v1000g_vcf
	file mills_vcf
	file dbsnp_vcf_index
	file v1000g_vcf_index
	file mills_vcf_index
	file intervals
	file read from read1

	output:
	file("recalibration_report.grp") into BaseRecalibrator_output

	script:
	def intervals = intervals.collect{ "-L $it" }.join(" ")
	"java -Dsnappy.disable=true -XX:ParallelGCThreads=4 \
	  -Dsamjdk.use_async_io=false -Xmx13G \
	  -jar $gatk4 \
	  BaseRecalibrator \
	  --reference $fasta_ref \
	  --input $read \
	  -O recalibration_report.grp \
	  --knownSites $dbsnp_vcf \
	  --knownSites $v1000g_vcf \
	  --knownSites $mills_vcf \
	  $intervals"
}