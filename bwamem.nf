#! /usr/bin/env nextflow

//params.fastqs = "/home/oskar/01-workspace/01-data/fastq/kek/NA12878*R*.fastq"
params.gatk4 = "/home/oskar/wdl_pipeline_bak/tools/GATK-4.jar"
params.gatk3 = "/home/oskar/wdl_pipeline_bak/tools/GenomeAnalysisTK.jar"
params.read1 = "/home/oskar/wdl_pipeline_bak/data/test_R1.fastq"
params.read2 = "/home/oskar/wdl_pipeline_bak/data/test_R2.fastq"
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
params.input_tsv = "/home/oskar/01-workspace/00-temp/wdl_pipeline/intervals/template_sample_manifest_na12878.tsv"
params.contigs = "/home/oskar/01-workspace/00-temp/wdl_pipeline/intervals/groups.list"

gatk4 = file(params.gatk4)
gatk3 = file(params.gatk3)
read1 = file(params.read1)
read2 = file(params.read2)
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
  .fromFilePairs( "/home/oskar/01-workspace/01-data/fastq/kek/test-L*_R{1,2}.fastq", flat: true)
  .into { reads; reads2 }

Channel
  .fromPath('/home/oskar/01-workspace/00-temp/nextflow-testing/template_sample_manifest_na12878.tsv')
  .splitCsv(sep:'\t') 
  .map { cols -> tuple(file(cols[0]), file(cols[1]), file(cols[2]), file(cols[3]), file(cols[4]), file(cols[5]), file(cols[6])) }
  .into { tsv_ch1; tsv_ch2 }

process BwaMem {
	tag { sample_name }
	maxForks = 1
    
    input:
	file reads from tsv_ch1
    file fasta_ref
    file fasta_ref_fai
    file fasta_ref_sa
    file fasta_ref_bwt
    file fasta_ref_ann
    file fasta_ref_amb
    file fasta_ref_pac

    output:
    file("bwamem.sam") into BwaMem_output

    script:
    """
    bwa mem -t 2 $fasta_ref \
      -R '@RG\\tID:${reads[1]}\\tSM:${reads[0]}\\tLB:${reads[5]}\\tPL:${reads[6]}\\tPU:NotDefined' \
      -M ${reads[3]} ${reads[4]} > bwamem.sam
    """
}

process FastqToSam {

	input:
	file reads from tsv_ch2
	file gatk4

	output:
	file("FastqToSam.bam") into FastqToSam_output

	script:
	"""
	java -Xmx16G -Dsnappy.disable=true -XX:ParallelGCThreads=4 -Djava.io.tmpdir=`pwd`/tmp -jar \
      $gatk4 \
      FastqToSam \
      --FASTQ ${reads[3]} \
      --FASTQ2 ${reads[4]} \
      -O FastqToSam.bam \
      --SAMPLE_NAME ${reads[0]} \
      --READ_GROUP_NAME ${reads[1]} \
      --LIBRARY_NAME ${reads[5]} \
      --PLATFORM ${reads[6]} \
      --SORT_ORDER coordinate \
      2> info
	"""
}

sam_and_bam_ch = FastqToSam_output.phase(BwaMem_output).map { left, right -> tuple(left[0], left[1], right[1]) }

process MergeBamAlignment {

	input:
	file gatk4
	file fasta_ref
	file fasta_ref_fai
	file fasta_ref_dict
	set pair_id, file(bam), file(sam) from sam_and_bam_ch

	output:
	set pair_id, file("mergebam.fastqtosam.bwa.bam") into MergeBamAlignment_output

	script:
	"""
	echo ${pair_id} ${bam} ${sam} && \
	java -Dsnappy.disable=true -Xmx16G -XX:ParallelGCThreads=16 -Djava.io.tmpdir=`pwd`/tmp -jar \
      $gatk4 \
      MergeBamAlignment \
      --VALIDATION_STRINGENCY SILENT \
      --EXPECTED_ORIENTATIONS FR \
      --ATTRIBUTES_TO_RETAIN X0 \
      --ALIGNED_BAM $sam \
      --UNMAPPED_BAM $bam \
      -O mergebam.fastqtosam.bwa.bam \
      --reference $fasta_ref \
      --SORT_ORDER coordinate \
      --IS_BISULFITE_SEQUENCE false \
      --ALIGNED_READS_ONLY false \
      --CLIP_ADAPTERS false \
      --MAX_RECORDS_IN_RAM 200000 \
      --ADD_MATE_CIGAR true \
      --MAX_INSERTIONS_OR_DELETIONS -1 \
      --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
      --PROGRAM_RECORD_ID "bwamem" \
      --PROGRAM_GROUP_VERSION "0.7.12-r1039" \
      --PROGRAM_GROUP_COMMAND_LINE "bwa mem -t 18 -R -M Input1 Input2 > output.sam" \
      --PROGRAM_GROUP_NAME "bwamem"
    """
}

process MarkDup {

	input:
	file gatk4
	file bam_files from MergeBamAlignment_output.collect()

	output:
	file("markduplicates.mergebam.fastqtosam.bwa.bam") into MarkDup_bamoutput
	file("markduplicates.mergebam.fastqtosam.bwa.bai") into MarkDup_baioutput
	file("markduplicates.metrics")

	script:
	def input_args = bam_files.collect{ "-I $it" }.join(" ")
	"java -Dsnappy.disable=true -Xmx16G -XX:ParallelGCThreads=16 -Djava.io.tmpdir=`pwd`/tmp -jar \
	  $gatk4 \
	  MarkDuplicates \
	  $input_args \
	  -O markduplicates.mergebam.fastqtosam.bwa.bam \
	  --VALIDATION_STRINGENCY LENIENT \
	  --METRICS_FILE markduplicates.metrics \
	  --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 200000 \
	  --CREATE_INDEX true"
}

process BaseRecalibrator {
	maxForks = 2

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
	each intervals from file(params.contigs).readLines().findAll{ it }
	file bam from MarkDup_bamoutput
	file bai from MarkDup_baioutput

	output:
	file("recalibration_report.grp") into BaseRecalibrator_output

	script:
    def interval = intervals.tokenize('\t').collect{ "-L $it" }.join(" ")
	"java -Dsnappy.disable=true -XX:ParallelGCThreads=4 \
	  -Dsamjdk.use_async_io=false -Xmx13G \
	  -jar $gatk4 \
	  BaseRecalibrator \
	  --reference $fasta_ref \
	  --input $bam \
	  -O recalibration_report.grp \
	  --knownSites $dbsnp_vcf \
	  --knownSites $v1000g_vcf \
	  --knownSites $mills_vcf \
	  $interval 2> info"
}

process GatherBqsrReports {
	
	input:
	file gatk4
	file input_bqsr_reports from BaseRecalibrator_output.collect()

	output:
	file("gathered_recalibration_reports.grp") into GatherBqsrReports_output

	script:
	def input_reports = input_bqsr_reports.collect{ "--input $it" }.join(" ")
	"java -Dsnappy.disable=true -Xmx16G -jar \
      $gatk4 \
      GatherBQSRReports \
      $input_reports \
      -O gathered_recalibration_reports.grp"
}

process ApplyBQSR {

	input:
	file gatk4
	file read from read1
	file recalibration_report from GatherBqsrReports_output
	each intervals from file(params.contigs).readLines().findAll{ it }
	file fasta_ref
	file fasta_ref_fai
	file fasta_ref_dict

	output:
	file("applybqsr.markdup.mergebam.fastqtosam.bwa.bam") into ApplyBQSR_bamoutput
	file("applybqsr.markdup.mergebam.fastqtosam.bwa.bai") into ApplyBQSR_baioutput

	shell:
    def interval = intervals.tokenize('\t').collect{ "-L $it" }.join(" ")
	"""
	java -Dsnappy.disable=true -Xmx13G -XX:ParallelGCThreads=4 \
	  -jar $gatk4 \
	  ApplyBQSR \
	  --reference $fasta_ref \
	  --input $read \
	  -O applybqsr.markdup.mergebam.fastqtosam.bwa.bam \
	  --createOutputBamIndex true \
	  -bqsr $recalibration_report \
	  $interval
	"""
}

process GatherBamFiles {
	input:
	file gatk4
	file input_bams from ApplyBQSR_bamoutput.collect()

	output:
	file("gatheredbams.applybqsr.baserecal.markdup.fastqtosam.bwa.bam") into GatherBamFiles_bamoutput
	file("gatheredbams.applybqsr.baserecal.markdup.fastqtosam.bwa.bai") into GatherBamFiles_baioutput

	script:
	"java -Dsnappy.disable=true -Xmx16G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=`pwd`/tmp -jar \
      $gatk4 \
      GatherBamFiles \
      --input $input_bams \
      -O gatheredbams.applybqsr.baserecal.markdup.fastqtosam.bwa.bam \
      --CREATE_INDEX true"
}

process HaplotypeCaller {
	input:
	file gatk3
	file bam from GatherBamFiles_bamoutput
	file bam_index from GatherBamFiles_baioutput
	file fasta_ref_dict
	file fasta_ref
	file fasta_ref_fai

	output:
	file("haplotypecaller.applybqsr.markdup.mergebam.fastqtosam.bwa.g.vcf") into HaplotypeCaller_vcf
	file("haplotypecaller.applybqsr.markdup.mergebam.fastqtosam.bwa.g.vcf.idx") into HaplotypeCaller_idx

	shell:
	"""
	java -XX:ParallelGCThreads=4 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx13G \
	  -jar $gatk3 \
	  -T HaplotypeCaller \
	  -nct 4 \
	  -R $fasta_ref \
	  -o haplotypecaller.applybqsr.markdup.mergebam.fastqtosam.bwa.g.vcf \
	  -I $bam \
	  -ERC GVCF
	"""
}

process GenotypeGVCF {
	input:
	file gatk3
	file vcf from HaplotypeCaller_vcf
	file vcf_index from HaplotypeCaller_idx
	file fasta_ref
	file fasta_ref_fai
	file fasta_ref_dict

	output:
	file("genotype.haplotypecaller.applybqsr.markdup.mergebam.fastqtosam.bwa.g.vcf") into GenotypeGVCF_vcfoutput
	file("genotype.haplotypecaller.applybqsr.markdup.mergebam.fastqtosam.bwa.g.vcf.idx") into GenotypeGVCF_idxoutput

	shell:
	"""
	java -Xmx16G -XX:ParallelGCThreads=18 -jar \
	$gatk3 \
	-T GenotypeGVCFs \
	-nt 4 \
	-R $fasta_ref \
	-o genotype.haplotypecaller.applybqsr.markdup.mergebam.fastqtosam.bwa.g.vcf \
	--variant $vcf
	"""
}

process VariantRecalibratorSNP {
	input:
	file gatk3
	file vcf from GenotypeGVCF_vcfoutput
	file idx from GenotypeGVCF_idxoutput
	file v1000g_vcf
	file omni_vcf
	file dbsnp_vcf
	file hapmap_vcf
	file v1000g_vcf_index
	file omni_vcf_index
	file dbsnp_vcf_index
	file hapmap_vcf_index
	file fasta_ref
	file fasta_ref_fai
	file fasta_ref_dict

	output:
	file("GVCF.SNP.recal") into VariantRecalibratorSNP_recal
	file("GVCF.SNP.tranches") into VariantRecalibratorSNP_tranches

	shell:
	"""
	java -Xmx16G -XX:ParallelGCThreads=18 -jar \
	$gatk3 \
	-T VariantRecalibrator \
	-nt 4 \
	-R $fasta_ref \
	-input $vcf \
	-mode SNP \
	-resource:v1000G,known=false,training=true,truth=false,prior=10.0 $v1000g_vcf \
	-resource:omni,known=false,training=true,truth=true,prior=12.0 $omni_vcf \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp_vcf \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap_vcf \
	-an QD -an MQ -an DP -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 \
	-tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 \
	-tranche 97.0 -tranche 90.0 \
	-recalFile GVCF.SNP.recal \
	-tranchesFile GVCF.SNP.tranches
	"""
}

process VariantRecalibratorINDEL {
	input:
	file gatk3
	file vcf from GenotypeGVCF_vcfoutput
	file idx from GenotypeGVCF_idxoutput
	file dbsnp_vcf
	file dbsnp_vcf_index
	file fasta_ref
	file fasta_ref_fai
	file fasta_ref_dict

	output:
	file("GVCF.INDEL.recal") into VariantRecalibratorINDEL_recal
	file("GVCF.INDEL.tranches") into VariantRecalibratorINDEL_tranches

	shell:
	"""
	java -Xmx16G -XX:ParallelGCThreads=18 -jar \
	$gatk3 \
	-T VariantRecalibrator \
	-nt 4 \
	-R $fasta_ref \
	-input $vcf \
	-mode INDEL \
	-resource:mills,known=false,training=true,truth=true,prior=12.0 $mills_vcf \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp_vcf \
	-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 \
	-tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 \
	-tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
	-recalFile GVCF.INDEL.recal \
	-tranchesFile GVCF.INDEL.tranches
	"""
}

process ApplyVQSRSNP {
	input:
	file gatk3
	file fasta_ref
	file fasta_ref_fai
	file fasta_ref_dict
	file vcf from HaplotypeCaller_vcf
	file recal from VariantRecalibratorSNP_recal
	file tranches from VariantRecalibratorSNP_tranches

	output:
	file("GVCF.SNP.g.vcf")
	file("GVCF.SNP.g.vcf.idx")

	shell:
	"""
	java -XX:ParallelGCThreads=18 -jar -Xmx16G \
	$gatk3 \
	-T ApplyRecalibration \
	-nt 4 \
	-input $vcf \
	-R $fasta_ref \
	-mode SNP \
	--ts_filter_level 99.0 \
	-tranchesFile $tranches \
	-recalFile $recal \
	-o GVCF.SNP.g.vcf
	"""
}

process ApplyVQSRINDEL {
	input:
	file gatk3
	file fasta_ref
	file fasta_ref_fai
	file fasta_ref_dict
	file vcf from HaplotypeCaller_vcf
	file recal from VariantRecalibratorINDEL_recal
	file tranches from VariantRecalibratorINDEL_tranches

	output:
	file("GVCF.INDEL.g.vcf")
	file("GVCF.INDEL.g.vcf.idx")

	shell:
	"""
	java -XX:ParallelGCThreads=18 -jar -Xmx16G \
	$gatk3 \
	-T ApplyRecalibration \
	-nt 4 \
	-input $vcf \
	-R $fasta_ref \
	-mode INDEL \
	--ts_filter_level 95.0 \
	-tranchesFile $tranches \
	-recalFile $recal \
	-o GVCF.INDEL.g.vcf
	"""
}