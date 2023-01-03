#!/bin/bash
# refactor WES analysis pipline from paper 
# A Somatic Activating NRAS Variant Associated with Kaposiform Lymphangiomatosis
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6565516/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6565516/bin/NIHMS1521476-supplement-Supplementary_File.pdf
# dependencies path setting
picard_jar_path="/media/whuang022/DATA/KLA_datapipline/picard-2.8.3.jar"
bwa_path='/media/whuang022/DATA/KLA_datapipline/bwa-0.7.15'
gatk_path='/media/whuang022/DATA/KLA_datapipline/gatk-4.0.12.0/gatk'
GenomeAnalysisTK_jar_path="/media/whuang022/DATA/KLA_datapipline/GenomeAnalysisTK.jar"
# ref data path setting
ref_genome_path='/media/whuang022/DATA/KLA_datapipline/refactor_pipline/ref_data_hg38/Homo_sapiens_assembly38.fasta'
S04380110_Covered='/media/whuang022/DATA/KLA_datapipline/refactor_pipline/ref_data_hg38/S04380110_Covered.bed'
dbsnp_138_hg38_vcf='/media/whuang022/DATA/KLA_datapipline/refactor_pipline/ref_data_hg38/dbsnp_138.hg38.vcf.gz'
Mills_and_1000G_gold_standard_indels_hg38_vcf='/media/whuang022/DATA/KLA_datapipline/refactor_pipline/ref_data_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
af_only_gnomad_vcf_path="/media/whuang022/DATA/KLA_datapipline/refactor_pipline/ref_data_hg38/af-only-gnomad.hg38.vcf.gz"
S_1000g_pon_hg38_vcf="/media/whuang022/DATA/KLA_datapipline/refactor_pipline/ref_data_hg38/1000g_pon.hg38.vcf.gz"
year=$(date +%Y)
#sample_ID="father"
#sample_path="/media/whuang022/DATA/KLA_datapipline/wes_hg38/father"
sample_ID="your-ID"
sample_path="your/path"
tumor_n=$sample_ID
out_path=$sample_path # set output_path same as sample path
# folder setting
raw_folder_name='raw'
tmp_folder_name='tmp'
bam__folder_name='bam'
pdf_folder_name='pdf'
stat_folder_name='qc'
vcf_folder_name='vcf'
# path prefix setting
bam_prefix=$out_path/$bam__folder_name
pdf_prefix=$out_path/$pdf_folder_name
stat_prefix=$out_path/$stat_folder_name
raw_prefix=$out_path/$raw_folder_name
vcf_prefix=$out_path/$vcf_folder_name
tmp_path=$out_path/$tmp_folder_name
# raw data file name setting 
R1_postfix='_R1.fq'
R2_postfix='_R2.fq'
read_group_name='739.2'
platfoorm='illumina'
library_name='test'
sequencing_center='test'
# machine setting
java_max_mem="-Xmx50G"
n_cpu=16
min_base_quality_score=20
mbqs=$min_base_quality_score

echo "==========WES analysis =========="
echo "pipeline start ......"
echo "sample-ID = $sample_ID"
 

## make tmp 
mkdir $tmp_path

## check sample path
if [ -d "$sample_path" ];
then
    echo "current sample path  $sample_path "
else
    echo "sample path $sample_path not exist , please check it ! "
    exit 1
fi
## check raw fastq 

FASTQ_R1_path=$raw_prefix/$sample_ID$R1_postfix
FASTQ_R2_path=$raw_prefix/$sample_ID$R2_postfix 

if [[ -f $FASTQ_R1_path && $FASTQ_R2_path ]];
then
    echo "decated raw fastq data:"
    echo "R1 = $FASTQ_R1_path"
    echo "R2 = $FASTQ_R1_path"
else
    echo "raw fastq not exist! please check it!"
    exit 1
fi
## creat QC file path
if [ -d "$stat_prefix" ];
then
    echo "$stat_prefix  QC stat file paths already exists! please delete it and re run.( make sure the files we do not need any more)."
    exit 1
else
    mkdir $stat_prefix
fi

## creat chart pdf file path
if [ -d "$pdf_prefix" ];
then
    echo "$pdf_prefix chart path already exists! please delete it and re run.( make sure the files we do not need any more)."
    exit 1
else
    mkdir $pdf_prefix
fi
#
## creat bam file path
if [ -d "$bam_prefix" ];
then
    echo "$bam_prefix BAM path already exists! please delete it and re run.( make sure the files we do not need any more)."
    exit 1
else
    mkdir $bam_prefix
fi

## creat vcf file path
if [ -d "$vcf_prefix" ];
then
    echo "$vcf_prefix VCF path already exists! please delete it and re run.( make sure the files we do not need any more)."
    exit 1
else
    mkdir $vcf_prefix
fi


## step1 fastq to unmapped bam
echo "step1 fastq to unmapped bam (using PICARD)"
java $java_max_mem -jar $picard_jar_path FastqToSam \
FASTQ=$FASTQ_R1_path \
FASTQ2=$FASTQ_R2_path \
OUTPUT=$bam_prefix/$sample_ID.unmapped.bam \
READ_GROUP_NAME=$read_group_name \
SAMPLE_NAME=$sample_ID \
LIBRARY_NAME=$library_name \
PLATFORM=$platfoorm \
SEQUENCING_CENTER=$sequencing_center \
RUN_DATE="$year" \
TMP_DIR=$tmp_path
## step2 mark illumina adaptors
echo "step2 unmapped bam mark illumina adaptors "
java $java_max_mem -jar $picard_jar_path MarkIlluminaAdapters  \
I=$bam_prefix/$sample_ID.unmapped.bam  \
O=$bam_prefix/$sample_ID.markilluminaadapters.unmapped.bam  \
M=$stat_prefix/$sample_ID.markilluminaadapters_metrics.txt  \
TMP_DIR=$tmp_path
## step3 use gatk pipeline (with bwa) mapping
echo "step3 use gatk pipeline (with bwa) to generate clean, mapped bam"
java $java_max_mem -jar $picard_jar_path SamToFastq \
I=$bam_prefix/$sample_ID.markilluminaadapters.unmapped.bam \
FASTQ=/dev/stdout \
CLIPPING_ATTRIBUTE=XT \
CLIPPING_ACTION=2 \
INTERLEAVE=true \
NON_PF=true \
TMP_DIR=$tmp_path \
| \
$bwa_path mem \
-M \
-t $n_cpu \
-p $ref_genome_path \
/dev/stdin \
| \
java $java_max_mem -jar $picard_jar_path MergeBamAlignment \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM=$bam_prefix/$sample_ID.markilluminaadapters.unmapped.bam \
OUTPUT=$bam_prefix/$sample_ID.bwa.bam \
R=$ref_genome_path \
CREATE_INDEX=true \
ADD_MATE_CIGAR=true \
CLIP_ADAPTERS=false \
CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true \
MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
ATTRIBUTES_TO_RETAIN=XS \
TMP_DIR=$tmp_path
## step4 markduplicates
echo "step4 use gatk pipeline markduplicates reads"
$gatk_path --java-options $java_max_mem SortSam \
--INPUT $bam_prefix/$sample_ID.bwa.bam \
--OUTPUT $bam_prefix/$sample_ID.bwa.query_sort.bam \
--SORT_ORDER "queryname" \
--CREATE_INDEX true \
--CREATE_MD5_FILE true \
--TMP_DIR $tmp_path
$gatk_path --java-options $java_max_mem MarkDuplicates \
-I $bam_prefix/$sample_ID.bwa.query_sort.bam \
-M $stat_prefix/$sample_ID.Duplicate_Metrics \
-O $bam_prefix/$sample_ID.bwa.dups_marked.bam \
--VALIDATION_STRINGENCY SILENT \
--ASSUME_SORT_ORDER "queryname" \
--CREATE_MD5_FILE true \
--CREATE_INDEX true \
--TMP_DIR $tmp_path
$gatk_path --java-options $java_max_mem SortSam \
--INPUT $bam_prefix/$sample_ID.bwa.dups_marked.bam \
--OUTPUT $bam_prefix/$sample_ID.bwa.dups_marked.sorted.bam \
--SORT_ORDER "coordinate" \
--CREATE_INDEX true \
--CREATE_MD5_FILE true \
--TMP_DIR $tmp_path
## step5 QC by BQSR ( Base Quality Score Recalibration) score
# Realign base quality scores
echo "step5 qc BQSR"
$gatk_path --java-options $java_max_mem SortSam \
--INPUT $bam_prefix/$sample_ID.bwa.bam \
--OUTPUT $bam_prefix/$sample_ID.bwa.query_sort.bam \
--SORT_ORDER "queryname" \
--CREATE_INDEX true \
--CREATE_MD5_FILE true\
--TMP_DIR $tmp_path
$gatk_path --java-options $java_max_mem MarkDuplicates \
-I $bam_prefix/$sample_ID.bwa.query_sort.bam \
-M $stat_prefix/$sample_ID.Duplicate_Metrics \
-O $bam_prefix/$sample_ID.bwa.dups_marked.bam \
--VALIDATION_STRINGENCY SILENT \
--ASSUME_SORT_ORDER "queryname" \
--CREATE_MD5_FILE true \
--CREATE_INDEX true \
--TMP_DIR $tmp_path
$gatk_path --java-options $java_max_mem SortSam \
--INPUT $bam_prefix/$sample_ID.bwa.dups_marked.bam \
--OUTPUT $bam_prefix/$sample_ID.bwa.dups_marked.sorted.bam \
--SORT_ORDER "coordinate" \
--CREATE_INDEX true \
--CREATE_MD5_FILE true\
--TMP_DIR $tmp_path
$gatk_path --java-options $java_max_mem BaseRecalibrator \
-R $ref_genome_path \
-I $bam_prefix/$sample_ID.bwa.dups_marked.sorted.bam \
--use-original-qualities \
-O $stat_prefix/$sample_ID.bwa.recal_data.table \
--known-sites $dbsnp_138_hg38_vcf \
--known-sites $Mills_and_1000G_gold_standard_indels_hg38_vcf \
-L $S04380110_Covered \
--use-original-qualities
$gatk_path --java-options $java_max_mem ApplyBQSR \
-R $ref_genome_path \
-I $bam_prefix/$sample_ID.bwa.dups_marked.sorted.bam \
-O $bam_prefix/$sample_ID.bwa.recal.bam \
-L $S04380110_Covered \
-bqsr $stat_prefix/$sample_ID.bwa.recal_data.table \
--add-output-sam-program-record \
--create-output-bam-md5 \
--use-original-qualities \
--create-output-bam-index
## step6 QC by QSD ( Quality Score Distribution)
# Determine quality score distribution (QSD) before and after BQSR using qualityscoredistribution;
# Determine mismatch rate using collectalignmentsummarymetrics; Calculate depth of coverage
# using depthofcoverage (GATK)
echo "step6 QC by QSD & depthofcoverage"
$gatk_path --java-options $java_max_mem QualityScoreDistribution \
-I $bam_prefix/$sample_ID.bwa.dups_marked.sorted.bam \
-CHART $pdf_prefix/$sample_ID.before.Chart.pdf \
-O $stat_prefix/$sample_ID.before.QSD.txt \
--ALIGNED_READS_ONLY true
$gatk_path --java-options $java_max_mem QualityScoreDistribution \
-I $bam_prefix/$sample_ID.bwa.recal.bam \
-CHART $pdf_prefix/$sample_ID.recal.Chart.pdf \
-O $stat_prefix/$sample_ID.recal.QSD.txt \
--ALIGNED_READS_ONLY true
$gatk_path --java-options $java_max_mem CollectAlignmentSummaryMetrics \
-I $bam_prefix/$sample_ID.bwa.recal.bam \
-O $stat_prefix/$sample_ID.bwa.recal.metrics.txt \
-R $ref_genome_path
java -jar $java_max_mem $GenomeAnalysisTK_jar_path \
-T DepthOfCoverage \
-R $ref_genome_path \
-o $bam_prefix/$sample_ID.exome_depth_Q20 \
-I $bam_prefix/$sample_ID.bwa.recal.bam \
-L $S04380110_Covered \
-mbq $mbqs\
-omitBaseOutput \
-ct 4 \
-ct 8 \
-ct 15 \
-ct 20 \
-ct 40 \
-ct 50 \
-ct 80 \
-ct 100 \
-ct 120 \
-ct 150 \
-ct 200 \
-ct 250 \
-ct 300
## step7 Mutect2 call vcf
# Tumour-only exomes Mutect2 call vcf
echo "step7 Mutect2 call vcf Tumour-only mode"
$gatk_path --java-options $java_max_mem Mutect2 \
-R $ref_genome_path \
-I $bam_prefix/$sample_ID.bwa.recal.bam  \
-tumor $tumor_n \
-pon $S_1000g_pon_hg38_vcf \
--germline-resource $af_only_gnomad_vcf_path \
--af-of-alleles-not-in-resource 0.0000025 \
-L $S04380110_Covered \
-O $vcf_prefix/$sample_ID.m2.vcf.gz \
-bamout $bam_prefix/$sample_ID.m2.bam

echo "pipeline done !"