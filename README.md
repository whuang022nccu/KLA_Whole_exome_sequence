# Kaposiform Lymphangiomatosis  Whole Exome s equence Analysis pipeline

WES pipeline for rerun a Kaposiform Lymphangiomatosis reserarch paper
:

[A Somatic Activating NRAS Variant Associated with Kaposiform Lymphangiomatosis](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6565516/)

the source code was refactor from the [supplementary file](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6565516/bin/NIHMS1521476-supplement-Supplementary_File.pdf)

### Setup

```
git clone https://github.com/whuang022nccu/KLA_Whole_exome_sequence

```
and change your PICARD , BWA and GATK paths :

* picard_jar_path="your_path_to/picard-2.8.3.jar"
* bwa_path='your_path_to/bwa-0.7.15'
* gatk_path='your_path_to/gatk-4.0.12.0/gatk'
* GenomeAnalysisTK_jar_path="your_path_to/GenomeAnalysisTK.jar"

and the data related :

* ref_genome_path='your_path_to/Homo_sapiens_assembly38.fasta'
* S04380110_Covered='your_path_to/S04380110_Covered.bed'
* dbsnp_138_hg38_vcf='your_path_to/dbsnp_138.hg38.vcf.gz'
* Mills_and_1000G_gold_standard_indels_h g38_vcf='your_path_to/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
* af_only_gnomad_vcf_path="your_path_to/af-only-gnomad.hg38.vcf.gz"
* S_1000g_pon_hg38_vcf="your_path_to/1000g_pon.hg38.vcf.gz"


### Useage

change the sample path  and sample ID :

* sample_ID="your-ID"
* sample_path="your/path"

move your raw fastq data in your/path/raw

and run :
```
./wes_hg38_kla.sh

```