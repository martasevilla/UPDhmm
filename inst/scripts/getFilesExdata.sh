############################################################
############################################################
#1. Creation of test_het_mat.vcf.gz
#File test_het_mat.vcf.gz is comes from GIB cohort data
############################################################
############################################################

#1. Download chromosome files from ftp site:

wget -r ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/integrated_call_sets/


#2. Extract from the files only the samples from the trio "NA19675,NA19679,NA19678"

sh extracTrioNA19675.sh


#3. Concatenate the chromosome files from the trio:

ls NA19675_*vcf.gz > file_list_NA19675

bcftools concat -f file_list_NA19675 -o trio_NA19675.vcf.gz -Oz


#4. Perform a serie of filtering steps in order to reduce the size of the file:

#Extract only exome calling
bcftools view -i 'SNPSOURCE="EXOME"' -Oz -o exome_trio_NA19675.vcf.gz trio_NA19675.vcf.gz

#Exract only biallelic sites
bcftools view -m2 -M2 -Oz -o exome_trio_NA19675_def.vcf.gz exome_trio_NA19675.vcf.gz

#Normalize vcf
bcftools norm -d all -o sort_trio_NA19675.vcf.gz exome_trio_NA19675_def.vcf.gz

#Remove commons deletions and duplications, as well as HLA and KIR regions
bcftools view -T cnps.bed -Oz -o filter_cnp_trio_NA19675.vcf.gz sort_trio_NA19675.vcf.gz

#Remove GT 0/0 ,AF > 0.1 and chrX to minimize file's size
bcftools view -i 'GT!="0|0"' -Oz -o filter_hom_trio_NA19675.vcf.gz filter_cnp_trio_NA19675.vcf.gz
bcftools view -i 'EUR_AF>0.1' -Oz -o filter_hom_AF_trio_NA19675.vcf.gz filter_hom_trio_NA19675.vcf.gz
bcftools view -i 'CHROM!="X"' -Oz -o filter_hom_AF_autosomal_trio_NA19675.vcf.gz filter_hom_AF_trio_NA19675.vcf.gz


# * Note cnps.bed is a bed file downloaded from GenomeBrowser with commons duplications and deletions, #is also provided in inst/scripts

#5. Once the filtering is done, a simulation in chromosome 3q is done. For that, the simulation script #is the following: 

#The scripts is called simulation_3q.R and is stored in inst/script


Rscript simulation_3q.R --input filter_hom_AF_autosomal_trio_NA19675.vcf.gz --region 3 --pos1 96468334 --pos2 198022430 --type heterodisomy --parent mother --out test_het_mat.vcf

#6. Finally, the file is compressed

bgzip test_het_mat.vcf


############################################################
############################################################
2. Creation of test.gz
############################################################
############################################################

#Extract only some positions for the performing of tests function

bcftools view -r 6:32489853-32489876 test_het_mat.vcf.gz -Oz -o test.vcf.gz