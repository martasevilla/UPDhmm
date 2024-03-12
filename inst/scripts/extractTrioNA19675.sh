#!/bin/bash
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrx")

sample_list="NA19675,NA19679,NA19678"

# Iterate over chromosomes
for chromosome in "${chromosomes[@]}"; do
    # Input VCF file
    input_vcf="ALL.${chromosome}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz"

    # Output VCF file
    output_vcf="NA19675_${chromosome}.vcf.gz"

    # Run bcftools command
    bcftools view -s "$sample_list" -Oz -o "$output_vcf" "$input_vcf"
done