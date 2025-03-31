#!/bin/bash

main() {
    set -exuo pipefail  # Strict mode to catch errors

    echo "Inputs specified:"
    echo "Value of input bam: '$input_bam'"
    echo "Value of input bai: '$input_bai'"
    echo "Value of reference fasta: '$reference_fasta'"
    echo "Value of reference fai: '$reference_fai'"
    echo "Value of region list: '$region_list'"


    echo "Downloading input files..."
    dx-download-all-inputs --parallel

    # Extract filename and sample name
    input_bam_name=$(basename "$input_bam_path")  # Get the actual filename
    input_bam_prefix="${input_bam_name%.bam}"  # Remove .bam extension
    sample_name="${input_bam_prefix%%_*}"  # Extract sample prefix
    find ~/in/ -type f -name "*" -print0 | xargs -0 -I {} mv {} /home/dnanexus

    echo "Extracted Sample Name: $sample_name"

    if [ -n "$senteion_vcf" ]; then
        dx download "$senteion_vcf" -o senteion.vcf.gz
    fi

    if [ -n "$senteion_tbi" ]; then
        dx download "$senteion_tbi" -o sentieon.tbi
    fi

    if [ -n "$region_list" ]; then
        dx download "$region_list" -o regions.txt
        region_list=$(cat regions.txt)
    fi

    output_vcfs=()
    for region in $region_list; do
        echo "Processing region: $region"
        output_vcf="output_${region//:/_}.vcf.gz"

        bcftools mpileup -d 8000 -f "$reference_fasta_name" "$input_bam_name" -r "$region" -a FORMAT/AD,FORMAT/DP -Ou | \
        bcftools call -mv -Oz -o "$output_vcf"

        tabix -p vcf "$output_vcf"
        output_vcfs+=("$output_vcf")
    done
    # merge VCFs if multiple regions
    final_vcf="${input_bam_prefix}.vcf.gz"
    if [ ${#output_vcfs[@]} -gt 1 ]; then
        echo "Merging VCFs..."
        bcftools concat "${output_vcfs[@]}" -Oz -o "$final_vcf"
        # Normalise VCF
        bcftools norm "$final_vcf" -f "$reference_fasta_name" -m -any --keep-sum AD -Oz -o "$final_vcf"
        tabix -p vcf "$final_vcf"
    else
        # Normalise VCF
        bcftools norm "$output_vcf" -f "$reference_fasta_name" -m -any --keep-sum AD -Oz -o "$final_vcf"
    fi


    # Merge with senteion VCF
    # if [ -n "$senteion_vcf" ]; then
    #     bcftools merge -m none -Oz -o merged.vcf.gz "$final_vcf" "$senteion_vcf"
    #     final_vcf="${sample_name}_merged.vcf.gz"
    # fi

    # Upload final VCF
    uploaded_vcf=$(dx upload "$final_vcf" --brief)

    dx-jobutil-add-output output_vcf "$uploaded_vcf" --class=file
}