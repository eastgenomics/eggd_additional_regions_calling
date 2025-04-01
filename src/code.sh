#!/bin/bash
# TODO: Change to modular functions

main() {
    set -exuo pipefail  # Strict mode to catch errors

    echo "Downloading input files..."
    dx-download-all-inputs --parallel

    # Extract filename and sample name
    input_bam_prefix="${input_bam_name%.bam}"  # Remove .bam extension
    sample_name="${input_bam_prefix%%_*}"  # Extract sample prefix
    find ~/in/ -type f -name "*" -print0 | xargs -0 -I {} mv {} /home/dnanexus

    echo "Extracted Sample Name: $sample_name"

    if [[ ! -f "$region_list_name" ]]; then
    echo "Error: File $region_list_name does not exist!" >&2
    exit 1
    fi

    output_vcfs=()
    region_list=$(cat "$region_list_name")
    for region in $region_list; do
        echo "Processing region: $region"
        output_vcf="output_${region//:/_}.vcf.gz"

        bcftools mpileup -d 8000 -f "$reference_fasta_name" "$input_bam_name" -r "$region" -a FORMAT/AD,FORMAT/DP -Ou | \
        bcftools call -mv -Oz -o "$output_vcf"

        tabix -p vcf "$output_vcf"
        output_vcfs+=("$output_vcf")
    done

    # Merge all region VCFs
    merged_vcf="${input_bam_prefix}_additional.vcf.gz"
    norm_vcf="${input_bam_prefix}_additional_normalized.vcf.gz"
    if [ ${#output_vcfs[@]} -gt 0 ]; then
        if [ ${#output_vcfs[@]} -gt 1 ]; then
            echo "Merging VCFs..."
            bcftools concat "${output_vcfs[@]}" -Oz -o "$merged_vcf"
        else
            # If only one VCF file, use it as final VCF
            merged_vcf="${output_vcfs[0]}"
            # Rename final VCF
            mv "$merged_vcf" "${input_bam_prefix}_additional.vcf.gz"
            merged_vcf="${input_bam_prefix}_additional.vcf.gz"
        fi

        # Normalize VCF
        bcftools norm "$merged_vcf" -f "$reference_fasta_name" -m -any --keep-sum AD -Oz -o "$norm_vcf"
        tabix -p vcf "$norm_vcf"
    else
        echo "Error: No VCF files were generated." >&2
        exit 1
    fi

    # Merge with senteion VCF, TODO: add concat for sention VCF
    # IF output_combined set as true
    final_vcf="${sample_name}_additional_normalised_combined.vcf.gz"
    if [ "$output_combined" = true ]; then
        bcftools concat -a "$norm_vcf" "$senteion_vcf_name" -Oz -o "${final_vcf}"
    else
        echo "Skipping merging with senteion SNV VCF"
        mv "$norm_vcf" "$final_vcf"
    fi

    # Upload final VCF
    uploaded_vcf=$(dx upload "$final_vcf" --brief)

    dx-jobutil-add-output output_vcf "$uploaded_vcf" --class=file
}