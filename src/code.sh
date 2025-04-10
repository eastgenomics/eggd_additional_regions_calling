#!/bin/bash
# TODO: Change to modular functions

download_inputs() {
    echo "Downloading input files..."
    dx-download-all-inputs --parallel
}

extract_sample_info() {
    # Extract filename and sample name
    input_bam_prefix="${input_bam_name%.bam}"  # Remove .bam extension
    sample_name="${input_bam_prefix%%_*}"      # Extract sample prefix
    find ~/in/ -type f -name "*" -print0 | xargs -0 -I {} mv {} /home/dnanexus

    echo "Extracted Sample Name: $sample_name"
}

check_region_list_exists() {
    if [[ ! -f "$region_list_name" ]]; then
        echo "Error: File $region_list_name does not exist!" >&2
        exit 1
    fi
}


generate_region_vcfs() {
    output_vcfs=()
    # Expect the region_list file to have lines: chromPos  REF  ALT  threshold
    while IFS=$'\t' read -r chromPos knownRef knownAlt thresh; do
        echo "Processing region: $chromPos with REF=$knownRef ALT=$knownAlt threshold=$thresh"
        output_vcf="output_${chromPos//:/_}.vcf.gz"

        bcftools mpileup -d 8000 -f "$reference_fasta_name" "$input_bam_name" \
            -r "$chromPos" -a FORMAT/AD,FORMAT/DP -Ou | \
        bcftools call -mv -Oz -o "$output_vcf"

        tabix -p vcf "$output_vcf"

        # Check if at least one variant in this region meets the threshold
        # The code checks if REF/ALT match what we expect (knownRef, knownAlt)
        # and then calculates ratio = ALT/(REF+0.0001). If ratio >= threshold,
        # keep (pass=true). Otherwise, ignore.
        pass=false
        while read -r chrom pos ref alt ad; do
            # Split AD
            IFS=',' read -r refCount altCount <<< "$ad"
            refCount="${refCount:-0}"
            altCount="${altCount:-0}"

            if [[ "$ref" == "$knownRef" && "$alt" == "$knownAlt" ]]; then
                # Attempt ratio; if refCount=0, return an error line.
                ratio=$(awk -v r="$refCount" -v a="$altCount" '
                BEGIN {
                  if (r == 0) {
                    print "ERROR_DIV_ZERO"
                    exit
                  }
                  print a/r
                }')

                # If AWK printed ERROR_DIV_ZERO, skip this variant or handle how you wish
                if [[ "$ratio" == "ERROR_DIV_ZERO" ]]; then
                    echo "Warning: Division by zero for $chrom:$pos. Skipping variant." >&2
                    continue
                fi
                # Check if ratio meets threshold
                # Use awk to compare ratio and threshold
                # awkPass will be 1 if ratio >= threshold, else 0
                awkPass=$(awk -v val="$ratio" -v thr="$thresh" 'BEGIN{if(val>=thr) print 1; else print 0}')
                if [[ "$awkPass" == "1" ]]; then
                    pass=true
                    break
                fi
            fi
        done < <(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' "$output_vcf")

        # If any variant in this region passes the threshold, add it to the list
        if [[ "$pass" == true ]]; then
            echo "Region $chromPos meets threshold; keeping VCF."
            output_vcfs+=("$output_vcf")
        else
            echo "Region $chromPos below threshold; ignoring VCF."
            rm -f "$output_vcf" "$output_vcf.tbi"
        fi

    done < "$region_list_name"
}

merge_normalize_vcfs() {
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
}

merge_with_senteion_vcf() {
    # Merge with senteion VCF, TODO: add concat for sention VCF
    # IF output_combined set as true
    final_vcf="${sample_name}_additional_normalised_combined.vcf.gz"
    if [ "$output_combined" = true ]; then
        bcftools concat -a "$norm_vcf" "$senteion_vcf_name" -Oz -o "${final_vcf}"
    else
        echo "Skipping merging with senteion SNV VCF"
        mv "$norm_vcf" "$final_vcf"
    fi
}

upload_final_vcf() {
    # Upload final VCF
    uploaded_vcf=$(dx upload "$final_vcf" --brief)
    dx-jobutil-add-output output_vcf "$uploaded_vcf" --class=file
}

main() {
    set -exuo pipefail  # Strict mode to catch errors

    download_inputs
    extract_sample_info
    check_region_list_exists
    generate_region_vcfs
    merge_normalize_vcfs
    merge_with_senteion_vcf
    upload_final_vcf
}

main "$@"