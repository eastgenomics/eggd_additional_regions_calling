#!/bin/bash
# This script processes a BAM file and a list of genomic regions to generate VCF files.
# It uses bcftools to call variants and filters them based on a threshold.
# It merges the generated VCFs with a provided Sentieon VCF and uploads the final result.

_download_inputs() {
    mark-section "Downloading inputs"
    echo "Downloading input files..."
    dx-download-all-inputs --parallel
    # Move files to the correct directory
    find ~/in/ -type f -name "*" -print0 | xargs -0 -I {} mv {} /home/dnanexus
}

_extract_reference() {
    mark-section "Extracting reference"
    echo "Extracting reference archive..."
    tar -I pigz -xf "$fasta_tar_name"
    reference_fasta_name="genome.fa"
    echo "Reference extracted to genome.fa"
}

_extract_sample_info() {
    # Extract filename and sample name
    input_bam_prefix="${input_bam_name%.bam}" # Remove .bam extension
    sample_name="${input_bam_prefix%%_*}"     # Extract sample prefix

    echo "Extracted Sample Name: $sample_name"
}

_check_region_list_exists() {
    if [[ ! -f "$region_list_name" ]]; then
        echo "Error: File $region_list_name does not exist!" >&2
        exit 1
    fi
}

_index_vcf_if_missing() {
    local vcf_file="$1"
    if [[ "$vcf_file" == *.vcf.gz && ! -f "$vcf_file.tbi" ]]; then
        echo "Indexing VCF: $vcf_file"
        tabix -p vcf "$vcf_file"
    fi
}

_check_if_variant_exists_already() {
    local chromPos="$1"
    local knownRef="$2"
    local knownAlt="$3"
    local hits

    hits=$(bcftools view -r "$chromPos" "$sentieon_vcf_name" | bcftools query -f '%REF\t%ALT\n' | grep -P "^${knownRef}\t${knownAlt}$" || true)
    if [[ -n "$hits" ]]; then
        echo "Variant ${chromPos} with REF=${knownRef} and ALT=${knownAlt} already exists in the VCF."
        return 0
    else
        return 1
    fi
}

_generate_region_vcfs() {
    mark-section "Regional calling with bcftools"

    output_vcfs=()
    while IFS=$'\t' read -r chromPos knownRef knownAlt thresh; do
        echo "Processing region: ${chromPos} with REF=${knownRef} ALT=${knownAlt} threshold=${thresh}"
        output_vcf="${sample_name}_${chromPos//:/}_${knownRef}_${knownAlt}.vcf.gz"

        # Check if the variant already exists in the VCF
        if _check_if_variant_exists_already "$chromPos" "$knownRef" "$knownAlt"; then
            echo "Skipping region ${chromPos} as it already exists in the VCF."
            continue
        fi

        # Call variants using bcftools mpileup and filter based on the threshold
        # Use the -Ou option to output uncompressed BCF format
        # Use -a to specify the format fields to include in the output

        bcftools mpileup -d 8000 -f "$reference_fasta_name" "$input_bam_name" \
            -r "$chromPos" -a FORMAT/AD,FORMAT/DP -Ou | \
            bcftools call -mv -Oz -o "$output_vcf"

        bcftools view -i "REF=='${knownRef}' && ALT=='${knownAlt}'" -Oz -o "${output_vcf}_filtered.vcf.gz" "$output_vcf"
        mv "${output_vcf}_filtered.vcf.gz" "$output_vcf"
        _index_vcf_if_missing "$output_vcf"

        # Check if at least one variant in this region meets the threshold
        # The code checks if REF/ALT match what we expect (knownRef, knownAlt)
        # and then calculates ratio = ALT/REF. If ratio >= threshold,
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
                      if (a == 0) {
                        print "NA"     # No reads at all
                      } else {
                        print "INF"    # All reads are ALT
                      }
                    } else {
                      print a / r      # Normal case
                    }
                  }')

                # If AWK printed ERROR_DIV_ZERO, skip this variant or handle how you wish
                if [[ "$ratio" == "NA" ]]; then
                    echo "Warning: No reads at $chrom:$pos. Skipping variant." >&2
                    continue
                fi

                if [[ "$ratio" == "INF" ]]; then
                    pass=true
                    break
                fi
                # Check if ratio meets threshold
                # Use awk to compare ratio and threshold
                # awkPass will be 1 if ratio >= threshold, else 0
                awkPass=$(awk -v val="${ratio}" -v thr="${thresh}" 'BEGIN{if(val>=thr) print 1; else print 0}')
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

_merge_region_vcfs() {
    mark-section "Merging regional VCFs"
    # Merge all VCFs generated for additional regions into one file (merged_vcf).
    merged_vcf="${input_bam_prefix}_additional.vcf.gz"

    if [ ${#output_vcfs[@]} -gt 0 ]; then
        if [ ${#output_vcfs[@]} -gt 1 ]; then
            echo "Merging multiple additional region VCFs..."
            bcftools concat "${output_vcfs[@]}" -Oz -o "$merged_vcf"
            _index_vcf_if_missing "$merged_vcf"
        else
            # Only one VCF, rename it as the merged output
            merged_vcf="${output_vcfs[0]}"
            mv "$merged_vcf" "${input_bam_prefix}_additional.vcf.gz"
            merged_vcf="${input_bam_prefix}_additional.vcf.gz"
            _index_vcf_if_missing "$merged_vcf"
        fi
    else
        echo "No additional region VCF files were generated. Using sentieon VCF as final."
        merged_vcf="$sentieon_vcf_name"
    fi
}

_merge_with_sentieon_vcf() {
    mark-section "Merging regional VCF with Sentieon VCF"

    # Conditionally merge the merged_vcf with the sentieon VCF if output_combined is true.

    # If merged_vcf == sentieon_vcf_name, just reuse the original file
    if [[ "$merged_vcf" == "$sentieon_vcf_name" ]]; then
        echo "No region VCFs were added; reusing Sentieon VCF as final."
        final_vcf="$merged_vcf"
        return
    fi

    final_vcf="${sample_name}_additional_combined.vcf.gz"

    if [ "$output_combined" = true ]; then
        echo "Merging with sentieon VCF..."
        bcftools concat -a "$merged_vcf" "$sentieon_vcf_name" -Oz -o "$final_vcf"
    else
        echo "Skipping merging with sentieon VCF."
        final_vcf="$merged_vcf"
    fi
}

_normalize_vcf() {
    # Normalize the final VCF using bcftools norm.
    # The normalized file overwrites final_vcf to keep usage consistent.
    local temp_norm_vcf="${final_vcf%.vcf.gz}_normalized.vcf.gz"

    echo "Normalizing final VCF..."
    bcftools norm "$final_vcf" -f "$reference_fasta_name" -m -any --keep-sum AD -Oz -o "$temp_norm_vcf"
    tabix -p vcf "$temp_norm_vcf"

    final_vcf="$temp_norm_vcf"
}

_upload_final_vcf() {
    mark-section "Uploading outputs"

    # Upload final VCF to DNAnexus, record output
    uploaded_vcf=$(dx upload "$final_vcf" --brief)
    dx-jobutil-add-output output_vcf "$uploaded_vcf" --class=file

    # If merged_vcf differs from the original input VCF, upload it as an additional file
    if [ "$merged_vcf" != "$sentieon_vcf_name" ] && [ -f "$merged_vcf" ]; then
        uploaded_merged_vcf=$(dx upload "$merged_vcf" --brief)
        dx-jobutil-add-output merged_additional_regions_vcf "$uploaded_merged_vcf" --class=file
    else
        # Create an empty VCF to satisfy DNAnexus requirements
        echo -e "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > empty_merged.vcf
        bgzip empty_merged.vcf

        placeholder_vcf=$(dx upload empty_merged.vcf.gz --brief)
        dx-jobutil-add-output merged_additional_regions_vcf "$placeholder_vcf" --class=file
    fi
}

main() {
    set -exo pipefail # Strict mode to catch errors

    _download_inputs
    _extract_sample_info
    _extract_reference
    _check_region_list_exists
    _index_vcf_if_missing "$sentieon_vcf_name"
    _generate_region_vcfs
    _merge_region_vcfs
    _merge_with_sentieon_vcf
    _normalize_vcf
    _upload_final_vcf
    echo "Done!"
    mark-success
}
