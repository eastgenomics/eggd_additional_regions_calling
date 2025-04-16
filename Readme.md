# eggd_additional_region_calling (DNAnexus Platform App)

## What does this app do?

This app calls variants on specified low-coverage genomic regions with known REF/ALT alleles and merges the results with an existing Sentieon VCF. It now also checks if a variant already exists in the Sentieon VCF to avoid duplicating calls. Finally, it normalizes and indexes the resulting VCF.

## What are typical use cases for this app?

- Investigating known low-coverage sites in a BAM file to confirm important known variants.

## What data are required for this app to run?

- A BAM file (input_bam) and its index (input_bai) for the region calls.
- A reference FASTA tar (fasta_tar) that includes genome.fa and genome.fa.fai.
- A region listing (region_list) specifying chrom:pos, REF, ALT, and threshold.
- A Sentieon VCF (sentieon_vcf) for merging and duplicate-variant checks.

## What are the optional inputs for this app?

- A boolean flag (output_combined) to merge these new region calls with the Sentieon VCF (defaults to true).

## How does the app work?

1. Reads each region, skipping any variants already found in the Sentieon VCF.
2. Uses bcftools mpileup/call to generate VCF files for each region.
3. Checks REF, ALT, and calculates an ALT-to-REF ratio to filter on a user-defined threshold.
4. Merges passing variants into a single VCF.
5. Optionally concatenates the merged file with the Sentieon VCF.
6. Normalizes the final VCF and uploads it along with any additional region VCFs.

## What does this app output?

- A final normalized VCF (and index) containing all regions that exceeded the set threshold and any variants from the Sentieon VCF (if merging is enabled).
- A merged additional_regions VCF (if new regions were identified).

## Notes

- bcftools mpileup/call is used for variant calling on each region.
- The script automatically indexes new and merged VCFs if indexes are missing.
- No new call is made if the variant is already in the Sentieon VCF.