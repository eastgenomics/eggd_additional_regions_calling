<!-- dx-header -->
# eggd_additional_region_calling (DNAnexus Platform App)

## What does this app do?
This app performs additional variant calling on genomic regions that typically have low coverage but contain known variants needing closer evaluation. It works alongside an existing Sentieon variant call set by generating separate region-based VCFs, filtering them on a user-defined threshold, and merging everything into a final VCF.

## What are typical use cases for this app?
- Investigating low-coverage regions that are missed in a standard run to include known variants.

## What data are required for this app to run?
- A BAM file to call from (input_bam)
- A reference FASTA-index tar (fasta_tar)
- A list of regions, each with an expected REF, ALT, and threshold (region_list)
- A Sentieon VCF file (senteion_vcf) to merge final results with

## What are the optional inputs for this app?
- A boolean flag (output_combined) to specify whether to merge the new region calls with the Sentieon VCF. Defaults to true.

## How does the app work?
1. Reads the region list, each line specifying chrom:pos, known REF, known ALT, and a threshold.
2. Calls variants using bcftools mpileup/call on each region.
3. Evaluates whether REF and ALT in the called variants match the known REF/ALT.
4. Calculates the ratio of ALT to REF reads. If it meets or exceeds the threshold, the region passes.
5. Merges all passing region VCFs.
6. Optionally merges these with the Sentieon VCF.
7. Normalizes and indexes the final VCF.
8. Uploads the final VCF to the DNAnexus project.

## What does this app output?
- A final VCF (and index) containing the additional region calls that meet the defined threshold. When merging is enabled, it also includes variants from your Sentieon VCF.
- The merged additional regions VCF is uploaded to the DNAnexus project.

## Notes
- The script uses bcftools for variant calling, merging, and normalization.
- Shell scripts in this repository handle tasks like error checking and skipping regions that fail.

<!-- /dx-header -->