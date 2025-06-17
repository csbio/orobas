# Orobas

Orobas is an R package for scoring chemogenomic CRISPR screening data.

## Installation

Orobas is compatible with R version >= 3.6, and should work on any operating system capable of installing the packages listed in DESCRIPTION.

Install Orobas directly from Github with the `install_github("HenryWard/Orobas")` command from the devtools package. This should take no more than a minute or two.

## Usage

To learn how to use the package, please follow the vignette in the vignettes directory.

## Misc

Questions, comments or concerns can be directed to henry.neil.ward@gmail.com. Feedback regarding potential new features, bugs, or ease of use are especially welcome. 

## Expected Outputs

Output directory and files from single-screen scoring:


   .
    ├── <output>
    │   ├── <screen-batch-1>                    
    │   │   ├── ... 
    │   ├── <screen-batch-2>
    │   │   ├── qc                                       
    │   │   │   ├── essential_PR_QC.tsv
    │   │   │   ├── lfc_heatmap.png
    │   │   │   ├── replicate_cor.tsv
    │   │   │   ├── <screen-batch-2>_<condition-screen-1-replicate-A>_vs_<screen-batch-2>_<condition-screen-1-replicate-B>_replicate_comparison.png
    │   │   │   ├── ...
    │   │   │   ├── <screen-batch-2>_<control-screen-1-replicate-A>_vs_<screen-batch-2>_<control-screen-1-replicate-B>_replicate_comparison.png
    │   │   │   ├── ...
    │   │   │   ├── reads
    │   │   │   │   ├── total_reads.png
    │   │   │   │   ├── reads_heatmap.png
    │   │   │   │   ├── <screen-T0-1>_raw_reads_histogram.png
    │   │   │   │   ├── ...
    │   │   │   │   ├── <screen-batch-2>_<control-screen-1-replicate-A>_raw_reads_histogram.png
    │   │   │   │   ├── ...
    │   │   │   │   ├── <screen-batch-2>_<condition-screen-1-replicate-A>_raw_reads_histogram.png
    │   │   │   │   ├── ...
    │   │   ├── guide_dlfc
    │   │   │   ├── <screen-batch-2>_<condition-screen-1>_vs_<control-screen-1>_guide_dlfc_pre_jk.tsv
    │   │   │   ├── ...
    │   │   ├── plots
    │   │   │   ├── <screen-batch-2>_<condition-screen-1>_vs_<control-screen-1>_scatter.png
    │   │   │   ├── ...
    │   │   ├── condition_gene_calls.tsv
    │   │   ├── t0_normalized_screens_guide_level.tsv
    │   ├──  <screen-batch-3>
    │   │   ├── ... 
    │   ├── ... 
    │   ├── differential_LFC_scores.tsv      # outputs are created after running R scripts
    │   ├── fdr_scores.tsv
    

