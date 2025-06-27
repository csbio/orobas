# Orobas

Orobas is an R package for scoring chemogenomic CRISPR screening data.

## Overview

## Installation

## How to run

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
    │   │   │   ├── <screen-batch-2>_<condition-screen-1>_vs_<screen-batch-2>_<control-screen-1>_guide_dlfc_pre_jk.tsv
    │   │   │   ├── ...
    │   │   ├── plots
    │   │   │   ├── <screen-batch-2>_<condition-screen-1>_vs_<screen-batch-2>_<control-screen-1>_scatter.png
    │   │   │   ├── ...
    │   │   ├── condition_gene_calls.tsv
    │   │   ├── t0_normalized_screens_guide_level.tsv
    │   ├──  <screen-batch-3>
    │   │   ├── ... 
    │   ├── ... 
    │   ├── differential_LFC_scores.tsv
    │   ├── fdr_scores.tsv

Output directory organization and files from global-normalization:

    .
    ├── global_normalization
    │   ├── global_normalized_dLFC_scores.tsv                  
    │   ├── fdr_scores_all.tsv
    │   ├── scores_all.csv
    │   ├── wbc_scores.csv
    │   ├── sd_scale_table.tsv
    │   ├── control
    │   │   ├── control
    │   │   │   ├── control_effect_scores.tsv
    │   │   ├── control_control_map_table.tsv
    │   │   ├── control_replicates_map_table.tsv
    │   │   ├── replicate_cor.tsv
    │   ├── LDA_evaluation_plots
    │   │   ├── bc_lda_<component_number>_histogram.png
    │   │   ├── ...
    │   │   ├── bc_lda_<component_number>_roc.png
    │   │   ├── ...
    │   ├── plots
    │   │   ├── <screen-batch-1>_<condition-screen-1>_vs_<screen-batch-1>_<control-screen-1>_scatter.png
    │   │   ├── <screen-batch-1>_<condition-screen-2>_vs_<screen-batch-1>_<control-screen-2>_scatter.png
    │   │   ├── <screen-batch-2>_<condition-screen-1>_vs_<screen-batch-2>_<control-screen-1>_scatter.png
    │   │   ├── ...
