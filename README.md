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
    │   │   ├── qc    # screen replicate LFC scatter plots and other quality control files                                       
    │   │   │   ├── essential_PR_QC.tsv    # precision-recall AUC of essential-targeting guides from screen replicates
    │   │   │   ├── lfc_heatmap.png    # heatmap of Pearson Correlation among screen replicate LFCs
    │   │   │   ├── replicate_cor.tsv    # Pearson Correlation among screen replicate LFCs
    │   │   │   ├── <screen-batch-2>_<condition-screen-1-replicate-A>_vs_<screen-batch-2>_<condition-screen-1-replicate-B>_replicate_comparison.png
    │   │   │   ├── ...
    │   │   │   ├── <screen-batch-2>_<control-screen-1-replicate-A>_vs_<screen-batch-2>_<control-screen-1-replicate-B>_replicate_comparison.png
    │   │   │   ├── ...
    │   │   │   ├── reads    # raw read count histograms of screen replicates and other quality control files
    │   │   │   │   ├── total_reads.png    # bar-plot of raw read counts from all screen replicates   
    │   │   │   │   ├── reads_heatmap.png    # heatmap of Pearson Correlation among screen replicate raw read counts
    │   │   │   │   ├── <screen-batch-2_T0>_raw_reads_histogram.png
    │   │   │   │   ├── ...
    │   │   │   │   ├── <screen-batch-2>_<control-screen-1-replicate-A>_raw_reads_histogram.png
    │   │   │   │   ├── ...
    │   │   │   │   ├── <screen-batch-2>_<condition-screen-1-replicate-A>_raw_reads_histogram.png
    │   │   │   │   ├── ...
    │   │   ├── guide_dlfc    # guide-level replicate-level dLFC score file
    │   │   │   ├── <screen-batch-2>_<condition-screen-1>_vs_<screen-batch-2>_<control-screen-1>_guide_dlfc_pre_jk.tsv
    │   │   │   ├── ...
    │   │   ├── plots    # scatter plots of gene-level condition LFCs vs control LFCs with negative and positive interactions
    │   │   │   ├── <screen-batch-2>_<condition-screen-1>_vs_<screen-batch-2>_<control-screen-1>_scatter.png
    │   │   │   ├── ...
    │   │   ├── condition_gene_calls.tsv    # score file containing gene-level screen-level LFC, dLFC, FDR, significant hits and other values
    │   │   ├── t0_normalized_screens_guide_level.tsv    # guide-level replicate-level LFC score file 
    │   ├──  <screen-batch-3>
    │   │   ├── ... 
    │   ├── ... 
    │   ├── differential_LFC_scores.tsv    # gene-level dLFC scores from all screens from all screen-batches 
    │   ├── fdr_scores.tsv    # gene-level FDR scores from all screens from all screen-batches

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
