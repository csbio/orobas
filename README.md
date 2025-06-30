# Orobas

Orobas is an R package for scoring chemogenomic CRISPR screening data.

## Overview

## Installation

## How to run

## Expected Outputs

Output directory and files from single-screen scoring:

    .
    ├── <output>
    │   ├── <screen-batch-1>    # (Directory)                  
    │   │   ├── ... 
    │   ├── <screen-batch-2>    # (Directory)
    │   │   ├── qc    # (Directory) screen replicate LFC scatter plots and other quality control files                                       
    │   │   │   ├── essential_PR_QC.tsv    # precision-recall AUC of essential-targeting guides from screen replicates
    │   │   │   ├── lfc_heatmap.png    # heatmap of Pearson Correlation among screen replicate LFCs
    │   │   │   ├── replicate_cor.tsv    # Pearson Correlation among screen replicate LFCs
    │   │   │   ├── <screen-batch-2>_<condition-screen-1-replicate-A>_vs_<screen-batch-2>_<condition-screen-1-replicate-B>_replicate_comparison.png
    │   │   │   ├── ...
    │   │   │   ├── <screen-batch-2>_<control-screen-1-replicate-A>_vs_<screen-batch-2>_<control-screen-1-replicate-B>_replicate_comparison.png
    │   │   │   ├── ...
    │   │   │   ├── reads    # (Directory) raw read count histograms of screen replicates and other quality control files
    │   │   │   │   ├── total_reads.png    # bar-plot of raw read counts from all screen replicates   
    │   │   │   │   ├── reads_heatmap.png    # heatmap of Pearson Correlation among screen replicate raw read counts
    │   │   │   │   ├── <screen-batch-2_T0>_raw_reads_histogram.png
    │   │   │   │   ├── ...
    │   │   │   │   ├── <screen-batch-2>_<control-screen-1-replicate-A>_raw_reads_histogram.png
    │   │   │   │   ├── ...
    │   │   │   │   ├── <screen-batch-2>_<condition-screen-1-replicate-A>_raw_reads_histogram.png
    │   │   │   │   ├── ...
    │   │   ├── guide_dlfc    # (Directory) guide-level replicate-level dLFC score file
    │   │   │   ├── <screen-batch-2>_<condition-screen-1>_vs_<screen-batch-2>_<control-screen-1>_guide_dlfc_pre_jk.tsv
    │   │   │   ├── ...
    │   │   ├── plots    # (Directory) scatter plots of gene-level condition LFCs vs control LFCs with negative and positive interactions
    │   │   │   ├── <screen-batch-2>_<condition-screen-1>_vs_<screen-batch-2>_<control-screen-1>_scatter.png
    │   │   │   ├── ...
    │   │   ├── condition_gene_calls.tsv    # ***score file containing gene-level screen-level LFC, dLFC, FDR, significant hits and other values
    │   │   ├── t0_normalized_screens_guide_level.tsv    # guide-level replicate-level LFC score file 
    │   ├──  <screen-batch-3>    # (Directory)
    │   │   ├── ... 
    │   ├── ... 
    │   ├── differential_LFC_scores.tsv    # gene-level dLFC scores from all screens from all screen-batches 
    │   ├── fdr_scores.tsv    # gene-level FDR scores from all screens from all screen-batches

Output directory organization and files from global-normalization:

    .
    ├── global_normalization    # (Directory) output files generated after running global normalization
    │   ├── global_normalized_dLFC_scores.tsv    # ***file with normalized dLFC scores from all selected condition screens                      
    │   ├── fdr_scores_all.tsv    # file with FDR scores from all selected condition screens
    │   ├── scores_all.csv    # file with LFC, normalized dLFC, FDR scores, and updated significant hits from all selected condition screens
    │   ├── wbc_scores.csv    # file with within-between correlation scores after each normalization step
    │   ├── sd_scale_table.tsv    # standard deviation of dLFC scores before and after scaling step
    │   ├── control    # (Directory) control screen files
    │   │   ├── control    # (Directory) control dLFC score file
    │   │   │   ├── control_effect_scores.tsv    # file with dLFC scores from control screens
    │   │   ├── control_control_map_table.tsv    # control_control_map table used in generating control dLFC scores
    │   │   ├── control_replicates_map_table.tsv    # control_replicates_map table used in generating control dLFC scores
    │   │   ├── replicate_cor.tsv    # Pearson correlation among control screen replicates
    │   ├── LDA_evaluation_plots    # (Directory) global ROCAUC and per-screen ROCAUC histograms at each LDA component removal step
    │   │   ├── bc_lda_<component_number>_histogram.png
    │   │   ├── ...
    │   │   ├── bc_lda_<component_number>_roc.png
    │   │   ├── ...
    │   ├── plots    # (Directory) scatter plots of gene-level condition LFCs vs control LFCs with negative and positive interactions
    │   │   ├── <screen-batch-1>_<condition-screen-1>_vs_<screen-batch-1>_<control-screen-1>_scatter.png
    │   │   ├── <screen-batch-1>_<condition-screen-2>_vs_<screen-batch-1>_<control-screen-2>_scatter.png
    │   │   ├── <screen-batch-2>_<condition-screen-1>_vs_<screen-batch-2>_<control-screen-1>_scatter.png
    │   │   ├── ...
