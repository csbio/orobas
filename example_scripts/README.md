# Example scripts for scoring chemical–genetic screens

This directory contains example scripts for scoring raw reads data from chemical–genetic screens produced and published by Olivieri et al. (2020) [1].

## Input Data
Download the input data and metadata files from the [Zenodo repository](https://doi.org/10.5281/zenodo.16935293).

## Example run instructions
1. Create a directory called 'run_orobas'. 
2. Put the input data in a directory (name it 'data') inside the 'run_orobas' directory. 
3. Download the orobas software from the Github repo. The dowloaded package should be in a directory called 'orobas-main'. Put the 'orobas-main' directory inside the 'run_orobas' directory.
4. Download and put the example scripts (score_single_screens_wrapper_olivier.R, global_normalization_wrapper_olivier.R) in the 'run_orobas' directory. 

The figure below summarizes the directory tree structure required to run the example files if you want to run them without modifying the file or directory paths.
```bash
    run_orobas # (Directory)
    ├── score_single_screens_wrapper_olivier.R
    ├── global_normalization_wrapper_olivier.R
    ├── data # (Directory)
    │   ├── Dataset_all_readcounts.txt  
    │   ├── screen_replicate_map_table.tsv 
    │   ├── condition_control_map_table.tsv 
    ├── orobas-main # (Directory) contains all orobas software files
    │   ├── python 
    │   │   ├── ... # all python scripts 
    │   ├── R
    │   │   ├── ... # all R scripts
    │   ├──  ... # other orobas files/directories
  ```

- **Score single screens independently**  
  Run:
  ```bash
  Rscript score_single_screens_wrapper_olivier.R
  ```

- **Normalize scored dLFCs across multiple screens**  
  Run:
  ```bash
  Rscript global_normalization_wrapper_olivier.R
  ```

## Reference
[1] Olivieri, M., Cho, T., Álvarez-Quilón, A., Li, K., Schellenberg, M.J., Zimmermann, M., Hustedt, N., Rossi, S.E., Adam, S., Melo, H., et al. (2020). A Genetic Map of the Response to DNA Damage in Human Cells. Cell 182, 481-496.e21. https://doi.org/10.1016/j.cell.2020.05.040.
