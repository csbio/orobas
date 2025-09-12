# Example scripts for scoring chemical–genetic screens

This folder contains example scripts for scoring raw reads data from chemical–genetic screens produced and published by Olivieri et al. (2020) [1].

## Input Data
Download the input data and metadata files from the [Zenodo repository](https://doi.org/10.5281/zenodo.16935293).

## Usage
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
