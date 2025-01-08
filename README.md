# scPANDA
DrPet is a tool designed to infer cell types in blood single-cell RNA sequencing (scRNA-seq) datasets using a three-layer annotation approach.

## Usage
To run DrPet, use the following command:

python -m drpet <file_path> <cell_type_column> [optional arguments]

## Required Arguments
file_path (str): The file path of the query dataset. cell_type_column (str): The name of the cell type column in the query dataset.

## Optional Arguments
--alpha (float, default=0.1): The 1st correlation threshold for Layer 1 inference (T cells or not). --beta (float, default=0): The 2nd correlation threshold for Layer 1 inference (T cells or not). --gamma (float, default=0.2): The correlation threshold for reannotation. -t, --t_cell_only (store_true): Consider only T cells for reannotation. -w, --write_adata (store_true): Write annotated data to output. -p, --plot (store_true): Generate plots. -pp, --preprocess (store_true): Preprocess the data before annotation.
