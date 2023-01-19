## SCPA v1.3.0

#### Major changes

- The main function of SCPA has been rewritten to be much more efficient with memory usage.
This means that RAM availability is likely not going to be limiting on a typical computer
i.e. 8GB RAM, even with a huge number of comparisons.
- The speed of the pathway comparisons is also slightly improved

#### Minor changes
- Added the option to filter gene set size in `compare_seurat()`

## SCPA v1.2.1

#### Minor changes

- Add functionality to choose assay in Seurat object

## SCPA v1.2

#### Major changes

- Add functionality to deal with SingleCellExperiment objects
- `sce_extract` can now extract data from SCE objects based on colData
- `compare_sce` can now compare pathways directly within an SCE object

## SCPA v1.1

#### Major changes

- Add functionality to visualise results from SCPA with
`plot_rank()` and `plot_heatmap()`

#### Minor changes
- Improve messaging with `compare_pathways()`

## SCPA v1.0.0

#### Minor changes

- Allow usage of multiple file types for pathways e.g. gmt, csv
- Fix CRAN archive with multicross and crossmatch

## SCPA v0.0.0.9000

- Release SCPA for multivariate testing of pathways in scRNA-seq data
