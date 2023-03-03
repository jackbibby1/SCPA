## SCPA v1.5.2

### Minor changes

Ability to specify multiple gmt filepaths in the pathway argument of `compare_pathways()`. For
exmaples see [here](https://jackbibby1.github.io/SCPA/articles/using_gene_sets.html#using-a-gmt-file)

## SCPA v1.5.1

### Minor changes

Fix downsample so it isn't forcing 500 cells [#34](https://github.com/jackbibby1/SCPA/issues/34)

## SCPA v1.5.0

#### Major changes

- Parallel processing is now built into into the `compare_pathways()`, `compare_seurat()`,
and `compare_sce()` functions. Users just
need to specify `parallel = TRUE` and `cores = x` to use parallel processing instead of calling
`compare_pathways_parallel()`

#### Minor changes

- Fix downsample [#31](https://github.com/jackbibby1/SCPA/pull/31/commits/da5b7bf3a11abbf071ca5e2a9c5743a3a9f320fb)

## SCPA v1.4.0

#### Major changes

- Parallel implementation for `compare_pathways()` via `compare_pathways_parallel()`.
Users can specify number of cores to run in the pathway comparison

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
