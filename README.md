# LAM_Spatial_Transcriptome

## PACKAGES REQUIRED:

1. shiny (version 1.9.1)
2. dplyr (version 1.1.4)
3. Seurat (version 5.1.0)
4. patchwork (version 1.2.0)
5. hdf5r (version 1.3.11)


## DESCRIPTION:

### Spatial Visualization Tool (lite)
The user picks an RDS file from a dropdown menu.

The main panel of this app consists of four quadrants.
- Upper-left quadrant: An image of the cell tissue from which the sample was drawn.
- Upper-right quadrant: A two-dimensional graph of the cells, reduced using UMAP. Clusters must be pre-defined in the RDS file before upload; this app does not compute clusters from scratch.
- Lower-left quadrant: Up to 4 violin plots of biomarker expression across clusters. This quadrant remains empty unless one or more biomarkers are selected in the sidebar panel.
- Lower-right quadrant: Up to 4 scatterplots of biomarker expression across cells. The higher the expression of the biomarker, the darker the cells will be colored. This quadrant remains empty unless one or more biomarkers are selected in the sidebar panel.

Inputs to this app must be QC’ed RDS files. The app does not do quality control, for the sake of space and time complexity.
If the names of pre-defined clusters are included, they will be displayed in the upper right quadrant. Otherwise, the clusters will be labeled 0, 1, 2, and so on.
Visualizations are created with the Seurat package in R.

I’m also considering adding a cellxgene visualization to the web portal. However, I don’t know how feasible that’s going to be. Those visualizations are created using input from Command Prompt (Windows) or Terminal (Mac).
While this Shiny app can handle RDS files, cellxgene requires that files be converted to H5AD first.

NOTE: The current RDS files I wish to upload to this repository are too large (177.4 MB, 177.4 MB, and 371.8 MB). GitHub Desktop only allows files under 100 MB to be uploaded.