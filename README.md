# CengenApp

Discovery and analysis of the C. elegans Neuronal Gene Expression Network – CeNGEN 

The CeNGEN app is based on Rshiny and uses multiple functions or modified functions from the Seurat package version 3.2.2 and others. This README provides descriptions of all functionalities presented in each tab.

1.- Gene expression by cell type.

This tab is designed to explore CeNGEN expression matrices under the different thresholds described in Taylor et al. from more stringent (“4”) to less stringent (“1”) or with no filter at all (“Unfiltered”). Gene expression values under the threshold have been converted to 0.

This tab offers three possibilities from left to right:
	
•	Cell type: Choose a cell type and retrieve all genes expressed in that cell type under the chosen threshold sorted by expression.

•	Gene: Choose a gene and retrieve all cell types expressing that gene under the chosen threshold sorted by expression.

•	Multiple genes: Retrieve the expression matrix with the chosen threshold for the subset of genes introduced in the box and separated by “,”. 

2.- Find markers based on percentage of expression.

This tab allows the user to query which genes are expressed in a given percentage of cells. The user needs to provide groups of cell types separated by “,” in either of the two boxes: group 1 or group 2. The first box retrieves genes whose percentage of cells expressing it is over the chosen threshold.  The second box retrieves genes whose percentage of cells expressing the is below the chosen threshold. The two boxes can be combined to retrieve genes over a threshold in a group of cells and under another threshold in another or the same group of cells at the same time. The result of each independent query and the combined query are displayed in columns 1-2 and 3, respectively. 

To calculate percentages of cells expressing each gene we consider a gene is expressed in a given cell if uncorrected counts in the Seurat’s RNA assay are greater than 0.

3.- Enriched genes by Cell Type

This function allows the user to browse the table of pre-calculated marker genes obtained from the function “FindMarkers” implemented in Seurat using default parameters on the SCT assay. The user can choose to browse two datasets: “All cell types2 and “Neurons only”. The first dataset includes non-neuronal cell types, affecting the background of cells against which the markers of each cell have been calculated.
 
4.- Find Differential Expression between Cell Types

In this tab the user can perform differential expression analysis between two cell types or groups of cell types, and also calculate enriched genes for a given cell type or group of cell types against the rest of cells in the datasets. The pairwise comparison can be achieved using the first column of options, selecting the two cell types to compare. The is also the option to calculate differential expression using a number of tests implemented in Seurat´s function “FindMarkers”. In the third column, the user can introduce groups of cell types separated by commas. If in the second group of cells the user writes “ALL”, then the comparison is calculated between group 1 and the rest of cells in the dataset. Be aware that calculation against ALL can take up to 5 minutes to be computed.

5.- Single Cell Plot

This tab displays multiple plots for gene expression and metadata visualization in single cells embedded in two UMAP dimensions. On the left menu the user can choose the entire CeNGEN dataset or only neurons. The second box allows to choose the plotting of metadata (i.e. Neuron Type, Tissue Type, etc.), single gene expression or two-gene co-expression.

•	Individual genes. Writing a gene name and clicking the plot button displays three panels.
•	Single cells colored by cell type.
•	Single cells colored by the scaled expression from the Seurat’s SCT data slot. 
•	Ridge plot of the density ex scaled expression sorted by cell type.

The first panel controls the zoom capabilities. Double click zooms in, double click resets zoom.

•	Colocalization. Filling the two boxes with Gene 1 and Gene 2 the user obtains two panels: 
		The top one shows the expression level of each gene in red in green and the colocalizations in shades of yellow. The degree of color blend can be customized using the blend bar on top of the panel from less stringent 0, to more stringent 1.
		Bottom panel shows cell types and control zoom. Again, double click zooms in, double click resets zoom.


6.- Heatmaps of gene expression.

This tab allows the user to introduce a list of genes separated by commas, spaces, newlines, and also mixing gene names, wormbase ID and “seqnames” to visualize a heatmap/dotplot of expression across cell types. 

By hovering and clicking with the mouse over dots the user can read the information of the values of each dot in terms of proportion of cells expressing the gene, the scaled TPM expression values, cell type and modality. 


More detailed documentation on [cengen.org](https://www.cengen.org/single-cell-rna-seq/).




# Updates and local installation

Download the app source code from [Github](https://github.com/cengenproject/CengenApp). Ensure all packages listed at the beginning of the script `ui.R` are installed.

Download the required dataset from the [cengen.org download page](https://www.cengen.org/downloads/) or by request to authors. Run the content of `validate_dataset.R` to ensure the files and format are correct.

For local use, ensure the configuration file is available and updated, e.g. `deployment/config_L4.R`. Importantly, make sure in this file `data_dir` points to the correct dataset.

If you have several config files, at the top of `global.R`, set the default file.

This should work for a local installation.

For a deployment to shinyapps.io, you additionally want to have a manifest file, e.g. `deployment/manifest_L4.txt`, that lists the files that are to be uploaded to the shinyapps.io server. Note that any file not listed in the manifest will not be uploaded. When the app is started, it looks for files named `deployment/config_xxx.R`, if it finds a single file, this describes the dataset to use.

Use the `rsconnect` calls below to upload the app, using only the manifest for one specific dataset.



## Deployment

The same App is used with several datasets and slightly different presentations. No `envVars` on Shinyapps.io, We include a file `deployment/config_xx.R` indicating we're using the app for `xx` (e.g. L4, L1, ...).


L4 app:
```
rsconnect::deployApp(appDir = ".", appFileManifest = "deployment/manifest_L4.txt", appName = "L4app", appTitle = "L4 app", upload = TRUE)
```

L1 app:
```
rsconnect::deployApp(appDir = ".", appFileManifest = "deployment/manifest_L1.txt", appName = "L1app", appTitle = "L1 app", upload = TRUE)
```

Adult app:
```
rsconnect::deployApp(appDir = ".", appFileManifest = "deployment/manifest_adult.txt", appName = "adult", appTitle = "adult app", upload = TRUE)
```








