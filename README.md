# cellxplorer: Interactive cell selection and marker exploration

 ## Example

See Selector tab on:

https://crem-bu.shinyapps.io/cellxplorer/

Use interactive Lasso or square tooltip to select a group of cells.

Allows the pairwise comparison between (one or more) groups for any of categorical variable in the metadata.

 ## Requirements

The ShinyCell app directory must have a file with Highly Variable Genes: var_features.rds

It assumes there is a column named `orig.ident` in the metadata

## Output

A table with the top 30 markers for the selected cells containing:
* gene name
* average lognorm expression in selected cells
* average lognorm expression in rest of cells
* Welch t-test statistic for selected vs rest of cells
* z-scores for selected cells (vs all)
* EnrichR results for 4 gene-set collections (only shown the top 5 gene-sets for each collection) using as query the top 200 markers of the selected cells


