# cellxplorer: Interactive cell selection and marker exploration

 ## Example

See Selector tab on:

https://crem-bu.shinyapps.io/cellxplorer/

Use interactive Lasso or square tooltip to select a group of cells.

Allows the pairwise comparison between (one or more) groups for any of categorical variable in the metadata.

## Requirements

It assumes there is a column named `orig.ident` in the metadata.
For marker discovery on the Selector tab, package `presto` must be installed.

## Output

A table with the top 30 markers for the selected cells containing:
* `feature` (gene name)
* `group`
* `avgExpr`
* `logFC`
* `statistic` (Wilcoxon U statistic)
* `auc`
* `pval`
* `padj`
* `pct_in`
* `pct_out`
* EnrichR results for 4 gene-set collections (only shown the top 5 gene-sets for each collection) using as query the top 200 genes ranked by `auc`

