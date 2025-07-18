# cellxplorer: Interactive cell selection and marker exploration

 ## Example

https://crem-bu.shinyapps.io/cellxplorer/


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

