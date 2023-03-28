# SCoR - An R package for Microbiome Correlation

## Installation
```
library(devtools)
devtools:install_github("AbhishakeL/SCoR")
```
## Pre-requisites
Install [FastSpar](https://github.com/scwatts/fastspar).  

## Pipeline

1. Create a phyloseq object. Possible examples can be found [here](https://cryptick-lab.github.io/NGS-Analysis/_site/R-PhyloseqIntro.html) or [here](https://joey711.github.io/phyloseq/import-data.html) .
2. Use `FSprep()` to create metadata attribute-wise partitioned OTU table.
3. Write the output of the previous command to a file.
4. Run FastSpar on the previous file.
5. ### Remove the `#` from the first cell in the correlation table and p-value table output from `FastSpar`. ###
6. Use `scor()` to generate the final output. Please follow example for details of parameter.
7. Visualise network using `igraph` or `Cytoscape'.
