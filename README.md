# TenXTools
Putting the fun in 10XGenomics gene expression analysis

Right now this is just a small collection of functions to facilitate common tasks in the course of analyzing and visualizing 10XGenomicx Chromium gene expression data.  Over time this will grow into a proper library.

Here is a quick example.

```
# load up the results from cellranger count
setwd("/media/secondary/someproject")
path <- "./H7AGG"
gbm <- load_cellranger_matrix(path)
analysis_results <- load_cellranger_analysis_results(path)

source("TenXTools.R")

# look up some genes and use them to annotate a TSNE plot
tsne_proj <- analysis_results$tsne
data <- tenXsym(gbm, c("TNNT2", "ACTN1", "SIRPA", "NKX2-5", "SHOX2"))
tenXoverlay(tsne_proj[c("TSNE.1","TSNE.2")], 
            size = data$SHOX2, 
            color = grepl("-1", colnames(gbm)), # library 1 vs others
            marker_range = c(1,6))
```
