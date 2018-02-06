library(cellrangerRkit)
library(reshape2)
library(Biobase)
library(org.Hs.eg.db)
library(GSEABase)

tenXplot <- function (gbm, gene_probes, projection, limits = c(0, 10), marker_size = 0.1, marker_scale = c(.5, 4),
          title = NULL) 
{
  gbm_trunc <- trunc_gbm_by_genes(gbm, gene_probes)
  gene_values <- t(as.matrix(exprs(gbm_trunc)))
  gene_values[gene_values < limits[1]] <- limits[1]
  gene_values[gene_values > limits[2]] <- limits[2]
  colnames(gene_values) <- gene_probes
  projection_names <- colnames(projection)
  colnames(projection) <- c("Component.1", "Component.2")
  
  if(marker_size %in% gene_probes) {
    size = gene_values[,marker_size]
  } else {
    size = rep(marker_size, nrow(gene_values))
  }
  proj_gene <- data.frame(cbind(projection, gene_values, size))
  proj_gene_melt <- melt(proj_gene, id.vars = c("Component.1", 
                                                "Component.2", "size"))
  p <- ggplot(proj_gene_melt, aes(Component.1, Component.2)) + 
    geom_point(aes(fill = value, size=size), color = "dimgray", stroke=.2, alpha=0.5, pch=21) + 
    facet_wrap(~variable) + 
    labs(x = projection_names[1], 
         y = projection_names[2]) +
    scale_size_continuous(range = marker_scale, guide="none") +
    scale_fill_gradient(low = "darkblue", high = "firebrick1")
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  p <- p + theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(p)
}

tenXoverlay <- function (projection, size=NULL, color=NULL, marker_range=c(0.5, 5),
                      title = NULL, size_lab="size", color_lab="color") 
{
  projection_names <- colnames(projection)
  colnames(projection) <- c("Component.1", "Component.2")
  
  if(is.null(size)) {
    size = rep(0.5, nrow(projection))
  } else if(length(size) != nrow(projection)) {
    stop("Size must either be null or of length nrow(projection)");
  }
  
  if(is.null(color)) {
    color = rep(1, nrow(projection))
  } else if(length(color) != nrow(projection)) {
    stop("Color must either be null or of length nrow(projection)");
  }
  proj_data <- data.frame(cbind(projection, color, size))
  proj_data_melt <- melt(proj_data, id.vars = c("Component.1", 
                                                "Component.2", "color", "size"))
  p <- ggplot(proj_data_melt, aes(Component.1, Component.2)) + 
    geom_point(aes(fill = color, size=size), color = "dimgray", stroke=.2, alpha=0.5, pch=21) + 
    labs(x = projection_names[1], 
         y = projection_names[2])
   if(!is.discrete(size)) {
     p <- p + scale_size_continuous(range = marker_range)
   }
  if(!is.discrete(color)) {
    p <- p + scale_fill_gradient(low = "darkblue", high = "firebrick1")
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  p <- p + theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(p)
}

tenXsym <- function(gbm, sym) {
  sym <- sym[which(sym %in% fData(gbm)$symbol)]
  ix <- match(sym, fData(gbm)$symbol)
  data <- exprs(gbm)[ix,]
  rownames(data) <- sym
  as.data.frame(t(as(data, "matrix")))
}

tenXpdf <- function(...) {
  f <- tempfile(fileext = ".pdf")
  pdf(file=f, height=4, width=8)
  print(tenXplot(...))
  dev.off()
  openPDF(f)
}

# extract expression data from gbm and set rownames
# rownames will be symbol by default. If symbol=FALSE 
# then will use ensemble gene id
tenXtoDf <- function(gbm, symbol=TRUE) {
  if(symbol) {
    rn <- fData(gbm)$symbol
  } else {
    rn <- fData(gbm)$id
  }
  df <- as.data.frame(as(exprs(gbm), "matrix"))
  ix <- which(duplicated(rn))
  df <- df[-ix,]
  rn <- rn[-ix]
  rownames(df) <- rn
  df
}

# extract expression data from gbm and set rownames
# rownames will be symbol by default. If symbol=FALSE 
# then will use ensemble gene id
tenXtoMat <- function(gbm, symbol=TRUE) {
  if(symbol) {
    rn <- fData(gbm)$symbol
  } else {
    rn <- fData(gbm)$id
  }
  df <- as(exprs(gbm), "matrix")
  ix <- which(duplicated(rn))
  if(length(ix)) {
    df <- df[-ix,]
    rn <- rn[-ix]
  }
  rownames(df) <- rn
  df
}

# fetch gene symbols associated with a GO id
goToSym <- function(term) {
  as.character(mget(mget(term, org.Hs.egGO2ALLEGS)[[1]], org.Hs.egSYMBOL))
}

# fetch gene symbols associated with a GO id
goToGeneSet <- function(term) {
  gg <- as.character(mget(mget(term, org.Hs.egGO2ALLEGS)[[1]], org.Hs.egSYMBOL))
  GeneSet(gg, setName=term)
}

