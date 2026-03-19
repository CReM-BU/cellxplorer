library(shiny) 
library(shinyhelper) 
library(data.table) 
library(Matrix) 
library(DT) 
library(magrittr) 
library(ggplot2) 
library(ggrepel) 
library(hdf5r) 
library(ggdendro) 
library(gridExtra) 
library(plotly)
library(httr)
library(jsonlite)
library(dplyr)
library(stringr)

sc1conf = readRDS("sc1conf.rds")
sc1def  = readRDS("sc1def.rds")
sc1gene = readRDS("sc1gene.rds")
sc1meta = readRDS("sc1meta.rds")

sc1meta$cellid <- if ("cellid" %in% names(sc1meta)) {
  as.character(sc1meta$cellid)
} else {
  rownames(sc1meta)
}

gene_idx <- as.integer(sc1gene)
genes_by_idx <- names(sc1gene)[order(gene_idx)]
cellid_lookup <- sc1meta$cellid

# Variable features for quick marker mode (pre-intersected with available genes)
var_features <- readRDS("var_features.rds")
var_features <- var_features[var_features %in% genes_by_idx]

lookupCellIndices <- function(ids) {
  idx <- match(as.character(ids), cellid_lookup)
  sort(unique(idx[!is.na(idx)]))
}

# Shared cache for the sparse expression matrix.
# Populated either by a background preload (session$onFlushed) or
# synchronously on first marker request — whichever comes first.
.expr_cache <- new.env(parent = emptyenv())
.expr_cache$mat <- NULL

buildExprSparse <- function() {
  h5file <- H5File$new("sc1gexpr.h5", mode = "r")
  on.exit(h5file$close_all(), add = TRUE)

  h5data <- h5file[["grp"]][["data"]]
  dims   <- h5data$dims
  nGenes <- dims[1]
  nCells <- dims[2]

  # Read in row-chunks and accumulate sparse triplets to avoid
  # materialising the full dense matrix (which doubles peak RAM).
  # Target ~32 MB per chunk, adapting to dataset width.
  chunk_size <- max(1L, as.integer(floor(32e6 / (nCells * 8L))))
  n_chunks   <- ceiling(nGenes / chunk_size)
  ti <- vector("list", n_chunks)
  tj <- vector("list", n_chunks)
  tx <- vector("list", n_chunks)

  for (k in seq_len(n_chunks)) {
    r1 <- (k - 1L) * chunk_size + 1L
    r2 <- min(k * chunk_size, nGenes)
    chunk <- h5data$read(args = list(r1:r2, quote(expr = )))
    chunk[is.na(chunk)] <- 0
    nz <- which(chunk != 0, arr.ind = TRUE)
    if (nrow(nz) > 0L) {
      ti[[k]] <- nz[, 1L] + (r1 - 1L)
      tj[[k]] <- nz[, 2L]
      tx[[k]] <- chunk[nz]
    }
  }

  Matrix::sparseMatrix(
    i    = unlist(ti),
    j    = unlist(tj),
    x    = as.double(unlist(tx)),
    dims = c(nGenes, nCells),
    dimnames = list(genes_by_idx, NULL)
  )
}

selectorExprSparse <- function() {
  if (!is.null(.expr_cache$mat)) return(.expr_cache$mat)
  .expr_cache$mat <- buildExprSparse()
  .expr_cache$mat
}

# Keep Enrichr network I/O off the main Shiny worker so marker tables render first.
if (inherits(future::plan(), "sequential")) {
  future::plan(future::multisession, workers = 2)
}
# Pre-warm one worker so subprocess startup doesn't compete with first renders.
tryCatch(future::future({ TRUE }), error = function(e) NULL)



### Useful stuff 
# Colour palette 
cList = list(c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84", 
               "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000"), 
             c("#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF", 
               "#FEE090","#FDAE61","#F46D43","#D73027")[c(1,1:9,9)], 
             c("#FDE725","#AADC32","#5DC863","#27AD81","#21908C", 
               "#2C728E","#3B528B","#472D7B","#440154")) 
names(cList) = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple") 

# Panel sizes 
pList = c("400px", "600px", "800px") 
names(pList) = c("Small", "Medium", "Large") 
pList2 = c("500px", "700px", "900px") 
names(pList2) = c("Small", "Medium", "Large") 
pList3 = c("600px", "800px", "1000px") 
names(pList3) = c("Small", "Medium", "Large") 
sList = c(18,24,30) 
names(sList) = c("Small", "Medium", "Large") 
lList = c(5,6,7) 
names(lList) = c("Small", "Medium", "Large") 

# Function to extract legend 
g_legend <- function(a.gplot){  
  tmp <- ggplot_gtable(ggplot_build(a.gplot))  
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")  
  legend <- tmp$grobs[[leg]]  
  legend 
}  

# Plot theme 
sctheme <- function(base_size = 24, XYval = TRUE, Xang = 0, XjusH = 0.5){ 
  oupTheme = theme( 
    text =             element_text(size = base_size, family = "Helvetica"), 
    panel.background = element_rect(fill = "white", colour = NA), 
    axis.line =   element_line(colour = "black"), 
    axis.ticks =  element_line(colour = "black", size = base_size / 20), 
    axis.title =  element_text(face = "bold"), 
    axis.text =   element_text(size = base_size), 
    axis.text.x = element_text(angle = Xang, hjust = XjusH), 
    legend.position = "bottom", 
    legend.key =      element_rect(colour = NA, fill = NA) 
  ) 
  if(!XYval){ 
    oupTheme = oupTheme + theme( 
      axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
      axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
  } 
  return(oupTheme) 
} 

### Common plotting functions 
# Plot cell information on dimred 
scDRcell <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt, inplab){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                       inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "Y", "val", "sub") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  if(inpord == "Max-1st"){ 
    ggData = ggData[order(val)] 
  } else if(inpord == "Min-1st"){ 
    ggData = ggData[order(-val)] 
  } else if(inpord == "Random"){ 
    ggData = ggData[sample(nrow(ggData))] 
  } 
  
  # Do factoring if required 
  if(!is.na(inpConf[UI == inp1]$fCL)){ 
    ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]] 
    names(ggCol) = levels(ggData$val) 
    ggLvl = levels(ggData$val)[levels(ggData$val) %in% unique(ggData$val)] 
    ggData$val = factor(ggData$val, levels = ggLvl) 
    ggCol = ggCol[ggLvl] 
  } 
  
  # Actual ggplot 
  ggOut = ggplot(ggData, aes(X, Y, color = val)) 
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + 
    geom_point(size = inpsiz, shape = 16) + xlab(inpdrX) + ylab(inpdrY) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt) 
  if(is.na(inpConf[UI == inp1]$fCL)){ 
    ggOut = ggOut + scale_color_gradientn("", colours = cList[[inpcol]]) + 
      guides(color = guide_colorbar(barwidth = 15)) 
  } else { 
    sListX = min(nchar(paste0(levels(ggData$val), collapse = "")), 200) 
    sListX = 0.75 * (sList - (1.5 * floor(sListX/50))) 
    ggOut = ggOut + scale_color_manual("", values = ggCol) + 
      guides(color = guide_legend(override.aes = list(size = 5),  
                                  nrow = inpConf[UI == inp1]$fRow)) + 
      theme(legend.text = element_text(size = sListX[inpfsz])) 
    if(inplab){ 
      ggData3 = ggData[, .(X = mean(X), Y = mean(Y)), by = "val"] 
      lListX = min(nchar(paste0(ggData3$val, collapse = "")), 200) 
      lListX = lList - (0.25 * floor(lListX/50)) 
      ggOut = ggOut + 
        geom_text_repel(data = ggData3, aes(X, Y, label = val), 
                        color = "grey10", bg.color = "grey95", bg.r = 0.15, 
                        size = lListX[inpfsz], seed = 42) 
    } 
  } 
  if(inpasp == "Square") { 
    ggOut = ggOut + coord_fixed(ratio = rat) 
  } else if(inpasp == "Fixed") { 
    ggOut = ggOut + coord_fixed() 
  } 
  return(ggOut) 
} 

scDRnum <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, 
                    inpH5, inpGene, inpsplt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), 
                   with = FALSE] 
  colnames(ggData) = c("group", "sub") 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Split inp1 if necessary 
  if(is.na(inpConf[UI == inp1]$fCL)){ 
    if(inpsplt == "Quartile"){nBk = 4} 
    if(inpsplt == "Decile"){nBk = 10} 
    ggData$group = cut(ggData$group, breaks = nBk) 
  } 
  
  # Actual data.table 
  ggData$express = FALSE 
  ggData[val2 > 0]$express = TRUE 
  ggData1 = ggData[express == TRUE, .(nExpress = .N), by = "group"] 
  ggData = ggData[, .(nCells = .N), by = "group"] 
  ggData = ggData1[ggData, on = "group"] 
  ggData = ggData[, c("group", "nCells", "nExpress"), with = FALSE] 
  ggData[is.na(nExpress)]$nExpress = 0 
  ggData$pctExpress = 100 * ggData$nExpress / ggData$nCells 
  ggData = ggData[order(group)] 
  colnames(ggData)[3] = paste0(colnames(ggData)[3], "_", inp2) 
  return(ggData) 
} 
# Plot gene expression on dimred 
scDRgene <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2, 
                     inpH5, inpGene, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "Y", "sub") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
  ggData[val < 0]$val = 0 
  h5file$close_all() 
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  if(inpord == "Max-1st"){ 
    ggData = ggData[order(val)] 
  } else if(inpord == "Min-1st"){ 
    ggData = ggData[order(-val)] 
  } else if(inpord == "Random"){ 
    ggData = ggData[sample(nrow(ggData))] 
  } 
  
  # Actual ggplot 
  ggOut = ggplot(ggData, aes(X, Y, color = val)) 
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + 
    geom_point(size = inpsiz, shape = 16) + xlab(inpdrX) + ylab(inpdrY) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt) +  
    scale_color_gradientn(inp1, colours = cList[[inpcol]]) + 
    guides(color = guide_colorbar(barwidth = 15)) 
  if(inpasp == "Square") { 
    ggOut = ggOut + coord_fixed(ratio = rat) 
  } else if(inpasp == "Fixed") { 
    ggOut = ggOut + coord_fixed() 
  } 
  return(ggOut) 
} 

# Plot gene coexpression on dimred 
bilinear <- function(x,y,xy,Q11,Q21,Q12,Q22){ 
  oup = (xy-x)*(xy-y)*Q11 + x*(xy-y)*Q21 + (xy-x)*y*Q12 + x*y*Q22 
  oup = oup / (xy*xy) 
  return(oup) 
} 
scDRcoex <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inp2, 
                     inpsub1, inpsub2, inpH5, inpGene, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "Y", "sub") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val1 = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
  ggData[val1 < 0]$val1 = 0 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Generate coex color palette 
  cInp = strsplit(inpcol, "; ")[[1]] 
  if(cInp[1] == "Red (Gene1)"){ 
    c10 = c(255,0,0) 
  } else if(cInp[1] == "Orange (Gene1)"){ 
    c10 = c(255,140,0) 
  } else { 
    c10 = c(0,255,0) 
  } 
  if(cInp[2] == "Green (Gene2)"){ 
    c01 = c(0,255,0) 
  } else { 
    c01 = c(0,0,255) 
  } 
  c00 = c(217,217,217) ; c11 = c10 + c01 
  nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2 
  gg = data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1))) 
  gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid 
  gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid 
  gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1]) 
  gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2]) 
  gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3]) 
  gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255) 
  gg = gg[, c("v1", "v2", "cMix")] 
  
  # Map colours 
  ggData$v1 = round(nTot * ggData$val1 / max(ggData$val1)) 
  ggData$v2 = round(nTot * ggData$val2 / max(ggData$val2)) 
  ggData$v0 = ggData$v1 + ggData$v2 
  ggData = gg[ggData, on = c("v1", "v2")] 
  if(inpord == "Max-1st"){ 
    ggData = ggData[order(v0)] 
  } else if(inpord == "Min-1st"){ 
    ggData = ggData[order(-v0)] 
  } else if(inpord == "Random"){ 
    ggData = ggData[sample(nrow(ggData))] 
  } 
  
  # Actual ggplot 
  ggOut = ggplot(ggData, aes(X, Y)) 
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + 
    geom_point(size = inpsiz, shape = 16, color = ggData$cMix) + 
    xlab(inpdrX) + ylab(inpdrY) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt) + 
    scale_color_gradientn(inp1, colours = cList[[1]]) + 
    guides(color = guide_colorbar(barwidth = 15)) 
  if(inpasp == "Square") { 
    ggOut = ggOut + coord_fixed(ratio = rat) 
  } else if(inpasp == "Fixed") { 
    ggOut = ggOut + coord_fixed() 
  } 
  return(ggOut) 
} 

scDRcoexLeg <- function(inp1, inp2, inpcol, inpfsz){ 
  # Generate coex color palette 
  cInp = strsplit(inpcol, "; ")[[1]] 
  if(cInp[1] == "Red (Gene1)"){ 
    c10 = c(255,0,0) 
  } else if(cInp[1] == "Orange (Gene1)"){ 
    c10 = c(255,140,0) 
  } else { 
    c10 = c(0,255,0) 
  } 
  if(cInp[2] == "Green (Gene2)"){ 
    c01 = c(0,255,0) 
  } else { 
    c01 = c(0,0,255) 
  } 
  c00 = c(217,217,217) ; c11 = c10 + c01 
  nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2 
  gg = data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1))) 
  gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid 
  gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid 
  gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1]) 
  gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2]) 
  gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3]) 
  gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255) 
  gg = gg[, c("v1", "v2", "cMix")] 
  
  # Actual ggplot 
  ggOut = ggplot(gg, aes(v1, v2)) + 
    geom_tile(fill = gg$cMix) + 
    xlab(inp1) + ylab(inp2) + coord_fixed(ratio = 1) + 
    scale_x_continuous(breaks = c(0, nTot), label = c("low", "high")) + 
    scale_y_continuous(breaks = c(0, nTot), label = c("low", "high")) + 
    sctheme(base_size = sList[inpfsz], XYval = TRUE) 
  return(ggOut) 
} 

scDRcoexNum <- function(inpConf, inpMeta, inp1, inp2, 
                        inpsub1, inpsub2, inpH5, inpGene){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpsub1]$ID), with = FALSE] 
  colnames(ggData) = c("sub") 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val1 = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
  ggData[val1 < 0]$val1 = 0 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Actual data.table 
  ggData$express = "none" 
  ggData[val1 > 0]$express = inp1 
  ggData[val2 > 0]$express = inp2 
  ggData[val1 > 0 & val2 > 0]$express = "both" 
  ggData$express = factor(ggData$express, levels = unique(c("both", inp1, inp2, "none"))) 
  ggData = ggData[, .(nCells = .N), by = "express"] 
  ggData$percent = 100 * ggData$nCells / sum(ggData$nCells) 
  ggData = ggData[order(express)] 
  colnames(ggData)[1] = "expression > 0" 
  return(ggData) 
} 

# Plot violin / boxplot 
scVioBox <- function(inpConf, inpMeta, inp1, inp1b, inp2,
                     inpsub1, inpsub2, inpH5, inpGene,
                     inptyp, inppts, inpsiz, inpfsz){
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]}

  # Check if secondary grouping is enabled
  useFill = !is.null(inp1b) && inp1b != "(none)" && inp1b %in% inpConf$UI

  # Prepare ggData
  if(useFill){
    ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inp1b]$ID, inpConf[UI == inpsub1]$ID),
                     with = FALSE]
    colnames(ggData) = c("X", "X2", "sub")
  } else {
    ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID),
                     with = FALSE]
    colnames(ggData) = c("X", "sub")
  }

  # Load in either cell meta or gene expr
  if(inp2 %in% inpConf$UI){
    ggData$val = inpMeta[[inpConf[UI == inp2]$ID]]
  } else {
    h5file <- H5File$new(inpH5, mode = "r")
    h5data <- h5file[["grp"]][["data"]]
    ggData$val = h5data$read(args = list(inpGene[inp2], quote(expr=)))
    ggData[val < 0]$val = 0
    set.seed(42)
    tmpNoise = rnorm(length(ggData$val)) * diff(range(ggData$val)) / 1000
    ggData$val = ggData$val + tmpNoise
    h5file$close_all()
  }
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){
    ggData = ggData[sub %in% inpsub2]
  }

  # Do factoring for X-axis
  ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]]
  names(ggCol) = levels(ggData$X)
  ggLvl = levels(ggData$X)[levels(ggData$X) %in% unique(ggData$X)]
  ggData$X = factor(ggData$X, levels = ggLvl)
  ggCol = ggCol[ggLvl]

  # Do factoring for fill variable if secondary grouping is enabled
  if(useFill){
    ggCol2 = strsplit(inpConf[UI == inp1b]$fCL, "\\|")[[1]]
    names(ggCol2) = levels(ggData$X2)
    ggLvl2 = levels(ggData$X2)[levels(ggData$X2) %in% unique(ggData$X2)]
    ggData$X2 = factor(ggData$X2, levels = ggLvl2)
    ggCol2 = ggCol2[ggLvl2]
  }

  # Actual ggplot
  if(useFill){
    if(inptyp == "violin"){
      ggOut = ggplot(ggData, aes(X, val, fill = X2)) + geom_violin(scale = "width")
    } else {
      ggOut = ggplot(ggData, aes(X, val, fill = X2)) + geom_boxplot()
    }
  } else {
    if(inptyp == "violin"){
      ggOut = ggplot(ggData, aes(X, val, fill = X)) + geom_violin(scale = "width")
    } else {
      ggOut = ggplot(ggData, aes(X, val, fill = X)) + geom_boxplot()
    }
  }
  if(inppts){
    ggOut = ggOut + geom_jitter(size = inpsiz, shape = 16)
  }
  ggOut = ggOut + xlab(inp1) + ylab(inp2) +
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1)

  if(useFill){
    ggOut = ggOut + scale_fill_manual(inp1b, values = ggCol2)
  } else {
    ggOut = ggOut + scale_fill_manual("", values = ggCol) +
      theme(legend.position = "none")
  }
  return(ggOut)
} 

# Plot proportion plot 
scProp <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, 
                   inptyp, inpflp, inpfsz){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inp2]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "grp", "sub") 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  ggData = ggData[, .(nCells = .N), by = c("X", "grp")] 
  ggData = ggData[, {tot = sum(nCells) 
  .SD[,.(pctCells = 100 * sum(nCells) / tot, 
         nCells = nCells), by = "grp"]}, by = "X"] 
  
  # Do factoring 
  ggCol = strsplit(inpConf[UI == inp2]$fCL, "\\|")[[1]] 
  names(ggCol) = levels(ggData$grp) 
  ggLvl = levels(ggData$grp)[levels(ggData$grp) %in% unique(ggData$grp)] 
  ggData$grp = factor(ggData$grp, levels = ggLvl) 
  ggCol = ggCol[ggLvl] 
  
  # Actual ggplot 
  if(inptyp == "Proportion"){ 
    ggOut = ggplot(ggData, aes(X, pctCells, fill = grp)) + 
      geom_col() + ylab("Cell Proportion (%)") 
  } else { 
    ggOut = ggplot(ggData, aes(X, nCells, fill = grp)) + 
      geom_col() + ylab("Number of Cells") 
  } 
  if(inpflp){ 
    ggOut = ggOut + coord_flip() 
  } 
  ggOut = ggOut + xlab(inp1) + 
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +  
    scale_fill_manual("", values = ggCol) + 
    theme(legend.position = "right") 
  return(ggOut) 
} 

# Get gene list 
scGeneList <- function(inp, inpGene){ 
  geneList = data.table(gene = unique(trimws(strsplit(inp, ",|;|
")[[1]])), 
                        present = TRUE) 
  geneList[!gene %in% names(inpGene)]$present = FALSE 
  return(geneList) 
} 

# Plot gene expression bubbleplot / heatmap 
scBubbHeat <- function(inpConf, inpMeta, inp, inpGrp, inpPlt, 
                       inpsub1, inpsub2, inpH5, inpGene, inpScl, inpRow, inpCol, 
                       inpcols, inpfsz, save = FALSE){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Identify genes that are in our dataset 
  geneList = scGeneList(inp, inpGene) 
  geneList = geneList[present == TRUE] 
  shiny::validate(need(nrow(geneList) <= 50, "More than 50 genes to plot! Please reduce the gene list!")) 
  shiny::validate(need(nrow(geneList) > 1, "Please input at least 2 genes to plot!")) 
  
  # Prepare ggData 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData = data.table() 
  for(iGene in geneList$gene){ 
    tmp = inpMeta[, c("sampleID", inpConf[UI == inpsub1]$ID), with = FALSE] 
    colnames(tmp) = c("sampleID", "sub") 
    tmp$grpBy = inpMeta[[inpConf[UI == inpGrp]$ID]] 
    tmp$geneName = iGene 
    tmp$val = h5data$read(args = list(inpGene[iGene], quote(expr=))) 
    ggData = rbindlist(list(ggData, tmp)) 
  } 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  shiny::validate(need(uniqueN(ggData$grpBy) > 1, "Only 1 group present, unable to plot!")) 
  
  # Aggregate 
  ggData$val = expm1(ggData$val) 
  ggData = ggData[, .(val = mean(val), prop = sum(val>0) / length(sampleID)), 
                  by = c("geneName", "grpBy")] 
  ggData$val = log1p(ggData$val) 
  
  # Scale if required 
  colRange = range(ggData$val) 
  if(inpScl){ 
    ggData[, val:= scale(val), keyby = "geneName"] 
    colRange = c(-max(abs(range(ggData$val))), max(abs(range(ggData$val)))) 
  } 
  
  # hclust row/col if necessary 
  ggMat = dcast.data.table(ggData, geneName~grpBy, value.var = "val") 
  tmp = ggMat$geneName 
  ggMat = as.matrix(ggMat[, -1]) 
  rownames(ggMat) = tmp 
  if(inpRow){ 
    hcRow = dendro_data(as.dendrogram(hclust(dist(ggMat)))) 
    ggRow = ggplot() + coord_flip() + 
      geom_segment(data = hcRow$segments, aes(x=x,y=y,xend=xend,yend=yend)) + 
      scale_y_continuous(breaks = rep(0, uniqueN(ggData$grpBy)), 
                         labels = unique(ggData$grpBy), expand = c(0, 0)) + 
      scale_x_continuous(breaks = seq_along(hcRow$labels$label), 
                         labels = hcRow$labels$label, expand = c(0, 0.5)) + 
      sctheme(base_size = sList[inpfsz]) + 
      theme(axis.title = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.y = element_blank(), 
            axis.text.x = element_text(color="white", angle = 45, hjust = 1)) 
    ggData$geneName = factor(ggData$geneName, levels = hcRow$labels$label) 
  } else { 
    ggData$geneName = factor(ggData$geneName, levels = rev(geneList$gene)) 
  } 
  if(inpCol){ 
    hcCol = dendro_data(as.dendrogram(hclust(dist(t(ggMat))))) 
    ggCol = ggplot() + 
      geom_segment(data = hcCol$segments, aes(x=x,y=y,xend=xend,yend=yend)) + 
      scale_x_continuous(breaks = seq_along(hcCol$labels$label), 
                         labels = hcCol$labels$label, expand = c(0.05, 0)) + 
      scale_y_continuous(breaks = rep(0, uniqueN(ggData$geneName)), 
                         labels = unique(ggData$geneName), expand=c(0,0)) + 
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) + 
      theme(axis.title = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.x = element_blank(), 
            axis.text.y = element_text(color = "white")) 
    ggData$grpBy = factor(ggData$grpBy, levels = hcCol$labels$label) 
  } 
  
  # Actual plot according to plottype 
  if(inpPlt == "Bubbleplot"){ 
    # Bubbleplot 
    ggOut = ggplot(ggData, aes(grpBy, geneName, color = val, size = prop)) + 
      geom_point() +  
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +  
      scale_x_discrete(expand = c(0.05, 0)) +  
      scale_y_discrete(expand = c(0, 0.5)) + 
      scale_size_continuous("proportion", range = c(0, 8), 
                            limits = c(0, 1), breaks = c(0.00,0.25,0.50,0.75,1.00)) + 
      scale_color_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) + 
      guides(color = guide_colorbar(barwidth = 15)) + 
      theme(axis.title = element_blank(), legend.box = "vertical") 
  } else { 
    # Heatmap 
    ggOut = ggplot(ggData, aes(grpBy, geneName, fill = val)) + 
      geom_tile() +  
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) + 
      scale_x_discrete(expand = c(0.05, 0)) +  
      scale_y_discrete(expand = c(0, 0.5)) + 
      scale_fill_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) + 
      guides(fill = guide_colorbar(barwidth = 15)) + 
      theme(axis.title = element_blank()) 
  } 
  
  # Final tidy 
  ggLeg = g_legend(ggOut) 
  ggOut = ggOut + theme(legend.position = "none") 
  if(!save){ 
    if(inpRow & inpCol){ggOut =  
      grid.arrange(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),  
                   layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))  
    } else if(inpRow){ggOut =  
      grid.arrange(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),  
                   layout_matrix = rbind(c(1,3),c(2,NA)))  
    } else if(inpCol){ggOut =  
      grid.arrange(ggOut, ggLeg, ggCol, heights = c(1,7,2),  
                   layout_matrix = rbind(c(3),c(1),c(2)))  
    } else {ggOut =  
      grid.arrange(ggOut, ggLeg, heights = c(7,2),  
                   layout_matrix = rbind(c(1),c(2)))  
    }  
  } else { 
    if(inpRow & inpCol){ggOut =  
      arrangeGrob(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),  
                  layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))  
    } else if(inpRow){ggOut =  
      arrangeGrob(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),  
                  layout_matrix = rbind(c(1,3),c(2,NA)))  
    } else if(inpCol){ggOut =  
      arrangeGrob(ggOut, ggLeg, ggCol, heights = c(1,7,2),  
                  layout_matrix = rbind(c(3),c(1),c(2)))  
    } else {ggOut =  
      arrangeGrob(ggOut, ggLeg, heights = c(7,2),  
                  layout_matrix = rbind(c(1),c(2)))  
    }  
  } 
  return(ggOut) 
} 




### Start server code
shinyServer(function(input, output, session) {
  ### For all tags and Server-side selectize
  observe_helpers()

  # Preload status: drives the live status line on the Selector tab.
  # States: "idle" → "loading" → "ready" (or "failed").
  preload_status_rv <- reactiveVal("idle")

  # Phase 1: Show the status banner and disable buttons the instant the user
  # opens the Selector tab — before any heavy work starts.
  observeEvent(input$main_nav, {
    if (identical(input$main_nav, "selector") &&
        preload_status_rv() == "idle" &&
        is.null(.expr_cache$mat)) {
      preload_status_rv("loading")
    }
  })

  # Phase 2: Start the actual H5 preload only after the Plotly UMAP has
  # rendered (signalled by the JS onRender callback).  This guarantees the
  # UMAP appears fast — the subprocess launch doesn't compete for CPU.
  # The worker is already warm (pre-warmed at startup), so the future
  # dispatches almost instantly.
  observeEvent(input$sel_umap_rendered, {
    if (is.null(.expr_cache$mat)) {
      if (preload_status_rv() != "loading") preload_status_rv("loading")
      promises::future_promise(
        { buildExprSparse() },
        packages = c("hdf5r", "Matrix")
      ) %>% promises::then(
        onFulfilled = function(mat) {
          if (is.null(.expr_cache$mat)) .expr_cache$mat <- mat
          preload_status_rv("ready")
        },
        onRejected = function(err) {
          warning("Background sparse-matrix preload failed: ",
                  conditionMessage(err))
          preload_status_rv("failed")
        }
      )
    }
  }, once = TRUE)

  output$preload_status <- renderUI({
    st <- preload_status_rv()
    if (st == "loading") {
      div(class = "sel-status-banner is-loading",
        div(class = "status-header",
          span(class = "status-icon", icon("cog", class = "fa-spin")),
          "Preparing marker engine"
        ),
        div(class = "status-detail",
          "Loading expression data in the background.",
          br(),
          "Draw a lasso selection on the UMAP while this completes."
        ),
        div(class = "sel-progress-track",
          div(class = "sel-progress-bar")
        )
      )
    } else if (st == "ready") {
      tagList(
        div(class = "sel-status-banner is-ready", id = "sel-status-ready",
          div(class = "status-header",
            icon("check-circle"), "Marker engine ready"
          ),
          div(class = "status-detail",
            "Select cells and click \u2018Find markers\u2019 to discover differentially expressed genes."
          )
        ),
        tags$script(HTML("
          setTimeout(function() {
            var el = document.getElementById('sel-status-ready');
            if (el) {
              el.classList.add('fade-out');
              setTimeout(function() { if (el) el.style.display = 'none'; }, 700);
            }
          }, 4000);
        "))
      )
    } else if (st == "failed") {
      div(class = "sel-status-banner is-failed",
        div(class = "status-header",
          icon("exclamation-triangle"), "Preload failed"
        ),
        div(class = "status-detail",
          "The expression matrix will load on demand when you click \u2018Find markers\u2019.",
          br(),
          "This may take a few extra seconds."
        )
      )
    } else {
      NULL
    }
  })

  # Disable "Find markers" buttons while engine is warming up,
  # re-enable once ready (or on failure, so synchronous fallback works).
  observe({
    st <- preload_status_rv()
    enabled <- st %in% c("idle", "ready", "failed")
    session$sendCustomMessage("sel_engine_state", list(enabled = enabled))
  })

  optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }"
  updateSelectizeInput(session, "sc1a1inp2", choices = names(sc1gene), server = TRUE, 
                       selected = sc1def$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "sc1a3inp1", choices = names(sc1gene), server = TRUE, 
                       selected = sc1def$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "sc1a3inp2", choices = names(sc1gene), server = TRUE, 
                       selected = sc1def$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "sc1b2inp1", choices = names(sc1gene), server = TRUE, 
                       selected = sc1def$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "sc1b2inp2", choices = names(sc1gene), server = TRUE, 
                       selected = sc1def$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "sc1c1inp2", server = TRUE, 
                       choices = c(sc1conf[is.na(fID)]$UI,names(sc1gene)), 
                       selected = sc1conf[is.na(fID)]$UI[1], options = list( 
                         maxOptions = length(sc1conf[is.na(fID)]$UI) + 3, 
                         create = TRUE, persist = TRUE, render = I(optCrt))) 
  
  
  ### Plots for tab a1 
  output$sc1a1sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1a1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1a1sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1a1sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1a1oup1 <- renderPlot({ 
    scDRcell(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp1,  
             input$sc1a1sub1, input$sc1a1sub2, 
             input$sc1a1siz, input$sc1a1col1, input$sc1a1ord1, 
             input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1) 
  }) 
  output$sc1a1oup1.ui <- renderUI({ 
    plotOutput("sc1a1oup1", height = pList[input$sc1a1psz]) 
  }) 
  output$sc1a1oup1.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                   input$sc1a1inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a1oup1.h, width = input$sc1a1oup1.w, useDingbats = FALSE, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp1,   
                      input$sc1a1sub1, input$sc1a1sub2, 
                      input$sc1a1siz, input$sc1a1col1, input$sc1a1ord1,  
                      input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1) ) 
    }) 
  output$sc1a1oup1.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                   input$sc1a1inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a1oup1.h, width = input$sc1a1oup1.w, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp1,   
                      input$sc1a1sub1, input$sc1a1sub2, 
                      input$sc1a1siz, input$sc1a1col1, input$sc1a1ord1,  
                      input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1) ) 
    }) 
  output$sc1a1.dt <- renderDataTable({ 
    ggData = scDRnum(sc1conf, sc1meta, input$sc1a1inp1, input$sc1a1inp2, 
                     input$sc1a1sub1, input$sc1a1sub2, 
                     "sc1gexpr.h5", sc1gene, input$sc1a1splt) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("pctExpress"), digits = 2) 
  }) 
  
  output$sc1a1oup2 <- renderPlot({ 
    scDRgene(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp2,  
             input$sc1a1sub1, input$sc1a1sub2, 
             "sc1gexpr.h5", sc1gene, 
             input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2, 
             input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt) 
  }) 
  output$sc1a1oup2.ui <- renderUI({ 
    plotOutput("sc1a1oup2", height = pList[input$sc1a1psz]) 
  }) 
  output$sc1a1oup2.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                   input$sc1a1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a1oup2.h, width = input$sc1a1oup2.w, useDingbats = FALSE, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp2,  
                      input$sc1a1sub1, input$sc1a1sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2, 
                      input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt) ) 
    }) 
  output$sc1a1oup2.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                   input$sc1a1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a1oup2.h, width = input$sc1a1oup2.w, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp2,  
                      input$sc1a1sub1, input$sc1a1sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2, 
                      input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt) ) 
    }) 
  
  
  ### Plots for tab a2 
  output$sc1a2sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1a2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1a2sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1a2sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1a2oup1 <- renderPlot({ 
    scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp1,  
             input$sc1a2sub1, input$sc1a2sub2, 
             input$sc1a2siz, input$sc1a2col1, input$sc1a2ord1, 
             input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab1) 
  }) 
  output$sc1a2oup1.ui <- renderUI({ 
    plotOutput("sc1a2oup1", height = pList[input$sc1a2psz]) 
  }) 
  output$sc1a2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                   input$sc1a2inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a2oup1.h, width = input$sc1a2oup1.w, useDingbats = FALSE, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp1,   
                      input$sc1a2sub1, input$sc1a2sub2, 
                      input$sc1a2siz, input$sc1a2col1, input$sc1a2ord1,  
                      input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab1) ) 
    }) 
  output$sc1a2oup1.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                   input$sc1a2inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a2oup1.h, width = input$sc1a2oup1.w, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp1,   
                      input$sc1a2sub1, input$sc1a2sub2, 
                      input$sc1a2siz, input$sc1a2col1, input$sc1a2ord1,  
                      input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab1) ) 
    }) 
  
  output$sc1a2oup2 <- renderPlot({ 
    scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp2,  
             input$sc1a2sub1, input$sc1a2sub2, 
             input$sc1a2siz, input$sc1a2col2, input$sc1a2ord2, 
             input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab2) 
  }) 
  output$sc1a2oup2.ui <- renderUI({ 
    plotOutput("sc1a2oup2", height = pList[input$sc1a2psz]) 
  }) 
  output$sc1a2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                   input$sc1a2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a2oup2.h, width = input$sc1a2oup2.w, useDingbats = FALSE, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp2,   
                      input$sc1a2sub1, input$sc1a2sub2, 
                      input$sc1a2siz, input$sc1a2col2, input$sc1a2ord2,  
                      input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab2) ) 
    }) 
  output$sc1a2oup2.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                   input$sc1a2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a2oup2.h, width = input$sc1a2oup2.w, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp2,   
                      input$sc1a2sub1, input$sc1a2sub2, 
                      input$sc1a2siz, input$sc1a2col2, input$sc1a2ord2,  
                      input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab2) ) 
    }) 
  
  
  ### Plots for tab a3 
  output$sc1a3sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1a3sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1a3sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1a3sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1a3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a3sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1a3sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1a3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a3sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1a3oup1 <- renderPlot({ 
    scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp1,  
             input$sc1a3sub1, input$sc1a3sub2, 
             "sc1gexpr.h5", sc1gene, 
             input$sc1a3siz, input$sc1a3col1, input$sc1a3ord1, 
             input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) 
  }) 
  output$sc1a3oup1.ui <- renderUI({ 
    plotOutput("sc1a3oup1", height = pList[input$sc1a3psz]) 
  }) 
  output$sc1a3oup1.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a3drX,"_",input$sc1a3drY,"_",  
                                   input$sc1a3inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a3oup1.h, width = input$sc1a3oup1.w, useDingbats = FALSE, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp1,  
                      input$sc1a3sub1, input$sc1a3sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a3siz, input$sc1a3col1, input$sc1a3ord1, 
                      input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) ) 
    }) 
  output$sc1a3oup1.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a3drX,"_",input$sc1a3drY,"_",  
                                   input$sc1a3inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a3oup1.h, width = input$sc1a3oup1.w, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp1,  
                      input$sc1a3sub1, input$sc1a3sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a3siz, input$sc1a3col1, input$sc1a3ord1, 
                      input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) ) 
    }) 
  
  output$sc1a3oup2 <- renderPlot({ 
    scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp2,  
             input$sc1a3sub1, input$sc1a3sub2, 
             "sc1gexpr.h5", sc1gene, 
             input$sc1a3siz, input$sc1a3col2, input$sc1a3ord2, 
             input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) 
  }) 
  output$sc1a3oup2.ui <- renderUI({ 
    plotOutput("sc1a3oup2", height = pList[input$sc1a3psz]) 
  }) 
  output$sc1a3oup2.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a3drX,"_",input$sc1a3drY,"_",  
                                   input$sc1a3inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a3oup2.h, width = input$sc1a3oup2.w, useDingbats = FALSE, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp2,  
                      input$sc1a3sub1, input$sc1a3sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a3siz, input$sc1a3col2, input$sc1a3ord2, 
                      input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) ) 
    }) 
  output$sc1a3oup2.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a3drX,"_",input$sc1a3drY,"_",  
                                   input$sc1a3inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a3oup2.h, width = input$sc1a3oup2.w, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp2,  
                      input$sc1a3sub1, input$sc1a3sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a3siz, input$sc1a3col2, input$sc1a3ord2, 
                      input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) ) 
    }) 
  
  
  ### Plots for tab b2 
  output$sc1b2sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1b2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1b2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1b2sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1b2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1b2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1b2sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1b2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1b2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1b2oup1 <- renderPlot({ 
    scDRcoex(sc1conf, sc1meta, input$sc1b2drX, input$sc1b2drY,   
             input$sc1b2inp1, input$sc1b2inp2, input$sc1b2sub1, input$sc1b2sub2, 
             "sc1gexpr.h5", sc1gene, 
             input$sc1b2siz, input$sc1b2col1, input$sc1b2ord1, 
             input$sc1b2fsz, input$sc1b2asp, input$sc1b2txt) 
  }) 
  output$sc1b2oup1.ui <- renderUI({ 
    plotOutput("sc1b2oup1", height = pList2[input$sc1b2psz]) 
  }) 
  output$sc1b2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                   input$sc1b2inp1,"_",input$sc1b2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1b2oup1.h, width = input$sc1b2oup1.w, useDingbats = FALSE, 
      plot = scDRcoex(sc1conf, sc1meta, input$sc1b2drX, input$sc1b2drY,  
                      input$sc1b2inp1, input$sc1b2inp2, input$sc1b2sub1, input$sc1b2sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1b2siz, input$sc1b2col1, input$sc1b2ord1, 
                      input$sc1b2fsz, input$sc1b2asp, input$sc1b2txt) ) 
    }) 
  output$sc1b2oup1.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                   input$sc1b2inp1,"_",input$sc1b2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1b2oup1.h, width = input$sc1b2oup1.w, 
      plot = scDRcoex(sc1conf, sc1meta, input$sc1b2drX, input$sc1b2drY,  
                      input$sc1b2inp1, input$sc1b2inp2, input$sc1b2sub1, input$sc1b2sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1b2siz, input$sc1b2col1, input$sc1b2ord1, 
                      input$sc1b2fsz, input$sc1b2asp, input$sc1b2txt) ) 
    }) 
  output$sc1b2oup2 <- renderPlot({ 
    scDRcoexLeg(input$sc1b2inp1, input$sc1b2inp2, input$sc1b2col1, input$sc1b2fsz) 
  }) 
  output$sc1b2oup2.ui <- renderUI({ 
    plotOutput("sc1b2oup2", height = "300px") 
  }) 
  output$sc1b2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                   input$sc1b2inp1,"_",input$sc1b2inp2,"_leg.pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = 3, width = 4, useDingbats = FALSE, 
      plot = scDRcoexLeg(input$sc1b2inp1, input$sc1b2inp2, input$sc1b2col1, input$sc1b2fsz) ) 
    }) 
  output$sc1b2oup2.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                   input$sc1b2inp1,"_",input$sc1b2inp2,"_leg.png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = 3, width = 4, 
      plot = scDRcoexLeg(input$sc1b2inp1, input$sc1b2inp2, input$sc1b2col1, input$sc1b2fsz) ) 
    }) 
  output$sc1b2.dt <- renderDataTable({ 
    ggData = scDRcoexNum(sc1conf, sc1meta, input$sc1b2inp1, input$sc1b2inp2, 
                         input$sc1b2sub1, input$sc1b2sub2, "sc1gexpr.h5", sc1gene) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("percent"), digits = 2) 
  }) 
  
  
  ### Plots for tab c1 
  output$sc1c1sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1c1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1c1sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1c1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1c1sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1c1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1c1oup <- renderPlot({
    scVioBox(sc1conf, sc1meta, input$sc1c1inp1, input$sc1c1inp1b, input$sc1c1inp2,
             input$sc1c1sub1, input$sc1c1sub2,
             "sc1gexpr.h5", sc1gene, input$sc1c1typ, input$sc1c1pts,
             input$sc1c1siz, input$sc1c1fsz)
  }) 
  output$sc1c1oup.ui <- renderUI({ 
    plotOutput("sc1c1oup", height = pList2[input$sc1c1psz]) 
  }) 
  output$sc1c1oup.pdf <- downloadHandler(
    filename = function() { paste0("sc1",input$sc1c1typ,"_",input$sc1c1inp1,"_",
                                   input$sc1c1inp2,".pdf") },
    content = function(file) { ggsave(
      file, device = "pdf", height = input$sc1c1oup.h, width = input$sc1c1oup.w, useDingbats = FALSE,
      plot = scVioBox(sc1conf, sc1meta, input$sc1c1inp1, input$sc1c1inp1b, input$sc1c1inp2,
                      input$sc1c1sub1, input$sc1c1sub2,
                      "sc1gexpr.h5", sc1gene, input$sc1c1typ, input$sc1c1pts,
                      input$sc1c1siz, input$sc1c1fsz) )
    }) 
  output$sc1c1oup.png <- downloadHandler(
    filename = function() { paste0("sc1",input$sc1c1typ,"_",input$sc1c1inp1,"_",
                                   input$sc1c1inp2,".png") },
    content = function(file) { ggsave(
      file, device = "png", height = input$sc1c1oup.h, width = input$sc1c1oup.w,
      plot = scVioBox(sc1conf, sc1meta, input$sc1c1inp1, input$sc1c1inp1b, input$sc1c1inp2,
                      input$sc1c1sub1, input$sc1c1sub2,
                      "sc1gexpr.h5", sc1gene, input$sc1c1typ, input$sc1c1pts,
                      input$sc1c1siz, input$sc1c1fsz) )
    }) 
  
  
  ### Plots for tab c2 
  output$sc1c2sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1c2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1c2sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1c2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1c2sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1c2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1c2oup <- renderPlot({ 
    scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,  
           input$sc1c2sub1, input$sc1c2sub2, 
           input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz) 
  }) 
  output$sc1c2oup.ui <- renderUI({ 
    plotOutput("sc1c2oup", height = pList2[input$sc1c2psz]) 
  }) 
  output$sc1c2oup.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1c2typ,"_",input$sc1c2inp1,"_",  
                                   input$sc1c2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1c2oup.h, width = input$sc1c2oup.w, useDingbats = FALSE, 
      plot = scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,  
                    input$sc1c2sub1, input$sc1c2sub2, 
                    input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz) ) 
    }) 
  output$sc1c2oup.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1c2typ,"_",input$sc1c2inp1,"_",  
                                   input$sc1c2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1c2oup.h, width = input$sc1c2oup.w, 
      plot = scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,  
                    input$sc1c2sub1, input$sc1c2sub2, 
                    input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz) ) 
    }) 
  
  
  ### Plots for tab d1 
  output$sc1d1sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1d1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1d1sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1d1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1d1sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1d1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1d1oupTxt <- renderUI({ 
    geneList = scGeneList(input$sc1d1inp, sc1gene) 
    if(nrow(geneList) > 50){ 
      HTML("More than 50 input genes! Please reduce the gene list!") 
    } else { 
      oup = paste0(nrow(geneList[present == TRUE]), " genes OK and will be plotted") 
      if(nrow(geneList[present == FALSE]) > 0){ 
        oup = paste0(oup, "<br/>", 
                     nrow(geneList[present == FALSE]), " genes not found (", 
                     paste0(geneList[present == FALSE]$gene, collapse = ", "), ")") 
      } 
      HTML(oup) 
    } 
  }) 
  output$sc1d1oup <- renderPlot({ 
    scBubbHeat(sc1conf, sc1meta, input$sc1d1inp, input$sc1d1grp, input$sc1d1plt, 
               input$sc1d1sub1, input$sc1d1sub2, "sc1gexpr.h5", sc1gene, 
               input$sc1d1scl, input$sc1d1row, input$sc1d1col, 
               input$sc1d1cols, input$sc1d1fsz) 
  }) 
  output$sc1d1oup.ui <- renderUI({ 
    plotOutput("sc1d1oup", height = pList3[input$sc1d1psz]) 
  }) 
  output$sc1d1oup.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1d1plt,"_",input$sc1d1grp,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1d1oup.h, width = input$sc1d1oup.w, 
      plot = scBubbHeat(sc1conf, sc1meta, input$sc1d1inp, input$sc1d1grp, input$sc1d1plt, 
                        input$sc1d1sub1, input$sc1d1sub2, "sc1gexpr.h5", sc1gene, 
                        input$sc1d1scl, input$sc1d1row, input$sc1d1col, 
                        input$sc1d1cols, input$sc1d1fsz, save = TRUE) ) 
    }) 
  output$sc1d1oup.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1d1plt,"_",input$sc1d1grp,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1d1oup.h, width = input$sc1d1oup.w, 
      plot = scBubbHeat(sc1conf, sc1meta, input$sc1d1inp, input$sc1d1grp, input$sc1d1plt, 
                        input$sc1d1sub1, input$sc1d1sub2, "sc1gexpr.h5", sc1gene, 
                        input$sc1d1scl, input$sc1d1row, input$sc1d1col, 
                        input$sc1d1cols, input$sc1d1fsz, save = TRUE) ) 
    }) 
  ### Plots for tab e1 
  
  selected_cells_rv <- reactive({
    ids <- input$sel_ingroup_keys
    if (is.null(ids) || length(ids) == 0) NULL else ids
  })
  outgroup_cells_rv <- reactive({
    ids <- input$sel_outgroup_keys
    if (is.null(ids) || length(ids) == 0) NULL else ids
  })
  marker_tbl_rv     <- reactiveVal(NULL)
  sel_gene_rv       <- reactiveVal(NULL)   # gene to overlay on UMAP (NULL = metadata mode)

  # Click a row in the marker table → show that gene's expression on UMAP.
  # Click the same row again (deselect) → revert to metadata coloring.
  observeEvent(input$sel_markers_tbl_rows_selected, {
    row_idx <- input$sel_markers_tbl_rows_selected
    if (length(row_idx) == 1 && !is.null(marker_tbl_rv())) {
      sel_gene_rv(as.character(marker_tbl_rv()$feature[row_idx]))
    } else {
      sel_gene_rv(NULL)
    }
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  #### Tab e1 specific functions

  get_palette_for_meta <- function(meta_col, sc1conf, sc1meta) {
    # find the row in sc1conf
    row <- sc1conf[UI == meta_col]
    fcl <- row$fCL   # e.g. "red|green|blue|orange"
    fid <- row$fID   # e.g. "Cluster1|Cluster2|Cluster3|Cluster4"
    
    # if there is no fixed‐colour list, return NULL
    if ( is.na(fcl) || !nzchar(fcl) ) return(NULL)
    
    # split into R vectors
    raw_cols  <- strsplit(fcl,  "\\|")[[1]]
    raw_levels<- strsplit(fid, "\\|")[[1]]
    
    # which of those levels actually appear in the data?
    present_levels <- intersect(raw_levels, unique(sc1meta[[ meta_col ]]))
    if (length(present_levels)==0) return(NULL)
    
    # subset and match colours to exactly the same positions
    sel_idx <- match(present_levels, raw_levels)
    pal_cols<- raw_cols[ sel_idx ]
    names(pal_cols) <- present_levels
    
    # also coerce the metadata column into a factor (most are already)
    sc1meta[[ meta_col ]] <<- factor(sc1meta[[ meta_col ]],
                                     levels = present_levels)
    
    pal_cols
  }
  
  runEnrichr <- function(genes, 
                        libraries = c("MSigDB_Hallmark_2020",
                                      "GO_Biological_Process_2025",
                                      "GO_Molecular_Function_2025",
                                      "KEGG_2021_Human"),
                        top_n     = 5) {
    # collapse your gene vector into the \n‐delimited string that Enrichr expects
    genes_str <- paste(genes, collapse = "\n")
    payload   <- list(list = genes_str,
                      description = "Shiny selector top 200")

    # 1) submit
    resp1 <- POST("https://maayanlab.cloud/Enrichr/addList",
                  body   = payload,
                  encode = "multipart")
    stop_for_status(resp1)
    ulid <- fromJSON(content(resp1, as = "text", encoding = "UTF-8"))$userListId

    # 2) fetch each library
    out <- lapply(libraries, function(lib) {
      url <- sprintf("https://maayanlab.cloud/Enrichr/enrich?userListId=%s&backgroundType=%s",
                    ulid, lib)
      resp2 <- GET(url)
      stop_for_status(resp2)
      js <- fromJSON(content(resp2, as = "text", encoding = "UTF-8"))
      res <- js[[lib]]
      if (is.null(res)) return(NULL)
      # each entry is a length‐>6 list: rank, term, pval, odds, combined, [gene vector]
      df <- data.frame(
        Rank          = sapply(res, `[[`, 1),
        Term          = sapply(res, `[[`, 2),
        Pvalue        = sapply(res, `[[`, 3),
        OddsRatio     = sapply(res, `[[`, 4),
        CombinedScore = sapply(res, `[[`, 5),
        Genes         = I(lapply(res, `[[`, 6)),
        stringsAsFactors = FALSE
      )
      df$Library <- lib
      # keep only the top_n by P‐value
      df %>%
        arrange(Pvalue) %>%
        slice_head(n = top_n)
    })

    # bind and return (dropping any NULL libraries)
    bind_rows(out)

  }

  runSelectorMarkers <- function(group1, group2, quick = FALSE) {
    shiny::validate(shiny::need(requireNamespace("presto", quietly = TRUE), "Package 'presto' is required for marker finding."))
    group1 <- sort(unique(as.integer(group1)))
    group2 <- sort(unique(as.integer(group2)))
    group2 <- setdiff(group2, group1)

    shiny::validate(shiny::need(length(group1) > 5, "Ingroup must have at least 5 cells"))
    shiny::validate(shiny::need(length(group2) > 5, "Outgroup must have at least 5 cells"))

    full_mat  <- selectorExprSparse()
    sel_cols  <- sort(unique(c(group1, group2)))
    all_cells <- length(sel_cols) == ncol(full_mat)

    # Skip the column copy when comparing ingroup vs all other cells
    expr <- if (all_cells) full_mat else full_mat[, sel_cols, drop = FALSE]

    # Quick mode: restrict to variable features (~7x fewer tests)
    if (quick && length(var_features) > 0L) {
      vf_mask <- rownames(expr) %in% var_features
      if (sum(vf_mask) > 0L) expr <- expr[vf_mask, , drop = FALSE]
    }

    # Drop genes with zero expression across all selected cells;
    # they produce uninformative test results and slow down presto.
    # tabulate on @i is O(nnz) and avoids materialising a logical matrix.
    nz_per_gene <- tabulate(expr@i + 1L, nbins = nrow(expr))
    keep <- nz_per_gene > 0L
    if (any(!keep)) expr <- expr[keep, , drop = FALSE]

    y <- rep("outgroup", length(sel_cols))
    y[match(group1, sel_cols)] <- "ingroup"
    y <- factor(y, levels = c("ingroup", "outgroup"))

    res <- data.table::as.data.table(presto::wilcoxauc(X = expr, y = y))
    res <- res[group == "ingroup", .(feature, group, avgExpr, logFC, statistic, auc, pval, padj, pct_in, pct_out)]
    res <- res[order(-auc, pval, padj, -avgExpr)]
    res
  }

  # Enrichr Reactive value:
  enrichr_res_rv <- reactiveVal(NULL)
  enrichr_status_rv <- reactiveVal(list(state = "idle", message = "Enrichr results will appear here."))
  enrichr_job_state <- new.env(parent = emptyenv())
  enrichr_job_state$request_id <- 0L
  enrichr_job_state$active <- TRUE

  session$onSessionEnded(function() {
    enrichr_job_state$active <- FALSE
  })

  beginEnrichrRequest <- function() {
    enrichr_job_state$request_id <- enrichr_job_state$request_id + 1L
    enrichr_res_rv(NULL)
    enrichr_status_rv(list(state = "idle", message = "Enrichr results will appear here."))
    enrichr_job_state$request_id
  }

  launchEnrichrAsync <- function(top_genes, request_id) {
    if (length(top_genes) <= 1) {
      if (enrichr_job_state$active && identical(request_id, enrichr_job_state$request_id)) {
        enrichr_res_rv(NULL)
        enrichr_status_rv(list(
          state = "empty",
          message = "Not enough ranked genes were available to run Enrichr for the latest selection."
        ))
      }
      return(invisible(NULL))
    }

    enrichr_status_rv(list(
      state = "loading",
      message = "Enrichr is updating in the background..."
    ))

    session$userData$enrichr_promise <- promises::then(
      promises::future_promise(
        {
          runEnrichr(top_genes)
        },
        packages = c("httr", "jsonlite", "dplyr")
      ),
      onFulfilled = function(enr) {
        if (!enrichr_job_state$active || !identical(request_id, enrichr_job_state$request_id)) {
          return(NULL)
        }

        if (is.null(enr) || nrow(enr) == 0) {
          enrichr_res_rv(NULL)
          enrichr_status_rv(list(
            state = "empty",
            message = "Enrichr returned no results for the latest selection."
          ))
        } else {
          enrichr_res_rv(enr)
          enrichr_status_rv(list(state = "ready", message = NULL))
        }
        NULL
      },
      onRejected = function(error) {
        if (!enrichr_job_state$active || !identical(request_id, enrichr_job_state$request_id)) {
          return(NULL)
        }

        err_msg <- conditionMessage(error)
        is_network <- grepl("resolve|connect|timeout|refused", err_msg, ignore.case = TRUE)
        user_msg <- if (is_network) {
          paste0("Could not reach Enrichr (network error: ", err_msg, "). Check your internet connection.")
        } else {
          paste0("Enrichr request failed: ", err_msg)
        }

        enrichr_res_rv(NULL)
        enrichr_status_rv(list(state = "error", message = user_msg))
        shiny::showNotification(user_msg, type = "error", duration = 8)
        NULL
      }
    )

    invisible(NULL)
  }

  # Resolve UMAP column names once using sc1conf (same source as all other tabs)
  sel_umap_x <- sc1conf[dimred == TRUE & grepl("UMAP", UI, ignore.case = TRUE)]$ID[1]
  sel_umap_y <- sc1conf[dimred == TRUE & grepl("UMAP", UI, ignore.case = TRUE)]$ID[2]

  output$sel_umap <- renderPlotly({
    req(input$sel_meta_col)

    gene     <- sel_gene_rv()
    meta_col <- input$sel_meta_col

    # Base coordinates shared by both modes
    plot_df <- data.frame(
      x      = sc1meta[[sel_umap_x]],
      y      = sc1meta[[sel_umap_y]],
      cellid = sc1meta$cellid,
      stringsAsFactors = FALSE
    )

    if (!is.null(gene) && gene %in% names(sc1gene)) {
      # --- Expression mode: color by gene expression ---
      if (!is.null(.expr_cache$mat) && gene %in% rownames(.expr_cache$mat)) {
        expr_vals <- as.numeric(.expr_cache$mat[gene, ])
      } else {
        gene_idx <- sc1gene[[gene]]
        h5file   <- H5File$new("sc1gexpr.h5", mode = "r")
        on.exit(h5file$close_all(), add = TRUE)
        expr_vals <- as.numeric(h5file[["grp"]][["data"]]$read(
          args = list(gene_idx, quote(expr = ))))
      }
      plot_df$expr <- pmax(expr_vals, 0)

      p <- plot_ly(
        data   = plot_df,
        x      = ~x,
        y      = ~y,
        type   = "scattergl",
        mode   = "markers",
        marker = list(size = 5, opacity = 0.7,
                      color     = ~expr,
                      colorscale = list(c(0, "#e0e0e0"), c(0.01, "#fee0d2"),
                                        c(0.25, "#fc9272"), c(0.5, "#de2d26"),
                                        c(1, "#67000d")),
                      colorbar  = list(title = gene, len = 0.5)),
        text      = ~paste0(gene, ": ", round(expr, 2)),
        hoverinfo = "text",
        key       = ~cellid,
        customdata = ~cellid,
        source    = "select_umap"
      )
      p_title <- paste0("UMAP \u2014 ", gene, " expression")
    } else {
      # --- Metadata mode: color by categorical metadata ---
      pal <- get_palette_for_meta(meta_col, sc1conf, sc1meta)
      plot_df$color <- sc1meta[[meta_col]]

      p <- plot_ly(
        data   = plot_df,
        x      = ~x,
        y      = ~y,
        type   = "scattergl",
        mode   = "markers",
        marker = list(size = 5, opacity = 0.6),
        color  = ~color,
        colors = pal,
        text     = ~paste0(meta_col, ": ", color),
        hoverinfo= "text",
        key      = ~cellid,
        customdata = ~cellid,
        source   = "select_umap"
      )
      p_title <- paste("UMAP coloured by", meta_col)
    }

    p %>%
      layout(
        title  = p_title,
        xaxis  = list(title       = sel_umap_x,
                      scaleanchor = "y",
                      scaleratio  = 1),
        yaxis  = list(title = sel_umap_y),
        showlegend = FALSE,
        dragmode = "lasso"
      ) %>%
      event_register("plotly_selected") %>%
      htmlwidgets::onRender("
function(el, x) {
  // Signal Shiny that the UMAP widget is on screen
  Shiny.setInputValue('sel_umap_rendered', true, {priority: 'event'});

  // Fade out the loading spinner overlay
  var spinner = document.getElementById('sel_umap_spinner');
  if (spinner) {
    spinner.style.transition = 'opacity 0.3s ease';
    spinner.style.opacity = '0';
    setTimeout(function() { spinner.style.display = 'none'; }, 300);
  }

  var shiftHeld = false;
  document.addEventListener('keydown', function(e) {
    if (e.key === 'Shift') shiftHeld = true;
  });
  document.addEventListener('keyup', function(e) {
    if (e.key === 'Shift') shiftHeld = false;
  });

  var ingroupKeys = [];
  var outgroupKeys = [];
  var ingroupShapeData = null;
  var outgroupShapes = [];
  var lastDrawnPath = null;

  Shiny.setInputValue('sel_ingroup_keys', null, {priority: 'event'});
  Shiny.setInputValue('sel_outgroup_keys', null, {priority: 'event'});

  function buildPathInfo(lassoPoints, range) {
    if (lassoPoints && lassoPoints.x && lassoPoints.x.length > 2) {
      var d = 'M' + lassoPoints.x[0] + ',' + lassoPoints.y[0];
      for (var i = 1; i < lassoPoints.x.length; i++) {
        d += 'L' + lassoPoints.x[i] + ',' + lassoPoints.y[i];
      }
      d += 'Z';
      var maxY = -Infinity, maxIdx = 0;
      for (var i = 0; i < lassoPoints.y.length; i++) {
        if (lassoPoints.y[i] > maxY) { maxY = lassoPoints.y[i]; maxIdx = i; }
      }
      return {svgPath: d, labelX: lassoPoints.x[maxIdx], labelY: maxY};
    } else if (range && range.x && range.y) {
      var x0 = range.x[0], x1 = range.x[1];
      var y0 = range.y[0], y1 = range.y[1];
      var d = 'M'+x0+','+y0+'L'+x1+','+y0+'L'+x1+','+y1+'L'+x0+','+y1+'Z';
      return {svgPath: d, labelX: (x0+x1)/2, labelY: Math.max(y0,y1)};
    }
    return null;
  }

  function updateVisuals() {
    var shapes = [];
    var annotations = [];
    if (ingroupShapeData) {
      shapes.push({
        type: 'path', path: ingroupShapeData.svgPath,
        fillcolor: 'rgba(205, 92, 92, 0.15)',
        line: {color: 'rgba(178, 34, 34, 0.6)', width: 2},
        xref: 'x', yref: 'y'
      });
      annotations.push({
        x: ingroupShapeData.labelX, y: ingroupShapeData.labelY,
        text: '<b>Ingroup</b>',
        showarrow: false,
        font: {color: 'rgb(178, 34, 34)', size: 14},
        xref: 'x', yref: 'y', yshift: 15
      });
    }
    for (var j = 0; j < outgroupShapes.length; j++) {
      var s = outgroupShapes[j];
      shapes.push({
        type: 'path', path: s.svgPath,
        fillcolor: 'rgba(100, 149, 237, 0.15)',
        line: {color: 'rgba(70, 130, 180, 0.6)', width: 2},
        xref: 'x', yref: 'y'
      });
      if (j === 0) {
        annotations.push({
          x: s.labelX, y: s.labelY,
          text: '<b>Outgroup</b>',
          showarrow: false,
          font: {color: 'rgb(70, 130, 180)', size: 14},
          xref: 'x', yref: 'y', yshift: 15
        });
      }
    }
    Plotly.relayout(el, {shapes: shapes, annotations: annotations});
  }

  // Capture the lasso/box path during each drag operation
  el.on('plotly_selecting', function(eventData) {
    if (eventData) {
      var p = buildPathInfo(eventData.lassoPoints, eventData.range);
      if (p) lastDrawnPath = p;
    }
  });

  el.on('plotly_selected', function(eventData) {
    if (!eventData || !eventData.points || eventData.points.length === 0) return;
    var allKeys = eventData.points.map(function(pt) { return pt.customdata; });
    // Use path captured during drag; fall back to event data
    var pathInfo = lastDrawnPath || buildPathInfo(eventData.lassoPoints, eventData.range);
    lastDrawnPath = null;

    if (shiftHeld && ingroupKeys.length > 0) {
      outgroupKeys = allKeys.filter(function(k) {
        return ingroupKeys.indexOf(k) < 0;
      });
      if (pathInfo) outgroupShapes.push(pathInfo);
    } else {
      ingroupKeys = allKeys;
      outgroupKeys = [];
      ingroupShapeData = pathInfo;
      outgroupShapes = [];
    }

    updateVisuals();
    Shiny.setInputValue('sel_ingroup_keys', ingroupKeys, {priority: 'event'});
    Shiny.setInputValue('sel_outgroup_keys', outgroupKeys, {priority: 'event'});
  });

  el.on('plotly_deselect', function() {
    ingroupKeys = [];
    outgroupKeys = [];
    ingroupShapeData = null;
    outgroupShapes = [];
    lastDrawnPath = null;
    updateVisuals();
    Shiny.setInputValue('sel_ingroup_keys', null, {priority: 'event'});
    Shiny.setInputValue('sel_outgroup_keys', null, {priority: 'event'});
  });
}
      ")
  })
  
  
  output$sel_ncells <- renderText({
    ids <- selected_cells_rv()
    out_ids <- outgroup_cells_rv()
    if (is.null(ids)) {
      "None selected"
    } else if (!is.null(out_ids) && length(out_ids) > 0) {
      paste0("Ingroup: ", length(ids), " cells | Outgroup: ", length(out_ids), " cells")
    } else {
      paste0(length(ids), " cells (vs all other cells)")
    }
  })
  
  observeEvent(input$do_marker, {
    sel_gene_rv(NULL)
    withProgress(message = "Processing markers...", value = 0, {
      enrichr_request_id <- beginEnrichrRequest()
      ids <- selected_cells_rv()
      shiny::validate(shiny::need(!is.null(ids) && length(ids) > 5, "Select at least 5 cells"))

      group1 <- lookupCellIndices(ids)
      outgroup_ids <- outgroup_cells_rv()
      if (!is.null(outgroup_ids) && length(outgroup_ids) > 0) {
        group2 <- lookupCellIndices(outgroup_ids)
        shiny::validate(shiny::need(length(group2) > 5, "Outgroup must have at least 5 cells"))
      } else {
        group2 <- setdiff(seq_len(nrow(sc1meta)), group1)
        req(length(group2) > 5)
      }

      if (preload_status_rv() == "idle") preload_status_rv("loading")
      incProgress(0.15, detail = "Loading expression matrix...")
      res <- runSelectorMarkers(group1, group2, quick = isTRUE(input$sel_quick_mode))
      if (preload_status_rv() != "ready") preload_status_rv("ready")
      incProgress(0.7, detail = "Computing wilcoxauc markers...")
      if (is.null(res) || nrow(res) == 0) {
        enrichr_res_rv(NULL)
        marker_tbl_rv(NULL)
        shiny::showNotification("No marker statistics could be computed for the selection.", type = "error")
        return(NULL)
      }
      # rank genes and run Enrichr on the top 200 genes (by AUC)
      res_sorted <- res[order(-auc, pval, padj, -avgExpr)]
      top200 <- head(res_sorted$feature, 200)
      marker_tbl_rv(res_sorted)
      res <- data.table::copy(head(res_sorted, 30))
      # Round numeric columns to 3 decimal places
      numeric_cols <- sapply(res, is.numeric)
      res[, (names(res)[numeric_cols]) := lapply(.SD, round, 3), .SDcols = numeric_cols]
      
      output$sel_markers_tbl <- DT::renderDT({
        DT::datatable(res, options = list(dom = "t", pageLength = 30),
                      rownames = FALSE, selection = "single")
      })

      launchEnrichrAsync(top200, enrichr_request_id)
    })
  })

  # --- Metadata-based selection ---

  # Optional subset filter: renderUI for filter-value multi-select
  output$sel_meta_filter_ui <- renderUI({
    filter_col <- input$sel_meta_filter_col
    if (is.null(filter_col) || filter_col == "") return(NULL)
    choices <- sort(unique(as.character(sc1meta[[filter_col]])))
    tagList(
      selectInput("sel_meta_filter_vals",
                  paste0("Include only (", filter_col, "):"),
                  choices = choices, selected = NULL, multiple = TRUE),
      uiOutput("sel_filter_count")
    )
  })

  # Cell indices surviving the optional filter
  filtered_cells_rv <- reactive({
    filter_col  <- input$sel_meta_filter_col
    filter_vals <- input$sel_meta_filter_vals
    if (is.null(filter_col) || filter_col == "" ||
        is.null(filter_vals) || length(filter_vals) == 0) {
      return(seq_len(nrow(sc1meta)))
    }
    which(as.character(sc1meta[[filter_col]]) %in% filter_vals)
  })

  # Show how many cells pass the filter
  output$sel_filter_count <- renderUI({
    n     <- length(filtered_cells_rv())
    total <- nrow(sc1meta)
    if (n < total) {
      tags$small(style = "color:#666; font-style:italic;",
                 sprintf("%s / %s cells in subset", format(n, big.mark = ","),
                         format(total, big.mark = ",")))
    }
  })

  meta_sel_col_rv <- reactive({ req(input$sel_meta_col); input$sel_meta_col })
  meta_sel_in_rv  <- reactive({ req(input$sel_meta_in);  input$sel_meta_in  })
  meta_sel_out_rv <- reactive({ req(input$sel_meta_out); input$sel_meta_out })

  # Ingroup choices — restricted to values present in filtered cells
  output$sel_meta_in_ui <- renderUI({
    req(input$sel_meta_col)
    cell_idx <- filtered_cells_rv()
    choices <- sort(unique(as.character(sc1meta[[input$sel_meta_col]][cell_idx])))
    selectInput("sel_meta_in", "In-group (one or more):",
                choices = choices, selected = NULL, multiple = TRUE)
  })

  # Outgroup choices — same restriction
  output$sel_meta_out_ui <- renderUI({
    req(input$sel_meta_col)
    cell_idx <- filtered_cells_rv()
    choices <- sort(unique(as.character(sc1meta[[input$sel_meta_col]][cell_idx])))
    selectInput("sel_meta_out", "Out-group (one or more):",
                choices = choices, selected = NULL, multiple = TRUE)
  })
  observeEvent(input$do_marker_meta, {
    sel_gene_rv(NULL)
    withProgress(message = "Processing markers (metadata selection)...", value = 0, {
      enrichr_request_id <- beginEnrichrRequest()
      meta_col <- meta_sel_col_rv()
      in_ids <- meta_sel_in_rv()
      out_ids <- meta_sel_out_rv()
      shiny::validate(shiny::need(length(in_ids) > 0, "Select at least one ingroup identity"))
      shiny::validate(shiny::need(length(out_ids) > 0, "Select at least one outgroup identity"))
      cell_idx <- filtered_cells_rv()
      group1 <- intersect(which(sc1meta[[meta_col]] %in% in_ids), cell_idx)
      group2 <- intersect(which(sc1meta[[meta_col]] %in% out_ids), cell_idx)
      shiny::validate(shiny::need(length(group1) > 5, "Ingroup must have at least 5 cells (check subset filter)"))
      shiny::validate(shiny::need(length(group2) > 5, "Outgroup must have at least 5 cells (check subset filter)"))

      if (preload_status_rv() == "idle") preload_status_rv("loading")
      incProgress(0.15, detail = "Loading expression matrix...")
      res <- runSelectorMarkers(group1, group2, quick = isTRUE(input$sel_quick_mode))
      if (preload_status_rv() != "ready") preload_status_rv("ready")
      incProgress(0.7, detail = "Computing wilcoxauc markers...")
      if (is.null(res) || nrow(res) == 0) {
        enrichr_res_rv(NULL)
        marker_tbl_rv(NULL)
        shiny::showNotification("No marker statistics could be computed for the metadata selection.", type = "error")
        return(NULL)
      }
      # rank genes by AUC before selecting top 200 for Enrichr
      res_sorted <- res[order(-auc, pval, padj, -avgExpr)]
      top200 <- head(res_sorted$feature, 200)
      marker_tbl_rv(res_sorted)
      res <- data.table::copy(head(res_sorted, 30))
      numeric_cols <- sapply(res, is.numeric)
      res[, (names(res)[numeric_cols]) := lapply(.SD, round, 3), .SDcols = numeric_cols]
      output$sel_markers_tbl <- DT::renderDT({
        DT::datatable(res, options = list(dom = "t", pageLength = 30),
                      rownames = FALSE, selection = "single")
      })

      launchEnrichrAsync(top200, enrichr_request_id)
    })
  })

  output$sel_marker_dl <- downloadHandler(
    filename = function() "selector_markers.csv",
    content = function(file) {
      data.table::fwrite(marker_tbl_rv(), file)
    }
  )   

  enrichrPlaceholder <- function() {
    status <- enrichr_status_rv()
    if (!is.null(status$message) && nzchar(status$message)) {
      status$message
    } else {
      "Enrichr results will appear here."
    }
  }
  
  # Enrichr reactive table and plot
  output$enrichr_table <- DT::renderDT({
    df <- enrichr_res_rv()
    required_cols <- c("Library", "Term", "Pvalue", "CombinedScore", "Genes")
    if (!is.data.frame(df) || nrow(df) == 0 || !all(required_cols %in% names(df))) {
      shiny::validate(shiny::need(FALSE, enrichrPlaceholder()))
    }
    # if you want a comma‐joined genes column:
    df$GeneList <- vapply(df$Genes, function(x) paste(unlist(x), collapse = ", "), character(1))
    out <- df[, c("Library", "Term", "Pvalue", "CombinedScore", "GeneList"), drop = FALSE]
    out$Pvalue <- as.numeric(formatC(out$Pvalue, format = "e", digits = 3))

    DT::datatable(
      out,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE)
    ) %>%
      DT::formatRound(
      columns = "CombinedScore",
      digits = 3
    )
  })

  output$enrichr_plot <- renderPlot({
    df <- enrichr_res_rv()
    required_cols <- c("Library", "Term", "CombinedScore")
    if (!is.data.frame(df) || nrow(df) == 0 || !all(required_cols %in% names(df))) {
      shiny::validate(shiny::need(FALSE, enrichrPlaceholder()))
    }
    df <- df %>%
      mutate(
        TermWrap  = stringr::str_wrap(Term, width = 40),
        LibShort  = recode(Library,
                          MSigDB_Hallmark_2020           = "Hallmarks",
                          GO_Biological_Process_2025     = "GO_BP",
                          GO_Molecular_Function_2025     = "GO_MF",
                          KEGG_2021_Human                 = "KEGG_2021")
      )
    ggplot(df, aes(x = CombinedScore, y = reorder(TermWrap, CombinedScore))) +
      geom_col(fill = "#c43c37", width = 0.7) +
      facet_wrap(~ LibShort, scales = "free", ncol = 2) +
      labs(x = "Combined Score", y = NULL) +
      theme_minimal(base_size = 14) +
      theme(
        strip.text    = element_text(face = "bold", size = 14),
        axis.text.y   = element_text(size = 10),
        panel.grid.major.y = element_blank()
      )
  })

})
