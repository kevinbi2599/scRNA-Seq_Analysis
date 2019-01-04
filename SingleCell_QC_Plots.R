library(MASS)
library(cowplot)
library(ggplot2)
library(viridis)
library(ggpubr)

#Define function get_density which returns a vector of 
get_density <- function(x, y, ...) {
dens <- MASS::kde2d(x, y, ...)
ix <- findInterval(x, dens$x)
iy <- findInterval(y, dens$y)
ii <- cbind(ix, iy)
return(dens$z[ii])
}

output_qcplots <- function(seurat_object){
  graphics.off()
  meta <- seurat_object@meta.data
  meta$log_nGene <- log(meta$nGene, 10)
  meta$log_nUMI <- log(meta$nUMI,10)
  meta$DENS_logUMI_lognGene <- get_density(meta$log_nUMI, meta$log_nGene, n = 500)
meta$DENS_logUMI_percent.mito <- get_density(meta$log_nUMI, meta$percent.mito, n = 500)
meta$DENS_logUMI_percent.rp <- get_density(meta$log_nUMI, meta$percent.rp, n = 500)
meta$DENS_percent.mito_percent.rp <- get_density(meta$percent.mito, meta$percent.rp, n = 500)
meta$DENS_lognGene_percent.rp <- get_density(meta$log_nGene, meta$percent.rp, n = 500)
percent.mito_filtercounter_cellslost <- length(subset(meta, percent.mito >= 0.2)$percent.mito)

my.formula <- log_nGene ~ log_nUMI
PLOT1_nGene_nUMI <- ggplot(meta) + geom_point(aes(log_nUMI, log_nGene, color = DENS_logUMI_lognGene), size = 0.1) + scale_color_viridis() + theme(legend.position="none") + labs(title = "Library Complexity", x = "log10(nUMI)", y = "log10(nGene)") + stat_cor(aes(log_nUMI, log_nGene),method = "pearson", label.x = 3.5, label.y = 2.5) + geom_smooth(aes(log_nUMI, log_nGene),method = "lm", se = FALSE, color = "red", linetype = "dashed", size = 0.5)

PLOT2_percentmito_nUMI <- ggplot(meta) + geom_point(aes(log_nUMI, percent.mito, color = DENS_logUMI_percent.mito), size = 0.1) + scale_color_viridis() + theme(legend.position="none") + geom_hline(yintercept=0.2, linetype="dashed", color = "red") + labs(title = paste("MT-RNA Filter:",percent.mito_filtercounter_cellslost,"cells lost"), x = "log10(nUMI)", y = "Fraction mitochondrial molecules")

PLOT3_percentrp_nUMI <- ggplot(meta) + geom_point(aes(log_nUMI, percent.rp, color = DENS_logUMI_percent.rp), size = 0.1) + scale_color_viridis() + theme(legend.position="none") + labs(title = paste("Ribosomal vs. Library Size"), x = "log10(nUMI)", y = "Fraction ribosomal molecules")

PLOT4_nGene_percentrp <- ggplot(meta) + geom_point(aes(log_nGene, percent.rp, color = DENS_lognGene_percent.rp), size = 0.1) + scale_color_viridis() + theme(legend.position="none") + labs(title = paste("Ribosomal vs. nGene"), x = "log(nGene)", y = "Fraction ribosomal molecules")

plot_grid(PLOT1_nGene_nUMI, PLOT2_percentmito_nUMI, PLOT3_percentrp_nUMI, PLOT4_nGene_percentrp, ncol = 2)
}


