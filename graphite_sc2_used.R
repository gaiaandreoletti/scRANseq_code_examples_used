library(Seurat)
library(assertthat)  # For when we want to make sanity checks
library(dplyr)       # For inline modification of matrices
library(cowplot)     # For pretty plots
library(ggplot2)     # For pretty plots
library(dittoSeq)    # For pretty, colorblind friendly plots
library(grid)        # For plotting multiple plots in one frame
library(gridExtra)   # For plotting multiple plots in one frame
library(scales)      # To access break formatting functions

setwd("/Users/gandreoletti/Library/CloudStorage/OneDrive-GraphiteBio/scRNAseq/")
setwd("/Users/gandreoletti/Downloads/GSE144568_RAW")
load("/Users/gandreoletti/Downloads/GSE144568_cell_information")
load("/Users/gandreoletti/Downloads/GSE144568_UMI_count")

# Creates generate_profile_plot function
generate_profile_plot <- function(sobj, feature1, feature2, feature1_binwidth=100,
                                  feature2_binwidth=100, visual_outlier_cutoff1=0.999,
                                  visual_outlier_cutoff2=0.999) {
  suppmsg <- assert_that(feature1 %in% colnames(sobj@meta.data), 
                         msg=paste0(feature1, " was not present in the metadata of sobj"))
  suppmsg <- assert_that(feature2 %in% colnames(sobj@meta.data), 
                         msg=paste0(feature2, " was not present in the metadata of sobj"))
  suppmsg <- assert_that(0 < visual_outlier_cutoff1 && visual_outlier_cutoff1 <=1.0, 
                         msg="visual_outlier_cutoff1 must be in the range (0,1]")
  suppmsg <- assert_that(0 < visual_outlier_cutoff2 && visual_outlier_cutoff2 <=1.0, 
                         msg="visual_outlier_cutoff2 must be in the range (0,1]")
  
  lay <- rbind(c(1,  1,2,2,2,2),
               c(1,  1,2,2,2,2),
               c(1,  1,2,2,2,2),
               c(NA,NA,3,3,3,3),
               c(NA,NA,3,3,3,3))
  
  lims = as.vector(
    c(quantile(sobj@meta.data[[feature1]], visual_outlier_cutoff1), 
      quantile(sobj@meta.data[[feature2]], visual_outlier_cutoff2)))
  xticks <- as.vector(quantile(sobj@meta.data[[feature1]], seq(0, max(0.9, visual_outlier_cutoff1), 0.1)))
  if (xticks[length(xticks)] != lims[1]) {
    xticks <- c(xticks, lims[1])
  }
  yticks <- as.vector(quantile(sobj@meta.data[[feature2]], seq(0, max(0.9, visual_outlier_cutoff2), 0.1)))
  if (yticks[length(yticks)] != lims[2]) {
    yticks <- c(yticks, lims[2])
  }
  
  main <- ggplot(sobj@meta.data, aes_string(x=feature1, y=feature2)) + 
    geom_point(aes(col="red"), size=0.5) + 
    xlim(NA, lims[1]) +
    ylim(NA, lims[2]) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()
    ) +
    NoLegend()
  
  y_hist <- ggplot(sobj@meta.data, aes_string(x=feature2)) + 
    geom_histogram(aes(col="red", fill="red"), binwidth=feature2_binwidth) +
    theme_bw() +
    theme(axis.title.x=element_blank(), 
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank()) +
    scale_x_continuous(limits=c(NA, lims[2]),
                       sec.axis = sec_axis(trans = ~., 
                                           breaks=yticks, labels=NULL)) +
    NoLegend() + 
    coord_flip() + 
    scale_y_reverse()
  
  
  x_hist <- ggplot(sobj@meta.data, aes_string(x=feature1)) + 
    geom_histogram(aes(col="red", fill="red"), binwidth=feature1_binwidth) +
    theme_bw() +
    theme(axis.title.y=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) +
    scale_x_continuous(limits=c(NA, lims[1]),
                       sec.axis = sec_axis(trans = ~., 
                                           breaks=xticks, labels=NULL)) +
    
    NoLegend() +
    scale_y_reverse()
  empty_plot <- ggplot() + theme_void()
  grid.arrange(grobs=list(y_hist, main, x_hist), layout_matrix = lay)
}

# A function to look at the distribution of a gene across the dataset
triplot <- function(sobj, features, reduction.use="umap", group.by="final_res") {
  plots <- list()
  layout = rbind(c(1,1,2,2),
                 c(1,1,2,2),
                 c(3,3,3,3),
                 c(3,3,3,3))
  for (f in features){
    tmp_plots <- list(p1=dittoDimPlot(group.by, 
                                      sobj, 
                                      reduction.use=reduction.use, 
                                      labels.repel=T, 
                                      do.label=F, 
                                      legend.show=F),
                      p2=dittoDimPlot(
                                      sobj, 
                                      f,
                                      reduction.use=reduction.use),
                      p3=dittoBoxPlot(
                                      sobj,
                                      f,
                                      group.by=group.by, 
                                      legend.show = F)
    )
    plots[[f]] <- grid.arrange(grobs=tmp_plots, layout_matrix=layout)
  }
  plots
}

# GLOBAL VARIABLES
PROJECT <- "graphite"
RIBO_GENES_FILE <- "/Users/gandreoletti/Library/CloudStorage/OneDrive-GraphiteBio/scRNAseq/auxiliary_files/ribo_genes.tsv"
CELL_CYCLE_GENS <- "/Users/gandreoletti/Library/CloudStorage/OneDrive-GraphiteBio/scRNAseq/auxiliary_files/cell_cycle_genes.tsv"

samples = read.table("/Users/gandreoletti/Downloads/GSE144568_RAW/balf_1 copy.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")

all.files <- samples$path
sample <- samples$sample
names(all.files) <- sample

i = 0
# for(sample_s in samples$sample){
#   i = i + 1
#   print(sample_s)
#   sample_i = samples %>% dplyr::filter(.,sample == sample_s)
#   print(samples$path[i])
#   data.i <- Read10X_h5(paste0("/Users/gandreoletti/Library/CloudStorage/OneDrive-GraphiteBio/scRNAseq/", samples$path[i]))
#   covid_sc <- CreateSeuratObject(counts = data.i , 
#                            project = sample_s, # sample name is best to put here. This become sobj@meta.data$orig.ident
#                            min.cells = 3, 
#                            min.features = 100)
# remove(data.i)
  
for(sample_s in samples$sample){
    i = i + 1
    print(sample_s)
    sample_i = samples %>% dplyr::filter(.,sample == sample_s)
    print(samples$path[i])
    data.i <- Read10X(paste0("/Users/gandreoletti/Downloads/GSE144568_RAW/", samples$path[i]))
    covid_sc <- CreateSeuratObject(counts = data.i , 
                                   project = sample_s, # sample name is best to put here. This become sobj@meta.data$orig.ident
                                   min.cells = 3, 
                                   min.features = 100)
remove(data.i)

covid_sc <- PercentageFeatureSet(covid_sc, 
                             pattern = "^Mt-",
                             col.name = "percent.mt")

ribo_genes <- read.table(RIBO_GENES_FILE, 
                         sep = "\t", header=TRUE, 
                         stringsAsFactors = FALSE)
ribo_genes <- ribo_genes[ribo_genes[['HUGO']] %in% rownames(covid_sc), ]
covid_sc <- PercentageFeatureSet(covid_sc,
                             features = ribo_genes[["HUGO"]],
                             col.name = "percent.ribo")
rm(ribo_genes)

cc_genes <- read.table(CELL_CYCLE_GENS, 
                       sep = "\t", 
                       header=T, 
                       stringsAsFactors = FALSE)

covid_sc <- CellCycleScoring(covid_sc,
                         s.features = cc_genes[cc_genes$stage=="G1-S", 'HUGO'],
                         g2m.features = cc_genes[cc_genes$stage=="G2-M", "HUGO"],
                         nbin = 12)
rm(cc_genes)


plots <- list(p1=generate_profile_plot(covid_sc,
                                       feature1 = "nCount_RNA",
                                       feature2 = "percent.mt",
                                       feature1_binwidth=100,
                                       feature2_binwidth=0.1),
              p2=generate_profile_plot(covid_sc,
                                       feature1 = "nCount_RNA",
                                       feature2 = "percent.ribo",
                                       feature1_binwidth=100,
                                       feature2_binwidth=0.1),
              p3=generate_profile_plot(covid_sc,
                                       feature1 = "nCount_RNA",
                                       feature2 = "nFeature_RNA",
                                       feature1_binwidth=100,
                                       feature2_binwidth=100),
              
              p4=generate_profile_plot(covid_sc,
                                       feature1 = "percent.ribo",
                                       feature2 = "percent.mt",
                                       feature1_binwidth=0.1,
                                       feature2_binwidth=0.1),
              
              p5=generate_profile_plot(covid_sc,
                                       feature1 = "nFeature_RNA",
                                       feature2 = "percent.mt",
                                       feature1_binwidth=100,
                                       feature2_binwidth=0.1),
              
              p6=generate_profile_plot(covid_sc,
                                       feature1 = "percent.ribo",
                                       feature2 = "nFeature_RNA",
                                       feature1_binwidth=0.1,
                                       feature2_binwidth=100)
)

df <- data.frame(
  cell_counts=seq(0, 1.01, 0.1)*dim(covid_sc@meta.data)[1],
  percent.mt=quantile(covid_sc@meta.data[["percent.mt"]], seq(0, 1.01, 0.1)),
  percent.ribo=quantile(covid_sc@meta.data[["percent.ribo"]], seq(0, 1.01, 0.1)),
  nFeature_RNA=quantile(covid_sc@meta.data[["nFeature_RNA"]], seq(0, 1.01, 0.1)),
  row.names=seq(0, 1.01, 0.1)
)
write.table(format(df, digits=2), file=paste0(sample_s,"covid_sc_dualscatter_pre.tsv"), row.names=T, 
            col.names=T, quote=F, sep="\t")

pdf(paste0(sample_s,"covid_sc_dualscatter_pre.pdf"), width = 21, height = 14)
print(CombinePlots(plots = plots, ncol=3))
dev.off()

keep_cells = colnames(covid_sc)[covid_sc$percent.mt<=10 & 
                              covid_sc$percent.ribo <=50 & 
                              covid_sc$nFeature_RNA > 300 &
                              covid_sc$nFeature_RNA < 5000 ]
covid_sc <- subset(covid_sc, cells = keep_cells)

plots <- list(p1=generate_profile_plot(covid_sc,
                                       feature1 = "nCount_RNA",
                                       feature2 = "percent.mt",
                                       feature1_binwidth=100,
                                       feature2_binwidth=0.1),
              p2=generate_profile_plot(covid_sc,
                                       feature1 = "nCount_RNA",
                                       feature2 = "percent.ribo",
                                       feature1_binwidth=100,
                                       feature2_binwidth=0.1),
              p3=generate_profile_plot(covid_sc,
                                       feature1 = "nCount_RNA",
                                       feature2 = "nFeature_RNA",
                                       feature1_binwidth=100,
                                       feature2_binwidth=100),
              
              p4=generate_profile_plot(covid_sc,
                                       feature1 = "percent.ribo",
                                       feature2 = "percent.mt",
                                       feature1_binwidth=0.1,
                                       feature2_binwidth=0.1),
              
              p5=generate_profile_plot(covid_sc,
                                       feature1 = "nFeature_RNA",
                                       feature2 = "percent.mt",
                                       feature1_binwidth=100,
                                       feature2_binwidth=0.1),
              
              p6=generate_profile_plot(covid_sc,
                                       feature1 = "percent.ribo",
                                       feature2 = "nFeature_RNA",
                                       feature1_binwidth=0.1,
                                       feature2_binwidth=100)
)

df <- data.frame(
  cell_counts=seq(0, 1.01, 0.1)*dim(covid_sc@meta.data)[1],
  percent.mt=quantile(covid_sc@meta.data[["percent.mt"]], seq(0, 1.01, 0.1)),
  percent.ribo=quantile(covid_sc@meta.data[["percent.ribo"]], seq(0, 1.01, 0.1)),
  nFeature_RNA=quantile(covid_sc@meta.data[["nFeature_RNA"]], seq(0, 1.01, 0.1)),
  row.names=seq(0, 1.01, 0.1)
)
write.table(format(df, digits=2), file=paste0(sample_s,"covid_sc_dualscatter_post.tsv"), row.names=T, 
            col.names=T, quote=F, sep="\t")

pdf(paste0(sample_s,"covid_sc_dualscatter_post.pdf"), width = 21, height = 14)
print(CombinePlots(plots = plots, ncol=3))
dev.off()

covid_sc <- SCTransform(covid_sc,
                    vars.to.regress = c("percent.mt", 
                                        "percent.ribo",
                                        "S.Score", 
                                        "G2M.Score"),
                    seed.use=21212,
                    verbose = FALSE)

saveRDS(covid_sc, file=paste0(sample_s,"covid_sc_SCTransformed.RDS"))

covid_sc <- RunPCA(covid_sc, 
               verbose = FALSE, 
               seed.use = 21212)

covid_sc <- RunUMAP(covid_sc,
                dims = 1:30,  # Num PCs to use
                # Default. Controls how UMAP balances local (low) 
                # versus global (large) structure in the data
                n.neighbors = 30,  
                # Default. Controls the size of the clusters. 
                # Should be smaller than spread
                min.dist = 0.3,   
                # Default. Controls the inter-cluster distances to some extent. 
                # Should be larger than min_dist
                spread = 1,  
                # Default. Can be used with b instead of using min.dist/spread
                a = NULL,  
                # Default. Can be used with a instead of using min.dist/spread
                b = NULL,  
                verbose = FALSE,
                seed.use = 21212)

covid_sc <- FindNeighbors(covid_sc,
                      dims = 1:30,  # Num PCs to use
                      k.param = 20,  # k for the knn algorithm
                      verbose = FALSE)
# Use the neighborhood graph to cluster the data
covid_sc <- FindClusters(covid_sc, 
                     verbose = TRUE,
                     # Use the Louvain algorithm. Leiden = 4
                     algorithm = 1,
                     # Clustering resolution. 
                     # Higher values (usually above 1 for Louvain) 
                     # lead to more clusters.
                     resolution = 0.8, 
                     random.seed = 21212)


# covid_sc_filtered <- subset(covid_sc, idents=c(1,4,6,7,8,9))
# saveRDS(covid_sc_filtered, file=paste0(sample_s,'covid_sc_filtered.RDS'))
saveRDS(covid_sc, file=paste0(sample_s,'HPSCCharac_sc.RDS'))
}

### merge

sobjs <- list()
i = 0
for(sample_s in samples$sample){
  i = i + 1
  print(sample_s)
sobjs[[sample_s]] <- readRDS(paste0(sample_s,"covid_sc_SCTransformed.RDS"))
}                             
                             
integration_features <- SelectIntegrationFeatures(object.list = sobjs, 
                                                  nfeatures = 3000)
options(future.globals.maxSize = 8000 * 1024^2)
sobjs <- PrepSCTIntegration(object.list = sobjs, 
                            anchor.features = integration_features, 
                            verbose = FALSE)

anchors <- FindIntegrationAnchors(object.list = sobjs, 
                                  normalization.method = "SCT", 
                                  anchor.features = integration_features,
                                  verbose = FALSE)
merged_data <- IntegrateData(anchorset = anchors, 
                             normalization.method = "SCT", 
                             verbose = FALSE)

# Run Umap, findneighbors and findclusters
saveRDS(merged_data, file="GSE144568_RAW-merged_data.RDS")
merged_data <- readRDS("GSE144568_RAW-merged_data.RDS")

merged_data <- FindVariableFeatures(merged_data, selection.method = "vst", nfeatures = 3000)

# Identify the 10 most highly variable genes
top50 <- head(VariableFeatures(merged_data), 50)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(merged_data)
plot2 <- LabelPoints(plot = plot1, points = top50, repel = TRUE)
# plot1
# plot2

all.genes <- rownames(merged_data)
pbmc <- ScaleData(merged_data, features = all.genes, vars.to.regress = c("percent.mt", "G2M.Score","S.Score","percent.ribo"))

cc_genes <- read.table("~/OneDrive - Graphite Bio/Barcoding/genes_c_trminus", stringsAsFactors=F, fill=T, header=F)
# 3 additional genes Beeke wanted me to test 2 Dec 2022
cc_genes <- read.table("~/OneDrive - Graphite Bio/Barcoding/genes_c_trminus_more_genes", stringsAsFactors=F, fill=T, header=F)

merged_data_ccgenes <- intersect(row.names(merged_data), cc_genes$V1)

merged_data <- RunPCA(merged_data, features = VariableFeatures(object = merged_data))

# Examine and visualize PCA results through text
print(merged_data[["pca"]], dims = 1:5, nfeatures = 5)

#graphically investigate the top 4 PCs and the genes which define them
VizDimLoadings(merged_data, dims = 1:4, reduction = "pca")

#graphically access the PCs on a biplot to help identify relational/spatial information 
DimPlot(merged_data, reduction = "pca")

#Plot a heatmap of the top 7 genes making up PC1
DimHeatmap(merged_data, dims = 1, cells = 500, balanced = TRUE)

#Use heatmap functionality to investigate the first 15 PCs and their top 7 defining genes
DimHeatmap(merged_data, dims = 1:15, cells = 500, balanced = TRUE)

# NOTE: This process can take a long time for big datasets or older processors. More approximate techniques such as those implemented in ElbowPlot can be used to reduce computation time

merged_data <- JackStraw(merged_data, num.replicate = 100)
merged_data <- ScoreJackStraw(merged_data, dims = 1:20)

### until PC13 are significant
JackStrawPlot(merged_data, dims = 1:20)
ElbowPlot(merged_data)
saveRDS(merged_data, file = "merged_data.rds")

## Step 7: Clustering
# merged_data <- FindNeighbors(merged_data, dims = 1:13)
merged_data <- FindNeighbors(merged_data,
                      dims = 1:13,  # Num PCs to use
                      k.param = 20,  # k for the knn algorithm
                      verbose = FALSE)
for (res in c(0.2, 0.4, 0.8, 1.2, 2.0)){
  # Use the neighborhood graph to cluster the data
  merged_data <- FindClusters(merged_data, 
                       verbose = TRUE,
                       # Use the Leiden algorithm. Leiden = 4
                       algorithm = 1,
                       # Clustering resolution. 
                       # Higher values (usually above 1 for Louvain) 
                       # lead to more clusters.
                       resolution = res, 
                       random.seed = 21212)
}


pdf("prova.pdf")
dittoDimPlot(merged_data, 'ident', reduction.use = "umap")
dittoDimPlot(merged_data, 'seurat_clusters', reduction.use = "umap")
dittoDimPlot(merged_data, 'integrated_snn_res.2', reduction.use = "umap")
dittoDimPlot(merged_data, 'integrated_snn_res.1.2',  reduction.use = "umap")
dittoDimPlot(merged_data, 'integrated_snn_res.0.8',  reduction.use = "umap")
dittoDimPlot(merged_data, 'integrated_snn_res.0.4',  reduction.use = "umap")
dittoDimPlot(merged_data, 'integrated_snn_res.0.2',  reduction.use = "umap")
dev.off()
# ## 4,6,8,1
# merged_data <- FindClusters(merged_data, resolution = 0.4)
# # Look at cluster IDs of the first 5 cells
head(Idents(merged_data), 5)

# UMAP now runs in entirely in R!
merged_data <- RunUMAP(merged_data, dims = 1:13)
# pbmc <- RunTSNE(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
# if umap is giving you trouble, replace "umap" with "tsne" for this workshop
# DimPlot(merged_data, 'integrated_snn_res.0.4', reduction = "umap")
pdf("umaps_diff_resolutions.pdf")
dittoDimPlot(merged_data, 'integrated_snn_res.2', reduction.use = "umap")
dittoDimPlot(merged_data, 'integrated_snn_res.1.2',  reduction.use = "umap")
dittoDimPlot(merged_data, 'integrated_snn_res.0.8',  reduction.use = "umap")
dittoDimPlot(merged_data, 'integrated_snn_res.0.4',  reduction.use = "umap")
dittoDimPlot(merged_data, 'integrated_snn_res.0.2',  reduction.use = "umap")
dev.off()

Idents(merged_data) <- merged_data@meta.data$integrated_snn_res.0.2

## add metadata
git_metadata <- read.table("meta.txt", header = T)
rownames(git_metadata) <- git_metadata$ID
git_metadata <- as.data.frame(git_metadata)

# Abbreviated snippets of AddMetaData
####add  sample info
# sample_info = as.data.frame(merged_data@meta.data)
sample_info = as.data.frame(colnames(merged_data))
colnames(sample_info) = c('ID')
rownames(sample_info) = sample_info$ID
sample_info$sample = merged_data@meta.data$orig.ident
sample_info = dplyr::left_join(sample_info,git_metadata, by = c("sample" = "GEO_IDs"))
rownames(sample_info) = sample_info$ID
merged_data = AddMetaData(object = merged_data, metadata = sample_info)

#Save
saveRDS(merged_data, "merged_data_metadata.rds")
merged_data <- readRDS("merged_data_metadata.rds")

pdf("umaps_diff_split.pdf")
DimPlot(merged_data, split.by = "Sample", na.value = element_blank())
DimPlot(merged_data, split.by = "Day", na.value = element_blank())
dev.off()
##### 0.8 ######
Idents(merged_data) <- merged_data@meta.data$integrated_snn_res.0.2
# Idents(merged_data) <- merged_data@meta.data$orig.ident
# pdf("umap_group_and_disease.pdf")
# DimPlot(merged_data, group.by  = "group")
# DimPlot(merged_data, group.by  = "disease")
# dev.off()

# find markers for every cluster compared to all remaining cells, report only the positive ones
# merged_data.markers <- FindAllMarkers(merged_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4, test.use= "MAST", latent.vars = "orig.ident" )
# merged_data.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

i = 0.8
merged_data_markers_i <- FindAllMarkers(merged_data, 
                               logfc.threshold=0.4, 
                               test.use="MAST",
                               min.pct = 0.1,
                               only.pos = T, 
                               random.seed = 21212)

write.table(merged_data_markers_i, 
            file= paste0("sobj_",i,"_markers_filtered.tsv"),
            sep='\t', 
            row.names = F, 
            col.names=T,
            quote=F)

##### 0.4 ######
Idents(merged_data) <- merged_data@meta.data$integrated_snn_res.0.4
# Idents(merged_data) <- merged_data@meta.data$orig.ident
# pdf("umap_group_and_disease.pdf")
# DimPlot(merged_data, group.by  = "group")
# DimPlot(merged_data, group.by  = "disease")
# dev.off()

# find markers for every cluster compared to all remaining cells, report only the positive ones
# merged_data.markers <- FindAllMarkers(merged_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4, test.use= "MAST", latent.vars = "orig.ident" )
# merged_data.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

i = 0.4
merged_data_markers_i <- FindAllMarkers(merged_data, 
                                        logfc.threshold=0.4, 
                                        test.use="MAST",
                                        min.pct = 0.1,
                                        only.pos = T, 
                                        random.seed = 21212)

##### 1.2 ######
Idents(merged_data) <- merged_data@meta.data$integrated_snn_res.1.2
# Idents(merged_data) <- merged_data@meta.data$orig.ident
# pdf("umap_group_and_disease.pdf")
# DimPlot(merged_data, group.by  = "group")
# DimPlot(merged_data, group.by  = "disease")
# dev.off()

# find markers for every cluster compared to all remaining cells, report only the positive ones
# merged_data.markers <- FindAllMarkers(merged_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4, test.use= "MAST", latent.vars = "orig.ident" )
# merged_data.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

i = 1.2
merged_data_markers_i <- FindAllMarkers(merged_data, 
                                        logfc.threshold=0.4, 
                                        test.use="MAST",
                                        min.pct = 0.1,
                                        only.pos = T, 
                                        random.seed = 21212)

write.table(merged_data_markers_i, 
            file= paste0("sobj_",i,"_markers_filtered.tsv"),
            sep='\t', 
            row.names = F, 
            col.names=T,
            quote=F)

write.table(merged_data_markers_i, 
            file= paste0("sobj_",i,"_markers_filtered.tsv"),
            sep='\t', 
            row.names = F, 
            col.names=T,
            quote=F)

##### 2.0 ######
Idents(merged_data) <- merged_data@meta.data$integrated_snn_res.2.0
# Idents(merged_data) <- merged_data@meta.data$orig.ident
# pdf("umap_group_and_disease.pdf")
# DimPlot(merged_data, group.by  = "group")
# DimPlot(merged_data, group.by  = "disease")
# dev.off()

# find markers for every cluster compared to all remaining cells, report only the positive ones
# merged_data.markers <- FindAllMarkers(merged_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4, test.use= "MAST", latent.vars = "orig.ident" )
# merged_data.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

i = 2.0
merged_data_markers_i <- FindAllMarkers(merged_data, 
                                        logfc.threshold=0.4, 
                                        test.use="MAST",
                                        min.pct = 0.1,
                                        only.pos = T, 
                                        random.seed = 21212)

write.table(merged_data_markers_i, 
            file= paste0("sobj_",i,"_markers_filtered.tsv"),
            sep='\t', 
            row.names = F, 
            col.names=T,
            quote=F)

write.table(merged_data_markers_i, 
            file= paste0("sobj_",i,"_markers_filtered.tsv"),
            sep='\t', 
            row.names = F, 
            col.names=T,
            quote=F)


pdf("umaps.pdf")
DimPlot(merged_data, 
        group.by='integrated_snn_res.0.8', 
        reduction = "umap", 
        label = T) + NoLegend()


DimPlot(merged_data, 
        group.by='integrated_snn_res.0.4', 
        reduction = "umap", 
        label = T) + NoLegend()
DimPlot(merged_data, 
        group.by='integrated_snn_res.1.2', 
        reduction = "umap", 
        label = T) + NoLegend()
DimPlot(merged_data, 
        group.by='integrated_snn_res.2', 
        reduction = "umap", 
        label = T) + NoLegend()

dev.off()
#####finding the optimal numbers of clusters #######
library(clustree)
clustree(merged_data, prefix = "integrated_snn_res.")
library(IKAP)
source("IKAP_Seurat3.R")
# run IKAP
sobj <- IKAP(merged_data, out.dir = "./IKAP")

# save the Seurat object with IKAP results
saveRDS(sobj,"./IKAP/sobj.rds")
################
merged_data_markers.8 <- read.table("sobj_0.8_markers_filtered.tsv", header = T, sep = "\t")
merged_data_markers1.2 <- read.table("sobj_1.2_markers_filtered.tsv", header = T, sep = "\t")

top5.08 <- merged_data_markers.8 %>% 
  arrange(cluster, -avg_logFC) %>% 
  group_by(cluster) %>% 
  top_n(5, wt=avg_logFC)

top5.1.2 <- merged_data_markers1.2 %>% 
  arrange(cluster, -avg_logFC) %>% 
  group_by(cluster) %>% 
  top_n(5, wt=avg_logFC)

png(file= "DoHeatmap_markers_0.8.png")
DoHeatmap(merged_data, features=top5.08$gene, size = 3)
# dev.off()
# png(file= "DoHeatmap_markers_1.2.png", width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
DoHeatmap(merged_data, features=top5.1.2$gene, size = 3)
dev.off()

################
merged_data@meta.data$final_res <- merged_data@meta.data$integrated_snn_res.0.2  # For convenience

### annotation
triplot <- function(sobj, features, reduction.use="umap", group.by="final_res") {
  plots <- list()
  layout = rbind(c(1,1,2,2),
                 c(1,1,2,2),
                 c(3,3,3,3),
                 c(3,3,3,3))
  for (f in features){
    tmp_plots <- list(p1=dittoDimPlot(sobj, group.by, 
                                      reduction.use=reduction.use, 
                                      labels.repel=T, 
                                      do.label=F, 
                                      legend.show=F),
                      p2=dittoDimPlot(sobj, f, 
                                      reduction.use=reduction.use,max = 5, min = -1),
                      p3=dittoBoxPlot(sobj, f, 
                                      group.by=group.by, 
                                      legend.show = F,   jitter.size = 0)
    )
    plots[[f]] <- grid.arrange(grobs=tmp_plots, layout_matrix=layout)
  }
  plots
}

triplot_noditto <- function(sobj, features, reduction.use="umap", group.by="final_res") {
  plots <- list()
  layout = rbind(c(1,1,2,2),
                 c(1,1,2,2),
                 c(3,3,3,3),
                 c(3,3,3,3))
  for (f in features){
    tmp_plots <- list(p1=DimPlot(sobj, 
                                 group.by=group.by, 
                                 label=F) + NoLegend(),
                      p2=FeaturePlot(sobj, 
                                     features = f,
                                     reduction=reduction.use, cols = c("grey","red"),min.cutoff = -1, max.cutoff = 5),
                      p3=VlnPlot(sobj,
                                 features = f, 
                                 group.by=group.by, assay = "SCT", pt.size = 0)  + NoLegend()
    )
    plots[[f]] <- grid.arrange(grobs=tmp_plots, layout_matrix=layout)
  }
  plots
}

# sobj <- merged_data
# marker_dimplot <- dittoDimPlot(merged_data, 'integrated_snn_res.0.8', 
#                                reduction.use = "umap", 
#                                labels.repel = T, 
#                                do.label = T, 
#                                legend.show = F)
# marker_featureplot <- dittoDimPlot(merged_data, 'PTPRC',  #CD45
#                                    reduction.use="umap", min = 0.1)

triplot(merged_data, "TUBA1B")  # CD45
triplot_noditto(merged_data, "PTPRC")  # CD45


triplot_noditto(merged_data, "TUBA1A")  # TUBA1A
triplot_noditto(merged_data, "TUBA1B")  # TUBA1B
triplot_noditto(merged_data, "GAPDH")  # GAPDH

pdf("FeaturePlot_markers_0.8.pdf")
FeaturePlot(merged_data, features ="CD3D", cols = c("grey","red"),min.cutoff = -1, max.cutoff = 5) ## Tcell
FeaturePlot(merged_data, features ="CD68", cols = c("grey","red"),min.cutoff = -1, max.cutoff = 5) ## Macrophages 
FeaturePlot(merged_data, features ="TPPP3", cols = c("grey","red"),min.cutoff = -1, max.cutoff = 5)## Epithelial-Ciliated 
FeaturePlot(merged_data, features ="HLA-DRB1", cols = c("grey","red"),min.cutoff = -1, max.cutoff = 5) ## Bcells, Myloyd cells
FeaturePlot(merged_data, features ="C1QC", cols = c("grey","red"),min.cutoff = -1, max.cutoff = 5) ##  macrophages
FeaturePlot(merged_data, features ="VCAN", cols = c("grey","red"),min.cutoff = -1, max.cutoff = 5) ## monocytes
FeaturePlot(merged_data, features ="G0S2", cols = c("grey","red"),min.cutoff = -1, max.cutoff = 5) ## Neutrophil and monocytes
FeaturePlot(merged_data, features ="TPPP3", cols = c("grey","red"),min.cutoff = -1, max.cutoff = 5) ## Epithelial-Ciliated 
FeaturePlot(merged_data, features ="KRT18", cols = c("grey","red"),min.cutoff = -1, max.cutoff = 5) ## Epithelial-Secretory
FeaturePlot(merged_data, features ="CD68", cols = c("grey","red"),min.cutoff = -1, max.cutoff = 5) ## Macrophages 
FeaturePlot(merged_data, features ="FCGR3B", cols = c("grey","red"),min.cutoff = -1, max.cutoff = 5) ## Neutrophils 
FeaturePlot(merged_data, features ="CD1C", cols = c("grey","red"),min.cutoff = -1, max.cutoff = 5) ## mDC 
FeaturePlot(merged_data, features ="LILRA4", cols = c("grey","red"),min.cutoff = -1, max.cutoff = 5) ## pDC 
FeaturePlot(merged_data, features ="TPSB2", cols = c("grey","red"),min.cutoff = -1, max.cutoff = 5) ## Mast cells 
FeaturePlot(merged_data, features ="CD3D", cols = c("grey","red"),min.cutoff = -1, max.cutoff = 5) ## T cells 
FeaturePlot(merged_data, features ="KLRD1", cols = c("grey","red"),min.cutoff = -1, max.cutoff = 5) ## NK cells 
FeaturePlot(merged_data, features ="MS4A1", cols = c("grey","red"),min.cutoff = -1, max.cutoff = 5) ## B cells 
FeaturePlot(merged_data, features ="IGHG4", cols = c("grey","red"),min.cutoff = -1, max.cutoff = 5) ## Plasma cells 
dev.off()

pdf("triplot_noditto_markers_0.8.pdf")
triplot_noditto(merged_data, "CD3D") ## Tcell
triplot_noditto(merged_data, "CD68") ## Macrophages 
triplot_noditto(merged_data, "TPPP3")## Epithelial-Ciliated 
triplot_noditto(merged_data, "HLA-DRB1") ## Bcells, Myloyd cells
triplot_noditto(merged_data, "C1QC") ##  macrophages
triplot_noditto(merged_data, "VCAN") ## monocytes
triplot_noditto(merged_data, "G0S2") ## 
triplot_noditto(merged_data, "TPPP3") ## Epithelial-Ciliated 
triplot_noditto(merged_data, "KRT18") ## Epithelial-Secretory
triplot_noditto(merged_data, "CD68") ## Macrophages 
triplot_noditto(merged_data, "FCGR3B") ## Neutrophils 
triplot_noditto(merged_data, "CD1C") ## mDC 
triplot_noditto(merged_data, "LILRA4") ## pDC 
triplot_noditto(merged_data, "TPSB2") ## Mast cells 
triplot_noditto(merged_data, "CD3D") ## T cells 
triplot_noditto(merged_data, "KLRD1") ## NK cells 
triplot_noditto(merged_data, "MS4A1") ## B cells 
triplot_noditto(merged_data, "IGHG4") ## Plasma cells 
dev.off()
# Epithelial-Ciliated (TPPP3)
# Epithelial-Secretory (KRT18)
# Macrophages (CD68)
# Neutrophils (FCGR3B)
# mDC (CD1C)
# pDC (LILRA4)
# Mast cells (TPSB2)
# T cells (CD3D)
# NK cells (KLRD1)
# B cells (MS4A1)
# Plasma cells (IGHG4)
# print(marker_dimplot)
# print(marker_featureplot)
# print(marker_boxplot)
# 
# marker_dimplot <- DimPlot(merged_data, 
#                           group.by='integrated_snn_res.0.8', 
#                           reduction = "umap", 
#                           label = T) + NoLegend()
# 
# marker_featureplot <- FeaturePlot(merged_data, 
#                                   features = 'PTPRC',
#                                   reduction="umap", min.cutoff = -1, max.cutoff = 5)
# marker_boxplot <- VlnPlot(merged_data, 
#                           features = 'PTPRC',
#                           group.by='final_res', pt.size = 0) + NoLegend()

# print(marker_dimplot)
# print(marker_featureplot)
# print(marker_boxplot)
# triplot_noditto(merged_data, features = 'DKK1', reduction.use = 'umap', group.by = 'final_res')

id.gene.list = c('IL7R', 'CCR7',        #Naive CD4+ T
                 'IL7R',  #Memory CD4+
                 'CD14', 'LYZ',     #CD14+ Mono
                 'MS4A1',                 #B
                 'CD8A',              #CD8+ T
                 'FCGR3A', 'MS4A7', #FCGR3A+ Mono
                 'GNLY', 'NKG7',        #NK
                 'FCER1A', 'CST3',  #DC
                 'PPBP')                  #Platelet

FeaturePlot(merged_data, features = id.gene.list, min.cutoff = -1, max.cutoff = 5)

id.gene.list = c("TUBA1B" ,   "HIST1H2BJ", "LAMP1",     "SMC1A"   ,  "ST6GAL1" ,  "FUS")
id.gene.list = c("TUBA1B" ,   "TUBA1A", "GAPDH")

library(ggpubr)
VlnPlot(object = merged_data, features =c(id.gene.list), assay = "SCT", pt.size = 0.1, group.by = 'final_res') #log = T 
VlnPlot(object = merged_data, features =c(id.gene.list), assay = "SCT", pt.size = 0.1, split.by  = 'Day') + stat_compare_means() #log = T 
VlnPlot(object = merged_data, features =c(id.gene.list), assay = "SCT", pt.size = 0.1, split.by  = 'Sample')  #log = T 

  VlnPlot(object = merged_data, features =c(id.gene.list), assay = "SCT", pt.size = 0.1, split.by  = 'final_res') + #log = T 
  
triplot_noditto(merged_data, features =c(id.gene.list))  # CD45

  
comparisons <- list(c("Day0", "Day5"))
p_case1(gene_signature = id.gene.list, file_name = "test", test_sign = comparisons)

TUBA1A_TUBA1B_GAPDH

VlnPlot(object = merged_data, features ="TUBA1B", assay = "SCT", pt.size = 0.1,group.by = 'Day')+ stat_compare_means(comparisons = list(c("Day0", "Day5")), label = "p.signif")
VlnPlot(object = merged_data, features ="TUBA1A", assay = "SCT", pt.size = 0.1,group.by = 'Day') + stat_compare_means()
VlnPlot(object = merged_data, features ="GAPDH", assay = "SCT", pt.size = 0.1,group.by = 'Day') + stat_compare_means()

my_comparisons <- list(c("Sample1", "Sample2"),c("Sample1", "Sample3"),c("Sample1", "Sample4"),c("Sample2", "Sample3"),c("Sample2", "Sample4"),c("Sample3", "Sample4"))
VlnPlot(object = merged_data, features ="TUBA1B", assay = "SCT", pt.size = 0.1,group.by = 'Sample', y.max = 10) + stat_compare_means(comparisons = my_comparisons, label = "p.signif")  + stat_compare_means(label.y = 10)
VlnPlot(object = merged_data, features ="TUBA1A", assay = "SCT", pt.size = 0.1,group.by = 'Sample' , y.max = 10) + stat_compare_means(comparisons = my_comparisons, label = "p.signif")  + stat_compare_means(label.y = 10)
VlnPlot(object = merged_data, features ="GAPDH", assay = "SCT", pt.size = 0.1,group.by = 'Sample',  y.max = 10) + stat_compare_means(comparisons = my_comparisons, label = "p.signif")  + stat_compare_means(label.y = 10)
VlnPlot(object = merged_data, features ="FUS", assay = "SCT", pt.size = 0.1,group.by = 'Sample',  y.max = 10) + stat_compare_means(comparisons = my_comparisons, label = "p.signif")  + stat_compare_means(label.y = 10)

VlnPlot(object = merged_data, features =c(id.gene.list), assay = "SCT", pt.size = 0.1, group.by  = 'Sample') + stat_compare_means() #log = T 

nonVargenes <-setdiff(x = rownames(x = merged_data), y = VariableFeatures(object = merged_data))
write.csv(nonVargenes, "non-var-genes.csv", row.names = FALSE)

####### annotations #######
annotations <- c(
  "Macrophages 1",#0
  "Monocytes 1",#1
  "Macrophages 2",#2
  "Macrophages 3",#3
  "Non classical monocytes",#4
  "Macrophages 4",#5
  "Macrophages 5",#6
  "Neutrophils",#7
  "Cd8 T cells ",#8
  "Macrophages (m2)",#9
  "Epithelial cells 2",#10
  "B cells- Plasma cells",#11
  "Basophils",#12
  "CD4 T cells",#13
  "Eosinophils",#14
  "T cell ",#15
  "Dividing T cell ",#16
  "Epithelial cells",#17
  "Epithelial cells ISG ",#18
  "DCs",#19
  "Lung epithelium",#20
  "Alveolar epithelial cells, mixture of cells"#21
)

Idents(merged_data) <- merged_data@meta.data$integrated_snn_res.0.8 
#DimPlot(merged_data, label = TRUE)+ NoLegend()
names(annotations) <- levels(merged_data)
merged_data <- RenameIdents(merged_data, annotations)
DimPlot(merged_data, label = TRUE) + NoLegend()
saveRDS(merged_data, "merged_data_annotated.rds")
# merged_data <- readRDS("merged_data_metadata.rds")
merged_data <- readRDS("merged_data_annotated.rds")

### Ansun genes ####
selected_genes <- read.table("selected_genes.txt", header = F)
selected_genes <-selected_genes$V1
selected_genes <- as.character(selected_genes)

DotPlot(merged_data, features = selected_genes, split.by = "group", cols = c("blue", "red", "green")) + RotatedAxis()

DotPlot(merged_data, features = rev(c(selected_genes)), cols = c("blue", "red", "green"), dot.scale = 8, 
        split.by = "group") + RotatedAxis()

DotPlot(merged_data, features = rev(c(selected_genes)), cols = c("blue", "red", "green"), dot.scale = 8, 
        split.by = "group") + RotatedAxis()


png(file="dot_plot_byGroup.png", width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
DotPlot(merged_data, features = selected_genes, assay = "SCT", split.by = "group",cols = c("blue", "red", "green")) + RotatedAxis()
dev.off()

DotPlot(merged_data, features = selected_genes, assay = "SCT", group.by  = "group") + RotatedAxis()


dpi = 300
png(file="dot_plot.png", width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
DotPlot(merged_data, features = selected_genes, assay = "SCT") + RotatedAxis()
dev.off()

png(file=paste('filter/',sample_s,"_qc.png",sep=''), width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')

DoHeatmap(merged_data, features = selected_genes, size = 3)

#### violin plots tests #####
pdf("violin.pdf")
VlnPlot(object = merged_data, features = "CD3D", assay = "SCT", pt.size = 0.1, split.by  = "disease", split.plot = TRUE) #log = T 

VlnPlot(object = merged_data, features = "CD3D", assay = "SCT", pt.size = 0.1, split.by  = "group") #log = T 


VlnPlot(object = merged_data, features = c(selected_genes), assay = "SCT", pt.size = 0.1, split.by  = "disease", split.plot = TRUE) #log = T 

VlnPlot(object = merged_data, features =c(selected_genes), assay = "SCT", pt.size = 0.1, split.by  = "group", split.plot = TRUE) #log = T 

dev.off()
VlnPlot(object = merged_data, features =c(selected_genes), assay = "SCT", pt.size = 0.1,  group.by = "disease") #log = T 

##### stack violin plots ######
Idents(merged_data) <- merged_data@meta.data$integrated_snn_res.0.8
merged_data <- RenameIdents(merged_data, annotations)

library(Seurat)
library(patchwork)
library(ggplot2)

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, assay = "SCT", split.by  = "disease", pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "left", 
          # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


## not available "SIGLECL1","SIGLEC" SIGLEC16???
Siglec_family_members_1 <- c( "CD22","SIGLEC1",
"CD33", "MAG",
"SIGLEC5", "SIGLEC6",
"SIGLEC7")
StackedVlnPlot(obj = merged_data, features = Siglec_family_members_1)

Siglec_family_members_2 <- c("SIGLEC8",
"SIGLEC9","SIGLEC10","SIGLEC11","SIGLEC12","SIGLEC14","SIGLEC15")
StackedVlnPlot(obj = merged_data, features = Siglec_family_members_2)

c_type_lectin_receptors <- c("MRC1",
"LY75",
"CLEC7A",
"CLEC6A",
"CLEC4E",
"CD209",
"CLEC10A",
"CLEC9A",
"CD207")
StackedVlnPlot(obj = merged_data, features = c_type_lectin_receptors)

TLR_family_members <- c("TLR1",
"TLR10",
"TLR2",
"TLR3",
"TLR5",
"TLR6",
"TLR7",
"TLR8",
"TLR9",
"TLR4")
StackedVlnPlot(obj = merged_data, features = TLR_family_members)

 remeining_members <- c("FCN1",
"PLA2G2D",
"CD24")
StackedVlnPlot(obj = merged_data, features = remeining_members)

##########
merged_data <- RenameIdents(merged_data, annotations)
DimPlot(merged_data, label = FALSE,split.by  = 'group') + NoLegend()
DimPlot(merged_data, label = FALSE,split.by  = 'disease') + NoLegend()

FeaturePlot(object = merged_data, features = c(selected_genes), cols = c("grey", "red"), reduction = "umap",min.cutoff = -1, max.cutoff = 5)
DoHeatmap(merged_data, features=c(selected_genes),group.by = 'group', assay = "SCT")
triplot_noditto(merged_data, c(selected_genes),group.by = 'final_res')  

merged_data@meta.data$final_res <- merged_data@meta.data$integrated_snn_res.0.8  # For convenience

pdf("Ansun_genes_triplot_noditto_final_res.pdf")
triplot_noditto(merged_data, "SIGLECL1",group.by = 'final_res')  
triplot_noditto(merged_data, "SIGLEC1",group.by = 'final_res')  
triplot_noditto(merged_data, "CD22",group.by = 'final_res')  
triplot_noditto(merged_data, "CD33",group.by = 'final_res')  
triplot_noditto(merged_data, "MAG",group.by = 'final_res')  
triplot_noditto(merged_data, "SIGLEC5",group.by = 'final_res')  
triplot_noditto(merged_data, "SIGLEC6",group.by = 'final_res')  
triplot_noditto(merged_data, "SIGLEC7",group.by = 'final_res')  
triplot_noditto(merged_data, "SIGLEC8",group.by = 'final_res')  
triplot_noditto(merged_data, "SIGLEC9",group.by = 'final_res')  
triplot_noditto(merged_data, "SIGLEC10",group.by = 'final_res')  
triplot_noditto(merged_data, "SIGLEC11",group.by = 'final_res')  
triplot_noditto(merged_data, "SIGLEC12",group.by = 'final_res')  
triplot_noditto(merged_data, "SIGLEC14",group.by = 'final_res')  
triplot_noditto(merged_data, "SIGLEC15",group.by = 'final_res')  
triplot_noditto(merged_data, "SIGLEC16",group.by = 'final_res')  
triplot_noditto(merged_data, "MRC1",group.by = 'final_res')  
triplot_noditto(merged_data, "LY75",group.by = 'final_res')  
triplot_noditto(merged_data, "CLEC7A",group.by = 'final_res')  
triplot_noditto(merged_data, "CLEC6A",group.by = 'final_res')  
triplot_noditto(merged_data, "CLEC4E",group.by = 'final_res')  
triplot_noditto(merged_data, "CD209",group.by = 'final_res')  
triplot_noditto(merged_data, "CLEC10A",group.by = 'final_res')  
triplot_noditto(merged_data, "CLEC9A",group.by = 'final_res')  
triplot_noditto(merged_data, "CD207",group.by = 'final_res')  
triplot_noditto(merged_data, "TLR1",group.by = 'final_res')  
triplot_noditto(merged_data, "TLR10",group.by = 'final_res')  
triplot_noditto(merged_data, "TLR2",group.by = 'final_res')  
triplot_noditto(merged_data, "TLR3",group.by = 'final_res')  
triplot_noditto(merged_data, "TLR5",group.by = 'final_res')  
triplot_noditto(merged_data, "TLR6",group.by = 'final_res')  
triplot_noditto(merged_data, "TLR7",group.by = 'final_res')  
triplot_noditto(merged_data, "TLR8",group.by = 'final_res')  
triplot_noditto(merged_data, "TLR9",group.by = 'final_res')  
triplot_noditto(merged_data, "TLR4",group.by = 'final_res')  
triplot_noditto(merged_data, "FCN1",group.by = 'final_res')  
triplot_noditto(merged_data, "PLA2G2D",group.by = 'final_res')  
triplot_noditto(merged_data, "CD24",group.by = 'final_res')  
dev.off()

pdf("Ansun_genes_triplot_noditto_group.pdf")
triplot_noditto(merged_data, "SIGLECL1",group.by = 'group')  
triplot_noditto(merged_data, "SIGLEC1",group.by = 'group')  
triplot_noditto(merged_data, "CD22",group.by = 'group')  
triplot_noditto(merged_data, "CD33",group.by = 'group')  
triplot_noditto(merged_data, "MAG",group.by = 'group')  
triplot_noditto(merged_data, "SIGLEC5",group.by = 'group')  
triplot_noditto(merged_data, "SIGLEC6",group.by = 'group')  
triplot_noditto(merged_data, "SIGLEC7",group.by = 'group')  
triplot_noditto(merged_data, "SIGLEC8",group.by = 'group')  
triplot_noditto(merged_data, "SIGLEC9",group.by = 'group')  
triplot_noditto(merged_data, "SIGLEC10",group.by = 'group')  
triplot_noditto(merged_data, "SIGLEC11",group.by = 'group')  
triplot_noditto(merged_data, "SIGLEC12",group.by = 'group')  
triplot_noditto(merged_data, "SIGLEC14",group.by = 'group')  
triplot_noditto(merged_data, "SIGLEC15",group.by = 'group')  
triplot_noditto(merged_data, "SIGLEC16",group.by = 'group')  
triplot_noditto(merged_data, "MRC1",group.by = 'group')  
triplot_noditto(merged_data, "LY75",group.by = 'group')  
triplot_noditto(merged_data, "CLEC7A",group.by = 'group')  
triplot_noditto(merged_data, "CLEC6A",group.by = 'group')  
triplot_noditto(merged_data, "CLEC4E",group.by = 'group')  
triplot_noditto(merged_data, "CD209",group.by = 'group')  
triplot_noditto(merged_data, "CLEC10A",group.by = 'group')  
triplot_noditto(merged_data, "CLEC9A",group.by = 'group')  
triplot_noditto(merged_data, "CD207",group.by = 'group')  
triplot_noditto(merged_data, "TLR1",group.by = 'group')  
triplot_noditto(merged_data, "TLR10",group.by = 'group')  
triplot_noditto(merged_data, "TLR2",group.by = 'group')  
triplot_noditto(merged_data, "TLR3",group.by = 'group')  
triplot_noditto(merged_data, "TLR5",group.by = 'group')  
triplot_noditto(merged_data, "TLR6",group.by = 'group')  
triplot_noditto(merged_data, "TLR7",group.by = 'group')  
triplot_noditto(merged_data, "TLR8",group.by = 'group')  
triplot_noditto(merged_data, "TLR9",group.by = 'group')  
triplot_noditto(merged_data, "TLR4",group.by = 'group')  
triplot_noditto(merged_data, "FCN1",group.by = 'group')  
triplot_noditto(merged_data, "PLA2G2D",group.by = 'group')  
triplot_noditto(merged_data, "CD24",group.by = 'group')  
dev.off()

pdf("Ansun_genes_triplot_noditto_disease.pdf")
triplot_noditto(merged_data, "SIGLECL1",group.by = 'disease')  
triplot_noditto(merged_data, "SIGLEC1",group.by = 'disease')  
triplot_noditto(merged_data, "CD22",group.by = 'disease')  
triplot_noditto(merged_data, "CD33",group.by = 'disease')  
triplot_noditto(merged_data, "MAG",group.by = 'disease')  
triplot_noditto(merged_data, "SIGLEC5",group.by = 'disease')  
triplot_noditto(merged_data, "SIGLEC6",group.by = 'disease')  
triplot_noditto(merged_data, "SIGLEC7",group.by = 'disease')  
triplot_noditto(merged_data, "SIGLEC8",group.by = 'disease')  
triplot_noditto(merged_data, "SIGLEC9",group.by = 'disease')  
triplot_noditto(merged_data, "SIGLEC10",group.by = 'disease')  
triplot_noditto(merged_data, "SIGLEC11",group.by = 'disease')  
triplot_noditto(merged_data, "SIGLEC12",group.by = 'disease')  
triplot_noditto(merged_data, "SIGLEC14",group.by = 'disease')  
triplot_noditto(merged_data, "SIGLEC15",group.by = 'disease')  
triplot_noditto(merged_data, "SIGLEC16",group.by = 'disease')  
triplot_noditto(merged_data, "MRC1",group.by = 'disease')  
triplot_noditto(merged_data, "LY75",group.by = 'disease')  
triplot_noditto(merged_data, "CLEC7A",group.by = 'disease')  
triplot_noditto(merged_data, "CLEC6A",group.by = 'disease')  
triplot_noditto(merged_data, "CLEC4E",group.by = 'disease')  
triplot_noditto(merged_data, "CD209",group.by = 'disease')  
triplot_noditto(merged_data, "CLEC10A",group.by = 'disease')  
triplot_noditto(merged_data, "CLEC9A",group.by = 'disease')  
triplot_noditto(merged_data, "CD207",group.by = 'disease')  
triplot_noditto(merged_data, "TLR1",group.by = 'disease')  
triplot_noditto(merged_data, "TLR10",group.by = 'disease')  
triplot_noditto(merged_data, "TLR2",group.by = 'disease')  
triplot_noditto(merged_data, "TLR3",group.by = 'disease')  
triplot_noditto(merged_data, "TLR5",group.by = 'disease')  
triplot_noditto(merged_data, "TLR6",group.by = 'disease')  
triplot_noditto(merged_data, "TLR7",group.by = 'disease')  
triplot_noditto(merged_data, "TLR8",group.by = 'disease')  
triplot_noditto(merged_data, "TLR9",group.by = 'disease')  
triplot_noditto(merged_data, "TLR4",group.by = 'disease')  
triplot_noditto(merged_data, "FCN1",group.by = 'disease')  
triplot_noditto(merged_data, "PLA2G2D",group.by = 'disease')  
triplot_noditto(merged_data, "CD24",group.by = 'disease')  
dev.off()

####### plot clusters by patients #########
ggplot(merged_data@meta.data) + geom_bar(aes(x=group , fill=seurat_clusters) , stat="count" , position="fill") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + facet_grid(~seurat_clusters) +     scale_y_continuous(labels=percent)

+ facet_wrap(vars(group))
ggplot(merged_data@meta.data) + geom_bar(aes(x=group , fill=seurat_clusters) , stat="count" , position="fill")
ggplot(merged_data@meta.data) + geom_bar(aes(x=group , fill=seurat_clusters) , stat="count" , position="fill") + theme(axis.text.x=element_text(angle=45, hjust = 1)) + facet_wrap(vars(group))

  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  merged_data@meta.data$dose<-as.factor(merged_data@meta.data$seurat_clusters)
  ggplot(merged_data@meta.data) +  geom_bar(aes(x=group, fill=seurat_clusters), stat="count", position=position_dodge())+
    scale_fill_brewer(palette="Paired")+
    theme_minimal()

  ####### plot clusters by patients #########
  # merged_data@meta.data$sample_new
  
  ggplot(merged_data@meta.data) + geom_bar(aes(x=seurat_clusters , fill=sample_new) , stat="count" , position="fill") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
  DimPlot(merged_data, label = FALSE,split.by  = 'sample_new') + NoLegend()
  DimPlot(merged_data, label = FALSE,group.by   = 'sample_new') + NoLegend()
  saveRDS(merged_data, file='merged_data.RDS')
  
  
  + facet_grid(~seurat_clusters) +     scale_y_continuous(labels=percent)
  seurat_clusters
  + facet_wrap(vars(group))
  ggplot(merged_data@meta.data) + geom_bar(aes(x=group , fill=seurat_clusters) , stat="count" , position="fill")
  ggplot(merged_data@meta.data) + geom_bar(aes(x=group , fill=seurat_clusters) , stat="count" , position="fill") + theme(axis.text.x=element_text(angle=45, hjust = 1)) + facet_wrap(vars(group))
  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  merged_data@meta.data$dose<-as.factor(merged_data@meta.data$seurat_clusters)
  ggplot(merged_data@meta.data) +  geom_bar(aes(x=group, fill=seurat_clusters), stat="count", position=position_dodge())+
    scale_fill_brewer(palette="Paired")+
    theme_minimal()
  
  
  ########
  library(ggplot2)
  library(scales)
  #
  # Displays bar heights as percents with percentages above bars
  #
  ggplot(merged_data@meta.data) + aes(x= test2,  group=test1)) +     geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count") + geom_text(aes( label = scales::percent(..prop..),y= ..prop.. ), stat= "count", vjust = -.5) +
    labs(y = "Percent", fill="test2") +
    facet_grid(~test1) +
    scale_y_continuous(labels=percent)
  #
  # Displays bar heights as percents with counts above bars
  #
  ggplot(test, aes(x= test2,  group=test1)) + 
    geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count") +
    geom_text(aes(label = ..count.., y= ..prop..), stat= "count", vjust = -.5) +
    labs(y = "Percent", fill="test2") +
    facet_grid(~test1) +
    scale_y_continuous(labels=percent)
#### DKK1 ######
FeaturePlot(merged_data, features = "DKK1", min.cutoff = -1, max.cutoff = 5,  assay = "SCT") 
VlnPlot(object = merged_data, features ="DKK1", assay = "SCT", pt.size = 0.1, group.by = 'final_res') #log = T 
VlnPlot(object = merged_data, features ="DKK1", assay = "SCT", pt.size = 0.1, group.by = 'group') #log = T 
VlnPlot(object = merged_data, features ="DKK1", assay = "SCT", pt.size = 0.1, group.by = 'disease') #log = T 
DotPlot(merged_data, features = "DKK1", assay = "SCT") + RotatedAxis()

####################### extra #########
# Immune Cells
triplot_noditto(merged_data, "PTPRC",group.by = 'final_res')  # CD45

# T Cell lineage 
triplot_noditto(sobj, "CD3E",group.by = 'final_res')  # CD3 Epsilon
d

# Cytotoxic T Cells
triplot_noditto(sobj, "CD8A",group.by = 'final_res')

# CD4 T Cells (And myeloid cells)
triplot_noditto(sobj, "CD4",group.by = 'final_res')
triplot_noditto(sobj, "IL7R",group.by = 'final_res')
# Memory CD4
triplot_noditto(sobj, "CCR7",group.by = 'final_res')
# Naive CD4
triplot_noditto(sobj, "S100A4",group.by = 'final_res')

#Tregs
triplot_noditto(sobj, "FOXP3",group.by = 'final_res')


# B Cells (CD19, CD20)
triplot_noditto(sobj, "CD19",group.by = 'final_res')
triplot_noditto(sobj, "MS4A1",group.by = 'final_res')  # CD20
triplot_noditto(sobj, "CD79A",group.by = 'final_res')

# Plasma Cells (These can be CD45+/-, CD19/20 +/- too)
triplot_noditto(sobj, "MZB1",group.by = 'final_res') # This also stains pDCs (useful for separating cDCs from pDCs)

# NK Cells
triplot_noditto(sobj, "GNLY",group.by = 'final_res')
triplot_noditto(sobj, "NKG7",group.by = 'final_res')
triplot_noditto(sobj, "KLRB1",group.by = 'final_res')  # NK1.1
triplot_noditto(sobj, "NCAM1",group.by = 'final_res') # CD56

# Myeloid Cells
triplot_noditto(sobj, "CD68",group.by = 'final_res')  # Marker
triplot_noditto(sobj, "CST3",group.by = 'final_res')

# Monocytes/Neutrophils
triplot_noditto(sobj, "S100A8",group.by = 'final_res')
triplot_noditto(sobj, "S100A9",group.by = 'final_res')

# Monocytes
triplot_noditto(sobj, "MS4A7",group.by = 'final_res')
triplot_noditto(sobj, "VCAN",group.by = 'final_res') # classical monocytes
triplot_noditto(sobj, "VMO1",group.by = 'final_res') # Non Classical Monocyte
triplot_noditto(sobj, "LYPD2",group.by = 'final_res') # Non Classical Monocyte


# Monocytes and NK Cells
triplot_noditto(sobj, "FCGR3A",group.by = 'final_res') # CD16

# Neutrophil
triplot_noditto(sobj, "NAMPT",group.by = 'final_res')
triplot_noditto(sobj, "FCGR3B",group.by = 'final_res')

# Macrophages
triplot_noditto(sobj, "APOE",group.by = 'final_res')
triplot_noditto(sobj, "C1QA",group.by = 'final_res')

#DCs
triplot_noditto(sobj, "FCER1A",group.by = 'final_res')
triplot_noditto(sobj, "FLT3",group.by = 'final_res')  # separates Macs from cDCs
triplot_noditto(sobj, "IL3RA",group.by = 'final_res') # pDCs, specifically
triplot_noditto(sobj, "CD1E",group.by = 'final_res') # myeloid DC
triplot_noditto(sobj, "CD1D",group.by = 'final_res') # myeloid DC
triplot_noditto(sobj, "PKIB",group.by = 'final_res')

# RBCs
triplot_noditto(sobj, "HBB",group.by = 'final_res')
triplot_noditto(sobj, "HBA1",group.by = 'final_res')
triplot_noditto(sobj, "HBA2",group.by = 'final_res')

#Platelets
triplot_noditto(sobj, "PPBP",group.by = 'final_res')
triplot_noditto(sobj, "PF4",group.by = 'final_res')

# MHCII
print(triplot_noditto(sobj, 'HLA-DRA',group.by = 'final_res'))

DoHeatmap(sobj, features=c("PTPRC", "CD3E", "CD8A", "CD4", "IL7R", 
                           "CCR7", "S100A4", "FOXP3", "CD19", "MS4A1", 
                           "CD79A", "MZB1", "GNLY", "NKG7", "KLRB1", "NCAM1",
                           "CD68", "CST3", "S100A8", "S100A9", "MS4A7", 
                           "VCAN", "FCGR3A", "NAMPT", "APOE", "FCER1A", 
                           "HBB", "HBA1", "HBA2", "PPBP", "PF4"))

annotations <- c(
  "Macrophages",   # 0
  "Macrophages",   # 1
  "Monocytes",         # 2
  "Monocytes - non classicals",   # 3
  "CD8 T",    # 4
  "CD4 T",               # 5
  "Dentritic Cells",      # 6
  "?",     # 7
  "Dentritic Cells",  # 8
  "CD4 T",   # 9
  "Epithelium",  # 10
  "Macrophage",          # 11
  "Epithelium",   # 12
  "Epithelium"    # 13
)
names(annotations) <- levels(sobj)
sobj <- RenameIdents(sobj, annotations)
```

Here's the finished product. Congraulations, you just analyzed your first Seurat object!

```{r final UMAP, message=FALSE}
dittoDimPlot('ident', sobj, reduction.use="umap")
```

P.S: Don't forget to save your processed object!
``` {r save_processed}
save(sobj, file=file.path(box_dir, 
                          "data",
                          "C148",
                          "processed",
                          "sobj_SCTransformed_processed.Robj"))
```
########
merged_data <- RunPCA(merged_data, 
               verbose = FALSE, 
               seed.use = 21212)


merged_data <- RunUMAP(merged_data,
                dims = 1:13,  # Num PCs to use
                # Default. Controls how UMAP balances local (low) 
                # versus global (large) structure in the data
                n.neighbors = 30,  
                # Default. Controls the size of the clusters. 
                # Should be smaller than spread
                min.dist = 0.3,   
                # Default. Controls the inter-cluster distances to some extent. 
                # Should be larger than min_dist
                spread = 1,  
                # Default. Can be used with b instead of using min.dist/spread
                a = NULL,  
                # Default. Can be used with a instead of using min.dist/spread
                b = NULL,  
                verbose = FALSE,
                seed.use = 21212)

merged_data <- FindNeighbors(merged_data,
                      dims = 1:30,  # Num PCs to use
                      k.param = 20,  # k for the knn algorithm
                      verbose = FALSE)
# Use the neighborhood graph to cluster the data
merged_data <- FindClusters(merged_data, 
                     verbose = TRUE,
                     # Use the Louvain algorithm. Leiden = 4
                     algorithm = 1,
                     # Clustering resolution. 
                     # Higher values (usually above 1 for Louvain) 
                     # lead to more clusters.
                     resolution = 0.8, 
                     random.seed = 21212)

saveRDS(merged_data, file='merged_data.RDS')

#to call the different objects use:
FeaturePlot(object = merged_data, reduction ="umap",features = "SIGLEC1")
