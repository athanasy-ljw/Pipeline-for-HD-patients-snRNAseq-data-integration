# Pipeline for HD patients snRNAseq data integration

```R
##Configure multithreading
library(future)
plan("multisession", workers = 15)
options(future.globals.maxSize = 100000 * 1024^5)


##Loading R packages
library(KernSmooth)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(harmony)
library(cowplot)
library(DoubletFinder)
library(clustree)


##Read source data

GSE152058 <- readRDS(file = "PATH/GSE152058.rds")

new_counts <- read.table(file="PATH/GSE180928_cell_counts.csv",sep = ",", row.names = 1, header = T)
GSE180928 <- CreateSeuratObject(counts = new_counts, min.cells = 3, project = "GSE180928")
meta_data <- read.table("PATH/GSE180928_metadata.csv",header=T,row.names=1,sep=",")
GSE180928 <- AddMetaData(GSE180928, metadata = meta_data)


GSE281069 <- readRDS(file = "PATH/GSE281069.rds")


##Merge objects and construct a Seurat object
HD_patients_merge_data <- merge(GSE152058, y = c(GSE281069,OLG_data), add.cell.ids = c("GSE152058", "GSE281069","GSE180928"), project = "HD_patients_data")


##Calculate the proportion of mitochondrial gene expression
HD_patients_merge_data[["percent.mt"]] <- PercentageFeatureSet(HD_patients_merge_data, pattern = "^MT-")


##Calculate the proportion of mitochondrial gene expression
HD_patients_merge_data <- subset(HD_patients_merge_data, subset = percent.mt < 20 & nFeature_RNA >= 200 & nFeature_RNA  < 10000)


##NormalizeData
HD_patients_merge_data <- NormalizeData(HD_patients_merge_data, normalization.method = "LogNormalize", scale.factor = 10000)


##FindVariableFeatures
HD_patients_merge_data <- FindVariableFeatures(HD_patients_merge_data, selection.method = "vst", nfeatures = 2000)


##ScaleData
HD_patients_merge_data <- ScaleData(HD_patients_merge_data, features = rownames(HD_patients_merge_data))


##RunPCA
HD_patients_merge_data <- RunPCA(HD_patients_merge_data, features = VariableFeatures(object = HD_patients_merge_data))

##Use Harmony to integrate data.
HD_patients_merge_data <- RunHarmony(HD_patients_merge_data, "Tissue_ID", project.dim = F)

##ElbowPlot
ElbowPlot(HD_patients_merge_data,ndims = 50)

##FindNeighbors
HD_patients_merge_data <- FindNeighbors(HD_patients_merge_data, dims = 1:20)



##run clustree
HD_patients_merge_data <- FindClusters(
  object = HD_patients_merge_data,
  resolution = c(seq(.1,2,.2))
)

pdf("clustree.pdf",width=15,height=25)
clustree(HD_patients_merge_data@meta.data, prefix = "RNA_snn_res.")
dev.off()


##Analyzing at a resolution of 0.9.
HD_patients_merge_data <- FindClusters(HD_patients_merge_data, resolution = 0.9)
head(HD_patients_merge_data@meta.data)



##RunTSNE
HD_patients_merge_data <- RunTSNE(HD_patients_merge_data, reduction = "harmony", dims = 1:20)

pdf("HD_patients_merge_data.tsne.pdf",width =12, height =10)
DimPlot(HD_patients_merge_data, reduction = "tsne", label = TRUE)
dev.off()


##FindAllMarkers
Idents(share_HD_patients) <- "seurat_clusters"  # 根据实际列名调整
cluster_markers <- FindAllMarkers(share_HD_patients, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)




##Cell type annotation transfer
prefix_vec <- gsub("^(GSE\\d+)_.*", "\\1", colnames(share_HD_patients))
stripped_bc <- sub("^GSE\\d+_", "", colnames(share_HD_patients))

share_HD_patients$source_annotation <- NA_character_

idx1 <- which(prefix_vec == "GSE152058")
share_HD_patients$source_annotation[idx1] <- GSE152058@meta.data[stripped_bc[idx1], "CellType"]

idx2 <- which(prefix_vec == "GSE281069")
share_HD_patients$source_annotation[idx2] <- GSE281069@meta.data[stripped_bc[idx2], "clusters_named"]

idx3 <- which(prefix_vec == "GSE180928")
share_HD_patients$source_annotation[idx3] <- OLG_data@meta.data[stripped_bc[idx3], "Lineage"]


##Expression of classic cellular marker genes

select_genes_Microglia <- c("RUNX1","ARHGAP24","CALCR","P2RY12","LRMDA","DOCK8","DOCK2","C1QB","ITGAM","CX3CR1","C1QC","APOE","FCER1G","TYROBP","C3","CD14","A2M","CD84","CSF1R","DHRS9","ENTPD1","IL13RA1","IRF8","ALOX5AP","RGS1","RGS10","SRGN","ZFP36L1","C1QA","HCK","LAPTM5","LTC4S","SYK")
pdf("str_Microglia.marker.dotplot.pdf", width=25,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_Microglia) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Neurons <- c("SYT1","ATP1B1","MAP1B","SNAP25","UCHL1","RTN1","CALM1","STMN2","VSNL1","THY1","AHI1","PCSK1N","RBFOX3")
pdf("str_Neurons.marker.dotplot.pdf", width=20,height=10)
DotPlot(object = HD_patients_merge_data, features = select_genes_Neurons) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Excitatory_neurons <- c("LDB2","TLE4","CAMK2A","SATB2","TBR1","CBLN2")
pdf("str_Excitatory_neurons.marker.dotplot.pdf", width=18,height=10)
DotPlot(object = HD_patients_merge_data, features = select_genes_Excitatory_neurons) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Inhibitory_neurons <- c("GAD1","GAD2","SLC6A1","LHFPL3")
pdf("str_Inhibitory_neurons.marker.dotplot.pdf", width=10,height=8)
DotPlot(object = HD_patients_merge_data, features = select_genes_Inhibitory_neurons) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_NSC <- c("DLX1")
pdf("str_NSC.marker.dotplot.pdf", width=10,height=10)
DotPlot(object = HD_patients_merge_data, features = select_genes_NSC) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Neuroblast <- c("TBR1","NEUROD2")
pdf("str_Neuroblast.marker.dotplot.pdf", width=10,height=10)
DotPlot(object = HD_patients_merge_data, features = select_genes_Neuroblast) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Astrocyte <- c("GFAP","SLC1A2","SLC1A3","AQP4","CLU","APOE","ATP1B2","NTM","CPE","GJA1","FGFR3","SFXN5","HIF3A","GRAMD1C","CD38","CLDN10")
pdf("str_Astrocyte.marker.dotplot.pdf", width=18,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_Astrocyte) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Oligodendrocyte <- c("MBP","ERBIN","CNP","MAG","PLP1","ST18","CRYAB","TTLL7","CLDN11","APLP1","SLC44A1","QDPR","TMEM144","TF","LAMP2","SCD","DOCK5","CNTN2","HAPLN2","SLC12A2","PPP1R14A","MOG","GJB1")
pdf("str_Oligodendrocyte.marker.dotplot.pdf", width=20,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_Oligodendrocyte) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_OPC <- c("VCAN","LHFPL3","NOVA1","BCAN","PTPRZ1","TNR","PDGFRA","NG2","CSPG4")
pdf("str_OPC.marker.dotplot.pdf", width=15,height=12)
DotPlot(object = HD_patients_merge_data, features = select_genes_OPC) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_NFO <- c("FYN")
pdf("str_NFO.marker.dotplot.pdf", width=10,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_NFO) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Endothelial <- c("COLEC12","SLC16A9","CLDN5","EBF1","ID3","IGFBP7","ITM2A","DUSP1","TMSB10","EGFL7","SLC2A1","FLT1")
pdf("str_Endothelial.marker.dotplot.pdf", width=18,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_Endothelial) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Myeloid_cell <- c("MNDA","CXCR4")
pdf("str_Myeloid_cell.marker.dotplot.pdf", width=10,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_Myeloid_cell) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Fibroblast <- c("FGFR1","COL1A1","COL3A1")
pdf("str_Fibroblast.marker.dotplot.pdf", width=15,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_Fibroblast) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Vascular <- c("FLT1")
pdf("str_Vascular.marker.dotplot.pdf", width=15,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_Vascular) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Secretory_ependymal <- c("TTR")
pdf("str_Secretory_ependymal.marker.dotplot.pdf", width=10,height=10)
DotPlot(object = HD_patients_merge_data, features = select_genes_Secretory_ependymal) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Ciliated_ependymal <- c("SNTN")
pdf("str_Ciliated_ependymal.marker.dotplot.pdf", width=15,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_Ciliated_ependymal) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Macrophages <- c("MRC1")
pdf("str_Macrophages.marker.dotplot.pdf", width=15,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_Macrophages) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_MSN <- c("PPP1R1B")
pdf("str_MSN.marker.dotplot.pdf", width=5,height=10)
DotPlot(object = HD_patients_merge_data, features = select_genes_MSN) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_dSPN <- c("TAC1","DRD1","PDYN","EBF1","SLC35D3","SFXN1","ASIC4","ARX","NRXN1","FOXP2","DLK1")
pdf("str_dSPN.marker.dotplot.pdf", width=18,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_dSPN) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_iSPN <- c("DRD2","ADORA2A","PENK","GPR52","SP9","CASZ1","OTOF","SIX3","P2RY1","FAM163B","NRXN2","FIG4","UNC5D","GABRA3","NT5E","NELL1","CHRM3","NTS","PTPRM","PLXDC1","NECAB1","GRIK3")
pdf("str_iSPN.marker.dotplot.pdf", width=20,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_iSPN) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Foxp2_Olfm3_expressing_neuron <- c("FOXP2","OLFM3")
pdf("str_Foxp2_Olfm3_expressing_neuron.marker.dotplot.pdf", width=10,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_Foxp2_Olfm3_expressing_neuron) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Cholinergic_interneurons <- c("CHAT","ACHE","SLC5A7","SLC18A3")
pdf("str_Cholinergic_interneurons1.marker.dotplot.pdf", width=15,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_Cholinergic_interneurons) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Serotonergic_neurons <- c("FEV","SLC6A4","TPH2")
pdf("str_Serotonergic_neurons.marker.dotplot.pdf", width=15,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_Serotonergic_neurons) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Dopaminergic_neurons <- c("DBH","FOXA2","KCNJ6","LMX1B","NR4A2","PRKN","SLC6A2","SLC6A3","TH","TOR1A")
pdf("str_Dopaminergic_neurons.marker.dotplot.pdf", width=15,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_Dopaminergic_neurons) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Glutamatergic_neurons <- c("FOLH1B","GOT1","GRIN1","GRIN2B","HDLBP","SLC17A8","GLS","SLC1A1","SLC1A6","SLC17A6","SLC17A7")
pdf("str_Glutamatergic_neurons.marker.dotplot.pdf", width=15,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_Glutamatergic_neurons) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_PENK_interneurons <- c("PENK")
pdf("str_PENK_interneurons.marker.dotplot.pdf", width=15,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_PENK_interneurons) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_GABAergic_neurons <- c("GAD2")
pdf("str_GABAergic_neurons.marker.dotplot.pdf", width=10,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_GABAergic_neurons) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Pvalb_expressing_GABAergic_interneurons <- c("HS3ST2","KIT","PVALB","CNTNAP4")
pdf("str_Pvalb_expressing_GABAergic_interneurons.marker.dotplot.pdf", width=10,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_Pvalb_expressing_GABAergic_interneurons) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Sst_Npy_expressing_GABAergic_interneurons <- c("NOS1","NPY","SST")
pdf("str_Sst_Npy_expressing_GABAergic_interneurons.marker.dotplot.pdf", width=10,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_Sst_Npy_expressing_GABAergic_interneurons) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_VIP_expressing_GABAergic_interneurons <- c("VIP")
pdf("str_VIP_expressing_GABAergic_interneurons.marker.dotplot.pdf", width=10,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_VIP_expressing_GABAergic_interneurons) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Mural_cells <- c("PDGFRB")
pdf("str_Mural_cells.marker.dotplot.pdf", width=10,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_Mural_cells) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Migrating_neuroblasts <- c("DCX")
pdf("str_Migrating_neuroblasts.marker.dotplot.pdf", width=10,height=10)
DotPlot(object = HD_patients_merge_data, features = select_genes_Migrating_neuroblasts) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Cycling_cell <- c("TOP2A")
pdf("str_Cycling_cell.marker.dotplot.pdf", width=10,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_Cycling_cell) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_T_cell <- c("CD3D","CD3E","CD8A","CD4")
pdf("str_T_cell.marker.dotplot.pdf", width=10,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_T_cell) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


select_genes_Granule_precursors <- c("ATOH1")
pdf("Granule_precursors_cell.marker.dotplot.pdf", width=10,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_Granule_precursors) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Granule_cell <- c("NHLH1","RELN","RBFOX3")
pdf("Granule_cell.marker.dotplot.pdf", width=10,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_Granule_cell) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_Purkinje_cell <- c("GAD1","CALB1","PCP2")
pdf("Purkinje_cell.marker.dotplot.pdf", width=10,height=15)
DotPlot(object = HD_patients_merge_data, features = select_genes_Purkinje_cell) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



select_genes_Ex_neurons <- c("SLC17A7","LDB2","TLE4","CAMK2A","SATB2","TBR1","CBLN2")
pdf("Ex_neurons.marker.dotplot.pdf", width=18,height=10)
DotPlot(object = HD_patients_merge_data, features = select_genes_Ex_neurons) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


select_genes_In_neurons <- c("GAD1","GAD2","SLC6A1","LHFPL3")
pdf("In_neurons.marker.dotplot.pdf", width=10,height=8)
DotPlot(object = HD_patients_merge_data, features = select_genes_In_neurons) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_peripheral_immune_cell <- c("CD2","CD8","CD8A")
pdf("peripheral_immune_cell.marker.dotplot.pdf", width=10,height=8)
DotPlot(object = HD_patients_merge_data, features = select_genes_peripheral_immune_cell) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_ependymal_cells <- c("SOX9","CFAP299","TTR","SNTN")
pdf("ependymal_cell.marker.dotplot.pdf", width=10,height=8)
DotPlot(object = HD_patients_merge_data, features = select_genes_ependymal_cells) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_ctx_layer_neuron <- c("LINC00507","RORB","POU3F1","HTR2C","NPFFR2","TWIST1")
pdf("ctx_layer_neuron.marker.dotplot.pdf", width=10,height=8)
DotPlot(object = HD_patients_merge_data, features = select_genes_ctx_layer_neuron) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


select_genes_cingulate_neuron <- c("PCDH17","PHGFC","CDH4","ADAMTS3")
pdf("cingulate_neuron.marker.dotplot.pdf", width=10,height=8)
DotPlot(object = HD_patients_merge_data, features = select_genes_cingulate_neuron) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

select_genes_LAMP5_neuron <- c("LAMP5")
pdf("LAMP5_neuron.marker.dotplot.pdf", width=10,height=8)
DotPlot(object = HD_patients_merge_data, features = select_genes_LAMP5_neuron, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



##Plotting relevant key indicators

DimPlot(object = HD_patients_merge_data, reduction = "pca", pt.size = .1, group.by = "GEO_number", raster=FALSE)

DimPlot(object = HD_patients_merge_data, reduction = "harmony", pt.size = .1, group.by = "GEO_number", raster=FALSE)


Idents(HD_patients_merge_data) <- "seurat_clusters"
DimPlot(HD_patients_merge_data, reduction = "tsne", label = TRUE, raster=FALSE)


Idents(HD_patients_merge_data) <- "Tissue"
DimPlot(HD_patients_merge_data, reduction = "tsne", label = F, raster=FALSE)


Idents(HD_patients_merge_data) <- "Tissue"
DimPlot(HD_patients_merge_data, reduction = "tsne", label = F, raster=FALSE)


Idents(HD_patients_merge_data) <- "Condition"
DimPlot(HD_patients_merge_data, reduction = "tsne", label = F, raster=FALSE)



## tSNE for various brain regions
Idents(human_str) <- "cell_type2"
jpeg("human_str.tsne.jpg", width = 10, height = 8, units = "in", res = 1200, quality = 100)
DimPlot(human_str, reduction = "tsne",label = F, raster=FALSE)
dev.off()

Idents(human_cortex) <- "cell_type2"
jpeg("human_cortex.tsne.jpg", width = 10, height = 8, units = "in", res = 1200, quality = 100)
DimPlot(human_cortex, reduction = "tsne",label = F, raster=FALSE)
dev.off()

Idents(human_Hippocampus) <- "cell_type2"
jpeg("human_Hippocampus.tsne.jpg", width = 10, height = 8, units = "in", res = 1200, quality = 100)
DimPlot(human_Hippocampus, reduction = "tsne",label = F, raster=FALSE)
dev.off()

Idents(human_Cerebellum) <- "cell_type2"
jpeg("human_Cerebellum.tsne.jpg", width = 10, height = 8, units = "in", res = 1200, quality = 100)
DimPlot(human_Cerebellum, reduction = "tsne",label = F, raster=FALSE)
dev.off()

Idents(human_Cingulate) <- "cell_type2"
jpeg("human_Cingulate.tsne.jpg", width = 10, height = 8, units = "in", res = 1200, quality = 100)
DimPlot(human_Cingulate, reduction = "tsne",label = F, raster=FALSE)
dev.off()

Idents(human_Accumbens) <- "cell_type2"
jpeg("human_Accumbens.tsne.jpg", width = 10, height = 8, units = "in", res = 1200, quality = 100)
DimPlot(human_Accumbens, reduction = "tsne",label = F, raster=FALSE)
dev.off()



##Ratio for various brain regions
Ratio <- HD_patients_merge_data@meta.data %>%
  group_by(merge_tissue, cell_type2) %>% # 分组
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))


library(ggalluvial)
library(ggplot2)

ggplot(Ratio, aes(x =merge_tissue, y= relative_freq, fill = cell_type2,
                  stratum=cell_type2, alluvium=cell_type2)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.4, knot.pos=0.5)+
  theme_classic() +
  labs(x='sample',y = 'Ratio')

ggsave("all_tissue_proportion.pdf", width =7, height =8)






```
