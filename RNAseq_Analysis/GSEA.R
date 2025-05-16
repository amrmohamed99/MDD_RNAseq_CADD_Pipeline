############### Gene set enrichment analysis using Clusterprofiler#
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(cowplot)
library(GOplot)
library(readxl)
library(dplyr)
library(GOSemSim)
library(pathview)
library(writexl)

# reading in data from deseq2
setwd("E:/grad")
df <-  read_excel("all_ordered_genes_DEA.xlsx")
row.names(df) = df$genes
eg = bitr(df$genes , fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
eg <- eg[!duplicated(eg$ENSEMBL), ]

row.names(eg) = eg$ENSEMBL

common_genes <- intersect(rownames(eg), rownames(df))
df2 = df[common_genes, ]
df2 = cbind(df2, eg$ENTREZID)

rownames(df) = df$genes
colnames(df)[7] ='ENTREZID'

df2 <- df2[!duplicated(df2$`eg$ENTREZID`), ]
rownames(df2) = df2$`eg$ENTREZID`
colnames(df2)[8] ='ENTREZID'

# ####### convert between gene annotation types
# eg = bitr(df$ENTREZID, fromType="ENTREZID", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
# 
# # Inner join based on symbols
# merged_df <- inner_join(df, eg, by = c("ENTREZID" = "ENTREZID"))
# 
# # Check if the SYMOL column contains any NAs
# if (any(is.na(merged_df$SYMBOL))) {
#   print("The column contains NAs.")
# } else {
#   print("The column does not contain NAs.")
# }
# 
# merged_df = merged_df[,-7]
# rownames(merged_df) = merged_df$SYMBOL

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$ENTREZID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

#######################################################
# we want the log2 fold change 
original_gene_list2 <- df2$log2FoldChange

# name the vector
names(original_gene_list2) <- df2$ENTREZID

# omit any NA values 
gene_list2<-na.omit(original_gene_list2)

# sort the list in decreasing order (required for clusterProfiler)
gene_list2 = sort(gene_list2, decreasing = TRUE)

organism = "org.Hs.eg.db"
keytypes(org.Hs.eg.db)

edo <- enrichDGN(gene_list2)

###################################################
############## Gene Set Enrichment Analysis #####
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

res_GO = as.data.frame(gse)

## Dot-plot:
require(DOSE)
A= dotplot(gse, showCategory=5, split=".sign", title = "Enriched Biological Processes") + facet_grid(.~.sign) 

## Encrichment Map:
kk = pairwise_termsim (gse) ## create similarity matrix
emapplot(kk, showCategory = 15)

## Category Netplot:
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 1)

## ridgeplot:
ridgeplot(gse) + labs(x = "enrichment distribution")

## GSEA Plot:
gseaplot(gse, by = "all", title = gse$Description[42], geneSetID = 1)


## PubMed trend of enriched terms
terms <- gse$Description["neurotransmitter secretion"]
pmcplot(terms, 2010:2024, proportion=FALSE)

#################################################
## KEGG Gene Set Enrichment Analysis (Biological Pathways)
kk2 <- gseKEGG(geneList     = gene_list2,
               organism     = "hsa",
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

kegg_pathways = as.data.frame(kk2)


## Dotplot:
B= dotplot(kk2, showCategory = 10, title = "Enriched KEGG Pathways" , split=".sign") + facet_grid(.~.sign)

## Encrichment Map:
k = pairwise_termsim (kk2) ## create similarity matrix
emapplot(k, showCategory = 50)

## Category Netplot:
cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)

## ridgeplot:
ridgeplot(kk2) + labs(x = "KEGG-Enrichment Distribution")

## GSEA Plot:
gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)


# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=gene_list, pathway.id="hsa04330", species = "hsa")

# Produce a different plot (PDF) (not displayed here)
dme <- pathview(gene.data=gene_list, pathway.id="hsa04330", species = "hsa", kegg.native = F)

knitr::include_graphics("hsa04330.pathview.png")


browseKEGG(kk2, 'hsa04330')

write_xlsx(kegg_pathways, "kegg_pathways.xlsx", col_names = TRUE,
           format_headers = TRUE)


C= plot_grid(A, B, nrow=2)

tiff("GSEA.tiff", 
     width = 15, height =18, units = "in",  # Width and height in inches
     res = 600, compression = "lzw")  # High resolution (600 DPI)
C
dev.off()
#####################################################
# Split gene IDs in each row and unlist them
all_genes <- unlist(strsplit(res_GO$core_enrichment, "/"))

# Count the occurrence of each gene across all enriched pathways
gene_counts <- table(all_genes)

# Sort genes by frequency
sorted_genes <- sort(gene_counts, decreasing = TRUE)

# Print top contributing genes
top_contributing_genes <- names(sorted_genes)[1:3]  # Select top 10 genes as an example
print(top_contributing_genes)

eg = bitr(top_contributing_genes, fromType="ENTREZID", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")











