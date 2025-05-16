##########################################
#   -   RNA-Seq data analysis            #
#                                        #
#   -   2023-8-2                         #
#   -   copyright : Amr Mohamed          #
#                                        #
##########################################
# R version 4.3.1 (2023-06-16 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# Matrix products: default
# locale:
# [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
# [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# time zone: Africa/Cairo
# tzcode source: internal



############ DESeq2 tutorial https://lashlock.github.io/compbio/R_presentation.html

#### used libraries ######
library(readr)
library(DESeq2)
library(ggfortify)
library(rgl)
library(plot3D)
library(plotly)
library(stats)
library(scatterplot3d)
library(genefilter)
library(matrixStats)
library(ComplexHeatmap)
library(pd.hg.u133.plus.2)
library(readxl)
library(EnhancedVolcano)
library(writexl)


##### load the mRNA-Seq data #####
setwd("E:/grad")

exp = read_csv("ML/merged_dataset.csv")
id = exp[,1]
exp = exp[,-1]
exp = as.matrix(exp)

exp = cbind(id , exp)
colnames(exp)[1] ='id'

######## check duplicates of gene symbol #######
sum(duplicated(exp$id))

###### aggregation step ####
# exp.data.agg= aggregate(exp, list(exp$id),FUN=mean)
# genes = exp.data.agg[1]
# exp.data.agg = exp.data.agg[,c(-1,-2)]
# 
# exp.data.agg = cbind(exp.data.agg , genes)
# rownames(exp.data.agg) = exp.data.agg$Group.1
# exp.data.agg = exp.data.agg[,-16]
# 
# exp.data.agg <- exp.data.agg[complete.cases(exp.data.agg), ] 

rownames(exp) = exp$id
exp = exp [,-1]


########## Data exploration ########
exp.data = apply(exp , 2, as.numeric)

hist(log2(exp.data+1),main = "RNA-Seq Histogram")

boxplot(log2(exp.data+1),main = "RNA-Seq Box Plot",col=seq(1:15))

########### load mRNA sample sheet #######
pheno = read_excel("ML/ML_Conditions.xlsx")

rownames(pheno) = pheno$Title
all(colnames(exp) %in% rownames(pheno))
all(colnames(exp) == rownames(pheno))

############### impute the missing values mean ##############
# Calculate the proportion of zeros in each row
# prop_zeros <- rowSums(exp == 0) / ncol(exp)
# 
# # Identify rows with less than 40% zeros
# rows_to_fill <- which(prop_zeros < 0.4)
# imputed_genes = rows_to_fill
# 
# # Calculate the row means for these rows
# row_means <- rowMeans(exp[rows_to_fill, ], na.rm = TRUE)
# 
# # Replace the zeros in these rows with the row means and remove others rows
# exp[rows_to_fill, ][exp[rows_to_fill, ] == 0] <- row_means
# exp = exp[imputed_genes,]
# 
# exp.count = round(exp)

############### Defferential expression analysis using DESeq2 #########
table(pheno$diagnosis)


dds <- DESeqDataSetFromMatrix(countData = exp,
                              colData = pheno,
                              design = ~ 0 + diagnosis)
dds

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds

# set the factor level
dds$diagnosis <- relevel(dds$diagnosis, ref = "Healthy_control")

dds.run = DESeq(dds)

res=results(dds.run)
res.df = as.data.frame(res)
summary(res)

resultsNames(dds.run)

########### make contrasts (comparisons between all conditins) ########
contrast1 = results(dds.run, contrast = c("diagnosis", "Healthy_control","MDD"))
contrast1 = contrast1[complete.cases(contrast1),]
contrast1 = as.data.frame(contrast1)
contrast1$genes = rownames(contrast1)

res.degs1 = contrast1[contrast1$padj<0.05 & abs(contrast1$log2FoldChange)>log2(2),]
degs.genes1= rownames(res.degs1)
write.table(degs.genes1,file = "final_DEGs.txt",row.names = F,col.names = F,quote = F)

contrast2 = results(dds.run, contrast = c("complex", "ctrlF","MDDM"))
contrast2 = contrast2[complete.cases(contrast2),]
contrast2 = as.data.frame(contrast2)
res.degs2 = contrast2[contrast2$padj<0.05 & abs(contrast2$log2FoldChange)>log2(1),]
degs.genes2= rownames(res.degs2)
write.table(degs.genes2,file = "DEGs2.txt",row.names = F,col.names = F,quote = F)
exp.degs2=exp[degs.genes2,]
contrast2$genes = rownames(contrast2)
write_xlsx(contrast2, "all_genes_contrast2.xlsx", col_names = TRUE,
           format_headers = TRUE)


contrast3 = results(dds.run, contrast = c("condition", "Healthy_controlNNo_SI","MDDNHigh_SI"))
contrast3 = contrast3[complete.cases(contrast3),]
contrast3 = as.data.frame(contrast3)
res.degs3 = contrast3[contrast3$padj<0.05 & abs(contrast3$log2FoldChange)>log2(1),]
degs.genes3= rownames(res.degs3)
write.table(degs.genes3,file = "DEGs3.txt",row.names = F,col.names = F,quote = F)

contrast4 = results(dds.run, contrast = c("condition", "Healthy_controlNNo_SI","MDDNLow_SI"))
contrast4 = contrast4[complete.cases(contrast4),]
contrast4 = as.data.frame(contrast4)
res.degs4 = contrast4[contrast4$padj<0.05 & abs(contrast4$log2FoldChange)>log2(1),]
degs.genes4= rownames(res.degs4)
write.table(degs.genes4,file = "DEGs4.txt",row.names = F,col.names = F,quote = F)

contrast5 = results(dds.run, contrast = c("condition", "Healthy_controlNNo_SI","MDDNNo_SI"))
contrast5 = contrast5[complete.cases(contrast5),]
contrast5 = as.data.frame(contrast5)
res.degs5 = contrast5[contrast4$padj<0.05 & abs(contrast5$log2FoldChange)>log2(1),]
degs.genes5= rownames(res.degs5)
write.table(degs.genes5,file = "DEGs5.txt",row.names = F,col.names = F,quote = F)

############# get the DEGs based on adj pval , LFC ############
res.df = as.data.frame(res.degs1)
res.df=res.df[complete.cases(res.degs1),]
res.df$genes = rownames(res.df)
write_xlsx(res.df, "all_genes_DEA1.xlsx", col_names = TRUE,
           format_headers = TRUE)

res.df.ordered=res.df[order(res.df[,6]) ,]
write_xlsx(res.df.ordered, "all_ordered_genes_DEA.xlsx", col_names = TRUE,
           format_headers = TRUE)

res.degs = res.df[res.df$padj<0.05 & abs(res.df$log2FoldChange)>log2(2),]
res.degs=res.degs[order(res.degs[,6]) ,]
degs.genes= rownames(res.degs)
exp.degs=exp[degs.genes,]
write.table(degs.genes,file = "DEGs.txt",row.names = F,col.names = F,quote = F)

write_xlsx(df, "filename.xlsx")

############### classify DEGs based on up or down regulation ###########
res.degs= as.data.frame(res.degs)
res.degs[,'REGULATION']='DOWN'
res.degs[res.degs$log2FoldChange > 0 ,]$REGULATION ="UP"
res.degs.up = res.degs[res.degs$REGULATION=='UP' ,]
res.degs.down = res.degs[res.degs$REGULATION=='DOWN' ,]

############# do normalization for all exp data to further analysis ########
ntd=normTransform(dds)
exp.norm= assay(ntd)

############### creating a heatmap for the top 100 DEG genes #####
exp.degs2=exp.norm[degs.genes,]
top100_DEGS = row.names(exp.degs2)[1:20]
exp100_DEGS = exp.degs2[top100_DEGS,]
exp21_DEGS = exp100_DEGS[, c(101:110,363:373)]

pheno2 = pheno[c(101:110,363:373),]
column_ha = HeatmapAnnotation(sample.type = pheno2$diagnosis)
A= Heatmap(exp21_DEGS,name = 'Exp', row_names_gp= gpar(fontsize=10) , column_names_gp = gpar(fontsize=5)
        , top_annotation = column_ha)


############### 2D PCA ############3
expression_t= t(exp100_DEGS)
expression.pca = prcomp(expression_t , center = TRUE , scale. = TRUE)
summary(expression.pca)
B= autoplot(expression.pca, data = pheno, colour = 'diagnosis')


C= EnhancedVolcano(contrast1,
                lab = contrast1$genes,
                x = 'log2FoldChange',
                y = 'padj',
                xlim =  c(-5, 5),
                ylim = c(0.5, 35),
                pCutoff = 0.05,
                pointSize = 2,
                FCcutoff = 1,
                pCutoffCol='padj',
                title = "(fold change cutoff = 1, padj cutoff = 0.05)"
)

# 1. boxplot
D= boxplot(exp21_DEGS,main = "processed genes Box Plot",col=seq(1:15))

# 2. histogram
cols=colnames(exp100_DEGS)
hist(exp100_DEGS,main = "processed genes Histogram")


# Save as TIFF with custom size and resolution
tiff("Box.tiff", 
     width = 9, height = 7, units = "in",  # Width and height in inches
     res = 600, compression = "lzw")  # High resolution (600 DPI)
D
dev.off()

ggsave("Box.png", plot = D, width = 10, height = 8, dpi = 600, device = "tiff")


