suppressPackageStartupMessages({
  library(GEOquery)
  library(oligo)
  library(limma)
  library(Biobase)
  library(splitstackshape)
  library(tidyverse)
  library(arrayQualityMetrics)
  library(affy)
  library(SummarizedExperiment)
})

raw.data <- ReadAffy(celfile.path = "E-GEOD-7621-CEL/")
raw.data
normalized.data <- affy::rma(raw.data)
normalized.expr <- as.data.frame(exprs(normalized.data))

gse <- getGEO("GSE7621", GSEMatrix = TRUE)
feature.data <- gse$GSE7621_series_matrix.txt.gz@featureData@data
feature.data <- feature.data[,c(1,11)]

normalized.expr<-normalized.expr %>%
  rownames_to_column(var = "ID") %>%  
  inner_join(., feature.data, by = 'ID')

counts <- normalized.expr
rownames(counts) <- counts$ID
counts <- counts[, -which(names(counts) == "ID")]
counts <- counts[, -which(names(counts) == "Gene Symbol")]

# Create design matrix for the main analysis
pdata <- pData(gse[[1]])
colnames(pdata)
pdata <- pdata %>%
  rename(Disease.state = characteristics_ch1)%>%
  rename(Sex = characteristics_ch1.1)
# Create factors for Disease state and Sex
pdata$Disease.state <- factor(pdata$`Disease.state`, levels = c("Old Control", "Parkinson's Disease"))
pdata$Sex <- factor(pdata$Sex, levels = c("male", "female"))

# Check the factors
table(pdata$Disease.state)
table(pdata$Sex)

# Create design matrix
design <- model.matrix(~ Disease.state + Sex, data = pdata)

group <- factor(pdata$`disease state`, levels = c("Old_Control", "Parkinsons_Disease"))
colnames(design) <- c("Intercept", "Parkinsons_Disease", "Sex")

# Create design matrix for analysis excluding 'Sex'
design_disease <- model.matrix(~ Disease.state, data = pdata)
colnames(design_disease) <- c("Intercept", "Parkinsons_Disease")

# Fit the linear model
fit <- lmFit(counts, design)
fit <- eBayes(fit)

# Define contrast to compare Parkinson's Disease vs Old Control while adjusting for sex
contrast_matrix <- makeContrasts(Parkinsons_Disease - Sex, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract top results
topTable_sex <- topTable(fit2, adjust.method = "fdr", number = Inf)
sig_genes_sex <- topTable_sex[topTable_sex$adj.P.Val < 0.05,]
sig_genes_sex <- sig_genes_sex[order(sig_genes_sex$adj.P.Val),]
top20_genes_sex <- head(sig_genes_sex, 20)

fit_disease <- lmFit(counts, design_disease)
fit_disease <- eBayes(fit_disease)

# Define contrast to compare Parkinson's Disease vs Old Control
contrast_matrix_disease <- makeContrasts(Parkinsons_Disease, levels = design_disease)
fit2_disease <- contrasts.fit(fit_disease, contrast_matrix_disease)
fit2_disease <- eBayes(fit2_disease)

# Extract top results
topTable_disease <- topTable(fit2_disease, adjust.method = "fdr", number = Inf)
sig_genes_disease <- topTable_disease[topTable_disease$adj.P.Val < 0.05,]
sig_genes_disease <- sig_genes_disease[order(sig_genes_disease$adj.P.Val),]
top20_genes_disease <- head(sig_genes_disease, 20)

# Filter the normalized counts for the top significant genes
norm_counts_sex <- counts[rownames(counts) %in% rownames(top20_genes_sex),]
norm_counts_disease <- counts[rownames(counts) %in% rownames(top20_genes_disease),]

# Merge the results
merged_df <- rbind(norm_counts_disease, norm_counts_sex)
merged_df<-merged_df %>%
  rownames_to_column(var = "ID") %>%  
  inner_join(., feature.data, by = 'ID')%>% 
  filter("Gene Symbol" != "")%>%
  rename(Gene_symbol = "Gene Symbol")%>%
  group_by(Gene_symbol) %>%
  summarize(across(where(is.numeric), mean, na.rm = TRUE))
merged_df <- as.data.frame(merged_df)

rownames(merged_df) <- make.names(merged_df$`Gene_symbol`, unique=TRUE)
merged_df <- merged_df[, -which(names(merged_df) == "Gene_symbol")]

# Convert merged_df to a matrix and calculate correlations
matrix_df <- as.matrix(t(merged_df))
correlation_matrix <- cor(matrix_df)

# Generate correlation plot
library(corrplot)
corrplot(correlation_matrix, order = 'hclust', tl.cex = 0.5)

top20_genes_sex<-top20_genes_sex %>%
  rownames_to_column(var = "ID") %>%  
  inner_join(., feature.data, by = 'ID')

top20_genes_disease<-top20_genes_disease %>%
  rownames_to_column(var = "ID") %>%  
  inner_join(., feature.data, by = 'ID')