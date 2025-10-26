library(DESeq2)
library(tibble)
library(VennDiagram)
library(ggplot2)
library(ggrepel)
library(ggalt)

bulk_sample_meta_path <- "./stadtmauer_bulk_sample_metadata.csv"
bulk_counts_path      <- "./bulkseq/stadtmauer_bulk_counts.csv"

bulk_sample_meta <- read.csv(bulk_sample_meta_path, row.names = 1)
bulk_counts      <- as.matrix(read.csv(bulk_counts_path, row.names = 1))

pseudobulk_sample_meta_path <- "./cite_fib_pseudobulk_sample_metadata.csv"
pseudobulk_counts_path      <- "./cite_fib_pseudobulk_counts.csv"
pseudobulk_gene_meta_path   <- "./cite_fib_pseudobulk_gene_metadata.csv"

pseudobulk_sample_meta <- read.csv(pseudobulk_sample_meta_path, row.names = 1)
pseudobulk_counts      <- as.matrix(read.csv(pseudobulk_counts_path, row.names = 1))

common_genes <- intersect(colnames(bulk_counts), colnames(pseudobulk_counts))
bulk_counts_filtered <- bulk_counts[,common_genes]
pseudobulk_counts_filtered <- pseudobulk_counts[,common_genes]


# create DESeqDataSet object for bulk RNA-seq data
dds_bulk <- DESeqDataSetFromMatrix(
  countData = t(bulk_counts),
  colData = bulk_sample_meta,
  design = ~ treatment  
)


# run the DESeq2 pipeline
dds_bulk <- DESeq(dds_bulk)

# compare PGE2 to non PGE2
res_bulk <- results(dds_bulk, contrast = c("treatment","PGE2MPA","MPA"))

summary(res_bulk)

res_bulk_df1 = as.data.frame(res_bulk) %>% rownames_to_column('gene')

res_bulk_df1 %>% write.csv('PGE2_MPA_results.csv')


# compare cAMP to non cAMP
res_bulk <- results(dds_bulk, contrast = c("treatment","CAMPMPA","MPA"))

summary(res_bulk)

res_bulk_df2 = as.data.frame(res_bulk) %>% rownames_to_column('gene')

res_bulk_df2 %>% write.csv('cAMP_MPA_results.csv')


# process single cell pseudobulk data
pseudobulk_sample_meta_db = pseudobulk_sample_meta %>% dplyr::filter(tissue == 'basalis')
pseudobulk_sample_meta_db = pseudobulk_sample_meta_db %>% dplyr::filter(donor_id != 'sc01')
pseudobulk_counts_db <- pseudobulk_counts[rownames(pseudobulk_sample_meta_db), ]

# create DESeqDataSet object for pseudobulk data
dds_pseudobulk <- DESeqDataSetFromMatrix(
  countData = t(pseudobulk_counts_db),
  colData = pseudobulk_sample_meta_db,
  design = ~ donor_id+cell_type  
)


# run the DESeq2 analysis
dds_pseudobulk <- DESeq(dds_pseudobulk)

# compare decFib to hpFig (former genFIB)
res_pseudobulk <- results(dds_pseudobulk, contrast = c("cell_type","decFIB","genFIB"))

summary(res_pseudobulk)

res_pseudobulk_df3 = as.data.frame(res_pseudobulk) %>% rownames_to_column('gene')

res_pseudobulk_df3 %>% write.csv('decFIB_hpFIB_results.csv')


# filter for upregulated genes in each DESeq2 result dataframe.
# adjust the criteria (padj < 0.05 and log2FoldChange > 0) as appropriate.

up_pseudobulk <- res_pseudobulk_df3$gene[
  !is.na(res_pseudobulk_df3$padj) & 
    res_pseudobulk_df3$padj < 0.05 & 
    res_pseudobulk_df3$log2FoldChange > 0 
]

up_bulk_df2 <- res_bulk_df2$gene[
  !is.na(res_bulk_df2$padj) & 
    res_bulk_df2$padj < 0.05 & 
    res_bulk_df2$log2FoldChange > 0
]

up_bulk_df <- res_bulk_df$gene[
  !is.na(res_bulk_df$padj) & 
    res_bulk_df$padj < 0.05 & 
    res_bulk_df$log2FoldChange > 0
]

# print number of upregulated genes in each dataset
cat("Upregulated genes in pseudobulk:", length(up_pseudobulk), "\n")
cat("Upregulated genes in bulk_df2:", length(up_bulk_df2), "\n")
cat("Upregulated genes in bulk_df:", length(up_bulk_df), "\n")

# calculate overlaps between all three datasets
overlap_all <- Reduce(intersect, list(up_pseudobulk, up_bulk_df2, up_bulk_df))
cat("Number of genes upregulated in all three:", length(overlap_all), "\n")

# calculate pairwise overlaps
overlap_pseudobulk_bulk_camp_db <- intersect(up_pseudobulk, up_bulk_df2)
overlap_pseudobulk_bulk_pge_db  <- intersect(up_pseudobulk, up_bulk_df)
overlap_bulk_df2_bulk_df_camp_pge    <- intersect(up_bulk_df2, up_bulk_df)

cat("Overlap (decFIB & cAMP):", length(overlap_pseudobulk_bulk_camp_db), "\n")
cat("Overlap (decFIB & PGE2):", length(overlap_pseudobulk_bulk_pge_db), "\n")
cat("Overlap (PGE2 & cAMP):", length(overlap_bulk_df2_bulk_df_camp_pge), "\n")

# calculate genes that are unique to each dataset
unique_pseudobulk <- setdiff(up_pseudobulk, union(up_bulk_df2, up_bulk_df))
unique_bulk_df2   <- setdiff(up_bulk_df2, union(up_pseudobulk, up_bulk_df))
unique_bulk_df    <- setdiff(up_bulk_df, union(up_pseudobulk, up_bulk_df2))

cat("Unique to decFIB:", length(unique_pseudobulk), "\n")
cat("Unique to bulk_cAMP:", length(unique_bulk_df2), "\n")
cat("Unique to bulk_PGE2:", length(unique_bulk_df), "\n")



gene_lists <- list(
  decFIB = up_pseudobulk,
  cAMP_Prog   = up_bulk_df2,
  PGE2_Prog    = up_bulk_df
)

# open pdf
pdf("./venn_diagram.pdf", width = 6, height = 6)

# generate the Venn diagram
venn.plot <- venn.diagram(
  x = gene_lists,
  filename = NULL,
  fill = c("red", "blue", "darkgreen"),
  alpha = 0.4,
  cex = 1.5,
  cat.cex = 1.5,
  main = "Upregulated Genes Overlap"
)

# draw the diagram onto the PDF
grid::grid.draw(venn.plot)

# close the PDF device
dev.off()



# set thresholds
padj_cutoff <- 0.05
log2fc_cutoff <- 1

res_df = res_pseudobulk_df3

# create capped versions of variables
res_df$capped_log2FoldChange <- with(res_df, ifelse(log2FoldChange > 10, 10,
                                                    ifelse(log2FoldChange < -10, -10, log2FoldChange)))
res_df$padj_fixed <- ifelse(res_df$padj == 0, 1e-300, res_df$padj)
res_df$logp <- -log10(res_df$padj_fixed)
res_df$capped_logp <- ifelse(res_df$logp > 100, 100, res_df$logp)

# define significance status
p_adj_cutoff <- 0.05
log2fc_cutoff <- 1
res_df$sig <- "Not Significant"
res_df$sig[!is.na(res_df$padj) & res_df$padj < p_adj_cutoff & res_df$log2FoldChange > log2fc_cutoff] <- "Upregulated"
res_df$sig[!is.na(res_df$padj) & res_df$padj < p_adj_cutoff & res_df$log2FoldChange < -log2fc_cutoff] <- "Downregulated"


# define the genes to label
selected_genes <- c("PRL",'ARG1','IGFBP1','CNR1','VEGFA','CD82','IL1R2','IL1RL1','WNT1','WNT6','LOX','SIGLEC7',
                    'MAOA','TREM1','CTSD','IGFBP3','MIF','EPOR','INHA','HSD11B1','IGF1','ACTA2','IGFBP3','IGFBP5',
                    'SEMA7A','SERPINB2','WNT1', 'WNT6', 'WNT10B', 'WNT3', 'WNT5A', 'WNT4', 'WNT2', 'WNT2B'
                    , 'WNT9A', 'SFRP4', 'SFRP1', 'CTNNB1', 'HOXA10', 'FOXO1', 'IL15','SIGLEC9')   

label_rows <- res_df$gene %in% selected_genes

# create the volcano plot with all points and labels
p <- ggplot(res_df, aes(x = capped_log2FoldChange, y = capped_logp, color = sig)) +
  geom_point(alpha = 0.8, size = 1.75) +
  geom_point(
    data = res_df[label_rows, ],
    aes(x = capped_log2FoldChange, y = capped_logp),
    shape = 21,              # shape 21 allows a colored fill and border was 21
    fill = "darkgreen",         # fill color for highlighted points
    color = "black",         # border color for highlighted points
    size = 2,                # larger size for emphasis
    stroke = 1               # thickness of the border
  ) +
  scale_color_manual(values = c("Downregulated" = "lightgrey",
                                "Not Significant" = "lightgrey",
                                "Upregulated" = "lightgrey"),
                     guide='none') +
  theme_classic() +
  labs(title = "Volcano Plot decFIB vs. genFIB",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value") +
  geom_text_repel(
    data = res_df[label_rows, ],
    aes(label = gene),
    size = 2,
    color = "black",          # label text color
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "black",  # connector line color
    segment.size = 0.5, max.overlaps = 1000
  )

ggsave("./volcano_plot_decFIB.pdf", plot = p, device = "pdf", width = 5, height = 5)

print(p)




res_df = res_bulk_df

# create capped versions of variables
res_df$capped_log2FoldChange <- with(res_df, ifelse(log2FoldChange > 10, 10,
                                                    ifelse(log2FoldChange < -10, -10, log2FoldChange)))

res_df$padj_fixed <- ifelse(res_df$padj == 0, 1e-300, res_df$padj)
res_df$logp <- -log10(res_df$padj_fixed)
res_df$capped_logp <- ifelse(res_df$logp > 100, 100, res_df$logp)

# define significance status
p_adj_cutoff <- 0.05
log2fc_cutoff <- 1
res_df$sig <- "Not Significant"
res_df$sig[!is.na(res_df$padj) & res_df$padj < p_adj_cutoff & res_df$log2FoldChange > log2fc_cutoff] <- "Upregulated"
res_df$sig[!is.na(res_df$padj) & res_df$padj < p_adj_cutoff & res_df$log2FoldChange < -log2fc_cutoff] <- "Downregulated"


label_rows <- res_df$gene %in% selected_genes

# create the volcano plot with all points and labels
p <- ggplot(res_df, aes(x = capped_log2FoldChange, y = capped_logp, color = sig)) +
  geom_point(alpha = 0.8, size = 1.75) +
  geom_point(
    data = res_df[label_rows, ],
    aes(x = capped_log2FoldChange, y = capped_logp),
    shape = 21,              # shape 21 allows a colored fill and border was 21
    fill = "darkgreen",         # fill color for highlighted points
    color = "black",         # border color for highlighted points
    size = 2,                # larger size for emphasis
    stroke = 1               # thickness of the border
  ) +
  scale_color_manual(values = c("Downregulated" = "lightgrey",
                                "Not Significant" = "lightgrey",
                                "Upregulated" = "lightgrey"),
                     guide='none') +
  theme_classic() +
  labs(title = "Volcano Plot PGE2+MPA vs. MPA",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value") +
  geom_text_repel(
    data = res_df[label_rows, ],
    aes(label = gene),
    size = 2,
    color = "black",          # label text color
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "black",  # connector line color
    segment.size = 0.5, max.overlaps = 1000
  )

ggsave("./volcano_plot_PGE2.pdf", plot = p, device = "pdf", width = 5, height = 5)

print(p)




res_df = res_bulk_df2

# create capped versions of variables
res_df$capped_log2FoldChange <- with(res_df, ifelse(log2FoldChange > 10, 10,
                                                    ifelse(log2FoldChange < -10, -10, log2FoldChange)))


res_df$padj_fixed <- ifelse(res_df$padj == 0, 1e-300, res_df$padj)
res_df$logp <- -log10(res_df$padj_fixed)
res_df$capped_logp <- ifelse(res_df$logp > 100, 100, res_df$logp)

# define significance status as before (adjust thresholds as needed)
p_adj_cutoff <- 0.05
log2fc_cutoff <- 1
res_df$sig <- "Not Significant"
res_df$sig[!is.na(res_df$padj) & res_df$padj < p_adj_cutoff & res_df$log2FoldChange > log2fc_cutoff] <- "Upregulated"
res_df$sig[!is.na(res_df$padj) & res_df$padj < p_adj_cutoff & res_df$log2FoldChange < -log2fc_cutoff] <- "Downregulated"


# create a logical vector for rows to label:
label_rows <- res_df$gene %in% selected_genes

# create the volcano plot with all points and labels
p <- ggplot(res_df, aes(x = capped_log2FoldChange, y = capped_logp, color = sig)) +
  geom_point(alpha = 0.8, size = 1.75) +
  # Highlight selected genes by plotting them with a different shape/size/color
  geom_point(
    data = res_df[label_rows, ],
    aes(x = capped_log2FoldChange, y = capped_logp),
    shape = 21,              # shape 21 allows a colored fill and border was 21
    fill = "darkgreen",         # fill color for highlighted points
    color = "black",         # border color for highlighted points
    size = 2,                # larger size for emphasis
    stroke = 1               # thickness of the border
  ) +
  scale_color_manual(values = c("Downregulated" = "lightgrey",
                                "Not Significant" = "lightgrey",
                                "Upregulated" = "lightgrey"),
                     guide='none') +
  theme_classic() +
  labs(title = "Volcano Plot cAMP+MPA vs. MPA",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value") +
  geom_text_repel(
    data = res_df[label_rows, ],
    aes(label = gene),
    size = 2,
    color = "black",          # label text color
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "black",  # connector line color
    segment.size = 0.5, max.overlaps = 1000
  )

ggsave("./volcano_plot_cAMP.pdf", plot = p, device = "pdf", width = 5, height = 5)

print(p)




rownames(res_bulk_df) = res_bulk_df$gene
rownames(res_pseudobulk_df3) = res_pseudobulk_df3$gene


# identify common genes
common_genes <- intersect(res_bulk_df$gene, res_pseudobulk_df3$gene)

# subset both results to the common genes
res1_common <- res_bulk_df[common_genes, ]
res2_common <- res_pseudobulk_df3[common_genes, ]

p_adj_cutoff <- 0.05
log2fc_cutoff <- 0
res1_common$sig <- 0
res1_common$sig[!is.na(res1_common$padj) & res1_common$padj < p_adj_cutoff & res1_common$log2FoldChange > log2fc_cutoff] <- 1
res1_common$sig[!is.na(res1_common$padj) & res1_common$padj < p_adj_cutoff & res1_common$log2FoldChange < -log2fc_cutoff] <- 1


p_adj_cutoff <- 0.05
log2fc_cutoff <- 0
res2_common$sig <- 0
res2_common$sig[!is.na(res2_common$padj) & res2_common$padj < p_adj_cutoff & res2_common$log2FoldChange > log2fc_cutoff] <- 1
res2_common$sig[!is.na(res2_common$padj) & res2_common$padj < p_adj_cutoff & res2_common$log2FoldChange < -log2fc_cutoff] <- 1


# create a data frame combining the log2 fold changes
df <- data.frame(
  gene = common_genes,
  log2FC1 = res1_common$log2FoldChange,
  log2FC2 = res2_common$log2FoldChange,
  sig1 = res1_common$sig,
  sig2 = res2_common$sig
)

df$sig = df$sig1 + df$sig2

df = df %>% dplyr::filter(.$sig == 2)

df$log2FC1 <- with(df, ifelse(log2FC1 > 10, 10,
                                                    ifelse(log2FC1 < -10, -10, log2FC1)))

df$sig <- "Not Significant"

p <- ggplot(df, aes(x = log2FC1, y = log2FC2, color = sig)) +
  geom_point(alpha = 0.8, size = 1.75) +
  geom_point(
    data = df%>%dplyr::filter(gene %in% selected_genes),
    aes(x = log2FC1, y = log2FC2),
    shape = 21,              # shape 21 allows a colored fill and border was 21
    fill = "darkgreen",         # fill color for highlighted points
    color = "black",         # border color for highlighted points
    size = 2,                # larger size for emphasis
    stroke = 1               # thickness of the border
  ) +
  scale_color_manual(values = c("Significant" = "blue",
                                "Not Significant" = "lightgrey"),
                     guide='none') +
  theme_classic() +
  labs(title = "Correlation of DESeq2 results (Log2FoldChanges)",
       x = "MPA + PGE2 vs. MPA",
       y = "decFIB vs. genFIB") +
  geom_text_repel(
    data = df%>%dplyr::filter(gene %in% selected_genes),
    aes(label = gene),
    size = 2,
    color = "black",          # label text color
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "black",  # connector line color
    segment.size = 0.5, max.overlaps = 1000
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")

ggsave("./corr_lfc_pge_decFIB.pdf", plot = p, device = "pdf", width = 5, height = 5)

print(p)


bulk_df = t(bulk_counts) %>% as.data.frame() %>% rownames_to_column('gene')

df_filt


df_filt %>% 
  dplyr::summarise(
    up_up = sum(log2FC1 > 0 & log2FC2 > 0),
    down_down = sum(log2FC1 < 0 & log2FC2 < 0),
    up_down = sum(log2FC1 > 0 & log2FC2 < 0),
    down_up = sum(log2FC1 < 0 & log2FC2 > 0)
  )

