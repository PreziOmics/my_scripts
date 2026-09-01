# ==============================================================================
# RNA-Seq & Functional Enrichment Analysis Utility Functions
# ==============================================================================

# Core Libraries
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(edgeR)

# ------------------------------------------------------------------------------
# 1. GSEA Analysis Function
# ------------------------------------------------------------------------------
GSEA_function <- function(dir_path,
                          df,
                          prefix,
                          pair_1,
                          pair_2,
                          ranking_by = "logFC",
                          species = "hs",
                          ontology_type = "BP",
                          padj_cutoff = 0.05,
                          num_categories = 10,
                          font_size = 10) {
  
  # Select Annotation Database
  if (species == "mm") {
    require(org.Mm.eg.db)
    organism <- org.Mm.eg.db
  } else if (species == "hs") {
    require(org.Hs.eg.db)
    organism <- org.Hs.eg.db
  } else {
    stop("Unsupported species. Use 'hs' or 'mm'.")
  }
  
  # Clean ENSEMBL IDs
  df$gene_id_clean <- sub("\\..*", "", df$gene_id)
  
  # Rank Genes
  if (ranking_by == "logFC") {
    geneList <- df$logFC
  } else {
    geneList <- df$FC_and_FDR
  }
  names(geneList) <- df$gene_id_clean
  
  geneList <- geneList[!is.na(geneList)]
  geneList <- geneList[!duplicated(names(geneList))]
  geneList <- sort(geneList, decreasing = TRUE)
  
  # Run GSEA
  gse <- gseGO(geneList = geneList,
               ont = ontology_type,
               keyType = "ENSEMBL",
               minGSSize = 10,
               maxGSSize = 2000,
               pvalueCutoff = padj_cutoff,
               pAdjustMethod = "BH",
               OrgDb = organism,
               verbose = FALSE,
               eps = 0,
               by = "fgsea")
  
  # Subset Up / Down
  gse_pos <- gse
  gse_pos@result <- gse@result[gse@result$NES > 0 & gse@result$p.adjust <= padj_cutoff, ]
  
  gse_neg <- gse
  gse_neg@result <- gse@result[gse@result$NES < 0 & gse@result$p.adjust <= padj_cutoff, ]
  
  # Save Table
  gse_df <- as.data.frame(gse)
  gse_df <- gse_df[order(gse_df$p.adjust), ]
  
  file_name <- paste0("gsea_", prefix, "_", pair_2, "_vs_", pair_1, "_", ontology_type, ".tsv")
  write.table(gse_df, file.path(dir_path, file_name), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Plot
  pdf_path <- file.path(dir_path, paste0("GSEA_dotplot_", prefix, "_", pair_2, "_vs_", pair_1, "_", ontology_type, ".pdf"))
  pdf(pdf_path, width = 6, height = 7)
  
  if (nrow(gse_pos@result) > 0) {
    print(dotplot(gse_pos, font.size = font_size, title = "Up-regulated", orderBy = "p.adjust", showCategory = num_categories))
  }
  if (nrow(gse_neg@result) > 0) {
    print(dotplot(gse_neg, font.size = font_size, title = "Down-regulated", orderBy = "p.adjust", showCategory = num_categories))
  }
  
  dev.off()
  return(gse)
}

# ------------------------------------------------------------------------------
# 2. EnrichGO Analysis Function
# ------------------------------------------------------------------------------
enrichGO_function <- function(dir_path, df, pair_1, pair_2, species = "hs", ontology_type = "BP", 
                              DE_genes_padj_cutoff = 0.05, GO_padj_cutoff = 0.05, 
                              num_categories = 10, font_size = 10) {
  
  if (species == "mm") {
    require(org.Mm.eg.db)
    organism <- "org.Mm.eg.db"
  } else if (species == "hs") {
    require(org.Hs.eg.db)
    organism <- "org.Hs.eg.db"
  }
  
  upregulated_df <- df[df$logFC > 0 & df$FDR < DE_genes_padj_cutoff, ]
  downregulated_df <- df[df$logFC < 0 & df$FDR < DE_genes_padj_cutoff, ]
  
  if (nrow(upregulated_df) == 0 && nrow(downregulated_df) == 0) {
    message("No significant genes found for enrichGO.")
    return(NULL)
  }
  
  upregulated <- sub("\\..*", "", upregulated_df$gene_id)
  downregulated <- sub("\\..*", "", downregulated_df$gene_id)
  
  # Enrich Down
  if (length(downregulated) > 0) {
    gse_down <- enrichGO(gene = downregulated, ont = ontology_type, keyType = "ENSEMBL",
                         OrgDb = organism, pAdjustMethod = "BH", pvalueCutoff = GO_padj_cutoff)
    write.table(as.data.frame(gse_down), file.path(dir_path, paste0("enrichGO_DOWN_", pair_2, "_vs_", pair_1, ".tsv")), 
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  # Enrich Up
  if (length(upregulated) > 0) {
    gse_up <- enrichGO(gene = upregulated, ont = ontology_type, keyType = "ENSEMBL",
                       OrgDb = organism, pAdjustMethod = "BH", pvalueCutoff = GO_padj_cutoff)
    write.table(as.data.frame(gse_up), file.path(dir_path, paste0("enrichGO_UP_", pair_2, "_vs_", pair_1, ".tsv")), 
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  # Plot
  pdf(file.path(dir_path, paste0("dotplot_enrichGO_", pair_2, "_vs_", pair_1, ".pdf")), width = 6, height = 7)
  if (exists("gse_up") && !is.null(gse_up)) print(dotplot(gse_up, showCategory = num_categories, title = "Up-regulated"))
  if (exists("gse_down") && !is.null(gse_down)) print(dotplot(gse_down, showCategory = num_categories, title = "Down-regulated"))
  dev.off()
}

# ------------------------------------------------------------------------------
# 3. Volcano Plot Function
# ------------------------------------------------------------------------------
volcano_plot_Genes_DiffExpr <- function(dir_path, df, pair_1, pair_2,
                                        FDR_cutoff = 0.05, log_FC = 1, 
                                        legend_name_up = "Up", legend_name_down = "Down", 
                                        labels = "yes", label_FDR_cutoff = 0.01, label_FC_cutoff = 2) {
  
  stopifnot(all(c("logFC", "FDR", "gene_name") %in% colnames(df)))
  
  df$FDR[df$FDR == 0] <- .Machine$double.xmin
  
  df$color <- "NS"
  df$color[df$FDR <= FDR_cutoff & df$logFC > log_FC] <- legend_name_up
  df$color[df$FDR <= FDR_cutoff & df$logFC < -log_FC] <- legend_name_down
  
  p <- ggplot(df, aes(x = logFC, y = -log10(FDR), color = color)) +
    geom_point(size = 1, alpha = 0.8) +
    scale_color_manual(values = setNames(c("red", "blue", "grey"), c(legend_name_up, legend_name_down, "NS"))) +
    theme_classic() +
    labs(title = paste0(pair_2, " vs ", pair_1), x = "log2 Fold Change", y = "-log10(FDR)")
  
  if (labels == "yes") {
    label_genes <- subset(df, FDR <= label_FDR_cutoff & abs(logFC) >= label_FC_cutoff)
    p <- p + geom_text_repel(data = label_genes, aes(label = gene_name), size = 2.5, max.overlaps = 15)
  }
  
  pdf(file.path(dir_path, paste0("volcano_", pair_2, "_vs_", pair_1, ".pdf")), width = 6, height = 6)
  print(p)
  dev.off()
  
  return(df)
}

# ------------------------------------------------------------------------------
# 4. Heatmap / Clustergram Function
# ------------------------------------------------------------------------------
clustergram <- function(significative_genes_cpm, dir_path, width = 6, height = 8,
                        file_name_prefix = "heatmap", plot_title_prefix = "Heatmap",
                        scale_type = "row", by_rows = TRUE, by_cols = TRUE,
                        size_row = 6, size_col = 8) {
  
  pdf(file.path(dir_path, paste0(file_name_prefix, ".pdf")), width = width, height = height)
  hm <- pheatmap(
    significative_genes_cpm, 
    main = plot_title_prefix,
    scale = scale_type,               
    clustering_method = "ward.D", 
    cluster_rows = by_rows,         
    cluster_cols = by_cols,         
    show_rownames = TRUE,        
    show_colnames = TRUE,        
    border_color = NA,
    color = colorRampPalette(c("darkblue", "lightskyblue", "white", "tomato", "darkred"))(250),
    fontsize_row = size_row,
    fontsize_col = size_col,
    angle_col = 45
  )
  print(hm)
  dev.off()
}

# ------------------------------------------------------------------------------
# 5. Differential Expression with edgeR
# ------------------------------------------------------------------------------
differential_expression_EdgeR <- function(df, dir_path, file_prefix, pair_1, pair_2) {
  yy <- estimateDisp(df)
  dge <- exactTest(yy, pair = c(pair_1, pair_2))
  merged_df <- as.data.frame(dge$table)
  merged_df$FDR <- p.adjust(merged_df$PValue, method = "BH")
  merged_df <- merged_df[order(merged_df$PValue), ]
  merged_df$gene_name <- rownames(merged_df)
  
  merged_df <- merged_df[, c("gene_name", "logFC", "logCPM", "PValue", "FDR")]
  write.table(merged_df, file.path(dir_path, paste0(file_prefix, "_", pair_2, "_vs_", pair_1, ".tsv")), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(merged_df)
}