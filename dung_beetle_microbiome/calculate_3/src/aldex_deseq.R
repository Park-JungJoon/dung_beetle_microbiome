library(phyloseq)
library(DESeq2)
library(ALDEx2)
library(ggplot2)
library(vegan)
library(ALDEx2)


analyze_microbiome <- function(input_path, output_dir) {
  # 파일명에서 종명 추출
  species_name <- gsub("_cli_formatted.tsv", "", basename(input_path))
  
  # 결과 저장 디렉토리 생성
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("Reading data: ", input_path)
  raw <- read.delim(input_path, row.names = 1, check.names = FALSE)
  
  # 메타데이터 추출
  meta <- data.frame(SiteType = as.character(unlist(raw["Class", ])))
  rownames(meta) <- colnames(raw)
  meta <- meta[-1, , drop = FALSE]
  
  # OTU table for phyloseq / DESeq2
  otu <- raw[-1, ]
  otu <- apply(otu, 2, as.numeric)
  rownames(otu) <- rownames(raw)[-1]
  otu_physeq <- otu_table(otu, taxa_are_rows = TRUE)
  ps <- phyloseq(otu_physeq, sample_data(meta))
  
  # DESeq2 분석
  dds <- phyloseq_to_deseq2(ps, ~ SiteType)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds)
  res <- results(dds)
  sig <- res[which(res$padj < 0.05), ]
  deseq_result <- sig[order(sig$log2FoldChange, decreasing = TRUE), ]
  
  # TSV 저장
  write.table(as.data.frame(deseq_result), 
              file = file.path(output_dir, paste0(species_name, "_DESeq2_result.tsv")),
              sep = "\t", quote = FALSE, col.names = NA)
  
  # ALDEx2 준비
  otu_mat <- apply(raw[-1, ], 2, as.numeric)
  rownames(otu_mat) <- rownames(raw)[-1]
  conds <- as.character(unlist(raw["Class", colnames(otu_mat)]))
  stopifnot(length(conds) == ncol(otu_mat))
  
  # ALDEx2 실행
  otu_clr <- aldex.clr(otu_mat, conds = conds, mc.samples = 128, denom = "all", verbose = TRUE)
  otu_tt <- aldex.ttest(otu_clr)
  otu_effect <- aldex.effect(otu_clr)
  otu_all <- data.frame(otu_tt, otu_effect)
  otu_sig <- subset(otu_all, we.eBH < 0.05)
  
  # TSV 저장
  write.table(otu_all, 
              file = file.path(output_dir, paste0(species_name, "_ALDEx2_all.tsv")),
              sep = "\t", quote = FALSE, col.names = NA)
  
  write.table(otu_sig, 
              file = file.path(output_dir, paste0(species_name, "_ALDEx2_significant.tsv")),
              sep = "\t", quote = FALSE, col.names = NA)
  
  message("Analysis complete for: ", species_name)
}

analyze_microbiome(
  input_path = "calculate_3/output/Onthophagus_nuchicornis_cli_formatted.tsv",
  output_dir = "calculate_3/output/aldex_deg"
)

analyze_microbiome(
  input_path = "calculate_3/output/Bodilus_sordescens_cli_formatted.tsv",
  output_dir = "calculate_3/output/aldex_deg"
)

analyze_microbiome(
  input_path="calculate_3/output/Onthophagus_gibbulus_cli_formatted.tsv",
  output_dir="calculate_3/output/aldex_deg"
)

analyze_microbiome(
  input_path="calculate_3/output/Gymnopleurus_mopsus_cli_formatted.tsv",
  output_dir="calculate_3/output/aldex_deg"
)
