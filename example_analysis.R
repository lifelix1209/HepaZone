# ==============================================================================
# HepaZone 使用示例
# 使用模拟的肝脏单细胞数据演示完整的空间转录组学重建流程
# ==============================================================================

library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)

# 加载 HepaZone 函数
source("R/data_io.R")
source("R/preprocessing.R")
source("R/zonation.R")
source("R/statistics.R")
source("R/visualization.R")
source("R/main.R")

cat("=============================================================\n")
cat("  HepaZone 肝脏空间转录组学重建 - 使用示例\n")
cat("=============================================================\n\n")

# ==============================================================================
# 1. 模拟肝脏单细胞数据
# ==============================================================================
cat(">>> 步骤1: 创建模拟的肝脏单细胞RNA-seq数据\n\n")

set.seed(42)
n_genes <- 500
n_cells <- 200

# 创建模拟表达矩阵
# 模拟肝小叶不同区域（中央静脉 vs 门静脉）的基因表达模式
counts <- matrix(rpois(n_genes * n_cells, lambda = 5), nrow = n_genes)

# 添加一些空间模式
# 中央静脉标志基因 (CV): 沿空间轴高表达
cv_genes <- paste0("Cyp", 1:20)
# 门静脉标志基因 (PN): 沿空间轴低表达
pn_genes <- paste0("Alb", 1:20)
# 其他基因
other_genes <- setdiff(1:n_genes, c(
  match(cv_genes, rownames(counts)),
  match(pn_genes, rownames(counts))
))
other_genes <- other_genes[1:(n_genes - 40)]
other_genes <- paste0("Gene", other_genes)

# 设置基因名
all_genes <- c(cv_genes, pn_genes, other_genes)
rownames(counts) <- all_genes
colnames(counts) <- paste0("Cell_", 1:n_cells)

# 添加空间位置信息（用于验证）
spatial_pos <- runif(n_cells)  # 0-1, 0=CV, 1=PV
cat("空间位置范围:", min(spatial_pos), "-", max(spatial_pos), "\n")

# 根据空间位置调整表达量
for (i in 1:length(cv_genes)) {
  if (cv_genes[i] %in% rownames(counts)) {
    # CV基因在CV区域高表达
    counts[cv_genes[i], ] <- counts[cv_genes[i], ] * (1 - spatial_pos) * 3 + 1
  }
}
for (i in 1:length(pn_genes)) {
  if (pn_genes[i] %in% rownames(counts)) {
    # PN基因在PV区域高表达
    counts[pn_genes[i], ] <- counts[pn_genes[i], ] * spatial_pos * 3 + 1
  }
}

# 创建Seurat对象
seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "Liver_Zonation_Example",
  min.cells = 3,
  min.features = 50
)

cat("创建Seurat对象:\n")
cat("  - 基因数:", nrow(seurat_obj), "\n")
cat("  - 细胞数:", ncol(seurat_obj), "\n\n")

# ==============================================================================
# 2. 数据预处理
# ==============================================================================
cat(">>> 步骤2: 数据预处理\n\n")

seurat_obj <- preprocess_zonation(
  seurat_obj,
  mt_pattern = "^Mt",
  remove_mup = TRUE,
  mup_pattern = "^Mup"
)

cat("\n")

# ==============================================================================
# 3. 计算空间位置分数
# ==============================================================================
cat(">>> 步骤3: 计算空间位置分数 (CL Score)\n\n")

# 使用CV和PN标志基因
cv_markers <- cv_genes[1:10]  # 使用前10个CV标志基因
pn_markers <- pn_genes[1:10]  # 使用前10个PN标志基因

seurat_obj <- calculate_spatial_position(
  seurat_obj,
  cv_markers = cv_markers,
  pn_markers = pn_markers
)

cat("\nCL Score 统计:\n")
cat("  - 平均值:", mean(seurat_obj$CL_score), "\n")
cat("  - 标准差:", sd(seurat_obj$CL_score), "\n")
cat("  - 范围:", min(seurat_obj$CL_score), "-", max(seurat_obj$CL_score), "\n\n")

# ==============================================================================
# 4. 映射细胞到空间层
# ==============================================================================
cat(">>> 步骤4: 映射细胞到空间层 (10个zone)\n\n")

n_zones <- 10
prob_matrix <- map_cells_to_layers(seurat_obj, n_zones = n_zones)

cat("概率矩阵维度:", nrow(prob_matrix), "x", ncol(prob_matrix), "\n")
cat("每层平均细胞数:\n")
zone_counts <- colSums(prob_matrix)
for (z in 1:n_zones) {
  cat(sprintf("  Zone %d: %.1f\n", z, zone_counts[z]))
}
cat("\n")

# ==============================================================================
# 5. 重建空间表达谱
# ==============================================================================
cat(">>> 步骤5: 重建基因空间表达谱\n\n")

spatial_result <- reconstruct_spatial_expression(seurat_obj, prob_matrix)

cat("空间表达矩阵维度:", nrow(spatial_result$mean_expression), "x",
    ncol(spatial_result$mean_expression), "\n")
cat("  (基因 x 空间层)\n\n")

# 显示一些关键基因的表达谱
key_genes <- c("Cyp1", "Cyp2", "Alb1", "Alb2")
key_genes <- key_genes[key_genes %in% rownames(spatial_result$mean_expression)]
if (length(key_genes) > 0) {
  cat("关键基因表达谱 (前5个zone):\n")
  for (g in key_genes[1:min(3, length(key_genes))]) {
    cat(sprintf("  %s: ", g))
    cat(paste0(round(spatial_result$mean_expression[g, 1:5], 3), " "), "\n")
  }
  cat("\n")
}

# ==============================================================================
# 6. Bootstrap 估计标准误差
# ==============================================================================
cat(">>> 步骤6: Bootstrap标准误差估计 (50次迭代)\n\n")

# 从Seurat对象获取归一化的表达矩阵
mat_norm <- .get_data(seurat_obj)

boot_result <- bootstrap_se(
  expression_matrix = mat_norm,
  prob_matrix = prob_matrix,
  n_bootstrap = 50,
  seed = 42
)

cat("Bootstrap完成!\n")
cat("SE矩阵维度:", nrow(boot_result$se), "x", ncol(boot_result$se), "\n\n")

# ==============================================================================
# 7. 置换检验
# ==============================================================================
cat(">>> 步骤7: 置换检验 (100次置换)\n\n")

perm_result <- permutation_test(
  expression_matrix = mat_norm,
  prob_matrix = prob_matrix,
  n_permutations = 100,
  seed = 123
)

cat("置换检验完成!\n")

# 计算q值
qvals <- calculate_qvalues(perm_result$p_values)

# 识别显著的空间变异基因
sig_threshold <- 0.05
n_sig <- sum(qvals < sig_threshold, na.rm = TRUE)
cat(sprintf("显著空间变异基因 (q < %s): %d\n\n", sig_threshold, n_sig))

# ==============================================================================
# 8. 运行完整流程 (hepa_zone_reconstruct)
# ==============================================================================
cat(">>> 步骤8: 使用 hepa_zone_reconstruct() 运行完整流程\n\n")

# 创建新的数据用于完整流程演示
set.seed(123)
counts_full <- matrix(rpois(100 * n_cells, lambda = 5), nrow = 100)
all_genes_full <- c(cv_genes[1:10], pn_genes[1:10], paste0("OtherGene", 1:80))
rownames(counts_full) <- all_genes_full
colnames(counts_full) <- paste0("Cell_", 1:n_cells)

seurat_full <- CreateSeuratObject(
  counts = counts_full,
  project = "Full_Example",
  min.cells = 3,
  min.features = 50
)

# 运行完整流程
result <- hepa_zone_reconstruct(
  seurat_full,
  cv_markers = cv_genes[1:5],
  pn_markers = pn_genes[1:5],
  n_zones = 5,
  n_bootstrap = 50,
  n_permutations = 100,
  seed = 42,
  verbose = FALSE
)

cat("完整分析完成!\n")
cat("结果对象类型:", class(result), "\n")
cat("基因数:", nrow(result$mean_expression), "\n")
cat("空间层数:", result$n_zones, "\n")
cat("显著SVG数:", sum(result$svg_results$significant), "\n\n")

# ==============================================================================
# 9. 可视化结果
# ==============================================================================
cat(">>> 步骤9: 生成可视化图表\n\n")

# 创建输出目录
output_dir <- "example_output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 9.1 CL分数分布直方图
cat("  - 生成 CL Score 分布图...\n")
p1 <- plot_cl_distribution(seurat_obj, prob_matrix = prob_matrix)
ggsave(file.path(output_dir, "cl_score_distribution.png"), p1, width = 8, height = 6)
cat("    保存到:", file.path(output_dir, "cl_score_distribution.png"), "\n")

# 9.2 空间表达热图 (选择高变异基因)
cat("  - 生成空间表达热图...\n")
top_genes <- head(result$svg_results$gene, 10)
if (length(top_genes) > 0) {
  p2 <- plot_spatial_heatmap(result, genes = top_genes, n_genes = 10)
  ggsave(file.path(output_dir, "spatial_expression_heatmap.png"), p2, width = 10, height = 8)
  cat("    保存到:", file.path(output_dir, "spatial_expression_heatmap.png"), "\n")
}

# 9.3 基因表达梯度图
cat("  - 生成基因表达梯度图...\n")
# 选择CV和PN标志基因
plot_genes <- c(cv_genes[1], pn_genes[1])
plot_genes <- plot_genes[plot_genes %in% rownames(result$mean_expression)]
if (length(plot_genes) > 0) {
  p3 <- plot_gradient(result, genes = plot_genes)
  ggsave(file.path(output_dir, "gene_expression_gradient.png"), p3, width = 8, height = 6)
  cat("    保存到:", file.path(output_dir, "gene_expression_gradient.png"), "\n")
}

# 9.4 CV vs PN 标志基因表达
cat("  - 生成 CV vs PN 标志基因表达图...\n")
p4 <- plot_marker_expression(seurat_obj, cv_genes[1:5], pn_genes[1:5])
ggsave(file.path(output_dir, "marker_expression.png"), p4, width = 8, height = 6)
cat("    保存到:", file.path(output_dir, "marker_expression.png"), "\n")

# 9.5 时间点比较 (模拟)
cat("  - 生成时间点比较图...\n")
# 创建模拟的多时间点结果
result_list <- list(
  ZT00 = result,
  ZT12 = result
)
# 稍微修改第二个结果以模拟差异
result_list$ZT12$mean_expression <- result$mean_expression * runif(length(result$mean_expression), 0.8, 1.2)

if (any(rownames(result$mean_expression) %in% c(cv_genes[1], pn_genes[1]))) {
  p5 <- plot_time_comparison(result_list, gene = cv_genes[1],
                             time_points = c("ZT00", "ZT12"))
  ggsave(file.path(output_dir, "time_comparison.png"), p5, width = 8, height = 6)
  cat("    保存到:", file.path(output_dir, "time_comparison.png"), "\n")
}

# ==============================================================================
# 10. 导出结果
# ==============================================================================
cat("\n>>> 步骤10: 导出分析结果\n\n")

# 导出完整结果
export_hepa_zone_results(result, output.dir = output_dir)

# 单独导出关键结果
cat("  - 导出空间表达矩阵...\n")
export_spatial_expression(
  result$mean_expression,
  format = "csv",
  filepath = file.path(output_dir, "spatial_expression_matrix")
)

cat("  - 导出q值矩阵...\n")
export_spatial_expression(
  result$qvalues,
  format = "csv",
  filepath = file.path(output_dir, "qvalues_matrix")
)

# ==============================================================================
# 结果摘要
# ==============================================================================
cat("\n=============================================================\n")
cat("  分析完成! 结果摘要\n")
cat("=============================================================\n\n")

cat("输出目录:", output_dir, "\n\n")

cat("生成的文件:\n")
cat("  1. cl_score_distribution.png   - CL分数分布直方图\n")
cat("  2. spatial_expression_heatmap.png - 空间表达热图\n")
cat("  3. gene_expression_gradient.png   - 基因表达梯度图\n")
cat("  4. marker_expression.png          - CV vs PN标志基因图\n")
cat("  5. time_comparison.png            - 时间点比较图\n")
cat("  6. mean_expression.csv            - 空间表达矩阵\n")
cat("  7. qvalues.csv                    - q值矩阵\n")
cat("  8. hepa_zone_full_results.rds     - 完整结果对象\n\n")

cat("关键统计:\n")
cat(sprintf("  - 分析细胞数: %d\n", ncol(seurat_obj)))
cat(sprintf("  - 分析基因数: %d\n", nrow(seurat_obj)))
cat(sprintf("  - 空间层数: %d\n", result$n_zones))
cat(sprintf("  - 显著空间变异基因数 (q<0.05): %d\n\n", sum(result$svg_results$significant)))

# 显示Top 10 SVG
cat("Top 10 空间变异基因:\n")
top10 <- head(result$svg_results, 10)
for (i in 1:nrow(top10)) {
  cat(sprintf("  %2d. %-15s q=%.4f\n", i, top10$gene[i], top10$q_value[i]))
}

cat("\n=============================================================\n")
cat("  示例分析完成!\n")
cat("=============================================================\n")
