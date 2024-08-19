# 加载必要的库
library(reshape2)
library(pheatmap)

# 读取数据
data <- read.table("data.txt", header = TRUE, sep = "\t")

# 修改列名以确保一致性
colnames(data) <- c("Type", "Data", "Name", "Completeness", "Contamination")


data_long <- melt(data, id.vars = c("Type","Name"), measure.vars = c("Completeness"))

# 检查数据格式
head(data_long)

# 重新组织数据以适应聚类分析
data_for_clustering <- dcast(data_long, Name ~ Type +  variable, value.var = "value")

# 检查数据格式
head(data_for_clustering)

# 去掉名称列以进行聚类分析
data_matrix <- as.matrix(data_for_clustering[, -1])

# 为行和列添加标签
rownames(data_matrix) <- data_for_clustering$Name

# 进行层次聚类分析
dist_matrix <- dist(scale(data_matrix))
hclust_result <- hclust(dist_matrix)

# 绘制聚类树状图
plot(hclust_result, labels = rownames(data_matrix), main = "Cluster Analysis of Microbiome, Remove, and Raw Data")

# 使用pheatmap包进行聚类可视化
pheatmap(data_matrix, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
         clustering_method = "complete", main = "Heatmap with MAG numbers")



# 定义自定义颜色
my_colors <- colorRampPalette(c("#67a9cf", "#ffffbf", "#ef8a62"))(100)

# 使用pheatmap包进行聚类可视化，设置比例尺为0-3
pheatmap(data_matrix, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete", 
         main = "Heatmap with MAG numbers",
         color = my_colors,
         breaks = seq(0, 3, length.out = 101)) # 设置比例尺范围为0到3




library(ggplot2)
library(car)
data <- read.table("data.txt", header = TRUE)

# 1) Plot Completeness

data$Type <- factor(data$Type, levels = c("Microbiome",  "Raw","Remove"))

ggplot(data, aes(x = Name, y = Completeness, fill = Type)) +
  geom_boxplot() +
  labs(x = "", y = "Completeness rate (%)", title = "", color = "black", fill = "Type") +
  theme_minimal() +
  theme_test(base_size = 20) +
  theme(
    panel.border = element_rect(size = 2, fill = 'transparent'),
    axis.text = element_text(color = 'black', size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better visibility
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5, vjust = 0.5),
    legend.title = element_text(size = 20, color = "black"),
    legend.text = element_text(size = 20, color = "black")
  )  +
  scale_fill_manual(values = c( "#eb746a","#2cad3f", "#6792cd"))


# 2) Plot Contamination

data$Type <- factor(data$Type, levels = c("Microbiome", "Remove", "Raw"))

ggplot(data, aes(x = Name, y = Contamination, fill = Type)) +
  geom_boxplot() +
  labs(x = "", y = "Contamination rate (%)", title = "", color = "black", fill = "Type") +
  theme_minimal() +
  theme_test(base_size = 20) +
  theme(
    panel.border = element_rect(size = 2, fill = 'transparent'),
    axis.text = element_text(color = 'black', size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better visibility
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5, vjust = 0.5),
    legend.title = element_text(size = 20, color = "black"),
    legend.text = element_text(size = 20, color = "black")
  )  +
  scale_fill_manual(values =  c( "#eb746a","#2cad3f", "#6792cd"))
