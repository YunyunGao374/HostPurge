data1=read.table('../eggnog/microbiome0.txt', sep = '', header = TRUE)
data2=read.table('../eggnog/raw0.txt', sep = '', header = TRUE)
data3=read.table('../eggnog/remove0.txt', sep = '', header = TRUE)
# 去掉每个表格的第一行
data1 <- data1[-1, ]
data2 <- data2[-1, ]
data3 <- data3[-1, ]
# 修改data2和data3的第二列的列名
colnames(data2)[2] <- 'raw0'
colnames(data3)[2] <- 'remove0'

df <- merge(data1, data2, by = 'KEGG_ko', all = TRUE)
df <- merge(df, data3, by = 'KEGG_ko', all = TRUE)
df[is.na(df)] <- 0


library(ggplot2)

# 计算microbiome0和raw0之间的回归并做图
df_raw <- df[, c("microbiome0", "raw0")]

cor_raw <- cor.test(df_raw$microbiome0, df_raw$raw0, method = "spearman")
m_raw <- lm(microbiome0 ~ raw0, data = df_raw)

p_raw <- ggplot(df_raw, aes(raw0, microbiome0)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = paste("rho = ", round(cor_raw$estimate, 3), 
                     ", P = ", signif(cor_raw$p.value, 3), 
                     ", R2 = ", round(summary(m_raw)$r.squared, 3), sep = "")) +
  theme_bw()

p_raw

ggsave("../microbiome0_vs_raw0.pdf", p_raw, width = 4, height = 2.5)

df_remove <- df[, c("microbiome0", "remove0")]

cor_remove <- cor.test(df_remove$microbiome0, df_remove$remove0, method = "spearman")
m_remove <- lm(microbiome0 ~ remove0, data = df_remove)

p_remove <- ggplot(df_remove, aes(remove0, microbiome0)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = paste("rho = ", round(cor_remove$estimate, 3), 
                     ", P = ", signif(cor_remove$p.value, 3), 
                     ", R2 = ", round(summary(m_remove)$r.squared, 3), sep = "")) +
  theme_bw()

p_remove
ggsave("../microbiome0_vs_remove0.pdf", p_remove, width = 4, height = 2.5)

###################################################################################
data1=read.table('../eggnog/microbiome1.txt', sep = '', header = TRUE)
data2=read.table('../eggnog/raw1.txt', sep = '', header = TRUE)
data3=read.table('../eggnog/remove1.txt', sep = '', header = TRUE)
# 去掉每个表格的第一行
data1 <- data1[-1, ]
data2 <- data2[-1, ]
data3 <- data3[-1, ]
# 修改data2和data3的第二列的列名
colnames(data2)[2] <- 'raw1'
colnames(data3)[2] <- 'remove1'

df <- merge(data1, data2, by = 'KEGG_ko', all = TRUE)
df <- merge(df, data3, by = 'KEGG_ko', all = TRUE)
df[is.na(df)] <- 0


library(ggplot2)

# 计算microbiome1和raw1之间的回归并做图
df_raw <- df[, c("microbiome1", "raw1")]

cor_raw <- cor.test(df_raw$microbiome1, df_raw$raw1, method = "spearman")
m_raw <- lm(microbiome1 ~ raw1, data = df_raw)

p_raw <- ggplot(df_raw, aes(raw1, microbiome1)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = paste("rho = ", round(cor_raw$estimate, 3), 
                     ", P = ", signif(cor_raw$p.value, 3), 
                     ", R2 = ", round(summary(m_raw)$r.squared, 3), sep = "")) +
  theme_bw()

p_raw

ggsave("../microbiome1_vs_raw1.pdf", p_raw, width = 6, height = 2.5)

df_remove <- df[, c("microbiome1", "remove1")]

cor_remove <- cor.test(df_remove$microbiome1, df_remove$remove1, method = "spearman")
m_remove <- lm(microbiome1 ~ remove1, data = df_remove)

p_remove <- ggplot(df_remove, aes(remove1, microbiome1)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = paste("rho = ", round(cor_remove$estimate, 3), 
                     ", P = ", signif(cor_remove$p.value, 3), 
                     ", R2 = ", round(summary(m_remove)$r.squared, 3), sep = "")) +
  theme_bw()

p_remove
ggsave("../microbiome1_vs_remove1.pdf", p_remove, width = 6, height = 2.5)

###################################################################################
data1=read.table('../eggnog/microbiome2.txt', sep = '', header = TRUE)
data2=read.table('../eggnog/raw2.txt', sep = '', header = TRUE)
data3=read.table('../eggnog/remove2.txt', sep = '', header = TRUE)
# 去掉每个表格的第一行
data1 <- data1[-1, ]
data2 <- data2[-1, ]
data3 <- data3[-1, ]
# 修改data2和data3的第二列的列名
colnames(data2)[2] <- 'raw2'
colnames(data3)[2] <- 'remove2'

df <- merge(data1, data2, by = 'KEGG_ko', all = TRUE)
df <- merge(df, data3, by = 'KEGG_ko', all = TRUE)
df[is.na(df)] <- 0


library(ggplot2)

# 计算microbiome2和raw2之间的回归并做图
df_raw <- df[, c("microbiome2", "raw2")]

cor_raw <- cor.test(df_raw$microbiome2, df_raw$raw2, method = "spearman")
m_raw <- lm(microbiome2 ~ raw2, data = df_raw)

p_raw <- ggplot(df_raw, aes(raw2, microbiome2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = paste("rho = ", round(cor_raw$estimate, 3), 
                     ", P = ", signif(cor_raw$p.value, 3), 
                     ", R2 = ", round(summary(m_raw)$r.squared, 3), sep = "")) +
  theme_bw()

p_raw

ggsave("../microbiome2_vs_raw2.pdf", p_raw, width = 6, height = 2.5)

df_remove <- df[, c("microbiome2", "remove2")]

cor_remove <- cor.test(df_remove$microbiome2, df_remove$remove2, method = "spearman")
m_remove <- lm(microbiome2 ~ remove2, data = df_remove)

p_remove <- ggplot(df_remove, aes(remove2, microbiome2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = paste("rho = ", round(cor_remove$estimate, 3), 
                     ", P = ", signif(cor_remove$p.value, 3), 
                     ", R2 = ", round(summary(m_remove)$r.squared, 3), sep = "")) +
  theme_bw()

p_remove
ggsave("../microbiome2_vs_remove2.pdf", p_remove, width = 6, height = 2.5)





