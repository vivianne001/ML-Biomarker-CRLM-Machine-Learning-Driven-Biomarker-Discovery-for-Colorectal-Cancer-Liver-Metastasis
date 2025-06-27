######################### Quantitative Analysis and Visualization of IHC Score ############################

rm(list = ls())

# 1. Load Data
library(readxl)
IHC = read_excel("Supplement Tables.xlsx", col_names = TRUE, sheet = 10)
IHC = as.data.frame(IHC)
colnames(IHC) = IHC[1, ]
IHC = IHC[-1, ]
IHC = IHC[, -1]
IHC = as.data.frame(IHC)

# Convert IHC Score to numeric
IHC$`IHC Score` = as.numeric(IHC$`IHC Score`)

# Check group information
table(IHC$Group)
table(IHC$Outcome)

# Assign group for outcome
IHC$group = ifelse(IHC$Outcome == "PD", "PD", "DC")
table(IHC$group)

# 2. Barplot: Total IHC Score by Group
total_scores <- IHC %>%
  group_by(Group) %>%
  summarise(Total_IHC_Score = sum(`IHC Score`, na.rm = TRUE))

t_test_result <- t.test(`IHC Score` ~ Group, data = IHC)
p_value <- t_test_result$p.value
signif_label <- ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", "ns")))

library(ggplot2)
library(ggsignif)

ggplot(total_scores, aes(x = Group, y = Total_IHC_Score, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("CRLM" = "red", "CRC AIS" = "lightblue")) +
  labs(title = "ACMSD IHC Score Comparison between Groups of CRLM & CRC AIS",
       x = "Group",
       y = "ACMSD IHC Score") +
  scale_x_discrete(labels = c("CRLM" = "CRLM", "CRC AIS" = "CRC AIS")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_signif(comparisons = list(c("CRLM", "CRC AIS")),
              annotations = signif_label,
              y_position = max(total_scores$Total_IHC_Score) * 1.1,
              tip_length = 0.03)

# 3. Violin and Boxplot: IHC Score by Group
IHC$Group <- factor(IHC$Group, levels = c("LCRC", "CRLM"))

ggplot(IHC, aes(x = Group, y = `IHC Score`, fill = Group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  scale_fill_manual(values = c("CRLM" = "red", "LCRC" = "lightblue")) +
  geom_signif(comparisons = list(c("LCRC", "CRLM")),
              test = "wilcox.test",
              map_signif_level = TRUE,
              y_position = 5.2) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 14)) +
  ylab("ACMSD IHC Score") +
  xlab("group: LCRC & CRLM")

# 4. Boxplot: IHC Score by Outcome Group (PD vs DC)
IHC$group <- factor(IHC$group, levels = c("PD", "DC"))

ggplot(IHC, aes(x = group, y = `IHC Score`, fill = group)) +
  geom_boxplot(width = 0.75, outlier.shape = NA) +
  scale_fill_manual(values = c("PD" = "red", "DC" = "lightblue")) +
  geom_signif(comparisons = list(c("PD", "DC")),
              test = "wilcox.test",
              map_signif_level = TRUE) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 14)) +
  ylab("ACMSD IHC Score") +
  xlab("group: DC & PD")