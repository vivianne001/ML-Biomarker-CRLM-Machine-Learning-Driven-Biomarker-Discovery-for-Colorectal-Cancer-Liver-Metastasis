# Visualization of Confusion Matrix for Classification Results

# Load required package
if (!require(ggplot2)) {
  install.packages("ggplot2")
}
library(ggplot2)

# Confusion matrix for validation set
confusion_validation <- matrix(c(65, 1, 14, 24), nrow = 2)
colnames(confusion_validation) <- c("0", "1")
rownames(confusion_validation) <- c("0", "1")
df_validation <- as.data.frame(as.table(confusion_validation))

p_validation <- ggplot(df_validation, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "#FBDFE2", high = "#B83945") +
  labs(x = "Predicted Label", y = "True Label", fill = "Count") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  ) +
  geom_text(aes(label = Freq), size = 8, color = "black")

print(p_validation)

# Confusion matrix for test set
confusion_test <- matrix(c(35, 4, 2, 92), nrow = 2)
colnames(confusion_test) <- c("0", "1")
rownames(confusion_test) <- c("0", "1")
df_test <- as.data.frame(as.table(confusion_test))

p_test <- ggplot(df_test, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "#FBDFE2", high = "#B83945") +
  labs(x = "Predicted Label", y = "True Label", fill = "Count") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  ) +
  geom_text(aes(label = Freq), size = 8, color = "black")

print(p_test)