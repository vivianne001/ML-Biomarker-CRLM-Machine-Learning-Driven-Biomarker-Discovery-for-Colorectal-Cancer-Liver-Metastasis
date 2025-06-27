# Optimal K Value Selection and Visualization for KNN Classification Model

library(e1071)
set.seed(123)

# Load training set features and labels
# Replace 'trainset' with your properly preprocessed data frame
x <- trainset[, -30]
y <- trainset[, 30]

# Parameter tuning using 10-fold cross-validation
tuning_result <- tune.knn(
  x, y,
  k = 1:10,
  tunecontrol = tune.control(sampling = "cross", cross = 10)
)

# Print and summarize tuning results
print(tuning_result)
summary(tuning_result)

# Visualization of tuning results
plot(tuning_result) # Plot error rates across K

# Add a red dashed line at the position of the optimal K (e.g., K = 7 as in result)
y_value_for_text <- 0.95 * max(tuning_result$results$error, na.rm = TRUE)

abline(v = 7, col = "red", lty = 2) # Mark K = 7 (or use tuning_result$best.parameters$k if generalized)
text(x = 7, y = y_value_for_text, labels = "the lowest error", col = "red", pos = 4)