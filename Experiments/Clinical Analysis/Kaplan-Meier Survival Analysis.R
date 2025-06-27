
############# Kaplan-Meier Survival Analysis for Gene Expression Data ############
################################################################################

# Section 1: Setup and Environment Preparation

rm(list = ls())

# Load required libraries
# install.packages(c("survival", "survminer")) # Uncomment if needed
library(survival)
library(survminer)

# Section 2: Data Loading and Initial Filtering

# Load the pre-processed RData file containing expression (exp) and phenotype (ph) data
load("GSE17536_Raw_Data.Rdata")

# Filter for patients with AJCC stage 2 or 3
ph <- ph[ph$`ajcc_stage:ch1` %in% c("2", "3"), ]

# Section 3: Data Preprocessing and Cleaning

# Harmonize samples between expression and phenotype data
common_samples <- intersect(rownames(ph), colnames(exp))
ph <- ph[common_samples, ]
exp <- exp[, common_samples]

# Ensure sample order matches
identical(rownames(ph), colnames(exp))

# Prepare survival data
ph <- ph[, c("overall_event (death from any cause):ch1", "overall survival follow-up time:ch1")]
colnames(ph) <- c("OS_event", "OS_time")
ph <- na.omit(ph)

# Convert event status to numeric (1 = event, 0 = censored)
ph$OS_event <- ifelse(ph$OS_event == "recurrence", 1, 0)
ph$OS_time <- as.numeric(ph$OS_time)
ph$OS_time <- ph$OS_time / 12 # Convert months to years
ph <- na.omit(ph)

# Filter by follow-up time (e.g., 5-year survival)
ph <- ph[ph$OS_time < 5, ]

# Final harmonization
final_common_samples <- intersect(colnames(exp), rownames(ph))
exp <- exp[, final_common_samples]
ph <- ph[final_common_samples, ]
identical(colnames(exp), rownames(ph))

# Merge data for analysis
exp_t <- as.data.frame(t(exp))
exp_t$ID <- rownames(exp_t)
ph$ID <- rownames(ph)
merged_data <- merge(ph, exp_t, by = "ID")
rownames(merged_data) <- merged_data$ID
merged_data <- merged_data[, -1]

# Section 4: Kaplan-Meier Analysis (Optimal Cutoff Method)

# Find optimal cutpoints for gene expression
res_cut <- surv_cutpoint(
  data = merged_data,
  time = "OS_time",
  event = "OS_event",
  variables = c("SERPINA3", "PTGIS", "ACMSD", "SLC2A14"),
  minprop = 0.1,
  progressbar = TRUE
)
summary(res_cut)
plot(res_cut, "ACMSD", palette = "npg")

# Categorize samples based on optimal cutpoints
res_cat <- surv_categorize(res_cut)

# Kaplan-Meier analysis and plot (optimal cutoff)
fit_optimal <- survfit(Surv(OS_time, OS_event) ~ ACMSD, data = res_cat)
km_plot_optimal_cutoff <- ggsurvplot(
  fit_optimal,
  data = res_cat,
  pval = TRUE,
  risk.table = TRUE,
  surv.median.line = "hv",
  legend.title = "ACMSD Expression",
  legend.labs = c("Low", "High"),
  title = "5-Year OS by ACMSD Expression Level (Optimal Cutoff)",
  ylab = "Disease-Free Survival Probability",
  xlab = "Time (Years)",
  censor.shape = 124,
  censor.size = 2,
  conf.int = FALSE,
  break.x.by = 1,
  risk.table.y.text = FALSE,
  ggtheme = theme_classic()
)
print(km_plot_optimal_cutoff)

# Section 5: Kaplan-Meier Analysis (Median Cutoff Method)

# Create groups based on median expression
median_value <- median(merged_data$ACMSD)
merged_data$group_median <- ifelse(merged_data$ACMSD > median_value, "high", "low")

# Kaplan-Meier analysis and plot (median cutoff)
fit_median <- survfit(Surv(OS_time, OS_event) ~ group_median, data = merged_data)
km_plot_median_cutoff <- ggsurvplot(
  fit_median,
  data = merged_data,
  pval = TRUE,
  risk.table = TRUE,
  legend.title = "ACMSD Expression",
  legend.labs = c("High", "Low"),
  title = "5-Year OS by ACMSD Expression Level (Median Cutoff)",
  ylab = "Disease-Free Survival Probability",
  xlab = "Time (Years)"
)
print(km_plot_median_cutoff)