###################### Visualization of Drug Response Proportions by Pathological Stage in GSE72970 #########################

rm(list = ls())

library(readxl)
library(dplyr)
library(ggplot2)
library(scales)

# Filter for CRLM samples
# load("GSE72970_processed.Rdata")

table(ph$`synchronous metastase:ch1`)
ph = ph[ph$`synchronous metastase:ch1` == "Yes", ]
k = intersect(colnames(exp), rownames(ph))
dat = exp[, k]
ph = ph[k, ]
identical(colnames(dat), rownames(ph))

# Group assignment based on treatment regimen
table(ph$characteristics_ch1.8)
table(ph$`response category:ch1`)
ph$group = ifelse(
  ph$characteristics_ch1.8 %in% c(
    "regimen: FOLFIRI+BEVACIZUMAB",
    "regimen: FOLFIRI+ERBITUX",
    "regimen: FOLFIRINOX+BEVACIZUMAB",
    "regimen: FOLFOX+BEVACIZUMAB"
  ),
  "target", "chemotherapy"
)
table(ph$group)

# Prepare data for plotting
dat = as.data.frame(t(dat))
dat$group = ph$group
dat$outcome = ph$`response category:ch1`

# Group by ACMSD expression median
median_score <- median(dat$ACMSD)
dat$Group <- ifelse(dat$ACMSD >= median_score, "high", "low")
table(dat$Group)

# Separate by treatment group
dat1 = dat[dat$group == "target", ]
dat2 = dat[dat$group == "chemotherapy", ]

# Plot for target therapy group
group_counts <- table(dat1$Group)
x_labels <- paste(names(group_counts), "(n=", group_counts, ")", sep = "")
dat_summary <- dat1 %>%
  group_by(Group, outcome) %>%
  summarise(count = n()) %>%
  mutate(prop = count / sum(count))
custom_colors <- c("PD" = "#BA3E45", "SD" = "#20AEDD", "PR" = "#B6D7E9")
dat_summary <- dat_summary %>%
  mutate(outcome = factor(outcome, levels = c("PD", "SD", "PR")))

p <- ggplot(dat_summary, aes(x = Group, y = prop, fill = outcome)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = percent(prop)), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Group", y = "Proportion", fill = "Treatment Response", title = "Targeted Therapy Response") +
  theme_minimal() +
  theme(axis.line = element_line()) +
  scale_x_discrete(labels = x_labels)
print(p)

# Plot for chemotherapy group
group_counts <- table(dat2$Group)
x_labels <- paste(names(group_counts), "(n=", group_counts, ")", sep = "")
dat_summary <- dat2 %>%
  group_by(Group, outcome) %>%
  summarise(count = n()) %>%
  mutate(prop = count / sum(count))
custom_colors <- c("PD" = "#BA3E45", "SD" = "#20AEDD", "PR" = "#B6D7E9", "CR" = "#E1F3FB")
dat_summary <- dat_summary %>%
  mutate(outcome = factor(outcome, levels = c("PD", "SD", "PR", "CR")))

p <- ggplot(dat_summary, aes(x = Group, y = prop, fill = outcome)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = percent(prop)), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Group", y = "Proportion", fill = "Treatment Response", title = "Chemotherapy Response") +
  theme_minimal() +
  theme(axis.line = element_line()) +
  scale_x_discrete(labels = x_labels)
print(p)