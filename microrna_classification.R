# MICRORNA-BASED CLASSIFICATION OF HODGKIN LYMPHOMA
# Exploratory analysis, group comparisons and binomial logistic regression

library(ggplot2)
library(pROC)

# OUTPUT DIRECTORY
# Save figures in the project folder.
dir.create("figures", showWarnings = FALSE)

# DATASET
# MicroRNA expression measurements from 52 subjects in three groups
# The file can be replaced with a comparable dataset by updating the
# variable names and group definitions below.
# Select the data file manually if it is not in the working directory

if (file.exists("linfomi.txt")) {
  data = read.delim("linfomi.txt")} else {data = read.delim(file.choose())}

head(data)
dim(data)
str(data)
summary(data)
colSums(is.na(data))

mirna = c("miR221", "miR222", "miR146b", "miR26a")
required_columns = c("Type", mirna)

if (!all(required_columns %in% names(data))) {
  stop("The dataset does not contain the required variables.")}

# Group definitions: 0 = Healthy, 1 = Hodgkin, 2 = Non-Hodgkin.
data$Type = factor(as.character(data$Type),
                   levels = c("0", "1", "2"),
                   labels = c("Healthy", "Hodgkin", "Non-Hodgkin"))

if (any(is.na(data$Type))) {
  stop("Unexpected or missing group labels. Check the Type variable.")}

table(data$Type)
prop.table(table(data$Type))

# Log transformations require strictly positive expression values
if (any(as.matrix(data[, mirna]) <= 0, na.rm = TRUE)) {
  stop("MicroRNA measurements must be positive for log transformation.")}

# GRAPHIC SETTINGS
colors = c("#E18B45", "#A9C5E0", "#72A6C6")

# EXPLORATORY ANALYSIS
summary(data[, mirna])
apply(data[, mirna], 2, sd, na.rm = TRUE)
apply(data[, mirna], 2, IQR, na.rm = TRUE)

# FIGURE 1: Original microRNA expression.
boxplot(data[, mirna], main = "MicroRNA expression",
        col = "#A9C5E0", ylab = "Expression level")
png("figures/01_mirna_expression.png",
    width = 1200, height = 800, res = 150)

boxplot(data[, mirna], main = "MicroRNA expression",
        col = "#A9C5E0", ylab = "Expression level")
dev.off()

# FIGURE 2: Log-transformed microRNA expression.
boxplot(log(data[, mirna]), main = "Log-transformed expression",
        col = "#A9C5E0", ylab = "Log expression")
png("figures/02_log_expression.png",
    width = 1200, height = 800, res = 150)

boxplot(log(data[, mirna]), main = "Log-transformed expression",
        col = "#A9C5E0", ylab = "Log expression")
dev.off()

# NORMALITY ASSESSMENT: 
# Shapiro-Wilk tests and Q-Q plots on log-transformed expression
par(mfrow = c(2, 2))

for (variable in mirna) {
  cat("\nNormality test:", variable, "\n")
  print(shapiro.test(na.omit(log(data[[variable]]))))
  qqnorm(log(data[[variable]]), main = variable)
  qqline(log(data[[variable]]), col = "#E18B45")}
par(mfrow = c(1, 1))

# CORRELATION ANALYSIS
# Spearman correlations among the four microRNA markers
cor(data[, mirna], method = "spearman", use = "complete.obs")
for (i in 1:(length(mirna) - 1)) {
  for (j in (i + 1):length(mirna)) {
    cat("\nCorrelation:", mirna[i], "vs", mirna[j], "\n")
    print(cor.test(data[[mirna[i]]], data[[mirna[j]]],
                   method = "spearman", exact = FALSE))}}

# FIGURE 3: MICRORNA EXPRESSION BY CLINICAL GROUP
# The figure contains four panels, one for each microRNA
# Each panel compares Healthy, Hodgkin and Non-Hodgkin subjects
par(mfrow = c(2, 2), mar = c(7, 4, 3, 1))
for (variable in mirna) {
  boxplot(log(data[[variable]]) ~ data$Type,
          main = variable,
          xlab = "Clinical group", ylab = "Log expression",
          col = colors, las = 2)}
par(mfrow = c(1, 1))

png("figures/03_expression_by_group.png",
    width = 1500, height = 1100, res = 150)

par(mfrow = c(2, 2), mar = c(7, 4, 3, 1))
for (variable in mirna) {
  boxplot(log(data[[variable]]) ~ data$Type,
          main = variable,
          xlab = "Clinical group", ylab = "Log expression",
          col = colors, las = 2)}
dev.off()
par(mfrow = c(1, 1))

# THREE-GROUP COMPARISONS
# Kruskal-Wallis tests and pairwise Wilcoxon tests with Holm adjustment
# These comparisons are exploratory; they do not assume normality
for (variable in mirna) {
  cat("\nThree-group comparison:", variable, "\n")
  print(kruskal.test(log(data[[variable]]) ~ data$Type))
  print(pairwise.wilcox.test(log(data[[variable]]), data$Type,
                             p.adjust.method = "holm", exact = FALSE))}

# HODGKIN VS NON-HODGKIN COMPARISONS
# Kruskal-Wallis tests and pairwise Wilcoxon tests with Holm adjustment
# These comparisons are exploratory and do not assume normality
for (variable in mirna) {
  cat("\nThree-group comparison:", variable, "\n")
  print(kruskal.test(log(data[[variable]]) ~ data$Type))
  print(pairwise.wilcox.test(log(data[[variable]]),
                             data$Type,
                             p.adjust.method = "holm",
                             exact = FALSE))}

# HODGKIN VS NON-HODGKIN COMPARISONS
# Healthy controls are excluded only from this two-group comparison

groups = droplevels(subset(data, Type != "Healthy"))
pvalues = numeric(length(mirna))
for (i in seq_along(mirna)) {
  variable = mirna[i]
  test = wilcox.test(log(groups[[variable]]) ~ groups$Type,
                     exact = FALSE)
  pvalues[i] = test$p.value}

comparison = data.frame(
  microRNA = mirna,
  p_value = pvalues,
  p_holm = p.adjust(pvalues, method = "holm"),
  p_bonferroni = p.adjust(pvalues, method = "bonferroni"))
print(comparison)

# LOGISTIC REGRESSION

# Binary outcome:
# 1 = Hodgkin lymphoma
# 0 = Healthy controls or Non-Hodgkin lymphoma
# All three clinical groups are included in the model

data$case = as.integer(data$Type == "Hodgkin")
table(data$case)
model_data = data[complete.cases(data[, c("case", mirna)]), ]

# Full model with all four log-transformed microRNAs
full_formula = reformulate(paste0("log(", mirna, ")"), response = "case")

fit_full = glm(full_formula, family = binomial,data = model_data)

summary(fit_full)

# BACKWARD VARIABLE SELECTION
# Backward selection based on likelihood-ratio tests
# At each step, the least significant predictor is removed if p > 0.05

fit_final = fit_full

repeat {tests = drop1(fit_final, test = "Chisq")
  print(tests)
  pvalues_backward = tests[-1, "Pr(>Chi)"]
  if (length(pvalues_backward) == 0 ||
      all(is.na(pvalues_backward))) break
  worst = which.max(replace(pvalues_backward,
                            is.na(pvalues_backward), -Inf))
  if (pvalues_backward[worst] <= 0.05) break
  predictor = rownames(tests)[-1][worst]
  cat("\nRemoving:", predictor, "\n")
  fit_final = update(fit_final, paste(". ~ . -", predictor))}

# Final selected model
summary(fit_final)
formula(fit_final)

# MODEL COMPARISON
# Compare the final model with the intercept-only model

fit_null = glm(case ~ 1,
               family = binomial,
               data = model_data)

anova(fit_null, fit_final, test = "Chisq")
AIC(fit_final)

# Odds ratios and 95% Wald confidence intervals
exp(coef(fit_final))
exp(confint.default(fit_final))

# IN-SAMPLE CLASSIFICATION
# Predictions are evaluated on the same observations used for model fitting.

predicted_prob = predict(fit_final, type = "response")
predicted_class = ifelse(predicted_prob >= 0.5, 1, 0)

# Confusion matrix
confusion = table(
  Observed = factor(model_data$case, levels = c(0, 1)),
  Predicted = factor(predicted_class, levels = c(0, 1)))
print(confusion)

# Classification performance
TP = sum(model_data$case == 1 & predicted_class == 1)
TN = sum(model_data$case == 0 & predicted_class == 0)
FP = sum(model_data$case == 0 & predicted_class == 1)
FN = sum(model_data$case == 1 & predicted_class == 0)

accuracy = (TP + TN) / nrow(model_data)
sensitivity = TP / (TP + FN)
specificity = TN / (TN + FP)

c(accuracy = accuracy,
  sensitivity = sensitivity,
  specificity = specificity)

# ROC CURVE
# ROC analysis using in-sample predicted probabilities
roc_fit = roc(model_data$case, predicted_prob,
              levels = c(0, 1), direction = "<", quiet = TRUE)
auc(roc_fit)

roc_data = data.frame(FPR = 1 - roc_fit$specificities,
  TPR = roc_fit$sensitivities)

# Create the ROC plot
plot_roc = ggplot(roc_data, aes(x = FPR, y = TPR)) +
  geom_line(color = "#2864A0", linewidth = 1) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed",
              color = "grey55") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = paste0("Hodgkin classification (AUC = ",
                      round(as.numeric(auc(roc_fit)), 3), ")"),
       x = "1 - Specificity",
       y = "Sensitivity") +
  theme_classic(base_size = 12)
print(plot_roc)

# Save the ROC curve in the figures folder
ggsave("figures/04_roc_curve.png",
       plot = plot_roc,
       width = 7,
       height = 5,
       dpi = 300)


# CONCLUSIONS
# The analysis compares four microRNA markers across healthy controls,
# Hodgkin lymphoma patients and Non-Hodgkin lymphoma patients.
# Logistic regression distinguishes Hodgkin lymphoma (case = 1)
# from the combined Healthy / Non-Hodgkin group (case = 0).
# Backward selection based on likelihood-ratio tests retains
# log(miR222) as the only predictor in the final model.
# Higher miR222 expression is associated with higher odds of
# belonging to the Hodgkin lymphoma group (OR = 24.12, 95% CI: 3.56-163.28).
# The model achieves an in-sample accuracy of 88.5%,
# sensitivity of 82.6%, specificity of 93.1% and ROC AUC of 0.955.
# These performance measures may be optimistic because the model
# was developed and evaluated on the same observations.
# Results should be interpreted cautiously given the small sample size,
# wide confidence intervals and data-driven variable selection.
# Independent validation is required before clinical application.

# CHECK SAVED FIGURES
# List all figures generated in the project folder.
list.files("figures")