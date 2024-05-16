library(readr)
library(readxl)
library(glmnet)
library(caret)

source("scripts/support_glmnet.R")

GDSC1_fitted_dose_response_24Jul22 <- as.data.frame(read_excel("data/GDSC1_fitted_dose_response_24Jul22.xlsx"))

rnaseq_tpm_20220624_logp1 <- as.data.frame(read_csv("data/rnaseq_tpm_20220624_logp1.csv"))
row.names(rnaseq_tpm_20220624_logp1) <- rnaseq_tpm_20220624_logp1$GENE_SYMBOLS

expr_matrix <- as.data.frame(t(rnaseq_tpm_20220624_logp1[,-1]))

LN_IC50s <- GDSC1_fitted_dose_response_24Jul22[GDSC1_fitted_dose_response_24Jul22$DRUG_NAME == "Pelitinib",c(6,16)]
row.names(LN_IC50s) <- LN_IC50s$SANGER_MODEL_ID
LN_IC50s <- LN_IC50s[,-1,drop = F]

expr_matrix <- merge(LN_IC50s,expr_matrix, by = "row.names")
row.names(expr_matrix) <- expr_matrix$Row.names
expr_matrix <- expr_matrix[,-1]

#ML happens here
set.seed(42)
X = model.matrix(LN_IC50 ~ ., expr_matrix)[, -1]
y = expr_matrix$LN_IC50


fit_ridge_cv = cv.glmnet(X, y, alpha = 0, nlambda = 100)
plot(fit_ridge_cv)
sqrt(fit_ridge_cv$cvm[fit_ridge_cv$lambda == fit_ridge_cv$lambda.1se])
coefs_ridge_cv <- get_coef_dataframe(fit_ridge_cv)

fit_to_plot <- fit_ridge_cv
data_to_model <- expr_matrix
best_lambda <- fit_to_plot$lambda.1se
response <- data_to_model$LN_IC50
predictions <- predict(fit_to_plot, s = best_lambda, newx = as.matrix(data_to_model[,-which(names(data_to_model) == "LN_IC50")]))
obs_vs_pred <- data.frame(Observed = response, Predicted = predictions)
names(obs_vs_pred) <- c("Observed","Predicted")

ggplot(obs_vs_pred, aes(x = Observed, y = Predicted)) +
  geom_point() +
  geom_smooth(method = "lm", col = "red") +
  ggtitle("Predicted vs Observed Values") +
  xlab("Observed Values") +
  ylab("Predicted Values") + theme_minimal()

## now with TFs
GDSC_TF_activities <- as.data.frame(read_csv("results/GDSC_TF_activities.csv"))
row.names(GDSC_TF_activities) <- GDSC_TF_activities$source

TF_matrix <- as.data.frame(t(GDSC_TF_activities[,-1]))

TF_matrix <- merge(LN_IC50s,TF_matrix, by = "row.names")
row.names(TF_matrix) <- TF_matrix$Row.names
TF_matrix <- TF_matrix[,-1]

#ML happens here
set.seed(42)
X = model.matrix(LN_IC50 ~ ., TF_matrix)[, -1]
y = TF_matrix$LN_IC50


fit_ridge_TF_cv = cv.glmnet(X, y, alpha = 0, nlambda = 100)
plot(fit_ridge_TF_cv)
sqrt(fit_ridge_TF_cv$cvm[fit_ridge_TF_cv$lambda == fit_ridge_TF_cv$lambda.1se])
coefs_ridge_TF_cv <- get_coef_dataframe(fit_ridge_TF_cv)

fit_to_plot <- fit_ridge_TF_cv
data_to_model <- TF_matrix
best_lambda <- fit_to_plot$lambda.1se
response <- data_to_model$LN_IC50
predictions <- predict(fit_to_plot, s = best_lambda, newx = as.matrix(data_to_model[,-which(names(data_to_model) == "LN_IC50")]))
obs_vs_pred <- data.frame(Observed = response, Predicted = predictions)
names(obs_vs_pred) <- c("Observed","Predicted")

ggplot(obs_vs_pred, aes(x = Observed, y = Predicted)) +
  geom_point() +
  geom_smooth(method = "lm", col = "red") +
  ggtitle("Predicted vs Observed Values") +
  xlab("Observed Values") +
  ylab("Predicted Values") + theme_minimal()
