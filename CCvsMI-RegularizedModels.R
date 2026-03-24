#==============================================================================#
# Title: Regularized Models with Multiple Imputation vs Complete Case Analysis:
#        Assessing Clinical Decision Utility via Decision Curve Analysis
# Author: Diego da Silva Santos
# Context: Analysis of predictive models for medical outcomes (e.g., intubation)
#==============================================================================#

# 1. DEPENDENCIES & SETUP --------------------------------------------------
# Loading necessary libraries for imputation, modeling, and visualization
libs <- c(
  "foreach", "ggplot2", "doParallel", "mice", "miselect", "data.table",
  "stargazer", "caret", "BaBooN", "ROCR", "precrec", "pROC", "glmnet",
  "parallelly", "parallel", "tidyverse", "dcurves", "patchwork", "knitr"
)

# Install missing packages and load all
lapply(libs, function(x) {
  if (!require(x, character.only = TRUE)) install.packages(x)
  library(x, character.only = TRUE)
})

# 2. DATA LOADING & PRE-PROCESSING -----------------------------------------
base_filtrada <- read.csv("dados_sinteticos_hc1.csv", header = TRUE, 
                          sep = ",", na.strings = c("NA",""), dec = ".")

# Handling potential string "NA" and cleaning (Adjusting to use base_filtrada)
for(i in 1:ncol(base_filtrada)){
  base_filtrada[,i] <- ifelse(base_filtrada[,i] == "NA", NA, base_filtrada[,i])
}

# 3. ANALYSIS I: COMPLETE-CASE (aLASSO & aENET) ----------------------------
# Focus: Models using only patients with no missing values

# Omitting the missing data
set.seed(101)
base_comp <- na.omit(base_filtrada)
nrow(base_comp
     )
cat("Original Proportion:\n")
print(prop.table(table(base_filtrada$intubation)))

cat("Proportion Data Omitted:\n")
print(prop.table(table(base_comp$intubation)))

# Creates an index maintaining the proportion of the response variable Y.
indice_teste <- createDataPartition(base_comp$intubation, p = 0.2, list = FALSE)
teste <- base_comp[indice_teste, ]
base_comp_treino <- base_comp[-indice_teste, ]

cat("\Proportion in Training:\n")
print(prop.table(table(base_comp_treino$intubation))) # Deve ser ~igual à original

cat("\Proportion in the Test:\n")
print(prop.table(table(teste$intubation)))

# Adjusting the adaptive LASSO model
start.time <- Sys.time()
fit_lasso <- cv.glmnet(as.matrix(base_comp_treino[,-5]),
                       as.matrix(base_comp_treino[,5]),family = "binomial",
                       alpha = 0)
initial_coefs <- coef(fit_lasso, s = "lambda.min")
initial_coefs_no_intercept <- as.vector(initial_coefs)[-1]
adaptive_weights <- 1 / (abs(initial_coefs_no_intercept) + 1e-10)

cv_adaptive_lasso <- cv.glmnet(
  as.matrix(base_comp_treino[,-5]), 
  as.matrix(base_comp_treino[,5]),
  family = "binomial",
  alpha = 1,                 # alpha = 1 for LASSO
  penalty.factor = adaptive_weights
)
end.time <- Sys.time()
time.taken <- end.time - start.time;time.taken
plot(cv_adaptive_lasso)
final_coefs <- coef(cv_adaptive_lasso, s = "lambda.min")

# adjusting the adaptive ENET model
alpha_grid <- seq(0.1, 0.9, by = 0.1)
best_cvm <- Inf
best_alpha <- NA
best_model <- NULL
start.time <- Sys.time()
for (alpha_val in alpha_grid) {
  
  # Perform cross-validation for the current alpha, with adaptive weights.
  cv_model <- cv.glmnet(
    as.matrix(base_comp_treino[,-5]), 
    as.matrix(base_comp_treino[,5]), 
    alpha = alpha_val,
    penalty.factor = adaptive_weights
  )
  
  # Check if this model is better than the best one found so far.
  if (min(cv_model$cvm) < best_cvm) {
    best_cvm <- min(cv_model$cvm)
    best_alpha <- alpha_val
    best_model <- cv_model
  }
}
end.time <- Sys.time()
time.taken <- end.time - start.time;time.taken

best_lambda <- best_model$lambda.min
final_coefs_enet <- coef(best_model, s = "lambda.min")

# Data for ranking power values
# estimated coefficients aLASSO and aENET
beta <- final_coefs
beta2 <- final_coefs_enet

# test covariates
newx <- base_comp[c(indice_teste),-5]

y.testlasso <- base_comp[c(indice_teste),5]
n_newx1 <- nrow(newx)

intercep1 <- rep(1, times=n_newx1)
beta_1 <- as.matrix(beta)
beta_2 <- as.matrix(beta2)

# Predicted probabilities
n_fitlasso <- predict(cv_adaptive_lasso, newx = as.matrix(newx), 
                      s = "lambda.min", type = "response")

n_fitenet <- predict(best_model, newx = as.matrix(newx), 
                     s = "lambda.min", type = "response")


# Obtaining AUC from models
roc_object_lasso <- roc(y.testlasso,n_fitlasso, auc=T, ci=T); roc_object_lasso
roc_object_enet <- roc(y.testlasso,n_fitenet, auc=T, ci=T); roc_object_enet

# Finding the best cut for the aLASSO
coords_lasso <- coords(roc_object_lasso, "best", best.method = "youden", 
                       ret = c("threshold", "sensitivity", "specificity"))
print(coords_lasso)

# threshold
corte_otimolasso <- coords_lasso$threshold
pred1_otimizadolasso <- ifelse(n_fitlasso >= corte_otimolasso, 1, 0)

metricas_precisas <- ci.coords(
  roc_object_lasso, 
  x = corte_otimolasso,      # threshold 
  input = "threshold",
  ret = c("sensitivity", "specificity", "accuracy", "ppv", "npv"),
  conf.level = 0.95          # Nível de confiança (padrão 95%)
)

print(metricas_precisas)

# Finding the best cut for aENET
coords_enet <- coords(roc_object_enet, "best", best.method = "youden", 
                      ret = c("threshold", "sensitivity", "specificity"))
print(coords_enet)

# threshold
corte_otimo2enet <- coords_enet$threshold
pred2_otimizadoenet <- ifelse(n_fitenet >= corte_otimo2enet, 1, 0)

metricas_precisas1 <- ci.coords(
  roc_object_enet, 
  x = corte_otimo2enet,      # threshold
  input = "threshold",       
  ret = c("sensitivity", "specificity", "accuracy", "ppv", "npv"),
  conf.level = 0.95          # Nível de confiança (padrão 95%)
)

print(metricas_precisas1)

# Auxiliary function to extract data from ROCR into a data frame.
extract_roc <- function(perf_obj, model_name) {
  data.frame(
    FPR = unlist(perf_obj@x.values),
    TPR = unlist(perf_obj@y.values),
    Model = model_name
  )
}

# Predicted values for the ROC object
preds <- prediction(as.numeric(n_fitlasso),as.factor(y.testlasso))
preds2 <- prediction(as.numeric(n_fitenet),as.factor(y.testlasso))
perflasso <- performance(preds, "tpr", "fpr" )
perfenet <- performance(preds2, "tpr", "fpr" )

# 4. ANALYSIS II: MULTIPLE IMPUTATION (GaLASSO & SaENET) --------------------
# Focus: Models using Bayesian Bootstrap Predictive Mean Matching (BBPMM)

# Creates an index maintaining the proportion of the response variable Y.
indice_teste1 <- createDataPartition(base_filtrada$intubation, p = 0.20, list = FALSE)

# Create the Training and Testing Bases
# Note: We use [-test_index, ] for training (80%) and [test_index, ] for testing (20%).
base_analise_t <- base_filtrada[-indice_teste1,] # Base de TREINO (Imputação)
base_teste     <- base_filtrada[indice_teste1,]  # Base de TESTE (Validação)

# confirming if the proportions are similar.

cat("Original Ratio base:\n")
print(prop.table(table(base_filtrada$intubation)))

cat("\nProportion in Training:\n")
print(prop.table(table(base_analise_t$intubation))) # Deve ser ~igual à original

cat("\nProportion in the Test:\n")
print(prop.table(table(base_teste$intubation)))

# Settings to avoid internal parallelism in BaBooN (which can cause conflicts) 
if (requireNamespace("rhpcBLASctl", quietly = TRUE)) {
  rhpcBLASctl::blas_set_num_threads(1)
  rhpcBLASctl::omp_set_num_threads(1)
}


# Imputation by Bayesian Bootstrap Predict Mean Meaning
run_parallel_BBPMM <- function(base_analysis_t,
                               M_total = 5,
                               ncores = 4,
                               BBPMM_args = list(stepmod = "stepAIC",
                                                 maxit.mult = 3,
                                                 maxit.glm = 25,
                                                 maxPerc = 0.98,
                                                 verbose = FALSE,
                                                 setSeed = NULL,
                                                 chainDiagnostics = TRUE)) {
  
  # number of workers
  ncores <- max(1, min(ncores, parallel::detectCores()))
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  # export packages/objects needed to workers
  # If BaBooN exports S3 functions, it's best to load the package in the workers.
  clusterEvalQ(cl, {
    library(BaBooN)           # BaBooN carries every worker
    if (requireNamespace("rhpcBLASctl", quietly = TRUE)) {
      rhpcBLASctl::blas_set_num_threads(1)
      rhpcBLASctl::omp_set_num_threads(1)
    }
    NULL
  })
  
  # Execute M_total executions with M=1 in parallel.
  results <- foreach(i = seq_len(M_total), .packages = "BaBooN", 
                     .export = ls(envir = globalenv())) %dopar% {
                       # define seed by repetition for reproducibility
                       if (!is.null(BBPMM_args$setSeed)) set.seed(BBPMM_args$setSeed + i)
                       
                       # Call: adapt the argument names if necessary.
                       res <- do.call(BBPMM, c(list(base_analysis_t, M = 1, nIter = 2,
                                                    outFile = NULL), BBPMM_args))
                       # returns only the object of imputation (adjustment according to structure)
                       res$impdata   # Assuming the structure you showed: imp$impdata
                     }
  
  stopCluster(cl)
  registerDoSEQ()
  
  # results is a list of length M_total containing the data from each run.
  # Add according to the expected format:
  # e.g., create a list with all the import data
  all_imputations <- results
  
  return(all_imputations)
}

# training data imputation
start.time <- Sys.time()
imp <- run_parallel_BBPMM(base_analise_t, M_total = 5, ncores = 4,
                          BBPMM_args = list(stepmod = "stepAIC",
                                            maxit.mult = 3, maxit.glm = 25,
                                            maxPerc = 0.98, verbose = FALSE,
                                            setSeed = 101))
end.time <- Sys.time()
time.taken <- end.time - start.time;time.taken

conjuntos_imputados <- imp

x <- list()
y <- list()
for (i in 1:5) {
  x[[i]] <- as.matrix(conjuntos_imputados[[i]][,-5])
  y[[i]] <- conjuntos_imputados[[i]]$intubation
}
pf       <- rep(1, 45)
adWeight <- rep(1, 45)
weights <- 1 - rowMeans(is.na(base_analise_t))

# Estimation of MI-LASSO (GaLASSO) logistic parameters
start.time <- Sys.time()
fit1 <- cv.galasso(x, y, pf, adWeight,family = "binomial")
end.time <- Sys.time()
time.taken <- end.time - start.time;time.taken

estimated.coefficientes1 <- coef(fit1, lambda = fit1$lambda.1se)
stargazer(estimated.coefficientes1,title = "Estimated coefficients",
          type="text",style = "aer",digits = 4,flip = TRUE)

# Estimation of MI-aENET (SaENET) logistic parameters
start.time <- Sys.time()
fit2 <- cv.saenet(x, y, pf, adWeight, weights, family = "binomial")
end.time <- Sys.time()
time.taken <- end.time - start.time;time.taken

estimated.coefficientes2 <- coef(fit2,lambda = fit2$lambda.min)
stargazer(estimated.coefficientes2,title = "Estimated coefficients",
          type="text",style = "aer",digits = 4,flip = TRUE)

# BBPMM imputation to the test set
start.time <- Sys.time()
imp1 <- run_parallel_BBPMM(base_teste, M_total = 5, ncores = 4,
                           BBPMM_args = list(stepmod = "stepAIC",
                                             maxit.mult = 3, maxit.glm = 25,
                                             maxPerc = 0.98, verbose = FALSE,
                                             setSeed = 101))
end.time <- Sys.time()
time.taken <- end.time - start.time;time.taken

conjuntos_imputados1 <- imp1

# Test sets 
x.test1 <- as.matrix(conjuntos_imputados1[[1]][,-5])
y.test1 <- conjuntos_imputados1[[1]]$intubation
newx1 = x.test1
newx1 <- data.matrix(newx1)

x.test2 <- as.matrix(conjuntos_imputados1[[2]][,-5])
y.test2 <- conjuntos_imputados1[[2]]$intubation
newx2 = x.test2
newx2 <- data.matrix(newx2)

x.test3 <- as.matrix(conjuntos_imputados1[[3]][,-5])
y.test3 <- conjuntos_imputados1[[3]]$intubation
newx3 = x.test3
newx3 <- data.matrix(newx3)

x.test4 <- as.matrix(conjuntos_imputados1[[4]][,-5])
y.test4 <- conjuntos_imputados1[[4]]$intubation
newx4 = x.test4
newx4 <- data.matrix(newx4)

x.test5 <- as.matrix(conjuntos_imputados1[[5]][,-5])
y.test5 <- conjuntos_imputados1[[5]]$intubation
newx5 = x.test5
newx5 <- data.matrix(newx5)

y.test <- (y.test1+ y.test2+y.test3+y.test4+y.test5)/5
x.test <- (x.test1+ x.test2+x.test3+x.test4+x.test5)/5
x.test <- data.frame(x.test)

# List of test matrices (each imputation)
newx_list <- list(newx1, newx2, newx3, newx4, newx5)
intubation <- data.frame(y.test)

#------------------------------------------------------------------------------#
# GaLASSO
# Prepare aligned coefficients.
coef_list <- lapply(estimated.coefficientes1, function(v) {
  vnum <- as.numeric(v)
  names(vnum) <- names(v)
  vnum
})
all_names <- sort(unique(unlist(lapply(coef_list, names))))
coef_list_aligned <- lapply(coef_list, function(v) {
  vt <- numeric(length(all_names))
  names(vt) <- all_names
  vt[names(v)] <- v
  vt
})

# Auxiliary function: calculates probabilities for a given matrix newx
calc_probs <- function(newx_mat, coef_list) {
  beta_names <- colnames(newx_mat)
  m <- length(coef_list)
  probs <- matrix(NA, nrow = nrow(newx_mat), ncol = m)
  for (i in seq_len(m)) {
    coef_i <- coef_list[[i]]
    intercept_i <- if ("(Intercept)" %in% names(coef_i)) coef_i["(Intercept)"] else 0
    # betas aligned with the columns of newx_mat
    beta_i <- numeric(length(beta_names))
    names(beta_i) <- beta_names
    beta_i[intersect(names(coef_i), beta_names)] <- coef_i[intersect(names(coef_i), beta_names)]
    eta_i <- intercept_i + as.numeric(newx_mat %*% beta_i)
    probs[, i] <- 1 / (1 + exp(-eta_i))
  }
  # average per row (pooling within each imputation)
  rowMeans(probs)
}

# Calculate probs_pooled for each imputation.
pooled_list <- lapply(newx_list, calc_probs, coef_list = coef_list)

# Calculate the final average between imputations.
probs_final <- Reduce("+", pooled_list) / length(pooled_list)

# Standard deviation between imputations
probs_mat_all <- do.call(cbind, pooled_list)
probs_final_sd <- apply(probs_mat_all, 1, sd)

# Final result
head(probs_final)     # final average probabilities
head(probs_final_sd)  # desvio padrão entre imputações

# Auxiliary function to calculate problems for a SaENET stacked imputation
calc_probs_saenet <- function(newx_mat, beta) {
  intercep <- rep(1, nrow(newx_mat))
  eta <- cbind(intercep, newx_mat) %*% beta
  as.vector(1 / (1 + exp(-eta)))
}

# Calculate probabilities for each imputation.
probs2_list <- lapply(newx_list, calc_probs_saenet, beta = beta)

# Average between imputations (final pooling)
probs2_final <- Reduce("+", probs2_list) / length(probs2_list)

# Standard deviation between imputations
probs2_mat <- do.call(cbind, probs2_list)
probs2_final_sd <- apply(probs2_mat, 1, sd)

# Final result
head(probs2_final)     # final average probability SaENET

# Obtaining AUCs from the models
roc_object_galasso <- roc(y.test,probs_final, auc=T, ci=T);roc_object_galasso
roc_object_saenet <- roc(y.test,probs2_final, auc=T, ci=T);roc_object_saenet

# DeLong test to compare roc1 and roc2
test_result <- roc.test(roc_object_galasso, roc_object_saenet, method="delong")
print(test_result)

# Finding the best cut for the GALASSO
coords_galasso <- coords(roc_object_galasso, "best", best.method = "youden", 
                         ret = c("threshold", "sensitivity", "specificity"))
print(coords_galasso)

# threshold
corte_otimo <- coords_galasso$threshold
pred1_otimizado <- ifelse(probs_final >= corte_otimo, 1, 0)

metricas_precisasga <- ci.coords(
  roc_object_galasso, 
  x = corte_otimo,      # threshold
  input = "threshold",  
  ret = c("sensitivity", "specificity", "accuracy", "ppv", "npv"), 
  conf.level = 0.95          # Confidence level (standard 95%)
)

print(metricas_precisasga)

# Finding the best cut for SaENET
coords_saenet <- coords(roc_object_saenet, "best", best.method = "youden", 
                        ret = c("threshold", "sensitivity", "specificity"))
print(coords_saenet)

# threshold
corte_otimo2 <- coords_saenet$threshold
pred2_otimizado <- ifelse(probs2_final >= corte_otimo2, 1, 0)

metricas_precisassa <- ci.coords(
  roc_object_saenet, 
  x = corte_otimo2,      # threshold
  input = "threshold",  
  ret = c("sensitivity", "specificity", "accuracy", "ppv", "npv"),
  conf.level = 0.95          # Confidence level (standard 95%)
)

print(metricas_precisassa)

# Performance objects for ROC curves
pred1s <- prediction(as.numeric(probs_final), as.factor(y.test))
pred2s <- prediction(as.numeric(probs2_final), as.factor(y.test))
perf1 <- performance(pred1s, "tpr", "fpr" )
perf2 <- performance(pred2s, "tpr", "fpr")

# 5. MODEL EVALUATION & VISUALIZATION --------------------------------------

# 5.1 ROC curves 
# Extracting to Graph (a)
df_cc <- rbind(
  extract_roc(perflasso, "aLASSO"),
  extract_roc(perfenet, "aENET")
)

# Extracting to Graph (b)
df_mi <- rbind(
  extract_roc(perf1, "GaLASSO"),
  extract_roc(perf2, "SaENET")
)

# Basic aesthetics for ROC graphs
plot_roc_theme <- function(df, title, x_label = "False Positive Rate", 
                           y_label = "True Positive Rate") {
  ggplot(df, aes(x = FPR, y = TPR, color = Model)) +
    geom_line(linewidth = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + 
    # Defining fixed colors to ensure consistency
    scale_color_manual(values = c("aLASSO" = "black", "aENET" = "red", 
                                  "GaLASSO" = "black", "SaENET" = "red")) +
    labs(title = title, x = x_label, y = y_label) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.spacing.x = unit(0.8, "cm"), 
      legend.text = element_text(margin = margin(r = 10)),
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 10)
    ) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE))
}

p1_roc <- plot_roc_theme(df_cc, "(a) ROC Curves: Complete-Case")
p2_roc <- plot_roc_theme(df_mi, "(b) ROC Curves: Multiple Imputation")

combined_roc <- p1_roc + p2_roc + 
  plot_layout(ncol = 2) & # Ensures side-by-side graphics
  theme(legend.position = "bottom")

# View the result for ROC curves
combined_roc


# 5.2 DCA curves
# Limiting probabilities
n_fitlasso <- pmin(pmax(n_fitlasso, 0), 1)
n_fitenet  <- pmin(pmax(n_fitenet, 0), 1)

# Create a clean data frame for the DCA function.
datapred1 <- data.frame(newx,
                        y.testlasso = as.numeric(y.testlasso),# Response variable
                        aLASSO = as.numeric(n_fitlasso),# LASSO Probabilities
                        aENET = as.numeric(n_fitenet)   # Probabilities aENET
)

# Create the clean MI data frame for the DCA function.
x.test$GaLASSO = as.numeric(probs_final)
x.test$SaENET = as.numeric(probs2_final)
GaLASSO <- data.frame(probs_final)
SaENET <- data.frame(probs2_final)
datapred <- data.frame(x.test,
                       as.numeric(y.test),
                       "GaLASSO"=GaLASSO, 
                       "SaENET"=SaENET)

# Prepare the data 
df_plot_cc <- dcurves::dca(y.testlasso ~ aLASSO + aENET, data = datapred1) %>% 
  as_tibble() %>%
  mutate(model_type = "Complete-Case Analysis")

df_plot_mi <- dcurves::dca(y.test ~ GaLASSO + SaENET, data = datapred) %>% 
  as_tibble() %>%
  mutate(model_type = "Multiple Imputation")

# Internal function to standardize aesthetics.
plot_dca_custom <- function(df, title) {
  ggplot(df, aes(x = threshold, y = net_benefit, color = label)) +
    stat_smooth(method = "loess", se = FALSE, span = 0.2, linewidth = 0.5) +
    coord_cartesian(ylim = c(-0.05, max(df$net_benefit, na.rm = TRUE) + 0.1)) +
    #scale_color_brewer(palette = "Set1") + # Cores distintas e sóbrias
    scale_color_manual(values = c("aLASSO" = "black", "aENET" = "red", 
                                  "GaLASSO" = "black", "SaENET" = "red", 
                                  "Treat All"="darkslategray4", "Treat None"="gray")) +
    labs(
      title = title,
      x = "Threshold Probability",
      y = "Net Benefit",
      color = "Strategy"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.spacing.x = unit(0.5, "cm"),   
      legend.text = element_text(margin = margin(r = 10)),  
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 10))+
    guides(color = guide_legend(nrow = 1, byrow = TRUE))
}

# Generate the two graphs individually.
p1 <- plot_dca_custom(df_plot_cc, "(a) DCA Curves: Complete-Case")
p2 <- plot_dca_custom(df_plot_mi, "(b) DCA Curves: Multiple Imputation")

# Combine using patchwork
combined_plot <- p1 + p2 + 
  #plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# View the result for DCA curves
combined_plot