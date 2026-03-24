Regularized Models with Multiple Imputation vs. Complete-Case Analysis

Overview

This repository contains the R implementation for comparing Regularized Logistic Regression models under two different data handling strategies:

    Complete-Case Analysis (CC): Using only observations with no missing values.

    Multiple Imputation (MI): Using Bayesian Bootstrap Predictive Mean Matching (BBPMM) to handle missingness.

The project evaluates the performance of Adaptive LASSO (aLASSO), Adaptive Elastic Net (aENET), and their MI-specific counterparts: GaLASSO and SaENET.

Scientific Context

In clinical datasets, missing data is a common challenge. Simply omitting rows (Complete-Case) can lead to bias and loss of statistical power. This study implements:

    Variable Selection: Identifying clinical predictors for patient intubation.

    Regularization: Penalization techniques to prevent overfitting.

    Pooling: Combining results from m=5 imputed datasets using Rubin's rules or stacked approaches.

Requirements

To run the analysis, you will need the following R packages:
R

install.packages(c("glmnet", "mice", "miselect", "BaBooN", "pROC", "dcurves", "patchwork", "tidyverse"))

How to Run

    Clone the repository:
    Bash

    git clone https://github.com/DiegoSsant/R-code-for-Regularized-Models-with-Multiple-Imputation-versus-Complete-Case-Analysis

    Place your dataset synthetic_data_hc.csv in the root folder.

    Open analysis_script.R in RStudio.

    Run the script to generate:

        Coefficient Tables: Regularized parameter estimates.

        ROC Curves: Comparison of AUC across models.

        DCA (Decision Curve Analysis): Clinical utility evaluation.

Key Results

The analysis generates comparative plots to visualize model performance:
1. ROC Curves

Comparison of predictive power between Complete-Case (aLASSO/aENET) and Multiple Imputation (GaLASSO/SaENET) approaches.
2. Decision Curve Analysis (DCA)

Evaluation of "Net Benefit" to determine if the models are clinically useful for decision-making across different probability thresholds.

File Structure

    CCvsMI-RegularizedModels.R: Main R script with the complete pipeline.

    synthetic_data_hc.csv: (Sample/Synthetic) Clinical data used for the study.

    README.md: Project documentation.

License

This project is licensed under the MIT License - see the LICENSE file for details.
