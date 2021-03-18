# Please update your R to latest version (It is required by bioconductor)
# Following codes are for synergyfinder installation.

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("synergyfinder")

# Please make sure you installed synergyfinder version 2.4.9

library(synergyfinder)
library(openxlsx)

# Two drug combination -----------------------------------------------
# I sent one example dataset "ReplExample.xlsx" together with this file. It has
# 3 replicates for each tested concentration combination points.

# Load data
data <- read.xlsx("/path/to/ReplExample.xlsx") # Example data with replicates

# Reshape data
res <- ReshapeData(data, data_type = "inhibition")

# Plot monotherapy dose-response curves and dose response matrix
PlotDoseResponse(res, statistic = "ci", block_ids = 1)
Plot2DrugSurface(res, plot_block = 1)

# Calculate and visualize synergy scores
res <- CalculateSynergy(res, method = c("ZIP"))
# 2D heatmap with 95% confidence intervals
Plot2DrugHeatmap(res,plot_value = "ZIP_synergy", plot_block = 1, statistic = "ci") 
# 2D heatmap with 95% confidence intervals and mean ZIP synergy for whole matrix
Plot2DrugHeatmap(res, plot_value = "ZIP_synergy", plot_block = 1,
                 statistic = "ci", summary_statistic = "mean")
# 3D landscape
Plot2DrugSurface(res, plot_value = "ZIP_synergy", plot_block = 1)

# Calculate CSS score
res <- CalculateSensitivity(res)
res$drug_pairs # The summarized sensitivity scores
res$sensitivity_scores_statistics # The statistics for sensitivity scores

# S-S (sensitivity-synergy) plot for all blocks in the dataset
PlotSensitiveSynergy(res)

# Barometer for selected concentration combination
PlotBarometer(res, plot_block = 1, plot_concs = c(100,10))

# The analysis work flow for combination without replicates is same.
# You could use build-in example data for testing
# data("mathews_screening_data")
# data <- mathews_screening_data

# Multiple drug combination -----------------------------------------------

# All the functions above are available to multiple drug combination. If you
# want to use Plot2DrugSurface, or Plot2DrugHearmap, please specify which two
# drugs in the combination you want to visualize by setting parameter "drugs"

# Load data
data("NCATS_10023_data") # Build-in 3 drug combination data without replicate.
res <- ReshapeData(NCATS_10023_data)
res <- CaluclateSynergy(res)

# Visualize dose response and synergy score for multiple drug combination.
PlotMultiDrugBar(res, plot_block = 1, plot_value = c("ZIP_synergy", "response"),
                 sort_by = "ZIP_synergy")

# Visualize synergy score surface in 2D concentraction space.
PlotMultiDrugSurface(res, plot_block = 1, plot_value = c("ZIP_synergy"))
