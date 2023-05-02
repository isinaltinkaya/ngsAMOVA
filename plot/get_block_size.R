#!/bin/env Rscript



###############################################################################
# get_block_size.R v0.0.0
#'
#' @author Isin Altinkaya, \email{isinaltinkaya@@gmail.com}
#'
#' @description This script estimates the block size that can be used to assume
#' independence between SNPs using the LD decay fit from ngsLD.
#'
#' @references \url{https://github.com/fgvieira/ngsLD/blob/master/scripts/fit_LDdecay.R}
#'
#' @usage Rscript get_block_size.R --ld_files ld_file_list --fit_level 3 --seed 4 --decay_threshold 99
#'
#' @note Current limitations: Does not work with multiple files
###############################################################################



# ______________________________________________________________________________
# LOAD LIBRARIES
library(optparse)
library(tools)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)


source("fit_LDdecay_functions.R")
# ______________________________________________________________________________


# ______________________________________________________________________________
# READ ARGUMENTS

option_list <- list(
  make_option(c("--ld_files"), action = "store", type = "character", default = NULL, help = "File with list of LD files to fit and plot (if ommited, can be read from STDIN)"),
  make_option(c("--header"), action = "store_true", type = "logical", default = FALSE, help = "Input file has header"),
  make_option(c("--col"), action = "store", type = "numeric", default = 3, help = "Which column is distance between sites? [%default]"),
  make_option(c("--ld"), action = "store", type = "character", default = "r2", help = "Which LD stats to plot (r2pear, D, Dp, r2) [%default]"),
  make_option(c("--n_ind"), action = "store", type = "numeric", default = 0, help = "Number of individuals (for 1-parameter r^2 fitting correction)?"),
  make_option(c("-r", "--use_recomb_rate"), action = "store_true", type = "logical", default = FALSE, help = "Assume constant recombination rate. [%default]"),
  make_option(c("--recomb_rate"), action = "store", type = "numeric", default = 1, help = "Recombination rate (or probability of recombination between adjacent sites in cM/Mb) to calculate genetic distances from physical distances. It is assumed to be constant throughout the whole dataset and, for human datasets, a common rule-of-thumb value is 1cM/Mb (1e-6). [%default]"),
  make_option(c("--max_kb_dist"), action = "store", type = "numeric", default = Inf, help = "Maximum distance between SNPs (in kb) to include in the fitting analysis. [%default]"),
  make_option(c("--fit_boot"), action = "store", type = "numeric", default = 0, help = "Number of bootstrap replicates for fitting CI. [%default]"),
  make_option(c("--fit_bin_size"), action = "store", type = "numeric", default = 250, help = "Bin data into fixed-sized windows for fitting. [default %default bps]"),
  make_option(c("--fit_level"), action = "store", type = "numeric", default = 1, help = "Fitting level 0) no fitting, best of 1) Nelder-Mead, 2) and BFGS, 3) and L-BFGS-B). [%default]"),
  make_option(c("--plot_group"), action = "store", type = "character", default = "File", help = "Group variable"),
  make_option(c("--plot_bin_size"), action = "store", type = "numeric", default = 0, help = "Bin data into fixed-sized windows for plotting. [default %default bps]"),
  make_option(c("--plot_x_lim"), action = "store", type = "numeric", default = NULL, help = "X-axis plot limit (in kb). [%default]"),
  make_option(c("--plot_y_lim"), action = "store", type = "numeric", default = NULL, help = "Y-axis plot limit. [%default]"),
  make_option(c("--plot_axis_scales"), action = "store", type = "character", default = "fixed", help = "Plot axis scales: fixed (default), free, free_x or free_y"),
  make_option(c("--plot_size"), action = "store", type = "character", default = "1,2", help = "Plot size (height,width). [%default]"),
  make_option(c("--plot_scale"), action = "store", type = "numeric", default = 1.5, help = "Plot scale. [%default]"),
  make_option(c("--plot_wrap"), action = "store", type = "numeric", default = 0, help = "Plot in WRAP with X columns (default in GRID)"),
  make_option(c("--plot_no_legend"), action = "store_true", type = "logical", default = FALSE, help = "Remove legend from plot"),
  make_option(c("--plot_shapes"), action = "store_true", type = "logical", default = FALSE, help = "Use also shapes (apart from colors)"),
  make_option(c("--plot_line_smooth"), action = "store", type = "numeric", default = 1000, help = "LD decay curve smoothness"),
  make_option(c("--bin_quant"), action = "store", type = "numeric", default = 0, help = "Quantile to represent the bins (e.g. 0 = mean, 50 = median). [%default]"),
  make_option(c("-f", "--plot_wrap_formula"), action = "store", type = "character", default = NULL, help = "Plot formula for WRAP. [%default]"),
  make_option(c("-o", "--out"), action = "store", type = "character", default = NULL, help = "Output file"),
  make_option(c("--seed"), action = "store", type = "numeric", default = NULL, help = "Seed for random number generator"),
  make_option(c("--decay_threshold"), action = "store", type = "numeric", default = "0.99", help = "Specify the percentage decay threshold"),
  make_option(c("--debug"), action = "store_true", type = "logical", default = FALSE, help = "Debug mode. Extra output...")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ______________________________________________________________________________


# Set random seed
if (is.null(opt$seed)) {
  opt$seed <- as.integer(runif(1, 1, 100000))
}
cat("Random seed:", opt$seed, fill = TRUE)
set.seed(opt$seed)



ld_files <- read.table(opt$ld_files, header = FALSE, stringsAsFactors = FALSE)
colnames(ld_files)[1] <- "File"
# Keep file list ordered
for (id in colnames(ld_files)) {
  ld_files[, id] <- factor(ld_files[, id], unique(ld_files[, id]), ordered = TRUE)
}
n_files <- nrow(ld_files)
# Parse input LD files
if (is.null(opt$ld_files)) {
  opt$ld_files <- file("stdin")
}
if (opt$debug) {
  print(ld_files)
}
# Parse LD stats to plot
if (!is.null(opt$ld)) {
  opt$ld <- unlist(strsplit(opt$ld, ","))
}
if (!all(opt$ld %in% c("r2pear", "D", "Dp", "r2"))) {
  stop("Invalid LD measure to plot", call. = opt$debug)
}
n_ld <- length(opt$ld)

# Check if number of individuals was specified
if (opt$n_ind < 0) {
  stop("Number of individuals must be greater than zero", call. = opt$debug)
}
if (opt$n_ind > 0 && length(unique(ld_files$File)) > 1) {
  warning("Sample size for correction will be assumed equal for all samples!", call. = opt$debug)
}
if (!any(opt$ld %in% c("r2", "r2pear")) && opt$n_ind) {
  stop("Number of individuals is only used for r^2 fitting", call. = opt$debug)
}
for (i in opt$ld) {
  cat("==> Fitting", i, "LD decay assuming a", ifelse(i == "Dp" || opt$n_ind == 0, "three (rate of decay, max LD and min LD)", "one (rate of decay)"), "parameter decay model", fill = TRUE)
}

# Check max_kb_dist parameter
if (opt$max_kb_dist < 50) {
  warning("Fitting of LD decay is highly unreliable at short distances (<50kb).", call. = opt$debug)
}
# Check binning sizes
if (opt$plot_bin_size > 0 && opt$plot_bin_size < opt$fit_bin_size) {
  stop("Ploting bin size must be greater than fiting bin size!", call. = opt$debug)
}
if (length(opt$plot_size) < 2) {
  opt$plot_size <- c(opt$plot_size, opt$plot_size)
}
# Plot formula
if (!is.null(opt$plot_wrap_formula)) {
  opt$plot_wrap_formula <- as.formula(opt$plot_wrap_formula)
}
# Set output file name (if not defined)
if (is.null(opt$out)) {
  if (is.null(opt$ld_files)) {
    stop("Output file name required, when reading LD files from STDIN")
  }
  # opt$out <- paste(basename(file_path_sans_ext(opt$ld_files)), ".pdf", sep = "")
  opt$out <- basename(file_path_sans_ext(opt$ld_files))
}



ld_data <- read_ld_files(n_files, ld_files, opt$col, opt$ld, opt$max_kb_dist, opt$use_recomb_rate, opt$recomb_rate, opt$fit_bin_size, opt$bin_quant)


if (is.null(opt$plot_x_lim)) {
  opt$plot_x_lim <- max(ld_data$Dist)
} else {
  opt$plot_x_lim <- opt$plot_x_lim * 1000
}
# Set maximum Y-axis
if (!is.null(opt$plot_y_lim)) {
  opt$plot_y_lim <- c(0, opt$plot_y_lim)
}


# Add extra info
ld_data <- merge(ld_files, ld_data, by = "File", sort = FALSE)
n_groups <- length(unique(ld_data[, opt$plot_group]))
n_plots <- n_files * n_ld / n_groups












fit_data <- fit_ld_decay_dist(opt$fit_level, opt$plot_x_lim, opt$plot_line_smooth, opt$fit_boot)

ggplot() +
  ylab("Linkage Disequilibrium") +
  xlab("Distance") +
  geom_point(data = ld_data, aes_string(x = "Dist", y = "value"), size = 0.05, alpha = 0.2) +
  geom_line(data = fit_data, aes_string(x = "Dist", y = "value")) +
  theme_bw()










# ______________________________________________________________________________
# PLOT
cat("==> Plotting data...", fill = TRUE)


# Add LD decay best fit
# Select fields that are variable
header <- names(which(lapply(lapply(fit_data, unique), length) > 1))
# Exclude non-relevant fields
grp <- header[!header %in% unique(c(as.character(opt$plot_wrap_formula), opt$plot_group, "Dist", "value", "File", "ci_l", "ci_u"))]
# Define line type
if (length(grp) == 0) grp <- NULL
if (length(grp) > 1) stop("invalid number of linetype groups!")



ggplot() +
  theme(panel.spacing = unit(1, "lines")) +
  coord_cartesian(xlim = c(0, opt$plot_x_lim), ylim = opt$plot_y_lim) +
  ylab("Linkage Disequilibrium") +
  xlab("Distance") +
  geom_point(data = ld_data, aes_string(x = "Dist", y = "value"), size = 0.05, alpha = 0.2) +
  geom_line(data = fit_data, aes_string(x = "Dist", y = "value")) +
  theme_bw()

out_file_id <- paste0(opt$out, "_ldDecay_fit")
ggsave(filename = paste0(out_file_id, ".png"), width = 6, height = 4, units = "in", dpi = 300)



first_value <- fit_data[which(min(fit_data$Dist) == fit_data$Dist), ]$value
last_value <- fit_data[which(max(fit_data$Dist) == fit_data$Dist), ]$value
total_decay <- first_value - last_value
fit_data$ld_pct_decay <- 100 - (((fit_data$value - last_value) / total_decay) * 100)


decay_thres <- opt$decay_threshold
thres_point <- fit_data[fit_data$ld_pct_decay > decay_thres, ][1, ]
ggplot() +
  geom_line(data = fit_data, aes(x = Dist, y = ld_pct_decay)) +
  geom_vline(xintercept = thres_point$Dist, linetype = "dashed", color = "red") +
  geom_hline(yintercept = thres_point$ld_pct_decay, linetype = "dashed", color = "red") +
  geom_text(aes(x = thres_point$Dist, y = 100, label = paste0(decay_thres, "%\n", round(thres_point$Dist))), hjust = 1, vjust = 5) +
  theme_bw() +
  labs(x = "Distance", y = "LD decay (%)")


out_file_id <- paste0(opt$out, "_ldDecay_thres", decay_thres)
cat(thres_point$Dist, file = paste0(out_file_id, ".txt"))

ggsave(filename = paste0(out_file_id, ".png"), width = 6, height = 4, units = "in", dpi = 300)
