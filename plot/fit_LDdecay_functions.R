###############################################################################
# The LD decay functions in this file are copied from the fit_LDdecay.R script
# in ngsLD: https://github.com/fgvieira/ngsLD/blob/master/scripts/fit_LDdecay.R
#
# The original script was written by Filipe G. Vieira and Emma Fox:
## # FileName: fit_LDdecay.R v1.0.9
##  Author: "Filipe G. Vieira (fgarrettvieira _at_ gmail [dot] com)"
##  Author: "Emma Fox (e.fox16 _at_ imperial [dot] ac [dot] uk)"
###############################################################################


### Fit decay
# Model function
ld_exp <- function(par, dist, ld_stat) {
    par <- as.numeric(par)
    if (ld_stat == "r2" || ld_stat == "r2pear") {
        C <- par[1] * dist
        r2h <- par[2]
        r2l <- par[3]
        if (opt$n_ind) {
            # LD decay curve adjusted for finite samples sizes
            ((10 + C) / ((2 + C) * (11 + C))) * (1 + ((3 + C) * (12 + 12 * C + C^2)) / (opt$n_ind * (2 + C) * (11 + C)))
        } else {
            # Theoretical expectation under to drift
            # 1 / (1 + C)
            # Theoretical expectation with r2_high and r2_low
            (r2h - r2l) / (1 + C) + r2l
        }
    } else if (ld_stat == "Dp") {
        D0 <- 1
        t <- par[1]
        Dh <- par[2]
        Dl <- par[3]
        Dl + (Dh - Dl) * D0 * (1 - dist * opt$recomb_rate / 1e6)^t
    }
}


# Evaluation function
fit_eval <- function(par, obs_data) {
    if (length(unique(obs_data$LD)) != 1) {
        stop("Invalid data.frame (several LD measures)", call. = opt$debug)
    }
    model <- ld_exp(par, obs_data$Dist, obs_data$LD[1])
    eval <- sum((model - obs_data$value)^2)
    if (opt$debug) {
        str(par)
        print(eval)
    }
    eval
}


# Fitting function
fit_func <- function(x, fit_level) {
    ld_stat <- x$LD[1]
    # There is no fitting model for D
    if (ld_stat == "D") {
        return(NULL)
    }

    optim_tmp <- list()
    n_iter <- ifelse(fit_level >= 10, fit_level, 1)
    for (iter in 1:n_iter) {
        # Fit LD model
        init_vals <- runif(3)
        if (ld_stat == "Dp") {
            init_vals[1] <- runif(1, 10, 20)
            par <- data.frame(init = init_vals, low_lim = c(0, 0, 0), up_lim = c(Inf, 1, 1))
        } else { # r2 and r2pear
            init_vals[1] <- runif(1, 0, 0.1)
            par <- data.frame(init = init_vals, low_lim = c(0, 0, 0), up_lim = c(1, 1, 1))
        }
        if (opt$debug) str(par)

        optim_tmp <- append(optim_tmp, list("BFGS" = optim(par$init, fit_eval, obs_data = x, method = "BFGS")))
        if (fit_level > 1) optim_tmp <- append(optim_tmp, list("Nelder-Mead" = optim(par$init, fit_eval, obs_data = x, method = "Nelder-Mead")))
        if (fit_level > 2) optim_tmp <- append(optim_tmp, list("L-BFGS-B" = optim(par$init, fit_eval, obs_data = x, method = "L-BFGS-B", lower = par$low_lim, upper = par$up_lim)))
    }
    if (opt$debug) str(optim_tmp)

    # If not using the theoretical r2 decay curve (with r2h and r2l)
    if (opt$n_ind > 0) {
        if (ld_stat != "Dp") {
            optim_tmp <- lapply(optim_tmp, function(x) {
                x$par[2] <- x$par[3] <- 0
                x
            })
        }
    }

    # Filter out runs that not-converged and/or with out-of-bound parameters
    # Columns stand for: par1 (rate), par2 (high), par3 (low), score, counts, convergence, message
    optim_tmp <- Filter(function(x) {
        x$convergence == 0 &
            x$par[1] >= par$low_lim[1] & x$par[1] <= par$up_lim[1] &
            x$par[2] >= par$low_lim[2] & x$par[2] <= par$up_lim[2] &
            x$par[3] >= par$low_lim[3] & x$par[3] <= par$up_lim[3] &
            x$par[2] >= x$par[3]
    }, optim_tmp)

    # Pick best run
    optim_tmp <- optim_tmp[order(sapply(optim_tmp, "[[", 2))[1]]
    if (length(optim_tmp[[1]]) == 0) stop("convergence analyses failed. Please try increasing the fit level (--fit_level)", call. = opt$debug)
    if (opt$debug) cat("Best fit for ", as.character(x$File[1]), " (", as.vector(ld_stat), "): ", names(optim_tmp), sep = "", fill = TRUE)
    optim_tmp[[1]]$par
}



read_ld_files <- function(n_files, ld_files, col, ld, max_kb_dist, use_recomb_rate, recomb_rate, fit_bin_size, bin_quant, header=c("Dist", "r2pear", "D", "Dp", "r2")){
    ld_data <- data.frame()
    for (i in 1:n_files) {
        ld_file <- as.character(ld_files$File[i])
        # Read point data
        tmp_data <- read.table(gzfile(ld_file), sep = "\t", quote = "\"", dec = ".")[-(1:(col - 1))]
        # Check if file is valid
        if (ncol(tmp_data) < 5) {
            stop("Invalid LD file format.\n")
        }
        # Add column labels
        colnames(tmp_data) <- header
        # Extract relevant columns
        tmp_data <- tmp_data[, which(names(tmp_data) %in% c("Dist", ld))]
        # Filter by minimum distance
        tmp_data <- tmp_data[which(tmp_data$Dist < max_kb_dist * 1000), ]
        # Convert all 'Inf' to NA
        tmp_data[mapply(is.infinite, tmp_data)] <- NA
        # Calculate genetic distances, according to Haldane's formula (assumes constant rate across all dataset)
        if (use_recomb_rate && !is.null(recomb_rate)) {
            tmp_data$Dist <- (1 - (1 - recomb_rate * 0.01 / 1e6)^(tmp_data$Dist)) / 2
        }
        # Bin data
        if (fit_bin_size > 1) {
            breaks <- seq(0, max(tmp_data$Dist) + fit_bin_size, fit_bin_size)
            tmp_data$Dist <- cut(tmp_data$Dist, breaks, head(breaks, -1))
            if (bin_quant > 0) {
                tmp_data <- aggregate(. ~ Dist, data = tmp_data, quantile, probs = bin_quant / 100, na.rm = TRUE)
            } else {
                tmp_data <- aggregate(. ~ Dist, data = tmp_data, mean, na.rm = TRUE)
            }
        }
        tmp_data$File <- ld_file
        ld_data <- rbind(ld_data, reshape2::melt(tmp_data, c("File", "Dist"), variable.name = "LD", na.rm = TRUE))
    }
    # Remove factor in Dist
    ld_data$Dist <- as.numeric(levels(ld_data$Dist))[ld_data$Dist]

    return(ld_data)
}

fit_ld_decay_dist <- function(fit_level, plot_x_lim, plot_line_smooth, fit_boot) {
    if (fit_level > 0) {
        # Define line "resolution"
        smooth <- seq(1, plot_x_lim, length = plot_line_smooth)
        # Full data
        optim_fit <- ddply(ld_data, .(File, LD), fit_func, fit_level = fit_level)

        # Bootstrap
        if (fit_boot > 0) {
            boot_rep_fit <- c()
            for (b in 1:fit_boot) {
                ld_bootdata <- ddply(ld_data, .(File, LD), function(x) {
                    x[sample(nrow(x), size = nrow(x), replace = TRUE), ]
                })
                boot_rep_fit <- rbind(boot_rep_fit, data.frame(ddply(ld_bootdata, .(File, LD), fit_func, fit_level = fit_level), Rep = b))
            }
            boot_fit <- as.data.frame(as.matrix(aggregate(cbind(V1, V2, V3) ~ File + LD, boot_rep_fit, quantile, probs = c(0.025, 0.975), names = FALSE)), stringsAsFactors = FALSE)
            optim_fit <- merge(optim_fit, boot_fit, sort = FALSE)
            optim_data <- ddply(optim_fit, .(File, LD), function(x) data.frame(LD = x[, "LD"], Dist = smooth, value = ld_exp(x[, c("V1", "V2", "V3")], smooth, x[, "LD"]), ci_l = ld_exp(x[, c("V1.2", "V2.1", "V3.1")], smooth, x[, "LD"]), ci_u = ld_exp(x[, c("V1.1", "V2.2", "V3.2")], smooth, x[, "LD"])))
        } else {
            optim_data <- ddply(optim_fit, .(File, LD), function(x) data.frame(LD = x[, "LD"], Dist = smooth, value = ld_exp(x[, c("V1", "V2", "V3")], smooth, x[, "LD"])))
        }
        # Print best FIT parameters
        optim_fit <- plyr::rename(optim_fit, c("V1" = "DecayRate", "V2" = "LDmax", "V3" = "LDmin", "V1.1" = "DecayRate_CI.u", "V1.2" = "DecayRate_CI.l", "V2.1" = "LDmax_CI.l", "V2.2" = "LDmax_CI.u", "V3.1" = "LDmin_CI.l", "V3.2" = "LDmin_CI.u"), warn_missing = FALSE)
        print(optim_fit)

        # Merge data together with extra info from input
        fit_data <- merge(ld_files, optim_data, sort = FALSE)
        return(fit_data)
    }else{
        return(NULL)
    }
}
