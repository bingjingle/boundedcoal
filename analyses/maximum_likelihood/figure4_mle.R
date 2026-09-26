#!/usr/bin/env Rscript
# Figure 4: the MLE portions of legacy/fig4_syn{1,2,3}_ntip100.R.
# Preserve the original data (including their eight-decimal rounding), bounds,
# and estimation routines; do not run the posterior chains in those scripts.

args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(file_arg) != 1L) stop("Run this entrypoint with Rscript.")
script_file <- gsub("~+~", " ", sub("^--file=", "", file_arg), fixed = TRUE)
script_dir <- dirname(normalizePath(script_file))
data_dir <- file.path(script_dir, "..", "synthetic", "data")
output_dir <- if (length(args)) args[[1]] else file.path(script_dir, "generated")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (pkg in c("ape", "phylodyn")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install required package: ", pkg)
}
if (!exists("bounded_skyline_ascent", asNamespace("phylodyn"))) {
  stop("Use the JuliaPalacios/phylodyn version described in README.md.")
}

bounds <- c(1, 0.7, 0.71)
# One-based rows; exactly match the embedded legacy vectors after rounding.
dataset_rows <- c(25L, 6L, 13L)
trajectories <- list(function(t) 1 + t - t,
                     function(t) 3 * exp(-t),
                     function(t) 25 * exp(-5 * t))
titles <- c("Ne(t) = 1", "Ne(t) = 3 exp(-t)", "Ne(t) = 25 exp(-5t)")
fits <- vector("list", 3L)

# Topology is immaterial to skyline coalescent times, but fix its random draw.
set.seed(1)
for (i in seq_len(3L)) {
  input <- new.env(parent = emptyenv())
  load(file.path(data_dir, sprintf("syn%d_bounded_ntip100.rda", i)), input)
  times <- round(as.numeric(input$coal_time[dataset_rows[i], ]), 8)
  stopifnot(length(times) == 99L, all(diff(times) > 0), max(times) < bounds[i])
  dat <- list(coal_times = times, lineages = 100:2,
              intercoal_times = c(times[1], diff(times)),
              samp_times = 0, n_sampled = 100)
  # The upstream optimizer prints its entire iteration trace; retain it on disk.
  capture.output(
    bounded <- phylodyn:::bounded_skyline_ascent(dat, bound = bounds[i]),
    file = file.path(output_dir, sprintf("figure4_syn%d_optimizer.log", i)))
  tree <- phylodyn:::generate_newick(dat)$newick
  standard <- ape::skyline(tree)
  x <- seq(0, bounds[i] + 1e-4, length.out = 101)
  stopifnot(all(is.finite(bounded$Ne)), all(bounded$Ne > 0))
  fits[[i]] <- list(standard = standard, bounded = bounded, times = times,
                     x = x, truth = trajectories[[i]](x))
  write.csv(data.frame(time = bounded$grid, Ne = bounded$Ne),
            file.path(output_dir, sprintf("figure4_syn%d_bc.csv", i)), row.names = FALSE)
  write.csv(data.frame(time = standard$time, Ne = standard$population.size),
            file.path(output_dir, sprintf("figure4_syn%d_sc.csv", i)), row.names = FALSE)
  write.csv(data.frame(time = x, Ne = trajectories[[i]](x)),
            file.path(output_dir, sprintf("figure4_syn%d_truth.csv", i)), row.names = FALSE)
}

pdf(file.path(output_dir, "figure4_mle.pdf"), width = 13, height = 4.5)
par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 1))
for (i in seq_len(3L)) {
  fit <- fits[[i]]
  standard <- fit$standard
  tick_height <- c(0.025, 0.003, 0.005)[i]
  yrange <- range(c(standard$population.size, fit$bounded$Ne,
                    fit$truth, tick_height)) * c(0.8, 1.2)
  # This is ape's skyline step convention, with explicit labels and full limits.
  plot(c(0, standard$time),
       c(standard$population.size, tail(standard$population.size, 1)),
       type = "s", xlim = c(max(standard$time), 0), ylim = yrange, log = "y",
       main = titles[i], xlab = "Time before sampling", ylab = expression(N[e](t)))
  lines(fit$x, fit$truth, col = "red", lwd = 3)
  lines(fit$bounded$grid, fit$bounded$Ne, col = "blue", lwd = 2, type = "s")
  legend("topleft", c("Truth", "SC-MLE", "BC-MLE"),
         col = c("red", "black", "blue"), lwd = c(3, 1.5, 2), bty = "n")
  points(fit$times, rep(tick_height, length(fit$times)),
         pch = 124, col = "black", cex = 0.6)
}
invisible(dev.off())
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Figure 4 MLE curves written to ", normalizePath(output_dir))
