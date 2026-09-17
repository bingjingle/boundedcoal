#!/usr/bin/env Rscript
# Redraw every figure in base R, from the CSVs in ../coords/ only.
#
# Nothing here reads a sampler output or an .npz -- the coordinates are the whole
# input.  So to change colours, line widths, axis limits or the credible-band
# convention, edit the few constants below and re-run:
#
#     Rscript replot.R                # -> figures_R/
#
# Usage from an interactive session, if you would rather drive it yourself:
#     source("replot.R"); panel("syn2", "RI_BM")
#
# ---- the only things you are likely to want to change ----------------------
COL_BC  <- "#d62728"     # bounded coalescent
COL_SC  <- "#1f77b4"    # standard coalescent
ALPHA   <- 0.25      # band fill opacity, as adjustcolor(., alpha.f=) in covid_plot.R
BAND    <- "eq"      # "eq" = equal-tailed 2.5/97.5% (published); "hpd" = minimal width
LWD_MED <- 2; LWD_BAND <- 1
# ---------------------------------------------------------------------------

here    <- tryCatch(dirname(normalizePath(sub("^--file=", "", grep("^--file=",
             commandArgs(FALSE), value = TRUE)[1]))), error = function(e) ".")
if (is.na(here) || !nzchar(here)) here <- "."

here<-"~/Documents/boundedcoal/boundedcoal/reproduce/"
COORD <- file.path(here, "coords")
OUTD  <- file.path(here, "figures_R")
dir.create(OUTD, showWarnings = FALSE)

ROWS <- list(syn1 = list(lab = expression(N[e1](t) == 1),           tau = 1.00, ylim = 5),
             syn2 = list(lab = expression(N[e2](t) == 3 * e^-t),    tau = 0.70, ylim = 14),
             syn3 = list(lab = expression(N[e3](t) == 25 * e^(-5*t)), tau = 0.71, ylim = 120))
            # syn4 = list(lab = expression(N[e4](t) == 25 * e^(-5*t)), tau = 0.55, ylim = 120))
#COLS <- list(RI_BM = "RI-BM", RI_SE = "RI-SE", DIS = "Discrete (phylodyn)")
COLS <- list(RI_BM = "RI", DIS = "Discrete")
T0   <- as.Date("2020-06-08")      # covid_plot.R: date of sampling

band_cols <- function(d, lik) {
  lo <- paste0(lik, "_low_", BAND); hi <- paste0(lik, "_hi_", BAND)
  if (!lo %in% names(d)) { lo <- paste0(lik, "_low_eq"); hi <- paste0(lik, "_hi_eq") }
  list(lo = d[[lo]], hi = d[[hi]])
}

draw_one <- function(d, lik, col, step) {
  if (!paste0(lik, "_med") %in% names(d)) return(invisible())
  b <- band_cols(d, lik); x <- d$x
  ty <- if (step) "S" else "l"
  polygon(c(x, rev(x)), c(b$lo, rev(b$hi)),
          col = adjustcolor(col, alpha.f = ALPHA), border = NA)
  lines(x, b$lo, col = col, lwd = LWD_BAND, lty = 2, type = ty)
  lines(x, b$hi, col = col, lwd = LWD_BAND, lty = 2, type = ty)
  lines(x, d[[paste0(lik, "_med")]], col = col, lwd = LWD_MED, type = ty)
}

panel <- function(tag, key, ylim = NULL, main = NULL) {
  f <- file.path(COORD, sprintf("panel_%s_%s.csv", tag, key))
  if (!file.exists(f)) { message("missing: ", basename(f)); return(invisible(FALSE)) }
  d <- read.csv(f); step <- key == "DIS"
  cfg <- ROWS[[tag]]
  if (is.null(ylim)) ylim <- if (!is.null(cfg)) cfg$ylim else
    stats::quantile(unlist(d[grep("_hi_", names(d))]), 0.97, na.rm = TRUE)
  tau <- if (!is.null(cfg)) cfg$tau else max(d$x)
  covid <- tag == "covid"
  plot(NA, xlim = c(tau, 0), ylim = if (tag=="covid") c(0,10) else c(0, ylim),
       xlab = if (tag=="syn3") "Time (past to present)" else "", ylab = if (key=="RI_BM") expression(N[e](t)) else "",
       main = if (is.null(main)) COLS[[key]] else main,
       xaxt = if (covid) "n" else "s", bty = "l")
  if (covid) {
    at <- seq(0, tau, by = 0.1)
    axis(1, at = at, labels = paste0(format(T0 - round(at * 365), "%b"), "'",
                                     format(T0 - round(at * 365), "%y")), cex.axis = 0.8)
  }
  draw_one(d, "sc", COL_SC, step)          # BC drawn last, on top
  draw_one(d, "bc", COL_BC, step)
  if ("truth" %in% names(d)) lines(d$x, d$truth, col = "black", lwd = LWD_MED)
  ev <- file.path(COORD, sprintf("events_%s.csv", tag))
  if (file.exists(ev)) rug(read.csv(ev)$coal_time, ticksize = 0.03, lwd = 0.7)
  invisible(TRUE)
}

# ---- the 4 x 3 simulation grid ---------------------------------------------
pdf(file.path(OUTD, "panels.pdf"), width = 10.5, height = 12)
pdf(file.path(OUTD, "panels.pdf"))
#op <- par(mfrow = c(4, 3), mar = c(4, 4.2, 2.4, 1), oma = c(3, 2, 0, 0))
op <- par(mfrow = c(3, 2), mar = c(4, 4.2, 2.4, 1), oma = c(3, 2, 0, 0),cex.main = 1.3,
          cex.lab  = 1.2,
          cex.axis = 1.2)
for (tag in names(ROWS)) for (key in names(COLS))
  panel(tag, key, main = if (tag == "syn1") COLS[[key]] else "")
par(op)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", horiz = TRUE, bty = "n",
       legend = c("Truth", "BC (bounded)", "SC (standard)"),
       col = c("black", COL_BC, COL_SC), lwd = LWD_MED)
invisible(dev.off())

# ---- individual panels, as main2.tex includes them -------------------------
for (tag in names(ROWS)) for (key in names(COLS)) {
  pdf(file.path(OUTD, sprintf("%s_%s.pdf", tag, key)), width = 3.6, height = 3.2)
  par(mar = c(4, 4.2, 2.4, 1)); panel(tag, key); invisible(dev.off())
}


df <- data.frame(
  month = factor(c("Jan","Feb","Mar","Apr","May","Jun"),
                 levels = c("Jan","Feb","Mar","Apr","May","Jun"),
                 ordered = TRUE),
  cases = c(20, 44, 6947, 8698, 7315, 13662)
)


# ---- the three COVID figures -----------------------------------------------
pdf(file.path(OUTD, sprintf("covid_2_%s.pdf", key)), width = 8.6, height = 3.6)
par(mfrow=c(1,3),mar = c(4, 4.2, 2.4, 1));  

for (key in names(COLS)) {
  f <- file.path(COORD, sprintf("panel_covid_%s.csv", key))
  if (!file.exists(f)) next
  panel("covid", key)
  if (key=="RI_BM"){legend("top", horiz = TRUE, bty = "n",
                        legend = c("BC (bounded)", "SC (standard)"),
                        col = c(COL_BC, COL_SC), lwd = LWD_MED)}
}


  barplot(df$cases, names.arg=df$month, col="steelblue", border=NA,
          xlab="Month of 2020", ylab="Cases", cex.names=0.8)
  
  invisible(dev.off())
}
cat("wrote", length(list.files(OUTD)), "files to figures_R/\n")
