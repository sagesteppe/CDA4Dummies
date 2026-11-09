## FALSIFICATION TEST of the mechanism "summary-vs-exact gap is driven by the distance between a
## model's EVALUATION prevalence and the TARGET prevalence". We sweep the TARGET prevalence and,
## for each P:A ratio, track the gap. Prediction: each ratio's gap bottoms out when the target
## equals that ratio's natural evaluation prevalence (1/(1+ratio)); flip the target to 50% and the
## 1:1 models should thrive while 1:30 suffer -- the mirror of the 2.5% result. Strong form: rescale
## the x-axis to (target / natural prevalence) and every ratio should collapse onto one curve.
suppressPackageStartupMessages({library(data.table); library(dplyr); library(ggplot2)})
HERE <- "/home/sagesteppe/Documents/assoRted/CDA4Dummies"

ap_pointwise <- function(truth, prob, pi) {
  n1 <- sum(truth == 1); n0 <- sum(truth == 0); if (n1 == 0 || n0 == 0) return(NA_real_)
  w <- ifelse(truth == 1, pi / n1, (1 - pi) / n0)
  o <- order(prob, decreasing = TRUE); truth <- truth[o]; w <- w[o]
  tpw <- cumsum(w * (truth == 1)); fpw <- cumsum(w * (truth == 0))
  sum(c(tpw[1] / pi, diff(tpw / pi)) * (tpw / (tpw + fpw)))
}
ap_binormal <- function(mu, pi, grid = 3000) {
  f <- seq(1e-6, 1 - 1e-6, length.out = grid); tpr <- pnorm(qnorm(f) + mu)
  prec <- (pi * tpr) / (pi * tpr + (1 - pi) * f)
  sum(diff(c(0, tpr)) * (head(c(1, prec), -1) + tail(c(1, prec), -1)) / 2)
}
exact_nb <- function(truth, prob, pi, pt) {
  n1 <- sum(truth == 1); n0 <- sum(truth == 0)
  w <- ifelse(truth == 1, pi / n1, (1 - pi) / n0); n <- sum(w)
  vapply(pt, function(t) { pp <- prob >= t
    sum(w[pp & truth == 1]) / n - sum(w[pp & truth == 0]) / n * (t / (1 - t)) }, numeric(1))
}
binorm_nb <- function(mu, pi, pt) vapply(pt, function(t) {
  m <- ((1 - pi) / pi) * (t / (1 - t)); cc <- (log(m) + mu^2 / 2) / mu
  pi * pnorm(mu - cc) - (1 - pi) * pnorm(-cc) * (t / (1 - t)) }, numeric(1))

## mean gap over a prevalence-RELATIVE threshold window (same relative tolerance at every target)
mean_gap <- function(truth, prob, mu, pi_target) {
  pt <- pmin(pi_target * seq(0.25, 10, by = 0.25), 0.95)
  mean(binorm_nb(mu, pi_target, pt) - exact_nb(truth, prob, pi_target, pt))
}

d <- fread(file.path(HERE, "data", "holdout_point_predictions.csv"))
d <- d[calibration == "raw" & design == "spatial_mapping"]
d <- d[!(seed %in% c(5001, 5002, 5003))]           # drop the 3-seed replicate experiment; keep 5101-5801 reps
d[, pa_ratio := round(as.numeric(sub(".*PA1:([0-9.]+)DO.*", "\\1", model)))]  # collapse ~1:1 cluster to 1

## fit mu once per model (anchored at its own evaluation prevalence + PR-AUC)
mods <- d[, {
  if (uniqueN(Presence) < 2 || .N < 30) .SD[0] else {
    pih <- mean(Presence); pr <- ap_pointwise(Presence, prob, pih)
    mu  <- tryCatch(uniroot(function(m) ap_binormal(m, pih) - pr, c(1e-3, 8))$root, error = function(e) NA_real_)
    .(pa_ratio = pa_ratio[1], eval_prev = pih, mu = mu)
  }
}, by = model]
pt_list <- split(d[, .(Presence, prob)], d$model)

targets <- 10^seq(log10(0.01), log10(0.50), length.out = 13)
long <- rbindlist(lapply(seq_len(nrow(mods)), function(i) {
  m <- mods[i]; if (is.na(m$mu)) return(NULL)
  pp <- pt_list[[m$model]]
  data.table(model = m$model, pa_ratio = m$pa_ratio, eval_prev = m$eval_prev,
             target_prev = targets,
             gap = vapply(targets, function(tp) mean_gap(pp$Presence, pp$prob, m$mu, tp), numeric(1)))
}))
long[, nat_prev := 1 / (1 + pa_ratio)]

agg <- long[, .(gap = median(gap), n = .N), by = .(pa_ratio, nat_prev, target_prev)]
fwrite(agg, file.path(HERE, "data", "gap_vs_target_prevalence.csv"))

## Fig A: gap vs target prevalence, one line per ratio; dots mark each ratio's natural prevalence
markA <- long[, .(target_prev = nat_prev[1], gap = median(gap[which.min(abs(target_prev - nat_prev))])),
              by = pa_ratio]
pA <- ggplot(agg, aes(target_prev, gap, colour = nat_prev, group = pa_ratio)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(linewidth = 0.6) +
  geom_point(data = markA, aes(target_prev, gap), colour = "black", size = 1.3, inherit.aes = FALSE) +
  scale_x_log10(labels = scales::percent) +
  scale_colour_viridis_c(name = "natural\nprevalence", trans = "log", breaks = c(0.02, 0.05, 0.1, 0.25, 0.5)) +
  labs(x = "Target landscape prevalence", y = "Median net-benefit gap (summary - exact)",
       title = "Each P:A ratio's gap bottoms out at its own evaluation prevalence",
       subtitle = "Dots = each ratio's natural prevalence 1/(1+ratio). Above 0 = summary route optimistic.") +
  theme_bw(base_size = 11)
ggsave(file.path(HERE, "figures", "gap_vs_target_prevalence.png"), pA, width = 8, height = 5, dpi = 150)

## Fig B: strong form -- rescale x to target / natural prevalence; curves should collapse, cross 0 at 1
pB <- ggplot(agg, aes(target_prev / nat_prev, gap, colour = nat_prev, group = pa_ratio)) +
  geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 1, linetype = "dotted", colour = "grey55") +
  geom_line(linewidth = 0.6) +
  scale_x_log10() +
  scale_colour_viridis_c(name = "natural\nprevalence", trans = "log", breaks = c(0.02, 0.05, 0.1, 0.25, 0.5)) +
  labs(x = "Target prevalence / evaluation prevalence", y = "Median net-benefit gap (summary - exact)",
       title = "Collapse test: the gap is a function of target / evaluation prevalence",
       subtitle = "If the mechanism holds, every ratio falls on one curve crossing 0 at x = 1 (dotted).") +
  theme_bw(base_size = 11)
ggsave(file.path(HERE, "figures", "gap_collapse_test.png"), pB, width = 8, height = 5, dpi = 150)

cat("targets swept:", paste(round(targets, 3), collapse = " "), "\n")
cat("sign check (median gap, summary - exact):\n")
dcast(agg[target_prev %in% targets[c(1, 7, 13)]], pa_ratio + nat_prev ~ round(target_prev, 3),
      value.var = "gap")[order(pa_ratio)] |> print()
cat("\nwrote data/gap_vs_target_prevalence.csv, figures/gap_vs_target_prevalence.png, gap_collapse_test.png\n")
