## Decision-curve usefulness by AOA domain across assumed landscape prevalence.
## Self-contained: reads data/groundtruth_point_predictions.csv (no models / HDD).
## Net benefit at DECISION-RELEVANT thresholds (rare-species survey tolerates false visits),
## reweighted to a target prevalence, for raw-fixed vs prevalence-recalibrated probabilities.
suppressPackageStartupMessages(library(tidyverse))
HERE <- "/home/sagesteppe/Documents/assoRted/CDA4Dummies"
for (f in list.files(file.path(HERE, "R"), "^fn_.*\\.R$", full.names = TRUE)) source(f)
d <- read.csv(file.path(HERE,"data","groundtruth_point_predictions.csv"))
## RAW probabilities for all three designs (uniform; classic/results_classicsplit has no Platt
## rows, and the Note's question is PREVALENCE recalibration, not Platt). The 'recalibrated'
## scenario below is a logit shift to the target prevalence applied on top of raw.
d <- d[d$calibration=="raw", ]

pi_grid <- c(0.339,0.2,0.1,0.05,0.02,0.01,0.005,0.002)
pt_grid <- c(0.01,0.05,0.10)                 # decision-relevant survey thresholds
## helpers: nb_exact(), nb_treat_all() (fn_net_benefit.R); recalibrate_to_prevalence() (fn_calibration.R)
## best reference strategy = max(treat-all, treat-none=0) at the reweighted prevalence pi
best_ref <- function(pi, pt) max(nb_treat_all(pi, pt), 0)

cells <- d |> group_by(design,domain,model) |> group_split()
res <- map_dfr(cells, function(cc){
  if(length(unique(cc$Presence))<2 || nrow(cc)<8) return(NULL)
  p_obs <- mean(cc$Presence)
  expand_grid(pi=pi_grid, pt=pt_grid, scenario=c("raw (fixed)","recalibrated")) |>
    pmap_dfr(function(pi,pt,scenario){
      pr <- if(scenario=="recalibrated") recalibrate_to_prevalence(cc$prob,p_obs,pi) else cc$prob
      nb_m <- nb_exact(cc$Presence, pr, pi, pt)
      tibble(design=cc$design[1], domain=cc$domain[1], pi=pi, pt=pt, scenario=scenario,
             nb_model=nb_m, excess=nb_m-best_ref(pi,pt))
    })
})
agg <- res |> group_by(design,domain,pi,pt,scenario) |>
  summarise(nb_model=median(nb_model), excess=median(excess),
            frac_useful=mean(excess>0), n_models=n(), .groups="drop")
agg$domain <- factor(agg$domain, levels=c("Interp <=270m","Interp >270m","Extrapolation"))
write.csv(agg, file.path(HERE,"data","dca_prevalence_sweep.csv"), row.names=FALSE)

options(width=170)
cat("=== median net benefit at pt=0.05 (>0 and above best reference = useful) ===\n")
agg |> filter(pt==0.05) |> transmute(design,domain,scenario,pi,excess=round(excess,4)) |>
  pivot_wider(names_from=pi, values_from=excess) |> arrange(domain,design,scenario) |>
  as.data.frame() |> print(row.names=FALSE)

## plot: excess net benefit vs prevalence at pt=0.05, by domain, design colour, scenario linetype
p <- ggplot(filter(agg,pt==0.05), aes(pi, excess, colour=design, linetype=scenario)) +
  geom_hline(yintercept=0, colour="grey60", linewidth=0.3) +
  geom_line(linewidth=0.7) + geom_point(size=1.2) +
  facet_wrap(~domain) + scale_x_log10() +
  labs(x="Assumed landscape prevalence (log)", y="Excess net benefit vs best reference (pt=0.05)",
       title="Is an SDM useful at pixel prevalence? Net benefit at a survey-relevant threshold",
       subtitle="Above 0 = beats treat-none/treat-all. Faceted by AOA domain.") +
  theme_bw(base_size=11) + theme(legend.position="bottom")
ggsave(file.path(HERE,"figures","dca_prevalence_sweep_pt05.png"), p, width=10, height=4.5, dpi=150)
cat("\nwrote data/dca_prevalence_sweep.csv and figures/dca_prevalence_sweep_pt05.png\n")
