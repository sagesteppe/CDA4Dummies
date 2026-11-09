## Decision-curve usefulness by AOA domain across assumed landscape prevalence.
## Self-contained: reads data/groundtruth_point_predictions.csv (no models / HDD).
## Net benefit at DECISION-RELEVANT thresholds (rare-species survey tolerates false visits),
## reweighted to a target prevalence, for raw-fixed vs prevalence-recalibrated probabilities.
suppressPackageStartupMessages(library(tidyverse))
HERE <- "/home/sagesteppe/Documents/assoRted/CDA4Dummies"
d <- read.csv(file.path(HERE,"data","groundtruth_point_predictions.csv"))
## RAW probabilities for all three designs (uniform; classic/results_classicsplit has no Platt
## rows, and the Note's question is PREVALENCE recalibration, not Platt). The 'recalibrated'
## scenario below is a logit shift to the target prevalence applied on top of raw.
d <- d[d$calibration=="raw", ]

pi_grid <- c(0.339,0.2,0.1,0.05,0.02,0.01,0.005,0.002)
pt_grid <- c(0.01,0.05,0.10)                 # decision-relevant survey thresholds
clamp   <- function(p) pmin(pmax(p,1e-6),1-1e-6)
recal   <- function(prob,p_obs,pi) plogis(qlogis(clamp(prob)) - qlogis(clamp(p_obs)) + qlogis(pi))

## weighted net benefit of the model and of treat-all at threshold pt, prevalence set by w
nb <- function(truth, prob, w, pt){
  n  <- sum(w); prev <- sum(w*truth)/n; pp <- prob >= pt
  nb_m <- sum(w[pp & truth==1])/n - sum(w[pp & truth==0])/n * (pt/(1-pt))
  nb_a <- prev - (1-prev)*(pt/(1-pt))
  c(model=nb_m, best_ref=max(nb_a,0))
}
reweight <- function(truth,pi){ n1<-sum(truth==1); n0<-sum(truth==0); ifelse(truth==1, pi/n1, (1-pi)/n0) }

cells <- d |> group_by(design,domain,model) |> group_split()
res <- map_dfr(cells, function(cc){
  if(length(unique(cc$Presence))<2 || nrow(cc)<8) return(NULL)
  p_obs <- mean(cc$Presence)
  expand_grid(pi=pi_grid, pt=pt_grid, scenario=c("raw (fixed)","recalibrated")) |>
    pmap_dfr(function(pi,pt,scenario){
      w  <- reweight(cc$Presence,pi)
      pr <- if(scenario=="recalibrated") recal(cc$prob,p_obs,pi) else cc$prob
      v  <- nb(cc$Presence, pr, w, pt)
      tibble(design=cc$design[1], domain=cc$domain[1], pi=pi, pt=pt, scenario=scenario,
             nb_model=v["model"], excess=v["model"]-v["best_ref"])
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
