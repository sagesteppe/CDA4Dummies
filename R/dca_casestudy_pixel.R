## Case study, PIXEL (naive) formulation: a classic decision curve from the real
## E. coloradense ground-truth predictions. Each ground-truth POINT is an independent
## decision unit (survey it / don't). Self-contained: reads the per-point writeout only.
##
## Net benefit vs threshold probability p_t, faceted by AOA domain, with the two reference
## strategies (treat-all = survey everything, treat-none = survey nothing = 0). The model
## curve is the median across the raw models in each domain.
##
## PREVALENCE: the field sample is enriched and intentionally "difficult", so we DON'T score
## at the sample prevalence. We reweight every point to a target LANDSCAPE prevalence PI with
## the same case weights as dca_prevalence_sweep.R (presence = pi/n1, absence = (1-pi)/n0).
PI <- 0.05                                # target landscape prevalence
suppressPackageStartupMessages(library(tidyverse))
HERE <- "/home/sagesteppe/Documents/assoRted/CDA4Dummies"

d <- read.csv(file.path(HERE, "data", "groundtruth_point_predictions.csv"))
d <- d[d$calibration == "raw", ]
d$domain <- factor(d$domain, levels = c("Interp <=270m", "Interp >270m", "Extrapolation"))

pt_grid <- seq(0.005, 0.15, by = 0.0025) # decision-relevant thresholds (rare-species survey)

reweight <- function(truth, pi) {
  n1 <- sum(truth == 1); n0 <- sum(truth == 0)
  ifelse(truth == 1, pi / n1, (1 - pi) / n0)
}

## prevalence-reweighted net benefit of a single model's predictions across the pt grid
nb_model_curve <- function(truth, prob, pt_grid, pi) {
  w <- reweight(truth, pi); n <- sum(w)
  vapply(pt_grid, function(pt) {
    pp <- prob >= pt
    sum(w[pp & truth == 1]) / n - sum(w[pp & truth == 0]) / n * (pt / (1 - pt))
  }, numeric(1))
}

## --- model curves: net benefit per (domain, model), then median across models ---
model_cells <- d |> group_by(domain, design, model) |> group_split()
model_nb <- map_dfr(model_cells, function(cc) {
  if (length(unique(cc$Presence)) < 2 || nrow(cc) < 8) return(NULL)
  tibble(domain = cc$domain[1], model = cc$model[1],
         pt = pt_grid, nb = nb_model_curve(cc$Presence, cc$prob, pt_grid, PI))
})
model_curve <- model_nb |>
  group_by(domain, pt) |>
  summarise(net_benefit = median(nb), n_models = n_distinct(model), .groups = "drop") |>
  mutate(strategy = "SDM (pixel)")

## --- reference strategies at the target landscape prevalence PI (same in every domain) ---
ref_curve <- tibble(pt = pt_grid) |>
  mutate(treat_all  = PI - (1 - PI) * (pt / (1 - pt)),
         treat_none = 0) |>
  crossing(domain = levels(d$domain)) |>
  mutate(domain = factor(domain, levels = levels(d$domain))) |>
  pivot_longer(c(treat_all, treat_none), names_to = "strategy", values_to = "net_benefit") |>
  mutate(strategy = recode(strategy, treat_all = "Survey all (treat-all)",
                           treat_none = "Survey none (treat-none)")) |>
  select(domain, pt, strategy, net_benefit)

curves <- bind_rows(
  model_curve |> select(domain, pt, strategy, net_benefit),
  ref_curve
) |>
  mutate(strategy = factor(strategy, levels = c(
    "SDM (pixel)", "Survey all (treat-all)", "Survey none (treat-none)")))

write.csv(curves, file.path(HERE, "data", "dca_casestudy_pixel.csv"), row.names = FALSE)

## --- plot ---
pal <- c("SDM (pixel)" = "#2a78d6",
         "Survey all (treat-all)" = "#d1495b",
         "Survey none (treat-none)" = "grey45")

is_ref <- curves$strategy != "SDM (pixel)"
p <- ggplot(mapping = aes(pt, net_benefit, colour = strategy)) +
  ## bold-dashed references underneath a thinner SDM line, so where the SDM coincides with
  ## treat-all (in-domain, low thresholds) the red dashes still show through
  geom_line(data = curves[is_ref, ], aes(linetype = strategy), linewidth = 1.1) +
  geom_line(data = curves[!is_ref, ], linewidth = 0.6) +
  facet_wrap(~domain) +
  scale_linetype_manual(values = c("Survey all (treat-all)" = "12",
                                   "Survey none (treat-none)" = "solid"), guide = "none") +
  scale_colour_manual(values = pal, name = NULL) +
  coord_cartesian(ylim = c(-0.02, PI * 1.05)) +
  labs(x = expression("Threshold probability  " * p[t]),
       y = "Net benefit",
       title = sprintf("Pixel formulation at %.0f%% landscape prevalence, by AOA domain", 100 * PI),
       subtitle = "Each ground-truth point is one decision unit. Curve above both references = the model adds net benefit.") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(HERE, "figures", "dca_casestudy_pixel.png"), p,
       width = 10, height = 4.2, dpi = 150)

cat(sprintf("target landscape prevalence PI = %.3f\n", PI))
cat("where SDM beats both references (max p_t with SDM net benefit > treat-all and > 0):\n")
model_curve |>
  left_join(tibble(pt = pt_grid, ta = PI - (1 - PI) * (pt_grid / (1 - pt_grid))), by = "pt") |>
  group_by(domain) |>
  summarise(max_pt_useful = { ok <- net_benefit > pmax(ta, 0); if (any(ok)) max(pt[ok]) else NA_real_ },
            .groups = "drop") |>
  as.data.frame() |> print(row.names = FALSE)
cat("\nwrote data/dca_casestudy_pixel.csv and figures/dca_casestudy_pixel.png\n")
