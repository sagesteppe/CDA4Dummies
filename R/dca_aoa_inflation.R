## A priori AOA downgrade for the example model: the naive decision curve (scored on the whole
## hold-out) vs the AOA-honest curve (scored only on in-domain points — presences vs the HARD,
## environmentally-similar absences). Both reweighted to the SAME prevalence, so the gap between
## them is purely the inflation from trivially rejecting easy, far-from-suitable extrapolation
## absences. No ground truth used (hold-out predictions + presences-only DI).
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
HERE <- "/home/sagesteppe/Documents/assoRted/CDA4Dummies"
for (f in list.files(file.path(HERE, "R"), "^fn_.*\\.R$", full.names = TRUE)) source(f)
EX   <- "1arc-Iteration1-PA1:30DO:0-Seed5101"
PI   <- 0.05                                      # in-AOA landscape prevalence (both curves)

p <- fread(file.path(HERE,"data","holdout_point_predictions.csv"))[calibration=="raw" & design=="spatial_mapping"]
a <- fread(file.path(HERE,"data","holdout_point_aoa.csv"))
m <- merge(p, a[,.(model,point_id,DI,DI_thresh)], by=c("model","point_id")); m[, extrap := DI>DI_thresh]

## net benefit nb_exact(), treat-all nb_treat_all() (fn_net_benefit.R); roc_auc() (fn_roc.R)
e <- m[model==EX]; ind <- e[extrap==FALSE]
pt <- seq(0.001, 0.25, by=0.002)
curves <- rbind(
  data.table(pt=pt, nb=nb_exact(e$Presence,   e$prob,   PI, pt), strategy="SDM — naive (full hold-out)"),
  data.table(pt=pt, nb=nb_exact(ind$Presence, ind$prob, PI, pt), strategy="SDM — AOA-honest (in-domain)"),
  data.table(pt=pt, nb=nb_treat_all(PI, pt),                     strategy="Survey all (treat-all)"),
  data.table(pt=pt, nb=0,                                        strategy="Survey none (treat-none)"))
lv <- c("SDM — naive (full hold-out)","SDM — AOA-honest (in-domain)","Survey all (treat-all)","Survey none (treat-none)")
curves[, strategy := factor(strategy, levels=lv)]
is_ref <- grepl("treat", curves$strategy)

p1 <- ggplot(mapping=aes(pt, nb, colour=strategy)) +
  geom_line(data=curves[is_ref], aes(linetype=strategy), linewidth=1.0) +
  geom_line(data=curves[!is_ref], linewidth=0.9) +
  scale_colour_manual(values=c("SDM — naive (full hold-out)"="#e8963a","SDM — AOA-honest (in-domain)"="#2a78d6",
                               "Survey all (treat-all)"="#d1495b","Survey none (treat-none)"="grey45"), name=NULL) +
  scale_linetype_manual(values=c("Survey all (treat-all)"="12","Survey none (treat-none)"="solid"), guide="none") +
  scale_x_continuous(labels=scales::percent) +
  coord_cartesian(ylim=c(-0.03, PI*1.08)) +
  labs(x=expression("Threshold probability  "*p[t]), y="Net benefit",
       title=sprintf("The AOA downgrade: naive vs in-domain decision curve (prevalence %.1f%%)", 100*PI),
       subtitle="Both curves at the same prevalence; the gap is the inflation from trivially rejecting easy absences.") +
  theme_bw(base_size=11) + theme(legend.position="bottom")
ggsave(file.path(HERE,"figures","dca_aoa_inflation.png"), p1, width=9, height=5, dpi=150)

## ensemble ROC-AUC (prevalence-invariant), full vs in-domain -- roc_auc() in fn_roc.R
ens <- m[, .(roc_full=roc_auc(Presence,prob), roc_in=roc_auc(Presence[!extrap],prob[!extrap])), by=model]
fwrite(ens, file.path(HERE,"data","aoa_roc_drop_by_model.csv"))
cat(sprintf("example in-domain n=%d (%d pres); ROC full->in: %.3f -> %.3f\n",
  nrow(ind), sum(ind$Presence), roc_auc(e$Presence,e$prob), roc_auc(ind$Presence,ind$prob)))
cat(sprintf("ensemble ROC-AUC full median=%.3f, in-domain median=%.3f; dropped in %.0f%% of %d models\n",
  median(ens$roc_full,na.rm=T), median(ens$roc_in,na.rm=T), 100*mean(ens$roc_in<ens$roc_full,na.rm=T), nrow(ens)))
cat("wrote figures/dca_aoa_inflation.png, data/aoa_roc_drop_by_model.csv\n")
