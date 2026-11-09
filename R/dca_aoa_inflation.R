## A priori AOA downgrade for the example model: the naive decision curve (scored on the whole
## hold-out) vs the AOA-honest curve (scored only on in-domain points — presences vs the HARD,
## environmentally-similar absences). Both reweighted to the SAME prevalence, so the gap between
## them is purely the inflation from trivially rejecting easy, far-from-suitable extrapolation
## absences. No ground truth used (hold-out predictions + presences-only DI).
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
HERE <- "/home/sagesteppe/Documents/assoRted/CDA4Dummies"
EX   <- "1arc-Iteration1-PA1:30DO:0-Seed5101"
PI   <- 0.05                                      # in-AOA landscape prevalence (both curves)

p <- fread(file.path(HERE,"data","holdout_point_predictions.csv"))[calibration=="raw" & design=="spatial_mapping"]
a <- fread(file.path(HERE,"data","holdout_point_aoa.csv"))
m <- merge(p, a[,.(model,point_id,DI,DI_thresh)], by=c("model","point_id")); m[, extrap := DI>DI_thresh]

nb_sdm <- function(truth, prob, pi, pt){
  n1<-sum(truth==1); n0<-sum(truth==0); w<-ifelse(truth==1, pi/n1, (1-pi)/n0); n<-sum(w)
  vapply(pt, function(t){ pp<-prob>=t
    sum(w[pp & truth==1])/n - sum(w[pp & truth==0])/n*(t/(1-t)) }, numeric(1))
}
e <- m[model==EX]; ind <- e[extrap==FALSE]
pt <- seq(0.001, 0.25, by=0.002)
curves <- rbind(
  data.table(pt=pt, nb=nb_sdm(e$Presence,   e$prob,   PI, pt), strategy="SDM — naive (full hold-out)"),
  data.table(pt=pt, nb=nb_sdm(ind$Presence, ind$prob, PI, pt), strategy="SDM — AOA-honest (in-domain)"),
  data.table(pt=pt, nb=PI-(1-PI)*(pt/(1-pt)),                   strategy="Survey all (treat-all)"),
  data.table(pt=pt, nb=0,                                      strategy="Survey none (treat-none)"))
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

## ensemble ROC-AUC (prevalence-invariant), full vs in-domain
rocauc <- function(truth,prob){ if(length(unique(truth))<2) return(NA_real_)
  r<-rank(prob); n1<-sum(truth==1); n0<-sum(truth==0); (sum(r[truth==1])-n1*(n1+1)/2)/(n1*n0) }
ens <- m[, .(roc_full=rocauc(Presence,prob), roc_in=rocauc(Presence[!extrap],prob[!extrap])), by=model]
fwrite(ens, file.path(HERE,"data","aoa_roc_drop_by_model.csv"))
cat(sprintf("example in-domain n=%d (%d pres); ROC full->in: %.3f -> %.3f\n",
  nrow(ind), sum(ind$Presence), rocauc(e$Presence,e$prob), rocauc(ind$Presence,ind$prob)))
cat(sprintf("ensemble ROC-AUC full median=%.3f, in-domain median=%.3f; dropped in %.0f%% of %d models\n",
  median(ens$roc_full,na.rm=T), median(ens$roc_in,na.rm=T), 100*mean(ens$roc_in<ens$roc_full,na.rm=T), nrow(ens)))
cat("wrote figures/dca_aoa_inflation.png, data/aoa_roc_drop_by_model.csv\n")
