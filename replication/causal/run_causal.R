# Executed S12 baseline-treatment workflow. All generated subjects are synthetic.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
here <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")) else normalizePath(".", winslash = "/")
lib <- Sys.getenv("COHERENTSURV_LIBRARY", "")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
library(coherentSurv)
out <- file.path(here, "results"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(20261123)
n <- 1600L; rows <- list(); index <- 0L
for (id in seq_len(n)) {
  x <- rbinom(1, 1, .5); w <- rbinom(1, 1, .3 + .4 * x)
  t <- 0; counts <- c(0, 0)
  while (t < 1) {
    a <- min(1, counts[1]); b <- min(1, counts[2])
    rates <- c(.8, .35, .4) * exp(c(.25*x-.35*w+.15*a+.10*b,
        .15*x+.10*w+.10*a+.20*b, .30*x-.25*w+.25*a+.20*b))
    ch <- .2 * exp(.2*x+.15*w)
    end <- min(1, t + rexp(1, sum(rates) + ch))
    mark <- if (end < 1) sample.int(4, 1, prob = c(rates, ch)) else 0L
    index <- index + 1L
    rows[[index]] <- data.frame(id=id,start=t,stop=end,event_type=if(mark %in% 1:2)mark else 0L,
        terminal=as.integer(mark==3),censored=as.integer(mark==4),x=x,treatment=w)
    if (mark %in% 1:2) counts[mark] <- counts[mark] + 1L
    t <- end
    if (mark %in% 3:4) break
  }
}
d <- do.call(rbind, rows)
write.csv(d,file.path(out,"synthetic_causal_histories.csv"),row.names=FALSE)
fit <- cs_causal(~x,d,"treatment",1,"censored","years",event_types=c("1","2"),
  spec=list(history_cap=1L,grid=64L),seed=81)
fine <- cs_causal(~x,d,"treatment",1,"censored","years",event_types=c("1","2"),
  folds=fit$folds,spec=list(history_cap=1L,grid=128L),seed=81)
effects <- cs_causal_contrast(fit)
simultaneous <- cs_causal_contrast(fit,interval="multiplier",draws=4000,seed=91)
ratios <- cs_causal_contrast(fit,"ratio")
oracle <- read.csv(file.path(out,"dense_gformula_oracle.csv"))
effects$oracle0 <- oracle$oracle[oracle$arm==0][match(effects$target,oracle$target[oracle$arm==0])]
effects$oracle1 <- oracle$oracle[oracle$arm==1][match(effects$target,oracle$target[oracle$arm==1])]
effects$oracle_difference <- effects$oracle1-effects$oracle0
write.csv(effects,file.path(out,"causal_contrasts.csv"),row.names=FALSE)
write.csv(simultaneous,file.path(out,"causal_simultaneous.csv"),row.names=FALSE)
write.csv(ratios,file.path(out,"causal_ratios.csv"),row.names=FALSE)
write.csv(fit$joint_EIF,file.path(out,"causal_individual_eif.csv"),row.names=FALSE)
write.csv(fit$joint_covariance,file.path(out,"causal_joint_covariance.csv"))
write.csv(data.frame(target=fit$target_names,grid64=unname(fit$estimate),grid128=unname(fine$estimate),
  absolute_difference=abs(unname(fit$estimate-fine$estimate))),file.path(out,"grid_refinement.csv"),row.names=FALSE)
diagnostics <- data.frame(arm=0:1,score=vapply(fit$arms,function(z)z$score,numeric(1)),
  iterations=vapply(fit$arms,function(z)z$iterations,numeric(1)),
  min_propensity=vapply(fit$arms,function(z)min(z$propensity),numeric(1)),
  min_censor_survival=vapply(fit$arms,function(z)exp(-max(z$censor_rates)),numeric(1)),
  interpolation_error=vapply(fit$arms,function(z)z$interpolation_error,numeric(1)))
write.csv(diagnostics,file.path(out,"causal_diagnostics.csv"),row.names=FALSE)
cs_audit_report(fit,format="json",file=file.path(out,"causal_audit.json"))
stopifnot(fit$convergence$score_solved,fine$convergence$score_solved,
  max(abs(fit$estimate-fine$estimate))<1e-4,
  max(abs(fit$covariance-cov(fit$EIF)/n))<1e-12)
print(effects);print(diagnostics)

# Independent numerical remainder check in S12's one-absorbing-event example.
# These calculations do not invoke the targeting engine.
set.seed(20261124); N <- 1000000L; lambda0 <- .9; censor0 <- .7; e0 <- .4
W <- rbinom(N,1,e0); T <- rexp(N,lambda0); C <- rexp(N,censor0)
end <- pmin(T,C,1); delta <- as.integer(T <= C & T <= 1)
truth <- -expm1(-lambda0)
cases <- data.frame(label=c("all_correct","event_correct","weight_product_correct","all_wrong"),
  lambda=c(.9,.9,1.3,1.3),censor=c(.7,.2,.7,.2),propensity=c(.4,.65,.4,.65))
integral_exp <- function(a,t) if(abs(a)<1e-10)t else expm1(a*t)/a
checks <- lapply(seq_len(nrow(cases)),function(i){
  row <- cases[i,]; l <- row$lambda; c <- row$censor; e <- row$propensity
  v <- -expm1(-l)
  residual <- exp(-l)*(delta*exp((l+c)*end)-l*integral_exp(l+c,end))
  phi <- v+W/e*residual
  analytic_bias <- (lambda0-l)*exp(-l)*
      (e0/e*integral_exp(l-lambda0+c-censor0,1)-integral_exp(l-lambda0,1))
  data.frame(setting=row$label,n=N,estimate=mean(phi),truth=truth,bias=mean(phi)-truth,
      analytic_remainder=analytic_bias,mcse=sd(phi)/sqrt(N),
      discrepancy=mean(phi)-truth-analytic_bias)
})
checks <- do.call(rbind,checks)
stopifnot(all(abs(checks$discrepancy)<5*checks$mcse))
write.csv(checks,file.path(out,"s12_remainder_check.csv"),row.names=FALSE)
print(checks)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
writeLines(c("Synthetic baseline-treatment experiment, not a clinical treatment-effect claim.",
  "1,600 subjects; two recurrent marks and one absorbing terminal event; binary observed baseline confounder.",
  "Independent dense matrix-exponential population oracle supplied by causal_oracle.py.",
  "Cross-fitted treatment probabilities, event rates and constant censor hazards; known generating model is within the declared state class.",
  "All tables are freshly calculated. Joint arm EIF covariance is retained. Grid refinement is a numerical check, not sampling coverage.",
  "The one-million-subject remainder check uses a separate one-terminal-event design and analytic integration.",
  "No new empirical causal coverage, universal efficiency or finite-sample guarantee is inferred from these checks."),file.path(out,"INTERPRETATION.txt"))
