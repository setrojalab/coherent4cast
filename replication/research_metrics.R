# Companion reporting functions, sourced explicitly; not coherent4cast exports.
# One independent validation unit per row at a prespecified landmark/horizon.
# The package fit and censoring assumptions remain necessary.

rc_research_report <- function(predicted, time, status, horizon,
                               thresholds = c(0.1, 0.2, 0.5),
                               censor_survival = NULL, bootstrap_reps = 0L,
                               seed = 20261123L, level = 0.95,
                               false_positive_cost = NULL,
                               false_negative_cost = NULL,
                               loss_unit = NULL) {
  if (!requireNamespace("coherent4cast", quietly = TRUE)) stop("Install coherent4cast first.")
  if (!is.numeric(thresholds) || !length(thresholds) || anyNA(thresholds) ||
      any(!is.finite(thresholds)) || any(thresholds <= 0 | thresholds >= 1) ||
      anyDuplicated(thresholds)) stop("Use unique thresholds strictly inside (0,1).")
  if (length(bootstrap_reps) != 1L || !is.numeric(bootstrap_reps) ||
      !is.finite(bootstrap_reps) || bootstrap_reps < 0 ||
      bootstrap_reps != floor(bootstrap_reps)) stop("bootstrap_reps must be a nonnegative integer.")
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must lie in (0,1).")
  costs <- c(false_positive_cost, false_negative_cost)
  has_cost <- !is.null(false_positive_cost) || !is.null(false_negative_cost)
  if (has_cost && (length(false_positive_cost) != 1L || length(false_negative_cost) != 1L ||
      !is.numeric(costs) || any(!is.finite(costs)) || any(costs < 0) ||
      length(loss_unit) != 1L || is.na(loss_unit) || !nzchar(loss_unit)))
    stop("Supply both nonnegative finite costs and a nonempty loss_unit.")

  calculate <- function(p, tt, ss, gg) {
    ans <- coherent4cast::validate(p, tt, ss, horizon, censor_survival = gg,
                              threshold = thresholds[1L])
    w <- ans$weights; y <- as.numeric(ss == 1 & tt <= horizon)
    if(any(!is.finite(w)))stop("Inverse censoring weights exceed the finite numerical range; revise the horizon or censoring model.")
    if(!any(w>0))stop("No observed cases or horizon survivors are available for validation.")
    wc <- sum(w * y); wn <- sum(w * (1-y)); n <- length(p)
    rows <- which(w > 0)
    ordered <- rows[order(p[rows], decreasing=TRUE)]
    groups <- if(length(ordered))split(ordered,cumsum(c(TRUE,diff(p[ordered])!=0))) else list()
    cum_tp <- cum_fp <- 0; previous_recall <- 0; ap <- 0
    curve <- data.frame(threshold=numeric(), precision=numeric(), recall=numeric())
    for (g in groups) {
      cum_tp <- cum_tp + sum(w[g] * y[g]); cum_fp <- cum_fp + sum(w[g] * (1-y[g]))
      recall <- if (wc > 0) cum_tp / wc else NA_real_
      precision <- cum_tp / (cum_tp + cum_fp)
      if (wc > 0) ap <- ap + (recall - previous_recall) * precision
      previous_recall <- recall
      curve <- rbind(curve, data.frame(threshold=p[g[1]], precision=precision, recall=recall))
    }
    # Require events and non-events for an informative discrimination summary.
    if (wc <= 0 || wn <= 0) ap <- NA_real_
    threshold_rows <- lapply(thresholds, function(cut) {
      a <- p >= cut
      tp <- sum(w*y*a); fn <- sum(w*y*(!a)); fp <- sum(w*(1-y)*a); tn <- sum(w*(1-y)*(!a))
      ratio <- function(x,z) if (z > 0) x/z else NA_real_
      nb <- (tp - fp*cut/(1-cut))/n
      all_nb <- (wc - wn*cut/(1-cut))/n
      loss <- if(has_cost)(fp/n)*false_positive_cost+(fn/n)*false_negative_cost else NA_real_
      if(has_cost && !is.finite(loss))stop("Requested fixed loss exceeds the finite numerical range; rescale loss units.")
      data.frame(threshold=cut, TP=tp, FN=fn, FP=fp, TN=tn,
        sensitivity=ratio(tp,tp+fn), specificity=ratio(tn,tn+fp),
        PPV=ratio(tp,tp+fp), NPV=ratio(tn,tn+fn),
        false_negative_rate=ratio(fn,tp+fn), false_positive_rate=ratio(fp,fp+tn),
        alert_fraction=mean(a), net_benefit=nb, net_benefit_all=all_nb,
        net_benefit_none=0, gain_over_all=nb-all_nb, gain_over_none=nb,
        fixed_loss=loss)
    })
    threshold_table <- do.call(rbind, threshold_rows)
    native <- ans$metrics
    native$reason <- ifelse(native$defined, "defined", "insufficient outcomes or denominator")
    if (!identical(ans$calibration_status, "estimated")) {
      bad <- native$metric %in% c("CITL", "calibration_slope")
      native$estimate[bad] <- NA_real_; native$defined[bad] <- FALSE
      native$reason[bad] <- ans$calibration_status
    }
    metric <- rbind(native, data.frame(metric="average_precision_IPCW", estimate=ap,
      defined=is.finite(ap), reason=if(is.finite(ap))"defined" else "requires events and non-events"))
    list(metrics=metric, thresholds=threshold_table, precision_recall=curve,
         native=ans, weighted_event_fraction=if(wc+wn>0)wc/(wc+wn)else NA_real_)
  }
  out <- calculate(predicted, time, status, censor_survival)
  out$uncertainty <- "No validation-metric confidence interval requested."
  if (bootstrap_reps > 0L) {
    if (length(seed) != 1L || !is.finite(seed)) stop("Supply one finite bootstrap seed.")
    old_kind <- RNGkind(); had_seed <- exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
    if(had_seed) old_seed <- get(".Random.seed", envir=.GlobalEnv)
    on.exit({do.call(RNGkind, as.list(old_kind)); if(had_seed)assign(".Random.seed",old_seed,envir=.GlobalEnv)
             else if(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv)}, add=TRUE)
    set.seed(seed)
    flatten <- function(a) c(stats::setNames(a$metrics$estimate,a$metrics$metric),
      stats::setNames(a$thresholds$net_benefit,paste0("NB_at_",thresholds)),
      stats::setNames(a$thresholds$fixed_loss,paste0("fixed_loss_at_",thresholds)))
    original <- flatten(out)
    boot <- matrix(NA_real_,bootstrap_reps,length(original),dimnames=list(NULL,names(original)))
    errors <- character(bootstrap_reps)
    for (b in seq_len(bootstrap_reps)) {
      ix <- sample.int(length(predicted), replace=TRUE)
      gg <- if(is.null(censor_survival))NULL else lapply(censor_survival,function(v)v[ix])
      boot[b,] <- tryCatch(flatten(calculate(predicted[ix],time[ix],status[ix],gg)),
        error=function(e){errors[b]<<-conditionMessage(e); rep(NA_real_,length(original))})
    }
    ok <- colSums(is.finite(boot))
    quant <- t(vapply(seq_along(original), function(j) {
      # Failed or undefined bootstrap draws are retained, never silently discarded.
      if(ok[j] == bootstrap_reps && bootstrap_reps >= 100L)
        stats::quantile(boot[,j],c((1-level)/2,1-(1-level)/2),names=FALSE) else c(NA_real_,NA_real_)
    },numeric(2)))
    out$bootstrap <- data.frame(metric=names(original),estimate=unname(original),
      lower=quant[,1],upper=quant[,2],finite_replicates=ok,total_replicates=bootstrap_reps,
      status=if(bootstrap_reps<100L)rep("insufficient_bootstrap_replicates",length(ok)) else
                    ifelse(ok==bootstrap_reps,"percentile_approximation","undefined_draws_interval_withheld"))
    out$bootstrap_draws <- boot; out$bootstrap_errors <- errors
    out$uncertainty <- paste("Independent-unit percentile bootstrap conditional on a fixed fitted prediction model;",
      if(is.null(censor_survival))"marginal KM censoring distribution is re-estimated per draw." else
      "supplied censoring probabilities are held fixed; their estimation uncertainty is excluded.",
      "Not simultaneous, finite-sample exact, or a substitute for retraining/internal validation.")
  }
  out$definitions <- c("Average precision = sum of recall increments times precision, with tied risks grouped; no trapezoidal interpolation.",
    "AUC and AP refer to event by horizon versus observed event-free beyond horizon, with IPCW; no competing-event recoding is automatic.",
    "Confusion cells are IPCW totals, not necessarily integer patient counts.",
    "Net benefit and fixed loss are threshold-weighted prognostic utilities, not causal treatment effects or demonstrated savings.",
    "Cost summaries require externally specified false-positive/false-negative losses; uncertainty in those losses is excluded.",
    "One row per independent eligible validation unit. Repeated landmarks or clustered units need a design-specific resampling scheme.")
  out$loss_unit <- if(has_cost)loss_unit else NA_character_
  out
}

rc_value_burden <- function(predictions, unit_cost, cost_unit) {
  required <- c("unit","landmark","horizon","target","estimate","time_unit")
  if (!is.data.frame(predictions) || !all(required %in% names(predictions)))
    stop("Use a forecast table containing unit, landmark, horizon, target, estimate and time_unit.")
  d <- as.data.frame(predictions[grepl("^burden_",predictions$target),required,drop=FALSE])
  if (!nrow(d) || any(!is.finite(d$estimate)) || any(d$estimate<0) || anyNA(d))
    stop("Require nonnegative finite burden forecasts and nonmissing keys.")
  if (anyDuplicated(d[c("unit","landmark","horizon","target")])) stop("Duplicate burden rows would double count costs.")
  if (!is.numeric(unit_cost) || is.null(names(unit_cost)) || anyDuplicated(names(unit_cost)) ||
      any(!is.finite(unit_cost)) || any(unit_cost<0) || !setequal(names(unit_cost),unique(d$target)))
    stop("Supply one named, finite, nonnegative cost for each and only each burden target.")
  if(length(cost_unit)!=1L || is.na(cost_unit) || !nzchar(cost_unit))stop("Declare cost_unit, including currency/year where applicable.")
  d$unit_cost <- unname(unit_cost[d$target]); d$expected_cost <- d$estimate*d$unit_cost
  if(any(!is.finite(d$expected_cost)))stop("Expected-cost multiplication exceeds the finite numerical range; rescale cost units.")
  total <- stats::aggregate(expected_cost~unit+landmark+horizon+time_unit,d,sum)
  if(any(!is.finite(total$expected_cost)))stop("Expected-cost aggregation exceeds the finite numerical range; rescale cost units.")
  total$cost_unit <- cost_unit
  list(components=d,total=total,
       interpretation="Expected additive event cost = sum of fixed per-event costs times expected counts, by linearity of expectation. No independence between event types is required.",
       assumptions="Externally fixed, non-overlapping cost components; no discounting, severity-dependent prices or capacity effects. Not observed expenditure, savings or a causal intervention effect. Joint forecast draws are required for joint uncertainty; interval widths and standard errors are not summed.")
}
