#' Forecast by retaining a mixture of pre-landmark complete histories
#' @export
cs_predict_paths <- function(object, paths, owner, log_weights, landmark, horizon,
                             type = c("risk", "burden", "survival", "rmst", "all"),
                             quadrature = object$control$quadrature,...) {
  type <- match.arg(type)
  if (!inherits(object, "coherentSurv")) .cs_error("fit", "object must be a coherentSurv fit.")
  if (length(landmark) != 1L || !is.finite(landmark) || landmark <= 0)
    .cs_error("landmark", "A path mixture requires a positive finite landmark.")
  s <- object$schema
  d <- cs_as_trajectory(paths, s$unit, s$start, s$stop, s$event_type, s$terminal,
    s$observed, s$event_types, s$time_unit, "coarsened")
  schema <- attr(d, "cs_schema"); n <- length(unique(d$.cs_unit))
  if (length(schema$gaps) || length(schema$delayed_entry) || any(!d$.cs_observed))
    .cs_error("path_history", "Each candidate complete history must cover [0,landmark] without gaps or delayed entry.")
  if (any(d$.cs_stop > landmark) || any(d$.cs_terminal != 0L))
    .cs_error("path_future", "Paths must contain only information through the landmark and must be terminal-event-free there.")
  ends <- vapply(split(d$.cs_stop, d$.cs_unit), max, numeric(1))
  if (any(abs(ends - landmark) > 1e-12 * max(1, landmark)))
    .cs_error("path_history", "Each complete path must end at the declared landmark.")
  if (!is.numeric(owner) || length(owner) != n || any(!is.finite(owner)) ||
      any(owner < 1 | owner != as.integer(owner)) ||
      !identical(sort(unique(as.integer(owner))), seq_len(length(unique(owner)))))
    .cs_error("path_owner", "owner must map every path to contiguous positive observed-subject integers.")
  if (!is.numeric(log_weights) || length(log_weights) != n || any(!is.finite(log_weights)))
    .cs_error("path_weights", "Supply one finite log integration mass per pre-landmark path.")
  dd <- .cs_design(d, object$formula, object$breaks, object$control$history_cap,
    object$nonrecurrent, object$design$terms, object$design$xlevels, object$design$contrasts)
  X <- if (is.null(object$design$active)) dd$design else dd$design[, object$design$active, drop = FALSE]
  dd$design<-X
  z <- .cs_fit_likelihood(object$coefficients,dd,object$control,
    if(is.null(object$frailty))"gamma"else object$frailty,
    if(is.null(object$frailty_spec))2L else object$frailty_spec$points,information=FALSE)
  posterior <- numeric(n)
  for (rows in split(seq_len(n), owner)) {
    lw <- z$subject_loglik[rows] + log_weights[rows]
    w <- exp(lw - max(lw)); posterior[rows] <- w / sum(w)
  }
  # Conditional forecasts preserve each path's own posterior frailty law.
  # Average the forecasts, not the distribution parameters.
  exact_object <- object; exact_object$schema$observation <- "exact"
  pred <- predict(exact_object, newdata = paths, landmark = landmark,
    horizon = horizon, type = type, quadrature = quadrature,...)
  ids <- unique(as.character(d[[s$unit]])); index <- match(pred$unit, ids)
  if (anyNA(index)) .cs_error("path_history", "Prediction lost a declared candidate path.")
  pred$owner <- owner[index]; pred$path_weight <- posterior[index]
  groups <- interaction(pred$owner, pred$horizon, pred$target, drop = TRUE)
  draw_values<-attr(pred,"covariate_draw_forecasts")
  out <- lapply(split(seq_len(nrow(pred)), groups), function(rows) {
    ans <- pred[rows[1L], setdiff(names(pred), c("owner", "path_weight")), drop = FALSE]
    ans$unit <- as.character(pred$owner[rows[1L]])
    ans$estimate <- sum(pred$path_weight[rows] * pred$estimate[rows])
    if(!is.null(draw_values)) {
      joint_draw<-as.vector(draw_values[,rows,drop=FALSE]%*%pred$path_weight[rows])
      ans$covariate_mcse<-stats::sd(joint_draw)/sqrt(length(joint_draw))
    }
    for (v in intersect(c("series_error_bound", "probability_roundoff_bound","frailty_compression_bound"), names(ans)))
      ans[[v]] <- sum(pred$path_weight[rows] * pred[[v]][rows])
    ans
  })
  ans <- do.call(rbind, out); rownames(ans) <- NULL
  class(ans) <- c("coherentSurv_prediction", "cs_forecast", "data.frame")
  attr(ans, "path_weights") <- data.frame(path = ids, owner = owner, posterior_weight = posterior)
  family <- if (is.null(object$frailty)) "gamma" else object$frailty
  attr(ans, "assumptions") <- unique(c(
    "candidate paths and integration weights describe only the observed history through landmark",
    "terminal-event-free at landmark",
    paste("posterior path mixture retains path-specific frailty distributions; fitted family:", family),
    attr(pred, "assumptions"),
    "observation-fibre quadrature error is separate from the uniformization series bound"))
  ans
}
