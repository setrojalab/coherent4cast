# Static scientific plots from an executed cross_field_usability.R run.
out <- Sys.getenv("COHERENT4CAST_CROSS_FIELD_OUTPUT", "replication/results/cross_field_usability")
report <- readRDS(file.path(out, "companion_report.rds"))
if (!is.null(report$bootstrap_draws)) write.csv(data.frame(replicate = seq_len(nrow(report$bootstrap_draws)),
  report$bootstrap_draws, check.names = FALSE), file.path(out, "companion_bootstrap_draws.csv"), row.names = FALSE)
if (!is.null(report$bootstrap_errors)) write.csv(data.frame(replicate = seq_along(report$bootstrap_errors),
  error = report$bootstrap_errors), file.path(out, "companion_bootstrap_errors.csv"), row.names = FALSE)
curve <- report$precision_recall
decision <- report$thresholds[order(report$thresholds$threshold), ]
export <- function(name, draw) {
  grDevices::pdf(file.path(out, paste0(name, ".pdf")), width = 7, height = 5)
  draw(); grDevices::dev.off()
  grDevices::svg(file.path(out, paste0(name, ".svg")), width = 7, height = 5)
  draw(); grDevices::dev.off()
  grDevices::png(file.path(out, paste0(name, ".png")), width = 2100, height = 1500, res = 300)
  draw(); grDevices::dev.off()
}
export("precision_recall", function() {
  par(mar = c(5.1, 4.7, 3.2, 1.2), las = 1, family = "sans")
  plot(curve$recall, curve$precision, type = "S", col = "#24506A", lwd = 2,
    xlim = c(0, 1), ylim = c(0, 1), xlab = "Sensitivity (recall)", ylab = "Positive predictive value (precision)",
    main = "Terminal-event prediction: precision and recall", bty = "l")
  abline(h = report$weighted_event_fraction, lty = 2, col = "#9C5A34", lwd = 1.5)
  ap <- report$metrics$estimate[report$metrics$metric == "average_precision_IPCW"]
  legend("topright", c(sprintf("IPCW average precision = %.3f", ap),
    sprintf("Weighted event fraction = %.3f", report$weighted_event_fraction)),
    col = c("#24506A", "#9C5A34"), lty = c(1, 2), lwd = c(2, 1.5), bty = "n", cex = .9)
  mtext("One synthetic cohort; landmark 1 year, horizon 1 year; fixed model", side = 1, line = 3.6, cex = .78)
})
export("decision_curve", function() {
  par(mar = c(5.1, 4.7, 3.2, 1.2), las = 1, family = "sans")
  yy <- range(c(decision$net_benefit, decision$net_benefit_all, 0))
  plot(decision$threshold, decision$net_benefit, type = "b", pch = 16, col = "#24506A", lwd = 2,
    ylim = yy, xlab = "Decision threshold probability", ylab = "Net benefit per eligible unit",
    main = "Terminal-event prediction: decision curves", bty = "l")
  lines(decision$threshold, decision$net_benefit_all, type = "b", pch = 1, lty = 2, col = "#9C5A34", lwd = 1.5)
  abline(h = 0, lty = 3, col = "#555555", lwd = 1.5)
  legend("bottomleft", c("Model flags", "Flag all", "Flag none"),
    col = c("#24506A", "#9C5A34", "#555555"), lty = c(1, 2, 3), pch = c(16, 1, NA), bty = "n", cex = .9)
  mtext("Prespecified grid; IPCW totals; prognostic utility, not an intervention effect", side = 1, line = 3.6, cex = .76)
})
writeLines(c("Both figures use one synthetic held-out cohort and the same frozen fitted predictor.",
  "IPCW assumes marginal independent censoring in this generating experiment.",
  "Average precision is the step integral sum(delta recall * precision), not trapezoidal PR-AUC.",
  "Decision-curve lines join evaluated threshold points for readability; values between points were not evaluated.",
  "Net benefit does not establish prospective impact, causal benefit or cost savings.",
  "No confidence ribbon is drawn; bootstrap uncertainty is reported in companion_bootstrap.csv."),
  file.path(out, "FIGURE_NOTES.txt"))
