"""Run the installed R package from Python; this is not a standalone Python package.

Example:
    python run_forecast.py
    python run_forecast.py --input events.csv --subject-id 1 --output forecasts.csv

R and coherentSurv must already be installed. Supply --rscript if Rscript is
not on PATH. The example uses exact interval histories, formula ~ x, recurrent
marks 1/2, months, landmark 0 and horizons 0.25/0.5/1. Adapt the R_CODE analysis
specification for your study; it is not an automatic analysis of arbitrary data.
"""
import argparse
import csv
import pathlib
import shutil
import subprocess
import tempfile

R_CODE = r'''
args <- commandArgs(trailingOnly = TRUE)
if (!requireNamespace("coherentSurv", quietly = TRUE))
  stop("Install the coherentSurv R package before running this example.")
library(coherentSurv)
if (args[1] != "--demo") {
  d <- read.csv(args[1], stringsAsFactors = FALSE, check.names = FALSE)
} else {
  data("sim_perinatal", package = "coherentSurv")
  d <- sim_perinatal
}
needed <- c("id", "start", "stop", "event_type", "terminal", "x")
if (!all(needed %in% names(d))) stop("Required columns: ", paste(needed, collapse = ", "))
subject <- if (args[2] != "--first") args[2] else as.character(d$id[1])
one <- d[as.character(d$id) == subject, , drop = FALSE]
if (!nrow(one)) stop("The requested subject ID is absent.")

# Study-specific choices: adapt these together with the event dictionary.
fit <- cs_fit(~ x, data = d, event_types = c("1", "2"),
                     time_unit = "months")
print(fit$convergence)
prediction <- predict(fit, newdata = one, landmark = 0,
                      horizon = c(0.25, 0.5, 1), type = "all")
result <- prediction[, c("unit", "horizon", "target", "estimate")]
stopifnot(nrow(result) > 0L, all(is.finite(result$estimate)))
write.csv(result, args[3], row.names = FALSE)
invisible(cs_audit_report(fit, predictions = prediction, format = "json",
                          file = paste0(args[3], ".audit.json")))
cat("R package version:", as.character(packageVersion("coherentSurv")), "\n")
cat("Fit status:", fit$convergence$code, "\n")
cat("Forecast CSV:", args[3], "\n")
cat("This worked example is not held-out predictive validation.\n")
'''

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=pathlib.Path, help="Optional interval-history CSV; default: bundled synthetic data.")
    parser.add_argument("--subject-id", default="", help="Subject to forecast; default: first subject.")
    parser.add_argument("--output", type=pathlib.Path, default=pathlib.Path("forecasts.csv"))
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable name or full path.")
    args = parser.parse_args()
    rscript = shutil.which(args.rscript)
    if rscript is None:
        parser.error("Rscript was not found. Install R and supply --rscript with its full executable path.")
    if args.input is not None and not args.input.is_file():
        parser.error("The input CSV does not exist.")
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    # A script file avoids platform-dependent quoting of multiline R code.
    with tempfile.TemporaryDirectory(prefix="coherentSurv_") as temporary:
        script = pathlib.Path(temporary) / "forecast.R"
        script.write_text(R_CODE, encoding="utf-8")
        subprocess.run(
            [rscript, "--vanilla", str(script),
             str(args.input.resolve()) if args.input else "--demo",
             args.subject_id or "--first", str(output)],
            check=True,
        )
    with output.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError("R returned an empty forecast table.")
    print(f"Loaded {len(rows)} forecast rows into Python.")
    for row in rows[:5]:
        print(row)
    print(f"Audit report: {output}.audit.json")

if __name__ == "__main__":
    main()
