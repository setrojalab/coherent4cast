# Small, wholly synthetic mapping examples; not certified/native adapters.
# Install recurcast and dependencies first. No external data are downloaded.
args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- if (length(args)) normalizePath(file.path(dirname(sub("^--file=", "", args[1])), ".."), winslash = "/") else normalizePath(".", winslash = "/")
lib <- Sys.getenv("RECURCAST_LIBRARY", "")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
library(recurcast)
out <- Sys.getenv("RECURCAST_MAPPING_OUTPUT", file.path(root, "replication/results/standards_mapping"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checks <- data.frame(check = character(), passed = logical())
check <- function(name, value) {
  ok <- isTRUE(value)
  checks <<- rbind(checks, data.frame(check = name, passed = ok))
  if (!ok) stop(name, call. = FALSE)
}
save_csv <- function(x, name) write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE, na = "")

# ADTTE-style example: preserve the selected parameter's derivation.
# OS has 0=death, 1=administrative censoring, 2=other censoring.
# Reverse follow-up has a DIFFERENT endpoint definition and is not pooled.
origin <- as.Date("2020-01-01")
adtte <- data.frame(USUBJID = rep(paste0("SYN", 1:4), 2),
  PARAMCD = rep(c("OS", "FUPTIMOS"), each = 4),
  STARTDT = origin, AVAL = rep(c(5, 10, 7, 12), 2),
  CNSR = c(0L, 1L, 2L, 0L, 1L, 0L, 0L, 1L),
  stringsAsFactors = FALSE)
adtte$ADT <- adtte$STARTDT + adtte$AVAL - 1
selected <- adtte[adtte$PARAMCD == "OS", ]
check("ADTTE parameter selection retains four unique OS rows", nrow(selected) == 4L && !anyDuplicated(selected$USUBJID))
check("ADTTE declared inclusive-day derivation is retained", all(selected$AVAL == as.numeric(selected$ADT - selected$STARTDT) + 1))
check("ADTTE source censor codes are explicitly declared", all(selected$CNSR %in% 0:2))
adtte_intervals <- data.frame(id = selected$USUBJID, start = 0,
  stop = selected$AVAL, event_type = "0", terminal = as.integer(selected$CNSR == 0),
  source_parameter = selected$PARAMCD, source_censor_code = selected$CNSR)
adtte_trajectory <- as_trajectory(adtte_intervals, time_unit = "day", event_types = character())
check("ADTTE positive censor codes map to zero, never negative status", identical(adtte_trajectory$terminal, c(1L, 0L, 0L, 1L)))
check("ADTTE two observed deaths and two censored subjects retained", sum(adtte_trajectory$.rc_terminal) == 2L && nrow(adtte_trajectory) == 4L)
check("ADTTE keeps new and legacy trajectory classes", inherits(adtte_trajectory, "rc_trajectory") && inherits(adtte_trajectory, "surv_trajectory"))
save_csv(adtte, "synthetic_adtte_source")
save_csv(adtte_intervals, "adtte_os_intervals")

# OMOP-style source tables: code 0 means unmapped, not a fabricated standard
# concept. EXAMPLE_A/B are artificial local source codes for this toy only.
# Real analyses must supply reviewed concepts, episode and ascertainment rules.
person <- data.frame(person_id = 1:3)
observation_period <- data.frame(person_id = 1:3,
  observation_period_start_date = origin,
  observation_period_end_date = origin + c(12, 15, 10))
condition_occurrence <- data.frame(condition_occurrence_id = 1:6,
  person_id = c(1, 1, 1, 2, 2, 3),
  condition_concept_id = 0L,
  condition_source_value = c("EXAMPLE_A", "EXAMPLE_A", "EXAMPLE_B", "EXAMPLE_A", "EXAMPLE_A", "EXAMPLE_B"),
  condition_start_date = origin + c(3, 3, 8, 4, 9, 2), stringsAsFactors = FALSE)
death <- data.frame(person_id = 1L, death_date = origin + 10)
check("OMOP toy person linkage is one-to-one at baseline", !anyDuplicated(person$person_id) && !anyDuplicated(observation_period$person_id) && all(condition_occurrence$person_id %in% person$person_id))
local_code_map <- c(EXAMPLE_A = "A", EXAMPLE_B = "B")
check("OMOP toy codes are unmapped artificial local codes", all(condition_occurrence$condition_concept_id == 0L) && all(condition_occurrence$condition_source_value %in% names(local_code_map)))
episode_key <- c("person_id", "condition_source_value", "condition_start_date")
episodes <- condition_occurrence[!duplicated(condition_occurrence[episode_key]), ]
check("OMOP declared same-day same-code episode rule removes one duplicate", nrow(episodes) == 5L)
rows <- list()
for (id in person$person_id) {
  op <- observation_period[observation_period$person_id == id, ]
  dd <- death$death_date[death$person_id == id]
  end_date <- if (length(dd)) min(op$observation_period_end_date, dd) else op$observation_period_end_date
  e <- episodes[episodes$person_id == id, ]
  e <- e[order(e$condition_start_date), ]
  check(paste0("OMOP subject ", id, " has strictly ordered observed episodes before exit"),
    all(e$condition_start_date > op$observation_period_start_date) &&
      all(e$condition_start_date < end_date) && !anyDuplicated(e$condition_start_date))
  times <- as.numeric(e$condition_start_date - op$observation_period_start_date)
  stop_time <- as.numeric(end_date - op$observation_period_start_date)
  terminal <- as.integer(length(dd) == 1L && dd <= op$observation_period_end_date)
  rows[[length(rows) + 1L]] <- data.frame(id = id,
    start = c(0, times), stop = c(times, stop_time),
    event_type = c(unname(local_code_map[e$condition_source_value]), "0"),
    terminal = c(rep(0L, length(times)), terminal), observed = TRUE,
    terminal_cause = c(rep(NA_character_, length(times)), if (terminal) "death" else NA_character_))
}
omop_intervals <- do.call(rbind, rows)
omop_trajectory <- as_trajectory(omop_intervals, event_types = c("A", "B"),
  observed = "observed", terminal_cause = "terminal_cause", time_unit = "day")
check("OMOP transformed histories contain five episodes and one death", sum(omop_trajectory$.rc_mark > 0L) == 5L && sum(omop_trajectory$.rc_terminal) == 1L)
check("OMOP person-time stops at death or declared observation end", sum(omop_trajectory$stop - omop_trajectory$start) == 35)
check("OMOP exact synthetic monitoring has no gaps or delayed entry", isTRUE(attr(omop_trajectory, "observation_process_details")$continuous_observation_from_zero))
check("OMOP observed pattern does not infer a scientific model", identical(attr(omop_trajectory, "rc_schema")$regime_details$model_inferred, FALSE))
save_csv(person, "synthetic_omop_person")
save_csv(observation_period, "synthetic_omop_observation_period")
save_csv(condition_occurrence, "synthetic_omop_condition_occurrence")
save_csv(death, "synthetic_omop_death")
save_csv(episodes, "omop_toy_episodes")
save_csv(omop_intervals, "omop_toy_intervals")
save_csv(checks, "mapping_checks")
rc_audit_report(omop_trajectory, source_files = file.path(root, "replication/standards_mapping.R"),
  format = "json", file = file.path(out, "mapping_audit.json"))
writeLines(c(
  "All source rows are constructed synthetic examples. No patient data are present.",
  "These scripts illustrate analyst-controlled mappings, not native, certified or generally validated CDISC/OMOP adapters.",
  "ADTTE: select OS; preserve AVAL and its example-specific inclusive-day derivation; use CNSR==0, not 1-CNSR. Do not mix reverse follow-up parameters.",
  "OMOP: artificial unmapped local codes EXAMPLE_A/B (concept_id=0); reviewed standard concept sets are required in real work.",
  "OMOP toy episode rule: collapse same person/source-code/date. Other episode rules need clinical justification; simultaneous distinct marks require a declared observation model.",
  "OMOP elapsed days are date differences from the observation origin; death ends follow-up. Table presence alone does not establish complete death or event ascertainment.",
  "Continuous exact monitoring is true only by construction in this synthetic fixture; real recording gaps must not be encoded as observed zeros.",
  "No likelihood, treatment effect or regulatory compliance is inferred from these mappings."), file.path(out, "INTERPRETATION.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(checks, row.names = FALSE)
cat("All", nrow(checks), "synthetic mapping checks passed. Outputs:", out, "\n")
