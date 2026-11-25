# Run from the project root: Rscript replication/reproduce_cases.R run_name
args<-commandArgs(trailingOnly=TRUE)
name<-if(length(args))args[1]else"fresh_run"
if(!grepl("^[A-Za-z0-9_-]+$",name))stop("Use an alphanumeric run name.")
base<-file.path("replication","fresh_runs",name)
if(dir.exists(base))stop("Fresh-run destination already exists; choose a new name.")
dir.create(base,recursive=TRUE)
Sys.setenv(COHERENT4CAST_CASE_OUTPUT=file.path(base,"cases"),COHERENT4CAST_CASE_CACHE=file.path("build","fresh_case_fits",name),COHERENT4CAST_REUSE_CASE_FITS="0")
source("replication/run_case_studies.R",local=new.env(parent=globalenv()))
writeLines("Completed from fresh fitting; no cached model objects reused.",file.path(base,"COMPLETED.txt"))
