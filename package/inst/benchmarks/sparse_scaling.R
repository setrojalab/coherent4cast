# Deterministic scaling benchmark. Run after installing coherentSurv.
# No benchmark is used to tune or select its own verification tolerances.
library(coherentSurv)
sizes <- c(1000L, 10000L, 50000L)
results <- list()
for (J in sizes) {
  M <- min(J, 1000L); k <- seq_len(J)-1L; map <- k %% M+1L
  p <- tabulate(map, M)/J
  denominator <- .2+.8*((k*104729) %% 50023)/50023
  numerator <- denominator*(.1+.8*((k*7919) %% 50021)/50021)
  elapsed <- system.time(z <- cs_bounds_sparse(p, numerator, denominator,
      kappa=3, recording_fraction=.6, counts=tabulate(map,M)*20L,
      observation_map=map))["elapsed"]
  results[[length(results)+1L]] <- data.frame(module="bounds_and_projection", states_or_atoms=J,
    fibres=M, elapsed_seconds=unname(elapsed), input_bytes=sum(vapply(list(p,numerator,denominator,map),function(x) as.numeric(object.size(x)),numeric(1))),
    native_workspace_estimate=z$plugin$workspace_bytes,
    endpoint_or_probability=z$upper, certificate=max(z$plugin$gap,z$outer$gap))
  q <- cs_generator_csr(J, seq_len(J-1L), seq_len(J-1L)+1L, rep(.1,J-1L))
  probability <- numeric(J); probability[1] <- 1
  elapsed <- system.time(y <- cs_propagate_sparse(probability,q,10,rep(1,J)))["elapsed"]
  stopifnot(abs(y$probability[1]-exp(-1))<1e-12,abs(y$reward-10)<1e-10)
  results[[length(results)+1L]] <- data.frame(module="CSR_uniformization",states_or_atoms=J,
    fibres=NA_integer_,elapsed_seconds=unname(elapsed),input_bytes=as.numeric(object.size(q))+as.numeric(object.size(probability)),
    native_workspace_estimate=y$workspace_bytes,
    endpoint_or_probability=y$probability[1],certificate=y$error_bound)
}
out <- do.call(rbind,results)
print(out, row.names=FALSE)
args <- commandArgs(trailingOnly=TRUE)
if(length(args))write.csv(out,args[1],row.names=FALSE)
cat("\nRuntime and memory are measured on this machine. Native workspace excludes R copies.\n")
print(sessionInfo())
