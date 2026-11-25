library(coherent4cast)
stopifnot(as.character(packageVersion('coherent4cast')) == '0.10.2')
stopifnot(identical(coherent4cast::coherent4cast, coherent4cast::recurcast))
stopifnot(identical(coherent4cast::coherent4cast_control, coherent4cast::recurcast_control))
stopifnot('coherent4cast' %in% names(getLoadedDLLs()))
stopifnot(requireNamespace('recurcast', quietly=TRUE))
old_d <- recurcast::rc_simulate(80, horizon=1, rates=c(.9,.5,.4), dispersion=0, seed=701)
new_d <- coherent4cast::rc_simulate(80, horizon=1, rates=c(.9,.5,.4), dispersion=0, seed=701)
stopifnot(identical(old_d,new_d))
old_fit <- recurcast::recurcast(~x, old_d, time_unit='years', frailty='none', event_types=c('1','2'))
new_fit <- coherent4cast::coherent4cast(~x, new_d, time_unit='years', frailty='none', event_types=c('1','2'))
coefficient_error <- max(abs(coef(old_fit)-coef(new_fit)))
covariance_error <- max(abs(vcov(old_fit)-vcov(new_fit)))
old_prediction <- predict(old_fit, newdata=old_d[old_d$id==1,], landmark=0, horizon=c(.25,.5,1), type='all')
new_prediction <- predict(new_fit, newdata=new_d[new_d$id==1,], landmark=0, horizon=c(.25,.5,1), type='all')
numeric_columns <- names(old_prediction)[vapply(old_prediction,is.numeric,logical(1))]
stopifnot(identical(names(old_prediction), names(new_prediction)))
forecast_error <- max(abs(as.matrix(old_prediction[numeric_columns])-as.matrix(new_prediction[numeric_columns])),na.rm=TRUE)
stopifnot(coefficient_error<1e-12,covariance_error<1e-12,forecast_error<1e-12)
stopifnot(inherits(new_fit,'coherent4cast'))
report <- list(previous_version=as.character(packageVersion('recurcast')),
  current_version=as.character(packageVersion('coherent4cast')),
  generator_identical=TRUE, compatibility_aliases_identical=TRUE,
  coefficient_max_absolute_difference=coefficient_error,
  covariance_max_absolute_difference=covariance_error,
  forecast_max_absolute_difference=forecast_error,
  new_namespace_and_DLL_loaded=TRUE)
jsonlite::write_json(report,'rename_parity.json',pretty=TRUE,auto_unbox=TRUE)
print(report)



