test_that("AUROC preserves distinct representable risk scores", {
  high <- .5 + 2*.Machine$double.eps
  extract_auc <- function(p, time, status) {
    x <- validate(p,time,status,horizon=1)
    x$metrics$estimate[x$metrics$metric=="AUROC"]
  }
  expect_equal(extract_auc(c(high,.5),c(.5,2),c(1,0)),1)
  expect_equal(extract_auc(c(.5,high),c(.5,2),c(1,0)),0)
  expect_equal(extract_auc(c(.5,.5),c(.5,2),c(1,0)),.5)
  x <- validate(c(high,.5),c(.5,2),c(1,0),horizon=1,
                censor_survival=list(case=c(.5,1),horizon=c(1,.25)))
  expect_equal(x$metrics$estimate[x$metrics$metric=="AUROC"],1)
})

test_that("no horizon outcome support cannot appear as perfect prediction", {
  expect_error(validate(c(.2,.3),c(.2,.4),c(0,0),horizon=1),
               "No observed cases or horizon survivors",class="recurcast_error")
})
