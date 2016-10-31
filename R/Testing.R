#
# require(Rcpp)
# sourceCpp('~/Dropbox (2.0)/Work/Software/choicepp/src/environments.cpp')
# sourceCpp('~/Dropbox (2.0)/Work/Software/choicepp/src/CPTcpp.cpp')
# sourceCpp('~/Dropbox (2.0)/Work/Software/choicepp/src/sample.cpp')
#
#
#
# prob = pgen_rnd(c(1000,1000,1000),2,2,ecological = FALSE,stdom1_test = TRUE, stdom2_test = TRUE,d_crit = .1)
# rprb = p_arrange(prob)
#
# ss = sampl_n(rprb,1000)
# exp = edit_exp(ss,rprb)
#
# rprb[[1]] - a[[1]]
#
# cpt_rndchoice(c(1,1,10),rprb,000)
#
#

#' @useDynLib choicepp
#' @importFrom Rcpp sourceCpp
NULL


#require(roxygen2)

