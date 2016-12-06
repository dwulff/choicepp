
# require(Rcpp)
# # sourceCpp('~/Dropbox (2.0)/Work/Software/choicepp/src/environments.cpp')
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


# smpl_sws(c(1,0,1,0),c(2,-10,.8,.2),1,10,100)
# smpl_sure(c(2,1,.3,.7),c(2,-10,.99,.01),.1,.01)
#
#
# require(choicepp)
#
# p   = pgen_rnd(c(0,0,100),2,2,lower=-10,upper=10,ecological = T)
#
# p = matrix(rep(c(32,0,.025,.975,3,0,.25,.75),10000),ncol=8,byrow=T)
#
# par = p_arrange(p,2)
# sure = sampl_sure(par,.1,.01)
# swe  = sampl_swe(par,1)
#
#
# d = lapply(sure,function(x) (which(x[[1]] == 32)-1) / (length(x[[1]]-1)))
# hist(unlist(d),breaks=40)
#
# d = lapply(sure,function(x) mean(x[[1]]==32))
# mean(unlist(d))
#
# d = lapply(swe,function(x) mean(x[[1]]==32))
# mean(unlist(d))
#
# mean(unlist(d)==32)
#
# d = lapply(swe,function(x) x[[1]])
# mean(unlist(d)==32)
#
#
#
