plot_lotteries = function(a,b) {
  oa = c(a[1:(length(a)/2)])
  pa = c(a[(length(a)/2+1):length(a)])[order(oa)]
  ob = c(b[1:(length(b)/2)])
  pb = c(b[(length(b)/2+1):length(b)])[order(ob)]
  print(sum(pb))
  oa = oa[order(oa)]
  ob = ob[order(ob)]
  cpa = cumsum(pa)
  cpb = cumsum(pb)
  plot.new();plot.window(ylim=c(0,1),xlim=c(min(c(oa,ob)),max(c(oa,ob))))
  lines(c(oa[1],oa[1]),c(0,cpa[1]),col='red',lwd=3)
  for(i in 1:(length(a)-1)){
    lines(c(oa[i+1],oa[i+1]),c(cpa[i],cpa[i+1]),col='red',lwd=3)
    lines(c(oa[i],oa[i+1]),c(cpa[i],cpa[i]),col='red',lwd=3)
  }
  lines(c(ob[1],ob[1]),c(0,cpb[1]),col='blue',lwd=3)
  for(i in 1:(length(b)-1)){
    lines(c(ob[i+1],ob[i+1]),c(cpb[i],cpb[i+1]),col='blue',lwd=3)
    lines(c(ob[i],ob[i+1]),c(cpb[i],cpb[i]),col='blue',lwd=3)
  }
  points(oa,cpa,pch=15)
  points(ob,cpb,pch=15)
  axis(1);axis(2)
}
