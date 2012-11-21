# anoles.R -- a unit test showing wrightscape performing the same fit as in the OUCH example

rm(list=ls())
require(ouch)

data(bimac)
tree <- with(bimac,ouchtree(node,ancestor,time/max(time),species))
ouch <- hansen(log(bimac['size']),tree,bimac['OU.LP'],sqrt.alpha=1,sigma=1)

require(wrightscape)
ws <- multiTypeOU(log(bimac[['size']]),tree,bimac[['OU.LP']], model_spec=list(alpha="global", sigma="global", theta="indep"))

