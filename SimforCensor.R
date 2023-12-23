#===================================================================
# Demo: Simulation to define cs.ita (censor time ~ Uniform(0,cs.ita))
#-------------------------------------------------------------------
# Written by:
#              Lican Kang (kanglican@whu.edu.cn)
#              Ning  SU (suning@whu.edu.cn)
# This version: Feb 17, 2023 
#===================================================================
gc(rm(list=ls())) 
# parameter of generating data 
T0 = 15
n = 300
p = 5000
rho = 0.2
sigma = 1
R = 10
seednum = 123
kind = 1
# parameter of simulation to determine cs.ita
sim = 100
#------------------------------------------------------------------------
# Generate beta
source("GenData.R")
beta = get.beta(p, T0, R, seednum)
betatrue = beta[[1]]
supptrue = beta[[2]]

# Simulation to determine parameter cs.ita
cs.ita = 2e+6
c.r = censor.rate(cs.ita,n,p,betatrue,sigma,rho,kind,sim)
c.r




