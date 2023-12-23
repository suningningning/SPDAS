#===================================================================
# Demo: solving AFT L0 regularizatio problem via SPDAS 
#-------------------------------------------------------------------
# Written by:
#              Lican Kang (kanglican@whu.edu.cn)
#              Ning  SU (suning@whu.edu.cn)
# This version: Dec 1, 2022 
#===================================================================
gc(rm(list=ls())) 
# parameter of generating data 
T0 = 15
n = 300
p = 5000
rho = 0.2
cs.ita = 2e+6 # for 0.25 average censor rate
sigma = 1
R = 10
seednum = 123
kind = 1
# parameter of SPDAS Algorithm
alpha = 0.85
m = 40
K = 20
#------------------------------------------------------------------------
# Generate beta
source("GenData.R")
beta = get.beta(p, T0, R, seednum)
betatrue = beta[[1]]
supptrue = beta[[2]]

# Generate data 
res = get.data(n, p, betatrue, sigma, rho, seednum, kind, cs.ita)
x0 = res[[1]]
y0 = res[[2]]
c0 = res[[3]]


# data processing
data = xy.process(x0, y0, c0, n)
x.tilde = data[[1]]
y.tilde = data[[2]]
D = data[[3]]


# SPDAS computing
source("S_PDAS.R")
t.start = Sys.time()
result = SPDAS(K, x.tilde, y.tilde, D, n, p, m, alpha)
t.end = Sys.time() 

betahatsupp = which(result[[1]]!=0)
AE = round(max(abs(result[[1]] - betatrue)),4)
RE = round(sqrt(sum((result[[1]] - betatrue)^2)/sum((betatrue)^2)),4)
MS = length(betahatsupp)
RP = all(sort(betahatsupp) == sort(supptrue))
time = t.end - t.start

AE 
RE
MS
RP 
time 

