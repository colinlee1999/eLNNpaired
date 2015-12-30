library("Biobase")
# library("gtools")
library("clues")

source('gen_data_by_eLNNpaired.r')
source('eLNNpaired.r')
source('rough_analyze.r')

G = 10000
n = 100
times = 1

delta_1 = 1
delta_2 = 1
k = 2
lambda = 3
nu = 1

psi = c(delta_1, delta_2, k, lambda, nu)
b = c(2,2,2)
t_pi = c(0.05, 0.05, 0.90)

es = gen_data_by_eLNNpaired(G,n,psi,t_pi)
result = eLNNpaired(es, verbose = 1, is_sim =1)