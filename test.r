library("Biobase")
# library("gtools")
library("clues")

source('gen_data_by_eLNNpaired.r')
source('eLNNpaired.r')
source('rough_analyze.r')

G = 10000
n = 100
times = 100

delta_1 = 1
k_1 = 2
lambda_1 = 3
nu_1 = 1

delta_2 = 1
k_2 = 2
lambda_2 = 3
nu_2 = 1

k_3 = 2
lambda_3 = 3
nu_3 = 1

psi = c(delta_1, k_1, lambda_1, nu_1,
        delta_2, k_2, lambda_2, nu_2,
        0,       k_3, lambda_3, nu_3)
b = c(2,2,2)
t_pi = c(0.05, 0.05, 0.90)


es = gen_data_by_eLNNpaired(G,n,psi,t_pi)
result = eLNNpaired(es, verbose = 1, is_sim =1, max_iteration_num_in_optim = 10000)