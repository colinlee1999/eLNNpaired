library("Biobase")
# library("gtools")
library("clues")

source('gen_data_by_eLNNpaired_cluster_wise_limma_prior.r')
source('eLNNpaired_cluster_wise_limma_prior.r')
source('rough_analyze.r')

G = 1000
n = 100
times = 1

delta_1 = 1
k_1 = exp(1)
lambda_1 = 3
nu_1 = 1

delta_2 = 1
k_2 = exp(1)
lambda_2 = 3
nu_2 = 1

k_3 = exp(1)
lambda_3 = 3
nu_3 = 1

psi = c(delta_1, k_1, lambda_1, nu_1,
        delta_2, k_2, lambda_2, nu_2,
        0,       k_3, lambda_3, nu_3)
b = c(2,2,2)
t_pi = c(0.05, 0.05, 0.90)


es = gen_data_by_eLNNpaired_cluster_wise_limma_prior(G,n,psi,t_pi)
result = eLNNpaired_cluster_wise_limma_prior(es, verbose = 1, is_sim =1)

print(table(fData(result$E_Set)$est_cluster))