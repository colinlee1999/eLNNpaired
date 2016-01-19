rm(list = ls())

library("Biobase")
# library("gtools")
library("clues")

source('gen_data_by_eLNNpaired_cluster_wise_limma_prior.r')
source('eLNNpaired_cluster_wise_limma_prior.r')
#source('rough_analyze.r')
library('iCheck')

start_time = proc.time()

# G = 10000
# n = 30
# times = 1

# delta_1 = 0
# xi_1 = 0
# lambda_1 = 3
# nu_1 = -1

# delta_2 = 0
# xi_2 = 0
# lambda_2 = 3
# nu_2 = -1

# xi_3 = -3
# lambda_3 = 3
# nu_3 = -1

G = 10000
n = 30
times = 1

delta_1 = -1.25
xi_1 = -0.79
lambda_1 = log(0.5)
nu_1 = -1.88

delta_2 = -1.26
xi_2 = -0.48
lambda_2 = log(0.5)
nu_2 = -2.37

xi_3 = -1.09
lambda_3 = 0.22
nu_3 = -2.70

psi = c(delta_1, xi_1, lambda_1, nu_1,
        delta_2, xi_2, lambda_2, nu_2,
        0,       xi_3, lambda_3, nu_3)
b = c(2,2,2)
t_pi = c(0.05, 0.05, 0.90)

c1 = qnorm(0.95)
c2 = qnorm(0.05)

generate = 1
if (generate)
{
	E_Set = gen_data_by_eLNNpaired_cluster_wise_limma_prior(G,n,psi,t_pi, c1, c2)
	saveRDS(E_Set, file = 'E_Set.Rdata')
}else
{
	E_Set = readRDS('E_Set.Rdata')
}

result_limma = lmFitPaired(
		E_Set, 
		probeID.var = "true_cluster", 
		gene.var = "est_cluster", 
		chr.var = "flag",
		verbose = FALSE)

result_limma$memGenes[which(result_limma$memGenes==2)] = 4
result_limma$memGenes[which(result_limma$memGenes==3)] = 2
result_limma$memGenes[which(result_limma$memGenes==4)] = 3

result = eLNNpaired_cluster_wise_limma_prior(E_Set, verbose = 1, is_sim = 1, c1 = c1, c2 = c2, max_repeated_times = 500)
result$time = proc.time() - start_time

print(table(fData(result$E_Set)$true_cluster,fData(result$E_Set)$est_cluster))
print(table(fData(result$E_Set)$true_cluster,result_limma$memGenes))
print(result$time)

source('draw_mu_g_tau_g.r')

# sub_script = intersect(which(result_limma$memGenes == 1),which(fData(result$E_Set)$true_cluster == 3))[1]

# dat = exprs(E_Set)
# hist(dat[sub_script,])

# fDat = fData(result$E_Set)
# postscript(file = 'temp.ps')
# if (length(which(fDat$flag==FALSE))>0)
# {
# 	temp_subscript = which(fDat$flag==FALSE)[1]
# 	hist(exprs(E_Set)[temp_subscript,])
# 	print(fDat[temp_subscript,])
# }
# dev.off()
