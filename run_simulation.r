# load(file = 'data.Rdata')

rm(list = ls())

start_time = proc.time()

library("Biobase")
library("clues")

source('gen_data_by_eLNNpaired_cluster_wise_limma_prior.r')
source('eLNNpaired_cluster_wise_limma_prior.r')
source('rough_analyze.r')

G = 10000
n = 30
times = 1

delta_1 = 0
xi_1 = 0
lambda_1 = 3
nu_1 = -1

delta_2 = 0
xi_2 = 0
lambda_2 = 3
nu_2 = -1

xi_3 = -3
lambda_3 = 3
nu_3 = -1

psi = c(delta_1, xi_1, lambda_1, nu_1,
        delta_2, xi_2, lambda_2, nu_2,
        0,       xi_3, lambda_3, nu_3)
b = c(2,2,2)
t_pi = c(0.05, 0.05, 0.90)

c1 = qnorm(0.95)
c2 = qnorm(0.05)



file_to_save_psi = 'psi.txt'
file_to_save_t_pi = 't_pi.txt'
file_to_save_accuracy = 'accuracy.txt'
file_to_save_notes = 'notes.txt'
file_to_save_adjusted_rand = 'adjusted_rand.txt'

# write to file_to_save_psi, file_to_save_t_pi, file_to_save_accuracy
write(toString(Sys.time()), file = file_to_save_t_pi)
write(toString(Sys.time()), file = file_to_save_psi)
write(toString(Sys.time()), file = file_to_save_accuracy)
write(c("total","npositive","npositive_recognize","npositive_not_recognize","false_positive","nnegative","nnegative_recognize","nnegative_not_recognize","false_negative"), file = file_to_save_accuracy, ncolumns = 9, append = TRUE, sep = '\t')



# write to file_to_save_notes
write(toString(Sys.time()), file = file_to_save_notes)
write(c('k_prior=', k_prior), file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('pi_prior=', t_pi_prior), file = file_to_save_notes, append = TRUE, sep = '\t')

write('', file = file_to_save_notes, append = TRUE, sep = '\t')
write('TRUE VALUES:', file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('delta_1=', delta_1), file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('delta_2=', delta_2), file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('k=', k), file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('lambda=', lambda), file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('nu=', nu), file = file_to_save_notes, append = TRUE, sep = '\t')

write('', file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('t_pi=', t_pi), file = file_to_save_notes, append = TRUE, sep = '\t')

write('', file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('G=', G), file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('n=', n), file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('times=', times), file = file_to_save_notes, append = TRUE, sep = '\t')

write('', file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('b=', b), file = file_to_save_notes, append = TRUE, sep = '\t')



# write to file_to_save_adjusted_rand
write(toString(Sys.time()), file = file_to_save_adjusted_rand)
write(c('Rand','HA','MA','FM','Jaccard'), file = file_to_save_adjusted_rand, ncolumns = 5, append = TRUE, sep = '\t')

for (i in 1:times)
{
  print(c("now running loop:", i))

  E_Set = gen_data(G, n, psi, t_pi)
  save.image(file = 'data.Rdata')  

  result = lnn(E_Set, t_pi_prior, k_prior, b, 0)

  write(c(result[[2]]$convergence,result[[2]]$par), file = file_to_save_psi, ncolumns = 6, append = TRUE, sep = '\t')
  write(result[[3]], file = file_to_save_t_pi, ncolumns = 3, append = TRUE, sep = '\t')

  rough_analysis = rough_analyze(as.matrix(fData(result[[1]])))
  write(rough_analysis, file = file_to_save_accuracy, ncolumns = length(rough_analysis), append = TRUE, sep = '\t')

  write(unname(adjustedRand(as.matrix(fData(result[[1]]))[,1],as.matrix(fData(result[[1]]))[,2])), file = file_to_save_adjusted_rand, ncolumns = 5, append = TRUE, sep = '\t')


  save.image(file = paste(c('simulation_',i,'.Rdata'), collapse = ''))
  
  print(proc.time()-start_time)
  #print(fData(result[[1]])[3098,])
  #print(sum(fData(result[[1]])[,3]))
}

end_time = proc.time()
print(end_time-start_time)
