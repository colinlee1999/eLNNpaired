# generate an ExpressionSet object
# assayData[['exprs']] is an G by n matrix
# featureData describes if gene i is differentially expressed

gen_data_by_eLNNpaired <- function(G, n, psi, t_pi)
{

  delta_1 = psi[1]
  delta_2 = psi[2]
  k = psi[3]
  lambda = psi[4]
  nu = psi[5]

  alpha = exp(psi[4])
  beta = exp(psi[5])

  data = matrix(, nrow = G, ncol = n)
  category_info = matrix(rep(0,G*3),G,3)

  colnames(category_info) = c("true_cluster","est_cluster","flag")

  for (row in 1:G)
  {
    # determine which category, using a uniform distribution to compare with t_pi
  	temp = runif(1)
  	if (temp<t_pi[1]) category = 1
  	else if(temp<t_pi[1]+t_pi[2]) category = 2
  	else category = 3

    category_info[row,1] = category
  	
  	tau_g = rgamma(1, alpha, beta)

  	if (category==1) mu_g = rnorm(1, mean = exp(delta_1), sd = sqrt(k/tau_g))
  	else if (category==2) mu_g = rnorm(1, mean = -exp(delta_2), sd = sqrt(k/tau_g))
  	else if (category==3) mu_g = rnorm(1, mean = 0, sd = sqrt(k/tau_g))

  	data[row,] = rnorm(n,mean = mu_g, sd = sqrt(1/tau_g))
  }

  #return (data)
  return (ExpressionSet(assayData = data, featureData = new("AnnotatedDataFrame", data = data.frame(category_info))))
  
}