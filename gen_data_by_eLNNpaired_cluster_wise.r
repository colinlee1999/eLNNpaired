# *** this version is for cluster-wise eLNN paired ***

# generate an ExpressionSet object
# assayData[['exprs']] is an G by n matrix
# featureData describes if gene i is differentially expressed

gen_data_by_eLNNpaired_cluster_wise <- function(G, n, psi, t_pi)
{
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

    mu_0 = exp(psi[category*4-3])
    k = psi[category*4-2]
    alpha = exp(psi[category*4-1])
    beta = exp(psi[category*4])

    if (category == 1)
    {
      mu_0 = mu_0
    }
    else if (category ==2)
    {
      mu_0 = - mu_0
    }
    else
    {
      mu_0 = 0
    }
    
    tau_g = rgamma(1, alpha, beta)

    mu_g = rnorm(1, mean = mu_0, sd = sqrt(k/tau_g))

    data[row,] = rnorm(n,mean = mu_g, sd = sqrt(1/tau_g))
  }

  #return (data)
  return (ExpressionSet(assayData = data, featureData = new("AnnotatedDataFrame", data = data.frame(category_info)))) 
}