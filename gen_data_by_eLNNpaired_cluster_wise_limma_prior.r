# *** this version is for cluster-wise eLNN paired ***

# generate an ExpressionSet object
# assayData[['exprs']] is an G by n matrix
# featureData describes if gene i is differentially expressed

gen_data_by_eLNNpaired_cluster_wise_limma_prior <- function(G, n, psi, t_pi, c1, c2)
{
  data = matrix(, nrow = G, ncol = n)
  # t_matrix <<- matrix(, nrow = G, ncol = 1)
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

    if (category == 1)
      {
        mu_0 = exp(psi[category*4-3])
        k = pnorm(psi[category*4-2])
        beta = exp(psi[category*4])
        alpha = exp(psi[category*4-1]) + 1 + beta * c1^2 * (1+sqrt(k))^2/mu_0^2
      }
      else if (category == 2)
      {
        mu_0 = -exp(psi[category*4-3])
        k = pnorm(psi[category*4-2])
        beta = exp(psi[category*4])
        alpha = exp(psi[category*4-1]) + 1 + beta * c2^2 * (1+sqrt(k))^2/mu_0^2
      }
      else
      {
        mu_0 = 0
        k = pnorm(psi[category*4-2])
        beta = exp(psi[category*4])
        alpha = exp(psi[category*4-1])
      }

    tau_g = rgamma(1, alpha, beta)

    mu_g = rnorm(1, mean = mu_0, sd = sqrt(k/tau_g))

    # t_matrix[row,] <<- mu_g * sqrt(tau_g)

    data[row,] = rnorm(n,mean = mu_g, sd = sqrt(1/tau_g))
  }

  #return (data)
  return (ExpressionSet(assayData = data, featureData = new("AnnotatedDataFrame", data = data.frame(category_info)))) 
}