# returns an ExpressionSet, containing information about classification in featureData
# 'E_set' is an ExpressionSet
# 't_pi_prior' is the initial value for 't_pi'
# 'k_prior' is the intial value for 'k'
# 'b' is the concentration parameter of Dirichlet distribution
# 'plot' is a boolean indicating if plot frequency and approximated gamma distribution

# G is the number of genes,
# n is the number of test samples for every genes

lnn <- function(E_Set, t_pi_prior = c(0.33, 0.33, 0.34), k_prior = 0.1, b = c(2,2,2), plot = 0, is_sim = 0, verbose = 0)
{
  G = nrow(assayData(E_Set)[['exprs']])
  n = ncol(assayData(E_Set)[['exprs']])

  data = matrix(assayData(E_Set)[['exprs']], nrow = G, ncol = n)

  # G = nrow(data)
  # n = ncol(data)

  infinity = 1e100

  # 'sum_dgl_by_l' is an G * 1 matrix, the summation result of every row of 'data'
  sum_dgl_by_l = apply(data,1,sum)

  # 'sum_dgl_square_by_l' is an G * 1 matrix, the summation of every squared elements of 'data' by row
  sum_dgl_square_by_l = apply(data^2,1,sum)

  # lf123 returns log likelihood log p_1 log p_2 log p_3,
  # specific return value depends on the last parameter 'category'
  # 'psi' is the parameters, 'psi' = (\delta_1, \delta_2, k, \lambda, \nu)
  # 'G_index' is a vector, it should be 1, 2, 3, ..., G
  lf123 <- function(psi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, category)
  {
    delta_1 = psi[1]
    delta_2 = psi[2]
    k = psi[3]
    lambda = psi[4]
    nu = psi[5]

    alpha = exp(psi[4])
    beta = exp(psi[5])

    A = n/(2*(n*k+1))
    B_g = sum_dgl_by_l[G_index]/n
    C_g = beta + sum_dgl_square_by_l[G_index]/2 - (sum_dgl_by_l[G_index])^2/(2*n)
    D = alpha * log(beta) + lgamma(n/2+alpha) - lgamma(alpha) - n/2 * log(2*pi) - log(n*k+1)/2

    if (category == 1)
    {
      return (D - (n/2+alpha) * log (C_g + A * (exp(delta_1) - B_g)^2))
    } else if (category == 2)
    {
      return (D - (n/2+alpha) * log (C_g + A * (exp(delta_2) + B_g)^2))
    } else
    {
      return (D - (n/2+alpha) * log (C_g + A * (B_g)^2))
    }
  }

  # it returns the expectation of \z_{g1}, \z_{g2}, \z_{g3}
  # this function is used in E-step
  # 'psi' is the parameters, 'psi' = (\delta_1, \delta_2, k, \lambda, \nu)
  # 't_pi' = (\pi_1, \pi_2, \pi_3)
  # 'G_index' is a vector, it should be 1, 2, 3, ..., G
  get_tilde_z <- function(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n)
  {
    t1 = t_pi[1] * exp(lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, 1))
    t2 = t_pi[2] * exp(lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, 2))
    t3 = t_pi[3] * exp(lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, 3))
    total = t1 + t2 + t3

    return (matrix(c(t1/total, t2/total, t3/total), nrow = length(G_index), ncol = 3))
  }

  # returns log likelihood in question
  # this is the objective function to be maximized
  l_c <- function(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)
  {
    if (t_pi[1]<1e-6 || t_pi[2]<1e-6 || t_pi[3]<1e-6) return (-infinity)
    result = sum(tilde_z[,1] * lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, 1)) + sum(tilde_z[,2] * lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, 2)) + sum(tilde_z[,3] * lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, 3)) + sum(tilde_z[,1] * log(t_pi[1])) + sum(tilde_z[,2] * log(t_pi[2])) + sum(tilde_z[,3] * log(t_pi[3]))

    #with dirichlet distribution
    result = result + lgamma(b[1]+b[2]+b[3]) - lgamma(b[1]) - lgamma(b[2]) - lgamma(b[3])
    result = result + sum((b-1) * log(t_pi))

    return (result)
  }

  # returns negative log likelihood in question
  # maximizing over l_c is equivalent to minimizing over negative_l_c
  negative_l_c <- function(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)
  {
    return (-l_c(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z))
  }

  # returns gradient w.r.t. 'psi' = (\delta_1, \delta_2, k, \lambda, \nu)
  gradient_l_c <- function(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)
  {
    delta_1 = psi[1]
    delta_2 = psi[2]
    k = psi[3]
    lambda = psi[4]
    nu = psi[5]

    alpha = exp(psi[4])
    beta = exp(psi[5])

    G = length(G_index)

    A = n/(2*(n*k+1))
    B_g = sum_dgl_by_l[G_index]/n
    C_g = beta + sum_dgl_square_by_l[G_index]/2 - (sum_dgl_by_l[G_index])^2/(2*n)
    D = alpha * log(beta) + lgamma(n/2+alpha) - lgamma(alpha) -n/2 * log(2*pi) - log(n*k+1)/2

    d_delta_1 = - sum( tilde_z[,1] * (2*A*(exp(delta_1) - B_g))/(A*(exp(delta_1) - B_g)^2 + C_g)) * (n/2+alpha) * exp(delta_1)

    d_delta_2 = - sum( tilde_z[,2] * (2*A*(exp(delta_2) + B_g))/(A*(exp(delta_2) + B_g)^2 + C_g)) * (n/2+alpha) * exp(delta_2)

    d_k = - (n * G)/(2*(n*k+1)) + (n^2)/(2*(n*k+1)^2) * (n/2 + alpha) * (sum(tilde_z[,1] * (exp(delta_1) - B_g)^2 / (A * (exp(delta_1) - B_g)^2 + C_g)) + sum(tilde_z[,2] * (exp(delta_2) + B_g)^2 / (A * (exp(delta_2) + B_g)^2 + C_g)) + sum(tilde_z[,3] * (B_g)^2 / (A * (B_g)^2 + C_g)))

    d_lambda = G * exp(lambda) * (nu + digamma(n/2+alpha) - digamma(alpha)) - exp(lambda) * ( sum(tilde_z[,1] * log(A * (exp(delta_1) - B_g)^2 + C_g)) + sum(tilde_z[,2] * log(A * (exp(delta_2) + B_g)^2 + C_g)) + sum(tilde_z[,3] * log(A * (B_g)^2 + C_g)))

    d_nu = G * exp(lambda) - exp(nu) * (n/2 + exp(lambda)) * (sum(tilde_z[,1] / (A * (exp(delta_1) - B_g)^2 + C_g)) + sum(tilde_z[,2] / (A * (exp(delta_2) + B_g)^2 + C_g)) + sum(tilde_z[,3] / (A * (B_g)^2 + C_g)))

    return (c(d_delta_1, d_delta_2, d_k, d_lambda, d_nu))
  }

  # returns gradient of nagative log likelihood function
  gradient_negative_l_c <- function(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)
  {
    return (-gradient_l_c(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z))
  }

  get_t_pi <- function(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)
  {
    #without Dirichlet
    # t1 = sum(tilde_z[,1]) / length(G_index)
    # t2 = sum(tilde_z[,2]) / length(G_index)
    # t3 = sum(tilde_z[,3]) / length(G_index)

    #with Dirichlet distribution
    denominator = length(G_index) + sum(b) - 3

    # print("===========================")
    # print(sum(tilde_z[,1]))
    # print(sum(tilde_z[,2]))
    # print(sum(tilde_z[,3]))
    
    t1 = (sum(tilde_z[,1]) + b[1] - 1) / denominator
    t2 = (sum(tilde_z[,2]) + b[2] - 1) / denominator
    t3 = (sum(tilde_z[,3]) + b[3] - 1) / denominator

    return (c(t1,t2,t3))
  }

  # limitations of parameters
  delta_1_min = -6
  delta_1_max = 15

  delta_2_min = delta_1_min
  delta_2_max = 15

  k_min = 0.001
  k_max = 1/k_min

  lambda_min = -6
  lambda_max = 6

  nu_min = -6
  nu_max = 6

  # initial values of 'psi' = (\delta_1, \delta_2, k, \lambda, \nu)

  median_dgl_by_l = apply(data, 1, median)

  temp = median(sort(median_dgl_by_l)[(ceiling(G * 0.95)):G])
  if (temp>0) delta_1 = log(temp)
  else delta_1 = delta_1_min

  temp = median(sort(median_dgl_by_l)[1:(trunc(G * 0.05))])
  if (temp<0) delta_2 = log(-temp)
  else delta_2 = delta_2_min

  temp_tau = 1 / (apply(data, 1, mad)^2)
  omega = mad(median_dgl_by_l)^2
  k_prior = omega * median(temp_tau)

  k = k_prior


  lambda = log((median(temp_tau)^2)/(mad(temp_tau)^2))
  nu = log(median(temp_tau)/(mad(temp_tau)^2))

  if (plot)
  {
    hist(log10(temp_tau))
  }

  psi = c(delta_1, delta_2, k, lambda, nu)

  # print(psi)
  t_pi = t_pi_prior
  G_index = c(1:G)

  tilde_z = get_tilde_z(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n)
  t_pi = get_t_pi(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)

  ###############################################################################
  ####-----------For Test Use Only-------------------------------------------####
  ###############################################################################

  # focus = 5
  # h = 0.0001

  # psi[focus] = psi[focus] - h
  # f1 = l_c(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)
  # psi[focus] = psi[focus] + h
  # f2 = l_c(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)
  # psi[focus] = psi[focus] + h
  # f3 = l_c(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)
  # psi[focus] = psi[focus] - h

  # print(c("f1,f2,f3",f1,f2,f3))

  # derivative_approx = (f3 - f1)/(2*h)
  # derivative_accurate = gradient_l_c(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)[focus]

  # print(c("difference:", derivative_accurate - derivative_approx))

  # return(0)

  ###############################################################################
  ####----END----For Test Use Only----END------------------------------------####
  ###############################################################################

  mleinfo = optim(par = psi, fn = negative_l_c, gr = gradient_negative_l_c, 
    t_pi = t_pi, sum_dgl_by_l = sum_dgl_by_l, sum_dgl_square_by_l = sum_dgl_square_by_l, 
    G_index = G_index, n = n, tilde_z = tilde_z, method = 'L-BFGS-B', 
    lower=c(delta_1_min, delta_2_min, k_min, lambda_min, nu_min), 
    upper=c(delta_1_max, delta_2_max, k_max, lambda_max, nu_max),
    control = list(maxit = 10000))

  repeated_times = 0
  repeated_times_max = 1000

  while (repeated_times<repeated_times_max)
  {
    repeated_times = repeated_times + 1
    if (verbose)
    {
      print(c("repeated times:", repeated_times))
    }

    psi = mleinfo$par
    
    tilde_z = get_tilde_z(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n)

    t_pi = get_t_pi(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)

    if (verbose)
    {
      print(psi)
      print(t_pi)
      # print(lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, 1)[17366])
      # print(lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, 2)[17366])
      # print(lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, 3)[17366])
      # print(tilde_z[17366,])
    #   print(negative_l_c(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z))
    }

    last_mleinfo = mleinfo
    mleinfo = optim(par = psi, fn = negative_l_c, gr = gradient_negative_l_c, 
      t_pi = t_pi, sum_dgl_by_l = sum_dgl_by_l, sum_dgl_square_by_l = sum_dgl_square_by_l, 
      G_index = G_index, n = n, tilde_z = tilde_z, method = 'L-BFGS-B', 
      lower=c(delta_1_min, delta_2_min, k_min, lambda_min, nu_min), 
      upper=c(delta_1_max, delta_2_max, k_max, lambda_max, nu_max),
      control = list(maxit = 10000))

    if (abs(last_mleinfo$value - mleinfo$value)<1e-6) break
  }

  if (plot)
  {
    x = seq(0.01, max(temp_tau), length=1000)
    lines(x, dgamma(x,exp(mleinfo$par[4]),exp(mleinfo$par[5]))*G, col = 'green')
  }

  tilde_z = get_tilde_z(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n)

  # 3 means non differential expression
  # 1 means over expression
  # 2 means under expression

  fData(E_Set)[,2] = apply(tilde_z,1,which.max)
  if (is_sim)
  {
    fData(E_Set)[,3] = (fData(E_Set)[,2] == fData(E_Set)[,1])
  }
  # write.csv(tilde_z, file = "tilde_z.csv")

  return (list(E_Set, mleinfo, t_pi))

}
