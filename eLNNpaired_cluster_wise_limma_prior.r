# *** this version is for cluster-wise eLNN paired ***

# returns an ExpressionSet, containing information about classification in featureData
# 'E_set' is an ExpressionSet
# 't_pi_prior' is the initial value for 't_pi'
# 'k_prior' is the intial value for 'k'
# 'b' is the concentration parameter of Dirichlet distribution
# 'plot' is a boolean indicating if plot frequency and approximated gamma distribution

# G is the number of genes,
# n is the number of test samples for every genes

library(iCheck)

eLNNpaired_cluster_wise_limma_prior <- function(
  E_Set,
  b = c(2,2,2), 
  t_pi_prior = c(0.05, 0.05, 0.90),
  plot = 0, 
  is_sim = 0, 
  verbose = 0,
  infinity = 1e100,
  converge_threshold = 1e-6,
  param_limit_min = c(-6,0.01,-6,-6,-6,0.01,-6,-6,0,0.01,-6,-6),
  param_limit_max = c(6,100,6,6,6,100,6,6,0,100,6,6),
  max_iteration_num_in_optim = 100,
  max_repeated_times = 1000
  )
{
  G = nrow(E_Set)
  n = ncol(E_Set)

  data_matrix_of_E_Set = exprs(E_Set)

  # 'sum_dgl_by_l' is an G * 1 matrix, the summation result of every row of 'data_matrix_of_E_Set'
  sum_dgl_by_l = apply(data_matrix_of_E_Set,1,sum, na.rm=TRUE)

  # 'sum_dgl_square_by_l' is an G * 1 matrix, the summation of every squared elements of 'data_matrix_of_E_Set' by row
  sum_dgl_square_by_l = apply(data_matrix_of_E_Set^2,1,sum, na.rm = TRUE)

  # initial values of 'psi' = (\delta_1, \delta_2, k, \lambda, \nu)
  # lf123 returns log likelihood log p_1 log p_2 log p_3,
  # specific return value depends on the last parameter 'category'
  # 'psi' is the parameters, 'psi' = (\delta_1, \delta_2, k, \lambda, \nu)
  lf123 <- function(
    psi, 
    sum_dgl_by_l, 
    sum_dgl_square_by_l, 
    n)
  {
    result = matrix(rep(0,length(sum_dgl_by_l)*3),length(sum_dgl_by_l),3)
    colnames(result) = c("logf1", "logf2", "logf3")

    for (category in 1:3)
    {
      mu_0 = exp(psi[category*4-3])

      if (category == 2) mu_0 = -mu_0
      if (category == 3) mu_0 = 0

      k = exp(psi[category*4-2])
      alpha = exp(psi[category*4-1])
      beta = exp(psi[category*4])
    
      A = n/(2*(n*k+1))
      B_g = sum_dgl_by_l/n
      C_g = beta + sum_dgl_square_by_l/2 - (sum_dgl_by_l)^2/(2*n)
      D = alpha * log(beta) + lgamma(n/2+alpha) - lgamma(alpha) - n/2 * log(2*pi) - log(n*k+1)/2
    
      result[,category] = D - (n/2+alpha) * log (C_g + A * (mu_0 - B_g)^2)
    }

    return (result)
  }

  # it returns the expectation of \z_{g1}, \z_{g2}, \z_{g3}
  # this function is used in E-step
  # 'psi' is the parameters, 'psi' = (\delta_1, \delta_2, k, \lambda, \nu)
  # 't_pi' = (\pi_1, \pi_2, \pi_3)
  get_tilde_z <- function(
    psi, 
    t_pi, 
    sum_dgl_by_l, 
    sum_dgl_square_by_l, 
    n)
  {
    logf = lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, n)
    max_logf = apply(logf, 1, max, na.rm = TRUE)
    t1 = t_pi[1] * exp(logf[,1] - max_logf)
    t2 = t_pi[2] * exp(logf[,2] - max_logf)
    t3 = t_pi[3] * exp(logf[,3] - max_logf)
    total = t1 + t2 + t3

    result = cbind(t1, t2, t3)/total
    return(result)
  }

  # returns expected log likelihood in question
  # this is the objective function to be maximized
  l_c <- function(
    psi, 
    t_pi, 
    sum_dgl_by_l, 
    sum_dgl_square_by_l, 
    n, 
    tilde_z, 
    b, 
    converge_threshold, 
    infinity)
  {
    if (t_pi[1]<converge_threshold && t_pi[2]<converge_threshold) return (-infinity)
    if (t_pi[1]<converge_threshold && t_pi[3]<converge_threshold) return (-infinity)
    if (t_pi[2]<converge_threshold && t_pi[3]<converge_threshold) return (-infinity)

    logf = lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, n)

    result = 0
    result = result + sum(tilde_z[,1] * logf[,1], na.rm = TRUE)
    result = result + sum(tilde_z[,2] * logf[,2], na.rm = TRUE)
    result = result + sum(tilde_z[,3] * logf[,3], na.rm = TRUE)
    result = result + sum(tilde_z[,1] * log(t_pi[1]), na.rm = TRUE)
    result = result + sum(tilde_z[,2] * log(t_pi[2]), na.rm = TRUE)
    result = result + sum(tilde_z[,3] * log(t_pi[3]), na.rm = TRUE)

    #with dirichlet distribution
    result = result + lgamma(b[1]+b[2]+b[3]) - lgamma(b[1]) - lgamma(b[2]) - lgamma(b[3])
    result = result + sum((b-1) * log(t_pi), na.rm = TRUE)

    return (result)
  }

  # returns negative log likelihood in question
  # maximizing over l_c is equivalent to minimizing over negative_l_c
  negative_l_c <- function(
    psi, 
    t_pi, 
    sum_dgl_by_l, 
    sum_dgl_square_by_l, 
    n, 
    tilde_z, 
    b, 
    converge_threshold, 
    infinity)
  {
    return (-l_c(
      psi, 
      t_pi, 
      sum_dgl_by_l, 
      sum_dgl_square_by_l, 
      n, 
      tilde_z, 
      b, 
      converge_threshold, 
      infinity))
  }

  # returns gradient w.r.t. 'psi' = (\delta_1, \delta_2, k, \lambda, \nu)
  gradient_l_c <- function(
    psi, 
    t_pi, 
    sum_dgl_by_l, 
    sum_dgl_square_by_l, 
    n, 
    tilde_z, 
    b, 
    converge_threshold, 
    infinity)
  {
    G = length(sum_dgl_by_l)
    gradient_temp = matrix(rep(0,12),4,3)

    for (category in 1:3)
    {
      mu_0 = exp(psi[category*4-3])
        k = exp(psi[category*4-2])
        
        lambda = psi[category*4-1]
        nu = psi[category*4]

        alpha = exp(psi[category*4-1])
        beta = exp(psi[category*4])


        A = n/(2*(n*k+1))
        B_g = sum_dgl_by_l/n
        C_g = beta + sum_dgl_square_by_l/2 - (sum_dgl_by_l)^2/(2*n)
        D = alpha * log(beta) + lgamma(n/2+alpha) - lgamma(alpha) -n/2 * log(2*pi) - log(n*k+1)/2
        
        if (category == 2) mu_0 = - mu_0
        if (category == 3) mu_0 = 0

        d_delta = - sum(tilde_z[,category] * (2*A*(mu_0 - B_g))/(A*(mu_0 - B_g)^2 + C_g) * (n/2 + alpha) * mu_0)

        if (category == 3)
        {
          d_delta = 0
        }

        d_xi = k * ( - (n * sum(tilde_z[,category])) / (2*(n*k+1)) + (n^2)/(2*(n*k+1)^2) * (n/2+alpha) *
                sum(tilde_z[,category]*(mu_0-B_g)^2/(A*(mu_0-B_g)^2+C_g)))
        
        d_lambda = sum(tilde_z[,category] * alpha * (nu + digamma(n/2+alpha) - digamma(alpha))) - alpha * ( sum(tilde_z[,category] * log(A * (mu_0 - B_g)^2 + C_g)) )
        
        d_nu = sum(tilde_z[,category]) * alpha - beta * (n/2 + alpha) * (sum(tilde_z[,category] / (A * (mu_0 - B_g)^2 + C_g)))
        
        gradient_temp[,category] = c(d_delta, d_xi, d_lambda, d_nu)
    }

    result = c(gradient_temp)

    return (result)
  }

  # returns gradient of nagative log likelihood function
  gradient_negative_l_c <- function(
    psi, 
    t_pi, 
    sum_dgl_by_l, 
    sum_dgl_square_by_l, 
    n, 
    tilde_z, 
    b, 
    converge_threshold, 
    infinity)
  {
    return (-gradient_l_c(
      psi, 
      t_pi, 
      sum_dgl_by_l, 
      sum_dgl_square_by_l, 
      n, 
      tilde_z, 
      b, 
      converge_threshold, 
      infinity))
  }

  get_t_pi <- function(
    psi, 
    t_pi, 
    sum_dgl_by_l, 
    sum_dgl_square_by_l, 
    n, 
    tilde_z, 
    b,
    converge_threshold,
    infinity)
  {
    #with Dirichlet distribution
    denominator = length(sum_dgl_by_l) + sum(b, na.rm=TRUE) - 3

    t1 = (sum(tilde_z[,1], na.rm=TRUE) + b[1] - 1) / denominator
    t2 = (sum(tilde_z[,2], na.rm=TRUE) + b[2] - 1) / denominator
    t3 = (sum(tilde_z[,3], na.rm=TRUE) + b[3] - 1) / denominator

    return (c(t1,t2,t3))
  }

  column_names = colnames(fData(E_Set))

  result_limma = lmFitPaired(
    E_Set,
    probeID.var = column_names[1],
    gene.var = column_names[2],
    chr.var = column_names[3],
    verbose = FALSE)

  frame_unsorted = result_limma$frame.unsorted
  over_expressed_sub_script = frame_unsorted$pos[which(frame_unsorted$stats > 0 & frame_unsorted$p.adj < 0.05)]
  under_expressed_sub_script = frame_unsorted$pos[which(frame_unsorted$stats < 0 & frame_unsorted$p.adj < 0.05)]
  non_diff_sub_script = frame_unsorted$pos[which(frame_unsorted$p.adj >= 0.05)]

  t_pi_prior = c(
    length(under_expressed_sub_script)/G,
    length(over_expressed_sub_script)/G,
    0)
  t_pi_prior[3] = 1 - t_pi_prior[1] - t_pi_prior[2]

  over_expressed_E_Set = E_Set[over_expressed_sub_script]
  under_expressed_E_Set = E_Set[under_expressed_sub_script]
  non_diff_E_Set = E_Set[non_diff_sub_script]
  
  #########################################
  # cluster 1
  data_matrix_of_E_Set = exprs(over_expressed_E_Set)

  median_dgl_by_l = apply(data_matrix_of_E_Set, 1, median, na.rm=TRUE)

  sorted_median_dgl_by_l = sort(median_dgl_by_l)

  temp = median(sorted_median_dgl_by_l, na.rm=TRUE)
  if (temp>0) delta_1 = log(temp)
  else delta_1 = param_limit_min[1]

  temp_tau = 1 / (apply(data_matrix_of_E_Set, 1, mad, na.rm=TRUE)^2)
  omega = mad(median_dgl_by_l)^2
  k_prior = omega * median(temp_tau, na.rm=TRUE)

  xi_1 = log(k_prior)

  lambda_1 = log((median(temp_tau, na.rm=TRUE)^2)/(mad(temp_tau, na.rm=TRUE)^2))
  nu_1 = log(median(temp_tau, na.rm=TRUE)/(mad(temp_tau, na.rm=TRUE)^2))
  # end of cluster 1
  #########################################

  #########################################
  # cluster 2
  data_matrix_of_E_Set = exprs(under_expressed_E_Set)

  median_dgl_by_l = apply(data_matrix_of_E_Set, 1, median, na.rm=TRUE)

  sorted_median_dgl_by_l = sort(median_dgl_by_l)

  temp = median(sorted_median_dgl_by_l, na.rm=TRUE)
  if (temp>0) delta_2 = log(temp)
  else delta_2 = param_limit_min[5]

  temp_tau = 1 / (apply(data_matrix_of_E_Set, 1, mad, na.rm=TRUE)^2)
  omega = mad(median_dgl_by_l)^2
  k_prior = omega * median(temp_tau, na.rm=TRUE)

  xi_2 = log(k_prior)

  lambda_2 = log((median(temp_tau, na.rm=TRUE)^2)/(mad(temp_tau, na.rm=TRUE)^2))
  nu_2 = log(median(temp_tau, na.rm=TRUE)/(mad(temp_tau, na.rm=TRUE)^2))
  # end of cluster 2
  #########################################

  #########################################
  # cluster 3
  data_matrix_of_E_Set = exprs(non_diff_E_Set)

  median_dgl_by_l = apply(data_matrix_of_E_Set, 1, median, na.rm=TRUE)

  sorted_median_dgl_by_l = sort(median_dgl_by_l)

  temp_tau = 1 / (apply(data_matrix_of_E_Set, 1, mad, na.rm=TRUE)^2)
  omega = mad(median_dgl_by_l)^2
  k_prior = omega * median(temp_tau, na.rm=TRUE)

  xi_3 = log(k_prior)

  lambda_3 = log((median(temp_tau, na.rm=TRUE)^2)/(mad(temp_tau, na.rm=TRUE)^2))
  nu_3 = log(median(temp_tau, na.rm=TRUE)/(mad(temp_tau, na.rm=TRUE)^2))
  # end of cluster 3
  #########################################

  data_matrix_of_E_Set = exprs(E_Set)

  if (plot)
  {
    hist(log10(temp_tau))
  }

  psi = c(delta_1, xi_1, lambda_1, nu_1,
          delta_2, xi_2, lambda_2, nu_2,
          0,       xi_3, lambda_3, nu_3)

  t_pi = t_pi_prior

  tilde_z = get_tilde_z(
    psi, 
    t_pi, 
    sum_dgl_by_l, 
    sum_dgl_square_by_l,
    n)
  
  mleinfo = optim(par = psi, fn = negative_l_c, gr = gradient_negative_l_c, 
    t_pi = t_pi, sum_dgl_by_l = sum_dgl_by_l, sum_dgl_square_by_l = sum_dgl_square_by_l, 
    n = n, tilde_z = tilde_z, method = 'L-BFGS-B', 
    b = b, converge_threshold = converge_threshold, infinity = infinity,
    # lower=param_limit_min, 
    # upper=param_limit_max,
    control = list(maxit = max_iteration_num_in_optim))

  t_pi = get_t_pi(
    psi, 
    t_pi, 
    sum_dgl_by_l, 
    sum_dgl_square_by_l, 
    n, 
    tilde_z,
    b,
    converge_threshold,
    infinity)

  repeated_times = 0  

  while (repeated_times<max_repeated_times)
  {
    repeated_times = repeated_times + 1
    if (verbose)
    {
      print(c("repeated times:", repeated_times))
    }

    psi = mleinfo$par
    
    tilde_z = get_tilde_z(
      psi, 
      t_pi, 
      sum_dgl_by_l, 
      sum_dgl_square_by_l, 
      n)

    if (verbose)
    {
      print(psi)
      print(t_pi)
    }

    last_mleinfo = mleinfo
    mleinfo = optim(par = psi, fn = negative_l_c, gr = gradient_negative_l_c, 
      t_pi = t_pi, sum_dgl_by_l = sum_dgl_by_l, sum_dgl_square_by_l = sum_dgl_square_by_l, 
      n = n, tilde_z = tilde_z, method = 'L-BFGS-B', 
      b = b, converge_threshold = converge_threshold, infinity = infinity,
      # lower=param_limit_min, 
      # upper=param_limit_max,
      control = list(maxit = max_iteration_num_in_optim))

    t_pi = get_t_pi(
      psi, 
      t_pi, 
      sum_dgl_by_l, 
      sum_dgl_square_by_l, 
      n, 
      tilde_z,
      b,
      converge_threshold,
      infinity)

    #if (abs(last_mleinfo$value - mleinfo$value)<converge_threshold) break
    if (sum(abs(last_mleinfo$par - mleinfo$par))<converge_threshold) break
  }

  tilde_z = get_tilde_z(
    psi, 
    t_pi, 
    sum_dgl_by_l, 
    sum_dgl_square_by_l, 
    n)

  fData(E_Set)$est_cluster = apply(tilde_z,1,which.max)
  if (is_sim)
  {
    fData(E_Set)$flag = (fData(E_Set)$est_cluster == fData(E_Set)$true_cluster)
  }

  result = list(
    E_Set = E_Set, 
    mleinfo = mleinfo, 
    t_pi = t_pi)
  invisible(result)
}
