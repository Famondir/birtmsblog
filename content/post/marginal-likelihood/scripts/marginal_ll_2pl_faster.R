# uses matrix objects instead of tibbles to store data and vectorisation to calculate loglikelihood
mll_parallel_brms_2pl_secondapproach <- function(fit, MFUN, n_nodes = 11, best_only = FALSE, cores = 4) {

  # ----- create a temporary logging file ----
  tempDir <- tempfile()
  dir.create(tempDir)
  logFile <- file.path(tempDir, "log.txt")
  viewer <- getOption("viewer")
  viewer(logFile)

  # ----- initialise multiple workers ----
  cl <- parallel::makeCluster(cores, outfile = logFile)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl)) # terminate workes when finished

  # we get the dataset from the brms fit and add row and person numbers instead of passing it to the function
  data_list2 <- fit$data %>% mutate(resp_number = row_number(),
                                    person_number = as.numeric(factor(person, levels = unique(person))))

  draws2 <- list(sd = tidybayes::spread_draws(fit, sd_person__theta_Intercept) %>%
                   rename(sd_person = sd_person__theta_Intercept) %>%
                   select(!starts_with('.')) %>% as.matrix(),
                 theta = tidybayes::spread_draws(fit, r_person__theta[person,]) %>%
                   pivot_wider(values_from = r_person__theta, names_from = person) %>%
                   select(!starts_with('.')) %>% as.matrix()
  )

  n_iter2 <- nsamples(fit)

  post_means2 <- map(draws2, ~matrix(colMeans(.), nrow = 1))

  # Seperate out draws for residuals and their SD
  resid2 <- ranef(fit)$person[,1,1] # mean thetas
  stddev2 <- ranef(fit)$person[,2,1] # sd of thetas

  n_persons <- length(resid2)

  # Get standard quadrature points
  std_quad <- statmod::gauss.quad.prob(n_nodes, "normal", mu = 0, sigma = 1)
  std_log_weights <- log(std_quad$weights)

  linear_terms <- fitted(fit, scale = 'linear', summary = FALSE)
  linear_terms_mean <- matrix(colMeans(linear_terms), nrow = 1)

  itemparams <- fit %>% tidybayes::spread_draws(r_item__beta[item,], r_item__logalpha[item,], b_logalpha_Intercept) %>%
    mutate(r_item__alpha = exp(b_logalpha_Intercept + r_item__logalpha)) %>% as.data.frame()

  alphas <- matrix(itemparams$r_item__alpha, nrow = n_iter2)
  alphas <- rbind(alphas, colMeans(alphas))

  start = 1
  if (best_only) start = n_iter2

  # Extra iteration is to evaluate marginal log-likelihood at parameter means.
  ll <- foreach::`%dopar%`(foreach::foreach(i = start:(n_iter2 + 1), .combine = rbind,
                                            .packages = "matrixStats"
  ),
  {
    my_options <- options(digits.secs = 3)
    on.exit(options(my_options))

    if(i %% 100 == 0 ) {
      print(paste(i, "/", n_iter2, ":", strptime(Sys.time(), "%Y-%m-%d %H:%M:%OS") ))
    }

    ll_j <- matrix(NA, nrow = 1, ncol = n_persons)

    for(j in 1:n_persons) {

      # Set up adaptive quadrature using SD for residuals either from draws or
      # posterior mean (for best_ll).
      sd_i <- ifelse(i <= n_iter2, draws2$sd[[i]], post_means2$sd[[1]])
      adapt_nodes <- resid2[[j]] + stddev2[[j]] * std_quad$nodes
      log_weights <- log(sqrt(2*pi)) + log(stddev2[[j]]) + std_quad$nodes^2/2 +
        dnorm(adapt_nodes, sd = sd_i, log = TRUE) + std_log_weights

      # Evaluate mll with adaptive quadrature. If at n_iter + 1, evaluate
      # marginal likelihood at posterior means.
      if(i <= n_iter2) {
        loglik_by_node <- MFUN(adapt_nodes,  person = j, iter = i,
                               data_list = data_list2, draws = draws2, linear_terms = linear_terms,
                               alphas = alphas)

        weighted_loglik_by_node <- loglik_by_node + log_weights
        ll_j[1,j] <- matrixStats::logSumExp(weighted_loglik_by_node)
      } else {
        loglik_by_node <- MFUN(adapt_nodes,  person = j, iter = 1,
                               data_list = data_list2, draws = post_means2, linear_terms = linear_terms_mean,
                               alphas = alphas)
        weighted_loglik_by_node <- loglik_by_node + log_weights
        ll_j[1,j] <- matrixStats::logSumExp(weighted_loglik_by_node)
      }

    }

    ll_j

  })

  if(best_only) {
    return(ll[nrow(ll), ])
  } else {
    return(list(ll = ll[-nrow(ll), ], best_ll = ll[nrow(ll), ]))
  }

}

MFUN_2pl_secondapproach <- function(node, person, iter, data_list2, draws2, linear_terms, alphas) {
  #browser()

  resp_numbers <- data_list2$resp_number[data_list2$person_number == person]
  y <- data_list2$dich[resp_numbers]
  base_term <- linear_terms[iter, resp_numbers] - alphas[iter,]*draws2$theta[[iter, person]]

  p2 <- brms::inv_logit_scaled(matrix(rep(base_term, length(node)), nrow = length(node), byrow = TRUE) + node %*% t(alphas[iter,]))
  rowSums(dbinom(matrix(rep(y, length(node)), nrow = length(node), byrow = TRUE), 1, p2, log = TRUE))
}

ll_marg_brms_2pl2 <- mll_parallel_brms_2pl_secondapproach(fit_2pl, MFUN_2pl_secondapproach)

chain_brms_2pl2 <- fit_2pl %>% tidybayes::spread_draws(sd_person__theta_Intercept) %>% pull(.chain)
loo_ll_marg_brms_2pl2 <- loo(ll_marg_brms_2pl2$ll, r_eff = relative_eff(ll_marg_brms_2pl2$ll, chain_brms_2pl2))
plot(loo_ll_marg_brms3)
print(loo_ll_marg_brms_2pl2)
loo(fit_2pl)
