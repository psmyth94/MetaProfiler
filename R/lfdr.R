start_arg_default <- function (x, distr)
{
  if (distr == "norm") {
    n <- length(x)
    sd0 <- sqrt((n - 1)/n) * sd(x)
    mx <- mean(x)
    start <- list(mean = mx, sd = sd0)
  }
  else if (distr == "lnorm") {
    if (any(x <= 0)) 
      stop("values must be positive to fit a lognormal distribution")
    n <- length(x)
    lx <- log(x)
    sd0 <- sqrt((n - 1)/n) * sd(lx)
    ml <- mean(lx)
    start <- list(meanlog = ml, sdlog = sd0)
  }
  else if (distr == "exp") {
    if (any(x < 0)) 
      stop("values must be positive to fit an exponential  distribution")
    start <- list(rate = 1/mean(x))
  }
  else if (distr == "gamma") {
    if (any(x < 0)) 
      stop("values must be positive to fit an gamma  distribution")
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n * var(x)
    start <- list(shape = m^2/v, rate = m/v)
  }
  else if (distr == "beta") {
    if (any(x < 0) | any(x > 1)) 
      stop("values must be in [0-1] to fit a beta distribution")
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n * var(x)
    aux <- m * (1 - m)/v - 1
    start <- list(shape1 = m * aux, shape2 = (1 - m) * aux)
  }
  else if (distr == "weibull") {
    if (any(x < 0)) 
      stop("values must be positive to fit an Weibull  distribution")
    m <- mean(log(x))
    v <- var(log(x))
    shape <- 1.2/sqrt(v)
    scale <- exp(m + 0.572/shape)
    start <- list(shape = shape, scale = scale)
  }
  else if (distr == "logis") {
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n * var(x)
    start <- list(location = m, scale = sqrt(3 * v)/pi)
  }
  else if (distr == "cauchy") {
    start <- list(location = median(x), scale = IQR(x)/2)
  }
  else if (distr == "unif") {
    start <- list(min = 0, max = 1)
  }
  else if (distr == "invgamma") {
    if (any(x < 0)) 
      stop("values must be positive to fit an inverse gamma  distribution")
    m1 <- mean(x)
    m2 <- mean(x^2)
    shape <- (2 * m2 - m1^2)/(m2 - m1^2)
    scale <- m1 * m2/(m2 - m1^2)
    start <- list(shape = shape, scale = scale)
  }
  else if (distr == "llogis") {
    if (any(x < 0)) 
      stop("values must be positive to fit a log-logistic  distribution")
    p25 <- as.numeric(quantile(x, 0.25))
    p75 <- as.numeric(quantile(x, 0.75))
    shape <- 2 * log(3)/(log(p75) - log(p25))
    scale <- exp(log(p75) + log(p25))/2
    start <- list(shape = shape, scale = scale)
  }
  else if (distr == "invweibull") {
    if (any(x < 0)) 
      stop("values must be positive to fit an inverse Weibull distribution")
    g <- log(log(4))/(log(log(4/3)))
    p25 <- as.numeric(quantile(x, 0.25))
    p75 <- as.numeric(quantile(x, 0.75))
    shape <- exp((g * log(p75) - log(p25))/(g - 1))
    scale <- log(log(4))/(log(shape) - log(p25))
    start <- list(shape = shape, scale = max(scale, 1e-09))
  }
  else if (distr == "pareto1") {
    if (any(x < 0)) 
      stop("values must be positive to fit a Pareto distribution")
    x1 <- min(x)
    m1 <- mean(x)
    n <- length(x)
    shape <- (n * m1 - x1)/(n * (m1 - x1))
    min <- x1 * (n * shape - 1)/(n * shape)
    start <- list(shape = shape, min = min)
  }
  else if (distr == "pareto") {
    if (any(x < 0)) 
      stop("values must be positive to fit a Pareto distribution")
    m1 <- mean(x)
    m2 <- mean(x^2)
    scale <- (m1 * m2)/(m2 - 2 * m1^2)
    shape <- 2 * (m2 - m1^2)/(m2 - 2 * m1^2)
    start <- list(shape = shape, scale = scale)
  }
  else if (distr == "lgamma") {
    if (any(x < 0)) 
      stop("values must be positive to fit a log-gamma distribution")
    m1 <- mean(log(x))
    m2 <- mean(log(x)^2)
    alpha <- m1^2/(m2 - m1^2)
    lambda <- m1/(m2 - m1^2)
    start <- list(shapelog = alpha, ratelog = lambda)
  }
  else if (distr == "trgamma") {
    if (any(x < 0)) 
      stop("values must be positive to fit an trans-gamma  distribution")
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n * var(x)
    start <- list(shape1 = m^2/v, shape2 = 1, rate = m/v)
  }
  else if (distr == "invtrgamma") {
    if (any(x < 0)) 
      stop("values must be positive to fit an inverse trans-gamma  distribution")
    n <- length(1/x)
    m <- mean(1/x)
    v <- (n - 1)/n * var(1/x)
    start <- list(shape1 = m^2/v, shape2 = 1, rate = m/v)
  }
  else stop(paste0("Unknown starting values for distribution ", distr, "."))
  return(start)
}

.fit_true_distribution <- function(data, kde, kde0, distr, times, start_prior_prob) {
  start_arg = mapply(function(dname, dat) {
    start_arg = start_arg_default(dat, dname)
  }, distr, data, SIMPLIFY = F)
  x = seq(0, 100, length.out = 512)
  y = seq(0, 100, length.out = 512)
  tbl = setNames(expand.grid(x,y), c('x', 'y'))
  prob = ks::dkde(tbl, kde)
  prob[prob < 0] = 0
  i = sample(seq_len(nrow(tbl)), 10000, replace = T, prob = prob)
  data1 = tbl[i,]
  u1 = VineCopula::pobs(data1)
  cop = copula::normalCopula(dim = 2)
  fit <- copula::fitCopula(cop, u1)
  estimates = fit@estimate
  vstart = c(estimates, unlist(start_arg), start_prior_prob)
  nestim = length(estimates)
  split_by = rep(1:2, lengths(start_arg))
  arg_names = unlist(lapply(start_arg, names))
  obs = as.matrix(data1)
  nparams = length(vstart) - 1
  p0 = ks::dkde(obs, kde0)
  mle = function(x) {
    args = setNames(x[(nestim+1):nparams], arg_names)
    args = split(setNames(as.list(args), arg_names), split_by)
    cop@parameters = x[1:nestim]
    mv_cop = try(suppressMessages(copula::mvdc(cop, distr, args)), silent = T)
    if(inherits(mv_cop, "try-error")) return(.Machine$double.xmax)
    prior_prob = x[nparams + 1]
    if (prior_prob < 0 || prior_prob > 1) return(.Machine$double.xmax)
    p1 = copula::dMvdc(obs, mv_cop)
    p = p0 * prior_prob + p1 * (1 - prior_prob)
    val = -sum(log(p))
    val
  }
  opt <- optim(par = vstart, f = mle, method = "Nelder-Mead")
  x = opt$par
  prior_prob = x[nparams + 1]
  args = setNames(x[(nestim+1):nparams], arg_names)
  args = split(setNames(as.list(args), arg_names), split_by)
  cop@parameters = x[1:nestim]
  mv_cop = try(suppressMessages(copula::mvdc(cop, distr, args)), silent = T)
  list(copula = cop, mvdc = mv_cop, prior_prob = prior_prob)
}

lfdr <- function(Object, by = Object@time_unit, estimate_prior_probability = c("mle", "analytically"),
                 distr = rep("weibull", 2), trace = T, ...) {
  observations = c(Object@incorporation_name, Object@labeling_ratio_name)
  estimate_prior_probability = match.arg(estimate_prior_probability, c("mle", "analytically"))
  dname = distr
  ddname = paste0("d", dname)
  viable_ddname = exists(ddname, mode = "function")
  if(!viable_ddname) {
    n = sum(!viable_ddname)
    msg = combine_words(ddname[!viable_ddname])
    stop("The function ", ddname, " is not defined.")
  }
  condition_names <- colnames(Object@design)[colnames(Object@design) != Object@time_unit]
  if (length(condition_names))
    Object@data[,condition_names] <- Object@data[,lapply(.SD, as.character), .SDcols = condition_names]
  design <- unique(Object@design[get(Object@time_unit) != Object@time_zero, ..by])
  design <- design[order(get(Object@time_unit))]
  false_discoveries <- Object@data[get(Object@time_unit) == Object@time_zero]
  false <- false_discoveries
  false_observation_table <- false[,..observations]
  x0 = as.matrix(false_observation_table[, lapply(.SD, function(x) rep(x, false$N))])
  H0 = ks::Hpi(x0)
  kde0 = ks::kde(x0)
  Object@data$id <- seq_len(nrow(Object@data))
  Object@data$LFDR <- NA_real_
  
  for(i in 1:nrow(design)) {
    mixture <- Object@data[design[i,], , on = by]
    mixture_observation_table <- mixture[, ..observations]
    if(length(condition_names)){
      tmp = false_discoveries[unique(mixture[,..condition_names]), , on = condition_names]
      if(!is.logical(all.equal(false, tmp, ignore.row.order = T)))
      {
        false = tmp
        false_observation_table = false[,..observations]
        x0 = as.matrix(false_observation_table[, lapply(.SD, function(x) rep(x, false$N))])
        H0 = ks::Hpi(x0)
        kde0 = ks::kde(x0)
      }
    }
    prior_prob = sum(false$N)/sum(mixture$N)
    x = as.matrix(mixture_observation_table[, lapply(.SD, function(x) rep(x, mixture$N))])
    H = ks::Hpi(x)
    kde = ks::kde(x, H)
    
    if(estimate_prior_probability == "mle") {
      if(!length(distr)) {
        stop("mle was chosen for estimating prior probability. However, the parametric distribution for true discoveries is not specified in `distr`.")
      }
      if (trace) {
        cat(crayon::blue("Fitting bivariate Gaussian copula to", combine_words(observations), "for true discoveries at", 
                         tolower(Object@time_unit), unique(mixture[,get(Object@time_unit)])))
        cat(crayon::blue("..."))
      }
      fit <- .fit_true_distribution(mixture_observation_table, kde, kde0, distr, mixture$N, prior_prob)
      if (trace) {
        cat(crayon::blue(" done.\n"))
      }
      x = as.matrix(mixture_observation_table)
      prior_prob = fit$prior_prob
      LTDR = (1 - prior_prob) * (copula::dMvdc(x, fit$mvdc)/ks::dkde(x, kde))
      LTDR[LTDR > 1] = 1
      LFDR = 1 - LTDR
    } else {
      fit = NULL
      LFDR = prior_prob * ks::dkde(x, kde0)/ks::dkde(x, kde)
    }
    Object@data$LFDR[mixture$id] = LFDR
    Object@pdf[paste(by, design[i,], collapse = ", ")] = list(c(fit, list(data = mixture_observation_table, kde = kde, kde0 = kde0)))
  }
  Object
}

filter_data <- function(Object, peptide_centric = TRUE, LFDR_threshold = c("strong", "substantial"),
                        min_nb_timepoints = 0, min_labeling_ratio = 0, min_nb_sample = NULL,
                        score_threshold = NULL, higher_score_better = T, trace = T) {
  data <- data.table::copy(Object@data)
  if(peptide_centric) {
    by = c(Object@peptide_column_PTMs, Object@peptide_column_no_PTMs, Object@accession_column)
  } else {
    by = Object@accession_column
  }
  data = data[get(Object@time_unit) != Object@time_zero]
  if(LFDR_threshold[1] == "substantial") {
    LFDR_threshold <- 1-1/(1+1/3)
  } else if (LFDR_threshold[1] == "strong") {
    LFDR_threshold <- 1-1/(1+1/10)
  }
  if(!is.null(score_threshold)) {
    if(higher_score_better) {
      data <- data[(LFDR <= LFDR_threshold) &  (get(Object@score_name) >= score_threshold) & (get(Object@score_columns[1]) >= score_threshold) & (LR >= min_labeling_ratio)]
    } else {
      data <- data[(LFDR <= LFDR_threshold) &  (get(Object@score_name) <= score_threshold) & (get(Object@score_columns[1]) <= score_threshold) & (LR >= min_labeling_ratio)]
    }
  } else {
    data <- data[(LFDR <= LFDR_threshold) & (get(Object@labeling_ratio_name) >= min_labeling_ratio)]
  }
  condition_names <- colnames(Object@design)[colnames(Object@design) != Object@time_unit]
  if(!is.null(min_nb_sample) && length(condition_names)) {
    data <- data[,test := (length(unique(get(Object@time_unit))) >= min_nb_timepoints) & (nrow(unique(.SD)) >= min_nb_sample), by = by, .SDcols = condition_names]
  } else {
    data <- data[,test := (length(unique(get(Object@time_unit))) >= min_nb_timepoints), by = by]
  }
  data <- data[(test),-"test"]
  if(trace) {
    message(paste(length(unique(data[[Object@peptide_column_PTMs]])), "peptides left"))
  }
  # unequal_var_shrink <- function (stat, vars) 
  # {
  #   n = length(stat)
  #   if (!n == length(vars)) 
  #     stop("stat and vars must be the same length.")
  #   vars[is.na(vars)] = 0
  #   mle_result = optimize(function(a_hat) {
  #     likelihood = 1
  #     for (i in 1:length(stat)) {
  #       if(vars[i] != 0) {
  #         likelihood = likelihood - (((stat[i]^2)/vars[i]) * (vars[i]/(vars[i] + a_hat))) + log(vars[i]/(vars[i] + a_hat))/2
  #       }
  #     }
  #     likelihood
  #   }, interval = c(0, 1e+06), maximum = TRUE)
  #   B_hat_mle = numeric(length = n)
  #   for (i in 1:n) {
  #     B_hat_mle[i] = vars[i]/(vars[i] + mle_result$maximum)
  #   }
  #   mu_hat_mle = numeric(length = n)
  #   for (i in 1:n) {
  #     mu_hat_mle[i] = (1 - B_hat_mle[i]) * stat[i]
  #   }
  #   res_df = data.frame(stat = stat, variance = vars, b_hat = B_hat_mle, 
  #                       shrink_stat = mu_hat_mle)
  #   return(res_df)
  # }
  col <- c(Object@incorporation_columns[1], Object@incorporation_name, Object@intensity_name, Object@intensity_columns[1], Object@score_name, Object@score_columns[1], Object@labeling_ratio_name, "LFDR")
  LR <- data[,get(Object@labeling_ratio_name)]
  if(.check_if_percantage_or_fraction(LR)) {
    LR <- LR/100
  }
  
  data[,Object@intensity_name] <- LR * data[[Object@intensity_columns[1]]]/(1 - LR)
  
  # X <- copy(data)
  # X[,col] <- X[,..col]/100
  # X[[Object@time_unit]] <- as.character(X[[Object@time_unit]])
  # X_mean <- X[, lapply(.SD, mean, na.rm = T), by = group_by, .SDcols = col]
  # X[,paste0(col, "_scaled")] <- NA_real_
  # X <- X[, c(paste0(col, "_mean")) := lapply(.SD, mean, na.rm = T), by = c(Object@time_unit), .SDcols = col]
  # X <- X[, c(paste0(col, "_sd")) := lapply(.SD, sd, na.rm = T), by = c(Object@time_unit), .SDcols = col]
  # X <- X[, c(paste0(col, "_scaled")) := lapply(.SD, function(x) {
  #   (x - mean(x, na.rm = T))/sd(x, na.rm = T)
  # }), by = c(Object@time_unit), .SDcols = col]
  # X <- X[,c(lapply(.SD, var, na.rm = T), lapply(.SD, mean, na.rm = T)), by = group_by, .SDcols = c(col, paste0(col, "_scaled"), paste0(col, "_sd"), paste0(col, "_mean"))]
  # colnames(X)[seq_along(c(col, paste0(col, "_scaled"), paste0(col, "_sd"), paste0(col, "_mean"))) + length(group_by)] <- paste0(c(col, paste0(col, "_scaled"), paste0(col, "_sd"), paste0(col, "_mean")), "_var")
  # shrink_stat <- X[, .(Peptides = Peptides,
  #                      RIA = (unequal_var_shrink(RIA_scaled, RIA_scaled_var)$shrink_stat * RIA_sd + RIA_mean)*100,
  #                      LR = (unequal_var_shrink(LR_scaled, LR_scaled_var)$shrink_stat * LR_sd + LR_mean)*100), by = c(Object@time_unit)]
  dt <- data[, lapply(.SD, median, na.rm = T), by = c(by, Object@time_unit), .SDcols = col]
  
  by2 = by
  by2[!(by2 %in% make.names(by2))] <- paste0('`', by2[!(by2 %in% make.names(by2))], '`')
  formula = eval(parse(text = paste(paste0(by2, collapse = "+"), " ~ ", Object@time_unit)))
  dt <- data.table::dcast(dt, formula = formula, value.var = col)
  timepoints <- c(Object@time_zero, Object@timepoints)
  col <- c(by, paste0(rep(col, each = length(timepoints)), "_", rep(timepoints, length(col))))
  col_add = setdiff(col, colnames(dt))
  dt[,col_add] <- NA
  dt <- dt[,..col]
  colnames(dt)[-(seq_along(by))] <- Object@get_cols(c(paste("Light", Object@incorporation_name),
                                                      paste("Heavy", Object@incorporation_name),
                                                      paste("Light", Object@intensity_name),
                                                      paste("Heavy", Object@intensity_name),
                                                      paste("Light", Object@score_name),
                                                      paste("Heavy", Object@score_name),
                                                      Object@labeling_ratio_name, "LFDR"), timepoints)
  Object@master_tbl <- dt
  if(length(Object@function_columns))
    Object <- getFunc(Object)
  if(length(Object@taxon_columns))
    Object <- getTaxon(Object)
  Object
}






