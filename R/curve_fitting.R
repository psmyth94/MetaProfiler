curve_fitting <- function(Object,
                          var,
                          use_knn = T,
                          lower = c(0,0,0,0),
                          upper = c(100,100,100,100),
                          control = DEoptim.control(trace = F),
                          equation = 3,
                          progress = T,
                          seed = 362436069) {
  # if(is.null(fn)) {
  #   fn = cppFunction(
  #     'double three_exponential_function(const Rcpp::NumericVector & x, const Rcpp::NumericVector & xr, const Rcpp::NumericVector & yr) {
  #       double u = ((x[3] + x[1] + x[3]) - std::sqrt(std::pow((x[3] + x[1] + x[3]), 2) - (4 * x[1]*x[3]))) / 2;
  #       double v = ((x[3] + x[1] + x[3]) + std::sqrt(std::pow((x[3] + x[1] + x[3]), 2) - (4 * x[1]*x[3]))) / 2;
  #       double yu = x[1] * x[2]*(u - x[3]) / ((u - v)*(u - x[2])*u);
  #       double yv = x[1] * x[2]*(v - x[3]) / ((v - u)*(v - x[2])*v);
  #       double ykbi = x[1] * (x[2] - x[3]) / ((u - x[2])*(v - x[2]));
  #       double sum = 0;
  #       for(int i; i < xr.size(); ++i) {
  #         double value = (1 + yu * exp(-u * xr[i]) + yv * exp(-v * xr[i]) + ykbi * exp(-x[2] * xr[i]))*100;
  #         sum += std::pow(value - yr[i], 2);
  #       }
  #       return sum;
  #     }'
  #   )
  # }
  ##fn1  <- function(par) fn(par, ...)
  if (length(lower) != length(upper))
    stop("'lower' and 'upper' are not of same length")
  if (!is.vector(lower))
    lower <- as.vector(lower)
  if (!is.vector(upper))
    upper <- as.vector(upper)
  if (any(lower > upper))
    stop("'lower' > 'upper'")
  if (any(lower == "Inf"))
    warning("you set a component of 'lower' to 'Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (any(lower == "-Inf"))
    warning("you set a component of 'lower' to '-Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (any(upper == "Inf"))
    warning("you set a component of 'upper' to 'Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (any(upper == "-Inf"))
    warning("you set a component of 'upper' to '-Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (!is.null(names(lower)))
    consts <- names(lower)
  else if (!is.null(names(upper)) & is.null(names(lower)))
    consts <- names(upper)
  else
    consts <- paste("par", 1:length(lower), sep = "")
  consts <- c(consts, "score")
  if(Object@peptide_centric) {
    group <- Object@peptide_column_PTMs
  } else {
    group <- Object@protein_column
  }
  col <- get_cols(var)
  mat <- Object@master_tbl[,c(group, col),with=F]
  if(use_knn)
  {
    shush(mat[,colnames(mat)[-seq_along(group)]] <- as.data.table(
      impute.knn(as.matrix(mat[,-..group]), rowmax = 1, colmax = 1)$data
    ))
  }
  if(is.null(equation)) {
    mat <- melt(mat, id.vars = group)
    p <- dplyr::progress_estimated(length(unique(mat[[group]])))
    f <- function(x, y) {
      k <- cobs::cobs(x, y,
                      pointwise = matrix(c(0,0,0), 1, 3),
                      print.warn = F, print.mesg = F, trace = F)
      p$tick()$print()
      list(k)
    }
    consts <- mat[, .(k = f(Object@timepoints, value)), by = group]
    consts$equation <- 4
    Object@consts[var] <- list(consts)
    return(Object)
  }
  
  
  mat <- cbind(mat[,..group], data.table(
    xr = list(Object@timepoints),
    yr = unlist(apply(as.matrix(mat[,-..group]), 1, function(x) list(x)), recursive = F)
  ))
  # makeEnv <- function(...) {
  #   list2env(list(...))
  # }
  # env <- mapply(function(xr, yr){
  #   makeEnv(xr = xr, yr = yr)
  # }, mat$xr, mat$yr)
  ctrl <- do.call(DEoptim.control, as.list(control))
  ctrl$npar <- length(lower)
  if (ctrl$NP < 4) {
    warning("'NP' < 4; set to default value 50\n", immediate. = TRUE)
    ctrl$NP <- 50
  }
  if (ctrl$NP < 10*length(lower))
    warning("For many problems it is best to set 'NP' (in 'control') to be at least ten times the length of the parameter vector. \n", immediate. = TRUE)
  if (!is.null(ctrl$initialpop)) {
    ctrl$specinitialpop <- TRUE
    if(!identical(as.numeric(dim(ctrl$initialpop)), c(ctrl$NP, ctrl$npar)))
      stop("Initial population is not a matrix with dim. NP x length(upper).")
  } else {
    ctrl$specinitialpop <- FALSE
    ctrl$initialpop <- matrix(0,1,1)    # dummy matrix
  }
  ##
  ctrl$trace <- as.numeric(ctrl$trace)
  ctrl$specinitialpop <- as.numeric(ctrl$specinitialpop)
  set.seed(seed)
  mat[, consts] <-as.data.table(t(curve_fitting_c(mat$xr,
                                                  mat$yr,
                                                  minbound = lower,
                                                  maxbound = upper,
                                                  equation = equation,
                                                  control = ctrl,
                                                  verbose = progress)))
  mat$equation <- equation
  Object@consts[var] <- list(mat[,c(group, consts, "equation"), with = F])
  Object
}

curve_fitting_test <- function(Object, var,
                               use_knn = T,
                               lower = c(0,0,0,0),
                               upper = c(100,100,100,100),
                               control = DEoptim.control(trace = F),
                               progress = T,
                               equation = 3,
                               seed = 362436069) {
  # if(is.null(fn)) {
  #   fn = cppFunction(
  #     'double three_exponential_function(const Rcpp::NumericVector & x, const Rcpp::NumericVector & xr, const Rcpp::NumericVector & yr) {
  #       double u = ((x[3] + x[1] + x[3]) - std::sqrt(std::pow((x[3] + x[1] + x[3]), 2) - (4 * x[1]*x[3]))) / 2;
  #       double v = ((x[3] + x[1] + x[3]) + std::sqrt(std::pow((x[3] + x[1] + x[3]), 2) - (4 * x[1]*x[3]))) / 2;
  #       double yu = x[1] * x[2]*(u - x[3]) / ((u - v)*(u - x[2])*u);
  #       double yv = x[1] * x[2]*(v - x[3]) / ((v - u)*(v - x[2])*v);
  #       double ykbi = x[1] * (x[2] - x[3]) / ((u - x[2])*(v - x[2]));
  #       double sum = 0;
  #       for(int i; i < xr.size(); ++i) {
  #         double value = (1 + yu * exp(-u * xr[i]) + yv * exp(-v * xr[i]) + ykbi * exp(-x[2] * xr[i]))*100;
  #         sum += std::pow(value - yr[i], 2);
  #       }
  #       return sum;
  #     }'
  #   )
  # }
  ##fn1  <- function(par) fn(par, ...)
  if (length(lower) != length(upper))
    stop("'lower' and 'upper' are not of same length")
  if (!is.vector(lower))
    lower <- as.vector(lower)
  if (!is.vector(upper))
    upper <- as.vector(upper)
  if (any(lower > upper))
    stop("'lower' > 'upper'")
  if (any(lower == "Inf"))
    warning("you set a component of 'lower' to 'Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (any(lower == "-Inf"))
    warning("you set a component of 'lower' to '-Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (any(upper == "Inf"))
    warning("you set a component of 'upper' to 'Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (any(upper == "-Inf"))
    warning("you set a component of 'upper' to '-Inf'. May imply 'NaN' results", immediate. = TRUE)
  group <- colnames(Object@master_tbl)[1]
  col <- get_cols(var)
  mat <- Object@master_tbl[,c(group, col),with=F]
  if(use_knn) {
    shush(mat[,colnames(mat)[-seq_along(group)]] <- as.data.table(
      impute.knn(as.matrix(mat[,-..group]), rowmax = 1, colmax = 1)$data
    ))
  }
  mat <- cbind(mat[,..group], data.table(
    xr = list(Object@timepoints),
    yr = unlist(apply(as.matrix(mat[,-..group]), 1, function(x) list(x)), recursive = F)
  ))
  # makeEnv <- function(...) {
  #   list2env(list(...))
  # }
  # env <- mapply(function(xr, yr){
  #   makeEnv(xr = xr, yr = yr)
  # }, mat$xr, mat$yr)
  ctrl <- do.call(DEoptim.control, as.list(control))
  ctrl$npar <- length(lower)
  if (ctrl$NP < 4) {
    warning("'NP' < 4; set to default value 50\n", immediate. = TRUE)
    ctrl$NP <- 50
  }
  if (ctrl$NP < 10*length(lower))
    warning("For many problems it is best to set 'NP' (in 'control') to be at least ten times the length of the parameter vector. \n", immediate. = TRUE)
  if (!is.null(ctrl$initialpop)) {
    ctrl$specinitialpop <- TRUE
    if(!identical(as.numeric(dim(ctrl$initialpop)), c(ctrl$NP, ctrl$npar)))
      stop("Initial population is not a matrix with dim. NP x length(upper).")
  } else {
    ctrl$specinitialpop <- FALSE
    ctrl$initialpop <- matrix(0,1,1)    # dummy matrix
  }
  ##
  ctrl$trace <- as.numeric(ctrl$trace)
  ctrl$specinitialpop <- as.numeric(ctrl$specinitialpop)
  set.seed(seed)
  ctrl$timepoints <- Object@timepoints
  SE <-curve_fitting_test_c(mat$xr,
                            mat$yr,
                            minbound = lower,
                            maxbound = upper,
                            equation = equation,
                            control = ctrl,
                            verbose = progress)
  SE[t(is.na(Object@master_tbl[,..col]))] <- NA
  RMSE <- sqrt(colMeans(SE, na.rm = T))
}


impute <- function(Object, vars, timepoints = Object@timepoints)
{
  for(var in vars) {
    value = Object@model(var = var, x = timepoints)
    ori = Object@master_tbl[, get_cols(var, timepoints), with = F]
    ori[is.na(ori)] <- value[is.na(ori)]
    Object@master_tbl[,get_cols(var, timepoints)] <- ori
  }
  Object
}


impute2 <- function(Object, vars, timepoints = Object@timepoints)
{
  for(var in impute) {
    value = Object@model(var = var, x = timepoints)
    Object@master_tbl[,get_cols(var, timepoints)] <- as.data.table(value)
  }
  Object
}


