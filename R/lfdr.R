.get_pdf <- function(data,
                     design,
                     time_unit,
                     time_zero,
                     by,
                     features,
                     bandwidth) {
  condition_names <- colnames(design)[colnames(design) != time_unit]
  design <- unique(design[get(time_unit) != time_zero, ..by])
  design <- design[order(get(time_unit))]
  false_discoveries <- data[get(time_unit) == time_zero]
  pdf_list <- list()
  for(i in 1:nrow(design)) {
    mixture <- data[design[i,], , on = by]
    mixture_feature_table <- mixture[, ..features]
    
    false <- false_discoveries[unique(mixture[,..condition_names]), , on = condition_names]
    false_feature_table <- false[,..features]
    
    d0 <- lapply(false_feature_table, function(x) {
      density(rep(x, false$N), bw = bandwidth, n = 1e4)
    })
    g0 <- setNames(lapply(d0, approxfun, yleft = 0, yright = 0), features)
    
    nt <- sum(mixture$N)
    n0 <- sum(false$N)
    
    prior_prob <- n0/nt
    prior_odds <- n0/(nt - n0)
    
    g0 <- setNames(lapply(d0, function(den) {
      approxfun(den$x, den$y, yleft = 0, yright = 0)
    }), features)
    
    d <- lapply(mixture_feature_table, function(x) {
      density(rep(x, mixture$N), bw = bandwidth, n = 1e4)
    })
    g <- setNames(lapply(d, approxfun, yleft = 0, yright = 0), features)

    g1 <- setNames(mapply(function(a, b, c) {
      rg <- range(c)
      x <- seq(rg[1], rg[2], length.out = 1e4)
      d <- b(x) * prior_prob
      e <- a(x)
      d[d > e] <- e[d > e]
      f = (e-d)
      approxfun(x, f, yleft = 0, yright = 0)
    }, g, g0, mixture_feature_table, SIMPLIFY = F), features)
    
    pdf_list[paste0(design[i,], collapse = ",")] = list(list(pdf0 = g0, pdf1 = g1, pdf = g, prior_prob = prior_prob))
  }
  pdf_list
}


calc_LFDR <- function(data,
                      design,
                      time_unit,
                      time_zero,
                      by,
                      features,
                      bandwidth = "nrd0") {
  pdf_list = .get_pdf(data, design, time_unit, time_zero, by, features, bandwidth)
  design <- unique(design[get(time_unit) != time_zero, ..by])
  design <- design[order(get(Object@time_unit))]
  data$id <- seq_len(nrow(data))
  data$LFDR <- NA_real_
  for(i in 1:nrow(design)) {
    
    id <- data[design[i,], id, on = by]
    
    x <- as.list(data[design[i,], ..features, on = by])
    
    pdfs = pdf_list[paste0(design[i,], collapse = ",")]
    
    pdf0 <- function(x) {
      mapply(function(a,b) a(b), pdfs$pdf0, x)
    }
    
    pdf <- function(x) {
      mapply(function(a,b) a(b), pdfs$pdf, x)
    }
    
    BF <- pdf0(x)/pdf(x)
    n <- nrow(BF)
    LFDR <- numeric(n)
    prior_prob = pfds$prior_prob
    for (j in seq_len(n)) {
      LFDR[j] <- prior_prob * prod(BF[j,], na.rm = T)
    }
    data[id,"LFDR"] <- LFDR
  }
  data[,-"id"]
}

filter_data <- function(Object,
                        LFDR_threshold = c("strong", "substantial"),
                        min_nb_timepoints = 0,
                        min_labeling_ratio = 0,
                        min_nb_sample = NULL,
                        score_threshold = NULL,
                        higher_score_better = T,
                        trace = T) {
  data <- copy(Object@data)
  if(Object@peptide_centric) {
    by = c(Object@peptide_column_PTMs, Object@peptide_column_no_PTMs, Object@protein_column)
  } else {
    by = Object@protein_column
  }
  data = data[get(Object@time_unit) != Object@time_zero]
  if(LFDR_threshold[1] == "substantial") {
    LFDR_threshold <- 1-1/(1+1/3)
  } else if (LFDR_threshold[1] == "strong") {
    LFDR_threshold <- 1-1/(1+1/10)
  }
  score_columns <- .check_columns(Object@score_name, colnames(data))
  if(!is.null(score_threshold)) {
    if(higher_score_better) {
      data <- data[(LFDR <= LFDR_threshold) &  (get(Object@score_name) >= score_threshold) & (get(score_columns[1]) >= score_threshold) & (LR >= min_labeling_ratio)]
    } else {
      data <- data[(LFDR <= LFDR_threshold) &  (get(Object@score_name) <= score_threshold) & (get(score_columns[1]) <= score_threshold) & (LR >= min_labeling_ratio)]
    }
  } else {
    data <- data[(LFDR <= LFDR_threshold) & (get(Object@labeling_ratio_name) >= min_labeling_ratio)]
  }
  if(!is.null(min_nb_sample)) {
    condition_names <- colnames(Object@design)[colnames(Object@design) != Object@time_unit]
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
  incorporation_columns <- .get_columns(Object@incorporation_name, colnames(data))
  intensity_columns <- .get_columns(Object@intensity_name, colnames(data))
  score_columns <- .get_columns(Object@score_name, colnames(data))
  col <- c(incorporation_columns[1], Object@incorporation_name, Object@intensity_name, intensity_columns[1], Object@score_name, score_columns[1], "LR", "LFDR")
  max_LR <- max(data[,get(Object@labeling_ratio_name)])
  if(max_LR > 1 && max_LR <= 100) {
    LR <- data[,get(Object@labeling_ratio_name)]/100
  } else if (max_LR > 0 && max_LR <= 1) {
    LR <- data[,get(Object@labeling_ratio_name)]
  }
  
  data[,Object@intensity_name] <- LR * data[[intensity_columns[1]]]/(1 - LR)
  
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
  # plot(density(shrink_stat$RIA[X$day == 8]))
  dt <- data[, lapply(.SD, median, na.rm = T), by = c(by ,Object@time_unit), .SDcols = col]
  by2 = by
  by2[!(by2 %in% make.names(by2))] <- paste0('`', by2[!(by2 %in% make.names(by2))], '`')
  formula = eval(parse(text = paste(paste0(by2, collapse = "+"), " ~ ", Object@time_unit)))
  dt <- dcast(dt, formula = formula, value.var = col)
  timepoints <- sort(unique(data$day))
  colnames(dt)[-(seq_along(by))] <- get_cols(c(paste("Light", Object@incorporation_name),
                                               paste("Heavy", Object@incorporation_name),
                                               paste("Light", Object@intensity_name),
                                               paste("Heavy", Object@intensity_name),
                                               paste("Light", Object@score_name),
                                               paste("Heavy", Object@score_name),
                                               "LR", "LFDR"), timepoints)
  Object@master_tbl <- dt
  Object <- getFunc(Object)
  Object <- getTaxon(Object)
  Object
}





