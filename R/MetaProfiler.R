setClass("MetaProfiler",
         representation(design = "data.table",
                        data = "data.table",
                        incorporation_name = "character",
                        intensity_name = "character",
                        labeling_ratio_name = "character",
                        score_name = "character",
                        peptide_column_no_PTMs = "character",
                        peptide_column_PTMs = "character",
                        protein_column = "character",
                        peptide_centric = "logical",
                        annotate_by_peptide = "character",
                        master_tbl = "data.table",
                        timepoints="numeric",
                        time_zero = "numeric",
                        time_unit = "character",
                        annotation_table = "data.table",
                        taxon_columns = "character",
                        function_columns = "character",
                        deconvolve = "function",
                        consts = "list",
                        model = "function",
                        pdf = "list"))

MetaProFiler <- function(design,
                         time_unit,
                         time_zero = 0,  
                         peptide_centric = T,
                         
                         incorporation_name = "RIA",
                         intensity_name = "INT",
                         labeling_ratio_name = NULL,
                         score_name = "Cor.",
                         incorporation_columns = "auto",
                         intensity_columns = "auto",
                         labeling_ratio_columns = "auto",
                         score_columns = "auto",
                         as_percentage = T,
                         
                         ...,
                         
                         data = NULL,
                         data_peptide_column_PTMs = "guess",
                         data_peptide_column_no_PTMs = "guess",
                         data_accession_column = "guess",
                         
                         pep2pro = NULL,
                         pep2pro_peptide_column = "guess",
                         pep2pro_accession_column = "guess",
                         pep2taxon = NULL,
                         pep2taxon_peptide_column = "guess",
                         rank_columns = "guess",
                         
                         compute_razor_protein = F,
                         
                         pro2func = NULL,
                         pro2func_accession_column = "guess",
                         pro2func_function_columns = "guess",
                         
                         annotate_by_peptide = c("unmodified", "um", "modified", "m"),
                         
                         keep_feature_type = c("feature", "id"),
                         cluster = T, 
                         cluster_by = c("RIA", "LR"),
                         radius = rep(0.1, length(cluster_by)),
                         distance_method = rep("relative", length(cluster_by)),
                         element_wise = T,
                         group_by = peptide_centric,
                         
                         LFDR_features = NULL,
                         LFDR_by = NULL,
                         bandwidth = "nrd0",
                         
                         trace = T,
                         progress = T) {
  
  if(is.null(data) || is.character(data)) {
    data <- load_results(design,
                         data,
                         time_unit,
                         incorporation_name,
                         intensity_name,
                         labeling_ratio_name,
                         score_name,
                         incorporation_columns,
                         intensity_columns,
                         labeling_ratio_columns,
                         score_columns,
                         as_percentage,
                         progress,
                         ...)
    check = unlist(design[, lapply(.SD, function(x) {
      test <- try(all(file.exists(x)), silent = T)
      ifelse(inherits(test, "try-error"), F, test)
    })])
    design <- design[,!check, with = F]
  }
  cols_intersect <- colnames(data)[colnames(data) %in% colnames(design)]
  cols_diff <- setdiff( colnames(design), colnames(data))
  if(trace && (length(cols_diff) > 0)) {
    message(paste0(ifelse(length(cols_diff) > 1,
                          paste0("Columns c(", paste0("\"", cols_diff, "\"", collapse = ", "),") in `design` are"),
                          paste0("Column \"", cols_diff, "\" in `design` is")),
                   " not found in `data`. Removing from table `design`."))
  }
  design <- design[,..cols_intersect]
  reorder_cols = c(time_unit, setdiff(colnames(design), time_unit))
  design <- design[.mixedorder(design[,..reorder_cols]),]
  
  if(is.character(pep2pro)) {
    stopifnot(file.exists(pep2pro))
    pep2pro <- fread(pep2pro)
  }
  
  if(is.character(pep2taxon)) {
    stopifnot(file.exists(pep2taxon))
    pep2taxon <- fread(pep2taxon)
  }
  
  if(is.character(pro2func)) {
    stopifnot(file.exists(pro2func))
    pro2func <- fread(pro2func)
  }
  
  args = .guess_columns(data,
                        data_peptide_column_PTMs,
                        data_peptide_column_no_PTMs,
                        data_accession_column,
                        pep2pro,
                        pep2pro_peptide_column,
                        pep2pro_accession_column,
                        pep2taxon,
                        pep2taxon_peptide_column,
                        rank_columns,
                        pro2func,
                        pro2func_accession_column,
                        pro2func_function_columns,
                        annotate_by_peptide,
                        trace)
  data = args$data
  data_peptide_column_PTMs = args$data_peptide_column_PTMs
  data_peptide_column_no_PTMs = args$data_peptide_column_no_PTMs
  data_accession_column = args$data_accession_column
  pep2pro = args$pep2pro
  pep2pro_peptide_column = args$pep2pro_peptide_column
  pep2pro_accession_column = args$pep2pro_accession_column
  pep2taxon = args$pep2taxon
  pep2taxon_peptide_column = args$pep2taxon_peptide_column
  rank_columns = args$rank_columns
  pep2taxon_columns = args$pep2taxon_columns
  pro2func = args$pro2func
  pro2func_accession_column = args$pro2func_accession_column
  pro2func_function_columns = args$pro2func_function_columns
  
  if(!is.numeric(design[[time_unit]])) {
    design[[time_unit]] <- as.numeric(design[[time_unit]])
    if(all(is.na(design[[time_unit]]))) {
      stop(paste0("Could not convert column `", time_unit, "` to numeric. Are there any non-numeric characters in the column?"))
    }
  }
  
  annotation_table <- make_annotation_table(time_unit,
                                            data,
                                            data_peptide_column_PTMs,
                                            data_peptide_column_no_PTMs,
                                            data_accession_column,
                                            annotate_by_peptide,
                                            pep2pro,
                                            pep2pro_peptide_column,
                                            pep2pro_accession_column,
                                            compute_razor_protein,
                                            accession_delimiter,
                                            pep2taxon,
                                            pep2taxon_peptide_column,
                                            rank_columns,
                                            pro2func,
                                            pro2func_accession_column,
                                            pro2func_function_columns,
                                            trace,
                                            progress)

  which_peptide = ifelse(any(tolower(annotate_by_peptide) %in% c("unmodified", "um")),
                         data_peptide_column_no_PTMs,
                         data_peptide_column_PTMs)
  data[,data_accession_column] <- annotation_table[data[[which_peptide]], Proteins, on = "Peptides"]
  data <- get_signals(data,
                      data_peptide_column_no_PTMs,
                      incorporation_name,
                      intensity_name,
                      score_name,
                      labeling_ratio_name,
                      incorporation_columns,
                      intensity_columns,
                      labeling_ratio_columns,
                      score_columns,
                      as_percentage,
                      keep_feature_type,
                      trace)
  
  if(is.null(labeling_ratio_name)) {
    labeling_ratio_name = "LR"
  }
  
  if(cluster) {
    if (is.logical(group_by) && group_by) {
      if (peptide_centric) {
        group_by = c(data_accession_column, data_peptide_column_PTMs, data_peptide_column_no_PTMs, colnames(design), colnames(data)[tolower(colnames(data)) %like% "charge"])
      } else {
        group_by = c(data_accession_column, colnames(design))
      }
    }
    
    if(trace) {
      message("Feature clustering:")
      cat(crayon::silver(paste0("Grouping on c(", paste0("\"", group_by, "\"", collapse = ", "), ").\nClustering each group by c(", paste0("\"", cluster_by, "\"", collapse = ", "), ") with radius c(", paste0(radius, collapse = ", "), ")\n")))
    }
    
    data <- cluster_features(data,
                                  incorporation_name,
                                  intensity_name,
                                  labeling_ratio_name,
                                  score_name,
                                  incorporation_columns,
                                  intensity_columns,
                                  labeling_ratio_columns,
                                  score_columns,
                                  group_by,
                                  cluster_by,
                                  radius,
                                  distance_method,
                                  element_wise,
                                  trace,
                                  progress)
  }
  
  if(is.null(LFDR_by)) {
    LFDR_by = time_unit 
  }
  
  if(is.null(LFDR_features)) {
    LFDR_features = c(incorporation_name, labeling_ratio_name)
  }
  
  pdf_list = .get_pdf(data, design, time_unit, time_zero, LFDR_by, LFDR_features, bandwidth)
  d <- unique(design[get(time_unit) != time_zero, ..LFDR_by])
  d <- d[order(get(time_unit))]
  data$id <- seq_len(nrow(data))
  data$LFDR <- NA_real_
  
  for(i in 1:nrow(d)) {
    
    row <- data[d[i,], id, on = LFDR_by]
    
    x <- as.list(data[d[i,], ..LFDR_features, on = LFDR_by])
    
    pdfs = pdf_list[[paste0(d[i,], collapse = ",")]]
    
    pdf0 <- function(x) {
      mapply(function(a,b) a(b), pdfs$pdf0, x)
    }
    
    pdf <- function(x) {
      mapply(function(a,b) a(b), pdfs$pdf, x)
    }
    
    BF <- pdf0(x)/pdf(x)
    n <- nrow(BF)
    LFDR <- numeric(n)
    prior_prob = pdfs$prior_prob
    for (j in seq_len(n)) {
      LFDR[j] <- prior_prob * prod(BF[j,], na.rm = T)
    }
    data[row,"LFDR"] <- LFDR
  }
  data = data[,-"id"]
  
  timepoints = sort(unique(design[[time_unit]]))
  timepoints = timepoints[timepoints != time_zero]
  
  
  model <- function(x, var, take_mean = F, subset = NULL) {
    if(is.null(Object@consts[[var]])) {
      stop(paste0("Curve fitting was not performed on variable `", var, "`."))
    }
    eq <- unique(Object@consts[[var]]$equation)
    if(eq == 4) {
      consts <- Object@consts[[var]]
      value <- lapply(consts$k, predict, x)
      value <- as.data.table(do.call(rbind, lapply(consts$k, predict, x))) 
      value$Peptides <- rep(consts[, Peptides], each = length(x))
      value <- as.matrix(dcast(value, Peptides ~ z, value.var = "fit"), rownames = T)
    } else {
      if(is.null(subset)) {
        subset <- 1:nrow(Object@consts[[var]])
      }
      if(length(x) > 1 & !is.matrix(x))
      {
        x <- matrix(x, nrow = nrow(Object@consts[[var]][subset,]), ncol = length(x), byrow = T)
      }
      if(eq == 2) {
        kbi = Object@consts[[var]]$kbi[subset]
        k0a = Object@consts[[var]]$k0a[subset]
        conc = Object@consts[[var]]$conc[subset]
        v = kbi
        yv = kbi/(k0a - kbi)
        ykbi = k0a/(kbi - k0a)
        value = (1 + yv * exp(-v * x) + ykbi * exp(-kbi * x)) * conc
      } else if (eq == 3) {
        kbi = Object@consts[[var]]$kbi[subset]
        k0a = Object@consts[[var]]$k0a[subset]
        kst = Object@consts[[var]]$kst[subset]
        kbt = Object@consts[[var]]$kbt[subset]
        conc = Object@consts[[var]]$conc[subset]
        u = ((kst + k0a + kbt) - sqrt((kst + k0a + kbt)^2 - (4 * k0a*kbt))) / 2
        v = ((kst + k0a + kbt) + sqrt((kst + k0a + kbt)^2 - (4 * k0a*kbt))) / 2
        yu = k0a * kbi*(u - kbt) / ((u - v)*(u - kbi)*u)
        yv = k0a * kbi*(v - kbt) / ((v - u)*(v - kbi)*v)
        ykbi = k0a * (kbi - kbt) / ((u - kbi)*(v - kbi))
        value = (1 + yu * exp(-u * x) + yv * exp(-v * x) + ykbi * exp(-kbi * x))*conc
      }
    }
    if(take_mean) {
      value = colMeans(as.matrix(value))
    }
    return(value)
  }
  
  
  
  deconvolve <- function(features = c(Object@incorporation_name, Object@labeling_ratio_name),
                                xlim = NULL,
                                ylim = NULL,
                                filename = NULL,
                                plot = T,
                                width = NA,
                                height = NA)
  {
    p = list()
    i = 1
    for(f in features) {
      dt = rbindlist(lapply(Object@pdf, function(pdfs) {
        v = Object@data[,get(f)]
        rg = range(v)
        x = seq(rg[1], rg[2], length.out = 512)
        y0 = pdfs$pdf0[[f]](x) * pdfs$prior_prob
        y1 = pdfs$pdf1[[f]](x)
        y = pdfs$pdf[[f]](x)
        data.table(x = rep(x, 3), y = c(y1,y0,y), Distribution = rep(c("true discoveries", "false discoveries", "mixture density"), each = 512))
      }), id = "variable")
      tag = LETTERS[i]
      i=1+i
      dt$Distribution <- factor(dt$Distribution, levels = c("true discoveries", "false discoveries", "mixture density"))
      dt$variable <- factor(dt$variable, levels = names(Object@pdf))
      p[f] = list(ggplot(dt, aes(x = x, y = y, color = Distribution)) +
                    theme_bw() + xlab(paste(f, "(%)")) + labs(tag = tag) +
                    ylab("density") + geom_line(stat = "identity", size = 0.05) +
                    scale_color_manual(values = c("blue", "red", "black")) +
                    guides(colour = guide_legend(override.aes = list(size = 1))) +
                    scale_x_continuous(expand = c(0,0), limits = xlim[[f]]) +
                    scale_y_continuous(expand = expand_scale(mult = c(0, .1)), limits = ylim[[f]]) +
                    facet_wrap(vars(variable)))
    }
    if(length(p) > 1) {
      library(ggplotify)
      g <- ggplotGrob(p[[1]] + theme(legend.position="right"))$grobs
      legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
      lwidth <- sum(legend$width)
      p <- arrangeGrob(
        do.call(arrangeGrob, c(lapply(p, function(x)
          x + theme(legend.position="none")), list(nrow = 1))),
        legend,
        nrow = 1,
        widths = unit.c(unit(1, "npc") - lwidth, lwidth))
    } else {
      p <- p[[1]]
    }
    if(!is.null(filename)) {
      save_plot(filename = filename, plot = p, width = width, height = height)
    }
    if (plot) {
      grid.newpage()
      grid.draw(p)
    }
    p
  }
  
  
  Object <- new("MetaProfiler",
      design = design,
      data = data,
      incorporation_name = incorporation_name,
      intensity_name = intensity_name,
      labeling_ratio_name = labeling_ratio_name,
      score_name = ifelse(is.null(score_name), character(), score_name),
      peptide_column_no_PTMs = ifelse(is.null(data_peptide_column_no_PTMs), character(), data_peptide_column_no_PTMs),
      peptide_column_PTMs = ifelse(is.null(data_peptide_column_PTMs), character(), data_peptide_column_PTMs),
      protein_column = data_accession_column,
      peptide_centric = peptide_centric,
      annotate_by_peptide = which_peptide,
      master_tbl = data.table(),
      timepoints = timepoints,
      time_zero = time_zero,
      time_unit = time_unit,
      annotation_table = annotation_table,
      taxon_columns = pep2taxon_columns,
      function_columns = pro2func_function_columns,
      deconvolve = deconvolve,
      consts = list(),
      model = model,
      pdf = pdf_list)
}
  
load_results <- function(design,
                         data,
                         time_unit,
                         incorporation_name,
                         intensity_name,
                         labeling_ratio_name = NULL,
                         score_name = NULL,
                         incorporation_columns = "auto",
                         intensity_columns = "auto",
                         labeling_ratio_columns = NULL,
                         score_columns = NULL,
                         as_percentage = T,
                         progress = T,
                         ...) {
  design <- as.data.table(design)
  if(!any(colnames(design) == time_unit)) {
    stop(paste0("No column named `", time_unit, "` in the table `design`. Please add/rename the column for the timepoints or change variable `time_unit` so that it matches the timepoint column in table `design`."))
  }
  if(!is.numeric(design[[time_unit]])) {
    design[[time_unit]] <- as.numeric(design[[time_unit]])
    if(all(is.na(design[[time_unit]]))) {
      stop(paste0("Could not convert column `", time_unit, "` to numeric. Are there any non-numeric characters in the column?"))
    }
  }
  check = unlist(design[, lapply(.SD, function(x) {
    test <- try(all(file.exists(x)), silent = T)
    ifelse(inherits(test, "try-error"), F, test)
  })])
  if((sum(check) == 1) && !is.null(data)) {
    if(nrow(design) != length(data)) {
      stop("Make sure to provide the experimental conditions (e.g. sample names, timepoints) for all the files listed in `data`")
    }
    design <- cbind(data, design)
  } else {
    idx <- order(check, decreasing = T)
    design <- design[,idx,with=F]
  }
  reorder_cols = c(time_unit, setdiff(colnames(design)[-1], time_unit))
  design <- design[.mixedorder(design[,..reorder_cols]),]
  readr_params <- list(...)
  options("readr.show_progress" = F)
  f <- function(...)
  {
    args <- list(...)
    if(!file.exists(args[[1]])) {
      stop("Either the first column in `design` does not contain the filenames or that file \"", x, "\" does not exist.")
    }
    if(!any(names(readr_params) == "delim")) {
      readr_params$delim = "\t"
    }
    if(progress) {
      writeLines(paste0("reading ", gsub("[\\S\\s]+\\/", "", args[[1]], perl = T), " (", paste0(colnames(design)[-1], " = ", unlist(args)[-1], collapse = ", "),")..."))
    }
    con <- file(args[[1]], "r")
    headers <- readLines(con,n=1)
    close(con)
    headers = ifelse(grepl("\t", headers), strsplit(headers, "\t"),
                     ifelse(grepl(",", headers), strsplit(headers, ","),
                            strsplit(headers, " ")))[[1]]
    if(length(incorporation_columns) == 1 && incorporation_columns == "auto") {
      incorporation_columns = .check_columns(incorporation_name, headers)
    }
    if(!is.null(intensity_columns) && length(intensity_columns) == 1 && intensity_columns == "auto") {
      intensity_columns = .check_columns(intensity_name, headers)
    }
    if(!is.null(score_columns) && length(score_columns) == 1 && score_columns == "auto") {
      score_columns = .check_columns(score_name, headers)
    }
    if(!is.null(labeling_ratio_columns) && length(labeling_ratio_columns) == 1 && labeling_ratio_columns == "auto") {
      labeling_ratio_columns = .check_columns(labeling_ratio_name, headers)
    }
    if(is.null(intensity_columns) && is.null(labeling_ratio_columns)) {
      stop("Both the intensity and labeling ratio columns are NULL. At least one of these variables need to be specified.")
    }
    group = c(incorporation_columns, intensity_columns, score_columns, labeling_ratio_columns)
    group[!(group %in% make.names(group))] <- paste0('`', group[!(group %in% make.names(group))], '`')
    col_types <- eval(parse(text = paste("cols(", paste0(group, " = col_double()", collapse = ","), ")")))
    shush(data <- do.call(readr::read_delim, c(list(file = args[[1]], col_types = col_types), readr_params)))
    data <- as.data.table(setNames(cbind(args[-1], data), c(colnames(design)[-1], colnames(data))))
    data[, c(incorporation_columns, labeling_ratio_columns)] = data[, lapply(.SD, function(x) {
      if(.check_if_percantage(x)) {
        if(as_percentage) {
          return(x)
        } else {
          return(x * 100)
        }
      } else {
        if(as_percentage) {
          return(x * 100)
        } else {
          return(x)
        }
      }
    }), .SDcols = c(incorporation_columns, labeling_ratio_columns)]
    data
  }
  args <- c(list(FUN = f), as.list(design), list(SIMPLIFY = F))
  data <- rbindlist(do.call(mapply, args))
  options("readr.show_progress" = T)
  data
}

get_signals <- function(data,
                        incorporation_name,
                        intensity_name,
                        score_name,
                        labeling_ratio_name,
                        peptide_column = "guess",
                        incorporation_columns = "auto",
                        intensity_columns = "auto",
                        labeling_ratio_columns = "auto",
                        score_columns = "auto",
                        as_percentage = T,
                        feature = c("feature", "id"),
                        spike_in = T,
                        trace = T) {
  data = as.data.table(data)
  if(peptide_column == "guess") {
    peptide_column = .guess_unmodified_peptide_column(data, trace)
    if(is.list(peptide_column)) {
      data = peptide_column[[1]]
      peptide_column = peptide_column[[2]]
    }
  }
  if(!is.null(data$Feature)) {
    feature <- match.arg(feature, c("feature", "id"), several.ok = T)
    data <- data[feature, , on = "Feature"]
  }
  headers = colnames(data)
  if(length(incorporation_columns) == 1 && incorporation_columns == "auto") {
    incorporation_columns = .check_columns(incorporation_name, headers, F)
  }
  if(!is.null(intensity_columns) && length(intensity_columns) == 1 && intensity_columns == "auto") {
    intensity_columns = .check_columns(intensity_name, headers, F)
  }
  if(!is.null(score_columns) && length(score_columns) == 1 && score_columns == "auto") {
    score_columns = .check_columns(score_name, headers, F)
  }
  if(!is.null(labeling_ratio_columns) && length(labeling_ratio_columns) == 1 && labeling_ratio_columns == "auto") {
    labeling_ratio_columns = .check_columns(labeling_ratio_name, headers, F)
  }
  seq <- unique(data[[peptide_column]])
  seq <- setNames(strsplit(seq, "", fixed = T), seq)
  seq <- lapply(seq, function(x) setNames(as.list(x), as.character(seq_along(x))))
  seq <- rbindlist(seq, fill = T, idcol = "pep")
  seq[,colnames(seq)[-1]] <- seq[,lapply(.SD[,-1], function(x) AA_tbl[x, N, on = "AA"])]
  seq <- seq[, sum(unlist(.SD), na.rm = T), by = pep]
  data$isotopes <- seq[data[[peptide_column]], V1, on = "pep"]
  incorporation_value <- data[,..incorporation_columns]
  if(.check_if_percantage(incorporation_value)) {
    incorporation_value = incorporation_value/100
  }
  light_cols <- c(incorporation_columns[1], intensity_columns[1], score_columns[1], labeling_ratio_columns[1])
  heavy_cols <- list(incorporation_columns[-1], intensity_columns[-1], score_columns[-1], labeling_ratio_columns[-1])
  keep = !sapply(heavy_cols, function(x) {
    length(x) == 0
  })
  heavy_cols = heavy_cols[keep]
  value.name = unlist(list(incorporation_name, intensity_name, score_name, labeling_ratio_name)[keep])
  test <- rowSums(incorporation_value > (1/data$isotopes), na.rm = T) > 0
  if(spike_in) {
    test <- test & (rowSums(incorporation_value[,1] <= (1/data$isotopes), na.rm = T) > 0)
  }
  data <- data[test]
  
  if(is.null(intensity_columns) && is.null(labeling_ratio_columns)) {
    stop("Both the intensity and labeling ratio columns are NULL. At least one of these variables need to be specified.")
  }
  group <- colnames(data)[!(colnames(data) %in% c(incorporation_columns, intensity_columns, score_columns, labeling_ratio_columns))]
  col <-  colnames(data)[colnames(data) %in% c(group, light_cols)]
  data <- melt(data, id.vars = col, measure = heavy_cols, value.name = value.name, na.rm = T)[,-"variable"]
  if(length(labeling_ratio_columns) == 1) {
    colnames(data)[colnames(data) == labeling_ratio_columns] = labeling_ratio_name
  }
  if(is.null(labeling_ratio_columns)) {
    if(trace) {
      cat(crayon::blue(paste0("Computing labeling ratio and inserting it into column `LR`.\n")))
    }
    LR = data[[intensity_name]]/(data[[intensity_name]] + data[[intensity_columns[1]]])
    if(as_percentage) {
      LR = LR * 100
    }
    data[,"LR"] <- LR
  }
  data$N <- 1
  data
}

cluster_features <- function(data,
                             incorporation_name,
                             intensity_name,
                             labeling_ratio_name,
                             score_name,
                             incorporation_columns = "auto",
                             intensity_columns = "auto",
                             labeling_ratio_columns = NULL,
                             score_columns = NULL,
                             group_by,
                             cluster_by = c("RIA", "LR", "RT"),
                             radius = c(1,1,0.5),
                             distance_method = c("absolute", "absolute", "absolute"),
                             element_wise = T,
                             trace = T,
                             progress = T) {
  data = as.data.table(data)
  headers = colnames(data)
  if(length(incorporation_columns) == 1 && incorporation_columns == "auto") {
    incorporation_columns = .check_columns(incorporation_name, headers)
  }
  if(!is.null(intensity_columns) && length(intensity_columns) == 1 && intensity_columns == "auto") {
    intensity_columns = .check_columns(intensity_name, headers)
  }
  if(!is.null(score_columns) && length(score_columns) == 1 && score_columns == "auto") {
    score_columns = .check_columns(score_name, headers)
  }
  if(is.null(intensity_columns) && is.null(labeling_ratio_columns)) {
    stop("Both the intensity and labeling ratio columns are NULL. At least one of these variables need to be specified.")
  }
  if(progress) {
    writeLines("clustering features...")
  }
  clust <- qtclust(data,
                   cluster_by = cluster_by,
                   radius = radius,
                   distance_method = distance_method,
                   centers = unique(c(cluster_by, incorporation_name, incorporation_columns[1],
                                      score_name, score_columns[1], intensity_name, intensity_columns[1], labeling_ratio_name)),
                   group_by = group_by,
                   element_wise = element_wise, 
                   progress = progress)
  if(trace) {
    message(paste0("clustered ", nrow(data), " features to ", nrow(clust$centers), "."))
  }
  clust$centers  
}




