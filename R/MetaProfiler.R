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
                         
                         observations = NULL,
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
  
  args = guess_columns(data,
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
  data <- flatten_data(data,
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
  
  if(is.null(observations)) {
    observations = c(incorporation_name, labeling_ratio_name)
  }
  
  pdf_list = get_pdf(data, design, time_unit, time_zero, LFDR_by, observations, bandwidth)
  data <- calc_LFDR(data, design, time_unit, time_zero, LFDR_by, observations, pdf_list, bandwidth)
  
  timepoints = sort(unique(design[[time_unit]]))
  timepoints = timepoints[timepoints != time_zero]
  
  
  Object <- new(
    "MetaProfiler",
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
    pdf = pdf_list
  )
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
  

flatten_data <- function(data,
                         incorporation_name,
                         intensity_name,
                         score_name,
                         labeling_ratio_name,
                         peptide_column = "guess",
                         incorporation_columns = guess_measurements(incorporation_name, colnames(data)),
                         intensity_columns = guess_measurements(intensity_name, colnames(data)),
                         labeling_ratio_columns = guess_measurements(labeling_ratio_name, colnames(data)),
                         score_columns = guess_measurements(score_name, colnames(data)),
                         heavy_columns = list(incorporation_columns[-1],
                                              intensity_columns[-1],
                                              score_columns[-1],
                                              labeling_ratio_columns[-1]),
                         LR_as_percentage = T,
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
    if(LR_as_percentage) {
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
                             cluster_by = c(incorporation_name, labeling_ratio_name),
                             radius = c(0.1,0.1),
                             distance_method = c("relative", "relative"),
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

