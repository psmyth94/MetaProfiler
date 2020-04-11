setClass(
  "MetaProfiler",
  slots = c(
    design = "data.table",
    data = "data.table",
    timepoints="numeric",
    time_zero = "numeric",
    time_unit = "character",
    isotope = "character",
    incorporation_name = "character",
    incorporation_columns = "character",
    intensity_name = "character",
    intensity_columns = "character",
    labeling_ratio_name = "character",
    labeling_ratio_columns = "character",
    score_name = "character",
    score_columns = "character",
    higher_score_better = "logical",
    peptide_column_no_PTMs = "character",
    peptide_column_PTMs = "character",
    accession_column = "character",
    peptide_centric = "logical",
    annotate_with = "character",
    master_tbl = "data.table",
    annotation_table = "data.table",
    taxon_columns = "character",
    function_columns = "character",
    deconvolve = "function",
    consts = "list",
    model = "function",
    pdf = "list",
    get_cols = "function"
  ),
  prototype = list(
    design = data.table::data.table(),
    data = data.table::data.table(),
    timepoints=numeric(),
    time_zero = numeric(),
    time_unit = character(),
    isotope = character(),
    incorporation_name = character(),
    incorporation_columns = character(),
    intensity_name = character(),
    intensity_columns = character(),
    labeling_ratio_name = character(),
    labeling_ratio_columns = character(),
    score_name = character(),
    score_columns = character(),
    peptide_column_no_PTMs = character(),
    peptide_column_PTMs = character(),
    accession_column = character(),
    peptide_centric = logical(),
    annotate_with = character(),
    master_tbl = data.table::data.table(),
    annotation_table = data.table::data.table(),
    taxon_columns = character(),
    function_columns = character(),
    deconvolve = function(observations = c(Object@incorporation_name, Object@labeling_ratio_name),
                          xlim = NULL,
                          ylim = NULL,
                          filename = NULL,
                          plot = T,
                          width = NA,
                          height = NA)
    {
      data = Object@data[,..observations]
      v = data[,lapply(.SD, function(x) {
        seq(min(x), max(x), length.out = 512)
      })]
      for(pdfs in Object@pdf) {
        k <- c(mapply(function(x, d, p1, p0, o) {
          ggplotify::as.grob(function() {
            h<-hist(rep(d, mixture$N), breaks=20, xlab = o, main = NULL)
            xfit<- x
            yfit0 <- p0(x)
            yfit0 <- yfit0*diff(h$mids[1:2])* sum(data$N) * pdfs$prior_prob
            yfit1 <- p1(x)
            yfit1 <- yfit1*diff(h$mids[1:2])* sum(data$N) * (1 - pdfs$prior_prob)
            yfite <- yfit1 + yfit0
            lines(xfit, yfit1, col="blue", lty = "dashed", lwd=2)
            lines(xfit, yfit0, col="red", lty = "dashed", lwd=2)
            lines(xfit, yfite, col="black", lwd=2)
          })
        },v, data, pdfs$pdf1, pdfs$pdf0, observations, SIMPLIFY = F),
        list(
          ggplotify::as.grob(function() {
            plot.new()
            legend("center", c("true", "false", "mixture"), lty = c(2,2,1), col = c("blue", "red", "black"))
          })
        ))
      }
      if(length(filename)) {
        save_plot(filename = filename, plot = p, width = width, height = height)
      }
      if (plot) {
        grid::grid.newpage()
        grid::grid.draw(p)
      }
      p
    },
    consts = list(),
    model = function(x, var, take_mean = F, subset = NULL) {
      if(!length(Object@consts[[var]])) {
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
        if(!length(subset)) {
          subset <- 1:nrow(Object@consts[[var]])
        }
        if(length(x) > 1 && !is.matrix(x))
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
    },
    pdf = list(),
    get_cols = function(var, x = Object@timepoints, time_unit = Object@time_unit) {
      nt = length(x)
      nv = length(var)
      paste0(rep(var, each = nt), " (", time_unit, " ", rep(x, nv), ")")
    }
  )
)

MetaProfiler <- function(design, data, time_unit, time_zero, isotope = "N", peptide_centric = T,
                         incorporation_name, intensity_name, labeling_ratio_name = NULL, score_name = NULL, higher_score_better = T,
                         incorporation_columns = "auto", intensity_columns = "auto", labeling_ratio_columns = "auto", score_columns = "auto",
                         peptide_column_PTMs = "guess", peptide_column_no_PTMs = "guess", accession_column = "guess", annotate_by_peptide = c("unmodified", "um", "modified", "m"),
                         feature_type_column = NULL, feature_type = c("feature", "id"),
                         as_percentage = T, progress = T, trace = T,
                         ...) {
  annotate_by_peptide = match.arg(annotate_by_peptide)
  Object <- new("MetaProfiler")
  Object@design <- data.table::as.data.table(design)
  Object@time_unit <- time_unit
  Object@time_zero <- time_zero
  Object@higher_score_better = higher_score_better
  Object@peptide_centric = peptide_centric
  Object@isotope = isotope
  if(length(incorporation_name))
    Object@incorporation_name = incorporation_name
  if(length(intensity_name))
    Object@intensity_name = intensity_name
  if(length(labeling_ratio_name))
    Object@labeling_ratio_name = labeling_ratio_name
  if(length(score_name))
    Object@score_name = score_name
  if(length(incorporation_columns) && (incorporation_columns != "auto"))
    Object@incorporation_columns = incorporation_columns
  if(length(intensity_columns) && (intensity_columns != "auto"))
    Object@intensity_columns = intensity_columns
  if(length(labeling_ratio_columns) && (labeling_ratio_columns != "auto"))
    Object@labeling_ratio_columns = labeling_ratio_columns
  if(length(score_columns) && (score_columns != "auto"))
    Object@score_columns = score_columns
  
  if(!any(colnames(Object@design) == Object@time_unit)) {
    stop(paste0("No column named `", Object@time_unit, "` in the table `design`. Please add/rename the column for the timepoints or change variable `time_unit` so that it matches the timepoint column in table `design`."))
  }
  if(!is.numeric(Object@design[[Object@time_unit]])) {
    Object@design[[Object@time_unit]] <- as.numeric(Object@design[[Object@time_unit]])
    if(all(is.na(Object@design[[Object@time_unit]]))) {
      stop(paste0("Could not convert column `", Object@time_unit, "` to numeric. Are there any non-numeric characters in the column?"))
    }
  }
  check = unlist(Object@design[, lapply(.SD, function(x) {
    test <- try(all(file.exists(x)), silent = T)
    ifelse(inherits(test, "try-error"), F, test)
  })])
  if(!inherits(data, "data.frame")) {
    if(length(data) && is.character(data)) {
      if (nrow(Object@design) != length(data)) {
        stop("Make sure to provide the experimental conditions (e.g. sample names, timepoints) for all the files listed in `data`")
      } else {
        Object@design <- cbind(data, Object@design)
      }
    } else if (sum(check) == 1) {
      idx <- order(check, decreasing = T)
      Object@design <- Object@design[,idx,with=F]
    } else if (any(check)) {
      stop("More than one column in table `design` contains filepaths. Initialize `data` with the appropriate filepaths instead.")
    }
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
        writeLines(paste0("reading ", gsub("[\\S\\s]+\\/", "", args[[1]], perl = T), " (", paste0(colnames(Object@design)[-1], " = ", unlist(args)[-1], collapse = ", "),")..."))
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
      if(length(intensity_columns) == 1 && intensity_columns == "auto") {
        intensity_columns = .check_columns(intensity_name, headers)
      }
      if(length(score_columns) == 1 && score_columns == "auto") {
        score_columns = .check_columns(score_name, headers)
      }
      if(length(labeling_ratio_columns) == 1 && labeling_ratio_columns == "auto") {
        labeling_ratio_columns = .check_columns(labeling_ratio_name, headers)
      }
      if(!length(intensity_columns) && !length(labeling_ratio_columns)) {
        stop("Both the intensity and labeling ratio columns are NULL. At least one of these variables need to be specified.")
      }
      group = c(incorporation_columns, intensity_columns, score_columns, labeling_ratio_columns)
      group[!(group %in% make.names(group))] <- paste0('`', group[!(group %in% make.names(group))], '`')
      col_types <- eval(parse(text = paste("readr::cols(", paste0(group, " = readr::col_double()", collapse = ","), ")")))
      shush(data <- do.call(readr::read_delim, c(list(file = args[[1]], col_types = col_types), readr_params)))
      data <- data.table::as.data.table(setNames(cbind(args[-1], data), c(colnames(design)[-1], colnames(data))))
      data[, c(incorporation_columns, labeling_ratio_columns)] = data[, lapply(.SD, function(x) {
        is_percentage <- .check_if_percantage_or_fraction(x)
        if(!length(is_percentage)) {
          return(as.numeric(x))
        } else if (is_percentage) {
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
    args <- c(list(FUN = f), as.list(Object@design), list(SIMPLIFY = F))
    data <- data.table::rbindlist(do.call(mapply, args))
  }
  res <- .guess_columns(data, data_peptide_column_PTMs = peptide_column_PTMs,
                        data_peptide_column_no_PTMs = peptide_column_no_PTMs,
                        data_accession_column = accession_column, trace = trace)
  Object@data = res$data
  if(length(res$data_peptide_column_PTMs))
    Object@peptide_column_PTMs = res$data_peptide_column_PTMs
  if(length(res$data_peptide_column_no_PTMs))
    Object@peptide_column_no_PTMs = res$data_peptide_column_no_PTMs
  if(length(res$data_accession_column))
    Object@accession_column = res$data_accession_column
  
  which_peptide = ifelse(
    any(tolower(annotate_by_peptide) %in% c("unmodified", "um")),
    Object@peptide_column_no_PTMs,
    Object@peptide_column_PTMs
  )
  Object@annotate_with = which_peptide
  Object@design <- Object@design[,!check, with = F]
  
  headers = colnames(data)
  if(length(incorporation_columns) == 1 && incorporation_columns == "auto") {
    Object@incorporation_columns = .check_columns(Object@incorporation_name, headers)
  }
  if(length(intensity_columns) == 1 && intensity_columns == "auto") {
    Object@intensity_columns = .check_columns(Object@intensity_name, headers)
  }
  if(length(score_columns) == 1 && score_columns == "auto") {
    Object@score_columns = .check_columns(Object@score_name, headers)
  }
  if(length(labeling_ratio_columns) == 1 && labeling_ratio_columns == "auto") {
    Object@labeling_ratio_columns = .check_columns(Object@labeling_ratio_name, headers)
  }
  
  if(length(feature_type_column)) {
    Object@data <- Object@data[feature_type, , on = feature_type_column]
  }
  
  Object@timepoints = sort(unique(Object@data[,get(Object@time_unit)]))
  Object@timepoints = Object@timepoints[Object@timepoints != Object@time_zero]
  options("readr.show_progress" = T)
  Object
}

get_data <- function(Object, labeling_ratio_as_percentage = T, light_peptide = T,
                         trace = T, score_threshold = NA, labeling_ratio_threshold = NA) {
  Object@data = data.table::copy(Object@data)
  light_cols <- c(Object@incorporation_columns[1], Object@intensity_columns[1], Object@score_columns[1], Object@labeling_ratio_columns[1])
  heavy_cols <- list(Object@incorporation_columns[-1], Object@intensity_columns[-1], Object@score_columns[-1], Object@labeling_ratio_columns[-1])
  keep = !sapply(heavy_cols, function(x) {
    length(x) == 0
  })
  heavy_cols = heavy_cols[keep]
  value.name = unlist(list(Object@incorporation_name, Object@intensity_name, Object@score_name, Object@labeling_ratio_name)[keep])
  seq <- unique(Object@data[[Object@annotate_with]])
  seq <- setNames(strsplit(seq, "", fixed = T), seq)
  seq <- lapply(seq, function(x) setNames(as.list(x), as.character(seq_along(x))))
  seq <- data.table::rbindlist(seq, fill = T, idcol = "pep")
  seq[,colnames(seq)[-1]] <- seq[,lapply(.SD[,-1], function(x) AA_tbl[x, , on = "AA"][[Object@isotope]])]
  seq <- seq[, sum(unlist(.SD), na.rm = T), by = pep]
  Object@data$isotopes <- seq[Object@data[[Object@annotate_with]], V1, on = "pep"]
  incorporation_value <- Object@data[,Object@incorporation_columns,with=F]
  is_percentage = .check_if_percantage_or_fraction(incorporation_value)
  if(is_percentage) {
    incorporation_value <- incorporation_value/100
  }
  test <- rowSums(incorporation_value > (1/Object@data$isotopes), na.rm = T) > 0
  if(!light_peptide) {
    swap <- is.na(incorporation_value[[2]]) & test
    incorporation_value[[2]][swap] = incorporation_value[[1]][swap]
    incorporation_value[[1]][swap] = 0
    if(is_percentage)
      Object@data[,Object@incorporation_columns] <- incorporation_value * 100
    else
      Object@data[,Object@incorporation_columns] <- incorporation_value
    if(length(Object@intensity_columns)){
      Object@data[[Object@intensity_columns[2]]][swap] = Object@data[[Object@intensity_columns[1]]][swap]
      Object@data[[Object@intensity_columns[1]]][swap] = 0
    }
    if(length(Object@score_columns)){
      Object@data[[Object@score_columns[2]]][swap] = Object@data[[Object@score_columns[1]]][swap]
      Object@data[[Object@score_columns[1]]][swap] = 0
    }
    if(length(Object@labeling_ratio_columns)){
      Object@data[[Object@labeling_ratio_columns[2]]][swap] = Object@data[[Object@labeling_ratio_columns[1]]][swap]
      Object@data[[Object@labeling_ratio_columns[1]]][swap] = 0
    }
    if(length(Object@score_columns)){
      Object@data[[Object@score_columns[2]]][swap] = Object@data[[Object@score_columns[1]]][swap]
      if(Object@higher_score_better) {
        Object@data[[Object@score_columns[1]]][swap] = max(Object@data[,Object@score_columns,with=F],na.rm = T)
      } else {
        Object@data[[Object@score_columns[1]]][swap] = min(Object@data[,Object@score_columns,with=F],na.rm = T)
      }
    }
  }
  test <- test & (rowSums(incorporation_value[,1] <= (1/Object@data$isotopes), na.rm = T) > 0)
  Object@data <- Object@data[test]
  measurements <- c(Object@incorporation_columns, Object@intensity_columns, Object@score_columns, Object@labeling_ratio_columns)
  group <- colnames(Object@data)[!(colnames(Object@data) %in% c(Object@incorporation_columns, Object@intensity_columns, Object@score_columns, Object@labeling_ratio_columns))]
  col <-  colnames(Object@data)[colnames(Object@data) %in% c(group, light_cols)]
  Object@data[, measurements] <- Object@data[, lapply(.SD, as.numeric), .SDcols = measurements]
  Object@data <- data.table::melt(Object@data, id.vars = col, measure = heavy_cols, value.name = value.name, na.rm = T)[,-"variable"]
  if(length(Object@labeling_ratio_columns) == 1) {
    if(length(Object@labeling_ratio_name)) {
      colnames(Object@data)[colnames(Object@data) == Object@labeling_ratio_columns] = Object@labeling_ratio_name
    }
    else {
      colnames(Object@data)[colnames(Object@data) == Object@labeling_ratio_columns] = "LR"
      Object@labeling_ratio_name = "LR"
    }
  }
  if(!length(Object@labeling_ratio_columns)) {
    if(trace) {
      cat(crayon::blue(paste0("Computing labeling ratio and inserting it into column `LR`.\n")))
    }
    LR = Object@data[[Object@intensity_name]]/(Object@data[[Object@intensity_name]] + Object@data[[Object@intensity_columns[1]]])
    if(labeling_ratio_as_percentage) {
      LR = LR * 100
    }
    Object@data[,"LR"] <- LR
    Object@labeling_ratio_name = "LR"
  }
  
  if(!is.na(score_threshold)) {
    if(Object@higher_score_better) {
      Object@data = Object@data[(get(Object@score_name) >= score_threshold) & (get(Object@score_columns[1]) >= score_threshold)]
    } else {
      Object@data = Object@data[(get(Object@score_name) <= score_threshold) & (get(Object@score_columns[1]) <= score_threshold)]
    }
  }
  if(!is.na(labeling_ratio_threshold)) {
    Object@data = Object@data[get(Object@labeling_ratio_name) >= labeling_ratio_threshold]
  }
  Object@data$N <- 1
  Object
}





