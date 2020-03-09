create_experimental_design <- function(filenames, results_directory, ...) {
  if(missing(filenames)) {
    filenames <- dir(result_directory, pattern = c("tsv|TSV|txt|TXT|csv|CSV"), full.names = T)
  }
  args <- list(...)
  data.table(filenames = filenames, as.data.table(lapply(args, function(x) {
    stringi::stri_extract_last_regex(filenames, x)
  })))
}
get_cols <- function(var, x = Object@timepoints, time_unit = Object@time_unit) {
  nt = length(x)
  nv = length(var)
  paste0(rep(var, each = nt), " (", time_unit, " ", rep(x, nv), ")")
}

get_ranks <- function(x) {
  colnames(x)[tolower(colnames(x)) %in% 
                c("superkingdom", "kingdom", "subkingdom",
                  "superphylum", "phylum", "subphylum",
                  "superclass", "class", "subclass", "infraclass",
                  "superorder", "order", "suborder", "infraorder", "parvorder",
                  "superfamily", "family", "subfamily",
                  "tribe", "subtribe",
                  "genus", "subgenus",
                  "species group", "species subgroup", "species", "subspecies",
                  "varietas", "forma")]
}

shush <- function(expr, all = TRUE) {
  if (Sys.info()['sysname'] == "Windows") {
    file <- "NUL"
  } else {
    file <- "/dev/null"
  }
  
  if (all) {
    suppressWarnings(suppressMessages(suppressPackageStartupMessages(
      utils::capture.output(expr, file = file) 
    )))
  } else {
    utils::capture.output(expr, file = file)
  }
}


save_plot <- function (filename, plot = last_plot(), device = NULL, path = NULL, 
                       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
                       dpi = 300, limitsize = TRUE, ...) {
  dpi <- ggplot2:::parse_dpi(dpi)
  dev <- ggplot2:::plot_dev(device, filename, dpi = dpi)
  dim <- ggplot2:::plot_dim(c(width, height), scale = scale, units = units, 
                            limitsize = limitsize)
  if (!is.null(path)) {
    filename <- file.path(path, filename)
  }
  old_dev <- grDevices::dev.cur()
  dev(filename = filename, width = dim[1], height = dim[2])
  on.exit(utils::capture.output({
    grDevices::dev.off()
    if (old_dev > 1) grDevices::dev.set(old_dev)
  }))
  if(inherits(plot, c("Heatmap", "AnnotationFunction", "Legends", "HeatmapAnnotation", "SingleAnnotation", "HeatmapList"))) {
    draw(plot, ...)
  } else if (inherits(plot, "igraph")) {
    plot(plot, ...)
  } else {
    grid.draw(plot)
  }
  invisible()
}

.guess_accession_column <- function(d, trace = T, pattern = "[-_]") {
  cols <- toupper(colnames(d))
  suppressWarnings(idx <- which(cols %like% "ACC|PRO|QUERY|NAME|GENE|LEAD|RAZOR|GROUP"))
  if(!all(is.finite(idx))) {
    warning(paste0("Could not guess the name of the accession/protein column for `", deparse(substitute(d)), "`. Please provide the name of this column."))
    return(NULL)
  }
  lookahead <- NULL
  for(i in idx) {
    c = d[,get(colnames(d)[i])]
    if(is.character(c)) {
      lookahead <- c(lookahead, sum(stringi::stri_count_regex(d[[colnames(d)[i]]], pattern)/nchar(d[[colnames(d)[i]]]), na.rm = T))
    } else {
      lookahead <- c(lookahead, -Inf)
    }
  }
  idx <- idx[which.max(lookahead)]
  if(any(sapply(lookahead, is.finite))) {
    if(trace) {
      cat(crayon::silver(paste0("Using `", colnames(d)[idx], "` as the accession/protein column.\n")))
    }
    return(colnames(d)[idx])
  } else {
    warning(paste0("Could not guess the name of the accession/protein column for `", deparse(substitute(d)), "`. Please provide the name of this column."))
    return(NULL)
  }
}

.convert_modified_peptide_column <- function(d, col, trace = T) {
  if(trace) {
    cat(crayon::silver(paste0("Looks like `", col, "` is a peptide column, but appears to contain post translational modifications.\n")))
  }
  new_column = paste(col, "(no PTMs)")
  if(trace) {
    cat(crayon::blue(paste0("Creating column `", new_column,"` with sequences from column `", col, "`.\n")))
  }
  d[,new_column] = d[[col]]
  if(trace) {
    cat(crayon::blue(paste0("Removing PTMs in column `", new_column, "`...\n")))
  }
  convert <- unique(d[,new_column,with=F])
  convert[,"seq"] <- .remove_modifications(convert[[new_column]])
  
  if(all(!grepl("[^ILKMFTWVRHANDCEQGPSY]", convert$seq, perl = T)))
  {
    d[,new_column] <- convert[d, seq, on = new_column]
  } else {
    cat(crayon::black(paste0("Failed to remove PTMs, returning NULL...\n")))
    return(NULL)
  }
  if(trace) {
    cat(crayon::blue("done.\n"))
  }
  list(d, new_column)
}

.guess_unmodified_peptide_column <- function(d, trace = T) {
  cols <- toupper(colnames(d))
  suppressWarnings(idx <- which(cols %like% "PEP|SEQ"))
  if(!all(is.finite(idx))) {
    warning(paste0("Could not guess the name of the unmodified peptide column for `", deparse(substitute(d)), "`. Please provide the name of this column."))
    return(NULL)
  }
  lookahead <- NULL
  for(i in idx) {
    c = d[,get(colnames(d)[i])]
    if(is.character(c)) {
      lookahead <- c(lookahead, sum(!grepl("[^ILKMFTWVRHANDCEQGPSY]", c, perl = T)))
    } else {
      lookahead <- c(lookahead, 0)
    }
  }
  idx2 <- idx[lookahead == nrow(d)]
  if(length(idx2) > 0) {
    idx <- idx[which.max(lookahead > 0)]
    if(trace) {
      cat(crayon::silver(paste0("Using `", colnames(d)[idx], "` as the peptide column.\n")))
    }
    return(colnames(d)[idx])
  } else {
    idx2 <- idx[lookahead > nrow(d) * 0.5]
    if(length(idx2) > 0) {
      idx <- idx[which.max(lookahead > 0)]
      return(.convert_modified_peptide_column(d, colnames(d)[idx], trace))
    }
    warning(paste0("Could not guess the name of the unmodified peptide column for `", deparse(substitute(d)), "`. Please provide the name of this column."))
    return(NULL)
  }
}

.guess_modified_peptide_column <- function(d, trace = T) {
  cols <- toupper(colnames(d))
  suppressWarnings(idx <- which(cols %like% "PEP|SEQ"))
  if(!all(is.finite(idx))) {
    warning(paste0("Could not guess the name of the modified peptide column for `", deparse(substitute(d)), "`. Please provide the name of this column."))
    return(NULL)
  }
  
  lookahead <- NULL
  for(i in idx) {
    c = d[,get(colnames(d)[i])]
    if(is.character(c)) {
      c2 <- .remove_modifications(c)
      lookahead <- c(lookahead, list(c(sum(!grepl("[^ILKMFTWVRHANDCEQGPSY]", c, perl = T)), sum(grepl("[\\(\\)\\.\\{\\}\\<\\>\\[\\]\\_]", c, perl = T)))))
    } else {
      lookahead <- c(lookahead,list(c(-Inf, -Inf)))
    }
  }
  idx <- idx[sapply(lookahead, function(x) (x[1] > 0) && (x[2] > 0))]
  if(length(idx) > 0) {
    if(trace) {
      cat(crayon::silver(paste0("Using `", colnames(d)[idx][1], "` as the peptide column.\n")))
    }
    colnames(d)[idx][1]
  } else {
    warning(paste0("Could not guess the name of the modified peptide column for `", deparse(substitute(d)), "`. Please provide the name of this column."))
    return(NULL)
  }
}

.get_columns <- function(x, cols) {
  x = gsub('([\\[!@#$%^&*(),.?":{}|<>\\]])', "\\\\\\1", x, perl = T)
  cols[grepl(paste0("^", x,"\\s?\\S+"), cols)]
}

.check_columns <- function(x, cols) {
  if(length(x) == 0 || is.na(x)) return(NULL)
  cols <- .get_columns(x, cols)
  if(length(unique(cols)) == 0)
  {
    stop(paste0("No valid columns that contains the variable `", x,"`. Please add columns where the first word is the variable's name, followed by a unique identifier (e.g. ",
                paste0(x, " ", 1:3, ", ", collapse = ""), "..., or ", paste0(x, " ", c("light", "heavy"), ", ", collapse = ""),")."))
  }
  cols
}

# Removes PTMs by deleting characters within brackets and by deleting every character that is not a one-letter code for an amino acid...
.remove_modifications <- function(x) { 
  stri_replace_all_regex(x, "\\([\\S\\s]+?\\)|\\[[\\S\\s]+?\\]|\\{[\\S\\s]+?\\}|\\<[\\S\\s]+?\\>|_[\\S\\s]+?_|[^ILKMFTWVRHANDCEQGPSY]", "")
}

.check_if_percantage <- function(x) {
  shush(v <- max(x, na.rm = T))
  if(v > 1 && v <= 100) {
    return(T)
  }
  F
}

.mixedorder <- function (data, decreasing = FALSE, na.last = TRUE, blank.last = FALSE, 
                         numeric.type = c("decimal", "roman"), roman.case = c("upper", 
                                                                              "lower", "both")) 
{
  numeric.type <- match.arg(numeric.type)
  roman.case <- match.arg(roman.case)
  data = as.list(data)
  order.frame <- as.list(do.call(cbind, lapply(data, function(x) {
    if (length(x) < 1) 
      return(NULL)
    else if (length(x) == 1) 
      return(1)
    if (!is.character(x)) 
      return(order(x, decreasing = decreasing, na.last = na.last))
    delim = "\\$\\@\\$"
    if (numeric.type == "decimal") {
      regex <- "((?:(?i)(?:[-+]?)(?:(?=[.]?[0123456789])(?:[0123456789]*)(?:(?:[.])(?:[0123456789]{0,}))?)(?:(?:[eE])(?:(?:[-+]?)(?:[0123456789]+))|)))"
      numeric <- function(x) as.numeric(x)
    }
    else if (numeric.type == "roman") {
      regex <- switch(roman.case, both = "([IVXCLDMivxcldm]+)", 
                      upper = "([IVXCLDM]+)", lower = "([ivxcldm]+)")
      numeric <- function(x) roman2int(x)
    }
    else stop("Unknown value for numeric.type: ", numeric.type)
    nonnumeric <- function(x) {
      ifelse(is.na(numeric(x)), toupper(x), NA)
    }
    x <- as.character(x)
    which.nas <- which(is.na(x))
    which.blanks <- which(x == "")
    delimited <- gsub(regex, paste(delim, "\\1", delim, sep = ""), 
                      x, perl = TRUE)
    step1 <- strsplit(delimited, delim)
    step1 <- lapply(step1, function(x) x[x > ""])
    suppressWarnings(step1.numeric <- lapply(step1, numeric))
    suppressWarnings(step1.character <- lapply(step1, nonnumeric))
    maxelem <- max(sapply(step1, length))
    step1.numeric.t <- lapply(1:maxelem, function(i) sapply(step1.numeric, 
                                                            function(x) x[i]))
    step1.character.t <- lapply(1:maxelem, function(i) sapply(step1.character, 
                                                              function(x) x[i]))
    rank.numeric <- sapply(step1.numeric.t, rank)
    rank.character <- sapply(step1.character.t, function(x) as.numeric(factor(x)))
    rank.numeric[!is.na(rank.character)] <- 0
    rank.character <- t(t(rank.character) + apply(matrix(rank.numeric), 
                                                  2, max, na.rm = TRUE))
    rank.overall <- ifelse(is.na(rank.character), rank.numeric, 
                           rank.character)
    order.frame <- as.data.frame(rank.overall)
    if (length(which.nas) > 0) 
      if (is.na(na.last)) 
        order.frame[which.nas, ] <- NA
    else if (na.last) 
      order.frame[which.nas, ] <- Inf
    else order.frame[which.nas, ] <- -Inf
    if (length(which.blanks) > 0) 
      if (is.na(blank.last)) 
        order.frame[which.blanks, ] <- NA
    else if (blank.last) 
      order.frame[which.blanks, ] <- 1e+99
    else order.frame[which.blanks, ] <- -1e+99
    order.frame
  })))
  order.frame$decreasing <- decreasing
  order.frame$na.last <- NA
  retval <- do.call("order", order.frame)
  return(retval)
}

.index_column <- function(data, idx, default, trace) {
  col = colnames(data)[idx]
  if(is.null(col)) {
    if(trace) {
      cat(crayon::blue(paste0("There is no column name at ", deparse(substitute(data)), "[", idx, "]. Renaming this column as ", default, ".")))
    }
    return(default)
  }
  return(col)
} 

.guess_columns <- function(data,
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
                           trace) {
  if(!is.null(data)){
    if(any(c(data_peptide_column_PTMs,
             data_peptide_column_no_PTMs,
             data_accession_column) == "guess") && trace) {
      message("data:")
    }
    if(!is.null(data_peptide_column_PTMs) && data_peptide_column_PTMs == "guess") {
      data_peptide_column_PTMs <- .guess_modified_peptide_column(data, trace)
    } else if (is.numeric(data_peptide_column_PTMs)) {
      tmp = .index_column(data, data_peptide_column_PTMs, "Peptides (PTMs)", trace)
      colnames(data)[data_peptide_column_PTMs] = tmp
      data_peptide_column_PTMs = tmp
    }
    if (!is.null(data_peptide_column_no_PTMs) && data_peptide_column_no_PTMs == "guess") {
      data_peptide_column_no_PTMs = .guess_unmodified_peptide_column(data, trace)
      if(is.list(data_peptide_column_no_PTMs)) {
        data = data_peptide_column_no_PTMs[[1]]
        data_peptide_column_no_PTMs = data_peptide_column_no_PTMs[[2]]
      }
    } else if (is.numeric(data_peptide_column_no_PTMs)) {
      tmp = .index_column(data, data_peptide_column_no_PTMs, "Peptides (no PTMs)", trace)
      colnames(data)[data_peptide_column_no_PTMs] = tmp
      data_peptide_column_no_PTMs = tmp
    }
    if(is.null(data_peptide_column_PTMs) & is.null(data_peptide_column_no_PTMs)) {
      stop("Both variable `data_peptide_column_PTMs` and `data_peptide_column_no_PTMs` cannot be NULL")
    }
    if (is.null(data_peptide_column_PTMs) && !is.null(data_peptide_column_no_PTMs)) {
      data_peptide_column_PTMs <- paste(data_peptide_column_no_PTMs, "(PTMs)")
      if(trace) {
        cat(crayon::blue(paste0("Copying and pasting sequences in column `", data_peptide_column_no_PTMs, "` to column `", data_peptide_column_PTMs, "`.\n")))
      }
      data[data_peptide_column_PTMs] = data[[data_peptide_column_no_PTMs]]
    }
    if (!is.null(data_peptide_column_no_PTMs) && is.null(data_peptide_column_no_PTMs)) {
      data_peptide_column_no_PTMs <- paste(data_peptide_column_PTMs, "(no PTMs)")
      if(trace) {
        cat(crayon::blue(paste0("Copying and pasting sequences in column `", data_peptide_column_PTMs, "` to column `", data_peptide_column_no_PTMs, "`.\n")))
      }
      data[,data_peptide_column_no_PTMs] = data[[data_peptide_column_PTMs]]
      if(trace) {
        cat(crayon::blue(paste0("Removing PTMs in column `", data_peptide_column_no_PTMs, "`...\n")))
      }
      convert <- unique(data[,data_peptide_column_no_PTMs,with=F])
      convert[,"seq"] <- .remove_modifications(convert[[data_peptide_column_no_PTMs]])
      data[,data_peptide_column_no_PTMs] <- convert[data, seq, on = data_peptide_column_no_PTMs]
      if(trace) {
        cat(crayon::blue("done.\n"))
      }
    }
    
    rows_with_PTMs = data[[data_peptide_column_no_PTMs]] %like% "[^ILKMFTWVRHANDCEQGPSY]"
    if(any(rows_with_PTMs)) {
      res = .convert_modified_peptide_column(data, data_peptide_column_no_PTMs)
      if(is.list(res)) {
        data = res[[1]]
        data_peptide_column_no_PTMs = res[[2]]
      }
    }
    
    if(is.null(data_accession_column)) {
      data_accession_column <- "Proteins"
    } else if(data_accession_column == "guess") {
      data_accession_column <- .guess_accession_column(data, trace)
    } else if (is.numeric(data_accession_column)) {
      tmp = .index_column(data, data_accession_column, "Proteins", trace)
      colnames(data)[data_accession_column] = tmp
      data_accession_column = tmp
    }
  }
  
  if(!is.null(pep2pro)) {
    if(any(c(pep2pro_peptide_column,
             pep2pro_accession_column) == "guess") && trace) {
      message("pep2pro:")
    }
    if(pep2pro_accession_column == "guess") {
      pep2pro_accession_column <- .guess_accession_column(pep2pro, trace)
      stopifnot(!is.null(pep2pro_accession_column))
    } else if (is.numeric(pep2pro_accession_column)) {
      tmp = .index_column(pep2pro, pep2pro_accession_column, "Proteins", trace)
      colnames(pep2pro)[pep2pro_accession_column] = tmp
      pep2pro_accession_column = tmp
    }
    if(any(tolower(annotate_by_peptide) %in% c("unmodified", "um"))){
      if(pep2pro_peptide_column == "guess") {
        pep2pro_peptide_column <- .guess_unmodified_peptide_column(pep2pro, trace)
        if(is.list(pep2pro_peptide_column)) {
          pep2pro = pep2pro_peptide_column[[1]]
          pep2pro_peptide_column = pep2pro_peptide_column[[2]]
        }
        stopifnot(!is.null(pep2pro_peptide_column))
      } else if (is.numeric(pep2pro_peptide_column)) {
        tmp = .index_column(pep2pro, pep2pro_peptide_column, "Peptides (no PTMs)", trace)
        colnames(pep2pro)[pep2pro_peptide_column] = tmp
        pep2pro_peptide_column = tmp
      }
      rows_with_PTMs = pep2pro[[pep2pro_peptide_column]] %like% "[^ILKMFTWVRHANDCEQGPSY]"
      if(any(rows_with_PTMs)) {
        res = .convert_modified_peptide_column(pep2pro, pep2pro_peptide_column)
        if(is.list(res)) {
          pep2pro = res[[1]]
          pep2pro_peptide_column = res[[2]]
        }
      }
    } else if(pep2pro_peptide_column == "guess") {
      pep2pro_peptide_column <- .guess_modified_peptide_column(pep2pro, trace)
    } else if (is.numeric(pep2pro_peptide_column)) {
      tmp = .index_column(pep2pro, pep2pro_peptide_column, "Peptides (PTMs)", trace)
      colnames(pep2pro)[pep2pro_peptide_column] = tmp
      pep2pro_peptide_column = tmp
    }
    rows_with_aa = pep2pro[[pep2pro_peptide_column]] %like% "[ILKMFTWVRHANDCEQGPSY]"
    if(!all(rows_with_aa)) {
      stop(paste0("Column `", pep2pro_peptide_column, "` in `pep2pro` does not seem to be a peptide column because some rows do not contain a single amino acid character."))
    }
  }
  if(!is.null(pep2taxon)) {
    if(any(c(pep2taxon_peptide_column,
             rank_columns) == "guess") && trace) {
      message("pep2taxon:")
    }
    if(any(tolower(annotate_by_peptide) %in% c("unmodified", "um"))){
      if(pep2taxon_peptide_column == "guess") {
        pep2taxon_peptide_column <- .guess_unmodified_peptide_column(pep2taxon, trace)
        if(is.list(pep2taxon_peptide_column)) {
          pep2taxon = pep2taxon_peptide_column[[1]]
          pep2taxon_peptide_column = pep2taxon_peptide_column[[2]]
        }
      } else if (is.numeric(pep2taxon_peptide_column)) {
        tmp = .index_column(pep2taxon, pep2taxon_peptide_column, "Peptides (no PTMs)", trace)
        colnames(pep2taxon)[pep2taxon_peptide_column] = tmp
        pep2taxon_peptide_column = tmp
      }
      rows_with_PTMs = pep2taxon[[pep2taxon_peptide_column]] %like% "[^ILKMFTWVRHANDCEQGPSY]"
      if(any(rows_with_PTMs)) {
        res = .convert_modified_peptide_column(pep2taxon, pep2taxon_peptide_column)
        if(is.list(res)) {
          pep2taxon = res[[1]]
          pep2taxon_peptide_column = res[[2]]
        }
      }
    } else if(pep2taxon_peptide_column == "guess") {
      pep2taxon_peptide_column <- .guess_modified_peptide_column(pep2taxon, trace)
    } else if (is.numeric(pep2taxon_peptide_column)) {
      tmp = .index_column(pep2taxon, pep2taxon_peptide_column, "Peptides (PTMs)", trace)
      colnames(pep2taxon)[pep2taxon_peptide_column] = tmp
      pep2taxon_peptide_column = tmp
    }
    rows_with_aa = pep2taxon[[pep2taxon_peptide_column]] %like% "[ILKMFTWVRHANDCEQGPSY]"
    if(!all(rows_with_aa)) {
      stop(paste0("Column `", pep2taxon_peptide_column, "` in `pep2taxon` does not seem to be a peptide column because some rows do not contain a single amino acid character."))
    }
    if((length(rank_columns) == 1) && (rank_columns == "guess")) {
      rank_columns = get_ranks(pep2taxon)
    }
    pep2taxon_columns = colnames(pep2taxon)
    lca_column = pep2taxon_columns[tolower(pep2taxon_columns) %like% "lca"]
    lca_rank_column = pep2taxon_columns[tolower(pep2taxon_columns) %like% "rank"]
    pep2taxon[pep2taxon == ""] <- NA
    if(length(lca_rank_column) == 0) {
      pep2taxon <- pep2taxon[, rank := colnames(.SD)[max(which(!is.na(unlist(.SD))))], by = 1:nrow(pep2taxon), .SDcols = ranks]
      lca_rank_column = "rank"
    }
    if(length(lca_column) == 0) {
      pep2taxon <- pep2taxon[, lca := .SD[, max(which(!is.na(unlist(.SD)))), with = F], by = 1:nrow(pep2taxon), .SDcols = ranks]
      lca_column = "lca"
    }
    pep2taxon_columns <- c(lca_column, lca_rank_column, rank_columns)
    
  }
  
  if(!is.null(pro2func)) {
    if(any(c(pro2func_accession_column, pro2func_function_columns) == "guess") && trace) {
      message("pro2func:")
    }
    if((length(pro2func_function_columns) == 1) && 
       (pro2func_function_columns == "guess")) {
      pro2func_function_columns <- colnames(pro2func)[colnames(pro2func) %like% "COG|NOG|KEGG|GO|BRITE|REACTOME"]
      if(trace) {
        cat(crayon::silver(paste0("Using c(", paste0('"', pro2func_function_columns, '"', collapse = ", "), ") as the funtional annotation columns.\n")))
      }
    }
    if(!is.null(pro2func_accession_column) && (pro2func_accession_column == "guess")) {
      pro2func_accession_column <- .guess_accession_column(pro2func, trace)
    } else if (is.numeric(pro2func_accession_column)) {
      tmp = .index_column(pro2func, pro2func_accession_column, "Proteins", trace)
      colnames(pro2func)[pro2func_accession_column] = tmp
      pro2func_accession_column = tmp
    }
  }
  list(data = data,
       data_peptide_column_PTMs = data_peptide_column_PTMs,
       data_peptide_column_no_PTMs = data_peptide_column_no_PTMs,
       data_accession_column = data_accession_column,
       pep2pro = pep2pro,
       pep2pro_peptide_column = pep2pro_peptide_column,
       pep2pro_accession_column = pep2pro_accession_column,
       pep2taxon = pep2taxon,
       pep2taxon_peptide_column = pep2taxon_peptide_column,
       rank_columns = rank_columns,
       pep2taxon_columns = pep2taxon_columns,
       pro2func = pro2func,
       pro2func_accession_column = pro2func_accession_column,
       pro2func_function_columns = pro2func_function_columns)
  
}
