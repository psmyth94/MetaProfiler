create_experimental_design <- function(filenames, results_directory, ...) {
  if(missing(filenames)) {
    filenames <- dir(result_directory, pattern = c("tsv|TSV|txt|TXT|csv|CSV"), full.names = T)
  }
  args <- list(...)
  d <- data.table::data.table(filenames = filenames, data.table::as.data.table(lapply(args, function(x) {
    stringi::stri_extract_last_regex(filenames, x)
  })))
  shush(d <- d[,lapply(.SD, function(x) {
    x_num <- as.numeric(x)
    if(all(is.na(x_num))) {
      return(x)
    }
    x_num
  })])
  d
}

get_cols <- function(var, x = Object@timepoints, time_unit = Object@time_unit) {
  nt = length(x)
  nv = length(var)
  paste0(rep(var, each = nt), " (", time_unit, " ", rep(x, nv), ")")
}

plural <- function (singular, plural, n = NULL)
{
  return (ifelse(n>1, plural, singular))
}

combine_words <- function (words, sep = ", ", and = " and ", before = "`", after = before) 
{
  n = length(words)
  rs = xfun::raw_string
  if (n == 0) 
    return(words)
  words = paste0(before, words, after)
  if (n == 1) 
    return(rs(words))
  if (n == 2) 
    return(rs(paste(words, collapse = and)))
  if (grepl("^ ", and) && grepl(" $", sep)) 
    and = gsub("^ ", "", and)
  words[n] = paste0(and, words[n])
  rs(paste(words, collapse = sep))
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

rescale <- function(x, to = c(0, 1), from = range(x, na.rm = TRUE, finite = TRUE)) {
  (x - from[1]) / diff(from) * diff(to) + to[1]
}

save_plot <- function (filename, plot = last_plot(), device = NULL, path = NULL, 
                       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
                       dpi = 300, limitsize = TRUE, ...) {
  dpi <- ggplot2:::parse_dpi(dpi)
  dev <- ggplot2:::plot_dev(device, filename, dpi = dpi)
  dim <- ggplot2:::plot_dim(c(width, height), scale = scale, units = units, 
                            limitsize = limitsize)
  if (length(path)) {
    filename <- file.path(path, filename)
  }
  old_dev <- grDevices::dev.cur()
  dev(filename = filename, width = dim[1], height = dim[2])
  on.exit(utils::capture.output({
    grDevices::dev.off()
    if (old_dev > 1) grDevices::dev.set(old_dev)
  }))
  if(inherits(plot, c("Heatmap", "AnnotationFunction", "Legends", "HeatmapAnnotation", "SingleAnnotation", "HeatmapList"))) {
    ComplexHeatmap::draw(plot, ...)
  } else if (inherits(plot, "igraph")) {
    plot(plot, ...)
  } else if (inherits(plot, "gtable")) {
    gridExtra::grid.arrange(plot, ...)
  } else {
    grid::grid.draw(plot)
  }
  invisible()
}

.guess_accession_column <- function(d, trace = T, pattern = "[-_]") {
  cols <- toupper(colnames(d))
  suppressWarnings(idx <- which(grepl("ACC|PRO|QUERY|NAME|GENE|LEAD|RAZOR|GROUP", cols, perl = T)))
  if(length(idx) == 0) {
    message(paste0("Could not guess the name of the accession/protein column for `", deparse(substitute(d)), "`. Please provide the name of this column."))
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
    message(paste0("Could not guess the name of the accession/protein column for `", deparse(substitute(d)), "`. Please provide the name of this column."))
    return(NULL)
  }
}

.convert_modified_peptide_column <- function(d, col, trace = T) {
  if(trace) {
    cat(crayon::blue(paste0("Looks like `", col, "` is a peptide column, but appears to contain post translational modifications.\n")))
  }
  new_column = paste(col, "(no PTMs)")
  if(trace) {
    cat(crayon::blue(paste0("Creating column `", new_column,"` with sequences from column `", col, "`.\n")))
  }
  d[,new_column] = d[[col]]
  if(trace) {
    cat(crayon::blue(paste0("Removing PTMs...")))
  }
  convert <- unique(d[,new_column,with=F])
  convert[,"seq"] <- .remove_modifications(convert[[new_column]])
  
  if(all(!grepl("[^ILKMFTWVRHANDCEQGPSY]", convert$seq, perl = T)))
  {
    d[,new_column] <- convert[d, seq, on = new_column]
  } else {
    cat(crayon::black(paste0("Failed to remove PTMs, returning NULL.\n")))
    return(NULL)
  }
  if(trace) {
    cat(crayon::blue("done.\n"))
  }
  list(d, new_column)
}

.guess_unmodified_peptide_column <- function(d, trace = T) {
  cols <- toupper(colnames(d))
  suppressWarnings(idx <- which(grepl("PEP|SEQ", cols)))
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
    idx <- idx[which.max(lookahead == nrow(d))]
    if(trace) {
      cat(crayon::silver(paste0("Using `", colnames(d)[idx], "` as the unmodified peptide column.\n")))
    }
    return(colnames(d)[idx])
  } else {
    idx2 <- idx[lookahead > (nrow(d) * 0.5)]
    if(length(idx2) > 0) {
      idx <- idx[which.max(lookahead > (nrow(d) * 0.5))]
      res <- .convert_modified_peptide_column(d, colnames(d)[idx], trace)
      cat(crayon::silver(paste0("Using `", res[[2]], "` as the unmodified peptide column.\n")))
      return(res)
    }
    warning(paste0("Could not guess the name of the unmodified peptide column for `", deparse(substitute(d)), "`. Please provide the name of this column."))
    return(NULL)
  }
}

.guess_modified_peptide_column <- function(d, trace = T) {
  cols <- toupper(colnames(d))
  suppressWarnings(idx <- which(grepl("PEP|SEQ", cols)))
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
      cat(crayon::silver(paste0("Using `", colnames(d)[idx][1], "` as the modified peptide column.\n")))
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

guess_measurements <- function(x, cols) {
  if(length(x) == 0 || is.na(x)) return(NULL)
  cols <- .get_columns(x, cols)
  if(length(unique(cols)) == 0)
  {
    warning(paste0("No valid columns that contains the variable `", x,"`. Please add columns where the first word is the variable's name, followed by a unique identifier (e.g. ",
                   paste0(x, " ", 1:3, ", ", collapse = ""), "..., or ", paste0(x, " ", c("light", "heavy"), ", ", collapse = ""),"). Returning NULL."))
    return(NULL)
  }
  cols
}

.check_columns <- function(x, cols) {
  if(!length(x) || is.na(x)) return(character())
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
  stringi::stri_replace_all_regex(stringi::stri_replace_all_regex(x, "\\([\\S\\s]+?\\)|\\[[\\S\\s]+?\\]|\\{[\\S\\s]+?\\}|\\<[\\S\\s]+?\\>|_[\\S\\s]+?_", ""), "[^ILKMFTWVRHANDCEQGPSY]", "")
}

.check_if_percantage_or_fraction <- function(x) {
  shush(v <- max(x, na.rm = T))
  if(v > 1 && v <= 100) {
    return(T)
  } else if(v > 0 && v <= 1) {
    return(F)
  } else {
    return(NULL)
  }
}

.index_column <- function(data, idx, default, trace) {
  col = colnames(data)[idx]
  if(!length(col)) {
    if(trace) {
      cat(crayon::blue(paste0("There is no column name at ", deparse(substitute(data)), "[", idx, "]. Renaming this column as ", default, ".")))
    }
    return(default)
  }
  return(col)
}



.guess_columns <- function(data,
                           data_peptide_column_PTMs = "guess",
                           data_peptide_column_no_PTMs = "guess",
                           data_accession_column = "guess",
                           pep2pro = NULL,
                           pep2pro_peptide_column = "guess",
                           pep2pro_accession_column = "guess",
                           pep2taxon = NULL,
                           pep2taxon_peptide_column = "guess",
                           rank_columns = "guess",
                           pro2func = NULL,
                           pro2func_accession_column = "guess",
                           function_columns = "guess",
                           annotate_with = NULL,
                           accession_pattern = "[-_]",
                           trace = T) {
  if(length(data)){
    if(any(c(data_peptide_column_PTMs,
             data_peptide_column_no_PTMs,
             data_accession_column) == "guess") && trace) {
      message("data:")
    }
    if(length(data_peptide_column_PTMs) && data_peptide_column_PTMs == "guess") {
      data_peptide_column_PTMs <- .guess_modified_peptide_column(data, trace)
    } else if (is.numeric(data_peptide_column_PTMs)) {
      tmp = .index_column(data, data_peptide_column_PTMs, "Peptides (PTMs)", trace)
      colnames(data)[data_peptide_column_PTMs] = tmp
      data_peptide_column_PTMs = tmp
    }
    if (length(data_peptide_column_no_PTMs) && data_peptide_column_no_PTMs == "guess") {
      data_peptide_column_no_PTMs = .guess_unmodified_peptide_column(data, T)
      if(is.list(data_peptide_column_no_PTMs)) {
        data = data_peptide_column_no_PTMs[[1]]
        data_peptide_column_no_PTMs = data_peptide_column_no_PTMs[[2]]
      }
    } else if (is.numeric(data_peptide_column_no_PTMs)) {
      tmp = .index_column(data, data_peptide_column_no_PTMs, "Peptides (no PTMs)", trace)
      colnames(data)[data_peptide_column_no_PTMs] = tmp
      data_peptide_column_no_PTMs = tmp
    }
    if(!length(data_peptide_column_PTMs) && !length(data_peptide_column_no_PTMs)) {
      stop("Both variable `data_peptide_column_PTMs` and `data_peptide_column_no_PTMs` cannot be NULL")
    }
    if (!length(data_peptide_column_PTMs) && length(data_peptide_column_no_PTMs)) {
      data_peptide_column_PTMs <- paste(data_peptide_column_no_PTMs, "(PTMs)")
      if(trace) {
        cat(crayon::blue(paste0("Copying and pasting sequences in column `", data_peptide_column_no_PTMs, "` to column `", data_peptide_column_PTMs, "`.\n")))
      }
      data[data_peptide_column_PTMs] = data[[data_peptide_column_no_PTMs]]
    }
    if (length(data_peptide_column_no_PTMs) && !length(data_peptide_column_no_PTMs)) {
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
    
    rows_with_PTMs = grepl("[^ILKMFTWVRHANDCEQGPSY]", data[[data_peptide_column_no_PTMs]])
    if(any(rows_with_PTMs)) {
      res = .convert_modified_peptide_column(data, data_peptide_column_no_PTMs)
      if(is.list(res)) {
        data = res[[1]]
        data_peptide_column_no_PTMs = res[[2]]
      }
    }
    if(length(data_accession_column) && data_accession_column == "guess") {
      data_accession_column <- .guess_accession_column(data, trace, accession_pattern)
    } else if (is.numeric(data_accession_column)) {
      tmp = .index_column(data, data_accession_column, "Proteins", trace)
      colnames(data)[data_accession_column] = tmp
      data_accession_column = tmp
    } 
    
    if(!length(data_accession_column)) {
      data_accession_column <- "Proteins"
    }
  }
  
  if(length(pep2pro)) {
    if(any(c(pep2pro_peptide_column,
             pep2pro_accession_column) == "guess") && trace) {
      message("pep2pro:")
    }
    if(pep2pro_accession_column == "guess") {
      pep2pro_accession_column <- .guess_accession_column(pep2pro, trace, accession_pattern)
      stopifnot(!is.null(pep2pro_accession_column))
    } else if (is.numeric(pep2pro_accession_column)) {
      tmp = .index_column(pep2pro, pep2pro_accession_column, "Proteins", trace)
      colnames(pep2pro)[pep2pro_accession_column] = tmp
      pep2pro_accession_column = tmp
    }
    if(annotate_with == data_peptide_column_no_PTMs){
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
      rows_with_PTMs = grepl("[^ILKMFTWVRHANDCEQGPSY]", pep2pro[[pep2pro_peptide_column]])
      if(any(rows_with_PTMs)) {
        res = .convert_modified_peptide_column(pep2pro, pep2pro_peptide_column)
        if(is.list(res)) {
          pep2pro = res[[1]]
          pep2pro_peptide_column = res[[2]]
        }
      }
    } else if(annotate_with == data_peptide_column_PTMs && pep2pro_peptide_column == "guess") {
      pep2pro_peptide_column <- .guess_modified_peptide_column(pep2pro, trace)
    } else if (is.numeric(pep2pro_peptide_column)) {
      tmp = .index_column(pep2pro, pep2pro_peptide_column, "Peptides (PTMs)", trace)
      colnames(pep2pro)[pep2pro_peptide_column] = tmp
      pep2pro_peptide_column = tmp
    }
    rows_with_aa = grepl("[ILKMFTWVRHANDCEQGPSY]", pep2pro[[pep2pro_peptide_column]])
    if(!all(rows_with_aa)) {
      stop(paste0("Row(s) [", ifelse(sum(!rows_with_aa) > 5, 
                                     paste0(paste0(which(!rows_with_aa)[1:5], collapse = ","), "..."),
                                     paste0(which(!rows_with_aa), collapse = ",")),
                  "] in column `", pep2pro_peptide_column, "` for `pep2pro` seem to contain characters that are not single letter amino acid codes"))
    }
  }
  pep2taxon_columns = NULL
  if(length(pep2taxon)) {
    if(any(c(pep2taxon_peptide_column,
             rank_columns) == "guess") && trace) {
      message("pep2taxon:")
    }
    if(annotate_with == data_peptide_column_no_PTMs){
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
      rows_with_PTMs = grepl("[^ILKMFTWVRHANDCEQGPSY]", pep2taxon[[pep2taxon_peptide_column]])
      if(any(rows_with_PTMs)) {
        res = .convert_modified_peptide_column(pep2taxon, pep2taxon_peptide_column)
        if(is.list(res)) {
          pep2taxon = res[[1]]
          pep2taxon_peptide_column = res[[2]]
        }
      }
      rows_with_aa = grepl("[ILKMFTWVRHANDCEQGPSY]", pep2taxon[[pep2taxon_peptide_column]])
      if(!all(rows_with_aa)) {
        stop(paste0("Row(s) [", ifelse(sum(!rows_with_aa) > 5, 
                                       paste0(paste0(which(!rows_with_aa)[1:5], collapse = ","), "..."),
                                       paste0(which(!rows_with_aa), collapse = ",")),"] in column `", pep2taxon_peptide_column, "` for `pep2taxon` seem to contain characters that are not single letter amino acid codes"))
      }
    } else if(annotate_with == data_peptide_column_PTMs && pep2taxon_peptide_column == "guess") {
      pep2taxon_peptide_column <- .guess_modified_peptide_column(pep2taxon, trace)
    } else if (is.numeric(pep2taxon_peptide_column)) {
      tmp = .index_column(pep2taxon, pep2taxon_peptide_column, "Peptides (PTMs)", trace)
      colnames(pep2taxon)[pep2taxon_peptide_column] = tmp
      pep2taxon_peptide_column = tmp
    }
    if((length(rank_columns) == 1) && (rank_columns == "guess")) {
      rank_columns = get_ranks(pep2taxon)
    }
    pep2taxon_columns = colnames(pep2taxon)
    lca_column = pep2taxon_columns[grepl("lca", tolower(pep2taxon_columns))]
    lca_rank_column = pep2taxon_columns[grepl("rank", tolower(pep2taxon_columns))]
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
  if(length(pro2func)) {
    if(any(c(pro2func_accession_column, function_columns) == "guess") && trace) {
      message("pro2func:")
    }
    if((length(function_columns) == 1) && 
       (function_columns == "guess")) {
      function_columns <- colnames(pro2func)[grepl("COG|NOG|KEGG|GO|BRITE|REACTOME|PRO[\\s\\S]+?NAME", colnames(pro2func), perl = T)]
      if(trace) {
        cat(crayon::silver(paste0("Using c(", paste0('"', function_columns, '"', collapse = ", "), ") as the functional annotation columns.\n")))
      }
    }
    if(length(pro2func_accession_column) && (pro2func_accession_column == "guess")) {
      pro2func_accession_column <- .guess_accession_column(pro2func, trace, accession_pattern)
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
       function_columns = function_columns)
  
}


store <- function(mz,RT,sequence,charge,aa_before,aa_after,start,end,ids,scan,
                  pro2id, output, db, Experiment, document_id, used_target_decoy = T, charges = "1+, 4+",
                  fixed_modifications = c("Carbamidomethyl (C)"),
                  variable_modifications = c("Acetyl (N-term)", "Oxidation (M)"),
                  taxnomomy = "", enzyme = "trypsin",
                  mass_type = "monoisotopic", missed_cleavages = 2, precursor_peak_tolerance = 0,
                  precursor_peak_tolerance_ppm = F, peak_mass_tolerance = 0, peak_mass_tolerance_ppm = F,
                  protein_score_type = "", protein_charges = "")
{
  require(dplyr)
  os <- paste(
    # IdXML header
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n",
    "<?xml-stylesheet type=\"text/xsl\" href=\"https://www.openms.de/xml-stylesheet/IdXML.xsl\" ?>\n",
    "<IdXML version=\"1.5\"",
    " id=\"", document_id, "\"",
    " xsi:noNamespaceSchemaLocation=\"https://www.openms.de/xml-schema/IdXML_1_5.xsd\" ",
    "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n",
    
    # IdXML search parameters.
    "\t<SearchParameters charges=\"", charges, "\"",
    " id=\"SP_0\"",
    " db=\"", paste0(db, collapse = " "), "\"",
    " db_version=\"0\"",
    " taxnomomy=\"", taxnomomy, "\"",
    " mass_type=\"", mass_type, "\"",
    " enzyme=\"", enzyme, "\"",
    " missed_cleavages=\"", missed_cleavages, "\"",
    " precursor_peak_tolerance=\"", precursor_peak_tolerance, "\"",
    " precursor_peak_tolerance_ppm=\"", ifelse(precursor_peak_tolerance_ppm, "true", "false"), "\"",
    " peak_mass_tolerance=\"", peak_mass_tolerance, "\"", "",
    " peak_mass_tolerance_ppm=\"", ifelse(peak_mass_tolerance_ppm, "true", "false"), "\" >\n",
    paste0("\t\t<FixedModification name=\"", fixed_modifications, "\" />", collapse = "\n"),
    "\n",
    paste0("\t\t<VariableModification name=\"", variable_modifications ,"\" />", collapse = "\n"),
    "\n",
    "\t\t\t<UserParam type=\"string\" name=\"TargetDecoyApproach\" value=\"",
    ifelse(used_target_decoy, "true", "false"), "\"/>\n",
    "\t</SearchParameters>\n",
    
    "\t<IdentificationRun",
    " date=\"0000-00-00T00:00:00\"",
    " search_engine=\"MetaLab\"",
    " search_engine_version=\"2.0.0\"",
    " search_parameters_ref=\"SP_0\"",
    " >\n",
    
    "\t\t<ProteinIdentification",
    " score_type=\"\"",
    " charges=\"0\"",
    " higher_score_better=\"false\"",
    " significance_threshold=\"0\" >\n",
    paste0("\t\t\t<ProteinHit ",
           "id=\"PH_" ,  pro2id$id,
           "\" accession=\"", pro2id$protein,
           "\" score=\"0\" ",
           "sequence=\"", pro2id$sequence, "\" >\n",
           "\t\t\t\t<UserParam type=\"string\" name=\"target_decoy\" value=\"",
           case_when(grepl("DECOY_", pro2id$protein) ~ "decoy", T ~ "target"), "\"/>\n",
           "\t\t\t\t<UserParam type=\"string\" name=\"Description\" value=\"", pro2id$description, "\"/>\n",
           "\t\t\t</ProteinHit>",
           collapse = "\n"),
    "\n\t\t</ProteinIdentification>\n",
    
    paste0("\t\t<PeptideIdentification ",
           "score_type=\"PEP\" ",
           "higher_score_better=\"false\" ",
           "significance_threshold=\"0\" ",
           "MZ=\"" ,  mz ,  "\" ",
           "RT=\"" ,  RT,  "\" >\n",
           "\t\t\t<PeptideHit",
           " score=\"", score, "\"",
           " sequence=\"" ,  sequence,  "\"",
           " charge=\"" , charge,  "\"",
           " aa_before=\"", aa_before, "\"",
           " aa_after=\"", aa_after, "\"",
           " start=\"", start, "\"",
           " end=\"", end, "\"",
           " protein_refs=\"", ids, "\"",
           " spectrum_reference=\"scan=",scan,"\"",
           " >\n",
           "\t\t\t</PeptideHit>\n",
           "\t\t</PeptideIdentification>",
           collapse = "\n"),
    
    "\n\t</IdentificationRun>\n",
    "</IdXML>\n", sep = ""
  )
  writeLines(os, con = file.path(output, paste(Experiment, ".idXML", sep = "")))
}

MQ2idXML <- function(MQ.output, db)
{
 
  evi = data.table::fread(file.path(MQ.output, "evidence.txt"))
  
  evi = evi[evi$`Potential contaminant` != "+" & evi$Proteins != "",]
  
  evi$`Modified sequence` <- gsub("\\(ac\\)", "(Acetyl)", evi$`Modified sequence`)
  evi$`Modified sequence` <- gsub("\\(ox\\)", "(Oxidation)", evi$`Modified sequence`)
  evi$`Modified sequence` <- gsub("_\\(", ".(", evi$`Modified sequence`)
  evi$`Modified sequence` <- gsub("_", "", evi$`Modified sequence`)
  evi$`Modified sequence` <- gsub("C", "C(Carbamidomethyl)", evi$`Modified sequence`)
  
  annotation_table <- unique(evi[, .(Sequence, Proteins)])
  annotation_table$Proteins <- stringi::stri_split_fixed(annotation_table$Proteins, ";")
  annotation_table <- annotation_table[, .(Proteins = unlist(Proteins)), by = Sequence]
  
  proteins <- unique(annotation_table$Proteins)
  fasta <- data.table::as.data.table(read_fasta(db, proteins))
  fasta$sequence <- gsub("\n", "", fasta$sequence)
  annotation_table$`Protein Sequence` <- fasta[annotation_table$Proteins,sequence, on = "accession"]
  annotation_table$`Protein Descriptions` <- fasta[annotation_table$Proteins,description, on = "accession"]
  annotation_table <- na.omit(annotation_table, "Protein Sequence")
  annotation_table$pos <- stringi::stri_locate_all_fixed(annotation_table$`Protein Sequence`, annotation_table$Sequence)
  
  annotation_table$aa <- stringi::stri_match_all_regex(
    annotation_table$`Protein Sequence`, paste0("(\\w?)", annotation_table$Sequence, "(\\w?)")
  )
  
  
  f <- function(pos, aa) {
    mat <- cbind(pos[[1]], aa[[1]])[,-3,drop=F]
    data.table::setnames(data.table::data.table(mat), c("start", "end", "aa_before", "aa_after"))
  }
  annotation_table <- annotation_table[, f(pos, aa), by = .(Sequence, Proteins, `Protein Sequence`, `Protein Descriptions`)]
  annotation_table$aa_before[annotation_table$aa_before == ""] <- "["
  annotation_table$aa_after[annotation_table$aa_after == ""] <- "]"
  annotation_table$start <- annotation_table$start - 1
  annotation_table$end <- annotation_table$end - 1
  col.names <- colnames(annotation_table)
  annotation_table <- annotation_table[, c(lapply(.SD[, col.names[2:4], with = F], list), lapply(.SD[, col.names[5:8], with = F], paste0, collapse = " ")), by = Sequence]
  setkey(annotation_table, Sequence)
  evi[,col.names] <- annotation_table[evi$Sequence,]
  
  f <- function(x)
  {
    if(is.character(x))
    {
      if(all(!grepl(";", x)))
      {
        return(x)
      }
      return(stringi::stri_split_fixed(x, ";"))
    }
    return(x)
  }
  evi <- evi[, lapply(.SD, f)]
  evi[is.na(evi)] <- 0
  for (e in unique(evi$Experiment))
  {
    x = evi[Experiment == e,]
    pro2id <- data.table(protein = unlist(x$Proteins), sequence = unlist(x$`Protein Sequences`),
                         description = unlist(x$`Protein Descriptions`))
    pro2id <- unique(pro2id)
    pro2id$id <- 1:nrow(pro2id)
    setkey(pro2id, protein)
    map2pep <- rep(1:nrow(x), lengths(x$Proteins))
    ids <- pro2id[unlist(x$Proteins), id]
    ids <- split(ids, map2pep)
    ids <- sapply(ids, function(x) paste0("PH_", x, collapse = " "))
    x$Protein_ids <- ids
    store(x, pro2id, db, e, e)
  }
}
Scaffold2idXML <- function(file, db)
{
  csv = data.table::fread(file)
  csv <- csv[-nrow(csv)]
  csv[csv == ""] <- NA
  csv$Sequence = Sequence = gsub("m", "M", csv$`Peptide sequence`)
  annotation_table <- unique(csv[, .(Sequence, `Protein accession numbers`)])
  csv$`Peptide sequence` <- gsub("m", "M(Oxidation)", csv$`Peptide sequence`)
  annotation_table <- na.omit(annotation_table)
  annotation_table$`Protein accession numbers` <- stringi::stri_split_fixed(annotation_table$`Protein accession numbers`, ",")
  annotation_table <- unique(annotation_table[, .(Proteins = unlist(`Protein accession numbers`)), by = Sequence])
  
  proteins <- unique(annotation_table$Proteins)
  fasta <- unique(data.table::as.data.table(read_fasta(db, proteins)))
  fasta$sequence <- gsub("\n", "", fasta$sequence)
  annotation_table$`Protein Sequence` <- fasta[annotation_table$Proteins,sequence, on="accession"]
  annotation_table$`Protein Descriptions` <- fasta[annotation_table$Proteins,description, on="accession"]
  annotation_table$pos <- stringi::stri_locate_all_fixed(annotation_table$`Protein Sequence`, annotation_table$Sequence)
  
  annotation_table$aa <- stringi::stri_match_all_regex(
    annotation_table$`Protein Sequence`, paste0("(\\w?)", annotation_table$Sequence, "(\\w?)")
  )
  
  f <- function(pos, aa) {
    mat <- cbind(pos[[1]], aa[[1]])[,-3,drop=F]
    data.table::setnames(data.table::data.table(mat), c("start", "end", "aa_before", "aa_after"))
  }
  annotation_table <- annotation_table[, f(pos, aa), by = .(Sequence, Proteins, `Protein Sequence`, `Protein Descriptions`)]
  annotation_table <- na.omit(annotation_table, c("Protein Sequence"))
  annotation_table$aa_before[annotation_table$aa_before == ""] <- "["
  annotation_table$aa_after[annotation_table$aa_after == ""] <- "]"
  annotation_table$start <- as.numeric(annotation_table$start)
  annotation_table$end <- as.numeric(annotation_table$end)
  annotation_table$start <- annotation_table$start - 1
  annotation_table$end <- annotation_table$end - 1
  col.names <- colnames(annotation_table)
  annotation_table <- annotation_table[, c(lapply(.SD[, col.names[2:4], with = F], list), lapply(.SD[, col.names[5:8], with = F], paste0, collapse = " ")), by = Sequence]
  csv <- csv[csv$Sequence %in% annotation_table$Sequence]
  csv[,col.names] <- annotation_table[csv$Sequence,,on="Sequence"]
  csv$`Observed m/z` <- as.numeric(gsub(",", "", csv$`Observed m/z`))
  f <- function(x)
  {
    if(is.character(x))
    {
      if(all(!grepl(",", x)))
      {
        return(x)
      }
      return(stringi::stri_split_regex(x, ",\\s*"))
    }
    return(x)
  }
  csv <- csv[, lapply(.SD, f)]
  csv[is.na(csv)] <- 0
  for (e in unique(csv$`MS/MS sample name`))
  {
    x = csv[`MS/MS sample name` == e]
    pro2id <- data.table::data.table(protein = unlist(x$Proteins), sequence = unlist(x$`Protein Sequence`),
                         description = unlist(x$`Protein Descriptions`))
    pro2id <- unique(pro2id)
    pro2id$id <- 1:nrow(pro2id)
    map2pep <- rep(1:nrow(x), lengths(x$Proteins))
    ids <- pro2id[unlist(x$Proteins), id, on = "protein"]
    ids <- split(ids, map2pep)
    ids <- sapply(ids, function(x) paste0("PH_", x, collapse = " "))
    x$Protein_ids <- ids
    mz = x$`Observed m/z`
    # RT,sequence,charge,aa_before,aa_after,start,end,ids,scan
    # mz,RT,sequence,charge,aa_before,aa_after,start,end,ids,scan
    store(x, pro2id, "C:/Users/Patrick Smyth/Documents/workspace", db, e, e, fixed_modifications = NULL, variable_modifications = "Oxidation (M)")
    
  }
}

