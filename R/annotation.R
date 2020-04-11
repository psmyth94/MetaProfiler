.accession_razor <- function(data, protein, peptide, accession_delimiter, progress) {
  x <- data.table::copy(data)
  x <- x[, unique := as.numeric(!grepl(accession_delimiter, get(protein))) , by = protein]
  x[[protein]] <- stringi::stri_split_fixed(x[[protein]], accession_delimiter)
  x <- x[, lapply(.SD, unlist), by=1:nrow(x)][,-1]
  x <- x[, c("count", "unique") := .(.N, sum(unique)), by = protein]
  x <- x[order(count, unique, decreasing = T)]
  x <- x[, .(paste0(get(protein), collapse = ","), get(protein)[[1]], paste0(count, collapse = ","), paste0(unique, collapse = ",")), by = peptide]
  x <- data.table::setnames(x, c(peptide, protein, paste(protein, "(razor)"), "count", "unique"))
  x
  # groups <- data.table::as.data.table(razor(x[[peptide]], x[[protein]], x[["unique"]]))
  # groups <- groups[, c(lapply(.SD, function(.sd) sapply(.sd, "[", 1)), list(x = x)), by = 1:nrow(groups), .SDcols = c("id", "count", "unique")]
  # groups <- groups[, .(x = unlist(x)), by = .(id, count, unique)]
  # groups <- groups[, order()]
  # groups[,peptide] <- x[[peptide]]
  # groups$Proteins
}
make_annotation_table <- function(Object,
                                  
                                  pep2pro = NULL,
                                  pep2pro_peptide_column = "guess",
                                  pep2pro_accession_column = "guess",
                                  
                                  compute_razor_protein = F,
                                  accession_delimiter = "[;\t,]",
                                  function_delimiter = accession_delimiter,
                                  accession_pattern = "[-_]",
                                  
                                  pep2taxon = NULL,
                                  pep2taxon_peptide_column = "guess",
                                  rank_columns = "guess",
                                  pro2func = NULL,
                                  pro2func_accession_column = "guess",
                                  function_columns = "guess",
                                  trace = T,
                                  progress = T) {
  
  if(!any(colnames(Object@data) == Object@time_unit)) {
    stop(paste0("No column named `", Object@time_unit, "` in slot `data`. Please add/rename the column for the timepoints or change slot `time_unit` so that it matches the timepoint column in slot `data`."))
  }
  
  if(is.character(pep2pro)) {
    stopifnot(file.exists(pep2pro))
    pep2pro <- data.table::fread(pep2pro)
  }
  
  if(is.character(pep2taxon)) {
    stopifnot(file.exists(pep2taxon))
    pep2taxon <- data.table::fread(pep2taxon)
  }
  
  if(is.character(pro2func)) {
    stopifnot(file.exists(pro2func))
    pro2func <- data.table::fread(pro2func)
  }
  
  args = .guess_columns(Object@data,
                       Object@peptide_column_PTMs,
                       Object@peptide_column_no_PTMs,
                       Object@accession_column,
                       pep2pro,
                       pep2pro_peptide_column,
                       pep2pro_accession_column,
                       pep2taxon,
                       pep2taxon_peptide_column,
                       rank_columns,
                       pro2func,
                       pro2func_accession_column,
                       function_columns,
                       Object@annotate_with,
                       accession_pattern,
                       trace)
  
  Object@data = args$data
  Object@peptide_column_PTMs = args$data_peptide_column_PTMs
  Object@peptide_column_no_PTMs = args$data_peptide_column_no_PTMs
  Object@accession_column = args$data_accession_column
  pep2pro = args$pep2pro
  pep2pro_peptide_column = args$pep2pro_peptide_column
  pep2pro_accession_column = args$pep2pro_accession_column
  pep2taxon = args$pep2taxon
  pep2taxon_peptide_column = args$pep2taxon_peptide_column
  rank_columns = args$rank_columns
  if(length(args$pep2taxon_columns))
    Object@taxon_columns = args$pep2taxon_columns
  pro2func = args$pro2func
  pro2func_accession_column = args$pro2func_accession_column
  if(length(args$function_columns))
    Object@function_columns = args$function_columns
  if(!is.null(pep2pro)) {
    if(compute_razor_protein) {
      if(trace) {
        cat(crayon::blue("Computing razor protein."))
      }
      delimiter <- table(unlist(stringi::stri_extract_all_regex(pep2pro[[pep2pro_accession_column]], accession_delimiter)))
      delimiter <- names(delimiter)[which.max(delimiter)]
      if(!length(delimiter)) delimiter = ";"
      pep2pro = data.table::setnames(pep2pro[,paste0(unique(get(pep2pro_accession_column)), collapse = delimiter), by = pep2pro_peptide_column], c(pep2pro_peptide_column,pep2pro_accession_column))
      pep2pro = .accession_razor(pep2pro, pep2pro_accession_column, pep2pro_peptide_column, delimiter, progress)      
    }
    peptides = unique(pep2pro[Object@data[[Object@annotate_with]], , on = pep2pro_peptide_column] )
    annotation_table = data.table::data.table(Proteins = peptides[[pep2pro_accession_column]], Peptides = peptides[[pep2pro_peptide_column]])
    annotation_table <- na.omit(annotation_table, "Proteins")
  } else {
    peptides = unique(Object@data[, c(Object@accession_column, Object@annotate_with), with = F])
    annotation_table = data.table::data.table(Proteins = peptides[[Object@accession_column]], Peptides = peptides[[Object@annotate_with]])
    if(compute_razor_protein) {
      if(progress) {
        writeLines("Computing razor protein...")
      }
      annotation_table$Proteins = .accession_razor(annotation_table, accession_delimiter, progress)      
    }
    annotation_table <- na.omit(annotation_table, "Proteins")
  }
  if(!nrow(annotation_table)) {
    stop("It looks like the accession column is empty. Please provide the Peptides file that contains the assigned proteins or specify the accession column in the result file.")
  }
  if(!length(Object@accession_column)) {
    Object@data[,Object@accession_column] <- annotation_table[Object@data[[Object@annotate_with]], Proteins, on = "Peptides"]
  } else {
    Object@data[,"Proteins"] <- annotation_table[Object@data[[Object@annotate_with]], Proteins, on = "Peptides"]
    Object@accession_column = "Proteins"
  }
  if(!is.null(pep2taxon)) {
    tax <- pep2taxon[annotation_table$Peptides, c(Object@taxon_columns), with = F, on = c(pep2taxon_peptide_column)]
    add <- rowSums(!is.na(tax)) > 0
    annotation_table[add,Object@taxon_columns] <- tax[add,]
  }
  if(!is.null(pro2func)) {
    pro2func[pro2func == ""] <- NA
    if(!length(Object@function_columns)) {
      stop("Names for the functional annotation columns are empty. Please specify the columns.")
    }
    delimiter <- table(unlist(stringi::stri_extract_all_regex(pro2func[[Object@function_columns[[1]]]], function_delimiter)))
    delimiter <- names(delimiter)[which.max(delimiter)]
    if(!length(delimiter)) delimiter = ";"
    pro2func = data.table::setnames(pro2func[,lapply(.SD, function(x) paste0(unique(x), collapse = accession_delimiter)), by = pro2func_accession_column, .SDcols = Object@function_columns], c(pro2func_accession_column,Object@function_columns))
    annot <- pro2func[annotation_table$Proteins, c(Object@function_columns), with = F, on = c(pro2func_accession_column)]
    add <- rowSums(!is.na(annot)) > 0
    annotation_table[add,Object@function_columns] <- annot[add,]
  }
  Object@annotation_table = annotation_table
  Object
}

getTaxon <- function(Object) {
  taxon_columns <- Object@taxon_columns
  Object@master_tbl[,taxon_columns] <- Object@annotation_table[Object@master_tbl[[Object@annotate_with]], ..taxon_columns, on = "Peptides"]
  Object
}

getFunc <- function(Object) {
  function_columns <- Object@function_columns
  Object@master_tbl[,function_columns] <- Object@annotation_table[Object@master_tbl[[Object@annotate_with]], ..function_columns, on = "Peptides"]
  Object
}




