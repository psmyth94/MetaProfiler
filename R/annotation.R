.accession_razor <- function(annotation_table, accession_delimiter, progress) {
  x <- copy(annotation_table)
  x <- x[, unique := as.numeric(!grepl(accession_delimiter, Proteins)) , by = Proteins]
  x$Proteins <- stringi::stri_split_fixed(x$Proteins, accession_delimiter)
  x <- x[, lapply(.SD, unlist), by=1:nrow(x)]
  y <- x[, .(Peptides = list(Peptides),
             freq = length(Peptides),
             unique = sum(unique)), by = Proteins]
  y <- y[with(y, order(freq, unique, decreasing = T))]
  
  groups <- as.data.table(razor(y$Peptides, y$Proteins, progress))
  groups$id <- sapply(groups$id, "[", 1)
  groups$freq <- lengths(groups$x)
  groups <- groups[, lapply(.SD, unlist), by=1:nrow(groups), .SDcols = c("x", "id", "freq")]
  groups$unique <- y[groups$id, unique, on = "Proteins"]
  groups <- groups[with(groups, order(freq, unique, decreasing = T))]
  groups <- groups[, .(Proteins = list(id), unique = list(unique), freq = list(freq)), by = x]
  groups <- groups[annotation_table$Peptides, , on = "x"]
  groups$Proteins
}

make_annotation_table <- function(time_unit,
                                  
                                  data,
                                  data_peptide_column_PTMs = "guess",
                                  data_peptide_column_no_PTMs = "guess",
                                  data_accession_column = "guess",
                                  
                                  annotate_by_peptide = c("unmodified", "um", "modified", "m"),
                                  pep2pro = NULL,
                                  pep2pro_peptide_column = "guess",
                                  pep2pro_accession_column = "guess",
                                  
                                  compute_razor_protein = F,
                                  accession_delimiter = ";",
                                  
                                  pep2taxon = NULL,
                                  pep2taxon_peptide_column = "guess",
                                  rank_columns = "guess",
                                  pro2func = NULL,
                                  pro2func_accession_column = "guess",
                                  pro2func_function_columns = "guess",
                                  trace = T,
                                  progress = T) {
  
  annotate_by_peptide <- match.arg(annotate_by_peptide)
  data <- as.data.table(data)
  if(!any(colnames(data) == time_unit)) {
    stop(paste0("No column named `", time_unit, "` in the table `data`. Please add/rename the column for the timepoints or change variable `time_unit` so that it matches the timepoint column in table `data`."))
  }
  
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
  
  if(!is.null(pep2pro)) {
    which_peptide = ifelse(any(tolower(annotate_by_peptide) %in% c("unmodified", "um")),
                           data_peptide_column_no_PTMs,
                           data_peptide_column_PTMs)
    peptides = unique(pep2pro[data[[which_peptide]], , on = pep2pro_peptide_column] )
    annotation_table = data.table(Proteins = peptides[[pep2pro_accession_column]], Peptides = peptides[[pep2pro_peptide_column]])
    if(compute_razor_protein) {
      if(progress) {
        writeLines("Computing razor protein...")
      }
      annotation_table$Proteins = .accession_razor(annotation_table, accession_delimiter, progress)      
    }
    annotation_table <- na.omit(annotation_table, "Proteins")
  } else {
    which_peptide = ifelse(any(tolower(annotate_by_peptide) %in% c("unmodified", "um")),
                           data_peptide_column_no_PTMs,
                           data_peptide_column_PTMs)
    peptides = unique(data[, c(data_accession_column, which_peptide), with = F])
    annotation_table = data.table(Proteins = peptides[[data_accession_column]], Peptides = peptides[[which_peptide]])
    if(compute_razor_protein) {
      if(progress) {
        writeLines("Computing razor protein...")
      }
      annotation_table$Proteins = .accession_razor(annotation_table, accession_delimiter, progress)      
    }
    annotation_table <- na.omit(annotation_table, "Proteins")
  }
  if(all(is.na(annotation_table[[data_accession_column]]))) {
    stop("It looks like the accession column is empty. Please provide the Peptides file that contains the assigned proteins or specify the accession column in the result file.")
  }
  
  if(!is.null(pep2taxon)) {
    tax <- pep2taxon[annotation_table$Peptides, pep2taxon_columns, with = F, on = c(pep2taxon_peptide_column)]
    add <- rowSums(!is.na(tax)) > 0
    annotation_table[add,pep2taxon_columns] <- tax[add,]
  }
  if(!is.null(pro2func)) {
    pro2func[pro2func == ""] <- NA
    annot <- pro2func[annotation_table$Proteins, pro2func_function_columns, with = F, on = c(pro2func_accession_column)]
    add <- rowSums(!is.na(annot)) > 0
    annotation_table[add,pro2func_function_columns] <- annot[add,]
  }
  annotation_table
}

getTaxon <- function(Object) {
  taxon_columns <- Object@taxon_columns
  Object@master_tbl[,taxon_columns] <- Object@annotation_table[Object@master_tbl[[Object@annotate_by_peptide]], ..taxon_columns, on = "Peptides"]
  Object
}

getFunc <- function(Object) {
  function_columns <- Object@function_columns
  Object@master_tbl[,function_columns] <- Object@annotation_table[Object@master_tbl[[Object@annotate_by_peptide]], ..function_columns, on = "Peptides"]
  Object
}




