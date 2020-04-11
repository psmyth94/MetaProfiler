# Contents of this file

 * Introduction
 * Citation
 * Installation
 * Examples
 * License

# Introduction
![alt text](https://github.com/psmyth94/MetaProfiler/blob/master/man/logo/logo.png)

This repository contains the source code of the `MetaProfiler` software package.
It provides calculations for local false discovery rates of protein-based stable isotopic probing (SIP) results and performs taxonomic, functional, phylogenetic, and time series analysis of microbiome dynamics.

`MetaProfiler` has only been tested on MetaProSIP results from OpenMS, but it is designed to work with multiple tools that extract heavy peptide features from light peptide identifications. More tools such as ProteinTurnover will be tested in the future.

# Citation

If you use `MetaProfiler` in your projects, please cite the preprint

Patrick Smyth, Xu Zhang, Zhibin Ning, Janice Mayne, Jasmine I Moore, Krystal Walker, Mathieu Lavallée-Adam and Daniel Figeys 2020, *Studying the dynamics of the gut microbiota using metabolically stable isotopic labeling and metaproteomics* [doi:10.1101/982884](https://doi.org/10.1101/2020.03.09.982884)

# Installation

Make sure to have `R >= 3.5.0` installed. Paste the following lines into your `R` session.

```{R}
# instal devtools, if you do not have it.
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# install MetaProfiler via devtools
library(devtools)
install_github("psmyth94/MetaProfiler")
```

# Examples

This an example R script that creates the MetaProfiler class object.

```{R}
  
# the units for the time measurements.
time_unit = "Day"
# name of the incorporation measurements. MetaProSIP calls it RIA.
incorporation_name = "RIA"
# name of the intensity measurements. MetaProSIP calls it INT.
intensity_name = "INT"
# name of the score values. MetaProSIP uses correlation scores.
score_name = "Cor."
# name of the labeling ratio measurements. We will not specify it yet as we will not use the LR values from MetaProSIP.
labeling_ratio_name = NULL
# automatically create the experimental design table using the names in the result directory.
design <- create_experimental_design(
  results_directory = "./Protein_SIP_results", # the directory with the files containing the information about the heavy peptide features.
  Sample = "Ref\\d+", # the sample names
  Day = "(?<=_D)\\d+|(?<=_D\\-)\\d+" # the day the stool sample was collected. Be sure to use the same name as specified in variable time_unit.
)
# lets sort by day and then by sample name.
design <- design[MetaProfiler:::.mixedorder(design[,.(Day,Sample)])]

# load the results to data.
data <- load_results(
  design,
  # The file names containing the information about the heavy peptide features.
  # Does not need to be specified if the `design` contain a column with the filenames.
  data = NULL, 
  time_unit = time_unit,
  # the names of the variables.
  incorporation_name = "RIA",
  intensity_name = "INT",
  labeling_ratio_name = NULL, # we will use the Global Peptide LR column for this value. It will be named LR by default.
  score_name = "Cor.",
  # the columns containing the information.
  # Does not need to be specified as long as the first word is the corresponding variable's name, followed by a unique identifier
  # e.g. [RIA 1, RIA 2, RIA 3, etc] or [RIA light, RIA heavy].
  incorporation_columns = NULL,
  intensity_columns = NULL,
  labeling_ratio_columns = "Global Peptide LR",
  score_columns = NULL,
  as_percentage = TRUE, # whether the incorporation and labeling ratio value should be a decimal or a percentage.
  progress = TRUE # track the progress.
)

# no longer needed
design <- design[,-"filenames"]

# all files were generated using MetaLab (http://dashboard.imetalab.ca/#/).
# The pep2pro, pep2taxon, and pro2func variable can also be a matrix, data.frame, or data.table.
# The function will try to guess the accession, peptide, taxon, and function columns of these tables.
# The peptide, taxon, and (somewhat) function columns are easy to guess, but the accession column can be a little tricky.
# if it fails to guess, you will need to specify the column names.
annotation_table <- make_annotation_table(
  time_unit = time_unit, 
  data = data, # the table containing the heavy peptide information.
  pep2pro = "Examples/closed_search/peptides.txt", # an optional peptide to protein file. However, MetaProSIP does not report razor proteins, thus the use of this parameter. 
  pep2taxon = "Examples/taxonomy_analysis/BuiltIn.pepTaxa.csv", # a peptide to taxon file.
  pro2func = "Examples/functional_analysis/functions.csv" # a peptide to taxon file.
)

# the list of peptide sequence to search for in the annotation table.
seq <- data[["Peptide Sequence (no PTMs)"]]
# search for the sequence in the peptide column of the annotation table and return their corresponding protein.
data[,"Proteins"] <- annotation_table[seq, Proteins, on = "Peptides"]

# lets extract the heavy peptide features
data <- flatten_data(
  data = data,
  incorporation_name = incorporation_name,
  intensity_name = intensity_name,
  score_name = score_name,
  labeling_ratio_name = labeling_ratio_name, # if set to NULL. It will calculate LR from intensity of the heavy and light peptides.
)
labeling_ratio_name = "LR"

data <- cluster_features(
  data = data,
  incorporation_name = incorporation_name,
  intensity_name = intensity_name,
  score_name = score_name,
  labeling_ratio_name = labeling_ratio_name,
  radius = c(1,1),
  distance_method = c("absolute", "absolute"),
  cluster_by = c(incorporation_name, labeling_ratio_name)
)

pdf_list <- get_pdf(
  data = data,
  design = design, 
  time_unit = time_unit,
  time_zero = 0,
  by = time_unit, # Calculate distribution at each time point. Must be specified in table design.
  observations = c(incorporation_name, labeling_ratio_name),
  bandwidth = 5 # size of the bandwidth for the kernel densitiy estimation.
)

data <- calc_LFDR(
  data = data,
  design = design, 
  time_unit = time_unit,
  time_zero = 0,
  by = time_unit, # Calculate distribution at each time point. Must be specified in table design.
  observations = c(incorporation_name, labeling_ratio_name),
  bandwidth = 5 # size of the bandwidth for the kernel densitiy estimation.
)

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
# make the master table where most statistics and analysis are done on. Filter at an LFDR of 10%.
make_master_table(Object, LFDR_threshold = 0.1)
# Verify that the density for false discovery fits well.
vars = c(incorporation_name, labeling_ratio_name)
Object@deconvolve(vars = vars, xlim = setNames(list(c(0,100), c(0,100)), vars))
# Check the distribution of LR.
plot(Object, "LR")

# Check the distribution of RIA.
plot(Object, "RIA")										  
```

# Licensing and contributions
`MetaProfiler` is licensed under the GPL (>= 2) license. Contributions are welcome.