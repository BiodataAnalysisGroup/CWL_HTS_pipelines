# Load libraries ----------------------------------------------------------

.libPaths("/home/analysis/custom_R_installation_dir")
library(data.table)
library(ballgown)

# Extract arguments -------------------------------------------------------

# get command-line arguments
args = commandArgs(trailingOnly = FALSE)

# Stringtie directories
# multiple paths (strings) as input
# must be the final argument!
stringtie_dirs = grep(
  pattern = "--stringtie-output",
  x = args
)
stringtie_dirs = args[(stringtie_dirs+1):length(args)]

# TSv/CSV file with phenotype(s)
phenotypeFile = grep(
  pattern = "--metadata",
  x = args
)
phenotypeFile = args[phenotypeFile+1]

# phenotype of interest
phenotype = grep(
  pattern = "--phenotype",
  x = args
)
phenotype = args[phenotype+1]

# sample names/IDs
samples = grep(
  pattern = "--samples",
  x = args
)
samples = args[samples+1]

# timecourse
timecourse = grep(
  pattern = "--timecourse",
  x = args
)
# timecourse = args[timecourse+1]

# feature
feature = grep(
  pattern = "--feature",
  x = args
)
feature = args[feature+1]

# measure
measure = grep(
  pattern = "--measure",
  x = args
)
measure = args[measure+1]

# confounders
# take one comma-separated string as input,
# which possibly contains multiple confounders
# matching colnames (variables) in phenotypeFile!
confounders = grep(
  pattern = "--confounders",
  x = args
)
if(length(confounders)!=0){
  confounders = unlist(strsplit(args[confounders+1], split = ","))
} else {
  confounders = NULL
}

# custom model comparison
custom_model = grep(
  pattern = "--custom-model",
  x = args
)

# mod
# formula with variable(s) of interest
# take one comma-separated string as input,
# which possibly contains multiple variables
# matching colnames in phenotypeFile!
mod = grep(
  pattern = "--mod",
  x = args
)
mod = unlist(strsplit(args[mod+1], split = ","))

# mod0
# formula without variable(s) of interest
# take one comma-separated string as input,
# which possibly contains multiple variables
# matching colnames in phenotypeFile!
mod0 = grep(
  pattern = "--mod0",
  x = args
)
mod0 = unlist(strsplit(args[mod0+1], split = ","))

# Test acceptable values for arguments ------------------------------------

# test acceptable values for `feature` argument
if(!feature %in% c("exon", "intron", "transcript")){
  stop('--feature should be set to one of the following values: "exon", "intron", "transcript"')
}
# test acceptable values for `measure` argument
if(!measure %in% c("cov", "FPKM", "rcount", "ucount", "mrcount", "mcov")){
  stop('--measure should be set to one of the following values: "cov", "FPKM", "rcount", "ucount", "mrcount", "mcov"')
}
# test custom model design
if(length(custom_model) != 0 & (length(mod) == 0 | length(mod0) == 0) ){
  stop('--mod and --mod0 should be set with the appropriate variables if --custom-model is used')
}

# # test acceptable values for `timecourse` argument
# if(!is.logical(timecourse)){
#   stop('--timecourse should be set to either TRUE (for comparisons of continuous covariates that are not necessarily time) or FALSE. If TRUE, --phenotype should be the column name of a continuous variable.')
# }

# Create ballgown object --------------------------------------------------

bg = ballgown(
  samples = stringtie_dirs,
  verbose = TRUE,
  meas = "all"
)

# load phenotype data.frame
phenotype_df = as.data.frame(fread(phenotypeFile))
# match sample IDs
idx = match(sampleNames(bg), phenotype_df[,samples])
# reorder
phenotype_df = phenotype_df[idx,]

# add phenotype information
pData(bg) = phenotype_df

# save ballgown class object
save(bg, file = "ballgown.rda")

# Differential expression analysis ----------------------------------------

# check for comparisons with continuous covariate or not

if(length(timecourse) != 0){
  # statistical tests for differential expression in ballgown
  ballgown_de_test = stattest(
    bg, 
    feature=feature, 
    meas=measure,
    covariate = phenotype,
    adjustvars = confounders,
    timecourse=TRUE
  )
} else {
  # statistical tests for differential expression in ballgown
  ballgown_de_test = stattest(
    bg, 
    feature=feature, 
    meas=measure,
    covariate = phenotype,
    adjustvars = confounders,
    timecourse=FALSE
  )
}

# extract feature annotations
if(feature == "transcript") {
  # get annotation fd for trans
  anno_df = bg@expr$trans
} else if(feature == "exon") {
  # get annotation fd for exon
  anno_df = bg@expr$exon
} else if(feature == "intron") {
  # get annotation fd for intron
  anno_df = bg@expr$intron
}

# merge annotation with stattest results
ballgown_de_test = merge(
  x = ballgown_de_test,
  y = anno_df,
  by.x = "id",
  by.y = colnames(anno_df)[1]
)

# save in TSV
fwrite(
  ballgown_de_test,
  "ballgown_DE_results.tsv",
  sep = '\t',
  row.names = FALSE,
  quote = FALSE,
  col.names = TRUE
)

# (Optional) Perform additional custom model comparison
if(length(custom_model) != 0){
  
  # mod formula
  mod_formula = as.formula(paste('~',paste0(mod, collapse = " + ")))
  # mod formula
  mod_zero_formula = as.formula(paste('~',paste0(mod0, collapse = " + ")))
  
  # create design matrices:
  mod_mt = model.matrix(mod_formula, data = pData(bg))
  mod_zero_mt = model.matrix(mod_zero_formula, data = pData(bg))
  
  # check for comparisons with continuous covariate or not
  
  if(length(timecourse) != 0){
    # statistical tests for differential expression in ballgown
    adjusted_results = stattest(
      bg, 
      feature=feature, 
      meas=measure, 
      mod0=mod_zero_mt, 
      mod=mod_mt,
      timecourse = TRUE)
  } else {
    # statistical tests for differential expression in ballgown
    adjusted_results = stattest(
      bg, 
      feature=feature, 
      meas=measure, 
      mod0=mod_zero_mt, 
      mod=mod_mt,
      timecourse = FALSE)
  }
  
  # merge annotation with stattest results
  adjusted_results = merge(
    x = adjusted_results,
    y = anno_df,
    by.x = "id",
    by.y = colnames(anno_df)[1]
  )
  
  # save in TSV
  fwrite(
    adjusted_results,
    'ballgown_CM_DE_results.tsv',
    sep = '\t',
    row.names = FALSE,
    quote = FALSE,
    col.names = TRUE
  )
}
