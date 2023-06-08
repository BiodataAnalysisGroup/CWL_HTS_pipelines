# Load libraries ----------------------------------------------------------

.libPaths("/home/analysis/custom_R_installation_dir")
library(DESeq2)
library(dplyr)
library(apeglm)
library(data.table)

# Extract arguments -------------------------------------------------------

# get command-line arguments
args = commandArgs(trailingOnly = FALSE)

# count matrix
count_matrix_file = grep(
  pattern = "--count-matrix",
  x = args
)
count_matrix_file = args[count_matrix_file+1]

# phenotype + potential batch effects
metadata = grep(
  pattern = "--metadata",
  x = args
)
metadata = args[metadata+1]

# sample names/IDs
samples = grep(
  pattern = "--samples",
  x = args
)
samples = args[samples+1]

# pass formula as design to the DESeqDataSet object
formula_design = grep(
  pattern = "--design",
  x = args
)
formula_design = args[formula_design+1]

# minimum sum of reads for pre-filtering
min_sum_of_reads = grep(
  pattern = "--min-sum-of-reads",
  x = args
)
min_sum_of_reads = as.integer(args[min_sum_of_reads+1])

# reference level for factor of phenotype of interest
reference_level = grep(
  pattern = "--reference-level",
  x = args
)
reference_level = args[reference_level+1]

# phenotype of interest (colname)
phenotype = grep(
  pattern = "--phenotype",
  x = args
)
phenotype = args[phenotype+1]

# contrast
de_contrast = grep(
  pattern = "--contrast",
  x = args
)

# numerator
numerator = grep(
  pattern = "--numerator",
  x = args
)
numerator = args[numerator+1]

# denominator
denominator = grep(
  pattern = "--denominator",
  x = args
)
denominator = args[denominator+1]

# lfcThreshold
lfcThreshold = grep(
  pattern = "--lfcThreshold",
  x = args
)
lfcThreshold = as.numeric(args[lfcThreshold+1])

# pAdjustMethod
pAdjustMethod = grep(
  pattern = "--pAdjustMethod",
  x = args
)
pAdjustMethod = args[pAdjustMethod+1]

# alpha
alpha = grep(
  pattern = "--alpha",
  x = args
)
alpha = as.numeric(args[alpha+1])

# parallelization
parallelization = grep(
  pattern = "--parallelization",
  x = args
)

# cores
cores = grep(
  pattern = "--cores",
  x = args
)
cores = as.integer(args[cores+1])

# transformation
transformation = grep(
  pattern = "--transformation",
  x = args
)
transformation = args[transformation+1]

# blind
blind = grep(
  pattern = "--blind",
  x = args
)
if(length(blind)!=0){
  blind = TRUE
} else {
  blind = FALSE
}

# hypothesis test
hypothesis_test = grep(
  pattern = "--hypothesis",
  x = args
)
hypothesis_test = args[hypothesis_test+1]

# reduced variable
reduced = grep(
  pattern = "--reduced",
  x = args
)
reduced = args[reduced+1]

# check if performing analysis for hidden batches
hidden_batches = grep(
  pattern = "--hidden-batch-effects",
  x = args
) 

# hidden_batches_rm
hidden_batches_rm = grep(
  pattern = "--hidden-batch-rowMeans",
  x = args
)
hidden_batches_rm = as.integer(args[hidden_batches_rm+1])

# method for hidden batches analysis (SVA or RUV)
hidden_batches_method = grep(
  pattern = "--hidden-batch-method",
  x = args
)
hidden_batches_method = args[hidden_batches_method+1]

# number of variables to be calculated by hidden batches analysis
variables = grep(
  pattern = "--hidden-batch-variables",
  x = args
)
variables = as.integer(args[variables+1])

# Pre-process count matrix and metadata -----------------------------------

# load count matrix
count_matrix = as.data.frame(fread(count_matrix_file, skip = 1))

# remove rows with missing gene ids
count_matrix = subset(count_matrix, Geneid != "")

# set gene ids as rownames
rownames(count_matrix) = count_matrix$Geneid

# keep sample columns
idx_samples = grep(colnames(count_matrix), pattern = "_sorted.bam$")
count_matrix = count_matrix[,idx_samples]

# edit colnames
colnames(count_matrix) = basename(colnames(count_matrix)) %>% 
  gsub(pattern = "_sorted.bam$", replacement = "")

# load metadata
coldata = as.data.frame(fread(metadata))

# set sample ids as rownames
rownames(coldata) = coldata[,samples]

# Ensure same order between:
# 1. count matrix colnames 
# 2. coldata rownames
idx_order = match(colnames(count_matrix), rownames(coldata))
coldata = coldata[idx_order,]

# convert coldata columns to class factor 
for(i in 2:ncol(coldata)){
  coldata[,i] = as.factor(coldata[,i])
}

# Set parameters for DE comparisons ---------------------------------------

# Create design (formula) if missing, with the last column of the metadata table
if(length(phenotype) == 0){
  phenotype = tail(colnames(coldata),1)
}
if(length(formula_design) == 0){
  formula_design = paste0("~",phenotype)
}

# create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = as.formula(formula_design)
)

# minimum sum of reads per gene for all samples
keep <- rowSums(counts(dds)) >= min_sum_of_reads
dds <- dds[keep,]

# # (Optional) set reference level to phenotype of interest
# if(length(reference_level) != 0){
#   # set reference level
#   colData(dds)[,phenotype] <- relevel(colData(dds)[,phenotype], ref = reference_level)  
# }

# drop levels of condition without samples
colData(dds)[,phenotype] <- droplevels(colData(dds)[,phenotype])

# Differential gene expression analysis -----------------------------------

# Check for parallelization
if(length(parallelization) != 0){
  # set parallel
  idx_parallel = TRUE
  # load BiocParallel
  library(BiocParallel)
  # register cores
  register(MulticoreParam(cores)) # default: 2
} else {
  idx_parallel = FALSE
}

# Check for hypothesis test
if(hypothesis_test == "Wald"){
  # Wald test
  dds = DESeq(dds, test = "Wald", parallel = idx_parallel)
} else if(hypothesis_test == "LRT"){
  # check for reduced formula
  if(length(reduced) != 0) {
    # Likelihood ratio test test
    dds = DESeq(dds, test = "LRT", parallel = idx_parallel, reduced = as.formula(paste0("~", reduced)))
  } else {
    stop("Set --reduced to the variable(s) of interest, which will be removed from the full formula design.")
  }
} else {
  stop("--hypothesis should be set to either 'Wald' or 'LRT'.")
}

# Perform DE analysis -----------------------------------------------------
# (Optional) Removing hidden batch effects --------------------------------

# check for selection of SVA for hidden batch effects correction
if(length(hidden_batches)!=0){
  # check methods
  if(hidden_batches %in% c("SVA","RUV")){
    # Hidden batch effects - SVA
    if(hidden_batches_method == "SVA"){
      # SVA
      library(sva)
      # extract normalized counts
      dat  <- counts(dds, normalized = TRUE)
      # average count across samples is larger than 10
      idx  <- rowMeans(dat) > hidden_batches_rm
      dat  <- dat[idx, ]
      mod  <- model.matrix(as.formula(paste0("~",phenotype)), colData(dds))
      mod0 <- model.matrix(~1, colData(dds))
      svseq <- svaseq(dat, mod, mod0, n.sv = variables)
      # add to colData and edit design formula
      svseq_df = as.data.frame(svseq$sv)
      colnames(svseq_df) = paste0("SV", as.character(seq(variables)))
      colData(dds) = cbind(colData(dds), svseq_df)
      design(dds) = as.formula(
        paste0(
          "~",
          paste(
            c( colnames(svseq_df), phenotype), collapse = "+"
          )
        )
      )
      
      # Check for hypothesis test
      if(hypothesis_test == "Wald"){
        # Wald test
        dds = DESeq(dds, test = "Wald", parallel = idx_parallel)
      } else if(hypothesis_test == "LRT"){
        # check for reduced formula
        if(length(reduced) != 0) {
          # Likelihood ratio test test
          dds = DESeq(dds, test = "LRT", parallel = idx_parallel, reduced = as.formula(paste0("~", reduced)))
        } else {
          stop("Set --reduced to the variable(s) of interest, which will be removed from the full formula design.")
        }
      } else {
        stop("--hypothesis should be set to either 'Wald' or 'LRT'.")
      }
      
      # Specify contrasts (depends on phenotype) --------------------------------
      idx_args = length(de_contrast) != 0 & length(denominator) != 0 & length(numerator) != 0
      #
      if(idx_args){
        # check levels of phenotype
        idx_label = numerator %in% levels(colData(dds)[,phenotype]) | denominator %in% levels(colData(dds)[,phenotype])
        #
        if(!idx_label){
          stop("--numerator and --denominator must be set to values of labels existing in the phenotype (either manually-defined or last metadata column)")
        }
        # set contrast vector
        de_contrast = c(
          phenotype,
          numerator,
          denominator
        )
        # execute DE analysis
        res <- results(
          dds, 
          contrast = de_contrast,
          lfcThreshold = lfcThreshold,
          pAdjustMethod = pAdjustMethod,
          alpha = alpha,
          parallel = idx_parallel
        )
      } else {
        # execute DE analysis
        res <- results(
          dds, 
          # contrast = de_contrast,
          lfcThreshold = lfcThreshold,
          pAdjustMethod = pAdjustMethod,
          alpha = alpha,
          parallel = idx_parallel
        )
      }
      
      # Hidden batch effects - RUV
    } else if(hidden_batches_method == "RUV"){
      
      # Specify contrasts (depends on phenotype) --------------------------------
      idx_args = length(de_contrast) != 0 & length(denominator) != 0 & length(numerator) != 0
      #
      if(idx_args){
        # check levels of phenotype
        idx_label = numerator %in% levels(colData(dds)[,phenotype]) | denominator %in% levels(colData(dds)[,phenotype])
        #
        if(!idx_label){
          stop("--numerator and --denominator must be set to values of labels existing in the phenotype (either manually-defined or last metadata column)")
        }
        # set contrast vector
        de_contrast = c(
          phenotype,
          numerator,
          denominator
        )
        # execute DE analysis
        res <- results(
          dds, 
          contrast = de_contrast,
          lfcThreshold = lfcThreshold,
          pAdjustMethod = pAdjustMethod,
          alpha = alpha,
          parallel = idx_parallel
        )
      } else {
        # execute DE analysis
        res <- results(
          dds, 
          # contrast = de_contrast,
          lfcThreshold = lfcThreshold,
          pAdjustMethod = pAdjustMethod,
          alpha = alpha,
          parallel = idx_parallel
        )
      }
      
      # RUV
      library(RUVSeq)
      # extract counts in newSeqExpressionSet
      set <- newSeqExpressionSet(counts(dds))
      # idx  <- rowSums(counts(set) > 5) >= 2 # at least 2 samples with >5 read counts
      idx <- rowMeans(counts(set)) > hidden_batches_rm
      set  <- set[idx, ]
      set <- betweenLaneNormalization(set, which="upper")
      not.sig <- rownames(res)[which(res$pvalue > .1)]
      empirical <- rownames(set)[ rownames(set) %in% not.sig ]
      set <- RUVg(set, empirical, k= variables)
      # add to colData and edit design formula
      ruv_df = pData(set)
      colnames(ruv_df) = paste0("W_", as.character(seq(variables))) # redundant...
      colData(dds) = cbind(colData(dds), ruv_df)
      design(dds) = as.formula(
        paste0(
          "~",
          paste(
            c( colnames(ruv_df), phenotype), collapse = "+"
          )
        )
      )
      # Re-run DE analysis:
      # Check for hypothesis test
      if(hypothesis_test == "Wald"){
        # Wald test
        dds = DESeq(dds, test = "Wald", parallel = idx_parallel)
      } else if(hypothesis_test == "LRT"){
        # check for reduced formula
        if(length(reduced) != 0) {
          # Likelihood ratio test test
          dds = DESeq(dds, test = "LRT", parallel = idx_parallel, reduced = as.formula(paste0("~", reduced)))
        } else {
          stop("Set --reduced to the variable(s) of interest, which will be removed from the full formula design.")
        }
      } else {
        stop("--hypothesis should be set to either 'Wald' or 'LRT'.")
      }
      
      # Specify contrasts (depends on phenotype) --------------------------------
      idx_args = length(de_contrast) != 0 & length(denominator) != 0 & length(numerator) != 0
      #
      if(idx_args){
        # check levels of phenotype
        idx_label = numerator %in% levels(colData(dds)[,phenotype]) | denominator %in% levels(colData(dds)[,phenotype])
        #
        if(!idx_label){
          stop("--numerator and --denominator must be set to values of labels existing in the phenotype (either manually-defined or last metadata column)")
        }
        # set contrast vector
        de_contrast = c(
          phenotype,
          numerator,
          denominator
        )
        # execute DE analysis
        res <- results(
          dds, 
          contrast = de_contrast,
          lfcThreshold = lfcThreshold,
          pAdjustMethod = pAdjustMethod,
          alpha = alpha,
          parallel = idx_parallel
        )
      } else {
        # execute DE analysis
        res <- results(
          dds, 
          # contrast = de_contrast,
          lfcThreshold = lfcThreshold,
          pAdjustMethod = pAdjustMethod,
          alpha = alpha,
          parallel = idx_parallel
        )
      }
      
    }
  } else {
    stop("--hidden-batches-method should be set to 'SVA' or 'RUV'")
  }
} else {
  # Run DE analysis without any hidden batch effect correction analysis
  # Specify contrasts (depends on phenotype) --------------------------------
  idx_args = length(de_contrast) != 0 & length(denominator) != 0 & length(numerator) != 0
  #
  if(idx_args){
    # check levels of phenotype
    idx_label = numerator %in% levels(colData(dds)[,phenotype]) | denominator %in% levels(colData(dds)[,phenotype])
    #
    if(!idx_label){
      stop("--numerator and --denominator must be set to values of labels existing in the phenotype (either manually-defined or last metadata column)")
    }
    # set contrast vector
    de_contrast = c(
      phenotype,
      numerator,
      denominator
    )
    # execute DE analysis
    res <- results(
      dds, 
      contrast = de_contrast,
      lfcThreshold = lfcThreshold,
      pAdjustMethod = pAdjustMethod,
      alpha = alpha,
      parallel = idx_parallel
    )
  } else {
    # execute DE analysis
    res <- results(
      dds, 
      # contrast = de_contrast,
      lfcThreshold = lfcThreshold,
      pAdjustMethod = pAdjustMethod,
      alpha = alpha,
      parallel = idx_parallel
    )
  }
}

# order p-values and adjusted p-values
resOrdered <- res[order(res$pvalue, decreasing = FALSE),]

# Log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef=tail(resultsNames(dds),1), type="apeglm")

# Transformation of read counts for visualization -------------------------

# (Optional) check and apply transformation
if(length(transformation) != 0){
  # check type
  if(transformation == "vst"){
    # vst
    transformed_counts <- varianceStabilizingTransformation(dds, blind=blind)  
  } else if(transformation == "rlog"){
    # rlog
    transformed_counts <- rlog(dds, blind=blind)
  } else {
    # stop
    stop("Set --transformation to either 'vst' or 'rlog'.")
  }
}

# Store DESeqDataSet objects and results ----------------------------------

# save DESeqDataSet R object
save(
  dds,
  file = "DESeq2_DESeqDataSet.rda"
)

# save lfcShrink as R object
save(
  resLFC, 
  file = "DESeq2_lfcShrink.rda"
)

# save transformed R object
if(length(transformation) != 0){
  save(
    transformed_counts,
    file = "DESeq2_DESeqTransform.rda"
  )
}

# save DESeq2 results
res_df = cbind(
  data.frame(genes = rownames(resOrdered)),
  as.data.frame(resOrdered)
)
# TSV
fwrite(
  res_df,
  file = "DESeq2_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
