# Load libraries ----------------------------------------------------------

.libPaths("/home/analysis/custom_R_installation_dir")
library(DiffBind)
library(rtracklayer)
library(stringr)
library(data.table)
library(DESeq2)

# -------------------------------------------------------------------------
# Extract arguments -------------------------------------------------------
# -------------------------------------------------------------------------

# get command-line arguments
args = commandArgs(trailingOnly = FALSE)

# treatment ΒΑΜ files 
treatmentBAM_idx = grep(
  pattern = "--treatment-bam",
  x = args
)
# add values
treatmentBAM = c()
# iterate
for(i in args[(treatmentBAM_idx+1):length(args)]){
  # check for trailing
  if(!str_detect(string = basename(i), pattern = "^--")){
    treatmentBAM = c(treatmentBAM, i)
  } else {
    break
  }
}

# control BAM files
controlBAM_idx = grep(
  pattern = "--control-bam",
  x = args
)
if(length(controlBAM_idx) != 0){
  # add values
  controlBAM = c()
  # iterate
  for(i in args[(controlBAM_idx+1):length(args)]){
    # check for trailing
    if(!str_detect(string = basename(i), pattern = "^--")){
      controlBAM = c(controlBAM, i)
    } else {
      break
    }
  }
} else {
  controlBAM = NULL
}

# files with peak information (.bed, .narrowPeak, etc.)
peakFiles_idx = grep(
  pattern = "--peak-files",
  x = args
)
# add values
peakFiles = c()
# iterate
for(i in args[(peakFiles_idx+1):length(args)]){
  # check for trailing
  if(!str_detect(string = basename(i), pattern = "^--")){
    peakFiles = c(peakFiles, i)
  } else {
    break
  }
}

# PeakCaller
peakCaller_idx = grep(
  pattern = "--peak-caller",
  x = args
)
peakCaller = args[peakCaller_idx+1]

### Metadata file ###

# Contains additional info on treatment samples including (colnames):
# SampleID (important for matching respective files!!!)
# Tissue
# Factor
# Condition
# Treatment
# Replicate
# ControlID

metadataFile_idx = grep(
  pattern = "--metadata",
  x = args
)
metadataFile = args[metadataFile_idx+1]

### DiffBind options ###

# consensus
consensus_idx = grep(
  pattern = "--consensus",
  x = args
)
if(length(consensus_idx) != 0) {
  # add values
  consensus = c()
  # iterate
  for (i in args[(consensus_idx + 1):length(args)]) {
    # check for trailing
    if (!str_detect(string = basename(i), pattern = "^--")) {
      consensus = c(consensus, i)
    } else {
      break
    }
  }
}

# minOverlap
minOverlap_idx = grep(
  pattern = "--minOverlap",
  x = args
)
minOverlap = as.numeric(args[minOverlap_idx+1])

# summits
summits_idx = grep(
  pattern = "--summits",
  x = args
)
summits = as.integer(args[summits_idx+1])

# blacklist
blacklist_idx = grep(
  pattern = "--blacklist",
  x = args
)
if(length(blacklist_idx) !=0){
  # add values
  blacklistFile = c()
  # iterate
  for (i in args[(blacklist_idx + 1):length(args)]) {
    # check for trailing
    if (!str_detect(string = basename(i), pattern = "^--")) {
      blacklistFile = c(blacklistFile, i)
    } else {
      break
    }
  }
}

# greylist
greylist_idx = grep(
  pattern = "--greylist",
  x = args
)
if(length(greylist_idx) !=0){
  # add values
  greylistFile = c()
  # iterate
  for (i in args[(greylist_idx + 1):length(args)]) {
    # check for trailing
    if (!str_detect(string = basename(i), pattern = "^--")) {
      greylistFile = c(greylistFile, i)
    } else {
      break
    }
  }
}

# cores
cores_idx = grep(
  pattern = "--cores",
  x = args
)
cores = as.integer(args[cores_idx+1])

# fragmentSize
fragmentSize_idx = grep(
  pattern = "--fragmentSize",
  x = args
)
fragmentSize = as.integer(args[fragmentSize_idx+1])

# bParallel
bParallel = grep(
  pattern = "--bParallel",
  x = args
)

# normalization method
norm_idx = grep(
  pattern = "--normalization",
  x = args
)
normalizationMethod = get(args[norm_idx+1])

# library calculation  method
lib_idx = grep(
  pattern = "--library",
  x = args
)
libraryMethod = get(args[lib_idx+1])

# background normalization 
background = grep(
  pattern = "--background",
  x = args
)

# spikein normalization (parallel factor)
spikein_idx = grep(
  pattern = "--spikein",
  x = args
)
if(length(spikein_idx) != 0){
  spikein = args[spikein_idx+1]
}
if(length(spikein_idx) !=0){
  # add values
  spikein = c()
  # iterate
  for (i in args[(spikein_idx + 1):length(args)]) {
    # check for trailing
    if (!str_detect(string = basename(i), pattern = "^--")) {
      spikein = c(spikein, i)
    } else {
      break
    }
  }
}

# formula design
formula_design_idx = grep(
  pattern = "--design",
  x = args
)
if(length(formula_design_idx) != 0){
  formula_design = args[formula_design_idx+1]
}

# contrast
contrast_idx = grep(
  pattern = "--contrast",
  x = args
)
if(length(contrast_idx) != 0) {
  # add values
  contrast = c()
  # iterate
  for (i in args[(contrast_idx + 1):length(args)]) {
    # check for trailing
    if (!str_detect(string = basename(i), pattern = "^--")) {
      contrast = c(contrast, i)
    } else {
      break
    }
  }
}

# reorderMeta_factor
reorderMeta_factor_idx = grep(
  pattern = "--reorderMeta-factor",
  x = args
)
if(length(reorderMeta_factor_idx) != 0) {
  # add values
  reorderMeta_factor = c()
  # iterate
  for (i in args[(reorderMeta_factor_idx + 1):length(args)]) {
    # check for trailing
    if (!str_detect(string = basename(i), pattern = "^--")) {
      reorderMeta_factor = c(reorderMeta_factor, i)
    } else {
      break
    }
  }
}

# reorderMeta_value
reorderMeta_value_idx = grep(
  pattern = "--reorderMeta-value",
  x = args
)
if(length(reorderMeta_value_idx) != 0) {
  # add values
  reorderMeta_value = c()
  # iterate
  for (i in args[(reorderMeta_value_idx + 1):length(args)]) {
    # check for trailing
    if (!str_detect(string = basename(i), pattern = "^--")) {
      reorderMeta_value = c(reorderMeta_value, i)
    } else {
      break
    }
  }
}

# retrieve or not consensus peaks used in the analysis (BED format)
retrieve_consensus_idx = grep(
  pattern = "--retrieve-consensus",
  x = args
)
if(length(retrieve_consensus_idx) != 0){
  retrieve_consensus = TRUE
} else {
  retrieve_consensus = FALSE
}

# DiffBind results
diffbind_results_idx = grep(
  pattern = "--diffbind-results",
  x = args
)
diffbind_results = args[diffbind_results_idx+1]

# DiffBind correlation heatmap
correlation_heatmap_idx = grep(
  pattern = "--correlation-heatmap",
  x = args
)
correlation_heatmap = args[correlation_heatmap_idx+1]

# DiffBind consensus
diffbind_consensus_idx = grep(
  pattern = "--diffbind-consensus",
  x = args
)
diffbind_consensus = args[diffbind_consensus_idx+1]

# DiffBind low read count filter (for dba.count function)
low_read_count_filter_idx = grep(
  pattern = "--low-read-count-filter",
  x = args
)
low_read_count_filter = args[low_read_count_filter_idx+1]

# DiffBind low read count filter method (for dba.count function)
diffbind_filterFun_idx = grep(
  pattern = "--diffbind-filterFun",
  x = args
)
diffbind_filterFun = args[diffbind_filterFun_idx+1]

# -------------------------------------------------------------------------
# Construct metadata table for DBA object ---------------------------------
# -------------------------------------------------------------------------

# Peak and BAM files will be automatically added to the metadata table,
# based on the "SampleID" and "ControlID" columns!

# load metadata table
metadata_dt = fread(metadataFile)

# add treatment peak files
treatment_peakFiles = sapply(
  X = metadata_dt$SampleID,
  FUN = function(j){
    # find respective file
    peakFile_idx = grep(pattern = j, x = basename(peakFiles))
    # return file
    return(peakFiles[peakFile_idx])
  }
)
all(metadata_dt$SampleID == names(treatment_peakFiles))
# TRUE
#
metadata_dt$Peaks = treatment_peakFiles

# add treatment BAM files
treatmentBAM_reordered = sapply(
  X = metadata_dt$SampleID,
  FUN = function(j){
    # find respective file
    treatmentBAM_idx = grep(pattern = j, x = basename(treatmentBAM))
    # return file
    return(treatmentBAM[treatmentBAM_idx])
  }
)
all(metadata_dt$SampleID == names(treatmentBAM_reordered))
# TRUE
#
metadata_dt$bamReads = treatmentBAM_reordered

# add control BAM files
if(!is.null(controlBAM)){
  controlBAM_reordered = sapply(
    X = metadata_dt$ControlID,
    FUN = function(j){
      # find respective file
      controlBAM_idx = grep(pattern = j, x = basename(controlBAM))[1]
      # return file
      return(controlBAM[controlBAM_idx])
    }
  )
  all(metadata_dt$ControlID == names(controlBAM_reordered))
  # TRUE
  #
  metadata_dt$bamControl = controlBAM_reordered
}

# add PeakCaller column
metadata_dt$PeakCaller = peakCaller

# # fix for acceptable values of data.table with the make.names function
# metadata_dt_fix = apply(
#   X = metadata_dt[,!colnames(metadata_dt) %in% c("bamReads", "Peaks", "bamControl", "Replicate"), with = FALSE],
#   MARGIN = 2,
#   FUN = make.names
# )
# metadata_dt = cbind(
#   metadata_dt_fix,
#   metadata_dt[,colnames(metadata_dt) %in% c("bamReads", "Peaks", "bamControl", "Replicate"), with = FALSE]
# )

# reorder data.frame input
expected_column_order = c(
  'SampleID',
  'Tissue',
  'Factor',
  'Condition',
  'Treatment',
  'Replicate',
  'bamReads',
  'bamControl',
  'Spikein',
  'ControlID',
  'Peaks',
  'PeakCaller',
  'PeakFormat',
  'ScoreCol',
  'LowerBetter'
)
expected_column_order = expected_column_order[expected_column_order %in% colnames(metadata_dt)]
# reorder data.table
setcolorder(metadata_dt, expected_column_order)

# -------------------------------------------------------------------------
# Perform differential binding analysis -----------------------------------
# -------------------------------------------------------------------------

# Load the data -----------------------------------------------------------

### Create initial DBA object ###

# default peak consensus with: minOverlap=2
# check for consensus option and minOverlap:

# if only minOverlap (different from default) and no consensus specified
if(minOverlap != 2 & length(consensus_idx) == 0){
  dba_object = dba(sampleSheet = metadata_dt, minOverlap = minOverlap)
} else {
  dba_object = dba(sampleSheet = metadata_dt, minOverlap = 2)
}

### Blacklists and greylists ###

# For a blacklist:
# a. Set to TRUE for auto-detection
# b. Assign a GRanges object for custom blacklist from a BED file
# c. Specify ready-to-use ENCODE blacklists from within DiffBind:

# DBA_BLACKLIST_HG19: Homo sapiens 19 (chromosomes have "chr")
# DBA_BLACKLIST_HG38: Homo sapiens 38 (chromosomes have "chr")
# DBA_BLACKLIST_GRCH37: Homo sapiens 37 (chromosomes are numbers
# DBA_BLACKLIST_GRCH38: Homo sapiens 38 (chromosomes are numbers)
# DBA_BLACKLIST_MM9: Mus musculus 9
# DBA_BLACKLIST_MM10: Mus musculus 10
# DBA_BLACKLIST_CE10: C. elegans 10
# DBA_BLACKLIST_CE11: C. elegans 11
# DBA_BLACKLIST_DM3: Drosophila melanogaster 3
# DBA_BLACKLIST_DM6: Drosophila melanogaster 6

# Configure blacklist value
if(length(blacklist_idx) !=0){
  
  if(is.null(blacklistFile)){
    # auto-detect
    blacklistFile = TRUE
    
  } else {
    
    if(str_detect(string = blacklistFile, pattern = "^DBA_")){
      # DiffBind package - DBA constant variable
      blacklistFile = get(blacklistFile)
    } else {
      # GRanges object from BED file
      blacklistFile = rtracklayer::import(blacklistFile)
    }
  }
} else {
  blacklistFile = FALSE
}

# For a greylist: 
# If equal to TRUE, the control bam files (if present), 
# will be examined to determine an appropriate reference genome. 
# Genomes associated with a valid BSgenome can be detected. 
# If successful, this genome will be used to generate greylists for each 
# available control (eg specified as bamControls in the sample sheet.)).

# Configure greylist value
if(length(greylist_idx) !=0){
  
  if(is.null(greylistFile)){
    # auto-detect
    greylistFile = TRUE
    
  } else {
    
    if(str_detect(string = greylistFile, pattern = "^DBA_")){
      # DiffBind package - DBA constant variable
      greylistFile = get(greylistFile)
    } else {
      # GRanges object from BED file
      greylistFile = rtracklayer::import(greylistFile)
    }
  }
} else {
  greylistFile = FALSE
}

# Assign blacklist and greylist regions to DBA object:
# If blacklist equals to TRUE or an ENCODE blacklist option then
# blacklist is used to determine genome type/version before greylist and 
# overwrites greylist if the options differ.

dba_object <- dba.blacklist(dba_object,
                            blacklist = blacklistFile,
                            greylist = greylistFile,
                            cores = cores)

# Counting reads (affinity scores) ----------------------------------------

# check for parallel execution
if(length(bParallel) != 0){
  bParallel = TRUE
} else {
  bParallel = FALSE
}

### Calculate consensus ###

# 1. minOverlap parameter for simple selection of number peakSets required to 
# overlap a peak to be selected
### OR ###
# 2. consensus parameter for custom consensus peaks design 
# (e.g., Tissue & Condition masks)

if(length(consensus_idx) != 0){
  # calculate consensus
  dba_object_consensus <- dba.peakset(
    dba_object,
    consensus = sapply(X = consensus, FUN = base::get),
    minOverlap = minOverlap
  )
  # create DBA object with new consensus
  dba_object_consensus <- dba(dba_object_consensus,
                              mask=dba_object_consensus$masks$Consensus,
                              minOverlap=1)
  # retrieve consensus
  consensus_peaks <- dba.peakset(dba_object_consensus, bRetrieve=TRUE)
  
  # Counting reads
  dba_object = dba.count(
    dba_object, 
    summits = summits, 
    fragmentSize = fragmentSize,
    peaks = consensus_peaks,
    bParallel = bParallel,
    filter = low_read_count_filter,
    filterFun = diffbind_filterFun)

} else {
  
  # Counting reads
  dba_object = dba.count(
    dba_object, 
    summits = summits, 
    fragmentSize = fragmentSize,
    bParallel = bParallel,
    filter = low_read_count_filter,
    filterFun = diffbind_filterFun)
}

# # Export Proportion of Retained Sites
# rate.filterFun <- dba.count(
#   dba_object, 
#   summits = FALSE,
#   peaks = NULL,
#   filter = 1:250,
#   filterFun = diffbind_filterFun)
# 
# # export TIF 
# tiff(
#   filename = "Proportion_of_Retained_Sites_per_Filter_value.tif",
#   compression = "lzw",
#   res = 300,
#   units = "in",
#   width = 8,
#   height = 8
# )
# plot(
#   0:250,
#   rate.filterFun / rate.filterFun[1],
#   type = 'l',
#   xlab = "Filter value",
#   ylab = "Proportion Retained Sites"
# )
# abline(v = low_read_count_filter, col = "blue", lwd=2, lty=2)
# abline(
#   h = (rate.filterFun / rate.filterFun[1])[low_read_count_filter],
#   col = "red",
#   lwd = 2,
#   lty = 2
# )
# # Add a legend
# legend(x = "topright", 
#        legend=c(
#          paste0("Filter value: ", as.character(low_read_count_filter)), 
#          paste0("Proportion of Retained Sites: ", as.character(round((rate.filterFun/rate.filterFun[1])[low_read_count_filter], 2)))
#        ),
#        col=c("blue", "red"), 
#        lty=c(2,2), 
#        cex=0.8)
# dev.off()

# Normalizing the data ----------------------------------------------------

# assign background normalization
if(length(background) != 0){
  background = TRUE
} else {
  background = FALSE
}

# assign spikein for parallel factor normalization 
# application of spike-in with foreign genome alignment (e.g., from Drosophila) 
# will be implemented in the future (requires additional steps in the CWL workflow...)
if(length(spikein_idx) != 0){
  if(length(spikein) == 0){
    # use spike-in tracks from the metadata table (one for each sample)
    spikein = TRUE
  } else {
    # load from BED file
    spikein = rtracklayer::import(spikein)
  }
} else {
  spikein = FALSE
}

# apply normalization
dba_object <- dba.normalize(
  # DBA object
  dba_object,
  # core normalization methods
  normalize = normalizationMethod,
  # library size normalization methods
  library = libraryMethod,
  # background normalization with any methods
  # Note that computing background reads requires access to the full sequencing data 
  # (BAM files)
  background = background,
  # spike-in or parallel factor normalization
  spikein = spikein
)

# Establishing a model design and contrast --------------------------------

# Utilizing GLM functionality included in DESeq2 and edgeR 
# so that the confounding factor can be modeled, enabling more sensitive results:

# create reorderMeta list with named items
condition1 = !(is.null(reorderMeta_value) & is.null(reorderMeta_factor))
condition2 = length(reorderMeta_value) == length(reorderMeta_factor)

# check reorderMeta
# with reorderMeta
if(condition1 & condition2){
  # check contrast
  # with specific contrast
  if( length(contrast_idx) != 0 ){

    # name and assign to list
    names(reorderMeta_value) = reorderMeta_factor
    reorderMeta_list = as.list(reorderMeta_value)
    
    # contrasts
    dba_object <- dba.contrast(
      dba_object, 
      # valid formula design
      design = formula_design,
      # contrasts of character vector of length three
      contrast = contrast, 
      reorderMeta = reorderMeta_list
    )
    
    # without specific contrast
  } else {
    
    # name and assign to list
    names(reorderMeta_value) = reorderMeta_factor
    reorderMeta_list = as.list(reorderMeta_value)
    
    # contrasts
    dba_object <- dba.contrast(
      dba_object, 
      # valid formula design
      design = formula_design,
      reorderMeta = reorderMeta_list
    )
  }
# without reorderMeta
} else {
  # check contrast
  # with specific contrast
  if( length(contrast_idx) != 0 ){
    
    # contrasts
    dba_object <- dba.contrast(
      dba_object, 
      # valid formula design
      design = formula_design,
      # contrasts of character vector of length three
      contrast = contrast
    )
    
    # without specific contrast
  } else {
    # contrasts
    dba_object <- dba.contrast(
      dba_object, 
      # valid formula design
      design = formula_design
    )
  }
}

# 6. Performing the differential analysis ---------------------------------
if(length(blacklist_idx) != 0){
  bBlacklist = TRUE
}
if(length(greylist_idx) != 0){
  bGreylist = TRUE
}
#
dba_object <- dba.analyze(
  dba_object, 
  bParallel = bParallel,
  bBlacklist = bBlacklist,
  bGreylist = bGreylist
)

# # show results and plot 
# dba.show(dba_object, bContrasts=TRUE)
# plot(dba_object, contrast=1)

# 7. Retrieving the differentially bound sites ----------------------------
dba_object.DB <- dba.report(dba_object)

# store results in TSV file

# extract differentially bound peaks
df_peaks = data.frame(
  seqnames = seqnames(dba_object.DB),
  start = start(dba_object.DB),
  end = end(dba_object.DB),
  strand = strand(dba_object.DB)
)
# extract metrics
df_metrics = as.data.frame(mcols(dba_object.DB))

# save in TSV file
fwrite(
  x = cbind(df_peaks, df_metrics),
  file = diffbind_results, #"DiffBind_diff_peaks.tsv",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

# save correlation heatmap
tiff(
  filename = correlation_heatmap, #"DiffBind_diff_peaks_corr_heatmap.tif",
  units = "in",
  width = 8,
  height = 8,
  compression = 'lzw',
  res = 300
)
plot(dba_object, contrast = 1)
dev.off()

# (Optional) save consensus peaks in BED file
if (retrieve_consensus) {
  if(length(consensus_idx) != 0){
    # retrieve consensus
    rtracklayer::export.bed(consensus_peaks, diffbind_consensus) #"DiffBind_diff_consensus_peaks.bed")
  } else {
    # retrieve consensus
    consensus_peaks <- dba.peakset(dba_object, bRetrieve=TRUE)
    rtracklayer::export.bed(consensus_peaks, diffbind_consensus) #"DiffBind_diff_consensus_peaks.bed")
  }
}

# (Optional) retrieve table with normalized counts:
normCounts <- dba.peakset(dba_object, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
fwrite(
  x = normCounts,
  file = "DiffBind_normalized_counts.csv",
  quote = FALSE,
  sep = ","
)

# (Optional) retrieve DBA object
save(
  dba_object,
  file = "DiffBind_DBA_object.rda"
)
