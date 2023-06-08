# Load libraries ----------------------------------------------------------

.libPaths("/home/analysis/custom_R_installation_dir")
library(ChIPQC)
library(rtracklayer)
library(stringr)
library(data.table)

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

### Metadata file ##

# Contains additional info on treatment samples including (colnames):
# SampleID (important for matching respective files)
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

### ChIPQC options ###

# consensus
consensus = grep(
  pattern = "--consensus",
  x = args
)

# bCount
bCount = grep(
  pattern = "--bCount",
  x = args
)

# summits
summits_idx = grep(
  pattern = "--summits",
  x = args
)
summits = as.integer(args[summits_idx+1])

# annotation
annotation_idx = grep(
  pattern = "--annotation",
  x = args
)
if(length(annotation_idx)!=0){
  annotation_genome = args[annotation_idx+1]
  # load the 'TxDb.Hsapiens.UCSC.hg38.knownGene' library if the annotation genome points to the hg38 version
  if(annotation_genome == "hg38"){
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  }
}

# chromosomes
chromosomes_idx = grep(
  pattern = "--chromosomes",
  x = args
)
# add values
chromosomes = c()
if(length(chromosomes_idx) != 0){
  # iterate
  for(i in args[(chromosomes_idx+1):length(args)]){
    # check for trailing
    if(!str_detect(string = basename(i), pattern = "^--")){
      chromosomes = c(chromosomes, i)
    } else {
      break
    }
  }
}

# blacklist
blacklist_idx = grep(
  pattern = "--blacklist",
  x = args
)
blacklistFile = args[blacklist_idx+1]

# mapQCth
mapQCth_idx = grep(
  pattern = "--mapQCth",
  x = args
)
mapQCth = as.integer(args[mapQCth_idx+1])

# profileWin
profileWin_idx = grep(
  pattern = "--profileWin",
  x = args
)
profileWin = as.integer(args[profileWin_idx+1])

# fragmentLength
fragmentLength_idx = grep(
  pattern = "--fragmentLength",
  x = args
)
fragmentLength = as.integer(args[fragmentLength_idx+1])

# bParallel
bParallel = grep(
  pattern = "--bParallel",
  x = args
)

# facetBy for built-in HTML report
facetBy_idx = grep(
  pattern = "--facetBy",
  x = args
)
if(length(facetBy_idx) != 0){
  # add values
  facetBy = c()
  # iterate
  for(i in args[(facetBy_idx+1):length(args)]){
    # check for trailing
    if(!str_detect(string = basename(i), pattern = "^--")){
      facetBy = c(facetBy, i)
    } else {
      break
    }
  }
}

# reportFolder
reportFolder_idx = grep(
  pattern = "--reportFolder",
  x = args
)
reportFolder = args[reportFolder_idx+1]

# ChIPQC_experiment
chipqc_experiment_idx = grep(
  pattern = "--chipqc-experiment",
  x = args
)
chipqc_experiment = args[chipqc_experiment_idx+1]

# ChIPQC_report
chipqc_report_idx = grep(
  pattern = "--chipqc-report",
  x = args
)
chipqc_report = args[chipqc_report_idx+1]

# -------------------------------------------------------------------------
# Construct metadata table for ChIPQCexperiment object --------------------
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

# BiocParallel fix
BiocParallel::register(BiocParallel::SerialParam())

# Construct ChIPQCexperiment object ---------------------------------------

# set values for ChIPQC options:

# consensus
if(length(consensus) == 0){
  consensus = FALSE
} else {
  consensus = TRUE
}

# bCount
if(length(bCount) == 0){
  bCount = FALSE
} else {
  bCount = TRUE
}

# bParallel
if(length(bParallel) == 0){
  bParallel = FALSE
} else {
  bParallel = TRUE
}
# blacklist
if(length(blacklistFile) == 0){
  blacklist = NULL
} else {
  blacklist = rtracklayer::import(blacklistFile)
}

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

# Execute ChIPQC
if(length(annotation_idx) == 0){
  
  # without genomic feature analysis
  ChIPQC_result = ChIPQC(
    experiment = metadata_dt,
    consensus = consensus,
    bCount = bCount,
    summits = summits,
    blacklist = blacklist,
    bParallel = bParallel,
    mapQCth = mapQCth,
    profileWin = profileWin,
    fragmentLength = fragmentLength
  )
} else {
  
  # with genomic feature analysis
  # check for chromosomes
  if(length(chromosomes) == 0){
    # default to 1st chromosome
    ChIPQC_result = ChIPQC(
      experiment = metadata_dt,
      consensus = consensus,
      bCount = bCount,
      summits = summits,
      annotation = annotation_genome, 
      blacklist = blacklist,
      bParallel = bParallel,
      mapQCth = mapQCth,
      profileWin = profileWin,
      fragmentLength = fragmentLength
    )
    
  } else {
    
    # selected or all chromosomes
    if(chromosomes == "all"){
      # use all chromosomes
      ChIPQC_result = ChIPQC(
        experiment = metadata_dt,
        consensus = consensus,
        bCount = bCount,
        summits = summits,
        annotation = annotation_genome, 
        chromosomes = NULL,
        blacklist = blacklist,
        bParallel = bParallel,
        mapQCth = mapQCth,
        profileWin = profileWin,
        fragmentLength = fragmentLength
      )
      
    } else if(length(chromosomes) > 1){
      
      # convert chromosome index to integer
      chromosomes = as.integer(chromosomes)
      # user selected chromosome
      ChIPQC_result = ChIPQC(
        experiment = metadata_dt,
        consensus = consensus,
        bCount = bCount,
        summits = summits,
        annotation = annotation_genome, 
        chromosomes = chromosomes,
        blacklist = blacklist,
        bParallel = bParallel,
        mapQCth = mapQCth,
        profileWin = profileWin,
        fragmentLength = fragmentLength
      )
      
    } else {
      
      # run with user selected chromosome
      ChIPQC_result = ChIPQC(
        experiment = metadata_dt,
        consensus = consensus,
        bCount = bCount,
        summits = summits,
        annotation = annotation_genome, 
        chromosomes = chromosomes,
        blacklist = blacklist,
        bParallel = bParallel,
        mapQCth = mapQCth,
        profileWin = profileWin,
        fragmentLength = fragmentLength
      )
      
    }
    
  }
}

# Custom fix for ChIPQC::QCmetrics error:

# Error in names(res) <- c("Reads", "Map%", "Filt%", "Dup%", "ReadL", "FragL",  : 
#                            'names' attribute [9] must be the same length as the vector [7]

QCmetrics_manual <- function (object) {
  res = c(
    ifelse(length(reads(object, FALSE)) != 0, reads(object, FALSE), 0), # 1
    
    ifelse(length(signif((mapped(object) / reads(object,
                                                 FALSE)) * 100, 3)) != 0,
           signif((mapped(object) / reads(object,
                                          FALSE)) * 100, 3), 
           0), # 2
    
    ifelse(length(signif((
      1 - reads(object, TRUE) / reads(object,
                                      FALSE)
    ) * 100, 3)) != 0, signif((
      1 - reads(object, TRUE) / reads(object,
                                      FALSE)
    ) * 100, 3), 0) , # 3
    
    ifelse(length(signif(duplicateRate(object) * 100, 3)) != 0, 
           signif(duplicateRate(object) * 100, 3), 0) , # 4
    
    ifelse(length(readlength(object)) != 0, readlength(object), 0) , # 5
    
    ifelse(
      length(fragmentlength(object, width = readlength(object))) != 0, 
      fragmentlength(object, width = readlength(object)), 
      0), # 6
    
    ifelse(length(signif(RelativeCrossCoverage(object), 3)) != 0, 
           signif(RelativeCrossCoverage(object), 3), 
           0), # 7
    
    ifelse(length(signif(ssd(object), 3)) != 0, signif(ssd(object), 3), 0), # 8
    
    ifelse(
      length(signif(frip(object) * 100, 3)) != 0, 
      signif(frip(object) * 100, 3), 
      0) # 9
  )
  names(res) = c("Reads",
                 "Map%",
                 "Filt%",
                 "Dup%",
                 "ReadL",
                 "FragL",
                 "RelCC",
                 "SSD",
                 "RiP%")
  blk = ribl(object)
  if (!is.na(blk)[1]) {
    names(blk) <- "RiBL%"
    blk = signif(blk / res[1] * 100, 3)
    res = c(res, blk)
  }
  return(res)
}

# QCmetrics
QCmetrics_list = lapply(
  X = ChIPQC_result@Samples,
  FUN = QCmetrics_manual
)
# matrix
QCmetrics_mt = t(do.call(cbind, QCmetrics_list))

# 1. Save ChIPQCexperiment object
save(ChIPQC_result, file = chipqc_experiment)

# Convert ChIPQC report to data.frame
ChIPQC_report = cbind(
  QCmetadata(ChIPQC_result),
  QCmetrics_mt # QCmetrics(ChIPQC_result)
)

# 2. Save table report in CSV
fwrite(
  ChIPQC_report,
  file = chipqc_report,
  sep = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# 3. Save build-in HTML report
try(
  ChIPQCreport(
    ChIPQC_result,
    facetBy=facetBy,
    reportFolder=reportFolder
  )
)
