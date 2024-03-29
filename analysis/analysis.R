###################################
### Smad translocation analysis ###
###################################

# Analysis pipeline for detecting nuclear translocation of Smad2/3 proteins into the nucleus
# (c) 2014 Till Dettmering

# Before starting, set working directory with setwd("DIR") to directory containing CellProfiler output files

# Read CellProfiler results

if (!exists("img")) img <- read.csv("Image.csv")
if (!exists("nuc")) nuc <- read.csv("Nuclei.csv")
if (!exists("cytopl")) cytopl <- read.csv("Cytoplasm.csv")

z.threshold <- 2 # How many SDs must the nuclear signal be away from the cytoplasm signal to count as translocated?

# Load sources

library(devtools) # needed for https source from github

source_url('https://raw.githubusercontent.com/tdett/r-helpers/master/generateList.R')
source_url('https://raw.githubusercontent.com/dettmering/r-helpers/master/dfBackup.R')

# Set classifiers; the results will be separated by these columns

classifiers <- c(
  'Metadata_Dose',
  'Metadata_Treatment',
  'Metadata_Time',
  'Metadata_Celltype',
  'Metadata_Experiment',
  'Metadata_Replicate'
)

image_area_cm2 <- 635.34 * 474.57 / 10^8 # Image size in cm^2. These values are likely different for your microscope.
img$cells_per_cm2 <- img$Count_Nuclei / image_area_cm2

# Add Cytoplasm data

nuc$Cytoplasm_Mean <- cytopl$Intensity_MeanIntensity_OrigPOI
nuc$Cytoplasm_SD <- cytopl$Intensity_StdIntensity_OrigPOI

# Calculate Cytoplasm-to-Nucleus ratio

nuc$Ratio <- nuc$Intensity_MeanIntensity_OrigPOI / nuc$Cytoplasm_Mean

# Calculate z-score

nuc$z.score <- (nuc$Intensity_MeanIntensity_OrigPOI - nuc$Cytoplasm_Mean) / nuc$Cytoplasm_SD

# Determine translocation status

nuc$Translocated <- nuc$z.score > z.threshold # Binary operator: Is Smad translocated?

# Area

nuc$Cells_Area <- nuc$AreaShape_Area + cytopl$AreaShape_Area
nuc$AreaRatio <- nuc$AreaShape_Area / nuc$Cells_Area * 100

# QC

nuc <- subset(nuc, AreaRatio < 100) # Exclude nuclei without cytoplasm

##########################################################
# Calculate percentage translocated and other parameters #
##########################################################

summary <- generateList(nuc, classifiers)

for (i in 1:length(summary$Metadata_Dose)) {
  temp <- merge(nuc, summary[i, classifiers])
  temp_img <- merge(img, summary[i, classifiers])
  
  summary[i, 'Translocated'] <- sum(temp$Translocated, na.rm = TRUE)
  summary[i, 'Mean.Intensity'] <- mean(temp$Intensity_MeanIntensity_OrigPOI, na.rm = TRUE)
  summary[i, 'Median.Intensity'] <- median(temp$Intensity_MeanIntensity_OrigPOI, na.rm = TRUE)
  summary[i, 'SD.Intensity'] <- sd(temp$Intensity_MeanIntensity_OrigPOI, na.rm = TRUE)
  summary[i, 'Mean.Ratio'] <- mean(temp$Ratio, na.rm = TRUE)
  summary[i, 'SD.Ratio'] <- sd(temp$Ratio, na.rm = TRUE)
  
  summary[i, 'Mean.cells_per_cm2'] <- mean(temp_img$cells_per_cm2)
}

summary$CV.Intensity <- summary$SD.Intensity / abs(summary$Mean.Intensity) * 100
summary$Translocated.Percent <- summary$Translocated / summary$n * 100

################################################
# Summarize translocation over all experiments #
################################################

transloc_classifiers <- c(
  'Metadata_Dose',
  'Metadata_Treatment',
  'Metadata_Time',
  'Metadata_Celltype'
)

transloc <- generateList(summary, transloc_classifiers)

for (i in 1:length(transloc$Metadata_Dose)) {
  temp <- merge(summary, transloc[i, transloc_classifiers])
  
  transloc[i, 'Mean.Translocated'] <- mean(temp$Translocated.Percent)
  transloc[i, 'SD.Translocated'] <- sd(temp$Translocated.Percent)
}

transloc$SEM.Translocated <- transloc$SD.Translocated / sqrt(transloc$n)

# Export raw data and summary table to csv (working directory)

dfBackup(c('img', 'summary', 'transloc'))