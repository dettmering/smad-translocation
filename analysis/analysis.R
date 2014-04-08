###################################
### Smad translocation analysis ###
###################################

# Read CellProfiler results

# Before starting, set working directory with setwd("DIR") to directory containing CellProfiler output files

img <- read.csv("Image.csv")
cells <- read.csv("Cells.csv")
nuc <- read.csv("Nuclei.csv")
cytopl <- read.csv("Cytoplasm.csv")

# Load sources

library(devtools) # needed for https source from github

source_url('https://raw.github.com/tdett/r-helpers/master/generateList.R')

# Set classifiers

classifiers <- c(
  'Metadata_Dose',
  'Metadata_Treatment',
  'Metadata_Time',
  'Metadata_Celltype'
)

# Parse information from file name: Celltype-Dosis-Treatment-Time-Sample-_***.tif

image_area_cm2 <- 635.34 * 474.57 / 10^8
img$cells_per_cm2 <- img$Count_Nuclei / image_area_cm2

# Write required values from img data frame to nuc data frame

for (i in img$ImageNumber) {
  nuc[nuc$ImageNumber == i,'Intensity_Background'] <- img[img$ImageNumber == i,'Intensity_MedianIntensity_MaskBackground']
}

# Add Cytoplasm data

nuc$Cytoplasm_Mean <- cytopl$Intensity_MeanIntensity_OrigPOI
nuc$Cytoplasm_SD <- cytopl$Intensity_StdIntensity_OrigPOI

nuc$Cytoplasm_Mean_Corr <- nuc$Cytoplasm_Mean - nuc$Intensity_Background
nuc$Nucleus_Mean_Corr <- nuc$Intensity_MeanIntensity_OrigPOI - nuc$Intensity_Background

# Calculate Cytoplasm-to-Nucleus ratio

nuc$Ratio <- nuc$Intensity_MeanIntensity_OrigPOI / nuc$Cytoplasm_Mean
nuc$Ratio_Corr <- nuc$Nucleus_Mean_Corr / nuc$Cytoplasm_Mean_Corr

# Determine translocation status

nuc$Translocated <- nuc$Intensity_MeanIntensity_OrigPOI > nuc$Cytoplasm_Mean + (2 * nuc$Cytoplasm_SD) # Binary operator: Is Smad translocated?
nuc$Translocated_Corr <- nuc$Nucleus_Mean_Corr > nuc$Cytoplasm_Mean_Corr + (2 * nuc$Cytoplasm_SD) # Binary operator: Is Smad translocated?

# Area

nuc$Cells_Area <- cells$AreaShape_Area
nuc$AreaRatio <- nuc$AreaShape_Area / nuc$Cells_Area * 100

# QC

nuc <- subset(nuc, AreaRatio < 100) # Exclude nuclei without cytoplasm
nuc <- subset(nuc, AreaRatio < (mean(nuc$AreaRatio) + 3 * sd(nuc$AreaRatio))) # Exclude outliers in ratio, e.g. cells with very small cytoplasm

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

# Export raw data and summary table to csv (working directory)

write.csv(img, paste0(format(Sys.time(), "%Y-%m-%d"), "_rawdata-img.csv"), row.names = F)
if (length(nuc$Metadata_Dose) < 50000) write.csv(nuc, paste0(format(Sys.time(), "%Y-%m-%d"), "_rawdata.csv"), row.names = F)
write.csv(summary, paste0(format(Sys.time(), "%Y-%m-%d"), "_results.csv"), row.names = F)

###########
# Figures #
###########

library(ggplot2)

pdf(paste0(format(Sys.time(), "%Y-%m-%d"), "_results.pdf"), width = 5.83, height = 4.13)

  ggplot(nuc, aes(x = Ratio)) + geom_histogram() + facet_grid(Metadata_Treatment ~ Metadata_Dose) + xlim(c(0,10))
  
  ggplot(nuc, aes(x = Metadata_Dose, y = Ratio)) +
    geom_boxplot(aes(fill = Metadata_Treatment)) +
    geom_hline(yintercept = 1, col = "red") +
    xlab("Dose (Gy)") +
    ylab("Ratio nucleus/cytoplasm") +
    scale_fill_discrete(name = "Condition") +
    facet_grid(. ~ Metadata_Celltype) +
    theme_bw()

  ggplot(nuc, aes(x = Metadata_Dose, y = Intensity_MedianIntensity_OrigPOI)) +
    geom_boxplot(aes(fill = Metadata_Treatment)) +
    xlab("Dose (Gy)") +
    ylab("Median intensity in nucleus (AU)") +
    scale_fill_discrete(name = "Condition") +
    facet_grid(. ~ Metadata_Celltype) +
    theme_bw()

  ggplot(summary, aes(x = Metadata_Dose, y = Translocated.Percent)) +
    geom_bar(aes(fill = Metadata_Treatment), position = position_dodge(width = 0.9)) +
    xlab("Dose (Gy)") +
    ylab("Cells with translocated protein (%)") +
    scale_fill_discrete(name = "Condition") +
    facet_grid(. ~ Metadata_Celltype) +
    theme_bw()

dev.off()