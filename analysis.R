###################################
### Smad translocation analysis ###
###################################

# Read CellProfiler results

setwd(wdir)

img <- read.csv("DefaultOUT_Image.csv")
cells <- read.csv("DefaultOUT_Cells.csv")
nuc <- read.csv("DefaultOUT_Nuclei.csv")
cytopl <- read.csv("DefaultOUT_Cytoplasm.csv")

# Load sources

library(devtools) # needed for https source from github

source_url('https://raw.github.com/tdett/r-helpers/master/generateList.R')

# Set classifiers

classifiers <- c(
  'Dose',
  'Treatment',
  'Time.after.irr',
  'Celltype'
)

# Parse information from file name: Celltype-Dosis-Treatment-Time-Sample-_***.tif

image_area_cm2 <- 635.34 * 474.57 / 10^8
img$cells_per_cm2 <- img$Count_Nuclei / image_area_cm2

res <- strsplit(as.character(img$FileName_OrigDNA), '-')
res <- do.call(rbind, res)

img$Celltype <- res[,1]
img$Dose <- res[,2]
img$Treatment <- res[,3]
img$Time.after.irr <- res[,4]
img$Sample <- res[,5]

# Write parsed information into nuclei data frame

for (i in img$ImageNumber) {
  nuc[nuc$ImageNumber == i,'Dose'] <- img[img$ImageNumber == i,'Dose']
  nuc[nuc$ImageNumber == i,'Treatment'] <- img[img$ImageNumber == i,'Treatment']
  nuc[nuc$ImageNumber == i,'Celltype'] <- img[img$ImageNumber == i,'Celltype']
  nuc[nuc$ImageNumber == i,'Time.after.irr'] <- img[img$ImageNumber == i,'Time.after.irr']
  nuc[nuc$ImageNumber == i,'Sample'] <- img[img$ImageNumber == i,'Sample']
  nuc[nuc$ImageNumber == i,'Intensity_Background'] <- img[img$ImageNumber == i,'Intensity_MedianIntensity_MaskBackground']
}

#nuc$Dose <- factor(nuc$Dose, levels = c('2ab', '0 Gy', '0.5 Gy', '6 Gy', 'TGFb'))

# Add Cytoplasm data

nuc$Cytoplasm_Mean <- cytopl$Intensity_MeanIntensity_OrigPOI
nuc$Cytoplasm_SD <- cytopl$Intensity_StdIntensity_OrigPOI

nuc$Cytoplasm_Mean_Corr <- nuc$Cytoplasm_Mean - nuc$Intensity_Background
nuc$Nucleus_Mean_Corr <- nuc$Intensity_MeanIntensity_OrigPOI - nuc$Intensity_Background

nuc$Ratio <- nuc$Intensity_MeanIntensity_OrigPOI / nuc$Cytoplasm_Mean
nuc$Ratio_Corr <- nuc$Nucleus_Mean_Corr / nuc$Cytoplasm_Mean_Corr
nuc[nuc$Ratio_Corr < 0, 'Ratio_Corr'] <- NA
nuc[nuc$Ratio_Corr > 50, 'Ratio_Corr'] <- NA

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

for (i in 1:length(summary$Dose)) {
  temp <- merge(nuc, summary[i, classifiers])
  temp_img <- merge(img, summary[i, classifiers])
  
  summary[i, 'Translocated'] <- sum(temp$Translocated)
  summary[i, 'Mean.Intensity'] <- mean(temp$Intensity_MeanIntensity_OrigPOI)
  summary[i, 'Median.Intensity'] <- median(temp$Intensity_MeanIntensity_OrigPOI)
  summary[i, 'SD.Intensity'] <- sd(temp$Intensity_MeanIntensity_OrigPOI)
  summary[i, 'Mean.Ratio'] <- mean(temp$Ratio)
  summary[i, 'SD.Ratio'] <- sd(temp$Ratio)
  
  summary[i, 'Mean.cells_per_cm2'] <- mean(temp_img$cells_per_cm2)
}

summary$CV.Intensity <- summary$SD.Intensity / abs(summary$Mean.Intensity) * 100
summary$Translocated.Percent <- summary$Translocated / summary$n * 100

###########
# Figures #
###########

library(ggplot2)

ggplot(nuc, aes(x = Ratio)) + geom_histogram() + facet_grid(Treatment ~ Dose) + xlim(c(0,10))

ggplot(nuc, aes(x = Dose, y = Ratio)) +
  geom_boxplot(aes(fill = Treatment)) +
  annotate("segment", x = 0, xend = 10, y = 1, yend = 1, colour = "red") +  
  facet_grid(. ~ Celltype) +
  theme_bw()

ggplot(summary, aes(x = Dose, y = Translocated.Percent)) +
  geom_bar(aes(fill = Treatment), position = position_dodge(width=0.9)) +
  facet_grid(. ~ Celltype) + 
  theme_bw()