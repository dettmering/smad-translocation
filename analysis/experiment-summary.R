classifiers <- c(
  'Metadata_Dose',
  'Metadata_Treatment',
  'Metadata_Time',
  'Metadata_Celltype'
)

transloc <- generateList(summary, classifiers)

for (i in 1:length(transloc$Metadata_Dose)) {
  temp <- merge(summary, transloc[i, classifiers])
  
  transloc[i, 'Mean.Translocated'] <- mean(temp$Translocated.Percent)
  transloc[i, 'SD.Translocated'] <- sd(temp$Translocated.Percent)
}

transloc$SEM.Translocated <- transloc$SD.Translocated / sqrt(transloc$n)

write.csv(img, paste0(format(Sys.time(), "%Y-%m-%d"), "_exp-summary.csv"), row.names = F)

library(ggplot2)

ggplot(transloc, aes(x = Metadata_Dose, y = Mean.Translocated)) +
  geom_bar(aes(fill = Metadata_Treatment), position = position_dodge(width = 0.9), stat = "identity") +
  geom_errorbar(aes(ymin = Mean.Translocated - SEM.Translocated, ymax = Mean.Translocated + SEM.Translocated, group = Metadata_Treatment), position = position_dodge(width = 0.9), stat = "identity", width = 0.1) +
  geom_text(aes(y = 0, label = n, group = Metadata_Treatment), position = position_dodge(width = 0.9), size = 2, vjust = -1) +
  xlab("Dose (Gy)") +
  ylab("Cells with translocated protein (%)") +
  scale_fill_discrete(name = "Treatment") +
  facet_grid(Metadata_Celltype ~ Metadata_Time)