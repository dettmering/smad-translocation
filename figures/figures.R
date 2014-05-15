###########
# Figures #
###########

library(ggplot2)

pdf(paste0(format(Sys.time(), "%Y-%m-%d"), "_results.pdf"), width = 8.27, height = 5.83)

ggplot(nuc, aes(x = Ratio)) + geom_histogram() + facet_grid(Metadata_Treatment ~ Metadata_Dose) + xlim(c(0,10))

ggplot(nuc, aes(x = Metadata_Dose, y = Ratio)) +
  geom_boxplot(aes(fill = Metadata_Treatment)) +
  geom_hline(yintercept = 1, col = "red") +
  xlab("Dose (Gy)") +
  ylab("Ratio nucleus/cytoplasm") +
  scale_fill_discrete(name = "Treatment") +
  facet_grid(Metadata_Celltype ~ Metadata_Time) +
  theme_bw()

ggplot(nuc, aes(x = Metadata_Dose, y = z.score)) +
  geom_boxplot(aes(fill = Metadata_Treatment)) +
  geom_hline(yintercept = z.threshold, col = "red") +
  xlab("Dose (Gy)") +
  ylab("z-score") +
  scale_fill_discrete(name = "Treatment") +
  facet_grid(Metadata_Celltype ~ Metadata_Time)

ggplot(nuc, aes(x = Metadata_Dose, y = Nucleus_Mean_Corr * 2^16)) +
  geom_boxplot(aes(fill = Metadata_Treatment)) +
  xlab("Dose (Gy)") +
  ylab("Mean intensity in nucleus (Grey value)") +
  scale_fill_discrete(name = "Treatment") +
  facet_grid(Metadata_Celltype ~ Metadata_Time) +
  theme_bw()

ggplot(transloc, aes(x = Metadata_Dose, y = Mean.Translocated)) +
  geom_bar(aes(fill = Metadata_Treatment), position = position_dodge(width = 0.9), stat = "identity") +
  geom_errorbar(aes(ymin = Mean.Translocated - SEM.Translocated, ymax = Mean.Translocated + SEM.Translocated, group = Metadata_Treatment), position = position_dodge(width = 0.9), stat = "identity", width = 0.1) +
  geom_text(aes(y = 0, label = n, group = Metadata_Treatment), position = position_dodge(width = 0.9), size = 2, vjust = -1) +
  xlab("Dose (Gy)") +
  ylab("Cells with translocated protein (%)") +
  ylim(c(0, 100)) +
  scale_fill_discrete(name = "Treatment") +
  facet_grid(Metadata_Celltype ~ Metadata_Time)

dev.off()