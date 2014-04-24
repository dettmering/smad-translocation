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

ggplot(summary, aes(x = Metadata_Dose, y = Translocated.Percent)) +
  geom_bar(aes(fill = Metadata_Treatment), position = position_dodge(width = 0.9), stat = "identity") +
  xlab("Dose (Gy)") +
  ylab("Cells with translocated protein (%)") +
  scale_fill_discrete(name = "Treatment") +
  facet_grid(Metadata_Celltype ~ Metadata_Time) +
  theme_bw()

dev.off()