
library(readxl)
library(ggplot2)
library(dplyr)
library(readr)

df<-read_excel("/Users/hanli/Desktop/FYP/sulfer_genes/plant_soil_metadata_for_plotting.xlsx",sheet = "Final",col_names = T,skip = 0)

df$sphere<-as.factor(df$SPHERE)
cols<-c("#3CB371","#FFD700","#B8860B")

box<-ggplot(df,aes(x= SPHERE,y=CONTAMINATION,fill=SPHERE))+
  geom_boxplot() +
  geom_jitter(color="black", size=0.2, alpha=0.2) +
  scale_fill_manual(values=cols)+
  theme_classic() +
  theme(
    legend.position="none",
    plot.title = element_text(size=20)
  ) +
  ggtitle("") +
  ylab("Completeness")+
  xlab("")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12))
box
ggsave(
  "boxplot_genome_completeness.tiff",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6,
  height = 6,
  units = "in",
  dpi = 400,
)
dev.off()

