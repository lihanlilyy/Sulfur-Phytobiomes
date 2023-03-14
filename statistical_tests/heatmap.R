library(devtools)
library(readxl)
library(gplots)
library(ggplot2)
library(pheatmap)
install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)


df<-read_excel("/Users/hanli/Desktop/FYP/stats/phyloglm/plyloglm_for_plotting.xlsx",sheet = "phyllo",col_names = T,skip = 0)
df<- as.data.frame(df)
df
for(i in 1:nrow(df)) {
  for (j in 1:ncol(df)) {
    if (is.na(df[i,j])) {
      df[i,j]<- 0
    }
  }
}
df
rownames(df) <- df$...1
df <- df[,-1]
df<-as.matrix(df)
# heatmap.2(df, scale = "none", col = bluered(100), 
#           trace = "none", density.info = "none")
mycols <- colorRamp2(breaks = c(-2, 0, 6), 
                     colors = c("blue", "white", "red"))

tiff(file="phyllo_heatmap.tiff")
Heatmap(df, 
        name = "Estimate", #title of legend
        column_title = "Orders", row_title = "Genes",
        row_names_gp = gpar(fontsize = 4.5), # Text size for row names
        column_names_gp = gpar(fontsize =8),
        column_names_centered = TRUE,
        column_title_rot = 0,
        column_names_rot = 0,
        col = mycols
)

dev.off()



        