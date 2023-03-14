
library(readxl)
library(ggplot2)
library(dplyr)
library(readr)
order_list <- list( 'Actinomycetales',
                    'Pseudomonadales', 
                    'Streptomycetales', 
                    'Bacillales', 
                    'Enterobacterales',
                    'Xanthomonadales', 
                    'Sphingomonadales', 
                    'Flavobacteriales',
                    'Paenibacillales', 
                    'Sphingobacteriales', 
                    'Propionibacteriales',
                    'Chitinophagales', 
                    'Caulobacterales',
                    'Mycobacteriales',
                    'Rhizobiales',
                    'Burkholderiales')
Res <- NULL
df1<-read_excel("/Users/hanli/Desktop/FYP/stats/phyloglm/combined_phyloglm_results.xlsx",sheet = i,col_names = T,skip = 0)
df1$order <- "random"
Res <- df1[1,]
Res <- Res[-1,]
Res
for (i in order_list) {
  df<-read_excel("/Users/hanli/Desktop/FYP/stats/phyloglm/combined_phyloglm_results.xlsx",sheet = i,col_names = T,skip = 0)
  #raw_rhizo
  # raw_rhizo <- NULL
  # raw_phyllo <- NULL
  # binary_rhizo <- NULL
  # binary_phyllo <- NULL
  # for (i in 1:nrow(df)) {
  #   if (df[i, "comparison_type"] == "raw" && df[i, "niche"] == "rhizo") {
  #     print(df[i,])
  #     raw_rhizo <- rbind(raw_rhizo, df[i,])
  #   }
  #   if (df[i, "comparison_type"] == "raw" && df[i, "niche"] == "phyllo") {
  #     print(df[i,])
  #     raw_phyllo <- rbind(raw_phyllo, df[i,])
  #   }
  #   if (df[i, "comparison_type"] == "binary" && df[i, "niche"] == "rhizo") {
  #     print(df[i,])
  #     binary_rhizo <- rbind(binary_rhizo, df[i,])
  #   }
  #   if (df[i, "comparison_type"] == "binary" && df[i, "niche"] == "phyllo") {
  #     print(df[i,])
  #     binary_phyllo <- rbind(binary_phyllo, df[i,])
  #   }
  # }
  # 
  # list.dfs <- list(raw_phyllo,raw_rhizo,binary_rhizo,binary_phyllo)
  # for (num in list.dfs) {
  #  print(num)
  # }
  # for (df in 1:length(list.dfs)) {
  # for (i in 1:nrow(df)) {
  #   if (df[i, "p.value"] == "NA") {
  #     print(df[i,])
  #     df[i, "p.value"] <- "1"
  #     print(df[i,])
  #   }
  # 
  #   df$adj_p.value<-p.adjust(df$p.value,method = "fdr")
  df$order <- i
  for (ii in 1:nrow(df)) {
    if (df[ii,"adj_p.value"] != "NA") {
      if(as.numeric(df[ii,"adj_p.value"]) <= 0.05) {
        print(df[ii,])
        Res <- rbind(Res,df[ii,])
      }
    }
  }
  # }

}
Res
write.csv(Res, "/Users/hanli/Desktop/FYP/stats/phyloglm/filter_results_all.csv", row.names=FALSE)
