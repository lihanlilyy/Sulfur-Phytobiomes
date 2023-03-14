library(readxl)
library(ggplot2)
library(dplyr)
library(readr)
#Actinomycetales
Actinomycetales<-read_excel("/Users/hanli/Desktop/FYP/sulfer_genes/final_matrix_with_metadata.xlsx",sheet = "Actinomycetales",col_names = T,skip = 0)
#head(Actinomycetales)
map_Actinomycetales<-read_excel("/Users/hanli/Desktop/FYP/sulfer_genes/final_matrix_with_metadata.xlsx",sheet = "sphere_Actinomycetales",col_names = T,skip = 0)
#head(map_Actinomycetales)
sum(Actinomycetales$Gene==map_Actinomycetales$Gene)

Actinomycetales<-as.data.frame(Actinomycetales)
map_Actinomycetales<-as.data.frame(map_Actinomycetales)
head(Actinomycetales)
head(map_Actinomycetales)
rownames(Actinomycetales)<-Actinomycetales$Gene
rownames(map_Actinomycetales)<-map_Actinomycetales$Gene
Actinomycetales<-Actinomycetales[,-1]
Actinomycetales<-Actinomycetales[,-1]
#map_Actinomycetales<-map_Actinomycetales[,-1]
#head(Actinomycetales)
#rownames(map_Actinomycetales)
#head(map_Actinomycetales)
N.soil<-sum(Actinomycetales[map_Actinomycetales$sphere=="Soil",])
N.phyllo<-sum(Actinomycetales[map_Actinomycetales$sphere=="Phyllosphere",])




allgenes_hyp <- function(genes,Map){
  #tree <- Dat$tree
  #Map <- Dat$Map
  #genes <- Dat$genes >=1
  
  N.pgyllo <- sum(genes[ Map$sphere == "Phyllosphere", ])
  N.soil <- sum(genes[ Map$sphere == "Soil", ])
  #N.genomes <- nrow(Map)
  Res <- NULL
  for(gene in 1:ncol(genes)){
    print (gene)
    gene.id <- colnames(genes)[gene]
    Map$gene <- genes[,gene]
    
    #genes.per.class <- aggregate(gene ~ Classification, data = Map,FUN = sum)
    #genomes.per.class <- aggregate(gene ~ Classification, data = Map,FUN = length)
    #N.genes <- sum(Map$gene)
    
    #ftable(gene ~ Classification , Map)
    
    #dhyper(x = 0, m = 2, n = 2, k = 2)
    #dhyper(x = 1, m = 2, n = 2, k = 2)
    #dhyper(x = 2, m = 2, n = 2, k = 2)
    #phyper(q = 0, m = 2, n = 2, k = 2)
    #phyper(q = 1, m = 2, n = 2, k = 2)
    #phyper(q = 2, m = 2, n = 2, k = 2)
    # Following is wrong, counts genomes for zero and genes for 1+
    # pval <- phyper(q = sum(subset(Map, Classification == "PA" & gene > 0)$gene) - 1,
    #                m = sum(Map$gene),
    #                n = nrow(subset(Map,gene == 0)),
    #                k = sum(subset(Map, Classification == "PA")$gene))
    #pval <- 1 - pval
    
    # Binary version
    pval <- phyper(q = nrow(subset(Map, sphere == "Phyllosphere" & gene > 0)) - 1,
                   m = sum(Map$gene > 0),
                   n = nrow(subset(Map,gene == 0)),
                   k = nrow(subset(Map, sphere == "Phyllosphere")),lower.tail=FALSE)
    #pval <- 1 - pval
    #pval[ pval == 0 ] <- 1e-16
    score <- -log10(pval)
    
    #raw counts
    pval2 <- phyper(q = sum(subset(Map, sphere == "Phyllosphere")$gene) - 1, 
                    m = N.phyllo, n = N.soil, k = sum(Map$gene),lower.tail=FALSE)
    #pval2 <- 1 - pval2
    #pval2[ pval2 == 0 ] <- 1e-16 
    
    
    #Test for depletion binary version 
    pval_dep <- phyper(q = nrow(subset(Map, sphere == "Phyllosphere" & gene > 0)),
                       m = sum(Map$gene > 0),
                       n = nrow(subset(Map,gene == 0)),
                       k = nrow(subset(Map, sphere == "Phyllisphere")),
                       lower.tail=TRUE)
    #pval_dep[ pval_dep == 0 ] <- 1e-16
    score_dep<--log10(pval_dep)
    
    #Test for depletion Raw counts version 
    pval2_dep <- phyper(q = sum(subset(Map, sphere == "Phullosphere")$gene),
                        m = N.phyllo, n = N.soil, k = sum(Map$gene),lower.tail=TRUE)
    #pval2_dep[ pval2_dep == 0 ] <- 1e-16
    score2_dep<--log10(pval2_dep)
    
    
    
    #res <- data.frame(gene.id = gene.id, score = score, p.value = pval,
    #                  full.score = -log10(pval2), full.p.value = pval2)
    res <- data.frame(gene.id = gene.id, score_enriched_binary = score, p.value_enriched_binary = pval,
                      score_depletion_binary=score_dep,p.value_depletion_binary=pval_dep,
                      score_enriched_rawcounts = -log10(pval2), p.value_enriched_rawcounts = pval2,
                      score_depletion_rawcounts=score2_dep,p.value_depletion_rawcounts=pval2_dep)
    
    Map$gene <- NULL
    Res <- rbind(Res,res)
  }
  Res<-data.frame(gene.id=Res$gene.id,score_enriched_binary=Res$score_enriched_binary,
                  z.score_enriched_binary= (Res$score_enriched_binary - mean(Res$score_enriched_binary)) / sd(Res$score_enriched_binary),
                  p.value_enriched_binary=Res$p.value_enriched_binary,
                  score_depletion_binary=Res$score_depletion_binary,
                  z.score_depletion_binary=(Res$score_depletion_binary - mean(Res$score_depletion_binary)) / sd(Res$score_depletion_binary),
                  p.value_depletion_binary=Res$p.value_depletion_binary,
                  score_enriched_rawcounts=Res$score_enriched_rawcounts,
                  z.score_enriched_rawcounts=(Res$score_enriched_rawcounts - mean(Res$score_enriched_rawcounts)) / sd(Res$score_enriched_rawcounts),
                  p.value_enriched_rawcounts=Res$p.value_enriched_rawcounts,
                  score_depletion_rawcounts=Res$score_depletion_rawcounts,
                  z.score_depletion_rawcounts=(Res$score_depletion_rawcounts - mean(Res$score_depletion_rawcounts)) / sd(Res$score_depletion_rawcounts),
                  p.value_depletion_rawcounts=Res$p.value_depletion_rawcounts)
  
  return(Res)
}

df<-allgenes_hyp(Actinomycetales,map_Actinomycetales)
#df

#correct them for pvalue corrections
df$adj_p.value_enriched_binary<-p.adjust(df$p.value_enriched_binary,method = "fdr")
df$adj_p.value_depletion_binary<-p.adjust(df$p.value_depletion_binary,method = "fdr")
df$adj_p.value_enriched_rawcounts<-p.adjust(df$p.value_enriched_rawcounts,method = "fdr")
df$adj_p.value_depletion_rawcounts<-p.adjust(df$p.value_depletion_rawcounts,method = "fdr")

write.table(df,"/Users/hanli/Desktop/FYP/stats/hyperg/phyllo_Actinomycetales_final_res_hyperg.tsv",sep = "\t")

