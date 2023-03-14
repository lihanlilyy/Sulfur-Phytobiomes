# (C) Copyright 2017 Isai Salas Gonzalez
#
#    This file is part of gfobap.
#
#    gfobap is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    gfobap is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with gfobap.  If not, see <http://www.gnu.org/licenses/>.

# args=commandArgs(trailingOnly=TRUE)
# # test if there is at least one argument: if not, return an error
# if (length(args)<1) {
#   stop("One argument must be supplied a matrix file", call.=FALSE)
# }
# mat<-args[1]
library(ape)
library(phytools)
library(readxl)
rhizo_list <- list( 'Actinomycetales',
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
#Sphingobacteriales changed to Sphingobacteriale due to length
phyllo_list <- list('Rhizobiales',
                    'Burkholderiales',
                    'Actinomycetales',
                    'Pseudomonadales',
                    'Bacillales',
                    'Enterobacterales',
                    'Xanthomonadales',
                    'Sphingomonadales',
                    'Flavobacteriales',
                    'Sphingobacteriale',
                    'Mycobacteriales')
i <- "Mycobacteriales"
#for (i in rhizo_list) {
  traits <- paste("traits_scoary_raw_rhizo_soil_", paste(i, ".csv", sep =""), sep = "")
  matrix <- paste("matrix_scoary_raw_rhizo_soil_", paste(i, ".csv", sep =""), sep = "")
  print(i)
  ######Rhizo versus Soil######
#Map<-read.table("/pine/scr/i/s/isai/3837_march_2017/metadata_3837_genomes_june2016_complete.tsv",header=T,sep="\t")
  Tab<-read_excel("/Users/hanli/Desktop/FYP/sulfer_genes/final_matrix_with_metadata.xlsx",sheet = i,col_names = T,skip = 0)
  head(Tab)
  # Remove the sphere column, which is the second column
  Tab <- Tab[,-2]
  Map<-read_excel("/Users/hanli/Desktop/FYP/sulfer_genes/final_matrix_with_metadata.xlsx",sheet = paste("sphere_", i, sep = ""),col_names = T,skip = 0)
  #Map
  Tab[, c('Gene')]
  for(i in 1:nrow(Tab)) {       # for-loop over rows
    Tab[i,'Gene' ] <-sub('.', '', Tab[i,'Gene' ])
  }
  
  for(i in 1:nrow(Map)) {       # for-loop over rows
    Map[i,'Gene' ] <-sub('.', '', Map[i,'Gene' ])
  }
  Map<-as.data.frame(Map)
  Tab<-as.data.frame(Tab)

  Map
  Tab
  #----------------------
  # Gene <- Tab[,1]
  # part2 <- Tab[,2:ncol(Tab)]
  # part2
  # for(i in 1:ncol(part2)) {       # for-loop over rows
  #   part2[,i ] <- ifelse(part2[,i] >=1, 1, 0)
  # }
  # Tab <- cbind(Gene, part2)
  #----------------------
  
  
  
  
  # Tab
  # foreach cell in Tab[, c('Gene')] {
  #   cell = sub('.', '', cell)
  # }
  #Map<-read.table("metadata.tsv",header=T,sep="\t")
  #Tab<-read.table("metadata_full_v2.csv",header=T,sep=",",row.names=1,check.names=F)
  # no need, since already transposed
  #Tab<-t(Tab)
  #The logic is to subset first using the columns of the matrix
  rownames(Map) <- Map$Gene
  Map
  #Tab$Gene
  #rownames(Map)
  #Map
  Map_soil<-as.character(Map$Gene[which(Map$sphere=="Soil")])
  Map_ra<-as.character(Map$Gene[which(Map$sphere=="Rhizosphere")])
  myids<-c(Map_ra,Map_soil)
  #myids
  Map<-Map[match(myids,Map$Gene),]
  Map<-droplevels(Map)
  Map
  head(Tab)
  colnames(Tab)
  #Map$Gene
  #Tab<-Tab[match(rownames(Map),Tab$Gene),]
  #Remove cols that sum to zero
  Tab
  rownames(Tab) <- Tab$Gene
  Tab <- Tab[ ,-1]
  # Tab<-Tab[,which(colSums(Tab)!=0)]
  Tab<-t(Tab)
  #Tab
  colnames(Tab)
  
  Map$sphere
  
  trait<-c(rep(1,length(Map_ra)), rep(0,length(Map_soil)))
  df_traits<-data.frame(TraitY=trait)
  rownames(df_traits)<-Map$Gene
  
  df_Tab<-data.frame(genes=rownames(Tab))
  df_Tab_vals<-as.data.frame(Tab)
  df_Tab<-cbind(df_Tab,df_Tab_vals)
  
  #Read the general phylogenetic tree
  #tree<-read.tree("/pine/scr/i/s/isai/3837_march_2017/3837_genomes_31scg_june2016.newick")
  # tree<-read.tree("isolates_file.tree")
  # tree$tip.label
  #Root the tree first using MRCA
  #Use the firmicutes
  #Paenibacillus Isolate 2517572151
  #Bacillus Isolate 2623620997
  #root_node<-phytools::findMRCA(tree=tree,tips=c("2517572151","2623620997"))
  #tree<-reroot(tree,node.number=root_node)
  #Map$Gene
  #subtree<-drop.tip(phy=tree,tip=which(!(tree$tip.label%in%Map$Gene)))

  print(traits)
  print(matrix)
  
  write.table(df_traits,traits,row.names=T,col.names=NA,sep=",",append=F,quote=F)
  write.table(df_Tab,matrix ,row.names=F,col.names=T,sep=",",append=F,quote=F)
  #write.tree(subtree,"tree_scoary_rhizo_soil_Rhizobiales.newick")
#}
  
  
for (i in rhizo_list) {
  print(i)
  Tab<-read.csv(paste("/Users/hanli/Desktop/FYP/scoary/code/results/rhizo/",paste(i, ".csv", sep=""), sep =""), header = TRUE, sep=",")
  Tab$adj_p.value<-p.adjust(Tab$Empirical_p,method = "fdr")
  write.csv(Tab, paste("/Users/hanli/Desktop/FYP/scoary/code/results/rhizo/adj_",paste(i, ".csv", sep=""), sep =""), row.names=FALSE)
}
  
  for (i in rhizo_list) {
    print(i)
    Tab<-read.csv(paste("/Users/hanli/Desktop/FYP/scoary/code/results/rhizo/binary_",paste(i, ".csv", sep=""), sep =""), header = TRUE, sep=",")
    Tab$adj_p.value<-p.adjust(Tab$Empirical_p,method = "fdr")
    write.csv(Tab, paste("/Users/hanli/Desktop/FYP/scoary/code/results/rhizo/adj_binary",paste(i, ".csv", sep=""), sep =""), row.names=FALSE)
  }
  
  for (i in phyllo_list) {
    print(i)
    Tab<-read.csv(paste("/Users/hanli/Desktop/FYP/scoary/code/results/phyllo/",paste(i, ".csv", sep=""), sep =""), header = TRUE, sep=",")
    Tab$adj_p.value<-p.adjust(Tab$Empirical_p,method = "fdr")
    write.csv(Tab, paste("/Users/hanli/Desktop/FYP/scoary/code/results/phyllo/adj_",paste(i, ".csv", sep=""), sep =""), row.names=FALSE)
  }
  
  for (i in phyllo_list) {
    print(i)
    Tab<-read.csv(paste("/Users/hanli/Desktop/FYP/scoary/code/results/phyllo/binary_",paste(i, ".csv", sep=""), sep =""), header = TRUE, sep=",")
    Tab$adj_p.value<-p.adjust(Tab$Empirical_p,method = "fdr")
    write.csv(Tab, paste("/Users/hanli/Desktop/FYP/scoary/code/results/phyllo/adj_binary",paste(i, ".csv", sep=""), sep =""), row.names=FALSE)
  }
  