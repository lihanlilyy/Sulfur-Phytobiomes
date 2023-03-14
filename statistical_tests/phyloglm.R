


##Trying out phyloGLM for BGC MAGs data
library(ape)
library(phylolm)
library(readxl)
library(xlsx)

rhizo_list <- list('Actinomycetales',
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
i <- "Enterobacterales"
for (i in rhizo_list) {
  #df<-read.table('/Users/hanli/Desktop/FYP/stats/phyloglm/MAGs_final_tree.txt',sep = "\t")
  filee <- paste("/Users/hanli/Desktop/FYP/stats/phyloglm/rhizo_binary_", paste(i, ".tsv", sep =""), sep = "")
  df<-read_excel("/Users/hanli/Desktop/FYP/sulfer_genes/final_matrix_with_metadata.xlsx",sheet = i ,col_names = T,skip = 0)
  df<-as.data.frame(df)
  head(df)
  rownames(df)
  df$niche <- ifelse(df$sphere == 'Soil', 0, 1)
  # colnames(df)
  # colnames(df)<-df[1,]
  for(i in 1:nrow(df)) {       # for-loop over rows
    df[i,'Gene' ] <-sub('.', '', df[i,'Gene' ])
  }
  
  rownames(df) <-df[,1]
  
  # df <- df[-1,]
  df <- df[,-1]
  df <- df[,-1]
  df
  #***************
  for(i in 1:nrow(df)) {       # for-loop over rows
    df[i, ] <- ifelse(df[i, ] >=1, 1, 0)
  }
  #***************
  
  colnames(df)
  # rownames(df) <- df$Assembly
  # names(df)
  # df<-df[,-12]
  df
  tree<-read.tree('/Users/hanli/Desktop/FYP/stats/phyloglm/isolates_file.tree')
  tree$tip.label
  dim(df)
  # df$niche<-rep(c(0,1),c(414,414)) # soil is 0
  # remake : assembly should be nivhe (assembly)
  #
  df$niche
  y<-df$niche
  df[,12]
  names(y)<-rownames(df)
  Res<-NULL
  rownames(df)
  tree$tip.label
  df
  subtree<-drop.tip(phy = tree,tip = which(!(tree$tip.label%in%rownames(df))))
  head(df)
  subtree
  # x<-sample(c(0,1),828,replace = TRUE)# a dummy x value
  x <- df
  x <- as.data.frame(x)
  # names(x)<-rownames(df)
  
  dat<-as.data.frame(cbind(y,x))
  dat
  for (column in 1:ncol(df)) {
    x<-df[,column]
    orthogroup<-colnames(df)[column]
    names(x)<-rownames(df)
    dat<-as.data.frame(cbind(y,x))
    m1<-tryCatch(phyloglm(formula = y~x,data = dat,phy = subtree,method = "logistic_IG10"),error = function(e) list(coefficients = NA))
    if(is.na(coef(m1)[1])){
      res<-data.frame(orthogroup.id = orthogroup,
                      Estimate = NA,
                      SE = NA,
                      z.value = NA,
                      p.value = NA)
      
    }else{
      m1.sum<-summary(m1)
      res<-data.frame(orthogroup.id = orthogroup,
                      Estimate = m1.sum$coefficients["x",1],
                      SE = m1.sum$coefficients["x",2],
                      z.value = m1.sum$coefficients["x",3],
                      p.value = m1.sum$coefficients["x",4])
    }
    rm(m1,m1.sum)
    Res<-rbind(Res,res)
  }
  Res$adj_p.value<-p.adjust(Res$p.value,method = "fdr")
  
  #filee <- paste("/Users/hanli/Desktop/FYP/stats/phyloglm/", paste(i, ".tsv", sep =""), sep = "")
  
  write.csv(Res, toString(filee), row.names=FALSE)
  print('************************')
  print(i)
  print(filee)
  print('************************')

  
  
  #The following works,  but fails to converge for dummy data
  #m<-phyloglm(x~y,data = dat,phy = subtree,method = "logistic_IG10")
}

for (i in rhizo_list) {
  #df<-read.table('/Users/hanli/Desktop/FYP/stats/phyloglm/MAGs_final_tree.txt',sep = "\t")
  filee <- paste("/Users/hanli/Desktop/FYP/stats/phyloglm/rhizo_", paste(i, ".tsv", sep =""), sep = "")
  df<-read_excel("/Users/hanli/Desktop/FYP/sulfer_genes/final_matrix_with_metadata.xlsx",sheet = i ,col_names = T,skip = 0)
  df<-as.data.frame(df)
  head(df)
  rownames(df)
  df$niche <- ifelse(df$sphere == 'Soil', 0, 1)
  # colnames(df)
  # colnames(df)<-df[1,]
  for(i in 1:nrow(df)) {       # for-loop over rows
    df[i,'Gene' ] <-sub('.', '', df[i,'Gene' ])
  }
  
  rownames(df) <-df[,1]
  
  # df <- df[-1,]
  df <- df[,-1]
  df <- df[,-1]
  df
  # for(i in 1:nrow(df)) {       # for-loop over rows
  #   df[i, ] <- ifelse(df[i, ] >=1, 1, 0)
  # }
  colnames(df)
  # rownames(df) <- df$Assembly
  # names(df)
  # df<-df[,-12]
  df
  tree<-read.tree('/Users/hanli/Desktop/FYP/stats/phyloglm/isolates_file.tree')
  tree$tip.label
  dim(df)
  # df$niche<-rep(c(0,1),c(414,414)) # soil is 0
  # remake : assembly should be nivhe (assembly)
  #
  df$niche
  y<-df$niche
  df[,12]
  names(y)<-rownames(df)
  Res<-NULL
  rownames(df)
  tree$tip.label
  df
  subtree<-drop.tip(phy = tree,tip = which(!(tree$tip.label%in%rownames(df))))
  head(df)
  subtree
  # x<-sample(c(0,1),828,replace = TRUE)# a dummy x value
  x <- df
  x <- as.data.frame(x)
  # names(x)<-rownames(df)
  
  dat<-as.data.frame(cbind(y,x))
  dat
  for (column in 1:ncol(df)) {
    x<-df[,column]
    orthogroup<-colnames(df)[column]
    names(x)<-rownames(df)
    dat<-as.data.frame(cbind(y,x))
    m1<-tryCatch(phyloglm(formula = y~x,data = dat,phy = subtree,method = "logistic_IG10"),error = function(e) list(coefficients = NA))
    if(is.na(coef(m1)[1])){
      res<-data.frame(orthogroup.id = orthogroup,
                      Estimate = NA,
                      SE = NA,
                      z.value = NA,
                      p.value = NA)
      
    }else{
      m1.sum<-summary(m1)
      res<-data.frame(orthogroup.id = orthogroup,
                      Estimate = m1.sum$coefficients["x",1],
                      SE = m1.sum$coefficients["x",2],
                      z.value = m1.sum$coefficients["x",3],
                      p.value = m1.sum$coefficients["x",4])
    }
    rm(m1,m1.sum)
    Res<-rbind(Res,res)
  }
  Res$adj_p.value<-p.adjust(Res$p.value,method = "fdr")
  
  #filee <- paste("/Users/hanli/Desktop/FYP/stats/phyloglm/", paste(i, ".tsv", sep =""), sep = "")
  
  write.csv(Res, toString(filee), row.names=FALSE)
  print('************************')
  print(i)
  print(filee)
  print('************************')
  
  
  
  #The following works,  but fails to converge for dummy data
  #m<-phyloglm(x~y,data = dat,phy = subtree,method = "logistic_IG10")
}

