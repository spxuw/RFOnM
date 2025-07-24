library(igraph) # 4.2.0
library(ggplot2)
library(ggpubr)
library(reshape2)
library(NetSci)

rm(list = ls())

setwd("path_for_code")
thres = c(0.85,0.86,0.87,0.88,0.89,0.91,0.92,0.93,0.94,0.95)
upper_size = c(500,600,700,800,900,1000)

results <- data.frame()

for (disease in c('Asthma','Alzheimer','COPD','Diabetes', 'GS-BRCA', 'GS-COAD', 'GS-GBM', 'GS-LGG', 'GS-OV', 'C9', 'MAPT', 'GRN')){
  print(disease)
  theta = read.csv(paste('../results/RFOnM/State_',disease,'.csv',sep = ""),sep=',',header = F)
  gene_name = read.csv(paste('../data/features_',disease,'_head.csv',sep = ""), sep=',')
  adj = read.csv(paste('../data/graph_',disease,'.csv',sep = ""), sep=',',header = F)
  
  sin_values = sin(theta)
  cos_values = cos(theta)
  
  mean_zz = c()
  max_thre = c()
  for (uppers in upper_size){
    mean_z = c()
    for (thre in thres){
      counts = rowSums(sin_values>thre | cos_values>thre)
      Hc = which.min(abs(counts - uppers))
      disease_gene <- which((sin_values[Hc,] > thre) | (cos_values[Hc,] > thre), arr.ind = TRUE)
      g <- graph_from_adjacency_matrix(as.matrix(adj), mode = "undirected")
      subg <- induced_subgraph(g, disease_gene[,2])
      clu <- igraph::components(subg)
      largest_component <- which(clu$csize==max(clu$csize))
      node_ids_largest_component <- which(as.vector(clu$membership)%in%largest_component)
      node_ids_largest_component <- as.numeric(sub("V", "", names(clu$membership)[node_ids_largest_component]))
      mean_z = c(mean_z, max(sum(gene_name$v1[node_ids_largest_component])/length(node_ids_largest_component),  
                             sum(gene_name$v2[node_ids_largest_component])/length(node_ids_largest_component)))
    }
    mean_zz = c(mean_zz,mean(mean_z))
    max_thre = c(max_thre,thres[which.max(mean_z)])
  }
  
  optimal_size = upper_size[which.max(mean_zz)]
  thre = max_thre[which.max(mean_zz)]
  counts = rowSums(sin_values>thre | cos_values>thre)
  Hc = which.min(abs(counts - optimal_size))
  
  disease_gene <- which((sin_values[Hc,] > thre) | (cos_values[Hc,] > thre), arr.ind = TRUE)
  g <- graph_from_adjacency_matrix(as.matrix(adj), mode = "undirected")
  subg <- induced_subgraph(g, disease_gene[,2])
  clu <- igraph::components(subg)
  largest_component <- which(clu$csize==max(clu$csize))
  node_ids_largest_component <- which(as.vector(clu$membership)%in%largest_component)
  node_ids_largest_component <- as.numeric(sub("V", "", names(clu$membership)[node_ids_largest_component]))
  disease_gene <- node_ids_largest_component

  disease_genes_expr <- order(-gene_name$v1)[1:length(disease_gene)]
  disease_genes_wgs <- order(-gene_name$v2)[1:length(disease_gene)]
  
  disease_robust0 <- read.csv(paste0('../results/robust/0_', disease, '.csv'), sep = ',', header = TRUE)
  disease_robust1 <- read.csv(paste0('../results/robust/1_', disease, '.csv'), sep = ',', header = TRUE)
  disease_robust0 <- disease_robust0$vertex+1
  disease_robust1 <- disease_robust1$vertex+1
  
  diease_diamond0 <- read.csv(paste0('../results/diamond/first_100_added_nodes_weight_', disease, '_0_1.txt'), sep = '\t', header = TRUE)
  diease_diamond1 <- read.csv(paste0('../results/diamond/first_100_added_nodes_weight_', disease, '_1_1.txt'), sep = '\t', header = TRUE)
  diease_diamond0 <- diease_diamond0$DIAMOnD_node+1
  diease_diamond1 <- diease_diamond1$DIAMOnD_node+1
  
  if (disease %in% c("GS-LGG")){
    diease_DOMINO0 <- read.csv(paste0('../results/DOMINO/String/seed_',disease,'_1/','seed_',disease,'_1/', 'modules.out'), header = FALSE,sep = ',',na.strings=c("","NA"))
  } else{
    diease_DOMINO0 <- read.csv(paste0('../results/DOMINO/String/seed_',disease,'_0/','seed_',disease,'_0/', 'modules.out'), header = FALSE,sep = ',',na.strings=c("","NA"))
  }
  if (disease %in% c("COPD","GS-COAD","GRN")){
    diease_DOMINO1 = diease_DOMINO0
  } else {
    diease_DOMINO1 <- read.csv(paste0('../results/DOMINO/String/seed_',disease,'_1/','seed_',disease,'_1/', 'modules.out'), header = FALSE,sep = ',',na.strings=c("","NA"))
  }
  
  
  diease_DOMINO0 = melt(as.matrix(diease_DOMINO0))
  diease_DOMINO0 = diease_DOMINO0[complete.cases(diease_DOMINO0),]
  diease_DOMINO0 = data.frame(V1 = diease_DOMINO0$value)
  diease_DOMINO0$V1 = gsub("\\[","",diease_DOMINO0$V1)
  diease_DOMINO0$V1 = gsub("\\]","",diease_DOMINO0$V1)
  diease_DOMINO0 = gsub(" ","",diease_DOMINO0$V1)
  
  diease_DOMINO1 = melt(as.matrix(diease_DOMINO1))
  diease_DOMINO1 = diease_DOMINO1[complete.cases(diease_DOMINO1),]
  diease_DOMINO1 = data.frame(V1 = diease_DOMINO1$value)
  diease_DOMINO1$V1 = gsub("\\[","",diease_DOMINO1$V1)
  diease_DOMINO1$V1 = gsub("\\]","",diease_DOMINO1$V1)
  diease_DOMINO1 = gsub(" ","",diease_DOMINO1$V1)
  
  diease_DOMINO0 <- which(rownames(gene_name) %in% diease_DOMINO0==TRUE)+1
  diease_DOMINO1 <- which(rownames(gene_name) %in% diease_DOMINO1==TRUE)+1
  
  G <- graph_from_adjacency_matrix(as.matrix(adj), mode='undirected')
  
  # rfonm ========================================================
  LCC_i = LCC_Significance(N=1000,Targets =  paste("V",disease_gene,sep = ""),G=G)
  results <- rbind(results, data.frame(Dataset = disease, Method = "RFOnM", ZScore = (LCC_i$Z)))
  
  # disease_genes_expr ========================================================
  LCC_i = LCC_Significance(N=1000,Targets =  paste("V",disease_genes_expr,sep = ""),G=G)
  results <- rbind(results, data.frame(Dataset = disease, Method = "Expression", ZScore = (LCC_i$Z)))
  
  # disease_genes_gws ========================================================
  LCC_i = LCC_Significance(N=1000,Targets =  paste("V",disease_genes_wgs,sep = ""),G=G)
  results <- rbind(results, data.frame(Dataset = disease, Method = "GWAS", ZScore = (LCC_i$Z)))
  
  # disease_robust0 ========================================================
  LCC_i = LCC_Significance(N=1000,Targets =  paste("V",disease_robust0,sep = ""),G=G)
  results <- rbind(results, data.frame(Dataset = disease, Method = "ROBUST (expression)", ZScore = (LCC_i$Z)))
  
  # disease_robust1 ========================================================                                                    
  LCC_i = LCC_Significance(N=1000,Targets =  paste("V",disease_robust1,sep = ""),G=G)
  results <- rbind(results, data.frame(Dataset = disease, Method = "ROBUST (GWAS)", ZScore = (LCC_i$Z)))
  
  LCC_i = LCC_Significance(N=1000,Targets =  paste("V",diease_DOMINO0,sep = ""),G=G)
  results <- rbind(results, data.frame(Dataset = disease, Method = "DOMINO (GWAS)", ZScore = (LCC_i$Z)))
  
  LCC_i = LCC_Significance(N=1000,Targets =  paste("V",diease_DOMINO1,sep = ""),G=G)
  results <- rbind(results, data.frame(Dataset = disease, Method = "DOMINO (expression)", ZScore = (LCC_i$Z)))
  
  LCC_i = LCC_Significance(N=1000,Targets =  paste("V",diease_diamond0,sep = ""),G=G)
  results <- rbind(results, data.frame(Dataset = disease, Method = "DIAMonD (GWAS)", ZScore = (LCC_i$Z)))
  
  LCC_i = LCC_Significance(N=1000,Targets =  paste("V",diease_diamond1,sep = ""),G=G)
  results <- rbind(results, data.frame(Dataset = disease, Method = "DIAMonD (expression)", ZScore = (LCC_i$Z)))
  
}

results$Method = factor(results$Method,levels = c("RFOnM","Expression","GWAS","ROBUST (expression)","ROBUST (GWAS)",
                                            "DIAMonD (expression)","DIAMonD (GWAS)",
                                            "DOMINO (expression)","DOMINO (GWAS)"))
results$Dataset = factor(results$Dataset,levels = c("Asthma","Alzheimer","COPD","Diabetes","GS-BRCA","GS-COAD","GS-GBM","GS-LGG","GS-OV","C9","MAPT","GRN"))


g1 = ggplot(results, aes(x = Method, y = ZScore, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("#5E81ACFF", "#8FA87AFF", "#BF616AFF", "#E7D202FF", "#7D5329FF", "#F49538FF", "#66CDAAFF", "#D070B9FF", "#98FB98FF", "#FCA3B7FF"))+
  theme_bw()+xlab("")+ylab("Z-score")+  facet_wrap(~Dataset, scales = "free", nrow = 2)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 10),
        axis.text.x = element_text(angle = 40, hjust=1),
        strip.background = element_blank(),
        legend.position = "none") 


#####################################################################################################################################
set_to_set_percentile_distance <- function(graph, setA, setB) {
  # Check inputs
  all_nodes <- V(graph)$name
  
  # Filter sets to valid nodes
  setA <- intersect(setA, all_nodes)
  setB <- intersect(setB, all_nodes)
  
  if (length(setA) == 0 || length(setB) == 0) {
    warning("One or both gene sets are empty after filtering to graph nodes.")
    return(NA)
  }
  
  # Function to convert a vector of distances to percentile ranks
  rank_row <- function(row) {
    finite <- is.finite(row)
    ranks <- rep(NA, length(row))
    ranks[finite] <- rank(row[finite], ties.method = "average") / sum(finite)
    return(ranks)
  }
  
  # Compute A → B
  dists_AB <- distances(graph, v = setA, to = setB, mode = "all")
  percentiles_AB <- t(apply(dists_AB, 1, rank_row))
  values <- as.vector(percentiles_AB)
  
  summary_stat = "median"
  # Return summary statistic
  if (summary_stat == "mean") {
    return(mean(values, na.rm = TRUE))
  } else {
    return(median(values, na.rm = TRUE))
  }
}


IDs = c("OT-MONDO_0004979-associated-targets-4_28_2025-v25_03_asthma.tsv","OT-MONDO_0004975-associated-targets-4_28_2025-v25_03_Alzheimer.tsv","OT-EFO_0000341-associated-targets-4_28_2025-v25_03_COPD.tsv",
        "OT-EFO_0000400-associated-targets-4_28_2025-v25_03_Diabetes.tsv","OT-MONDO_0007254-associated-targets-4_28_2025-v25_03_BRCA.tsv",
        "OT-EFO_1001949-associated-targets-4_28_2025-v25_03_COAD.tsv","OT-EFO_0000519-associated-targets-4_28_2025-v25_03_GBM.tsv","OT-EFO_0005543-associated-targets-4_28_2025-v25_03_LGG.tsv",
        "OT-MONDO_0008170-associated-targets-4_28_2025-v25_03_OV.tsv","OT-MONDO_0017276-associated-targets-4_28_2025-v25_03_FtD.tsv",
        "OT-MONDO_0017276-associated-targets-4_28_2025-v25_03_FtD.tsv","OT-MONDO_0017276-associated-targets-4_28_2025-v25_03_FtD.tsv")
index = 1
distance_1 = c()
meth = c()
diseases_all = c()

for (disease in c('Asthma','Alzheimer','COPD','Diabetes', 'GS-BRCA', 'GS-COAD', 'GS-GBM', 'GS-LGG', 'GS-OV', 'C9', 'MAPT', 'GRN')){
  print(disease)
  theta = read.csv(paste('../results/RFOnM/State_',disease,'.csv',sep = ""),sep=',',header = F)
  gene_name = read.csv(paste('../data/features_',disease,'_head.csv',sep = ""), sep=',')
  adj = read.csv(paste('../data/graph_',disease,'.csv',sep = ""), sep=',',header = F)
  disgnet = read.csv(file = paste('../data/',IDs[index],sep = ""),header = T, sep = "\t")
  disgnet = disgnet[order(disgnet$globalScore,decreasing = T),]
  disgnet = disgnet[1:100,]
  
  sin_values = sin(theta)
  cos_values = cos(theta)
  
  mean_zz = c()
  max_thre = c()
  for (uppers in upper_size){
    mean_z = c()
    for (thre in thres){
      counts = rowSums(sin_values>thre | cos_values>thre)
      Hc = which.min(abs(counts - uppers))
      disease_gene <- which((sin_values[Hc,] > thre) | (cos_values[Hc,] > thre), arr.ind = TRUE)
      g <- graph_from_adjacency_matrix(as.matrix(adj), mode = "undirected")
      subg <- induced_subgraph(g, disease_gene[,2])
      clu <- igraph::components(subg)
      largest_component <- which(clu$csize==max(clu$csize))
      node_ids_largest_component <- which(as.vector(clu$membership)%in%largest_component)
      node_ids_largest_component <- as.numeric(sub("V", "", names(clu$membership)[node_ids_largest_component]))
      mean_z = c(mean_z, max(sum(gene_name$v1[node_ids_largest_component])/length(node_ids_largest_component),  
                             sum(gene_name$v2[node_ids_largest_component])/length(node_ids_largest_component)))
    }
    mean_zz = c(mean_zz,mean(mean_z))
    max_thre = c(max_thre,thres[which.max(mean_z)])
  }
  
  optimal_size = upper_size[which.max(mean_zz)]
  thre = max_thre[which.max(mean_zz)]
  counts = rowSums(sin_values>thre | cos_values>thre)
  Hc = which.min(abs(counts - optimal_size))
  
  disease_gene <- which((sin_values[Hc,] > thre) | (cos_values[Hc,] > thre), arr.ind = TRUE)
  g <- graph_from_adjacency_matrix(as.matrix(adj), mode = "undirected")
  subg <- induced_subgraph(g, disease_gene[,2])
  clu <- igraph::components(subg)
  largest_component <- which(clu$csize==max(clu$csize))
  node_ids_largest_component <- which(as.vector(clu$membership)%in%largest_component)
  node_ids_largest_component <- as.numeric(sub("V", "", names(clu$membership)[node_ids_largest_component]))
  disease_gene <- rownames(gene_name)[node_ids_largest_component]
  
  disease_genes_expr <- order(-gene_name$v1)[1:length(disease_gene)]
  disease_genes_wgs <- order(-gene_name$v2)[1:length(disease_gene)]
  disease_genes_expr = rownames(gene_name)[disease_genes_expr]
  disease_genes_wgs = rownames(gene_name)[disease_genes_wgs]
  
  disease_robust0 <- read.csv(paste0('../results/robust/0_', disease, '.csv'), sep = ',', header = TRUE)
  disease_robust1 <- read.csv(paste0('../results/robust/1_', disease, '.csv'), sep = ',', header = TRUE)
  disease_robust0 <- disease_robust0$vertex+1
  disease_robust1 <- disease_robust1$vertex+1
  disease_robust0 = rownames(gene_name)[disease_robust0]
  disease_robust1 = rownames(gene_name)[disease_robust1]
  
  diease_diamond0 <- read.csv(paste0('../results/diamond/first_100_added_nodes_weight_', disease, '_0_1.txt'), sep = '\t', header = TRUE)
  diease_diamond1 <- read.csv(paste0('../results/diamond/first_100_added_nodes_weight_', disease, '_1_1.txt'), sep = '\t', header = TRUE)
  diease_diamond0 <- diease_diamond0$DIAMOnD_node+1
  diease_diamond1 <- diease_diamond1$DIAMOnD_node+1
  
  diease_diamond0 = rownames(gene_name)[diease_diamond0]
  diease_diamond1 = rownames(gene_name)[diease_diamond1]
  
  if (disease %in% c("GS-LGG")){
    diease_DOMINO0 <- read.csv(paste0('../results/DOMINO/String/seed_',disease,'_1/','seed_',disease,'_1/', 'modules.out'), header = FALSE,sep = ',',na.strings=c("","NA"))
  } else{
    diease_DOMINO0 <- read.csv(paste0('../results/DOMINO/String/seed_',disease,'_0/','seed_',disease,'_0/', 'modules.out'), header = FALSE,sep = ',',na.strings=c("","NA"))
  }
  if (disease %in% c("COPD","GS-COAD","GRN")){
    diease_DOMINO1 = diease_DOMINO0
  } else {
    diease_DOMINO1 <- read.csv(paste0('../results/DOMINO/String/seed_',disease,'_1/','seed_',disease,'_1/', 'modules.out'), header = FALSE,sep = ',',na.strings=c("","NA"))
  }
  
  
  diease_DOMINO0 = melt(as.matrix(diease_DOMINO0))
  diease_DOMINO0 = diease_DOMINO0[complete.cases(diease_DOMINO0),]
  diease_DOMINO0 = data.frame(V1 = diease_DOMINO0$value)
  diease_DOMINO0$V1 = gsub("\\[","",diease_DOMINO0$V1)
  diease_DOMINO0$V1 = gsub("\\]","",diease_DOMINO0$V1)
  diease_DOMINO0 = gsub(" ","",diease_DOMINO0$V1)
  
  diease_DOMINO1 = melt(as.matrix(diease_DOMINO1))
  diease_DOMINO1 = diease_DOMINO1[complete.cases(diease_DOMINO1),]
  diease_DOMINO1 = data.frame(V1 = diease_DOMINO1$value)
  diease_DOMINO1$V1 = gsub("\\[","",diease_DOMINO1$V1)
  diease_DOMINO1$V1 = gsub("\\]","",diease_DOMINO1$V1)
  diease_DOMINO1 = gsub(" ","",diease_DOMINO1$V1)
  
  
  otp_high_conf = disgnet$symbol
  g <- graph_from_adjacency_matrix(as.matrix(adj), mode = "undirected")
  V(g)$name <- rownames(gene_name)
  
  disease_robust0 <- disease_robust0[disease_robust0%in% V(g)$name]
  disease_robust1 <- disease_robust1[disease_robust1%in% V(g)$name]
  diease_diamond0 <- diease_diamond0[diease_diamond0%in% V(g)$name]
  disease_robust1 <- disease_robust1[disease_robust1%in% V(g)$name]
  diease_DOMINO0 <- diease_DOMINO0[diease_DOMINO0%in% V(g)$name]
  diease_DOMINO1 <- diease_DOMINO1[diease_DOMINO1%in% V(g)$name]
  
  otp_high_conf <- otp_high_conf[otp_high_conf%in%V(g)$name]
  
  observed_dist1 <- set_to_set_percentile_distance(g, disease_gene, otp_high_conf)
  observed_dist2 <- set_to_set_percentile_distance(g, disease_genes_expr, otp_high_conf)
  observed_dist3 <- set_to_set_percentile_distance(g, disease_genes_wgs, otp_high_conf)
  observed_dist4 <- set_to_set_percentile_distance(g, disease_robust0, otp_high_conf)
  observed_dist5 <- set_to_set_percentile_distance(g, disease_robust1, otp_high_conf)
  observed_dist6 <- set_to_set_percentile_distance(g, diease_diamond0, otp_high_conf)
  observed_dist7 <- set_to_set_percentile_distance(g, diease_diamond1, otp_high_conf)
  observed_dist8 <- set_to_set_percentile_distance(g, diease_DOMINO0, otp_high_conf)
  observed_dist9 <- set_to_set_percentile_distance(g, diease_DOMINO1, otp_high_conf)

  
  distance_1 = c(distance_1,observed_dist1)
  meth = c(meth,rep("RFOnM (with)",sum(disease_gene%in%disgnet$symbol)))
  diseases_all = c(diseases_all,rep(disease,sum(disease_gene%in%disgnet$symbol)))
  
  distance_1 = c(distance_1,observed_dist2)
  meth = c(meth,rep("Expression",sum(disease_genes_expr%in%disgnet$symbol)))
  diseases_all = c(diseases_all,rep(disease,sum(disease_genes_expr%in%disgnet$symbol)))
  
  distance_1 = c(distance_1,observed_dist3)
  meth = c(meth,rep("GWAS",sum(disease_genes_wgs%in%disgnet$symbol)))
  diseases_all = c(diseases_all,rep(disease,sum(disease_genes_wgs%in%disgnet$symbol)))
  
  distance_1 = c(distance_1,observed_dist4)
  meth = c(meth,rep("ROBUST (expression)",sum(disease_robust0%in%disgnet$symbol)))
  diseases_all = c(diseases_all,rep(disease,sum(disease_robust0%in%disgnet$symbol)))
  
  distance_1 = c(distance_1,observed_dist5)
  meth = c(meth,rep("ROBUST (GWAS)",sum(disease_robust1%in%disgnet$symbol)))
  diseases_all = c(diseases_all,rep(disease,sum(disease_robust1%in%disgnet$symbol)))
  
  distance_1 = c(distance_1,observed_dist6)
  meth = c(meth,rep("DIAMonD (expression)",sum(diease_diamond0%in%disgnet$symbol)))
  diseases_all = c(diseases_all,rep(disease,sum(diease_diamond0%in%disgnet$symbol)))
  
  distance_1 = c(distance_1,observed_dist7)
  meth = c(meth,rep("DIAMonD (GWAS)",sum(diease_diamond1%in%disgnet$symbol)))
  diseases_all = c(diseases_all,rep(disease,sum(diease_diamond1%in%disgnet$symbol)))
  
  distance_1 = c(distance_1,observed_dist8)
  meth = c(meth,rep("DOMINO (expression)",sum(diease_DOMINO0%in%disgnet$symbol)))
  diseases_all = c(diseases_all,rep(disease,sum(diease_DOMINO0%in%disgnet$symbol)))
  
  distance_1 = c(distance_1,observed_dist9)
  meth = c(meth,rep("DOMINO (GWAS)",sum(diease_DOMINO1%in%disgnet$symbol)))
  diseases_all = c(diseases_all,rep(disease,sum(diease_DOMINO1%in%disgnet$symbol)))
  
  
  index = index + 1
}

dat1 = data.frame(Acccuracy = distance_1, method = rep(c("RFOnM","Expression","GWAS","ROBUST (expression)","ROBUST (GWAS)",
                                                "DIAMonD (expression)","DIAMonD (GWAS)","DOMINO (expression)","DOMINO (GWAS)"),12),
                  dataset = rep(c("Asthma","Alzheimer","COPD","Diabetes","GS-BRCA","GS-COAD","GS-GBM","GS-LGG","GS-OV","C9","MAPT","GRN"),each=9)) 

dat1$method = factor(dat1$method,levels = c("RFOnM","Expression","GWAS","ROBUST (expression)","ROBUST (GWAS)",
                                            "DIAMonD (expression)","DIAMonD (GWAS)",
                                            "DOMINO (expression)","DOMINO (GWAS)"))
dat1$dataset = factor(dat1$dataset,levels = c("Asthma","Alzheimer","COPD","Diabetes","GS-BRCA","GS-COAD","GS-GBM","GS-LGG","GS-OV","C9","MAPT","GRN"))


g2 = ggplot(data = dat1,aes(method,Acccuracy,color=method))+geom_point(size=3)+
  scale_color_manual(values = c("#5E81ACFF","#8FA87AFF", "#BF616AFF", "#E7D202FF", "#7D5329FF", "#F49538FF", "#66CDAAFF", "#D070B9FF", "#98FB98FF", "#FCA3B7FF"))+
  theme_bw()+xlab("")+ylab("Distance")+facet_wrap(~dataset,scales = "free",nrow=2)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 10),
        axis.text.x = element_text(angle = 40, hjust=1),
        strip.background = element_blank(),
        legend.position = "none") 

################################################################################################################
library(clusterProfiler) # R 4.2.0
library(org.Hs.eg.db)
library(enrichplot)


index = 1
ascore = c()
meth = c()
diseases_all = c()
module_size = c()

for (disease in c('Asthma','Alzheimer','COPD','Diabetes', 'GS-BRCA', 'GS-COAD', 'GS-GBM', 'GS-LGG', 'GS-OV', 'C9', 'MAPT', 'GRN')){
  print(disease)
  theta = read.csv(paste('../results/RFOnM/State_',disease,'.csv',sep = ""),sep=',',header = F)
  gene_name = read.csv(paste('../data/features_',disease,'_head.csv',sep = ""), sep=',')
  adj = read.csv(paste('../data/graph_',disease,'.csv',sep = ""), sep=',',header = F)

  sin_values = sin(theta)
  cos_values = cos(theta)
  
  mean_zz = c()
  max_thre = c()
  for (uppers in upper_size){
    mean_z = c()
    for (thre in thres){
      counts = rowSums(sin_values>thre | cos_values>thre)
      Hc = which.min(abs(counts - uppers))
      disease_gene <- which((sin_values[Hc,] > thre) | (cos_values[Hc,] > thre), arr.ind = TRUE)
      g <- graph_from_adjacency_matrix(as.matrix(adj), mode = "undirected")
      subg <- induced_subgraph(g, disease_gene[,2])
      clu <- igraph::components(subg)
      largest_component <- which(clu$csize==max(clu$csize))
      node_ids_largest_component <- which(as.vector(clu$membership)%in%largest_component)
      node_ids_largest_component <- as.numeric(sub("V", "", names(clu$membership)[node_ids_largest_component]))
      mean_z = c(mean_z, max(sum(gene_name$v1[node_ids_largest_component])/length(node_ids_largest_component),  
                   sum(gene_name$v2[node_ids_largest_component])/length(node_ids_largest_component)))
    }
    mean_zz = c(mean_zz,mean(mean_z))
    max_thre = c(max_thre,thres[which.max(mean_z)])
  }
  
  optimal_size = upper_size[which.max(mean_zz)]
  thre = max_thre[which.max(mean_zz)]
  counts = rowSums(sin_values>thre | cos_values>thre)
  Hc = which.min(abs(counts - optimal_size))
  
  disease_gene <- which((sin_values[Hc,] > thre) | (cos_values[Hc,] > thre), arr.ind = TRUE)
  g <- graph_from_adjacency_matrix(as.matrix(adj), mode = "undirected")
  subg <- induced_subgraph(g, disease_gene[,2])
  clu <- igraph::components(subg)
  largest_component <- which(clu$csize==max(clu$csize))
  node_ids_largest_component <- which(as.vector(clu$membership)%in%largest_component)
  node_ids_largest_component <- as.numeric(sub("V", "", names(clu$membership)[node_ids_largest_component]))
  disease_gene <- rownames(gene_name)[node_ids_largest_component]
  
  module_size = c(module_size, length(disease_gene))

  disease_genes_expr <- order(-gene_name$v1)[1:length(disease_gene)]
  disease_genes_wgs <- order(-gene_name$v2)[1:length(disease_gene)]
  disease_genes_expr = rownames(gene_name)[disease_genes_expr]
  disease_genes_wgs = rownames(gene_name)[disease_genes_wgs]
  
  disease_robust0 <- read.csv(paste0('../results/robust/0_', disease, '.csv'), sep = ',', header = TRUE)
  disease_robust1 <- read.csv(paste0('../results/robust/1_', disease, '.csv'), sep = ',', header = TRUE)
  disease_robust0 <- disease_robust0$vertex+1
  disease_robust1 <- disease_robust1$vertex+1
  disease_robust0 = rownames(gene_name)[disease_robust0]
  disease_robust1 = rownames(gene_name)[disease_robust1]
  
  diease_diamond0 <- read.csv(paste0('../results/diamond/first_100_added_nodes_weight_', disease, '_0_1.txt'), sep = '\t', header = TRUE)
  diease_diamond1 <- read.csv(paste0('../results/diamond/first_100_added_nodes_weight_', disease, '_1_1.txt'), sep = '\t', header = TRUE)
  diease_diamond0 <- diease_diamond0$DIAMOnD_node+1
  diease_diamond1 <- diease_diamond1$DIAMOnD_node+1
  
  diease_diamond0 = rownames(gene_name)[diease_diamond0]
  diease_diamond1 = rownames(gene_name)[diease_diamond1]
  
  if (disease %in% c("GS-LGG")){
    diease_DOMINO0 <- read.csv(paste0('../results/DOMINO/String/seed_',disease,'_1/','seed_',disease,'_1/', 'modules.out'), header = FALSE,sep = ',',na.strings=c("","NA"))
  } else{
    diease_DOMINO0 <- read.csv(paste0('../results/DOMINO/String/seed_',disease,'_0/','seed_',disease,'_0/', 'modules.out'), header = FALSE,sep = ',',na.strings=c("","NA"))
  }
  if (disease %in% c("COPD","GS-COAD","GRN")){
    diease_DOMINO1 = diease_DOMINO0
  } else {
    diease_DOMINO1 <- read.csv(paste0('../results/DOMINO/String/seed_',disease,'_1/','seed_',disease,'_1/', 'modules.out'), header = FALSE,sep = ',',na.strings=c("","NA"))
  }
  
  
  diease_DOMINO0 = melt(as.matrix(diease_DOMINO0))
  diease_DOMINO0 = diease_DOMINO0[complete.cases(diease_DOMINO0),]
  diease_DOMINO0 = data.frame(V1 = diease_DOMINO0$value)
  diease_DOMINO0$V1 = gsub("\\[","",diease_DOMINO0$V1)
  diease_DOMINO0$V1 = gsub("\\]","",diease_DOMINO0$V1)
  diease_DOMINO0 = gsub(" ","",diease_DOMINO0$V1)
  
  diease_DOMINO1 = melt(as.matrix(diease_DOMINO1))
  diease_DOMINO1 = diease_DOMINO1[complete.cases(diease_DOMINO1),]
  diease_DOMINO1 = data.frame(V1 = diease_DOMINO1$value)
  diease_DOMINO1$V1 = gsub("\\[","",diease_DOMINO1$V1)
  diease_DOMINO1$V1 = gsub("\\]","",diease_DOMINO1$V1)
  diease_DOMINO1 = gsub(" ","",diease_DOMINO1$V1)

  gene_list_converted_rfon <- bitr(disease_gene, fromType = "SYMBOL", 
                                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  gene_list_converted_expr <- bitr(disease_genes_expr, fromType = "SYMBOL", 
                                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  gene_list_converted_wgs <- bitr(disease_genes_wgs, fromType = "SYMBOL", 
                                  toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  gene_list_converted_robust0 <- bitr(disease_robust0, fromType = "SYMBOL", 
                                      toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  gene_list_converted_robust1 <- bitr(disease_robust1, fromType = "SYMBOL", 
                                      toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  gene_list_converted_diamond0 <- bitr(diease_diamond0, fromType = "SYMBOL", 
                                      toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  gene_list_converted_diamond1 <- bitr(diease_diamond1, fromType = "SYMBOL", 
                                      toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  gene_list_converted_DOMINO0 <- bitr(diease_DOMINO0, fromType = "SYMBOL", 
                                       toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  gene_list_converted_DOMINO1 <- bitr(diease_DOMINO1, fromType = "SYMBOL", 
                                       toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  kegg_result <- enrichKEGG(gene = gene_list_converted_rfon$ENTREZID, 
                            organism = 'hsa', 
                            keyType = 'kegg',
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05)
  kegg_result_expr <- enrichKEGG(gene = gene_list_converted_expr$ENTREZID, 
                                 organism = 'hsa', 
                                 keyType = 'kegg',
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 0.05)
  kegg_result_wgs <- enrichKEGG(gene = gene_list_converted_wgs$ENTREZID, 
                                organism = 'hsa', 
                                keyType = 'kegg',
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05)
  kegg_result_robust0 <- enrichKEGG(gene = gene_list_converted_robust0$ENTREZID, 
                                    organism = 'hsa', 
                                    keyType = 'kegg',
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.05)
  kegg_result_robust1 <- enrichKEGG(gene = gene_list_converted_robust1$ENTREZID, 
                                    organism = 'hsa', 
                                    keyType = 'kegg',
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.05)
  
  kegg_result_DOMINO0 <- enrichKEGG(gene = gene_list_converted_DOMINO0$ENTREZID, 
                                    organism = 'hsa', 
                                    keyType = 'kegg',
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.05)
  kegg_result_DOMINO1 <- enrichKEGG(gene = gene_list_converted_DOMINO1$ENTREZID, 
                                    organism = 'hsa', 
                                    keyType = 'kegg',
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.05)
  
  kegg_result_diamond0 <- enrichKEGG(gene = gene_list_converted_diamond0$ENTREZID, 
                                    organism = 'hsa', 
                                    keyType = 'kegg',
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.05)
  kegg_result_diamond1 <- enrichKEGG(gene = gene_list_converted_diamond1$ENTREZID, 
                                    organism = 'hsa', 
                                    keyType = 'kegg',
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.05)
  
  ascore = c(ascore,kegg_result@result$p.adjust,kegg_result_expr@result$p.adjust,kegg_result_wgs@result$p.adjust,kegg_result_robust0@result$p.adjust,kegg_result_robust1@result$p.adjust,
             kegg_result_diamond0@result$p.adjust,kegg_result_diamond1@result$p.adjust,kegg_result_DOMINO0@result$p.adjust,kegg_result_DOMINO1@result$p.adjust)
  meth = c(meth,rep("RFOnM",length(kegg_result@result$p.adjust)),rep("Expression",length(kegg_result_expr@result$p.adjust)),rep("GWAS",length(kegg_result_wgs@result$p.adjust)),
           rep("ROBUST (expression)",length(kegg_result_robust0@result$p.adjust)),rep("ROBUST (GWAS)",length(kegg_result_robust1@result$p.adjust)),
           rep("DIAMonD (expression)",length(kegg_result_diamond0@result$p.adjust)),rep("DIAMonD (GWAS)",length(kegg_result_diamond1@result$p.adjust)),
           rep("DOMINO (expression)",length(kegg_result_DOMINO0@result$p.adjust)),rep("DOMINO (GWAS)",length(kegg_result_DOMINO1@result$p.adjust)))
  diseases_all = c(diseases_all,rep(disease,length(kegg_result@result$p.adjust)+length(kegg_result_expr@result$p.adjust)+length(kegg_result_wgs@result$p.adjust)+length(kegg_result_robust0@result$p.adjust)+length(kegg_result_robust1@result$p.adjust)+
                                    length(kegg_result_diamond0@result$p.adjust)+length(kegg_result_diamond1@result$p.adjust)+
                                      length(kegg_result_DOMINO0@result$p.adjust)+length(kegg_result_DOMINO1@result$p.adjust)))
  
  index = index + 1
}

dat2 = data.frame(pvalue=ascore,method=meth,dataset=diseases_all)
dat2$method = factor(dat2$method,levels = c("RFOnM","Expression","GWAS","ROBUST (expression)","ROBUST (GWAS)",
                                            "DIAMonD (expression)","DIAMonD (GWAS)",
                                            "DOMINO (expression)","DOMINO (GWAS)"))
dat2$dataset = factor(dat2$dataset,levels = c("Asthma","Alzheimer","COPD","Diabetes","GS-BRCA","GS-COAD","GS-GBM","GS-LGG","GS-OV","C9","MAPT","GRN"))

my_comparisons = list(c("RFOnM","DOMINO (expression)"),c("RFOnM","DOMINO (GWAS)"))

g3 = ggplot(data = dat2,aes(method,-log10(pvalue),fill=method))+
  stat_compare_means(comparisons = my_comparisons,size=2,method = "wilcox.test")+
  geom_boxplot(outlier.size = 0.2,size=0.2,alpha=0.8)+
  scale_fill_manual(values = c("#5E81ACFF", "#8FA87AFF", "#BF616AFF", "#E7D202FF", "#7D5329FF", "#F49538FF", "#66CDAAFF", "#D070B9FF", "#98FB98FF", "#FCA3B7FF"))+
  theme_bw()+xlab("")+ylab("-log10(P-adjusted)")+facet_wrap(~dataset,scales = "free",nrow=2)+
  scale_y_sqrt()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 10),
        axis.text.x = element_text(angle = 40, hjust=1),
        strip.background = element_blank(),
        legend.position = "none") 


ggsave(g1,file="../figures/LCC.pdf",width=14.5, height=7.6,scale = 0.8,dpi = 500)
ggsave(g2,file="../figures/Overlap.pdf",width=15, height=7.6,scale = 0.8,dpi = 500)
ggsave(g3,file="../figures/KEGG.pdf",width=15, height=7.7,scale = 0.8,dpi = 500)


  
