library(igraph)
library(ggplot2)

setwd("/udd/spxuw/RFON/code")
thre = 0.9

results <- data.frame()

for (disease in c('Asthma','Alzheimer','COPD','Diabetes')){
  theta = read.csv(paste('../results/GEO_v4/State_',disease,'.csv',sep = ""),sep=',',header = F)
  gene_name = read.csv(paste('../data/GEO/features_',disease,'_head.csv',sep = ""), sep=',')
  adj = read.csv(paste('../data/GEO/graph_',disease,'.csv',sep = ""), sep=',',header = F)

  sin_values = sin(theta)
  cos_values = cos(theta)
  
  counts = rowSums(sin_values>thre | cos_values>thre)
  Hc = which.min(abs(counts - 800))
  
  lcc_size = c()
  for (i in 1:nrow(sin_values)){
    disease_gene <- which((sin_values[i,] > thre) | (cos_values[i,] > thre), arr.ind = TRUE)
    lcc_size = c(lcc_size,length(disease_gene))
  }
  Hc = which.min(abs(lcc_size - 500))
  disease_gene <- which((sin_values[Hc,] > thre) | (cos_values[Hc,] > thre), arr.ind = TRUE)
  
  g <- graph_from_adjacency_matrix(as.matrix(adj), mode = "undirected")
  subg <- induced_subgraph(g, disease_gene[,2])
  clu <- components(subg)
  #grouped_nodes <- groups(clu)
  largest_component <- which.max(clu$csize)
  node_ids_largest_component <- which(clu$membership==largest_component)
  node_ids_largest_component <- as.numeric(sub("V", "", names(node_ids_largest_component)))
  disease_gene1 <- node_ids_largest_component
  
  disease_genes_expr <- order(-gene_name$v1)[1:length(disease_gene)]
  disease_genes_wgs <- order(-gene_name$v2)[1:length(disease_gene)]
  
  disease_robust0 <- read.csv(paste0('../results/robust/0_', disease, '.csv'), sep = ',', header = TRUE)
  disease_robust1 <- read.csv(paste0('../results/robust/1_', disease, '.csv'), sep = ',', header = TRUE)
  disease_robust0 <- disease_robust0$vertex
  disease_robust1 <- disease_robust1$vertex
  
  for (evaluations in c("Meandegree","Zscore")){
    if (evaluations == "Zscore"){
      results1 = read.table(file="../results/GEO_v4/LCC.csv",sep=",",row.names = 1,header = T)
      results1$evaluation=evaluations
      results <- rbind(results, results1)

    }
    if (evaluations == "Meandegree"){
      G <- graph_from_adjacency_matrix(as.matrix(adj), mode='undirected')
      subg <- induced_subgraph(G, paste("V",disease_gene,sep = ""))
      degrees <- degree(subg)
      mean_degree <- mean(degrees)
      results <- rbind(results, data.frame(Dataset = disease, Method = "RFOnM", ZScore = mean_degree,evaluation=evaluations))

      # disease_genes_expr ========================================================
      subg <- induced_subgraph(G, paste("V",disease_genes_expr,sep = ""))
      degrees <- degree(subg)
      mean_degree <- mean(degrees)
      results <- rbind(results, data.frame(Dataset = disease, Method = "Expression", ZScore = mean_degree,evaluation=evaluations))
      
      # disease_genes_gws ========================================================
      subg <- induced_subgraph(G, paste("V",disease_genes_wgs,sep = ""))
      degrees <- degree(subg)
      mean_degree <- mean(degrees)
      results <- rbind(results, data.frame(Dataset = disease, Method = "GWAS", ZScore = mean_degree,evaluation=evaluations))
      
      # disease_robust0 ========================================================
      subg <- induced_subgraph(G, paste("V",disease_robust0,sep = ""))
      degrees <- degree(subg)
      mean_degree <- mean(degrees)
      results <- rbind(results, data.frame(Dataset = disease, Method = "Robust (expression)", ZScore = mean_degree,evaluation=evaluations))
      
      # disease_robust1 ========================================================                                                    
      subg <- induced_subgraph(G, paste("V",disease_robust1,sep = ""))
      degrees <- degree(subg)
      mean_degree <- mean(degrees)
      results <- rbind(results, data.frame(Dataset = disease, Method = "Robust (GWAS)", ZScore = mean_degree,evaluation=evaluations))
    }
  }
}

results$Method = factor(results$Method,levels = c("RFOnM","Expression","GWAS","Robust (expression)","Robust (GWAS)"))
write.csv(results,file = "../results/GEO_v4/Structural.csv")

results$facet_label <- paste(results$evaluation, results$Dataset, sep = " | ")


g1 = ggplot(results, aes(x = Method, y = ZScore, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("#6497B1FF", "#6A359CFF", "#FFB04FFF", "#679C35FF", "#CD1076FF"))+
  theme_bw()+xlab("")+ylab("Mean degree")+  facet_wrap(~facet_label, scales = "free", nrow = 2)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 10),
        axis.text.x = element_text(angle = 40, hjust=1),
        strip.background = element_blank(),
        legend.position = "none") 
#####################################################################################################################################
library(caret)
library(igraph)
library(gridExtra)
library(ggpubr)


IDs = c("C0004096_disease_gda_summary_Asthma.tsv","C0002395_disease_gda_summary_Alzheimer.tsv","C0024117_disease_gda_summary_COPD.tsv","C0011849_disease_gda_summary_Diabetes.tsv")
index = 1
F_1 = c()
ascore = c()
meth = c()
diseases_all = c()

for (disease in c('Asthma','Alzheimer','COPD','Diabetes')){
  theta = read.csv(paste('../results/GEO_v4/State_',disease,'.csv',sep = ""),sep=',',header = F)
  gene_name = read.csv(paste('../data/GEO/features_',disease,'_head.csv',sep = ""), sep=',')
  adj = read.csv(paste('../data/GEO/graph_',disease,'.csv',sep = ""), sep=',',header = F)
  disgnet = read.csv(paste('../data/DisGENT/',IDs[index],sep = ""), sep='\t')
  
  sin_values = sin(theta)
  cos_values = cos(theta)
  
  counts = rowSums(sin_values>thre | cos_values>thre)
  Hc = which.min(abs(counts - 800))
  
  lcc_size = c()
  for (i in 1:nrow(sin_values)){
    disease_gene <- which((sin_values[i,] > thre) | (cos_values[i,] > thre), arr.ind = TRUE)
    lcc_size = c(lcc_size,length(disease_gene))
  }
  Hc = which.min(abs(lcc_size - 500))
  disease_gene <- which((sin_values[Hc,] > thre) | (cos_values[Hc,] > thre), arr.ind = TRUE)
  # g <- graph_from_adjacency_matrix(as.matrix(adj), mode = "undirected")
  # subg <- induced_subgraph(g, disease_gene[,2])
  # clu <- components(subg)
  # largest_component <- which.max(clu$csize)
  # node_ids_largest_component <- which(clu$membership==largest_component)
  # node_ids_largest_component <- as.numeric(sub("V", "", names(node_ids_largest_component)))
  disease_gene <- rownames(gene_name)[disease_gene]
  
  F_1 = c(F_1,sum(disease_gene%in%disgnet$Gene)/length(disease_gene))
  ascore = c(ascore,disgnet$Score_gda[match(disease_gene[disease_gene%in%disgnet$Gene],disgnet$Gene)])
  meth = c(meth,rep("RFOnM",sum(disease_gene%in%disgnet$Gene)))
  diseases_all = c(diseases_all,rep(disease,sum(disease_gene%in%disgnet$Gene)))
  
  disease_genes_expr <- order(-gene_name$v1)[1:length(disease_gene)]
  disease_genes_wgs <- order(-gene_name$v2)[1:length(disease_gene)]
  disease_genes_expr = rownames(gene_name)[disease_genes_expr]
  disease_genes_wgs = rownames(gene_name)[disease_genes_wgs]
  
  disease_robust0 <- read.csv(paste0('../results/robust/0_', disease, '.csv'), sep = ',', header = TRUE)
  disease_robust1 <- read.csv(paste0('../results/robust/1_', disease, '.csv'), sep = ',', header = TRUE)
  disease_robust0 <- disease_robust0$vertex
  disease_robust1 <- disease_robust1$vertex
  disease_robust0 = rownames(gene_name)[disease_robust0]
  disease_robust1 = rownames(gene_name)[disease_robust1]
  
  
  F_1 = c(F_1,sum(disease_genes_expr%in%disgnet$Gene)/length(disease_genes_expr))
  ascore = c(ascore,disgnet$Score_gda[match(disease_genes_expr[disease_genes_expr%in%disgnet$Gene],disgnet$Gene)])
  meth = c(meth,rep("Expression",sum(disease_genes_expr%in%disgnet$Gene)))
  diseases_all = c(diseases_all,rep(disease,sum(disease_genes_expr%in%disgnet$Gene)))
  
  F_1 = c(F_1,sum(disease_genes_wgs%in%disgnet$Gene)/length(disease_genes_wgs))
  ascore = c(ascore,disgnet$Score_gda[match(disease_genes_wgs[disease_genes_wgs%in%disgnet$Gene],disgnet$Gene)])
  meth = c(meth,rep("GWAS",sum(disease_genes_wgs%in%disgnet$Gene)))
  diseases_all = c(diseases_all,rep(disease,sum(disease_genes_wgs%in%disgnet$Gene)))
  
  F_1 = c(F_1,sum(disease_robust0%in%disgnet$Gene)/length(disease_robust0))
  ascore = c(ascore,disgnet$Score_gda[match(disease_robust0[disease_robust0%in%disgnet$Gene],disgnet$Gene)])
  meth = c(meth,rep("Robust (expression)",sum(disease_robust0%in%disgnet$Gene)))
  diseases_all = c(diseases_all,rep(disease,sum(disease_robust0%in%disgnet$Gene)))
  
  F_1 = c(F_1,sum(disease_robust1%in%disgnet$Gene)/length(disease_robust1))
  ascore = c(ascore,disgnet$Score_gda[match(disease_robust1[disease_robust1%in%disgnet$Gene],disgnet$Gene)])
  meth = c(meth,rep("Robust (GWAS)",sum(disease_robust1%in%disgnet$Gene)))
  diseases_all = c(diseases_all,rep(disease,sum(disease_robust1%in%disgnet$Gene)))
  
  index = index + 1
}

dat1 = data.frame(Acccuracy = F_1, method = rep(c("RFOnM","Expression","GWAS","Robust (expression)","Robust (GWAS)"),4),
                  dataset = rep(c('Asthma','Alzheimer','COPD','Diabetes'),each=5)) 

dat1$method = factor(dat1$method,levels = c("RFOnM","Expression","GWAS","Robust (expression)","Robust (GWAS)"))

dat2 = data.frame(score=ascore,method=meth,dataset=diseases_all)
dat2$method = factor(dat2$method,levels = c("RFOnM","Expression","GWAS","Robust (expression)","Robust (GWAS)"))

write.csv(dat1,file = "../results/GEO_v4/disgnet1.csv")
write.csv(dat2,file = "../results/GEO_v4/disgnet2.csv")



g2 = ggplot(data = dat1,aes(method,Acccuracy,fill=method))+geom_bar(stat="identity")+
  scale_fill_manual(values = c("#6497B1FF", "#6A359CFF", "#FFB04FFF", "#679C35FF", "#CD1076FF"))+
  theme_bw()+xlab("")+ylab("Overlap ratio")+facet_wrap(~dataset,scales = "free",nrow=1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 10),
        axis.text.x = element_text(angle = 40, hjust=1),
        strip.background = element_blank(),
        legend.position = "none") 

g3 = ggplot(data = dat2,aes(method,(score),color=method))+geom_jitter(size=1.5,alpha=0.7)+
  scale_color_manual(values = c("#6497B1FF", "#6A359CFF", "#FFB04FFF", "#679C35FF", "#CD1076FF"))+
  theme_bw()+xlab("")+ylab("GDA score")+facet_wrap(~dataset,scales = "free",nrow=1)+
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

IDs = c("C0004096_disease_gda_summary_Asthma.tsv","C0002395_disease_gda_summary_Alzheimer.tsv","C0024117_disease_gda_summary_COPD.tsv","C0011849_disease_gda_summary_Diabetes.tsv")
index = 1
ascore = c()
meth = c()
diseases_all = c()

for (disease in c('Asthma','Alzheimer','COPD','Diabetes')){
  theta = read.csv(paste('../results/GEO_v4/State_',disease,'.csv',sep = ""),sep=',',header = F)
  gene_name = read.csv(paste('../data/GEO/features_',disease,'_head.csv',sep = ""), sep=',')
  adj = read.csv(paste('../data/GEO/graph_',disease,'.csv',sep = ""), sep=',',header = F)
  disgnet = read.csv(paste('../data/DisGENT/',IDs[index],sep = ""), sep='\t')
  
  sin_values = sin(theta)
  cos_values = cos(theta)
  
  counts = rowSums(sin_values>thre | cos_values>thre)
  Hc = which.min(abs(counts - 800))
  
  lcc_size = c()
  for (i in 1:nrow(sin_values)){
    disease_gene <- which((sin_values[i,] > thre) | (cos_values[i,] > thre), arr.ind = TRUE)
    lcc_size = c(lcc_size,length(disease_gene))
  }
  Hc = which.min(abs(lcc_size - 500))
  disease_gene <- which((sin_values[Hc,] > thre) | (cos_values[Hc,] > thre), arr.ind = TRUE)
  # g <- graph_from_adjacency_matrix(as.matrix(adj), mode = "undirected")
  # subg <- induced_subgraph(g, disease_gene[,2])
  # clu <- components(subg)
  # largest_component <- which.max(clu$csize)
  # node_ids_largest_component <- which(clu$membership==largest_component)
  # node_ids_largest_component <- as.numeric(sub("V", "", names(node_ids_largest_component)))
  disease_gene <- rownames(gene_name)[disease_gene]
  
  disease_genes_expr <- order(-gene_name$v1)[1:length(disease_gene)]
  disease_genes_wgs <- order(-gene_name$v2)[1:length(disease_gene)]
  disease_genes_expr = rownames(gene_name)[disease_genes_expr]
  disease_genes_wgs = rownames(gene_name)[disease_genes_wgs]
  
  disease_robust0 <- read.csv(paste0('../results/robust/0_', disease, '.csv'), sep = ',', header = TRUE)
  disease_robust1 <- read.csv(paste0('../results/robust/1_', disease, '.csv'), sep = ',', header = TRUE)
  disease_robust0 <- disease_robust0$vertex
  disease_robust1 <- disease_robust1$vertex
  disease_robust0 = rownames(gene_name)[disease_robust0]
  disease_robust1 = rownames(gene_name)[disease_robust1]
  
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
  
  ascore = c(ascore,kegg_result@result$p.adjust,kegg_result_expr@result$p.adjust,kegg_result_wgs@result$p.adjust,kegg_result_robust0@result$p.adjust,kegg_result_robust1@result$p.adjust)
  meth = c(meth,rep("RFOnM",length(kegg_result@result$p.adjust)),rep("Expression",length(kegg_result_expr@result$p.adjust)),rep("GWAS",length(kegg_result_wgs@result$p.adjust)),
           rep("Robust (expression)",length(kegg_result_robust0@result$p.adjust)),rep("Robust (GWAS)",length(kegg_result_robust1@result$p.adjust)))
  diseases_all = c(diseases_all,rep(disease,length(kegg_result@result$p.adjust)+length(kegg_result_expr@result$p.adjust)+length(kegg_result_wgs@result$p.adjust)+length(kegg_result_robust0@result$p.adjust)+length(kegg_result_robust1@result$p.adjust)))
  
  index = index + 1
}

write.csv(dat2,file = "../results/GEO_v4/KEGG.csv")

dat2 = data.frame(pvalue=ascore,method=meth,dataset=diseases_all)
dat2$method = factor(dat2$method,levels = c("RFOnM","Expression","GWAS","Robust (expression)","Robust (GWAS)"))

my_comparisons = list(c("RFOnM","Robust (GWAS)"),c("RFOnM","Robust (expression)"))

g4 = ggplot(data = dat2,aes(method,-log10(pvalue),fill=method))+geom_boxplot(outlier.size = 1,size=0.2)+
  scale_fill_manual(values = c("#6497B1FF", "#6A359CFF", "#FFB04FFF", "#679C35FF", "#CD1076FF"))+
  stat_compare_means(comparisons = my_comparisons,size=2,method = "t.test")+
  theme_bw()+xlab("")+ylab("-log10(P-adjusted)")+facet_wrap(~dataset,scales = "free",nrow=1)+
  scale_y_sqrt()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 10),
        axis.text.x = element_text(angle = 40, hjust=1),
        strip.background = element_blank(),
        legend.position = "none") 


p1 = ggarrange(g1,g2,g3,g4,ncol = 1, nrow = 4,align = "v",labels = c("a","b","c","d"),heights = c(0.95,0.5,0.5,0.5))
ggsave(p1,file="../figures/combined.pdf",width=8, height=15.8,scale = 0.8,dpi = 500)


  