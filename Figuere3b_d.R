library(caret)
library(igraph)
library(gridExtra)
library(ggpubr)
library(dplyr)

setwd("/udd/spxuw/RFON/code")

thre = 0.9

disease = "COPD"

theta_COPD = read.csv(paste('../results/GEO_v4/State_',disease,'.csv',sep = ""),sep=',',header = F)
gene_name_COPD = read.csv(paste('../data/GEO/features_',disease,'_head.csv',sep = ""), sep=',')
adj = read.csv(paste('../data/GEO/graph_',disease,'.csv',sep = ""), sep=',',header = F)

sin_values = sin(theta_COPD)
cos_values = cos(theta_COPD)

counts = rowSums(sin_values>thre | cos_values>thre)
Hc = which.min(abs(counts - 800))

lcc_size = c()
for (i in 1:nrow(sin_values)){
  disease_gene <- which((sin_values[i,] > thre) | (cos_values[i,] > thre), arr.ind = TRUE)
  # g <- graph_from_adjacency_matrix(as.matrix(adj), mode = "undirected")
  # subg <- induced_subgraph(g, disease_gene[,2])
  # clu <- components(subg)
  # largest_component <- which.max(clu$csize)
  # node_ids_largest_component <- which(clu$membership==largest_component)
  # node_ids_largest_component <- as.numeric(sub("V", "", names(node_ids_largest_component)))
  # #disease_gene <- rownames(gene_name)[disease_gene[,2]]
  # node_ids_largest_component <- as.numeric(sub("V", "", node_ids_largest_component))
  # disease_gene <- rownames(gene_name)[node_ids_largest_component]
  lcc_size = c(lcc_size,length(disease_gene))
}
Hc = which.min(abs(lcc_size - 500))
disease_gene <- which((sin_values[Hc,] > thre) | (cos_values[Hc,] > thre), arr.ind = TRUE)
disease_gene_COPD <- rownames(gene_name_COPD)[disease_gene]

###########################################################################################
disease = "COPD_meta"

theta_ICGC = read.csv(paste('../results/GEO_v4/State_',disease,'.csv',sep = ""),sep=',',header = F)
gene_name_ICGC = read.csv(paste('../data/GEO/features_',disease,'_head.csv',sep = ""), sep=',')
adj = read.csv(paste('../data/GEO/graph_',disease,'.csv',sep = ""), sep=',',header = F)

sin_values = sin(theta_ICGC)
cos_values = cos(theta_ICGC)

counts = rowSums(sin_values>thre | cos_values>thre)
Hc = which.min(abs(counts - 800))

lcc_size = c()
for (i in 1:nrow(sin_values)){
  disease_gene <- which((sin_values[i,] > thre) | (cos_values[i,] > thre), arr.ind = TRUE)
  # g <- graph_from_adjacency_matrix(as.matrix(adj), mode = "undirected")
  # subg <- induced_subgraph(g, disease_gene[,2])
  # clu <- components(subg)
  # largest_component <- which.max(clu$csize)
  # node_ids_largest_component <- which(clu$membership==largest_component)
  # node_ids_largest_component <- as.numeric(sub("V", "", names(node_ids_largest_component)))
  # #disease_gene <- rownames(gene_name)[disease_gene[,2]]
  # node_ids_largest_component <- as.numeric(sub("V", "", node_ids_largest_component))
  # disease_gene <- rownames(gene_name)[node_ids_largest_component]
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
disease_gene_ICGC <- rownames(gene_name_ICGC)[disease_gene]

overlap_gene = intersect(rownames(gene_name_COPD),rownames(gene_name_ICGC))
gene_name_COPD_sub = gene_name_COPD[overlap_gene,]
gene_name_ICGC_sub = gene_name_ICGC[overlap_gene,]

unions = unique(c(disease_gene_COPD,disease_gene_ICGC))

dat1 = data.frame(z1=gene_name_COPD_sub$v1,z2=gene_name_COPD_sub$v2,z3=gene_name_ICGC_sub$v1,z4=gene_name_ICGC_sub$v2)
dat1['C1'] = 'Out-module'
rownames(dat1) = overlap_gene
dat1$C1[rownames(dat1)%in%disease_gene_COPD] = "In-module"


g5 = ggplot(data = dat1,aes(C1,z3,fill=C1)) +
  geom_violin(trim=FALSE,scale = "width",lwd=0.2)+
  geom_boxplot(width=0.1, fill="white",outlier.size = 0.5,outlier.shape = 21,lwd=0.2)+
  stat_compare_means() +
  scale_fill_manual(values = c("#004488FF", "#DDAA33FF")) +
  xlab("") + ylab("Z-score") +theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 10),
        axis.text.x = element_text(angle = 40, vjust = 0.5, hjust=0.5),
        strip.background = element_blank(),
        legend.position = "none")

g6 = ggplot(data = dat1,aes(C1,z4,fill=C1)) +
  geom_violin(trim=FALSE,scale = "width",lwd=0.2)+
  geom_boxplot(width=0.1, fill="white",outlier.size = 0.5,outlier.shape = 21,lwd=0.2)+
  stat_compare_means() +
  scale_fill_manual(values = c("#004488FF", "#DDAA33FF")) +
  xlab("") + ylab("Z-score") + theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 10),
        axis.text.x = element_text(angle = 40, vjust = 0.5, hjust=0.5),
        strip.background = element_blank(),
        legend.position = "none")


dat1_sub = dat1
dat2_sub = dat1
dat1_sub = dat1_sub[dat1_sub$z1<0,]
dat2_sub = dat2_sub[dat2_sub$z2<0,]

dat1_sub['maxs'] = pmax(dat1_sub$z3,dat1_sub$z4)


g7 = ggplot(data = dat1_sub,aes(C1,maxs,fill=C1)) +
  geom_violin(trim=FALSE,scale = "width",lwd=0.2)+
  geom_boxplot(width=0.1, fill="white",outlier.size = 0.5,outlier.shape = 21,lwd=0.2)+
  stat_compare_means() +
  scale_fill_manual(values = c("#004488FF", "#DDAA33FF")) +
  xlab("") + ylab("Z-score") +theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 10),
        axis.text.x = element_text(angle = 40, vjust = 0.5, hjust=0.5),
        strip.background = element_blank(),
        legend.position = "none")

# g8 = ggplot(data = dat2_sub,aes(C1,maxs,fill=C1)) +
#   geom_violin(trim=FALSE,scale = "width",lwd=0.2)+
#   geom_boxplot(width=0.1, fill="white",outlier.size = 0.5,outlier.shape = 21,lwd=0.2)+
#   stat_compare_means() +
#   scale_fill_manual(values = c("#004488FF", "#DDAA33FF")) +
#   xlab("") + ylab("Z-score") + theme_bw()+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.text = element_text(face = "bold.italic", size = 10),
#         axis.text.x = element_text(angle = 40, vjust = 0.5, hjust=0.5),
#         strip.background = element_blank(),
#         legend.position = "none")

p1 = ggarrange(g5,g6,g7,ncol = 3, nrow = 1,align = "hv",labels = c("a","b","c"))
ggsave(p1,file="../figures/copd_com.pdf",width=9, height=4.5,scale = 0.7,dpi = 500)

####################################################################################
edges = read.csv("../data/STRING/9606.protein.physical.links.v12.0.txt",sep = " ")
annotation = read.csv("../data/STRING/9606.protein.info.v12.0.txt",sep = "\t")
edges = edges[edges$combined_score>=400,]
edges$combined_score = 1

# convert ENSP to gene symbol
edges$protein1 = annotation$preferred_name[match(edges$protein1,annotation$X.string_protein_id)]
edges$protein2 = annotation$preferred_name[match(edges$protein2,annotation$X.string_protein_id)]

edges_temp = edges

edges_temp <- edges_temp %>%
  filter((protein1 %in% disease_gene_ICGC)) %>%
  filter((protein2 %in% disease_gene_ICGC))

write.table(edges_temp,file = "copd_module_list.csv",sep = ",")
