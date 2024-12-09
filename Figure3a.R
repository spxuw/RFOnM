library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(disgenet2r)
library(igraph)
library(gridExtra)
library(ggpubr)

setwd("/udd/spxuw/RFON/code")

for (disease in c('ICGC')){
  theta = read.csv(paste('../results/GEO_v4/State_',disease,'.csv',sep = ""),sep=',',header = F)
  gene_name = read.csv(paste('../data/GEO/features_',disease,'_head.csv',sep = ""), sep=',')
  adj = read.csv(paste('../data/GEO/graph_',disease,'.csv',sep = ""), sep=',',header = F)

  sin_values = sin(theta)
  cos_values = cos(theta)
  
  counts = rowSums(sin_values>0.9 | cos_values>0.9)
  Hc = which.min(abs(counts - 800))
  
  lcc_size = c()
  for (i in 1:nrow(sin_values)){
    disease_gene <- which((sin_values[i,] > 0.9) | (cos_values[i,] > 0.9), arr.ind = TRUE)
    # g <- graph_from_adjacency_matrix(as.matrix(adj), mode = "undirected")
    # subg <- induced_subgraph(g, disease_gene[,2])
    # clu <- components(subg)
    # largest_component <- which.max(clu$csize)
    # node_ids_largest_component <- which(clu$membership==largest_component)
    # node_ids_largest_component <- as.numeric(sub("V", "", names(node_ids_largest_component)))
    # node_ids_largest_component <- as.numeric(sub("V", "", node_ids_largest_component))
    # disease_gene <- rownames(gene_name)[node_ids_largest_component]
    lcc_size = c(lcc_size,length(disease_gene))
  }
  Hc = which.min(abs(lcc_size - 500))
  disease_gene <- which((sin_values[Hc,] > 0.9) | (cos_values[Hc,] > 0.9), arr.ind = TRUE)
  # g <- graph_from_adjacency_matrix(as.matrix(adj), mode = "undirected")
  # subg <- induced_subgraph(g, disease_gene[,2])
  # clu <- components(subg)
  # largest_component <- which.max(clu$csize)
  # node_ids_largest_component <- which(clu$membership==largest_component)
  # node_ids_largest_component <- as.numeric(sub("V", "", names(node_ids_largest_component)))
  # node_ids_largest_component <- as.numeric(sub("V", "", node_ids_largest_component))
  disease_gene <- rownames(gene_name)[disease_gene]
  
  gene_list_converted_rfon <- bitr(disease_gene, fromType = "SYMBOL", 
                                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)

  
  kegg_result <- enrichKEGG(gene = gene_list_converted_rfon$ENTREZID, 
                            organism = 'hsa', 
                            keyType = 'kegg',
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05)

}

# all_genes <- keys(org.Hs.eg.db, keytype = "SYMBOL")
# 
# ego <- enrichGO(gene = gene_list_converted_rfon$ENTREZID,
#                 OrgDb = org.Hs.eg.db,
#                 keyType = "ENTREZID",  # Make sure to set the correct keyType
#                 ont = "BP",
#                 pAdjustMethod = "BH",
#                 qvalueCutoff = 0.05,
#                 readable = TRUE)

custom_colors <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1")

# pdf("../figures/kegg_enrichment_results.pdf")
# dotplot(ego, showCategory = 10, color = custom_colors) +
#   scale_y_discrete(labels = names) + # Set custom labels for the y-axis
#   theme(axis.text.y = element_text(size = 8)) 
# dev.off()
# 
# edox <- setReadable(kegg_result, 'org.Hs.eg.db', 'ENTREZID')
# p1 <- cnetplot(edox,showCategory = 20)
reactome_result <- enrichPathway(
  gene = gene_list_converted_rfon$ENTREZID, 
  organism = 'human',  # Specify the organism
  pAdjustMethod = "BH", 
  pvalueCutoff = 0.05,
  readable = TRUE      # Converts the gene IDs to gene symbols if possible
)

filtered_results <- reactome_result@result[
  reactome_result@result$Count <= 25,]

reactome_result_rep = reactome_result
reactome_result_rep@result <- filtered_results

dp <- dotplot(reactome_result_rep,showCategory=20,  color="qvalue")+
  theme_bw()+coord_flip()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 6),
        axis.text.x = element_text(angle = 45, hjust=1,size = 6),
        axis.text.y = element_text(size = 6),
        strip.background = element_blank()) 
  
ggsave(dp,file="../figures/KEGG_COPD.pdf",width=12, height=5,scale = 0.82,dpi = 500)
