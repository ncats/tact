# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences
#

library(dplyr)
library(reshape2)
library(dendextend)
library(mclust)
library(cluster)
library(corrplot)
library(maptree)
library(ggplot2)
library(tidyr)
library(stringr)


# Ref: http://www.sthda.com/english/articles/28-hierarchical-clustering-essentials/91-comparing-dendrograms-essentials/
# Ref: https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
# Ref: https://stat.ethz.ch/R-manual/R-devel/library/cluster/html/agnes.object.html
# Ref: https://www.rdocumentation.org/packages/cluster/versions/2.0.6/topics/agnes.object
# Ref: https://www.rdocumentation.org/packages/maptree/versions/1.4-7/topics/kgs
# Ref: https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
# Ref: https://www.rdocumentation.org/packages/stats/versions/3.5.1/topics/prcomp
# Ref: https://stackoverflow.com/questions/17499013/how-do-i-make-a-list-of-data-frames
# Ref: https://www.rdocumentation.org/packages/maptree/versions/1.4-7/topics/kgs
# Ref: http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
# Ref: https://www.datamentor.io/r-programming/saving-plot/
# Ref: http://rfunction.com/archives/812
# Ref: https://stackoverflow.com/questions/10302364/how-to-control-font-size-in-png
# Ref: https://stackoverflow.com/questions/10156417/subscripts-in-plots-in-r
# Ref: https://stackoverflow.com/questions/14290364/heatmap-with-values-ggplot2
# Ref: https://uc-r.github.io/tidyr
# Ref: http://www.sthda.com/english/wiki/ggplot2-title-main-axis-and-legend-titles
# Ref: https://statisticsglobe.com/change-font-size-of-ggplot2-plot-in-r-axis-text-main-title-legend
# Ref: http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
# Ref: https://stackoverflow.com/questions/14942681/change-size-of-axes-title-and-labels-in-ggplot2
# Ref: https://stackoverflow.com/questions/14487188/increase-distance-between-text-and-title-on-the-y-axis
# Ref: https://stackoverflow.com/questions/35839408/r-dplyr-drop-multiple-columns
# Ref: https://stackoverflow.com/questions/26838005/putting-x-axis-at-top-of-ggplot2-chart
# Ref: https://stackoverflow.com/questions/21187603/replace-characters-from-a-column-of-a-data-frame-r
# Ref: https://cran.r-project.org/web/packages/dendroextras/dendroextras.pdf
# Ref: https://stat.ethz.ch/R-manual/R-patched/library/cluster/html/agnes.object.html
# Ref: Tanglegram analysis was introduced to me by Rajarshi Guha, PhD when working together on an earlier and unrelated study at NCATS.


clean.up <- function (D) {
  row.names(D) <- D$CellLine
  D <- D %>% dplyr::select(-CellLine)
  
  return (D)
}

replace.cl.names <- function (D) {
  D$CellLine <- str_replace_all(D$CellLine, '163-KH2-024', '163-KH2-021')
  
  return(D)
}




compute.Adjusted.RandIndex = function (REFERENCE_CLUSTERING, CLUSTERING) { 
  #
  # REFERENCE_CLUSTERING, CLUSTERING: dataframes of : origID, clusterID
  #
  
  #install.packages("mclust")
  
  REFERENCE_CLUSTERING = REFERENCE_CLUSTERING[order(REFERENCE_CLUSTERING$origID),]
  origID_order = REFERENCE_CLUSTERING$origID
  ORDERED_CLUSTERING = CLUSTERING[match(origID_order, CLUSTERING$origID),]
  
  return (adjustedRandIndex(REFERENCE_CLUSTERING$clusterID, ORDERED_CLUSTERING$clusterID))
  
}



perform.clustering <- function (M, nr.clusters, dist.method, clust.method) {
  M_dist = as.matrix(M %>% dist(method=dist.method))
  HC_M = agnes(M_dist, diss = TRUE, method = clust.method, keep.diss = TRUE)
  
  return (HC_M)
}

extract.clustering <- function (clustering, nr.clusters) {
  Clusters_M = cutree(as.hclust(clustering), k = nr.clusters)
  HCM_M <- as.data.frame (Clusters_M)
  HCM_M$origID <- rownames(HCM_M)
  colnames(HCM_M) <- c('clusterID', 'origID')
  HCM_M <- HCM_M %>% select (origID, clusterID)
  #HCM_M = cbind (clustering$order.lab, Clusters_M)
  #HCM_M = as.data.frame (HCM_M)
  #colnames (HCM_M) = c ('origID','clusterID')
  
  return (HCM_M)
}

compute.KelleyIndex <- function (clustering) {
  KGS <- kgs (clustering, clustering$diss)
  KGS <- as.data.frame(KGS)
  KGS$ClusterNr <- rownames(KGS)
  KGS <- subset (KGS, select = c('ClusterNr','KGS'))
  colnames(KGS) <- c('ClusterNr','Kelley.index')
  KGS$ClusterNr <- as.numeric(KGS$ClusterNr)
  
  return (KGS)
}


pairwise.Adjusted.RandIndex <- function (clusterings, names) {
  
  first = TRUE
  for (i in 1:length(clusterings)) {
    clustering.A <- clusterings[[i]]
    
    for (j in 1:length(clusterings)) {
      
      
      clustering.B <- clusterings[[j]]
      
      RI <-compute.Adjusted.RandIndex (REFERENCE_CLUSTERING = clustering.A, CLUSTERING = clustering.B)
      
      RI <- data.frame ("cluster.A" = names[i], "cluster.B" = names[j], "RI" = RI)
      
      if (first) {
        R <- RI
        first = FALSE
      }
      else {
        R <- rbind (R, RI)
      }
      
    }
  }
  
  return (R)
  
}

do.pca <- function (D) {
  D.pca <- prcomp(D, scale = TRUE, center = TRUE)
  plot(D.pca$sdev)
  
  return (D.pca)
}





D <- read.csv (file = '../input/transformed.sliding.window.input.tab', header = TRUE, sep = '\t', comment.char = '')
D <- D %>% dplyr::rename ('CellLine' = 'cell.line')
saveRDS(D, file = '../analysis/CellLines_SlidingWindow.Rds')


triple.slidingwindow <- readRDS ('../analysis/CellLines_SlidingWindow.Rds')
triple.ic50 <- readRDS('../analysis/CellLines_IC50.Rds')
triple.top <- readRDS('../analysis/CellLines_Top.Rds')
triple.middle <- readRDS('../analysis/CellLines_Middle.Rds')
triple.bottom <- readRDS('../analysis/CellLines_Bottom.Rds')

triple.slidingwindow$CellLine = as.character(triple.slidingwindow$CellLine)
triple.ic50$CellLine = as.character(triple.ic50$CellLine)
triple.top$CellLine = as.character(triple.top$CellLine)
triple.middle$CellLine = as.character(triple.middle$CellLine)
triple.bottom$CellLine = as.character(triple.bottom$CellLine)


# Ref: https://stackoverflow.com/questions/2261079/how-to-trim-leading-and-trailing-whitespace

triple.slidingwindow$CellLine = trimws(triple.slidingwindow$CellLine)
triple.ic50$CellLine = trimws(triple.ic50$CellLine)
triple.top$CellLine = trimws(triple.top$CellLine)
triple.middle$CellLine = trimws(triple.middle$CellLine)
triple.bottom$CellLine = trimws(triple.bottom$CellLine)

triple.slidingwindow = replace.cl.names(triple.slidingwindow)
triple.ic50 = replace.cl.names(triple.ic50)
triple.top = replace.cl.names(triple.top)
triple.middle = replace.cl.names(triple.middle)
triple.bottom = replace.cl.names(triple.bottom)




triple.ic50 <- clean.up (triple.ic50)
triple.top <- clean.up (triple.top)
triple.middle <- clean.up (triple.middle)
triple.bottom <- clean.up (triple.bottom)
triple.slidingwindow <- clean.up (triple.slidingwindow)


### PCA



triple.slidingwindow.pca <- do.pca(triple.slidingwindow)
# Based on scree-plot, #princ. components to retain = 4
triple.slidingwindow <- triple.slidingwindow.pca$x[,1:4]

triple.ic50.pca <- do.pca(triple.ic50)
# Based on scree-plot, #princ. components to retain = 5
triple.ic50 <- triple.ic50.pca$x[,1:5]




triple.top.pca <- do.pca(triple.top)
# Based on scree-plot, #princ. components to retain = 3
triple.top <- triple.top.pca$x[,1:3]


triple.middle.pca <- do.pca(triple.middle)
# Based on scree-plot, #princ. components to retain = 4
triple.middle <- triple.middle.pca$x[,1:4]


triple.bottom.pca <- do.pca(triple.bottom)
# Based on scree-plot, #princ. components to retain = 4
triple.bottom <- triple.bottom.pca$x[,1:4]


### Cluster number determined by Kelley-index

clustering.slidingwindow <- perform.clustering (triple.slidingwindow, dist.method = 'euclidean', clust.method = 'ward')
kelley.slidingwindow <- compute.KelleyIndex(clustering = clustering.slidingwindow)
ggplot(kelley.slidingwindow, aes (x = ClusterNr, y = Kelley.index)) + geom_point()
clustering.slidingwindow <- extract.clustering (clustering.slidingwindow, nr.clusters = 4)

clustering.ic50 <- perform.clustering (triple.ic50, dist.method = 'euclidean', clust.method = 'ward')
kelley.ic50 <- compute.KelleyIndex(clustering = clustering.ic50)
ggplot(kelley.ic50, aes (x = ClusterNr, y = Kelley.index)) + geom_point()
clustering.ic50 <- extract.clustering (clustering.ic50, nr.clusters = 4)


clustering.top <- perform.clustering (triple.top, dist.method = 'euclidean', clust.method = 'ward')
kelley.top <- compute.KelleyIndex(clustering = clustering.top)
ggplot(kelley.top, aes (x = ClusterNr, y = Kelley.index)) + geom_point()
clustering.top <- extract.clustering (clustering.top, nr.clusters = 5)



clustering.middle <- perform.clustering (triple.middle, dist.method = 'euclidean', clust.method = 'ward')
kelley.middle <- compute.KelleyIndex(clustering = clustering.middle)
ggplot(kelley.middle, aes (x = ClusterNr, y = Kelley.index)) + geom_point()
clustering.middle <- extract.clustering (clustering.middle, nr.clusters = 4)



clustering.bottom <- perform.clustering (triple.bottom, dist.method = 'euclidean', clust.method = 'ward')
kelley.bottom <- compute.KelleyIndex(clustering = clustering.bottom)
ggplot(kelley.bottom, aes (x = ClusterNr, y = Kelley.index)) + geom_point()
clustering.bottom <- extract.clustering (clustering.bottom, nr.clusters = 5)



slidingwindow.dend <- triple.slidingwindow %>% dist(method='euclidean') %>% hclust(method = "ward.D") %>% as.dendrogram
ic50.dend <- triple.ic50 %>% dist(method='euclidean') %>% hclust(method = "ward.D") %>% as.dendrogram
top.dend <- triple.top %>% dist(method='euclidean') %>% hclust(method = "ward.D") %>% as.dendrogram
middle.dend <- triple.middle %>% dist(method='euclidean') %>% hclust(method = "ward.D") %>% as.dendrogram
bottom.dend <- triple.bottom %>% dist(method='euclidean') %>% hclust(method = "ward.D") %>% as.dendrogram



#dendrograms <- dendlist("Sliding.Window" = slidingwindow.dend, "AC50" = ic50.dend, "Top" = top.dend, "Middle" = middle.dend, "Bottom" = bottom.dend)
# Top: low concentration
# Bottom: high-conentration

dendrograms <- dendlist("Sliding Window" = slidingwindow.dend, 'AC50' = ic50.dend, "Low Conc." = top.dend, "Middle Conc." = middle.dend, "High Conc." = bottom.dend)


corrplot(cor.dendlist(dendrograms), "pie", "lower")
corrplot(cor.dendlist(dendrograms), "number", "lower")


cor.dendlist(dendrograms, method = "cophenetic")
cor.dendlist(dendrograms, method = "baker")


png("../plots/corrplot_cophenetic_pie.png", 5000, 5000, pointsize=15, res = 600)
corrplot(cor.dendlist(dendrograms, method = "cophenetic"), "pie", "lower")
dev.off()


png("../plots/corrplot_cophenetic_number.png", 5000, 5000, pointsize=15, res = 600)
corrplot(cor.dendlist(dendrograms, method = "cophenetic"), "number", "lower")
dev.off()


png("../plots/corrplot_baker_pie.png", 5000, 5000, pointsize=15, res = 600)
corrplot(cor.dendlist(dendrograms, method = "baker"), "pie", "lower")
dev.off()


png("../plots/corrplot_baker_number.png", 5000, 5000, pointsize=15, res = 600)
corrplot(cor.dendlist(dendrograms, method = "baker"), "number", "lower")
dev.off()


df.all.entanglement <- data.frame()

# Top: Low concentration,bottom: high concentration (confisung naming, but this it's correct, has historical reasons in the project)
png("../plots/tanglegram_ic50_low_conc.png", 5000, 5000, pointsize=15, res = 600)
dendrograms.ic50.top <- dendlist("AC50" = ic50.dend, "Low Conc." = top.dend)
tanglegram(dendrograms.ic50.top, sort = TRUE, common_subtrees_color_lines = TRUE,
            highlight_distinct_edges  = TRUE, highlight_branches_lwd = FALSE, margin_inner = 8, margin_outer = 3,
           main_left = expression('AC'[50]), main_right = expression('Low Conc.'))
dev.off()
df.entanglement <- data.frame ('dendrogram_pair' = 'AC50_LowConc', 'entaglement' = entanglement(dendrograms.ic50.top, L = 2))
df.all.entanglement <- rbind (df.all.entanglement, df.entanglement)


png("../plots/tanglegram_ic50_middle_conc.png", 5000, 5000, pointsize=15, res = 600)
dendrograms.ic50.middle <- dendlist("AC50" = ic50.dend, "Middle Conc." = middle.dend)
tanglegram(dendrograms.ic50.middle, sort = TRUE, common_subtrees_color_lines = TRUE,
            highlight_distinct_edges  = TRUE, highlight_branches_lwd = FALSE, margin_inner = 8, margin_outer = 3,
           main_left = expression('AC'[50]), main_right = expression('Middle Conc.'))
dev.off()
df.entanglement <- data.frame ('dendrogram_pair' = 'AC50_MiddleConc', 'entaglement' = entanglement(dendrograms.ic50.middle, L = 2))
df.all.entanglement <- rbind (df.all.entanglement, df.entanglement)

png("../plots/tanglegram_ic50_high_conc.png", 5000, 5000, pointsize=15, res = 600)
#dendrograms.ic50.bottom <- dendlist("AC50" = ic50.dend, "Bottom" = bottom.dend)
dendrograms.ic50.bottom <- dendlist("AC50" = ic50.dend, "High Conc." = bottom.dend)
tanglegram(dendrograms.ic50.bottom, sort = TRUE, common_subtrees_color_lines = TRUE,
            highlight_distinct_edges  = TRUE, highlight_branches_lwd = FALSE, margin_inner = 8, margin_outer = 3,
           main_left = expression('AC'[50]), main_right = expression('High Conc.'))
dev.off()
df.entanglement <- data.frame ('dendrogram_pair' = 'AC50_HighConc', 'entaglement' = entanglement(dendrograms.ic50.bottom, L = 2))
df.all.entanglement <- rbind (df.all.entanglement, df.entanglement)

png("../plots/tanglegram_ic50_sw.png", 5000, 5000, pointsize=15, res = 600)
dendrograms.ic50.sw <- dendlist("AC50" = ic50.dend, "Sliding Window" = 
                                  slidingwindow.dend)
tanglegram(dendrograms.ic50.sw, sort = TRUE, common_subtrees_color_lines = TRUE,
            highlight_distinct_edges  = TRUE, highlight_branches_lwd = FALSE, margin_inner = 8, margin_outer = 3,
           main_left = expression('AC'[50]), main_right = expression('Sliding Window'))
dev.off()
df.entanglement <- data.frame ('dendrogram_pair' = 'AC50_SW', 'entaglement' = entanglement(dendrograms.ic50.sw, L = 2))
df.all.entanglement <- rbind (df.all.entanglement, df.entanglement)

df.all.entanglement

write.table(df.all.entanglement, file = '../analysis/entaglement_values.tab', sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
saveRDS (df.all.entanglement, file = '../analysis/entaglement_values.Rds')


all.clusterings <- list(clustering.ic50, clustering.top, clustering.middle, clustering.bottom, clustering.slidingwindow)
#a<- all.clusterings[[1]]
#a
#all.names <- c('AC50', 'Top', 'Middle', 'Bottom', 'SW')
all.names <- c('AC50', 'Low Conc.', 'Middle Conc.', 'High Conc.', 'Sliding Window')

RIs <- pairwise.Adjusted.RandIndex(clusterings = all.clusterings, names = all.names)
pRIs <- dcast (RIs, cluster.A ~ cluster.B, value.var = 'RI')
pRIs


#pRIs.long <- melt(pRIs, value.name = 'rand_index', id.vars = 'cluster.A') %>% rename ('cluster.B' = 'variable') %>% rename ('method1' = 'cluster.A', 'method2' = 'cluster.B')

#pRIs.long

### For heatmap-reordering
rownames(pRIs) <- pRIs$cluster.A
M <- pRIs %>% select (-cluster.A)
M <- as.matrix(M)

hm.dd <- as.dist (M)
hm.hc <- hclust(hm.dd)

M.reordered <- M[hm.hc$order, hm.hc$order]

get_lower_tri<-function(dist.matrix){
  dist.matrix[upper.tri(dist.matrix)] <- NA
  return(dist.matrix)
}

get_upper_tri <- function(dist.matrix){
  dist.matrix[lower.tri(dist.matrix)]<- NA
  return(dist.matrix)
}

M.lower.triangle <- get_lower_tri(M.reordered)
# Melt the correlation matrix
M.long <- melt(M.lower.triangle, na.rm = TRUE)
M.long <- M.long %>% rename ('method1' = 'Var1', 'method2' = 'Var2', 'RI' = 'value')
RIs.ordered.halfmartix <- M.long

###

#pRIs.long <- pRIs %>% gather (cluster.B, RI, AC50:SW) %>% rename ('method1' = 'cluster.A', 'method2' = 'cluster.B')
#RIs <- RIs %>% rename ('method1' = 'cluster.A', 'method2' = 'cluster.B')

saveRDS(RIs.ordered.halfmartix, file = '../analysis/rand_indices.Rds')
write.table(RIs.ordered.halfmartix, file = '../analysis/rand_indices.tab', sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)



png("../plots/hm_rand_indices.png", 5000, 5000, pointsize=15, res = 600)

ggheatmap <- ggplot(RIs.ordered.halfmartix, aes(method1, method2, fill = RI)) +
    geom_tile(color = 'white') +
    geom_text(aes(label = round(RI, 2)), size = 6) +
    scale_fill_gradient2(low = "white", high = "orange", 
                       midpoint = 0, limit = c(0,1), space = "Lab", 
                       name="Adjusted\nRand-Index") +
     xlab("Clustering Method 1") +
     ylab("Clustering Method 2") +
  
    theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 16, hjust = 1),
        axis.text.y = element_text(size = 16),
        #axis.title=element_text(),
        axis.title.x = element_text(margin = margin(t = 20, r = 20, b = 20, l = 20), size=18,face="bold"),
        axis.title.y = element_text(margin = margin(t = 20, r = 20, b = 20, l = 20), size=18,face="bold"),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.4, 0.85),
        legend.direction = "horizontal"
        
        )+
  scale_y_discrete(position = "right")+
  
  coord_fixed()
  
print(ggheatmap)


dev.off()
