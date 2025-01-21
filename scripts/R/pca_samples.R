library("fst")
library(tidyverse)
library(factoextra)
library(data.table)
library(ggbiplot)
library(ggbiplot)

#maybe just load dseq and plot matix
# args(dcast)
fst <- read_fst("resultsOverlapsStats_5mCpG_5hmCpG_sortedBed_minCov10.fst", as.data.table = TRUE)


dt_activeL1 <- dcast(fst[variable=="Freq_5mCpG.geom_mean",],  ... ~ samples, value.var = "value") 
dt_activeL1 <- dt_activeL1[,!c("id", "RepeatID", "variable")] 

dim(dt_activeL1)


setnafill(dt_activeL1, fill=0)
res.pca <- prcomp(transpose(dt_activeL1), scale = TRUE)
# res.pca <- prcomp(dt_activeL1, scale = TRUE)

##get proportion of variance explaiined
## importance of componenets  
summary(res.pca)
print(res.pca)

# plotPCA(vsd, intgroup=c("condition", "type"))

#install_github("vqv/ggbiplot")
args(ggbiplot)
g <- ggbiplot(res.pca,
              obs.scale = 1,
              var.scale = 1,
            #   groups = training$Species,
            # pc.biplot = FALSE,
              ellipse = TRUE,
              circle = TRUE,
              ellipse.prob = 0.68)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')

ggsave(g, filename = "pca_active_line1_samples.png")


fviz_pca_ind(iris.pca, label="none", habillage=iris$Species,
             addEllipses=TRUE, ellipse.level=0.95)