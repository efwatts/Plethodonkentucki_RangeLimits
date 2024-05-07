library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")
library("hierfstat")
library("ade4")
library("apex")
library("mmod")
library("poppr")

# Creating DNAbin objects
k_data <- read.structure("kentuckionly_labeledindiv.str", n.ind=21, n.loc=6803, onerowperind=FALSE,
                         col.lab=1, col.pop=2, row.marknames=0,
                         NA.char="-9", ask=TRUE, quiet=FALSE)
k_data@pop <- as.factor(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","16","16","17","18","19"))



# Run DAPC
dapc_result <- dapc(k_data, n.pca = 4, n.da = 2, center = TRUE, scale = FALSE,
                    var.contrib = FALSE, var.loadings = FALSE, pca.info = TRUE,
                    pca.select = c("nbEig", "percVar"), perc.pca = 5)

# Visualize DAPC results with scatterplot
scatter(dapc_result, geom = "scatter", posi.leg = "topright", bandwidth = 0.1)

# Find clusters
grp <- find.clusters(k_data, max.n.clust=20)
# k = 2 is the best fit 

table(pop(k_data), grp$grp)
table.value(table(pop(k_data), grp$grp), col.lab=paste("inf", 1:6),
            row.lab=paste("ori", 1:6))


# DAPC
dapc1 <- dapc(k_data, grp$grp, n.pca = 4, n.da = 2)
dapc1
scatter(dapc1)
