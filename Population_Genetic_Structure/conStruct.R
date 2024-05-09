library(fossil)
library(conStruct)
library(readr)

#formatting data
vignette(topic="format-data",package="conStruct")

#read in locality data
localities <- read.csv("kentucki_only_localities.csv")
pop.loc <- cbind.data.frame(localities$Longitude,localities$Latitude,localities$Locality)
colnames(pop.loc) <- c("long", "lat", "pop")
geodist <- earth.dist(pop.loc[c("long", "lat")], dist = FALSE)
geodist


#structure to construct
conStruct.final <- structure2conStruct(infile = "kentucki_Poplable.str",
                                      onerowperind = FALSE,
                                      start.loci = 3,
                                      start.samples = 1,
                                      missing.datum = -9,
                                      outfile = "kentuckiConStruct_final")

str(conStruct.final)

#read in coordinates
coords <- read_table2("kentucki.coord", col_names = FALSE)

DF <- data.frame(coords)
DF
matrix <- apply(as.matrix.noquote(DF),2,as.numeric)
matrix

#how to run a conStruct analysis
vignette(topic="run-conStruct",package="conStruct")

# run a conStruct analysis

#   you have to specify:
#       the number of layers (K)
#       the allele frequency data (freqs)
#       the geographic distance matrix (geoDist)
#       the sampling coordinates (coords)

my.run <- conStruct(spatial = TRUE, 
                    K = 4, 
                    freqs = conStruct.data,
                    geoDist = geodist, 
                    coords = matrix,
                    prefix = "spK4")

full.run <- conStruct(spatial = TRUE, 
                    K = 21, 
                    freqs = conStruct.data,
                    geoDist = geodist, 
                    coords = matrix,
                    prefix = "spK21")



#########actual construct run
K1 <- conStruct(spatial = TRUE, 
                    K = 1, 
                    freqs = conStruct.data,
                    geoDist = geodist, 
                    coords = matrix,
                n.chains = 10,
                n.iter = 50000,
                    prefix = "K1")

K2 <- conStruct(spatial = TRUE, 
                K = 2, 
                freqs = conStruct.final,
                geoDist = geodist, 
                coords = matrix,
                n.chains = 10,
                n.iter = 50000,
                prefix = "K2_final")

K3 <- conStruct(spatial = TRUE, 
                K = 3, 
                freqs = conStruct.data,
                geoDist = geodist, 
                coords = matrix,
                n.chains = 10,
                n.iter = 50000,
                prefix = "K3")

K4 <- conStruct(spatial = TRUE, 
                K = 4, 
                freqs = conStruct.data,
                geoDist = geodist, 
                coords = matrix,
                n.chains = 10,
                n.iter = 50000,
                prefix = "K4")

K5 <- conStruct(spatial = TRUE, 
                K = 5, 
                freqs = conStruct.data,
                geoDist = geodist, 
                coords = matrix,
                n.chains = 10,
                n.iter = 50000,
                prefix = "K5")

K6 <- conStruct(spatial = TRUE, 
                K = 6, 
                freqs = conStruct.data,
                geoDist = geodist, 
                coords = matrix,
                n.chains = 10,
                n.iter = 50000,
                prefix = "K6")

K7 <- conStruct(spatial = TRUE, 
                K = 7, 
                freqs = conStruct.data,
                geoDist = geodist, 
                coords = matrix,
                n.chains = 10,
                n.iter = 50000,
                prefix = "K7")

K8 <- conStruct(spatial = TRUE, 
                K = 8, 
                freqs = conStruct.data,
                geoDist = geodist, 
                coords = matrix,
                n.chains = 10,
                n.iter = 50000,
                prefix = "K8")

K9 <- conStruct(spatial = TRUE, 
                K = 9, 
                freqs = conStruct.data,
                geoDist = geodist, 
                coords = matrix,
                n.chains = 10,
                n.iter = 50000,
                prefix = "K9")

K10 <- conStruct(spatial = TRUE, 
                K = 10, 
                freqs = conStruct.data,
                geoDist = geodist, 
                coords = matrix,
                n.chains = 10,
                n.iter = 50000,
                prefix = "K10")

K11 <- conStruct(spatial = TRUE, 
                K = 11, 
                freqs = conStruct.data,
                geoDist = geodist, 
                coords = matrix,
                n.chains = 10,
                n.iter = 50000,
                prefix = "K11")

K12 <- conStruct(spatial = TRUE, 
                K = 12, 
                freqs = conStruct.data,
                geoDist = geodist, 
                coords = matrix,
                n.chains = 10,
                n.iter = 50000,
                prefix = "K12")

K13 <- conStruct(spatial = TRUE, 
                K = 13, 
                freqs = conStruct.data,
                geoDist = geodist, 
                coords = matrix,
                n.chains = 10,
                n.iter = 50000,
                prefix = "K13")

K14 <- conStruct(spatial = TRUE, 
                K = 14, 
                freqs = conStruct.data,
                geoDist = geodist, 
                coords = matrix,
                n.chains = 10,
                n.iter = 50000,
                prefix = "K14")

K15 <- conStruct(spatial = TRUE, 
                K = 15, 
                freqs = conStruct.data,
                geoDist = geodist, 
                coords = matrix,
                n.chains = 10,
                n.iter = 50000,
                prefix = "K15")

K16 <- conStruct(spatial = TRUE, 
                K = 16, 
                freqs = conStruct.data,
                geoDist = geodist, 
                coords = matrix,
                n.chains = 10,
                n.iter = 50000,
                prefix = "K16")

K17 <- conStruct(spatial = TRUE, 
                K = 17, 
                freqs = conStruct.data,
                geoDist = geodist, 
                coords = matrix,
                n.chains = 10,
                n.iter = 50000,
                prefix = "K17")

K18 <- conStruct(spatial = TRUE, 
                K = 18, 
                freqs = conStruct.data,
                geoDist = geodist, 
                coords = matrix,
                n.chains = 10,
                n.iter = 50000,
                prefix = "K18")

K19 <- conStruct(spatial = TRUE, 
                 K = 19, 
                 freqs = conStruct.data,
                 geoDist = geodist, 
                 coords = matrix,
                 n.chains = 10,
                 n.iter = 50000,
                 prefix = "K19")



# how to visualize the output of a conStruct model
vignette(topic="visualize-results",package="conStruct")

load("K2_conStruct.results.Robj")

# assign the MAP admixture proportions from 
#   the first MCMC chain to a variable 
#   with a new name

admix.props <- conStruct.results$chain_1$MAP$admix.proportions

# make a STRUCTURE plot using the 
#   maximum a posteriori (MAP) estimates
#   from the first chain of a conStruct run

make.structure.plot(admix.proportions = admix.props)

# order by membership in layer 1
make.structure.plot(admix.proportions = admix.props,
                    sort.by = 1,
                    sample.names = row.names(pop.loc))

# re-order the stacking order of the layers
make.structure.plot(admix.proportions = admix.props,
                    layer.order = c(2,1,3,4),
                    sort.by = 2)
# add sample names
make.structure.plot(admix.proportions = admix.props,
                    sample.names = row.names(coords),
                    mar = c(4.5,4,2,2))


##compare funs
# load output files from a run with 
#   the spatial model and K=4
load("spK4_conStruct.results.Robj")
load("spK4_data.block.Robj")

# assign to new variable names
spK4_cr <- conStruct.results
spK4_db <- data.block

# load output files from a run with 
#   the spatial model and K=21
load("spK21_conStruct.results.Robj")
load("spK21_data.block.Robj")

# assign to new variable names
spK21_cr <- conStruct.results
spK21_db <- data.block

# compare the two runs
compare.two.runs(conStruct.results1=spK4_cr,
                 data.block1=spK4_db,
                 conStruct.results2=spK21_cr,
                 data.block2=spK21_db,
                 prefix="spK21_vs_spK4")

# generates a bunch of pdf figures

# how to compare and select between different conStruct models
vignette(topic="model-comparison",package="conStruct")


my.xvals <- x.validation(train.prop = 0.9,
                         n.reps = 10,
                         K = 1:10,
                         freqs = conStruct.kentuck,
                         data.partitions = NULL,
                         geoDist = geodist,
                         coords = matrix,
                         prefix = "fullrun",
                         n.iter = 1e3,
                         make.figs = TRUE,
                         save.files = FALSE,
                         parallel = FALSE,
                         n.nodes = NULL)

sp.results <- as.matrix(
  read.table("fullrun_sp_xval_results.txt",
             header = TRUE,
             stringsAsFactors = FALSE)
)
nsp.results <- as.matrix(
  read.table("fullrun_nsp_xval_results.txt",
             header = TRUE,
             stringsAsFactors = FALSE)
)

# or, format results from the output list
sp.results <- Reduce("cbind",lapply(my.xvals,function(x){unlist(x$sp)}),init=NULL)
nsp.results <- Reduce("cbind",lapply(my.xvals,function(x){unlist(x$nsp)}),init=NULL)

# first, get the 95% confidence intervals for the spatial and nonspatial
#   models over values of K (mean +/- 1.96 the standard error)

sp.CIs <- apply(sp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
nsp.CIs <- apply(nsp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})

# then, plot cross-validation results for K=1:3 with 8 replicates

par(mfrow=c(1,2))
plot(rowMeans(sp.results),
     pch=19,col="blue",
     ylab="predictive accuracy",xlab="values of K",
     ylim=range(sp.results,nsp.results),
     main="cross-validation results")
points(rowMeans(nsp.results),col="green",pch=19)

# finally, visualize results for the spatial model
#   separately with its confidence interval bars
#
# note that you could do the same with the spatial model, 
#   but the confidence intervals don't really show up 
#   because the differences between predictive accuracies
#   across values of K are so large.

plot(rowMeans(sp.results),
     pch=19,col="blue",
     ylab="predictive accuracy",xlab="values of K",
     ylim=range(sp.CIs),
     main="spatial cross-validation results")
segments(x0 = 1:nrow(sp.results),
         y0 = sp.CIs[1,],
         x1 = 1:nrow(sp.results),
         y1 = sp.CIs[2,],
         col = "blue",lwd=2)

###Reordered for publication 
load("K2_conStruct.results.Robj")

admix.props <- conStruct.results$chain_1$MAP$admix.proportions

reordered.names = c("3","2","4","5","10","11", 
                    "12","6","7","1","8",
                    "13","14","14","14",
                    "15","16","17","18","9",
                    "19")
reorder = c(10,2,1,3,4,
            8,9,11,20,5,
            6,7,12,13,14,15,
            16,17,18,19,21)

make.structure.plot(admix.props, 
                    sample.order = reorder,
                    sample.names = reordered.names)
