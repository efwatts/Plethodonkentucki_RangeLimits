
## running bedassle

library(adegenet)
library(BEDASSLE)
library(dplyr)
library(ggplot2)
library(broom)
library(data.table)
library(vcfR)
library(fossil)
library(viridis)
library(patchwork)
library(sp)
library(raster)
library(ecodist)
#library(radiator)
library(wesanderson)

# set wd and read in locality data
setwd("~/Library/CloudStorage/GoogleDrive-efw24@vt.edu/My Drive/OUwork/Dissertation/Chapter2/BEDASSLE")
localities <- read.csv("kentucki_only_localities_19.csv")
pop.loc <- cbind.data.frame(localities$Longitude,localities$Latitude,localities$Locality)
pop.loc <- cbind.data.frame(unique(localities$Longitude),unique(localities$Latitude),unique(localities$Locality), unique(localities$Elevation..m.))
colnames(pop.loc) <- c("long", "lat", "pop", "elevation")
geodist <- earth.dist(pop.loc[c("long", "lat")], dist = FALSE)

### d ky

# read in .vcf, turn to genpop
#ky.vcf <- read.vcfR("populations.snps.vcf", convertNA=TRUE)
#ky.gen <- vcfR2genind(ky.vcf) 

#trying converting STRUCTURE file instead
ky.gen <- read.structure("kentuckionly_labeledpops.str", n.ind=21, n.loc=6803, onerowperind=FALSE,
                         col.lab=1, col.pop=2, row.marknames=0,
                         NA.char="-9", ask=TRUE, quiet=FALSE)
md <- propTyped(ky.gen, by=c("both")) # confirm we have expected amount of missing data
table(md)[1]/table(md)[2] # proportion of missing sites
table(ky.gen@tab)
#ky.pops <-  gsub( "_.*$", "", rownames(ky.gen@tab)) # add pops

ky.2pops <- as.factor(c("1","1","1","1","2","2","2","1","1","1","1","2","2","2","2","2","2","1","2","1","2")) #these are LEA pops
ky.19pops <- as.factor(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","14","14","15","16","17","18","19")) #these are each as individuals
ky.b <- genind2genpop(ky.gen, ky.19pops) # turn into gen pop

# convert to BEDASSLE format
ky.ac <- as.matrix(ky.b@tab)
del <- seq(2, ncol(ky.ac), 2) # sequence of integers to drop non-ref allele
ky.ac <- ky.ac[,-del] 

# create equally sized matrix for sample sizes
ky.n <- matrix(nrow=nrow(ky.ac), ncol=ncol(ky.ac))

# name our rows the same thing
rownames(ky.n) <- rownames(ky.ac)

# get sample size per population
sample.n <- table(ky.gen@pop) 

# turn this into a vector
sample.sizes <- as.vector(sample.n)

# populate each row of matrix with sample sizes for pops
for(i in 1:nrow(ky.n)){
  ky.n[i,] <- sample.sizes[i]
}
ky.n <- ky.n*2 # adjust to account for loss of one allele

# calculate pairwise Fst
ky.p.fst.all <- calculate.all.pairwise.Fst(ky.ac, ky.n)
ky.p.fst.all

# look at global Fst
ky.p.fst <- calculate.pairwise.Fst(ky.ac, ky.n)
ky.p.fst

# drop levels and calc distance
#not dropping pops
drop.pop <- c("13","14")
pop.loc.sat <- pop.loc[!pop.loc$pop %in% drop.pop,]
droplevels(pop.loc.sat)
ky.geo <- earth.dist(pop.loc[c("long", "lat")], dist = FALSE)

# standardize geographic distances
geo.effect <- sd(c(ky.geo)) # save effect size
ky.geo <- ky.geo/sd(c(ky.geo)) #trying without this to see if it helps things run

# turn to vectors
sat.dist <- as.vector(ky.geo)
sat.gen <- as.vector(ky.p.fst.all)

# make data frame with these variables
sat.df <- cbind.data.frame(sat.dist, sat.gen)
sat.df$species <- rep("plethodon_ky",nrow(sat.df))
colnames(sat.df) <- c("distance", "fst", "species")

# calculate environmental distance

worldclim <- getData("worldclim",var="bio",res=0.5,lon=localities$Longitude[1],lat=localities$Latitude[1])
proj <- as.character(worldclim[[2]]@crs)
worldclim <- worldclim[[c(1,3,7,8,9,12,15)]]
names(worldclim) <- c("temp","isotherm", "temp_ann_range", "mean_temp_wet", "mean_temp_dry", "precip", "precip_season")
sp1 <- SpatialPoints(pop.loc[,c('long', 'lat')], proj4string=CRS(proj))
values <- extract(worldclim, sp1)
loc.master <- cbind.data.frame(pop.loc, values)
loc.uniq <- cbind.data.frame(loc.master$pop, loc.master$elevation, loc.master$temp,
                             loc.master$isotherm, loc.master$temp_ann_range, loc.master$mean_temp_wet,
                             loc.master$mean_temp_dry, loc.master$precip, loc.master$precip_season)
colnames(loc.uniq) <- c("population", "elevation", "temp", "isotherm", "temp_ann_range", 
                        "mean_temp_wet", "mean_temp_dry", "precip", "precip_season")
rownames(loc.uniq) <- loc.uniq$population
loc.uniq <- as.data.frame(loc.uniq)
loc.uniq <- loc.uniq[,-which(names(loc.uniq) %in% c("population"))]
env.dist.all <- dist(loc.uniq, diag = TRUE, upper = TRUE)

# prep dist matrix
ky.env <- as.matrix(env.dist.all)
colnames(ky.env) <- NULL
rownames(ky.env) <- NULL

# standardize env distances
env.effect <- sd(c(ky.env)) # save effect size
ky.env <- ky.env/sd(c(ky.env)) #trying not standardizing here too


#########################################
#########################################
########################################
#######################################
#The value of delta may set off warnings, 
#so temporarily disable warnings. 
op <- options("warn") 
options(warn =-1)

#Testing regular MCMC
MCMC(   
  counts = ky.ac,
  sample_sizes = ky.n,
  D = ky.geo,  # geographic distances #this could be an issue (may need euclidean distances)
  E = ky.env,  # environmental distances #this only has 19 entries, since it used unique localities (fixing first)
  k = nrow(ky.ac), 
  loci = ncol(ky.ac),  # dimensions of the data
  #delta = 0.000000000238,  # a small, positive, number #trying even smaller, since some of these are very close localities 
  #aD_stp = 0.075,   # step sizes for the MCMC
  #aE_stp = 0.05,
  #a2_stp = 0.025,
  #thetas_stp = 0.2,
  #mu_stp = 0.35,
  #delta = 0.0000000007,  
  delta = 0.000001, #default delta for both analyses
  aD_stp = 0.01, 
  aE_stp = 0.1, 
  a2_stp = 0.03, 
  thetas_stp = 0.1,
  mu_stp = 0.3, 
  ngen = 10000,  
  printfreq = 100,
  savefreq = 10000,
  samplefreq = 100,     # record current state for posterior (2000)
  prefix = "kentucki_19_pops_test_",# filename prefix
  continue = FALSE)

## run below code on local machine
plot_all_acceptance_rates('kentucki_19_pops_test_MCMC_output1.Robj') #shown in viewer

#make a vector population names
pop_vector <- as.vector(rownames(ky.n))
plot_all_trace('kentucki_19_pops_test_MCMC_output1.Robj', percent.burnin = 0, thinning = 1, population.names = pop_vector)

posterior.predictive.sample(MCMC.output = "kentucki_19_pops_test_MCMC_output1.Robj", posterior.predictive.sample.size = 1000, output.file = "test_post1000", prefix = "test_")
plot_posterior_predictive_samples("test_test_post1000.Robj", save.figure = TRUE, figure.name = "test_r1_post1000_1.jpg")

#use a burn-in of the first 8 million
plot_all_trace('kentucki_19_pops_test_MCMC_output1.Robj', percent.burnin = 80, thinning = 1, population.names = pop_vector)

show(load("kentucki_19_pops_test_MCMC_output1.Robj"))

##testing BB model
MCMC_BB(
  counts = ky.ac,
  sample_sizes = ky.n,
  D = ky.geo,  # geographic distances #this could be an issue (may need euclidean distances)
  E = ky.env,  # environmental distances #this only has 19 entries, since it used unique localities (fixing first)
  k = nrow(ky.ac), 
  loci = ncol(ky.ac),  # dimensions of the data
  #delta = 0.000000000238,  # a small, positive, number #trying even smaller, since some of these are very close localities 
  #aD_stp = 0.075,   # step sizes for the MCMC
  #aE_stp = 0.05,
  #a2_stp = 0.025,
  #thetas_stp = 0.2,
  #mu_stp = 0.35,
  #delta = 0.0000000007,  
  delta = 0.000001, #default delta for both analyses
  aD_stp = 0.5, 
  aE_stp = 0.9, #increase
  a2_stp = 0.01, 
  thetas_stp = 0.3,
  phi_stp = 110, #increase
  mu_stp = 2, 
  ngen = 10000,  
  printfreq = 100,
  savefreq = 10000,
  samplefreq = 100,     # record current state for posterior (2000)
  prefix = "kentucki_19_pops_BBtest_",# filename prefix
  continue = FALSE
)

## run below code on local machine
plot_all_acceptance_rates('kentucki_19_pops_BBtest_MCMC_output1.Robj') #shown in viewer


#make a vector population names
pop_vector <- as.vector(rownames(ky.n))
plot_all_trace('kentucki_19_pops_BBtest_MCMC_output1.Robj', percent.burnin = 0, thinning = 1, population.names = pop_vector)

posterior.predictive.sample(MCMC.output = "kentucki_19_pops_BBtest_MCMC_output1.Robj", posterior.predictive.sample.size = 1000, output.file = "BBtest_post1000", prefix = "BB_test_ky1_")
plot_posterior_predictive_samples("BB_test_ky1_BBtest_post1000.Robj", save.figure = TRUE, figure.name = "BBtest_r1_post1000_1.jpg")

#use a burn-in of the first 8 million
plot_all_trace('kentucki_19_pops_BBtest_MCMC_output1.Robj', percent.burnin = 80, thinning = 1, population.names = pop_vector)

show(load("kentucki_19_pops_BBtest_MCMC_output1.Robj"))


###THE BB_MCMC was the better fit, but I do still need to reduce the environmental factors
# calculate environmental distance
MCMC_BB(
  counts = ky.ac,
  sample_sizes = ky.n,
  D = ky.geo,  # geographic distances #this could be an issue (may need euclidean distances)
  E = ky.env,  # environmental distances #this only has 19 entries, since it used unique localities (fixing first)
  k = nrow(ky.ac), 
  loci = ncol(ky.ac),  # dimensions of the data
  delta = 0.000001, #default delta for both analyses
  aD_stp = 0.5, 
  aE_stp = 0.9, #increase
  a2_stp = 0.01, 
  thetas_stp = 0.3,
  phi_stp = 110,
  mu_stp = 2, 
  ngen = 1000000,  #1 million
  printfreq = 100000,
  savefreq = 10000,
  samplefreq = 100,     # record current state for posterior (2000)
  prefix = "kentucki_19_pops_BB_",# filename prefix
  continue = FALSE
)

## run below code on local machine
plot_all_acceptance_rates('kentucki_19_pops_BB_MCMC_output1.Robj') #shown in viewer

#make a vector population names
pop_vector <- as.vector(rownames(ky.n))
plot_all_trace('kentucki_19_pops_BB_MCMC_output1.Robj', percent.burnin = 0, thinning = 1, population.names = pop_vector)

posterior.predictive.sample(MCMC.output = "kentucki_19_pops_BB_MCMC_output1.Robj", posterior.predictive.sample.size = 1000, output.file = "BB_post1000", prefix = "BB_ky19_")
plot_posterior_predictive_samples("BB_ky19_BB_post1000.Robj", save.figure = TRUE, figure.name = "BB_r1_post1000_1.jpg")

#use a burn-in of the first 8 million
plot_all_trace('kentucki_19_pops_BB_MCMC_output1.Robj', percent.burnin = 80, thinning = 1, population.names = pop_vector)

show(load("kentucki_19_pops_BB_MCMC_output1.Robj"))




#to correct the aE : aD ratios, multiply by sd(geom) and sd(genetic_fst)
#geom <- read.csv(geodist, row.names = 1, header = TRUE)
#geom <- as.matrix(geom)

#genetic_fst <- read.csv("population_FST.csv", row.names = 1, header = TRUE)
#genetic_fst <- as.matrix(genetic_fst)

geodist
ky.p.fst.all

#correct the aE : aD ratios
#1 is geology
geology_bray_curtis_ratio <- ((as.vector(aE[1,])* sd(geodist))/(aD))
#2 is soil pc1
soil_pca_pc1_ratio <- ((as.vector(aE[2,])* sd(geodist))/(aD)) #this one isn't working, and I'm not sure why
#3 is genetic fst
genetic_fst_ratio <- ((as.vector(aE[3,])* sd(geodist))/(aD * sd(ky.p.fst.all)))

#calculate output after 80% burn-in
library(coda)
all.param = as.mcmc(cbind(geology_bray_curtis_ratio[-(1:80000)], genetic_fst_ratio[-(1:80000)]))
s1 <- summary(all.param)
s2 <- as.data.frame(s1$statistics)
s2 <- s2[1:2,1:2] #changed from 1:3, 1:2 because no soil
s3 <- as.data.frame(HPDinterval(all.param, prob=0.95))
s3 <- cbind(s2, s3)
rownames(s3) <- c("geology_bray_curtis_ratio", "genetic_fst_ratio") #removed soil from the list
s3 <- round(s3, digits = 5)
s3

#this ends the code from second BEDASSLE

# check shit out
show(load("kentucki_19_pops_BB_MCMC_output1.Robj"))
layout(t(1:2))
plot(aD, xlab="MCMC generation", ylab="value", main="aD")
plot((aD_accept/aD_moves), xlab="MCMC generation", ylab="", main="aD acceptance", ylim=c(0,1))
plot(as.vector(aE), xlab="MCMC generation", ylab="value", main="aE")
plot((aE_accept/aE_moves), xlab="MCMC generation", ylab="", main="aE acceptance", ylim=c(0,1))
plot(a2, xlab="MCMC generation", ylab="value", main="a2")
plot((a2_accept/a2_moves), xlab="MCMC generation", ylab="", main="a2 acceptance", ylim=c(0,1))
plot(Prob, xlab="MCMC generation", main="log likelihood")
plot((mu_accept/mu_moves), xlab="MCMC generation", ylab="", main="mu acceptance", ylim=c(0,1) )
plot((thetas_accept/thetas_moves)[-(1:40)], xlab="MCMC generation", ylab="", main="thetas acceptance", ylim=c(0,1) )
hist((aE/aD),breaks=100,main="posterior of aE/aD ratio")

# back transform to get  per km effect
#these are all coming back NaN
aD <- aD/geo.effect
aE <- aE/env.effect
mean(aE[-c(1:2000)]/aD[-c(1:2000)]) 
sd(aE[-c(1:2000)]/aD[-c(1:2000)]) 


##From a tutorial I FINALLY found
layout(matrix(1:6,nrow=3))
cols <- adjustcolor(rainbow(64)[cut(seq_along(a0),breaks=64)],0.5)

plot(as.vector(aE), as.vector(aD), col=cols, pch=20, xlab=aE, ylab=aD) 
plot(as.vector(a2), as.vector(aE), col=cols, pch=20, xlab=a2, ylab=aE)
plot(as.vector(a2), as.vector(aD), col=cols, pch=20, xlab=a2, ylab=aD)
plot(as.vector(a0), as.vector(aE), col=cols, pch=20, xlab=a0, ylab=aE)
plot(as.vector(a0), as.vector(aD), col=cols, pch=20, xlab=a0, ylab=aD)
plot(as.vector(a0), as.vector(a2), col=cols, pch=20, xlab=a0, ylab=a2)

layout(1:1)

hist(aE/aD,breaks=100,main="posterior of aE/aD ratio")


################################################################################
################################################################################
################################################################################
################################################################################

###################testing IBD vs IBE in east group#############################

# set wd and read in locality data
setwd("~/Library/CloudStorage/GoogleDrive-efw24@vt.edu/My Drive/OUwork/Dissertation/Chapter2/BEDASSLE/east")
localities.east <- read.csv("east_localities.csv")
pop.loc.east <- cbind.data.frame(localities.east$Longitude,localities.east$Latitude,localities.east$Locality)
pop.loc.east <- cbind.data.frame(unique(localities.east$Longitude),unique(localities.east$Latitude),unique(localities.east$Locality), unique(localities.east$Elevation..m.))
colnames(pop.loc.east) <- c("long", "lat", "pop", "elevation")
geodist.east <- earth.dist(pop.loc.east[c("long", "lat")], dist = FALSE)

#convert STRUCTURE file instead
ky.gen.east <- read.structure("kentucki_east.str", n.ind=11, n.loc=6803, onerowperind=FALSE,
                         col.lab=1, col.pop=2, row.marknames=0,
                         NA.char="-9", ask=TRUE, quiet=FALSE)
md.east <- propTyped(ky.gen.east, by=c("both")) # confirm we have expected amount of missing data
table(md.east)[1]/table(md.east)[2] # proportion of missing sites
table(ky.gen.east@tab)
#ky.pops <-  gsub( "_.*$", "", rownames(ky.gen@tab)) # add pops

ky.pops.east <- as.factor(c("1","2","3","4","5","5","5","6","7","8","9")) #these are each as individuals
ky.b.east <- genind2genpop(ky.gen.east, ky.pops.east) # turn into gen pop

# convert to BEDASSLE format
ky.ac.east <- as.matrix(ky.b.east@tab)
del.east <- seq(2, ncol(ky.ac.east), 2) # sequence of integers to drop non-ref allele
ky.ac.east <- ky.ac.east[,-del.east] 

# create equally sized matrix for sample sizes
ky.n.east <- matrix(nrow=nrow(ky.ac.east), ncol=ncol(ky.ac.east))

# name our rows the same thing
rownames(ky.n.east) <- rownames(ky.ac.east)

# get sample size per population
sample.n.east <- table(ky.gen.east@pop) 

# turn this into a vector
sample.sizes.east <- as.vector(sample.n.east)

# populate each row of matrix with sample sizes for pops
for(i in 1:nrow(ky.n.east)){
  ky.n.east[i,] <- sample.sizes.east[i]
}
ky.n.east <- ky.n.east*2 # adjust to account for loss of one allele

# calculate pairwise Fst
ky.p.fst.all.east <- calculate.all.pairwise.Fst(ky.ac.east, ky.n.east) 
ky.p.fst.all.east
######done till here; looking good so far!
# look at global Fst
ky.p.fst.east <- calculate.pairwise.Fst(ky.ac.east, ky.n.east)
ky.p.fst.east

#calc distance
ky.geo.east <- earth.dist(pop.loc.east[c("long", "lat")], dist = FALSE)

# standardize geographic distances
geo.effect.east <- sd(c(ky.geo.east)) # save effect size
ky.geo.east <- ky.geo.east/sd(c(ky.geo.east)) #trying without this to see if it helps things run

# turn to vectors
sat.dist.east <- as.vector(ky.geo.east)
sat.gen.east <- as.vector(ky.p.fst.all.east)

# make data frame with these variables
sat.df.east <- cbind.data.frame(sat.dist.east, sat.gen.east)
sat.df.east$species <- rep("plethodon_ky",nrow(sat.df.east))
colnames(sat.df.east) <- c("distance", "fst", "species")

# calculate environmental distance
worldclim.east <- getData("worldclim",var="bio",res=0.5,lon=localities.east$Longitude[1],lat=localities.east$Latitude[1])
proj.east <- as.character(worldclim.east[[2]]@crs)
worldclim.east <- worldclim.east[[c(1,3,7,8,9,12,15)]]
names(worldclim.east) <- c("temp","isotherm", "temp_ann_range", "mean_temp_wet", "mean_temp_dry", "precip", "precip_season")
sp1 <- SpatialPoints(pop.loc.east[,c('long', 'lat')], proj4string=CRS(proj))
values <- extract(worldclim.east, sp1)
loc.master.east <- cbind.data.frame(pop.loc.east, values)
loc.uniq.east <- cbind.data.frame(loc.master.east$pop, loc.master.east$elevation, loc.master.east$temp,
                             loc.master.east$isotherm, loc.master.east$temp_ann_range, loc.master.east$mean_temp_wet,
                             loc.master.east$mean_temp_dry, loc.master.east$precip, loc.master.east$precip_season)
colnames(loc.uniq.east) <- c("population", "elevation", "temp", "isotherm", "temp_ann_range", 
                        "mean_temp_wet", "mean_temp_dry", "precip", "precip_season")
rownames(loc.uniq.east) <- loc.uniq.east$population
loc.uniq.east <- as.data.frame(loc.uniq.east)
loc.uniq.east <- loc.uniq[,-which(names(loc.uniq.east) %in% c("population"))]
env.dist.all.east <- dist(loc.uniq.east, diag = TRUE, upper = TRUE)

# prep dist matrix
ky.env.east <- as.matrix(env.dist.all.east)
colnames(ky.env.east) <- NULL
rownames(ky.env.east) <- NULL

# standardize env distances
env.effect.east <- sd(c(ky.env.east)) # save effect size
ky.env.east <- ky.env.east/sd(c(ky.env.east)) #trying not standardizing here too


#########################################
#########################################
########################################
#######################################
#The value of delta may set off warnings, 
#so temporarily disable warnings. 
op <- options("warn") 
options(warn =-1)

##testing BB model
MCMC_BB(
  counts = ky.ac.east,
  sample_sizes = ky.n.east,
  D = ky.geo.east,  # geographic distances #this could be an issue (may need euclidean distances)
  E = ky.env.east,  # environmental distances #this only has 19 entries, since it used unique localities (fixing first)
  k = nrow(ky.ac.east), 
  loci = ncol(ky.ac.east),  # dimensions of the data
  #delta = 0.000000000238,  # a small, positive, number #trying even smaller, since some of these are very close localities 
  #aD_stp = 0.075,   # step sizes for the MCMC
  #aE_stp = 0.05,
  #a2_stp = 0.025,
  #thetas_stp = 0.2,
  #mu_stp = 0.35,
  #delta = 0.0000000007,  
  delta = 0.000001, #default delta for both analyses
  aD_stp = 0.8,
  aE_stp = 0.5, 
  a2_stp = 0.01, 
  thetas_stp = 0.3,
  phi_stp = 80, 
  mu_stp = 2, 
  ngen = 10000,  
  printfreq = 100,
  savefreq = 10000,
  samplefreq = 100,     # record current state for posterior (2000)
  prefix = "kentucki_east_BBtest_",# filename prefix
  continue = FALSE
)

## run below code on local machine
plot_all_acceptance_rates('kentucki_east_BBtest_MCMC_output1.Robj') #shown in viewer

#make a vector population names
pop_vector <- as.vector(rownames(ky.n.east))
plot_all_trace('kentucki_east_BBtest_MCMC_output1.Robj', percent.burnin = 0, thinning = 1, population.names = pop_vector)

posterior.predictive.sample(MCMC.output = "kentucki_east_BBtest_MCMC_output1.Robj", posterior.predictive.sample.size = 1000, output.file = "east_BB_test_post1000", prefix = "east_BB_test_ky1_")
plot_posterior_predictive_samples("east_BB_test_ky1_BBtest_post1000.Robj", save.figure = TRUE, figure.name = "east_BBtest_r1_post1000_1.jpg")

#use a burn-in of the first 8 million
plot_all_trace('kentucki_east_BBtest_MCMC_output1.Robj', percent.burnin = 80, thinning = 1, population.names = pop_vector)

show(load("kentucki_east_BBtest_MCMC_output1.Robj"))


###THE BB_MCMC was the better fit, but I do still need to reduce the environmental factors
# calculate environmental distance
MCMC_BB(
  counts = ky.ac.east,
  sample_sizes = ky.n.east,
  D = ky.geo.east,  # geographic distances #this could be an issue (may need euclidean distances)
  E = ky.env.east,  # environmental distances #this only has 19 entries, since it used unique localities (fixing first)
  k = nrow(ky.ac.east), 
  loci = ncol(ky.ac.east),  # dimensions of the data
  #delta = 0.000000000238,  # a small, positive, number #trying even smaller, since some of these are very close localities 
  #aD_stp = 0.075,   # step sizes for the MCMC
  #aE_stp = 0.05,
  #a2_stp = 0.025,
  #thetas_stp = 0.2,
  #mu_stp = 0.35,
  #delta = 0.0000000007,  
  delta = 0.000001, #default delta for both analyses
  aD_stp = 0.8,
  aE_stp = 0.5, 
  a2_stp = 0.01, 
  thetas_stp = 0.3,
  phi_stp = 80, 
  mu_stp = 2, 
  ngen = 10000000,  
  printfreq = 10000,
  savefreq = 10000,
  samplefreq = 100,     # record current state for posterior (2000)
  prefix = "kentucki_19_pops_BB_",# filename prefix
  continue = FALSE
)

## run below code on local machine
plot_all_acceptance_rates('kentucki_east_BB_MCMC_output1.Robj') #shown in viewer

#make a vector population names
pop_vector <- as.vector(rownames(ky.n.east))
plot_all_trace('kentucki_east_BB_MCMC_output1.Robj', percent.burnin = 0, thinning = 1, population.names = pop_vector)

posterior.predictive.sample(MCMC.output = "kentucki_east_BB_MCMC_output1.Robj", posterior.predictive.sample.size = 1000, output.file = "east_BB_post1000", prefix = "BB_east_")
plot_posterior_predictive_samples("BB_east_east_BB_post1000.Robj", save.figure = TRUE, figure.name = "east_BB_post1000_1.jpg")

#use a burn-in of the first 8 million
plot_all_trace('kentucki_east_BB_MCMC_output1.Robj', percent.burnin = 80, thinning = 1, population.names = pop_vector)

show(load("kentucki_east_BB_MCMC_output1.Robj"))




#to correct the aE : aD ratios, multiply by sd(geom) and sd(genetic_fst)
#geom <- read.csv(geodist, row.names = 1, header = TRUE)
#geom <- as.matrix(geom)

#genetic_fst <- read.csv("population_FST.csv", row.names = 1, header = TRUE)
#genetic_fst <- as.matrix(genetic_fst)

geodist.east
ky.p.fst.all.east

#correct the aE : aD ratios
#1 is geology
geology_bray_curtis_ratio <- ((as.vector(aE[1,])* sd(geodist.east))/(aD))
#2 is soil pc1
soil_pca_pc1_ratio <- ((as.vector(aE[2,])* sd(geodist.east))/(aD)) #this one isn't working, and I'm not sure why
#3 is genetic fst
genetic_fst_ratio <- ((as.vector(aE[3,])* sd(geodist.east))/(aD * sd(ky.p.fst.all.east)))

#calculate output after 80% burn-in
library(coda)
all.param.east = as.mcmc(cbind(geology_bray_curtis_ratio[-(1:80000)], genetic_fst_ratio[-(1:80000)]))
s1.east <- summary(all.param.east)
s2.east <- as.data.frame(s1.east$statistics)
s2.east <- s2.east[1:2,1:2] #changed from 1:3, 1:2 because no soil
s3.east <- as.data.frame(HPDinterval(all.param.east, prob=0.95))
s3.east <- cbind(s2.east, s3.east)
rownames(s3.east) <- c("geology_bray_curtis_ratio", "genetic_fst_ratio") #removed soil from the list
s3.east <- round(s3.east, digits = 5)
s3.east

#this ends the code from second BEDASSLE

# check shit out
show(load("kentucki_east_BB_MCMC_output1.Robj"))
layout(t(1:2))
plot(aD, xlab="MCMC generation", ylab="value", main="aD")
plot((aD_accept/aD_moves), xlab="MCMC generation", ylab="", main="aD acceptance", ylim=c(0,1))
plot(as.vector(aE), xlab="MCMC generation", ylab="value", main="aE")
plot((aE_accept/aE_moves), xlab="MCMC generation", ylab="", main="aE acceptance", ylim=c(0,1))
plot(a2, xlab="MCMC generation", ylab="value", main="a2")
plot((a2_accept/a2_moves), xlab="MCMC generation", ylab="", main="a2 acceptance", ylim=c(0,1))
plot(Prob, xlab="MCMC generation", main="log likelihood")
plot((mu_accept/mu_moves), xlab="MCMC generation", ylab="", main="mu acceptance", ylim=c(0,1) )
plot((thetas_accept/thetas_moves)[-(1:40)], xlab="MCMC generation", ylab="", main="thetas acceptance", ylim=c(0,1) )
hist((aE/aD),breaks=100,main="posterior of aE/aD ratio")

# back transform to get  per km effect
#these are all coming back NaN
aD <- aD/geo.effect.east
aE <- aE/env.effect.east
mean(aE[-c(1:2000)]/aD[-c(1:2000)]) 
sd(aE[-c(1:2000)]/aD[-c(1:2000)]) 


##From a tutorial I FINALLY found
layout(matrix(1:6,nrow=3))
cols <- adjustcolor(rainbow(64)[cut(seq_along(a0),breaks=64)],0.5)

plot(as.vector(aE), as.vector(aD), col=cols, pch=20, xlab=aE, ylab=aD) 
plot(as.vector(a2), as.vector(aE), col=cols, pch=20, xlab=a2, ylab=aE)
plot(as.vector(a2), as.vector(aD), col=cols, pch=20, xlab=a2, ylab=aD)
plot(as.vector(a0), as.vector(aE), col=cols, pch=20, xlab=a0, ylab=aE)
plot(as.vector(a0), as.vector(aD), col=cols, pch=20, xlab=a0, ylab=aD)
plot(as.vector(a0), as.vector(a2), col=cols, pch=20, xlab=a0, ylab=a2)

layout(1:1)

hist(aE/aD,breaks=100,main="posterior of aE/aD ratio")

################################################################################
################################################################################
################################################################################
################################################################################

###################testing IBD vs IBE in west group#############################

