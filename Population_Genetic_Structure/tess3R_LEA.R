
##first, tess3r tutorial to be sure it all works

#these are additional functions to import input files from the structure format, 
#and to display nice geographic representations of ancestry coefficients with maps.
library(LEA)
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")


input.file = "http://membres-timc.imag.fr/Olivier.Francois/secondary_contact.str" 

###allelic markers
struct2geno(file = input.file, TESS = TRUE, diploid = TRUE, FORMAT = 2,
            extra.row = 0, extra.col = 0, output = "secondary_contact.geno")

obj.snmf = snmf("secondary_contact.geno", K = 3, alpha = 100, project = "new") 
qmatrix = Q(obj.snmf, K = 3)

barplot(t(qmatrix), col = c("orange","violet","lightgreen"), border = NA, space = 0, xlab = "Individuals", ylab = "Admixture coefficients")
coord = read.table("coordinates.coord")
pop = rep(1:60, each = 10)

K = 3 
Npop = length(unique(pop)) 
qpop = matrix(NA, ncol = K, nrow = Npop) 
coord.pop = matrix(NA, ncol = 2, nrow = Npop) 
for (i in unique(pop)){ 
  qpop[i,] = apply(qmatrix[pop == i,], 2, mean) 
  coord.pop[i,] = apply(coord[pop == i,], 2, mean)}

library(mapplots)
library(maps)
plot(coord, xlab = "Longitude", ylab = "Latitude", type = "n") 
map(add = T, col = "grey90", fill = TRUE) 
for (i in 1:Npop){ 
  add.pie(z = qpop[i,], x = coord.pop[i,1], y = coord.pop[i,2], labels = "", 
          col = c("orange","violet","lightgreen"))}

pop = scan("mypop.txt") #for my own data, I'll need to add a list of populations

###choosing clusters
obj.snmf = snmf("secondary_contact.geno", K = 1:8, ploidy = 2, entropy = T, alpha = 100, project = "new") 
plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19) #choose the value of K for which the cross-entropy curve exhibits a plateau, here K=3


###SNP data
url = "http://membres-timc.imag.fr/Olivier.Francois/Arabidopsis/A_thaliana_chr1.geno" 
download.file(url = url, destfile = "./A_thaliana_chr1.geno") 
url = "http://membres-timc.imag.fr/Olivier.Francois/Arabidopsis/at_coord.coord" 
download.file(url = url, destfile = "./at_coord.coord")

obj.at = snmf("./A_thaliana_chr1.geno", K = 1:10, ploidy = 1, entropy = T, CPU = 1, project = "new") 
plot(obj.at, col = "blue4", cex = 1.4, pch = 19)

qmatrix = Q(obj.at, K =5)

asc.raster="http://membres-timc.imag.fr/Olivier.Francois/RasterMaps/Europe.asc" 
grid=createGridFromAsciiRaster(asc.raster) 
constraints=getConstraintsFromAsciiRaster(asc.raster, cell_value_min=0) 
coord.at = read.table("at_coord.coord")

maps(matrix = qmatrix, coord.at, grid, constraints, method = "max", 
     main = "Ancestry coefficients", xlab = "Longitude", ylab = "Latitude", cex = .5) 
map(add = T, interior = F) #this didn't work, but I prefer the pie charts from the last example anyway



###loading my data into LEA
setwd("~/Library/CloudStorage/GoogleDrive-efw24@vt.edu/My Drive/OUwork/Dissertation/Chapter2/tess3r/glutinosus_id")

library(LEA)
##this part is all glutinosus, but I'm not really worried about the conclusions there. Rerunning just for better K figures
struct2geno("allpops_labeledindiv.str", ploidy=2, FORMAT = 2, extra.col = 2) 
##identify glutinosus in the data
#struct2geno(file = "kentucki_glutremoved.str", TESS = FALSE, diploid = TRUE, FORMAT = 2, extra.row = 0, extra.col = 2, output = "./kentucki_glut_genotype.geno")

#choose K
obj.snmf = snmf("allpops_labeledindiv.str.geno", K = 1:20, ploidy = 2, entropy = T, alpha = 100, project = "new", repetitions = 20, iterations = 500000) #going to pick 2 anyway because looking for glutinosus vs kentucki
plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19) #choose the value of K for which the cross-entropy curve exhibits a plateau, here K=2 or K=3


##using K =2 for glutinosus
obj.snmf = snmf("allpops_labeledindiv.str.geno", K = 2, alpha = 100, project = "new",  iterations = 500000) 
qmatrix = Q(obj.snmf, K = 2)

barplot(t(qmatrix), col = c("cornflowerblue","lightcoral"), border = NA, space = 0, xlab = "Individuals", ylab = "Admixture coefficients")
#grabbed this for loop...should give sample site names


coord = read.table("all_coordinates.coord")
pop = rep(1:35, each = 1)

K = 2 
Npop = length(unique(pop)) 
qpop = matrix(NA, ncol = K, nrow = Npop) 
coord.pop = matrix(NA, ncol = 2, nrow = Npop) 

for (i in 1:length(pop)){ 
  qpop[i,] = qmatrix[pop == i,] #qpop[i,] = apply(qmatrix[pop == i,], 2, mean)
  coord.pop[i,] = apply(coord[pop == i,], 2, min)}
qpop

library(mapplots)
library(maps)
library(raster)
model.extent<-extent(-84,-81,36,40) #numbers are for map of continental US (lat, long: bottom, right, top, left)

plot(model.extent, xlab = "Longitude", ylab = "Latitude", type ="n")
map('state',xlim = c(-84,-81), ylim = c(36,40), fill=TRUE, add=TRUE, col="grey90")
for (i in 1:Npop){ 
  add.pie(z = qpop[i,], x = coord.pop[i,1], y = coord.pop[i,2], labels = "", 
          col = c("cornflowerblue","lightcoral"), radius = 0.1)}

##now using only kentucki data
setwd("~/Library/CloudStorage/GoogleDrive-efw24@vt.edu/My Drive/OUwork/Dissertation/Chapter2/tess3r/kentucki_only")

library(LEA)
struct2geno("reroder_kentuckionly_labeledindiv.str", ploidy=2, FORMAT = 2, extra.col = 2) 
#struct2geno(file = "kentucki_glutremoved.str", TESS = FALSE, diploid = TRUE, FORMAT = 2, extra.row = 0, extra.col = 2, output = "./kentucki_glut_genotype.geno")

#choose K
obj.snmf = snmf("reroder_kentuckionly_labeledindiv.str.geno", K = 1:20, ploidy = 2, entropy = T, alpha = 100, project = "new", repetitions = 20, iterations = 500000) #using 20 here
plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19) #choose the value of K for which the cross-entropy curve exhibits a plateau, here K=2 or K=3

obj.snmf = snmf("kentuckionly_labeledindiv.str.geno", K = 1:10, ploidy = 2, entropy = T, alpha = 100, project = "new", repetitions = 20, iterations = 500000) #using 10 here to get a better look
plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19) #choose the value of K for which the cross-entropy curve exhibits a plateau, here K=2 or K=3...2 might be a little better


##using K=2 because it was the best fit
####here, using reordered only for barplot (can redo if using localities)
obj.snmf = snmf("reroder_kentuckionly_labeledindiv.str.geno", K = 2, alpha = 100, project = "new", iterations = 500000) 
qmatrix = Q(obj.snmf, K = 2)

barplot(t(qmatrix), col = c("tomato", "lightblue"), border = NA, space = 0, xlab = "Individuals", ylab = "Admixture coefficients")

coord = read.table("kentucki.coord")
pop = rep(1:21, each = 1)

K = 2 
Npop = length(unique(pop)) 
qpop = matrix(NA, ncol = K, nrow = Npop) 
coord.pop = matrix(NA, ncol = 2, nrow = Npop) 

for (i in 1:length(pop)){ 
  qpop[i,] = qmatrix[pop == i,] #qpop[i,] = apply(qmatrix[pop == i,], 2, mean)
  coord.pop[i,] = apply(coord[pop == i,], 2, min)}
qpop

library(mapplots)
library(maps)
library(raster)
model.extent<-extent(-84,-81,36,40) #numbers are for map of continental US (lat, long: bottom, right, top, left)

plot(model.extent, xlab = "Longitude", ylab = "Latitude", type ="n")
map('state',xlim = c(-84,-81), ylim = c(36,40), fill=TRUE, add=TRUE, col="grey90")
for (i in 1:Npop){ 
  add.pie(z = qpop[i,], x = coord.pop[i,1], y = coord.pop[i,2], labels = "", 
          col = c("cornflowerblue","lightcoral"), radius = 0.1)}

##now trying K=3 because it was also a good fit
####here, using reordered only for barplot (can redo if using localities)
obj.snmf = snmf("reroder_kentuckionly_labeledindiv.str.geno", K = 3, alpha = 100, project = "new", iterations = 500000) 
qmatrix = Q(obj.snmf, K = 3)

barplot(t(qmatrix), col = c("tomato", "orange","lightblue"), border = NA, space = 0, xlab = "Individuals", ylab = "Admixture coefficients")

coord = read.table("kentucki.coord")
pop = rep(1:21, each = 1)

K = 3 
Npop = length(unique(pop)) 
qpop = matrix(NA, ncol = K, nrow = Npop) 
coord.pop = matrix(NA, ncol = 2, nrow = Npop) 

for (i in 1:length(pop)){ 
  qpop[i,] = qmatrix[pop == i,] #qpop[i,] = apply(qmatrix[pop == i,], 2, mean)
  coord.pop[i,] = apply(coord[pop == i,], 2, min)}
qpop

library(mapplots)
library(maps)
library(raster)
model.extent<-extent(-84,-81,36,40) #numbers are for map of continental US (lat, long: bottom, right, top, left)

plot(model.extent, xlab = "Longitude", ylab = "Latitude", type ="n")
map('state',xlim = c(-84,-81), ylim = c(36,40), fill=TRUE, add=TRUE, col="grey90")
for (i in 1:Npop){ 
  add.pie(z = qpop[i,], x = coord.pop[i,1], y = coord.pop[i,2], labels = "", 
          col = c("cornflowerblue","lightcoral","seagreen"), radius = 0.1)}

##now trying K=4 because of my last paper
####here, using reordered only for barplot (can redo if using localities)
obj.snmf = snmf("reroder_kentuckionly_labeledindiv.str.geno", K = 4, alpha = 100, project = "new", iterations = 500000) 
qmatrix = Q(obj.snmf, K = 4)

barplot(t(qmatrix), col = c("plum4","lightblue","tomato","orange"), border = NA, space = 0, xlab = "Individuals", ylab = "Admixture coefficients")

coord = read.table("kentucki.coord")
pop = rep(1:21, each = 1)

K = 4 
Npop = length(unique(pop)) 
qpop = matrix(NA, ncol = K, nrow = Npop) 
coord.pop = matrix(NA, ncol = 2, nrow = Npop) 

for (i in 1:length(pop)){ 
  qpop[i,] = qmatrix[pop == i,] #qpop[i,] = apply(qmatrix[pop == i,], 2, mean)
  coord.pop[i,] = apply(coord[pop == i,], 2, min)}
qpop

library(mapplots)
library(maps)
library(raster)
model.extent<-extent(-84,-81,36,40) #numbers are for map of continental US (lat, long: bottom, right, top, left)

plot(model.extent, xlab = "Longitude", ylab = "Latitude", type ="n")
map('state',xlim = c(-84,-81), ylim = c(36,40), fill=TRUE, add=TRUE, col="grey90")
for (i in 1:Npop){ 
  add.pie(z = qpop[i,], x = coord.pop[i,1], y = coord.pop[i,2], labels = "", 
          col = c("cornflowerblue","lightcoral","seagreen","lightblue"), radius = 0.1)}



################################################################################
##############now trying actual tess3r##########################################
################################################################################
setwd("~/Library/CloudStorage/GoogleDrive-efw24@vt.edu/My Drive/OUwork/Dissertation/Chapter2/tess3r/tess3r_kentucki")
library(tess3r)
library(rworldmap)

#determine individuals who are glutinosus with tess3r too
genotype.glut <-read.geno("allpops_labeledindiv.str.geno")
coord.glut <- read.table("all_coordinates.coord")
glut.matrix <- as.matrix(coord.glut)

# Running the tess3 function
tess3.obj <- tess3(X = genotype.glut, coord = glut.matrix, K = 1:20,
                   method = "projected.ls",
                   ploidy = 2, rep = 20, max.iteration = 500000)

# Plot error
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")

# Retrieve the Q-matrix for K = 16 clusters is the best fit
q.matrix <- qmatrix(tess3.obj, K = 2)

## STRUCTURE-like barplot for the Q-matrix
barplot(q.matrix, border = NA, space = 0,
        xlab = "Individuals", ylab = "Ancestry proportions",
        main = "Ancestry matrix") -> bp
axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4)


# load data set
genotype.ky <- read.geno("kentuckionly_labeledindiv.str.geno")
coord.ky <- read.table("kentucki.coord")
coord.matrix <- as.matrix(coord.ky)

# Running the tess3 function
tess3.obj <- tess3(X = genotype.ky, coord = coord.matrix, K = 1:20,
                   method = "projected.ls",
                   ploidy = 2, rep = 20, max.iteration = 500000)

# Plot error
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")

# Retrieve the Q-matrix for K = 3 clusters 
q.matrix <- qmatrix(tess3.obj, K = 3)

## STRUCTURE-like barplot for the Q-matrix
barplot(q.matrix, border = NA, space = 0,
        xlab = "Individuals", ylab = "Ancestry proportions",
        main = "Ancestry matrix", col = my.palette) -> bp
axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4)

## Spatial interpolation of ancestry coefficient
my.colors <- c("orange", "lightblue", "tomato")
my.palette <- CreatePalette(my.colors, 9)

plot(q.matrix, coord.ky, method = "map.max", interpol = FieldsKrigModel(10),
     main = "Ancestry coefficients",
     xlab = "Longitude", ylab = "Latitude",
     resolution = c(500,500), cex = .4,
     col.palette = my.palette)
map('state',xlim = c(-84,-81), ylim = c(36,40), fill=FALSE, add=TRUE)

## Genome scan p-values for K = 3
p.values <- pvalue(tess3.obj, K = 3)
hist(p.values, col = "lightblue")

## Manhatan plot
plot(p.values, main = "Manhattan plot",
     xlab = "Locus id",
     ylab = "-log10(P-values)",
     cex = .3, col = "grey")

### Retrieve the Q-matrix for K = 2 clusters 
q.matrix <- qmatrix(tess3.obj, K = 2)

## STRUCTURE-like barplot for the Q-matrix
barplot(q.matrix, border = NA, space = 0,
        xlab = "Individuals", ylab = "Ancestry proportions",
        main = "Ancestry matrix", col = my.palette) -> bp
axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4)

## Spatial interpolation of ancestry coefficient
my.colors <- c("tomato", "lightblue")
my.palette <- CreatePalette(my.colors, 9)

plot(q.matrix, coord.ky, method = "map.max", interpol = FieldsKrigModel(10),
     main = "Ancestry coefficients",
     xlab = "Longitude", ylab = "Latitude",
     resolution = c(500,500), cex = .4,
     col.palette = my.palette)
map('state',xlim = c(-84,-81), ylim = c(36,40), fill=FALSE, add=TRUE)

## Genome scan p-values for K = 2
p.values <- pvalue(tess3.obj, K = 2)
hist(p.values, col = "lightblue")

## Manhattan plot
plot(p.values, main = "Manhattan plot",
     xlab = "Locus id",
     ylab = "-log10(P-values)",
     cex = .3, col = "grey")


### Retrieve the Q-matrix for K = 4 clusters 
q.matrix <- qmatrix(tess3.obj, K = 4)

## STRUCTURE-like barplot for the Q-matrix
my.colors <- c("orange", "tomato","plum4","lightblue")
my.palette <- CreatePalette(my.colors, 9)
barplot(q.matrix, border = NA, space = 0,
        xlab = "Individuals", ylab = "Ancestry proportions",
        main = "Ancestry matrix", col = my.palette) -> bp
axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4)

## Spatial interpolation of ancestry coefficient


plot(q.matrix, coord.ky, method = "map.max", interpol = FieldsKrigModel(10),
     main = "Ancestry coefficients",
     xlab = "Longitude", ylab = "Latitude",
     resolution = c(500,500), cex = .4,
     col.palette = my.palette)
map('state',xlim = c(-84,-81), ylim = c(36,40), fill=FALSE, add=TRUE)

## Genome scan p-values for K = 4
p.values <- pvalue(tess3.obj, K = 4)
hist(p.values, col = "lightblue")

## Manhattan plot
plot(p.values, main = "Manhattan plot",
     xlab = "Locus id",
     ylab = "-log10(P-values)",
     cex = .3, col = "grey")

