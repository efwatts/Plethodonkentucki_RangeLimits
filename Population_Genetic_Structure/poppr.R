##I'm doing this to measure Ho, He, and theta

library(poppr)
library(adegenet)
library(ape)
library(ade4)

#convert STRUCTURE object to genind
ky.gen <- read.structure("kentuckionly_labeledpops.str", n.ind=21, n.loc=6803, onerowperind=FALSE,
                         col.lab=1, col.pop=2, row.marknames=0,
                         NA.char="-9", ask=TRUE, quiet=FALSE)
md <- propTyped(ky.gen, by=c("both")) # confirm we have expected amount of missing data
table(md)[1]/table(md)[2] # proportion of missing sites
table(ky.gen@tab)

genind2genalex(ky.gen, filename = "ky_genalex.csv", overwrite = FALSE, quiet = FALSE, 
               pop = NULL, allstrata = TRUE, geo = FALSE, geodf = "xy", 
               sep = ",", sequence = FALSE)
ky.glex <- read.genalex("ky_genalex.csv", ploidy = 2, geo = FALSE, region = FALSE, 
             genclone = TRUE, sep = ",", recode = FALSE)

popdata <- poppr(ky.glex)
popdata
