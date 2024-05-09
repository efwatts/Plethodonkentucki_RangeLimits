##SDMTools for paper

#These are the parameters I used for the paper
##-80.5 – -84.5° longitude and 36 – 39.5° latitude
##19 bioclimatic variables from the WorldClim (http://www.worldclim.org) database at 30 arc-seconds resolution 

library(ggplot2)    # To plot locations
library(maps)       # To access useful maps
library(rasterVis)  # To plot raster objects
library(dismo)
library(raster)

#get climate data
# Set the path to the directory containing WorldClim data for North America
files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd",
                    full.names = TRUE)

#convert to raster object
predictors <- terra::rast(files)

k_folder <- "/Users/emily/Downloads/wc2.1_30s_bio"
k_raster_files <- list.files(path = k_folder, pattern = "\\.tif$", full.names = TRUE)
k_predictors <- terra::rast(k_raster_files)

# Assuming 'k_predictors' is your original SpatRaster object
# Set the new extent
new_extent <- extent(-84.5, -80.5, 36, 39.5)

# Crop the SpatRaster object to the new extent
k_predictors_cropped <- crop(k_predictors, new_extent)

# Print the updated SpatRaster object
print(k_predictors_cropped)


#see the environmental variables
names(predictors)

#prepare presence and background points
library(SDMtune)

#this is an example species
help(virtualSp)
p_coords <- virtualSp$presence
bg_coords <- virtualSp$background

library(readr)
kentucki <- read_csv("onlykentucki_rarefied_points.csv")
pk_coords <- as.data.frame(kentucki[, c("LONG", "LAT")])
colnames(pk_coords) <- c("x", "y")
print(pk_coords)

# Set the number of points
num_points <- 10000
# Set the latitude and longitude limits
latitude_limits <- c(36, 39.5)
longitude_limits <- c(-84.5, -80.5)
# Generate random latitude and longitude points
random_lat <- runif(num_points, min = latitude_limits[1], max = latitude_limits[2])
random_lon <- runif(num_points, min = longitude_limits[1], max = longitude_limits[2])
# Create a data frame with the random points
random_points <- data.frame(x = random_lon, y = random_lat)
# Display the first few rows of the resulting data frame
head(random_points)
bk_coords <- random_points
print(bk_coords)

#plot with presence locations (change view accordingly)
ggplot(map_data("world"), aes(long, lat)) +
  geom_polygon(aes(group = group), fill = "grey95", color = "gray40", size = 0.2) +
  geom_jitter(data = pk_coords, aes(x = x, y = y), color = "red",
              alpha = 0.4, size = 1) +
  labs(x = "longitude", y = "latitude") +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_fixed() +
  scale_x_continuous(limits = c(-84.5, -80.5)) +
  scale_y_continuous(limits = c(36, 39.5))

#plot background locations
ggplot(map_data("world"), aes(long, lat)) +
  geom_polygon(aes(group = group), fill = "grey95", color = "gray40", size = 0.2) +
  geom_jitter(data = as.data.frame(bk_coords), aes(x = x, y = y),
              color = "blue", alpha = 0.4, size = 0.5) +
  labs(x = "longitude", y = "latitude") +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_fixed() +
  scale_x_continuous(limits = c(-84.5, -80.5)) +
  scale_y_continuous(limits = c(36, 39.5))

#create SWD object to use for model training 
data <- prepareSWD(species = "kentucki", 
                   p = pk_coords, 
                   a = bk_coords,
                   env = k_predictors_cropped, 
                   categorical = "wc2.1_30s_bio_1")
data

#visualize data 
head(data@data)
#visualize coordinates
head(data@coords)
#visualize name of species 
data@species 

#save all data
swd2csv(data, 
        file_name = "final_data.csv")

#save data separately with presence and background
swd2csv(data, 
        file_name = c("final_presence.csv", "final_background.csv"))

##train the model 
#options are articifcial neural networks (ann), boosted regression trees (brt), MaxEnt using dismo or maxnet, or random forest RF...code should be the same across
#using maxent
library(rJava)
library(maxnet)

#these are the defaut settings 
model <- train(method = "Maxent", 
               data = data)
model
slotNames(model)
slotNames(model@model)

#here we create the class of models:
#lqpht: linear, quadratic, product, hinge, and threshold features
lqpht.1 <- train(method = "Maxent", 
               data = data, 
               fc = "lqpht", 
               reg = 1, 
               iter = 100)
lqpht.1.5 <- train(method = "Maxent", 
                data = data, 
                fc = "lqpht", 
                reg = 1.5, 
                iter = 100)
lqpht.2 <- train(method = "Maxent", 
                  data = data, 
                  fc = "lqpht", 
                  reg = 2, 
                  iter = 100)
lqpht.2.5 <- train(method = "Maxent", 
                  data = data, 
                  fc = "lqpht", 
                  reg = 2.5, 
                  iter = 100)
lqpht.3 <- train(method = "Maxent", 
                  data = data, 
                  fc = "lqpht", 
                  reg = 3, 
                  iter = 100)
lqpht.3.5 <- train(method = "Maxent", 
                  data = data, 
                  fc = "lqpht", 
                  reg = 3.5, 
                  iter = 100)
lqpht.4 <- train(method = "Maxent", 
                  data = data, 
                  fc = "lqpht", 
                  reg = 4, 
                  iter = 100)
lqpht.4.5 <- train(method = "Maxent", 
                  data = data, 
                  fc = "lqpht", 
                  reg = 4.5, 
                  iter = 100)
lqpht.5 <- train(method = "Maxent", 
                  data = data, 
                  fc = "lqpht", 
                  reg = 5, 
                  iter = 100)

#lqh: linear, quadratic, and hinge features
lqh.1 <- train(method = "Maxent", 
                 data = data, 
                 fc = "lqh", 
                 reg = 1, 
                 iter = 100)
lqh.1.5 <- train(method = "Maxent", 
                   data = data, 
                   fc = "lqh", 
                   reg = 1.5, 
                   iter = 100)
lqh.2 <- train(method = "Maxent", 
                 data = data, 
                 fc = "lqh", 
                 reg = 2, 
                 iter = 100)
lqh.2.5 <- train(method = "Maxent", 
                   data = data, 
                   fc = "lqh", 
                   reg = 2.5, 
                   iter = 100)
lqh.3 <- train(method = "Maxent", 
                 data = data, 
                 fc = "lqh", 
                 reg = 3, 
                 iter = 100)
lqh.3.5 <- train(method = "Maxent", 
                   data = data, 
                   fc = "lqh", 
                   reg = 3.5, 
                   iter = 100)
lqh.4 <- train(method = "Maxent", 
                 data = data, 
                 fc = "lqh", 
                 reg = 4, 
                 iter = 100)
lqh.4.5 <- train(method = "Maxent", 
                   data = data, 
                   fc = "lqh", 
                   reg = 4.5, 
                   iter = 100)
lqh.5 <- train(method = "Maxent", 
                 data = data, 
                 fc = "lqh", 
                 reg = 5, 
                 iter = 100)

#lq: linear and quadratic features
lq.1 <- train(method = "Maxent", 
               data = data, 
               fc = "lq", 
               reg = 1, 
               iter = 100)
lq.1.5 <- train(method = "Maxent", 
                 data = data, 
                 fc = "lq", 
                 reg = 1.5, 
                 iter = 100)
lq.2 <- train(method = "Maxent", 
               data = data, 
               fc = "lq", 
               reg = 2, 
               iter = 100)
lq.2.5 <- train(method = "Maxent", 
                 data = data, 
                 fc = "lq", 
                 reg = 2.5, 
                 iter = 100)
lq.3 <- train(method = "Maxent", 
               data = data, 
               fc = "lq", 
               reg = 3, 
               iter = 100)
lq.3.5 <- train(method = "Maxent", 
                 data = data, 
                 fc = "lq", 
                 reg = 3.5, 
                 iter = 100)
lq.4 <- train(method = "Maxent", 
               data = data, 
               fc = "lq", 
               reg = 4, 
               iter = 100)
lq.4.5 <- train(method = "Maxent", 
                 data = data, 
                 fc = "lq", 
                 reg = 4.5, 
                 iter = 100)
lq.5 <- train(method = "Maxent", 
               data = data, 
               fc = "lq", 
               reg = 5, 
               iter = 100)

#l: linear features
l.1 <- train(method = "Maxent", 
              data = data, 
              fc = "l", 
              reg = 1, 
              iter = 100)
l.1.5 <- train(method = "Maxent", 
                data = data, 
                fc = "l", 
                reg = 1.5, 
                iter = 100)
l.2 <- train(method = "Maxent",
              data = data, 
              fc = "l", 
              reg = 2, 
              iter = 100)
l.2.5 <- train(method = "Maxent", 
                data = data, 
                fc = "l", 
                reg = 2.5, 
                iter = 100)
l.3 <- train(method = "Maxent",
              data = data, 
              fc = "l", 
              reg = 3, 
              iter = 100)
l.3.5 <- train(method = "Maxent", 
                data = data, 
                fc = "l", 
                reg = 3.5, 
                iter = 100)
l.4 <- train(method = "Maxent",
              data = data, 
              fc = "l", 
              reg = 4, 
              iter = 100)
l.4.5 <- train(method = "Maxent", 
                data = data, 
                fc = "l", 
                reg = 4.5, 
                iter = 100)
l.5 <- train(method = "Maxent",
              data = data, 
              fc = "l", 
              reg = 5, 
              iter = 100)

#h: hinge features
h.1 <- train(method = "Maxent", 
              data = data, 
              fc = "h", 
              reg = 1, 
              iter = 100)
h.1.5 <- train(method = "Maxent", 
                data = data, 
                fc = "h", 
                reg = 1.5, 
                iter = 100)
h.2 <- train(method = "Maxent",
              data = data, 
              fc = "h", 
              reg = 2, 
              iter = 100)
h.2.5 <- train(method = "Maxent", 
                data = data, 
                fc = "h", 
                reg = 2.5, 
                iter = 100)
h.3 <- train(method = "Maxent",
              data = data, 
              fc = "h", 
              reg = 3, 
              iter = 100)
h.3.5 <- train(method = "Maxent", 
                data = data, 
                fc = "h", 
                reg = 3.5, 
                iter = 100)
h.4 <- train(method = "Maxent",
              data = data, 
              fc = "h", 
              reg = 4, 
              iter = 100)
h.4.5 <- train(method = "Maxent", 
                data = data, 
                fc = "h", 
                reg = 4.5, 
                iter = 100)
h.5 <- train(method = "Maxent",
              data = data, 
              fc = "h", 
              reg = 5, 
              iter = 100)

#evaluate model performance via prediction 
# Assuming models is a list of your trained models
model_names <- c("lqpht.1", "lqpht.1.5", "lqpht.2", "lqpht.2.5", "lqpht.3", "lqpht.3.5", "lqpht.4", "lqpht.4.5", "lqpht.5",
                 "lqh.1", "lqh.1.5", "lqh.2", "lqh.2.5", "lqh.3", "lqh.3.5", "lqh.4", "lqh.4.5", "lqh.5",
                 "lq.1", "lq.1.5", "lq.2", "lq.2.5", "lq.3", "lq.3.5", "lq.4", "lq.4.5", "lq.5",
                 "l.1", "l.1.5", "l.2", "l.2.5", "l.3", "l.3.5", "l.4", "l.4.5", "l.5",
                 "h.1", "h.1.5", "h.2", "h.2.5", "h.3", "h.3.5", "h.4", "h.4.5", "h.5")

models <- mget(model_names)

# Initialize a data frame to store the results
results <- data.frame(Model = character(),
                      OmissionRate = numeric(),
                      AUC = numeric(),
                      AICc = numeric(),  # Add this line
                      stringsAsFactors = FALSE)

# Now run the loop
for (i in seq_along(models)) {
  # Calculate the AUC
  auc_value <- auc(models[[i]])
  
  # Calculate the TSS
  tss_value <- tss(models[[i]])
  
  # Calculate the omission rate
  omission_rate <- 1 - tss_value
  
  # Calculate the AICc
  aicc_value <- aicc(models[[i]], env = k_predictors)  # replace k_predictors with the correct number of predictors
  
  # Print the values
  print(paste("Model:", names(models)[i]))
  print(paste("AUC:", auc_value))
  print(paste("TSS:", tss_value))
  print(paste("Omission_Rate:", omission_rate))
  print(paste("AICc:", aicc_value))
  
  # Append the result to the results data frame
  results <- rbind(results, data.frame(Model = names(models)[i],
                                       OmissionRate = omission_rate,
                                       AUC = auc_value,
                                       AICc = aicc_value))
}

# Rank the models based on lowest omission rate, highest AUC, and lowest AICc
results$Rank <- with(results, rank(OmissionRate, ties.method = "min") +
                       rank(-AUC, ties.method = "min") +
                       rank(AICc, ties.method = "min"))


# Order the results by the rank
results <- results[order(results$Rank), ]

# Print the results
print(results)

