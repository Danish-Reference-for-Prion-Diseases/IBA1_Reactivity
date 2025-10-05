library(here)
library("caroline")
library(tidyverse)
library(viridis)
library(ggplot2)

#Datasets
file_path <- here("Data", "Cellprofiler (morphology)", "MoBIE")
cp_data <- read.delim(here(file_path,"Microglia_total.txt"))# A Raw datafile of the total IBA1+ Objects from CP
cluster_id <- read.delim(here(file_path, "Microglia_clusterID.txt")) #Clust ID From the Cluster analysis
geomx_data <- read.csv(here("Data", "GeoMx (protein)", "3. HK normalized Log2 Transformed GeoMx Data.csv"))#Processed GeoMx Data

#Combine cluster_id with data
data_with_cluster <- merge(cp_data, cluster_id, by = c("ImageNumber", "ObjectNumber"), all.y=TRUE, all.x = TRUE)

#Imputate Nans with zero
data_with_cluster <- data_with_cluster %>% mutate(Cluster = ifelse(is.na(Cluster), 0, Cluster))

#Transform cluster to string
data_with_cluster$Cluster <- as.character(data_with_cluster$Cluster)
data_with_cluster <- transform(data_with_cluster, Cluster = sprintf("Cluster%s", Cluster))

#Add color scheme
data_with_cluster <- data_with_cluster %>% mutate(figuresColorScheme=case_when(
  Cluster=="Cluster0"~"224-224-224-220",
  Cluster=="Cluster1"~"0-123-255-220",
  Cluster=="Cluster2"~"220-20-60-220",
  Cluster=="Cluster3"~"225-165-0-220"))

#Rename columns before merging geomx data and cellprofiler data
geomx_data <- geomx_data %>% rename("ROI_Label" = "ROI..Label.")
data_with_cluster <- data_with_cluster %>% rename("ROI_Label" = "Metadata_FileLocation")
data_with_cluster$ROI_Label <- gsub("mask-", "", data_with_cluster$FileName_ROI_Mask)#Create ROI Label in morph data
data_with_cluster$ROI_Label <- gsub(".tif", "", data_with_cluster$ROI_Label)#Create ROI Label in morph data

#Merge
data_with_geomx <- merge(data_with_cluster, geomx_data, by = "ROI_Label", all.y=TRUE, all.x = TRUE)

#Add label ID
label_id <- 1:nrow(data_with_geomx)
data_with_id<- cbind(label_id, data_with_geomx)

#Replace NA with zero
data_with_id[is.na(data_with_id)] <- 0

#Export 
write.delim(data_with_id, here(file_path, "Microglia_total_withID.txt"),
            row.names = FALSE, sep = "\t")

sessionInfo()
