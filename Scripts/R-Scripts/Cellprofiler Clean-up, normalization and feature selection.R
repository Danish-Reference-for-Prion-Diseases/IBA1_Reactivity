###Libraries
library(here)
library(tidyverse)
library(VIM)
library(sjmisc)
library(caret)
library('collapse')
library("ComplexHeatmap")
library("circlize")



###Custom FUNCTIONS----------------------------------------------------------------------
#Rename Metadata to correct later format
recode_grouping_columns <- function(data){
  data <- data %>%
    mutate(Metadata_Cortex_region = recode(Metadata_Cortex_region, "OC" = "OCC"),
           Metadata_Layer = recode(Metadata_Layer, "1" = "sGM", "2" = "dGM", "3" = "WM"),
           Metadata_Pt_type = recode(Metadata_Pt_type, "NB"="CTRL"))
  return(data)
}


#Make aggregated profiles from the single cell dataset, including weights for each image 
weighted_aggregate <- function(merged_data){
  #get relevant annotations
  annotations <- merged_data[,c(1,2,3,5,7,8,10,12,17)] %>% #Annotations/metadata that you wish to save #IMPORTANT TO INCLUDE FileName_Mask FOR MOBIE!!!!!!!!!
    group_by(ImageNumber) %>%
    summarise(across(everything(), ~paste0(unique(.), collapse = "|")))
  
  #add weight for each cell that will be aggregated to get summed weights
  merged_data$weight <- 1/nrow(merged_data)
  
  #aggregate the features
  data_aggregated <- collap(merged_data[,(feature_column_start+2):ncol(merged_data)], by = merged_data$ImageNumber, w = merged_data$weight)
  
  #remove original single weight column
  merged_data = merged_data[,-ncol(merged_data)]
  data_aggregated = data_aggregated[,-ncol(data_aggregated)]
  
  #merge aggregated features with annotation
  data_aggregated <- merge(annotations, data_aggregated, by = "ImageNumber")
  
  #change FileName_Mask column for later incorporation with GeoMx data
  data_aggregated$FileName_ROI_Mask <- gsub("mask-", "", data_aggregated$FileName_ROI_Mask)
  data_aggregated$FileName_ROI_Mask <- gsub(".tif", "", data_aggregated$FileName_ROI_Mask)
  colnames(data_aggregated)[9] <- "ROI (Label)"
  
  return(data_aggregated)
}

#Normalize (robust z-score scaling) using Median Absolute Deviation
normalize_MAD <- function(data){
  data_MAD <- std(data, robust = "mad", append = FALSE)
  return(data_MAD)
}

#Remove noisy features (either with mean sd or minimum sd for each disease group)
remove_noise <- function(data, labels, filter_method, cutoff){
  data_for_noise <- cbind(Metadata_Pt_type = labels$Metadata_Pt_type, data)
  grouped_sd <- plyr::ddply(data_for_noise,plyr::.(Metadata_Pt_type), plyr::colwise(sd))
  grouped_sd <- as.data.frame.list(grouped_sd[2:ncol(grouped_sd)])
  if(filter_method == "mean"){
    mean_sd <- grouped_sd %>% select_if(colMeans(grouped_sd) < as.numeric(cutoff))
    filtered_noisy_features <- data[,c(colnames(mean_sd))]
    return(filtered_noisy_features)
  }
  else if(filter_method == "min"){
    minimum_sd <- grouped_sd %>% select_if(apply(grouped_sd,2,min) < as.numeric(cutoff))
    filtered_noisy_features <- data[,c(colnames(minimum_sd))]
    return(filtered_noisy_features)
  }
  else{
    return(print("Wrong filter_method"))
  }
}


eval_heatmap <- function(Data_Matrix, column_title, row_title){
  features <- t(Data_Matrix[,21:ncol(Data_Matrix)])
  
  #Colors
  heatmap_color_u <- colorRamp2(breaks = c(2, 0, -2), colors = c("yellow", "black", "blue"))
  
  #MAIN HEATMAP BASED 
  heat <- Heatmap(features, name="Z-score", 
                     column_title = column_title, row_title = row_title,
                     row_title_gp = gpar(fontsize = 10), 
                     column_title_gp = gpar(fontsize = 10),  
                     col = heatmap_color_u,
                     column_dend_reorder = FALSE,
                     clustering_method_columns = "ward.D2",
                     clustering_method_rows="ward.D2",
                     clustering_distance_columns = "spearman",
                     clustering_distance_rows = "spearman",
                     show_row_names = F,
                  heatmap_legend_param = list(legend_height = unit(60, "cm"),  # Adjust the height of the legend
                                              legend_direction = "horizontal",  # Place the legend horizontally
                                              title_position = "topcenter"))
  return(heat)
}

#-------------------------------------------------------------------------------
###IMPORT FILES
#main directory
main_dir <- here("Data", "Cellprofiler (morphology)")


#Microglia Soma and Total Cell Datasets from Cellprofiler
microglia_soma <- recode_grouping_columns(read.csv(here(main_dir, "Filtered_Microglia_Soma.csv"))) #PrimaryObject (Filtered from the SecondaryObject module)
microglia_total <- recode_grouping_columns(read.csv(here(main_dir, "Microglia_total.csv"))) #SecondaryObject



#----------------------------------------------------------------------------------------------
#Remove rendundant features and variables (Blocklisting)
microglia_soma = subset(microglia_soma, select = -c(Metadata_Frame, Metadata_FileLocation,
                                                    AreaShape_BoundingBoxMaximum_X,
                                                    AreaShape_BoundingBoxMaximum_Y,
                                                    AreaShape_BoundingBoxMinimum_X,
                                                    AreaShape_BoundingBoxMinimum_Y,
                                                    AreaShape_Center_X, AreaShape_Center_Y,
                                                    AreaShape_EulerNumber, Children_Microglia_total_Count,
                                                    Location_Center_X, Location_Center_Y, Location_Center_Z,
                                                    Neighbors_FirstClosestObjectNumber_150,
                                                    Neighbors_PercentTouching_150,
                                                    Neighbors_SecondClosestObjectNumber_150,
                                                    Number_Object_Number, Parent_Microglia_soma))


microglia_total = subset(microglia_total, select = -c(Metadata_Frame, Metadata_FileLocation,
                                                      AreaShape_BoundingBoxMaximum_X,
                                                      AreaShape_BoundingBoxMaximum_Y,
                                                      AreaShape_BoundingBoxMinimum_X,
                                                      AreaShape_BoundingBoxMinimum_Y,
                                                      AreaShape_Center_X, AreaShape_Center_Y,
                                                      AreaShape_EulerNumber, Parent_Microglia_soma,
                                                      Parent_Filtered_Microglia_Soma,
                                                      Children_Masked_Microglia_Count,
                                                      Location_Center_X, Location_Center_Y, Location_Center_Z,
                                                      Number_Object_Number))




#-------------------------------------------------------------------------------------------------------------------------
#First column with features
feature_column_start <- 21

#Add Soma and Total prefixes before merging
colnames(microglia_soma)[feature_column_start:ncol(microglia_soma)] <- paste("Soma", colnames(microglia_soma)[feature_column_start:ncol(microglia_soma)], sep = '_') #Mind the metadata columns
colnames(microglia_total)[feature_column_start:ncol(microglia_total)] <- paste("Total", colnames(microglia_total)[feature_column_start:ncol(microglia_total)], sep = '_') #Mind the metadata columns

#Merge Soma and Total Cell datasets
data_merged <- merge(microglia_soma, microglia_total, by=1:(feature_column_start-1))

#Metadata
labels <- data_merged[,1:(feature_column_start-1)] 
#Features 
features <- data_merged[feature_column_start:ncol(data_merged)] 


#-------------------------------------------------------------------------------------------------------------------------
#####Feature Filtering 
#1. Drop NA and inf columns (if any)
features <- features[, apply(features, 2, function(x) !all(is.infinite(x)))] # Remove columns with only infs
features <- Filter(function(x)!all(is.na(x)), features) #Remove columns with only NaNs


#2. Missing data:
# Evaluate missing data....
pMiss <- function(x){sum(is.na(x))/length(x)*100} # function
sort(apply(features, 2, pMiss)) # Evaluating missing data percentage in features (columns)
sort(apply(features, 1, pMiss), decreasing = TRUE) #Evaluating missing data percentage in samples (rows)
aggr_plot <- aggr(features, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, 
                  labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern")) #Summarizing plot and pattern plot



#3.NORMALIZE - (Robust standardization with Median Absolute Devation) 
normalized_features <- normalize_MAD(features)
colnames(normalized_features) <- substr(colnames(normalized_features), 1, nchar(colnames(normalized_features)) - 2) #Remove the z-suffix


#4. Variance thresholding (remove features with low variance)
nzv <- nearZeroVar(normalized_features, freqCut = 95/5, uniqueCut = 10)
dim(normalized_features)
if(length(nzv)!=0){
  normalized_features <- normalized_features[, -nzv]
}
dim(normalized_features) #Evaluate



#5. Remove noisy features
# Arguments (1. The dataset, 2. labels, 3. filter_method - either "mean" or "min", 4. standard deviation cutoff)
normalized_features <- remove_noise(normalized_features, labels, "min", 1)


#6. Correlation thresholding (filter highly correlated features)
excluded <- findCorrelation(cor(normalized_features), cutoff = 0.9, verbose = FALSE)
normalized_features <- normalized_features[-excluded]




#7. Object Outliers detection and filtering
#Find multivariate outliers using Mahalanobis Distance
mahalnobis_merged <- cbind(labels, normalized_features)#Combine metadata and features for this one
mahalnobis_scores <- mahalanobis(mahalnobis_merged[,feature_column_start:ncol(mahalnobis_merged)], 
                                 colMeans(mahalnobis_merged[,feature_column_start:ncol(mahalnobis_merged)]), 
                                 cov = cov(mahalnobis_merged[,feature_column_start:ncol(mahalnobis_merged)]), tol=1e-20) #Mahalanobis Score

cutoff <- qchisq(p = 0.9999, df = ncol(mahalnobis_merged[feature_column_start:ncol(mahalnobis_merged)])) #appropriate cutoff point chosen here (based on  later visual inspection in MoBIE)

outliers <- mahalnobis_merged[mahalnobis_scores > cutoff ,] #Identify outliers

Filtered_data <- anti_join(mahalnobis_merged, outliers) #Remove the outliers




#--------------------------------------------------
###Feature Selection Evaluation
#Unfiltered Features
data_merged[,21:ncol(data_merged)] <- normalize_MAD(data_merged[,21:ncol(data_merged)]) #MAD normalize
unfiltered_heat <- eval_heatmap(data_merged, "Unfiltered Single IBA1+ Objects (n = 7733)", "Unfiltered Morphometric Features (n = 208)")
unfiltered_heat
#png(here("Figures and Tables", "Figure 2", "Extra", "heatmap_unfiltered_features.png"), width = 5, height=4, units = "in", res = 300)
#draw(unfiltered_heat,heatmap_legend_side="bottom")
#dev.off()

#Filtered Features
filtered_heat <- eval_heatmap(Filtered_data, "Filtered Single IBA1+ Objects (n = 7559)", "Filtered Morphometric Features (n = 73)")
filtered_heat
#png(here("Figures and Tables", "Figure 2", "Extra", "heatmap_filtered_features.png"), width = 5, height=4, units = "in", res = 300)
#draw(filtered_heat,heatmap_legend_side="bottom")
#dev.off()



#---------------------------------------------------
#Add Patient ID and Slide Number
#Patient ID (n = 12)
Patient <- Filtered_data %>% mutate(patient=case_when(
  Metadata_Sample_number=="1"~"Patient 1",
  Metadata_Sample_number=="2"~"Patient 1",
  Metadata_Sample_number=="5"~"Patient 2",
  Metadata_Sample_number=="6"~"Patient 2",
  Metadata_Sample_number=="11"~"Patient 3",
  Metadata_Sample_number=="12"~"Patient 3",
  Metadata_Sample_number=="13"~"Patient 4",
  Metadata_Sample_number=="14"~"Patient 4",
  Metadata_Sample_number=="15"~"Patient 5",
  Metadata_Sample_number=="16"~"Patient 5",
  Metadata_Sample_number=="17"~"Patient 6",
  Metadata_Sample_number=="18"~"Patient 6",
  Metadata_Sample_number=="19"~"Patient 7",
  Metadata_Sample_number=="20"~"Patient 7",
  Metadata_Sample_number=="21"~"Patient 8",
  Metadata_Sample_number=="22"~"Patient 8",
  Metadata_Sample_number=="23"~"Patient 9",
  Metadata_Sample_number=="24"~"Patient 9",
  Metadata_Sample_number=="25"~"Patient 10",
  Metadata_Sample_number=="26"~"Patient 10",
  Metadata_Sample_number=="27"~"Patient 11",
  Metadata_Sample_number=="28"~"Patient 11",
  Metadata_Sample_number=="29"~"Patient 12",
  Metadata_Sample_number=="30"~"Patient 12",
)) %>% pull(patient)

ready_singlecell_data <- cbind(Patient, Filtered_data)

#Slide Number (n = 12)
Slide <- Filtered_data %>% mutate(slide=case_when(
  Metadata_Sample_number=="1"~"Slide 4",
  Metadata_Sample_number=="2"~"Slide 5",
  Metadata_Sample_number=="5"~"Slide 6",
  Metadata_Sample_number=="6"~"Slide 7",
  Metadata_Sample_number=="11"~"Slide 1",
  Metadata_Sample_number=="12"~"Slide 8",
  Metadata_Sample_number=="13"~"Slide 2",
  Metadata_Sample_number=="14"~"Slide 9",
  Metadata_Sample_number=="15"~"Slide 3",
  Metadata_Sample_number=="16"~"Slide 10",
  Metadata_Sample_number=="17"~"Slide 6",
  Metadata_Sample_number=="18"~"Slide 11",
  Metadata_Sample_number=="19"~"Slide 7",
  Metadata_Sample_number=="20"~"Slide 12",
  Metadata_Sample_number=="21"~"Slide 1",
  Metadata_Sample_number=="22"~"Slide 8",
  Metadata_Sample_number=="23"~"Slide 2",
  Metadata_Sample_number=="24"~"Slide 9",
  Metadata_Sample_number=="25"~"Slide 3",
  Metadata_Sample_number=="26"~"Slide 10",
  Metadata_Sample_number=="27"~"Slide 4",
  Metadata_Sample_number=="28"~"Slide 11",
  Metadata_Sample_number=="29"~"Slide 5",
  Metadata_Sample_number=="30"~"Slide 12",
)) %>% pull(slide)

ready_singlecell_data <- cbind(Slide, ready_singlecell_data)





#######AGGREGATED DATASET (MORPHOLOGICAL FEATURES FOR AGGREGATED OBJECTS PER IMAGE)
aggregated_data <- weighted_aggregate(ready_singlecell_data)


#Lasly Export for further downstream analysis:
#......Single cell profiles
write.csv(ready_singlecell_data, file=here(main_dir, "Feature selected profiles", "SingleMicroglia_profile.csv"),
          row.names = FALSE) 

#......Aggregated cell profiles
write.csv(aggregated_data, file=here(main_dir, "Feature selected profiles", "AggregatedMicroglia_profile.csv"),
          row.names = FALSE)


sessionInfo()
