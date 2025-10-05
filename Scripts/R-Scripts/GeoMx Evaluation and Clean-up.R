#Context for the figures: https://nanostring.com/wp-content/uploads/MK2593_GeoMx_Normalization-Protein.pdf

#Libraries
library(here)
library(tidyverse)
library(readxl)
library(gridExtra)
library(TidyDensity)
library(lme4)
library(emmeans)
library(ggh4x)

###Colors
pt_colors <- c("Patient 1"  = "#b15928", "Patient 2"  = "#ffff99",
"Patient 3"  = "#6a3d9a", "Patient 4"  = "#cab2d6",
"Patient 5"  = "#ff7f00", "Patient 6"  = "#fdbf6f",
"Patient 7"  = "#e31a1c", "Patient 8"  = "#fb9a99",
"Patient 9"  = "#34a02c", "Patient 10" = "#b2df8a",
"Patient 11" = "#1f78b4", "Patient 12" = "#a6cee3")

disease_colors <-c("CTRL" = "#390F6E", "AD"   = "#D6390E", "CJD"  = "#C6D300")



#Custom Functions
process_data <- function(dataset, column){
  #Transpose (inverse) rows and columns
  matrix_transposed = t(dataset)
  
  #Split, Rename columns, remove redundant rows and merge again
  matrix_transposed_split1 <- matrix_transposed[,1:column] # Split
  colnames(matrix_transposed_split1) <- matrix_transposed_split1[1,] # Rename
  matrix_transposed_split1 <- matrix_transposed_split1[-c(1,2,3,4),] # Remove rendundant rows
  
  matrix_transposed_split2 <- matrix_transposed[,(column+1):ncol(matrix_transposed)]
  colnames(matrix_transposed_split2) <- matrix_transposed_split2[4,]
  matrix_transposed_split2 <- matrix_transposed_split2[-c(1,2,3,4),]
  
  Data_Matrix <- bind_cols(matrix_transposed_split1, matrix_transposed_split2) # merge back
  
  #Transform variables to correct type (factorize and order incl.)
  Data_Matrix[,(column+2):ncol(matrix_transposed)] <- sapply(Data_Matrix[,32:ncol(matrix_transposed)], as.numeric)
  
  #rename columns
  colnames(Data_Matrix)[4] <- "Celltype"
  
  #convert to tibble
  Data_Matrix <- as_tibble(Data_Matrix)
  
  #Log2 Transform!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Data_Matrix[32:ncol(Data_Matrix)] <- log2(Data_Matrix[32:ncol(Data_Matrix)])
  
  #Add Patient ID (n = 12)
  Patient <- Data_Matrix %>% mutate(patient=case_when(
    PatientID=="1"~"Patient 1",
    PatientID=="2"~"Patient 1",
    PatientID=="5"~"Patient 2",
    PatientID=="6"~"Patient 2",
    PatientID=="11"~"Patient 3",
    PatientID=="12"~"Patient 3",
    PatientID=="13"~"Patient 4",
    PatientID=="14"~"Patient 4",
    PatientID=="15"~"Patient 5",
    PatientID=="16"~"Patient 5",
    PatientID=="17"~"Patient 6",
    PatientID=="18"~"Patient 6",
    PatientID=="19"~"Patient 7",
    PatientID=="20"~"Patient 7",
    PatientID=="21"~"Patient 8",
    PatientID=="22"~"Patient 8",
    PatientID=="23"~"Patient 9",
    PatientID=="24"~"Patient 9",
    PatientID=="25"~"Patient 10",
    PatientID=="26"~"Patient 10",
    PatientID=="27"~"Patient 11",
    PatientID=="28"~"Patient 11",
    PatientID=="29"~"Patient 12",
    PatientID=="30"~"Patient 12",
  )) %>% pull(patient)
  
  Data_Matrix <- cbind(Patient, Data_Matrix)
  Data_Matrix$Patient <- ordered(as.factor(Data_Matrix$Patient), levels = c("Patient 1", "Patient 2", "Patient 3", "Patient 4", "Patient 5", "Patient 6",
                                                                          "Patient 7", "Patient 8","Patient 9","Patient 10", "Patient 11", "Patient 12")) 
  
  #Rename slide.name to slides numbers
  Data_Matrix <- Data_Matrix %>% mutate(Slide=case_when(
    Slide == "4085a382-1f65-40a8-87a1-f9c228a89342"~"Slide 1",
    Slide == "27d61329-b3bb-4730-a6fe-7dd740ed660c"~"Slide 2",
    Slide == "f894f63c-821c-42c8-8e03-d10de60b419b"~"Slide 3",
    Slide == "ec9f36f2-25b4-491c-b96c-0237329057c3"~"Slide 4",
    Slide == "65104fd9-1c11-44f6-9e0b-d7bf6637ec05"~"Slide 5",
    Slide == "c4edbbd4-3bed-4d6b-bb6a-78140c67e8ac"~"Slide 6",
    Slide == "014a54bd-a2e1-4e0b-a305-ab87be868cf0"~"Slide 7",
    Slide == "fe82390f-9b3d-451e-8ac4-a8262d7c3228"~"Slide 8",
    Slide == "a3a16611-8d0e-4d68-a628-a13f903170fe"~"Slide 9",
    Slide == "325461f4-eb15-49b3-a3b8-7cad77c5aeb1"~"Slide 10",
    Slide == "0f9a2903-c87a-45aa-96ef-ed8a94389587"~"Slide 11",
    Slide == "1772d5b7-2964-4f36-826b-7c0e7ac499c3"~"Slide 12"
  ))
  
  
  #Add Age
  Age <- Data_Matrix %>% mutate(age=case_when(
    Patient=="Patient 1"~53,
    Patient=="Patient 2"~56,
    Patient=="Patient 3"~92,
    Patient=="Patient 4"~64,
    Patient=="Patient 5"~57,
    Patient=="Patient 6"~77,
    Patient=="Patient 7"~66,
    Patient=="Patient 8"~78,
    Patient=="Patient 9"~83,
    Patient=="Patient 10"~85,
    Patient=="Patient 11"~59,
    Patient=="Patient 12"~83
  )) %>% pull(age)
  
  Data_Matrix <- cbind(Age, Data_Matrix)
  Data_Matrix$Age <- as.numeric(Data_Matrix$Age)
  
  #Add Sex
  Sex <- Data_Matrix %>% mutate(sex=case_when(
    Patient=="Patient 1"~"Female",
    Patient=="Patient 2"~"Female",
    Patient=="Patient 3"~"Female",
    Patient=="Patient 4"~"Female",
    Patient=="Patient 5"~"Male",
    Patient=="Patient 6"~"Female",
    Patient=="Patient 7"~"Male",
    Patient=="Patient 8"~"Female",
    Patient=="Patient 9"~"Female",
    Patient=="Patient 10"~"Male",
    Patient=="Patient 11"~"Male",
    Patient=="Patient 12"~"Male"
  )) %>% pull(sex)
  
  Data_Matrix <- cbind(Sex, Data_Matrix)
  
  #Add post-mortem delay
  Postdelay <- Data_Matrix %>% mutate(postdelay=case_when(
    Patient=="Patient 1"~3,
    Patient=="Patient 2"~1,
    Patient=="Patient 3"~NA,
    Patient=="Patient 4"~1,
    Patient=="Patient 5"~1,
    Patient=="Patient 6"~14,
    Patient=="Patient 7"~NA,
    Patient=="Patient 8"~1,
    Patient=="Patient 9"~1,
    Patient=="Patient 10"~4,
    Patient=="Patient 11"~1,
    Patient=="Patient 12"~2
  ))%>% pull(postdelay)
  
  Data_Matrix <- cbind(Postdelay, Data_Matrix)
  Data_Matrix$Postdelay <- as.numeric(Data_Matrix$Postdelay)
  
  Formalinfixed <- Data_Matrix %>% mutate(formalinfixed=case_when(
    Patient=="Patient 1"~34,
    Patient=="Patient 2"~35,
    Patient=="Patient 3"~NA,
    Patient=="Patient 4"~43,
    Patient=="Patient 5"~27,
    Patient=="Patient 6"~275,
    Patient=="Patient 7"~182,
    Patient=="Patient 8"~133,
    Patient=="Patient 9"~86,
    Patient=="Patient 10"~119,
    Patient=="Patient 11"~77,
    Patient=="Patient 12"~35
  ))%>% pull(formalinfixed)
  
  Data_Matrix <- cbind(Formalinfixed, Data_Matrix)
  Data_Matrix$Formalinfixed <- as.numeric(Data_Matrix$Formalinfixed)
  
  #Change NB to CTRL
  Data_Matrix <- Data_Matrix %>% mutate(Disease=case_when(
    Disease=="NB"~"CTRL",
    Disease=="AD"~"AD",
    Disease=="CJD"~"CJD")) 
  
  
  #Change TM and BM to sGM and dGM respectively
  Data_Matrix <- Data_Matrix %>% mutate(CortexLayer=case_when(
    CortexLayer=="TM"~"sGM",
    CortexLayer=="BM"~"dGM",
    CortexLayer=="WM"~"WM")) 
  
  # Add proper Sample_ID
  Sample_ID <- substr(Data_Matrix$ROI_ID,1,nchar(Data_Matrix$ROI_ID)-6)
  Data_Matrix <- cbind(Sample_ID, Data_Matrix)
  
  #Factorize Slide for Later
  Data_Matrix$Slide <- ordered(as.factor(Data_Matrix$Slide), levels = c("Slide 1", "Slide 2", "Slide 3",
                                                                                          "Slide 4", "Slide 5", "Slide 6",
                                                                                          "Slide 7", "Slide 8", "Slide 9",
                                                                                          "Slide 10", "Slide 11", "Slide 12"))
    
  return(Data_Matrix)
}


normalization_visualization <- function(feature) {
  ggplot(Data_Matrix_evaluate, aes(x = Disease, y = !!sym(feature))) +  
    facet_wrap(BrainRegion ~ CortexLayer) +
    geom_boxplot(size = 1, aes(color = Disease),
                 position = position_dodge(width = 0.8),
                 outlier.shape = NA) +
    geom_jitter(size = 2, alpha = 0.75, stroke = 0.8,
                aes(x = Disease, fill = Patient),
                position = position_dodge(width = 0.8),
                shape = 21) +
    scale_color_manual(values = disease_colors) +
    scale_fill_manual(values = pt_colors) +
    labs(
      y = "Log2 Protein Expression",
      x = "Disease",
      title = feature
    ) +
    theme(
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 14),  
      axis.text.y   = element_text(size = 14),  
      axis.title.x  = element_text(size = 16),  
      axis.title.y  = element_text(size = 16),  
      plot.title    = element_text(size = 18, hjust = 0),  
      legend.title  = element_text(size = 14),  
      legend.text   = element_text(size = 12),   
      panel.grid.major = element_line(size = 0.5),  
      panel.grid.minor = element_line(size = 0.25),  
      plot.margin   = margin(10, 10, 10, 10)   
    ) +
    guides(fill = "none")
}

pivot_data <- function(data, columns) {
  data %>%
    pivot_longer(
      cols = 38:columns,                    
      names_to = "Protein",            
      values_to = "count"              
    )
}

ROI_boxplots <- function(dataset, name){
  ggplot(dataset, aes(x = Sample_ID, y = count, fill = Slide)) +
    geom_boxplot(size = 1) +  
    labs(x = "Sample ID",             
         y = "Log2 Protein Count",       
         title = name) +     
    scale_fill_brewer(palette = "Paired") +  
    theme_minimal() +                    
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  
      axis.text.y = element_text(size = 14),  
      axis.title.x = element_text(size = 16), 
      axis.title.y = element_text(size = 16),  
      plot.title = element_text(size = 18, hjust = 0.5),  
      legend.title = element_text(size = 14),  
      legend.text = element_text(size = 12),   
      panel.grid.major = element_line(size = 0.5),  
      panel.grid.minor = element_line(size = 0.25),  
      plot.margin = margin(10, 10, 10, 10)   
    )
}


ROI_density <- function(dataset, name){
  ggplot(dataset, aes(x = count, fill = Sample_ID)) +  
    geom_density(alpha = 0.25) +  
    labs(x = "Log2 Protein Count",             
         y = "Density",                        
         title = name) +          
    theme_minimal() +                         
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  
      axis.text.y = element_text(size = 14), 
      axis.title.x = element_text(size = 16),  
      axis.title.y = element_text(size = 16),  
      plot.title = element_text(size = 18, hjust = 0.5),  
      legend.title = element_text(size = 14),  
      legend.text = element_text(size = 12),   
      panel.grid.major = element_line(size = 0.5),  
      panel.grid.minor = element_line(size = 0.25),  
      plot.margin = margin(10, 10, 10, 10)   
    )
}


#########DATA##############################################################################################
#Working Directory
datadir <- here("Data", "GeoMx (protein)")

evaluateData <- read_excel(here(datadir, "1. IBA1 Raw Geomx Dataset.xlsx")) #Raw GeoMx Data
matrixHK <- read_excel(here(datadir, "2. Target Pruned and HK Normalized Dataset.xlsx")) #Housekeeper Normalized Data
MatrixNegative <- read_excel(here(datadir, "Other Normalization Methods", "NegativeRbIgG+MsIgG1_Normalized.xlsx")) #Negative Control Normalized Data
MatrixArea <- read_excel(here(datadir, "Other Normalization Methods", "Area_Normalized.xlsx")) #Area Normalized Data


####Operation - Metadata completion and Log2 Transformation for all previous data files
data_list <- list(evaluateData, matrixHK, MatrixNegative, MatrixArea)
data_list <- lapply(data_list, process_data, column=30)
list2env(setNames(data_list, c("Data_Matrix_evaluate", "processed_HKData", "processed_negativeData", "processed_areaData")), envir = .GlobalEnv)
############################################################################################################


#####EVALUATION#######################################################

#################################Evaluation of AOIs
#Calculate geomeans of HK proteins
Data_Matrix_evaluate$geomean_hk <- rowMeans(Data_Matrix_evaluate[, c("S6", "Histone H3", "GAPDH")], na.rm = TRUE)

#Calculate geomeans of negative control IgGs
Data_Matrix_evaluate$geomean_negative <- rowMeans(Data_Matrix_evaluate[, c("Rb IgG", "Ms IgG1", "Ms IgG2a")], na.rm = TRUE)


###Histograms
#House keepers
hist(Data_Matrix_evaluate$geomean_hk, xlab ="AOIs' log2 HK Geomean (Counts)", main=NULL,  xlim = c(0, 14))


#Negative Controls
hist(Data_Matrix_evaluate$geomean_negative, xlab ="AOIs' log2 IgG Geomean (Counts)", main=NULL, xlim = c(0,7))


###Geomean per AOI
ggplot(Data_Matrix_evaluate, aes(x = reorder(ROI, geomean_hk), y = geomean_hk, fill = Sample_ID)) + 
  geom_bar(stat = "identity") +
  xlab("Individuals AOIs") +
  ylab("Log2 HK Geomean (Counts)") +
  ggtitle("Geomean Housekeepers per AOI") +
  theme_minimal()+
  theme(
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank()  
  )

ggplot(Data_Matrix_evaluate, aes(x = reorder(ROI, geomean_negative), y = geomean_negative, fill = Sample_ID)) + 
  geom_bar(stat = "identity") +
  xlab("Individuals AOIs") +
  ylab("Log2 IgG Geomean (Counts)") +
  ggtitle("Geomean Negative IgGs per AOI") +
  theme_minimal()+
  theme(
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank()  
  )


#HK vs Geomean
ggplot(Data_Matrix_evaluate, aes(x = geomean_hk, y = geomean_negative, color = Sample_ID)) + 
  geom_point(size = 2.5) +
  xlab("Log2 HK Geomean (Counts)") +
  ylab("Log2 Negative Geomean (Counts)") +
  ggtitle("Geomean Housekeepers vs Geomean Negative IgGs") +
  theme_minimal()




#################################Evaluation of Normalization options
#HK genes
normalization_visualization('S6')
normalization_visualization('Histone H3')
normalization_visualization('GAPDH')


#Negative Control IgGs
normalization_visualization('Rb IgG')
normalization_visualization('Ms IgG2a')
normalization_visualization('Ms IgG1')


#Area & Nuclei
Data_Matrix_evaluate$`AOI surface area` <- as.numeric(Data_Matrix_evaluate$`AOI surface area`)
normalization_visualization('AOI surface area')
Data_Matrix_evaluate$`AOI nuclei count` <- as.numeric(Data_Matrix_evaluate$`AOI nuclei count`)
normalization_visualization("AOI nuclei count")


###############Evaluate After Normalization
#ALL To long format
#Evaluation Data
Data_Matrix_evaluate_long <- pivot_data(Data_Matrix_evaluate, 75) #number of total columns of dataset as argument
#HK data
HKData_long <- pivot_data(processed_HKData, 68) #number of total columns of dataset as argument
#Negative Igg Data
Negative_long <- pivot_data(processed_negativeData, 68) #number of total columns of dataset as argument
#Area data
Area_long <- pivot_data(processed_areaData, 68) #number of total columns of dataset as argument

###Boxplots
#Boxplot - Raw Data
raw <- ROI_boxplots(Data_Matrix_evaluate_long, "Raw GeoMx Data")
raw

#Boxplot - Housekeeper Normalized Data
hk_normalized <- ROI_boxplots(HKData_long, "Housekeeper Normalized GeoMx Data")
hk_normalized

#Boxplot - Negative Control Normalized Data
negative_normalized <- ROI_boxplots(Negative_long, "Negative IgG Normalized GeoMx Data")
negative_normalized

#Boxplot - Area Normalized Data
area_normalized <- ROI_boxplots(Area_long, "Area Normalized GeoMx Data")
area_normalized

grid.arrange(raw, hk_normalized, negative_normalized, area_normalized, nrow = 2, ncol = 2)  



###Density Plots
#Density Plot - Raw Data
raw <- ROI_density(Data_Matrix_evaluate_long, "Raw GeoMx Data")
raw

#Density Plot - Normalized Data
hk_normalized <- ROI_density(HKData_long, "Housekeeper Normalized GeoMx Data")
hk_normalized

#Density Plot - Negative Control Data
negative_normalized <- ROI_density(Negative_long, "Negative IgG Normalized GeoMx Data")
negative_normalized

#Density Plot - Area Normalized
area_normalized <- ROI_density(Area_long, "Area Normalized GeoMx Data")
area_normalized

grid.arrange(raw, hk_normalized, negative_normalized, area_normalized, nrow = 2, ncol = 2)  


#######
#From All Evaluation, Housekeeper normalization was chosen to be proceeded with
write.csv(processed_HKData, here(datadir, "3. HK normalized Log2 Transformed GeoMx Data.csv"),
          row.names = FALSE)









#####Extra
#Interaction plots of non-normalized housekeeper proteins
interaction_plot <- function(feature_name) {
  fit_interact <- lmer(Data_Matrix_evaluate[[feature_name]] ~ Disease * BrainRegion * CortexLayer + (1 | Patient/Sample_ID), data = Data_Matrix_evaluate)
  
  print("!!!!!!!!!!!!!!!!!Model Specifications!!!!!!!!!!!!!!")
  print(fit_interact)
  
  windows(width = 15, height = 15)
  print(performance::check_model(fit_interact))
  
  # Reorder for plotting
  interaction_plot_data <- emmip(fit_interact, Disease ~ CortexLayer | BrainRegion, plotit = FALSE, CIs = TRUE)
  interaction_plot_data$Disease <- ordered(as.factor(interaction_plot_data$Disease), levels = c("CTRL", "AD", "CJD"))
  interaction_plot_data$CortexLayer <- ordered(as.factor(interaction_plot_data$CortexLayer), levels = c("sGM", "dGM", "WM"))
  
  # Disease Contrasts
  emm_disease <- emmeans(fit_interact, pairwise ~ Disease | CortexLayer + BrainRegion)
  print("!!!!!!!!!!!DISEASE GROUP CONTRASTS!!!!!!!!!!!!!!!")
  print(emm_disease$contrasts)
  
  # Brain Region Contrasts
  emm_brainregion <- emmeans(fit_interact, pairwise ~ BrainRegion | CortexLayer + Disease)
  print("!!!!!!!!!!!BRAIN REGION CONTRASTS!!!!!!!!!!!!!!!")
  print(emm_brainregion$contrasts)
  
  # Cortex Layer Contrasts
  emm_cortexlayer <- emmeans(fit_interact, pairwise ~ CortexLayer | BrainRegion + Disease)
  print("!!!!!!!!!!!CORTEX LAYERS CONTRASTS!!!!!!!!!!!!!!!")
  print(emm_cortexlayer$contrasts)
  
  # Prepare data for plotting
  data <- Data_Matrix_evaluate
  data$Disease <- ordered(as.factor(data$Disease), levels = c("CTRL", "AD", "CJD"))
  data$CortexLayer <- ordered(as.factor(data$CortexLayer), levels = c("sGM", "dGM", "WM"))
  
  
  pd <- position_dodge(.5)
  p1 <- ggplot(data = interaction_plot_data, aes(x = CortexLayer, y = yvar, 
                                                 group = Disease, linetype = Disease, shape = Disease)) +
    facet_nested(~ "Brain Region" + BrainRegion) +
    geom_point(aes(x = CortexLayer, 
                   y = .data[[feature_name]], fill = as.factor(Patient), color = as.factor(Disease)), 
               data = data, size = 1.5, alpha = 0.8, stroke = 0.5,
               position = position_jitterdodge(dodge.width = 0.5, jitter.width = 1)) +
    geom_point(size = 2.5, fill = "black", color = "black", alpha = 0.6, position = pd) +
    geom_line(alpha = 0.6, color = "black", position = pd) +
    theme_grey() +
    scale_fill_manual(values = pt_colors) +
    scale_color_manual(values = disease_colors,
                       guide = guide_legend(override.aes = list(shape = c(21, 24, 22), alpha = 1, size = 1.5))) +
    scale_shape_manual(values = c("CTRL" = 21, "AD" = 24, "CJD" = 22)) +
    ggtitle(feature_name) +
    labs(x = "Cortex Layer", y = "Log2 Protein Expression (non-normalized)", shape = "Estimated Mean", color = "IBA1 segments", linetype = NULL) +
    guides(fill = "none", shape = guide_legend(order = 2), linetype = guide_legend(order = 3)) +
    theme(text = element_text(size = 12), plot.title = element_text(size = 20),
          strip.text.x = element_text(face = "bold", size = 10), axis.title.y = element_text(size = 8),
          legend.title = element_text(size = 8.5), legend.text = element_text(size = 7), legend.key.size = unit(1, "lines"), 
          legend.spacing = unit(1.05, "cm"), legend.margin = margin(t = -20))
  
  all_contrasts_df <- bind_rows(
    as.data.frame(emm_disease$contrasts) %>% mutate(contrast_type = "Disease"),
    as.data.frame(emm_brainregion$contrasts) %>% mutate(contrast_type = "BrainRegion"),
    as.data.frame(emm_cortexlayer$contrasts) %>% mutate(contrast_type = "CortexLayer")
  ) %>%
    select(contrast_type, contrast, Disease, everything()) %>%
    mutate(significant = ifelse(p.value < 0.05, "Yes", "No"))
  
  return(list(
    plot_interaction = p1,
    contrasts_df = all_contrasts_df
  ))
}

interaction_plot("S6")
interaction_plot("GAPDH")
interaction_plot("Histone H3")





sessionInfo()
