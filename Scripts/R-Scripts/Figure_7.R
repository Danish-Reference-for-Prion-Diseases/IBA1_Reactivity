library(here)
library(tidyverse)
library(factoextra)
library(plotly)
library(RColorBrewer)
library(circlize)
library(GGally)
library(gridExtra)
library(ggrepel)
library("ggh4x")
library(cowplot)
library(SCORPIUS)
library(ComplexHeatmap)
library(viridis)
library(MOFA2)

#Colors
pt_colors <- c("Patient 1"="#b15928", "Patient 2" = "#ffff99", "Patient 3" = "#6a3d9a",
               "Patient 4" = "#cab2d6", "Patient 5" = "#ff7f00", "Patient 6" = "#fdbf6f",
               "Patient 7" = "#e31a1c", "Patient 8" = "#fb9a99", "Patient 9" = "#34a02c",
               "Patient 10"="#b2df8a", "Patient 11" = "#1f78b4", "Patient 12" = "#a6cee3")

disease_color <- c("CJD" = "#C6D300", "AD" = "#D6390E", "CTRL" = "#390F6E")

layer_color <- c("WM" = "#56B4E9", "dGM" = "#009E73", "sGM"= "#dbc300")

region_color <- c("FC" = "orange", "OCC" = "#000000") 


###########DATA##############
#Aggregated Cellprofiler Data
Morph_data <- read.csv(here("Data", "Cellprofiler (morphology)", "Feature selected profiles", "AggregatedMicroglia_profile.csv"))
feature_column_start <- 15
cell_features <- Morph_data[,feature_column_start:ncol(Morph_data)]

##Geomx Protein Data
Protein_data <- read.csv(here("Data", "GeoMx (protein)", "3. HK normalized Log2 Transformed GeoMx Data.csv"))
colnames(Protein_data) <- gsub("\\.", " ", colnames(Protein_data))
Protein_data <- Protein_data %>% select(-S6, -GAPDH, -`Histone H3`, -`Myelin basic protein`, -NeuN, -MAP2, -C4B, -Vimentin, -EMP1, 
                                        -`Neurofilament light`, -GFAP, -Synaptophysin, -S100B, -Olig2) # Discard Housekeeper proteins

#Rename for later joining
Morph_data <- Morph_data %>% rename("ROI (Label)" = "ROI..Label.")
Protein_data <- Protein_data %>% rename("ROI (Label)" = "ROI  Label ")


#################MOFA#########################
#Combine datasets (left-join to preserve excess morph data)
combined_data <- left_join(Morph_data, Protein_data, by = c("ROI (Label)", "Patient", "Slide", "Sex", "Age", "Formalinfixed", "Postdelay"))
#Metadata
meta_data <- combined_data[1:c(feature_column_start-1)]

#Prepare View 1 - Proteins
Protein <- combined_data[,118:ncol(combined_data)] # Protein Expression dat 
Protein <- scale(Protein, scale = T, center = T) #Standardize and center
Protein <- as.matrix(t(Protein)) # Transpose

#Prepare View 2 - Morphometrics
Morphology <- combined_data[,c(feature_column_start:87)]
Morphology <- as.matrix(t(Morphology)) #Transpose

#Matrix list
data <- list(Protein,Morphology)

#Create MOFA-object 
MOFAobject <- create_mofa(data)
plot_data_overview(MOFAobject)

#Training-parameters - getting and keeping the default options
data_opts <- get_default_data_options(MOFAobject)
head(data_opts)
data_opts <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
head(model_opts)
train_opts <- get_default_training_options(MOFAobject)
head(train_opts)
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

#Train and save actual MOFA Object
outfile = here("Data", "MOFAmodelmicrominimal.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)

#Load and confirm MOFA Object
model <- load_model(outfile)
plot_data_overview(model)

#Add, rename and reorder actual metadata
meta_data_mofa <- cbind(MOFAobject@samples_metadata[["sample"]],meta_data)
meta_data_mofa <- meta_data_mofa %>% rename("sample" = "MOFAobject@samples_metadata[[\"sample\"]]", 
                                            "Brain Region" = "Metadata_Cortex_region",
                                            "Cortex Layer" = "Metadata_Layer",
                                            "Disease Group" = "Metadata_Pt_type")
meta_data_mofa$`Cortex Layer` <- factor(meta_data_mofa$`Cortex Layer`, levels=c("sGM", "dGM", "WM"))
meta_data_mofa$`Disease Group` <- factor(meta_data_mofa$`Disease Group`, levels=c("CTRL", "AD", "CJD"))
samples_metadata(model) <- meta_data_mofa #save

#VARIANCE DECOMPOSITION
head(model@cache$variance_explained$r2_total[[1]]) # Total R^2
plot_variance_explained(model, x ="view", y ="factor", plot_total = T)#Overview
head(model@cache$variance_explained$r2_per_factor[[1]]) #R^2 per factor
plot_variance_explained(model, x="view", y="factor") + #Per factor
  scale_x_discrete(labels = c(
    expression("Protein ("*R^2*" = 58.3%)"),
    expression("Morphometrics ("*R^2*" = 72.2%)")
  ))


#FEATURE LOADINGS
views <- c("view_1", "view_2")
plots <- list()
# Generate plots and store them in a list
for(i in 1:2) {
  for(j in 1:3) {
    p <- plot_top_weights(model, view = views[i], factor = j, nfeatures = 5
    ) + theme(text = element_text(size = 15))
    plots[[length(plots) + 1]] <- p
  }
}
#Display loadings
dev.new()
grid.arrange(grobs = plots, nrow = 2, ncol = 3)


#Generate new list isolating factor 1
plots <- list()
for (i in 1:2) {
  p <- plot_top_weights(model, view = views[i], factor = 1, nfeatures = 5) +
    theme(text = element_text(size = 12),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"))
  plots[[length(plots) + 1]] <- p
}
#Rename morphometric features before saving
levels(plots[[2]][["data"]][["feature_id"]]) <- c(
  "Soma - Texture Correlation(.3)",
  "Total - Maximum Radius",
  "Total - Mean Radius",
  "Total - Texture Correlation(.2)",
  "Total - Texture Correlation(.3)"
)
# Align the plots to ensure the same width
aligned_plots <- align_plots(plots[[1]], plots[[2]], align = "v", axis = "lr")
grid.arrange(grobs = aligned_plots, nrow = 2, ncol = 1)


#Extra - Check (linear) regression estimate between (top) features and factor values
plot_data_scatter(model, view = "view_1", #Change view
                  factor=1, #Change Factor
                  features=5,
                  add_lm=T,
                  text_size = 0,
                  color_by="Disease Group",
                  shape_by = "Cortex Layer")


plot_data_scatter(model, view = "view_2", #Change view
                  factor=1, #Change Factor
                  features=5,
                  add_lm=T,
                  color_by="Disease Group",
                  shape_by = "Cortex Layer")


#VISUALISATION OF FACTORS
#Overall overview
p <- plot_factor(model, 
                 factors = c(1:3),
                 group_by = "Brain Region",
                 color_by = "Disease Group",
                 shape_by = "Cortex Layer", #Cortex Layer or Brain Region
                 dot_size = 2,        # change dot size
                 dodge = T,           # dodge points with different colors
                 legend = T,          # remove legend
                 add_violin = T,      # add violin plots,
                 violin_alpha = 0.25,  # transparency of violin plots
                 show_missing = T)
p

#Proper Visualisation of Factor 1 in isolation
factors <- get_factors(model, factors = 1, as.data.frame = T) #Extract factor 1
factor_1 <- merge(factors, meta_data_mofa, by="sample") # Add metadata

ggplot(data = factor_1, aes(x = `Cortex Layer`, y = value, color=as.factor(`Disease Group`))) +
  facet_nested(~ "Brain Region" + `Brain Region`)+
  geom_violin(alpha=0.75, position =position_dodge(width=0.7), scale="width", width=0.5)+
  geom_point(aes(x=`Cortex Layer`, y = value, shape = `Disease Group`, group= `Disease Group`, fill=as.factor(Patient)), 
             data = factor_1, size = 1.5, alpha=0.8, stroke = 0.8,
             position = position_jitterdodge(dodge.width = 0.7, jitter.width = 1)) +
  theme_grey()+
  scale_fill_manual(values = pt_colors)+
  scale_color_manual(values = disease_color,
                     guide = guide_legend(override.aes = list(shape = c(21, 24, 22),alpha = 1, size = 1.5)))+
  scale_shape_manual(values = c("CTRL" = 21,  
                                "AD" = 24,    
                                "CJD" = 22))+ 
  labs(x = "Cortex Layer", y = "Value of Factor 1", color = "Disease Group", fill=NULL, shape=NULL) +
  guides(fill = "none",
         shape = "none")+
  theme(text = element_text(size = 12), plot.title = element_text(size = 20),
        strip.text.x = element_text(face = "bold", size = 10),
        legend.title = element_text(size = 8.5), legend.text = element_text(size = 7), legend.key.size = unit(1, "lines"), 
        legend.spacing = unit(1.05, "cm"),legend.margin = margin(t = -20))



#############SCORPIUS - Microglial Reactivity Grading############
##Prepare Data
#Standardize Protein Data
Geomx_Matrix <- Protein_data %>%
  select(1:36) %>%                           # Select the first 36 columns (metadata)
  bind_cols(Protein_data %>%                  # Bind columns with scaled protein expression
              select(38:ncol(Protein_data)) %>%               # Select columns 38 to end for protein expression data
              scale(center=T, scale=T)) 

#Combine Datasets
IBA1_matrix <- merge(Geomx_Matrix, Morph_data, by = "ROI (Label)") #SCORPIUS only works on complete data, thus merging instead of left_join

#Subset Grey matter ROIs
IBA1_matrix <- subset(IBA1_matrix, Metadata_Layer!="WM")

##Isolate feature sets
expression <- as.matrix(cbind(IBA1_matrix[,37:53], IBA1_matrix[,67:ncol(IBA1_matrix)]))#Including protein and Morp features

##Actual SCORPIUS analysis
space <- reduce_dimensionality(expression, "pearson")#Pearson as distance metric, spearman produces stochastic results
Interest_tract <- factor(IBA1_matrix$Disease, levels = c("CTRL", "AD", "CJD"))
BrainRegion <- as.factor(IBA1_matrix$BrainRegion)
Patient <- as.factor(IBA1_matrix$Patient.x)
traj <- infer_trajectory(space, maxit=1000, approx_points = 2) #Keep the trajectory straight

plot_data_combined <- data.frame(
  Dim1 = space[, "Comp1"],  # Component 1 and 2
  Dim2 = space[, "Comp2"],
  Interest_tract = Interest_tract,
  BrainRegion = BrainRegion,
  Patient = Patient
) # Prepare Plot data

draw_trajectory_plot(space, Interest_tract, traj$path, contour = TRUE, progression_group_palette = disease_color)+
  geom_point(data = plot_data_combined, aes(x = Dim1, y = Dim2, shape = BrainRegion, color = Interest_tract, fill=Interest_tract), size = 2.6) +
  scale_shape_manual(values = c(21, 25))+
  guides(shape = guide_legend(title = "Brain Region"), 
         color = guide_legend(title = "Disease Group"),
         fill="none") # Plot

#Feature importance
gimp <- gene_importances(
  expression, 
  traj$time, 
  num_permutations = 10, 
  num_threads = 8, 
  ntree = 10000,
  ntree_perm = 1000)
print(gimp, n=20)

# Select significant features
gene_sel <- gimp$gene[gimp$`pvalue` < .05]
expr_sel <- scale_quantile(expression[, gene_sel]) #Quantile Scaling

# Combine Disease (Interest_tract) and other annotations (metadata)
combined_annotation <- data.frame(
  Disease = Interest_tract,
  Sex = IBA1_matrix$Sex.x,
  Age = IBA1_matrix$Age.x,
  BrainRegion = IBA1_matrix$BrainRegion,
  Layer = IBA1_matrix$CortexLayer,
  Patient = IBA1_matrix$Patient.x,
  Reactivity = traj$time,
  Identifier = IBA1_matrix$`ROI (Label)`)
combined_annotation$Disease <- ordered(as.factor(combined_annotation$Disease), levels = c("CTRL", "AD", "CJD"))
combined_annotation$Layer <- ordered(as.factor(combined_annotation$Layer), levels=c("sGM", "dGM"))

# Rowname and transpose
rownames(combined_annotation) <- rownames(IBA1_matrix)
expr_sel <- t(expr_sel)

# Reorder combined_annotation by Reactivity
combined_annotation <- combined_annotation[order(combined_annotation$Reactivity), ]

# Reorder expr_sel matrix according to the new order of combined_annotation
expr_sel <- expr_sel[, rownames(combined_annotation)]

#change colnames
combined_annotation <- rename(combined_annotation,
                              `Brain Region` = BrainRegion,
                              `Reactivity Score` = Reactivity)

# Colors
annotation_colors <- list(
  "Disease Group" = disease_color,
  "Brain Region" = region_color,
  "Patient ID" = pt_colors,
  "Cortex Layer" = layer_color,
  "Reactivity Score" = colorRamp2(c(min(combined_annotation$Reactivity), max(combined_annotation$Reactivity)), c("white", "red")),
  "Sex" = c("Female"="violet", "Male"="#63C5DA"),
  "Age" = colorRamp2(c(50, 100), c("#99d8c9","#005824")))

# Create Heatmap annotations for columns
col_annotation <- HeatmapAnnotation(
  "Patient ID" = combined_annotation$Patient,
  Sex = combined_annotation$Sex,
  Age = combined_annotation$Age,
  "Disease Group" = combined_annotation$Disease,
  "Brain Region" = combined_annotation$`Brain Region`,
  "Cortex Layer" = combined_annotation$Layer,
  "Reactivity Score" = combined_annotation$`Reactivity Score`,
  col = annotation_colors, show_legend = c(FALSE,TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  annotation_name_gp = gpar(fontface="bold", fontsize = 10))

#Rename Morph RowNames 
rownames(expr_sel)[rownames(expr_sel) == "Soma_Texture_DifferenceEntropy_IBA1_3_02_256"] <- "Soma - Texture Difference Entropy(.2)"
rownames(expr_sel)[rownames(expr_sel) == "Soma_Texture_DifferenceEntropy_IBA1_3_00_256"] <- "Soma - Texture Difference Entropy(.0)"
rownames(expr_sel)[rownames(expr_sel) == "Soma_Texture_InfoMeas1_IBA1_3_03_256"] <- "Soma - Texture Info Meas1(.3)"
rownames(expr_sel)[rownames(expr_sel) == "Total_AreaShape_MeanRadius"] <- "Total - Mean Radius"

#Set Proteins as italics
custom_row_names <- rownames(expr_sel)
italic_rows <- c("IBA1", "CSF1R", "CD11b", "P2ry12", "TMEM119", "CD39", "CD68", "CD11c", "HLA DR")
custom_row_names[rownames(expr_sel) %in% italic_rows] <- 
  sapply(custom_row_names[rownames(expr_sel) %in% italic_rows], 
         function(x) as.expression(bquote(italic(.(x)))))

# Create a heatmap with ComplexHeatmap
ht <- Heatmap(
  expr_sel,                                  
  name = "Rescaled Feature Value",           
  column_title = "ROIs", row_title = "Features",
  top_annotation = col_annotation,          
  cluster_rows = TRUE,                      
  cluster_columns = FALSE,                  
  show_row_names = TRUE,                    
  row_names_gp = gpar(fontsize = 8.5),
  row_labels = custom_row_names, 
  show_column_names = FALSE,                
  col = colorRamp2(breaks = c(0, 0.5, 1), colors = c("blue", "black", "yellow")),               
  heatmap_legend_param = list(title = "Scaled Feature Value",
                              legend_height = unit(60, "cm"),  
                              legend_direction = "horizontal",  
                              title_position = "topcenter") 
)

# Draw the heatmap
draw(ht, heatmap_legend_side="bottom", annotation_legend_side="right", legend_grouping = "original")

#Reverse Trajector
traj <- reverse_trajectory(traj)#Optional!!!!!: Only if the Trajectory needs to be reversed for correct heatmap plot (run it if the heatmap is reversed)


sessionInfo()
