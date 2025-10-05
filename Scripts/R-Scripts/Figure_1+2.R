#Libraries
library(here)
library(tidyverse)
library(factoextra)
library(plotly)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(GGally)
library(gridExtra)
library("variancePartition")
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
library(ggrepel)
library("ggh4x")
library("glmmSeq")
library("kableExtra")
library("openxlsx")


#Colors
pt_colors <- c("Patient 1"="#b15928", "Patient 2" = "#ffff99", "Patient 3" = "#6a3d9a",
               "Patient 4" = "#cab2d6", "Patient 5" = "#ff7f00", "Patient 6" = "#fdbf6f",
               "Patient 7" = "#e31a1c", "Patient 8" = "#fb9a99", "Patient 9" = "#34a02c",
               "Patient 10"="#b2df8a", "Patient 11" = "#1f78b4", "Patient 12" = "#a6cee3")

disease_color <- c("CJD" = "#C6D300", "AD" = "#D6390E", "CTRL" = "#390F6E")

layer_color <- c("WM" = "#56B4E9", "dGM" = "#009E73", "sGM"= "#dbc300")

region_color <- c("FC" = "orange", "OCC" = "#000000") 


####Protein Data - GeoMx#############
#Datasets
Data_Matrix_Full <- read.csv(here("Data", "GeoMx (protein)", "3. HK normalized Log2 Transformed GeoMx Data.csv"))
colnames(Data_Matrix_Full) <- gsub("\\.", " ", colnames(Data_Matrix_Full))

#Discard Housekeeper proteins and non-specifc microglia proteins
Data_Matrix <- Data_Matrix_Full %>% select(-S6, -GAPDH, -`Histone H3`, -`Myelin basic protein`, -NeuN, -MAP2, -C4B, -Vimentin, -EMP1, 
                                           -`Neurofilament light`, -GFAP, -Synaptophysin, -S100B, -Olig2)

#Features and Metadata
features <- Data_Matrix[,38:ncol(Data_Matrix)]#features
labels <- c("Celltype", "Disease", "BrainRegion", "CortexLayer", "Patient", "Slide", "Sex", "Age", "Formalinfixed", "Postdelay")

#############################


#####Variance check
variance_per_protein <- sort(apply(features, 2, function(x) sd(x) / mean(x)), decreasing = T)
variance_per_protein


######Variance partitioning
#Features and metadata
t_features <- t(features) # Transposed Features
metadata <- Data_Matrix[1:36] #metadata

#####ACTUAL
#Model
reduced <- ~ (1|BrainRegion) + (1|CortexLayer) + (1|Disease)+ (1|Disease:BrainRegion) + (1|Disease:CortexLayer) + (1|Sex) + Age

#Partioning
varPart <- fitExtractVarPartModel(t_features, reduced, metadata)

#Rename and reorder
colnames(varPart) <- c("Brain Region", "Cortex Layer", "Disease", 
                       "Disease:Brain Region", "Disease:Cortex Layer", "Sex", "Age", "Residuals")
varPart <- varPart[order((varPart$Residuals), decreasing = F), ]
colVar <- c("#0044CC", "#FFD700", "#FF4500","#8B0000", "orange", "#9467BD", "#87CEFA", 
            "lightgrey")

#Variance Partioning Variables Overview
variables_var <- plotVarPart(varPart, col =colVar) + ylab("Variance Explained (%)")
variables_var

#Variance Partioning of all proteins
proteins_var <- plotPercentBars(varPart[1:17, ], 
                                col =colVar)+theme(
                                  legend.position = "right", 
                                  legend.box = "horizontal",
                                  legend.title = element_blank(),
                                  legend.key.size = unit(1, "lines"))+coord_flip()
proteins_var


#Collinearity Check
form <- ~ BrainRegion + CortexLayer + Disease + Sex + Age
C <- canCorPairs(form, metadata)
plotCorrMatrix(C)


####PCA ANALYSIS
###Compute PCA
prin_comp <- prcomp(features, scale = TRUE)
explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]
explained_variance_ratio <- 100 * explained_variance_ratio
components <- prin_comp[["x"]]
components <- data.frame(components)
components <- cbind(components,Data_Matrix[,labels]) #add labels

###Initial Exploration
#scree plot
fviz_eig(prin_comp, ncp = 40)

#Scatter matrix
interest <- components$CortexLayer #Variable of Interest
ggpairs(components[, 1:5], 
        legend = 1,
        aes(color = as.factor(interest), fill = as.factor(interest)), 
        upper = list(continuous = "points"),  
        lower = list(continuous = "points"),  
        diag = list(continuous = "barDiag")) +  
  scale_color_manual(values = c("#0072B2", "#E69F00", "#999999")) +  
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#999999")) +  
  theme(legend.position = "bottom")  

#Plot of feature loading
loading <- fviz_pca_var(prin_comp, axes = c(1, 2), title="", xlim = c(-0.4,1.2), ylim = c(-1.5,1.5), col.circle = rgb(0, 0, 0, alpha = 0),
                        col.var = "contrib", # Color by contributions to the PC
                        gradient.cols = c("#00AFBB", "#FC4E07"),
                        repel = TRUE , geom = "arrow") + 
  labs(color="Loading") + xlab(paste('PC 1 (',toString(round(explained_variance_ratio[1],1)),'%)',sep = '')) + 
  ylab(paste('PC 2 (',toString(round(explained_variance_ratio[2],1)),'%)',sep = ''))+
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10))


#Inclusion of Labels
label_df <- data.frame(x = loading$data$x, y = loading$data$y, label = loading$data$name)
label_df$contrib <- loading$data$contrib
loading_with_labels <- loading + 
  geom_text_repel(data = label_df, 
                  aes(x = x, y = y, label = label, colour = contrib), 
                  size = log(abs(label_df$contrib))+4, direction = "both", nudge_x = label_df$x * 0.2,  # Nudge labels away from the x-center
                  nudge_y = label_df$y * 0.3,
                  segment.linetype = "dotted",
                  show.legend = FALSE, 
                  max.overlaps = 10) +
  scale_colour_gradientn(colours = c("#00AFBB", "#FC4E07"), 
                         limits = c(0, 10), 
                         breaks = c(0, 5, 10)) +
  scale_x_continuous(limits = c(-0.5, 1.2))
loading_with_labels


####2D PCA Proper
#Reorder factors for plotting
components$Disease <- ordered(as.factor(components$Disease), levels = c("CJD", "AD", "CTRL"))
components$BrainRegion <- ordered(as.factor(components$BrainRegion), levels = c("OCC", "FC")) 
components$CortexLayer <- ordered(as.factor(components$CortexLayer), levels = c("WM", "dGM", "sGM")) 
components$Patient <- ordered(as.factor(components$Patient), levels = c("Patient 12", "Patient 11", "Patient 10", "Patient 9", "Patient 8", "Patient 7",
                                                                          "Patient 6", "Patient 5","Patient 4","Patient 3", "Patient 2", "Patient 1")) 
components$Slide <- ordered(as.factor(components$Slide), levels = c("Slide 12", "Slide 11", "Slide 10", "Slide 9", "Slide 8", "Slide 7",
                                                                      "Slide 6", "Slide 5","Slide 4","Slide 3", "Slide 2", "Slide 1"))
m <- list(
  l = 100,
  r = 150,
  t = 100,
  b = 100)
tit = ""

# Define axis settings to have a fixed range from -3 to 3 and ensure square aspect ratio
axis_settings <- list(
  range = c(-7.5, 7.5),                
  zerolinecolor = "#ffff",         
  zerolinewidth = 2,               
  gridcolor = '#ffff',             
  fixedrange = TRUE,
  dtick = 2,
  tickfont = list(size=16))

# Define the aspect ratio for square plots
aspect_ratio_settings <- list(x = 1, y = 1)

# PCA2D - Patient ID
fig_pca2d_pt <- plot_ly(components, type = "scatter", x = ~PC1, y = ~PC2, color = ~as.factor(Patient), 
                       colors = pt_colors,
                       mode = "markers", legendgroup = "Patient ID", marker = list(size = 5.5, symbol = "circle"), 
                       showlegend = FALSE) %>% 
  layout(xaxis = axis_settings, 
         yaxis = axis_settings, 
         autosize = FALSE, 
         width = 300, height = 300,   
         scene = list(aspectratio = aspect_ratio_settings)) %>%  
  add_annotations(text = "<b>Patient ID</b>", x = 0.5, y = 1.1, yref = "paper", xref = "paper", 
                  xanchor = "center", yanchor = "top", showarrow = FALSE, font = list(size = 20))

# PCA2D - Disease
fig_pca2d_disease <- plot_ly(components, type = "scatter", x = ~PC1, y = ~PC2, color = ~as.factor(Disease), 
                       colors = disease_color, mode = "markers", 
                       legendgroup = "Disease", marker = list(size = 5.5)) %>% 
  layout(xaxis = axis_settings, 
         yaxis = axis_settings, 
         autosize = FALSE, 
         width = 300, height = 300,   
         scene = list(aspectratio = aspect_ratio_settings)) %>%  
  add_annotations(text = "<b>Disease Group</b>", x = 0.5, y = 1.1, yref = "paper", xref = "paper", 
                  xanchor = "center", yanchor = "top", showarrow = FALSE, font = list(size = 20))

# PCA2D - Brain Region 
fig_pca2d_brainregion <- plot_ly(components, type = "scatter", x = ~PC1, y = ~PC2, color = ~as.factor(BrainRegion), 
                       colors = region_color, mode = "markers", legendgroup = "Brain Region", 
                       marker = list(size = 5.5)) %>%
  layout(xaxis = axis_settings, 
         yaxis = axis_settings, 
         autosize = FALSE, 
         width = 300, height = 300,   
         scene = list(aspectratio = aspect_ratio_settings)) %>%  
  add_annotations(text = "<b>Brain Region</b>", x = 0.5, y = 1.1, yref = "paper", xref = "paper", 
                  xanchor = "center", yanchor = "top", showarrow = FALSE, font = list(size = 20))

# PCA2D - Cortex Layer
fig_pca2d_cortexlayer <- plot_ly(components, type = "scatter", x = ~PC1, y = ~PC2, color = ~as.factor(CortexLayer), 
                       colors = layer_color, mode = "markers", 
                       legendgroup = "Cortex Layer", marker = list(size = 5.5)) %>%
  layout(xaxis = axis_settings, 
         yaxis = axis_settings, 
         autosize = FALSE, 
         width = 300, height = 300,   
         scene = list(aspectratio = aspect_ratio_settings)) %>%  
  add_annotations(text = "<b>Cortex Layer</b>", x = 0.5, y = 1.1, yref = "paper", xref = "paper", 
                  xanchor = "center", yanchor = "top", showarrow = FALSE, font = list(size = 20))

# Main 2D PCA PLOT
fig_pca2d <- subplot(
  fig_pca2d_pt,
  fig_pca2d_disease,
  fig_pca2d_brainregion,
  fig_pca2d_cortexlayer,
  nrows = 2, margin = c(0.05,0.05,0.05,0.05)
) %>% 
  layout(
    title = list(text = tit, font = list(size = 24), x=0.47),  
    plot_bgcolor = '#e5ecf6',
    margin = m,
    autosize = FALSE,
    xaxis = axis_settings,   
    yaxis = axis_settings,
    width = 850, height = 850,  
    legend = list(tracegroupgap = 40, y = 0.63, yanchor = "top", itemsizing = 'constant', font = list(size = 13)),
    scene = list(aspectratio = aspect_ratio_settings)  
  ) %>%
  config(
    toImageButtonOptions = list(
      format = "svg",
      filename = "myplot",
      width = NULL,
      height = NULL
    )
  ) %>%
  add_annotations(text = "<b>Disease Group<b>", xref = "paper", yref = "paper",
                  x = 1.02, xanchor = "left",
                  y = 0.622, yanchor = "bottom",
                  legendtitle = TRUE, font = list(size=15), showarrow = FALSE) %>%
  add_annotations(text = "<b>Brain Region<b>", xref = "paper", yref = "paper",
                  x = 1.02, xanchor = "left",
                  y = 0.47, yanchor = "bottom",
                  legendtitle = TRUE, font = list(size=15), showarrow = FALSE) %>%
  add_annotations(text = "<b>Cortex Layer<b>", xref = "paper", yref = "paper",
                  x = 1.02, xanchor = "left",
                  y = 0.35, yanchor = "bottom",
                  legendtitle = TRUE, font = list(size=15), showarrow = FALSE) %>%
  add_annotations(text = paste('PC 1 (', toString(round(explained_variance_ratio[1], 1)), '%)', sep = ''), 
                  xref = 'paper', yref = 'paper', x = 0.5, y = -0.1, font = list(color = "black", size = 24),
                  showarrow = F) %>%
  add_annotations(text = paste('PC 2 (', toString(round(explained_variance_ratio[2], 1)), '%)', sep = ''), 
                  xref = 'paper', yref = 'paper', x = -0.11, y = 0.5, font = list(color = "black", size = 24),
                  textangle = 270, showarrow = F)
fig_pca2d


#Extra PCAs
#PCA2D - Sex
fig_pca2d_sex <- plot_ly(components, type = "scatter", x = ~PC1, y = ~PC2, color = ~as.factor(Sex), 
                       colors = c("violet", "#63C5DA"),
                       mode = "markers", legendgroup = "Sex", marker = list(size = 5.5, symbol = "circle")) %>% 
  layout(xaxis = axis_settings, 
         yaxis = axis_settings, 
         autosize = FALSE, 
         width = 300, height = 300, 
         scene = list(aspectratio = aspect_ratio_settings)) %>%  
  add_annotations(text = "<b>Sex</b>", x = 0.5, y = 1.1, yref = "paper", xref = "paper", 
                  xanchor = "center", yanchor = "top", showarrow = FALSE, font = list(size = 20))

# PCA2D - Age
fig_pca2d_age <- plot_ly(components, type = "scatter", x = ~PC1, y = ~PC2, 
                      mode = "markers", showlegend = FALSE,
                       legendgroup = "Age", marker = list(size = 5.5, color = ~Age, colorscale = "Viridis",
                                                          colorbar = list(title = "Age", y = 0.6, len = 0.2))) %>% 
  layout(xaxis = axis_settings, 
         yaxis = axis_settings, 
         autosize = FALSE, 
         width = 300, height = 300,   
         scene = list(aspectratio = aspect_ratio_settings)) %>%  
  add_annotations(text = "<b>Age</b>", x = 0.5, y = 1.1, yref = "paper", xref = "paper", 
                  xanchor = "center", yanchor = "top", showarrow = FALSE, font = list(size = 20))

# PCA2D - Formalin fixation duration
formalin_non_na <- components[!is.na(components$Formalinfixed), ]
formalin_na     <- components[is.na(components$Formalinfixed), ]
fig_pca2d_formalin <- plot_ly() %>%
  add_trace(
    data = formalin_non_na,
    type = "scatter",
    x = ~PC1,
    y = ~PC2,
    mode = "markers",
    legendgroup = "Formalin Fixation",
    showlegend = FALSE,
    marker = list(
      size = 5.5, 
      color = ~Formalinfixed, 
      colorscale = "Viridis",
      colorbar = list(
        title = "Formalin",
        y = 0.4,        
        len = 0.2))) %>%
  add_trace(
    data = formalin_na,
    type = "scatter",
    x = ~PC1,
    y = ~PC2,
    showlegend = FALSE,
    mode = "markers",
    marker = list(size = 5.5, color = "darkgrey", symbol = "x"),
    name = "Missing Formalin Fixation") %>%
  layout(
    xaxis = axis_settings,
    yaxis = axis_settings,
    autosize = FALSE,
    width = 300, height = 300,
    scene = list(aspectratio = aspect_ratio_settings)) %>%
  add_annotations(text = "<b>Formalin Fixation Duration</b>", x = 0.5, y = 1.1, yref = "paper", xref = "paper",
                  xanchor = "center", yanchor = "top", showarrow = FALSE, font = list(size = 20))

# PCA2D - Postmortem delay
postmortem_non_na <- components[!is.na(components$Postdelay), ]
postmortem_na     <- components[is.na(components$Postdelay), ]
fig_pca2d_postmortem <- plot_ly() %>%
  add_trace(
    data = postmortem_non_na,
    type = "scatter",
    x = ~PC1,
    y = ~PC2,
    mode = "markers",
    showlegend = FALSE,
    marker = list(
      size = 5.5, 
      color = ~Postdelay,
      colorscale = "Viridis",  
      colorbar = list(
        title = "Postmortem",
        y = 0.2,        # middle position between Age and Formalin
        len = 0.2))) %>%
  add_trace(
    data = postmortem_na,
    type = "scatter",
    x = ~PC1,
    y = ~PC2,
    mode = "markers",
    showlegend = FALSE,
    marker = list(size = 5.5, color = "darkgrey", symbol = "x"),
    name = "Missing Postmortem Delay") %>%
  layout(
    xaxis = axis_settings,
    yaxis = axis_settings,
    autosize = FALSE,
    width = 300, height = 300,
    scene = list(aspectratio = aspect_ratio_settings)) %>%
  add_annotations(text = "<b>Postmortem Delay</b>", x = 0.5, y = 1.1, yref = "paper", xref = "paper",
                  xanchor = "center", yanchor = "top", showarrow = FALSE, font = list(size = 20))

# Main 2D PCA PLOT
fig_pca2d_extra <- subplot(
  fig_pca2d_sex,
  fig_pca2d_age,
  fig_pca2d_formalin,
  fig_pca2d_postmortem,
  nrows = 2, margin = c(0.05,0.05,0.05,0.05)
) %>% 
  layout(
    title = list(text = tit, font = list(size = 24), x=0.47),  
    plot_bgcolor = '#e5ecf6',
    margin = m,
    autosize = FALSE,
    xaxis = axis_settings,   
    yaxis = axis_settings,
    width = 850, height = 850,  
    legend = list(tracegroupgap = 40, y = 0.8, yanchor = "top", itemsizing = 'constant', font = list(size = 13)),
    scene = list(aspectratio = aspect_ratio_settings)  
  ) %>%
  config(
    toImageButtonOptions = list(
      format = "png",
      filename = "myplot",
      width = NULL,
      height = NULL
    )
  ) %>%
  add_annotations(text = paste('PC 1 (', toString(round(explained_variance_ratio[1], 1)), '%)', sep = ''), 
                  xref = 'paper', yref = 'paper', x = 0.5, y = -0.1, font = list(color = "black", size = 24),
                  showarrow = F) %>%
  add_annotations(text = paste('PC 2 (', toString(round(explained_variance_ratio[2], 1)), '%)', sep = ''), 
                  xref = 'paper', yref = 'paper', x = -0.11, y = 0.5, font = list(color = "black", size = 24),
                  textangle = 270, showarrow = F)%>%
  add_annotations(text = "Sex", xref = "paper", yref = "paper",
                  x = 1.03, xanchor = "left",
                  y = 0.8, yanchor = "bottom",
                  legendtitle = TRUE, font = list(size=14), showarrow = FALSE)
fig_pca2d_extra


##########Hierarchical Clustering and Heatmap
#Prepare Data
#Reorder Metadata
heat_data <- Data_Matrix
heat_data$Disease <- ordered(as.factor(heat_data$Disease), levels = c("CTRL", "AD", "CJD"))
heat_data$BrainRegion <- ordered(as.factor(heat_data$BrainRegion), levels = c("FC", "OCC"))
heat_data$CortexLayer <- ordered(as.factor(heat_data$CortexLayer), levels = c("sGM", "dGM", "WM")) 
heat_data$Patient <- ordered(as.factor(heat_data$Patient), levels = c("Patient 1", "Patient 2", "Patient 3",
                                                                          "Patient 4", "Patient 5", "Patient 6",
                                                                          "Patient 7", "Patient 8", "Patient 9",
                                                                          "Patient 10", "Patient 11", "Patient 12")) 

#Scale and Transpose the Protein Data
features_heat <- scale(heat_data[,38:ncol(heat_data)], center = TRUE, scale = TRUE)
features_heat <- t(features_heat)

#Defining labels and color palettes (annotation)
annotation_colors <- list("Patient ID" = pt_colors,
                          "Sex" = c("Female"="violet", "Male"="#63C5DA"),
                          "Age" = colorRamp2(c(50, 100), c("#99d8c9","#005824")),
                          "Disease Group" = disease_color,
                          "Brain Region" = region_color,
                          "Cortex Layer" = layer_color)

column_anno <- HeatmapAnnotation(
                        "Patient ID" = heat_data$Patient, "Sex" = heat_data$Sex, "Age" = heat_data$Age, "Disease Group" = heat_data$Disease, 
                        "Brain Region" = heat_data$BrainRegion, "Cortex Layer" = heat_data$CortexLayer, 
                        col = annotation_colors, show_legend = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE), annotation_name_gp = gpar(fontface="bold"))

heatmap_color_u <- colorRamp2(breaks = c(2.1, 0, -2.1), colors = c("yellow", "black", "blue"))

#Heatmap
htmp <- Heatmap(features_heat, name="z-score", 
                column_title = "ROIs", row_title = "Proteins",
                row_names_gp = gpar(fontsize = 11),
                top_annotation = column_anno,
                left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c(3,5)),
                                                                 labels = c("Inflammatory Protein Cluster", "Homeostatic Protein Cluster"), 
                                                                 labels_gp = gpar(col = "white", fontsize = 12, fontface = "bold"))),
                col = heatmap_color_u,
                row_split = 2,
                clustering_method_rows = "ward.D2",
                clustering_method_columns = "ward.D2",
                clustering_distance_rows = "pearson",
                clustering_distance_columns = "pearson",
                heatmap_legend_param = list(legend_height = unit(60, "cm"),  # Adjust the height of the legend
                                            legend_direction = "horizontal",  # Place the legend horizontally
                                            title_position = "topcenter"))
htmp


#######Interaction plots of single proteins
#Plot function
interaction_plot <- function(feature_name) {
  fit_interact <- lmer(Data_Matrix_Full[[feature_name]] ~ Disease * BrainRegion * CortexLayer + (1 | Patient/Sample_ID), data = Data_Matrix_Full)
  
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
  data <- Data_Matrix_Full
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
    scale_color_manual(values = disease_color,
                       guide = guide_legend(override.aes = list(shape = c(21, 24, 22), alpha = 1, size = 1.5))) +
    scale_shape_manual(values = c("CTRL" = 21, "AD" = 24, "CJD" = 22)) +
    ggtitle(feature_name) +
    labs(x = "Cortex Layer", y = "Log2 Protein Expression", shape = "Estimated Mean", color = "IBA1 segments", linetype = NULL) +
    guides(fill = "none", shape = guide_legend(order = 2), linetype = guide_legend(order = 3)) +
    theme(text = element_text(size = 12), plot.title = element_text(size = 20),
          strip.text.x = element_text(face = "bold", size = 10),
          legend.title = element_text(size = 8.5), legend.text = element_text(size = 7), legend.key.size = unit(1, "lines"), 
          legend.spacing = unit(1.05, "cm"), legend.margin = margin(t = -20))
  
  # Plot 2 - Bar plot with CTRL as reference
  emm <- emmeans(fit_interact, ~ Disease | CortexLayer + BrainRegion)
  contrast_ctrl <- contrast(emm, method = "trt.vs.ctrl", ref = "CTRL", infer = c(TRUE, TRUE))
  contrast_df <- as.data.frame(contrast_ctrl) %>%
    filter(contrast %in% c("AD - CTRL", "CJD - CTRL")) %>%
    mutate(Disease = ifelse(contrast == "AD - CTRL", "AD", "CJD"),
           Disease = factor(Disease, levels = c("AD", "CJD")),
           CortexLayer = factor(CortexLayer, levels = c("sGM", "dGM", "WM")))
  
  y_mid_diff <- median(contrast_df$estimate, na.rm = TRUE)
  if (y_mid_diff < 0 || feature_name == "GPNMB") {
    y_lower_diff <- -5
    y_upper_diff <- 2
  } else {
    y_lower_diff <- -2
    y_upper_diff <- 4.5
  }
  
  contrast_df <- contrast_df %>% mutate(TopLevel = "Brain Region")#Dummy column
  
  p2 <- ggplot(contrast_df, aes(x = CortexLayer, y = estimate, fill = Disease)) +
    facet_nested(~ TopLevel + BrainRegion) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", alpha = 0.8) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), position = position_dodge(width = 0.9), width = 0.25) +
    scale_fill_manual(values = disease_color) +
    theme_grey() +
    geom_hline(yintercept = 0, linewidth = 0.5, alpha=0.7, color = "darkred") +
    theme(text = element_text(size = 12),
          plot.title = element_text(size = 20),
          strip.text.x = element_text(face = "bold", size = 10),
          legend.title = element_text(size = 8.5),
          legend.text = element_text(size = 7),
          legend.key.size = unit(1, "lines"),
          legend.spacing = unit(1.05, "cm"),
          legend.margin = margin(t = -20)) +
    labs(x = "Cortex Layer", y = "Difference from CTRL (log2)", fill = "Disease") +
    ggtitle(feature_name) +
    scale_y_continuous(limits = c(y_lower_diff, y_upper_diff), breaks = seq(y_lower_diff, y_upper_diff, by = 1))
  
  all_contrasts_df <- bind_rows(
    as.data.frame(emm_disease$contrasts) %>% mutate(contrast_type = "Disease"),
    as.data.frame(emm_brainregion$contrasts) %>% mutate(contrast_type = "BrainRegion"),
    as.data.frame(emm_cortexlayer$contrasts) %>% mutate(contrast_type = "CortexLayer")
  ) %>%
    select(contrast_type, contrast, Disease, everything()) %>%
    mutate(significant = ifelse(p.value < 0.05, "Yes", "No"))
  
  return(list(
    plot_interaction = p1,
    plot_contrast = p2,
    contrasts_df = all_contrasts_df
  ))
}
  
#Plot - single proteins
interact <- interaction_plot("GPNMB")
print(interact$plot_interaction)
print(interact$plot_contrast)

####Mass save
#protein_list <- colnames(Data_Matrix_Full[,38:ncol(Data_Matrix_Full)])
#contrast_dfs_list <- list()
#for (protein in protein_list) {
  #result <- interaction_plot(protein)
  
  #ggsave(filename = paste0("Graphs", protein, "_interaction_plot.png"), plot = result$plot_interaction, width = 5, height = 3.5)#Save interaction plot
  #ggsave(filename = paste0("Graphs",protein, "_contrast_plot.png"), plot = result$plot_contrast, width = 5, height = 3.5)#Save contrast plot
  
  #df_with_protein <- result$contrasts_df %>%
    #mutate(Protein = protein) %>%
    #select(Protein, everything())
  #contrast_dfs_list[[protein]] <- df_with_protein
#}
#wb <- createWorkbook()
#for (protein in names(contrast_dfs_list)) {
  #addWorksheet(wb, protein)
  #writeData(wb, sheet = protein, contrast_dfs_list[[protein]])
#}
#saveWorkbook(wb, file = "protein_contrasts_all.xlsx", overwrite = TRUE)
#####


###Supplementary
features_full <- Data_Matrix_Full[,38:ncol(Data_Matrix_Full)]

######Variance partitioning - With Postmortem Delay and formalin fixation duration
#Model
full <- ~ (1|BrainRegion) + (1|CortexLayer) + (1|Disease)+ (1|Disease:BrainRegion) + (1|Disease:CortexLayer) + (1|Sex) + Age + Formalinfixed + Postdelay

#Partioning
varPart <- fitExtractVarPartModel(t(features_full), full, metadata)

#Rename and reorder
colnames(varPart) <- c("Brain Region", "Cortex Layer", "Disease", 
                       "Disease:Brain Region", "Disease:Cortex Layer", "Sex", "Age", "Formalin Fixation Duration", "Postmortem Delay", "Residuals")
varPart <- varPart[order((varPart$Residuals), decreasing = F), ]
colVar <- c("#0044CC", "#FFD700", "#FF4500","#8B0000", "orange", "#9467BD", "#87CEFA", "#F5F5F5", "black", "lightgrey")

#Variance Partioning Variables Overview
variables_var <- plotVarPart(varPart, col =colVar) + ylab("Variance Explained (%)")
variables_var

#Variance Partioning of all proteins
proteins_var <- plotPercentBars(varPart[1:31, ], 
                                col =colVar)+theme(
                                  legend.position = "right", 
                                  legend.box = "horizontal",
                                  legend.title = element_blank(),
                                  legend.key.size = unit(1, "lines"))+coord_flip()
proteins_var

#Collinearity Check
form <- ~ BrainRegion + CortexLayer + Disease + Sex + Age + Formalinfixed + Postdelay
C <- canCorPairs(form, metadata)
plotCorrMatrix(C)


###Evaluating effects of Age, Sex, Formalin Fixation Duration, Postmortem Delay on protein expression
#Age
results_age <- glmmQvals(lmmSeq(~ Age + (1|Patient) + (1|Sample_ID),
                                maindata = t(features_full),
                                metadata = Data_Matrix_Full[,1:37],
                                progress = TRUE))
stats <- summary(results_age)
kable(stats[order(stats[, 10]), ]) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "400px")

#Sex
results_sex <- glmmQvals(lmmSeq(~ Sex + (1|Patient) + (1|Sample_ID),
                                maindata = t(features_full),
                                metadata = Data_Matrix_Full[,1:37],
                                progress = TRUE))
stats <- summary(results_sex)
kable(stats[order(stats[, 10]), ]) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "400px")

#Formalin Fixation duration
combined_df <- cbind(Data_Matrix_Full, features_full)
combined_df_filtered <- combined_df[!is.na(combined_df$Formalinfixed), ]
metadata_clean <- combined_df_filtered[, colnames(Data_Matrix_Full)]
features_clean <- combined_df_filtered[, colnames(features_full)]
results_fixation <- glmmQvals(lmmSeq(~ Formalinfixed + (1|Patient) + (1|Sample_ID),
                                maindata = t(features_clean),
                                metadata = metadata_clean[,1:36],
                                progress = TRUE))
stats <- summary(results_fixation)
kable(stats[order(stats[, 10]), ]) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "400px")


#Postmortem Delay
combined_df <- cbind(Data_Matrix_Full, features_full)
combined_df_filtered <- combined_df[!is.na(combined_df$Postdelay), ]
metadata_clean <- combined_df_filtered[, colnames(Data_Matrix_Full)]
features_clean <- combined_df_filtered[, colnames(features_full)]
results_postmortem <- glmmQvals(lmmSeq(~ Postdelay + (1|Patient) + (1|Sample_ID),
                                     maindata = t(features_clean),
                                     metadata = metadata_clean[,1:36],
                                     progress = TRUE))
stats <- summary(results_postmortem)
kable(stats[order(stats[, 10]), ]) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "400px")



sessionInfo()
