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



#Colors
pt_colors <- c("Patient 1"="#b15928", "Patient 2" = "#ffff99", "Patient 3" = "#6a3d9a",
               "Patient 4" = "#cab2d6", "Patient 5" = "#ff7f00", "Patient 6" = "#fdbf6f",
               "Patient 7" = "#e31a1c", "Patient 8" = "#fb9a99", "Patient 9" = "#34a02c",
               "Patient 10"="#b2df8a", "Patient 11" = "#1f78b4", "Patient 12" = "#a6cee3")

disease_color <- c("CJD" = "#C6D300", "AD" = "#D6390E", "CTRL" = "#390F6E")

layer_color <- c("WM" = "#56B4E9", "dGM" = "#009E73", "sGM"= "#dbc300")

region_color <- c("FC" = "orange", "OCC" = "#000000")  


###########DATA##############
#Working Directory
#setwd(here("Figures and Tables", "Figure 1"))


#Datasets
Data_Matrix <- read.csv(here("Data", "GeoMx (protein)", "3. HK normalized Log2 Transformed GeoMx Data.csv"))
colnames(Data_Matrix) <- gsub("\\.", " ", colnames(Data_Matrix))

#Discard Housekeeper proteins
Data_Matrix <- Data_Matrix %>% select(-S6, -GAPDH, -`Histone H3`)

#Features and Metadata
features <- Data_Matrix[,36:ncol(Data_Matrix)]
labels <- c("Celltype", "Disease", "BrainRegion", "CortexLayer", "Patient", "Slide", "Sex", "Age")

#############################


#####Variance check
variance_per_protein <- sort(apply(features, 2, function(x) sd(x) / mean(x)), decreasing = T)
variance_per_protein


######Variance Partioning
#Features and metadata
t_features <- t(Data_Matrix[36:ncol(Data_Matrix)])
metadata <- Data_Matrix[1:33]

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

#ggsave("Extra/Varplotoverview.png", height=4, width = 8, device="png", dpi=300)

#Variance Partioning of all proteins
proteins_var <- plotPercentBars(varPart[1:17, ], 
                                col =colVar)+theme(
                                  legend.position = "right", 
                                  legend.box = "horizontal",
                                  legend.title = element_blank(),
                                  legend.key.size = unit(1, "lines")
                                )+coord_flip()
proteins_var
#ggsave("Varplotsingleproteins.png", height=3, width = 5, device="png", dpi=300)

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

#Plot of variables loading
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
    legend.text = element_text(size = 10)
  )


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
#ggsave("loading1+2PCA.png", width = 7, height = 7, device='png', dpi=300)



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
  r = 136,
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
  tickfont = list(size=16)
)

# Define the aspect ratio for square plots
aspect_ratio_settings <- list(
  x = 1, y = 1                     
)

# PCA2D - Patient ID
fig_pca2d_1 <- plot_ly(components, type = "scatter", x = ~PC1, y = ~PC2, color = ~as.factor(Patient), 
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
fig_pca2d_2 <- plot_ly(components, type = "scatter", x = ~PC1, y = ~PC2, color = ~as.factor(Disease), 
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
fig_pca2d_3 <- plot_ly(components, type = "scatter", x = ~PC1, y = ~PC2, color = ~as.factor(BrainRegion), 
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
fig_pca2d_4 <- plot_ly(components, type = "scatter", x = ~PC1, y = ~PC2, color = ~as.factor(CortexLayer), 
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
  fig_pca2d_1,
  fig_pca2d_2,
  fig_pca2d_3,
  fig_pca2d_4,
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
features <- scale(heat_data[,36:ncol(heat_data)], center = TRUE, scale = TRUE)
features <- t(features)

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
htmp <- Heatmap(features, name="z-score", 
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
#png("protein_heatmap.png", width = 18, height=8.5, units = "in", res = 300)
#draw(htmp, heatmap_legend_side="bottom", annotation_legend_side="right", legend_grouping = "original")
#dev.off()



#######Interaction Plots of single Proteins
#Plot function
interaction_plot <- function(feature_name) {
  fit_interact <- lmer(Data_Matrix[[feature_name]] ~ Disease * BrainRegion * CortexLayer + (1 | Patient/Sample_ID), data = Data_Matrix)
  
  print("!!!!!!!!!!!!!!!!!Model Specifications!!!!!!!!!!!!!!")
  print(fit_interact)
  
  windows(width = 15, height = 15)
  print(performance::check_model(fit_interact))
  
  #Reorder
  interaction_plot <- emmip(fit_interact, Disease ~ CortexLayer | BrainRegion, plotit = FALSE, CIs = TRUE)
  interaction_plot$Disease <- ordered(as.factor(interaction_plot$Disease), levels = c("CTRL", "AD", "CJD"))
  interaction_plot$CortexLayer <- ordered(as.factor(interaction_plot$CortexLayer), levels = c("sGM", "dGM", "WM"))
  
  #Disease Contrasts
  emm_disease <- emmeans(fit_interact, pairwise ~ Disease | CortexLayer + BrainRegion)
  print("!!!!!!!!!!!DISEASE CONTRASTS!!!!!!!!!!!!!!!")
  print(emm_disease$contrasts)
  
  #Brain Region Contrasts
  emm_brainregion <- emmeans(fit_interact, pairwise ~ BrainRegion | CortexLayer + Disease)
  print("!!!!!!!!!!!BRAIN REGION CONTRASTS!!!!!!!!!!!!!!!")
  print(emm_brainregion$contrasts)
  
  #Cortex Layer Contrasts
  emm_cortexlayer <- emmeans(fit_interact, pairwise ~ CortexLayer | BrainRegion + Disease)
  print("!!!!!!!!!!!CORTEX LAYERS CONTRASTS!!!!!!!!!!!!!!!")
  print(emm_cortexlayer$contrasts)
  
  
  #Data for single observations
  data <- Data_Matrix
  data$Disease <- ordered(as.factor(data$Disease), levels = c("CTRL", "AD", "CJD"))
  data$CortexLayer <- ordered(as.factor(data$CortexLayer), levels = c("sGM", "dGM", "WM"))
  
  
  #Graph settings
  y_limit <- max(data[[feature_name]], na.rm = TRUE)
  
  dev.new()
  
  #Create ggplot
  pd <- position_dodge(.5)
  ggplot(data = interaction_plot, aes(x = CortexLayer, y = yvar, 
                                      group = Disease, linetype = Disease, shape = Disease)) +
    facet_nested(~ "Brain Region" + BrainRegion)+
    geom_point(aes(x = CortexLayer, 
                   y = .data[[feature_name]], fill=as.factor(Patient), color=as.factor(Disease)), 
               data = data, size = 1.5, alpha=0.8, stroke = 0.5,
               position = position_jitterdodge(dodge.width = 0.5, jitter.width = 1)) +
    geom_point(size = 2.5, fill ="black", color="black", alpha = 0.6, position = pd) +
    geom_line(alpha=0.6, color = "black", position = pd) +
    theme_grey()+
    scale_fill_manual(values = pt_colors)+
    scale_color_manual(values = disease_color,
                       guide = guide_legend(override.aes = list(shape = c(21, 24, 22),alpha = 1, size = 1.5)))+
    scale_shape_discrete(name = "Estimated Mean",
                         breaks = c("CTRL", "AD", "CJD"),
                         labels = c("CTRL", "AD", "CJD")) +
    scale_shape_manual(values = c("CTRL" = 21,  
                                  "AD" = 24,    
                                  "CJD" = 22))+ 
    ggtitle(feature_name) +
    ylim(NA, y_limit) +  # Set y-axis limit dynamically based on the protein data
    labs(x = "Cortex Layer", y = "Log2 Protein Expression", shape = "Estimated Mean", color = "IBA1 segments", linetype=NULL) +
    guides(fill = "none", 
           shape = guide_legend(order = 2),
           linetype = guide_legend(order = 3))+
    theme(text = element_text(size = 12), plot.title = element_text(size = 20),
          strip.text.x = element_text(face = "bold", size = 10),
          legend.title = element_text(size = 8.5), legend.text = element_text(size = 7), legend.key.size = unit(1, "lines"), 
          legend.spacing = unit(1.05, "cm"),legend.margin = margin(t = -20))
}

#Plot - single proteins
interaction_plot("CD11c")

######Mass save
#protein_list = c("P2ry12", "Ki 67", "TMEM119", "CD11c", "CD39", "IBA1", "CD45", "CD68", "CD40", "HLA DR", "CD11b", "CSF1R", "CD31", "CTSD", "CD9", "MERTK", "GPNMB")
#for (protein in protein_list) {
#  interaction_plot(protein)
#  ggsave(paste0(protein, ".png"), path = here("Figures and Tables", "Figure 1", "Single Proteins"), 
#         width = 4.5, height = 3.5, device = 'png', dpi = 300)
#}


sessionInfo()