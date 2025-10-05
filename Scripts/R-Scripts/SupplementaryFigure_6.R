#Libraries
library(here)
library(tidyverse)
library(factoextra)
library(plotly)
library(RColorBrewer)
library(circlize)
library(GGally)
library(gridExtra)
library("variancePartition")
library(lme4)
library(emmeans)
library(glmmSeq)
library(performance)
library(kableExtra)
library(ggrepel)
library("ggh4x")
library(cowplot)

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


#####AGGREGATED MORPH FEATURES EXPLORATION####################
#Variance Check
variance_per_feature <- sort(apply(cell_features, 2, var), decreasing = T)
variance_per_feature

####PCA ANALYSIS
###Compute PCA
prin_comp <- prcomp(cell_features, scale = FALSE)
explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]
explained_variance_ratio <- 100 * explained_variance_ratio
components <- prin_comp[["x"]]
components <- data.frame(components)
components <- cbind(components,Morph_data[,1:(feature_column_start-1)])


###Initial Exploration
#scree plot
fviz_eig(prin_comp, ncp = 40)

#Scatter matrix
interest <- components$Metadata_Pt_type #Change for different variable
ggpairs(components[, 1:4], 
        legend = 1,
        aes(color = as.factor(interest), fill = as.factor(interest)),
        upper = list(continuous = "points"),  
        lower = list(continuous = "points"),  
        diag = list(continuous = "barDiag")) +  
  scale_color_manual(values = c("#0072B2", "#E69F00", "#999999")) +  
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#999999")) + 
  theme(legend.position = "bottom")  

#Plot of variables loading
loading <- fviz_pca_var(prin_comp, axes = c(1, 3), title="", xlim = c(-0.2,0.4), ylim = c(-0.3,0.3),
                        col.var = "contrib", 
                        gradient.cols = c("#00AFBB", "#FC4E07"),
                        repel = TRUE , geom = "arrow"
                        
) + labs(color="Loading") + xlab(paste('PC 1 (',toString(round(explained_variance_ratio[1],1)),'%)',sep = '')) + 
  ylab(paste('PC 3 (',toString(round(explained_variance_ratio[3],1)),'%)',sep = ''))+
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
loading

#Inclusion of Labels
label_df <- data.frame(x = loading$data$x, y = loading$data$y, label = loading$data$name)
label_df$contrib <- loading$data$contrib

#Isolate the 12 most contributing features
top_label_df <- label_df[order(-label_df$contrib), ][1:12, ]
#top_label_df <- top_label_df[-11,] #For better visualization in publication

#Rename
top_label_df <- top_label_df %>% mutate(label = recode(label,
                "Total_AreaShape_MaximumRadius" = "Total - Maximum Radius",
                "Total_AreaShape_Zernike_2_0" = "Total - Zernike 2.0",
                "Total_AreaShape_MeanRadius" = "Total - Mean Radius",
                "Total_AreaShape_EquivalentDiameter" = "Total - Equivalent Diameter",
                "Total_AreaShape_Extent" = "Total - Extent",
                "Total_Texture_SumEntropy_IBA1_3_02_256" = "Total - Texture Sum Entropy(.2)",
                "Total_AreaShape_MajorAxisLength" = "Total - Major Axis Length",
                "Total_Texture_Correlation_IBA1_3_00_256" = "Total - Texture Correlation(.0)",
                "Total_Texture_Correlation_IBA1_3_01_256" = "Total - Texture Correlation(.1)",
                "Total_Texture_Correlation_IBA1_3_03_256" = "Total - Texture Correlation(.3)",
                "Total_Texture_Correlation_IBA1_3_02_256" = "Total - Texture Correlation(.2)",
                "Soma_Texture_Correlation_IBA1_3_02_256" = "Soma - Texture Correlation(.2)"
))

loading_with_labels <- loading + 
  geom_text_repel(data = top_label_df, 
                  aes(x = x, y = y, label = label, colour = contrib), 
                  size = log(abs(top_label_df$contrib))+2.75, direction = "both",
                  segment.linetype = "dotted", nudge_y = ifelse(top_label_df$y>0, +0.05, -0.025), nudge_x = 0.02,
                  show.legend = FALSE) +
  scale_colour_gradientn(colours = c("#00AFBB", "#FC4E07"), 
                         limits = c(0, 10), 
                         breaks = c(0, 5, 10))
loading_with_labels


##2D PCA Proper
#Reorder factors for plotting
components$Metadata_Pt_type <- ordered(as.factor(components$Metadata_Pt_type), levels = c("CJD", "AD", "CTRL"))
components$Metadata_Cortex_region <- ordered(as.factor(components$Metadata_Cortex_region), levels = c("OCC", "FC")) 
components$Metadata_Layer <- ordered(as.factor(components$Metadata_Layer), levels = c("WM", "dGM", "sGM")) 
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
  range = c(-3.99, 3.99),                
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
fig_pca2d_ptID <- plot_ly(components, type = "scatter", x = ~PC1, y = ~PC3, color = ~as.factor(Patient), 
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
fig_pca2d_disease <- plot_ly(components, type = "scatter", x = ~PC1, y = ~PC3, color = ~as.factor(Metadata_Pt_type), 
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
fig_pca2d_brainregion <- plot_ly(components, type = "scatter", x = ~PC1, y = ~PC3, color = ~as.factor(Metadata_Cortex_region), 
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
fig_pca2d_cortexlayer <- plot_ly(components, type = "scatter", x = ~PC1, y = ~PC3, color = ~as.factor(Metadata_Layer), 
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
  fig_pca2d_ptID,
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
  add_annotations(text = paste('PC 3 (', toString(round(explained_variance_ratio[3], 1)), '%)', sep = ''), 
                  xref = 'paper', yref = 'paper', x = -0.11, y = 0.5, font = list(color = "black", size = 24),
                  textangle = 270, showarrow = F)
fig_pca2d



############Interaction plots of features of interest##############
interaction_plot <- function(feature_name, title) {
  fit_interact <- lmer(Morph_data[[feature_name]] ~ Metadata_Pt_type * Metadata_Cortex_region * Metadata_Layer + (1 | Patient/Metadata_Sample_number), data = Morph_data)
  
  print("!!!!!!!!!!!!!!!!!Model Specifications!!!!!!!!!!!!!!")
  print(fit_interact)
  
  windows(width = 15, height = 15)
  print(performance::check_model(fit_interact))
  
  #Reorder
  interaction_plot <- emmip(fit_interact, Metadata_Pt_type ~ Metadata_Layer | Metadata_Cortex_region, plotit = FALSE, CIs = TRUE)
  interaction_plot$Metadata_Pt_type <- ordered(as.factor(interaction_plot$Metadata_Pt_type), levels = c("CTRL", "AD", "CJD"))
  interaction_plot$Metadata_Layer <- ordered(as.factor(interaction_plot$Metadata_Layer), levels = c("sGM", "dGM", "WM"))
  
  #Metadata_Pt_type Contrasts
  emm_Metadata_Pt_type <- emmeans(fit_interact, pairwise ~ Metadata_Pt_type | Metadata_Layer + Metadata_Cortex_region)
  print("!!!!!!!!!!!DISEASE GROUP CONTRASTS!!!!!!!!!!!!!!!")
  print(emm_Metadata_Pt_type$contrasts)
  
  #Brain Region Contrasts
  emm_Metadata_Cortex_region <- emmeans(fit_interact, pairwise ~ Metadata_Cortex_region | Metadata_Layer + Metadata_Pt_type)
  print("!!!!!!!!!!!BRAIN REGION CONTRASTS!!!!!!!!!!!!!!!")
  print(emm_Metadata_Cortex_region$contrasts)
  
  #Cortex Layer Contrasts
  emm_Metadata_Cortex_region <- emmeans(fit_interact, pairwise ~ Metadata_Layer | Metadata_Cortex_region + Metadata_Pt_type)
  print("!!!!!!!!!!!CORTEX LAYERS CONTRASTS!!!!!!!!!!!!!!!")
  print(emm_Metadata_Cortex_region$contrasts)
  
  
  #Data for single observations
  data <- Morph_data
  data$Metadata_Pt_type <- ordered(as.factor(data$Metadata_Pt_type), levels = c("CTRL", "AD", "CJD"))
  data$Metadata_Layer <- ordered(as.factor(data$Metadata_Layer), levels = c("sGM", "dGM", "WM"))
  
  
  dev.new()
  #Create ggplot
  pd <- position_dodge(.5)
  ggplot(data = interaction_plot, aes(x = Metadata_Layer, y = yvar, 
                                      group = Metadata_Pt_type, linetype = Metadata_Pt_type, shape = Metadata_Pt_type)) +
    facet_nested(~ "Brain Region" + Metadata_Cortex_region)+
    geom_point(aes(x = Metadata_Layer, 
                   y = .data[[feature_name]], fill=as.factor(Patient), color=as.factor(Metadata_Pt_type)), 
               data = data, size = 1.5, alpha=0.8, stroke = 0.5,
               position = position_jitterdodge(dodge.width = 0.5, jitter.width = 1)) +
    geom_point(size = 2.5, fill ="black", color="black", alpha = 0.6, position = pd) +
    geom_line(alpha=0.6, color = "black", position = pd) +
    theme_grey()+
    scale_fill_manual(values = c("Patient 1"="#b15928", "Patient 2" = "#ffff99", "Patient 3" = "#6a3d9a",
                                 "Patient 4" = "#cab2d6", "Patient 5" = "#ff7f00", "Patient 6" = "#fdbf6f",
                                 "Patient 7" = "#e31a1c", "Patient 8" = "#fb9a99", "Patient 9" = "#34a02c",
                                 "Patient 10"="#b2df8a", "Patient 11" = "#1f78b4", "Patient 12" = "#a6cee3"))+
    scale_color_manual(values = c("CJD" = "#C6D300", "AD" = "#D6390E", "CTRL" = "#390F6E"),
                       guide = guide_legend(override.aes = list(shape = c(21, 24, 22),alpha = 1, size = 1.5)))+
    scale_shape_discrete(name = "Estimated Mean",
                         breaks = c("CTRL", "AD", "CJD"),
                         labels = c("CTRL", "AD", "CJD")) +
    scale_shape_manual(values = c("CTRL" = 21,  
                                  "AD" = 24,    
                                  "CJD" = 22))+ 
    ggtitle(title) +
    labs(x = "Cortex Layer", y = "Morphometry Feature Z-score", shape = "Estimated Mean", color = "ROIs", linetype=NULL) +
    guides(fill = "none", 
           shape = guide_legend(order = 2),
           linetype = guide_legend(order = 3))+
    theme(text = element_text(size = 12), plot.title = element_text(size = 16),
          strip.text.x = element_text(face = "bold", size = 10),
          legend.title = element_text(size = 8.5), legend.text = element_text(size = 7), legend.key.size = unit(1, "lines"), 
          legend.spacing = unit(1.05, "cm"),legend.margin = margin(t = -20))}

#Plot - Feature Maximum Radius
interaction_plot("Total_AreaShape_MaximumRadius", "Total - Maximum Radius")

#Plot - Feature Mean Radius
interaction_plot("Total_AreaShape_MeanRadius", "Total - Mean Radius")

#Extra Plot - Feature Texture Correlation(.2)
interaction_plot("Total_Texture_Correlation_IBA1_3_02_256", "Total - Texture Correlation(.2)")


###Supplementary
###Extra - PCA visualized with Covariates
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
      format = "svg",
      filename = "myplot",
      width = NULL,
      height = NULL
    )
  ) %>%
  add_annotations(text = paste('PC 1 (', toString(round(explained_variance_ratio[1], 1)), '%)', sep = ''), 
                  xref = 'paper', yref = 'paper', x = 0.5, y = -0.1, font = list(color = "black", size = 24),
                  showarrow = F) %>%
  add_annotations(text = paste('PC 3 (', toString(round(explained_variance_ratio[3], 1)), '%)', sep = ''), 
                  xref = 'paper', yref = 'paper', x = -0.11, y = 0.5, font = list(color = "black", size = 24),
                  textangle = 270, showarrow = F)%>%
  add_annotations(text = "Sex", xref = "paper", yref = "paper",
                  x = 1.03, xanchor = "left",
                  y = 0.8, yanchor = "bottom",
                  legendtitle = TRUE, font = list(size=14), showarrow = FALSE)
fig_pca2d_extra

###Extra - Evaluating effects of Age, Sex, Formalin Fixation Duration, Postmortem Delay on morphometrics
#Age
results_age <- glmmQvals(lmmSeq(~ Age + (1|Patient) + (1|Metadata_Sample_number),
                                maindata = t(cell_features),
                                metadata = Morph_data[,1:14],
                                progress = TRUE))
stats <- summary(results_age)
kable(stats[order(stats[, 10]), ]) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "400px")

#Sex
results_sex <- glmmQvals(lmmSeq(~ Sex + (1|Patient) + (1|Metadata_Sample_number),
                                maindata = t(cell_features),
                                metadata = Morph_data[,1:14],
                                progress = TRUE))
stats <- summary(results_sex)
kable(stats[order(stats[, 10]), ]) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "400px")

#Formalin Fixation duration
combined_df <- cbind(Morph_data, cell_features)
combined_df_filtered <- combined_df[!is.na(combined_df$Formalinfixed), ]
metadata_clean <- combined_df_filtered[, colnames(Morph_data)]
features_clean <- combined_df_filtered[, colnames(cell_features)]
results_fixation <- glmmQvals(lmmSeq(~ Formalinfixed + (1|Patient) + (1|Metadata_Sample_number),
                                     maindata = t(features_clean),
                                     metadata = metadata_clean[,1:14],
                                     progress = TRUE))
stats <- summary(results_fixation)
kable(stats[order(stats[, 10]), ]) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "400px")


#Postmortem Delay
combined_df <- cbind(Morph_data, cell_features)
combined_df_filtered <- combined_df[!is.na(combined_df$Postdelay), ]
metadata_clean <- combined_df_filtered[, colnames(Morph_data)]
features_clean <- combined_df_filtered[, colnames(cell_features)]
results_postmortem <- glmmQvals(lmmSeq(~ Postdelay + (1|Patient) + (1|Metadata_Sample_number),
                                       maindata = t(features_clean),
                                       metadata = metadata_clean[,1:14],
                                       progress = TRUE))
stats <- summary(results_postmortem)
kable(stats[order(stats[, 10]), ]) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "400px")


##### Extra - Variance Partitioning###############
t_cell_features <- t(cell_features)
metadata <- Morph_data[,1:(feature_column_start-1)]

#reduced Model
reduced <- ~ (1|Metadata_Cortex_region) + (1|Metadata_Layer) + (1|Metadata_Pt_type)+ (1|Metadata_Pt_type:Metadata_Cortex_region) + (1|Metadata_Pt_type:Metadata_Layer) + (1|Sex) + Age

#Partioning
varPart <- fitExtractVarPartModel(t_cell_features, reduced, metadata)

#Rename and reorder
colnames(varPart) <- c("Brain Region", "Cortex Layer", "Disease", 
                       "Disease:Brain Region", "Disease:Cortex Layer", "Sex", "Age", "Residuals")
varPart <- varPart[order((varPart$Residuals), decreasing = F), ]
varCol <- c("#0044CC", "#FFD700", "#FF4500","#8B0000", "orange", "#9467BD", "#87CEFA",
            "lightgrey")

#Variance Partioning Variables Overview
variables_var <- plotVarPart(varPart, col =varCol) + ylab("Variance Explained (%)")
variables_var

#Variance Partioning of all morphometrical features
morph_var <- plotPercentBars(varPart[1:73, ], 
                                col =varCol) + theme(
                                         legend.position = "right", 
                                         legend.box = "horizontal",
                                         legend.title = element_blank(),
                                         legend.key.size = unit(1, "lines"))
morph_var


#full model with Formalix fixation duration and postmortem delay
full <- ~ (1|Metadata_Cortex_region) + (1|Metadata_Layer) + (1|Metadata_Pt_type)+ (1|Metadata_Pt_type:Metadata_Cortex_region) + (1|Metadata_Pt_type:Metadata_Layer) + (1|Sex) + Age + Formalinfixed + Postdelay
varPart <- fitExtractVarPartModel(t_cell_features, full, metadata)
varPart <- varPart[order((varPart$Residuals), decreasing = F), ]
colnames(varPart) <- c("Brain Region", "Cortex Layer", "Disease", 
  "Disease:Brain Region", "Disease:Cortex Layer", "Sex", "Age", "Formalin Fixation Duration", "Postmortem Delay", "Residuals")
colVar <- c("#0044CC", "#FFD700", "#FF4500","#8B0000", "orange", "#9467BD", "#87CEFA", "#F5F5F5", "black", "lightgrey")
plotVarPart(varPart, col =colVar) + ylab("Variance Explained (%)")
plotPercentBars(varPart[1:73, ], 
                col =colVar) + theme(
                  legend.position = "right", 
                  legend.box = "horizontal",
                  legend.title = element_blank(),
                  legend.key.size = unit(1, "lines"))


sessionInfo()
