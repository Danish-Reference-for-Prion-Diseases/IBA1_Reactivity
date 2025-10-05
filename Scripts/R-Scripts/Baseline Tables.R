#Libraries
library(here)
library(tidyverse)
library("kableExtra")
library("table1")
library("ggalluvial")

#Data
data_dir <- here("Data", "GeoMx (protein)")
Full_Metadata_Geomx <- read.csv(here(data_dir, "Full_GeoMx_Metadata.csv"))

###################Sankey plot
#Prepare Data
Full_Metadata_Geomx <- Full_Metadata_Geomx %>% 
  mutate(CortexLayer = recode(CortexLayer,
                              "TM"="sGM",
                              "BM" = "dGM"))


grouped <- Full_Metadata_Geomx %>% group_by(Patient, Disease, BrainRegion, Slide, CortexLayer, Celltype) %>% count()

grouped$Celltype <- ordered(as.factor(grouped$Celltype), levels = c("Missing", "IBA1")) 
grouped$Disease <- ordered(as.factor(grouped$Disease), levels = c("CTRL", "AD", "CJD"))
grouped$BrainRegion <- ordered(as.factor(grouped$BrainRegion), levels = c("FC", "OCC")) 
grouped$Slide <- ordered(as.factor(grouped$Slide), levels = c("Slide 1", "Slide 2", "Slide 3", "Slide 4", "Slide 5", "Slide 6", "Slide 7", "Slide 8", "Slide 9", "Slide 10", "Slide 11", "Slide 12")) #Change this
grouped$CortexLayer <- ordered(as.factor(grouped$CortexLayer), levels = c("sGM", "dGM", "WM")) 
grouped$Patient <- ordered(as.factor(grouped$Patient), levels = c("Patient 1", "Patient 2", "Patient 3", "Patient 4", "Patient 5", "Patient 6",
                                                                  "Patient 7", "Patient 8","Patient 9","Patient 10", "Patient 11", "Patient 12")) 

levels(grouped$Slide) <- c("1","2","3","4","5","6","7","8","9","10","11","12")
levels(grouped$Patient) <- c("1","2","3","4","5","6","7","8","9","10","11","12")

# Create the plot
ggplot(as.data.frame(grouped),
       aes(y = n, axis1 = Disease, axis2 = Patient, axis3 = BrainRegion, axis4 = Slide, axis5 = CortexLayer, axis6 = Celltype)) +
  geom_alluvium(aes(fill = Patient), width = 1/14) +
  geom_curve(aes(x = 6.05, y = 180, xend = 6.65, yend = 180),
             inherit.aes = FALSE,
             curvature = -0.3,
             arrow = arrow(length = unit(0.2,"cm")))+
  geom_stratum(width = 1/12, color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Patient Group", "Patient ID", "Brain Region", "Slide Number", "Cortex Layer", "Segments (GeoMx)"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Paired", direction = -1) +
  labs(y = "ROIs (n = 186)", x = "") +
  ggtitle("") +
  theme_minimal() +  
  theme(
    legend.position = "none",   
    text = element_text(size = 16),   
    plot.margin = margin(5, 120, 5, 5),  
    axis.text = element_text(size = 12)  
  ) +
  annotate("text",x=6.4,y=93,label="IBA1 Segments (n = 176)", size = 6, fontface = "plain", angle=90)+
  coord_cartesian(ylim=c(0,190),xlim = c(1,6),clip="off")+
  geom_rect(aes(xmin = 6.7, xmax = 7.49,
                ymin = 124, ymax = 186),
            inherit.aes = FALSE,
            fill = "white", colour = "black") +
  annotate(
    "text",
    x = 6.76, y = 155,
    label = paste(
      "Missing segments",
      "AD4-FC-sGM: 4",
      "AD3-FC-dGM: 1",
      "AD3-FC-WM: 1",
      "AD4-OCC-sGM: 1",
      "AD7-OCC-dGM: 1",
      "CJD10-OCC-dGM: 1",
      "CJD11-OCC-WM: 1",
      sep = "\n"
    ),
    hjust = 0,
    lineheight = 1,
    size = 2.6)
  

sessionInfo()
