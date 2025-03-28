#Libraries
library(here)
library(tidyverse)
library("kableExtra")
library("table1")
library("ggalluvial")
library(readxl)
library(gt)
library(gtsummary)
library(grid)


#Data
data_dir <- here("Data", "GeoMx (protein)")
Full_Metadata_Geomx <- read.csv(here(data_dir, "Full_GeoMx_Metadata.csv"))
protein_list <- read_excel(here("Data", "GeoMx (protein)", "Protein_panels_list.xlsx"))

#Baseline tables of patientes and proteins
collapsed <- cbind(Full_Metadata_Geomx[,1:3], 'Disease Group' = Full_Metadata_Geomx[,31])
collapsed <- collapsed %>% distinct(Patient, .keep_all = TRUE)
collapsed <- collapsed[,c(3,4,1,2)]
collapsed$Patient <- factor(collapsed$Patient, ordered=TRUE, levels=c("Patient 1", "Patient 2", "Patient 3", "Patient 4", "Patient 5", "Patient 6",
                                                                      "Patient 7", "Patient 8", "Patient 9", "Patient 10", "Patient 11", "Patient 12"))
collapsed <- collapsed %>% arrange(Patient)
colnames(collapsed)[4] <- "Age (years at death)" 

collapsed$'Disease Group' <- factor(collapsed$'Disease Group', ordered=TRUE, levels=c("CTRL", "AD", "CJD"))


#####Patient Data
gt_tbl <- collapsed %>%
  gt() %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_body(columns = where(~ is.character(.) | is.factor(.)))  
  ) %>%
  tab_style(
    style = cell_text(align = "left"), 
    locations = cells_column_labels(columns = where(~ is.character(.) | is.factor(.)))  
  ) %>%
  
  tab_style(
    style = cell_text(align = "right"),  
    locations = cells_body(columns = where(~ is.numeric(.))) 
  ) %>%
  tab_style(
    style = cell_text(align = "right"),  
    locations = cells_column_labels(columns = where(~ is.numeric(.)))  
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),  
    locations = cells_column_labels(columns = everything())  
  ) %>%
 
  tab_options(
    table.width = pct(100), 
    table.layout = "auto"  
  ) %>%
  
  tab_footnote(
    footnote = "CTRL = Controls without neuropathological findings;",
    locations = cells_body(columns = "Disease Group", rows = 1:2)
  ) %>%
  tab_footnote(
    footnote = "AD = Alzheimer’s Disease (All Braak Stage VI)",
    locations = cells_body(columns = "Disease Group", rows = 3:7)
  ) %>%
  tab_footnote(
    footnote = "CJD = Creutzfeldt Jacob Disease (All sporadic, subtype MM1)",
    locations = cells_body(columns = "Disease Group", rows = 8:12)
  )

gt_tbl




#Protein list
# Replace NA values with empty strings
protein_list[is.na(protein_list)] <- ""

protein_list %>%
  gt() %>%
  tab_stubhead(label = "Protein Assay Group") %>%
  tab_style(
    style = cell_text(weight = "bold", size = px(14)),
    locations = cells_stubhead()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_options(
    table.width = pct(100),
    data_row.padding = px(8) 
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = 1) 
  )


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
  geom_stratum(width = 1/12, color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Patient Group", "Patient ID", "Brain Region", "Slide Number", "Cortex Layer", "Segments (GeoMx)"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Paired", direction = -1) +
  labs(y = "ROIs (n = 186)", x = "") +
  ggtitle("Subsampling Overview") +
  theme_minimal() +  
  
  theme(
    legend.position = "none",   
    text = element_text(size = 16),   
    plot.margin = margin(5, 55, 5, 5),  
    axis.text = element_text(size = 12)  
  ) +
  annotate("text",x=6.4,y=93,label="IBA1 Segments (n = 176)", size = 6, fontface = "plain", angle=90)+
  coord_cartesian(ylim=c(0,190),xlim = c(1,6),clip="off")



sessionInfo()
