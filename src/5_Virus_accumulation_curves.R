# manuscript: Johnson et al. 2026. Large-scale surveillance at 
#             wildlife-human interfaces reveals virus spillover  
#             risk in wildlife markets and trade supply chains. 
#             Nature Microbiology.
# analysis: Infection prevalence among bats by preferred roost habitat

###############################
## Table of Contents       ####
###############################

## 1.  LOAD AND FILTER DATA
## 2.  SCATTERPLOTS FOR # ANIMALS POSITIVE AND # VIRUSES DETECTED ~ # ANIMALS SAMPLED [FIG 3]
## 3.  FAMILY-LEVEL VIRUS ACCUMULATION CURVES [EXT DATA FIG 3]
## 4.  SPECIES-LEVEL VIRUS ACCUMULATION CURVES (ALL VIRUSES AND BY VIRUS FAMILY) [EXT DATA FIG 4-7]
## 5.  EXTRAPOLATED FAMILY-LEVEL CURVES USING iNEXT (USING HILL NUMBERS) [EXT DATA FIG 8]
## 6.  PARAMETRIC & NONPARAMETRIC ESTIMATORS OF ASYMPTOTES [SUPPL TABLE 5]
## 7.  MICHELIS-MENTEN EXTRAPOLATED VIRUS ACCUMULATION CURVES [FIG 4]

###############################
## 1. LOAD AND FILTER DATA ####
###############################

library(tidyverse)
library(vegan); library(iNEXT); 
library(data.table)
library(mgcv); library(scales)
library(minpack.lm); library(boot)
library(patchwork)
library(grid)
graphics.off()

output_dir <- "plots/Final Figs/"

text_to_binary <- function(x) {ifelse(is.na(x), 0, 1)}

pred <- read.csv("data3_virus_accumulation_curves.csv")

## CREATE DATASETS FOR COVs and PMXs
pred_original <- pred
pred_covs <- pred[pred$VirusFamilyTested=="Coronaviridae",]
pred_pmxs <- pred[pred$VirusFamilyTested=="Paramyxoviridae",]

##############################################################
## 2. SAMPLING EFFORT AND VIRUS DETECTION BY HOST SPECIES ####
## Figures 3a and 3b, 2014-2019 data                        ##
##############################################################

pred_strat <- pred_original %>% filter(Timeframe == "2014-2019")

number_positive <- pred_strat %>% group_by(Taxa,ScientificName,individual_id) %>% 
    mutate(positive_any = ifelse(any(TestResult == "Positive"), 1, 0)) %>% 
    summarize(positive_any = ifelse(any(positive_any == 1), 1, 0))

unique_detections <- number_positive %>% group_by(Taxa,ScientificName) %>% 
    summarise(Sample_Size = n_distinct(individual_id), Unique_detections = sum(positive_any))
unique_detections <- unique_detections[unique_detections$ScientificName!="",]
unique_detections$face <- "italic"
unique_detections <- unique_detections[order(unique_detections$Sample_Size), ]
unique_detections$alpha <- seq(0.6,1,length=length(unique_detections$ScientificName))

gam_model <- gam(Unique_detections~ s(Sample_Size), data = unique_detections)
unique_detections$newdata <- predict(gam_model)

#### figure labels -----------

labels_fig3a = data.frame(Taxa = c('bats', 'bats', 'bats', 'bats', 'bats', 
                                   'rodents & shrews', 'rodents & shrews', 
                                   'bats', 'rodents & shrews', 'bats', 
                                   'rodents & shrews', 'bats', 'bats', 
                                   'bats', 'rodents & shrews', 'bats'),
                          ScientificName = c('Myotis laniger', 'Myotis ricketti', 'Miniopterus schreibersii', 
                                             'Rhinopoma hardwickii', 'Mops condylurus', 'Rattus norvegicus', 
                                             'Rattus exulans', 'Hipposideros ruber', 'Rattus argentiventer', 
                                             'Pteropus giganteus', 'Atherurus africanus', 'Pteropus lylei', 
                                             'Scotophilus kuhlii', 'Rousettus aegyptiacus', 'Rattus rattus', 
                                             'Eidolon helvum'))

labels_fig3a_rt = c("Rattus exulans", "Rousettus aegyptiacus", "Pteropus giganteus",
                    "Rattus rattus", "Atherurus africanus")

labels_fig3b = data.frame(Taxa = c('bats', 'bats', 'bats', 'bats', 'bats',
                                   'rodents & shrews', 'rodents & shrews',
                                   'bats', 'bats', 'bats', 'bats',
                                   'rodents & shrews', 'bats', 'bats'),
                          ScientificName = c('Hipposideros cervinus', 'Rhinolophus affinis', 
                                             'Miniopterus schreibersii', 'Rhinolophus creaghi',
                                             'Hipposideros ruber', 'Rattus argentiventer', 
                                             'Atherurus africanus', 'Pteropus alecto', 
                                             'Pteropus lylei', 'Scotophilus kuhlii', 
                                             'Rousettus aegyptiacus', 'Rattus rattus', 
                                             'Eidolon helvum', 'Myotis laniger'))

labels_fig3b_rt = c("Rattus argentiventer","Atherurus africanus","Scotophilus kuhlii",
                    "Rattus rattus","Eidolon helvum")

#### FIG 3a -- individuals positive by host species ----

fig3a = ggplot(unique_detections, aes(x = Sample_Size, y = Unique_detections, 
                                      family='sans-serif', color=Taxa)) +
    geom_point(aes(shape=factor(Taxa)),alpha=1,size=1) + 
    scale_color_manual(values = c("#8E24AA","#219985")) +  
    labs(x = "Number of individuals sampled within a species", 
         y = "Number of individuals detected positive", title = "") + 
    theme_bw() + 
    scale_x_log10(breaks = c(1, 10, 100, 1000, 10000), expand = expansion(mult = c(0, 0.1))) + 
    theme(legend.position = 'none', 
          panel.border = element_blank(), 
          panel.grid = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          axis.title = element_text(size = 7, color = "black"), 
          axis.text = element_text(size = 6, color = "black")) +  
    geom_smooth(aes(color=NULL),method = "gam", formula = y ~ s(x),
                col='gray28', linewidth = 0.5) +
    
    
    
    # labels left of smooth line
    ggrepel::geom_text_repel(
        data = unique_detections %>% 
            mutate(ScientificName = ifelse(ScientificName %in% labels_fig3a$ScientificName &
                                               !ScientificName %in% labels_fig3a_rt,
                                           ScientificName, "")),
        aes(x = Sample_Size, y = Unique_detections, label = ScientificName),
        size = 1.8, 
        nudge_x = -0.15, nudge_y = 0.75, 
        fontface = "italic",
        segment.color = "black", 
        direction = "both",
        seed = 100,
        point.padding = 0.3, box.padding = 0.1,
        force = 0.2, force_pull = 10, 
        max.overlaps = Inf, 
        min.segment.length = 20,
        lineheight = 0.9) +
    
    # labels right of smooth line
    ggrepel::geom_text_repel(
        data = unique_detections %>% 
            mutate(ScientificName = ifelse(ScientificName %in% labels_fig3a_rt,
                                           ScientificName, ""),
                   ScientificName = gsub("Rousettus aegyptiacus","Rousettus\naegyptiacus",ScientificName)),
        aes(x = Sample_Size, y = Unique_detections, label = ScientificName),
        size = 1.8, 
        nudge_x = 0.1, nudge_y = -0.05,
        fontface = "italic",
        segment.color = "black", 
        direction = "both",
        seed = 100,
        point.padding = 0.3, box.padding = 0.1,
        force = 0.4, force_pull = 10, 
        max.overlaps = Inf,
        min.segment.length = 10,
        lineheight = 0.9) + 
    
    # Fig label "a"
    annotate(
        "text",
        x = 1,  
        y = 270,  
        label = "a",
        fontface = "bold",
        hjust = 0,    
        vjust = 1,     
        size = 2.46 
    )

gam_model <- gam(Unique_detections~ s(Sample_Size), data = unique_detections)
summary(gam_model)

#### FIG 3b -- unique viruses detected by host species ----

unique_counts <- pred_strat %>% group_by(Taxa,ScientificName) %>% 
    summarise(Sample_Size = n_distinct(individual_id), Unique_Virus_Count = n_distinct(VirusName))
unique_counts <- unique_counts[unique_counts$ScientificName!="",]
unique_counts$face <- "italic"
unique_counts <- unique_counts[order(unique_counts$Sample_Size), ]
unique_counts$alpha <- seq(0.6,1,length=length(unique_counts$ScientificName))

fig3b = ggplot(unique_counts, aes(x = Sample_Size, y = Unique_Virus_Count,family='sans-serif',color=Taxa)) +
    geom_point(aes(shape=factor(Taxa)),alpha=1,size=1) +
    scale_color_manual(values = c("#8E24AA","#219985")) +
    labs(x = "Number of individuals sampled within a species", 
         y = "Number of unique viruses detected", title = "") +   
    theme_bw() + 
    scale_x_log10(breaks = c(1, 10, 100, 1000, 10000), expand = expansion(mult = c(0, 0.1))) +   
    theme(legend.position = 'none', 
          panel.border = element_blank(), 
          panel.grid = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          axis.title = element_text(size = 7, color = "black"), 
          axis.text = element_text(size = 6, color = "black")) + 
    geom_smooth(aes(color=NULL), method = "gam", formula = y ~ s(x),
                col='gray28', linewidth = 0.5) +
    
    # labels left of smooth line
    ggrepel::geom_text_repel(
        data = unique_counts %>% 
            mutate(ScientificName = ifelse(ScientificName %in% labels_fig3b$ScientificName &
                                               !ScientificName %in% labels_fig3b_rt,
                                           ScientificName, ""),
                   ScientificName = gsub("Miniopterus schreibersii","Miniopterus\nschreibersii",ScientificName),
                   ScientificName = gsub("Rhinolophus sinicus","Rhinolophus\nsinicus",ScientificName)),
        aes(label = ScientificName), 
        size = 1.8,
        nudge_x = -0.15, nudge_y = 0.75,
        fontface = "italic",
        segment.color = "black", 
        direction = "both",
        seed = 100,
        point.padding = 0.3, box.padding = 0.1, 
        force = 0.2, force_pull = 10, 
        max.overlaps = Inf,
        min.segment.length = 0,
        lineheight = 0.9) +
    
    # labels right of smooth line
    ggrepel::geom_text_repel(
        data = unique_counts %>% 
            mutate(ScientificName = ifelse(ScientificName %in% labels_fig3b_rt,
                                           ScientificName, "")),
        aes(label = ScientificName),
        size = 1.8,
        nudge_x = 0.15, nudge_y = -0.1,
        fontface = "italic",
        segment.color = "black", 
        direction = "both",
        seed = 100,
        point.padding = 0.3, box.padding = 0.1,
        force = 0.2, force_pull = 10, 
        max.overlaps = Inf,
        min.segment.length = 0,
        lineheight = 0.9) +
    
    # Fig label "b"
    annotate(
        "text",
        x = 1,  
        y = 27,  
        label = "b",
        fontface = "bold",
        hjust = 0,      
        vjust = 1,      
        size = 2.46
    )

gam_model <- gam(Unique_Virus_Count~ s(Sample_Size), data = unique_counts)
summary(gam_model)

#### FIG 3, combined ----

combined_plot_fig3 = fig3a + fig3b +
    plot_layout(ncol = 2)

ggsave(paste0(output_dir, "Fig3ab_2026_FINAL.pdf"),
       combined_plot_fig3,
       dpi = 600, 
       width = 180,   
       height = 90,  
       units = "mm",
       device = cairo_pdf)

##############################
## 3. FAMILY-LEVEL CURVES ####
## Extended Data Figure 3   ##
##############################

pred_strat <- pred_original  
unique_counts <- pred_strat %>% group_by(Family) %>%
    summarise(Sample_Size = n_distinct(individual_id), 
              Unique_Virus_Count = n_distinct(VirusName))
unique_counts <- unique_counts[-1,]
list_excl <- unique_counts$Family[unique_counts$Sample_Size<=10 | 
                                      unique_counts$Unique_Virus_Count<=3]
pred_fam <- pred_strat[!(pred_strat$Family %in% list_excl), ]
fam_list <- unique(pred_fam$Family); fam_list <- sort(fam_list)[-1]

#### EXT DATA FIG 3 ----

pdf(paste0(output_dir,"FigExt3_2026_FINAL.pdf"), 
    width = 180/25.4, height = 93.75/25.4) 
par(mfrow = c(4,6)) 
par(mar = c(1.2, 1.2, 1.2, 0.4)) 
par(mgp = c(2, 0.3, 0))
par(tcl = -0.3)
par(font.main = 1) 
par(oma = c(1.5, 1, 0, 0))

for (i in fam_list) {
    pred_spp <- pred_fam[pred_fam$Family==i,]
    pred_spp <- pred_spp[,c("SamplingDate","individual_id","specimen_id","TestResult","VirusName")]
    pred_spp_byanimal <- pred_spp %>% filter(TestResult == "Positive" & VirusName!="") %>% 
        distinct(individual_id,specimen_id,VirusName) %>% group_by(individual_id) %>% 
        summarise(NumberDistinctViruses = length(unique(VirusName)), 
                  ListDistinctViruses = knitr::combine_words(unique(VirusName),sep=";",and=""))
    pred_spp2 <- merge(pred_spp, pred_spp_byanimal, by="individual_id", all.x=TRUE)
    pred_spp2$NumberDistinctViruses[is.na(pred_spp2$NumberDistinctViruses)] <- 0
    pred_spp3 <- pred_spp2 %>% 
        distinct(individual_id,SamplingDate,NumberDistinctViruses,ListDistinctViruses)
    
    df <- separate_rows(pred_spp3, ListDistinctViruses, sep = ";\\s*")
    df2 <- spread(df, ListDistinctViruses, ListDistinctViruses)
    df2[, -(1:3)] <- lapply(df2[, -(1:3)], text_to_binary)
    df2 <- df2 %>% arrange(SamplingDate)
    df2 <- df2[,-c(1:3)]
    sp1 <- specaccum(df2, method = "exact")
    
    plot(sp1, col="#D1D2E9", 
         xlab = "N", ylab = "VR", 
         main = i, cex.main = 0.75, cex.axis = 0.625, cex.lab = 0.625) 
    lines(sp1$richness,lwd = 2, col="#1A1A57") 
    box()
} 

grid.text(
    "Sample size (n)",
    x = unit(0.5, "npc"),
    y = unit(0.02, "npc"),
    just = c("center", "bottom"),
    gp = gpar(
        fontsize = 7
    )
)

grid.text(
    "Number of unique viruses",
    x = unit(0.01, "npc"),
    y = unit(0.5, "npc"),
    just = c("center", "center"),
    rot = 90,
    gp = gpar(
        fontsize = 7
    )
)

dev.off()



############################################
## 4. VIRUS DISCOVERY CURVES BY SPECIES ####
## Extended Data Figures 4-7              ##
############################################

#### EXT DATA FIG 4, 6, 7 -- INDIVIDUAL SPECIES, ALL VIRUSES ----

vac_type <- c("all_vac", "cov_vac", "pmx_vac")

for(vac in vac_type){
    
    if(vac == "all_vac"){
        color_1 <- "#D1D2E9"
        color_2 <- "#1A1A57"
        pred_strat <- pred_original
    }
    
    if(vac == "cov_vac"){
        color_1 <- "#A6D6D6"
        color_2 <- "#2D708E"
        pred_strat <- pred_covs
    }
    
    if(vac == "pmx_vac"){
        color_1 <- adjustcolor("#B8DE29", alpha.f = 0.6)
        color_2 <- "#4D780E" 
        pred_strat <- pred_pmxs
    }
    
    unique_counts <- pred_strat %>% 
        group_by(Taxa,ScientificName) %>% 
        summarise(Sample_Size = n_distinct(individual_id), 
                  Unique_Virus_Count = n_distinct(VirusName))
    unique_counts <- unique_counts[unique_counts$ScientificName!="",]
    
    if(vac == "pmx_vac"){
        list_excl <- unique_counts$ScientificName[unique_counts$Sample_Size<=10 | 
                                                      unique_counts$Unique_Virus_Count<3]
    } else{
        list_excl <- unique_counts$ScientificName[unique_counts$Sample_Size<=10 | 
                                                      unique_counts$Unique_Virus_Count<=3] 
    }
    
    pred_strat <- pred_strat[!(pred_strat$ScientificName %in% list_excl), ]
    spp_list <- unique(pred_strat$ScientificName); spp_list <- sort(spp_list)[-1]
    
    
    
    pdf(paste0(output_dir,vac,"_plot.pdf"), 
        width = 180/25.4, height = 240/25.4) 
    par(mfrow = c(12,6)) 
    par(mar = c(1.2, 1.2, 1.2, 0.4)) 
    par(mgp = c(2, 0.3, 0))
    par(tcl = -0.3)
    par(font.main = 3) 
    par(oma = c(1.5, 1, 0, 0))
    
    for (i in spp_list) {
        pred_spp <- pred_strat[pred_strat$ScientificName==i,]
        pred_spp <- pred_spp[,c("SamplingDate","individual_id", 
                                "specimen_id","TestResult","VirusName")]
        pred_spp_byanimal <- pred_spp %>% 
            filter(TestResult == "Positive" & VirusName!="") %>% 
            distinct(individual_id,specimen_id,VirusName) %>% 
            group_by(individual_id) %>% 
            summarise(NumberDistinctViruses = length(unique(VirusName)), 
                      ListDistinctViruses = knitr::combine_words(unique(VirusName),sep=";",and=""))
        
        pred_spp2 <- merge(pred_spp, pred_spp_byanimal, by="individual_id", all.x=TRUE)
        pred_spp2$NumberDistinctViruses[is.na(pred_spp2$NumberDistinctViruses)] <- 0
        pred_spp3 <- pred_spp2 %>% 
            distinct(individual_id, SamplingDate, NumberDistinctViruses, ListDistinctViruses)
        df <- separate_rows(pred_spp3, ListDistinctViruses, sep = ";\\s*")
        df2 <- spread(df, ListDistinctViruses, ListDistinctViruses)
        df2[, -(1:3)] <- lapply(df2[, -(1:3)], text_to_binary)
        df2 <- df2 %>% arrange(SamplingDate)
        df2 <- df2[,-c(1:3)]
        sp1 <- specaccum(df2, method = "exact")
        
        y_max <- ifelse(max(sp1$richness)<=4,
                        max(sp1$richness),
                        ifelse(max(sp1$richness)>=20,
                               ceiling(max(sp1$richness)/5)*5,
                               ceiling(max(sp1$richness)/4)*4))
        
        y_div <- ifelse(max(sp1$richness)<4,
                        max(sp1$richness), 
                        ifelse(max(sp1$richness)>=20, 5, 4))
        
        plot(sp1,col=color_1, xlab = "N", 
             ci.type = "polygon", ci.lty = 0, 
             yaxp =  c(0, y_max, y_div), 
             ylab = "VR", main = i, cex.main=0.75, cex.axis=0.625, cex.lab=0.625) 
        lines(sp1$richness,lwd = 2, col= color_2) 
        
        box()
    } 
    
    grid.text(
        "Sample size (n)",
        x = unit(0.5, "npc"),
        y = unit(0.01, "npc"),
        just = c("center", "bottom"),
        gp = gpar(
            fontsize = 7
        )
    )
    
    grid.text(
        "Number of unique viruses",
        x = unit(0.01, "npc"),
        y = unit(0.5, "npc"),
        just = c("center", "center"),
        rot = 90,
        gp = gpar(
            fontsize = 7
        )
    )
    
    dev.off()
}


#### EXT DATA FIG 5 -- OVERLAPPING COV/PMX VAC CURVES ----

## Coronavirus stratified dataset
pred_strat <- pred_covs 
unique_counts <- pred_strat %>% 
    group_by(ScientificName) %>% 
    summarise(Sample_Size = n_distinct(individual_id), 
              Unique_Virus_Count = n_distinct(VirusName)) 
unique_counts <- unique_counts[-1,]
list_excl <- unique_counts$ScientificName[unique_counts$Sample_Size<=10 | unique_counts$Unique_Virus_Count<=3]
pred_strat_covs <- pred_strat[!(pred_strat$ScientificName %in% list_excl), ]
spp_list_covs <- unique(pred_strat_covs$ScientificName); spp_list_covs <- sort(spp_list_covs)[-1]

## Paramyxovirus stratified dataset
pred_strat <- pred_pmxs 
unique_counts <- pred_strat %>% 
    group_by(ScientificName) %>% 
    summarise(Sample_Size = n_distinct(individual_id), 
              Unique_Virus_Count = n_distinct(VirusName))
unique_counts <- unique_counts[-1,]
list_excl <- unique_counts$ScientificName[unique_counts$Sample_Size<=10 | unique_counts$Unique_Virus_Count<3]
pred_strat_pmxs <- pred_strat[!(pred_strat$ScientificName %in% list_excl), ]
spp_list_pmxs <- unique(pred_strat_pmxs$ScientificName); spp_list_pmxs <- sort(spp_list_pmxs)[-1]

pred_virus_list <- list(pred_strat_covs,pred_strat_pmxs)
spp_list_viruses <- intersect(spp_list_covs, spp_list_pmxs)

# Loop through all species and both datasets

pdf(paste0(output_dir,"FigExt5_2026_FINAL.pdf"), 
    width = 180/25.4, height = 93.75/25.4) # width = 180/25.4, height = 93.75/25.4) 
par(mfrow = c(4,6)) 
par(mar = c(1.2, 1.2, 1.2, 0.4)) 
par(mgp = c(2, 0.3, 0))
par(tcl = -0.3)
par(font.main = 3)
par(oma = c(1.5, 1, 0, 0))

for (i in spp_list_viruses) {
    pred_strat <- pred_strat_covs
    pred_spp <- pred_strat[pred_strat$ScientificName==i,]
    pred_spp <- pred_spp[,c("SamplingDate","individual_id","specimen_id","TestResult","VirusName")]
    pred_spp_byanimal <- pred_spp %>% 
        filter(TestResult == "Positive" & VirusName!="") %>% 
        distinct(individual_id,specimen_id,VirusName) %>% 
        group_by(individual_id) %>% 
        summarise(NumberDistinctViruses = length(unique(VirusName)), 
                  ListDistinctViruses = knitr::combine_words(unique(VirusName),sep=";",and=""))
    pred_spp2 <- merge(pred_spp, pred_spp_byanimal, by="individual_id", all.x=TRUE)
    pred_spp2$NumberDistinctViruses[is.na(pred_spp2$NumberDistinctViruses)] <- 0
    pred_spp3 <- pred_spp2 %>% distinct(individual_id,SamplingDate,NumberDistinctViruses,ListDistinctViruses)
    df <- separate_rows(pred_spp3, ListDistinctViruses, sep = ";\\s*")
    df2 <- spread(df, ListDistinctViruses, ListDistinctViruses)
    df2[, -(1:3)] <- lapply(df2[, -(1:3)], text_to_binary)
    df2 <- df2 %>% arrange(SamplingDate)
    df2 <- df2[,-c(1:3)]
    sp1 <- specaccum(df2, method = "exact")
    
    pred_strat <- pred_strat_pmxs
    pred_spp <- pred_strat[pred_strat$ScientificName==i,]
    pred_spp <- pred_spp[,c("SamplingDate","individual_id","specimen_id","TestResult","VirusName")]
    pred_spp_byanimal <- pred_spp %>% 
        filter(TestResult == "Positive" & VirusName!="") %>% 
        distinct(individual_id,specimen_id,VirusName) %>% 
        group_by(individual_id) %>% 
        summarise(NumberDistinctViruses = length(unique(VirusName)), 
                  ListDistinctViruses = knitr::combine_words(unique(VirusName),sep=";",and=""))
    pred_spp2 <- merge(pred_spp, pred_spp_byanimal, by="individual_id", all.x=TRUE)
    pred_spp2$NumberDistinctViruses[is.na(pred_spp2$NumberDistinctViruses)] <- 0
    pred_spp3 <- pred_spp2 %>% 
        distinct(individual_id,SamplingDate,NumberDistinctViruses,ListDistinctViruses)
    df <- separate_rows(pred_spp3, ListDistinctViruses, sep = ";\\s*")
    df2 <- spread(df, ListDistinctViruses, ListDistinctViruses)
    df2[, -(1:3)] <- lapply(df2[, -(1:3)], text_to_binary)
    df2 <- df2 %>% arrange(SamplingDate)
    df2 <- df2[,-c(1:3)]
    sp2 <- specaccum(df2, method = "exact")
    
    max_samples <- max(max(sp1$sites), max(sp2$sites))
    max_richness <- max(max(sp1$richness), max(sp2$richness))
    
    plot(sp1,col = "#A6D6D6",
         ci.type = "polygon", ci.lty = 0,
         xlim = c(0, max_samples), ylim = c(0, max_richness+1),
         xlab = "N",  ylab = "VR", 
         main = i, cex.main = 0.75, cex.axis = 0.625, cex.lab = 0.625) ### default is 8pt = 1
    lines(sp1$richness,lwd = 2, col = '#2D708E') 
    lines(sp2,col = adjustcolor("#B8DE29", alpha.f = 0.6), ci.type = "polygon", ci.lty = 0)
    lines(sp2$richness,col="#4D780E",lwd = 2)
    box()
}

legend(
    "bottomright",
    inset = c(0, -0.6),
    legend = c("Coronaviruses", "Paramyxoviruses"),
    col = c("#2D708E", "#4D780E"),
    fill = c("#A6D6D6", adjustcolor("#B8DE29", alpha.f = 0.6)),
    lwd = 1.25,
    border = NA,
    bty = "n",
    horiz = TRUE,
    xpd = NA,
    cex = 0.625
)

grid.text(
    "Sample size (n)",
    x = unit(0.5, "npc"),
    y = unit(0.02, "npc"),
    just = c("center", "bottom"),
    gp = gpar(
        fontsize = 7
    )
)

grid.text(
    "Number of unique viruses",
    x = unit(0.01, "npc"),
    y = unit(0.5, "npc"),
    just = c("center", "center"),
    rot = 90,
    gp = gpar(
        fontsize = 7
    )
)

dev.off()


####################################
## 5. FAMILY CURVES USING iNEXT ####
## Extended Data Figure 8         ##
####################################

pred_strat <- pred_original 
unique_counts <- pred_strat %>% group_by(Family) %>%
    summarise(Sample_Size = n_distinct(individual_id), Unique_Virus_Count = n_distinct(VirusName))
unique_counts <- unique_counts[-1,]
list_excl <- unique_counts$Family[unique_counts$Sample_Size<=10 | unique_counts$Unique_Virus_Count<=3]
pred_fam <- pred_strat[!(pred_strat$Family %in% list_excl), ]
fam_list <- unique(pred_fam$Family)
fam_list <- sort(fam_list)[-1]
unique_counts <- pred_fam %>% 
    group_by(Family) %>%
    summarise(Sample_Size = n_distinct(individual_id), 
              Unique_Virus_Count = n_distinct(VirusName))
unique_counts <- unique_counts[-1,]

unique_counts2 <- pred_fam %>% group_by(Family) %>%
    summarise(No_Spp = n_distinct(ScientificName), 
              Sample_Size = n_distinct(individual_id), 
              Unique_Virus_Count = n_distinct(VirusName))
unique_counts2 <- unique_counts2[-1,]
unique_counts2$FamilySppRichness <- c(47,80,11,5,20,100,870,186,106,4,280,300,36,318) ## Manually enter based on dataset
unique_counts2$Taxa <- c("Bats","Bats","Rodents","Bats","Bats","Bats","Rodents", "Bats","Bats","Bats","Rodents","Rodents","Rodents","Bats") ## Manually enter based on dataset
temp <- unique_counts2


fam_list_lowN <- sort(unique(unique_counts$Family[unique_counts$Sample_Size<1500]))
fam_list_largeN <- sort(unique(unique_counts$Family[unique_counts$Sample_Size>=1500]))
fam_datasets <- list(fam_list_lowN, fam_list_largeN)
out = list()

for(j in 1:2){
    fam_data = fam_datasets[[j]]
    matrix_list <- list()
    for (i in fam_data) {
        pred_spp <- pred_fam[pred_fam$Family==i,]
        pred_spp <- pred_spp[,c("SamplingDate","individual_id","specimen_id","TestResult","VirusName")]
        pred_spp_byanimal <- pred_spp %>% 
            filter(TestResult == "Positive" & VirusName!="") %>% 
            distinct(individual_id,specimen_id,VirusName) %>% 
            group_by(individual_id) %>% 
            summarise(NumberDistinctViruses = length(unique(VirusName)), 
                      ListDistinctViruses = knitr::combine_words(unique(VirusName),sep=";",and=""))
        pred_spp2 <- merge(pred_spp, pred_spp_byanimal, by="individual_id", all.x=TRUE)
        pred_spp2$NumberDistinctViruses[is.na(pred_spp2$NumberDistinctViruses)] <- 0
        pred_spp3 <- pred_spp2 %>% 
            distinct(individual_id,SamplingDate,NumberDistinctViruses,ListDistinctViruses)
        df <- separate_rows(pred_spp3, ListDistinctViruses, sep = ";\\s*")
        df2 <- spread(df, ListDistinctViruses, ListDistinctViruses)
        df2[, -(1:3)] <- lapply(df2[, -(1:3)], text_to_binary)
        df2 <- df2 %>% arrange(SamplingDate)
        df3 <- df2[,-c(1:3)]
        df4 <- as.matrix(df3)
        rownames(df4) <- df2$individual_id
        df5 <- t(df4)
        matrix_list[[i]] <- df5
    } 
    out[[j]] <- iNEXT(matrix_list, q=0, datatype="incidence_raw")
}

#### EXT FIG 8 ----

for(j in 1:2){
    
    if(j==1){s=160}
    if(j==2){s=100}
    
    if(j==1){ll="a" 
    ly = 48}
    if(j==2){ll="b" 
    ly = 150}
    
    fam_datasets[[j+2]] =  ggiNEXT(out[[j]], type=1) +   
        theme_bw() +
        xlab("Sample size (n)") +
        ylab("Number of unique viruses") +
        annotate(
            "text",
            x = 0,
            y = ly,
            label = ll,
            fontface = "bold",
            hjust = 0,
            vjust = 1,
            size = 2.46) +
        theme(legend.position = "none",
              panel.grid = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line(color = 'black'),
              axis.ticks = element_line(color = "black", linewidth = 0.5),
              axis.text = element_text(size = 6, color = "black"),
              axis.title = element_text(size = 7, color = "black"))
    
    fam_datasets[[j+2]]$layers = Filter(function(layer) !inherits(layer$geom, "GeomPoint"),
                                        fam_datasets[[j+2]]$layers)
    
    fam_labels = ggplot_build(fam_datasets[[j+2]])$data[[1]] %>% 
        filter(linetype=="22") %>%
        group_by(group) %>%
        filter(x == max(x)) 
    
    fam_labels$label = fam_datasets[[j]]
    
    fam_datasets[[j+2]] = fam_datasets[[j+2]] +     
        ggrepel::geom_text_repel(
            data = fam_labels,
            aes(x = x, y = y, label = label),
            seed = s,
            size = 2.11, 
            color = "black",
            force = 5,
            force_pull = 1,
            hjust = 0,
            nudge_y = 0, 
            max.overlaps = 3,
            direction = "y",
            segment.size = 0.2) 
    
    if(j==1){
        fam_datasets[[j+2]] = fam_datasets[[j+2]] +
            scale_color_discrete(breaks = fam_datasets[[j]],
                                 type = c("EMBALLONURIDAE" = "#8E24AA",
                                          "MEGADERMATIDAE" = "#8E24AA",
                                          "MINIOPTERIDAE" = "#8E24AA",
                                          "RHINOPOMATIDAE" = "#8E24AA",
                                          "SCIURIDAE" = "#219985",
                                          "SPALACIDAE" = "#219985")) +
            scale_fill_discrete(breaks = fam_datasets[[j]],
                                type = c("EMBALLONURIDAE" = "#8E24AA",
                                         "MEGADERMATIDAE" = "#8E24AA",
                                         "MINIOPTERIDAE" = "#8E24AA",
                                         "RHINOPOMATIDAE" = "#8E24AA",
                                         "SCIURIDAE" = "#219985",
                                         "SPALACIDAE" = "#219985"))
    }
    
    if(j==2){
        fam_datasets[[j+2]] = fam_datasets[[j+2]] +
            scale_color_discrete(breaks = fam_datasets[[j]],
                                 type = c("HIPPOSIDERIDAE" = "#8E24AA",
                                          "HYSTRICIDAE" = "#219985",
                                          "MOLOSSIDAE" = "#8E24AA",
                                          "MURIDAE" = "#219985",
                                          "PTEROPODIDAE" = "#8E24AA",
                                          "RHINOLOPHIDAE" = "#8E24AA",
                                          "SORICIDAE" = "#219985",
                                          "VESPERTILIONIDAE" = "#8E24AA")) +
            scale_fill_discrete(breaks = fam_datasets[[j]],
                                type = c("HIPPOSIDERIDAE" = "#8E24AA",
                                         "HYSTRICIDAE" = "#219985",
                                         "MOLOSSIDAE" = "#8E24AA",
                                         "MURIDAE" = "#219985",
                                         "PTEROPODIDAE" = "#8E24AA",
                                         "RHINOLOPHIDAE" = "#8E24AA",
                                         "SORICIDAE" = "#219985",
                                         "VESPERTILIONIDAE" = "#8E24AA"))
    }
}


fam_datasets[[3]]$layers[[1]]$aes_params$linewidth <- 0.6
fam_datasets[[4]]$layers[[1]]$aes_params$linewidth <- 0.6

combined_plot_ext8 = fam_datasets[[3]] + fam_datasets[[4]] +
    plot_layout(ncol = 2)

ggsave(paste0(output_dir, "FigExt8ab_2026_FINAL.pdf"),
       combined_plot_ext8,
       dpi = 600, 
       width = 180,   
       height = 76,   
       units = "mm",
       device = cairo_pdf)


#####################################
## 6. PARA & NON-PARA ESTIMATORS ####
## Supplementary Table 5           ##
#####################################

pred_strat <- pred_original 
unique_counts <- pred_strat %>% group_by(Taxa,ScientificName) %>% 
    summarise(Sample_Size = n_distinct(individual_id), Unique_Virus_Count = n_distinct(VirusName))
unique_counts <- unique_counts[unique_counts$ScientificName!="",]
list_excl <- unique_counts$ScientificName[unique_counts$Sample_Size<=10 | unique_counts$Unique_Virus_Count<=3]
pred_strat <- pred_strat[!(pred_strat$ScientificName %in% list_excl), ]
spp_list <- unique(pred_strat$ScientificName)
spp_list <- sort(spp_list)[-1]

## these species do not converge 
elements_to_remove <- c("Rhinolophus ferrumequinum", "Rhinolophus lepidus", "Taphozous melanopogon")
spp_list <- setdiff(spp_list, elements_to_remove)

## empty data frame for loop 
results_df <- data.frame(ScientificName = character(),
                         K = numeric(), b = numeric(), sample_needed = numeric(),
                         Hill_est=numeric(), Hill_LL=numeric(), Hill_UL=numeric(),
                         chao_est=numeric(),chao_LL=numeric(),chao_UL=numeric(),
                         jack1_est=numeric(),jack1_LL=numeric(),jack1_UL=numeric(),
                         jack2_est=numeric(),jack2_LL=numeric(),jack2_UL=numeric(),
                         boot_est=numeric(),boot_LL=numeric(),boot_UL=numeric(), stringsAsFactors = FALSE)

for (i in spp_list) {
    
    ## Extract information for each species 
    pred_spp <- pred_strat[pred_strat$ScientificName==i,]
    spp_name_i <- unique(pred_strat$ScientificName[pred_strat$ScientificName==i])
    pred_spp <- pred_spp[,c("SamplingDate","individual_id","TestResult","VirusName")]
    pred_spp_byanimal <- pred_spp %>% 
        filter(TestResult == "Positive" & VirusName!="") %>% 
        distinct(individual_id,VirusName) %>% 
        group_by(individual_id) %>% 
        summarise(NumberDistinctViruses = length(unique(VirusName)), 
                  ListDistinctViruses = knitr::combine_words(unique(VirusName),sep=";",and=""))
    pred_spp2 <- merge(pred_spp, pred_spp_byanimal, by="individual_id", all.x=TRUE)
    pred_spp2$NumberDistinctViruses[is.na(pred_spp2$NumberDistinctViruses)] <- 0
    pred_spp3 <- pred_spp2 %>% 
        distinct(individual_id,SamplingDate,NumberDistinctViruses,ListDistinctViruses)
    df <- separate_rows(pred_spp3, ListDistinctViruses, sep = ";\\s*")
    df2 <- spread(df, ListDistinctViruses, ListDistinctViruses)
    df2[, -(1:3)] <- lapply(df2[, -(1:3)], text_to_binary)
    df2 <- df2 %>% arrange(SamplingDate)
    df3 <- df2[,-c(1:3)]
    
    ## Extract MM parameters 
    accum_curve <- specaccum(df3, method = "exact")
    richness <- accum_curve$richness
    sample_size <- accum_curve$sites
    model <- nlsLM(richness ~ K * sample_size / (b + sample_size), 
                   data = data.frame(richness, sample_size), 
                   start = list(K = max(richness), b = 1))
    K_i <- coef(model)["K"];   b_i <- coef(model)["b"] 
    extrapolated_sample_size <- seq(max(sample_size), max(sample_size) + 10000, by = 10)
    extrapolated_richness <- predict(model, newdata = data.frame(sample_size = extrapolated_sample_size))
    asymptote_90 <- max(extrapolated_richness)*.9
    combined_sample_size <- c(sample_size, extrapolated_sample_size)
    combined_richness <- c(richness, extrapolated_richness)
    sample_size_90_i <- combined_sample_size[which.min(abs(combined_richness - asymptote_90))]
    
    ## Extract Hill parameters 
    matrix_list <- list()
    df4 <- as.matrix(df3)
    rownames(df4) <- df2$individual_id
    df5 <- t(df4)
    matrix_list[[1]] <- df5
    out <- iNEXT(matrix_list, q=0, datatype="incidence_raw")
    Hill_est_i <- out$AsyEst[1,4]
    Hill_LL_i <- out$AsyEst[1,6]
    Hill_UL_i <- out$AsyEst[1,7]
    
    ## Extract nonparametric estimators
    sp1_pool <- poolaccum(df2, permutations = 5000)
    chao_est_i <- max(sp1_pool$chao)
    chao_LL_i <- max(sp1_pool$chao)- 1.96*sd(sp1_pool$chao)
    chao_UL_i <- max(sp1_pool$chao)+ 1.96*sd(sp1_pool$chao)
    jack1_est_i <- max(sp1_pool$jack1)
    jack1_LL_i <- max(sp1_pool$jack1)- 1.96*sd(sp1_pool$jack1)
    jack1_UL_i <- max(sp1_pool$jack1)+ 1.96*sd(sp1_pool$jack1)
    jack2_est_i <- max(sp1_pool$jack2)
    jack2_LL_i <- max(sp1_pool$jack2)- 1.96*sd(sp1_pool$jack2)
    jack2_UL_i <- max(sp1_pool$jack2)+ 1.96*sd(sp1_pool$jack2)
    boot_est_i <- max(sp1_pool$boot)
    boot_LL_i <- max(sp1_pool$boot)- 1.96*sd(sp1_pool$boot)
    boot_UL_i <- max(sp1_pool$boot)+ 1.96*sd(sp1_pool$boot)
    
    ## Add all results to data frame
    results_df <- rbind(results_df, data.frame(
        ScientificName = spp_name_i, 
        K=K_i,b=b_i,sample_needed=sample_size_90_i, 
        Hill_est=Hill_est_i, Hill_LL=Hill_LL_i, Hill_UL=Hill_UL_i,
        chao_est=chao_est_i, chao_LL=chao_LL_i, chao_UL=chao_UL_i,
        jack1_est=jack1_est_i,jack1_LL=jack1_LL_i,jack1_UL=jack1_UL_i,
        jack2_est=jack2_est_i,jack2_LL=jack2_LL_i,jack2_UL=jack2_UL_i,
        boot_est=boot_est_i,boot_LL=boot_LL_i,boot_UL=boot_UL_i ))
} 

results_df <- results_df %>% 
    mutate(across(where(is.numeric), ~round(.x, digits = 0)))

#### Supplementary Table 5 ----
write.csv(results_df, paste0(output_dir, "supplement_section6_table.csv"), row.names = F)

View(pred %>% filter(ScientificName %in% results_df$ScientificName) %>% 
         filter(VirusName!="") %>% 
         group_by(VirusFamilyTested) %>% 
         summarise(n_virus = n_distinct(VirusName), n_species = n_distinct(ScientificName)))


####################################################
## 7. MICHELIS-MENTEN APPROACH & VISUALIZATIONS ####
## Figures 4a-h                                   ##
####################################################

pred_strat <- pred_original 
unique_counts <- pred_strat %>% group_by(Taxa,ScientificName) %>% 
    summarise(Sample_Size = n_distinct(individual_id), Unique_Virus_Count = n_distinct(VirusName))
unique_counts <- unique_counts[unique_counts$ScientificName!="",]

list_excl <- unique_counts$ScientificName[unique_counts$Sample_Size<=10 | unique_counts$Unique_Virus_Count<=3]

pred_strat <- pred_strat[!(pred_strat$ScientificName %in% list_excl), ]
spp_list <- unique(pred_strat$ScientificName)
spp_list <- sort(spp_list)[-1]

## these species do not converge 
elements_to_remove <- c("Rhinolophus ferrumequinum", "Rhinolophus lepidus", "Taphozous melanopogon")
spp_list <- setdiff(spp_list, elements_to_remove)

## empty data frame for loop 
results_df_mm <- data.frame(ScientificName = character(),asymptote = numeric(), 
                            asymptote90 = numeric(), sample_needed = numeric(), 
                            K = numeric(), b = numeric(),stringsAsFactors = FALSE)

for (i in spp_list) {
    pred_spp <- pred_strat[pred_strat$ScientificName==i,]
    pred_spp <- pred_spp[,c("SamplingDate","individual_id","TestResult","VirusName")]
    pred_spp_byanimal <- pred_spp %>% filter(TestResult == "Positive" & VirusName!="") %>% 
        distinct(individual_id,VirusName) %>% group_by(individual_id) %>% 
        summarise(NumberDistinctViruses = length(unique(VirusName)), 
                  ListDistinctViruses = knitr::combine_words(unique(VirusName),sep=";",and=""))
    pred_spp2 <- merge(pred_spp, pred_spp_byanimal, by="individual_id", all.x=TRUE)
    pred_spp2$NumberDistinctViruses[is.na(pred_spp2$NumberDistinctViruses)] <- 0
    pred_spp3 <- pred_spp2 %>% 
        distinct(individual_id,SamplingDate,NumberDistinctViruses,ListDistinctViruses)
    df <- separate_rows(pred_spp3, ListDistinctViruses, sep = ";\\s*")
    df2 <- spread(df, ListDistinctViruses, ListDistinctViruses)
    text_to_binary <- function(x) {ifelse(is.na(x), 0, 1)}
    df2[, -(1:3)] <- lapply(df2[, -(1:3)], text_to_binary)
    df2 <- df2 %>% arrange(SamplingDate)
    df3 <- df2[,-c(1:3)]
    
    accum_curve <- specaccum(df3, method = "exact")
    richness <- accum_curve$richness
    sample_size <- accum_curve$sites
    model <- nlsLM(richness ~ K * sample_size / (b + sample_size), data = data.frame(richness, sample_size), start = list(K = max(richness), b = 1))
    K_i <- coef(model)["K"]
    b_i <- coef(model)["b"] 
    extrapolated_sample_size <- seq(max(sample_size), max(sample_size) + 10000, by = 10)
    extrapolated_richness <- predict(model, newdata = data.frame(sample_size = extrapolated_sample_size))
    asymptote_90 <- max(extrapolated_richness)*.9  
    combined_sample_size <- c(sample_size, extrapolated_sample_size)
    combined_richness <- c(richness, extrapolated_richness)
    sample_size_90_i <- combined_sample_size[which.min(abs(combined_richness - asymptote_90))]
    maxRichness_i <- max(extrapolated_richness)
    spp_name_i <- unique(pred_strat$ScientificName[pred_strat$ScientificName==i])
    
    results_df_mm <- rbind(results_df_mm, 
                           data.frame(ScientificName = spp_name_i, 
                                      asymptote = maxRichness_i, 
                                      asymptote90 = asymptote_90, 
                                      sample_needed = sample_size_90_i,
                                      K=K_i,b=b_i))
} 

unique_counts <- pred_strat %>% group_by(Taxa, Family, ScientificName) %>% 
    summarise(PRED_N = n_distinct(individual_id), PRED_VR = n_distinct(VirusName))
unique_counts <- unique_counts[-1,]
results_df_mm2 <- merge(results_df_mm,unique_counts,by="ScientificName")
results_df_mm2$PropVR <- results_df_mm2$PRED_VR/results_df_mm2$asymptote
results_df_mm2$PropVR[results_df_mm2$PropVR>1] <- 1


#### FIG 4a -- Create scatterplot for VR/N with each species as one data point ------

results_df_mm2$face <- "italic"

fig4a = ggplot(results_df_mm2, aes(x = sample_needed, y = asymptote90, color=Taxa)) + 
    stat_smooth(aes(color=NULL, ymin = after_stat(pmax(ymin, 0))), 
                method = "gam", formula = y ~ s(x), se = T, color='gray28', linewidth = 0.5,) + 
    geom_point(aes(shape=factor(Taxa)),size=1) +  
    labs(x = "Sample size (n, log scale)", 
         y = "Number of unique viruses", title = "") + 
    theme_bw() + 
    scale_x_log10(breaks = c(100, 300, 1000, 3000, 10000)) + 
    scale_y_continuous(breaks=c(0,10,20,30,40)) + 
    scale_color_manual(values = c("#8E24AA","#219985")) +  
    theme(legend.position = 'none', 
          panel.grid = element_blank(), 
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          axis.line = element_blank(),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black")) +  
    coord_cartesian(ylim = c(0, NA))

gam_model <- gam(asymptote90~ s(sample_needed), data = results_df_mm2)
summary(gam_model)


labels_fig4a_a = c("Rousettus aegyptiacus", "Pteropus giganteus", "Lissonycteris angolensis", 
                   "Pteropus alecto", "Rhinolophus sinicus", "Miniopterus schreibersii", 
                   "Hipposideros lekaguli", "Hipposideros cervinus", "Triaenops persicus")
labels_fig4a_b = c("Mastomys natalensis", "Rattus rattus", "Eidolon helvum", "Hipposideros ruber", 
                   "Rousettus leschenaultii", "Rattus tanezumi", "Cynopterus sphinx", 
                   "Rhizomys pruinosus")

#### add labels --------

fig4a_text = fig4a +
    
    # above labels
    ggrepel::geom_text_repel(
        data = results_df_mm2 %>% filter(ScientificName %in% labels_fig4a_a) %>% 
            mutate(ScientificName = gsub("Rousettus aegyptiacus",
                                         "Rousettus\naegyptiacus",ScientificName),
                   ScientificName = gsub("Miniopterus schreibersii",
                                         "Miniopterus\nschreibersii",ScientificName)),
        aes(label = ScientificName), 
        size = 1.76,
        nudge_x = -0.25, nudge_y = 1,
        fontface = "italic",
        segment.color = "black", 
        direction = "both",
        seed = 100,
        point.padding = 0.2, box.padding = 0.1,
        force = 0.3, force_pull = 15, 
        max.overlaps = Inf,
        min.segment.length = 0.1,
        segment.size = 0.2,
        lineheight = 0.9) +
    
    # below labels
    ggrepel::geom_label_repel(
        data = results_df_mm2 %>% filter(ScientificName %in% labels_fig4a_b) %>% 
            mutate(ScientificName = gsub("Rousettus leschenaultii",
                                         "Rousettus\nleschenaultii",ScientificName),
                   ScientificName = gsub("Eidolon helvum",
                                         "Eidolon\nhelvum",ScientificName),
                   ScientificName = gsub("Cynopterus sphinx",
                                         "Cynopterus\nsphinx",ScientificName),
                   ScientificName = gsub("Rhizomys pruinosus",
                                         "Rhizomys\npruinosus",ScientificName),
                   ScientificName = gsub("Rattus rattus",
                                         "Rattus\nrattus",ScientificName)), 
        aes(label = ScientificName), 
        label.size = NA, fill = alpha("white", 0.5), label.padding = 0.05,
        size = 1.76,
        nudge_x = 0.16, 
        nudge_y = -0.01,
        fontface = "italic",
        segment.color = "black", 
        direction = "both",
        seed = 100,
        point.padding = 0.2, box.padding = 0.1,
        force = 0.4, force_pull = 0.5, 
        max.overlaps = Inf,
        min.segment.length = 0.01,
        segment.size = 0.2,
        lineheight = 0.9) +
    
    # plot label "a"
    annotate(
        "text",
        x = 84,  
        y = 43,  
        label = "a",
        fontface = "bold",
        size = 2.46)+
    theme(axis.text = element_text(size = 6),
          axis.title = element_text(size = 7))


#### FIG. 4b-h -- PLOT RESULTS AS INDIVIDUAL VAC CURVES -----

asymptotic_curve <- function(K, b, x) {return(K * x / (b + x))}
unique_families <- unique(results_df_mm2$Family)
n_families <- length(unique_families)
letter <- letters[which(letters == "b"):which(letters == "h")]

pdf(paste0(output_dir,"Fig4ah_2026_FINAL.pdf"), 
    width = 180/25.4, height = 150/25.4) 

mat <- matrix(
    c(0, 1, 2,
      0, 3, 4,
      5, 6, 7),
    nrow = 3, byrow = TRUE
)

layout(mat)

par(mar = c(1.5, 1.5, 1.5, 0.5)) 
par(mgp = c(2, 0.5, 0))
par(tcl = -0.3)
par(las = 1)
par(oma = c(2, 2, 0, 0))

filtered_families <- unique_families[sapply(unique_families, function(family) {
    n_unique_species <- length(unique(results_df_mm2$ScientificName[results_df_mm2$Family == family]))
    n_unique_species >= 3})]
filtered_families <- sort(filtered_families)
x_vals <- c(5500,2750,2250,8300,8800,3500,2500) # manually set x-axis maximum
y_vals <- c(45,18,13,9,30,21,10.5) # manually set y-axis maximum

up = c("Hipposideros ruber", "Miniopterus schreibersii", "Chaerephon pumilus", 
       "Mops condylurus", "Mastomys natalensis", "Rattus norvegicus", 
       "Lissonycteris angolensis", "Rhinolophus affinis", "Rhinolophus creaghi", 
       "Pipistrellus coromandra", "Scotophilus kuhlii")

rt = c("Aselliscus stoliczkanus", "Hipposideros cervinus", "Hipposideros lekaguli", 
       "Miniopterus inflatus", "Miniopterus magnater", "Miniopterus pusillus", 
       "Chaerephon plicatus", "Bandicota bengalensis", "Rattus argentiventer", 
       "Rattus tanezumi", "Eidolon helvum", "Pteropus alecto", "Pteropus giganteus", 
       "Pteropus lylei", "Rousettus leschenaultii", "Rhinolophus clivosus", "Myotis siligorensis",
       "Rhinolophus pusillus", "Rhinolophus sinicus",  
       "Myotis ricketti", "Tylonycteris pachypus", "Vespertilio sinensis")

dn = c("Hipposideros armiger", "Hipposideros galeritus", "Mus musculus", "Myotis horsfieldii",
       "Rattus rattus", "Acerodon celebensis", "Rousettus aegyptiacus")

all_labs = c(up, rt, dn)

not_labeled = c("Crocidura olivieri", "Cynopterus brachyotis", "Cynopterus sphinx", 
                "Eonycteris spelaea", "Epomophorus gambianus", "Epomops franqueti", 
                "Hipposideros caffer", "Hipposideros diadema", "Hipposideros gigas", 
                "Hipposideros larvatus", "Hystrix brachyura", "Megaderma lyra", 
                "Megaloglossus woermanni", "Micropteropus pusillus", "Myotis laniger", 
                "Pipistrellus pipistrellus", "Rattus exulans",
                "Rattus losea", "Rhinolophus mehelyi", "Rhinopoma hardwickii", 
                "Rhizomys pruinosus", "Rousettus amplexicaudatus", "Scotophilus leucogaster", 
                "Suncus murinus", "Triaenops persicus")

for (f in 1:length(filtered_families)) {
    
    family_data <- results_df_mm2[results_df_mm2$Family == filtered_families[f], ]
    plot(NULL,xlim=c(0,x_vals[f]),ylim=c(0,y_vals[f]),xlab="Sample Size", 
         ylab='Predicted Total Virus Richness',main=filtered_families[f],
         cex.main = 0.875, cex.axis = 0.75)
    
    for (i in 1:length(family_data$ScientificName)) {
        x <- seq(0, family_data$PRED_N[i], length.out = 100)
        y <- asymptotic_curve(family_data$K[i], family_data$b[i],x)
        species_family <- family_data$Family[family_data$ScientificName == family_data$ScientificName[i]]
        
        if(species_family=="MURIDAE" & family_data$ScientificName[i] %in% all_labs)
        {family_color <- adjustcolor("#219985", alpha.f = 0.7)}
        if(species_family=="MURIDAE" & !(family_data$ScientificName[i] %in% all_labs))
        {family_color <- adjustcolor("#219985", alpha.f = 0.3)}
        if(species_family!="MURIDAE" & family_data$ScientificName[i] %in% all_labs)
        {family_color <- adjustcolor("#7A168F", alpha.f = 0.7)} 
        if(species_family!="MURIDAE" & !(family_data$ScientificName[i] %in% all_labs))
        {family_color <- adjustcolor("#7A168F", alpha.f = 0.3)}
        
        lines(x,y,type='l',col=family_color,lwd=1)}
    
    for (i in 1:length(family_data$ScientificName)) {
        if(family_data$PRED_N[i]<100){
            x <- seq(family_data$PRED_N[i]+10, family_data$sample_needed[i], length.out = 100)}
        if(family_data$PRED_N[i]>=100){
            x <- seq(family_data$PRED_N[i]+50, family_data$sample_needed[i], length.out = 100)}
        y <- asymptotic_curve(family_data$K[i], family_data$b[i],x)
        species_family <- family_data$Family[family_data$ScientificName == family_data$ScientificName[i]]
        
        if(species_family=="MURIDAE" & family_data$ScientificName[i] %in% all_labs)
        {family_color <- adjustcolor("#219985", alpha.f = 0.7)}
        if(species_family=="MURIDAE" & !(family_data$ScientificName[i] %in% all_labs))
        {family_color <- adjustcolor("#219985", alpha.f = 0.3)}
        if(species_family!="MURIDAE" & family_data$ScientificName[i] %in% all_labs)
        {family_color <- adjustcolor("#7A168F", alpha.f = 0.7)}
        if(species_family!="MURIDAE" & !(family_data$ScientificName[i] %in% all_labs))
        {family_color <- adjustcolor("#7A168F", alpha.f = 0.3)}
        
        if(!family_data$ScientificName[i] %in% c("Scotophilus kuhlii", "Chaerephon plicatus")){ 
            lines(x,y,type='l',col=family_color,lty=2,lwd=1) 
            # IF statement removes EXTRAPOLATION when we sampled past 90% of the asymptote
        }
        
        endpoint_x <- max(x)
        endpoint_y <- max(y)
        
        if(species_family=="MURIDAE"){family_color <- "#219985"}
        if(species_family!="MURIDAE"){family_color <- "#8E24AA"} 
        
        ## spp. labels
        # UP
        calibrate::textxy(endpoint_x-x_vals[f]*0.01, endpoint_y, 
                          labs = ifelse(family_data$ScientificName[i] %in% up, 
                                        family_data$ScientificName[i], ""),
                          pos = 3, cex=0.625,
                          col = family_color, font = 3) 
        # RIGHT
        calibrate::textxy(endpoint_x, endpoint_y,
                          labs = ifelse(family_data$ScientificName[i] %in% rt, 
                                        family_data$ScientificName[i], ""),
                          pos = 4, cex=0.625,
                          col = family_color, font = 3)
        # DOWN
        calibrate::textxy(endpoint_x+x_vals[f]*0.08, endpoint_y,
                          labs = ifelse(family_data$ScientificName[i] %in% dn, 
                                        family_data$ScientificName[i], ""),
                          pos = 1, cex=0.625,
                          col = family_color, font = 3)
        
        lims <- par("usr")
        xpad <- diff(lims[1:2]) * 0.04 
        ypad <- diff(lims[3:4]) * 0.05  
        text(
            x = lims[1] + xpad,
            y = lims[4] - ypad,
            labels = letter[f],
            adj = c(0, 1),   
            cex = 0.875,     
            font = 1         
        )
    } 
}

#### FIG 4, combined ----

g <- ggplotGrob(fig4a_text)

pushViewport(
    viewport(
        width = 1/3 * 1.03 + 0.0015,
        height = 2/3, 
        x = 0.012, 
        y = 1/3 + 0.0225, 
        just = c("left", "bottom")
    ))

grid.draw(g)

popViewport()

pushViewport(
    viewport(
        x = 0,
        y = 0,
        width = 1,
        height = 1,
        just = c("left", "bottom")
    )
)

grid.text(
    "Sample size (n)",
    x = unit(0.53, "npc"),
    y = unit(0.02, "npc"),
    just = c("center", "bottom"),
    gp = gpar(
        fontsize = 7
    )
)

grid.text(
    "Number of unique viruses",
    x = unit(0.0275, "npc"),
    y = unit(0.11, "npc"),
    just = c("left", "center"),
    rot = 90,
    gp = gpar(
        fontsize = 7
    )
)

popViewport()

dev.off()

