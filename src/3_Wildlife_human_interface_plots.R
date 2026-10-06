# manuscript: Johnson et al. 2026. Large-scale surveillance at 
#             wildlife-human interfaces reveals virus spillover  
#             risk in wildlife markets and trade supply chains. 
#             Nature Microbiology.
# analysis: Fig. 1 Plots

# Packages ---------------------------------------------------------------------

library(tidyverse)

# Load Data --------------------------------------------------------------------

fig1_data = list()
fig1_plots = list()

fig1_data$data1 = read.csv("data1_viruses_bats_rodents.csv", stringsAsFactors = F)

# Wildlife-human interfaces with potential for close (direct) contact 
# between humans and wildlife are labelled in bold in figures
fig1_data$bold = c("Hunted",
                   "For consumption",
                   "Private sale",
                   "In transit in trade supply chain",
                   "For sale in small market",
                   "For sale in medium market",
                   "For sale in large market",
                   "Sanctuary (or zoo)")

# model output for animal-human interfaces shown in Supplementary Tables 1.1-1.4
fig1_data$fig1b_interfaces = data.frame(
    suppl_table = c(1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 
                    1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 
                    1.3, 1.3, 1.3, 
                    1.4, 1.4, 1.4, 1.4),
    model_virus = c("All Viruses", "All Viruses", "All Viruses", "All Viruses", "All Viruses", 
                    "All Viruses", "All Viruses", "All Viruses", "All Viruses", "All Viruses", 
                    "All Viruses", "All Viruses", "All Viruses", "All Viruses", 
                    "Coronaviruses", "Coronaviruses", "Coronaviruses", "Coronaviruses", 
                    "Coronaviruses", "Coronaviruses", "Coronaviruses", "Coronaviruses", 
                    "Paramyxoviruses", "Paramyxoviruses", "Paramyxoviruses", 
                    "Coronaviruses", "Coronaviruses", "Coronaviruses", "Coronaviruses"),
    model_taxa = c("bats, rodents", "bats, rodents", "bats, rodents", "bats, rodents", "bats, rodents", 
                   "bats, rodents", "bats, rodents", "bats, rodents", "bats, rodents", "bats, rodents", 
                   "bats, rodents", "bats, rodents", "bats, rodents", "bats, rodents", 
                   "bats", "bats", "bats", "bats", "bats", "bats", "bats", "bats", "bats", "bats", "bats", 
                   "rodents", "rodents", "rodents", "rodents"),
    term = c("Dwellings", "Raiding crops", "Raiding markets", "Tourism site", "Wildlife management site", 
             "Guano farm", "Hunted", "For consumption", "Private sale", "In transit in trade supply chain", 
             "For sale in small market", "For sale in medium market", "For sale in large market", "Sanctuary (or zoo)", 
             "Dwellings", "Raiding crops", "Guano farm", "Private sale", "For sale in small market", 
             "For sale in large market", "In transit in trade supply chain", "Sanctuary (or zoo)", 
             "Dwellings", "For sale in large market", "Sanctuary (or zoo)", 
             "Dwellings", "For sale in medium market", "For sale in large market", "In transit in trade supply chain"),
    odds.ratio = c(0.857, 0.436, 0.293, 0.721, 0.649, 0.274, 0.559, 0.319, 0.343, 3.346, 2.895, 1.424, 3.355, 2.013, 
                   0.738, 0.506, 0.442, 0.225, 2.274, 11.025, 1.704, 1.468, 1.465, 3.682, 4.635, 3.206, 6.867, 2.947, 4.894),
    conf.low = c(0.764, 0.361, 0.072, 0.604, 0.439, 0.208, 0.462, 0.24, 0.231, 2.623, 1.9, 1.087, 2.582, 1.536, 
                 0.644, 0.399, 0.309, 0.107, 1.346, 7.188, 1.055, 1.109, 0.98, 1.396, 2.289, 1.985, 2.444, 1.341, 2.3),
    conf.high = c(0.962, 0.525, 1.195, 0.86, 0.958, 0.361, 0.677, 0.425, 0.511, 4.269, 4.411, 1.865, 4.359, 2.637, 
                  0.844, 0.643, 0.633, 0.473, 3.842, 16.91, 2.755, 1.943, 2.189, 9.706, 9.385, 5.179, 19.301, 6.476, 10.411),
    p.value = c("0.009", "< 0.001", "0.087", "< 0.001", "0.029", "< 0.001", "< 0.001", "< 0.001", "< 0.001", "< 0.001", 
                "< 0.001", "0.01", "< 0.001", "< 0.001", 
                "< 0.001", "< 0.001", "< 0.001", "< 0.001", "0.002", "< 0.001", "0.029", "0.007", 
                "0.062", "0.008", "< 0.001", "< 0.001", "< 0.001", "0.007", "< 0.001"),
    p.signif = c("*", "*", "", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", 
                 "*", "*", "*", "*", "*", "*", "*", "", "*", "*", "*", "*", "*", "*"))

# model output for specimen type shown in Supplementary Tables 1.1-1.4
fig1_data$fig1c_specimen_type = data.frame(  
    suppl_table = c(1.1, 1.1, 1.1, 1.1, 1.2, 1.2, 1.2, 1.2, 1.3, 1.3, 1.3, 1.3, 1.4, 1.4),
    model_virus = c("All Viruses", "All Viruses", "All Viruses", "All Viruses", 
                    "Coronaviruses", "Coronaviruses", "Coronaviruses", "Coronaviruses", 
                    "Paramyxoviruses", "Paramyxoviruses", "Paramyxoviruses", "Paramyxoviruses", 
                    "Coronaviruses", "Coronaviruses"),
    model_taxa = c("bats, rodents", "bats, rodents", "bats, rodents", "bats, rodents", 
                   "bats", "bats", "bats", "bats", 
                   "bats", "bats", "bats", "bats", 
                   "rodents", "rodents"),
    term = c("Feces", "Guano", "Rectal swab", "Urine", 
             "Feces", "Guano", "Rectal swab", "Urine", 
             "Guano", "Rectal swab", "Urine", "Urogenital swab", 
             "Feces", "Rectal swab"),
    odds.ratio = c(1.828, 2.112, 2.524, 7.413, 
                   2.174, 3.169, 4.932, 1.732, 
                   8.381, 2.933, 41.362, 14.221, 
                   2.052, 0.728),
    conf.low = c(1.407, 1.671, 2.313, 5.471, 
                 1.504, 2.387, 4.342, 1.004, 
                 3.696, 2.129, 22.9, 7.09, 
                 1.126, 0.611),
    conf.high = c(2.374, 2.668, 2.756, 10.044, 
                  3.143, 4.207, 5.602, 2.987, 
                  19.003, 4.039, 74.711, 28.525, 
                  3.739, 0.866),
    p.value = c("< 0.001", "< 0.001", "< 0.001", "< 0.001", 
                "< 0.001", "< 0.001", "< 0.001", "0.048", 
                "< 0.001", "< 0.001", "< 0.001", "< 0.001", 
                "0.019", "< 0.001"),
    p.signif = c("*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*"))

# model output for human analyses shown in Supplementary Tables 3.2
# load output directly generated with Human_models.R script
fig1_data$data2_model_ORs = read.csv("human_models_OR.csv", stringsAsFactors = F)


# Fig 1a Plot Data Prep ---------------------------------------------------

taxa_subtitle = paste0("</b><br><span style = 'font-size:6pt'>", "bats, rodents","</span>")

fig1_data$fig1a_plot_prep = fig1_data$data1 %>% 
    select(specimen_id, taxa_group, hunted:familytest_paramyxoviridae, -cr_flaviviridae, -familytest_flaviviridae) %>% 
    gather(key = AnimalHumanInterfaces, value = "InterfacePresent", c(hunted:raiding_livestock_food)) %>% 
    filter(InterfacePresent==1) %>% 
    mutate(cr_coronaviridae = ifelse(familytest_coronaviridae==0, NA, cr_coronaviridae),
           cr_orthomyxoviridae = ifelse(familytest_orthomyxoviridae==0, NA, cr_orthomyxoviridae),
           cr_paramyxoviridae = ifelse(familytest_paramyxoviridae==0, NA, cr_paramyxoviridae)) %>% 
    select(-c(familytest_coronaviridae, familytest_orthomyxoviridae, familytest_paramyxoviridae, InterfacePresent)) %>% 
    pivot_longer(cols = c(cr_coronaviridae, cr_orthomyxoviridae, cr_paramyxoviridae), 
                 names_to = "virus_family", values_to = "test_result") %>% 
    janitor::clean_names() %>% 
    filter(!animal_human_interfaces %in% c("crop_production", "raiding_livestock_food")) %>% 
    mutate(virus_family = gsub("cr_","", virus_family),
           virus_family = str_to_sentence(gsub("idae", "uses", virus_family)),
           virus_family = gsub("Orthomyxoviruses", "Influenza Viruses", virus_family), 
           animal_human_interfaces = str_to_sentence(gsub("_"," ",animal_human_interfaces)),
           animal_human_interfaces = gsub("Large markets","For sale in large market",animal_human_interfaces),
           animal_human_interfaces = gsub("Consumption","For consumption",animal_human_interfaces),
           animal_human_interfaces = gsub("Tourism","Tourism site",animal_human_interfaces),
           animal_human_interfaces = gsub("Wildlife management","Wildlife management site",animal_human_interfaces),
           animal_human_interfaces = gsub("Medium market","For sale in medium market",animal_human_interfaces),
           animal_human_interfaces = gsub("Zoo sanctuary","Sanctuary (or zoo)",animal_human_interfaces),
           animal_human_interfaces = gsub("Transit along valuechain","In transit in trade supply chain",animal_human_interfaces)) %>%
    mutate(animal_human_interfaces = ifelse(animal_human_interfaces %in% fig1_data$bold, 
                                            paste0("**", animal_human_interfaces, "**"), 
                                            as.character(animal_human_interfaces))) %>% 
    mutate(animal_human_interfaces = factor(animal_human_interfaces, 
                                            levels = c("Dwellings",
                                                       "Raiding crops",
                                                       "Raiding markets",
                                                       "Tourism site",
                                                       "Wildlife management site",
                                                       "Guano farm",
                                                       "**Hunted**",
                                                       "**For consumption**",
                                                       "**Private sale**",
                                                       "**In transit in trade supply chain**",
                                                       "**For sale in small market**",
                                                       "**For sale in medium market**",
                                                       "**For sale in large market**",
                                                       "**Sanctuary (or zoo)**"))) %>%
    group_by(animal_human_interfaces, taxa_group, virus_family) %>%
    summarise(test_result = mean(test_result, na.rm = T), n=n(), .groups = "drop") %>% 
    mutate(taxa_group = factor(taxa_group, levels = c("rodents & shrews", "bats")),
           virus_family = paste0("<b>", virus_family, taxa_subtitle),
           virus_family = factor(virus_family, levels = c(paste0("<b>", "Coronaviruses", taxa_subtitle),
                                                          paste0("<b>", "Paramyxoviruses", taxa_subtitle), 
                                                                 paste0("<b>", "Influenza Viruses", taxa_subtitle))),
           test_result = ifelse(test_result<0.001 & test_result!=0, 0.001, test_result)) # makes small non-zero percentages visible

# Fig 1a Plot -------------------------------------------------------------

fig1_plots$plot_1a = 
    ggplot(data = fig1_data$fig1a_plot_prep, 
           aes(x=animal_human_interfaces, y=test_result, fill=taxa_group)) + 
    geom_col(position = "dodge", width = 0.8) +
    scale_x_discrete(limits=rev) +
    scale_y_continuous(labels = scales::label_percent()) +
    coord_flip() +
    ggforce::facet_row(vars(virus_family), scales = 'fixed', space = 'free') +
    scale_fill_manual(values=c("rodents & shrews" = "#219985", 
                               "bats" = "#8E24AA")) + 
    ylab("Frequency of Virus Detection") + 
    xlab("") +
    theme(panel.grid = element_blank(), 
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25),
          panel.background = element_rect(fill = NA, color = NA),
          strip.background = element_rect(fill = NA, color = NA),
          legend.position = "none", 
          axis.ticks = element_line(color = "black", linewidth = 0.25),
          strip.text.x = ggtext::element_markdown(color = "black", size = 7),
          axis.title.x = element_text(color = "black", size = 7),
          axis.text = ggtext::element_markdown(color = "black", size = 5.5), 
          axis.text.y = ggtext::element_markdown(margin = margin(r = 4)),
          text = element_text(size = 7,  family = "sans"))

fig1_plots$plot_1a

# Fig 1b Plot Data Prep ---------------------------------------------------

fig1_data$fig1b_plot_prep = fig1_data$fig1b_interfaces %>% 
    mutate_at(vars(odds.ratio, conf.high, conf.low), log10) %>% 
    mutate(fill_OR = ifelse(grepl("*", p.signif, fixed=T), 
                            paste(model_taxa,"sig"), paste(model_taxa,"not sig")), 
           model = paste0("<b>", model_virus, "</b><br><span style = 'font-size:6pt'>", 
                          model_taxa,"</span>"),
           term2 = ifelse(term %in% fig1_data$bold, paste0("**",term,"**"), term)) %>% 
    mutate(model = factor(model, levels = unique(model)), 
           term2 = factor(term2, levels = unique(term2)))

# Fig 1b Plot -------------------------------------------------------------

fig1_plots$plot_1b = ggplot(fig1_data$fig1b_plot_prep, aes(y = term2, x = odds.ratio)) +
    geom_vline(lty=2, aes(xintercept=0), colour = 'gray50', linewidth=0.25)+
    geom_point(aes(fill = factor(fill_OR), color = factor(fill_OR)), size = 2, shape = 21) + 
    geom_segment(aes(x = conf.low, xend = conf.high, y = term2, yend = term2), 
                 colour = "black", linewidth=0.25) +
    scale_fill_manual(values = c("bats, rodents sig" = "#3487BB",
                                 "bats sig" = "#8E24AA","rodents sig" = "#219985", 
                                 "bats, rodents not sig" = "white", "bats not sig" = "white",
                                 "rodents not sig" = "white"))+
    scale_color_manual(values = c("bats, rodents sig" = "#3487BB","bats sig" = "#8E24AA",
                                  "rodents sig" = "#219985",
                                  "bats, rodents not sig" = "#3487BB", "bats not sig" = "#8E24AA",
                                  "rodents not sig" = "#219985"))+ 
    scale_y_discrete(limits = rev, breaks = fig1_data$fig1b_plot_prep$term2) + 
    facet_grid(~model, scales = "fixed", space = "fixed")+
    scale_x_continuous(breaks = c(-1, 0, 1, 2, 3), 
                       labels = c(0.1, 1, 10, 100, 1000), 
                       limits = c(-1.3, 2.3)) +
    ylab("") + 
    xlab("Odds Ratio (log scale)") +
    theme(panel.grid = element_blank(), 
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25),
          panel.background = element_rect(fill = NA, color = NA),
          strip.background = element_rect(fill = NA, color = NA),
          legend.position="none", 
          strip.text.x = ggtext::element_markdown(color = "black", size = 7),
          axis.title.x = element_text(colour = "black", size = 7),
          axis.ticks.x = element_line(color = "black", linewidth = 0.25),
          axis.text.x = element_text(color = "black", size = 5.5),
          axis.ticks.y = element_line(color = "black", linewidth = 0.25),
          axis.text.y = ggtext::element_markdown(color = "black", 
                                                 size = 5.5, margin = margin(r = 4))) +
    coord_cartesian(xlim = c(-1.3, 2.3),  ylim = c(1, 14), clip = "off")

fig1_plots$plot_1b

# Fig 1c Plot Data Prep ---------------------------------------------------

fig1_data$fig1c_plot_prep = fig1_data$fig1c_specimen_type %>% 
    mutate_at(vars(odds.ratio, conf.high, conf.low), log10) %>% 
    mutate(fill_OR = ifelse(grepl("*", p.signif, fixed=T), 
                            paste(model_taxa,"sig"), paste(model_taxa,"not sig"))) %>% 
    mutate(model = paste0("<b>", model_virus, "</b><br><span style = 'font-size:6pt'>", 
                          model_taxa,"</span>")) %>% 
    mutate(model = factor(model, levels = unique(model))) %>% 
    mutate(term = factor(term, levels = unique(term))) %>% 
    mutate(conf.h = ifelse(conf.high>1.3, 1.3, NA))

# Fig 1c Plot -------------------------------------------------------------

fig1_plots$plot_1c = ggplot(fig1_data$fig1c_plot_prep, aes(y = term, x = odds.ratio)) +
    geom_vline(lty=2, aes(xintercept=0), color = 'gray50', linewidth=0.25) +
    geom_point(aes(fill = factor(fill_OR), color = factor(fill_OR)), size = 2, shape = 21) + 
    geom_segment(aes(x = conf.low, xend = conf.high, y = term, yend = term), 
                 color = "black", linewidth=0.25) +
    scale_fill_manual(values = c("bats, rodents sig" = "#3487BB","bats sig" = "#8E24AA",
                                 "rodents sig" = "#219985",
                                 "bats, rodents not sig" = "white", "bats not sig" = "white",
                                 "rodents not sig" = "white"))+
    scale_color_manual(values = c("bats, rodents sig" = "#3487BB","bats sig" = "#8E24AA",
                                  "rodents sig" = "#219985",
                                  "bats, rodents not sig" = "#3487BB", "bats not sig" = "#8E24AA",
                                  "rodents not sig" = "#219985"))+ 
    scale_y_discrete(limits=rev, breaks = fig1_data$fig1c_plot_prep$term) + 
    facet_grid(~model, scales = "fixed", space="fixed") +
    scale_x_continuous(breaks = c(-1, 0, 1, 2, 3), 
                       labels = c(0.1, 1, 10, 100, 1000), 
                       limits = c(-1.3, 2.3)) +
    xlab("Odds Ratio (log scale)") + 
    ylab("") +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25),
          panel.background = element_rect(fill = NA, color = NA),
          strip.background = element_rect(fill = NA, color = NA),
          legend.position="none", 
          axis.title.x = element_text(color = "black", size = 7),
          axis.text.x = element_text(color = "black", size = 5.5),
          axis.ticks.x = element_line(color = "black", linewidth = 0.25),
          axis.ticks.y = element_line(color = "black", linewidth = 0.25),
          strip.text.x = ggtext::element_markdown(color = "black", size = 7),
          axis.text.y = element_text(color = "black", size = 5.5, margin = margin(r = 4))) +
    coord_cartesian(xlim = c(-1.3, 2.3),  ylim = c(1, 5), clip = "off")

fig1_plots$plot_1c


# Fig 1d Plot Data Prep ---------------------------------------------------

fig1_data$fig1d_plot_prep = fig1_data$data2_model_ORs %>% 
    mutate(model = factor(model, levels = c("Coronaviruses", "Paramyxoviruses", "Influenza Viruses", "Flaviviruses"))) %>% 
    mutate(model2 = paste0("<b>", model, "</b><br><span style = 'font-size:6pt'>", 
                           "humans","</span>")) %>% 
    mutate(model2 = factor(model2, levels = c(
        "<b>Coronaviruses</b><br><span style = 'font-size:6pt'>humans</span>",
        "<b>Paramyxoviruses</b><br><span style = 'font-size:6pt'>humans</span>",
        "<b>Influenza Viruses</b><br><span style = 'font-size:6pt'>humans</span>",
        "<b>Flaviviruses</b><br><span style = 'font-size:6pt'>humans</span>"))) %>% 
    mutate_at(vars(estimate, conf.high, conf.low), log10) %>% 
    mutate(fill_OR = ifelse(as.numeric(p.value)<0.05, "sig", "not sig")) %>% 
    mutate(term = gsub("_"," ",term)) %>% 
    mutate(term = gsub("acute fever","Acute fever",term)) %>% 
    mutate(term = gsub("age","Age ",term)) %>% 
    mutate(term = gsub("y","",term)) %>% 
    mutate(term = gsub("healthcare worker","Healthcare worker",term)) %>% 
    mutate(term = gsub("oral nasal","Oral/nasal",term)) %>% 
    mutate(term = gsub("blood serum","Blood/serum",term)) %>% 
    mutate(term = gsub("region asia","Asia",term)) %>% 
    mutate(term = factor(term, levels = c('Asia', 
                                          'Age 2-8', 'Age 9-17', 'Age 41+', 
                                          'Acute fever', 
                                          'Oral/nasal swab', 'Blood/serum', 
                                          'Healthcare worker')))%>% 
    mutate(conf.h = ifelse(conf.high>2.3, 2.3, NA))

# Fig 1d Plot -------------------------------------------------------------

fig1_plots$plot_1d = ggplot(fig1_data$fig1d_plot_prep, aes(y = term, x = estimate)) +
    geom_vline(lty=2, aes(xintercept=0), color = 'gray50', linewidth=0.25) +
    geom_point(aes(fill = factor(fill_OR), color = factor(fill_OR)), size = 2, shape = 21) + 
    geom_segment(aes(x = conf.low, xend = conf.high, y = term, yend = term), 
                 colour = "black", linewidth=0.25) +
    geom_segment(aes(x = conf.low, xend = conf.h, y = term, yend = term), 
                 colour = "black", linewidth=0.25, 
                 arrow = arrow(angle = 30, length = unit(1, "mm"), ends = "last", type = "open")) +
    scale_fill_manual(values = c("sig" = "#FA8775", "not sig" = "white")) +  
    scale_color_manual(values = c("sig" = "#FA8775", "not sig" = "#FA8775")) + 
    scale_y_discrete(limits = rev) +
    facet_grid(~model2, scales = "fixed", space="fixed") +
    scale_x_continuous(breaks = c(-1, 0, 1, 2), 
                       labels = c(0.1, 1, 10, 100), 
                       limits = c(-1.3, 2.3)) +
    xlab("Odds Ratio (log scale)") + 
    ylab("") +
    theme(panel.grid = element_blank(), 
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25),
          panel.background = element_rect(fill = NA, color = NA),
          strip.background = element_rect(fill = NA, color = NA),
          legend.position="none", 
          axis.ticks = element_line(color = "black", linewidth = 0.25),
          strip.text.x = ggtext::element_markdown(color = "black", size = 7),
          axis.title.x = element_text(color = "black", size = 7),
          axis.text = element_text(color = "black", size = 5.5), 
          axis.text.y = element_text(margin = margin(r = 4)),
          text = element_text(size = 7,  family="sans"))

fig1_plots$plot_1d

# Layout Packages ---------------------------------------------------------

library(ggh4x)
library(patchwork)
library(grid)

# Combine all Fig 1 plots -------------------------------------------------

file_name = "plots/Final Figs/Fig1_2026_FINAL.pdf"

w = 180 # 180 mm full page width, 88 mm half width
h = 188 # Aim for 185 mm height
d = 34
r = 2.84

plot_1a_std = fig1_plots$plot_1a + force_panelsizes(cols = unit(d+11.95, "mm"), # 3 cols only, whereas rest are 4 cols
                                         rows = unit(14*r, "mm")) # 14 y-axis labels x2 bars

plot_1b_std = fig1_plots$plot_1b + force_panelsizes(cols = unit(d, "mm"),
                                         rows = unit(14*r, "mm")) # 14 y-axis labels

plot_1c_std = fig1_plots$plot_1c + force_panelsizes(cols = unit(d, "mm"),
                                         rows = unit(5*r, "mm")) # 5 y-axis labels

plot_1d_std = fig1_plots$plot_1d + force_panelsizes(cols = unit(d, "mm"),
                                         rows = unit(8*r, "mm")) # 8 y-axis labels

combined_fig1 = plot_1a_std / plot_1b_std / plot_1c_std / plot_1d_std 


# PDF and annotate --------------------------------------------------------

pdf(file_name, width = w/25.4, height = h/25.4)

grid.draw(patchworkGrob(combined_fig1))

grid.text(
    "a",
    x = unit(4, "mm"),
    y = unit(184, "mm"),
    just = c("right", "bottom"),
    gp = gpar(fontsize = 7, col = "black", fontface = "bold")
)

grid.text(
    "ANIMAL-HUMAN INTERFACE",
    x = unit(35.6, "mm"),
    y = unit(180.5, "mm"),
    just = c("right", "bottom"),
    gp = gpar(fontsize = 6, col = "black", fontface = "bold")
)

grid.text(
    "b",
    x = unit(4, "mm"),
    y = unit(126, "mm"),
    just = c("right", "bottom"),
    gp = gpar(fontsize = 7, col = "black", fontface = "bold")
)

grid.text(
    "ANIMAL-HUMAN INTERFACE",
    x = unit(35.6, "mm"),
    y = unit(122.5, "mm"),
    just = c("right", "bottom"),
    gp = gpar(fontsize = 6, col = "black", fontface = "bold")
)

grid.text(
    "c",
    x = unit(4, "mm"),
    y = unit(67.5, "mm"),
    just = c("right", "bottom"),
    gp = gpar(fontsize = 7, col = "black", fontface = "bold")
)

grid.text(
    "SPECIMEN TYPE",
    x = unit(35.6, "mm"),
    y = unit(64, "mm"),
    just = c("right", "bottom"),
    gp = gpar(fontsize = 6, col = "black", fontface = "bold")
)

grid.text(
    "d",
    x = unit(4, "mm"),
    y = unit(34.5, "mm"),
    just = c("right", "bottom"),
    gp = gpar(fontsize = 7, col = "black", fontface = "bold")
)

grid.text(
    "HUMAN RISK FACTORS",
    x = unit(35.6, "mm"),,
    y = unit(31, "mm"),
    just = c("right", "bottom"),
    gp = gpar(fontsize = 6, col = "black", fontface = "bold")
)


dev.off()
