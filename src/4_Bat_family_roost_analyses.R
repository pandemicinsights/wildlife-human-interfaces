# manuscript: Johnson et al. 2026. Large-scale surveillance at 
#             wildlife-human interfaces reveals virus spillover  
#             risk in wildlife markets and trade supply chains. 
#             Nature Microbiology.
# analysis: Infection prevalence among bats by preferred roost habitat

# IUCN Redlist Data (v.2022-02) was used to determine bat cave use.
# The use of caves and similar habitats "cave use" was identified in records 
# coded with habitat categories 7.1 and 7.2. As not all records have habitat 
# codes, habitat notes were reviewed for indications that cave or mine use was 
# known, likely, or possible. Cave use included the use of cave entrances and 
# overhangs. A species was considered to use caves even if use was limited to 
# specific life cycle events or seasons. Variations in habitat use across 
# regions or by sex was not considered. In the absence of habitat codes 7.1 and
# 7.2, and any habitat notes indicating cave or similar habitat use, 
# a bat species was assumed to not use caves.
#
# 7.1 - Caves and Subterranean Habitats (non-aquatic) - Caves
# 7.2 - Caves and Subterranean Habitats (non-aquatic) - Other Subterranean Habitats

# Packages ---------------------------------------------------------------------

library(dplyr)

# Load Data --------------------------------------------------------------------

bat_data = read.csv("data1_viruses_bats_rodents.csv", stringsAsFactors = F) %>% 
    filter(!is.na(bat_cave_use)) %>% 
    mutate(family_group = ifelse(family == "PTEROPODIDAE", "Pteropodidae", "Not_Pteropodidae")) %>% 
    mutate(flaviviridae = ifelse(familytest_flaviviridae==0, NA, cr_flaviviridae)) %>%
    mutate(paramyxoviridae = ifelse(familytest_paramyxoviridae==0, NA, cr_paramyxoviridae)) %>% 
    mutate(orthomyxoviridae = ifelse(familytest_orthomyxoviridae==0, NA, cr_orthomyxoviridae)) %>%
    select(specimen_id, region, family_group, bat_cave_use, 
           flaviviridae, orthomyxoviridae, paramyxoviridae, 
           alphacoronavirus, betacoronavirus)

# Variable Prep ----------------------------------------------------------------

Var1 = c("Pteropodidae", "Not_Pteropodidae")
Var2 = c("flaviviridae", "orthomyxoviridae", "paramyxoviridae", 
         "alphacoronavirus", "betacoronavirus")

# Bat Family & Cave Use --------------------------------------------------------

Cave_Family_FisherExact = 
    data.frame(
        FisherExact_Pvalue = as.numeric(), 
        FisherExact_OR = as.numeric(), 
        VirusFamily = as.character(), 
        FamilyGroup = as.character(), 
        stringsAsFactors = F)

for(i in Var1){
    
    df = bat_data %>% filter(family_group==i)    
    
    print(i)
    
    for(j in Var2){
        
        print(j)      
        
        mytable=table(df$bat_cave_use, df[,j], useNA = "no")
        
        CT = suppressWarnings(
            gmodels::CrossTable(mytable, expected=T, prop.r=F, prop.c=F, 
                                prop.t=F, prop.chisq=F, fisher=T))    
        
        if(CT$prop.row[[1]]!=1){ 
            
            Cave_Family_FisherExact[nrow(Cave_Family_FisherExact)+1,] = 
                data.frame(
                CT$fisher.ts$p.value, 
                CT$fisher.ts$estimate[["odds ratio"]], 
                paste(names(df[j])), 
                paste(i),
                stringsAsFactors = F)
        }
        
    }}

Cave_Family_FisherExact = Cave_Family_FisherExact %>% 
    mutate(P.Adjust = p.adjust(FisherExact_Pvalue, 
                               method = "bonferroni", 
                               n = nrow(.)))

View(Cave_Family_FisherExact)

# Bat Family -------------------------------------------------------------------

Family_FisherExact = 
    data.frame(
        FisherExact_Pvalue = as.numeric(), 
        FisherExact_OR = as.numeric(), 
        VirusFamily = as.character(), 
        stringsAsFactors = F)

df = bat_data

for(j in Var2){
    
    print(j)      
    
    mytable=table(df$family_group, df[,j], useNA = "no")
    
    CT = suppressWarnings(
        gmodels::CrossTable(mytable, expected=T, prop.r=F, prop.c=F, 
                            prop.t=F, prop.chisq=F, fisher=T))    
    
    if(CT$prop.row[[1]]!=1){ 
        
        Family_FisherExact[nrow(Family_FisherExact)+1,] = data.frame(
            CT$fisher.ts$p.value, 
            CT$fisher.ts$estimate[["odds ratio"]], 
            paste(names(df[j])), 
            stringsAsFactors = F)
    }}


Family_FisherExact = Family_FisherExact %>% 
    filter(FisherExact_OR != 0) %>% 
    mutate(P.Adjust = p.adjust(FisherExact_Pvalue, 
                               method = "bonferroni", 
                               n = 3))

View(Family_FisherExact)
