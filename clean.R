library(tidyverse)
library(dplyr)
library(readxl)

# clean seine data --------------------------------------------------------------

# Bodega seine surveys 2019 & 2021
BB_seine_all <-  read_csv(file = "raw/Master Fish list.csv")

# subset to 2019
BB_seine_2019 <- BB_seine_all[1:809,1:11]

# rename column
BB_seine_2019 <- rename(BB_seine_2019, Length_mm = 'Length (mm)')

# create column SiteCode
BB_seine_2019$SiteCode <- BB_seine_2019$Subsite
BB_seine_2019$SiteCode <- as.factor(BB_seine_2019$SiteCode)
# recode with EWD SiteCode levels
BB_seine_2019 <- BB_seine_2019 %>% mutate(SiteCode=recode(SiteCode, 
                                        `WP` = c("A"), `MM` = c("B"), `NC` = c("C"), `CC` = c("D"), `MP` = c("E"), `BL` = c("F"), `SB` = c("G"), `SL` = c("H"))) 

# rename column Species to Common_name
BB_seine_2019 <- rename(BB_seine_2019, Common_name = Species)
# create column Species to fill with genus species names
BB_seine_2019$Species <- BB_seine_2019$Common_name
# recode with genus species names
BB_seine_2019 <- BB_seine_2019 %>% mutate(Species=recode(Species, 
                `Arrow goby` = c("Clevelandia ios"), `Black perch` = c("Embiotoca jacksoni"), `English sole` = c("Parophrys vetulus"), 
                `Kelpfish` = c("Gibbonsia metzi"), `Pile perch` = c("Phanerodon vacca"), `Pipefish` = c("Syngnathus leptorhynchus"),
                `Saddleback gunnel` = c("Pholis ornata"), `Shiner perch` = c("Cymatogaster aggregata"), `Staghorn sculpin` = c("Leptocottus armatus"),
                `Stickleback` = c("Gasterosteus aculeatus"), `Striped perch` = c("Embiotoca lateralis"), `Tubesnout` = c("Aulorhynchus flavidus"),
                `Black rockfish` = c("Sebastes melanops"), `Gopher rockfish` = c("Sebastes carnatus"), `Kelp perch` = c("Brachyistius frenatus"),
                `Lingcod` = c("Ophiodon elongatus"), `Penpoint gunnel` = c("Apodichthys flavidus"), `California halibut` = c("Paralichthys californicus"),
                `Topsmelt` = c("Atherinops affinis"), `Dwarf perch` = c("Micrometrus minimus"), `Fluffy sculpin` = c("Oligocottus snyderi"),
                `Midshipman` = c("Porichthys notatus"), `Yellowfin goby` = c("Acanthogobius flavimanus"), `Cabezon` = c("Scorpaenichthys marmoratus"), `Rock prickleback` = c("Xiphister mucosus"))) 

# correct common names according to fishbase website
# pipefish = bay pipefish; kelpfish = striped kelpfish; Staghorn sculpin = Pacific staghorn sculpin; Stickleback = three-spined stickleback
# Striped perch = Striped seaperch, Tubesnout = Tube-snout; topsmelt = topsmelt silverside; midshipman = plainfin midshipman 

write_csv(BB_seine_2019, file = "clean/BB_seine_2019.csv")

# site metadata for seine samples
BB_seine_2019_site_metadata <- BB_seine_2019 %>% 
  select(Sample, SiteCode) %>% 
  distinct(Sample, SiteCode) %>% 
  as.data.frame()

write_csv(BB_seine_2019_site_metadata, file = "clean/BB_seine_2019_site_metadata.csv")

# site metadata for seine samples with EWD metadata added
metadata_EWD_2019 <- read.table("metadata/fish_metadata_2019.txt", sep = "\t", header = TRUE) %>% 
  mutate(Year = as.integer("2019")) %>% 
  select(-sampleid,-LinkerPrimerSequence,-BarcodeSequence,-SampleType,-SiteCodeBBExtra,-DateCollected,-TotalWaterVolumeFiltered) %>% 
  filter(Region == "BB") %>% 
  distinct_all()

BB_seine_2019_EWD_metadata <- BB_seine_2019 %>% 
  select(Sample, SiteCode) %>% 
  distinct(Sample, SiteCode) %>% 
  left_join(metadata_EWD_2019, by = "SiteCode") %>% 
  as.data.frame()

write_csv(BB_seine_2019_EWD_metadata, file = "clean/BB_seine_2019_EWD_metadata.csv")
