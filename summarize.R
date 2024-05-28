library(tidyverse)
library(dplyr)

# import clean seine data
BB_seine_2019 <- read_csv(file = "clean/BB_seine_2019.csv")

# convert seine to character
BB_seine_2019$Seine <- as.character(BB_seine_2019$Seine)

# westside park has seven seines, all others have six seines
# remove 7th seine data from westside park site
BB_seine_2019 <- filter(BB_seine_2019, Seine != '7')

# create summaries of seine data for barplots -------------------------------

# calculate relative abundances of each Species per Site
# loop through a list of unique Site names to subset the df
# group_by Species and sum abundances of each Species, calculate sum of counts for each Site (N), 
# calculate relative abundances of each Species, check that sum of relative abundances per Site equal 100% 

# subset by columns of interest
BB_seine_2019_subset <- subset(BB_seine_2019, select = c(SiteCode, Species, Count))

# obtain unique list of seine names
site_unique <- unique(BB_seine_2019_subset$SiteCode)

# create list
seine_fish_table <- list()

for(i in seq_along(site_unique)) {
  site_subset=subset(BB_seine_2019_subset, SiteCode == site_unique[i])
  temp = site_subset %>% group_by(Species) %>% summarise_each(funs(sum (Count)))
  temp$SiteCode = site_unique[i]
  temp$N =  sum(temp$Count) 
  temp$RelAbund <- ( (temp$Count / temp$N) * 100 )
  temp$SumRelAbund <- sum(temp$RelAbund) 
  seine_fish_table[[i]] <- temp
  
  print(temp)
}

# bind items in seine_fish_table list with bind_rows
seine_fish_table <- bind_rows(seine_fish_table)

write.csv(seine_fish_table, file = "summary/BB_seine_2019_summary.csv",
          row.names = FALSE, quote = FALSE)

# create summaries of seine data for beta div -----------------------------

# species counts across seines 
BB_seine_2019_spp_total_count_summary <- BB_seine_2019 %>% 
  select(Sample, Species, Count) %>% 
  group_by(Species) %>% 
  summarize(total = sum(Count)) 

# summarize species counts per seine (multiple rows for measured individuals will be aggregated) 
BB_seine_2019_count_summary <- BB_seine_2019 %>% 
  select(Sample, Species, Count) %>% 
  group_by(Sample, Species) %>% 
  summarise(Count_seine = sum(Count))

write_csv(BB_seine_2019_count_summary, file = "summary/BB_seine_2019_count_summary.csv")

# repeat with only EWD sites
# summarize species counts per seine (multiple rows for measured individuals will be aggregated) 
additional_sites <- c("G", "H")
BB_seine_2019 <- filter(BB_seine_2019, !SiteCode %in% additional_sites)

BB_seine_2019_count_summary_EWD_sites <- BB_seine_2019 %>% 
  select(Sample, Species, Count) %>% 
  group_by(Sample, Species) %>% 
  summarise(Count_seine = sum(Count))

write_csv(BB_seine_2019_count_summary_EWD_sites, file = "summary/BB_seine_2019_count_summary_EWD_sites.csv")
