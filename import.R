# updated to handle cases of Order NA for taxonomy at 0.90-0.99 assignments

library(qiime2R)
library(phyloseq)
library(dplyr)
library(tidyverse)

# import metadata, add and recode factors for consistency across dataframes
metadata_2019 <- read.table("metadata/fish_metadata_2019.txt", sep = "\t", header = TRUE) %>% 
  mutate(Year = as.integer("2019"))
metadata_2020 <- read.table("metadata/BB_water_metadata_2020.txt", sep = "\t", header = TRUE)
metadata_2020$SampleType <- recode_factor(metadata_2020$SampleType, 
                                negative_control = "neg_control") 
metadata_2020$Year <- as.integer(c("2020"))

# create shared minimal metadata for filtering
m19 <- metadata_2019 %>% 
  select(c(sampleid, SampleType, Year))

m20 <- metadata_2020 %>% 
  select(c(sampleid, SampleType, Year))

meta_19_20 <- rbind(m19, m20)
write.csv(meta_19_20, "metadata/meta_19_20_sample.csv")
write.table(meta_19_20, "metadata/meta_19_20_sample.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# create shared maximum metadata for plotting
meta_all <- bind_rows(metadata_2019, metadata_2020)
write.table(meta_all, "metadata/meta_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# taxonomy ----------------------------------------------------------------

# read in taxonomy files
tax_99 <-  read_qza("output/rep_seqs_2019_2020_dada2_de_novo_97_taxonomy_12S_16S_18S_99.qza")
head(tax_99$data)
tax_99 <- parse_taxonomy(tax_99$data)

tax_97 <-  read_qza("output/rep_seqs_2019_2020_dada2_de_novo_97_taxonomy_12S_16S_18S_97.qza")
head(tax_97$data)
tax_97 <- parse_taxonomy(tax_97$data)

tax_95 <-  read_qza("output/rep_seqs_2019_2020_dada2_de_novo_97_taxonomy_12S_16S_18S_95.qza")
head(tax_95$data)
tax_95 <- parse_taxonomy(tax_95$data)

tax_90 <-  read_qza("output/rep_seqs_2019_2020_dada2_de_novo_97_taxonomy_12S_16S_18S_90.qza")
head(tax_90$data)
tax_90 <- parse_taxonomy(tax_90$data)

tax_80 <-  read_qza("output/rep_seqs_2019_2020_dada2_de_novo_97_taxonomy_12S_16S_18S_80.qza")
head(tax_80$data)
tax_80 <- parse_taxonomy(tax_80$data)

# remove rows with unassigned and add column with percent similarity blastn
tax_99_red <- filter(tax_99, Kingdom != "Unassigned")
tax_99_red$Percent_similarity <- c("0.99")
tax_99_red <- filter(tax_99_red, !is.na(Order))

tax_97_red <- filter(tax_97, Kingdom != "Unassigned")
tax_97_red$Percent_similarity <- c("0.97")
tax_97_red <- filter(tax_97_red, !is.na(Order))

tax_95_red <- filter(tax_95, Kingdom != "Unassigned")
tax_95_red$Percent_similarity <- c("0.95")
tax_95_red <- filter(tax_95_red, !is.na(Order))

tax_90_red <- filter(tax_90, Kingdom != "Unassigned")
tax_90_red$Percent_similarity <- c("0.90")
tax_90_red <- filter(tax_90_red, !is.na(Order))

tax_80_red <- filter(tax_80, Kingdom != "Unassigned")
tax_80_red$Percent_similarity <- c("0.80")

# get hash id for each table
hash_99 <- row.names(tax_99_red)
hash_97 <- row.names(tax_97_red)
hash_95 <- row.names(tax_95_red)
hash_90 <- row.names(tax_90_red)
hash_80 <- row.names(tax_80_red)

# find intersection of 97 to 99 hash ids
hash_intersect_97 <- intersect(hash_99, hash_97)
length(hash_intersect_97)

# remove intersect hash ids from 97 df
tax_97_fill <- tax_97_red %>% filter(!row.names(tax_97_red) %in% hash_intersect_97)
nrow(tax_97_fill)
nrow(tax_99_red)
# bind rows
tax_99_bind <- bind_rows(tax_99_red, tax_97_fill)
nrow(tax_99_bind)

# hash id from tax_99_bind 
hash_id_bind <- row.names(tax_99_bind)
length(hash_id_bind)

# find intersection 95 to hash_id_bind 
hash_intersect_95 <- intersect(hash_id_bind, hash_95)
length(hash_intersect_95) 

# remove intersect hash ids from 95 df
tax_95_fill <- tax_95_red %>% filter(!row.names(tax_95_red) %in% hash_intersect_95)
nrow(tax_95_fill) 
nrow(tax_99_bind)
# bind rows
tax_99_bind <- bind_rows(tax_99_bind, tax_95_fill)
nrow(tax_99_bind)

# hash id from tax_99_bind
hash_id_bind <- row.names(tax_99_bind)
length(hash_id_bind)

# find intersection 90 to hash_id_bind 
hash_intersect_90 <- intersect(hash_id_bind, hash_90)
length(hash_intersect_90)

# remove intersect hash ids from 90 df
tax_90_fill <- tax_90_red %>% filter(!row.names(tax_90_red) %in% hash_intersect_90)
nrow(tax_90_fill)
nrow(tax_99_bind)
# bind rows
tax_99_bind <- bind_rows(tax_99_bind, tax_90_fill)
nrow(tax_99_bind)

# hash id from tax_99_bind
hash_id_bind <- row.names(tax_99_bind)
length(hash_id_bind)

# find intersection 80 to hash_id_bind 
hash_intersect_80 <- intersect(hash_id_bind, hash_80)
length(hash_intersect_80)

# remove intersect hash ids from 80 df
tax_80_fill <- tax_80_red %>% filter(!row.names(tax_80_red) %in% hash_intersect_80)
nrow(tax_80_fill)
nrow(tax_99_bind)
# bind rows
tax_99_bind <- bind_rows(tax_99_bind, tax_80_fill)
nrow(tax_99_bind)

# change column name to Domain
colnames(tax_99_bind)[1] <- "Domain"
# remove domain underscore from levels of Domain
tax_99_bind$Domain <- gsub(tax_99_bind$Domain, pattern = "d__", replacement = "")

df <- tax_99_bind

# if lower taxonomy is NA fill with higher taxonomy
for(row in 1:nrow(df)) {
  Species <- df[row, "Species"][[1]]
  Genus <- df[row, "Genus"][[1]]
  Family <- df[row, "Family"][[1]]
  Order <- df[row, "Order"][[1]]
  Class <- df[row, "Class"][[1]]
  if (is.na(Species)) {
    df[row, "Species"] <- Genus
    Species <- df[row, "Species"][[1]]
    if (is.na(Species)) {
      df[row, "Species"] <- Family
      Species <- df[row, "Species"][[1]]
      if (is.na(Species)) {
        df[row, "Species"] <- Order
        Species <- df[row, "Species"][[1]]
        if (is.na(Species)) {
          df[row, "Species"] <- Class
        }
      }
    }
  }
}

# create column with genus species and percent similarity
df$Species_percent_sim <- paste(df$Species, " ", df$Percent_similarity) 

tax_99_bind <- df

# filter to 0.97 or higher similarity 
tax_only_97_sim_or_higher <- tax_99_bind %>% filter(Percent_similarity %in% c("0.97", "0.99"))

# convert df to matrix
tax_matrix <- as.matrix(tax_99_bind, rownames = TRUE, colnames = TRUE)

# import qza to phyloseq --------------------------------------------------

#create phyloseq object (rarefied table to 20000)
Fish_phyloseq <- qza_to_phyloseq(
  features = "output/table_2019_2020_dada2_de_novo_97_unassigned_80_removed_r20000.qza",
  metadata = "metadata/meta_all.txt"
)

# not rarefied
Fish_phyloseq_non_rarefied <- qza_to_phyloseq(
  features = "output/table_2019_2020_dada2_de_novo_97_unassigned_80_removed.qza",
  metadata = "metadata/meta_all.txt"
)

# add tax_table to phyloseq object
tax_table(Fish_phyloseq) <- tax_matrix
view(tax_table(Fish_phyloseq))
Fish_phyloseq # 175 samples
tax_table(Fish_phyloseq_non_rarefied) <- tax_matrix
view(tax_table(Fish_phyloseq_non_rarefied))
Fish_phyloseq_non_rarefied # 186 samples

# filter by year
fish_2019 <- subset_samples(Fish_phyloseq, Year == "2019") # 91 samples (6 samples lost to rarefying, 2 neg control 0 counts, 1 drop out sample)
fish_2020 <- subset_samples(Fish_phyloseq, Year == "2020") # 84 samples (5 samples lost to rarefying, 2 neg control 0 counts, 1 drop out sample)
# prune OTUs from phyloseq object with zero counts across all samples within each dataset
fish_2019_sp <- filter_taxa(fish_2019, function(x) sum(x) > 0, TRUE)
fish_2020_sp <- filter_taxa(fish_2020, function(x) sum(x) > 0, TRUE)

# filter by year non-rarefied
fish_2019_not_rare <- subset_samples(Fish_phyloseq_non_rarefied, Year == "2019") 
fish_2020_not_rare <- subset_samples(Fish_phyloseq_non_rarefied, Year == "2020")
# prune OTUs from phyloseq object with zero counts across all samples within each dataset
fish_2019_not_rare_sp <- filter_taxa(fish_2019_not_rare, function(x) sum(x) > 0, TRUE)
fish_2020_not_rare_sp <- filter_taxa(fish_2020_not_rare, function(x) sum(x) > 0, TRUE)

#add metadata to phyloseq object
metadata_2019 <- metadata_2019 %>% column_to_rownames("sampleid")
sample_data(fish_2019_sp) <- metadata_2019
sample_data(fish_2019_not_rare_sp) <- metadata_2019
metadata_2020 <- metadata_2020 %>% column_to_rownames("sampleid")
sample_data(fish_2020_sp) <- metadata_2020
sample_data(fish_2020_not_rare_sp) <- metadata_2020

saveRDS(fish_2019_sp, "output/phyloseq_object/fish_2019_phyloseq_v2.rds")
saveRDS(fish_2020_sp, "output/phyloseq_object/fish_2020_phyloseq_v2.rds")
saveRDS(fish_2019_not_rare_sp, "output/phyloseq_object/fish_2019_nonrarefied_phyloseq_v2.rds")
saveRDS(fish_2020_not_rare_sp, "output/phyloseq_object/fish_2020_nonrarefied_phyloseq_v2.rds")

Fish_phyloseq <- qza_to_phyloseq(
  tree = "output/rep_seqs_2019_2020_dada2_de_novo_97_unassigned_80_removed_alignment_mask_rooted_tree.qza",
  features = "output/table_2019_2020_dada2_de_novo_97_unassigned_80_removed_r20000.qza",
  metadata = "metadata/meta_all.txt"
)
# add tax_table to phyloseq object
tax_table(Fish_phyloseq) <- tax_matrix

saveRDS(Fish_phyloseq, "output/phyloseq_object/fish_phyloseq_v2.rds")


# sequence variant import -------------------------------------------------

# create phyloseq object on non-clustered and non-rarefied data with 97% similarity
Fish_phyloseq_nonrarefied_nonclustered <- qza_to_phyloseq(
  features = "output/table_2019_2020_dada2.qza",
  taxonomy = "output/rep_seqs_2019_2020_dada2_taxonomy_12S_16S_18S_97.qza",
  metadata = "metadata/meta_all.txt"
)
saveRDS(Fish_phyloseq_nonrarefied_nonclustered, file = "output/phyloseq_object/Fish_phyloseq_nonrarefied_nonclustered.rds")
view(tax_table(Fish_phyloseq_nonrarefied_nonclustered))

# create phyloseq object on non-clustered but rarefied to 20000 with taxonomic assignments 97% similarity
Fish_phyloseq_nonclustered_rarefied <- qza_to_phyloseq(
  features = "output/table_2019_2020_dada2_r20000.qza",
  taxonomy = "output/rep_seqs_2019_2020_dada2_taxonomy_12S_16S_18S_97.qza",
  metadata = "metadata/meta_all.txt"
)
saveRDS(Fish_phyloseq_nonclustered_rarefied, file = "output/phyloseq_object/Fish_phyloseq_nonclustered_rarefied.rds")
view(tax_table(Fish_phyloseq_nonclustered_rarefied))