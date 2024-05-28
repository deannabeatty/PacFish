conda activate qiime2-2020.11

# dada2 on 2019 data
qiime dada2 denoise-paired \
--i-demultiplexed-seqs output/2019/qza_files/import_trimmed.qza \
--p-trunc-len-f 200 \
--p-trunc-len-r 200 \
--o-representative-sequences output/2019/qza_files/rep_seqs_2019_dada2.qza \
--o-table output/2019/qza_files/table_2019_dada2.qza \
--o-denoising-stats output/2019/qza_files/stats_2019_dada2.qza

# dada2 on 2020 data
qiime dada2 denoise-paired \
--i-demultiplexed-seqs output/2020/qza_files/import_trimmed.qza \
--p-trunc-len-f 200 \
--p-trunc-len-r 200 \
--o-representative-sequences output/2020/qza_files/rep_seqs_2020_dada2.qza \
--o-table output/2020/qza_files/table_2020_dada2.qza \
--o-denoising-stats output/2020/qza_files/stats_2020_dada2.qza

# create qzv of summary stats
 qiime metadata tabulate \
 --m-input-file output/2019/qza_files/stats_2019_dada2.qza \
 --o-visualization output/2019/qza_files/stats_2019_dada2.qzv 
 qiime metadata tabulate \
 --m-input-file output/2020/qza_files/stats_2020_dada2.qza \
 --o-visualization output/2020/qza_files/stats_2020_dada2.qzv 
# create qzv of table
 qiime feature-table summarize \
 --i-table output/2019/qza_files/table_2019_dada2.qza \
 --o-visualization output/2019/qza_files/table_2019_dada2.qzv
 qiime feature-table summarize \
 --i-table output/2020/qza_files/table_2020_dada2.qza \
 --o-visualization output/2020/qza_files/table_2020_dada2.qzv
# create qzv of seqs
 qiime feature-table tabulate-seqs \
 --i-data output/2019/qza_files/rep_seqs_2019_dada2.qza \
 --o-visualization output/2019/qza_files/rep_seqs_2019_dada2.qzv 
 qiime feature-table tabulate-seqs \
 --i-data output/2020/qza_files/rep_seqs_2020_dada2.qza \
 --o-visualization output/2020/qza_files/rep_seqs_2020_dada2.qzv

# merge datasets
# merge tables
qiime feature-table merge \
--i-tables output/2019/qza_files/table_2019_dada2.qza output/2020/qza_files/table_2020_dada2.qza \
--o-merged-table output/table_2019_2020_dada2.qza
# merge seqs
qiime feature-table merge-seqs \
--i-data output/2019/qza_files/rep_seqs_2019_dada2.qza output/2020/qza_files/rep_seqs_2020_dada2.qza \
--o-merged-data output/rep_seqs_2019_2020_dada2.qza

# create qzv of merged table
qiime feature-table summarize \
--i-table output/table_2019_2020_dada2.qza \
--o-visualization output/table_2019_2020_dada2.qzv
# create qzv of merged seqs
qiime feature-table tabulate-seqs \
--i-data output/rep_seqs_2019_2020_dada2.qza \
--o-visualization output/rep_seqs_2019_2020_dada2.qzv

# classify sequences with blastn against ref database at 
# 0.97
qiime feature-classifier classify-consensus-blast \
--i-query output/rep_seqs_2019_2020_dada2.qza \
--i-reference-reads ref_database/Nov_2022_release/12S-16S-18S-seqs.qza \
--i-reference-taxonomy ref_database/Nov_2022_release/12S-16S-18S-tax.qza \
--o-classification output/rep_seqs_2019_2020_dada2_taxonomy_12S_16S_18S_97.qza \
--p-perc-identity 0.97 \
--p-evalue 0.00001

# denovo clustering at 97% similarity
qiime vsearch cluster-features-de-novo \
--i-table output/table_2019_2020_dada2.qza \
--i-sequences output/rep_seqs_2019_2020_dada2.qza \
--p-perc-identity 0.97 \
--o-clustered-table output/table_2019_2020_dada2_de_novo_97.qza \
--o-clustered-sequences output/rep_seqs_2019_2020_dada2_de_novo_97.qza

# create qzv of clustered table
qiime feature-table summarize \
--i-table output/table_2019_2020_dada2_de_novo_97.qza \
--o-visualization output/table_2019_2020_dada2_de_novo_97.qzv

# classify OTUs with blastn against ref database at 
# 0.99
qiime feature-classifier classify-consensus-blast \
--i-query output/rep_seqs_2019_2020_dada2_de_novo_97.qza \
--i-reference-reads ref_database/Nov_2022_release/12S-16S-18S-seqs.qza \
--i-reference-taxonomy ref_database/Nov_2022_release/12S-16S-18S-tax.qza \
--o-classification output/rep_seqs_2019_2020_dada2_de_novo_97_taxonomy_12S_16S_18S_99.qza \
--p-perc-identity 0.99 \
--p-evalue 0.00001

# create qzv of seqs with classification labels at 0.99 percent similarity
qiime metadata tabulate \
--m-input-file output/rep_seqs_2019_2020_dada2_de_novo_97_taxonomy_12S_16S_18S_99.qza \
--o-visualization output/rep_seqs_2019_2020_dada2_de_novo_97_taxonomy_12S_16S_18S_99.qzv

# 0.97
qiime feature-classifier classify-consensus-blast \
--i-query output/rep_seqs_2019_2020_dada2_de_novo_97.qza \
--i-reference-reads ref_database/Nov_2022_release/12S-16S-18S-seqs.qza \
--i-reference-taxonomy ref_database/Nov_2022_release/12S-16S-18S-tax.qza \
--o-classification output/rep_seqs_2019_2020_dada2_de_novo_97_taxonomy_12S_16S_18S_97.qza \
--p-perc-identity 0.97 \
--p-evalue 0.00001

# 0.95
qiime feature-classifier classify-consensus-blast \
--i-query output/rep_seqs_2019_2020_dada2_de_novo_97.qza \
--i-reference-reads ref_database/Nov_2022_release/12S-16S-18S-seqs.qza \
--i-reference-taxonomy ref_database/Nov_2022_release/12S-16S-18S-tax.qza \
--o-classification output/rep_seqs_2019_2020_dada2_de_novo_97_taxonomy_12S_16S_18S_95.qza \
--p-perc-identity 0.95 \
--p-evalue 0.00001

# 0.90
qiime feature-classifier classify-consensus-blast \
--i-query output/rep_seqs_2019_2020_dada2_de_novo_97.qza \
--i-reference-reads ref_database/Nov_2022_release/12S-16S-18S-seqs.qza \
--i-reference-taxonomy ref_database/Nov_2022_release/12S-16S-18S-tax.qza \
--o-classification output/rep_seqs_2019_2020_dada2_de_novo_97_taxonomy_12S_16S_18S_90.qza \
--p-perc-identity 0.90 \
--p-evalue 0.00001

# 0.80
qiime feature-classifier classify-consensus-blast \
--i-query output/rep_seqs_2019_2020_dada2_de_novo_97.qza \
--i-reference-reads ref_database/Nov_2022_release/12S-16S-18S-seqs.qza \
--i-reference-taxonomy ref_database/Nov_2022_release/12S-16S-18S-tax.qza \
--o-classification output/rep_seqs_2019_2020_dada2_de_novo_97_taxonomy_12S_16S_18S_80.qza \
--p-perc-identity 0.80 \
--p-evalue 0.00001

# create qzv of seqs with classification labels at 0.80 percent similarity
qiime metadata tabulate \
--m-input-file output/rep_seqs_2019_2020_dada2_de_novo_97_taxonomy_12S_16S_18S_80.qza \
--o-visualization output/rep_seqs_2019_2020_dada2_de_novo_97_taxonomy_12S_16S_18S_80.qzv

# remove taxa unassigned at 0.80 
# from the count table 
qiime taxa filter-table \
--i-table output/table_2019_2020_dada2_de_novo_97.qza \
--i-taxonomy output/rep_seqs_2019_2020_dada2_de_novo_97_taxonomy_12S_16S_18S_80.qza \
--p-exclude Unassigned,Bacteria  \
--o-filtered-table output/table_2019_2020_dada2_de_novo_97_unassigned_80_removed.qza
# from the sequence table
qiime taxa filter-seqs \
--i-sequences output/rep_seqs_2019_2020_dada2_de_novo_97.qza \
--i-taxonomy output/rep_seqs_2019_2020_dada2_de_novo_97_taxonomy_12S_16S_18S_80.qza \
--p-exclude Unassigned,Bacteria \
--o-filtered-sequences output/rep_seqs_2019_2020_dada2_de_novo_97_unassigned_80_removed.qza 

# create qzv of filtered table
qiime feature-table summarize \
--i-table output/table_2019_2020_dada2_de_novo_97_unassigned_80_removed.qza \
--o-visualization output/table_2019_2020_dada2_de_novo_97_unassigned_80_removed.qzv
# create qzv of filtered seqs with summary statistics
qiime feature-table tabulate-seqs \
--i-data output/rep_seqs_2019_2020_dada2_de_novo_97_unassigned_80_removed.qza \
--o-visualization output/rep_seqs_2019_2020_dada2_de_novo_97_unassigned_80_removed_tabulate.qzv
# create qzv of filtered seqs 
qiime metadata tabulate \
--m-input-file output/rep_seqs_2019_2020_dada2_de_novo_97_unassigned_80_removed.qza \
--o-visualization output/rep_seqs_2019_2020_dada2_de_novo_97_unassigned_80_removed.qzv

# mafft alignment and fasttree with midpoint root for phylogenetic tree
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences output/rep_seqs_2019_2020_dada2_de_novo_97_unassigned_80_removed.qza \
--o-alignment output/rep_seqs_2019_2020_dada2_de_novo_97_unassigned_80_removed_alignment.qza \
--o-masked-alignment output/rep_seqs_2019_2020_dada2_de_novo_97_unassigned_80_removed_alignment_mask.qza \
--o-tree output/rep_seqs_2019_2020_dada2_de_novo_97_unassigned_80_removed_alignment_mask_unrooted_tree.qza \
--o-rooted-tree output/rep_seqs_2019_2020_dada2_de_novo_97_unassigned_80_removed_alignment_mask_rooted_tree.qza

# rarefy to even sampling depth
qiime feature-table rarefy \
--i-table output/table_2019_2020_dada2_de_novo_97_unassigned_80_removed.qza \
--p-sampling-depth 20000 \
--o-rarefied-table output/table_2019_2020_dada2_de_novo_97_unassigned_80_removed_r20000.qza
