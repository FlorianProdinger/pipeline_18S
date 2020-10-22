#!/bin/tcsh

#This is an 18S analysis pipeline.
#This script needs:
# 1)  Qiime2 (2018 November version)
# 2)  qsub
#   *) an error and output dir for Qsub (in qsub commands section)
# 3) a directory defined as "OUT_DIR"
#  files that have to be in the out dir:
#   *)metadata (sample info, manually created and named 'METADATA_FILE')
#   *)quiime2 manifest (created manually or with the  python script)
#  oprional: 4) an INPUT_DIR
#   *)this dir will be passed to a script to generate the manifest needed for qiime2
#
##!! DO NOT forget the "/" at the end of all the paths!!!##

#####################################################################################################################
#QSUB COMMANDS 
#PBS -q cdb 
#PBS -N Qiime2_18S_ASV
#PBS -l mem=50gb  
#PBS -l ncpus=8 
#PBS -o qsub_out.txt
#PBS -e qsub_error.txt
#PBS -m ea
#PBS -M florian.prodinger@gmx.net
#####################################################################################################################


#In the OUT_DIR all files of the pipeline will be generated.
set OUT_DIR="/lustre1/aptmp/florian/uranouchi/18S_uranouchi/UU_euk-18S/qiime/20200527_18S_analysis/pipeline_v3_ASVs_dada_20200625/"
#"/lustre1/aptmp/florian/uranouchi/18S_uranouchi/UU_euk-18S/qiime/pipeline_v3_ASV_86_files_97search_20200501/"
#"/lustre1/aptmp/florian/uranouchi/18S_uranouchi/UU_euk-18S/qiime/pipeline_v3_86_files_20191121_97search/"
#"/lustre1/aptmp/florian/uranouchi/18S_uranouchi/UU_euk-18S/qiime/pipeline_v3_86_files_20190820/"
#"/lustre1/aptmp/florian/uranouchi/18S_uranouchi/UU_euk-18S/qiime/pipeline_v3_ALL_files/"

#this file can be gerenated with a script (0 MANIFEST), but manuall creation is recommended.
# tutorial: https://docs.qiime2.org/2018.11/tutorials/importing/
set MANIFEST="qiime_manifest.txt"

#This file is needed for a barplot and other plots and has to be created manually
# tutorial: https://docs.qiime2.org/2018.11/tutorials/metadata/
set METADATA_FILE="metadata_sheet.tsv"

#Loading necessary modules fot his sript
eval `/usr/bin/modulecmd tcsh load qiime2/2020.2`
eval `/usr/bin/modulecmd tcsh load Python/3.6.5`
#source /etc/profile.d/modules.sh
#module qiime2/2020.2
#module Python/3.6.5




########## PIPELINE STARTS ##########


##  1  IMPORT
#This script imports the files specified in the "manifest" 
#echo "import using $MANIFEST"
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path ${OUT_DIR}$MANIFEST \
--output-path ${OUT_DIR}out1-import.qza \
--input-format PairedEndFastqManifestPhred33

#### DADA2 includes: 
# 2  PRIMER REMOVAL !!!check if length is okay!!!
# 3  MERGING
# 4  QUALITY FILTER
# 5  DEREPLICATION
# 7  CHIMERA CHECK

qiime dada2 denoise-paired \
 --i-demultiplexed-seqs ${OUT_DIR}out1-import.qza \
 --p-trunc-len-f 240 \
 --p-trunc-len-r 240 \
 --p-trim-left-f 21  \
 --p-trim-left-r 20  \
 --p-n-threads 8 \
 --o-table ${OUT_DIR}out5-derep-table.qza \
 --o-representative-sequences ${OUT_DIR}out5-derep-sequences.qza \
 --o-denoising-stats ${OUT_DIR}out_dada2_denoising-stats



# 6
############### SILVA IMPORT ################
echo "SILVA Sequence and taxonomy import"

#6.1) importing 97% sequence data
set SILVA_FILE_SEQ_97="/lustre1/aptmp/florian/uranouchi/18S_uranouchi/UU_euk-18S/qiime/Silva/SILVA_132_QIIME_release/rep_set/rep_set_18S_only/97/silva_132_97_18S.fna"
qiime tools import --type 'FeatureData[Sequence]' \
--input-path $SILVA_FILE_SEQ_97  \
--output-path ${OUT_DIR}out6.1-silva_sequence_0.97_18Sonly.qza

#6.2) import 97% taxonomy data
set SILVA_FILE_97="/lustre1/aptmp/florian/uranouchi/18S_uranouchi/UU_euk-18S/qiime/Silva/SILVA_132_QIIME_release/taxonomy/18S_only/97/majority_taxonomy_all_levels.txt"
qiime tools import  --type 'FeatureData[Taxonomy]' \
--input-path $SILVA_FILE_97 \
--input-format HeaderlessTSVTaxonomyFormat \
--output-path ${OUT_DIR}out6.2-silva_0.97_taxonomy_majority_out.qza.qza

#6.3) importing 99% sequence data
set SILVA_FILE_SEQ_99="/lustre1/aptmp/florian/uranouchi/18S_uranouchi/UU_euk-18S/qiime/Silva/SILVA_132_QIIME_release/rep_set/rep_set_18S_only/99/silva_132_99_18S.fna"
qiime tools import --type 'FeatureData[Sequence]' \
--input-path $SILVA_FILE_SEQ_99  \
--output-path ${OUT_DIR}out6.1-silva_sequence_0.99_18Sonly.qza

#6.4) import 99% taxonomy data
set SILVA_FILE_99="/lustre1/aptmp/florian/uranouchi/18S_uranouchi/UU_euk-18S/qiime/Silva/SILVA_132_QIIME_release/taxonomy/18S_only/99/majority_taxonomy_all_levels.txt"
qiime tools import  --type 'FeatureData[Taxonomy]' \
--input-path $SILVA_FILE_99 \
--input-format HeaderlessTSVTaxonomyFormat \
--output-path ${OUT_DIR}out6.2-silva_0.99_taxonomy_majority_out.qza.qza
################## SILVA IMPORT ################


# 8.3) tree generation
#   ALIGNMENT & TREE
echo "alignment & tree generation"
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences ${OUT_DIR}out5-derep-sequences.qza \
--o-alignment ${OUT_DIR}out8-seqs-nonchimeras_clustered_ASV-aligned.qza \
--o-masked-alignment ${OUT_DIR}out8-seqs-nonchimeras_clustered_ASV-masked-aligned.qza \
--o-tree ${OUT_DIR}out8-seqs-nonchimeras_clustered_ASV_unrooted-tree.qza \
--o-rooted-tree ${OUT_DIR}out8-seqs-nonchimeras_clustered_ASV_rooted-tree.qza

#  9  TAXONOMIC ANOTATION
#This commands annotates reads
echo "vsearch feature classifier 97% similarity"
qiime feature-classifier classify-consensus-vsearch \
--i-query ${OUT_DIR}out5-derep-sequences.qza \
--i-reference-taxonomy ${OUT_DIR}out6.2-silva_0.97_taxonomy_majority_out.qza.qza \
--i-reference-reads ${OUT_DIR}out6.1-silva_sequence_0.97_18Sonly.qza \
--o-classification ${OUT_DIR}out9-classify-vsearch_ASVOTUs_0.97SILVA.qza \
--p-perc-identity 0.90 \
--p-maxaccepts 1 \
--p-threads 8 

#  9.1 SINGELTON FILTER
echo "filtering singelton OTUs"
qiime feature-table filter-features \
--i-table ${OUT_DIR}out5-derep-table.qza \
--p-min-frequency 2 \
--p-where "Taxon NOT LIKE '%Unassigned%'" \
--m-metadata-file ${OUT_DIR}out9-classify-vsearch_ASVOTUs_0.97SILVA.qza \
--o-filtered-table ${OUT_DIR}out9-derep-clustered_ASV-final_table_no_singelton.qza
#--o-filtered-table ${OUT_DIR}out9-derep-clustered_ASV_no_singelton.qza 
#-p-min-frequency <- can be used to remove OTUs with less than 2 reads
# --p-min-samples 10 <- can be used to remove OTUs that were found in less than 10 samples



# ###\\\///   PLOTTING   \\\///### #

# 8.4)  TABLES for debbuging
echo "visualizing several tables only for PLOTTING!"
##    a) ASV - table
qiime feature-table summarize \
--i-table ${OUT_DIR}out9-derep-clustered_ASV-final_table_no_singelton.qza\
--o-visualization ${OUT_DIR}out9-derep-clustered_ASV-final_table_no_singelton_OTU.qzv\
--m-sample-metadata-file ${OUT_DIR}$METADATA_FILE
##    b) ASV - sequence list
qiime feature-table tabulate-seqs \
--i-data ${OUT_DIR}out5-derep-sequences.qza \
--o-visualization ${OUT_DIR}out8.4-seqs-nonchimeras_clustered_ASV_OTU_Seqs.qzv

##    c) taxonomic read table as tsv and   ##   f) OTU table with reads as tsv
  mkdir -p ${OUT_DIR}export_for_OTU_table
  foreach LEVEL (2 3 4 5 6 7 8 9 10)
  qiime taxa collapse \
  --i-table ${OUT_DIR}out9-derep-clustered_ASV-final_table_no_singelton.qza \
  --i-taxonomy ${OUT_DIR}out9-classify-vsearch_ASVOTUs_0.97SILVA.qza\
  --p-level $LEVEL\
  --o-collapsed-table ${OUT_DIR}export_for_OTU_table/out9-derep-clustered_ASV-final_table_no_singelton_taxonomy_lvl{$LEVEL}.qza
  qiime tools export --input-path ${OUT_DIR}export_for_OTU_table/out9-derep-clustered_ASV-final_table_no_singelton_taxonomy_lvl{$LEVEL}.qza\
  --output-path ${OUT_DIR}export_for_OTU_table
  /bin/mv -f ${OUT_DIR}export_for_OTU_table/feature-table.biom ${OUT_DIR}export_for_OTU_table/out9-derep-clustered_ASV-final_table_no_singelton_taxonomy_lvl{$LEVEL}_feature-table.biom
  biom convert -i ${OUT_DIR}export_for_OTU_table/out9-derep-clustered_ASV-final_table_no_singelton_taxonomy_lvl{$LEVEL}_feature-table.biom\
  -o ${OUT_DIR}export_for_OTU_table/out9-derep-clustered_ASV-final_table_no_singelton_taxonomy_lvl{$LEVEL}.tsv\
  --to-tsv
  end

mkdir -p ${OUT_DIR}export_for_OTU_table
qiime tools export \
 --input-path ${OUT_DIR}out9-derep-clustered_ASV-final_table_no_singelton.qza \
 --output-path ${OUT_DIR}export_for_OTU_table

/bin/mv -f ${OUT_DIR}export_for_OTU_table/feature-table.biom \
${OUT_DIR}export_for_OTU_table/out9-derep-clustered_ASV-final_table_no_singelton_feature-table.biom

biom convert \
 -i ${OUT_DIR}export_for_OTU_table/out9-derep-clustered_ASV-final_table_no_singelton_feature-table.biom \
 -o ${OUT_DIR}export_for_OTU_table/out9-derep-clustered_ASV-final_table_no_singelton_OTUs.tsv \
 --to-tsv



##   d) OTU table with assigned taxonomy

echo "exporting classification (OTU - tax list) from qiime..."
qiime tools export \
 --input-path ${OUT_DIR}out9-classify-vsearch_ASVOTUs_0.97SILVA.qza \
 --output-path ${OUT_DIR}export_tax_OTU_list
echo "exporting OTU table without annotation from qiime..." #this command fails in qsub
qiime tools export \
 --input-path ${OUT_DIR}out5-derep-table.qza\
 --output-path ${OUT_DIR}export_tax_OTU_list

#mkdir -p ${OUT_DIR}export_tax_OTU_list

echo "using 'command_rename_first_line_taxonomy_file.R' to rename the classification (OTU - tax) list header"
echo `ls ${OUT_DIR}export_tax_OTU_list/taxonomy.tsv`
#replace fist line of taxonomy.tsv with "#OTUID       taxonomy        confidence"
Rscript /lustre1/aptmp/florian/uranouchi/18S_uranouchi/UU_euk-18S/qiime/pipeline_commands/command_rename_first_line_taxonomy_file.R ${OUT_DIR}export_tax_OTU_list/taxonomy.tsv ASV


echo "adding metadata (taxonomy) to the OTU table (feature-table.biom)"
biom add-metadata \
 -i ${OUT_DIR}export_tax_OTU_list/feature-table.biom \
 --observation-metadata-fp ${OUT_DIR}export_tax_OTU_list/taxonomy.tsvASV_renamed.tsv \
 --sc-separated taxonomy \
 -o ${OUT_DIR}export_tax_OTU_list/table-with-taxonomy_ASV.biom

echo "exporting the OTU table (feature-table.biom) WITH the added classification (tax)"
biom convert \
 -i ${OUT_DIR}export_tax_OTU_list/table-with-taxonomy_ASV.biom \
 --to-tsv \
 --header-key taxonomy \
 -o ${OUT_DIR}export_tax_OTU_list/out8-derep-filtered_table_clustered_ASV_with_tax.tsv
#rm -f ${OUT_DIR}export_tax_OTU_list/taxonomy.tsv
#rm -f ${OUT_DIR}export_tax_OTU_list/feature-table.biom

echo "file: ${OUT_DIR}export_tax_OTU_list/out8-derep-filtered_table_clustered_with_tax.tsv"


  ##   h) OTU table with all reads (singeltons, unassigned) for Shannon and Chao1
  mkdir -p ${OUT_DIR}export_for_OTU_table_with_singeltons
  qiime tools export --input-path ${OUT_DIR}out8-derep-filtered_table_clustered_ASV.qza \
  --output-path ${OUT_DIR}export_for_OTU_table_with_singeltons
  /bin/mv -f ${OUT_DIR}export_for_OTU_table_with_singeltons/feature-table.biom ${OUT_DIR}export_for_OTU_table_with_singeltons/out8-derep-filtered_table_clustered_ASV-final_table_with_singelton_feature-table.biom
  biom convert -i ${OUT_DIR}export_for_OTU_table_with_singeltons/out8-derep-filtered_table_clustered_ASV-final_table_with_singelton_feature-table.biom \
  -o ${OUT_DIR}export_for_OTU_table_with_singeltons/out8-derep-filtered_table_clustered_ASV-final_table_with_singelton_OTUs.tsv \
  --to-tsv


  
  #creating the barplot
  echo "creating barplot - PLOTTING"
  qiime taxa barplot \
  --i-table ${OUT_DIR}out9-derep-clustered_ASV-final_table_no_singelton.qza \
  --i-taxonomy ${OUT_DIR}out9-classify-vsearch_ASVOTUs_0.97SILVA.qza \
  --m-metadata-file ${OUT_DIR}$METADATA_FILE \
  --o-visualization ${OUT_DIR}out9_classify-table_barplot_ASV.qzv 
  
  #make some alphararefaction curves
  foreach SUBSAMPLING (1000 5000 3000 10000 30000 50000)
   echo "creating rare faction with $SUBSAMPLING subsampled reads of the ASV dataset..."
   qiime diversity alpha-rarefaction \
   --i-table ${OUT_DIR}out9-derep-clustered_ASV-final_table_no_singelton.qza\
   --p-max-depth {$SUBSAMPLING} \
   --m-metadata-file ${OUT_DIR}$METADATA_FILE\
   --o-visualization ${OUT_DIR}out9-derep-culstered_ASV-final_${SUBSAMPLING}.qzv \
   --p-steps 60
   
   echo "creating MDS type data with $SUBSAMPLING subsampled reads of the ASV dataset..."
   set OUT_FOLDER_NAME="core-metrics-phylogenetic-output_${SUBSAMPLING}_depth_ASV"
   /bin/rm -r -f $OUT_FOLDER_NAME
   qiime diversity core-metrics-phylogenetic \
   --i-table ${OUT_DIR}out9-derep-clustered_ASV-final_table_no_singelton.qza \
   --i-phylogeny ${OUT_DIR}out8-seqs-nonchimeras_clustered_ASV_rooted-tree.qza\
   --p-sampling-depth $SUBSAMPLING \
   --m-metadata-file ${OUT_DIR}$METADATA_FILE \
   --output-dir ${OUT_DIR}$OUT_FOLDER_NAME
   ####--p-n-jobs 16 \ #this does not work > "segmentation fault" error
   end



echo "pipeline has ended"
