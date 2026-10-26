# HomFam
python tools/summarize.py -i learnMSA_pub_homfam competition_pub_homfam \
--average-replicas --all --all-struct-scores \
--tools learnMSA learnMSA_struct clustalo famsa mafft muscle_super5 t-coffee kalign \
--barplots --barplot_title "HomFam Benchmark" --barplot_fontsize 12 \
--barplot_tool_names tools/tool_names.txt \

# LearnMSA embeddings vs 3Di
python tools/summarize.py -i learnMSA_pub_homfam \
--average-replicas --all --all-struct-scores \
--barplots --barplot_title "Embeddings vs 3Di (HomFam)" --barplot_fontsize 12 \
--barplot_tool_names tools/tool_names.txt \

# Homfam 3Di
python tools/summarize.py -i learnMSA_pub_homfam competition_pub_homfam \
--average-replicas --all --all-struct-scores \
--tools learnMSA learnMSA_struct mafft mafft_3Di \
--barplots --barplot_title "HomFam 3Di Benchmark" --barplot_fontsize 12 \
--barplot_tool_names tools/tool_names.txt \
--barplot_output plots/homfam_3di.pdf

# AFDB Clusters
python tools/summarize.py -i learnMSA_pub_afdb_clusters competition_pub_afdb_clusters \
--average-replicas --all --all-struct-scores \
--barplot_lddt --barplot_title "AFDB Clusters LDDT" --barplot_fontsize 12 \
--barplot_tool_names tools/tool_names_pdb.txt

# AFDB Max Cluster
python tools/summarize.py -i learnMSA_pub_afdb_max_cluster competition_pub_afdb_max_cluster \
--average-replicas --all --all-struct-scores \
--barplot_lddt --barplot_title "AFDB Max Cluster (35k) LDDT" --barplot_fontsize 12 \
--barplot_tool_names tools/tool_names_pdb.txt
