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
--tools learnMSA learnMSA_struct famsa foldmason muscle_3d muscle_3d_super7 \
--scatterplot_lddt --barplot_title "AFDB Clusters LDDT" --barplot_fontsize 12 \
--barplot_tool_names tools/tool_names_pdb.txt

# AFDB Max Cluster
python tools/summarize.py -i learnMSA_pub_afdb_max_cluster competition_pub_afdb_max_cluster \
--average-replicas --all --all-struct-scores \
--scatterplot_lddt --barplot_title "AFDB Max Cluster (35k) LDDT" --barplot_fontsize 12 \
--barplot_tool_names tools/tool_names_pdb.txt

# Homologous embeddings
python tools/summarize.py -i reuse_learnMSA_struct learnMSA_pub_homfam \
--average-replicas --all --fixed_order --all \
--tools learnMSA_struct_reuse050 learnMSA_struct_reuse070 learnMSA_struct_reuse080 \
learnMSA_struct_reuse090 learnMSA_struct_reuse095 learnMSA_struct \
--runtime_line \
--scatterplot_lddt --barplot_title "Homologs share embeddings" --barplot_fontsize 12 \
--barplot_tool_names tools/tool_names_homolog_emb.txt

# learnMSA size scaling (alignment depth)
for score in SP-Score TC lddt_muscle; do
python tools/summarize_size_scaling.py --run learnMSA_scale \
--tools learnMSA learnMSA_struct --score $score \
--plot --title "learnMSA Size Scaling" --fontsize 12 \
--tool_names tools/tool_names.txt \
--output plots/learnMSA_scale_$score.pdf
done
