python tools/summarize.py -i famsa_homfam learnMSA_homfam foldmason_homstradfam \
--all --average-replicas --barplots --barplot_title "HomFam benchmark" \
--barplot_fontsize 12 --barplot_tool_names tools/tool_names.txt


python tools/summarize.py -i famsa_homfam learnMSA_homfam foldmason_homstradfam \
--all --average-replicas --scatterplot_lddt lddt_muscle lddt \
--barplot_title "HomFam LDDT (Muscle vs FoldMason)" \
--barplot_fontsize 12 --barplot_tool_names tools/tool_names.txt


python tools/summarize.py -i famsa_homfam learnMSA_homfam foldmason_homstradfam \
--all --average-replicas --scatterplot_lddt lddt_core dali_z \
--barplot_title "HomFam LDDT (Core) vs Dali Z score" \
--barplot_fontsize 12 --barplot_tool_names tools/tool_names.txt
