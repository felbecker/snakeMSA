python tools/summarize.py -i famsa_homfam learnMSA_homfam foldmason_homstradfam \
--all --average-replicas --barplots --barplot_title "HomFam benchmark" \
--barplot_fontsize 12 --barplot_tool_names tools/tool_names.txt


python tools/summarize.py -i famsa_homfam learnMSA_homfam foldmason_homstradfam \
--all --average-replicas --barplot_lddt lddt --barplot_title "HomFam LDDT (FoldMason)" \
--barplot_fontsize 12 --barplot_tool_names tools/tool_names.txt


python tools/summarize.py -i famsa_homfam learnMSA_homfam foldmason_homstradfam \
--all --average-replicas --barplot_lddt lddt_core --barplot_title "HomFam LDDT (Core)" \
--barplot_fontsize 12 --barplot_tool_names tools/tool_names.txt


python tools/summarize.py -i famsa_homfam learnMSA_homfam foldmason_homstradfam \
--all --average-replicas --barplot_lddt lddt_muscle --barplot_title "HomFam LDDT (Muscle)" \
--barplot_fontsize 12 --barplot_tool_names tools/tool_names.txt


python tools/summarize.py -i famsa_homfam learnMSA_homfam foldmason_homstradfam \
--all --average-replicas --barplot_lddt dali_z --barplot_title "HomFam Dali Z score" \
--barplot_fontsize 12 --barplot_tool_names tools/tool_names.txt
