data = "data"
#DIRS, SAMPLES = glob_wildcards(data+"/{dir}/train/{sample}")
TOOLS = ["learnMSA_no_weights"] #["learnMSA", "learnMSA_language2", "learnMSA_paper", "famsa", "t_coffee", "muscle"]
#to run all tools only of HomFam, uncomment the next 3 lines
SAMPLES, = glob_wildcards(data+"/homfam/train/{sample}")
DIRS = ["homfam"]*len(SAMPLES)
#TOOLS = ["learnMSA", "learnMSA_language2", "learnMSA_paper", "mafft_sparsecore", "famsa", "t_coffee", "clustalo", "muscle", "magus"]
ruleorder: learnMSA_no_weights > learnMSA > learnMSA_language > learnMSA_language2 > learnMSA_paper > learnMSA_language_cpu > msa

#compute one output table per tool
rule all:
    input:
        expand("results/{tool}.out", tool=TOOLS)
        
        
rule msa:
    input:  
         data+"/{dir}/train/{sample}"
    output:
        "outputs/{tool}/alignments/{dir}/{sample}"
    threads: 32
    resources:
        mem_mb = 250000,
        runtime = "3d",
        msa_load = 1,
        partition = "pinky,vision"
    log:
        "outputs/{tool}/logs/{dir}/{sample}.log"
    benchmark:
        "outputs/{tool}/benchmarks/{dir}/{sample}.txt"
    params:
        tree = "outputs/{tool}/trees/{dir}/{sample}.mbed",
        tree_dir = "outputs/{tool}/trees/{dir}/",
        magus_tmp = "outputs/magus/tmp/{dir}/{sample}",
        tool = "{tool}",
        dir = "{dir}"
    run:
        import subprocess
        if params.tool == "famsa":
            cmd = f"famsa -t {threads} {input} {output} > {log}"
        if params.tool == "t_coffee":
            cmd = f"""mkdir -p {params.tree_dir} ; \
            clustalo --threads={threads} -i {input} --guidetree-out {params.tree} --force -o /dev/null ; \
            export MAX_N_PID_4_TCOFFEE=$(cat /proc/sys/kernel/pid_max) ; \
            t_coffee -reg -reg_method famsa_msa -reg_tree {params.tree} -seq {input} \
            -reg_nseq 1000 -reg_thread {threads} -outfile {output} > {log}"""
        if params.tool == "t_coffee_sparsecore": 
            #unable to run this, may try again in the future
            cmd = f"""mkdir -p {params.tree_dir} ; \
            clustalo --threads={threads} -i {input} --guidetree-out {params.tree} --force -o /dev/null ; \
            export MAX_N_PID_4_TCOFFEE=$(cat /proc/sys/kernel/pid_max) ; \
            t_coffee -reg -reg_method mafftsparsecore_msa -reg_tree {params.tree} -seq {input} \
            -reg_nseq 1000 -reg_thread {threads} -outfile {output} > {log}"""
        if params.tool == "kalign":
            cmd = f"~/mambaforge/envs/snakeMSA/bin/kalign -n {threads} -i {input} -o {output} > {log}"
        if params.tool == "clustalo":
            cmd = f"clustalo --threads={threads} -i {input} -o {output} > {log}"
        if params.tool == "magus":
            cmd = f"magus -d {params.magus_tmp} --numprocs {threads} -i {input} -o {output} > {log}"
        if params.tool == "muscle":
            if params.dir == "homstrad":
                cmd = f"muscle -align {input} -output {output} > {log}"
            else:
                cmd = f"muscle -super5 {input} -output {output} > {log}"
        if params.tool == "mafft":
            if "ext" in params.dir:
                cmd = f"mafft --thread {threads} --anysymbol --quiet --dpparttree {input} > {output}"
            else:
                cmd = f"mafft --thread {threads} --anysymbol --quiet {input} > {output}"
        if params.tool == "mafft_sparsecore":
            cmd = f"mkdir -p tmp/{input} && python3 util/replace_non_standard_aa.py {input} > tmp/{input}/data && mafft-sparsecore.rb -i tmp/{input}/data > {output} && rm -r tmp/{input}"

        p = subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
        stdout, stderr= p.communicate()
        
        print(stdout.decode(), file=sys.stdout)

        if not os.path.exists(output[0]):
            with open(output[0], "w") as fout:
                fout.write("FAILED")
            print(stderr.decode(), file=sys.stderr)
            
            
         
rule learnMSA:
    input:  
         data+"/{dir}/train/{sample}"
    output:
        "outputs/learnMSA/alignments/{dir}/{sample}"
    threads: 4
    resources:
        mem_mb = 250000, 
        gpu = 1,
        runtime = "3d",
        learnMSA_load = 1,
        partition = "vision"
    log:
        "outputs/learnMSA/logs/{dir}/{sample}.log"
    benchmark:
        "outputs/learnMSA/benchmarks/{dir}/{sample}.txt"
    params:
        cluster_dir="outputs/learnMSA/clustering/{dir}"
    shell:
        "python3 ../tmp_work/learnMSA/learnMSA.py -i {input} -o {output} -n 10 "
        "--sequence_weights --cluster_dir {params.cluster_dir} > {log}"


rule learnMSA_no_weights:
    input:  
         data+"/{dir}/train/{sample}"
    output:
        "outputs/learnMSA_no_weights/alignments/{dir}/{sample}"
    threads: 4
    resources:
        mem_mb = 250000, 
        gpu = 1,
        runtime = "1h",
        learnMSA_load = 1,
        partition = "vision"
    log:
        "outputs/learnMSA_no_weights/logs/{dir}/{sample}.log"
    benchmark:
        "outputs/learnMSA_no_weights/benchmarks/{dir}/{sample}.txt"
    shell:
        "python3 ../tmp_work/learnMSA/learnMSA.py -i {input} -o {output} -n 10 > {log}"
        
        

rule learnMSA_paper:
    input:  
         data+"/{dir}/train/{sample}"
    output:
        "outputs/learnMSA_paper/alignments/{dir}/{sample}"
    threads: 4
    resources:
        mem_mb = 250000, 
        gpu = 1,
        runtime = "3d",
        learnMSA_load = 1,
        partition = "vision"
    log:
        "outputs/learnMSA_paper/logs/{dir}/{sample}.log"
    benchmark:
        "outputs/learnMSA_paper/benchmarks/{dir}/{sample}.txt"
    params:
        cluster_dir="outputs/learnMSA_paper/clustering/{dir}"
    shell:
        "python3 ../tmp_work/learnMSA/learnMSA.py -i {input} -o {output} --unaligned_insertions -n 10 > {log}"
        


rule learnMSA_language:
    input:  
         data+"/{dir}/train/{sample}"
    output:
        "outputs/learnMSA_language/alignments/{dir}/{sample}"
    threads: 4
    resources:
        mem_mb = 250000, 
        gpu = 1,
        runtime = "3d",
        learnMSA_load = 1,
        partition = "vision"
    log:
        "outputs/learnMSA_language/logs/{dir}/{sample}.log"
    benchmark:
        "outputs/learnMSA_language/benchmarks/{dir}/{sample}.txt"
    params:
        cluster_dir="outputs/learnMSA_language/clustering/{dir}"
    shell:
        "python3 ../tmp_work/learnMSA/learnMSA.py -i {input} -o {output} -n 4 \
        --sequence_weights --cluster_dir {params.cluster_dir} --use_language_model --language_model protT5 \
        --scoring_model_activation softmax --scoring_model_dim 32 \
        --embedding_prior_components 10 --frozen_insertions > {log}"



rule learnMSA_language2:
    input:  
         data+"/{dir}/train/{sample}"
    output:
        "outputs/learnMSA_language2/alignments/{dir}/{sample}"
    threads: 32
    resources:
        mem_mb = 250000, 
        gpu = 1,
        runtime = "3d",
        learnMSA_load = 1,
        partition = "vision"
    log:
        "outputs/learnMSA_language2/logs/{dir}/{sample}.log"
    benchmark:
        "outputs/learnMSA_language2/benchmarks/{dir}/{sample}.txt"
    params:
        cluster_dir="outputs/learnMSA_language2/clustering/{dir}"
    shell:
        "python3 ../tmp_work/learnMSA/learnMSA.py -i {input} -o {output} -n 5 \
        --sequence_weights --cluster_dir {params.cluster_dir} --use_language_model --language_model protT5 \
        --scoring_model_activation sigmoid --scoring_model_dim 16 \
        --embedding_prior_components 32 --frozen_insertions > {log}"



rule learnMSA_language_cpu:
    input:  
         data+"/{dir}/train/{sample}"
    output:
        "outputs/learnMSA_language_cpu/alignments/{dir}/{sample}"
    threads: 32
    resources:
        mem_mb = 250000,
        runtime = "3d",
        msa_load = 1,
        partition = "pinky"
    log:
        "outputs/learnMSA_language_cpu/logs/{dir}/{sample}.log"
    benchmark:
        "outputs/learnMSA_language_cpu/benchmarks/{dir}/{sample}.txt"
    params:
        cluster_dir="outputs/learnMSA_language_cpu/clustering/{dir}"
    shell:
        "python3 ../tmp_work/learnMSA/learnMSA.py -i {input} -o {output} -n 5 \
        --sequence_weights --cluster_dir {params.cluster_dir} --use_language_model --language_model protT5 \
        --scoring_model_activation sigmoid --scoring_model_dim 16 \
        --embedding_prior_components 32 --frozen_insertions > {log}"
        


## Derive the sub-MSA of the reference sequences from the full MSA
rule project_references:
    input:
        msa = "outputs/{tool}/alignments/{dir}/{sample}",
        ref = data+"/{dir}/refs/{sample}"
    output:
        "outputs/{tool}/projections/{dir}/{sample}"
    threads: 32
    resources:
        #some aligners might produce extremely large fasta files
        #this allocates a huge chunk of memory, but the jobs are short
        mem_mb = 500000,
        partition = "pinky,vision",
        runtime = "1h"
    run:
        import os
        if os.path.getsize(input.msa) > int(2.5e11): #output MSA too large (>250GB) too handle, assume failed
             with open(output[0], "w") as fout:
                fout.write("FAILED")
        else:
            with open(input.msa, "r") as msa_file:
                exit_status = msa_file.readline().strip()
                if exit_status == "FAILED":
                    # pass failed samples through
                     with open(output[0], "w") as fout:
                        fout.write("FAILED")
                else:
                    shell("""id_list=$(sed -n '/^>/p' {input.ref} | sed 's/^.//') ; \
                            export MAX_N_PID_4_TCOFFEE=$(cat /proc/sys/kernel/pid_max) ; \
                            export MAX_N_PID_4_TCOFFEE=$(cat /proc/sys/kernel/pid_max) ; \
                            t_coffee -other_pg seq_reformat -in {input.msa} -action +extract_seq_list ${{id_list[@]}} +rm_gap \
                            > {output}""")
        
       
## Compute a number of some general statistics and scoring metrics and store them in per-sample files
rule score_msa:
    input:
        msa = "outputs/{tool}/projections/{dir}/{sample}",
        ref = data+"/{dir}/refs/{sample}",
        all = data+"/{dir}/train/{sample}"
    output:
        "outputs/{tool}/scores/{dir}/{sample}"
    threads: 8
    resources:
        mem_mb = 16000,
        partition = "batch",
        runtime = "1h"
    run:
        with open(input.msa, "r") as msa_file:
            exit_status = msa_file.readline().strip()
            if exit_status == "FAILED":
                # insert a score of 0 for failed files
                shell("""nseq=$(grep -c '>' {input.all}) ; \
                        avg_len=$(awk '{{/>/&&++a||b+=length()}}END{{print b/a}}' {input.all}) ;\
                        avg_ref_len=$(awk '{{/>/&&++a||b+=length()}}END{{print b/a}}' {input.ref}) ;\
                        echo {wildcards.dir} {wildcards.sample} $nseq 0 $avg_len $avg_ref_len 0 0 0 0 0 >> {output}""")
            else:
                shell("""nseq=$(grep -c '>' {input.all}) ; \
                        avg_len=$(awk '{{/>/&&++a||b+=length()}}END{{print b/a}}' {input.all}) ;\
                        avg_ref_len=$(awk '{{/>/&&++a||b+=length()}}END{{print b/a}}' {input.ref}) ;\
                        export MAX_N_PID_4_TCOFFEE=$(cat /proc/sys/kernel/pid_max) ; \
                        sp_out=$(t_coffee -other_pg aln_compare -al1 {input.ref} -al2 {input.msa} -compare_mode sp | \
                        grep -v 'seq1' | grep -v '*') ; \
                        nseq_ref=$(echo $sp_out | awk '{{ print $2}}') ; \
                        sim=$(echo $sp_out | awk '{{ print $3}}') ; \
                        sp=$(echo $sp_out | awk '{{ print $4}}') ; \
                        modeler=$(t_coffee -other_pg aln_compare -al1 {input.msa}  -al2 {input.ref} -compare_mode sp | \
                        grep -v 'seq1' | grep -v '*' | awk '{{ print $4}}') ; \
                        tc=$(t_coffee -other_pg aln_compare -al1 {input.ref} -al2 {input.msa} -compare_mode tc | \
                        grep -v 'seq1' | grep -v '*' | awk '{{ print $4}}') ; \
                        col=$(t_coffee -other_pg aln_compare -al1 {input.ref} -al2 {input.msa} -compare_mode column | \
                        grep -v 'seq1' | grep -v '*' | awk '{{ print $4}}') ; \
                        echo {wildcards.dir} {wildcards.sample} $nseq $nseq_ref $avg_len $avg_ref_len $sim $sp $modeler $tc $col >> {output}""")
        
        
rule concat_scores:
    input:
        scores = expand("outputs/{{tool}}/scores/{dir}/{sample}", zip, dir=DIRS, sample=SAMPLES),
        benchmarks = expand("outputs/{{tool}}/benchmarks/{dir}/{sample}.txt", zip, dir=DIRS, sample=SAMPLES)
    output:
        "results/{tool}.out"
    threads: 1
    resources:
        mem_mb = 1000,
        partition = "batch",
        runtime = "1h"
    run:
        shell("echo -n "" > {output}_tmp")
        for s,b in zip(input.scores, input.benchmarks):
            #row by row for each sample, concatenate the scores and the second line of the benchmark file
            shell("echo $(cat {s}) $(tail -n 1 {b}) >> {output}_tmp")
        shell("python3 util/fix_timestamps.py {output}_tmp > {output}")
        shell("rm {output}_tmp")
