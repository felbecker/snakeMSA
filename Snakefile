data = "data"
DIRS, SAMPLES = glob_wildcards(data+"/{dir}/train/{sample}")
#SAMPLES, = glob_wildcards(data+"/homfam/train/{sample}")
#DIRS = ["homfam"]*len(SAMPLES)
TOOLS = ["learnMSA_paper", "learnMSA", "learnMSA_language", "famsa", "t_coffee", "clustalo", "muscle", "mafft"]

ruleorder: learnMSA > learnMSA_language > learnMSA_paper > msa

#compute one output table per tool
rule all:
    input:
        expand("{tool}.tbl", tool=TOOLS)
        
        
rule msa:
    input:  
         data+"/{dir}/train/{sample}"
    output:
        "{tool}/alignments/{dir}/{sample}",
    threads: 32
    resources:
        mem_mb = 256000,
        runtime = "3d",
        msa_load = 1
    log:
        "{tool}/logs/{dir}/{sample}.log"
    benchmark:
        "{tool}/benchmarks/{dir}/{sample}.txt"
    params:
        tree = "{tool}/trees/{dir}/{sample}.mbed",
        tree_dir = "{tool}/trees/{dir}/",
        magus_tmp = "./magus/tmp/{sample}",
        tool = "{tool}"
    run:
        if params.tool == "famsa":
            cmd = f"famsa -t {threads} {input} {output} > {log}"
        if params.tool == "t_coffee":
            cmd = f"""mkdir -p {params.tree_dir} ; \
            clustalo --threads={threads} -i {input} --guidetree-out {params.tree} --force -o /dev/null ; \
            export MAX_N_PID_4_TCOFFEE=10000000 ; \
            t_coffee -reg -reg_method famsa_msa -reg_tree {params.tree} -seq {input} \
            -reg_nseq 1000 -reg_thread {threads} -outfile {output} > {log}"""
        if params.tool == "kalign":
            cmd = f"kalign -n {threads} -i {input} -o {output} > {log}"
        if params.tool == "clustalo":
            cmd = f"clustalo --threads={threads} -i {input} -o {output} > {log}"
        if params.tool == "magus":
            cmd = f"magus -d {params.magus_tmp} --numprocs {threads} -i {input} -o {output} > {log}"
        if params.tool == "muscle":
            cmd = f"muscle -super5 {input} -output {output} > {log}"
        if params.tool == "mafft":
            cmd = f"mafft --thread {threads} --anysymbol --quiet --dpparttree {input} > {output}"

        p = subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
        stdout, stderr= p.communicate()

        if p.returncode != 0:
            with open(output[0]) as fout:
                fout.write("FAILED")
            print(stderr.decode(), file=sys.stderr)
            
            
         
rule learnMSA:
    input:  
         data+"/{dir}/train/{sample}"
    output:
        "learnMSA/alignments/{dir}/{sample}"
    threads: 32
    resources:
        mem_mb = 240000, 
        gpu = 1,
        runtime = "3d",
        learnMSA_load = 1,
        partition = "vision"
    log:
        "learnMSA/logs/{dir}/{sample}.log"
    benchmark:
        "learnMSA/benchmarks/{dir}/{sample}.txt"
    params:
        slice_dir="learnMSA/slices/{dir}/{sample}",
        cluster_dir="learnMSA/clustering/{dir}"
    shell:
        "python3 ../learnMSA/learnMSA.py -i {input} -o {output} -n 10 "
        "--sequence_weights --cluster_dir {params.cluster_dir} > {log}"
        
        

rule learnMSA_paper:
    input:  
         data+"/{dir}/train/{sample}"
    output:
        "learnMSA_paper/alignments/{dir}/{sample}"
    threads: 32
    resources:
        mem_mb = 240000, 
        gpu = 1,
        runtime = "3d",
        learnMSA_load = 1,
        partition = "vision"
    log:
        "learnMSA_paper/logs/{dir}/{sample}.log"
    benchmark:
        "learnMSA_paper/benchmarks/{dir}/{sample}.txt"
    params:
        slice_dir="learnMSA_paper/slices/{dir}/{sample}",
        cluster_dir="learnMSA_paper/clustering/{dir}"
    shell:
        "python3 ../learnMSA/learnMSA.py -i {input} -o {output} --unaligned_insertions -n 10 > {log}"
        


rule learnMSA_language:
    input:  
         data+"/{dir}/train/{sample}"
    output:
        "learnMSA_language/alignments/{dir}/{sample}"
    threads: 32
    resources:
        mem_mb = 240000, 
        gpu = 1,
        runtime = "3d",
        learnMSA_load = 1,
        partition = "vision"
    log:
        "learnMSA_language/logs/{dir}/{sample}.log"
    benchmark:
        "learnMSA_language/benchmarks/{dir}/{sample}.txt"
    params:
        slice_dir="learnMSA_language/slices/{dir}/{sample}",
        cluster_dir="learnMSA_language/clustering/{dir}"
    shell:
        "python3 ../learnMSA/learnMSA.py -i {input} -o {output} -n 10 "
        "--sequence_weights --cluster_dir {params.cluster_dir} --use_language_model > {log}"
        
        
        
## Derive the sub-MSA of the reference sequences from the full MSA
rule project_references:
    input:
        msa = "{tool}/alignments/{dir}/{sample}",
        ref = data+"/{dir}/refs/{sample}"
    output:
        "{tool}/projections/{dir}/{sample}"
    threads: 32
    resources:
        #some aligners might produce extremely large fasta files
        #this allocates a huge chunk of memory, but the jobs are short
        mem_mb = 120000,
        partition = "batch"
    run:
        exit_status = open(input.msa).readline()[0].strip()
        if exit_status == "FAILED":
            # pass failed samples through
             with open(output[0]) as fout:
                fout.write("FAILED")
        else:
            shell("""id_list=$(sed -n '/^>/p' {input.ref} | sed 's/^.//') ; \
                    export MAX_N_PID_4_TCOFFEE=10000000 ; \
                    t_coffee -other_pg seq_reformat -in {input.msa} -action +extract_seq_list ${{id_list[@]}} +rm_gap \
                    > {output}""")
        
       
## Compute a number of some general statistics and scoring metrics and store them in per-sample files
rule score_msa:
    input:
        msa = "{tool}/projections/{dir}/{sample}",
        ref = data+"/{dir}/refs/{sample}",
        all = data+"/{dir}/train/{sample}"
    output:
        "{tool}/scores/{dir}/{sample}"
    threads: 8
    resources:
        mem_mb = 16000,
        partition = "batch"
    run:
        exit_status = open(input.msa).readline()[0].strip()
        if exit_status == "FAILED":
            # insert a score of 0 for failed files
            shell("""nseq=$(grep -c '>' {input.all}) ; \
                    avg_len=$(awk '{{/>/&&++a||b+=length()}}END{{print b/a}}' {input.all}) ;\
                    avg_ref_len=$(awk '{{/>/&&++a||b+=length()}}END{{print b/a}}' {input.ref}) ;\
                    echo {wildcards.dir} {wildcards.sample} $nseq $nseq_ref $avg_len $avg_ref_len $sim 0 0 0 0 >> {output}""")
        else:
            shell("""nseq=$(grep -c '>' {input.all}) ; \
                    avg_len=$(awk '{{/>/&&++a||b+=length()}}END{{print b/a}}' {input.all}) ;\
                    avg_ref_len=$(awk '{{/>/&&++a||b+=length()}}END{{print b/a}}' {input.ref}) ;\
                    export MAX_N_PID_4_TCOFFEE=10000000 ; \
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
        scores = expand("{{tool}}/scores/{dir}/{sample}", zip, dir=DIRS, sample=SAMPLES),
        benchmarks = expand("{{tool}}/benchmarks/{dir}/{sample}.txt", zip, dir=DIRS, sample=SAMPLES)
    output:
        "{tool}.tbl"
    threads: 1
    resources:
        mem_mb = 1000,
        partition = "batch"
    run:
        shell("echo -n "" > {output}")
        for s,b in zip(input.scores, input.benchmarks):
            #row by row for each sample, concatenate the scores and the second line of the benchmark file
            shell("echo $(cat {s}) $(tail -n 1 {b}) >> {output}")
