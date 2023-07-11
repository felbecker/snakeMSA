data = "data"
DIRS, SAMPLES = glob_wildcards(data+"/{dir}/train/{sample}")
TOOLS = ["learnMSA", "famsa"]

ruleorder: learnMSA > msa

#compute one output table per tool
rule all:
    input:
        expand("{tool}.tbl", tool=TOOLS)
        
        
rule msa:
    input:  
         data+"/{dir}/train/{sample}"
    output:
        "{tool}/alignments/{dir}/{sample}"
    threads: 32
    resources:
        mem_mb = 256000,
        runtime = "3d"
    log:
        "{tool}/logs/{dir}/{sample}.log"
    benchmark:
        "{tool}/benchmarks/{dir}/{sample}.txt"
    params:
        tree = "{tool}/trees/{dir}/{sample}.mbed",
        tool = "{tool}"
    run:
        if params.tool == "famsa":
            shell("famsa -t {threads} {input} {output} > {log}")
        if params.tool == "t_coffee":
            #according to Santus et al. 2023 Supplementary information (their own tool)
            shell("""clustalo --threads={threads} -i {input} --guidetree-out params.tree ; \
            export MAX_N_PID_4_TCOFFEE=10000000 ; \
            t_coffee -reg -reg_method clustalo_msa -reg_tree params.tree -seq {input} \
            -reg_nseq 1000 -reg_thread {threads} -outfile {output} > {log}""")
        if params.tool == "kalign":
            shell("kalign -n {threads} -i {input} -o {output} > {log}")
            
         
rule learnMSA:
    input:  
         data+"/{dir}/train/{sample}"
    output:
        "learnMSA/alignments/{dir}/{sample}"
    threads: 32
    resources:
        mem_mb = 1000000, #temporary to make vision work
        nvidia_gpu = 1,
        runtime = "3d"
    log:
        "learnMSA/logs/{dir}/{sample}.log"
    benchmark:
        "learnMSA/benchmarks/{dir}/{sample}.txt"
    params:
        slice_dir="learnMSA/slices/{dir}/{sample}",
        cluster_dir="learnMSA/clustering/{dir}"
    shell:
        "python3 ../learnMSA/learnMSA.py -i {input} -o {output} -n 10 -d 0 --align_insertions --insertion_slice_dir {params.slice_dir} "
        "--sequence_weights --cluster_dir {params.cluster_dir} > {log}"
        
        
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
        mem_mb = 500000
    shell:
        "id_list=$(sed -n '/^>/p' {input.ref} | sed 's/^.//') ; "
        "export MAX_N_PID_4_TCOFFEE=10000000 ; "
        "t_coffee -other_pg seq_reformat -in {input.msa} -action +extract_seq_list ${{id_list[@]}} +rm_gap"
        "> {output}"
        
       
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
        mem_mb = 16000
    shell:
        "nseq=$(grep -c '>' {input.all}) ; "
        "avg_len=$(awk '{{/>/&&++a||b+=length()}}END{{print b/a}}' {input.all}) ;"
        "avg_ref_len=$(awk '{{/>/&&++a||b+=length()}}END{{print b/a}}' {input.ref}) ;"
        "export MAX_N_PID_4_TCOFFEE=10000000 ; "
        "sp_out=$(t_coffee -other_pg aln_compare -al1 {input.ref} -al2 {input.msa} -compare_mode sp | "
        "grep -v 'seq1' | grep -v '*') ; "
        "nseq_ref=$(echo $sp_out | awk '{{ print $2}}') ; "
        "sim=$(echo $sp_out | awk '{{ print $3}}') ; "
        "sp=$(echo $sp_out | awk '{{ print $4}}') ; "
        "modeler=$(t_coffee -other_pg aln_compare -al1 {input.msa}  -al2 {input.ref} -compare_mode sp | "
        "grep -v 'seq1' | grep -v '*' | awk '{{ print $4}}') ; "
        "tc=$(t_coffee -other_pg aln_compare -al1 {input.ref} -al2 {input.msa} -compare_mode tc | "
        "grep -v 'seq1' | grep -v '*' | awk '{{ print $4}}') ; "
        "col=$(t_coffee -other_pg aln_compare -al1 {input.ref} -al2 {input.msa} -compare_mode column | "
        "grep -v 'seq1' | grep -v '*' | awk '{{ print $4}}') ; "
        "echo {wildcards.dir} {wildcards.sample} $nseq $nseq_ref $avg_len $avg_ref_len $sim $sp $modeler $tc $col >> {output}"
        
        
rule concat_scores:
    input:
        scores = expand("{{tool}}/scores/{dir}/{sample}", zip, dir=DIRS, sample=SAMPLES),
        benchmarks = expand("{{tool}}/benchmarks/{dir}/{sample}.txt", zip, dir=DIRS, sample=SAMPLES)
    output:
        "{tool}.tbl"
    threads: 1
    resources:
        mem_mb = 1000
    run:
        shell("echo -n "" > {output}")
        for s,b in zip(input.scores, input.benchmarks):
            #row by row for each sample, concatenate the scores and the second line of the benchmark file
            shell("echo $(cat {s}) $(tail -n 1 {b}) >> {output}")
