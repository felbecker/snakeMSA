# type: ignore

configfile: "configs/default.json"

TOOLS = config["tools"]
SAMPLES, = glob_wildcards(config["data_path"]+"/unaligned/{sample}")

rule all:
    input:
        expand("results_"+config["name"]+"/{tool}.out", 
               tool=TOOLS,
               sample=SAMPLES)


include: "rules/align"
include: "rules/score_msa"
