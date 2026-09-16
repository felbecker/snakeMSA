# type: ignore

TOOLS = config["tools"]


def data_file(*parts):
    """
    Resolves relative/absolute path. Relative paths do not start with `/`and
    are preceded by the `data_path` config variable.
    """
    return os.path.join(
        os.path.expanduser(config["data_path"]),
        *(os.path.expanduser(part) for part in parts)
    )


def resolve_pdb_path(wildcards):
    if "pdb" not in config:
        return ""
    return data_file(config["pdb"], os.path.splitext(wildcards.sample)[0])


SAMPLES, = glob_wildcards(data_file("unaligned", "{sample}"))

rule all:
    input:
        expand("results/"+config["name"]+"/{tool}.out",
               tool=TOOLS,
               sample=SAMPLES)


include: "rules/align"
include: "rules/score_msa"
include: "rules/evaluate"
