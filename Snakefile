# type: ignore

TOOLS = config["tools"]


def config_flag(key, default=False):
    """Boolean config value; strings from `--config key=false` are parsed."""
    value = config.get(key, default)
    if isinstance(value, str):
        return value.strip().lower() in ("true", "yes", "1")
    return bool(value)


# Fuse alignment and scoring into one job per sample and delete the MSA,
# projection and work dir once the scores are written (rules/align_and_score).
CLEANUP = config_flag("cleanup")


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
include: "rules/align_and_score"
include: "rules/evaluate"
