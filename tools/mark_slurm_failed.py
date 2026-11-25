"""Mark samples whose SLURM job ran out of memory or time as failed.

The align rule only catches a failing aligner process. When SLURM kills the
whole job (cgroup OOM killer, time limit), nothing writes the FAILED marker,
snakemake removes the outputs and the workflow fails. This script finds such
jobs in the snakemake logs and writes the outputs a caught failure would have
produced, so that the summary lists the samples as failed (success 0) instead
of the workflow failing again on the next run.

Run from the snakeMSA root. Samples that already have their output (e.g. a
later run succeeded) are left alone.

Example:
    python tools/mark_slurm_failed.py --dry
    python tools/mark_slurm_failed.py --all-logs
    python tools/mark_slurm_failed.py --all-logs -i competition_pub_*
"""

import argparse
import fnmatch
import glob
import os
import re
from typing import Dict, List, NamedTuple, Optional, Tuple


LOG_DIR = ".snakemake/log"

#: Columns of a snakemake benchmark file. A killed job has no record, so
#: every value is "NA", like snakemake writes for missing measurements.
BENCHMARK_COLS = [
    "s", "h:m:s", "max_rss", "max_vms", "max_uss", "max_pss",
    "io_in", "io_out", "mean_load", "cpu_time",
]

#: Number of score columns after "<sample> <success>" (7 FastSP scores,
#: lddt_muscle, lddt_muscle_L, lddt_muscle_conf, lddt_muscle_L_conf, coverage,
#: lddt_fm), see STRUCT_COLUMNS in rules/score_msa.
NUM_SCORES = 13

ERROR_RE = re.compile(r"^Error in rule (\S+):$")
SLURM_RE = re.compile(
    r"^\s+message: SLURM-job '(\d+)' failed, SLURM status is: '(\w+)'"
)
OUTPUT_RE = re.compile(r"^\s+output: (.*)$")
# outputs/<run>/<tool>/scores/<sample>.score or .../alignments/<sample>
PATH_RE = re.compile(r"^outputs/([^/]+)/([^/]+)/(scores|alignments)/(.+)$")


argparser = argparse.ArgumentParser(
    description="Report SLURM jobs that ran out of memory or time and mark "
                "their samples as failed."
)
argparser.add_argument(
    'logs', nargs='*', type=str,
    help="Snakemake log files to scan. Default: the most recent one in "
         f"{LOG_DIR}."
)
argparser.add_argument(
    '--all-logs', action='store_true', dest='all_logs',
    help=f"Scan every log in {LOG_DIR}. A sample whose outputs were deleted "
         "on purpose (summarize.py --delete_failed) is marked again if an "
         "older log has a failure for it."
)
argparser.add_argument(
    '-i', '--runs', nargs='+', type=str, metavar='RUN',
    help="Only mark samples of these runs (directories in outputs/, i.e. "
         "the config's name; supports wildcards with *). Default: all runs."
)
argparser.add_argument(
    '--states', nargs='+', type=str, default=["OUT_OF_MEMORY", "TIMEOUT"],
    help="SLURM states to treat as a failed sample. Default: OUT_OF_MEMORY "
         "TIMEOUT."
)
argparser.add_argument(
    '--dry', action='store_true',
    help="Only report the jobs and the files that would be written."
)


class KilledJob(NamedTuple):
    log: str
    rule: str
    slurm_id: str
    state: str
    output: str


def find_logs() -> List[str]:
    if args.logs:
        return args.logs
    # the timestamped names sort chronologically
    logs = sorted(glob.glob(os.path.join(LOG_DIR, "*.snakemake.log")))
    if not logs:
        raise SystemExit(f"No snakemake logs found in {LOG_DIR}.")
    return logs if args.all_logs else logs[-1:]


def parse_log(path: str) -> List[KilledJob]:
    """The failed SLURM jobs of a log with the first output of their rule.

    Snakemake reports each failure twice, once with the SLURM status and once
    with "message: None" in the final summary; only the former is used.
    """
    jobs = []
    rule = state = slurm_id = None
    with open(path) as log_file:
        for line in log_file:
            line = line.rstrip("\n")
            match = ERROR_RE.match(line)
            if match:
                rule, state = match.group(1), None
                continue
            if rule is None:
                continue
            match = SLURM_RE.match(line)
            if match:
                slurm_id, state = match.groups()
                continue
            match = OUTPUT_RE.match(line)
            if match:
                if state is not None:
                    output = match.group(1).split(", ")[0]
                    jobs.append(KilledJob(path, rule, slurm_id, state, output))
                rule = None
    return jobs


def failed_score_line(sample: str) -> str:
    # same sample name and layout as score_sample writes for a FAILED MSA
    sample_name = os.path.basename(sample).rsplit('.', 1)[0]
    return f"{sample_name} 0 " + " ".join(["0"] * NUM_SCORES)


def in_runs(job: KilledJob) -> bool:
    """Whether the job belongs to one of the runs given with --runs."""
    if not args.runs:
        return True
    match = PATH_RE.match(job.output)
    return match is not None and any(
        fnmatch.fnmatchcase(match.group(1), pattern) for pattern in args.runs
    )


def sample_of(job: KilledJob) -> Optional[Tuple[str, str]]:
    """(tool dir, sample) of a job of a handled rule, else None."""
    match = PATH_RE.match(job.output)
    if match is None:
        return None
    run_name, tool, kind, sample = match.groups()
    tool_dir = os.path.join("outputs", run_name, tool)
    if job.rule == "align_and_score" and kind == "scores":
        return tool_dir, sample[:-len(".score")]
    if job.rule == "align" and kind == "alignments":
        return tool_dir, sample
    return None


def files_to_write(job: KilledJob, tool_dir: str, sample: str
                   ) -> Dict[str, str]:
    """Path -> content of the files that mark the job's sample as failed."""
    if job.rule == "align_and_score":
        primary = failed_score_line(sample)
    else:
        # score_msa turns the marker into a failed score file
        primary = "FAILED"
    bench_stem = os.path.join(tool_dir, "benchmarks", sample)
    na_row = ["NA"] * len(BENCHMARK_COLS)
    return {
        job.output: primary,
        bench_stem + ".txt": "\t".join(BENCHMARK_COLS) + "\n"
                             + "\t".join(na_row) + "\n",
        # sidecars merged by concat_scores; -1 = failed / not monitored
        bench_stem + ".time": "-1",
        bench_stem + ".gpu": "-1",
    }


def write_files(job: KilledJob, tool_dir: str, sample: str) -> None:
    for path, content in files_to_write(job, tool_dir, sample).items():
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as out_file:
            out_file.write(content)
    # record the reason next to the aligner's own log
    err_log = os.path.join(tool_dir, "logs", sample + ".err")
    os.makedirs(os.path.dirname(err_log), exist_ok=True)
    with open(err_log, "a") as err_file:
        print(
            f"SLURM job {job.slurm_id} ended with {job.state}; "
            "marked as failed by tools/mark_slurm_failed.py",
            file=err_file,
        )


def main() -> None:
    logs = find_logs()
    # keep the latest failure per output
    jobs: Dict[str, KilledJob] = {}
    for log in logs:
        for job in parse_log(log):
            if job.state in args.states and in_runs(job):
                jobs[job.output] = job
    runs = f" of run(s) {' '.join(args.runs)}" if args.runs else ""
    print(f"Scanned {len(logs)} log(s), found {len(jobs)} job(s){runs} in "
          f"state {' or '.join(args.states)}.")
    if not jobs and args.runs and not (args.logs or args.all_logs):
        print("The newest log may be from another run, see --all-logs.")

    counts: Dict[str, int] = {}
    for output, job in sorted(jobs.items()):
        handled = sample_of(job)
        if handled is None:
            action = f"skipped, rule {job.rule} is not handled"
        elif os.path.exists(output):
            action = "skipped, output exists"
        elif args.dry:
            action = "would mark failed"
        else:
            write_files(job, *handled)
            action = "marked failed"
        print(f"{job.state:<14} {job.slurm_id:>10}  {output}  ({action})")
        counts[action] = counts.get(action, 0) + 1

    for action, n in sorted(counts.items()):
        print(f"{n} {action}")


if __name__ == '__main__':
    args = argparser.parse_args()
    main()
