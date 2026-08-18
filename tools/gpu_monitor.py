"""Out-of-process peak GPU memory tracking for benchmarked subprocesses.

Backend agnostic: per-process device memory is read from ``nvidia-smi``, so this
works for PyTorch, TensorFlow, JAX or any other CUDA program without the program
having to report anything itself.

Caveats on the reported number:

* It is the *allocator high-water mark*, not the sum of live tensors, and it
  includes the CUDA context (~300 MiB). PyTorch's caching allocator never
  returns memory to the driver. This is the right number for "how much GPU does
  this job need"; it is consistent between backends (measured 1322 MiB for torch
  vs 1282 MiB for TensorFlow for an identical 1 GiB allocation).
* TensorFlow preallocates the whole device by default, in which case the number
  is the pool size rather than the peak (a 1 GiB allocation reports 22402 MiB on
  a 24 GiB card). learnMSA is not affected -- it sets
  ``TF_GPU_ALLOCATOR=cuda_malloc_async`` in learnMSA/run/util.py, which does not
  preallocate -- but another TF-based tool would be. Prefix such a tool's
  command in the config with ``TF_FORCE_GPU_ALLOW_GROWTH=true`` to get a real
  number; note that TF wants the literal string ``true``, ``1`` does not work.
* On MIG partitions or in containers that hide PIDs, the per-process query can
  return no rows and the peak stays 0. A peak of 0 means "monitored, nothing
  attributed"; ``NULL_GPU_MEM`` (-1) means "not monitored at all".
"""

import shutil
import subprocess
import threading
import time

import psutil

#: Written to the benchmark when no GPU tracking took place (CPU job, or no
#: nvidia-smi on this machine).
NULL_GPU_MEM = -1

_QUERY = [
    "nvidia-smi",
    "--query-compute-apps=pid,used_gpu_memory",
    "--format=csv,noheader,nounits",
]

#: Poll fast at first so that short jobs are not missed, then back off. Mirrors
#: snakemake's own BENCHMARK_INTERVAL_SHORT / BENCHMARK_INTERVAL split.
FAST_INTERVAL = 0.2
FAST_PERIOD = 10.0
SLOW_INTERVAL = 1.0


def gpu_available():
    """Whether per-process GPU memory can be queried on this machine."""
    if shutil.which("nvidia-smi") is None:
        return False
    try:
        return subprocess.run(_QUERY, capture_output=True, timeout=30).returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False


def _query_used_memory():
    """Map pid -> used GPU memory in MiB, summed over all devices."""
    out = subprocess.run(
        _QUERY, capture_output=True, text=True, timeout=30
    ).stdout
    used = {}
    for line in out.strip().splitlines():
        fields = [f.strip() for f in line.split(",")]
        if len(fields) != 2:
            continue
        try:
            pid, mem = int(fields[0]), int(fields[1])
        except ValueError:
            # e.g. "[N/A]" for processes whose memory the driver cannot report
            continue
        used[pid] = used.get(pid, 0) + mem
    return used


class GPUMemoryMonitor:
    """Track the peak GPU memory of a process tree while it runs.

    The root pid is usually a shell (``subprocess.Popen(..., shell=True)``) and
    the process actually touching the GPU is a descendant -- ``conda run`` adds
    another hop on top -- so the whole tree is attributed to the job.
    """

    def __init__(self, pid):
        self.pid = int(pid)
        self.peak_mb = 0
        self._known = {self.pid}
        self._stop = threading.Event()
        self._thread = None

    def start(self):
        self._thread = threading.Thread(target=self._poll, daemon=True)
        self._thread.start()
        return self

    def stop(self):
        """Stop polling and return the peak in MiB."""
        if self._thread is not None:
            self._stop.set()
            self._thread.join(timeout=60)
            self._thread = None
        return self.peak_mb

    def __enter__(self):
        return self.start()

    def __exit__(self, *exc):
        self.stop()
        return False

    def _tick(self):
        try:
            self._known.update(
                p.pid for p in psutil.Process(self.pid).children(recursive=True)
            )
        except psutil.Error:
            # process tree gone or not inspectable; keep the pids seen so far
            pass
        used = _query_used_memory()
        total = sum(mem for pid, mem in used.items() if pid in self._known)
        self.peak_mb = max(self.peak_mb, total)

    def _poll(self):
        start = time.time()
        while not self._stop.is_set():
            try:
                self._tick()
            except Exception:
                # monitoring must never be able to fail the job it observes
                pass
            elapsed = time.time() - start
            interval = FAST_INTERVAL if elapsed < FAST_PERIOD else SLOW_INTERVAL
            self._stop.wait(interval)
