import argparse
import os
import shutil
from pathlib import Path

import numpy as np
from learnMSA.msa_hmm.SequenceDataset import SequenceDataset


def sample_subset(
    fasta : SequenceDataset,
    ref_fasta : SequenceDataset,
    old_indices : np.ndarray,
    subset_size : int,
) -> np.ndarray:
    """
    Samples a nested subset of given size from the old indices, ensuring all
    reference sequences are included.

    Args:
        fasta: SequenceDataset object containing the full sequence set.
        ref_fasta: SequenceDataset object containing the reference sequences.
        old_indices: Array of indices representing the current subset of the
            full set.
        subset_size: Desired size of the new subset.
    """
    assert subset_size < old_indices.size,\
        "The sampled subset has to be smaller than the old subset."
    ref_indices = np.array(
        [fasta.seq_ids.index(sid) for sid in ref_fasta.seq_ids]
    )
    for i in ref_indices:
        assert np.where(old_indices == i)[0].size == 1,\
            f"Ref. index {i} is either not at all or multiple times in the "\
            "old indices array."
    homolog_indices = np.setdiff1d(old_indices, ref_indices)
    homolog_subset = np.random.choice(
        homolog_indices, size=subset_size-ref_indices.size, replace=False
    )
    return np.concatenate([ref_indices, homolog_subset], axis=0)


def make_paths(
    fasta : SequenceDataset,
    n : int,
    k : int,
    prefix_path : Path,
) -> tuple[Path, Path]:
    """
    Utility function to create paths for storing the subset fasta and the
    reference fasta.
    """
    family = os.path.splitext(os.path.basename(fasta.filename))[0]
    subset_train_path = prefix_path / "unaligned" / f"{family}_{n}_{k}"
    subset_ref_path = prefix_path / "aligned" / f"{family}_{n}_{k}"
    return subset_train_path, subset_ref_path
    

def write_fasta(
    fasta : SequenceDataset,
    subset : np.ndarray,
    k : int,
    reference_path : Path,
    prefix_path : Path,
) -> None:
    """
    Writes the subset fasta and copies the reference fasta to the appropriate
    locations.
    Args:
        fasta: SequenceDataset object containing the full sequence set.
        subset: Array of indices representing the sampled subset.
        k: Index of the nested series.
        reference_path: Path to the reference fasta file.
    """
    subset_train_path, subset_ref_path = make_paths(
        fasta, subset.size, k, prefix_path
    )

    # Ensure output directories exist
    subset_train_path.parent.mkdir(parents=True, exist_ok=True)
    subset_ref_path.parent.mkdir(parents=True, exist_ok=True)

    if not os.path.exists(subset_train_path):
        # write the subset fasta
        with open(subset_train_path, "w") as train_file:
            for i in subset:
                train_file.write(
                    ">"+fasta.seq_ids[i]+"\n"+str(fasta.get_record(i).seq)+"\n"
                )
        # copy the reference msa
        shutil.copyfile(reference_path, subset_ref_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate nested subsets for size scaling experiments."
    )
    parser.add_argument(
        "--benchmark", type=Path, required=True,
        help="Path to benchmark collection containing 'aligned' and "\
            "'unaligned' subfolders."
    )
    parser.add_argument(
        "--lengths", nargs="+", type=int, required=True,
        help="List of sequence counts L."
    )
    parser.add_argument(
        "--n_series", type=int, default=3,
        help="Number of nested series per dataset."
    )
    parser.add_argument(
        "--prefix", type=Path, default=Path("data/size_scaling"),
        help="Base directory for storing subsets."
    )
    args = parser.parse_args()

    L = sorted(args.lengths)
    min_len, max_len = L[0], L[-1]

    unaligned_dir = args.benchmark / "unaligned"
    aligned_dir = args.benchmark / "aligned"

    c = 0

    for unaligned_file in unaligned_dir.iterdir():
        ref_file = aligned_dir / unaligned_file.name
        if not ref_file.exists():
            print(
                f"Skipping {unaligned_file.name}: no matching reference in "\
                "'aligned'."
            )
            continue

        fasta = SequenceDataset(unaligned_file)
        ref_fasta = SequenceDataset(ref_file)
        if fasta.num_seq < max_len:
            print(
                f"Skipping {unaligned_file.name} "\
                "(to include, reduce max(L))."
            )
            continue
        if ref_fasta.num_seq > min_len:
            print(
                f"Skipping {unaligned_file.name} "\
                "(to include, increase min(L))."
            )
            continue

        c += 1

        for k in range(args.n_series):
            old_indices = np.arange(fasta.num_seq)
            for size in reversed(L):
                subset = sample_subset(fasta, ref_fasta, old_indices, size)
                write_fasta(fasta, subset, k, ref_file, args.prefix)
                old_indices = subset

    print(f"Processed {c} families.")


if __name__ == "__main__":
    main()
