import argparse
import shutil
from pathlib import Path

import numpy as np
from learnMSA.util.dataset import Dataset
from learnMSA.util.embedding_cache import EmbeddingCache
from learnMSA.util.embedding_dataset import EmbeddingDataset
from learnMSA.util.sequence_dataset import SequenceDataset


def sample_subset(
    seq_ids : list[str],
    ref_ids : list[str],
    old_indices : np.ndarray,
    depth : int,
) -> np.ndarray:
    """
    Samples a nested subset of given depth from the old indices, ensuring all
    reference sequences are included. If the depth equals the depth of the old
    subset, the old subset is returned unchanged.

    Args:
        seq_ids: Sequence ids of the full sequence set.
        ref_ids: Sequence ids of the reference sequences.
        old_indices: Array of indices representing the current subset of the
            full set.
        depth: Desired depth (number of sequences) of the new subset.
    """
    assert depth <= old_indices.size,\
        "The sampled subset can not be larger than the old subset."
    if depth == old_indices.size:
        return old_indices
    index = {sid: i for i, sid in enumerate(seq_ids)}
    ref_indices = np.array([index[sid] for sid in ref_ids])
    for i in ref_indices:
        assert np.where(old_indices == i)[0].size == 1,\
            f"Ref. index {i} is either not at all or multiple times in the "\
            "old indices array."
    homolog_indices = np.setdiff1d(old_indices, ref_indices)
    homolog_subset = np.random.choice(
        homolog_indices, size=depth-ref_indices.size, replace=False
    )
    return np.concatenate([ref_indices, homolog_subset], axis=0)


def find_file(directory : Path, family : str) -> Path | None:
    """
    Finds the file belonging to a family in a directory. A file belongs to the
    family if its name is the family name or starts with the family name
    followed by a dot (e.g. 'fam.fasta', 'fam.3di.npz').
    """
    matches = [
        f for f in directory.iterdir()
        if f.is_file() and (f.name == family or f.name.startswith(family+"."))
    ]
    if len(matches) > 1:
        raise ValueError(
            f"Family {family} matches multiple files in {directory}: "\
            f"{sorted(f.name for f in matches)}."
        )
    return matches[0] if matches else None


def load_dataset(path : Path) -> Dataset:
    """
    Loads a fasta-like sequence file or an embedding file (.npz), e.g. as
    written by learnMSA-3di.
    """
    if path.suffix == ".npz":
        return EmbeddingDataset(path)
    return SequenceDataset(path, "fasta")


def output_path(
    directory : Path,
    source : Path,
    family : str,
    depth : int,
    k : int,
) -> Path:
    """
    The subset keeps everything in the name of the source file that follows the
    family name (the extension).
    """
    return directory / f"{family}_{depth}_{k}{source.name[len(family):]}"


def write_subset(dataset : Dataset, subset : np.ndarray, path : Path) -> None:
    """
    Writes the sequences with the given indices in the given order. Sequence
    files are written as fasta, embedding files as npz.
    """
    if isinstance(dataset, EmbeddingDataset):
        embeddings = [dataset.get_encoded_seq(int(i)) for i in subset]
        seq_lens = np.array([e.shape[0] for e in embeddings])
        cache = EmbeddingCache(
            seq_lens, embeddings[0].shape[1], cache=np.concatenate(embeddings)
        )
        subset_dataset = EmbeddingDataset(
            embedding_cache=cache,
            seq_ids=[dataset.seq_ids[i] for i in subset]
        )
        # EmbeddingDataset.write appends '.npz' if missing, so write to a
        # temporary name that has it and rename.
        tmp_path = path.with_name(path.name + ".tmp.npz")
        subset_dataset.write(tmp_path, metadata=dataset.metadata)
        tmp_path.rename(path)
    else:
        with open(path, "w") as file:
            for i in subset:
                file.write(
                    ">"+dataset.seq_ids[i]+"\n"+
                    str(dataset.get_record(i).seq)+"\n" # type: ignore
                )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate nested subsets for size scaling experiments."
    )
    parser.add_argument(
        "--ref", type=Path, required=True,
        help="Directory containing the reference alignments."
    )
    parser.add_argument(
        "--data", type=Path, nargs="+", required=True,
        help="One or more directories with per-family files that all hold "\
            "the same sequences: fasta files or .npz embedding files, e.g. "\
            "ProstT5 outputs written by learnMSA-3di. The order of "\
            "sequences may differ between directories. Files are matched "\
            "to each other and to --ref by family name, i.e. the file "\
            "name without extension. The families are taken from the first "\
            "directory, which also determines the sequence counts. Subsets "\
            "are written to subfolders of --prefix named like the input "\
            "directories."
    )
    parser.add_argument(
        "--depths", nargs="+", type=int, required=True,
        help="List of dataset depths (numbers of sequences) D."
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

    out_names = [d.name for d in args.data]
    if len(set(out_names)) != len(out_names) or "aligned" in out_names:
        parser.error(
            "The --data directories need distinct names that differ from "\
            "'aligned', which is used for the references."
        )

    D = sorted(args.depths)
    min_depth, max_depth = D[0], D[-1]

    c = 0

    for main_file in sorted(args.data[0].iterdir()):
        family = main_file.stem
        ref_file = find_file(args.ref, family)
        if ref_file is None:
            print(
                f"Skipping {main_file.name}: no matching reference in "\
                f"{args.ref}."
            )
            continue
        data_files = [main_file]
        for directory in args.data[1:]:
            data_file = find_file(directory, family)
            if data_file is None:
                break
            data_files.append(data_file)
        if len(data_files) < len(args.data):
            print(f"Skipping {main_file.name}: no matching file in {directory}.")
            continue

        datasets = [load_dataset(f) for f in data_files]
        ref_fasta = SequenceDataset(ref_file)
        fasta = datasets[0]
        if fasta.num_seq < max_depth:
            print(
                f"Skipping {main_file.name} (to include, reduce max(D))."
            )
            continue
        if ref_fasta.num_seq > min_depth:
            print(
                f"Skipping {main_file.name} (to include, increase min(D))."
            )
            continue
        try:
            # Bring all datasets into the order of the first one, such that
            # the same indices refer to the same sequences everywhere.
            for dataset in datasets[1:]:
                dataset.adapt_order(fasta)
        except ValueError:
            print(
                f"Skipping {main_file.name}: the sequence ids in the --data "\
                "directories do not match."
            )
            continue

        for directory in [args.prefix / n for n in out_names]\
                + [args.prefix / "aligned"]:
            directory.mkdir(parents=True, exist_ok=True)

        c += 1

        for k in range(args.n_series):
            old_indices = np.arange(fasta.num_seq)
            for depth in reversed(D):
                subset = sample_subset(
                    fasta.seq_ids, ref_fasta.seq_ids, old_indices, depth
                )
                for dataset, data_file, out_name in zip(
                    datasets, data_files, out_names
                ):
                    path = output_path(
                        args.prefix / out_name, data_file, family, depth, k
                    )
                    if not path.exists():
                        write_subset(dataset, subset, path)
                ref_path = output_path(
                    args.prefix / "aligned", ref_file, family, depth, k
                )
                if not ref_path.exists():
                    shutil.copyfile(ref_file, ref_path)
                old_indices = subset

    print(f"Processed {c} families.")


if __name__ == "__main__":
    main()
