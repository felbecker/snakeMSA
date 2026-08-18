#!/usr/bin/env python3
"""Propagate the gap pattern of an MSA onto matching unaligned sequences.

Given an aligned FASTA file (an MSA, possibly over one alphabet, e.g. amino
acids) and an unaligned FASTA file whose sequences use another alphabet
(e.g. 3Di) but share identifiers and per-sequence residue counts with the MSA,
this script induces an equivalent alignment on the unaligned sequences: every
non-gap column of the MSA is replaced by the corresponding character of the
matching unaligned sequence, while gaps are kept in place.
"""

import argparse
import sys
from typing import Dict, Iterator, List, Tuple, TextIO

GAP_CHARS: str = "-."


def read_fasta(handle: TextIO) -> Iterator[Tuple[str, str]]:
    """Yield ``(identifier, sequence)`` pairs from a FASTA stream.

    The identifier is the first whitespace-delimited token of the header line.
    """
    identifier: str = ""
    chunks: List[str] = []
    have_record: bool = False
    for line in handle:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if have_record:
                yield identifier, "".join(chunks)
            identifier = line[1:].split()[0] if line[1:].split() else ""
            chunks = []
            have_record = True
        elif have_record:
            chunks.append(line.strip())
    if have_record:
        yield identifier, "".join(chunks)


def read_fasta_dict(path: str) -> Dict[str, str]:
    """Read a FASTA file into an ``identifier -> sequence`` dictionary."""
    records: Dict[str, str] = {}
    with open(path, "r") as handle:
        for identifier, sequence in read_fasta(handle):
            if not identifier:
                raise ValueError(f"Empty identifier in FASTA file '{path}'.")
            if identifier in records:
                raise ValueError(
                    f"Duplicate identifier '{identifier}' in FASTA file "
                    f"'{path}'."
                )
            records[identifier] = sequence
    return records


def propagate_gaps(aligned: str, unaligned: str, identifier: str) -> str:
    """Replace non-gap characters of ``aligned`` by those of ``unaligned``.

    Gaps in ``aligned`` are preserved. The number of non-gap characters in
    ``aligned`` must equal the length of ``unaligned``.
    """
    residues: int = sum(1 for char in aligned if char not in GAP_CHARS)
    if residues != len(unaligned):
        raise ValueError(
            f"Length mismatch for '{identifier}': the MSA row has {residues} "
            f"non-gap characters but the unaligned sequence has "
            f"{len(unaligned)} characters."
        )
    result: List[str] = []
    pos: int = 0
    for char in aligned:
        if char in GAP_CHARS:
            result.append(char)
        else:
            result.append(unaligned[pos])
            pos += 1
    return "".join(result)


def write_fasta(
    handle: TextIO,
    records: Iterator[Tuple[str, str]],
    line_width: int = 60,
) -> None:
    """Write ``(identifier, sequence)`` pairs as FASTA to ``handle``."""
    for identifier, sequence in records:
        handle.write(f">{identifier}\n")
        if line_width > 0:
            for start in range(0, len(sequence), line_width):
                handle.write(sequence[start:start + line_width] + "\n")
        else:
            handle.write(sequence + "\n")


def propagate_msa(
    msa_path: str,
    unaligned_path: str,
    out_handle: TextIO,
    line_width: int = 60,
) -> None:
    """Induce the MSA gap pattern onto the matching unaligned sequences."""
    unaligned: Dict[str, str] = read_fasta_dict(unaligned_path)
    with open(msa_path, "r") as msa_handle:
        aligned_records = list(read_fasta(msa_handle))

    missing: List[str] = [
        identifier
        for identifier, _ in aligned_records
        if identifier not in unaligned
    ]
    if missing:
        raise ValueError(
            "No matching unaligned sequence for identifier(s): "
            + ", ".join(missing)
        )

    def induced() -> Iterator[Tuple[str, str]]:
        for identifier, aligned in aligned_records:
            yield identifier, propagate_gaps(
                aligned, unaligned[identifier], identifier
            )

    write_fasta(out_handle, induced(), line_width)


def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Propagate the gap pattern of an MSA onto matching unaligned "
            "sequences (possibly over a different alphabet), inducing an "
            "equivalent alignment."
        )
    )
    parser.add_argument(
        "msa",
        help="Aligned FASTA file (the MSA whose gap pattern is propagated).",
    )
    parser.add_argument(
        "unaligned",
        help=(
            "Unaligned FASTA file with matching identifiers and per-sequence "
            "residue counts (may use a different alphabet)."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        default="-",
        help="Output FASTA path ('-' for stdout, the default).",
    )
    parser.add_argument(
        "-w",
        "--line-width",
        type=int,
        default=60,
        help="Wrap sequences at this width (0 to disable wrapping).",
    )
    return parser.parse_args(argv)


def main(argv: List[str]) -> int:
    args = parse_args(argv)
    try:
        if args.output == "-":
            propagate_msa(
                args.msa, args.unaligned, sys.stdout, args.line_width
            )
        else:
            with open(args.output, "w") as out_handle:
                propagate_msa(
                    args.msa, args.unaligned, out_handle, args.line_width
                )
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
