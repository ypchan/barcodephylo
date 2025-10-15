#!/usr/bin/env python3

"""
Extract barcode-like sequences from genome assemblies in three steps:
  1) Cluster each barcode reference set with CD-HIT to produce representative sequences.
  2) BLAST those representatives against each genome assembly.
  3) Extract the best-hit region with user-defined flanking length and write to FASTA.

update: 2025-10-02
bugs: yanpengch@qq.com
"""

import argparse
import gzip
import io
import os
import sys
import tempfile
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

# =========================
# Utilities
# =========================

def require_cmd(cmd: str, human_name: str = None) -> None:
    """Ensure an external command exists in PATH; exit with a clear message if missing."""
    if shutil.which(cmd) is None:
        if human_name is None:
            human_name = cmd
        print(f"Error: `{human_name}` is required but not found in PATH. Please install it.", file=sys.stderr, flush=True)
        sys.exit(1)

def run(cmd: List[str], cwd: Path = None) -> None:
    """Run a command with subprocess.run, failing fast with helpful stdout/stderr on error."""
    try:
        subprocess.run(cmd, check=True, cwd=str(cwd) if cwd else None,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        msg = f"\n[Command failed]\n  CMD : {' '.join(cmd)}\n  CODE: {e.returncode}\n  STDOUT:\n{e.stdout}\n  STDERR:\n{e.stderr}\n"
        print(msg, file=sys.stderr, flush=True)
        sys.exit(1)

def open_maybe_gzip(path: Path, mode: str = "rt"):
    """Open a plain text file or .gz file transparently."""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)  # type: ignore
    return open(path, mode)

# =========================
# CLI parsing & input checks
# =========================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract barcode-like sequences from genome assemblies (cd-hit + blastn + flanking)."
    )
    parser.add_argument(
        "--sequence", "-s", required=True, action="append",
        help="Reference barcode sequences (multi-FASTA). May be given multiple times; pairs with --barcode_label."
    )
    parser.add_argument(
        "--barcode_label", "-b", required=True, action="append",
        help="Barcode labels (e.g., ITS TUB TEF RPB2). Must correspond 1:1 with --sequence."
    )
    parser.add_argument(
        "--flank_length", "-f", type=int, default=100,
        help="Number of bases to extract on both sides of the BLAST hit (default: 100)."
    )
    parser.add_argument(
        "--genome", "-g", required=True, action="append",
        help="Genome assemblies (.fa/.fasta/.fna/.gz). May be given multiple times; or a text file with one path per line."
    )
    parser.add_argument(
        "--prefix", "-p", type=str, default="",
        help="Prefix for output files (e.g., prefix_ITS.fasta)."
    )
    return parser.parse_args()

def check_input(seq_paths: List[str], labels: List[str]) -> None:
    """Ensure the number of sequence files equals the number of labels."""
    if len(seq_paths) != len(labels):
        print("Error: Each --sequence must have a corresponding --barcode_label.", file=sys.stderr, flush=True)
        sys.exit(1)

# =========================
# Input normalization
# =========================

def normalize_genome_list(g_args: List[str]) -> List[Path]:
    """Expand --genome/-g arguments into a deduplicated list of file paths.

    Supported forms:
      - FASTA files (plain or .gz)
      - A text list file (one path per line)

    Heuristics: if not a common FASTA suffix and the first char isn't '>', treat as list file.
    """
    genomes: List[Path] = []
    seen = set()

    for item in g_args:
        p = Path(item).expanduser().resolve()
        if not p.exists():
            print(f"Error: genome path not found: {p}", file=sys.stderr, flush=True)
            sys.exit(1)

        is_list = False
        if p.is_file() and p.suffix.lower() in {".txt", ".list"}:
            is_list = True
        else:
            if p.is_file() and p.suffix.lower() not in {".fa", ".fna", ".fasta", ".gz"}:
                is_list = True
            else:
                try:
                    with open_maybe_gzip(p, "rt") as fh:
                        ch = fh.read(1)
                        if ch != ">":
                            is_list = True
                except Exception:
                    is_list = False

        if is_list:
            with open(p, "rt") as fh:
                for line in fh:
                    path2 = Path(line.strip()).expanduser().resolve()
                    if not path2.exists():
                        print(f"Error: genome path in list not found: {path2}", file=sys.stderr, flush=True)
                        sys.exit(1)
                    if path2 not in seen:
                        seen.add(path2)
                        genomes.append(path2)
        else:
            if p not in seen:
                seen.add(p)
                genomes.append(p)

    if not genomes:
        print("Error: No valid genome paths provided.", file=sys.stderr, flush=True)
        sys.exit(1)

    return genomes

def read_genome_to_dict(genome: Path) -> Dict[str, str]:
    """Read a FASTA/FNA (optionally .gz) into a dict: {contig_id: sequence}."""
    seqs: Dict[str, List[str]] = {}
    current = None
    with open_maybe_gzip(genome, "rt") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                current = line.strip().split()[0][1:]  # drop ">", keep the first token
                if current in seqs:
                    print(f"Warning: duplicated contig id '{current}' in {genome}", file=sys.stderr)
                seqs[current] = []
            else:
                if current is None:
                    continue
                seqs[current].append(line.strip())

    return {k: "".join(v) for k, v in seqs.items()}

# =========================
# CD-HIT & BLAST
# =========================

def run_cdhit(seqs: List[str], labels: List[str], tmpdir: Path) -> Dict[str, Path]:
    """Cluster each reference set with CD-HIT, producing a representative FASTA per label.

    Returns: {label: path_to_representatives_fasta}
    """
    outdir = tmpdir / "temp_cdhit"
    outdir.mkdir(parents=True, exist_ok=True)

    label2rep: Dict[str, Path] = {}
    for fasta, label in zip(seqs, labels):
        in_fa = Path(fasta).expanduser().resolve()
        if not in_fa.exists():
            print(f"Error: sequence file not found: {in_fa}", file=sys.stderr, flush=True)
            sys.exit(1)

        # CD-HIT output prefix (it will create multiple files)
        prefix = outdir / f"{label}.cdhit.I0.9L0.8"
        cmd = [
            "cd-hit",
            "-i", str(in_fa),
            "-o", str(prefix),
            "-aL", "0.8",      # alignment coverage threshold on the longer sequence
            "-c", "0.90",      # sequence identity threshold (explicit for clarity)
        ]
        run(cmd)

        # Convert CD-HIT representative file to a standard FASTA
        src = prefix               # the rep FASTA is written to the prefix (no extension)
        dst = outdir / f"{label}.cdhit.rep.fasta"

        count = 0
        with open(src, "rt") as infh, open(dst, "wt") as outfh:
            for line in infh:
                if line.startswith(">"):
                    count += 1
                    outfh.write(f">{label}_{count}\n")
                else:
                    outfh.write(line)

        label2rep[label] = dst

    return label2rep

def run_blastn(label2rep: Dict[str, Path], genomes: List[Path], tmpdir: Path) -> Dict[Tuple[str, str], Path]:
    """BLAST each label's representative FASTA against each genome (subject mode).

    Uses outfmt 6 with explicit columns including bitscore, then writes TSV per label×genome.

    Returns: {(label, genome_basename): blast_tsv_path}
    """
    outdir = tmpdir / "temp_blastn"
    outdir.mkdir(parents=True, exist_ok=True)

    results: Dict[Tuple[str, str], Path] = {}
    for label, qfa in label2rep.items():
        for genome in genomes:
            gbase = genome.name
            out_tab = outdir / f"{label}_{gbase}.blastn.tsv"
            cmd = [
                "blastn",
                "-query", str(qfa),
                "-subject", str(genome),
                "-outfmt", "6 qseqid sseqid pident length evalue bitscore sstart send",
                "-max_target_seqs", "1000",
                "-max_hsps", "1"
            ]
            # Stream stdout to file directly for memory efficiency
            with open(out_tab, "wt") as ofh:
                try:
                    subprocess.run(cmd, check=True, stdout=ofh, stderr=subprocess.PIPE, text=True)
                except subprocess.CalledProcessError as e:
                    msg = f"\n[blastn failed]\n  CMD : {' '.join(cmd)}\n  CODE: {e.returncode}\n  STDERR:\n{e.stderr}\n"
                    print(msg, file=sys.stderr, flush=True)
                    sys.exit(1)

            results[(label, gbase)] = out_tab

    return results

# =========================
# Parse BLAST & extract windows
# =========================

def pick_best_hit(blast_tab: Path) -> Tuple[str, int, int]:
    """Pick the best BLAST hit by highest bitscore.

    Expects outfmt 6 with columns: qseqid sseqid pident length evalue bitscore sstart send.
    Returns (contig_id, sstart, send). Raises ValueError if no hits.
    """
    best = None  # (bitscore, contig, sstart, send)
    with open(blast_tab, "rt") as fh:
        for line in fh:
            if not line.strip():
                continue
            toks = line.strip().split("\t")
            if len(toks) < 8:
                continue
            contig = toks[1]
            try:
                bitscore = float(toks[5])
                sstart = int(toks[6])
                send = int(toks[7])
            except Exception:
                continue

            if best is None or bitscore > best[0]:
                best = (bitscore, contig, sstart, send)

    if best is None:
        raise ValueError("no hits")

    _, contig, sstart, send = best
    return contig, sstart, send

def slice_with_flank(seq: str, sstart: int, send: int, flank: int) -> str:
    """Convert 1-based inclusive BLAST coordinates to 0-based slice and extend by `flank` on both sides."""
    if sstart > send:
        sstart, send = send, sstart

    left = max(0, sstart - 1 - flank)
    right = min(len(seq), send + flank)  # send is inclusive in BLAST, so no -1 here

    return seq[left:right]

def extract_barcodes(labels: List[str],
                     genomes: List[Path],
                     blast_map: Dict[Tuple[str, str], Path],
                     flank_len: int) -> Dict[str, Dict[str, str]]:
    """For each label×genome, take the best hit and extract a flanked window.

    Returns: {label: {genome_basename_without_ext: SEQ}}
    Entries with no hits are omitted.
    """
    result: Dict[str, Dict[str, str]] = {lb: {} for lb in labels}

    # Cache genomes as dicts to avoid repeated I/O
    genome_cache: Dict[str, Dict[str, str]] = {}

    for (label, gbase), tab in blast_map.items():
        try:
            contig, sstart, send = pick_best_hit(tab)
        except ValueError:
            # No hits for this label×genome; skip
            continue

        # Retrieve the specific genome path by basename
        gpath_candidates = [g for g in genomes if g.name == gbase]
        if not gpath_candidates:
            continue
        gpath = gpath_candidates[0]

        if gbase not in genome_cache:
            genome_cache[gbase] = read_genome_to_dict(gpath)

        genome_dict = genome_cache[gbase]
        if contig not in genome_dict:
            # Rare case: sseqid differs from FASTA header tokenization
            # A more permissive matching (e.g., before first whitespace) could be added here.
            continue

        fragment = slice_with_flank(genome_dict[contig], sstart, send, flank_len)

        # Drop common FASTA extensions from the genome name
        name_noext = gbase
        for ext in (".fasta", ".fa", ".fna", ".fa.gz", ".fna.gz", ".fasta.gz"):
            if name_noext.endswith(ext):
                name_noext = name_noext[: -len(ext)]
                break

        result[label][name_noext] = fragment

    return result

# =========================
# Output
# =========================

def write_outputs(barcode_dict: Dict[str, Dict[str, str]], prefix: str = "") -> None:
    """Write one FASTA per label with only the genomes that had hits."""
    for label, items in barcode_dict.items():
        outf = f"{prefix + '_' if prefix else ''}{label}.fasta"
        with open(outf, "wt") as ofh:
            for genome_name, seq in items.items():
                if not seq:
                    continue
                ofh.write(f">{genome_name}\n{seq.upper()}\n")

# =========================
# Main
# =========================

def main():
    # 1) Parse CLI
    args = parse_args()
    check_input(args.sequence, args.barcode_label)

    # 2) External dependencies
    require_cmd("cd-hit", "cd-hit")
    require_cmd("blastn", "BLAST+ (blastn)")

    # 3) Normalize genomes
    genomes = normalize_genome_list(args.genome)

    # 4) Pipeline in a temporary working directory (auto-cleaned)
    with tempfile.TemporaryDirectory(prefix="temp_barcode_fisher_") as td:
        tmpdir = Path(td)
        # 4.1 CD-HIT: representatives per label
        label2rep = run_cdhit(args.sequence, args.barcode_label, tmpdir)
        # 4.2 BLAST: label×genome
        blast_map = run_blastn(label2rep, genomes, tmpdir)
        # 4.3 Extract flanked windows
        barcode_dict = extract_barcodes(
            args.barcode_label, genomes, blast_map, args.flank_length
        )

    # 5) Write outputs
    write_outputs(barcode_dict, args.prefix or "")

    # 6) Done
    sys.exit(0)


if __name__ == "__main__":
    main()
