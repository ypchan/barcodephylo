#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Find best substitution models per alignment with IQ-TREE ModelFinder
Optionally generating a MrBayes-ready NEXUS, and log all steps with a startup banner.

contact: yanpengch@qq.com
update: 2025-09-16
"""

from __future__ import annotations

import argparse
import datetime as dt
import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Iterable

# ====== Mapping: IQ-TREE models -> MrBayes settings ======
MODEL_MAP = {
    'GTR': 'nst=6',
    'GTR+I': 'nst=6 rates=propinv',
    'GTR+G': 'nst=6 rates=gamma',
    'GTR+I+G': 'nst=6 rates=invgamma',

    'SYM': ['nst=6', 'statefreqpr=fixed(equal)'],
    'SYM+I': ['nst=6 rates=propinv', 'statefreqpr=fixed(equal)'],
    'SYM+G': ['nst=6 rates=gamma', 'statefreqpr=fixed(equal)'],
    'SYM+I+G': ['nst=6 rates=invgamma', 'statefreqpr=fixed(equal)'],

    'HKY': 'nst=2',
    'HKY+I': 'nst=2 rates=propinv',
    'HKY+G': 'nst=2 rates=gamma',
    'HKY+I+G': 'nst=2 rates=invgamma',

    'K2P': ['nst=2', 'statefreqpr=fixed(equal)'],
    'K2P+I': ['nst=2 rates=propinv', 'statefreqpr=fixed(equal)'],
    'K2P+G': ['nst=2 rates=gamma', 'statefreqpr=fixed(equal)'],
    'K2P+I+G': ['nst=2 rates=invgamma', 'statefreqpr=fixed(equal)'],

    'F81': 'nst=1',
    'F81+I': 'nst=1 rates=propinv',
    'F81+G': 'nst=1 rates=gamma',
    'F81+I+G': 'nst=1 rates=invgamma',

    'JC': ['nst=1', 'statefreqpr=fixed(equal)'],
    'JC+I': ['nst=1 rates=propinv', 'statefreqpr=fixed(equal)'],
    'JC+G': ['nst=1 rates=gamma', 'statefreqpr=fixed(equal)'],
    'JC+I+G': ['nst=1 rates=invgamma', 'statefreqpr=fixed(equal)'],
}

# ====== Logging & Banner ======
def setup_logger(outdir: Path, level: str = "INFO") -> logging.Logger:
    """
    Configure logging to file and stdout. Generate a startup banner at the very beginning
    with script name, purpose, and last modified time.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    log_file = outdir / "iqtree_modelfinder.log"

    logger = logging.getLogger("iqtree_modelfinder")
    logger.setLevel(getattr(logging, level.upper(), logging.INFO))
    logger.handlers.clear()

    fmt = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s")

    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    # ---- Banner (first output lines) ----
    try:
        script_path = Path(__file__).resolve()
        script_name = script_path.name
        mtime = dt.datetime.fromtimestamp(script_path.stat().st_mtime).isoformat()
    except Exception:
        script_name = Path(sys.argv[0]).name if sys.argv else "iqtree_modelfinder.py"
        mtime = "unknown"

    purpose = (
        "Find best substitution models per alignment with IQ-TREE ModelFinder and "
        "optionally generate a MrBayes-ready NEXUS."
    )
    now = dt.datetime.now().isoformat()

    banner_lines = [
        "=" * 70,
        f" Script : {script_name}",
        f" Purpose: {purpose}",
        f" Updated: {mtime}",
        f" Started: {now}",
        "=" * 70,
    ]
    # Emit banner as the very first log records
    for line in banner_lines:
        logger.info(line)

    return logger

# ====== FASTA utilities ======
def read_fasta(path: Path) -> Dict[str, str]:
    """Read FASTA (multi-line sequences, ignore blank lines) -> {taxon: sequence}."""
    seqs: Dict[str, List[str]] = {}
    current = None
    with path.open("rt", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].strip()
                if current in seqs:
                    raise ValueError(f"Duplicate taxon '{current}' in {path}")
                seqs[current] = []
            else:
                if current is None:
                    raise ValueError(f"FASTA format error (sequence before header) in {path}")
                seqs[current].append(line)
    return {k: "".join(v) for k, v in seqs.items()}

def concatenate_msa(msa_files: Iterable[Path], outdir: Path, logger: logging.Logger) -> Dict[str, str]:
    """
    Concatenate multiple alignments; pad missing taxa with '?'.
    Write 'concatenated.fna' and return {taxon: concatenated_seq}.
    """
    logger.info("Concatenating %d alignments", len(list(msa_files)))
    all_dicts = {p: read_fasta(p) for p in msa_files}
    taxa = sorted({t for d in all_dicts.values() for t in d.keys()})

    # Determine per-alignment length (by first sequence)
    lengths: Dict[Path, int] = {}
    for p, d in all_dicts.items():
        if not d:
            raise ValueError(f"Empty alignment: {p}")
        lengths[p] = len(next(iter(d.values())))

    # Pad missing taxa with '?'
    for p, d in all_dicts.items():
        L = lengths[p]
        for t in taxa:
            if t not in d:
                d[t] = "?" * L

    # Concatenate in the given input order
    concatenated: Dict[str, str] = {}
    for t in taxa:
        concatenated[t] = "".join(all_dicts[p][t] for p in msa_files)

    # Write output
    out_path = outdir / "concatenated.fna"
    with out_path.open("wt", encoding="utf-8") as out:
        for t, seq in concatenated.items():
            out.write(f">{t}\n{seq}\n")
    logger.info("Wrote concatenated FASTA: %s", out_path)
    logger.info("File mtime: %s", dt.datetime.fromtimestamp(out_path.stat().st_mtime).isoformat())
    return concatenated

def get_partition_lengths(msa_files: Iterable[Path]) -> Dict[str, int]:
    """Return {partition_label: alignment_length} in the same order as inputs."""
    part_len: Dict[str, int] = {}
    for p in msa_files:
        label = p.stem
        d = read_fasta(p)
        if not d:
            raise ValueError(f"Empty alignment: {p}")
        length = len(next(iter(d.values())))
        part_len[label] = length
    return part_len

# ====== ModelFinder runner ======
def clean_iqtree_model(model: str) -> str:
    """
    Normalize IQ-TREE model token:
    - drop '+F'
    - map '+Gk' -> '+G' (k: number of discrete gamma categories)
    """
    m = model.replace("+F", "")
    m = re.sub(r"\+G\d+", "+G", m)
    return m

def run_modelfinder(
    msa_files: Iterable[Path],
    model_restriction: str,
    outdir: Path,
    threads: int,
    iqtree_bin: str,
    logger: logging.Logger,
) -> Dict[str, str]:
    """
    Run IQ-TREE ModelFinder for each alignment and return {partition_label: best_model}.
    Standard output and error are captured into '<prefix>.stdout.log' and '<prefix>.stderr.log'.
    """
    if shutil.which(iqtree_bin) is None:
        raise FileNotFoundError(f"Cannot find '{iqtree_bin}' in PATH")

    result: Dict[str, str] = {}

    for msa in msa_files:
        label = msa.stem
        prefix = outdir / f"{label}_modelfinder"
        if model_restriction == "mrbayes":
            cmd = [
                iqtree_bin, "-s", str(msa), "-T", str(threads),
                "--prefix", str(prefix),
                "-m", "TESTONLY",
                "--mset", "mrbayes",
                "--msub", "nuclear",
                "--redo",
            ]
        else:
            cmd = [
                iqtree_bin, "-s", str(msa), "-T", str(threads),
                "--prefix", str(prefix),
                "-m", "TESTONLY",
                "--msub", "nuclear",
                "--redo",
            ]

        logger.info("Running IQ-TREE ModelFinder for %s", msa.name)
        stdout_log = prefix.with_suffix(".stdout.log")
        stderr_log = prefix.with_suffix(".stderr.log")

        with stdout_log.open("wt", encoding="utf-8") as out_fh, stderr_log.open("wt", encoding="utf-8") as err_fh:
            try:
                subprocess.run(cmd, check=True, stdout=out_fh, stderr=err_fh)
            except subprocess.CalledProcessError as e:
                logger.error(
                    "IQ-TREE failed for %s (exit %s). See logs: %s / %s",
                    msa.name, e.returncode, stdout_log, stderr_log
                )
                raise

        # Parse best model from possible IQ-TREE outputs
        log_file = prefix.with_suffix(".log")
        best_model_raw = None
        candidates = [log_file, prefix.with_suffix(".iqtree"), stdout_log]
        for f in candidates:
            if not f.exists():
                continue
            text = f.read_text(encoding="utf-8", errors="ignore")
            m1 = re.search(r"Bayesian Information Criterion:\s*([A-Za-z0-9\+\-]+)", text)
            m2 = re.search(r"\bBIC(?:\s+best)?\s+model:\s*([A-Za-z0-9\+\-]+)", text, flags=re.I)
            if m1:
                best_model_raw = m1.group(1)
                break
            if m2:
                best_model_raw = m2.group(1)
                break

        if not best_model_raw:
            logger.error("Cannot parse best model for %s; checked: %s", msa.name, ", ".join(str(x) for x in candidates))
            raise RuntimeError(f"Cannot parse best model for {msa}")

        best_model = clean_iqtree_model(best_model_raw)
        if best_model not in MODEL_MAP:
            logger.error("Unknown/unsupported model after cleaning: '%s' (raw: '%s')", best_model, best_model_raw)
            raise ValueError(f"Unknown/unsupported model: {best_model} (raw: {best_model_raw})")

        result[label] = best_model
        logger.info("Best model for %s: %s (raw: %s)", msa.name, best_model, best_model_raw)

    return result

# ====== Writers ======
def write_best_scheme(part_len: Dict[str, int], part_model: Dict[str, str], outdir: Path, logger: logging.Logger) -> Path:
    """
    Write NEXUS 'sets' block to 'best_scheme.txt' and print a compatibility line on stdout.
    """
    out = outdir / "best_scheme.txt"
    with out.open("wt", encoding="utf-8") as fh:
        fh.write("#nexus\n")
        fh.write("begin sets;\n")

        start = 1
        for label, L in part_len.items():
            end = start + L - 1
            fh.write(f"charset {label} = {start}-{end};\n")
            start = end + 1

        charparts = ", ".join(f"{part_model[label]}:{label}" for label in part_model.keys())
        fh.write(f"charpartition ModelFinder = {charparts};\n")
        fh.write("end;\n")
    logger.info("Wrote best scheme: %s", out)
    logger.info("File mtime: %s", dt.datetime.fromtimestamp(out.stat().st_mtime).isoformat())
    # Keep original script's stdout line for compatibility
    print(f"Best scheme:{out}", file=sys.stdout, flush=True)
    return out

def write_mrbayes_nexus(
    part_model: Dict[str, str],
    part_len: Dict[str, int],
    concatenated: Dict[str, str],
    outdir: Path,
    outgroup: str | None,
    logger: logging.Logger,
) -> Path:
    """
    Generate 'run_mrbayes.nexus' ready to run in MrBayes.
    """
    mrbayes_file = outdir / "run_mrbayes.nexus"
    taxa = list(concatenated.keys())
    ntaxa = len(taxa)
    nchar = sum(part_len.values())
    longest = max(len(t) for t in taxa)

    with mrbayes_file.open("wt", encoding="utf-8") as ofh:
        ofh.write("#NEXUS\n")
        ofh.write("Begin data;\n")
        ofh.write(f"  DIMENSIONS NTAX={ntaxa} NCHAR={nchar};\n")
        ofh.write("  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;\n")
        ofh.write("  MATRIX\n")

        # Interleave into 80-column blocks
        sequences = {t: concatenated[t] for t in taxa}
        start = 0
        width = 80
        while start < nchar:
            end = min(start + width, nchar)
            for t in taxa:
                ofh.write(f"  {t:<{longest}} {sequences[t][start:end]}\n")
            ofh.write("\n")
            start = end

        ofh.write(";\nEND;\n\n")

        ofh.write("BEGIN MRBAYES;\n")
        ofh.write("  log start filename=run_mrbayes.log;\n")
        ofh.write("  [! non-interactive,no prompts, no warnings]\n")
        ofh.write("  set autoclose=yes nowarnings=yes;\n")
        if outgroup:
            ofh.write(f"  outgroup {outgroup};\n")
        ofh.write("\n")

        # Define charset blocks
        start = 1
        idx = 0
        for label, L in part_len.items():
            end = start + L - 1
            idx += 1
            ofh.write(f"  charset Subset{idx} = {start}-{end};\n")
            start = end + 1

        ofh.write(f"  partition ModelFinder = {idx}:" + ", ".join(f"Subset{i}" for i in range(1, idx+1)) + ";\n")
        ofh.write("  set partition=ModelFinder;\n\n")

        # Apply per-partition models
        i = 0
        for label in part_model.keys():
            i += 1
            model = part_model[label]
            mb = MODEL_MAP.get(model)
            if mb is None:
                raise ValueError(f"Unknown model for MrBayes: {model}")
            if isinstance(mb, list):
                ofh.write(f"  lset applyto=({i}) {mb[0]};\n")
                ofh.write(f"  prset applyto=({i}) {mb[1]};\n")
            else:
                ofh.write(f"  lset applyto=({i}) {mb};\n")

        ofh.write("\n")
        ofh.write("  prset applyto=(all) ratepr=variable;\n")
        ofh.write("  unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all) tratio=(all);\n\n")
        ofh.write("  mcmc ngen=20000000 Stoprule=yes Stopval=0.01;\n")
        ofh.write("  sump;\n  sumt;\n  log stop;\nEND;\n")

    logger.info("Wrote MrBayes NEXUS: %s", mrbayes_file)
    logger.info("File mtime: %s", dt.datetime.fromtimestamp(mrbayes_file.stat().st_mtime).isoformat())
    return mrbayes_file

# ====== CLI ======
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="iqtree_modelfinder.py",
        description="Find best substitution models per alignment with IQ-TREE ModelFinder; optionally generate MrBayes NEXUS."
    )
    p.add_argument("-i", "--input", metavar="alignment.fna", type=str, nargs="+", required=True,
                   help="Multiple sequence alignment files in FASTA format (one or more).")
    p.add_argument("-o", "--outdir", metavar="outdir", type=str, required=True, help="Output directory.")
    p.add_argument("--mrbayes_nexus", action="store_true", help="Generate a MrBayes-ready NEXUS file.")
    p.add_argument("--outgroup", metavar="Outgroup", type=str, default=None, help="Outgroup taxon for MrBayes (optional).")
    p.add_argument("-m", "--model_restriction", metavar="mrbayes|iqtree",
                   choices=["mrbayes", "iqtree"], default="mrbayes",
                   help="Restrict model search to a supported set (mrbayes|iqtree).")
    p.add_argument("--threads", type=int, default=4, help="Number of IQ-TREE threads (default: 4).")
    p.add_argument("--iqtree-bin", type=str, default="iqtree2", help="IQ-TREE executable name (default: iqtree2).")
    p.add_argument("--log-level", type=str, default="INFO", help="Logging level (DEBUG/INFO/WARN/ERROR).")
    return p.parse_args()

def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    logger = setup_logger(outdir, args.log_level)

    try:
        msa_files = [Path(x).resolve() for x in args.input]
        for p in msa_files:
            if not p.exists():
                raise FileNotFoundError(f"Input file not found: {p}")

        # Partition lengths (in input order)
        part_len = get_partition_lengths(msa_files)
        logger.info("Partitions (lengths): %s", ", ".join(f"{k}:{v}" for k, v in part_len.items()))

        # ModelFinder
        part_model = run_modelfinder(
            msa_files=msa_files,
            model_restriction=args.model_restriction,
            outdir=outdir,
            threads=args.threads,
            iqtree_bin=args.iqtree_bin,
            logger=logger,
        )

        # best_scheme.txt
        write_best_scheme(part_len, part_model, outdir, logger)

        # Concatenate and optionally write MrBayes NEXUS
        concatenated = concatenate_msa(msa_files, outdir, logger)
        if args.mrbayes_nexus:
            write_mrbayes_nexus(part_model, part_len, concatenated, outdir, args.outgroup, logger)

        logger.info("All done.")
        return 0
    except Exception as e:
        logger.exception("Failed: %s", e)
        return 1

if __name__ == "__main__":
    sys.exit(main())
