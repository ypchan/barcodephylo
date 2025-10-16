#!/usr/bin/env bash
# ============================================================
# barcodephylo.sh — Build a phylogenetic tree from DNA barcodes
# ============================================================
# Usage:
#   bash barcodephylo.sh <data_dir> <outgroup_label> "<marker_list>" <prefix>
#
# Arguments:
#   <data_dir>       Directory containing FASTA files, e.g. 01_data/
#   <outgroup_label> Taxon name used as outgroup (must match sequence header)
#   "<marker_list>"  Space-separated list of marker names, e.g. "ITS LSU TEF TUB HIS CAL"
#   <prefix>         File prefix shared by all FASTA files, e.g. Diaporthe
#
# Example:
#   bash barcodephylo.sh 01_data Diaporthella_corylina_CBS_121124 "ITS LSU TEF TUB HIS CAL" Diaporthe
#
# Description:
#   1. Align barcode sequences with MAFFT
#   2. Trim alignments with trimAl
#   3. Concatenate and find the best partition model via iqtree_modelfinder.py
#   4. Infer the phylogenetic tree with IQ-TREE2
# ============================================================

# ----------- Step 1. Parse input arguments -----------
DATA_DIR=$1      # Directory containing FASTA files
OUTGROUP=$2      # Outgroup taxon
BARCODE=$3       # Space-separated marker names
PREFIX=$4        # Prefix (genus name, e.g. Diaporthe)

# ----------- Step 2. Validate inputs -----------
if [ -z "$DATA_DIR" ] || [ -z "$OUTGROUP" ] || [ -z "$BARCODE" ] || [ -z "$PREFIX" ]; then
  echo "Usage: bash barcodephylo.sh <data_dir> <outgroup_label> \"<marker_list>\" <prefix>"
  echo
  echo "Example:"
  echo "  bash barcodephylo.sh 01_data Diaporthella_corylina_CBS_121124 \"ITS LSU TEF TUB HIS CAL\" Diaporthe"
  echo
  echo "Arguments:"
  echo "  data_dir       Folder containing FASTA files (e.g. 01_data/Diaporthe_ITS.fasta)"
  echo "  outgroup_label Sequence name to be used as outgroup"
  echo "  marker_list    Quoted list of markers, separated by spaces"
  echo "  prefix         Common prefix for filenames"
  exit 1
fi

if [ ! -d "$DATA_DIR" ]; then
  echo "[ERROR] Data directory '$DATA_DIR' not found."
  exit 1
fi

# ----------- Step 3. Create working directories -----------
mkdir -p 02_mafft 03_trimal 04_modelfinder 05_iqtree

# ----------- Step 4. Align sequences with MAFFT -----------
echo "[INFO] Running MAFFT alignments..."
for fasta in "$DATA_DIR"/*.fasta; do
  base=$(basename "$fasta" .fasta)
  echo "[INFO] Aligning $base..."
  mafft --localpair --thread 4 --adjustdirection "$fasta" > "02_mafft/${base}.mafft.fna"
done

# Remove reverse complement headers if added by MAFFT
sed -i 's/>_R_/>/' 02_mafft/*.mafft.fna

# ----------- Step 5. Trim alignments with trimAl -----------
echo "[INFO] Trimming alignments..."
for f in 02_mafft/*.mafft.fna; do
  base=$(basename "$f" .mafft.fna)
  trimal -in "$f" -gt 0.5 -out "03_trimal/${base}.mafft.trimal.fna"
done

# ----------- Step 6. Concatenate and find best model -----------
echo "[INFO] Running iqtree_modelfinder.py..."
iqtree_modelfinder.py \
  $(for marker in $BARCODE; do echo -n "03_trimal/${PREFIX}_${marker}.mafft.trimal.fna "; done) \
  -o 04_modelfinder \
  --mrbayes_nexus \
  --outgroup "$OUTGROUP"

# ----------- Step 7. Infer phylogeny with IQ-TREE2 -----------
echo "[INFO] Building phylogeny with IQ-TREE2..."
nohup iqtree2 \
  -s 04_modelfinder/concatenated.fna \            # concatenated alignment
  --seqtype DNA \                                 # DNA data
  -o "$OUTGROUP" \                                # outgroup label
  --prefix 05_iqtree/iqtree_ml \                  # output prefix
  -T AUTO \                                       # automatic CPU threads
  -p 04_modelfinder/best_scheme.txt \             # partition scheme
  --ufboot 1000 \                                 # ultrafast bootstrap replicates
  --alrt 1000 &                                   # SH-aLRT tests

echo "[INFO] IQ-TREE analysis started successfully."
echo "[INFO] Monitor progress using: tail -f nohup.out"
