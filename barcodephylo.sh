#!/usr/bin/env bash
# ============================================================
# barcodephylo.sh — Build a phylogenetic tree from DNA barcodes
# ============================================================
# Usage:
#   bash barcodephylo.sh <data.list> <outgroup_label> "<marker_list>" <prefix>
#
# Example:
#   bash barcodephylo.sh data.list Diaporthella_corylina_CBS_121124 "ITS LSU TEF TUB HIS CAL" Diaporthe
#
# Description:
#   1. Align barcode sequences with MAFFT
#   2. Trim alignments using trimAl
#   3. Concatenate and determine the best model via iqtree_modelfinder.py
#   4. Infer a phylogenetic tree with IQ-TREE2
# ============================================================

# ----------- Step 1. Parse input arguments -----------
DATA=$1            # Optional: a file list (for future use)
OUTGROUP=$2        # Outgroup taxon name
BARCODE=$3         # Space-separated list of marker names
PREFIX=$4          # Prefix for file naming (e.g., genus name)

# ----------- Step 2. Validate inputs -----------
if [ -z "$DATA" ] || [ -z "$OUTGROUP" ] || [ -z "$BARCODE" ] || [ -z "$PREFIX" ]; then
  echo "Usage: bash barcodephylo.sh <data.list> <outgroup_label> \"<marker_list>\" <prefix>"
  echo "Example: bash barcodephylo.sh data.list Diaporthella_corylina_CBS_121124 \"ITS LSU TEF TUB HIS CAL\" Diaporthe"
  exit 1
fi

# ----------- Step 3. Create working directories -----------
mkdir -p 02_mafft 03_trimal 04_modelfinder 05_iqtree

# ----------- Step 4. Align sequences with MAFFT -----------
echo "[INFO] Running MAFFT alignments..."
ls 01_data/ | grep fasta | while read a; do
  # Localpair + adjustdirection: recommended for barcode regions
  echo "mafft --localpair --thread 4 --adjustdirection 01_data/${a} > 02_mafft/${a%.fasta}.mafft.fna"
done > mafft.sh

bash mafft.sh

# Remove reverse complement indicators (e.g. “_R_”) added by MAFFT
sed -i 's/>_R_/>/' 02_mafft/*.mafft.fna

# ----------- Step 5. Trim alignments with trimAl -----------
echo "[INFO] Trimming alignments with trimAl..."
ls 01_data/ | grep fasta | sed 's/.fasta//' | while read a; do
  trimal -in 02_mafft/${a}.mafft.fna \
         -gt 0.5 \
         -out 03_trimal/${a}.mafft.trimal.fna
done

# ----------- Step 6. Concatenate and select models -----------
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
  --seqtype DNA \                                 # input is DNA
  -o "$OUTGROUP" \                                # define outgroup
  --prefix 05_iqtree/iqtree_ml \                  # output prefix
  -T AUTO \                                       # auto thread detection
  -p 04_modelfinder/best_scheme.txt \             # partition file
  --ufboot 1000 \                                 # 1000 ultrafast bootstraps
  --alrt 1000 &                                   # 1000 SH-aLRT replicates

echo "[INFO] IQ-TREE analysis started."
echo "[INFO] Monitor progress: tail -f nohup.out"