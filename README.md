# 🧬 barcodephylo

**`barcodephylo`** is a modular, reproducible bioinformatics pipeline for phylogenetic analysis based on barcode sequences (e.g., ITS, LSU, SSU, RPB2, TUB2). It streamlines the workflow from sequence collection to tree visualization, integrating tools such as **MAFFT**, **trimal**, **IQ-TREE**, **MrBayes**, and **ggtree**.

![WORKFLOW](https://img.shields.io/badge/WORKFLOW-8A2BE2)

### Workflow Overview

1. Quality control of in-house loci data from Sanger sequencing
2. Download public data using `read.GenBank` (APE package)
3. Multiple sequence alignment with MAFFT
4. Trimming MSA files
5. Model selection and MSA concatenation
6. Maximum Likelihood analysis (IQ-TREE/RAxML-NG)
7. Bayesian analysis (MrBayes)
8. Maximum Parsimony analysis (MPboot)
9. Time dating (MEGA)
10. Visualization (ggtree + Adobe Illustrator)

---

## 📦 Installation

### 🔧 Prerequisites

- Python ≥ 3.6  
- R ≥ 4.0  
- MAFFT, trimal, IQ-TREE3, MrBayes, MPI (optional)  
- R packages: `ape`, `seqinr`, `treedataverse`, etc.

```bash
git clone https://github.com/yourusername/barcodephylo.git
cd barcodephylo
bash setup.sh
```

---

## 🚀 Quick Start

### 1️⃣ Download Barcode Sequences

#### Download sequences from GenBank using `read.GenBank.R` (Windows)

> **Tip:** If you encounter errors, check:
> 1. Are accession numbers correct?
> 2. Are accessions public?
> 3. Is your internet connection stable?

```r
library(this.path)
library(openxlsx)
library(ape)
library(tidyverse)

setwd(this.dir())
taxa_tbl <- read.xlsx('../Boeremia_taxa_table_20250407.xlsx') %>% as_tibble()
species <- 'Boeremia'
marker_lst <- c('ITS', 'LSU', 'SSU', 'TEF', 'RPB2')

for (marker in marker_lst) {
  outfilename <- paste0(species, '_', marker, '.fasta')
  cat(outfilename, '\n')
  df_marker <- taxa_tbl %>% drop_na(all_of(marker))
  marker_obj <- read.GenBank(unlist(df_marker[[marker]]), seq.names = FALSE, quiet = FALSE, chunk.size = 20) # adjust the chunk.size to find the problematic accessions
  names(marker_obj) <- df_marker$longLabel
  write.FASTA(marker_obj, outfilename)
}
```

#### Add Your Own Barcode Sequences

Append your sequences to the relevant dataset files. Ensure headers and formats match existing data.

```bash
cat 01_data/00_sanger_raw_data/2021032807_ITS.fasta >> 01_data/Boeremia_ITS.fasta
cat 01_data/00_sanger_raw_data/2021032807_LROR_LR5.fasta | sed 's/_LROR_LR5//' >> 01_data/Boeremia_LSU.fasta
cat 01_data/00_sanger_raw_data/2021032807_PNS1_NS41.fasta | sed 's/_PNS1_NS41//' >> 01_data/Boeremia_SSU.fasta
cat 01_data/00_sanger_raw_data/2021032807_RPB2.fas | sed 's/_RPB2//' >> 01_data/Boeremia_RPB2.fasta
cat 01_data/00_sanger_raw_data/2021032807_TUB2.fas | sed 's/_TUB2//' >> 01_data/Boeremia_TUB2.fasta
```

---

### 2️⃣ Multiple Sequence Alignment (MAFFT)

```bash
mkdir -p 02_mafft
ls 01_data/ | grep fasta | while read a; do
  mafft --localpair --thread 4 --adjustdirection 01_data/${a} > 02_mafft/${a%.fasta}.mafft.fna
done

# Remove '_R_' prefix from headers if MAFFT adjusted reverse-complemented sequences
sed -i 's/>_R_/>/' 02_mafft/*.mafft.fna
```
> **Review alignments in AliView or similar software to remove abnormal sequences.**

---

### 3️⃣ Alignment Trimming (trimal)

```bash
mkdir -p 03_trimal
ls 01_data/ | grep fasta | sed 's/.fasta//' | while read a; do
  trimal -in 02_mafft/${a}.mafft.fna -gt 0.5 -out 03_trimal/${a}.mafft.trimal.fna
done
```

---

### 4️⃣ Concatenate MSAs & Model Selection

```bash
rm -rf 04_modelfinder
mkdir 04_modelfinder
outgroup_label="Phoma_herbarum_CBS_615.75"

iqtree_modelfinder.py -i $(ls 03_trimal/*) -o 04_modelfinder --mrbayes_nexus --outgroup ${outgroup_label}

# Specify file order if needed:
iqtree_modelfinder.py -i barcode1.mafft.trimal.fasta barcode3.mafft.trimal.fasta barcode2.mafft.trimal.fasta -o 04_modelfinder --mrbayes_nexus --outgroup ${outgroup_label}
```

---

### 5️⃣ Maximum Likelihood Tree (IQ-TREE3)

```bash
rm -rf 05_iqtree
mkdir 05_iqtree

nohup iqtree3 -s 04_modelfinder/concatenated.fna \
      --seqtype DNA -o ${outgroup_label} \
      --prefix 05_iqtree/iqtree_ml -T AUTO \
      -p 04_modelfinder/best_scheme.txt \
      --ufboot 1000 --alrt 1000 &
```

---

### 6️⃣ Bayesian Phylogenetics (MrBayes)

```bash
mkdir -p 06_mrbayes
nohup mpirun -n 4 mb < run_mrbayes.sh &
# Or:
# nohup bash mb < run_mrbayes &
```

---

### 7️⃣ Maximum Parsimony Analysis (MPboot)

```bash
mkdir -p 07_mpboot
nohup mpirun -n 4 mb < run_mrbayes.sh &
# Or:
# nohup bash mb < run_mrbayes &
```

---

### 8️⃣ Time Dating (MEGA)

> **See the tutorial:**  
> [Constructing a Timetree (ML)](https://www.megasoftware.net/web_help_10/Part_I_Getting_Started/A_Walk_Through_MEGA/Constructing_a_Timetree_(ML).htm)

---

## 📁 Project Folder Structure

```
project/
├── 00_in_house_data/
├── 01_data/              # Public loci data, concatenated with in-house data
├── 02_mafft/             # MAFFT alignments
├── 03_trimal/            # Trimmed alignments
├── 04_modelfinder/       # Concatenated alignments & model info
├── 05_iqtree/            # IQ-TREE analysis
├── 05_raxml_ng/          # RAxML-NG analysis
├── 06_mrbayes/           # MrBayes analysis
├── 07_mpboot/            # Maximum parsimony analysis (MPboot)
├── 08_realtime/          # Time dating
└── work.sh               # All command-lines
```

---

## 🌳 Tree Visualization (ggtree)

> **Troubleshooting:**  
> 1. Do tree labels match the Excel first column exactly?  
> 2. Are there special characters in labels?

### ML Tree (IQ-TREE3)

```r
library(this.path)
library(tidytree)
library(ggtree)
library(openxlsx)
library(treeio)
library(glue)
library(midrootiqtree)
library(export)

setwd(this.dir())
taxa_tbl <- grep('_taxa_table.xlsx$', list.files(), value = TRUE)
df <- read.xlsx(taxa_tbl[ifelse(length(taxa_tbl) == 2, 2, 1)])

outgroup_label <- 'Periconia_pseudodigitata_CBS_139699'
tree_string <- readLines('05_iqtree/iqtree_ml.treefile')
binary_tree_string <- mid_root_iqtree(tree_string, outgroup_label)
cat(binary_tree_string, file = 'iqtree_ml.binary.treefile')
tree_connection <- textConnection(binary_tree_string, open = "r")
tre <- read.iqtree(tree_connection)
close(tree_connection)

ml_tree <- ggtree(tre, size=0.25) %<+% df +
  geom_tiplab(aes(subset = (New == 1 & Type == 1),
                  label=glue("bolditalic({Genus})~bolditalic({Epithet})~bold({Collection})~bold({Number1})")),
              color = 'red', parse = TRUE, size = 2.5) +
  geom_tiplab(aes(subset = (New == 0 & Type == 1),
                  label=glue("bolditalic({Genus})~bolditalic({Epithet})~bold({Collection})~bold({Number1})")),
              color = 'black', parse = TRUE, size = 2.5) +
  geom_tiplab(aes(subset = (New == 0 & Type == 0),
                  label=paste0('italic(', Genus, ')~italic(', Epithet, ')~', Collection, '~', Number1)),
              color = 'black', parse = TRUE, size = 2.5) +
  geom_tiplab(aes(subset = (New == 1 & Type == 0),
                  label=paste0('italic(', Genus, ')~italic(', Epithet, ')~', Collection, '~', Number1)),
              color = 'red', parse = TRUE, size = 2.5) +
  geom_text2(aes(label=UFboot, subset = (UFboot >= 95)), vjust = -0.1, hjust = 1, size = 2) +
  xlim(NA, 0.35) + geom_treescale()
ml_tree
export::graph2ppt(ml_tree, file="ml.pptx", paper="a4", width=8.27, height=10, orient="landscape")
```

### MP Tree (MPboot)

```r
mpboot_tre <- read.iqtree('07_mpboot/mpboot.contree')

ggtree(mpboot_tre, size=0.25) + geom_tiplab() + xlim(NA, 120) +
  geom_text2(aes(label=UFboot, subset = (UFboot >= 95)), vjust = 0, hjust = 1, size = 2)

ggtree(mpboot_tre, size=0.25) %<+% df +
  geom_tiplab(aes(subset = (New == 1 & Type == 1),
                  label=glue("bolditalic({Genus})~bolditalic({Epithet})~bold({Collection})~bold({Number1})")),
              color = 'red', parse = TRUE, size = 2.5) +
  geom_tiplab(aes(subset = (New == 0 & Type == 1),
                  label=glue("bolditalic({Genus})~bolditalic({Epithet})~bold({Collection})~bold({Number1})")),
              color = 'black', parse = TRUE, size = 2.5) +
  geom_tiplab(aes(subset = (New == 0 & Type == 0),
                  label=paste0('italic(', Genus, ')~italic(', Epithet, ')~', Collection, '~', Number1)),
              color = 'black', parse = TRUE, size = 2.5) +
  geom_tiplab(aes(subset = (New == 1 & Type == 0),
                  label=paste0('italic(', Genus, ')~italic(', Epithet, ')~', Collection, '~', Number1)),
              color = 'red', parse = TRUE, size = 2.5) +
  geom_text2(aes(label=UFboot, subset = (UFboot >= 95)), vjust = 0, hjust = 1, size = 2) +
  xlim(NA, 2) + geom_treescale()
```

### Bayesian Tree (MrBayes)

```r
mrbayes_tre <- read.mrbayes('06_mrbayes/run_mrbayes.nexus.con.tre')
ggtree(mrbayes_tre, size=0.25) %<+% df +
  geom_tiplab(aes(subset = (New == 1 & Type == 1),
                  label=glue("bolditalic({Genus})~bolditalic({Epithet})~bold({Collection})~bold({Number1})")),
              color = 'red', parse = TRUE, size = 2.5) +
  geom_tiplab(aes(subset = (New == 0 & Type == 1),
                  label=glue("bolditalic({Genus})~bolditalic({Epithet})~bold({Collection})~bold({Number1})")),
              color = 'black', parse = TRUE, size = 2.5) +
  geom_tiplab(aes(subset = (New == 0 & Type == 0),
                  label=paste0('italic(', Genus, ')~italic(', Epithet, ')~', Collection, '~', Number1)),
              color = 'black', parse = TRUE, size = 2.5) +
  geom_tiplab(aes(subset = (New == 1 & Type == 0),
                  label=paste0('italic(', Genus, ')~italic(', Epithet, ')~', Collection, '~', Number1)),
              color = 'red', parse = TRUE, size = 2.5) +
  geom_text2(aes(label=round(as.numeric(prob), 2), subset = (prob >= 0.90 & !isTip)), vjust = 0, hjust = 1, size = 2) +
  xlim(NA, 0.25) + geom_treescale()
```

---

## 📚 Citation

If you use this pipeline, please cite:

> Yanpeng Chen. (2025). **barcodephylo: A Pipeline for Barcode-based Phylogenetics**. GitHub repository. [https://github.com/yourusername/barcodephylo](https://github.com/yourusername/barcodephylo)

---

## 🧠 Troubleshooting & Tips

- **Outgroup not rooted properly?** Double-check spelling and format of `outgroup_label`.
- **MrBayes hangs?** Try running without MPI or `mpirun`.

