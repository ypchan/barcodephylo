# 🧬 barcodephylo

**`barcodephylo`** is a bioinformatics pipeline for phylogenetic analysis based on barcode sequences (e.g., ITS, LSU, SSU, RPB2, TUB2). It offers a modular and reproducible workflow from sequence collection to tree construction using tools such as **MAFFT**, **trimal**, **IQ-TREE**, and **MrBayes**.

![WORKFLOW](https://img.shields.io/badge/WORKFLOW-8A2BE2)

1. QC in-house loci data genereated from Sanger sequencing
2. Downloading public data using read.GenBank function in APE package
3. Multiple sequence alignment using mafft
4. To trimm msa files
5. To determine the best evolution model and concatenate msa
6. ML analysis using iqtree/raxml-ng
7. MrBayes analysis
8. Minimum Parsimony analysis using mpboot
9. Time dating using mega
10. Vissualization using ggtree + Adobe Illustrator

---

## 📦 Installation

### 🔧 Prerequisites

- Python ≥ 3.6  
- R ≥ 4.0  
- MAFFT, trimal, IQ-TREE2, MrBayes, MPI (if needed)  
- R packages: `ape`, `seqinr`, etc.

### 🐍 Using `pip`

```bash
git clone https://github.com/yourusername/barcodephylo.git
cd barcodephylo
pip install .
```

---

## 🚀 Quick Start

### 1️⃣ Download Barcode Sequences

Download barcodes using the `read.GenBank.R` script. You may also append your Sanger sequences:

```bash
cat 01_data/00_sanger_raw_data/2021032807_ITS.fasta >> 01_data/Boeremia_ITS.fasta
cat 01_data/00_sanger_raw_data/2021032807_LROR_LR5.fasta | sed 's/_LROR_LR5//' >> 01_data/Boeremia_LSU.fasta
cat 01_data/00_sanger_raw_data/2021032807_PNS1_NS41.fasta | sed 's/_PNS1_NS41//' >> 01_data/Boeremia_SSU.fasta
cat 01_data/00_sanger_raw_data/2021032807_RPB2.fas | sed 's/_RPB2//' >> 01_data/Boeremia_RPB2.fasta
cat 01_data/00_sanger_raw_data/2021032807_TUB2.fas | sed 's/_TUB2//' >> 01_data/Boeremia_TUB2.fasta
```

---

### 2️⃣ Multiple Sequence Alignment with MAFFT

```bash
is_exist_folder() {
    folder_name=$1
    [[ ! -d ${folder_name} ]] && mkdir ${folder_name} && echo "${folder_name} created" || echo "${folder_name} already exists"
}

is_exist_folder 02_mafft

ls 01_data/ | grep fasta | while read a; do
  echo "mafft --localpair --thread 4 --adjustdirection 01_data/${a} > 02_mafft/${a%.fasta}.mafft.fna"
done > mafft.sh

bash mafft.sh
sed -i 's/>_R_/>/' 02_mafft/*.mafft.fna
```

---

### 3️⃣ Alignment Trimming with trimal

```bash
is_exist_folder 03_trimal

ls 01_data/ | grep fasta | sed 's/.fasta//' | while read a; do
  trimal -in 02_mafft/${a}.mafft.fna -gt 0.5 -out 03_trimal/${a}.mafft.trimal.fna
done
```

---

### 4️⃣ Concatenate MSAs & Select Models

```bash
is_exist_folder 04_modelfinder

outgroup_label=Phoma_herbarum_CBS_615.75
mafft_items=$(ls 03_trimal/*)

iqtree_modelfinder.py -i ${mafft_items} -o 04_modelfinder --mrbayes_nexus --outgroup ${outgroup_label}
```

---

### 5️⃣ Build Maximum Likelihood Tree with IQ-TREE

```bash
mkdir -p 05_iqtree

nohup iqtree2 -s 04_modelfinder/concatenated.fna \
            --seqtype DNA -o ${outgroup_label} \
            --prefix 05_iqtree/iqtree_ml -T AUTO \
            -p 04_modelfinder/best_scheme.txt \
            --ufboot 1000 --alrt 1000 &
```

---

### 6️⃣ Bayesian Phylogenetics with MrBayes

```bash
mkdir -p 06_mrbayes

nohup mpirun -n 4 mb < run_mrbayes.sh &  # Or:
# nohup bash mb < run_mrbayes &
```

---

## 📁 Project Folder
For a certain phylogenetic analysis project, the following is a example folder.

```
project/
├── 00_in_house_data/
├── 01_data/              # Public loci dada and concatenating them to in house dada
├── 02_mafft/             # MAFFT alignments
├── 03_trimal/            # Trimmed alignments
├── 04_modelfinder/       # Concatenated alignments and model info
├── 05_iqtree/            # IQ-TREE analysis
├── 05_raxml_ng/          # raxml-ng analysis
├── 06_mrbayes/           # MrBayes analysis
├── 07_mpboot/            # manimum parsimony analysis uing mpboot
├── 08_realtime/          # for time dating
└── work.sh               # all command-lines
```

---

## 📚 Citation

If you use this pipeline, please consider citing it:

> Yanpeng Chen. (2025). **barcodephylo: A Pipeline for Barcode-based Phylogenetics**. GitHub repository. [https://github.com/yourusername/barcodephylo](https://github.com/yourusername/barcodephylo)

---

## 🧠 Troubleshooting & Tips

- **Outgroup not rooted properly?** Check spelling and format of `outgroup_label`.
- **MrBayes hangs?** Try removing MPI or run without `mpirun`.