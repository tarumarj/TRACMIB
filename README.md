# TRACMIB Metagenome Analysis Pipeline

Bioinformatics workflow for analysing six sediment shotgun metagenomes. Covers quality control, metagenomic assembly, taxonomic classification, and functional annotation to characterise microbial community composition and metabolic potential.

---

## Main Repository Structure

```
TRACMIB/
├── 00_LOGS/          # Log files from all pipeline steps
├── 01_RAW/           # Raw sequencing reads (.fastq.gz)
├── 02_TRIMMED/       # Quality-trimmed reads
├── 03_ASSEMBLY/      # Metagenomic assemblies
├── 04_TAXONOMY/      # Taxonomic classification results
├── 05_ANNOTATION/    # Functional annotation results
└── 999_SCRIPTS/      # All scripts used in the pipeline
```

---

## Setup

### Create directory structure

```bash
mkdir -p 00_LOGS 01_RAW 02_TRIMMED 03_ASSEMBLY 04_TAXONOMY 05_ANNOTATION 999_SCRIPTS
```

### Rename raw reads

Rename files from sequencing center format to sample names:

```bash
cd 01_RAW/
for f in NG-A5412_Tracmib_*.fastq.gz; do
    newname=$(echo "$f" | sed 's/NG-A5412_\(Tracmib_[0-9]*\)_libLAP[0-9]*_\([12]\)/\1_\2/')
    mv "$f" "$newname"
done
```

**Result:** `NG-A5412_Tracmib_15_libLAP8065_1.fastq.gz` → `Tracmib_15_1.fastq.gz`

---

## Samples

| Sample | Forward read | Reverse read |
|--------|-------------|-------------|
| Tracmib_1 | Tracmib_1_1.fastq.gz | Tracmib_1_2.fastq.gz |
| Tracmib_5 | Tracmib_5_1.fastq.gz | Tracmib_5_2.fastq.gz |
| Tracmib_10 | Tracmib_10_1.fastq.gz | Tracmib_10_2.fastq.gz |
| Tracmib_15 | Tracmib_15_1.fastq.gz | Tracmib_15_2.fastq.gz |
| Tracmib_25 | Tracmib_25_1.fastq.gz | Tracmib_25_2.fastq.gz |
| Tracmib_30 | Tracmib_30_1.fastq.gz | Tracmib_30_2.fastq.gz |

---

## Pipeline Steps

### 1. Quality Control
Tools: FastQC, MultiQC, Fastp
- Assess raw read quality
- Trim adapters and low-quality bases
- Output: `02_TRIMMED/`

### 2. Metagenomic Assembly
Tools: MEGAHIT
- Assemble trimmed reads into contigs
- Output: `03_ASSEMBLY/`

### 3. Taxonomic Classification
Tools: Kraken2
- Classify reads/contigs against reference database
- Output: `04_TAXONOMY/`

### 4. Functional Annotation
- Predict and annotate genes on assembled contigs
- Output: `05_ANNOTATION/`

---

## Computing Environment

Analyses run on [CSC Puhti HPC cluster](https://docs.csc.fi/computing/systems-puhti/).
Data stored on [Allas object storage](https://docs.csc.fi/data/Allas/).
