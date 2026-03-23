# TRACMIB Metagenome Analysis Pipeline

Bioinformatics workflow for analysing six sediment shotgun metagenomes. Covers quality control, metagenomic assembly, taxonomic classification, and functional annotation to characterise microbial community composition and metabolic potential.

---

## Repository Structure

```
TRACMIB/
├── 00_LOGS/                  # Log files from all pipeline steps
├── 01_RAW/                   # Raw sequencing reads (.fastq.gz)
│   └── 01_QUALITY/           # FastQC + MultiQC reports (raw)
├── 02_TRIMMED/               # Quality-trimmed reads
│   ├── 01_QUALITY/           # fastp QC reports
│   └── 02_QUALITY/           # FastQC + MultiQC reports (trimmed)
├── 03_ASSEMBLY/              # Metagenomic assemblies
│   ├── 01_QUALITY/           # QUAST assembly QC
│   └── 02_CONTIGS/           # Final contig files
├── 04_TAXONOMY/              # Taxonomic classification results
│   ├── KRAKEN2_READ/         # Read-based Kraken2 + Bracken
│   │   ├── 01_OUTPUT/
│   │   ├── 02_REPORT/
│   │   └── 03_BRACKEN/
│   └── KRAKEN2_CONTIG/       # Contig-based Kraken2 + Bracken
│       ├── 01_OUTPUT/
│       ├── 02_REPORT/
│       └── 03_BRACKEN/
├── 05_MAPPING/               # Bowtie2 mapping results
│   ├── 01_INDEX/             # Bowtie2 indices
│   ├── 02_BAM/               # BAM files (all contigs)
│   ├── 03_BAM_FILT/          # BAM files (contigs >1000 bp)
│   └── RESFINDER/            # ARG mapping BAM files
├── 06_AUTHENTICATION/        # aDNA authentication
│   └── pyDamage/             # pyDamage results per sample
└── 999_SCRIPTS/              # All scripts used in the pipeline
```

---

## Setup

### Create directory structure

```bash
mkdir -p 00_LOGS \
         01_RAW/01_QUALITY \
         02_TRIMMED/01_QUALITY \
         02_TRIMMED/02_QUALITY \
         03_ASSEMBLY/01_QUALITY \
         03_ASSEMBLY/02_CONTIGS \
         04_TAXONOMY/KRAKEN2_READ/01_OUTPUT \
         04_TAXONOMY/KRAKEN2_READ/02_REPORT \
         04_TAXONOMY/KRAKEN2_READ/03_BRACKEN \
         04_TAXONOMY/KRAKEN2_CONTIG/01_OUTPUT \
         04_TAXONOMY/KRAKEN2_CONTIG/02_REPORT \
         04_TAXONOMY/KRAKEN2_CONTIG/03_BRACKEN \
         05_MAPPING/01_INDEX \
         05_MAPPING/02_BAM \
         05_MAPPING/03_BAM_FILT \
         05_MAPPING/RESFINDER \
         06_AUTHENTICATION/pyDamage \
         999_SCRIPTS
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

### Sample list

Create `sample_names.txt` (used by all array jobs):

```
Tracmib_1
Tracmib_5
Tracmib_10
Tracmib_15
Tracmib_25
Tracmib_30
```

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

## Pipeline

### Step 1 — Quality Control: Raw reads

**Tools:** FastQC v0.12.1, MultiQC v1.30

```bash
#!/bin/bash -l
#SBATCH --job-name=fastqc_raw
#SBATCH --account=project_2014298
#SBATCH --output=00_LOGS/fastqc_raw_out_%j.txt
#SBATCH --error=00_LOGS/fastqc_raw_err_%j.txt
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

module load biokit/11.3.0

fastqc -o 01_RAW/01_QUALITY/ -t 8 01_RAW/*.fastq.gz
```

Run MultiQC interactively after FastQC completes:

```bash
sinteractive --time=01:00:00 --mem=8G --cpus-per-task=4 --account=project_2014298
module load multiqc/1.30
cd 01_RAW/01_QUALITY/
multiqc . -p
```

---

### Step 2 — Trimming: fastp

**Tool:** fastp v0.24.0

```bash
#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --account=project_2014298
#SBATCH --output=00_LOGS/fastp_out_%A_%a.txt
#SBATCH --error=00_LOGS/fastp_err_%A_%a.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=15G
#SBATCH --time=02:00:00
#SBATCH --array=1-6

cd /scratch/project_2014298/tarumarj/TRACMIB

module load fastp/0.24.0

i=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_names.txt)

fastp \
    --in1 01_RAW/${i}_1.fastq.gz \
    --in2 01_RAW/${i}_2.fastq.gz \
    --out1 02_TRIMMED/${i}_R1.fastq.gz \
    --out2 02_TRIMMED/${i}_R2.fastq.gz \
    --trim_poly_g --detect_adapter_for_pe \
    --json 02_TRIMMED/01_QUALITY/${i}.fastp.json \
    --html 02_TRIMMED/01_QUALITY/${i}.fastp.html \
    --thread ${SLURM_CPUS_PER_TASK}
```

---

### Step 3 — Quality Control: Trimmed reads

**Tools:** FastQC v0.12.1, MultiQC v1.30

```bash
#!/bin/bash -l
#SBATCH --job-name=fastqc_trimmed
#SBATCH --account=project_2014298
#SBATCH --output=00_LOGS/fastqc_trimmed_out_%j.txt
#SBATCH --error=00_LOGS/fastqc_trimmed_err_%j.txt
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

module load biokit/11.3.0

fastqc -o 02_TRIMMED/02_QUALITY/ -t 8 02_TRIMMED/*.fastq.gz
```

Run MultiQC interactively after FastQC completes:

```bash
sinteractive --time=01:00:00 --mem=8G --cpus-per-task=4 --account=project_2014298
module load multiqc/1.30
cd 02_TRIMMED/02_QUALITY/
multiqc . -p
```

---

### Step 4 — Assembly: MEGAHIT

**Tool:** MEGAHIT v1.2.9  
Minimum contig length: 200 bp. Uses unmerged paired-end reads.

```bash
#!/bin/bash -l
#SBATCH --job-name=megahit
#SBATCH --account=project_2014298
#SBATCH --output=00_LOGS/megahit_out_%A_%a.txt
#SBATCH --error=00_LOGS/megahit_err_%A_%a.txt
#SBATCH --time=06:00:00
#SBATCH --mem=30G
#SBATCH --cpus-per-task=14
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-6
#SBATCH --gres=nvme:500

cd /scratch/project_2014298/tarumarj/TRACMIB

module purge
module load megahit/1.2.9

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_names.txt)

megahit \
    -1 02_TRIMMED/${SAMPLE}_R1.fastq.gz \
    -2 02_TRIMMED/${SAMPLE}_R2.fastq.gz \
    -o 03_ASSEMBLY/${SAMPLE} \
    -t 14 --min-contig-len 200 \
    --tmp-dir $LOCAL_SCRATCH \
    --out-prefix ${SAMPLE}
```

---

### Step 5 — Assembly QC: QUAST

**Tool:** QUAST v5.3.0

```bash
#!/bin/bash -l
#SBATCH --job-name=quast
#SBATCH --account=project_2014298
#SBATCH --output=00_LOGS/quast_out_%j.txt
#SBATCH --error=00_LOGS/quast_err_%j.txt
#SBATCH --time=01:00:00
#SBATCH --mem=10G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

cd /scratch/project_2014298/tarumarj/TRACMIB

module purge
module load quast/5.3.0

metaquast.py -o 03_ASSEMBLY/01_QUALITY/ \
    --max-ref-num 0 --fast --threads 4 \
    03_ASSEMBLY/02_CONTIGS/*.contigs.fa
```

---

### Step 6 — aDNA Authentication: pyDamage

Read-back mapping to assemblies to assess DNA damage patterns. Only contigs >1000 bp are retained. Results filtered at q-value ≤ 0.05 and predicted accuracy ≥ 0.5.

**Tools:** Bowtie2 v2.5.3, SAMtools v1.21, pyDamage

#### 6a. Map reads to filtered contigs (>1000 bp)

```bash
#!/bin/bash -l
#SBATCH --job-name=bowtie2_pydamage
#SBATCH --account=project_2014298
#SBATCH --output=00_LOGS/bowtie2_pydamage_out_%A_%a.txt
#SBATCH --error=00_LOGS/bowtie2_pydamage_err_%A_%a.txt
#SBATCH --time=08:00:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --array=1-6

cd /scratch/project_2014298/tarumarj/TRACMIB

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_names.txt)

module purge
module load bowtie2/2.5.3
module load samtools/1.21

# Filter contigs to >1000 bp
samtools faidx 03_ASSEMBLY/02_CONTIGS/${SAMPLE}.contigs.fa
awk '$2 > 1000 {print $1}' 03_ASSEMBLY/02_CONTIGS/${SAMPLE}.contigs.fa.fai \
    > 03_ASSEMBLY/02_CONTIGS/${SAMPLE}.contigs_gt1000.list
samtools faidx 03_ASSEMBLY/02_CONTIGS/${SAMPLE}.contigs.fa \
    -r 03_ASSEMBLY/02_CONTIGS/${SAMPLE}.contigs_gt1000.list \
    > 03_ASSEMBLY/02_CONTIGS/${SAMPLE}.contigs_gt1000.fa

# Build index and map reads
bowtie2-build --threads $SLURM_CPUS_PER_TASK \
    03_ASSEMBLY/02_CONTIGS/${SAMPLE}.contigs_gt1000.fa \
    05_MAPPING/01_INDEX/${SAMPLE}_gt1000_index

bowtie2 -x 05_MAPPING/01_INDEX/${SAMPLE}_gt1000_index \
    -1 02_TRIMMED/${SAMPLE}_R1.fastq.gz \
    -2 02_TRIMMED/${SAMPLE}_R2.fastq.gz \
    --very-sensitive -N 1 \
    --threads $SLURM_CPUS_PER_TASK \
| samtools view -Sb - \
| samtools sort -o 05_MAPPING/03_BAM_FILT/${SAMPLE}.sorted.bam

# Add MD tags and index
samtools calmd -b 05_MAPPING/03_BAM_FILT/${SAMPLE}.sorted.bam \
    03_ASSEMBLY/02_CONTIGS/${SAMPLE}.contigs_gt1000.fa \
    > 05_MAPPING/03_BAM_FILT/${SAMPLE}.sorted.calmd.bam

samtools index 05_MAPPING/03_BAM_FILT/${SAMPLE}.sorted.calmd.bam
```

#### 6b. Run pyDamage

```bash
#!/bin/bash
#SBATCH --job-name=pyDamage
#SBATCH --account=project_2014298
#SBATCH --time=02:00:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --error=00_LOGS/pydamage_err_%A_%a.txt
#SBATCH --output=00_LOGS/pydamage_out_%A_%a.txt
#SBATCH --array=1-6

cd /scratch/project_2014298/tarumarj/TRACMIB

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_names.txt)

export PATH="/projappl/project_2007691/pydamage/bin:$PATH"

pydamage --outdir 06_AUTHENTICATION/pyDamage/${SAMPLE} \
    analyze -w 30 \
    -p $SLURM_CPUS_PER_TASK \
    05_MAPPING/03_BAM_FILT/${SAMPLE}.sorted.calmd.bam
```

#### 6c. Filter pyDamage results

```bash
#!/bin/bash
#SBATCH --job-name=pyDamage_filter
#SBATCH --account=project_2014298
#SBATCH --time=00:15:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --error=00_LOGS/pydamage_filter_err_%A_%a.txt
#SBATCH --output=00_LOGS/pydamage_filter_out_%A_%a.txt
#SBATCH --array=1-6

cd /scratch/project_2014298/tarumarj/TRACMIB

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_names.txt)

export PATH="/projappl/project_2007691/pydamage/bin:$PATH"

pydamage --outdir 06_AUTHENTICATION/pyDamage/${SAMPLE} \
    filter 06_AUTHENTICATION/pyDamage/${SAMPLE}/pydamage_results.csv
```

---

### Step 7 — Taxonomic Classification: Kraken2 + Bracken

**Tools:** Kraken2 v2.17.0, Bracken v3.1  
**Database:** db_vistamilk_soil_metagenomic_kraken2_2024

#### 7a. Read-based classification

```bash
#!/bin/bash -l
#SBATCH --job-name=kraken2_reads
#SBATCH --output=00_LOGS/kraken2_reads_out_%A_%a.txt
#SBATCH --error=00_LOGS/kraken2_reads_err_%A_%a.txt
#SBATCH --time=04:00:00
#SBATCH --partition=hugemem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --account=project_2014298
#SBATCH --mem=370G
#SBATCH --array=1-6

module load kraken/2.17.0

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_names.txt)

kraken2 --db /scratch/project_2007691/00_DATABASES/db_vistamilk_soil_metagenomic_kraken2_2024 \
    --threads $SLURM_CPUS_PER_TASK \
    --paired \
    02_TRIMMED/${SAMPLE}_R1.fastq.gz \
    02_TRIMMED/${SAMPLE}_R2.fastq.gz \
    --output 04_TAXONOMY/KRAKEN2_READ/01_OUTPUT/${SAMPLE}.out \
    --report 04_TAXONOMY/KRAKEN2_READ/02_REPORT/${SAMPLE}.report \
    --use-names --report-minimizer-data
```

#### 7b. Contig-based classification

```bash
#!/bin/bash -l
#SBATCH --job-name=kraken2_contigs
#SBATCH --output=00_LOGS/kraken2_contigs_out_%A_%a.txt
#SBATCH --error=00_LOGS/kraken2_contigs_err_%A_%a.txt
#SBATCH --time=04:00:00
#SBATCH --partition=hugemem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --account=project_2014298
#SBATCH --mem=370G
#SBATCH --array=1-6

module load kraken/2.17.0

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_names.txt)

kraken2 --db /scratch/project_2007691/00_DATABASES/db_vistamilk_soil_metagenomic_kraken2_2024 \
    --threads $SLURM_CPUS_PER_TASK \
    --output 04_TAXONOMY/KRAKEN2_CONTIG/01_OUTPUT/${SAMPLE}.out \
    --report 04_TAXONOMY/KRAKEN2_CONTIG/02_REPORT/${SAMPLE}.report \
    --use-names --report-minimizer-data \
    03_ASSEMBLY/02_CONTIGS/${SAMPLE}.contigs.fa
```

#### 7c. Bracken abundance estimation

```bash
#!/bin/bash -l
#SBATCH --job-name=bracken
#SBATCH --account=project_2014298
#SBATCH --output=00_LOGS/bracken_out_%A_%a.txt
#SBATCH --error=00_LOGS/bracken_err_%A_%a.txt
#SBATCH --time=02:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --array=1-6

cd /scratch/project_2014298/tarumarj/TRACMIB

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_names.txt)

module load bracken/3.1

# Read-based
bracken -d /scratch/project_2007691/00_DATABASES/db_vistamilk_soil_metagenomic_kraken2_2024 \
    -i 04_TAXONOMY/KRAKEN2_READ/02_REPORT/${SAMPLE}.report \
    -o 04_TAXONOMY/KRAKEN2_READ/03_BRACKEN/${SAMPLE}.bracken.txt \
    -w 04_TAXONOMY/KRAKEN2_READ/03_BRACKEN/${SAMPLE}.breport.txt \
    -r 150 -l S

# Contig-based
bracken -d /scratch/project_2007691/00_DATABASES/db_vistamilk_soil_metagenomic_kraken2_2024 \
    -i 04_TAXONOMY/KRAKEN2_CONTIG/02_REPORT/${SAMPLE}.report \
    -o 04_TAXONOMY/KRAKEN2_CONTIG/03_BRACKEN/${SAMPLE}.bracken.txt \
    -w 04_TAXONOMY/KRAKEN2_CONTIG/03_BRACKEN/${SAMPLE}.breport.txt \
    -r 150 -l S
```

#### 7d. Combine Bracken reports (interactive)

```bash
# Convert to MPA format
for file in 04_TAXONOMY/KRAKEN2_READ/03_BRACKEN/*.breport.txt; do
    SAMPLE=$(basename "$file" .breport.txt)
    04_TAXONOMY/KrakenTools/kreport2mpa.py \
        -r "$file" --display-header \
        -o 04_TAXONOMY/KRAKEN2_READ/03_BRACKEN/${SAMPLE}_mpa.txt
done

# Combine all samples
04_TAXONOMY/KrakenTools/combine_mpa.py \
    -i 04_TAXONOMY/KRAKEN2_READ/03_BRACKEN/*_mpa.txt \
    -o 04_TAXONOMY/KRAKEN2_READ/03_BRACKEN/combined_mpa.txt

# Extract species level
awk '$1 ~ "clade_name" || $1 ~ "s__" {print $0}' \
    04_TAXONOMY/KRAKEN2_READ/03_BRACKEN/combined_mpa.txt \
    | grep -v "t__" > 04_TAXONOMY/KRAKEN2_READ/03_BRACKEN/mpa_species.txt
```

---

### Step 8 — ARG Mapping: ResFinder

Map trimmed reads against the ResFinder antibiotic resistance gene database.

**Tools:** Bowtie2 v2.5.3, SAMtools v1.21  
**Database:** ResFinder v2.6.0

```bash
#!/bin/bash -l
#SBATCH --job-name=bowtie2_ARG
#SBATCH --account=project_2014298
#SBATCH --output=00_LOGS/bowtie2_ARG_out_%A_%a.txt
#SBATCH --error=00_LOGS/bowtie2_ARG_err_%A_%a.txt
#SBATCH --time=12:00:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=15G
#SBATCH --array=1-6

cd /scratch/project_2014298/tarumarj/TRACMIB

name=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_names.txt)

module purge
module load bowtie2/2.5.3
module load samtools/1.21

bowtie2 -x /scratch/project_2014298/00_DATABASES/RESFINDER_2.6.0/resfinder99 \
    -1 02_TRIMMED/${name}_R1.fastq.gz \
    -2 02_TRIMMED/${name}_R2.fastq.gz \
    -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 \
    --threads "$SLURM_CPUS_PER_TASK" \
| samtools view -Sb - \
| samtools sort -o 05_MAPPING/RESFINDER/${name}.bam

samtools index 05_MAPPING/RESFINDER/${name}.bam
```

Build ARG count table (run interactively from `05_MAPPING/RESFINDER/`):

```bash
for bam in *.bam; do
    name=$(basename "$bam" .bam)
    echo -e "$name" > "${name}_counts"
    samtools idxstats "$bam" | grep -v "*" | cut -f3 >> "${name}_counts"
done

echo -e "gene" > gene_names
samtools idxstats Tracmib_1.bam | grep -v "*" | cut -f1 >> gene_names

paste gene_names *_counts > ARG_counts.txt
```

---

## Computing Environment

Analyses run on [CSC Puhti HPC cluster](https://docs.csc.fi/computing/systems-puhti/).  
Data stored on [Allas object storage](https://docs.csc.fi/data/Allas/).
