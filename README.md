# Enterobacterales-chromosomal-and-plasmid-clustering-pipeline
This repository offers description of the analytical steps needed for plasmid-encoded resistance dynamics investigation using Oxford Nanopore reads. This README outlines what to do and which tools to use. There are also automated softwares that provide similar analyses for non-experts, but the tools and parameters that they use are usually not fully visible. 

---

## 1. Goal & Scope

Identify and compare plasmids (and chromosomes) from relevant bacterial isolates sequenced with the RBK114 library preparation kit on a MinION flowcell:

- Reconstruct plasmids and annotate resistance/replicon markers
- Assess relatedness (SNP distances, plasmid clustering)
---

## 2. Workflow Overview

1. **Signal → Bases**: Basecall raw POD5 signals with Dorado
2. **Per‑sample separation**: Demultiplex barcoded reads (Dorado demux)
3. **Read cleaning**: Adapter trimming & QC filtering (*optional if basecaller model handles it*)
4. **Assembly**: De novo assembly (Flye)
5. **Consensus polishing**: Map & correct (Minimap2 + Medaka v2.0.1 with `--bacteria` flag)
6. **Plasmid reconstruction & typing**: MOB‑Suite (mob\_recon)
7. **AMR gene detection**: AMRFinderPlus
8. **Chromosomal relatedness**: SNP calling & core alignment (Snippy + snp-dists)
9. **Plasmid relatedness**: Mash first, then cluster (e.g. DCJ distance via Pling)
10. **Miscellaneous**: Manual SNP/Mash matrices & visualizations; species ID via Pathogenwatch or PubMLST

---

## 3. Inputs

- Raw POD5 or basecalled FASTQ (e.g., from MinKNOW)
- Barcode list (e.g. `01–24`)

---

## 4. Detailed Steps & Codes

> Replace placeholder paths, models, parameters as appropriate.

### 4.1 Basecalling

**Tool:** Dorado

```bash
DORADO_BIN="/path/to/dorado"
CONFIG_FILE="/path/to/dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
INPUT_DIR="/path/to/pod5_dir"
OUTDIR="/path/to/output"
KIT="SQK-RBK114-24"

mkdir -p "$OUTDIR"
"$DORADO_BIN" basecaller "$CONFIG_FILE" "$INPUT_DIR" \
  --emit-fastq --kit-name "$KIT" -r --no-trim \
  > "$OUTDIR/basecalled.fastq"
```

### 4.2 Demultiplexing

**Tool:** Dorado demux

```bash
"$DORADO_BIN" demux \
  --output-dir "$OUTDIR/demux" \
  --emit-fastq --kit-name "$KIT" \
  "$OUTDIR/basecalled.fastq"
# yields files like demux/barcode01.fastq
```

### 4.3 Trimming & Filtering

**Tools:** Porechop, NanoFilt, Seqkit

```bash
BCLIST=$(seq -w 01 24)
RAW_PREFIX="$OUTDIR/demux/barcode"
PC_DIR="$OUTDIR/porechop"
NF_DIR="$OUTDIR/filtered_reads"
STATS_DIR="$OUTDIR/stats"
QMIN=8; LMIN=200
mkdir -p "$PC_DIR" "$NF_DIR" "$STATS_DIR"

for bc in $BCLIST; do
  in="${RAW_PREFIX}${bc}.fastq"; [[ -s "$in" ]] || continue
  trim="$PC_DIR/trimmed_${bc}.fastq"
  out="$NF_DIR/filtered_barcode${bc}_passed.fastq"
  porechop -i "$in" -o "$trim"
  NanoFilt -q "$QMIN" -l "$LMIN" < "$trim" > "$out"
  rm -f "$trim"
done

for fq in "$NF_DIR"/*.fastq; do
  [[ -e "$fq" ]] || continue
  seqkit stats -T -a "$fq" > "$STATS_DIR/$(basename "${fq%.fastq}").txt"
done
```

### 4.4 Assembly

**Tool:** Flye

```bash
GENOME_SIZE="5m"; THREADS=8
FLY_DIR="$OUTDIR/flye"; mkdir -p "$FLY_DIR"
for bc in $BCLIST; do
  reads="$NF_DIR/filtered_barcode${bc}_passed.fastq"; [[ -s "$reads" ]] || continue
  flye --nano-hq "$reads" --out-dir "$FLY_DIR/barcode${bc}" \
       --genome-size "$GENOME_SIZE" --threads "$THREADS"
done
```

### 4.5 Consensus Polishing

**Tools:** Minimap2, Samtools, Medaka v2.0.1 (`--bacteria`)

```bash
MIN_DIR="$OUTDIR/minimap"; MEDAKA_DIR="$OUTDIR/medaka"
mkdir -p "$MIN_DIR" "$MEDAKA_DIR"
THREADS=8; MODEL="r1041_min_sup_g632"
for bc in $BCLIST; do
  asm="$FLY_DIR/barcode${bc}/assembly.fasta"
  reads="$NF_DIR/filtered_barcode${bc}_passed.fastq"
  [[ -s "$asm" && -s "$reads" ]] || continue
  bam="$MIN_DIR/aligned_barcode${bc}.bam"
  med_out="$MEDAKA_DIR/barcode${bc}"; mkdir -p "$med_out"

  minimap2 -t "$THREADS" -ax map-ont "$asm" "$reads" \
    | samtools view -b - \
    | samtools sort -@ "$THREADS" -o "$bam" -
  samtools index "$bam"
  medaka_consensus -b "$bam" -d "$asm" -o "$med_out" \
    -t "$THREADS" -m "$MODEL" --bacteria

done
```

### 4.6 Plasmid Reconstruction & Typing

**Tool:** MOB‑Suite

```bash
MOB_DIR="$OUTDIR/mob"; mkdir -p "$MOB_DIR"
for bc in $BCLIST; do
  cns="$MEDAKA_DIR/barcode${bc}/consensus.fasta"; [[ -s "$cns" ]] || continue
  mob_recon --infile "$cns" --outdir "$MOB_DIR/barcode${bc}"
done
```

### 4.7 AMR Gene Detection

**Tool:** AMRFinderPlus

```bash
AMR_DIR="$OUTDIR/amr"; mkdir -p "$AMR_DIR"
for bc in $BCLIST; do
  cns="$MEDAKA_DIR/barcode${bc}/consensus.fasta"; [[ -s "$cns" ]] || continue
  amrfinder -n "$cns" -o "$AMR_DIR/barcode${bc}.amr.txt" \
    --threads "$THREADS" --plus
done
```

### 4.8 Chromosomal SNP Analysis

**Tools:** Snippy, snippy-core, snp-dists

```bash
REFERENCE="/path/to/reference.fasta"
SNIPPY_DIR="$OUTDIR/snippy"; mkdir -p "$SNIPPY_DIR"
for SAMPLE in "$MEDAKA_DIR"/barcode*/consensus.fasta; do
  [[ -s "$SAMPLE" ]] || continue
  BN=$(basename "$(dirname "$SAMPLE")")
  snippy --cpus 4 --ref "$REFERENCE" --ctgs "$SAMPLE" --outdir "$SNIPPY_DIR/$BN"
done
snippy-core --ref "$REFERENCE" --prefix snippy_core "$SNIPPY_DIR"/*
snp-dists snippy_core.full.aln > "$SNIPPY_DIR/snp_distances.tsv"
```

### 4.9 Plasmid Clustering / Distances

**Tools:** Mash, Pling

```bash
# 1) Compute Mash distances
PLASMID_TYPE="IncN"
OUT_TXT="$OUTDIR/clustering/mash_${PLASMID_TYPE}.txt"; > "$OUT_TXT"
for i in {1..24}; do
  for j in {1..24}; do
    [[ $i -lt $j ]] || continue
    f1=$(printf "plasmid_%s_barcode%02d.fasta" "$PLASMID_TYPE" "$i")
    f2=$(printf "plasmid_%s_barcode%02d.fasta" "$PLASMID_TYPE" "$j")
    [[ -f "$f1" && -f "$f2" ]] && mash dist "$f1" "$f2" >> "$OUT_TXT" \
      || echo "Missing: $f1 or $f2"
  done

done

# 2) DCJ-based clustering with Pling (after selecting pairs below threshold)
pling align --containment_distance 0.3 --cores 8 --sourmash plasmid.txt "$OUTDIR/clustering/pling_out"
```

### 4.10 Miscellaneous & Reporting

- SNP/Mash distance matrices and visualizations were assembled manually
- Species identification via Pathogenwatch and/or PubMLST (see corresponding publication)
- Archive: command logs, software versions, parameter files, results summaries

---

## 5. Suggested Folder Structure

```text
project/
├─ raw/                     # POD5/FAST5 or raw FASTQs
├─ demux/
├─ filtered_reads/          # post-trim/filter
├─ flye/                    # raw assemblies
├─ medaka/                  # polished assemblies
├─ mob/                     # plasmid reconstruction
├─ amr/                     # AMR gene tables
├─ snippy/                  # chromosomal SNP results
├─ clustering/             # plasmid clustering & distances
└─ reports/                 # summary tables & figures
```

---

## 6. Software Checklist

- Dorado
- Porechop, NanoFilt, Seqkit
- Flye
- Minimap2, Samtools
- Medaka v2.0.1
- MOB-suite
- AMRFinderPlus
- Snippy, snippy-core, snp-dists
- Mash, Pling, Sourmash
