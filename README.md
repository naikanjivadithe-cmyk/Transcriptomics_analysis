# Transcriptomics_analysis
# Transcriptome and lncRNA Discovery Workflow

End-to-end pipeline: raw RNA-seq reads → QC → alignment → transcript assembly → novel lncRNA identification → differential expression → functional and downstream analysis.

> Edit the variables in Section 0 once. Every later command uses them, so no personal paths are hard-coded.

---

## Pipeline overview

```
Raw reads (SRA/ENA)
   │  fastp
Clean reads
   │  HISAT2 + samtools
Sorted BAM
   │  StringTie → merge → gffcompare
Merged transcripts with class codes
   │  keep i / u / x → BED → FASTA → length > 200 nt
Candidate transcripts
   │  ORF prediction + Pfam (RPS-BLAST) → CPC2 → CPAT → FEELnc
   │  intersection of the three predictors
High-confidence lncRNAs
   │  StringTie -e → prepDE → edgeR
Differentially expressed mRNAs and lncRNAs
   │  NONCODE / tRNA / miRNA / Rfam screening, GO/GSEA, correlation, Circos
Final results
```

---

## 0. Setup and configuration

### 0.1 Conda environment (`environment.yml`)

```yaml
name: transcriptome
channels: [conda-forge, bioconda, defaults]
dependencies:
  - sra-tools
  - fastp=0.21
  - hisat2=2.2.1
  - samtools=1.11
  - stringtie=2.1.4
  - gffcompare=0.11
  - bedtools
  - blast
  - python=3.9
  - pip
  - r-base
  - bioconductor-edger
  - r-ggplot2
  - pip:
      - CPAT
```

```bash
conda env create -f environment.yml
conda activate transcriptome
# CPC2, FEELnc, OrfPredictor, prepDE.py and Circos are installed separately (see Section 17)
```

**Why:** Pinning tool versions makes the results reproducible.

### 0.2 Config file (`config/paths.sh`)

```bash
#!/usr/bin/env bash
export PROJECT=$(pwd)
export THREADS=16

# Reference (download from NCBI Genome: FASTA, GTF, GFF)
export REF_FA=$PROJECT/reference/genome.fna
export REF_GFF=$PROJECT/reference/genome.gff
export REF_GTF=$PROJECT/reference/genome.gtf
export HISAT_INDEX=$PROJECT/reference/ht2/genome

# External databases
export NONCODE_FA=$PROJECT/db/NONCODE_species_lncRNA.fa
export CDS_FA=$PROJECT/db/cds_from_genomic.fna      # coding sequences from NCBI
export PFAM_DB=$PROJECT/db/pfam/Pfam
export RFAM_RIBOSWITCH_CM=$PROJECT/db/riboswitch_only.cm

# Output folders
mkdir -p $PROJECT/{raw,clean,qc,bam,stringtie,lncrna,counts,deg,downstream,logs}
```

### 0.3 Sample sheet (`config/samples.tsv`)

```
sample_id    condition    accession
S1           control      SRRxxxxxxx
S2           infected     SRRxxxxxxx
```

---

## 1. Download data

### 1.1 Raw reads from ENA/SRA

```bash
source config/paths.sh
cd $PROJECT/raw

# single file (-nc = "no clobber": do not overwrite existing files, safe to re-run)
wget -nc <ENA_FTP_LINK>

# many files: ENA gives a batch script
chmod +x ena-file-download-selected-files.sh     # +x = make executable
./ena-file-download-selected-files.sh
ls -l                                            # check permissions and sizes
```

### 1.2 Inspect files

```bash
ls *.fastq.gz
head <(zcat sample_1.fastq.gz)       # .gz is compressed, so zcat is needed to read it
gunzip -k -v *.fastq.gz              # -k keep the .gz, -v verbose
head sample_1.fastq
```

**FASTQ format:** 4 lines per read (ID, sequence, `+`, quality). Paired-end data has `_1` (forward) and `_2` (reverse) files.

### 1.3 Reference genome and annotation

Download the genome FASTA (`.fna`), GTF and GFF for your species from NCBI Genome into `reference/`.

- **GTF** = Gene Transfer Format
- **GFF** = General Feature Format

---

## 2. Quality control and trimming (fastp)

```bash
fastp \
  -i raw/${S}_1.fastq -I raw/${S}_2.fastq \
  -o clean/${S}_1_clean.fastq -O clean/${S}_2_clean.fastq \
  -h qc/${S}.html -j qc/${S}.json -w 8
```

| Flag | Meaning |
|---|---|
| `-i / -I` | read 1 / read 2 input |
| `-o / -O` | read 1 / read 2 output |
| `-h / -j` | HTML and JSON reports |

**Why:** removes adapters and low-quality bases, and corrects errors, so alignment is accurate.

**Batch loop** (`scripts/02_qc_trim.sh`):

```bash
#!/usr/bin/env bash
source config/paths.sh
tail -n +2 config/samples.tsv | while read S COND ACC; do
  fastp -i raw/${ACC}_1.fastq -I raw/${ACC}_2.fastq \
        -o clean/${ACC}_1_clean.fastq -O clean/${ACC}_2_clean.fastq \
        -h qc/${ACC}.html -j qc/${ACC}.json -w 8
done
```

---

## 3. Alignment (HISAT2) and BAM conversion (samtools)

### 3.1 Build genome index (once)

```bash
mkdir -p reference/ht2
hisat2-build -p $THREADS $REF_FA $HISAT_INDEX
```

### 3.2 Map reads

```bash
hisat2 -x $HISAT_INDEX \
  -1 clean/${ACC}_1_clean.fastq -2 clean/${ACC}_2_clean.fastq \
  -S bam/${ACC}.sam --summary-file bam/${ACC}.summary -p $THREADS
```

### 3.3 SAM → sorted BAM

```bash
samtools sort -@ 10 -O BAM -o bam/${ACC}.bam bam/${ACC}.sam
samtools index bam/${ACC}.bam
rm bam/${ACC}.sam          # SAM is huge, delete after conversion
```

**Why:** the SAM is text and very large. A sorted, indexed BAM is compact and required by StringTie.

Check `*.summary` for the overall alignment rate (aim above 80–90%).

---

## 4. Transcript assembly, merging and annotation comparison

### 4.1 Assemble each sample

```bash
stringtie bam/${ACC}.bam -G $REF_GFF -o stringtie/${ACC}.gtf -p 10 -A stringtie/${ACC}_abund.out
```

| Flag | Meaning |
|---|---|
| `-G` | reference annotation to guide assembly |
| `-o` | assembled transcripts (GTF) |
| `-p` | threads |
| `-A` | gene abundance table |

### 4.2 Merge all samples into one non-redundant set

```bash
ls stringtie/*[!_new].gtf > gtf_list.txt        # list of per-sample GTFs
stringtie --merge -G $REF_GFF -o stringtie/all_merged.gtf gtf_list.txt
```

### 4.3 Compare with reference and assign class codes

```bash
gffcompare -r $REF_GFF -s $REF_FA -T -V -i gtf_list.txt
```

Output: `gffcmp.combined.gtf` with a `class_code` for each transcript.

| Code | Meaning |
|---|---|
| `i` | fully within an intron |
| `u` | unknown, intergenic |
| `x` | exonic overlap on the opposite strand (antisense) |

---

## 5. Select candidate lncRNA classes

```bash
mkdir -p lncrna && cd lncrna

grep -c 'class_code "i"' ../gffcmp.combined.gtf      # count each class
grep -c 'class_code "u"' ../gffcmp.combined.gtf
grep -c 'class_code "x"' ../gffcmp.combined.gtf

grep  'class_code "i"' ../gffcmp.combined.gtf >  iux.gtf   # > creates/overwrites
grep  'class_code "u"' ../gffcmp.combined.gtf >> iux.gtf   # >> appends
grep  'class_code "x"' ../gffcmp.combined.gtf >> iux.gtf

wc -l iux.gtf
```

**Why:** lncRNAs are non-coding and do not overlap known coding exons on the same strand. Classes i, u and x capture novel candidates.

---

## 6. Extract sequences and length filter (> 200 nt)

### 6.1 Chromosome sizes

```bash
samtools faidx $REF_FA
cut -f1,2 ${REF_FA}.fai > chrom_sizes.txt
```

### 6.2 GTF → BED (transcript lines only)

```bash
awk '$3=="transcript"{gsub(/"|;/,"",$10); print $1"\t"($4-1)"\t"$5"\t"$10}' iux.gtf > iux.bed
```

GTF is 1-based and BED is 0-based, hence `$4-1`. Field `$10` is the `transcript_id`.

### 6.3 Add 50 bp flanks, then extract FASTA

```bash
bedtools slop -i iux.bed -g chrom_sizes.txt -b 50 > iux_slop.bed
bedtools getfasta -fi $REF_FA -bed iux_slop.bed -fo iux.fasta -nameOnly
```

> If you want exact transcript boundaries, skip `slop` and use `iux.bed` directly.

### 6.4 Length filter (`scripts/helpers/length_filter.py`)

```python
#!/usr/bin/env python3
"""Keep FASTA records longer than MIN_LEN nucleotides.
Usage: python3 length_filter.py in.fasta out.fasta [min_len]"""
import sys

def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name: yield name, "".join(seq)
                name, seq = line, []
            else:
                seq.append(line)
    if name: yield name, "".join(seq)

inp, out = sys.argv[1], sys.argv[2]
min_len = int(sys.argv[3]) if len(sys.argv) > 3 else 200

with open(out, "w") as o:
    for name, seq in read_fasta(inp):
        if len(seq) > min_len:
            o.write(f"{name}\n{seq}\n")
```

```bash
python3 ../scripts/helpers/length_filter.py iux.fasta length_filter.fasta 200
```

**Why:** the standard lncRNA definition is longer than 200 nt.

### 6.5 Reusable helper: keep FASTA records by ID list (`scripts/helpers/filter_fasta_by_ids.py`)

```python
#!/usr/bin/env python3
"""Keep (or drop) FASTA records whose header is in an ID list.
Usage: python3 filter_fasta_by_ids.py in.fasta ids.txt out.fasta [--exclude]
IDs in ids.txt may be with or without a leading '>'."""
import sys

fasta, idfile, out = sys.argv[1:4]
exclude = "--exclude" in sys.argv

ids = {l.strip().lstrip(">").split()[0] for l in open(idfile) if l.strip()}

def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name: yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line)
    if name: yield name, "".join(seq)

with open(out, "w") as o:
    for name, seq in read_fasta(fasta):
        if (name in ids) != exclude:
            o.write(f">{name}\n{seq}\n")
```

---

## 7. Coding-potential filtering

Four independent lines of evidence remove protein-coding transcripts: ORF length, protein-domain homology, and three coding-potential predictors.

### 7.1 ORF prediction (OrfPredictor)

```bash
mkdir -p orf && cd orf
touch empty_blast.txt        # no BLASTX input: ORFs predicted from sequence alone
OrfPredictor.pl ../length_filter.fasta empty_blast.txt 0 both 0.001 orf_pred.txt
```

| Argument | Meaning |
|---|---|
| input FASTA | nucleotide sequences |
| `empty_blast.txt` | optional BLASTX result (empty = none supplied) |
| `0` | no extra features |
| `both` | search both strands |
| `0.001` | E-value threshold (used only if BLASTX was supplied) |
| output | predicted peptides (frame, start, stop) |

Outputs: `noOrf.txt` (transcripts with no ORF), `ORF6frame.txt` (nucleotide ORFs), and the predicted protein FASTA.

**Get IDs of transcripts with no ORF:**

```bash
cut -f1 noOrf.txt | sed 's/^/>/' > noOrf_headers.txt
```

This keeps only the ID column (drops "No ORF found") and adds the `>` prefix.

### 7.2 Drop transcripts with long ORFs (≥ 100 aa = 300 nt) (`scripts/helpers/filter_long_orfs.py`)

```python
#!/usr/bin/env python3
"""Keep only peptides shorter than MAX_AA (default 100 aa).
Usage: python3 filter_long_orfs.py orf_proteins.fasta out.fasta [max_aa]"""
import sys

inp, out = sys.argv[1], sys.argv[2]
max_aa = int(sys.argv[3]) if len(sys.argv) > 3 else 100

name, seq = None, []
def flush(o):
    if name:
        s = "".join(seq).split("*")[0]
        if len(s) < max_aa:
            o.write(f"{name}\n{s}\n")

with open(out, "w") as o:
    for line in open(inp):
        line = line.strip()
        if line.startswith(">"):
            flush(o)
            name, seq = line.split("\t")[0], []
        else:
            seq.append(line)
    flush(o)
```

```bash
python3 ../scripts/helpers/filter_long_orfs.py orf_pred.fasta short_orf_proteins.fasta 100
```

**Why:** a transcript with a long ORF is likely to be protein-coding.

### 7.3 Protein-domain search (RPS-BLAST against Pfam)

```bash
rpsblast+ -db $PFAM_DB -query short_orf_proteins.fasta \
          -out rps.txt -evalue 0.001 -outfmt 7

# Query IDs with NO domain hit
grep -B 3 "0 hits found" rps.txt | head          # sanity check
grep -B 3 "0 hits found" rps.txt | grep "TCONS" | sed 's/# Query: /> /; s/^/>/' > zero_hits.txt
```

Combine with the no-ORF IDs, since both groups are non-coding candidates:

```bash
cat zero_hits.txt noOrf_headers.txt | sed 's/^>*/>/' | sort -u > zero_noOrf_ids.txt
python3 ../../scripts/helpers/filter_fasta_by_ids.py \
        ../length_filter.fasta zero_noOrf_ids.txt ../zero_noOrf.fasta
cd ..
```

**Why:** a peptide that matches a known protein family is probably coding. Transcripts with 0 hits or no ORF stay as lncRNA candidates.

> The `TCONS` pattern assumes StringTie/gffcompare IDs. Change it if your IDs differ.

### 7.4 CPC2

```bash
mkdir -p cpc2 && cd cpc2
CPC2.py -i ../zero_noOrf.fasta -o cpc2_out

# Keep only "noncoding" calls
grep -w "noncoding" cpc2_out.txt | cut -f1 | sed 's/^/>/' > cpc2_noncoding_ids.txt
python3 ../../scripts/helpers/filter_fasta_by_ids.py \
        ../zero_noOrf.fasta cpc2_noncoding_ids.txt cpc2_noncoding.fasta
cd ..
```

CPC2 output columns: ID, length, peptide length, Fickett score, pI, ORF integrity, coding probability, label (`coding`/`noncoding`).

### 7.5 CPAT (species-specific model)

CPAT needs a **trained model for your species**, because the default human or mouse models do not suit other organisms.

**Training data:** coding sequences (`cds_from_genomic.fna`, from NCBI) and known non-coding RNAs (NONCODE FASTA for your species).

```bash
mkdir -p cpat && cd cpat

# 1. Hexamer usage table
make_hexamer_tab -c $CDS_FA -n $NONCODE_FA > species_hexamer.tsv

# 2. Logistic regression model
make_logitmodel -x species_hexamer.tsv -c $CDS_FA -n $NONCODE_FA -o species

# 3. Predict
cpat -x species_hexamer.tsv -d species.logit.RData \
     --antisense --top-orf=1 \
     -g ../zero_noOrf.fasta -o cpat_out
```

Outputs from `make_logitmodel`: `species.feature.xls`, `species.logit.RData`, `species.make_logitModel.r`.

**Choose the cutoff:** keep transcripts below the coding-probability cutoff, which is calculated from the model's ROC curve for your species (not the human default of 0.364). Use the same cutoff for all samples.

```bash
CUT=0.431    # replace with your species-specific cutoff
# transcripts scored below the cutoff
awk -v c=$CUT 'NR>1 && $NF<c {print ">"$1}' cpat_out.ORF_prob.best.tsv > cpat_low.txt
# transcripts with no ORF at all
sed 's/^/>/; s/^>>/>/' cpat_out.no_ORF.txt > cpat_noorf.txt
cat cpat_low.txt cpat_noorf.txt | sort -u > cpat_noncoding_ids.txt
cd ..
```

> Check which column holds `coding_prob` in your CPAT version (`head cpat_out.ORF_prob.best.tsv`) and edit `$NF` if needed.

### 7.6 FEELnc

```bash
mkdir -p feelnc && cd feelnc
FEELnc_codpot.pl -i ../zero_noOrf.fasta -a $REF_GTF -g $REF_FA -l $NONCODE_FA
```

| Flag | Meaning |
|---|---|
| `-i` | candidate transcripts |
| `-a` | reference annotation (mRNA, for training) |
| `-g` | genome FASTA |
| `-l` | known lncRNAs (NONCODE), used as the lncRNA training set |

Outputs include `*.lncRNA.fa`, `*.mRNA.fa`, `*.noORF.fa`, and Random Forest QC plots and summaries.

```bash
# IDs predicted lncRNA or without ORF
grep -h ">" feelnc_codpot_out/*.lncRNA.fa feelnc_codpot_out/*.noORF.fa | sort -u > feelnc_noncoding_ids.txt
cd ..
```

---

## 8. Consensus set of high-confidence lncRNAs

Keep only transcripts called non-coding by **all three** predictors:

```bash
cd lncrna
for f in cpc2/cpc2_noncoding_ids.txt cpat/cpat_noncoding_ids.txt feelnc/feelnc_noncoding_ids.txt; do
  sed 's/^>//' $f | cut -f1 | sort -u > $(basename $f .txt).sorted
done

comm -12 cpc2_noncoding_ids.sorted cpat_noncoding_ids.sorted > common_cpc2_cpat.txt
comm -12 common_cpc2_cpat.txt feelnc_noncoding_ids.sorted     > common_all_three.txt
wc -l *.sorted common_cpc2_cpat.txt common_all_three.txt     # record the counts in the README
```

`comm -12` prints lines common to both (sorted) files. Each tool has different false positives, so their intersection gives the most reliable set.

### 8.1 Final lncRNA FASTA and GTF

```bash
# FASTA
python3 ../scripts/helpers/filter_fasta_by_ids.py \
        zero_noOrf.fasta common_all_three.txt final_lncRNA.fasta

# GTF: take all lines for these transcripts from the gffcompare output
grep -F -w -f common_all_three.txt ../gffcmp.combined.gtf > final_lncRNA.gtf
```

### 8.2 Known-database novelty check

```bash
makeblastdb -in $NONCODE_FA -dbtype nucl -out db/NONCODE_species

blastn -query final_lncRNA.fasta -db db/NONCODE_species \
       -out noncode_hits.txt -outfmt 7 \
       -evalue 0.001 -num_threads 24 \
       -qcov_hsp_perc 100 -perc_identity 100 -max_target_seqs 5
```

Transcripts with a 100% match are **known** lncRNAs; the rest are **novel**.

---

## 9. Quantification

### 9.1 Re-run StringTie against the merged annotation

For mRNAs use the merged GTF; for lncRNAs use `final_lncRNA.gtf`.

```bash
stringtie bam/${ACC}.bam -G stringtie/all_merged.gtf \
          -o counts/${ACC}_new.gtf -p 10 -e -A counts/${ACC}_new_abund.out
```

| Flag | Meaning |
|---|---|
| `-e` | estimate abundance only for the given transcripts and do not assemble new ones |
| `-A` | abundance table |

Run every sample with a loop:

```bash
tail -n +2 config/samples.tsv | while read S COND ACC; do
  stringtie bam/${ACC}.bam -G stringtie/all_merged.gtf \
            -o counts/${ACC}_new.gtf -p 10 -e -A counts/${ACC}_new_abund.out
done
```

Repeat with `-G lncrna/final_lncRNA.gtf` into a separate folder for lncRNA counts.

### 9.2 Count matrices (prepDE.py)

```bash
# prepDE_input.txt: two columns, sample ID and path to its GTF
awk 'NR>1{print $1"\tcounts/"$1"_new.gtf"}' config/samples.tsv > prepDE_input.txt

python3 prepDE.py -i prepDE_input.txt -v
```

Outputs:
- `gene_count_matrix.csv` for **mRNA / gene-level** DE
- `transcript_count_matrix.csv` for **lncRNA / transcript-level** DE (IDs are TCONS)

> Do not mix them up. lncRNA analysis must use the transcript matrix.

---

## 10. Differential expression (edgeR)

### 10.1 Install (once)

```bash
sudo apt update && sudo apt install -y r-base
```

```r
install.packages("BiocManager")
BiocManager::install("edgeR")
```

### 10.2 Script (`scripts/08_deg_analysis.R`)

```r
library(edgeR)

counts <- read.csv("gene_count_matrix.csv", row.names = 1)   # or transcript_count_matrix.csv
meta   <- read.delim("config/samples.tsv")
group  <- factor(meta$condition[match(colnames(counts), meta$sample_id)],
                 levels = c("control", "infected"))

y <- DGEList(counts = counts, group = group)

# 1. Remove low-expressed genes
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]

# 2. Normalize library size (TMM)
y <- calcNormFactors(y)

# 3. Model + dispersion
design <- model.matrix(~ group)
y   <- estimateDisp(y, design)

# 4. Test (quasi-likelihood F-test)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)

res <- topTags(qlf, n = Inf)$table
res$significant <- ifelse(res$FDR < 0.05 & abs(res$logFC) >= 1,
                          ifelse(res$logFC > 0, "UP", "DOWN"), "NS")

write.csv(res, "deg/DEG_infected_vs_control.csv")
write.csv(cpm(y, log = TRUE), "deg/logCPM_normalized.csv")

# QC plots
png("deg/MDS.png");  plotMDS(y, col = as.numeric(group)); dev.off()
png("deg/BCV.png");  plotBCV(y); dev.off()
```

Default thresholds: **FDR < 0.05** and **|log2FC| ≥ 1**. If you have multiple time points (e.g. 2 DPI, 6 DPI, 10 DPI), repeat the contrast for each and store them separately.

---

## 11. Downstream characterization of lncRNAs

### 11.1 Remove tRNA-derived and miRNA-precursor sequences

```bash
# tRNA database
makeblastdb -in tRNA_mature.fa -dbtype nucl -parse_seqids -out db/tRNA_db

blastn -query lncrna/final_lncRNA.fasta -db db/tRNA_db \
       -out downstream/lnc_vs_tRNA.tsv \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
       -evalue 1e-5 -word_size 7 -dust no -soft_masking false

# miRNA precursor database
makeblastdb -in miRBase_species.fa -dbtype nucl -parse_seqids \
            -out db/miRNA_precursor_db -title "miRBase miRNA precursors"

blastn -query lncrna/final_lncRNA.fasta -db db/miRNA_precursor_db \
       -out downstream/lnc_vs_miRNA.tsv -outfmt 6 -evalue 1e-5 -word_size 7
```

Short word size and relaxed E-value are used because these target sequences are very short. Transcripts with high coverage hits may be precursors rather than genuine lncRNAs.

### 11.2 Riboswitch and Rfam screen (Infernal)

```bash
cmscan --cpu 8 --fmt 2 \
       --tblout downstream/riboswitch.tbl \
       -o downstream/riboswitch.out \
       $RFAM_RIBOSWITCH_CM \
       lncrna/final_lncRNA.fasta
```

| Flag | Meaning |
|---|---|
| `cmscan` | scans sequences against covariance models (RNA family models) |
| `--cpu` | threads |
| `--fmt 2` | table format that includes clan information |
| `--tblout` | tabular hits |

**Convert the hit table to a clean TSV:**

```bash
awk 'BEGIN{OFS="\t"; print "idx","Riboswitch","Rfam_ID","Transcript_ID","Query_Start","Query_End","Strand","Score","E_value"}
     !/^#/ { split($4,a,"::"); print $1,$2,$3,a[1],$10,$11,$12,$17,$18 }' \
    downstream/riboswitch.tbl > downstream/riboswitch_results.tsv
```

> Column numbers depend on the `--fmt` used. Inspect `head downstream/riboswitch.tbl` and adjust.

---

## 12. Functional annotation of DE genes (OmicsBox)

### 12.1 Prepare gene BED and FASTA

```bash
# Combine DE gene lists from all contrasts and remove duplicates
cat deg/contrast1_DEG.txt deg/contrast2_DEG.txt deg/contrast3_DEG.txt \
  | awk -F'\t' '!seen[$1 FS $2]++' > deg/DEG_unique.txt
```

Convert to BED (chrom, start, end, gene ID; see Section 6.2 for the approach), then:

```bash
bedtools getfasta -fi $REF_FA -bed deg/DEG.bed -fo deg/DEG.fasta -nameOnly
```

### 12.2 OmicsBox steps (GUI)

1. **Functional Analysis → Load Sequences** and select the DEG FASTA. For duplicate IDs choose **"Add with new ID"**.
2. **BLAST:** use DIAMOND BLASTX, taxonomy filter = *your species*, E-value 1e-10.
3. **InterProScan:** run the scan and save.
4. **GO Mapping**, then **GO Annotation**; save.
5. **Merge InterPro GOs** and save.
6. **EggNOG-mapper:** run on the same FASTA and save; then **Merge EggNOG GOs**.
7. **EC-code mapping** from GO annotation.
8. **Combined Pathway Analysis:** set the top-priority taxon to your species.

---

## 13. GSEA

### 13.1 Ranked list

A two-column, tab-delimited file: **gene ID** and **log2 fold-change** (sorted descending). Use the gene/BED ID and the logFC from the edgeR table.

```bash
# Example: columns taken from the DEG table (edit column numbers)
awk -F, 'NR>1{print $1"\t"$2}' deg/DEG_infected_vs_control.csv | sort -k2,2gr > deg/ranked_list.rnk
```

Load this into OmicsBox/GSEA, with GO annotations as the gene sets.

---

## 14. Correlation between lncRNAs and mRNAs

Compute Spearman correlation on normalized expression (logCPM):

```r
lnc <- read.csv("deg/logCPM_lncRNA.csv", row.names = 1)
mrna <- read.csv("deg/logCPM_mRNA.csv", row.names = 1)

cor_mat <- cor(t(lnc), t(mrna), method = "spearman")
pairs <- as.data.frame(as.table(cor_mat))
colnames(pairs) <- c("lncRNA", "mRNA", "rho")
strong <- subset(pairs, abs(rho) > 0.9)       # strong positive or negative
write.csv(strong, "deg/lnc_mRNA_correlation_strong.csv", row.names = FALSE)
```

**Why:** highly correlated lncRNA–mRNA pairs are candidate regulatory relationships (co-expression).

---

## 15. Circos plot

### 15.1 Karyotype from the genome index

Circos needs a karyotype file, which a non-model genome does not come with, so build it from the FASTA index:

```bash
awk '{print "chr - "$1" "$1" 0 "$2" grey"}' ${REF_FA}.fai > circos/karyotype.txt
```

For readability, keep only the main chromosomes (drop small scaffolds).

### 15.2 Data tracks

Two tab-delimited tracks: **chromosome, start, end, log2FC** (one for DE genes, one for DE lncRNAs).

```
chr1   49876915   49877483   -9.94
chr1   92021419   92022645   -7.85
chr1   76319066   76321283    7.78
```

If the karyotype names differ (e.g. `NC_052532.1` vs `chr1`), rename consistently across all files.

### 15.3 `circos.conf`

```
karyotype = karyotype.txt

<ideogram>
  <spacing>
    default = 0.005r
  </spacing>
  radius    = 0.85r
  thickness = 20p
  fill      = yes
</ideogram>

<plots>
  <plot>
    type    = histogram
    file    = genes_track.txt
    r0      = 0.90r
    r1      = 0.98r
  </plot>
  <plot>
    type    = histogram
    file    = lncrna_track.txt
    r0      = 0.80r
    r1      = 0.88r
  </plot>
</plots>

<image>
  <<include etc/image.conf>>
</image>
<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
```

```bash
circos -conf circos/circos.conf
```

---

## 16. Repository layout and what to commit

```
project/
├── README.md
├── WORKFLOW.md                (this file)
├── environment.yml
├── .gitignore
├── config/        paths.sh, samples.tsv
├── scripts/       01_download.sh ... 08_deg_analysis.R
│   └── helpers/   length_filter.py, filter_fasta_by_ids.py, filter_long_orfs.py
├── results/       small tables only (DEG lists, count matrices, summaries)
├── figures/
└── docs/
```

**`.gitignore`**

```
raw/
clean/
bam/
reference/
db/
*.fastq
*.fastq.gz
*.sam
*.bam
*.bai
*.fna
*.fa
*.fasta
*.ht2
*.RData
```

**Commit:**

```bash
git add README.md WORKFLOW.md environment.yml .gitignore config scripts results figures docs
git commit -m "Add transcriptome and lncRNA analysis workflow"
git push origin main
```

**Do not commit:** raw reads, BAMs, genome indexes, server IPs, usernames, SSH commands or absolute personal paths. Large result files should go to Zenodo or Figshare, with a link in the README.

---

## 17. Tool versions and installation notes

| Tool | Version | Purpose |
|---|---|---|
| fastp | 0.21 | QC and trimming |
| HISAT2 | 2.2.1 | alignment |
| samtools | 1.11 | BAM handling |
| StringTie | 2.1.4 | assembly and quantification |
| gffcompare | 0.11 | class codes |
| bedtools | latest | BED/FASTA operations |
| OrfPredictor | latest | ORF prediction |
| BLAST+ / RPS-BLAST | 2.9.0+ | homology and domain search |
| CPC2 | latest | coding potential |
| CPAT | latest | coding potential (species-trained) |
| FEELnc | latest | lncRNA classification |
| Infernal (cmscan) | latest | Rfam / riboswitch |
| edgeR (R) | latest | differential expression |
| OmicsBox | latest | functional annotation |
| Circos | latest | circular plots |

CPC2, FEELnc, OrfPredictor.pl, prepDE.py and Circos are not on Conda, so install them from their GitHub or project pages and add them to your `PATH`. Record the exact versions you used in the README.
