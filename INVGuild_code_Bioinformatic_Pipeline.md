# **Bioinformatic pipeline used to analyse LUCAS data**

#### Primers used for the amplification
- **Euk575F**: `ASCYGYGGTAAYWCCAGC`
- **Euk895R**: `TCHNHGNATTTCACCNCT`

---

## **1. Set the directory containing paired-end sequences and the maximum number of available cores**
```bash
# raw data directory
RAWDIR=/path/to/rawdata/dir
# maximum number of processes that can be run simultaneously 
JOBS=32
```

---

## **2. Create reference database from SILVA**

#### 2.1. Import, reverse transcribe and clean sequences from [SILVA v.138.1](https://www.arb-silva.de/)
```bash
# activate qiime2 evironment
conda activate qiime2-amplicon-2024.10

# import sequences
qiime rescript get-silva-data \
  --p-version '138.1' \
  --p-target 'SSURef_NR99' \
  --o-silva-sequences silva-138.1-ssu-nr99-rna-seqs.qza \
  --o-silva-taxonomy silva-138.1-ssu-nr99-tax.qza \
  --p-no-rank-propagation \
  --p-include-species-labels \
  --parallel \
  --verbose

# reverse transcription
qiime rescript reverse-transcribe \
  --i-rna-sequences silva-138.1-ssu-nr99-rna-seqs.qza \
  --o-dna-sequences silva-138.1-ssu-nr99-seqs.qza

# remove low quality sequences (5+ degenerated bases, homopolymers with 8+ bases)
qiime rescript cull-seqs \
  --i-sequences silva-138.1-ssu-nr99-seqs.qza \
  --p-num-degenerates 5 \
  --p-homopolymer-length 8 \
  --o-clean-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza \
  --p-n-jobs $JOBS

# export taxonomy to tsv file
qiime tools export \
  --input-path silva-138.1-ssu-nr99-tax-cleaned.qza \
  --output-path exp
```

#### 2.2. Clean SILVA taxonomy
```r
R

library(dplyr)
library(stringr)

data<-read.table("exp/taxonomy.tsv", header=T, sep="\t")

# eliminate species-level ID that do not correspond to taxon names (= those starting wih lowercase) 
badIDs <- data %>%
  mutate(
    STR = str_extract(Taxon, "s__([a-z].+)")  # Extract string starting with lowercase after "s__"
  ) %>%
  filter(!is.na(STR)) %>%  # Remove NA values (cases without a match)
  pull(STR) %>%
  unique()
CleanData <- data %>%
  mutate(
    Taxon = case_when(
      # Case 1: First letter after "s__" is lowercase
      str_detect(Taxon, "s__[a-z]") ~ str_replace(Taxon, "s__[a-z].*", "s__"),
      # Case 2: First letter after "s__" is uppercase
      str_detect(Taxon, "s__[A-Z]") ~ str_replace(Taxon, 
                                                  "g__; s__([^_]+)",  # Capture "g__;" and genus portion after "s__"
                                                  "g__\\1; s__\\1"),  # Append genus after "g__" and keep it after "s__"
      # Default: Leave other rows unchanged
      TRUE ~ Taxon
    )
  )

# eliminate genus-level ID that do not correspond to taxon names (= those starting wih lowercase) 
badIDs <- CleanData %>%
  mutate(
    STR = str_extract(Taxon, "g__([a-z].+)")  # Extract string starting with lowercase after "g__"
  ) %>%
  filter(!is.na(STR)) %>%  # Remove NA values (cases without a match)
  pull(STR) %>%
  unique()  # Get unique values
# View unique genera
print(sort(badIDs))

### remove also genera wrong IDs ###
CleanData2 <- CleanData %>%
  mutate(
    Taxon = case_when(
      # If the first letter after "g__" is lowercase, remove everything after "g__" and keep "s__"
      str_detect(Taxon, "g__[a-z]") ~ str_replace(Taxon, "g__[a-z].*s__", "g__; s__"),
      
      # Otherwise, leave the string as it is
      TRUE ~ Taxon
    )
  )

write.table(CleanData2, "silva-138.1-ssu-nr99-tax-cleaned.tsv", sep="\t", quote = F, row.names = F)

q()
```

#### 2.3. Filter and dereplicate the reference database

```bash
# import the cleaned SILVA taxonomy in QIIME2
qiime tools import \
  --type FeatureData[Taxonomy] \
  --input-path silva-138.1-ssu-nr99-tax-cleaned.tsv \
  --input-format HeaderlessTSVTaxonomyFormat \
  --output-path silva-138.1-ssu-nr99-tax-cleaned.qza

# filter by length depending on taxon
qiime rescript filter-seqs-length-by-taxon \
  --i-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza \
  --i-taxonomy silva-138.1-ssu-nr99-tax-cleaned.qza \
  --p-labels Archaea Bacteria Eukaryota \
  --p-min-lens 900 1200 1400 \
  --o-filtered-seqs silva-138.1-ssu-nr99-seqs-filt.qza \
  --o-discarded-seqs silva-138.1-ssu-nr99-seqs-discard.qza 

# dereplicate sequences only if they have the same taxonomic ID
qiime rescript dereplicate \
  --i-sequences silva-138.1-ssu-nr99-seqs-filt.qza  \
  --i-taxa silva-138.1-ssu-nr99-tax-cleaned.qza \
  --p-mode 'uniq' \
  --o-dereplicated-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
  --o-dereplicated-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza \
  --p-threads $JOBS
```

#### 2.4. Extract the amplified region from the reference database using primer sequences
```bash
# extract the amplified region
qiime feature-classifier extract-reads \
  --i-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
  --p-f-primer ASCYGYGGTAAYWCCAGC \
  --p-r-primer TCHNHGNATTTCACCNCT \
  --p-identity 0.8 \
  --p-min-length 150 \
  --p-max-length 450 \
  --p-n-jobs $JOBS \
  --p-read-orientation forward \
  --o-reads silva-138.1-ssu-nr99-seqs_Euk575-895.qza 

# dereplicate sequences (again since some of these shorter sequences may be identical now)
qiime rescript dereplicate \
  --i-sequences silva-138.1-ssu-nr99-seqs_Euk575-895.qza  \
  --i-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza \
  --p-mode 'uniq' \
  --o-dereplicated-sequences silva-138.1-ssu-nr99-seqs_Euk575-895_derep-uniq.qza \
  --o-dereplicated-taxa silva-138.1-ssu-nr99-tax_Euk575-895_derep-uniq.qza \
  --p-threads $JOBS
```

---

## **3. Quality check of raw data**

```bash
# run fastQC
bash script_FastQC.sh ${RAWDIR} ${RAWDIR}/FastQC_results

# run multiqc to summarize fastQC results
multiqc ${RAWDIR}/FastQC_results
```
See *multiqc_report.html*

⚠️
There are illumina universal adapters in the reads that needs to be removed

---

## **4. Import sequences in QIIME2**
The downloaded paired-end reads are named following this scheme: "Lucas\<sampleID\>.\<sequencing-direction\>.fq.gz" (e.g., Lucas0001.1.fq.gz, Lucas0001.2.fq.gz, Lucas0002.1.fq.gz, Lucas0002.2.fq.gz, ...).

#### 4.1. Make a manifest file containing paths to raw reads

```bash
# get absolute path of reads
ls -1 ${RAWDIR}/*1.fq.gz > R1.txt
ls -1 ${RAWDIR}/*1.fq.gz > R2.txt
# get identifiers
ls ${RAWDIR}/*fq.gz | cut -c 1-9 | sort | uniq > ${RAWDIR}/ids.txt
# make manifest file
paste ${RAWDIR}/ids.txt ${RAWDIR}/R1.txt ${RAWDIR}/R2.txt > ${RAWDIR}/man.txt
echo sample-id$'\t'forward-absolute-filepath$'\t'reverse-absolute-filepath > ${RAWDIR}/manifest.tsv
cat ${RAWDIR}/man.txt >> ${RAWDIR}/manifest.tsv
rm ${RAWDIR}/*.txt
mv ${RAWDIR}/manifest.tsv manifest.tsv
```

#### 4.2. Import sequences and remove remaining Illumina adapters
```bash
# import sequences
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.tsv \
  --output-path seqs.qza \
  --input-format PairedEndFastqManifestPhred33V2
qiime demux summarize \
  --i-data seqs.qza \
  --o-visualization seqs.qzv

# remove illumina adapters
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences seqs.qza \
  --p-cores $JOBS \
  --p-adapter-f AGATCGGAAGAG \
  --p-adapter-r AGATCGGAAGAG \
  --o-trimmed-sequences seqs_trimmed.qza \
  --verbose
qiime demux summarize \
  --i-data seqs_trimmed.qza \
  --o-visualization seqs_trimmed.qzv
```

---

## **5. Denoising and sequence re-orientation**
Since barcodes and adapters were added after PCR amplification each fastq file contained both forward and reverse reads. So sequences were re-orientered using the reference database as guide.

#### 5.1. Denoising with DADA2
```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs seqs_trimmed.qza \
  --p-trim-left-f 18 \
  --p-trim-left-r 18 \
  --p-trunc-len-f 220 \
  --p-trunc-len-r 200 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-trunc-q 2 \
  --p-pooling-method 'pseudo' \
  --p-chimera-method 'pooled' \
  --p-n-reads-learn 30000000 \
  --p-n-threads $JOBS \
  --o-table table_MixedOrientation.qza \
  --o-representative-sequences rep-seqs_MixedOrientation.qza \
  --o-denoising-stats denoising-stats.qza \
  --verbose
qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv
qiime feature-table summarize \
  --i-table table_MixedOrientation.qza \
  --o-visualization table_MixedOrientation.qzv
qiime feature-table tabulate-seqs \
  --i-data rep-seqs_MixedOrientation.qza \
  --o-visualization rep-seqs_MixedOrientation.qzv
```
#### 5.2. Re-orientering sequences using the reference database
```bash
# re-orient sequences
qiime rescript orient-seqs \
  --i-sequences rep-seqs_MixedOrientation.qza \
  --i-reference-sequences silva-138.1-ssu-nr99-seqs_Euk575-895_derep-uniq.qza \
  --o-oriented-seqs rep-seqs.qza \
  --o-unmatched-seqs orientation_unmatched_sequences.qza \
  --p-threads $JOBS
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
# exclude not oriented sequences
qiime feature-table filter-features \
  --i-table table_MixedOrientation.qza \
  --m-metadata-file orientation_unmatched_sequences.qza \
  --p-exclude-ids \
  --o-filtered-table table.qza
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv
```

---

## **6. Taxonomic classification**

```bash
# fit a classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138.1-ssu-nr99-seqs_Euk575-895_derep-uniq.qza \
  --i-reference-taxonomy silva-138.1-ssu-nr99-tax_Euk575-895_derep-uniq.qza \
  --o-classifier NBclassifier.qza

# classify without a confidence threshold
qiime feature-classifier classify-sklearn \
  --i-classifier NBclassifier.qza \
  --i-reads rep-seqs.qza \
  --p-confidence=disable \
  --p-n-jobs $JOBS \
  --o-classification taxonomy_no_confidence.qza

# generate weights
qiime clawback generate-class-weights \
  --i-reference-taxonomy silva-138.1-ssu-nr99-tax_Euk575-895_derep-uniq.qza \
  --i-reference-sequences silva-138.1-ssu-nr99-seqs_Euk575-895_derep-uniq.qza \
  --i-samples table.qza \
  --i-taxonomy-classification taxonomy_no_confidence.qza \
  --o-class-weight class-weights.qza

# fit new classifier with weights
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138.1-ssu-nr99-seqs_Euk575-895_derep-uniq.qza \
  --i-reference-taxonomy silva-138.1-ssu-nr99-tax_Euk575-895_derep-uniq.qza \
  --i-class-weight Euk575-895_LUCAS_class-weights.qza \
  --o-classifier weightedNBclassifier.qza

# classify using a confidence threshold of 0.9
qiime feature-classifier classify-sklearn \
  --i-classifier weightedNBclassifier.qza \
  --i-reads rep-seqs.qza \
  --p-confidence 0.90 \
  --p-n-jobs $JOBS \
  --o-classification taxonomy.qza

```

## **7. Final filtering and exporting the results**

```bash
# remove singletons
qiime feature-table filter-features \
  --i-table table.qza \
  --p-min-frequency 2 \
  --o-filtered-table table.qza

# remove non-metazoan ASVs
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-include 'd__Eukaryota(Metazoa)' \
  --o-filtered-table table.qza

# export table (ASVs abundance per sample)
qiime tools export \
  --input-path table.qza \
  --output-path exp
biom convert -i exp/feature-table.biom -o table.tsv --to-tsv

# export taxonomy (ASVs taxonomic assignments)
qiime tools export \
  --input-path taxonomy.qza \
  --output-path exp
mv exp/taxonomy.tsv taxonomy.tsv

```

---
