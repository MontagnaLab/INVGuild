> [!IMPORTANT]
> Work in progress

## Denoising

The default values of the dada2 parameters are usually good for most applications, but you should experiment with the values of `--p-trunc-len-f` and `--p-trunc-len-r` and compare the results (in terms of model fitting, number of retained sequences per sample and sequences length) to choose the best values. 

Here an example of how to do it, assuming that q2-dada2 output files are named like this (<run.id> represent the variable part of the filename identifing a run using specific dada2 parameters):
- denoising-stats_<run.id>.qza
- table_<run.id>.qza
- rep-seqs_<run.id>.qza

```bash

# Export files
tmpdir="tmp_${TLF}_${TLR}"
mkdir -p "$tmpdir"

qiime tools export --input-path "denoising-stats_${TLF}_${TLR}.qza" --output-path "$tmpdir/exported_stats"
cp "${tmpdir}/exported_stats/stats.tsv" "denoising-stats_${TLF}_${TLR}.tsv"

qiime tools export --input-path "table_${TLF}_${TLR}.qza" --output-path "$tmpdir/exported_table"
biom convert -i "$tmpdir/exported_table/feature-table.biom" -o "table_${TLF}_${TLR}.tsv" --to-tsv

qiime tools export --input-path "rep-seqs_${TLF}_${TLR}.qza" --output-path "$tmpdir/exported_seqs"
cp "${tmpdir}/exported_seqs/dna-sequences.fasta" "rep-seqs_${TLF}_${TLR}.fasta"

rm -r "$tmpdir"


```

R script for comparing different dada2 runs.

```bash
library(tidyr)
library(dplyr)
library(ggpubr)
library(Biostrings)
library(plotly)

### PLOT 1 - for checking number of retained reads

# read data
cc=0
dat_list <- list()
for (fil in list.files(pattern = "denoising-stats_.*tsv")) {
  cc=cc+1
  xx<-read.csv(fil, sep="\t")
  xx<-xx[-1,]
  xx_long <- xx %>% 
    pivot_longer(
      cols      = -sample.id,
      names_to  = "metric",
      values_to = "value"
    )
  dat_list[[cc]]<-as.data.frame(xx_long)
}
names(dat_list) <- sub("^[^_]+_([^\\.]+)\\.tsv$", "\\1",
                       list.files(pattern = "denoising-stats_.*tsv"))

df <- bind_rows(dat_list, .id = "source_df")
df$value<-as.numeric(df$value)

# define comparisons and columns of interest
my_comparisons<-combn(unique(df$source_df), 2, simplify = FALSE)

ToPlot<-c("percentage.of.input.passed.filter",
          "percentage.of.input.merged",
          "percentage.of.input.non.chimeric")

# make the plot
plts<-list()
for (i in 1:length(ToPlot)) {
  dd<-df[df$metric==ToPlot[i],]
  plts[[i]]<-ggboxplot(dd, x = "source_df", y = "value",
                        add = "jitter",
                        title = ToPlot[i]) + xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

ggarrange(plotlist = plts, ncol=3)
#ggsave("Figure_xx.png", device="png", 
#       width = 29.7*1.5, height = 21, units = "cm", bg="white")

### PLOT 2 - for checking amplicon length

# read data
SeqsLengthTotalData<-list()
SeqsLengthPerSample<-list()
for (fil in list.files(pattern="^table_[0-9_]+\\.tsv$")) {
  id <- sub("^table_(.*)\\.tsv$", "\\1", fil)
  
  # read table
  counts <- read.delim(fil, skip = 1, check.names = FALSE) 
  rownames(counts) <- counts$"#OTU ID"
  counts <- counts[ , -1, drop = FALSE]
  
  # read sequences
  seqs <- readDNAStringSet(paste0("rep-seqs_",id,".fasta"))
  seq_lengths <- width(seqs)
  names(seq_lengths) <- names(seqs)
  
  # match only shared features (just for safety)
  shared_features <- intersect(rownames(counts), names(seq_lengths))
  counts <- counts[shared_features, , drop = FALSE]
  seq_lengths <- seq_lengths[shared_features]
  
  # save length of each seqs in total data
  seq_lengths2<-as.data.frame(seq_lengths)
  seq_lengths2$ASV<-row.names(seq_lengths2)
  SeqsLengthTotalData[[id]]<-seq_lengths2
  
  # compute average length per sample
  average_lengths <- sapply(counts, function(sample_counts) {
    total_reads <- sum(sample_counts)
    if (total_reads == 0) return(NA)
    sum(sample_counts * seq_lengths) / total_reads
  })
  
  # save average length per sample
  SeqsLengthPerSample[[id]]<-as.data.frame(average_lengths)
  SeqsLengthPerSample[[id]]$SampleID<-row.names(SeqsLengthPerSample[[id]])
  
}

# make the plot - static version
PLOTS_length<-list()
for (i in 1:length(SeqsLengthTotalData)) {
  df<-SeqsLengthTotalData[[i]]
  PLOTS_length[[i]] <- ggplot(df, aes(x = seq_lengths)) +
    geom_density(fill = "skyblue", alpha = 0.5) +
    labs(title = names(SeqsLengthTotalData)[i],
         x = "Sequence Length",
         y = "Density") +
    theme_pubr()
}

ggarrange(plotlist = PLOTS_length, nrow = 3, ncol = 3)
#ggsave("Figure_xx.png", device="png", 
#       width = 29.7*1.5, height = 21, units = "cm", bg="white")

# make the plot - interactive version
long_df <- bind_rows(
  lapply(names(SeqsLengthTotalData), function(name) {
    df <- SeqsLengthTotalData[[name]]
    df$method <- name
    return(df)
  })
)

p<-ggplot(long_df, aes(x = seq_lengths, color = method)) +
  geom_density(size = 1) +
  theme_minimal() +
  labs(
    x = "Sequence Length",
    y = "Density",
    color = "Denoising Method",
    title = "Sequence Length Distribution by Denoising Method"
  )

#ggsave("Figure_xx.png", device="png", 
#       width = 29.7*1.5, height = 21, units = "cm", bg="white")


ggplotly(p)


```

