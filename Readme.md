# R Analyses performed in the Bamboo study

## Description

The Bamboo study is a study in the microbial composition in infant feces after an intervention with milf fat or vegetable fat based infant formula. This study is a collaboration between the Center for Biomics Erasmus medical center, Horaizon NV and FrieslandCampina. 

## Data preparation

In this study, stool samples were obtained from infants. The parents and caretakers recorded a number of parameters which were summarized as a metadata table. From the stool DNA was extracted and sequencing libraries were generated using the Nextera XT method. These libraries were sequenced on Illumina HiSeq2500 sequencers. Paired-end reads were generated of 300 basepairs in lengths. These data have been submitted to the SRA in BioProject PRJNA798191.

In the downstream analysis, the Illumina adapter sequences were trimmed using [AdapterTrimmer](https://github.com/erasmus-center-for-biomics/AdapterTrimmer). The trimmed reads were analyzed using [MetaPhlan 3.0](https://huttenhower.sph.harvard.edu/metaphlan/) with the CHOCOPhlan library (mpa_v30_CHOCOPhlAn_201901). The `--add_viruses` and `-t rel_ab_w_read_stats` options were set as well. Paired reads were fed consecutively as MetaPhlan does not make use of those.

Trimming and alignment were performed using a Snakemake workflow. The equivalent bash code is shown below.

```bash
cat R1.fastq R2.fastq | \
  metaphlan \
  --input_type fastq \
  --samout output.sam \
  --bowtie2out output.bowtie \
  --output_file output.tsv \
  --sample_id_key output \
  --add_viruses \
  -t rel_ab_w_read_stats \
  --bowtie2db bowtie2_reference_directory \
  --index mpa_v30_CHOCOPhlAn_201901 \
  --bt2_ps sensitive-local \
  --nproc 12
```