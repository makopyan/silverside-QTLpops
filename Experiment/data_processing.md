# Silverside Experiment Data - Alignment to *Menidia menidia* genome

*Goal*: aligning whole genome sequencing data from Gen0 and Gen5 from the silverside experiment data to the Menidia menidia reference genome, to facilitate better comparison with the more recent data.

## Data processing

Most of this data processing workflow comes from following Nicolas' workflow for the GOSL cod project (<https://github.com/therkildsen-lab/gosl-cod/blob/master/markdowns/data_processing.md>).

### 0. File preparation

On the cbsunt246 server, I created a directory for these analyses (/workdir/jrick/mme_exp_analyses/), and created the recommended subdirectories within this directory. I copied the `NexteraPE_NT.fa` file to this directory, as well as the `Sample_Pop_Library_Overview.txt` file with all of the metadata for the fastq files. From this sample info file, I pulled only the individuals from Gen0 or Gen5, and made a list of the fastq files associated with these individuals (fastq files are in `/workdir/backup/silverside/Experiment_Geo_2015/Fastq/`), because these are the only ones that we're interested in. I then created the sample table for these individuals using R.

``` r
library(tidyverse)
fastq <- read_table("sample_lists/sample_fastq_list_Gen0_Gen5", col_names = "prefix")
info <- read_table2("../Sample_Pop_Library_Overview_Gen0_Gen5.txt")

new_info <- fastq %>%
    add_column(info) %>%
    relocate(prefix=prefix,lane_number=RunLane,seq_id=ID,sample_id=GNomEx_ID,population=Population) %>%
    select(prefix:population) %>%
    add_column(data_type="pe")

write_tsv(new_info, "sample_lists/sample_fastq_table.txt")
```

### 1. Adapter clipping

``` bash
sbatch --export=\
SERVER='cbsunt246 workdir/',\
BASEDIR=/fs/cbsunt246/workdir/jrick/mme_exp_analyses/,\
SAMPLELIST=/fs/cbsunt246/workdir/jrick/mme_exp_analyses/sample_lists/sample_fastq_list_Gen0_Gen5,\
SAMPLETABLE=/fs/cbsunt246/workdir/jrick/mme_exp_analyses/sample_lists/sample_fastq_table.txt,\
RAWFASTQDIR=/fs/cbsunt246/workdir/backup/silverside/Experiment_Geo_2015/Fastq/,\
RAWFASTQSUFFIX1=_1.txt.gz,\
RAWFASTQSUFFIX2=_2.txt.gz,\
ADAPTERS=/fs/cbsunt246/workdir/jrick/mme_exp_analyses/NexteraPE_NT.fa,\
TRIMMOMATIC=/programs/trimmomatic/trimmomatic-0.39.jar,\
THREADS=4,\
ARRAY_LENGTH=25,\
SCRIPT=/fs/cbsubscb16/storage/data-processing/scripts/adapter_clipping.sh \ 
--ntasks=4 \
--array=1-25 \
--mem=3000 \
--partition=short \
--output /home/jar573/slurm_logs/adapter_clipping_mme_exp.log \
/fs/cbsubscb16/storage/data-processing/scripts/adapter_clipping.slurm
```

### 2. Quality filtering

Skipped, as suggested by Nicolas.

### 3. Mapping to reference

I was having some issues with this step, so I copied the Mmenidia_refgenome_anchored.all_renamed_v2.fasta file from Arne's directory (/fs/cbsunt246/workdir/arne/ReferenceSeq/Mmenidia_refgenome_anchored.all_renamed_v2.fasta) to my directory, and then indexed it using bowtie with the command `/programs/bowtie2-2.4.3-linux-x86_64/bowtie2-build -f Mmenidia_refgenome_anchored_all_renamed_v2.fasta Mmenidia_refgenome_anchored_all_renamed_v2` (using v2.3.x made it so that the reverse indices were not created, which apparently is a known issue).

``` bash
## Tested on small number of files first, then did full run
## Each individual takes about 10 minutes
## 525 individuals total
## job id 4144646
## first ones finished ~1.5 hours
sbatch --export=\
SERVER='cbsubscb16 storage/',\
BASEDIR=/fs/cbsubscb16/storage/jrick/mme_exp_analyses/,\
FASTQDIR=/fs/cbsubscb16/storage/jrick/mme_exp_analyses/adapter_clipped/,\
SAMPLELIST=/fs/cbsubscb16/storage/jrick/mme_exp_analyses/sample_lists/sample_fastq_list_Gen0_Gen5,\
SAMPLETABLE=/fs/cbsubscb16/storage/jrick/mme_exp_analyses/sample_lists/sample_fastq_table.txt,\
FASTQSUFFIX1=_adapter_clipped_f_paired.fastq.gz,\
FASTQSUFFIX2=_adapter_clipped_r_paired.fastq.gz,\
MAPPINGPRESET=very-sensitive,\
REFERENCE=/fs/cbsubscb16/storage/jrick/mme_exp_analyses/reference/Mmenidia_refgenome_anchored_all_renamed_v2.fasta,\
REFNAME=Mmenidia_refgenome_anchored_all_renamed_v2,\
MINQ=0,\
BOWTIE=/programs/bowtie2-2.3.4.3/bowtie2,\
SAMTOOLS=samtools,\
THREADS=8,\
ARRAY_LENGTH=50,\
SCRIPT=/fs/cbsubscb16/storage/data-processing/scripts/low_coverage_mapping.sh \
--ntasks=8 \
--array=1-50 \
--mem=3000 \
--partition=regular \
--output=/home/jar573/slurm_logs/low_coverage_mapping_mme_exp.log \
/fs/cbsubscb16/storage/data-processing/scripts/low_coverage_mapping.slurm
```

After this was completed, I noticed that some of the array instances ran out of memory during the sorting step, so there are some samples that have only a `.bam` file, but no `.sorted.bam` file. To fix this in an inelegant manner, I did the following in an interactive session:

``` bash
for seq in `cat ../sample_lists/sample_fastq_list_Gen0_Gen5 | cut -f 1 -d'_'`; 
    do file=`ls *.bam | grep $seq | grep 'sorted'`; 
    if [ ! -z "$file" ]; 
            then echo "done"; 
            else echo "failed"; 
            ls *.bam | grep $seq >> orphaned_bams; 
    fi; 
done

for bam in `cat orphaned_bams`; 
    do echo $bam; 
    newbam=`echo $bam | sed 's/\.bam/_minq0_sorted\.bam/g'`; 
    samtools view -h -q 0 $bam | samtools view -@ 1 -buS - | samtools sort -@ 1 -o $newbam; 
done
```

### 4. Count reads and combine bam files from duplicate individuals

``` bash
nohup bash /fs/cbsunt246/workdir/data-processing/scripts/count_bam_unmerged.sh \
    sample_lists/sample_fastq_list_Gen0_Gen5 sample_lists/sample_fastq_table.txt \
    /fs/cbsubscb16/storage/jrick/mme_exp_analyses/ \
    Mmenidia_refgenome_anchored_all_renamed_v2 \
    samtools \
    >& /fs/cbsubscb16/storage/jrick/mme_exp_analyses/sample_lists/bam_count_unmerged_Gen0_Gen5.tsv 2> nohup.err < /dev/null &
```

(actually, that timed out, so I'll be submitting a job instead)

``` bash
sbatch \
  --nodelist cbsubscb16 \
  --partition short \
  --nodes=1 \
  --ntasks=4 \
  --mem=100 \
  --output /local/storage/jrick/mme_exp_analyses/sample_lists/bam_count_unmerged_Gen0_Gen5.tsv \
  --error /local/storage/jrick/mme_exp_analyses/sample_lists/bam_count_error.txt \
  /local/storage/data-processing/scripts/count_bam_unmerged.sh \
  /local/storage/jrick/mme_exp_analyses/sample_lists/sample_fastq_list_Gen0_Gen5 \
  /local/storage/jrick/mme_exp_analyses/sample_lists/sample_fastq_table.txt \
  /local/storage/jrick/mme_exp_analyses/ \
  Mmenidia_refgenome_anchored_all_renamed_v2 \
  samtools \
  4
```

**Create merged sample table**

``` r
library(tidyverse)

# Define base directory and reference name
basedir <- "/fs/cbsubscb16/storage/jrick/mme_exp_analyses/"
refname <- "Mmenidia_refgenome_anchored_all_renamed_v2"

# Read in unmerged sample table
info <- read_tsv("sample_lists/sample_fastq_table.txt") %>%
    mutate(sample_seq_id=paste(population,seq_id,lane_number,sep="_"))

# Create table with only one row for each unique sample
# and write to sample_lists directory
sample_table_merged <- group_by(info, seq_id) %>%
  summarise(population=unique(population), sample_id=ifelse(n()==1,sample_id, "merged"), lane_number=ifelse(length(unique(lane_number))==1,unique(lane_number), "merged"), data_type=paste0(sort(unique(data_type)), collapse = "")) %>%
  mutate(sample_seq_id=paste(population, seq_id, lane_number, data_type, sep = "_")) %>%
  select(sample_seq_id, lane_number, seq_id, sample_id, population, data_type)

write_tsv(sample_table_merged, "sample_lists/sample_table_merged.tsv")

 # Create bam lists as inputs for future steps
bam_list_merged <- paste0(basedir, "bam/", sample_table_merged$sample_seq_id, "_bt2_", refname, "_minq0_sorted.bam")
bam_list_dedup_overlapclipped <- transmute(sample_table_merged, suffix=ifelse(data_type=="se", paste0("_bt2_", refname, "_minq0_sorted_dedup.bam"), paste0("_bt2_", refname, "_minq0_sorted_dedup_overlapclipped.bam"))) %>%
  .$suffix %>%
  paste0(basedir, "bam/", sample_table_merged$sample_seq_id, .)
bam_list_realigned <- transmute(sample_table_merged, suffix=ifelse(data_type=="se", paste0("_bt2_", refname, "_minq0_sorted_dedup_realigned.bam"), paste0("_bt2_", refname, "_minq0_sorted_dedup_overlapclipped_realigned.bam"))) %>%
  .$suffix %>%
  paste0(basedir, "bam/", sample_table_merged$sample_seq_id, .)
write_lines(bam_list_merged, "sample_lists/bam_list_merged.txt")
write_lines(bam_list_dedup_overlapclipped, "sample_lists/bam_list_dedup_overlapclipped.txt")
write_lines(bam_list_realigned, "sample_lists/bam_list_realigned.txt")
```

Now, write script to merge bamfiles from the same individuals.

``` r
## Find all duplicated samples
duplicated_samples <- sample_table_merged %>%
  filter(str_detect(sample_seq_id, "merged")) %>%
  .$seq_id
duplicated_samples_seq_ids <-sample_table_merged %>%
  filter(str_detect(sample_seq_id, "merged")) %>%
  .$sample_seq_id
merging_script<-NULL

## Loop through all duplicated samples 
for (i in 1:length(duplicated_samples)){
  duplicated_sample <- duplicated_samples[i]
  duplicated_samples_seq_id <- duplicated_samples_seq_ids[i]
  ## Extract the bam file names from the unmerged sample table
  input <- filter(info,seq_id==duplicated_sample) %>% 
    mutate(unmerged_bam=paste(sample_id,seq_id,lane_number,data_type,"bt2",refname,"minq0_sorted.bam",sep="_")) %>% 
    .$unmerged_bam %>% 
    paste0(basedir,"bam/",.) %>% 
    paste(collapse=" ")
  
  ## Paste together the command line
  merging_script[i] <- paste0("samtools merge -@ 8 -o ", basedir, "bam/", duplicated_samples_seq_id, "_bt2_", refname, "_minq0_sorted.bam ", input)
}

# Write merging script
write_lines("#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --job-name=merge_bam\n#SBATCH --output=merge_bam.log", "scripts/merge_bam.sh")
write_lines("echo 'beginning bamfile merging' ","scripts/merge_bam.sh",append=TRUE)
write_lines(merging_script, "scripts/merge_bam.sh",append = TRUE)
```

Now, run the merging script (this takes a couple of minutes per individual)

``` bash
sbatch --nodelist cbsubscb16 \
--partition regular \
--nodes=1 \
--ntasks=8 \
--mem=1000 \
--output=/home/jar573/slurm_logs/merge_bams.log \
/fs/cbsubscb16/storage/jrick/mme_exp_analyses/scripts/merge_bam.sh
```

For this dataset, there ended up being 148 samples that were duplicated, with bamfiles that need to be merged. This leaves us with a total of 376 unique individuals.

These new, merged bamfiles now have names that are slightly different (and more helpful) than the non-merged bamfiles. So, I renamed the bamfiles (again, inelegantly) for non-merged samples as well to match them.

``` bash
mkdir newbams
for id in `cat ../sample_lists/sample_table_merged.tsv | grep -v merged | cut -f 4`
do echo $id
#id=`cat ../sample_lists/sample_table_merged.tsv | grep -v merged | head -n 2 | tail -n +2 | cut -f 4`
oldbam=`ls *sorted.bam | grep "${id}_"`
new=`cat ../sample_lists/sample_table_merged.tsv | grep -v merged | grep -P "${id}\t" | cut -f 1`
newbam="newbams/${new}_bt2_Mmenidia_refgenome_anchored_all_renamed_v2_minq0_sorted.bam"
if [[ -z $oldbam || -z $newbam ]]; then
echo "something wrong. $oldbam or $newbam missing"
else
echo "copying $oldbam to $newbam"
cp $oldbam $newbam
fi
done
```

It appears that I'm still missing a few bam files compared to my sample list: 

* RGen0_1085 
* RGen0_1144 
* R1Gen5_1359 
* RGen0_1391

Three of these have bamfiles, but not sorted bamfiles, so I can fix that by running the `samtools view | samtools view | samtools sort` step from the `low_coverage_mapping.sh` script. The fourth one (1144) needs to be merged, so I'll run the merging step on its individual bamfiles to create the merged file. Once that's done, I now have 376 bamfiles ready to go in my `newbams/` directory.

### 5. Deduplicate, overlap clip, and count merged bam files

``` bash
sbatch --export=\
SERVER='cbsubscb16 storage/',\
BAMLIST=/fs/cbsubscb16/storage/jrick/mme_exp_analyses/sample_lists/bam_list_merged.txt,\
SAMPLETABLE=/fs/cbsubscb16/storage/jrick/mme_exp_analyses/sample_lists/sample_table_merged.tsv,\
JAVA=java,\
PICARD=/programs/picard-tools-2.9.0/picard.jar,\
BAMUTIL=/programs/bamUtil/bam,\
ARRAY_LENGTH=376,\
SCRIPT=/fs/cbsubscb16/storage/data-processing/scripts/deduplicate_clipoverlap.sh \
--ntasks=1 \
--array=1-376 \
--mem=65G \
--partition=short \
--output=/home/jar573/slurm_logs/deduplicate_clipoverlap.log \
/fs/cbsubscb16/storage/data-processing/scripts/deduplicate_clipoverlap.slurm
```

``` bash
## Ran this on cbsubscb16
sbatch \
  --nodelist cbsubscb16 \
  --partition regular \
  --nodes=1 \
  --ntasks=4 \
  --mem=100 \
  --output /local/storage/jrick/mme_exp_analyses/sample_lists/bam_count_merged.tsv \
  --error /local/storage/jrick/mme_exp_analyses/sample_lists/bam_count_merged_error.txt \
  /local/storage/data-processing/scripts/count_bam_merged.sh \
  /local/storage/jrick/mme_exp_analyses/sample_lists/bam_list_merged.txt \
  /local/storage/jrick/mme_exp_analyses/sample_lists/sample_table_merged.tsv \
  samtools \
  20 \
  4
```

### 6. Realign around indels

``` bash
cd /local/storage/jrick/mme_exp_analyses/bam/newbams/
cp /local/storage/jrick/mme_exp_analyses/sample_lists/bam_list_dedup_overlapclipped.txt \
/local/storage/jrick/mme_exp_analyses/sample_lists/bam_list_dedup_overlapclipped.list

## Ran on cbsubscb16 (Job ID 4169412)
## note: this step requires both a .fai index and a .dict sequence dictionary for the ref genome
## these can be created using samtools faidx and picard CreateSequenceDictionary
## note #2: this step took ~5 days to run
sbatch \
  --nodelist=cbsubscb16 \
  --partition=long30 \
  --nodes=1 \
  --ntasks=1 \
  --mem=40G \
  --output=/local/storage/jrick/mme_exp_analyses/nohups/realign_indels.log \
  /local/storage/data-processing/scripts/realign_indels.sh \
  /local/storage/jrick/mme_exp_analyses/sample_lists/bam_list_dedup_overlapclipped.list \
  /local/storage/jrick/mme_exp_analyses/ \
  /local/storage/jrick/mme_exp_analyses/reference/Mmenidia_refgenome_anchored_all_renamed_v2.fasta \
  samtools \
  /usr/local/jdk1.8.0_121 \
  /programs/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
  1
```

### 7. Count read depth per position per individually

Here, I first count each bam file individually in a loop using samtools depth. Then R can be used to loop through these individual depth files and summarize the results. This information is important for specifying min/max read depths for ANGSD. I am doing this for all of the bamfiles that I'll be working with in ANGSD, so created a list of bamfiles (`all_bamlist`) that includes the files that I just aligned, as well as the ones that Arne has worked with previously (from his `/workdir/arne/sampleinfo/BamFiles_clean_list_Mme.txt` list).

``` bash
# note: here I am using a slightly modified version of the script in the data-processing
# repository, where I get to specify an output directory
sbatch --partition=short --nodes=1 --ntasks=48 --mem=50G --output=/fs/cbsubscb16/storage/jrick/mme_exp_analyses/nohups/count_depth_per_position_per_sample.log /fs/cbsubscb16/storage/jrick/mme_exp_analyses/scripts/count_depth_per_position_per_sample.sh /fs/cbsubscb16/storage/jrick/mme_exp_analyses/bam/all_bamlist 48 20 20 samtools /fs/cbsubscb16/storage/jrick/mme_exp_analyses/bam

nohup Rscript /fs/cbsunt246/workdir/data-processing/scripts/summarize_depth_per_position.R \
"/local/storage/jrick/mme_exp_analyses/sample_lists/bam_list_realigned_mincov_filtered.txt" \
"/local/storage/jrick/mme_exp_analyses/sample_lists/sample_table_merged_mincov_filtered.tsv" \
"/local/storage/jrick/mme_exp_analyses/" \
> /local/storage/jrick/mme_exp_analyses/nohups/summarize_depth_per_position.out 2>&1 &

## I was struggling to get the others to work, so I ended up doing this interactively:
parallel -j 8 --line-buffer "echo {} && samtools depth -aa -@ 2 {} -q 20 -Q 20 | cut -f 3 | gzip > {.}.depth.gz" ::: *realigned.bam
```

From this, I find that the mean depth per site per sample is 0.96x, with a mean of 47% of the reference covered per individual. The average total sequencing depth per site is 672.4. From these results, I'm going to set my variant calling parameters as: minimum combined sequencing depth of `0.33*nind=186`, maximum combined sequencing depth of `mean+2*sd = 672.4+(2*393.8) = 1460`, and minimum number of individuals at half (282).

### 7.5 Removing problematic individuals

After chatting with Nina, there are some individuals that need to be removed because of some sample ID mix ups, and then we'll also be paring the individual list down to just those from RGen0, D1Gen5, D2Gen5, U1Gen5, U2Gen5, R1Gen5, and R2Gen5. I've removed these inds from the bamlist and have a new bamlist that I'll be using from here for the purposes of Maria's paper, `all_bamlist_MixupsRemoved_expOnly`. This changes my filtering parameters slightly, and I'll now be using a minimum combined sequencing depth of 123, the same maximum combined sequencing depth (1460), and minInd of 186.

### 8. Global SNP calling

```sh
## Run on cbsubscb16
sbatch --nodelist=cbsubscb16 \
  --partition=long7 \
  --nodes=1 \
  --ntasks=16 \
  --mem=20G \
  --output=/local/storage/jrick/mme_exp_analyses/nohups/global_snp_calling_exp_wild.log \
  /fs/cbsunt246/workdir/genomic-data-analysis/scripts/angsd_global_snp_calling.sh \
  /local/storage/jrick/mme_exp_analyses/sample_lists/all_bamlist_MixupsRemoved_expOnly \
  /local/storage/jrick/mme_exp_analyses/ \
  /local/storage/jrick/mme_exp_analyses/reference/Mmenidia_refgenome_anchored_all_renamed_v2.fasta \
  123 1460 186 20 0 20 \
  /fs/cbsunt246/workdir/programs/angsd0.935/angsd/angsd \
  16
```

It looks like this produced a total of 15,504,517 SNPs in this group of individuals.

### 9. MAF for each population

```sh
# first, make a list of ind names and populations
for bam in `cat all_bamlist_MixupsRemoved_expOnly`
  do echo $bam
  ind=$(basename $bam | cut -f 1-2 -d'_')
  pop=$(basename $bam | cut -f 1 -d'_')
  echo "$ind,$pop" >> all_bamlist_MixupsRemoved_expOnly.poplist
done

# now, make individual bamlists for each population
POPLIST=all_bamlist_MixupsRemoved_expOnly.poplist
BAMLIST=all_bamlist_MixupsRemoved_expOnly.bamlist
OUTNAME=all_bamlist_MixupsRemoved_expOnly_

for POP in `tail -n +2 $POPLIST | cut -f 2 -d',' | sort | uniq`; do 
	echo $POP
	SAMPLESEQID=`awk  -F',' -v pop="$POP" 'BEGIN{OFS=FS} $2==pop' $POPLIST | cut -f 1 -d','` 
	grep -F "${SAMPLESEQID}" $BAMLIST > bam_list_per_pop/${OUTNAME}${POP}.txt
done

# and finally, use slurm script to get maf per population
## On the cluster 
sbatch --export=\
SERVER='cbsubscb16 storage/',\
BASEDIR=/fs/cbsubscb16/storage/jrick/mme_exp_analyses/,\
SAMPLETABLE=/fs/cbsubscb16/storage/jrick/mme_exp_analyses/sample_lists/all_bamlist_MixupsRemoved_expOnly.poplist,\
POPCOLUMN=2,\
BAMLISTPREFIX=all_bamlist_MixupsRemoved_expOnly_,\
REFERENCE=/fs/cbsubscb16/storage/jrick/mme_exp_analyses/reference/Mmenidia_refgenome_anchored_all_renamed_v2.fasta,\
SNPLIST=/fs/cbsubscb16/storage/jrick/mme_exp_analyses/angsd/global_snp_list_all_bamlist_MixupsRemoved_expOnly_mindp123_maxdp1460_minind186_minq20.txt,\
MINDP=5,\
MAXDP=1000,\
MININD=5,\
MINQ=20,\
MINMAPQ=20,\
ANGSD=/fs/cbsubscb16/storage/programs/angsd0.935/angsd/angsd,\
THREADS=16,\
EXTRA_ARG='-remove_bads 1 -only_proper_pairs 1 -C 50',\
SCRIPT=/fs/cbsubscb16/storage/jrick/mme_exp_analyses/scripts/get_maf_per_pop.sh \
--ntasks=16 \
--mem=20G \
--partition=long7 \
--array=1-11 \
--output=slurm_logs/get_maf_per_pop_mme_exp_wild.log \
--nodelist=cbsubscb16 \
/fs/cbsubscb16/storage/jrick/mme_exp_analyses/scripts/get_maf_per_pop.slurm
```

