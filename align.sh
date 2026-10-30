#!/bin/bash

#SBATCH --job-name=lis_smg                   
#SBATCH --cpus-per-task=12              
#SBATCH --mem=80G             
#SBATCH --account=your.hpc.account
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
#SBATCH --time=24:00:00


INPUT_DIR="path/to/basecalled/rawfastq/files"
FILE_PATTERN="*.fastq"

if [[ -z $SLURM_ARRAY_TASK_ID ]]; then
    FILES=($INPUT_DIR/$FILE_PATTERN)
    MAX_IDX=$(( ${#FILES[@]} - 1 ))
    mkdir -p logs
    sbatch --array=0-$MAX_IDX "$0"
    exit 0
fi


FILES=($INPUT_DIR/$FILE_PATTERN)
INPUT_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
BASENAME=$(basename "${INPUT_FILE}" .fastq)

ml apptainer

OUTPUT_DIR="/output/analysis"

SCRATCH_DIR="/your/path/to/tmp/folder/"

APPTAINER_PATH="/bioinfo/tools"

mkdir "$BASENAME"

cd  "$BASENAME"

#Filt 200bp

apptainer exec $APPTAINER_PATH/filtlong_v0.3.1.sif  filtlong --min_length 200 --keep_percent 90 "$INPUT_DIR"/"$BASENAME".fastq > "$OUTPUT_DIR"/"$BASENAME"_filt.fastq

#Nanostats

apptainer exec $APPTAINER_PATH/nanostat_v1.6.0.sif  NanoStat --fastq "$OUTPUT_DIR"/"$BASENAME"_filt.fastq --name "$BASENAME"_nanostat.tsv  --outdir "$OUTPUT_DIR" --tsv

#minimap against gtdb

apptainer exec "$APPTAINER_PATH"/minimap_v.2.31.sif minimap2 -ax map-ont -t 10 --secondary=no --split-prefix "$SCRATCH_DIR"/"$BASENAME"_  /shares/stephan.ils.uzh/bioinfo/databases/gtdb/gtdb_genomes_reps_r232/gtdb_v232.mmi "$BASENAME"_filt.fastq  > "$BASENAME".sam 

apptainer exec "$APPTAINER_PATH"/samtools_v1.24.sif samtools view -@ 10 -bS "$BASENAME".sam  > "$BASENAME".bam

apptainer exec "$APPTAINER_PATH"/samtools_v1.24.sif samtools sort -@ 10  "$BASENAME".bam > "$BASENAME"_sort.bam

apptainer exec "$APPTAINER_PATH"/samtools_v1.24.sif samtools index -@ 10  "$BASENAME"_sort.bam

apptainer exec "$APPTAINER_PATH"/samtools_v1.24.sif samtools coverage --ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY "$BASENAME"_sort.bam > "$BASENAME"_cov_cleaned.tsv

apptainer exec "$APPTAINER_PATH"/samtools_v1.24.sif samtools flagstat -@ 4  "$BASENAME"_sort.bam > "$BASENAME"_align_stats.txt

rm "$BASENAME".sam "$BASENAME".bam

#Extract listeria reads

apptainer exec "$APPTAINER_PATH"/samtools_v1.24.sif samtools view -@ 10 -H "$BASENAME"_sort.bam |  awk -F'\t' '$1=="@SQ"{sn="";ln=0;
    for(i=2;i<=NF;i++){if($i~/^SN:/)sn=substr($i,4); if($i~/^LN:/)ln=substr($i,4)}
    print sn"\t"ln}' > "$BASENAME"_refs_all.tsv

grep -F -f /bioinfo/databases/gtdb/gtdb_genomes_reps_r232/database/Listeria_genomes/Listeria_list_complete_gcf_ids_gtdb.txt "$BASENAME"_refs_all.tsv | awk -F'\t' -v OFS='\t' '{print $1,0,$2}' > "$BASENAME"_targets.bed

wc -l "$BASENAME"_targets.bed

apptainer exec "$APPTAINER_PATH"/samtools_v1.24.sif samtools view -@ 10 -b -M -L "$BASENAME"_targets.bed -o "$BASENAME"_lis_sort.bam "$BASENAME"_sort.bam

apptainer exec "$APPTAINER_PATH"/samtools_v1.24.sif samtools index -@ 10 "$BASENAME"_lis_sort.bam

apptainer exec "$APPTAINER_PATH"/samtools_v1.24.sif samtools fastq -@ 10 -n "$BASENAME"_lis_sort.bam  > "$BASENAME"_lis_sort.fastq

apptainer exec "$APPTAINER_PATH"/minimap_v.2.31.sif minimap2 -ax map-ont -t 10 -I 10G -K 2G --secondary=no --split-prefix $SCRATCH_DIR"/"$BASENAME"_   /bioinfo/databases/gtdb/gtdb_genomes_reps_r232/database/Listeria_genomes/Listeria_genomes_gtdbv232.fna "$BASENAME"_lis_sort.fastq > "$BASENAME"_lis.sam

apptainer exec "$APPTAINER_PATH"/samtools_v1.24.sif samtools view -@ 10 -bS "$BASENAME"_lis.sam  > "$BASENAME"_lis.bam

apptainer exec "$APPTAINER_PATH"/samtools_v1.24.sif samtools sort -@ 10  "$BASENAME"_lis.bam > "$BASENAME"_sortlis.bam

apptainer exec "$APPTAINER_PATH"/samtools_v1.24.sif samtools index -@ 10  "$BASENAME"_sortlis.bam

apptainer exec "$APPTAINER_PATH"/samtools_v1.24.sif samtools flagstat -@ 4 "$BASENAME"_sortlis.bam  > "$BASENAME"_align_lis_stats.txt

rm "$BASENAME"_lis.bam  "$BASENAME"_lis.sam

apptainer exec "$APPTAINER_PATH"/samtools_v1.24.sif samtools coverage --ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY "$BASENAME"_sortlis.bam > "$BASENAME"_cov_lis_cleaned.tsv
