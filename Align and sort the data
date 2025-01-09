# Index the genome
hisat2-build directory/references/SC.fa \
             directory/references/index/SC_hista2_index

# Align using Hisat2
for f in "directory/raw_data"/*
do
if [[ $f =~ * ]]; then
        # align paired-end reads respectively
        data1=""
        data2=""
        for data in "$f"/*
        do
            if [[ $data =~ *_1.fq.gz ]]; then
                data1=$data
            elif [[ $data =~ *_2.fq.gz ]]; then
                data2=$data
            fi
        done
        f_base="$(basename -- "$f")"
        directory/hisat2-2.1.0/hisat2 \
        -x directory/references/index/SC_hisat2_index \
        -1 "$data1" \
        -2 "$data2" \
        -S directory/bam/$f_base.sam
    fi
done

# The SAM output was then transformed to its binary version (BAM file)
for sam in "directory/bam"/*
do   
    if [[ $sam =~ *.sam ]]; then
        sam_base="$(basename -- "$sam")"
        sam_base_filename="${sam_base%.*}"
        samtools view -b "$sam" > directory/bam/$sam_base_filename.bam
    fi   
done

# Sort the BAM file so that the reads are in order with regard to chromosome number and position
for bam in "directory/bam"/*
do
    if [[ $bam == *.bam ]]; then
        bam_base="$(basename -- "$bam")"
        bam_base_filename="${bam_base%.*}"
        samtools sort "$bam" > directory/bam/$bam_base_filename.sorted.bam
    fi   
done

# Index the sorted BAM file to make accessing the data quicker for downstream tools
for sorted_bam in "directory/bam"/*
do
    if [[ $sorted_bam == *.sorted.bam ]]; then
        samtools index "$sorted_bam"
    fi
done

# QC of alignment
for sorted_bam in “directory/bam"/*
do
    if [[ $sorted_bam == *.sorted.bam ]]; then
        sorted_bam_base="$(basename -- "$sorted_bam")"
        sorted_bam_base_filename="${sorted_bam_base%.*.*}"

# Duplication metrics: Picard’s MarkDuplicates tool on the sorted bam file
java -jar directory/picard.jar MarkDuplicates \
     INPUT= "$sorted_bam" \
     OUTPUT= directory/bam/$sorted_bam_base_filename.mkdup.bam \
     METRICS_FILE= directory/bam/$sorted_bam_base_filename.mkdup_metrics.txt \
     CREATE_INDEX=true \
     USE_JDK_DEFLATER=true USE_JDK_INFLATER=true

# Alignment metrics
java -jar directory/picard.jar CollectAlignmentSummaryMetrics \
            INPUT= "$sorted_bam" \
            OUTPUT= directory/bam/$sorted_bam_base_filename.alignment_metrics.bam \
            REFERENCE_SEQUENCE= directory/SC.fa \
            USE_JDK_DEFLATER=true USE_JDK_INFLATER=true

# RNA alignment metrics
java -jar directory/picard.jar CollectRnaSeqMetrics \
            INPUT= "$sorted_bam" \
            REF_FLAT= directory/SC_gtf3.txt \
            OUTPUT= directory/bam/$sorted_bam_base_filename.RNA_metrics.txt \
            STRAND=NONE \
            USE_JDK_DEFLATER=true USE_JDK_INFLATER=true

    fi
done

# The count reads of each gene were calculated using the featureCounts function in the SubRead package
directory/subread-2.0.0-MacOS-x86_64/bin/featureCounts \
    -t exon \
    -g gene_id \
    --primary \
    -a directory/references/SC.gtf \
    -o directory/counts/WT_gtr2_pib2.featureCounts \
    directory/bam/*.sorted.bam
