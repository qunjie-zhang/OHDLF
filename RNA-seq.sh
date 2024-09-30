#!/bin/bash
# Usage example： RNA-seq.sh SRR18493739 SRR18493740
# Check if a sample ID has been provided as an argument.
LOG_FILE="RNA-seq_$(date +"%Y%m%d_%H%M%S").log"

if [ $# -eq 0 ]; then
    echo "没有提供样本ID。用法: $0 SAMPLE_ID [SAMPLE_ID...]"
    exit 1
fi

exec > >(tee -a "$LOG_FILE") 2>&1

# Iterate through each sample ID passed as an argument.
for sample_id in "$@"
do
    echo "开始处理样本: $sample_id"

    # Trimmomatic质量控制
    echo "开始 Trimmomatic 质量控制"
    START_TIME=$SECONDS
    trimmomatic PE \
        ${sample_id}_1.fastq.gz ${sample_id}_2.fastq.gz \
        ${sample_id}_1.trimed.fastq.gz ${sample_id}_1.unpaired.fastq.gz \
        ${sample_id}_2.trimed.fastq.gz ${sample_id}_2.unpaired.fastq.gz \
        ILLUMINACLIP:/home/hywang02/miniconda3/envs/RNA-seq/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10:1:true \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
    ELAPSED_TIME=$(($SECONDS - $START_TIME))
    echo "Trimmomatic 质量控制完成，耗时 $ELAPSED_TIME 秒"

    # Trinity组装
    echo "开始 Trinity 组装"
    START_TIME=$SECONDS
    Trinity --seqType fq --CPU 32 --max_memory 100G --output ${sample_id}_trinity \
        --full_cleanup --min_contig_length 200 --min_kmer_cov 3 \
        --left ${sample_id}_1.trimed.fastq.gz --right ${sample_id}_2.trimed.fastq.gz
    ELAPSED_TIME=$(($SECONDS - $START_TIME))
    echo "Trinity 组装完成，耗时 $ELAPSED_TIME 秒"

    # CD-HIT去除冗余，卡0.99 -n一定要填11 不然报错
    echo "开始 CD-HIT 去除冗余"
    START_TIME=$SECONDS
    cd-hit-est -i ${sample_id}_trinity.Trinity.fasta -o ${sample_id}_trinity_cdhit99.fasta -c 0.99 -n 11 -d 0 -T 8 -M 16000
    ELAPSED_TIME=$(($SECONDS - $START_TIME))
    echo "CD-HIT 去除冗余完成，耗时 $ELAPSED_TIME 秒"
	
	# 预测转录本
    echo "开始预测转录本"
    START_TIME=$SECONDS
    TransDecoder.LongOrfs -t ${sample_id}_trinity_cdhit99.fasta
    TransDecoder.Predict -t ${sample_id}_trinity_cdhit99.fasta
    ELAPSED_TIME=$(($SECONDS - $START_TIME))
    echo "转录本预测完成，耗时 $ELAPSED_TIME 秒"

    
done

echo "所有样本处理完毕。"


