#!/usr/bin/env nextflow

Channel.fromPath(params.reads)
    .splitCsv(sep: '\t', header: true)
    .map { row ->
        def all_R1 = [row.L001_R1, row.L002_R1].findAll { it }
        def all_R2 = [row.L001_R2, row.L002_R2].findAll { it }
        [row.id, all_R1, all_R2]
    }
    .set { reads_ch }

process fastp {
    tag "$sample_id"
    publishDir "${params.outdir}/fastp",
        mode: 'copy',
        pattern: "*.{html,json,fastq.gz}",
        saveAs: { filename ->
            if (filename.endsWith(".html")) "html/$filename"
            else if (filename.endsWith(".json")) "json/$filename"
            else if (filename.contains("trimmed")) "trimmed/$filename"
            else filename
            }

    input:
    tuple val(sample_id), path(reads_R1), path(reads_R2)

    output:
    tuple val(sample_id), path("*_R1_trimmed.fastq.gz"), path("*_R2_trimmed.fastq.gz"), emit: trimmed
    path("*fastp.html"), emit: html
    path("*fastp.json"), emit: json

    script:
    """
    ${params.fastp} \
        -i ${reads_R1[0]} \
        -I ${reads_R2[0]} \
        -o ${sample_id}_R1_trimmed.fastq.gz \
        -O ${sample_id}_R2_trimmed.fastq.gz \
        --html ${sample_id}_fastp.html \
        --json ${sample_id}_fastp.json \
        -p \
        -P 20 \
        --thread 25 \
        --detect_adapter_for_pe
    """
}

process bwa_mem {
    tag "$sample_id"
    publishDir "${params.outdir}/bwa_mem",
        mode: 'copy',
        pattern: "*.bam",
        saveAs: { filename -> "aligned/$filename" }

    input:
    tuple val(sample_id), path(reads_R1), path(reads_R2)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam

    script:
    """
    ${params.bwa} index ${reference}
    ${params.samtools} faidx ${reference}
    ${params.bwa} mem -t 20 ${reference} ${reads_R1} ${reads_R2} -o ${sample_id}.bam
    """
}

process sort_bam {
    tag "$sample_id"
    publishDir "${params.outdir}/sorted_bam",
        mode: 'copy',
        pattern: "*.sorted.bam",
        saveAs: { filename -> "sorted/$filename" }

    input:
    tuple val(sample_id), path(input_bam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: sorted_bam
    path("${sample_id}.sorted.bam.bai"), emit: bai

    script:
    """
    ${params.samtools} sort -@ 20 -o ${sample_id}.sorted.bam ${input_bam}
    ${params.samtools} index ${sample_id}.sorted.bam
    """
}

process mark_duplicates {
    memory '62.GB'
    tag "$sample_id"
    publishDir "${params.outdir}/mark_duplicates",
        mode: 'copy',
        pattern: "*.{bam,bai,metrics}",
        saveAs: { filename ->
            if (filename.endsWith(".bam")) "bam/$filename"
            else if (filename.endsWith(".bai")) "bai/$filename"
            else if (filename.endsWith(".metrics")) "metrics/$filename"
            else filename
            }

    input:
    tuple val(sample_id), path(input_bam)

    output:
    tuple val(sample_id), path("${sample_id}.markdup.bam"), emit: markdup_bam
    path("${sample_id}.markdup.bam.bai"), emit: markdup_bai
    path("${sample_id}.markdup.metrics"), emit: metrics

    script:
    """
    ${params.gatk} MarkDuplicates \
        -I ${input_bam} \
        -O ${sample_id}.markdup.bam \
        -M ${sample_id}.markdup.metrics \
        --CREATE_INDEX true \
        -R '@RG\tID:${sample_id}\tSM:${sample_id}\tPL:ILLUMINA\tLB:lib1\tPU:unit1'

    if [[ ! -f "${sample_id}.markdup.bam.bai" ]]; then
        ${params.samtools} index ${sample_id}.markdup.bam
    fi
    """
}

workflow {
    fastp_results = fastp(reads_ch)
    bwa_results = bwa_mem(fastp_results.trimmed, params.reference)
    sorted_results = sort_bam(bwa_results.bam)
    markdup_results = mark_duplicates(sorted_results.sorted_bam)
}
