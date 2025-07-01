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

workflow {


    fastp_results = fastp(reads_ch)
}
