#!/usr/bin/env nextflow

Channel.fromPath(params.reads)
    .splitCsv(sep: '\t', header: true)
    .map { row ->
        def all_R1 = [row.L001_R1, row.L002_R1].findAll { it }
        def all_R2 = [row.L001_R2, row.L002_R2].findAll { it }
        [row.id, all_R1, all_R2]
    }
    .set { reads_ch }

process fastqc {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith(".html")) "html/$filename"
            else if (filename.endsWith(".zip")) "zip/$filename"
            else filename
        }

    input:
    tuple val(sample_id), path(reads_R1), path(reads_R2)

    output:
    path("*.{zip,html}"), emit: reports
    path("."), emit: dir

    script:
    """
    ${params.fastqc} ${reads_R1.join(' ')} ${reads_R2.join(' ')} -o . -t ${task.cpus}
    """
}


workflow {
    fastqc(reads_ch)
}

workflow.onComplete {
    println "FastQC : ${params.outdir}/fastqc"
}
