process CONVERT_MOTIFS {
    tag "$meta.id"
    label "process_single"

    conda "bioconda:bioconductor-universalmotif==1.20.0--r43hf17093f_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-universalmotif:1.20.0--r43hf17093f_0':
        'biocontainers/bioconductor-universalmotif:1.20.0--r43hf17093f_0' }"

    input:
    tuple val(meta), path(in_file), val(in_type)
    val(out_type)

    output:
    tuple val(meta), path("${out_file}"), emit: converted
    path "versions.yml"                 , emit: versions

    script:
    out_file = "${meta.id}.converted.${out_type}"
    template "convert.R"
}
