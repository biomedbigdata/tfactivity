process COMBINE_COUNTS {
    tag "$meta.id"
    label "process_single"

    conda "conda-forge::mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6==fccb0c41a243c639e11dd1be7b74f563e624fcca-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6:fccb0c41a243c639e11dd1be7b74f563e624fcca-0':
        'biocontainers/mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6:fccb0c41a243c639e11dd1be7b74f563e624fcca-0' }"

    input:
    tuple val(meta), path(counts)
    tuple val(samples), path(extra_files)
    tuple val(meta2), path(gene_map)
    val(agg_method)

    output:
    tuple val(meta), path("*.clean.tsv"), emit: counts
    tuple val(meta), path("genes.txt")  , emit: genes
    path  "versions.yml"                , emit: versions

    script:
    template "combine.py"
}
