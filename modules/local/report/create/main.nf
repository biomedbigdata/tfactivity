import groovy.json.JsonOutput

process CREATE {
    label "process_medium"

    conda "bioconda::mulled-v2-ab48c38c3be93a696d7773767d9287b4a0d3bf19==e3c8a1ac0a27058d7922e8b6d02f303c30d93e3a-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ab48c38c3be93a696d7773767d9287b4a0d3bf19:e3c8a1ac0a27058d7922e8b6d02f303c30d93e3a-0':
        'biocontainers/mulled-v2-ab48c38c3be93a696d7773767d9287b4a0d3bf19:e3c8a1ac0a27058d7922e8b6d02f303c30d93e3a-0' }"

    cache false

    input:
    tuple val(meta), path(tf_ranking)
    tuple val(meta2), path(tg_ranking)
    tuple val(meta3), path(differential)
    val(params)
    path(schema)

    output:
    path("report")      , emit: report
    path("versions.yml"), emit: versions

    script:
    params_string = JsonOutput.toJson(params)
    template "build.py"
}
