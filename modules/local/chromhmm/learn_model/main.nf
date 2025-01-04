process LEARN_MODEL {
    tag "$meta.id"
    label "process_high"

    conda "bioconda::chromhmm=1.25"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chromhmm:1.25--hdfd78af_0' :
        'biocontainers/chromhmm:1.25--hdfd78af_0' }"

    input:
    tuple val(meta), path(binarized_bams, stageAs: "input")
    val states

    output:
    tuple val(meta), path("output/emissions_${states}.txt"), path("output/*_${states}_dense.bed"), emit: model
    path "versions.yml",                                                                           emit: versions

    script:
    """
    # Organism (PLACEHOLDER) only needed for downstream analysis of ChromHMM and therefore not supplied

    ChromHMM.sh LearnModel \\
        -p $task.cpus \\
        input \\
        output \\
        $states \\
        PLACEHOLDER \\
        -Xmx${task.memory.toMega()}M

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chromhmm: \$(ChromHMM.sh Version | cut -f4 -d" ")
    END_VERSIONS
    """
}
