#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/tfactivity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/tfactivity
    Website: https://nf-co.re/tfactivity
    Slack  : https://nfcore.slack.com/channels/tfactivity
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TFACTIVITY  } from './workflows/tfactivity'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_tfactivity_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_tfactivity_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_tfactivity_pipeline'
include { PREPARE_GENOME          } from './subworkflows/local/prepare_genome'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta     = getGenomeAttribute('fasta')
params.gtf       = getGenomeAttribute('gtf')
params.blacklist = getGenomeAttribute('blacklist')
params.taxon_id  = getGenomeAttribute('taxon_id')

if (!params.motifs && !params.taxon_id) {
    error "Please provide either a motifs file or a taxon ID"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_TFACTIVITY {

    take:
    samplesheet // channel: samplesheet read in from --input
    samplesheet_bam // channel: samplesheet read in from --input_bam
    counts_design // channel: counts design file read in from --counts_design

    main:

    ch_versions = Channel.empty()

    ch_fasta = Channel.value(file(params.fasta, checkIfExists: true))
    ch_gtf   = Channel.value(file(params.gtf, checkIfExists: true))
    ch_blacklist = params.blacklist ? Channel.value(file(params.blacklist, checkIfExists: true)) : Channel.value([])
    ch_motifs  = params.motifs ? Channel.value(file(params.motifs, checkIfExists: true)) : Channel.empty()
    ch_counts = Channel.value(file(params.counts, checkIfExists: true))
    ch_taxon_id = (!params.motifs && params.taxon_id) ? Channel.value(params.taxon_id) : Channel.empty()

    //
    // SUBWORKFLOW: Prepare genome
    //
    PREPARE_GENOME (
        ch_fasta,
        ch_gtf
    )

    ch_extra_counts = counts_design.filter{ meta, file -> file }

    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // WORKFLOW: Run pipeline
    //
    TFACTIVITY (
        samplesheet,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.gtf,
        ch_blacklist,
        ch_motifs,
        ch_taxon_id,
        PREPARE_GENOME.out.gene_lengths,
        PREPARE_GENOME.out.gene_map,
        PREPARE_GENOME.out.chrom_sizes,

        // ChromHMM
        samplesheet_bam,
        params.chromhmm_states,
        params.chromhmm_threshold,
        params.chromhmm_enhancer_marks.split(','),
        params.chromhmm_promoter_marks.split(','),

        // Peaks
        params.window_size,
        params.decay,
        params.merge_samples,
        params.affinity_aggregation,

        // Counts
        ch_counts,
        ch_extra_counts,
        Channel.value(file(params.counts_design, checkIfExists: true))
            .map{ design -> [[id: "design"], design]},
        params.min_count,
        params.min_tpm,
        params.expression_aggregation,
        params.min_count_tf,
        params.min_tpm_tf,

        // Dynamite
        params.dynamite_ofolds,
        params.dynamite_ifolds,
        params.dynamite_alpha,
        params.dynamite_randomize,

        // Ranking
        params.alpha,

        ch_versions
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_TFACTIVITY (
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.samplesheet_bam,
        PIPELINE_INITIALISATION.out.counts_design
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
