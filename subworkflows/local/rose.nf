include { GAWK as FILTER_CONVERT_GTF } from '../../modules/nf-core/gawk'
include { GNU_SORT as SORT_BED } from '../../modules/nf-core/gnu/sort'
include { GNU_SORT as SORT_CHROM_SIZES } from '../../modules/nf-core/gnu/sort'
include { BEDTOOLS_SLOP as CONSTRUCT_TSS } from '../../modules/nf-core/bedtools/slop'
include { BEDTOOLS_SUBTRACT as FILTER_PREDICTIONS } from '../../modules/nf-core/bedtools/subtract'
include { BEDTOOLS_COMPLEMENT as INVERT_TSS } from '../../modules/nf-core/bedtools/complement'
include { BEDTOOLS_MERGE as STITCHING } from '../../modules/nf-core/bedtools/merge'
include { BEDTOOLS_INTERSECT as TSS_OVERLAP } from '../../modules/nf-core/bedtools/intersect'
include { GAWK as FILTER_OVERLAPS } from '../../modules/nf-core/gawk'
include { BEDTOOLS_SUBTRACT as SUBTRACT_OVERLAPS } from '../../modules/nf-core/bedtools/subtract'
include { BEDTOOLS_INTERSECT as UNSTITCHED_REGIONS } from '../../modules/nf-core/bedtools/intersect'
include { GNU_SORT as CONCAT_AND_SORT } from '../../modules/nf-core/gnu/sort'

workflow ROSE {
    take:
    ch_bed
    ch_gtf
    chrom_sizes

    main:

    ch_versions = Channel.empty()

    // Convert GTF to BED format and collapse regions to a single base pair at their start positions
    FILTER_CONVERT_GTF(ch_gtf, [])

    // Downstream methods require sorted inputs
    SORT_BED(FILTER_CONVERT_GTF.out.output)

    // Sort chrom_sizes to have same ordering as bed file
    SORT_CHROM_SIZES(chrom_sizes)

    // Construct 2 * params.rose_tss_window bps window around transcription start site (TSS)
    CONSTRUCT_TSS(SORT_BED.out.sorted, SORT_CHROM_SIZES.out.sorted.map{meta, file -> file})

    INVERT_TSS(CONSTRUCT_TSS.out.bed, SORT_CHROM_SIZES.out.sorted.map{meta, file -> file})

    predicted_regions = ch_bed.branch{
        meta, file ->
        enhancers: meta.assay.contains('enhancers')
        promoters: meta.assay.contains('promoters')
    }

    ch_filter_predictions = Channel.empty()
        .mix(
            predicted_regions.enhancers.combine(CONSTRUCT_TSS.out.bed),
            predicted_regions.promoters.combine(INVERT_TSS.out.bed),
        )
        .map{meta1, pred, meta2, filtering -> [meta1, pred, filtering]}

    // Remove predictions contained within a TSS
    FILTER_PREDICTIONS(ch_filter_predictions)

    // Merge regions closer than params.rose_stichting_window bps from each other
    STITCHING(FILTER_PREDICTIONS.out.bed)

    // Get overlap counts of stitched regions with TSS
    ch_tss_overlap = STITCHING.out.bed
        .combine(CONSTRUCT_TSS.out.bed)
        .map{meta1, stitched, meta2, tss -> [meta1, stitched, tss]}

    TSS_OVERLAP(ch_tss_overlap, [[], []])

    // Filter regions that overlap at least 2 TSS
    FILTER_OVERLAPS(TSS_OVERLAP.out.intersect, [])

    // Remove regions that overlap at least 2 TSS from stitched regions
    ch_subtract_overlaps = STITCHING.out.bed
        .combine(FILTER_OVERLAPS.out.output)
        .filter{meta1, stitched, meta2, overlaps -> meta1.id == meta2.id}
        .map{meta1, stitched, meta2, overlaps -> [meta1, stitched, overlaps]}

    SUBTRACT_OVERLAPS(ch_subtract_overlaps)

    // Get original regions (before stitching) of stitched regions that overlap at least 2 TSS
    ch_unstitched_regions = FILTER_OVERLAPS.out.output
        .combine(ch_bed)
        .filter{meta1, overlaps, meta2, pred -> meta1.id == meta2.id}
        .map{meta1, overlaps, meta2, pred -> [meta1, overlaps, pred]}

    UNSTITCHED_REGIONS(ch_unstitched_regions, [[], []])

    // Combine correctly stitched (overlap with < 2 TSS) and original unstitched regions and sort
    ch_concat_and_sort = SUBTRACT_OVERLAPS.out.bed
        .combine(UNSTITCHED_REGIONS.out.intersect)
        .filter{meta1, stitched, meta2, unstitched -> meta1.id == meta2.id}
        .map{meta1, stitched, meta2, unstitched -> [meta1, [stitched, unstitched]]}

    CONCAT_AND_SORT(ch_concat_and_sort)

    ch_versions = ch_versions.mix(
        FILTER_CONVERT_GTF.out.versions,
        SORT_BED.out.versions,
        SORT_CHROM_SIZES.out.versions,
        CONSTRUCT_TSS.out.versions,
        INVERT_TSS.out.versions,
        FILTER_PREDICTIONS.out.versions,
        STITCHING.out.versions,
        TSS_OVERLAP.out.versions,
        FILTER_OVERLAPS.out.versions,
        SUBTRACT_OVERLAPS.out.versions,
        UNSTITCHED_REGIONS.out.versions,
        CONCAT_AND_SORT.out.versions,
    )

    emit:
    stitched = CONCAT_AND_SORT.out.sorted

    versions = ch_versions
}

