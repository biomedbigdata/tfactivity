workflow SNEEP {
    take:

    motifs // Coming from MOTIFS subworkflow
    snp_file // Coming from user/params/download
    fasta // Coming from pipeline
    scale_file // Coming from user/params/download
    motif_regions //Coming from FIMO

    main:

    motifs.view()
    // Subset this motifs file for motifs in FIMO output (transfac format)

    snp_file.view()
    // Check for file ending -> convert if necessary
    // Subset file by motif_regions

    motif_regions.view()
    // coming from FIMO as gff file
    // use gawk to get column 3, 4, 5 (format as tsv file)
    // use gnu sort to sort those regions
    // subset SNP file based on this regions (bedtools/bcftools)

    fasta.view()
    // Do nothing, already right

    scale_file.view()
    // Do nothing, maybe subset to motifs in FIMO output (motif_regions)
    // Motifs are in column 1 of the motif_regions file

}
