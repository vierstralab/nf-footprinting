#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"


process unstarch {
    tag "${id}"
    publishDir "${params.outdir}/${id}"
    conda params.conda

    input: 
        tuple val(id), path(hotspot_peak)

    output:
        tuple val(id), path(name)

    script:
    name = "${id}.intervals.bed"
    """
    unstarch ${hotspot_peak} > ${name}
    """
}

process learn_dm {
    tag "${id}"
    publishDir "${params.outdir}/${id}"
    conda params.conda
    label "high_cpu"

    input:
        tuple val(id), path(bam_file), path(bam_index), path(unstarch_file)
    
    output:
        tuple val(id), path(name)

    script:
    name = "${id}.dm.json"
    """
    ftd learn_dm \
        --bias_model_file ${params.bias} \
        ${unstarch_file} \
        ${bam_file} \
        ${params.genome_fasta_file} \
        --n_threads ${task.cpus} \
        > dm.json

    mv dm.json ${name}
    """
}

process plot_dm {
    tag "${id}"
    publishDir "${params.outdir}/${id}"
    conda params.conda

    input:
        tuple val(id), path(dispersion_model)

    output:
        tuple val(id), path(name)

    script:
    name = "${id}.dm.pdf"
    """
    ftd plot_dm ${dispersion_model} --histograms 5,50,100 > dm.pdf
    mv dm.pdf ${name}
    """
}

process detect_dm {
    tag "${id}"
    publishDir "${params.outdir}/${id}"
    conda params.conda
    label "high_cpu"

    input:
        tuple val(id), path(dispersion_model), path(bam_file), path(bam_index), path(unstarch_file)
    
    output:
        tuple val(id), path(name)

    script:
    name = "${id}.out.bedgraph"
    """
    ftd detect \
        --bias_model_file ${params.bias} \
        --dispersion_model_file ${dispersion_model} \
        ${unstarch_file} \
        ${bam_file} \
        --n_threads ${task.cpus} \
        ${params.genome_fasta_file}
    
    mv out.bedgraph ${name}
    """
}

process retrieve_dm { 
    tag "${id}"
    publishDir "${params.outdir}/${id}"
    conda params.conda

    input:
        tuple val(id), path(bedgraph), val(threshold)

    output:
        tuple val(id), val(threshold), path(name)

    script:
    name = "${id}.interval.all.fps.${threshold}.bed"
    """
    cat ${bedgraph} \
        | awk -v OFS="\t" -v cutoff=${threshold} '\$8 <= cutoff { print \$1, \$2-3, \$3+3; }' \
        |  bedops -m - \
        > ${name}
    """
 }
 
 process learn_bayes {
    tag "${id}"
    publishDir "${params.outdir}/${id}"
    conda params.conda
    scratch true

    input:
        tuple val(id), path(bedgraph)

    output:
        tuple val(id), path(name)

    script:
    name = "${id}.beta.txt"
    """
    ftd learn_beta \
        --fdr_cutoff 0.05 \
        --exp_cutoff 10 \
        ${bedgraph}
    mv beta.txt ${name}
    """
}


 workflow footprintsCalling {
    take:
        metadata_channel

    main:
        thresholds = Channel.from(0.1, 0.05, 0.01, 0.001, 0.0001)

        unstarch_channel = metadata_channel
            | map(it -> tuple(it[0], it[3]))
            | unstarch

        bams_channel = metadata_channel
            | map(it -> tuple(it[0], it[1], it[2]))
            | join(unstarch_channel)

        out = bams_channel
            | learn_dm
            | join(bams_channel)
            | detect_dm
            | combine(thresholds)
            | retrieve_dm
            | learn_bayes
        
        plot_dm(learn_dm.out)

    emit:
        out
}

 workflow {
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.ag_id,
             file(row.filtered_alignments_bam), 
             file(row?.bam_index ?: "${row.filtered_alignments_bam}.crai"),
             file(row.hotspot_peaks_point1per)))
        | footprintsCalling
 }
