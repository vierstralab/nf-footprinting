#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"
params.footprints_metadata = '/net/seq/data2/projects/ENCODE4Plus/REGULOME/footprints/fp_meta_filtered.tsv'
params.dhs_index = '/net/seq/data2/projects/ENCODE4Plus/REGULOME/footprints/cons_dhs.tsv'
params.d_min = 0.5 // call threshold


// Per DHS: stage-1 segmentation, footprints and differential p-values (footprint_tools differential_bayesian)
process diff_footprints {
    tag "${dhs_id}"
    conda params.conda
    publishDir "${params.outdir}/per_dhs/"
    label "diff_footprints"

    input:
        val dhs_id

    output:
        path name

    script:
    name = "diff_data.${dhs_id}.npz"
    """
    OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus} \
    python3 $moduleDir/bin/extract_fp_data.py \
        ${dhs_id} \
        ${params.dhs_index} \
        ${params.footprints_metadata} \
        ${name} \
        --d_min ${params.d_min}
    """
}


// Genome-wide: BH over all positions and over DHSs, calls, tables and the summary page
process diff_report {
    conda params.conda
    publishDir params.outdir
    label "diff_report"

    input:
        path diff_data

    output:
        path "diff.*"

    script:
    """
    python3 -m footprint_tools.stats.differential_bayesian.report . --prefix diff
    """
}


// Report only, from results already in outdir/per_dhs
workflow report {
    Channel.fromPath("${params.outdir}/per_dhs/diff_data.*.npz", checkIfExists: true)
        | collect
        | diff_report
}


workflow {
    Channel.fromPath(params.dhs_index)
        | splitCsv(header: true, sep: '\t')
        | map(it -> it.dhs_id)
        | diff_footprints
        | collect
        | diff_report
}
