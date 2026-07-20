#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"
params.footprints_metadata = '/net/seq/data2/projects/ENCODE4Plus/REGULOME/footprints/fp_meta_filtered.tsv'
params.dhs_index = '/net/seq/data2/projects/ENCODE4Plus/REGULOME/footprints/cons_dhs.tsv'


process diff_footprints {
    tag "${dhs_id}"
    conda "/home/sabramov/miniconda3/envs/pytorch"
    publishDir "${params.outdir}/per_dhs/"

    input: 
        val dhs_id

    output:
        tuple val(dhs_id), path(name)

    script:
    name = "diff_data.${dhs_id}.npz"
    """
    python3 $moduleDir/bin/extract_fp_data.py \
        ${dhs_id} \
        ${params.dhs_index} \
        ${params.footprints_metadata} \
        ${params.fp_index_w_hotspots_path} \
        ${name}
    """
}


process diff_summary {
    tag "${dhs_id}"
    conda "/home/sabramov/miniconda3/envs/pytorch"

    input: 
        tuple val(dhs_id), path(diff_data)

    output:
        tuple val(dhs_id), path(name)

    script:
    name = "diff_summary.${dhs_id}.tsv"
    """
    python3 $moduleDir/bin/generate_summary.py \
        ${dhs_id} \
        ${diff_data} \
        ${name}
    """
}


workflow {
    Channel.fromPath(params.dhs_index)
        | splitCsv(header: true, sep: '\t')
        | map(it -> it.dhs_id)
        | diff_footprints
        | diff_summary
        | map(it -> it[1])
        | collectFile(
            name: 'diff_test_summary.tsv',
            storeDir: params.outdir,
            skip: 1,
            keepHeader: true
        )
}
