#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"
params.footprints_metadata = '/net/seq/data2/projects/ENCODE4Plus/REGULOME/footprints/fp_meta_filtered.tsv'
params.dhs_index = '/net/seq/data2/projects/ENCODE4Plus/REGULOME/footprints/cons_dhs.tsv'


process diff_footprints {
    tag "${dhs_id}"
    conda "/home/sabramov/miniconda3/envs/pytorch"
    publishDir "${params.outdir}/${dhs_id}"

    input: 
        val dhs_id

    output:
        tuple val(dhs_id), path("*${dhs_id}.npz"), path("summary.${dhs_id}.tsv")

    script:
    """
    python3 $moduleDir/bin/extract_fp_data.py \
        ${dhs_id} \
        ${params.dhs_index} \
        ${params.footprints_metadata} \
        ${params.fp_index_w_hotspots_path} \
        ./
    """
}


workflow {
    Channel.fromPath(params.dhs_index)
        | splitCsv(header: true, sep: '\t')
        | map(it -> it.dhs_id)
        | diff_footprints
        | map(it -> it[2])
        | collectFile(
            name: 'summary.tsv',
            storeDir: params.outdir,
            skip: 1,
            keepHeader: true
        )
}
