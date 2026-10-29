#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"
params.footprints_metadata = '/net/seq/data2/projects/ENCODE4Plus/REGULOME/footprints/fp_meta_filtered.tsv'
params.dhs_index = '/net/seq/data2/projects/ENCODE4Plus/REGULOME/footprints/cons_dhs.tsv'
params.model_config = ''   // JSON from differential_bayesian.hyper; empty = package defaults
params.seqcorr = ''        // JSON from differential_bayesian.seqcorr; empty = 6-mer expected as is


process diff_footprints {
    tag "${dhs_id}"
    conda "/home/sabramov/miniconda3/envs/pytorch"
    publishDir "${params.outdir}/per_dhs/"
    label "diff_footprints"

    input: 
        val dhs_id

    output:
        tuple val(dhs_id), path(name), emit: results
        path stats, emit: stats

    script:
    name = "diff_data.${dhs_id}.npz"
    stats = "stats.${dhs_id}.npz"
    """
    python3 $moduleDir/bin/extract_fp_data.py \
        ${dhs_id} \
        ${params.dhs_index} \
        ${params.footprints_metadata} \
        ${params.genome_fasta_file} \
        ${name} \
        ${stats} \
        ${params.model_config ? "--model_config ${params.model_config}" : ""} \
        ${params.seqcorr ? "--seqcorr ${params.seqcorr}" : ""}
    """
}


// Global fits from the per-DHS statistics: sequence correction of the expected and model constants.
// Rerun the workflow with --seqcorr and --model_config set to these outputs.
process fit_params {
    conda "/home/sabramov/miniconda3/envs/pytorch"
    publishDir params.outdir

    input:
        path stats

    output:
        path "seqcorr.json"
        path "model_config.json"

    script:
    """
    python3 -m footprint_tools.stats.differential_bayesian.seqcorr stats.*.npz --out seqcorr.json
    python3 -m footprint_tools.stats.differential_bayesian.hyper stats.*.npz --out model_config.json \
        ${params.model_config ? "--config ${params.model_config}" : ""}
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
    echo summary_v1.1
    python3 $moduleDir/bin/generate_summary.py \
        ${dhs_id} \
        ${diff_data} \
        ${name}
    """
}


workflow summary_only {
    Channel.fromPath("${params.outdir}/per_dhs/diff_data.*.npz", checkIfExists: true)
        | map(it -> tuple(it.baseName - 'diff_data.', it))
        | diff_summary
        | map(it -> it[1])
        | collectFile(
            name: 'diff_test_summary.tsv',
            storeDir: params.outdir,
            skip: 1,
            keepHeader: true
        )
}


workflow {
    Channel.fromPath(params.dhs_index)
        | splitCsv(header: true, sep: '\t')
        | map(it -> it.dhs_id)
        | diff_footprints
    diff_footprints.out.stats | collect | fit_params
    diff_footprints.out.results
        | diff_summary
        | map(it -> it[1])
        | collectFile(
            name: 'diff_test_summary.tsv',
            storeDir: params.outdir,
            skip: 1,
            keepHeader: true
        )
}
