#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"
params.footprints_metadata = '/net/seq/data2/projects/ENCODE4Plus/REGULOME/footprints/fp_meta_filtered.tsv'
params.dhs_index = '/net/seq/data2/projects/ENCODE4Plus/REGULOME/footprints/cons_dhs.tsv'
params.sample_maps = ''      // TSV from -entry sample_maps; empty = k = 1, f = 0
params.model_config = ''     // JSON from fit_params; empty = package defaults
params.n_profile_dhs = 5000  // DHSs sampled for fitting the sample maps

// Run order:
//   1. -entry sample_maps                          -> sample_maps.tsv
//   2. --sample_maps sample_maps.tsv               -> model_config.json (and preliminary results)
//   3. --sample_maps ... --model_config model_config.json


def dhs_ids() {
    Channel.fromPath(params.dhs_index)
        | splitCsv(header: true, sep: '\t')
        | map(it -> it.dhs_id)
}

def extract(dhs_id, output, extra) {
    """
    python3 $moduleDir/bin/extract_fp_data.py \
        ${dhs_id} \
        ${params.dhs_index} \
        ${params.footprints_metadata} \
        ${output} \
        ${extra}
    """
}


process sample_profile {
    tag "${dhs_id}"
    conda "/home/sabramov/miniconda3/envs/pytorch"
    label "diff_footprints"

    input:
        val dhs_id

    output:
        path name

    script:
    name = "profile.${dhs_id}.npz"
    extract(dhs_id, name, "--profile")
}


process fit_sample_maps {
    conda "/home/sabramov/miniconda3/envs/pytorch"
    publishDir params.outdir

    input:
        path profiles

    output:
        path "sample_maps.tsv"

    script:
    """
    python3 -m footprint_tools.stats.differential_bayesian.samples profile.*.npz --out sample_maps.tsv
    """
}


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
    extract(dhs_id, name, "--stats ${stats}"
        + (params.sample_maps ? " --sample_maps ${params.sample_maps}" : "")
        + (params.model_config ? " --model_config ${params.model_config}" : ""))
}


// Global constants (tau2, b_var, omega_scale, group nuggets) from the per-DHS statistics.
// Rerun the workflow with --model_config set to this output.
process fit_params {
    conda "/home/sabramov/miniconda3/envs/pytorch"
    publishDir params.outdir

    input:
        path stats

    output:
        path "model_config.json"

    script:
    """
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
    python3 $moduleDir/bin/generate_summary.py \
        ${dhs_id} \
        ${diff_data} \
        ${name}
    """
}


workflow summarize {
    take:
        results
    main:
        results
            | diff_summary
            | map(it -> it[1])
            | collectFile(
                name: 'diff_test_summary.tsv',
                storeDir: params.outdir,
                skip: 1,
                keepHeader: true
            )
}


workflow sample_maps {
    dhs_ids()
        | randomSample(params.n_profile_dhs, 0)
        | sample_profile
        | collect
        | fit_sample_maps
}


workflow summary_only {
    Channel.fromPath("${params.outdir}/per_dhs/diff_data.*.npz", checkIfExists: true)
        | map(it -> tuple(it.baseName - 'diff_data.', it))
        | summarize
}


workflow {
    dhs_ids() | diff_footprints
    diff_footprints.out.stats | collect | fit_params
    diff_footprints.out.results | summarize
}
