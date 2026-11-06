#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"
params.footprints_metadata = '/net/seq/data2/projects/ENCODE4Plus/REGULOME/footprints/fp_meta_filtered.tsv'
params.dhs_index = '/net/seq/data2/projects/ENCODE4Plus/REGULOME/footprints/cons_dhs.tsv'
params.report_rule = 'base'         // calls shown in debug reports: 'base' (what the pipeline calls) or 'window'
params.spot_column = 'SPOT3_score'
params.report_ids = null            // optional file, one dhs_id per line: report only these DHSs


process diff_footprints {
    tag "${dhs_id}"
    conda "/home/sabramov/miniconda3/envs/pytorch"
    publishDir "${params.outdir}/per_dhs/"
    label "diff_footprints"

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


process diff_report {
    tag "${dhs_id}"
    conda "/home/sabramov/miniconda3/envs/pytorch"

    input:
        tuple val(dhs_id), path(diff_data)

    output:
        path pdf, emit: pdf
        path candidates, emit: candidates
        path summary, emit: summary

    script:
    pdf = "diff_report.${dhs_id}.pdf"
    candidates = "diff_report_candidates.${dhs_id}.tsv"
    summary = "diff_report_summary.${dhs_id}.tsv"
    """
    python3 $moduleDir/bin/plot_diff_report.py \
        ${dhs_id} \
        ${params.dhs_index} \
        ${params.footprints_metadata} \
        ${diff_data} \
        ${pdf} \
        ${candidates} \
        ${summary} \
        --rule ${params.report_rule} \
        --spot_column ${params.spot_column}
    """
}


process zip_reports {
    publishDir params.outdir

    input:
        path pdfs

    output:
        path name

    script:
    name = "diff_reports.zip"
    """
    python3 -c "import glob, zipfile; z = zipfile.ZipFile('${name}', 'w'); [z.write(f) for f in sorted(glob.glob('*.pdf'))]; z.close()"
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


workflow summary_only {
    Channel.fromPath("${params.outdir}/per_dhs/diff_data.*.npz", checkIfExists: true)
        | map(it -> tuple(it.baseName - 'diff_data.', it))
        | summarize
}


// Diagnostic reports for every DHS with saved results: one zip of PDFs, plus all candidates and per-DHS summaries as tables
workflow debug_reports {
    ids = params.report_ids ? file(params.report_ids).readLines()*.trim().findAll() as Set : null
    reports = Channel.fromPath("${params.outdir}/per_dhs/diff_data.*.npz", checkIfExists: true)
        | map(it -> tuple(it.baseName - 'diff_data.', it))
        | filter { ids == null || it[0] in ids }
        | diff_report
    reports.candidates
        | collectFile(name: 'diff_report_candidates.tsv', storeDir: params.outdir, skip: 1, keepHeader: true)
    reports.summary
        | collectFile(name: 'diff_report_summary.tsv', storeDir: params.outdir, skip: 1, keepHeader: true)
    reports.pdf
        | collect
        | zip_reports
}


workflow {
    Channel.fromPath(params.dhs_index)
        | splitCsv(header: true, sep: '\t')
        | map(it -> it.dhs_id)
        | diff_footprints
        | summarize
}
