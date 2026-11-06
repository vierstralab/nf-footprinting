#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"
params.footprints_metadata = '/net/seq/data2/projects/ENCODE4Plus/REGULOME/footprints/fp_meta_filtered.tsv'
params.dhs_index = '/net/seq/data2/projects/ENCODE4Plus/REGULOME/footprints/cons_dhs.tsv'
params.report_rule = 'base'         // calls shown in debug reports: 'base' (what the pipeline calls) or 'window'
params.spot_column = 'SPOT3_score'
params.report_ids = null            // optional file, one dhs_id per line: report only these DHSs
params.n_perm = 1000                // label permutations per DHS for the p-values
params.pv_batch = 20                // DHSs per p-value task (~15 min each at 1000 permutations)
params.pv_ids = null                // optional file, one dhs_id per line: p-values only for these DHSs
params.pv_q = 0.05                  // Benjamini-Hochberg level in the p-value report


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


process diff_pvalues {
    tag "${ids[0]} +${ids.size() - 1}"
    conda "/home/sabramov/miniconda3/envs/pytorch"
    publishDir "${params.outdir}/pvalues/", pattern: "*.npz"

    input:
        tuple val(ids), path(diff_data)

    output:
        path "diff_pvalues.*.npz", emit: npz
        path summary, emit: summary

    script:
    summary = "diff_pvalues_summary.${ids[0]}.tsv"
    """
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 $moduleDir/bin/diff_pvalues.py \
        ${ids.join(',')} \
        ${diff_data} \
        --n_perm ${params.n_perm} \
        --summary ${summary}
    """
}


process diff_pvalue_report {
    conda "/home/sabramov/miniconda3/envs/pytorch"
    publishDir params.outdir
    memory 64.GB

    input:
        path npz

    output:
        path "diff_pvalues.*"

    script:
    """
    python3 $moduleDir/bin/diff_pvalue_report.py . --q ${params.pv_q}
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


// Per-base and per-DHS p-values from saved stage-1 results, pooled into BH cutoffs, calls and histograms
workflow pvalues {
    ids = params.pv_ids ? file(params.pv_ids).readLines()*.trim().findAll() as Set : null
    pv = Channel.fromPath("${params.outdir}/per_dhs/diff_data.*.npz", checkIfExists: true)
        | map(it -> tuple(it.baseName - 'diff_data.', it))
        | filter { ids == null || it[0] in ids }
        | collate(params.pv_batch)
        | map(b -> tuple(b*.get(0), b*.get(1)))
        | diff_pvalues
    pv.summary
        | collectFile(name: 'diff_pvalues.per_dhs.tsv', storeDir: params.outdir, skip: 1, keepHeader: true)
    pv.npz
        | flatten
        | collect
        | diff_pvalue_report
}


// Report only, from p-values already in outdir/pvalues
workflow pvalue_report {
    Channel.fromPath("${params.outdir}/pvalues/diff_pvalues.*.npz", checkIfExists: true)
        | collect
        | diff_pvalue_report
}


workflow {
    Channel.fromPath(params.dhs_index)
        | splitCsv(header: true, sep: '\t')
        | map(it -> it.dhs_id)
        | diff_footprints
        | summarize
}
