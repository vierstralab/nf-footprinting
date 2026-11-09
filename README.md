# Footprint calling pipeline

Nextflow wrapper for [footprint-tools](https://github.com/vierstralab/footprint-tools) ([docs](https://footprint-tools.readthedocs.io/en/latest/))

## Requirements
- Nextflow (https://www.nextflow.io/)
- conda (https://docs.conda.io/en/latest/)



## Usage
 1) Create conda environment from `environment.yml` file or use the existing one
 2) Modify `nextflow.config` to computing enviroment specifications
 3) Fill in params paths in ```params.config```. You can also specify parameters in command line. Please find detailed explanation of the parameters in the [Config section](#config).
 4) Run the pipeline with `nextflow run main.nf -profile Altius`

## Config
There are two config files in the repository.
- ```nextflow.config``` - contains enviornment configuration. Detailed explanation can be found at https://www.nextflow.io/docs/latest/config.html. 
- ```params.config``` - specifies paths to input files.

Following parameters should be present in ```params.config```. Each option can be specified either in ```params.config``` file or with a command line.

- ```outdir``` - directory to save results into. May be a relative path.
- ```conda``` - path to an installed conda environment (from `environment.yml`). If not set, the environment is created from `environment.yml`.
- ```samples_file``` - tab-delimited file with metadata for samples. The file must contain a header and the following columns (other columns are permitted and ignored)
    - `sample_id` - unique identifier of the sample.
    - `cram_file` - path to a CRAM file with filtered alignments
    - `cram_index` - (optional) index of the CRAM file; if the column is absent, `${cram_file}.crai` is used
    - hotspots column (name set by `hotspots_column`, default `hotspots_file_0.05fdr`) - gzipped BED file of hotspots called with hotspot3; its first 3 columns are the regions to footprint
- `genome_fasta_file` - path to genome FASTA (with `.fai` index), the reference of the CRAM files
- `bias` - 6-mer DNase I cleavage bias model (default: `vierstra_et_al.6mer-model.txt` in this repository)

## Output

Written to `${outdir}/${sample_id}/`:
- `${sample_id}.dm.json`, `${sample_id}.dm.pdf` - dispersion model (`ftd learn_dm`) and its diagnostic plots (`ftd plot_dm`)
- `${sample_id}.out.bedgraph` - per-nucleotide statistics from `ftd detect`: `exp`, `obs`, `-log(pval)`, `-log(winpval)`, `fdr`
- `${sample_id}.out.bedgraph.gz` (+ `.tbi`) - the same, sorted, bgzipped and tabix-indexed
- `${sample_id}.interval.all.fps.${threshold}.bed` - footprints (merged positions with `fdr` <= threshold, extended by 3 bp) for thresholds 0.1, 0.05, 0.01, 0.001, 0.0001
- `${sample_id}.beta.txt` - parameters of the Beta prior (`ftd learn_beta`) used by `ftd posterior`

## Consensus footprints

`consensus_footprints.nf` jointly calls footprints across samples with `ftd posterior`, using the per-sample results of `main.nf`. Run it after `main.nf`, with the same `samples_file`:

```
nextflow run consensus_footprints.nf -profile Altius
```

It reads `sample_id` and the `hotspots_column` column from `samples_file`, and `${footprints_dir}/${sample_id}/${sample_id}.{dm.json,out.bedgraph.gz,beta.txt}` (`footprints_dir` defaults to the `outdir` of `main.nf`). Steps:
1. Merges the hotspots of all samples and splits them into chunks of 100 regions.
2. Writes `posterior_metadata.tsv` (`id`, `dm_file`, `tabix_file`, `beta_a`, `beta_b`) and runs `ftd posterior` on each chunk.
3. Merges the chunks into `all_samples.bedgraph.gz` (+ `.tbi`): per-nucleotide `-log(1 - posterior)` for every sample.
4. Writes `per_dataset/footprints.${sample_id}.0p01.bed.gz` for each sample: positions with `-log(1 - posterior)` > 2, extended by 3 bp and merged.

Additional parameters: `footprints_dir`, `threads` (CPUs per `ftd posterior` job, default 1). The merge steps use the `highmem` label (`bigmem` queue in the `Altius` profile).

Other entry points, both taking `--posterior_meta` (a `posterior_metadata.tsv`):
- `-entry fromInterval` - reuses the hotspot chunks from a previous run (`${outdir}/hotspot_chunks/`)
- `-entry retrieveFootprints` - only step 4, from existing merged posteriors (`tabix_file` column)

## Differential footprints

`differential_footprints.nf` finds footprints and differential footprints (positions where some groups of samples are more or less protected than the others) in each DHS, with the `differential_bayesian` module of footprint-tools ([methods](https://github.com/vierstralab/footprint-tools/blob/HEAD/footprint_tools/stats/differential_bayesian/methods.md)). It uses the per-sample results of `main.nf`.

```
nextflow run differential_footprints.nf -profile Altius \
    --dhs_index /path/to/dhs_index.tsv \
    --footprints_metadata /path/to/footprints_metadata.tsv
```

Parameters:
- `dhs_index` - tab-delimited file with a header: `dhs_id` and the DHS coordinates (`#chr` or `chrom`, `start`, `end`). One job per DHS.
- `footprints_metadata` - tab-delimited file with a header, one row per sample: `sample_id`, `dm_file` and `tabix_file` (`${sample_id}.dm.json` and `${sample_id}.out.bedgraph.gz` from `main.nf`), `extended_annotation` (group of the sample, e.g. tissue or cell type; differential footprints are tested between these groups).
- `d_min` - call threshold for the effect size `d` (default 0.5).

Output (in `outdir`):
- `per_dhs/diff_data.${dhs_id}.npz` - per-DHS results (segmentation, footprints, differential p-values); read them with `footprint_tools.stats.differential_bayesian.api.load_fit`
- `diff.{calls,regions,region_groups,dhs,cutoffs}.tsv` and `diff.summary.pdf` - genome-wide report: Benjamini-Hochberg over all positions and over DHSs, calls at the `d_min` cutoff, summary page

`-entry report` rebuilds only the genome-wide report from the files already in `${outdir}/per_dhs/`.
