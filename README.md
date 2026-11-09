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
