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
    - `hotspots_file` - gzipped BED file of hotspots called with hotspot3; its first 3 columns are the regions to footprint
- `genome_fasta_file` - path to genome FASTA (with `.fai` index), the reference of the CRAM files
- `bias` - 6-mer DNase I cleavage bias model (default: `vierstra_et_al.6mer-model.txt` in this repository)

## Output

Written to `${outdir}/${sample_id}/`:
- `${sample_id}.dm.json`, `${sample_id}.dm.pdf` - dispersion model (`ftd learn_dm`) and its diagnostic plots (`ftd plot_dm`)
- `${sample_id}.out.bedgraph` - per-nucleotide statistics from `ftd detect`: `exp`, `obs`, `-log(pval)`, `-log(winpval)`, `fdr`
- `${sample_id}.out.bedgraph.gz` (+ `.tbi`) - the same, sorted, bgzipped and tabix-indexed
- `${sample_id}.interval.all.fps.${threshold}.bed` - footprints (merged positions with `fdr` <= threshold, extended by 3 bp) for thresholds 0.1, 0.05, 0.01, 0.001, 0.0001
- `${sample_id}.beta.txt` - parameters of the Beta prior (`ftd learn_beta`) used by `ftd posterior`
