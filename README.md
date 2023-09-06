# Footprint calling pipeline

Nextflow wrapper for footprint-tools (https://footprint-tools.readthedocs.io/en/latest/)

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
- ```conda``` - path to installed conda (from environment.yml)
- ```samples_file``` - tab-delimited file with metadata for samples. The file must contain a header and the following columns (other columns are permitted and ignored)
    - ```ag_id``` - unique identifier of the sample.
    - `filtered_alignments_bam` - path to a bam/cram file
    - `hotspot_peaks_point1per` - path to peaks called with hotspot2

- `genome_fasta_file` - path to genome fasta


