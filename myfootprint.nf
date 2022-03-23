#!/usr/bin/env nextflow


params.genome_file = "/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.fa"
params.bias_file = "/home/jhutchinson/projects/encode4plus/vierstra_et_al.6mer-model.txt"
params.bam_file = "/net/seq/data/aggregations/LN2497/aggregation-4169/LN2497.GRCh38_no_alts.sorted.bam"
params.hotspots_file = "/net/seq/data/aggregations/LN2497/aggregation-4169/peaks_v2_1_1/LN2497.GRCh38_no_alts.uniques.sorted.hotspots.fdr0.05.starch"


process learn_dm {
	
  label "footprints"


  """
  module load bedops;
  unstarch $params.hotspots_file \
    | grep -v "_random" \
    | grep -v "chrUn" \
    | grep -v "chrM" \
    > intervals.bed;
  ftd learn_dm \
    --bias_model_file $params.bias_file \
    intervals.bed \
    $params.bam_file \
    $params.genome_file \
    > dm.json
  """ 

}
