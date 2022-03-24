#!/usr/bin/env nextflow


params.genome_file = "/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.fa"
params.bias_file = "/home/jhutchinson/projects/encode4plus/vierstra_et_al.6mer-model.txt"
params.bam_file = "/net/seq/data/aggregations/LN2497/aggregation-4169/LN2497.GRCh38_no_alts.sorted.bam"
params.hotspots_file = "/net/seq/data/aggregations/LN2497/aggregation-4169/peaks_v2_1_1/LN2497.GRCh38_no_alts.uniques.sorted.hotspots.fdr0.05.starch"



process unstarch {

 label "unstarch"

 input: 
 path hotspots from hotspots_file

 output:
 path(intervals) into intervals_ch

 """
   module load bedops;
  unstarch $params.hotspots_file \
    | grep -v "_random" \
    | grep -v "chrUn" \
    | grep -v "chrM" \
    > $intervals;
  """
}


process learn_dm {
	
  label "learn"

  input:
  path ref from genome_file
  path bias from bias_file
  path bam from bam_file
  path intervals from intervals_ch

  output:
  path model 


  """
  ftd learn_dm \
    --bias_model_file $bias \
    $intervals \
    $bam \
    $genome \
    > $model
  """ 

}


