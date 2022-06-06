#!/usr/bin/env nextflow

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//PARAMS


params.library = "LN2497"
params.aggregation = "4169"
params.genome  = "GRCh38"
params.genome_sub = "no_alts"
params.bias_file = "/home/jhutchinson/projects/encode4plus/footprinting_nf/vierstra_et_al.6mer-model.txt"


////////////////////////////////////////////////////////////////////////////////////////////////////////////
//FILE HANDLERS
genome = file("/net/seq/data/genomes/human/${params.genome}/noalts/${params.genome}_${params.genome_sub}.fa")
genome_index = file("${genome}.fai")
hotspots = file("/net/seq/data/aggregations/${params.library}/aggregation-${params.aggregation}/peaks_v2_1_1/${params.library}.GRCh38_no_alts.uniques.sorted.hotspots.fdr0.05.starch")
bam = file("/net/seq/data/aggregations/${params.library}/aggregation-${params.aggregation}/${params.library}.${params.genome}_${params.genome_sub}.sorted.bam")
bam_index = file("${bam}.bai")
bias = file(params.bias_file)

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//OUTPUT TO STDOUT
log.info """\
         footprinting-nf pipeline
         ===================================
         genome           : "${params.genome}_${params.genome_sub}"
         genome_file      : $genome
         genome_index     : $genome_index
         bias_file        : ${params.bias_file}
         bam_file         : $bam
         bam_file_index   : $bam_index
         hotspots         : $hotspots
         """.stripIndent()

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//publishing directories
outputDir = "./${params.library}/${params.aggregation}"

////////////////////////////////////////////////////////////////////////////////////////////////////////////
///SETUP CHANNELS
intervals_ch = Channel.create()
models_ch = Channel.create()
footprints_ch = Channel.create()

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//PROCESSES


process unstarch {
  label "unstarch"
  publishDir outputDir, mode: "copy", overwrite: true

  input: 
  file hotspots from hotspots

  output:
  file 'intervals.bed' into intervals_ch1, intervals_ch2

  beforeScript 'module load bedops'

  script:
  """
  unstarch $hotspots > intervals.bed
  """
}

process learn_dm {
	label "learn_dm"
  publishDir outputDir, mode: "copy", overwrite: true

  memory = '8 GB'
  cpus = 8

  input:
  file genome from genome
  file bias from bias
  file bam from bam
  file 'intervals.bed' from intervals_ch1

  output:
  file 'dm.json' into models_ch1, models_ch2

  script:
  """
  ln -s  $bam_index ./
  ftd learn_dm \
    --bias_model_file $bias \
    intervals.bed \
    $bam \
    $genome \
    --n_threads 8 \
    > dm.json
  """
}

process plot_dm {
  label "plot_dm"
  publishDir outputDir, mode: "copy", overwrite: true

  input:
  file 'dm.json' from models_ch1

  output:
  file 'dm.pdf'

  script:
  """
  ftd plot_dm dm.json --histograms 5,50,100 >dm.pdf
  """
 }


process detect_dm {

  label "detect_dm"
  publishDir outputDir, mode: "copy", overwrite: true

  memory = '8 GB'
  cpus = 8

   input:
   file bias from bias
   file 'dm.json' from models_ch2
   file 'intervals.bed' from intervals_ch2
   file bam from bam
   file genome from genome

   output:
   file('out*bed')
   file('out.bedgraph') into footprints_ch

   script:
   """
   ln -s  $bam_index ./
   ftd detect \
   --bias_model_file $bias \
   --dispersion_model_file dm.json \
   intervals.bed \
   $bam \
   --n_threads 8 \
   $genome 
   """
}


thresholds_ch = Channel.of(0.1, 0.05, 0.01, 0.001, 0.0001)
retrieve_params_ch = footprints_ch.combine(thresholds_ch)


process retrieve_dm { 
  label "retrieve_dm"
  publishDir outputDir, mode: "copy", overwrite: true


  input:
  set file('out.bedgraph'), val(threshold) from retrieve_params_ch
  //val threshold from thresholds_ch

  output:
  file "interval.all.fps.${threshold}.bed"

  beforeScript 'module load bedops'

  script:
  """
  cat out.bedgraph \
  | awk -v OFS="\t" -v cutoff=${threshold} '\$8 <= cutoff { print \$1, \$2-3, \$3+3; }' \
  |  bedops -m - \
  > interval.all.fps.${threshold}.bed
 """
 }








