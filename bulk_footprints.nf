#!/usr/bin/env nextflow


// TODO
// revamp processes to pull from metadata_ch

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//PARAMS
/// sample inputs
params.metadata = "/home/jhutchinson/projects/encode4plus/meta/test_metadata_for_footprints.csv"

// first input from file
Channel
  .fromPath(params.metadata)
  .splitCsv(header:true)
  .map{row -> tuple(row.id, row.ln_number, row.ds_number)}
  .set {id_ch}
id_ch.into {id_ch1; id_ch2; id_ch3 ; id_ch4; id_ch5}

Channel
  .fromPath(params.metadata)
  .splitCsv(header:true)
  .map{row -> tuple(file(row.all_alignments_bam), file(row.bam_index))}
  .set {alignments_ch}
alignments_ch.into {alignments_ch1; alignments_ch2}

Channel
  .fromPath(params.metadata)
  .splitCsv(header:true)
  .map{row -> tuple(file(row.hotspot_peaks))}
  .set {hotspots_ch}

/// common files for all samples
//// bias
params.bias_file = "/home/jhutchinson/projects/encode4plus/footprinting_nf/vierstra_et_al.6mer-model.txt"
bias = file(params.bias_file)
//// genome
params.genome  = "GRCh38"
params.genome_sub = "no_alts"
genome = file("/net/seq/data/genomes/human/${params.genome}/noalts/${params.genome}_${params.genome_sub}.fa")
genome_index = file("${genome}.fai")


////////////////////////////////////////////////////////////////////////////////////////////////////////////
///CHANNELS
intervals_ch = Channel.create()
models_ch = Channel.create()
footprints_ch = Channel.create()

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////PROCESSES


unstarch_ch = id_ch1.merge(hotspots_ch)
process unstarch {
  label "unstarch"
  publishDir "./DS$ds_number/LN$ln_number/AG$id/"

  input: 
  set(id, ln_number, ds_number, file(hotspot_peaks)) from unstarch_ch

  output:
  file 'intervals.bed' into intervals_ch1, intervals_ch2

  beforeScript 'module load bedops'

  script:
  """
  unstarch $hotspot_peaks > intervals.bed
  """
}


learn_dm_ch = id_ch2.merge(alignments_ch1).merge(intervals_ch1)
process learn_dm {
	label "learn_dm"
  publishDir "./DS$ds_number/LN$ln_number/AG$id/"

  memory = '8 GB'
  cpus = 8

  input:
  set(id, ln_number, ds_number, file(all_alignments_bam), file(bam_index), file(intervals)) from learn_dm_ch
  file genome from genome
  file bias from bias
  
  output:
  file 'dm.json' into models_ch1, models_ch2

  script:
  """
  ftd learn_dm \
    --bias_model_file $bias \
    intervals.bed \
    $all_alignments_bam \
    $genome \
    --n_threads 8 \
    > dm.json
  """
}


plot_dm_ch = id_ch3.merge(models_ch1)
process plot_dm {
  label "plot_dm"
  publishDir "./DS$ds_number/LN$ln_number/AG$id/"

  input:
  
  set(id, ln_number, ds_number, file('dm.json')) from plot_dm_ch

  output:
  file 'dm.pdf'

  script:
  """
  ftd plot_dm dm.json --histograms 5,50,100 >dm.pdf
  """
 }


detect_dm_chr=id_ch4.merge(intervals_ch2).merge(models_ch2).merge(alignments_ch2)
process detect_dm {
  label "detect_dm"
  publishDir "./DS$ds_number/LN$ln_number/AG$id/"

  memory = '8 GB'
  cpus = 8

   input:
   file bias from bias
   file genome from genome
   
    set(id, ln_number, ds_number, file(intervals), file('dm.json'), file(all_alignments_bam), file(bam_index)) from detect_dm_chr

  // file 'intervals.bed' from intervals_ch2
   //file 'dm.json' from models_ch2
   //set file(bam_index), file(all_alignments_bam) from metadata_ch

   output:
   file('out*bed')
   file('out.bedgraph') into footprints_ch

   script:
   """
   ftd detect \
   --bias_model_file $bias \
   --dispersion_model_file dm.json \
   $intervals \
   $all_alignments_bam \
   --n_threads 8 \
   $genome 
   """
}


thresholds_ch = Channel.of(0.1, 0.05, 0.01, 0.001, 0.0001)
temp_ch = id_ch5.merge(footprints_ch)
retrieve_params_ch = temp_ch.combine(thresholds_ch)
//retrieve_params_ch.view()


process retrieve_dm { 
  label "retrieve_dm"
  publishDir "./"

  input:
  set(id, ln_number, ds_number,  file('out.bedgraph'), val(threshold)) from retrieve_params_ch
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
 
 








