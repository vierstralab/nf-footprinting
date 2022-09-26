#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

outdir = params.outdir

//  params.conda = "$moduleDir/environment.yml"


process unstarch {
  tag "${id}"
  publishDir "${outdir}/AG${id}"

  conda params.conda

  input: 
    tuple val(id), val(hotspot_peak)

  output:
    tuple val(id), path(name)

  script:
  name = "intervals.${id}.bed"
  """
  unstarch ${hotspot_peak} > ${name}
  """
}

process learn_dm {
	tag "${id}"

  publishDir "${outdir}/AG$id"
  conda params.conda
  memory = '8 GB'
  cpus = 8

  input:
    tuple val(id), val(bam_file), path(unstarch_file)
  
  output:
    tuple val(id), path(name)

  script:
  name = "dm.${id}.json"
  """
  ftd learn_dm \
    --bias_model_file ${params.bias} \
    ${unstarch_file} \
    ${bam_file} \
    ${params.genome_fasta_file} \
    --n_threads ${task.cpus} \
    > dm.json
    mv dm.json ${name}
  """
}

process plot_dm {
  tag "${id}"
  conda params.conda
  publishDir "${outdir}/AG$id"

  input:
    tuple val(id), path(dispersion_model)

  output:
    tuple val(id), path(name)

  script:
  name = 'dm.pdf'
  """
  ftd plot_dm ${dispersion_model} --histograms 5,50,100 > ${name}
  """
 }

process detect_dm {
  tag "${id}"
  publishDir "${outdir}/AG$id"
  conda params.conda
  memory = '8 GB'
  cpus = 8

  input:
    tuple val(id), val(bam_file), path(unstarch_file), path(dispersion_model)
  
  output:
    tuple val(id), path('out.bedgraph')

  script:
  """
  ftd detect \
    --bias_model_file ${params.bias} \
    --dispersion_model_file ${dispersion_model} \
    ${unstarch_file} \
    ${bam_file} \
    --n_threads ${task.cpus} \
    ${params.genome_fasta_file}
  """
}

process retrieve_dm { 
  tag "${id}"
  publishDir "${outdir}/AG$id"
  conda params.conda

  input:
    tuple val(id), path(bedgraph), val(threshold)

  output:
    tuple val(id), val(threshold), path(name)

  script:
  name = "interval.all.fps.${threshold}.${id}.bed"
  """
  cat ${bedgraph} \
  | awk -v OFS="\t" -v cutoff=${threshold} '\$8 <= cutoff { print \$1, \$2-3, \$3+3; }' \
  |  bedops -m - \
  > ${name}
 """
 }
 
 workflow footprintsCalling {
  take:
    metadata_channel
  main:
    peaks_channel = metadata_channel.map(row -> tuple(row[0], row[2]))
    unstarch_channel = unstarch(peaks_channel)
  
    bams_channel = metadata_channel
      .map(row -> tuple(row[0], row[1]))
      .join(unstarch_channel)

    dm_models = learn_dm(bams_channel)
    plot_dm(dm_models)

    bedgraphs = detect_dm(bams_channel.join(dm_models))

    thresholds = Channel.of(0.1, 0.05, 0.01, 0.001, 0.0001)
    retrieve_dm(bedgraphs.combine(thresholds))
  
  emit:
    retrieve_dm.out

 }

 workflow {
  metadata_channel = Channel
    .fromPath(params.samples_file)
    .splitCsv(header:false)
    .view { row -> "${row[0]} - ${row[1]} - ${row[2]}" }
    //.map{row -> tuple(row.ag_id, row.bam_file, row.hotspots_file)}

  footprintsCalling(metadata_channel)
 }








