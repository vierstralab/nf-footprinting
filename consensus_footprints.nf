#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"
// per-sample outputs of main.nf (<footprints_dir>/<sample_id>/<sample_id>.{dm.json,out.bedgraph.gz,beta.txt})
params.footprints_dir = "${launchDir}/${params.outdir}"
params.threads = 1


process create_interval {
    label "highmem"
	publishDir "${params.outdir}"
	conda params.conda

	input:
	   path "hotspots/hotspots*.starch"
	
	output:
	   path "hotspot_chunks/chunk*"


	script:
	
	name = "interval_file.bed"
	"""

        #bedops -m \$(ls -d hotspots/* | sort) | sort-bed - > ${name}
	gunzip -c hotspots/* | grep -v '^#' | sort-bed - | bedops -m - > ${name}

	mkdir -p hotspot_chunks

	cat ${name} | split -l 100 -a 5 -d - hotspot_chunks/chunk.
	"""
}


process prepare_metadata {
	publishDir "${params.outdir}"
	conda params.conda

	input:
	   val aggregations

	output:
	   path ("${name}")

	script:
	name = "posterior_metadata.tsv"
	"""	
	echo ${aggregations} | sed 's/[][]//g'  | sed 's/ //g'  | tr "," "\n" > list_of_aggregations.txt
	footprint_path="${params.footprints_dir}"

	while read agg
	do
	
		sample_id="\${agg}"
		dm_file="\${footprint_path}/\${agg}/\${agg}.dm.json"
		tabix_file="\${footprint_path}/\${agg}/\${agg}.out.bedgraph.gz"
		beta_a=`cat \${footprint_path}/\${agg}/\${agg}.beta.txt | cut -f1`
		beta_b=`cat \${footprint_path}/\${agg}/\${agg}.beta.txt | cut -f2`

		echo -e "\${sample_id}\t\${dm_file}\t\${tabix_file}\t\${beta_a}\t\${beta_b}" >> tmp.${name}
	done<list_of_aggregations.txt

	
	echo -e "id\tdm_file\ttabix_file\tbeta_a\tbeta_b" > header.txt
	cat header.txt tmp.${name} > ${name}

	"""		
}

process posterior_probabilities {
    cpus params.threads
    // label "highmem"
    conda params.conda


    input:
    	tuple(path(posterior_metadata), path(chunkFile))

    output:
	path("${name}")

    script:
    name="out_${chunkFile}.bedgraph"
    """
  
    ftd posterior ${posterior_metadata} ${chunkFile} --n_threads ${task.cpus}
    mv out.bedgraph ${name}

    """
}

process merge_and_compress {
    label "highmem"
    conda params.conda
    publishDir "${params.outdir}"
    scratch true
   
    input:
    	path("out*bedgraph")

    output:
       	tuple path("${name}"), path("${name}.tbi")

    script:
    name="all_samples.bedgraph.gz"

    """
	cat out*bedgraph | grep -v '^#' | sort-bed - > tmp.${name}
	head -2 out1bedgraph | cat - tmp.${name} | bgzip -c > ${name}
	tabix -p bed ${name}

    """
}

process retrieve_footprints {


    conda params.conda
    publishDir "${params.outdir}/per_dataset/"


    input:
	tuple path(all_samples_bedgraph), path(all_samples_bedgraph_index), val(sample_id)

    output:
        tuple val(sample_id), path("${name}")

    script:
    threshold="0.01"
    threshold_str="0p01"
    name="footprints.${sample_id}.${threshold_str}.bed.gz"
    """

    col=`(zcat ${all_samples_bedgraph} | head -2 || true) | tail -n -1 | tr "\t" "\n" | grep -w -n ${sample_id} | cut -f1 -d':'`
	
    zcat ${all_samples_bedgraph} | tail -n +3 | cut -f1-3,\${col} |  awk -v OFS="\t" '{ print \$1, \$2, \$3, ".", \$4; }' > ${sample_id}_${threshold}.bed

    cat ${sample_id}_${threshold}.bed \
    |  awk -v thresh=${threshold} '{if(\$5 > 2) print}' \
    | bedops --range 3 -m - \
    | bedmap --delim "\t" --echo --max - ${sample_id}_${threshold}.bed \
    | awk -v OFS="\t" -v id=${sample_id} '{ midpoint=int((\$3-\$2)/2)+\$2; print \$1, \$2, \$3, id, \$4, midpoint, midpoint; }' \
    | tee  footprints.${sample_id}.${threshold_str}.bed | bgzip \
    > ${name}

    """


}

 workflow {
    chunks = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
	| map { row -> file(row[params.hotspots_column]) }
	| collect(sort: true)
        | create_interval
	| flatten()
 


    sample_ids = Channel.fromPath(params.samples_file)
    	| splitCsv(header:true, sep:'\t')
    	| map { row -> "${row.sample_id}" } 
    
    sample_ids 
        | collect(sort: true)
    	| prepare_metadata
    	| combine(chunks)
    	| posterior_probabilities
        | collect(sort: true)
        | merge_and_compress
	| combine(sample_ids)
        | retrieve_footprints
 }


workflow retrieveFootprints {

    sample_ids = Channel.fromPath(params.posterior_meta)
     	| splitCsv(header:true, sep:'\t')
        | map { row -> "${row.id}" }

    Channel.fromPath(params.posterior_meta)
	| splitCsv(header:true, sep:'\t')
	| map { it -> tuple(file("${it.tabix_file}"), file("${it.tabix_file}.tbi")) }
	| combine(sample_ids)
   	| retrieve_footprints
}


workflow fromInterval {


	sample_ids = Channel.fromPath(params.posterior_meta)
	    | splitCsv(header:true, sep:'\t')
            | map { row -> "${row.id}" }

	

	chunks = Channel.fromPath("${launchDir}/${params.outdir}/hotspot_chunks/*").flatten()
	
	sample_ids
            | collect(sort: true)
            | map { file(params.posterior_meta) }
            | combine(chunks)
            | posterior_probabilities
            | collect(sort: true)
            | merge_and_compress
            | combine(sample_ids)
	    | view()
            | retrieve_footprints
 }

