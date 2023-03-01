#!/usr/bin/env nextflow
/*
========================================================================================
                         MaestSi/nf-m6anet
========================================================================================
 MaestSi/nf-m6anet analysis pipeline.
 #### Homepage / Documentation
 https://github.com/MaestSi/nf-m6anet
----------------------------------------------------------------------------------------
*/

def helpMessage() {
		log.info"""
	Usage:
	nextflow -c nf-m6anet.conf run nf-m6anet.nf --samples="/path/to/samples.txt" --resultsDir="/path/to/resultsDir" -profile docker

	Mandatory argument:
	-profile                                                 Configuration profile to use. Available: docker, singularity
	Other mandatory arguments which may be specified in the nf-m6anet.conf file

	--samples                                                Path to the tab-separated sample file including sample name, condition, path to fast5 folder and path to fastq file
	--resultsDir                                             Path to a folder where to store results
	--transcriptome_fasta                                    Path to the transcriptome fasta file
	--gtf                                                    Path to genome annotation gtf file
	--prob_mod_thr                                           Probability modification threshold for calling a site as m6A+
	--postprocessingScript                                   Path to Transcript_to_genome.R script
	--bulkLevelScript                                        Path to Calculate_m6anet_bulk.R script

	""".stripIndent()
}

// Show help message
if (params.help) {
	helpMessage()
	exit 0
}

// Input of transcriptome fasta.
Channel
	.fromPath(params.transcriptome_fasta, checkIfExists:true)
	.into{transcriptome_fasta_minimap2;transcriptome_fasta_nanopolish}

// Input of sample names, conditions, FAST5s path and FASTQ.
Channel
	.fromPath( params.samples )
	.splitCsv(header: true, sep:'\t')
	.map{ row-> tuple(row.sample, row.condition, file(row.fast5_dir), file(row.fastq)) }
	.set{ samples_minimap2 }

// Transcriptome alignment.
process minimap2 {
	input:
	tuple val(sample), val(condition), val(fast5_dir), val(fastq) from samples_minimap2

	each file('transcriptome.fa') from transcriptome_fasta_minimap2
		
	output:
	tuple val(sample), val(condition), val(fast5_dir), val(fastq) into minimap2_nanopolish

    script:
    if(params.minimap2)
	"""
		mkdir -p ${params.resultsDir}/${condition}/${sample}/transcriptomeAlignment/
		minimap2 -x map-ont -k14 -t ${task.cpus} -a transcriptome.fa ${fastq} | samtools view -hSb | samtools sort -@ ${task.cpus} -o ${params.resultsDir}/${condition}/${sample}/transcriptomeAlignment/minimapT.bam
		samtools view ${params.resultsDir}/${condition}/${sample}/transcriptomeAlignment/minimapT.bam -bh -F 2324 | samtools sort -@ ${task.cpus} -o ${params.resultsDir}/${condition}/${sample}/transcriptomeAlignment/minimap.filt.sortT.bam	
		samtools index -@ ${task.cpus} ${params.resultsDir}/${condition}/${sample}/transcriptomeAlignment/minimap.filt.sortT.bam
		ln -s ${params.resultsDir}/${condition}/${sample}/transcriptomeAlignment/minimap.filt.sort.bam ./minimap.filt.sortT.bam
		ln -s ${params.resultsDir}/${condition}/${sample}/transcriptomeAlignment/minimap.filt.sort.bam.bai ./minimap.filt.sortT.bam.bai
	"""
	else
	"""
		ln -s ${params.resultsDir}/${condition}/${sample}/transcriptomeAlignment/minimap.filt.sort.bam ./minimap.filt.sortT.bam
		ln -s ${params.resultsDir}/${condition}/${sample}/transcriptomeAlignment/minimap.filt.sort.bam.bai ./minimap.filt.sortT.bam.bai
	"""
}


// Resquigling with nanopolish for each condition
process nanopolish {
	input:
	tuple val(sample), val(condition), val(fast5_dir), val(fastq) from minimap2_nanopolish
	each file('transcriptome.fa') from transcriptome_fasta_nanopolish

	output:
	tuple val(condition), val(sample) into nanopolish_m6anet1

    script:
    if(params.nanopolish)
	"""
		mkdir -p ${params.resultsDir}/${condition}/${sample}/nanopolish/
		seqsum_file=\$(dirname ${fast5_dir})/sequencing_summary.txt

        	if [[ -f "\$seqsum_file" ]]; then
			f5c index -d ${fast5_dir} ${fastq} -s \$seqsum_file
		else
			f5c index -d ${fast5_dir} ${fastq}
        	fi
		f5c eventalign --rna -r ${fastq} -b ${params.resultsDir}/${condition}/${sample}/transcriptomeAlignment/minimap.filt.sortT.bam -g transcriptome.fa --signal-index --scale-events --summary ${params.resultsDir}/${condition}/${sample}/nanopolish/summary.txt -t ${task.cpus} > ${params.resultsDir}/${condition}/${sample}/nanopolish/eventalign_readIndex.txt
	"""
	else
	"""
        	echo "Skipped"
	"""
}


// Data formatting for m6anet for each sample
process m6anet1 {
	input:
	tuple val(condition), val(sample) from nanopolish_m6anet1

	output:
	tuple val(condition), val(sample) into m6anet1_m6anet2


    script:
	if(params.m6anet1)
	"""
		mkdir -p ${params.resultsDir}/${condition}/${sample}/m6anet/
		m6anet dataprep --eventalign  ${params.resultsDir}/${condition}/${sample}/nanopolish/eventalign_readIndex.txt \
		--out_dir ${params.resultsDir}/${condition}/${sample}/m6anet --n_processes ${task.cpus}
	"""
	else
	"""
		ln -sf ${params.resultsDir}/${condition}/${sample}/m6anet m6anet
	"""
}


// From a single channel for all the alignments to one channel for each condition
m6anet1_m6anet2.groupTuple(by:0)
.set{ m6anet1_m6anet2_grouped_by_condition } 

// RNA modifications detection with m6anet
process m6anet2 {
	input:
	tuple val(condition), val(sample) from m6anet1_m6anet2_grouped_by_condition

	output:
    	val(condition) into m6anet_postprocessing
	script:
	if(params.m6anet2)
	"""
		mkdir -p ${params.resultsDir}/${condition}/m6anet
		preprocessing_dirs=\$(find ${params.resultsDir}/${condition} -maxdepth 2 -mindepth 2 -type d | grep "m6anet\$")
		m6anet inference --input_dir \$preprocessing_dirs --out_dir ${params.resultsDir}/${condition}/m6anet --n_processes ${task.cpus}
	"""
	else
	"""
		echo "Skipped"
	"""
}


// Processing of each output to obtain bed files
process postprocessing {
	input:
	val(condition) from m6anet_postprocessing

	output:

	script:
	if(params.postprocessing)
	"""
		mkdir -p ${params.resultsDir}/${condition}/m6anet_postprocessing

		Rscript ${params.postprocessingScript} \
			input_file=${params.resultsDir}/${condition}/m6anet/data.site_proba.csv \
			output_file=${params.resultsDir}/${condition}/m6anet_postprocessing/data.site_proba.thr${params.prob_mod_thr}.tsv \
			output_file_genome=${params.resultsDir}/${condition}/m6anet_postprocessing/data.site_proba.genome.thr${params.prob_mod_thr}.tsv \
			genome_gtf=${params.gtf} \
			mccores=${task.cpus} \
			prob_mod_thr=${params.prob_mod_thr}

		Rscript ${params.bulkLevelScript} \
			m6anet_output_file=${params.resultsDir}/${condition}/m6anet/data.site_proba.csv \
			report_file=${params.resultsDir}/${condition}/m6anet_postprocessing/m6A_bulk_level_estimate.txt \
			prob_mod_thr=${params.prob_mod_thr}
	"""
	else
	"""
		echo "Skipped"
	"""
}
