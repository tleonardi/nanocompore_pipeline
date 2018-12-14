#!/usr/bin/env nextflow

// Set up input channel from GTF annotation
if( params.gtf ){
    Channel
        .fromPath(params.gtf, checkIfExists:true)
        .set{ transcriptome_gtf }
}
else {
    exit 1, "No GTF annotation specified!"
}

// Set up input channel from Fasta file
if( params.fasta ){
    Channel
        .fromPath(params.fasta, checkIfExists:true)
        .set { genome_fasta }
}
else {
    exit 1, "No genome fasta file specified!"
}

// Setup input channel for target transcript list
if( params.target_trancripts){
	bed_filter = file(params.target_trancripts)
}
else{
	bed_filter = file("NO_FILE")
}

/* If the input paths are already basecalled
   define Albacore's output channels otherwise,
   execute Albacore
*/
if(params.input_is_basecalled){
  Channel
      .fromPath( params.samples )
      .splitCsv(header: true, sep:'\t')
      .map{ row-> tuple(row.SampleName, row.Condition, file(row.DataPath)) }
      .set{nanopolish_annot}

  Channel
      .fromPath( params.samples )
      .splitCsv(header: true, sep:'\t')
      .map{ row-> tuple(row.SampleName, file(row.DataPath)) }
      .into{albacore_outputs_pycoqc; albacore_outputs_minimap; albacore_outputs_nanopolish}
}
else{
  Channel
      .fromPath( params.samples )
      .splitCsv(header: true, sep:'\t')
      .map{ row-> tuple(row.SampleName, row.Condition, file(row.DataPath)) }
      .into{albacore_annot; nanopolish_annot}
  
  process albacore {
    publishDir "${params.resultsDir}/${sample}", mode: 'copy'
    input:
      set val(sample),val(condition),file(fast5) from albacore_annot
    output:
      set val("${sample}"), file("albacore") into albacore_outputs_pycoqc, albacore_outputs_minimap, albacore_outputs_nanopolish
    
    script:
      def outformat = params.keep_basecalled_fast5  ? "fastq,fast5" : "fastq"
    """
    read_fast5_basecaller.py -r -i ${fast5} -t ${task.cpus} -s albacore -f "FLO-MIN106" -k "SQK-RNA001" -o ${outformat} -q 0 --disable_pings --disable_filtering
    """
  }
}

// QC Albacore's output
process pycoQC {
  publishDir "${params.resultsDir}/${sample}", mode: 'copy'
  input:
    set val(sample),file(albacore_results) from albacore_outputs_pycoqc
  output:
    file "pycoqc.html" into pycoqc_outputs
  when:
    params.qc==true

  """
  echo pycoQC -f "${albacore_results}/sequencing_summary.txt" -o pycoqc.html --min_pass_qual 7 > pycoqc.html
  """
}

// Prepare BED and fasta annotation files
process prepare_annots {
  publishDir "${params.resultsDir}/references/", mode: 'copy'
  input:
    file transcriptome_gtf
    file genome_fasta
    file bed_filter
  output:
    file "reference_transcriptome.bed" into transcriptome_bed
    file "reference_transcriptome_fastaName.bed" into transcriptome_bed_faname
    file "reference_transcriptome.fa" into transcriptome_fasta_minimap, transcriptome_fasta_nanopolish
    file "reference_transcriptome.fa.fai" into transcriptome_fai_minimap

  script:
    def filter = bed_filter.name != 'NO_FILE' ? "| bedparse filter --annotation !{bed_filter}" : ''
  """
  bedparse gtf2bed --extraFields transcript_name ${transcriptome_gtf} ${filter} > reference_transcriptome.bed

  awk 'BEGIN{OFS=FS="\t"}{print \$1,\$2,\$3,\$4"::"\$1":"\$2"-"\$3"("\$6")",\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12}' reference_transcriptome.bed > reference_transcriptome_fastaName.bed
  bedtools getfasta -fi ${genome_fasta} -s -split -name -bed reference_transcriptome.bed > reference_transcriptome.fa
  samtools faidx reference_transcriptome.fa
 """
}

// Map the basecalled data to the reference with Minimap2
process minimap {
  publishDir "${params.resultsDir}/${sample}/", mode: 'copy'
  input:
    set val(sample),file(albacore_results) from albacore_outputs_minimap
    each file('transcriptome.fa') from transcriptome_fasta_minimap
    each file('transcriptome.fa.fai') from transcriptome_fai_minimap
  output:
    set val(sample), file("minimap.filt.sort.bam"), file("minimap.filt.sort.bam.bai") into minimap


"""
	minimap2 -x map-ont -t ${task.cpus} -a transcriptome.fa ${albacore_results}/workspace/*.fastq > minimap.sam
	samtools view minimap.sam -bh -t transcriptome.fa.fai -F 2324 | samtools sort -@ ${task.cpus} -m 15G -o minimap.filt.sort.bam
	samtools index minimap.filt.sort.bam minimap.filt.sort.bam.bai
"""  
}

// Run nanopolish and nanopolishComp
process nanopolish {
  publishDir "${params.resultsDir}/${sample}/", mode: 'copy'
  input:
    // The raw data file has to have a fixed name to avoid a collision with albacore_results filename
    set val(sample), file(albacore_results), val(label), file('raw_data'), file(bam_file), file(bam_index) from albacore_outputs_nanopolish.join(nanopolish_annot).join(minimap)
    each file(transcriptome_fasta) from transcriptome_fasta_nanopolish
  output: 
    file("reads_collapsed.tsv")
    file("reads_collapsed.tsv.idx")


script:
def cpus_each = (task.cpus/2).trunc(0)
"""
	nanopolish index -s ${albacore_results}/sequencing_summary.txt -d 'raw_data' ${albacore_results}/workspace/*.fastq
	nanopolish eventalign -t ${cpus_each} --reads ${albacore_results}/workspace/*.fastq --bam ${bam_file} --genome ${transcriptome_fasta} --samples --print-read-names --scale-events | NanopolishComp Eventalign_collapse -t ${cpus_each} -o reads_collapsed.tsv
"""
}


