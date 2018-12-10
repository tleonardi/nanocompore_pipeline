#!/usr/bin/env nextflow
 
params.samples = "$baseDir/samples.txt"
params.gtf = "/home/nanopore/references/Homo_sapiens.GRCh38.93.gtf"
params.fasta = "/home/nanopore/references/Homo_sapiens.GRCh38.dna.toplevel.fa"
params.target_trancripts = "targets"

params.minimap2="/home/nanopore/.local/src/minimap2"
params.samtools="/home/nanopore/.local/src/samtools-1.9/samtools" 

if( params.gtf ){
    Channel
        .fromPath(params.gtf, checkIfExists:true)
        .into { transcriptome_gtf }
}
else {
    exit 1, "No GTF annotation specified!"
}

if( params.fasta ){
    Channel
        .fromPath(params.fasta, checkIfExists:true)
        .into { genome_fasta }
}
else {
    exit 1, "No genome fasta file specified!"
}

if( params.target_trancripts ){
    Channel
        .fromPath(params.target_trancripts, checkIfExists:true)
        .into { bed_filter }
}
else {
    bed_filter = "NA"
}
Channel
    .fromPath( params.samples )
    .splitCsv(header: true, sep:'\t')
    .map{ row-> tuple(row.SampleName, row.Condition, file(row.DataPath)) }
    .into{albacore_annot; another_annot}

process albacore {
  publishDir "$baseDir/out/${sample}", mode: 'copy'
  input:
    set val(sample),val(condition),file(fast5) from albacore_annot
  output:
    set val("${sample}"), file("albacore") into albacore_outputs_pycoqc, albacore_outputs_minimap

  """
  VIRTUAL_ENV_DISABLE_PROMPT=true
  source /home/nanopore/DATA/nanopore_7SK_complete_analysis/virtualenv/bin/activate
  read_fast5_basecaller.py -r -i ${fast5} -t 12 -s albacore -f "FLO-MIN106" -k "SQK-RNA001" -o fastq,fast5 -q 0 --disable_pings --disable_filtering

  """
}

process pycoQC {
  publishDir "$baseDir/out/${sample}", mode: 'copy'
  input:
    set val(sample),file(albacore_results) from albacore_outputs_pycoqc
  output:
    file "pycoqc.html" into pycoqc_outputs

  """
  echo pycoQC -f "${albacore_results}/sequencing_summary.txt" -o pycoqc.html --min_pass_qual 7 > pycoqc.html
  """
}

process prepare_annots {
  publishDir "$baseDir/out/references/", mode: 'copy'
  input:
    file transcriptome_gtf
    file genome_fasta
    file bed_filter
  output:
    file "reference_transcriptome.bed" into transcriptome_bed
    file "reference_transcriptome_fastaName.bed" into transcriptome_bed_faname
    file "reference_transcriptome.fa" into transcriptome_fasta

  shell:
  '''
  if [[ -f !{bed_filter} ]]; then 
	  bedparse gtf2bed --extraFields transcript_name !{transcriptome_gtf} | bedparse filter --annotation !{bed_filter} > reference_transcriptome.bed
  else
	  bedparse gtf2bed --extraFields transcript_name !{transcriptome_gtf} > reference_transcriptome.bed
  fi

  awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,$4"::"$1":"$2"-"$3"("$6")",$5,$6,$7,$8,$9,$10,$11,$12}' reference_transcriptome.bed > reference_transcriptome_fastaName.bed
  bedtools getfasta -fi !{genome_fasta} -s -split -name -bed reference_transcriptome.bed > reference_transcriptome.fa
 ''' 
}

process map {
  publishDir "$baseDir/out/${sample}/", mode: 'copy'
  input:
    file transcriptome_fasta
    set val(sample),file(albacore_results) from albacore_outputs_minimap
  output:
    file "minimap.sam"
    file "minimap.filt.sort.bam"
    file "minimap.filt.sort.bam.bai"
    file "map"


"""
	echo "map" > map
	${params.minimap2} -ax map-ont ${transcriptome_fasta} ${albacore_results}/workspace/*.fastq > minimap.sam
	${params.samtools} view minimap.sam -bh -F 2324 | ${params.samtools} sort -o minimap.filt.sort.bam
	${params.samtools} index minimap.filt.sort.bam minimap.filt.sort.bam.bai
"""  
}
