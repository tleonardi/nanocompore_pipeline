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
        .set{ genome_fasta}
}
else {
    exit 1, "No genome fasta file specified!"
}

// Setup input channel for target transcript list
if( params.target_trancripts){
	bed_filter = file(params.target_trancripts)
}
else{
	bed_filter = file("$baseDir/assets/NO_FILE")
}

// Setup input channel for inverse filter 
if( params.exclude_trancripts){
        bed_invfilter = file(params.exclude_trancripts)
}
else{
        bed_invfilter = file("$baseDir/assets/NO_FILE2")
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
      .into{guppy_outputs_pycoqc; guppy_outputs_minimap; guppy_outputs_nanopolish}
}
else{
  Channel
      .fromPath( params.samples )
      .splitCsv(header: true, sep:'\t')
      .map{ row-> tuple(row.SampleName, row.Condition, file(row.DataPath)) }
      .into{guppy_annot; nanopolish_annot}
  
  process guppy {
    publishDir "${params.resultsDir}/${sample}", mode: 'copy'
    input:
      set val(sample),val(condition),file(fast5) from guppy_annot
    output:
      set val("${sample}"), file("guppy") into guppy_outputs_pycoqc, guppy_outputs_minimap, guppy_outputs_nanopolish
    
    script:
      def keep_fast5 = params.keep_basecalled_fast5  ? "--fast5_out" : ""
      def gpu_opts = ""
      if (params.GPU == "true") {
        gpu_opts = "-x 'cuda:0' --gpu_runners_per_device ${params.guppy_runners_per_device} --chunks_per_runner ${params.guppy_chunks_per_runner} --chunk_size ${params.guppy_chunk_size}"
      }
    """
    guppy_basecaller -i ${fast5} -s guppy  ${keep_fast5} ${gpu_opts} --recursive --num_callers ${task.cpus} --disable_pings --reverse_sequence true --u_substitution true --trim_strategy rna --flowcell ${params.flowcell} --kit ${params.kit}
    """
  }
}

// QC Albacore's output
process pycoQC {
  publishDir "${params.resultsDir}/${sample}", mode: 'copy'
  input:
    set val(sample),file(guppy_results) from guppy_outputs_pycoqc
  output:
    file "pycoqc.html" into pycoqc_outputs
  when:
    params.qc==true

  """
  pycoQC -f "${guppy_results}/sequencing_summary.txt" -o pycoqc.html --min_pass_qual 7 
  """
}

// Prepare BED and fasta annotation files
process prepare_annots {
  publishDir "${params.resultsDir}/references/", mode: 'copy'
  input:
    file transcriptome_gtf
    file genome_fasta
    file bed_filter
    file bed_invfilter
  output:
    file "reference_transcriptome.bed" into transcriptome_bed
    file "${genome_fasta}.fai" into genome_fai
    file "reference_transcriptome.fa" into transcriptome_fasta_minimap, transcriptome_fasta_nanopolish, transcriptome_fasta_nanocompore
    file "reference_transcriptome.fa.fai" into transcriptome_fai_minimap

  script:
    def filter = bed_filter.name != 'NO_FILE' ? "| bedparse filter --annotation !{bed_filter}" : ''
    def inv_filter = bed_invfilter.name != 'NO_FILE2' ? "| bedparse filter --annotation ${bed_invfilter} -v " : ''
  """
  bedparse gtf2bed ${transcriptome_gtf} ${filter} ${inv_filter} | awk 'BEGIN{OFS=FS="\t"}{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12}' > reference_transcriptome.bed
  bedtools getfasta -fi ${genome_fasta} -s -split -name -bed reference_transcriptome.bed -fo - | perl -pe 's/>(.+)::.+/>\$1/' > reference_transcriptome.fa
  samtools faidx reference_transcriptome.fa
 """
}

// Map the basecalled data to the reference with Minimap2
process minimap {
  publishDir "${params.resultsDir}/${sample}/", mode: 'copy'
  input:
    set val(sample),file(guppy_results) from guppy_outputs_minimap
    each file('transcriptome.fa') from transcriptome_fasta_minimap
    each file('transcriptome.fa.fai') from transcriptome_fai_minimap
  output:
    set val(sample), file("minimap.filt.sort.bam"), file("minimap.filt.sort.bam.bai") into minimap

script:
def mem = task.mem ? " -m ${(task.mem.toBytes()/1000000).trunc(0) - 1000}M" : ''
"""
	minimap2 -x map-ont -t ${task.cpus} -a transcriptome.fa ${guppy_results}/*.fastq > minimap.sam
	samtools view minimap.sam -bh -t transcriptome.fa.fai -F 2324 | samtools sort -@ ${task.cpus} ${mem} -o minimap.filt.sort.bam
	samtools index minimap.filt.sort.bam minimap.filt.sort.bam.bai
"""  
}

// Run nanopolish and nanopolishComp
process nanopolish {
  publishDir "${params.resultsDir}/${sample}/", mode: 'copy'
  input:
    // The raw data file has to have a fixed name to avoid a collision with guppy_results filename
    set val(sample), file(guppy_results), val(label), file('raw_data'), file(bam_file), file(bam_index) from guppy_outputs_nanopolish.join(nanopolish_annot).join(minimap)
    each file(transcriptome_fasta) from transcriptome_fasta_nanopolish
  output: 
    set val(sample), val(label), file("nanopolishcomp/out_eventalign_collapse.tsv"), file("nanopolishcomp/out_eventalign_collapse.tsv.idx") into nanopolishComp


script:
def cpus_each = task.cpus
"""
	cat ${guppy_results}/*.fastq > basecalled.fastq
	nanopolish index -s ${guppy_results}/sequencing_summary.txt -d 'raw_data' basecalled.fastq
	mkdir -p nanopolishcomp/
        nanopolish eventalign -t ${cpus_each} --reads basecalled.fastq --bam ${bam_file} --genome ${transcriptome_fasta} --samples --print-read-names --scale-events | NanopolishComp Eventalign_collapse -t ${cpus_each} -o nanopolishcomp/
"""
}



ni_ref=Channel.create()
ni_other=Channel.create()
nanopolishComp.groupTuple(by:1).choice( ni_ref, ni_other ) { a -> a[1] == params.reference_condition ? 0 : 1 }

process nanocompore {
  echo true
  publishDir "${params.resultsDir}/", mode: 'copy'
  input:
    set val(labels1), val(condition1), file('npcomp_ref*.tsv'), file('npcomp_ref*.tsv.idx') from ni_ref
    set val(labels2), val(condition2), file('npcomp_kd*.tsv'), file('npcomp_kd*.tsv.idx') from ni_other
    file transcriptome_fasta_nanocompore
    file transcriptome_bed
  output:
    file(nanocompore)

def downsample = params.downsample_high_cov ? " --downsample_high_coverage ${params.downsample_high_cov}" : " "
shell:
"""
IFS=','
f1=(npcomp_ref*.tsv)
f2=(npcomp_kd*.tsv)
nanocompore sampcomp --file_list1 "\${f1[*]}" --file_list2 "\${f2[*]}" \
 --label1 ${condition1} \
 --label2 ${condition2} \
 --fasta ${transcriptome_fasta_nanocompore} \
 --bed ${transcriptome_bed}\
 --outpath nanocompore \
 --sequence_context ${params.sequenceContext} \
 ${downsample} \
 --allow_warnings \
 --pvalue_thr ${params.pvalue_thr} \
 --min_coverage ${params.min_cov} \
 --logit \
 --nthreads ${task.cpus}
"""
}
