#!/usr/bin/env nextflow

//RETURNS ALIGNER NAMES/LABELS IF BOTH INDEXING AND ALIGNMENT TEMPLATES PRESENT
aligners = Channel.fromFilePairs("${workflow.projectDir}/templates/*_{index,align}.sh", maxDepth: 1, checkIfExists: true)
  .map { it[0] }
  .filter{ !params.debug || it.matches("(biokanga)\$") }


readPairs = Channel
  .fromFilePairs('split/*_R{1,2}.fastq.gz')
  // .println()


url = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-41/fasta/hordeum_vulgare/dna/Hordeum_vulgare.IBSC_v2.dna.toplevel.fa.gz'

process fetchRemoteReference {
  storeDir {executor == 'awsbatch' ? "${params.outdir}/downloaded" : "downloaded"} // storeDir "${workflow.workDir}/downloaded" put the datasets there and prevent generating cost to dataset creators through repeated downloads on re-runs
  scratch false
  tag{outname}
  label 'download'

  input:
    url

  output:
    file "${outname}" into references, refs4regions, refs4calls

  script:
    outname = url.replaceAll(".*/","").replaceAll("\\.gz\$","")
    if(url.endsWith(".gz")) {
      """
      curl $url | gunzip --stdout  > ${outname}
      """
    } else {
      """
      curl $url > ${outname}
      """
    }
}


process index {
  label 'index'
  container { this.config.process.get("withLabel:${tool}" as String).get("container") }
  tag("${tool} << ${ref}")

  input:
    set val (tool), file(ref) from aligners.combine(references)

  output:
    set val(meta), file("*") into indices

  script:
    meta = [tool: "${tool}", target: "${ref}"]
    template "${tool}_index.sh" //points to e.g. biokanga_index.sh in templates/
}

process align {
  label 'align'
  label 'out'

  container { this.config.process.get("withLabel:${idxmeta.tool}" as String).get("container") }
  tag("${meta}")

  input:
    set val(idxmeta), file("*"), val(sample), file(reads) from indices.combine(readPairs).take(10)

  output:
    set val(meta), file("*") into alignedDatasets
    file("*.bam") into BAMs

  script:
    meta = idxmeta.clone() + [sample: sample]
    r1 = reads[0]
    r2 = reads[1]
    template "${idxmeta.tool}_align.sh"  //points to e.g. biokanga_align.sh in templates/
}

BAMs
  .map { it.path }
  .collectFile(newLine: true)
  .set { BAMfofn }

// BAMfofn.subscribe {println "Entries are saved to file: $it"}


// refs4regions
//   .splitFasta( record: [id: true], decompress: true)
//   .subscribe { record -> println record.id }

process refIds {

  input:
    file ref from refs4regions

  output:
    stdout into refIds

  // exec:
  //   println ref
  //   println ref.name.endsWith(".gz")

  script:
  CAT = ref.name.endsWith(".gz") ? "zcat" : "cat"
  """
  ${CAT} ${ref} | grep '^>' | cut -d ' ' -f1 | sed 's/>//'
  """
}

// process recompressRef {
//   echo true
//   module 'samtools'
//   input:
//     file(ref) from refs4calls

//   output:
//     file('*') into bgzRef

//   script:
//   out = ref.name.replaceAll("\\.gz\$",".bgz")
//   """
//   gzip -dc ${ref} | bgzip -c > ${out}
//   """
// }



// // refIds
// //   .splitText()
// //   .subscribe { print it }
process pileupAndCall {
  label 'bcftools'
  // module 'bcftools/1.9.0'
  echo true
  tag "${id}"

  input:
    set val(id), file(bamFOFN), file(ref) from refIds.splitText().map{ it.trim() }.combine(BAMfofn).combine(refs4calls)


  // output:
  //   stdout into refIds

  // exec:
  //   // id = comb[0].trim()
  //   // // bamFiles = comb[(1..-2)].join("\")
  //   // ref = comb[-1]
  //   // print comb
  //   println id
  //   println bamFOFN
  //   println ref
  //   // println comb[(1..-2)]
  //   // println comb[-1]
  //   // println comb[(1..-1)]
  //   // println id.trim()
  //   // println meta
  //   // println files
  // //   println ref.name.endsWith(".gz")

  script:
  // """
  // ls -la
  """
  bcftools mpileup \
    --regions ${id}  \
    --bam-list ${bamFOFN} \
    --output-type z \
    --skip-indels \
    --annotate AD,DP \
    --fasta-ref ${ref} \
    --min-MQ 20 \
    --min-BQ 20  \
    --no-version \
    -o mpileup_${id}.vcf.gz;\
  """
  // #bcftools call --multiallelic-caller --variants-only --no-version Intermediate_files/3.mpileup/mpileup_{}.vcf.gz | sed -e 's|$(pwd)\/||g' -e 's/Intermediate_files\/2\.bam_alignments\///g' -e  's/\.R.\.fastq.sorted_bam//g'  > Intermediate_files/4.Raw_SNPs/raw_SNPs_{}.vcf;\

  // """
}
// parallel  --gnu --max-procs $NUM_CORES --keep-order "\

// bcftools mpileup --regions {} --output-type z --skip-indels --annotate AD,DP --fasta-ref $REF_GENOME --min-MQ 20 --min-BQ 20  --no-version -b Intermediate_files/2.bam_alignments/samples_list.txt -o Intermediate_files/3.mpileup/mpileup_{}.vcf.gz;\

// bcftools call --multiallelic-caller --variants-only --no-version Intermediate_files/3.mpileup/mpileup_{}.vcf.gz | sed -e 's|$(pwd)\/||g' -e 's/Intermediate_files\/2\.bam_alignments\///g' -e  's/\.R.\.fastq.sorted_bam//g'  > Intermediate_files/4.Raw_SNPs/raw_SNPs_{}.vcf;\

// " ::: `grep ">" $REF_GENOME | cut -d ' ' -f1 | sed 's/>//g'`
// process callVariants {

//   script:
//   """
//     Platypus.py callVariants --bamFiles="${listOfBAMs}" \
//       --bufferSize 1000000 \
// 	    --nCPU="${task.cpus}" \
// 	    --trimAdapter=0 --maxGOF=20 \
// 	    --minReads 5 --genIndels=1 --minFlank=5 \
// 	    --sbThreshold=0.01 --scThreshold=0.95 --hapScoreThreshold=15 \
// 	    --filterDuplicates=0 \

// 	    --refFile ${ref} \
// 	    --logFileName ${out}.txt \
// 	    --output ${out}.vcf
//   """
// }

// process align {
//   echo true
//   input:
//     set val(sample), file(reads) from readPairs

//   script:
//   """
//   ls -la
//   """
// }


// keys = Channel
//     .fromPath('meta/runs.csv')
//     .splitCsv(header:true)
//     .subscribe{ println it}

// reads = Channel
//     .fromFilePairs('data/Hordeum_run*/FASTQ_Generation_*/*/*R?_001.fastq.gz')
//     .subscribe { println it}
// reads =

// // runs = Channel.fromPath( 'meta/runs.tsv' )
// Channel.from(file('meta/runs.csv').text)
//   // .splitText()
//   // .splitCSV( header: true )
//   // .take(3)
//   .first()
//   .splitCSV( header:true )
//   // .subscribe { print it }
//   .subscribe { print it.getClass()}
//   // .subscribe { print it.getClass().methods}

// runs = Channel.from(1..6)
// indices = Channel.from(1..3)

// process splitReads {
//   input:
//     set val(run), val(index), val(key) from runs.combine(indices)

//   when:

//   exec:
//     println "run "+run+" index "+index
//   script:
//   """

//   """
// }

// /*
//   Generic method for extracting a string tag or a file basename from a metadata map
//  */
// def getDatasetTagFromMeta(meta, delim = '_') {
//   return meta.species+delim+meta.version+(params.trialLines == null ? "" : delim+params.trialLines+delim+"trialLines")
// }
