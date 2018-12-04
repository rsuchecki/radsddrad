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
    file "*" into references

  script:
    outname = url.replaceAll(".*/","")
    // if(params.trialLines == null) {
      """
      curl $url > ${outname}
      """
    // } else {
    //   outname = "${params.trialLines}_trialLines_${outname}"
    //   """
    //   curl $url | gunzip --stdout | head -n ${params.trialLines} | pigz -cp ${task.cpus} > ${outname}
    //   """
    // }
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
    set val(idxmeta), file("*"), val(sample), file(reads) from indices.combine(readPairs)

  output:
    set val(meta), file("*") into alignedDatasets

  script:
    meta = idxmeta.clone() + [sample: sample]
    r1 = reads[0]
    r2 = reads[1]
    template "${idxmeta.tool}_align.sh"  //points to e.g. biokanga_align.sh in templates/
}

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
