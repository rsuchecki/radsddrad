#!/usr/bin/env nextflow

//RETURNS ALIGNER NAMES/LABELS IF BOTH INDEXING AND ALIGNMENT TEMPLATES PRESENT
aligners = Channel.fromFilePairs("${workflow.projectDir}/templates/*_{index,align}.sh", maxDepth: 1, checkIfExists: true)
  .map { it[0] }
  .filter{ !params.debug || it.matches("(biokanga)\$") }


readPairs = Channel
  .fromFilePairs('split/*_R{1,2}.fastq.gz')
  // .println()

//INPUT PARAMS
trialLines = params.trialLines

url = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-41/fasta/hordeum_vulgare/dna/Hordeum_vulgare.IBSC_v2.dna.toplevel.fa.gz'


process fetchRemoteReference {
  storeDir {executor == 'awsbatch' ? "${params.outdir}/downloaded" : "downloaded"} // storeDir "${workflow.workDir}/downloaded" put the datasets there and prevent generating cost to dataset creators through repeated downloads on re-runs
  scratch false
  tag{outname}
  label 'download'

  input:
    url

  output:
    file ref // into references

  script:
    outname = (trialLines == null ? "" : trialLines+"_trialLines_")+url.replaceAll(".*/","").replaceAll("\\.gz\$","")
    ref = "${outname}"
    TLINES = params.trialLines == null ? "" : "| head -n ${trialLines}"
    GUNZIP = url.endsWith(".gz") ? "| gunzip --stdout" : ""
    """
    curl $url ${GUNZIP} ${TLINES}  > ${outname}
    """
}

//Duplicate refs channel for 4 different downstream processes
// refChannelsQueue = references.into(4) as Queue


process index {
  //Only use store dir with explicit fileNames otherwise exec skipped incorrectly
  //storeDir {executor == 'awsbatch' ? "${workDir}/indices/${tool}" : "${workDir}/indices/${tool}"}  //hmm... could go via params to allow user-defined location?
  label 'index'
  module { this.config.process.get("withLabel:${tool}" as String).get("module") }
  container { this.config.process.get("withLabel:${tool}" as String).get("container") }
  tag("${tool} << ${ref}")

  input:
    // set val (tool), file(ref) from aligners.combine(refChannelsQueue.poll())
    val tool from aligners
    file ref

  output:
    set val(meta), file("*") into indices

  script:
    meta = [tool: "${tool}", target: "${ref}"]
    template "${tool}_index.sh" //points to e.g. biokanga_index.sh in templates/
}

process align {
  label 'align'
  // label 'out'

  module { this.config.process.get("withLabel:${idxmeta.tool}" as String).get("module") }
  container { this.config.process.get("withLabel:${idxmeta.tool}" as String).get("container") }
  tag("${meta}")

  input:
    set val(idxmeta), file("*"), val(sample), file(reads) from indices.combine(readPairs).take(4)

  output:
    //set val(meta), file("*") into alignedDatasets //don't use '*' with scratch true
    set val(meta), file("*.?am") into BAMs //don't use '*' with scratch true
    // file("*.bam") into BAMs

  script:
    meta = idxmeta.clone() + [sample: sample, sorted: idxmeta.tool=='biokanga' ? true : false]
    r1 = reads[0]
    r2 = reads[1]
    template "${idxmeta.tool}_align.sh"  //points to e.g. biokanga_align.sh in templates/
}


BAMfofn = BAMs.map { meta, bam -> bam.path }.collectFile(newLine: true)
  //.toSortedList( { a,b -> a[0].sample <=> b[0].sample }) //.map { meta, bam -> bam.path }.collectFile(newLine: true)
  // .flatMap { meta,bam -> bam.path }
  // .subscribe { println it }
  // .subscribe { println it[0].sample }

// BAMfofn.subscribe {println "Entries are saved to file: $it"}


// refs4regions
//   .splitFasta( record: [id: true], decompress: true)
//   .subscribe { record -> println record.id }

process refFaidx {
  tag("${ref}")
  label 'samtools'
  input:
    file ref //from refChannelsQueue.poll()

  output:
    file fai //into faIndices

  script:
  fai = "${ref}.fai"
  """
  samtools faidx ${ref}
  """
}

process refIds {
  tag("${ref}")
  input:
    file ref //from refChannelsQueue.poll()

  output:
    file refIds

  // exec:
  //   println ref
  //   println ref.name.endsWith(".gz")

  script:
  cat = ref.name.endsWith(".gz") ? "zcat" : "cat"
  """
  ${cat} ${ref} | grep '^>' | cut -d ' ' -f1 | sed 's/>//'> refIds
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





// REFIDS.subscribe { println it }


process pileupAndCall {
  label 'bcftools'
  // module 'bcftools/1.9.0'
  echo true
  tag "${id}"

  input:
    file BAMfofn
    file fai
    file ref
    // set val(id), file(bamFOFN), file(ref), file(fai) from refIds.splitText().map{ it.trim() } \
    each id from refIds.splitText().map{ it.trim() } //(refIds.splitText().map{ it.trim() })
    // set val(id), file(ref), file(fai) from refIds.splitText().map{ it.trim() } \
                                                                  // .combine(BAMfofn) \
                                                                  // .combine(refChannelsQueue.poll()) \
                                                                  // .combine(faIndices)


  // output:
  //   file('*.vcf') into vcfs

  // exec:
  //   // id = comb[0].trim()
  //   // // bamFiles = comb[(1..-2)].join("\")
  //   // ref = comb[-1]
  //   // print comb
    // println id
    // println bamFOFN
    // println ref
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
    --targets ${id}  \
    --bam-list ${BAMfofn} \
    --skip-indels \
    --annotate AD,DP \
    --fasta-ref ${ref} \
    --min-MQ 20 \
    --min-BQ 20  \
    --no-version \
  | bcftools call \
    --multiallelic-caller \
    --variants-only \
    - \
  | awk 'BEGIN{FS=OFS="\\t"};{if(\$1=="#CHROM"){for(i=10;i<=NF;i++){gsub(/.*\\//,"",\$i);gsub(/\\.bam\$/,"",\$i)}};print}' \
    > ${id}.vcf
  """
  // #bcftools call --multiallelic-caller --variants-only --no-version Intermediate_files/3.mpileup/mpileup_{}.vcf.gz | sed -e 's|$(pwd)\/||g' -e 's/Intermediate_files\/2\.bam_alignments\///g' -e  's/\.R.\.fastq.sorted_bam//g'  > Intermediate_files/4.Raw_SNPs/raw_SNPs_{}.vcf;\
}

// process combineRegions {
//   label 'bcftools'
//   input:
//     file('*') from vcfs.collect()

//   output:
//     file combined

//   script:
//   combined = 'combined.vcf'
//   """
//   bcftools concat --no-version \$(ls -v *.vcf) > ${combined}
//   """
// }

// thresholds = Channel.from((1..50).findAll {it % 5 == 0 })
// process filterSNPs {
//   label 'vcftools'
//   label 'out'
//   tag("minDepth = ${minDepth}")

//   input:
//     file combined
//     val minDepth from thresholds

//   output:
//     file '*.vcf'

//   script:
//   """
//   vcftools \
//     --vcf ${combined} \
//     --minDP ${minDepth} \
//     --recode \
//     --stdout \
//     | awk '/#/ || /[0-9]\\/[0-9]/' \
//     > combined_minDepth_${minDepth}.vcf
//   """
// }

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
