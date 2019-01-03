#!/usr/bin/env nextflow

//  main.nf
//
//
//  Created by jianhong ou on 11/17/18.
//

import ParamsHelper


if (params.verbose)
echo true

today = new java.util.Date()

ParamsHelper.checkNonEmptyParam(params.dataDir, "dataDir");

log.info """\
$today
A T A C S E Q   P R O C E S S I N G   P I P E L I N E
=====================================================
dataDir : ${params.dataDir}
outdir  : ${params.outdir}
available CPU : ${params.cpus}
"""
.stripIndent()


// Open channel for left and right files and merge it into triples, the
// first entry is the LCS of the file names that can be used as a read
// pair identifier.
Channel
  .fromFilePairs("${params.dataDir}/*_{R,}{1,2}.f{ast,}q.gz")
  .into {readPairsFastQCOriginal; readPairsForTrim}

// --------------------------------------------------------------------------
// Step 1) Run FastQC
//
// - yields report
// --------------------------------------------------------------------------

process runFastQCOriginal {
cpus params.fastqc.cpus

input:
set runID, file(reads) from readPairsFastQCOriginal

storeDir "${params.outdir}/fastqc-original/${runID}"

output:
set file('*.zip'), file('*.html') into fastqcOutputOriginal

script:
"""
${params.fastqc.path} -t ${params.fastqc.cpus} -o . ${reads}
"""
}


// --------------------------------------------------------------------------
// Step 2) Run adapter trimming
//
// - yields trimmed read, used as downstream input
// --------------------------------------------------------------------------

process runTrimming {
cpus params.trim_galore.cpus

input:
set runID, file(reads) from readPairsForTrim

storeDir "${params.outdir}/trimmed/${runID}"

output:
set runID, file("${runID}_val_1.fq.gz"), file("${runID}_val_2.fq.gz") into readPairsTrimmed
set file('*.zip'), file('*.html') into trim_galoreOut

script:
readL=reads[0]
readR=reads[1]
"""
# call trim_galore
echo "PE: ${readL} ${readR}"
${params.trim_galore.path} -q 15 --paired --fastqc -o . \\
${readL} ${readR}
mv ${readL.simpleName}_val_1.fq.gz ${runID}_val_1.fq.gz
mv ${readR.simpleName}_val_2.fq.gz ${runID}_val_2.fq.gz
mv ${readL.simpleName}_val_1_fastqc.html ${runID}_val_1_fastqc.html
mv ${readL.simpleName}_val_1_fastqc.zip ${runID}_val_1_fastqc.zip
mv ${readR.simpleName}_val_2_fastqc.html ${runID}_val_2_fastqc.html
mv ${readR.simpleName}_val_2_fastqc.zip ${runID}_val_2_fastqc.zip
"""
}

readPairsTrimmed.into{readPairsTrimmed4BWA; readPairsTrimmed4blast}

// --------------------------------------------------------------------------
// Step 3) Align reads using BWA-MEM
//
// - align reads
// - sort
// - mark duplicates
// - yields alignment for downstream processing
// --------------------------------------------------------------------------

process runBWA {
  cpus params.bwa.cpus
  
  input:
    set runID, file(readL), file(readR) from readPairsTrimmed4BWA
  
  storeDir "${params.outdir}/bwa/${runID}"
  
  output:
    set runID, file("${runID}.bam"), file("${runID}.bam.bai") into mappedFiles
  
  script:
    """
    ${params.bwa.path} mem -M -t ${params.bwa.cpus} ${params.bwa.index} ${readL} ${readR} > ${runID}.sam
    ${params.samtools.path} sort -@ ${params.bwa.cpus} -n ${runID}.sam -o ${runID}.bam
    ${params.samtools.path} fixmate -@ ${params.bwa.cpus} -m ${runID}.bam ${runID}.fixout.bam
    ${params.samtools.path} sort ${runID}.fixout.bam -@ ${params.bwa.cpus} -o ${runID}.sort.bam
    ${params.samtools.path} markdup ${runID}.sort.bam ${runID}.rem.bam -@ ${params.bwa.cpus} -r
    mv ${runID}.rem.bam ${runID}.bam
    ${params.samtools.path} index ${runID}.bam
    """
}

mappedFiles.into{mappedFiles4MACS2; mappedFiles4ATACseqQC}

// --------------------------------------------------------------------------
// Step 4) Blastn reads to see if there is any cross-speices comtermination.
// --------------------------------------------------------------------------
if(false){
process runBlast {
  cpus params.blastn.cpus
  
  input:
    set runID, file(readL), file(readR) from readPairsTrimmed4blast
  
  storeDir "${params.outdir}/taxReport/${runID}"
  
  output:
    set runID, file("${runID}.blastn.txt"), file("${runID}.csv") into blastnFiles
  
  script:
    """
    zcat ${readL} | head -n ${params.blastn.numOfLine} > ${runID}.fq
    fq2sc -i ${runID}.fq -o ${runID}.txt
    s2c_2_fasta.pl --in ${runID}.txt --out ${runID}.fa --notfilterN
    blastn -num_threads ${params.blastn.cpus} -db nt -max_target_seqs 1 \\
    -outfmt '6 qaccver saccver staxid sblastname ssciname scomname pident length mismatch gapopen qstart qend sstart send evalue bitscore' \\
    -query ${runID}.fa -out ${runID}.blastn.txt
    cat <<EOT | R --vanilla
x <- read.delim("${runID}.blastn.txt", header=FALSE, stringsAsFactor=FALSE)
colnames(x) <- c("qaccver", "saccver", "staxid", "sblastname", "ssciname", "scomname", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
x <- x[x\$mismatch<=${params.blastn.mismatch}, ]
x <- x[order(x\$qaccver, -x\$pident, x\$evalue), ]
x <- x[!duplicated(x\$qaccver), ]
sc <- as.numeric(sub("SEQ_.*?_x", "", x\$qaccver))
species <- paste0(x\$ssciname, " (", x\$scomname, ")", " [", x\$sblastname, "]")
w <- rowsum(sc, species)
w <- w[order(-w[, 1]), ]
p <- w/${params.blastn.numOfLine}
write.table(p, '${runID}.percentage.txt', quote=FALSE, col.names=FALSE, sep='\t')
f <- dir(".", "percentage.txt", full.names = TRUE)
d <- lapply(f, read.delim, nrows=3, header=FALSE)
e <- sapply(d, function(.ele){
  if(nrow(.ele)>1){
    .ele[1, 2]/sum(.ele[-1, 2])
  }else{
    100
  }
})
e <- cut(e, breaks=c(0, 3, 5, 10, Inf), labels=c("contaminated", "warning", "minor", "clean"))
e <- rep(e, each=3)
fq <- sub(".percentage.txt", "", basename(f))
fq <- rep(fq, each=3)
d <- do.call(rbind, d)
d <- cbind(d, fastq=fq, cross.species.contamination=e)
colnames(d) <- c("species", "percentage", "fastq", "cross.species.contamination")
write.csv(d, "${runID}.csv", row.names=FALSE)
EOT
    """
}
}

// --------------------------------------------------------------------------
// Step 5) call peaks by macs2
// --------------------------------------------------------------------------

process runBWA {
  cpus params.macs2.cpus
  
  input:
    set runID, file(bam), file(bai) from mappedFiles4MACS2
  
  storeDir "${params.outdir}/macs2/${runID}"
  
  output:
    set runID, file("*") into macs2
    set runID, file("*_treat_pileup.bdg"), file{"${runID}.ChromInfo.txt"} into bedgraph
  
  script:
    """
    ${params.macs2.path} callpeak -t ${bam} \\
                 -f BAM -g ${params.macs2.genome}  -n ${runID} \\
                 --outdir . -q 0.05 \\
                 --nomodel --shift 37 --extsize 73 -B
    ${params.samtools.path} view -H ${bam} | grep @SQ | \\
        awk -F \$'\\t' 'BEGIN {OFS=FS} {gsub("[SL]N:", ""); print \$2,\$3}' \\
        > ${runID}.ChromInfo.txt
    """
}

// --------------------------------------------------------------------------
// Step 6) generate bigwig files
// --------------------------------------------------------------------------

process runBigWig {
  cpus params.bedGraphToBigWig.cpus
  
  input:
    set runID, file(bdg), file(size) from bedgraph
  
  storeDir "${params.outdir}/bw/${runID}"
  
  output:
    set runID, file("${runID}.bw") into bigwig
  
  script:
    """
    sort -k1,1 -k2,2n ${bdg} > ${bdg}.srt
    ${params.bedGraphToBigWig.path} ${bdg}.srt \\
                 ${size} \\
                 ${runID}.bw
    """
}

// --------------------------------------------------------------------------
// Step 7) ATACseqQC
// --------------------------------------------------------------------------

process runATACseqQC{
  cpus params.R.cpus
  
  input:
    set runID, file(bam), file(bai) from mappedFiles4ATACseqQC
  
  storeDir "${params.outdir}/ATACseqQC/${runID}"
  
  output:
    set runID, file("*.pdf") into ATACseqQCout
    set runID, file("*.bam"), file("*.bai") into ATACseqQCbam
  
  script:
    """
    cat <<EOF | ${params.R.path} --vanilla
    library(BiocManager)
    install("${params.R.BSgenome}")
    install("${params.R.TxDb}")
    library(ATACseqQC)
    bamfile <- "${bam}"
    bamfile.labels <- gsub(".bam", "", basename(bamfile))
    pdf("fragmentSizeDistribution.pdf")
    fragSize <- fragSizeDist(bamfile, bamfile.labels)
    dev.off()
    library(${params.R.TxDb})
    txs <- transcripts(${params.R.TxDb})
    library(${params.R.BSgenome})
    seqlev <- paste0("chr", c(1:21, "X", "Y"))
    which <- as(seqinfo(${params.R.BSgenome})[seqlev], "GRanges")
    gal <- readBamFile(bamfile, which=which, asMates=TRUE, bigFile=TRUE)
    gal1 <- shiftGAlignmentsList(gal, outbam="shifted.bam")
    pt <- PTscore(gal1, txs)
    pdf("PTscore.pdf")
    plot(mcols(pt)[, "log2meanCoverage"], mcols(pt)[, "PT_score"], 
     xlab="log2 mean coverage",
     ylab="Promoter vs Transcript")
    dev.off()
    nfr <- NFRscore(gal1, txs)
    pdf("NFRscore.pdf")
    plot(mcols(nfr)[, "log2meanCoverage"], mcols(nfr)[, "NFR_score"], 
     xlab="log2 mean coverage",
     ylab="Nucleosome Free Regions score",
     main="NFRscore for 200bp flanking TSSs",
     xlim=c(-10, 0), ylim=c(-5, 5))
    dev.off()
    tsse <- TSSEscore(gal1, txs)
    pdf("TSSEscore.pdf")
    hist(mcols(tsse)[, "TSS.enrichment.score"], breaks=100, 
    main="Transcription Start Site (TSS) Enrichment Score", 
    xlab="TSS enrichment score")
    dev.off()
    gc(reset=TRUE)
    genome <- ${params.R.BSgenome}
    objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome, outPath = ".")
    rm(gal1)
    gc(reset=TRUE)
    library(ChIPpeakAnno)
    bamfiles <- file.path(".",
                     c("NucleosomeFree.bam",
                     "mononucleosome.bam",
                     "dinucleosome.bam",
                     "trinucleosome.bam"))
    pdf("cumulativePercentage.pdf")
    cumulativePercentage(bamfiles[1:2], as(seqinfo(${params.R.BSgenome})[seqlev], "GRanges"))
    dev.off()
    TSS <- promoters(txs, upstream=0, downstream=1)
    TSS <- unique(TSS)
    librarySize <- estLibSize(bamfiles)
    NTILE <- 101
    dws <- ups <- 1010
    sigs <- enrichedFragments(gal=objs[c("NucleosomeFree", 
                                     "mononucleosome",
                                     "dinucleosome",
                                     "trinucleosome")], 
                          TSS=TSS,
                          librarySize=librarySize,
                          seqlev=seqlev,
                          TSS.filter=0.5,
                          n.tile = NTILE,
                          upstream = ups,
                          downstream = dws)
      sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
      pdf("featureAligndHeatmap.pdf")
      featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                      zeroAt=.5, n.tile=NTILE)
      dev.off()
      out <- featureAlignedDistribution(sigs, 
                                  reCenterPeaks(TSS, width=ups+dws),
                                  zeroAt=.5, n.tile=NTILE, type="l", 
                                  ylab="Averaged coverage")
    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    out <- apply(out, 2, range01)
    pdf("TSS.profile.pdf")
    matplot(out, type="l", xaxt="n", 
        xlab="Position (bp)", 
        ylab="Fraction of signal")
    axis(1, at=seq(0, 100, by=10)+1, 
     labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
    abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
    dev.off()
    library(MotifDb)
    CTCF <- query(MotifDb, c("CTCF"))
    CTCF <- as.list(CTCF)
    pdf("CTCF.footprint.pdf")
    sigs <- factorFootprints("shifted.bam", pfm=CTCF[[1]], 
                         genome=genome,
                         min.score="90%", seqlev=seqlev,
                         upstream=100, downstream=100)
    dev.off()
    pdf("CTCF.Vplot.pdf")
    vp <- vPlot("shifted.bam", pfm=CTCF[[1]], 
            genome=genome, min.score="90%", seqlev=seqlev,
            upstream=200, downstream=200, 
            ylim=c(30, 250), bandwidth=c(2, 1))
    dev.off()
    unlink("Rplots.pdf")
EOF
    """
}