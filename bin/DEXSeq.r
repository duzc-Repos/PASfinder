suppressPackageStartupMessages( library( "DEXSeq" ) )

inDir = commandArgs()[6]
tmp = list.files(inDir, pattern=".count$", full.names=TRUE)
controlFiles = list.files(inDir, pattern="^control", full.names=TRUE)
controlFiles = controlFiles[controlFiles %in% tmp]
treatmentFiles = list.files(inDir, pattern="^treatment", full.names=TRUE)
treatmentFiles = treatmentFiles[treatmentFiles %in% tmp]
gffFile = list.files(inDir, pattern=".gff$", full.names=TRUE)

test1 <- tmp
sampleTable = data.frame(
  row.names = unlist(strsplit(test1, split = '/'))[seq(4, length(test1)*4, 4)],
  condition = c(rep('control', length(controlFiles) ), rep('treatment', length(treatmentFiles)))
  )

dxd = DEXSeqDataSetFromHTSeq(
  test1,
  sampleData=sampleTable,
  design = ~ sample + exon + condition:exon,
  flattenedfile = gffFile
)

dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd )
dxd = testForDEU( dxd )
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")
dxr1 = as.data.frame( DEXSeqResults( dxd ) )

dxr2 = as.data.frame(dxr1)[c(3, 6:10)]
dxr2$pvalue[is.na(dxr2$pvalue)] <- 1
dxr2 = dxr2[order(dxr2$pvalue), ]
dxr2$padj[is.na(dxr2$padj)] <- 1
write.table(file = paste(inDir, 'DEXSeq_treatment_vs_control.txt', sep="/"), dxr2, sep = "\t", col.names = T, quote = F, row.names = T)


