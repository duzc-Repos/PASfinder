library('cleanUpdTSeq')

Args <- commandArgs()
input_file <- Args[6]
output_file <- Args[7]
#cat(Args)

testSet <- read.table(input_file, sep = "\t", header = T)
peaks <- BED2GRangesSeq(testSet, upstream.seq.ind=7, 
                        downstream.seq.ind=8, withSeq=TRUE)
testSet.NaiveBayes <- buildFeatureVector(peaks, 
                                         upstream=40, downstream=30, 
                                         wordSize=6, alphabet=c("ACGT"),
                                         sampleType="unknown", 
                                         replaceNAdistance=30, 
                                         method="NaiveBayes",
                                         ZeroBasedIndex=1, fetchSeq=FALSE)

data(classifier)
test_results <- predictTestSet(testSet.NaiveBayes=testSet.NaiveBayes,
                              classifier=classifier,
                              outputFile=NULL,
                              assignmentCutoff=0.5)

test_results <- test_results[test_results[2] < 0.5,  1]
write.table(test_results, output_file, quote = F, row.names = F, col.names = F)


