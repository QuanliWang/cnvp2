library(CNVPanelizer)
#prepare_reference.R
load("last.RData")

# Directory with the test data
sampleDirectory <- "/Users/quanliwang/Desktop/CNVtest/cnv"
# Vector with test filenames
sampleFilenames <- list.files(path = sampleDirectory, pattern = ".bam$",
                              full.names = TRUE)
# Read the sample of interest data set
sampleReadCounts <- ReadCountsFromBam(sampleFilenames, genomicRangesFromBed,
                                      sampleNames = sampleFilenames,
                                      ampliconNames = ampliconNames,
                                      removeDup = removePcrDuplicates)


normalizedReadCounts <- CombinedNormalizedCounts(sampleReadCounts, referenceReadCounts,
                                                 ampliconNames = ampliconNames)
# After normalization data sets need to be splitted again to perform bootstrap
samplesNormalizedReadCounts = normalizedReadCounts["samples"][[1]]
referenceNormalizedReadCounts = normalizedReadCounts["reference"][[1]]

# Number of bootstrap replicates to be used
replicates <- 10000

# Perform the bootstrap based analysis
bootList <- BootList(geneNames, samplesNormalizedReadCounts,
                     referenceNormalizedReadCounts,
                     replicates = replicates)

# Estimate the background noise left after normalization
backgroundNoise <- Background(geneNames, samplesNormalizedReadCounts,
                              referenceNormalizedReadCounts,
                              bootList,
                              replicates = replicates,
                              significanceLevel = 0.05,
                              sigene = 0.05,
                              robust = TRUE)

# Build report tables
reportTables <- ReportTables(geneNames, samplesNormalizedReadCounts,
                             referenceNormalizedReadCounts,
                             bootList,
                             backgroundNoise,
                             uppper_limit = 1.5,
                             lower_limit = 0.5)

PlotBootstrapDistributions(bootList, reportTables, save = TRUE)


#pdf(file="cnv.pdf",width=6.5,height=5)
#PlotBootstrapDistributions(bootList, reportTables)
#dev.off()



