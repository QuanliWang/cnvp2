library(CNVPanelizer)
bedFilepath <- "/Users/quanliwang/Desktop/CNVtest/DHS-005Z.roi.bed"

# The column number of the gene Names in the BED file.
amplColumnNumber <- 4

# Extract the information from a bed file
genomicRangesFromBed <- BedToGenomicRanges(bedFilepath,
                                           ampliconColumn = amplColumnNumber,
                                           split = "_")

metadataFromGenomicRanges <- elementMetadata(genomicRangesFromBed) 
geneNames = metadataFromGenomicRanges["geneNames"][, 1] 
ampliconNames = metadataFromGenomicRanges["ampliconNames"][, 1]

# Directory with the test data
sampleDirectory <- "/Users/quanliwang/Desktop/CNVtest/cnv"

# Directory with the reference data
referenceDirectory <- "/Users/quanliwang/Desktop/CNVtest/normal"

# Vector with test filenames
sampleFilenames <- list.files(path = sampleDirectory, pattern = ".bam$",
                              full.names = TRUE)

referenceFilenames <- list.files(path = referenceDirectory, pattern = ".bam$",
                                 full.names = TRUE)

# Should duplicated reads (same start, end site and strand) be removed
removePcrDuplicates <- FALSE # TRUE is only recommended for Ion Torrent data

# Read the Reference data set
referenceReadCounts <- ReadCountsFromBam(referenceFilenames, genomicRangesFromBed,
                                         sampleNames = referenceFilenames,
                                         ampliconNames = ampliconNames,
                                         removeDup = removePcrDuplicates)

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



pdf(file="cnv.pdf",width=6.5,height=5)
PlotBootstrapDistributions(bootList, reportTables)
dev.off()



