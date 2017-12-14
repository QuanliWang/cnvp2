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

# Should duplicated reads (same start, end site and strand) be removed
removePcrDuplicates <- FALSE # TRUE is only recommended for Ion Torrent data

# Directory with the reference data
referenceDirectory <- "/Users/quanliwang/Desktop/CNVtest/normal"
referenceFilenames <- list.files(path = referenceDirectory, pattern = ".bam$",
                                 full.names = TRUE)
# Read the Reference data set
referenceReadCounts <- ReadCountsFromBam(referenceFilenames, genomicRangesFromBed,
                                         sampleNames = referenceFilenames,
                                         ampliconNames = ampliconNames,
                                         removeDup = removePcrDuplicates)