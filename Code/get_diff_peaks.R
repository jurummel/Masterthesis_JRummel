library(DiffBind)

# Load sample sheets
samplesheet_BL6 <- read.csv("/projects/allelespecificchromatinmouse/work/MouseSeqData/results/DNA-Mapping/diff_peaks/BL6/samplesheet_BL6.csv", header=T, stringsAsFactors=F, sep=";")
samplesheet-CC <- read.csv("/projects/allelespecificchromatinmouse/work/MouseSeqData/results/DNA-Mapping/diff_peaks/CC/samplesheet_CC.csv", header=T, stringsAsFactors=F, sep=";")

# Create a DBA object
dba_obj_BL6 <- dba(sampleSheet = samplesheet_BL6, peakFormat="narrow")
dba_obj_CC <- dba(sampleSheet = samplesheet_CC, peakFormat="narrow")

# Count the number of reads in each peak for each sample
dba_count_BL6 <- dba.count(dba_obj_BL6, summits=75)
dba_count_CC <- dba.count(dba_obj_CC, summits=75)

# create contrast between both conditions
dba_contrast_BL6 <- dba.contrast(dba_count_BL6, categories=DBA_CONDITION, minMembers=2)
dba_contrast_CC <- dba.contrast(dba_count_CC, categories=DBA_CONDITION, minMembers=2)

# Perform differential peak calling
dba_res_BL6 <- dba.analyze(dba_contrast_BL6, method=DBA_DESEQ2)
dba_res_CC <- dba.analyze(dba_contrast_CC, method=DBA_DESEQ2)

# extract regions FDR < 0.05
dba_db_BL6 <- dba.report(dba_res_BL6)
dba_db_CC <- dba.report(dba_res_CC)
# extract all regions
dba_db_BL6 <- dba.report(dba_res_BL6, th=1, method=DBA_DESEQ2, bCounts=TRUE)
dba_db_CC <- dba.report(dba_res_CC, th=1, method=DBA_DESEQ2, bCounts=TRUE)

# Write results to file
write.table(dba_db, file="/projects/allelespecificchromatinmouse/work/MouseSeqData/results/DNA-Mapping/diff_peaks/BL6/differential_peaks_BL6_new.txt", sep="\t", row.names=T, col.names=T)
write.table(dba_db, file="/projects/allelespecificchromatinmouse/work/MouseSeqData/results/DNA-Mapping/diff_peaks/CC/differential_peaks_CC_new.txt", sep="\t", row.names=T, col.names=T)

#MA plot not normalized
pdf(file="/projects/allelespecificchromatinmouse/work/MouseSeqData/results/DNA-Mapping/diff_peaks/BL6/MA_BL6_not_Norm.pdf")
dba.plotMA(dba_res_BL6, bNormalized=FALSE)
dev.off()

pdf(file="/projects/allelespecificchromatinmouse/work/MouseSeqData/results/DNA-Mapping/diff_peaks/CC/MA_CC_not_Norm.pdf")
dba.plotMA(dba_res_CC, bNormalized=FALSE)
dev.off()

# MA plot normalized 
pdf(file="/projects/allelespecificchromatinmouse/work/MouseSeqData/results/DNA-Mapping/diff_peaks/BL6/MA_BL6_new.pdf")
dba.plotMA(dba_res_BL6)
dev.off()

pdf(file="/projects/allelespecificchromatinmouse/work/MouseSeqData/results/DNA-Mapping/diff_peaks/CC/MA_CC_new.pdf")
dba.plotMA(dba_res_CC)
dev.off()



