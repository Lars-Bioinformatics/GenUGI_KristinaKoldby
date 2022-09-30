# args = commandArgs(trailingOnly=TRUE)
# print(args)

library("sequenza")

### ECV2-4-biopsi-H1

data.file = "ECV2-4-biopsi-H1_tumor_tagseq-medexome-deep-seq.bins.seqz.gz"
seqzdata = sequenza.extract(data.file, chromosome.list = paste0("chr",c(1:22,"X","Y")))

# CP = sequenza.fit(seqzdata)

id = "ECV2-4-biopsi-H1_tumor_tagseq-medexome-deep-seq"
out.dir = paste0(id, "_ploidy3.7")
sequenza.results(sequenza.extract = seqzdata, sample.id = id, out.dir=out.dir, cellularity = 0.77, ploidy = 3.7)

### ECV2-4-biopsi-H2

# data.file = "ECV2-4-biopsi-H2_tumor_tagseq-medexome-deep-seq.bins.seqz.gz"
# seqzdata = sequenza.extract(data.file, chromosome.list = paste0("chr",c(1:22,"X","Y")))

# CP = sequenza.fit(seqzdata)

# id = "ECV2-4-biopsi-H2_tumor_tagseq-medexome-deep-seq"
# out.dir = paste0(id, "_ploidy3.6")
# sequenza.results(sequenza.extract = seqzdata, sample.id = id, out.dir=out.dir, cellularity = 0.74, ploidy = 3.6)


### Got this warning when runnings script for biopsi-H2:
# Warning message:
# In b[which(diff(b) == 0) + 1] <- b[diff(b) == 0] + offset :
#  number of items to replace is not a multiple of replacement length
