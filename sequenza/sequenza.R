args = commandArgs(trailingOnly=TRUE)
print(args)

library("sequenza")

data.file = args[1]
seqzdata = sequenza.extract(data.file, chromosome.list = paste0("chr",c(1:22,"X","Y")))

CP = sequenza.fit(seqzdata)

id = args[2]
out.dir = paste0("sequenza/", id)
sequenza.results(sequenza.extract = seqzdata, cp.table = CP, sample.id = id, out.dir=out.dir)
 
 
