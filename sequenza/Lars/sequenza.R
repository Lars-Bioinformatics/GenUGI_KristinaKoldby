args = commandArgs(trailingOnly=TRUE)
print(args)

library("sequenza")

#data.file = "G37-01226-06_G37-A085.seqz.bin500.noGL.noMT.gz"
data.file = args[1]
seqzdata = sequenza.extract(data.file, chromosome.list = c(1:22,"X","Y"))

CP = sequenza.fit(seqzdata)

id = paste0("Sequenza_", args[2])
out.dir = paste0(args[3], id, "_results")
sequenza.results(sequenza.extract = seqzdata, cp.table = CP, sample.id = id, out.dir=out.dir)
 
 