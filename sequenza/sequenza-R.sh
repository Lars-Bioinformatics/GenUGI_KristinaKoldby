library("sequenza")

data.file = "ECV2-4-biopsi-H1.bins.seqz.gz"

test2 = sequenza.extract(data.file)

CP = sequenza.fit(test2)

sequenza.results(sequenza.extract = test2,cp.table = CP, sample.id = "ECV2-4-biopsi-H1", out.dir="ECV2-4-biopsi-H1")