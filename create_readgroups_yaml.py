import glob, os

# input dirs containing fastq files
inputDirs = ["../../fastq_files_for_sanger"]

# Raw filetype
filetype = ".fastq.gz"

# Output yaml filename
def getYamlFileName(sm): 
	return("ReadGroups_"+sm+".yaml")

def writeYamlHeader(out, sm):
	out.write(
	"# abstracted as has to be same for all records\n" +
	"SM: "+ sm + "\n" +
	"# the actual readgroups\n" +
	"READGRPS:"
	)

def writeYamlSample(out, basename, sampleid, protocol, flowcell):
    # Very strict about using spaces and not tabs!
	out.write(
		"\n    " + basename + filetype + ":\n" +
		"        ID: " + sampleid + "_" + flowcell + "\n" +
		"        PL: ILLUMINA \n" +
		"        LB: " + protocol + "\n" +
		"        PU: " + flowcell + "_" + sampleid + "\n" +
		"        DS: WGS_ILLUMINA"
	)

for directory in inputDirs:
	#print directory
	for filename in glob.iglob(directory+"/*"+filetype):
		
		# Assume filename follows the following pattern:
		# SAMPLEID_PROTOCOL_FLOWCELLID_READ.filetype
		basename = filename.split("/")[-1].split(filetype)[0]
		sampleid, protocol, flowcell, lane, read = basename.split("_")
		
		#print sampleid, protocol, flowcell, read
		
		sm=sampleid+"_"+protocol+"_"+flowcell+"_"+lane
		
		# Get output yaml file name
		yaml_file = getYamlFileName(sm)
		
		# Write yaml file
		with open(yaml_file, "a") as out:
			# if file empty, create header
			if os.path.getsize(yaml_file) == 0:
				writeYamlHeader(out, sm)
			
			# Write read group information for sample
			writeYamlSample(out, basename, sampleid, protocol, flowcell)
			
