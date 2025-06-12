import csv
from Bio import SeqIO
import tools_m
import re


# Wrapper for countreadsizes
def countreadsizes_wrapper(txtinfile):
	f=open(txtinfile)
	inputparams=[line.strip() for line in f.readlines()] 
	f.close()
	
	variables={}
	nextlineisinput=0
	for inputline in inputparams:
		if inputline=="" or inputline[0]=="#":
			continue
		
		if nextlineisinput!=0:
			variables[nextlineisinput]=inputline
			nextlineisinput=0
		
		if inputline[0:5]=="INPUT":
		# The next line is an input line.
			nextlineisinput=inputline.split("_")[-1]
			
	print("Variables collected:")
	for key in (variables.keys()):
		print(key)
		print(variables[key])
		print("")
	
	outfile=variables["outfile"]
	lowest=int(variables["lowest"])
	highest=int(variables["highest"])
	inputfiles=list(variables["fastqnames"].split(","))
	writerfile=open(outfile+".csv", "w")	
	writer = csv.writer(writerfile)
	
	# Writer headers for csv file.
	writer.writerow(["size"]+list(range(highest+1))[lowest:])
	
	for fastqfile in inputfiles:	# Loop through the file list and count read sizes.
	
		samplename=re.split("[./]",fastqfile)
		# Take filename without path and extension, otherwise take whatever is there.
		if (len(samplename)>1):
			samplename=samplename[-2]	
		else:
			samplename=samplename[-1]
		print("Sample name = "+samplename)
		print("Motif = "+variables["motif"])
		
		# Get the histogram.
		sizes=countreadsizes(fastqfile,lowest,highest,variables["motif"])	

		writer.writerow([samplename+"_"+variables["motif"]]+sizes)		# write out the histogram.
	
	writerfile.close()				# Close the writer file.
	tools_m.transposecsv(outfile)	# Transpose csv for easy viewing.


# Simple program that reports read sizes.
# Filters for the terminal motif on 3' end.
# Putting in no motif will run on all reads.
def countreadsizes(f_in,lowest,highest,end3motif):
	# Setting counters to 0.
	skipped_size=0	
	skipped_motif=0
	withmotif=0
	totcount=0
	# Get length of motif.
	motiflen=len(end3motif)
	if (motiflen)>lowest:
		print("ERROR - filter motif longer than lower limit")
		exit()
	sizes=[0 for x in range(highest-lowest+1)]	# Establish counts table.
	
	recordsin=SeqIO.parse(f_in,"fastq")
	for record in recordsin:
		totcount+=1
		reclen=len(record)
		if reclen>highest or reclen<lowest:
			skipped_size+=1
			continue		
		if str(record[-motiflen:].seq)==end3motif or end3motif=="none":
			reclen=len(record)
			sizes[reclen-lowest]+=1
			withmotif+=1
		else:
			skipped_motif+=1
	
	
	print("Reads lacking motif = "+str(skipped_motif))
	print("Reads with motif = "+str(withmotif))
	print("")
	
	return(sizes)
	
	


# Wrapper for motif distribution
def endmotif_distribution_wrapper(txtinfile):
	f=open(txtinfile)
	inputparams=[line.strip() for line in f.readlines()] 
	f.close()
	
	variables={}
	nextlineisinput=0
	for inputline in inputparams:
		if inputline=="" or inputline[0]=="#":
			continue
		
		if nextlineisinput!=0:
			variables[nextlineisinput]=inputline
			nextlineisinput=0
		
		if inputline[0:5]=="INPUT":
		# The next line is an input line.
			nextlineisinput=inputline.split("_")[-1]
			
	print("Variables collected:")
	for key in (variables.keys()):
		print(key)
		print(variables[key])
		print("")
	
	outfile=variables["outfile"]
	motiflen=int(variables["motiflength"])
	whichend=int(variables["whichend"])
	inputfiles=list(variables["fastqnames"].split(","))
	writerfile=open(outfile+".csv", "w")	
	writer = csv.writer(writerfile)

	# Make an initialized motif dictionary.
	motif_dictionary_init={}
	i=1
	keylist=["A","C","G","T"]
	while i<motiflen:
		tempkeylist=list(keylist)
		keylist=[]
		for base2 in ["A","C","G","T"]:
			for base in tempkeylist:
				keylist.append(base+base2)		
		i+=1
			
	for key in keylist:
		motif_dictionary_init[key]=0
	
	temp_motifs_list=[]		
	for key in motif_dictionary_init.keys():
		temp_motifs_list.append(key)
	writer.writerow(["Motifs"]+temp_motifs_list)
			
			
	for fastqfile in inputfiles:	# Loop through the file list and count read sizes.
		motif_dictionary=dict(motif_dictionary_init) # Copy the empty dictionary for use.
		samplename=re.split("[./]",fastqfile)
		# Take filename without path and extension, otherwise take whatever is there.
		if (len(samplename)>1):
			samplename=samplename[-2]	
		else:
			samplename=samplename[-1]
		print("Sample name = "+samplename)
		print("Motif_length = "+str(motiflen))
		
		# Get the histogram.
		endmotif_distribution(fastqfile,whichend,motif_dictionary,motiflen)	

		temp_counts_list=[]
		for key in motif_dictionary.keys():
			temp_counts_list.append(motif_dictionary[key])
		writer.writerow([samplename+"_counts"]+temp_counts_list)
	
	writerfile.close()				# Close the writer file.
	tools_m.transposecsv(outfile)	# Transpose csv for easy viewing.


# Code to count particular nt motifs at the 5' or 3' end of reads in a fastq file.
def endmotif_distribution(f_in,whichend,motif_dictionary,motiflen):

	skipped_size=0	
	totcount=0
	skipped_N=0
	
	recordsin=SeqIO.parse(f_in,"fastq")
	for record in recordsin:
		totcount+=1
		reclen=len(record)
		if reclen<motiflen:
			skipped_size+=1
			continue
		if whichend==3:
			motif=str(record[-motiflen:].seq)
		elif whichend==5:
			motif=str(record[0:motiflen].seq)
		
		if motif in motif_dictionary.keys():
			motif_dictionary[motif]+=1
		else:
			skipped_N+=1
	
	print("Total reads = "+str(totcount))
	print("Reads outside range = "+str(skipped_size))	
	print("Reads skipped for N or other bases = "+str(skipped_N))
	print("")
	# There is no need to return the dictionary since dictionaries are mutable objects in Python.		
			
