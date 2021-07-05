from Bio import SeqIO
from Bio import Seq
import csv
import pickle
import gzip


	
# Put -1 for normalize to use all reads to compute rpm. Alternatively, put in total reads to divide by or 1E6 for no normalization.
# Smallsize and largesize are the limits (inclusive) of which read sizes to include. -
# Individual tracks for each sized reads are commented out because the files get unbearable (~20 GB) and tend to crash. Can be activated for specific needs.
# outfile is where the pickled density file will go.
# fasta in has to be "full name" fasta, with "|" separated record.id fields for each gene as described in code below.
def builddense_m(fasta_in,sam_in,outfile,normalize,smallsize,largesize):
	normalize=int(normalize)
	smallsize=int(smallsize)
	largesize=int(largesize)
	
	outputdata={}
	negativestrandreads=0
	fsam=open(sam_in)
	samgen=csv.reader(fsam,delimiter='	')
	mappedreads={}
	mappedreads["all"]=0	
	for readlen in list(range(smallsize,largesize+1)):
		mappedreads[readlen]=0
	
	
	# Add in sequences etc. 
	for record in SeqIO.parse(fasta_in, "fasta"):
		gene=record.id.split("|")[1]
		outputdata[gene]=["alias","sequence",{},"ORFstart","utr3start","transcriptlen"]
		outputdata[gene][0]=record.id.split("|")[5]	# alias
		outputdata[gene][1]=str(record.seq)	# Sequence of transcript
		outputdata[gene][3]=(record.id.split("|")[7]).split("-")[-1]	# ORF start
		outputdata[gene][4]=(record.id.split("|")[8]).split("-")[-1]	# UTR3 start
		outputdata[gene][5]=record.id.split("|")[6]	# Transcript length
		genelength=int(outputdata[gene][5])
		#for readlen in list(range(smallsize,largesize+1)):	# Put back for readsize specific work.
		#	outputdata[gene][2][str(readlen)+"_5"]=[0 for x in range(genelength)] # Put back for readsize specific work.
		#	outputdata[gene][2][str(readlen)+"_3"]=[0 for x in range(genelength)] # Put back for readsize specific work.
		outputdata[gene][2]["all_5"]=[0 for x in range(genelength)] 
		outputdata[gene][2]["all_3"]=[0 for x in range(genelength)] 
	
	print("Sequences added.")
	
	# Loop through the samfile.
	for read in samgen:
		if read[0][0] == '@':   # Ignore header lines.
			continue
		if read[1] == '4':      # A bowtie no match. 
			continue
        
		gene = read[2]             # gene identified for read in bowtie
		readid = read[0]            # read id
		startp = int(read[3]) - 1    # start position. Need to subtract 1 since genomic sequence starts at 1, but list structures for read summing start at 0.
		seq = Seq.Seq(read[9])      # sequence of the read
		length = len(seq)           # length of read
		genelength=int(outputdata[gene][5])
		
		if (length<smallsize or length > largesize):
			continue

		# Remove negative strands
		if (read[1] == '16'):
			negativestrandreads+=1
			continue

		if (read[1] == '0'):			
			#outputdata[gene][2][str(length)+"_5"][startp]+=1		# Put back for readsize specific work.
			#outputdata[gene][2][str(length)+"_3"][startp+length-1]+=1		# Put back for readsize specific work.
			#mappedreads[length]+=1		# Put back for readsize specific work.
			outputdata[gene][2]["all_5"][startp]+=1
			outputdata[gene][2]["all_3"][startp+length-1]+=1
			mappedreads["all"]+=1
	
		
	print(str(mappedreads["all"])+" reads mapped to transcriptome.")
	print(str(negativestrandreads)+" negative reads mapped to transcriptome.")


	# Normalize
	genecount=0
	if normalize==-1:
		for gene in outputdata.keys():
			genecount+=1
			if genecount%1000==0:
				print("Genes normalized="+str(genecount))
			#for readlen in list(range(smallsize,largesize+1)):
			#	outputdata[gene][2][str(readlen)+"_5"]=(normalizelist(outputdata[gene][2][str(readlen)+"_5"],(mappedreads[readlen]/float(1E6))))
			#	outputdata[gene][2][str(readlen)+"_3"]=(normalizelist(outputdata[gene][2][str(readlen)+"_3"],(mappedreads[readlen]/float(1E6))))
			outputdata[gene][2]["all_5"]=(normalizelist(outputdata[gene][2]["all_5"],(mappedreads["all"]/float(1E6))))
			outputdata[gene][2]["all_3"]=(normalizelist(outputdata[gene][2]["all_3"],(mappedreads["all"]/float(1E6))))
	else:
		for gene in outputdata.keys():
			genecount+=1
			if genecount%1000==0:
				print("Genes normalized="+str(genecount))
			#for readlen in list(range(smallsize,largesize+1)):
			#	outputdata[gene][2][str(readlen)+"_5"]=(normalizelist(outputdata[gene][2][str(readlen)+"_5"],(normalize/float(1E6))))
			#	outputdata[gene][2][str(readlen)+"_3"]=(normalizelist(outputdata[gene][2][str(readlen)+"_3"],(normalize/float(1E6))))
			outputdata[gene][2]["all_5"]=(normalizelist(outputdata[gene][2]["all_5"],(normalize/float(1E6))))
			outputdata[gene][2]["all_3"]=(normalizelist(outputdata[gene][2]["all_3"],(normalize/float(1E6))))

	
	# Pickle the dictionary
	f = gzip.open(outfile+'.pkl.gzip', 'wb')
	pickle.dump(outputdata,f)
	f.close()
	fsam.close()

	
def normalizelist(listin,normfactor):
	newlist=[]
	for element in listin:
		newlist.append(element/float(normfactor))	
	return newlist


	