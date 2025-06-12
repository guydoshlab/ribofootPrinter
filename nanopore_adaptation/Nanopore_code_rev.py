from Bio import SeqIO
from Bio import Seq
import csv

## This is the beginning of the mammalian_builddense.py code that is modified here for a new purpose.
## The goal of this code is to output the sequence in the transcriptome upstream of the 5' mapped end. It also provides an output of read lengths and read lengths normalized to the annotated length of the gene they map to.

def nanopore_motif_length(fasta_in,sam_in,outfile,smallsize,largesize):
	smallsize=int(smallsize)
	largesize=int(largesize)
	
	## text file with aligned read lengths for R analysis
	length_all=open(outfile+"_all.txt", "w") 
	
	## text file with aligned normalized read lengths for R analysis
	length_fraction=open(outfile+"_fraction.txt", "w") 
	
	## OUTPUT csv file with dimotifs at the -12-11 position from 5' end of the reads
	dimotiffile=open(outfile+"_dimotifs.csv", "w", newline="") 
	writer = csv.writer(dimotiffile, delimiter=',')
	
	## fasta file created for weblogo analysis, it takes 20 nt at the moment
	fastafile=open(outfile+".fa", "w") 
	
	outputdata={}
	negativestrandreads=0
	fsam=open(sam_in)
	samgen=csv.reader(fsam,delimiter='	')
	csv.field_size_limit(100000000) # This is needed to accomodate larger fields in the samfile. 
	mappedreads={}	# A dictionary that doesn't end up being used; for development.
	mappedreads["all"]=0	
	
	TT_motifs=0
	TA_motifs=0
	TG_motifs=0
	TC_motifs=0
	CC_motifs=0
	CA_motifs=0
	CG_motifs=0
	CT_motifs=0
	GG_motifs=0
	GA_motifs=0
	GT_motifs=0
	GC_motifs=0
	AT_motifs=0
	AA_motifs=0
	AG_motifs=0
	AC_motifs=0
	
	for readlen in list(range(smallsize,largesize+1)):
		mappedreads[readlen]=0

	# Add in sequences etc. 
	# Note that hashed out lines are eliminated from the original script.
	for record in SeqIO.parse(fasta_in, "fasta"):
		gene=record.id.split("|")[1]
		outputdata[gene]=["alias","sequence",{},"ORFstart","utr3start","transcriptlen"]
		#outputdata[gene][0]=record.id.split("|")[5]	# alias
		outputdata[gene][1]=str(record.seq)	# Sequence of transcript
		#outputdata[gene][3]=(record.id.split("|")[7]).split("-")[-1]	# ORF start
		#outputdata[gene][4]=(record.id.split("|")[8]).split("-")[-1]	# UTR3 start
		outputdata[gene][5]=record.id.split("|")[6]	# Transcript length
		genelength=int(outputdata[gene][5])
		#outputdata[gene][2]["all_5"]=[0 for x in range(genelength)] 
		#outputdata[gene][2]["all_3"]=[0 for x in range(genelength)] 
	
	print("Sequences added.")
	
	# Loop through the samfile. 
	for read in samgen:
		if read[0][0] == '@':   # Ignore header lines
			continue
		if read[1] == '4':      # A bowtie no match. 
			continue
        
		gene = read[2]             # gene identified for read in bowtie
		#readid = read[0]            # read id
		startp = int(read[3]) - 1    # start position. Need to subtract 1 since genomic sequence starts at 1, but list structures for read summing start at 0.
		seq = Seq.Seq(read[9])      # sequence of the read
		
		
		## First part of this code here writes out the lengths of all reads.
		length = len(seq)           # length of read
		genelength=int(outputdata[gene][5])
		length_all.write(str(length)+'\n')
		length_fraction.write(str(length/genelength)+'\n')
	
		## The remainder of this code now filters for size/strand and finds the motif.
		if (length<smallsize or length > largesize):
			continue

		# Remove negative strands
		if (read[1] == '16'):
			negativestrandreads+=1
			continue

		if (read[1] == '0'):
			if startp > genelength/3*2:	## Threshold for read fragments, <1/3 of total length.
				#this is a fasta file output for creating a weblogo, takes the first 20 
				window=outputdata[gene][1][startp-20:startp]
				window_rep=window.replace("T", "U")
				fastafile.write('>' +'\n'+ window_rep +'\n')
				#motif at 13 and 12 position 
				motif=outputdata[gene][1][startp-13:startp-11]
				if motif =="TT":
					TT_motifs+=1
				elif motif == "TA":
					TA_motifs+=1
				elif motif == "TG":
					TG_motifs+=1
				elif motif == "TC":
					TC_motifs+=1	
				elif motif == "CT":
					CT_motifs+=1
				elif motif == "CA":
					CA_motifs+=1
				elif motif == "CG":
					CG_motifs+=1
				elif motif == "CC":
					CC_motifs+=1
				elif motif == "GT":
					GT_motifs+=1
				elif motif == "GA":
					GA_motifs+=1
				elif motif == "GG":
					GG_motifs+=1
				elif motif == "GC":
					GC_motifs+=1
				elif motif == "AT":
					AT_motifs+=1
				elif motif == "AA":
					AA_motifs+=1
				elif motif == "AG":
					AG_motifs+=1
				elif motif == "AC":
					AC_motifs+=1
			
	writer.writerow(["UU"]+[str(TT_motifs)])
	writer.writerow(["UA"]+[str(TA_motifs)])
	writer.writerow(["UG"]+[str(TG_motifs)])
	writer.writerow(["UC"]+[str(TC_motifs)])
	
	writer.writerow(["CU"]+[str(CT_motifs)])
	writer.writerow(["CA"]+[str(CA_motifs)])
	writer.writerow(["CG"]+[str(CG_motifs)])
	writer.writerow(["CC"]+[str(CC_motifs)])
	
	writer.writerow(["GU"]+[str(GT_motifs)])
	writer.writerow(["GA"]+[str(GA_motifs)])
	writer.writerow(["GG"]+[str(GG_motifs)])
	writer.writerow(["GC"]+[str(GC_motifs)])
	
	writer.writerow(["AU"]+[str(AT_motifs)])
	writer.writerow(["AA"]+[str(AA_motifs)])
	writer.writerow(["AG"]+[str(AG_motifs)])
	writer.writerow(["AC"]+[str(AC_motifs)])
		
	length_all.close() 
	length_fraction.close()
	fastafile.close() 
	dimotiffile.close()
	fsam.close()
    