from Bio import SeqIO
from Bio import Seq
import csv
import pickle
import os
import struct
import re
import gzip



# Writegene2 wrapper function
def writegene2_m_wrapper(txtinfile):
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
	
	
	inputfiles=list(variables["picklenames"].split(","))
	genelist=list(variables["genenames"].split(","))
	
	# Handle Alias conversion:
	aliasdict={}
	f=gzip.open(inputfiles[0],"rb")
	footprints=pickle.load(f)
	f.close()
	for gene in footprints.keys():
		aliasdict[footprints[gene][0]]=gene
	for gnnum in range(len(genelist)):
		if genelist[gnnum][0:3]!="NM_":
			alias=genelist[gnnum]
			if alias in aliasdict:
				genelist[gnnum]=aliasdict[alias]
			else:
				print("Error - alias "+alias+" does not exist. Try the NM_... format name.")
				exit()

	for picklefile in inputfiles:
		writegene2_m(genelist,variables["endmode"],picklefile,variables["outfile"]+"_"+re.split("[./]",picklefile)[-3])


#### writegene2_m function.
### endmode is "all_3" or "all_5"
## This function is used to put the reads that map to a particular gene into a csv file. This can be used for making figures or other analysis.
# genenames is a list of gene names
# endmode is 5' or 3' end alignments
# picklefile is the ribosome density file
# outfile is the name of the output csv file
def writegene2_m(genenames,endmode,picklefile,outfile):
	f=gzip.open(picklefile,"rb")
	footprints=pickle.load(f)
	print("Data loaded.")
	f.close()
	
	writerfile=open(outfile+".csv", "w")	
	writer = csv.writer(writerfile)
		
	for gene in genenames:
		footprintcounts=list(footprints[gene][2][endmode])
		alias=footprints[gene][0]
		writer.writerow([gene]+[alias]+footprintcounts)
	
	writerfile.close()	
	transposecsv(outfile)	
	
def transposecsv(csvfile):
	writerfile=open(csvfile+"_transposed.csv", "w")
	writer = csv.writer(writerfile,delimiter=',')
	f=open(csvfile+".csv")
	readercsv=csv.reader(f)
	maxrow=0
	for row in readercsv:
		if len(row)>maxrow:
			maxrow=len(row)
	
	columnnum=maxrow
	for column in range(columnnum):
		f.close() 
		f=open(csvfile+".csv")
		readercsv=csv.reader(f)
		col=extractcolumn(readercsv,column)
		writer.writerow(col)
	f.close()
	writerfile.close()
	os.remove(csvfile+".csv")
	os.rename(csvfile+"_transposed.csv",csvfile+".csv")
	
# Extract the given column from a csv file.				
def extractcolumn(rowgen,colnum):
	column=[]
	for row in rowgen:
		if len(row)>colnum:
			column.append(row[colnum])
		else:
			column.append(float('nan'))
	return column
	
	

	

	
# Posavg wrapper function
# This handles the inputs and outputs the metacodon data to a csv file. 
# It allows use of "all" to look at every amino acid or nt combination; in this case pause scores computed from the average data are also written out.
def posavg_m_wrapper(txtinfile):
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
	
	
	
	inputfiles=list(variables["picklenames"].split(","))
	

	if variables["motif"]=="all":			# This will do all aas or triplets, and compute pause scores in addition to averages.
		pausescores=[]		
		filenames=[]
		motifused=[]
		
		if variables["kind"]=="0":
			motifs=["AAA","AAG","AAC","AAT","AGA","AGG","AGC","AGT","ACA","ACG","ACC","ACT","ATA","ATG","ATC","ATT","GAA","GAG","GAC","GAT","GGA","GGG","GGC","GGT","GCA","GCG","GCC","GCT","GTA","GTG","GTC","GTT","CAA","CAG","CAC","CAT","CGA","CGG","CGC","CGT","CCA","CCG","CCC","CCT","CTA","CTG","CTC","CTT","TAA","TAG","TAC","TAT","TGA","TGG","TGC","TGT","TCA","TCG","TCC","TCT","TTA","TTG","TTC","TTT"]
		elif variables["kind"]=="1":
			motifs=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
		else:
			exit()
			
	else:
		motifs=list(variables["motif"].split(","))
		
	
	# Do score and average data.
		
	writerfile=open(variables["outfile"]+"_avgdata.csv", "w")	# Note no "b" in python 3!
	writer = csv.writer(writerfile)
	
	for motif in motifs:
		
		for picklefile in inputfiles:
			avggene=posavg_m(variables["endmode"],picklefile,variables["outfile"]+"_"+motif+"_"+re.split("[./]",picklefile)[-3],motif,variables["kind"],variables["frame"],variables["bkndwindowthresh"],variables["bkndwindow"],variables["ORFnorm"],variables["UTRmode"],variables["shift"])
		
			#Compute pause score
			if variables["motif"]=="all":	
				bkndwindow=int(variables["bkndwindow"])
				num=sum(avggene[bkndwindow-1:bkndwindow+2])/3		# assume a peak of 3 nt.
				denom=sum(avggene)/len(avggene)	
				pausescores.append(num/denom)
				filenames.append(re.split("[./]",picklefile)[-3])
				motifused.append(motif)
			
			writer.writerow([re.split("[./]",picklefile)[-3]+"_"+motif]+avggene)
		
	writerfile.close()
	transposecsv(variables["outfile"]+"_avgdata")		
	
	
	# Score data not written out unless doing all. There's not really a reason to do just one score.
	if variables["motif"]=="all":		
		writerfile=open(variables["outfile"]+"_score.csv", "w")	
		writer = csv.writer(writerfile)
		writer.writerow(filenames)
		writer.writerow(motifused)
		writer.writerow(pausescores)
		writerfile.close()
		transposecsv(variables["outfile"]+"_score")	

	
	
	
	

	
	
	
	
#### Codon analysis position average analysis for 64 codons
# endmode is either "all_5" or "all_3".
# picklefile is the ribosome footprints file.
# outfile is the name of the output csv file.
# motif is the nt of aa of the motif that is being averaged. Note you can specify the first or last occurence as noted in the code.
# kind is whether the motif is an amino acid (1) or nucleotide sequence (0)
# frame is frame 0, 1, or 2 to be looking in. Put 3 for all 3 frames.
# bkndwindowthresh this is the minimal rpkm counts within the bkndwindow to be included in the average.
# bkndwindow - This is the size of the half-window around the codon of interest to be included in the average.
# UTRmode is whether doing UTR5,CDS,or UTR3 = 0,1,2, respectively.shift
# ORFnorm will normalize to the ORF for UTRs. 0 does not normalize. A positive value will ORF normalize and specify the rpkm ORF threshold for genes to include. This uses a hard-coded shift of +/-13.
# shift is nt to shift data for site of interest (P site, A site, etc).
# Pause score is computed from the average plot, uses a hardcoded window of +/- 1 nt around the peak.
# There is a hard-coded requirement that for ORFnorm, the UTR region cannot be >10x the ORF. 
def posavg_m(endmode,picklefile,outfile,motif,kind,frame,bkndwindowthresh,bkndwindow,ORFnorm,UTRmode,shift):
	
	
	f=gzip.open(picklefile,"rb")
	footprints=pickle.load(f)
	print("Data loaded.")
	f.close()


	kind=int(kind)
	frame=int(frame)
	bkndwindowthresh=int(bkndwindowthresh)
	bkndwindow=int(bkndwindow)
	ORFnorm=float(ORFnorm)
	UTRmode=int(UTRmode)
	shift=int(shift)

	if ORFnorm>0 and UTRmode==1:
		print("Invalid entries. Spike checking prevents iORF usage with ORF norm.")
		exit()



	firstlast=0		# Set to 1 if doing first only; 2 if last only. Note this will check all frames.
	if motif[0]=="f":
		motif=motif[1:]
		firstlast=1
	elif motif[0]=="l":
		motif=motif[1:]
		firstlast=2
	if firstlast>0 and kind!=0:
		print("Invalid use of firstlast on an amino acid motif. Nucleotide only.")
		exit()


	count=0
	i=0
	motiflen=len(motif)
	averagegene=[0 for x in range(2*bkndwindow)]
	
	
	if frame==3:
		framelist=[0,1,2]
	else:
		framelist=[frame]
			
	for frame in framelist:
		genecount=0
		print("Frame "+str(frame))

		for gene in footprints.keys():
			#if gene=="ENSG00000242485.5":
			#	continue	# Skip MRPL20 - weird stuff. This was for old gencode system.
			if gene=="NM_020791.2":
				continue
			if genecount%9000==0:
				print("Genes finished so far="+str(genecount))
			genecount+=1
			ORFstart=int(footprints[gene][3])
			UTR3start=int(footprints[gene][4])
		
		

			# For positive shifts:			
			if shift>=0:
				if UTRmode==0:
					counts=footprints[gene][2][endmode][0:ORFstart-shift]
					genesequence=footprints[gene][1][shift:ORFstart]
					if (ORFstart-shift)<0:		# Check for negative indexes: 
						continue
				elif UTRmode==1:
					counts=footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift]
					genesequence=footprints[gene][1][ORFstart:UTR3start]
					if (ORFstart-shift)<0 or (UTR3start-shift)<0:	# Check for negative indexes:
						continue
				elif UTRmode==2:
					if shift==0:
						counts=footprints[gene][2][endmode][UTR3start-shift:]
					else:
						counts=footprints[gene][2][endmode][UTR3start-shift:-shift]
					genesequence=footprints[gene][1][UTR3start:]
					if (UTR3start-shift)<0:	# Check for negative indexes:
						continue
		
			# FOR negative SHIFTS (3' end aligned):
			if shift<0:
				if UTRmode==0:
					counts=footprints[gene][2][endmode][-shift:ORFstart-shift]
					genesequence=footprints[gene][1][0:ORFstart]
					if (ORFstart-shift)>len(counts):		# Check for indexes off end:
						continue
				elif UTRmode==1:
					counts=footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift]
					genesequence=footprints[gene][1][ORFstart:UTR3start]
					if (ORFstart-shift)>len(counts) or (UTR3start-shift)>len(counts):	# Check for indexes off end:
						continue
				elif UTRmode==2:
					counts=footprints[gene][2][endmode][UTR3start-shift:]
					genesequence=footprints[gene][1][UTR3start:shift]
					if (UTR3start-shift)>len(counts):	# Check for indexes off end:
						continue

		
			if len(counts)==0:
				#Skip gene because shift isn't allowing us to get the full region.
				continue
			
			if shift<0:
				shift_loc=-13
			else:
				shift_loc=13
			ORFlevel=(sum(footprints[gene][2][endmode][ORFstart-shift_loc:UTR3start-shift_loc]))/(UTR3start-ORFstart) #avg rpm per len level						
		
			# Set a value that helps to eliminate poorly annotated UTRs.
			# Set the level at which ORFnorm no longer includes a 3'UTR. Typically this could be a UTR ORF that exceeds the main ORF.
			#spikecheck=(1*ORFlevel)
			spikecheck=(10*ORFlevel)
			#spikecheck=1000000000000000000   # Essentially skips this.
		
			# Take out genes if ORFnorm is being used and you're reading off ends of lists:
			if ORFnorm>0:
				if shift<0:
					if (ORFstart-shift_loc)>len(counts) or (UTR3start-shift_loc)>len(counts):
						continue
				if shift>=0:
					if (ORFstart-shift_loc)<0 or (UTR3start-shift_loc)<0:	# Check for negative indexes:
						continue
				if ((ORFlevel*1000)<ORFnorm): # rpkm conversion and check that ORF is over ORF thresh.
					continue
		
			genelen=len(genesequence)

			
		# Find motifs of interest
		
		
			i=frame
			
			if kind==0:
				while (i+motiflen)<genelen: 
					if motif==str(genesequence[i:i+motiflen]):
						
						if firstlast==1:	# Check if another instance before this starting just past first codon.
							if motif in str(genesequence[1:i]):
								i+=motiflen
								continue
						
						
						if firstlast==2:	# CHeck if another instance of motif is coming in any frame.
							if motif in str(genesequence[i+1:]):
								i+=motiflen
								continue
					
						if (i+bkndwindow)<genelen and (i-bkndwindow)>0:
							bknd=counts[i-bkndwindow:i+bkndwindow]
							
							if ((sum(bknd))/(len(bknd)/1000))<bkndwindowthresh:
								i+=motiflen
								continue
								
							if ORFnorm==0:
								normfactor=sum(bknd)/(len(bknd))	#rpm per length units.
							else:
								normfactor=ORFlevel
								if (sum(bknd))/(len(bknd))>=spikecheck:	# Spike check
									i+=motiflen
								#	print(gene) # For debugging.
									continue

							if normfactor>0:		# Note this will allow 0-read genes for ORFnorm but not nonORFnorm
								for j in range(len(averagegene)):
									averagegene[j]+=((counts[i-bkndwindow+j])/normfactor)	
								count+=1
								# For debugging
								#outputsite=[gene]+[i]+averagegene
								#writer.writerow(outputsite)
							else:
								i+=motiflen
								continue
						else:
							i+=motiflen
							continue
					else:
						i+=motiflen
						continue	
					i+=motiflen

			elif kind==1:
				i=0
				counts=counts[frame:]		# Note that this loss here of a base or two on the end eliminates around 100 genes. 
				genesequence=Seq.Seq(genesequence[frame:]).translate()
				genelen_seq=len(genesequence)
				genelen=len(counts)
	
				while (i+motiflen)<genelen_seq: 
					if motif==str(genesequence[i:i+motiflen]):
					
						if (i*3+bkndwindow)<genelen and (i*3-bkndwindow)>0:
							bknd=counts[i*3-bkndwindow:i*3+bkndwindow]
						
							if ((sum(bknd))/(len(bknd)/1000))<bkndwindowthresh:		
								i+=motiflen
								continue
						
							if ORFnorm==0:
								normfactor=sum(bknd)/(len(bknd))			#rpm per length units.
							else:
								normfactor=ORFlevel
								if (sum(bknd))/(len(bknd))>=spikecheck:	# Spike check 
									i+=motiflen
									continue
								
							if normfactor>0:	# Note this will allow 0-read genes for ORFnorm but not nonORFnorm
						
								for j in range(len(averagegene)):
									averagegene[j]+=((counts[i*3-bkndwindow+j])/normfactor)
								count+=1
							else:
								i+=motiflen
								continue

						else:
							i+=motiflen
							continue
					else:
						i+=motiflen
						continue	
					i+=motiflen
		

	if count==0:
		print("Error 0 count")
		quit()
		

	#Normalize:
	for i in range(len(averagegene)):
		averagegene[i]/=count

		
	print("Number of positions in average")
	print(count)
	return(averagegene)







# Posstats wrapper function
def posstats_m_wrapper(txtinfile):
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
	
	inputfiles=list(variables["picklenames"].split(","))
	for picklefile in inputfiles:
		posstats_m(variables["endmode"],picklefile,variables["outfile"]+"_"+re.split("[./]",picklefile)[-3],variables["motif"],variables["kind"],variables["frame"],variables["genethresh"],variables["pkwindow"],variables["shift"],variables["UTRmode"])
	
	

	
### Codon analysis posstats
# shift is a positive value for 5' end aligned, negative for 3' end aligned. It is amount to shift data before computing score (P site, A site, etc).
# UTRmode is whether doing UTR5,CDS,or UTR3 = 0,1,2, respectively.
# endmode is either "all_5" or "all_3".
# picklefile the input footprint file
# outfile is the name of the csv file for output.
# motif - aa or nt sequence to compute scores at
# kind is 0 for nt or 1 for aa.
# frame is 0, 1, or 2 to be looking for motifs in.
# genethresh - this is the minimal counts (rpkm) in the entire gene to be included in the analysis.
# pkwindow - how much on either side of peak to use in numerator of pause score.
# shift is nt to shift data for site of interest (P site, A site, etc).

def posstats_m(endmode,picklefile,outfile,motif,kind,frame,genethresh,pkwindow,shift,UTRmode):
	f=gzip.open(picklefile,"rb")
	footprints=pickle.load(f)
	print("Data loaded.")
	f.close()
	i=0
	motiflen=len(motif)
	
	kind=int(kind)
	frame=int(frame)
	genethresh=int(genethresh)
	pkwindow=int(pkwindow)
	UTRmode=int(UTRmode)
	shift=int(shift)
	
	writerfile=open(outfile+".csv", "w")	# Note no "b" in python 3!
	writer = csv.writer(writerfile)
	
	outputsite=["headers","alias","mRNA_position","localsequence_nt","localsequence_aa","numerator","denominator","pausescore"]
	writer.writerow(outputsite)
	for gene in footprints.keys():
		
		alias=footprints[gene][0]
		ORFstart=int(footprints[gene][3])
		UTR3start=int(footprints[gene][4])
		

		# For positive shifts:			#Note 0 shifts won't work.
		if shift>=0:
			if UTRmode==0:
				counts=footprints[gene][2][endmode][0:ORFstart-shift]
				genesequence=footprints[gene][1][shift:ORFstart]
				if (ORFstart-shift)<0:		# Check for negative indexes: 
					continue
			elif UTRmode==1:
				counts=footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift]
				genesequence=footprints[gene][1][ORFstart:UTR3start]
				if (ORFstart-shift)<0 or (UTR3start-shift)<0:	# Check for negative indexes:
					continue
			elif UTRmode==2:
				counts=footprints[gene][2][endmode][UTR3start-shift:-shift]
				genesequence=footprints[gene][1][UTR3start:]
				if (UTR3start-shift)<0:	# Check for negative indexes:
					continue
		
		# FOR negative SHIFTS (3' end aligned):
		if shift<0:
			if UTRmode==0:
				counts=footprints[gene][2][endmode][-shift:ORFstart-shift]
				genesequence=footprints[gene][1][0:ORFstart]
				if (ORFstart-shift)>len(counts):		# Check for indexes off end:
					continue
			elif UTRmode==1:
				counts=footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift]
				genesequence=footprints[gene][1][ORFstart:UTR3start]
				if (ORFstart-shift)>len(counts) or (UTR3start-shift)>len(counts):	# Check for indexes off end:
					continue
			elif UTRmode==2:
				counts=footprints[gene][2][endmode][UTR3start-shift:]
				genesequence=footprints[gene][1][UTR3start:shift]
				if (UTR3start-shift)>len(counts):	# Check for indexes off end:
					continue
					
					
					
		
		genelen=len(genesequence)
		# Find motifs of interest
		i=frame
		
		if kind==0:
			while (i+motiflen)<genelen: 

				if motif==str(genesequence[i:i+motiflen]):
					if (i+pkwindow)<genelen:
						numerator=sum(counts[i-pkwindow:i+pkwindow+1])/((2*pkwindow+1)/1000)		
						denominator=(sum(counts))/(genelen/1000)

						
						if denominator<genethresh:
							i+=motiflen
							continue
						if denominator>0:
							pausescore=numerator/denominator
						else:
							pausescore=-1

					else:
						i+=motiflen
						continue
					if (i-9)>=0 and (i+motiflen+9)<genelen:
						localsequence=genesequence[i-9:i+motiflen+9]
						localsequenceaa=Seq.Seq(genesequence[i-9:i+motiflen+9]).translate()
					else:
						localsequence="nolocseq"
						localsequenceaa="nolocseq"
						
					if UTRmode==0:
						localposition=i
					elif UTRmode==1:
						localposition=i+ORFstart
					elif UTRmode==2:
						localposition=i+UTR3start
					else:
						exit()
				else:
					i+=motiflen
					continue	
				outputsite=[gene+"_"+str(i),alias,str(localposition),localsequence,localsequenceaa,numerator,denominator,pausescore]
				writer.writerow(outputsite)
				i+=motiflen
		
		elif kind==1:
			i=0
			counts=counts[frame:]
			genesequencent=genesequence[frame:]
			genesequence=Seq.Seq(genesequence[frame:]).translate()
			genelen_seq=len(genesequence)
			genelen=len(counts)
	
			while (i+motiflen)<genelen_seq: 
				if motif==str(genesequence[i:i+motiflen]):
					
					if (i*3+pkwindow)<genelen:
						numerator=sum(counts[i*3-pkwindow:i*3+pkwindow+1])/((2*pkwindow+1)/1000)		
						denominator=(sum(counts))/(genelen/1000)
						if denominator<genethresh:
							i+=motiflen
							continue
						if denominator>0:
							pausescore=numerator/denominator
						else:
							pausescore=-1
					
					else:
						i+=motiflen
						continue
						
					if (i-9)>=0 and (i+motiflen+9)<genelen_seq:
						localsequenceaa=genesequence[i-9:i+motiflen+9]
						localsequence=genesequencent[i*3-18:i*3+motiflen*3+18]
												
						
					else:
						localsequenceaa="nolocseq"
						localsequence="nolocseq"
						
					if UTRmode==0:
						localposition=i
					elif UTRmode==1:
						localposition=i+ORFstart
					elif UTRmode==2:
						localposition=i+UTR3start
					else:
						exit()
				else:
					i+=motiflen
					continue
					
				outputsite=[gene+"_"+str(i),alias,str(localposition),localsequence,localsequenceaa,numerator,denominator,pausescore]
				writer.writerow(outputsite)
				i+=motiflen
	
	
# Metagenewrapper

def metagene_m_wrapper(txtinfile):
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
	inputfiles=list(variables["picklenames"].split(","))
	writerfile=open(outfile+".csv", "w")	
	writer = csv.writer(writerfile)
	for picklefile in inputfiles:
		metagene=metagene_m(variables["endmode"],picklefile,variables["kind"],variables["equalweighting"],variables["genethresh"],variables["range5"],variables["range3"])
		writer.writerow([re.split("[./]",picklefile)[-3]]+metagene)
	writerfile.close()	
	transposecsv(outfile)	
	
	
# Metagene
# endmode - must be "all_5". Right now endmode3 is not supported for thresholding since shift 13.
# picklefile - file where footprint density is stored
# kind is whether to do average around start codons (1) or stop codons (2)
# weighting is whether to equally weight all genes (1) or to make an unweighted average based on rpm values (0)
# genethresh is the minimal counts (rpkm) in the entire gene to be included in the analysis.
# range values are the 5' and 3' ranges around the position of interest.

def metagene_m(endmode,picklefile,kind,weighting,genethresh,range5,range3):
	if endmode!="all_5":
		print("ERROR - right now the code doesn't support anything except end5. This is easy to fix though with a little more programming.")
		exit()	
	
	# kind is either 1 or 2 (1 for start codons, 2 for stop codons; other features can come later).
	f=gzip.open(picklefile,"rb")
	footprints=pickle.load(f)
	print("Data loaded.")
	f.close()
	
	kind=int(kind)
	weighting=int(weighting)
	genethresh=int(genethresh)
	range5=int(range5)
	range3=int(range3)
	
	growingmetagene=[0 for x in range(range5+range3)]
	
	genecount=0
	
	for gene in footprints.keys():
		if gene=="NM_020791.2":
			continue	# Skip this gene, bad annotation in Refseq+MANE0.6
		tempmetagene=[0 for x in range(range5+range3)]
		
		ORFstart=int(footprints[gene][3])
		UTR3start=int(footprints[gene][4])
		
		# For thresholding and equalweighting - assumes no end effects (typically 30 nt) on ends of genes and that shift is 13 roughly.
		if (UTR3start-ORFstart)<=0 or ORFstart<13 or (len(footprints[gene][2][endmode])-UTR3start)<13:
			continue
		else:
			ORFcounts=(sum(footprints[gene][2][endmode][ORFstart+0-13:UTR3start-0-13]))/((UTR3start-ORFstart-0)/1000)
			if ORFcounts<genethresh:
				continue
		
		if kind==1:
			if ORFstart<range5 or (UTR3start-ORFstart)<range3:
				continue
			else:
				tempmetagene=footprints[gene][2][endmode][ORFstart-range5:ORFstart+range3]
				genecount+=1
			
		elif kind==2:
			if (len(footprints[gene][2][endmode])-UTR3start)<range3 or (UTR3start-ORFstart)<range5:
				continue
			else:
				tempmetagene=footprints[gene][2][endmode][UTR3start-range5:UTR3start+range3]
				genecount+=1
		
		else:
			print("error")
			
		for position in range(range5+range3):
			if weighting==0:
				growingmetagene[position]+=tempmetagene[position]
			else:
				growingmetagene[position]+=(tempmetagene[position]/(ORFcounts/1000))
			
	for position in range(range5+range3):
		growingmetagene[position]/=genecount
		
	print("Genes included = "+str(genecount))	
	
	return growingmetagene



	
# Wrapper function for counting reads on genes.
# Handles file I/O and calling the counter.
def genelist_m_wrapper(txtinfile):

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
	
	doextra=int(variables["doextra"])
	
	outfile=variables["outfile"]
	inputfiles=list(variables["picklenames"].split(","))
	
	
	firstwrite=0
	outputtype=int(variables["outputtype"])
	
	if outputtype==0:
		writerfile=open(outfile+".csv", "w")	
		writer = csv.writer(writerfile)
		if doextra==1:
			writerfile2=open(outfile+"_frame.csv", "w")	
			writer2 = csv.writer(writerfile2)

	else:
		print("Currently this functionality (direct counts output - output 1) is disabled until further notice. File sizes are unwieldy - 300 GB each. Transpose fails. Needs thresholding.")	
			#	writerfile=open(outfile+"_"+re.split("[./]",picklefile)[-3]+".csv", "w")	
			#	writer = csv.writer(writerfile)
			#	for gene in genedict.keys():
			#		writer.writerow([gene]+genedict[gene][2][variables["endmode"]])
			#	writerfile.close()
			#	#transposecsv(outfile+"_"+re.split("[./]",picklefile)[-3])		# Can't do big files with this.
			
	for picklefile in inputfiles:
		
		outdictionaries=genelist_m(variables["endmode"],picklefile,variables["shift"],doextra,outputtype)
		
		if outputtype==0:
			
			genedict=outdictionaries[0]
			
			if doextra==1:
				framedict=outdictionaries[1]
			else:
				framedict={}		# Just a placeholder.
			
			if firstwrite==0:
				growinggenedict=genedict
				growingframedict=framedict
				firstwrite=1
				continue
			else:
				for key in genedict.keys():
					growinggenedict[key]+=genedict[key][6:]
				for key in framedict.keys():
					growingframedict[key]+=framedict[key][1:]
		
		
	for key in genedict.keys():
		writer.writerow([key]+growinggenedict[key])					
	for key in framedict.keys():
		writer2.writerow([key]+growingframedict[key])		
			
	
	writerfile.close()		
	if doextra==1:
		writerfile2.close()
		

		
# This function counts reads on genes.
# endmode - whether you want 5' or 3' end aligned reads. "all_5" or "all_3"
# picklefile - input file of ribosome footprints
# shift - amount to shift data before counting. Typically do for P sites.
# doextra - setting this variable to 1 will output an extra file with CDS reads broken down by frame. Other extras could be added in future. 
# outputtype variable sets whether to count reads on genes (0) or just output the reads for every gene (1). 1 is currently disabled.
# NOTE: the dictionaries for gene and frame and pause are called "list" but they are dictionaries.
def genelist_m(endmode,picklefile,shift,doextra,outputtype):
	if endmode!="all_5":
		print("ERROR - right now the genelist code doesn't support anything except end5.")
	
	samplename=re.split("[./]",picklefile)[-3]
	
	# outputtype is either 0 or 1 (0 for traditional counts that are quantitated gene by gene, or 1 for just whole gene list).
	f=gzip.open(picklefile,"rb")
	footprints=pickle.load(f)
	print("Data loaded.")
	f.close()
	
	shift=int(shift)
	
	genecount=0
	outputtype=int(outputtype)
	if outputtype==1:
		
		return footprints
	else:
	
		genelist={}
		genelist["headers"]=["alias","UTR5len","CDSlen","UTR3len","transcriptlen","seq","UTR5_"+samplename,"CDS_"+samplename,"UTR3_"+samplename,"UTR5raw_"+samplename,"CDSraw_"+samplename,"UTR3raw_"+samplename]
		framelist={}
		framelist["headers"]=["alias","CDS0_"+samplename,"CDS1_"+samplename,"CDS2_"+samplename,"CDS0raw_"+samplename,"CDS1raw_"+samplename,"CDS2raw_"+samplename]
		
		# Code to find total mapped reads for picklefile.
		minoverall=1E6
		for gene in footprints.keys():
			tempgene = [i for i in footprints[gene][2][endmode] if i != 0]
			if len(tempgene)!=0:
				mingene=min(tempgene)
				if minoverall>mingene and mingene!=0:
					minoverall=mingene
		
		mappedreads=int(round(1E6/float(minoverall)))
				
		print("Assuming total mapped reads = "+str(mappedreads))
		
		for gene in footprints.keys():
			#if gene=="ENSG00000242485.5":
			#	continue	# Skip MRPL20 - weird stuff from other genome.
			
			if gene=="NM_020791.2": #This is a gene with known bad annotation.
				continue
			
			ORFstart=int(footprints[gene][3])
			UTR3start=int(footprints[gene][4])
			


			if ((ORFstart-shift)<=0) or (len(footprints[gene][2][endmode])-UTR3start)<=0:		# Gene lacks both UTRs.
				CDScount=""
				UTR5count=""
				UTR3count=""
				CDSraw=""
				UTR5raw=""
				UTR3raw=""
				CDScount_0=""
				CDScount_1=""
				CDScount_2=""
				CDScount_0raw=""
				CDScount_1raw=""
				CDScount_2raw=""
				
				continue		#Note this continue may not be ideal to compare transcriptomes.
			else:
				CDScount=sum(footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift])/((UTR3start-ORFstart)/1000)
				UTR5count=sum(footprints[gene][2][endmode][0:ORFstart-shift])/((ORFstart-shift)/1000)
				UTR3count=sum(footprints[gene][2][endmode][UTR3start-shift:-shift])/((len(footprints[gene][2][endmode])-UTR3start)/1000)
				CDSraw=int(round((CDScount*mappedreads/1E6)*((UTR3start-ORFstart)/1000)))
				UTR5raw=int(round((UTR5count*mappedreads/1E6)*((ORFstart-shift)/1000)))
				UTR3raw=int(round((UTR3count*mappedreads/1E6)*((len(footprints[gene][2][endmode])-UTR3start)/1000)))
				
				CDScount_0=sum(footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift][0::3])/((UTR3start-ORFstart)/1000)
				CDScount_1=sum(footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift][1::3])/((UTR3start-ORFstart)/1000)
				CDScount_2=sum(footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift][2::3])/((UTR3start-ORFstart)/1000)
				CDScount_0raw=int(round((CDScount_0*mappedreads/1E6)*((UTR3start-ORFstart)/1000)))
				CDScount_1raw=int(round((CDScount_1*mappedreads/1E6)*((UTR3start-ORFstart)/1000)))
				CDScount_2raw=int(round((CDScount_2*mappedreads/1E6)*((UTR3start-ORFstart)/1000)))
				
		
				
				genecount+=1
		
			genelist[gene]=[footprints[gene][0],ORFstart,UTR3start-ORFstart,len(footprints[gene][2][endmode])-UTR3start,len(footprints[gene][2][endmode]),footprints[gene][1],0,0,0,0,0,0]
			genelist[gene][6]=UTR5count
			genelist[gene][7]=CDScount
			genelist[gene][8]=UTR3count
			genelist[gene][9]=UTR5raw
			genelist[gene][10]=CDSraw
			genelist[gene][11]=UTR3raw
			
			if genelist[gene][4]>32767:
				genelist[gene][5]="toolongforExcel - lookup in fasta instead"
		
			
			# Output frame information		
			framelist[gene]=[footprints[gene][0],CDScount_0,CDScount_1,CDScount_2,CDScount_0raw,CDScount_1raw,CDScount_2raw]
		
			
		print("Genes included = "+str(genecount))	
	
		return [genelist,framelist]


# Function to find small ORFs in UTRs.
# endmode - whether you want 5' or 3' end aligned reads. "all_5" or "all_3"
# picklefile - input mapped ribosome density file
# shift - amount to shift data before analyzing it, typically for P sites.
# lengththresh - this is the threshold of amino acids. Must be longer than this. Negative is only this.
# mismatches - this is how many mismatches to allow in start codon - 1 or 0. For finding near cognates.
# This is whether to do 5'UTR or 3'UTR (5 or 3).

def dorflist_m(endmode,picklefile,shift,lengththresh,mismatches,UTR):

	if int(lengththresh)<0:		# Allow a neg length to be used absolutely.
		abslengththresh=int(lengththresh)*-1
		lengththresh=0
	
	if endmode!="all_5":
		print("ERROR - right now the genelist code doesn't support anything except end5. This is easy to fix though with a little more programming.")
	
	samplename=re.split("[./]",picklefile)[-3]
	
	f=gzip.open(picklefile,"rb")
	footprints=pickle.load(f)
	print("Data loaded "+samplename)
	f.close()
	
	shift=int(shift)
	genelist={}
	######
	genelist["headers"]=["alias","UTR5len","CDSlen","UTR3len","transcriptlen","smallorfseq","smallorfaa","sm_orfstartinUTR","sm_orfstopinUTR","sm_orf_"+samplename,"CDS_"+samplename]

	dorfcounttotal=0
	genecount=0
	for gene in footprints.keys():
		#if gene=="ENSG00000242485.5":
		#	continue	# Skip MRPL20 - weird stuff from other genome.
		
	#	if gene!="ENSG00000184009.11" and gene!="ENSG00000177606.6":		# Actin and Jun
	#		continue
		
		
		dorfcount=0
		ORFstart=int(footprints[gene][3])
		UTR3start=int(footprints[gene][4])
		
		
		if (ORFstart-shift)<0:	# Check for negative indexes:
				continue
		if (UTR3start-shift)<0:	# Check for negative indexes:
				continue
				
		if int(UTR)==3:
			
			counts=footprints[gene][2][endmode][UTR3start-shift:-shift]
			genesequence=footprints[gene][1][UTR3start:]
		else:	
			counts=footprints[gene][2][endmode][0:ORFstart-shift]
			genesequence=footprints[gene][1][shift:ORFstart]
		
		
		
		genelen=len(counts)

		i=0
		dorfcount=0		# Counting number of dorfs per gene.
		while (i+3)<genelen: 
			testmotif=str(genesequence[i:i+3])
			score=0
			if testmotif[0]!="A":
				score+=1
			if testmotif[1]!="T":
				score+=1
			if testmotif[2]!="G":
				score+=1
			if score>int(mismatches):		# How many mismatches to tolerate.
				i+=1
				continue
				
			# Find in frame stop codon.
				
			j=i+3
				
			while (j+3)<genelen:
				testmotif2=str(genesequence[j:j+3])
				if testmotif2=="TAA" or testmotif2=="TAG" or testmotif2=="TGA":
					endpos=j
					break
				j+=3
			else:
				i+=1
				continue # NO in frame stop.
				
			dorfseq=genesequence[i:j+3]
			dorfaa=Seq.Seq(dorfseq[:-3]).translate()
			dorflen=len(dorfaa)
			if dorflen<int(lengththresh) or (lengththresh==0 and dorflen!=abslengththresh):
				i+=1
				continue
				
			if ((ORFstart-shift)<=0):
				CDScount=-10
			else:
				CDScount=sum(footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift][0::3])/((UTR3start-ORFstart)/1000)
				
			# Note we don't shift since this is already shifted. Was error before.
			dorfreads=sum(counts[i:j+3][0::3])/(((dorflen+1)*3)/1000)
						
						
						
			dorfkey=gene+"_smorf_"+str(dorfcount)
		
			
			genelist[dorfkey]=[footprints[gene][0],ORFstart,UTR3start-ORFstart,len(footprints[gene][2][endmode])-UTR3start,len(footprints[gene][2][endmode]),dorfseq,dorfaa,i,j,0,0]
			genelist[dorfkey][9]=dorfreads
			genelist[dorfkey][10]=CDScount
			dorfcount+=1	
			dorfcounttotal+=1
			i+=1
	
		genecount+=1
		if genecount%1000==0:
			print("Genes finished so far="+str(genecount))		
	print("Total dorfs reported = "+str(dorfcounttotal))	
	return genelist


	

# Wrapper function for finding dORFs.	
def dorflist_m_wrapper(txtinfile):
	tossedgenes=0
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
	inputfiles=list(variables["picklenames"].split(","))
	
	firstwrite=0
	
	
	writerfile=open(outfile+".csv", "w")	
	writer = csv.writer(writerfile)
	
	
	for picklefile in inputfiles:
		dorfdict=dorflist_m(variables["endmode"],picklefile,variables["shift"],variables["lengththresh"],variables["mismatches"],variables["UTR"])		
		
		if firstwrite==0:
			growinggenedict=dorfdict
			firstwrite=1
			continue
		else:
			for key in dorfdict.keys():
				growinggenedict[key]+=dorfdict[key][9:]
		
	startkey="startkey"
	
	towrite="towrite"
	for key in growinggenedict.keys():
		tossedgenescounter=0
		# 0s for all dorfs are gone.
		for i in range(len(inputfiles)):
			if growinggenedict[key][9+i*2]!=0:
				tossedgenescounter=1
				break
		if tossedgenescounter==0:	
			tossedgenes+=1
			continue
		
		
		# CODE TO TAKE LONGEST ISOFORM ONLY. 
		keyrootnum=key.find("_smorf")
		testkey=key[0:keyrootnum]
		if testkey!=startkey:		# New gene (new key).
			startkey=testkey
			currentstop={}
			if towrite!="towrite":
				writer.writerow(towrite)		# Write out old stop codon.
			
			currentstop[growinggenedict[key][8]]=0		# 0 is just a placeholder.
			towrite=[key]+growinggenedict[key]		
			continue				
		if growinggenedict[key][8] in currentstop:	# Check if we are on a new stop codon. If yes, check if length is better and save it.
			if variables["smallest"]==1:		# Check if we are looking for smaller or longer ORFs. 
				if growinggenedict[key][7] > towrite[8]:		# > signals shortest
					towrite=[key]+growinggenedict[key]
				elif growinggenedict[key][7] == towrite[8]:	#Same length should not occur (not possible).
					print("Error")
					exit()
			else:
				if growinggenedict[key][7] < towrite[8]:		# < signals longest
					towrite=[key]+growinggenedict[key]	
				elif growinggenedict[key][7] == towrite[8]: #Same length should not occur (not possible).
					print("Error")
					exit()
	
		else:	
			writer.writerow(towrite)					# New stop codon, get ready to write out the best option and reset.
			currentstop[growinggenedict[key][8]]=0	
			towrite=[key]+growinggenedict[key]
			
	writer.writerow(towrite)							# Write out very last one.



	writerfile.close()
	print("Genes tossed since no reads in any sample.")
	print(tossedgenes)
	
	
	


