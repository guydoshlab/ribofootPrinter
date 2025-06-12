# ribofootPrinter
Python code for analysis of ribosome profiling data that has been aligned to a simplified transcriptome.
7/5/2021

Note that an earlier version of this code was placed in the ReducedTranscriptome repository.

The code in this Python3 repository can be used for analyzing ribosome profiling data from any organism with an annotated transcriptome. Footprint data needs to be first mapped to a set of transcripts (1 isoform per gene), such as those provided by the MANE project, with bowtie. An example of the fasta file header format to use is: 

>NM_002046.7|NM_002046.7|NM_002046.7|NM_002046.7|NM_002046.7|GAPDH|1285|UTR5:1-76|CDS:77-1084|UTR3:1085-1285

(Note that 1285 is the length of this mRNA isoform for the GAPDH gene).

The output SAM files can then be used by mammalian_builddense.py to create data structures in Python to hold the sequence and mapped footprint information for each gene (both 5' and 3' end assignments). These are then pickled in python and saved to disk. 

The code in tools_m.py can then be used with associated input text files to analyze this data for quantitation (counting reads), metagene, metacodon, or localized pause analysis. The code can also be used to determine dORFs and to write out footprint density from a spliced gene to file. 

Most of these tasks can be performed by using the associated "wrapper" function to read in the input file and pickled ribosome occupancy file.

Input files:

metagene.txt
Creates an average of reads from all genes after aligning spliced transcripts by start codon or stop codon. There are various options for weighting and thresholding available, as well as sizes of windows.

genelist.txt
Counter that adds reads up for various regions of transcripts (5'UTR, main ORF, 3'UTR), taking shift for P site, A site, etc into account.

posavg.txt
Multifaceted function that computes average reads in a region around particular sites of interest (such as AUG codons or Pro amino acids), or all permutations of a given length, such as 3 amino acids. Options available for normalization or thresholding.

posstats.txt
A function that computes pause scores (peak height/background) for particular sites of interest. Options available for choice of the windows around the peak.

writegene2.txt
Code to write out reads mapping to a spliced gene in a csv file, mainly for presentation purposes.

smorflist.txt
Experimental function that helps find small ORFs in UTR regions by reporting in-frame reads in the small ORF and main ORF. Relies on well-RNased ribosome footprint data and precise shift values.


Added 6/12/25

The nanopore_adaptation code is a version of mammalian_builddense that can be used with aligned nanopore data to ask what sequences exist in the transcriptome in locations near 5' ends. It also reports lengths of mapped reads (both absolute and normalized to the transcriptome).

The fastq_functions code includes seqtools_m.py. It offers the capability to: 

countreadsizes: 
Report the sizes of a fastq file and, if desired, filter for reads with a particular motif on the 3' end.  

endmotif_distribution: 
Reads in a fastq file and counts the number of sequence motifs (of specified length) at either 5' or 3' end.


