###############################################################################
####Title: getIsoforms.py
####Author: Michael Ward
####Purpose: To find unique isoforms from sequence reads.
####Last Modified: 5/2/2024
####
## This program takes a list of sequence reads, extracts the sequence, splices them and
## finds the unique ones. It outputs a file containing the unique isoforms, one containing
## isoforms that are likely subjects for nonsense mediated decay and a report listing the
## sequences read and list if any overlap with one another.
###############################################################################

############################## Import Libraries ###############################
from collections import Counter
from tabulate import tabulate
import itertools as it
import pandas as pd
import numpy as np
import re

############################## Declare Functions ##############################

# get a list of codons from the sequence #
# This function breaks the full sequence into codons(groups of three nucleotides) so they
# can be translated into their amino acid sequence.
def get_codons(gene):                       		# input sequence = gene
    codons = []										# format the list of codons
    start, stop, step = 0, 3, 1 					# set the steps
    while stop <= len(gene):
        codon = gene[slice(start, stop, step)]		# use slice function to get codon
        codons += [codon]							# append to list
        start += 3                            		# update position
        stop += 3
    return codons                             		#returns all codons as type list 

# translate codons into a peptide #
# This function takes a list of codons(groups of three nucleotides) and translates each
# into the appropriate amino acid using the standard (one letter) code. 
def translate(codons):                        		# takes codons from get_codons
    ucode = {                                 		# dictionary of AA one letter codes
    'ttt':'F', 'ttc':'F', 'tta':'L', 'ttg':'L', 'tct':'S',
    'tcc':'S', 'tca':'S', 'tcg':'S', 'tat':'Y', 'tac':'Y',
    'tgt':'C', 'tgc':'C', 'tgg':'W', 'ctt':'L', 'ctc':'L',
    'cta':'L', 'ctg':'L', 'cct':'P', 'ccc':'P', 'cca':'P',
    'ccg':'P', 'cat':'H', 'cac':'H', 'caa':'Q', 'cag':'Q',
    'cgt':'R', 'cgc':'R', 'cga':'R', 'cgg':'R', 'att':'I',
    'atc':'I', 'ata':'I', 'atg':'M', 'act':'T', 'acc':'T',
    'aca':'T', 'acg':'T', 'aat':'N', 'aac':'N', 'aaa':'K',
    'aag':'K', 'agt':'S', 'agc':'S', 'aga':'R', 'agg':'R',
    'gtt':'V', 'gtc':'V', 'gta':'V', 'gtg':'V', 'gct':'A',
    'gcc':'A', 'gca':'A', 'gcg':'A', 'gat':'D', 'gac':'D',
    'gaa':'E', 'gag':'E', 'ggt':'G', 'ggc':'G', 'gga':'G',
    'ggg':'G', 'taa':'', 'tag':'', 'tga':''}
    pep_seq = []									# format the list of peptides
    for codon in codons:
        if codon in ucode.keys():             		# look up AA code
            pep_seq += [ucode[codon]]				# append amino acid to list pep_seq
        else:
            pep_seq += 'o'                    		# unrecognized codon
                                              		# join them all together
    peptides = "".join([str(i.split()) for j,i in enumerate(pep_seq)])
    peptides = re.sub('[\[\]\']','',peptides)		# remove brackets
    
    return peptides                           		# return peptides as str type

# returns coordinates of splice sites #
# the splice sites listed in the sample file are in the format chr:start-stop. This
# function extracts the start and stop location of the splice.
def get_splice_sites(start,prob,up,down):
        PEx=re.findall(r'(\d+)-(\d+)',prob)[0]		# extract coordinate numbers
        sd=re.findall(r'(\d+)-(\d+)',up)[0]			# from the three different columns
        su=re.findall(r'(\d+)-(\d+)',down)[0]
        # create a list of the actual splice locations
        splics_sites=((int(PEx[0])-start,int(PEx[1])-start), \
                      (int(sd[0])-start,int(sd[1])-start), \
                        (int(su[0])-start,int(su[1])-start))
        return splics_sites                      	# return all splice sites
    
# splices the sequence based on splice sites given #
# this function takes the list of splice sites given and removes the nucleotides from the 
# full sequence based on the given splice locations.    
def splice_gene(code,splice_sites):                          
    isoform = ''.join(ch for i, ch in enumerate(code, 1) \
                      if not any(a <= i <= b for a, b in splice_sites))
    return isoform                                  # return spliced sequence

# find positions of stop codons #
# there are three known stop codons taa, tag, and tga this function takes a nucleotide
# sequence for a gene and finds all of these stop codons that appear in frame with the 
#start of the gene sequence.
def findStops(gene):
    codons = np.arange(0, len(gene), 3)		# make a list of codon start positions 
    stopPos = []							# create a list for stop codons
    for i in codons:
        codon = gene[i:i + 3]				# step through the list of codons
        if codon == 'taa' or codon == 'tag' or codon == 'tga': # identify if a stop
            stopPos.append(i + 1)			# add to list of stop codons
    return stopPos							# return the stop positions

# find last exon junction #
# this function uses the position of the splice sites to find all the splice junctions in 
# the spliced sequence and returns the last splice junction for NMD determination.
def getJunc(ssites):
    junc = [a for b in ssites for a in b]	# identify the junctions using splice sites
    ej = sorted(junc)						# add all junctions to a list and sort
    return ej[-1]							# return the last junction

# detect NMD in the current gene #
# Nonsense mediated decay is a translation-coupled mechanism that eliminates mRNAs 
# containing premature translation-termination codons (PTCs). This function looks for 
# two predictors of NMD. It looks for PTCs between 50 and 55 nuclotides upstream of the
# last splice junction. It also looks for sequences where there is a downstream region
# longer than 1000nt past the stop codon.
def detectNMD(segment, ej, stop):
    high,low = (ej-50),(ej-55)				# set region where to look for stop codons
    if any(low <= j <= high for j in stop):	# look for stop codons in this region
        return True							# True: this is a NMD candidate
    elif (len(segment) - stop[1]) > 1000:	# find long downstream regions
        return True							# True: this is a NMD candidate
    else:
        return False						# False: neither is true
    
# find unique isoforms #
# This function compares the different sequences after they have been spliced and returns
# the unique ones. If sequences are spliced differently it returns only the unique ones.
def findUnique(temp):
    unique = {}								# create a dictionary for the unique sequences
    for key, value in temp.items():			# read through the dictionary
        if value not in unique.values():	# identify if it is in the unique dictionary
            unique[key] = value				# if not add them
    return unique							# return a dictionary of unique sequences

# get full sequence without formatting #
# The full sequence comes from  a file that contains annotation and has a sequence that
# is formatted with line number and broken into small segments. This function extracts 
# the sequence from everything else. 
def get_sequence(df): 
    sequence = [] 							# create a list of sequences
    for i in df: 
        if i.startswith("ORIGIN"): 			# scan to the sequence after the annotation
            x = df.read() 					# read the sequence.
            no_num = re.sub("[^a-z\s+]",'', x, 0) 		#remove numbers from sequence
            sequence = re.sub(r'\s+', '', no_num, 0) 	#remove all spaces from sequence
    return sequence							# return unformatted sequence.

# write sequences to file #
# this function formats the output for the sequences to be written to file and then uses
# the write command to write them.
def write_out(count,dictionary,out):
    out.write('>Full Sequence' + '\n' + ori + '\n')		# write the first line
    for n,key,value in zip(count,dictionary.keys(),dictionary.values()): # go through list
        num = count[n]													 # write count
        out.write('>'+str(key[:-2]) +'  ['+str(num)+']  '+'\n'+str(value)[2:-2]+'\n')
        												# write each value to a file
        												
        												
################################ Main Program #################################

# The first stage of this program reads in three files. A file that describes the sequence
# reads, their start and stop positions, splice sites and chromosome. second, it need a 
# genome sequence for the organism. Finally, it need a complete gene sequence for the gene 
# being examined. This part of the code askes for the names of these files and opens them.
# It also names and opens the output files.

   
recount = (input("Please enter recount file name: "))	# give file containing the reads
genefile = (input("Please enter genome file name: "))	# give the file with full genome
original = (input("Please enter gene sequence file name: "))	# file with full sequence
out1,out2,out3 = 'unique_isoforms.fasta','NMD.fasta','report.txt'	# name output files
with (open(genefile, 'r') as seq,open(recount, 'r') as rec, \
      open(original, 'r') as ori,open(out1, 'w') as iso,open(out2, 'w') as nmd, \
      open(out3, 'w') as rep):                           # open all input and output files

# This section initializes some of the variables that will be used later. It sets the 
# count variable to 1 for the total count of processed sequences. It creates two tmp
# dictionaries that will store the sequences that are read. It inports the information on
# these reads into a dataframe where the necessary information can be extracted. The fline
# variable is the length of the first line in the genome file. This file is in fasta 
# format so it has a line that gives a name and other information about the sequence. 
# when 'seek' is used later it must account for this line and does this using this 
#variable. 
    count = 1                    						# initialize count
    tmp1,tmp2 = {},{}									# create temporary dictionaries
    data = pd.read_csv(rec)								# read data as dataframe
    fline = (len(seq.readline()))	# find length of the text at beginning of genome file 
    comp = {'a':'t', 't':'a', 'c':'g', 'g':'c', 'n':'n'}	# list of complement values

# This section starts a loop that will read several values from the file that gives 
# details of where the sequences were read. It reads The name of the read, the chromosome
# it is being read from, the start position, the length, if it is a forward or reverse
# read and the positions of the splice sites. It reads them from the dataframe that was 
# created from the input file (named recount). It also appends a count number to the name
# of each read. 
    for a,b,c,d,e,f,g,h in zip(data.iloc[:,1],data.iloc[:,2],data.iloc[:,3], \
            data.iloc[:,5],data.iloc[:,6],data.iloc[:,22],data.iloc[:,23], \
            data.iloc[:,44]):							# loop through values
        a = "%s:%d" % (a, count)						# add counter to name

# the first conditional that is decided is if this is a '+' or '-' strand. Plus strands
# can be read and processed as they are. Minus strands need to be converted to the 
# reverse complement before processing them further. For both the sequence is read and 
# automatically converted to the lower case. This is regardless of them being '+' or '-'
# the minus strands are converted to reverse complement. The next step for for both is to
# splice the sequences. For this the splice function is called which requires a list of 
# splice sites which it is given by nesting the function get_splice_sites and giving it
# the start position and the listed splice sites from the input file. The next part of 
# processing is NMD detection. For this the detectNMD function is called. This function
# requires splice sites which it gets as the previous function by calling the 
# get_splice_sites function as explained above. and feeds them into the getJunc function
# which identifies the position of the splice junctions on the spliced sequence. It also 
# requires information on stop codons which are identified by the findStops function which
# which needs the spliced sequence that has just been processed and assigned to r_comp for
# minus strands or gene for plus strands. If detectNMD returns 'True' the sequence is added
# to the first temporary dictionary. If it returns 'False' it is added to the second 
# temporary dictionary. The nucleotide sequences are translated as they are added to the 
# dictionary by calling the translate function nested with the get codons function.
# Finally, the count of total sequences processed is incremented.   
        if e == '-':									# if stand is minus
            seq.seek(int(c)+fline)						# adjust for fasta format
            gene = seq.read(int(d))						# read the sequence
            gene = gene.lower()							# make lower case
            r_comp = "".join(comp.get(nt, nt) for nt in reversed(gene))	# get reverse comp
            r_comp = splice_gene(r_comp,get_splice_sites(int(c),f,g,h))	# splie the gene
            if detectNMD(r_comp,getJunc(get_splice_sites(int(c),f,g,h)),findStops(r_comp)) == True:
                tmp1[a] = translate(get_codons(r_comp))	# check if NMD and add to NMD list
                count+=1								# increment count of processed 
            elif detectNMD(r_comp,getJunc(get_splice_sites(int(b),e,f,g)),findStops(r_comp)) == False:
                tmp2[a] = translate(get_codons(r_comp))	# check if NMD add to nonNMD list
                count+=1								# increment counter
                
        elif e == '+':									# if stand is plus
            seq.seek(int(c)+fline)						# adjust for fasta format
            gene = seq.read(int(d))						# read the sequence
            gene = gene.lower()							# make lower case
            gene = splice_gene(gene,get_splice_sites(int(c),f,g,h))		# splie the gene
            if detectNMD(gene,getJunc(get_splice_sites(int(c),f,g,h)),findStops(gene)) == True:
                tmp1[a] = translate(get_codons(gene))	# check if NMD and add to NMD list
                count+=1								# increment count of processed 
            elif detectNMD(gene,getJunc(get_splice_sites(int(c),f,g,h)),findStops(gene)) == False:
                tmp2[a] = translate(get_codons(gene))	# check if NMD add to nonNMD list
                count+=1								# increment counter


############################## Processing Output ##############################
# In the first section of output writing processes the dictionaries containing all of the 
# amino acid sequences processed above is used to create a dictionary of counts of each
# of the sequences. The counter is created at the same time the values are assigned to 
# the count dictionary. 
 
    tmp1count = Counter(tmp1.values())					# get count using Counter function
    tmp2count = Counter(tmp2.values())					# from the collections library
    
# Then the sequence of the original gene is assigned to a variable 'ori' as it is 
# translated.

    ori = translate(get_codons(get_sequence(ori)))		# translate original sequence
    
# the values of the temporary dictionaries are assigned to the isoform and NMD 
# dictionaries while being passed through the find unique function to determine if they
#  are unique values.
 
    isoforms = findUnique(tmp2).copy()					# copy only unique items 
    NMD = findUnique(tmp1).copy()
    
# Once the unique values have been transfered to the proper dictionaries both are written 
# to a separate file using the write_out function. Which adds the number of times the 
# unique sequence was observed. and the sequence of the full gene sequence of interest.

    write_out(tmp1count,NMD,nmd)						# write to output files
    write_out(tmp2count,isoforms,iso)
  
# Next the report is prepared by creating two tables. The first table shows information
# about the data that was processed including, name, chromosome, start position, length 
# and strand '+' or '-'. For this the dataframe output is created and the values are 
# assigned to specific columns.
       
    output = pd.DataFrame()								# create a dataframe 
    output['Name'],output['Chrom'],output['Start'],output['Length'],output['Strand']= \
        data.iloc[:,1],data.iloc[:,2],data.iloc[:,3],data.iloc[:,5],data.iloc[:,6]
        												# assign columns
        
# The first part, of course, declares some variables that will be filled out later.
# The three variables idS, chR and nuMs will contain a list of the names of the reads
# being compared, the chromosome number and the numbers that define the region  of 
# overlap observed. Finally, the variable grouped is a list containing the information
# listed above grouped by chromosome.
        
    idS = []											# what is being compared
    chR = []											# what chromosome
    nuMs = []											# overlap region
    grouped = [v for w, v in output.groupby('Chrom')]	# group by chrome
    
# In the next section the groups are iterated through. Therefore it loops through each 
# based on the chromosome it is on. At the beginning with each loop it assigns ranges as 
# an empty dictionary. Then the length of the group is checked. If a group has less than
# two entries there is nothing to compare. If there are 2 or more the range from start 
# position to end position is assigned to the ranges dictionary just created. This is 
# assigned as a list so every values where the two overlap. 

    for group in grouped:
        ranges = {}										# create ranges dictionary
        if len(group) < 2:								# less than two nothing to do
            pass
        elif len(group) >= 2:							# greater than two
            for a,b,c,d in zip(group.Name,group.Chrom,group.Start,group.Length):
                ranges[a] = list(range(c,d+1))			# expand the range of the region
                
# In this part I want to determine where these reads overlap. For example, if you
# find that two reads read the same sequence and overlap 100% but you get two different
# unique isoforms you know they are spliced diffwerently. I begin by getting a list of 
# keys from the ranges. Each key is the name of a specific read. Then using the 
# combinations function from the itertools library I creat a list of all possible 
# combinations then by making the ranges into sets I use the '&' operator to look for
# set intersection. If the intersection exists, the chromosome, start and stop positions 
# of the intersection and the names of the sequences being compared are assigned to lists.
# If there is no intersection the names of what was compared and chromosome are assigned
# to lists but the word 'None' is placed in the list for where the overlap occurs.
                
        keyList = list(ranges.keys())				# get names of keys
        for (x,y) in it.combinations(keyList, 2):	# go through all unique combinations
            L1,L2 = ranges[x],ranges[y]			
            inT = list(set(L1) & set(L2))			# look for intersection
            if len(inT) > 0:						# if there is an intersection
                chR += [str(b)]						# record chromosome names
                left,right = min(inT),max(inT)		# find begin and end of intersection
                idn = (str(x),',',str(y))			# record what is compared
                nums = (str(left),'-',str(right))	# create string of range
                idS += [''.join(idn)]				# convert to string
                nuMs += [''.join(nums)]
            elif len(inT) == 0:						# if no overlap
                chR += [str(b)]
                idn = (str(x),',',str(y))
                idS += [''.join(idn)]
                nuMs += ['None']					# record no intersection
          
# Prepare output table of the intersection data by making a dataframe and assigning the 
# above lists to columns
          
    output2 = pd.DataFrame()						# make dataframe
    output2['Chrom'],output2['Reads'],output2['Overlap'] = chR,idS,nuMs	# assign lists
    
# Write the first table to file using the tabulate function
    
    rep.write('Table of locations of the sequences being read.\n')	# write caption
    rep.write(tabulate(output, headers = 'keys', tablefmt = 'fancy_grid')+ '\n\n')
    																# write table
    
# Write the second table to file   
    
    rep.write('Table of overlaps in the sequences being read.\n')	# write caption
    rep.write(tabulate(output2, headers = 'keys', tablefmt = 'fancy_grid'))
    																# write table