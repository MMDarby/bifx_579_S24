#####################################################
####Title: Isoform_Finder
####Author: Michael Ward
####Purpose: This program searches through sequences and finds the unique
####          sequences. It also count how many times each unique one occurs.
####Last Modified: March 29,2024
####
## This program reads a fasta file and obtains the unique sequences, the number
## of each unique sequence and produces a new fasta file with only the unique 
## sequences that are annotated with the number of times that sequence occurred
####################################################

################# Declare Functions #################

# check type of sequence, returns True if protein sequence #
def checkType(sample):                        # expects a sample sequennce
    ntch = 'actgnACTGN'                       # letters in nt sequences
    aach = 'defhiklmpqrsvwxyDEFHIKLMPQRSVWXY' #letters in aa sequences
    for l in sample:
        if l in ntch:                         #check for matches
            return False
        elif l in aach:
            return True                       # return boolean

# get a list of codons from the sequence, breacks nt sequence into codons #
def get_codons(gene):                         # input sequence = gene
    codons = []
    start, stop, step = 0, 3, 1 
    while stop <= len(gene):
        codon = gene[slice(start, stop, step)]# use slice huntion to get codon
        codons += [codon]
        start += 3                            # update position
        stop += 3
    return codons                             #returns all codons as type list 

# translate codons into a peptide #
def translate(codons):                        # takes codons from get_codons
    ucode = {                                 # dictionary of AA one letter codes
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
    'ggg':'G', 'taa':'X', 'tag':'X', 'tga':'X'}
    pep_seq = []
    for codon in codons:
        if codon in ucode.keys():             # look up AA code
            pep_seq += [ucode[codon]]
        else:
            pep_seq += 'o'                    # unrecognized codon
                                              # join them all together
    peptides = "".join([str(i.split()) for j,i in enumerate(pep_seq)])
    return peptides                           # return peptides as str type

################# Main Program #################
temp = {}                             # a dictionary used to hold all sequences
isoforms = {}                  # a dictionary used to hold the unique sequences
reads = 'sample_data.nt.txt'                             # name of file to read                                        
with open(reads, 'r') as seq, open('out.txt', 'w') as out: # open above and out

## this section reads the sequence into a list called lines then it checks to
## see if it is a nucleotide or amino acid sequence. If a nucleotide sequence
## it translates it using the get_codons and translate functions. It then
## coonverts it to a sting and replaces the nt sequence in lines

    lines = [line.rstrip() for line in seq]
    for n, line in enumerate(lines):                #determine type of sequence
        if not (line.find('>')==0) and checkType(line) == False:
            trans = translate(get_codons(line))     # translate nt sequences
            trans = (trans.translate({ord(i): None for i in '[]\''}))
            lines[n] = trans                        # convert to str & replace

## since I dont know how many unique sequences there will be I create a 
## dictionary where new values can be added easily. Every time it finds a 
## unique sequence it adds it to the dictionary and starts counting.
    
    count = {i:lines.count(i) for i in lines[1::2]}    #create count dictionary

## this converts the list lines to a dictionary for further processing
## so values can be identified by keys.
            
    for u in range(0,len(lines)-1,2):               #create sequence dictionary
        temp[lines[u]] = lines[u+1:u+2]             # add each value in lines        

## this creates a seperate dictionary for the unique values called (isoforms)
## this will be the final output allong with count.               
        
    for key, value in temp.items():                      #find unique sequences
        if value not in isoforms.values():                 # check if its unique
            isoforms[key[-1:]] = value                   #add unique values
            
## this selects the values that will be written to a file and loops throught
## them until thay all have been written.         
            
    for n,key,value in zip(count,isoforms.keys(),isoforms.values()):    #output
        num = count[n]
        out.write('>'+str(key)+'  ['+str(num)+']  '+'\n'+str(value)[2:-2]+'\n')
        
## the outpus is formatted as:
## >Nmae [total number of matching sequences]
## sequence in AA format