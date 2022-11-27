#!/usr/bin/env python
##                      WORKSHOP BIOINFORMATICS M1 International Track - Silvia Fajardo and Andressa Leao 
##                                             X. Tropicalis Chromosome 2

###FIRST STEP - EXTRACT THE CHROMOSOME OF INTEREST 
### Extract fasta sequences from a fasta file based on a string in the name line
### python ExtractChromosome.py string (Chromosome 2) filename (X. Tropicalis.fa)

import re,os,sys
# I suppress the import of Bio module because it is not installed.

def ask():
    
    filename=raw_input('Seqfile name ? : ').split()
    string=raw_input('string ? : ')
    
    extract(string, filename)
    
def extract(string, filename):  
    seqfile=open(filename).read()
    seqlist=seqfile.split('>')[1:] 
    newlist=[]
    for x in seqlist:
        if string in x.split('\n')[0]:
            newlist.append(x) 
    with open('.'.join(filename.split('.')[:-1])+'_'+string+'.fa','w') as outfile:
       outfile.write('>'+'>'.join(newlist))
        
if __name__ == '__main__' : 
    if len(sys.argv)==3:
        # string, filename
        print(sys.argv[1],sys.argv[2])
        extract(sys.argv[1],sys.argv[2])
    else: ask()

##Output = X.tropicalis_Chr_2.fa
#edit the output manually - adding an ">" at the beginning of the fasta file 
#the script was edited by the professor - she added '>'+'
#Next step - trf :) 






###SECOND STEP - TRF
###Input = X.tropicalis_Chr_2.fa
## /Users/silviafajardo/Desktop/WS-repeat/Scripts/trf409.macosx

## Trf X.tropicalis_Chr_2.fa 2 5 7 80 10 50 200 -f -d -h

## ./Scripts/trf409.macosx Data/X.tropicalis_Chr_2.fa 2 5 7 80 10 50 200 -f -d -h

###Output=X.tropicalis_Chr_2.fa.2.5.7.80.10.50.200.dat 







###THIRD STEP - This step is totrf2fasta4.py
##WE USED THIS: Scripts/trf2fasta4.py Results /Users/silviafajardo/Desktop/WS-repeat/Results/X.tropicalis_Chr_2.fa.2.5.7.80.10.50.200.dat 30 

#!/usr/bin/env python
# Please ensure that the script can be found in a folder named Scripts
# Nicolas.Pollet@egce.cnrs-gif.fr , aurelie.hua-van@universite-paris-saclay.fr
# 
# Command line if you're in the script directory
# python3 trf2fasta4.py Datfile minsize

# 
# Command line if you're in the WS_repeat directory
# python3 Scripts/trf2fasta4.py Datfile minsize

# datfile = fastafile.dat
# minsize= for filtering
# The fasta file must be in your working directory because trf saves its output file in the working directory
import os, sys, argparse, subprocess

# check if the file can be found and return the absolute path
def Check_file(infile):
    if os.path.isfile(infile):
        return(os.path.abspath(infile))
    else: print('Cannot find the file, please try again')
# This function returns one satellite type depending of the size of the motif (period size)
# To what data correspond x?
def get_type(x):
    if int(x)<=6: return("microsatellite")
    elif 6<int(x)<=100: return("minisatellite")
    elif int(x)>100: return("satellite")
    
# This analyse all overlapping hits of a genomic region, and return the one with the best alignment score
# The argument is a list of list
def ProcessContainer(container):
    # print(container)
    
    # which hit among the overlapping ones have the best alignment score
    score=[int(i[7]) for i in container]
    best_score_index=score.index(max(score))
    
    # which hit among the overlapping ones have the smaller motif
    motif=[int(i[2]) for i in container]
    smaller_motif_index=motif.index(min(motif))
    
    # list of sat type overlapping
    type=set([i[-1] for i in container])
    
    # it returns the hit with the smaller motif
    return(container[smaller_motif_index])
    
    # it returns the hit with the best score
    #return container[best-score_index]
    
def read_trf_file(trffile,size=0):
    #open and read the file
    trf=open(trffile).read()
    # Create output files, for non filtered and non reduundant
    # fasta files and tab files
    nfoutfile=open(trffile+'.nf.fa','w')
    nfoutfile2=open(trffile+'.nf.tab','w')
   
    nroutfile=open(trffile+'.nr.tab','w')   
    
    foutfile=open(trffile+'.f.fa','w')
    foutfile2=open(trffile+'.f.tab','w')  

    ### Split the file content by sequences
    chromosomes=trf.split('Sequence: ')[1:]
    a,b,c=0,0,0  # some counters
    dict_chr={}
    
    ### For each sequence, split by hit
    for chr in chromosomes:
        # Extract the chromosome name
        chrname=chr.split('\n')[0].split(' ')[0]
        # Store the hits as a list (lines) of list (the values of each column)
        # each item 
        restrf=[line.split(' ')  for line in chr.split('\n')[7:-3]   ]
        # Discard empty lists, and calculate the number of hits
        if len(restrf) >0:
            print( 'trf found ',len(restrf),' tandem repeat sequences on ',chrname)
        else: continue
        rlen=0
        nrlen=0
        start,startf=0,0
        end,endf=0,0
        flag = True
        container=[]

        # for each hit (provided as a list of value)
        for j,i in enumerate(restrf):
            ## i[0], i[1]  = start and end of the sat on the chromosome
            ## i[2] :size of the motif (period)
            ## i[3] : number of repeats
            ## i[7] : the alignment score
            ## i[13] : the  sequence  of the motif
        
            # call the get_type function
            i.append(get_type(i[2]))
            #add the type to the list of information for the hit
            container.append(i)
            # save fasta and tab,  non filtered (nf.tab and. nf.fa)
            # chr, posX, posY, motif_size, repeat, trf_id, repeat_type
            nfoutfile.write('>'+chrname+'_trfid'+str(a)+'_'+i[0]+'-'+i[1]+'_'+i[2]+'_'+i[3]+'\n'+i[13]+'\n')
            nfoutfile2.write(chrname+'\t'+i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\ttrfid'+str(a)+'\t'+get_type(i[2])+'\n')
            
            
            # save fasta and tab, filtered if requested (argument size > 0)         
            if len(i[13])>=int(size) and int(size)>0: #size of  the motif above the threshold  
                foutfile.write('>'+chrname+'_trfid'+str(a)+'_'+i[0]+'-'+i[1]+'_'+i[2]+'_'+i[3]+'\n'+i[13]+'\n')
                foutfile2.write(chrname+'\t'+i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\ttrfid'+str(a)+'\t'+get_type(i[2])+'\n')  
                c+=1 # add 1 to the counter c (== nummber of hit after filtration by size)
                
            #  Count the nb of nt  and size          
            a+=1
            rlen += int(i[1])-int(i[0])
            keep=j
            if j+1==len(restrf): # it is the last line
                nrlen+=endf-startf
                continue   # go directly to the next for-loop (which does not exist since this is the end of the file!)
                
            else: # calculate some distances with the next hit  
                minmotifsize=min(int(i[2]),int(restrf[j+1][2]))
                overlap =  int(restrf[j+1][0])  -   int(i[1]) # overlap if negative
                shiftx = abs(int(i[0])-int(restrf[j+1][0]))
                shifty = abs(int(i[1])-int(restrf[j+1][1]))
                
            if overlap <0 and (shiftx < minmotifsize or shifty < minmotifsize): # large overlap with the next hit
                container.append(i)
                            
            elif  overlap >=0: # no overlap with the next hit -->  analyse thee previous set
                if len(container)>0: # if the previous hits are overlapping, choose the best hit to keep, save it, empty the container
                    keep=ProcessContainer(container)
                    nroutfile.write(chrname+'\t'+keep[0]+'\t'+keep[1]+'\t'+keep[2]+'\t'+keep[3]+'\ttrfid'+str(a)+'\t'+keep[-1]+'\n')  
                    b+=1
                    container=[]


        print(chrname,rlen, nrlen)           
    foutfile.close()
    foutfile2.close()
    nfoutfile.close()
    nfoutfile2.close()
    nroutfile.close()
    print()
    print( 'Found ',a,' tandem repeat sequences in  total')
    if int(size)>0:
        print( 'Found ',c,' tandem repeat sequences with motif >= ',size, 'in  total')
    print( 'Found ',b,' non redundant tandem repeat sequences in  total')
    return(trffile+'.fa')


infile=full_path_fasta_file=Check_file(sys.argv[1])
read_trf_file(infile,sys.argv[2])








###CHROMOSOME 2 LENGTH 
##WE USED THIS: ##python Scripts/Seqlen.py -g Data/X.tropicalis_Chr_2.fa 

import sys, argparse,os

def sequence_size(genome):
    genome_dict={}
    with open(genome, 'r') as F:
        for line in F:
            if line.startswith('>'):
                chrname=line.strip('>\n').split(' ')[0]
                genome_dict[chrname]=0
            else:
                genome_dict[chrname]+=len(line.rstrip())
    for record in genome_dict:
        seq_length=str(genome_dict[record])
        seq_id=record
        sys.stdout.write(''.join([seq_id,'\t',seq_length,'\n']))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="genome length statistics.")
    parser.add_argument("-g", "--genome", type=str, required=True, help="The reference genome.")
    Args = parser.parse_args()
    
    sequence_size(os.path.abspath(Args.genome))







###TO CREATE THE MAP FILE 
##WE USED: python Scripts/trf_count.py Data/X.tropicalis_Chr_2.chrmsize /Users/silviafajardo/Desktop/WS-repeat/Data/X.tropicalis_Chr_2.fa.2.5.7.80.10.50.200.dat.nr.tab 1000000 150000 


#!/usr/bin/env python
# Please ensure that perl scripts can be found in a folder named script
# Nicolas.Pollet@egce.cnrs-gif.fr , aurelie.hua-van@universite-paris-saclay.fr
# 
# Usage : python trf_count.py sequence_sizefile tabfile windowsize overlap

import os, sys, subprocess
from collections import defaultdict


# From the size file (size of the chromosomes), create a  dictioaary associating chromosome with its size
def SizeSeq(Sizefile):
    chrsize=open(Sizefile).readlines()
    dict_chr={}
    for chr in chrsize:
        #print(chr.rsplit())
        dict_chr[chr.split('\t')[0]]=int(chr.split('\t')[1])
    return(dict_chr, Sizefile)

            
# classify all hits found by trf according to the motif's size
# microsat: <=6
# minisat: between 7 and 100
# sat: over 100
# return three dictionaries  
def Categorize_repeat(datafile):
    microsat,minisat,sat=defaultdict(list),defaultdict(list),defaultdict(list)
    
    repeat=open(datafile).readlines()
    # Chr, trfid, x, y, motifsize, repeat
    for i in repeat:
        i=i.split('\t')
        if int(i[3])<=6: microsat[i[0]].append(i)
        elif 6 < int(i[3]) < 100: minisat[i[0]].append(i)
        elif  int(i[3]) >= 100: sat[i[0]].append(i)
    return(microsat,minisat,sat)
        
# Compute the dictionaries by  sliding windows
### ----------------
###     ----------------
###         ----------------
### default window size= 1000000
### default overlap size= 150000
def process_satellite(dict_chr, cat, sizefile, datafile, windowsize=1000000,overlap=150000):
    open(datafile+'.map', 'w').close()
    
    #create  a file with  the window and shift as a  bed format
    windowfile=os.path.abspath(sizefile)+'.'+str(windowsize)+'.'+str(overlap)+'.bed'
    open(windowfile, 'w').close()
    print("File ", datafile+'.map', " has been created")

    
    # format the 2 last argument as number   (integer)
    windowsize=int(windowsize)
    overlap=int(overlap)
    
    # for each chr
    for chr in dict_chr:
        print(chr, "is processed")
        coordx=0
        coordy=windowsize
        while coordy < dict_chr[chr]:
           write2(windowfile,coordx, coordy,  chr)
           micro,mini,sat='','',''
           m,n,o=0,0,0
           q,r,s=0,0,0
           # for each microsat on this chromosome
           for each in cat[0][chr]:
               if coordx  <= int(each[2]) < coordy:
                   m+=1
                   q+= int(each[2])-int(each[1])
           if m != 0   : write(datafile, coordx, coordy, chr, str(q), str(m), 'microsatellite')
                
           # for each minisat on this chromosome                 
           for each in cat[1][chr]:
               if coordx  <= int(each[2])<coordy:
                   n+=1
                   r+= int(each[2])-int(each[1])
           if n != 0   : write(datafile, coordx, coordy, chr, str(r), str(n), 'minisatellite')

           # for each sat on this chromosome                 
           for each in cat[2][chr]:
               if coordx  <= int(each[2])<coordy:
                   o+=1
                   s+= int(each[2])-int(each[1])                                
           if o != 0   : write(datafile, coordx, coordy, chr, str(s), str(o), 'satellite')
           
           
           coordx += overlap
           coordy += overlap
    
# Write output in thee same directory as the datafile
def write(datafile, coordx, coordy, chr, length, number, type):
    with open(datafile+'.map','a') as outfile:
        outfile.write('\t'.join(map(str,[chr, coordx, coordy, length, number, type]))+'\n')
        
def write2(newfile,coordx, coordy,  chr):
    with open(newfile,'a') as outfile2:
        outfile2.write('\t'.join(map(str,[chr, coordx, coordy]))+'\n')    
  
if __name__ == '__main__' : 
    dict_chr=SizeSeq(sys.argv[1])
    Categories=Categorize_repeat(sys.argv[2])
    process_satellite(dict_chr[0], Categories, sys.argv[1],sys.argv[2], sys.argv[3],sys.argv[4] )
    
#OK	
#nf.tab






