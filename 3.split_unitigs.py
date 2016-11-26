#!/usr/bin/python
from pyfasta import Fasta
import csv
import sys

def ReadTable(inputStr, header=0, delim='\t'):
    """
    Takes table in csv format, skips header line if needed and returns an array of arrays
    """
    rows = []
    with open(inputStr, 'rb') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=delim)
        if (header > 0):
            #print 'Skipping', header, 'lines of header from file', inputStr
	    for i in range(header):
                csvreader.next()
        for row in csvreader:
            rows.append(row)
    return rows

def main():
    """
    cuts fasta file at specific location
    """
    if len(sys.argv) == 2:
        prefix = sys.argv[1]
    else:
        print "Usage: python split.py <prefix>; assume that <prefix>_BspQI_key.txt <prefix>.fasta and <prefix>_cut_list.csv exist; output will be <prefix>_new.fasta; cut_list is <contigID>,<loc1>,[<loc2>] -- first line is scaling constant"
        return 0
     
    ren = ReadTable(prefix+'_BspQI_key.txt', 4, '\t') # 4 lines of header 
    #print ren
    cut = ReadTable(prefix+'_cut_list.csv', 0, ',') # no header (saved as MS-DOS csv via Excel)
    #print cut 

    # create a dictionary between contig id and FASTA id x[1] and FAST len x[2]
    renaming = {}
    for x in ren:
        renaming[int(x[0])]=(x[1],int(x[2])) # contigs names are converted into integers, as well as length
    #print renaming
  
    # collect the names of the contigs to be cut 
    location = {}
    scaling = float(cut[0][0])
    print 'scaling constant',scaling
    for x in cut[1:]:
        index = int(x[0]) # name of the contig to cut, convert contig name into into integer so we can match it
        if index in renaming:
	    if (len(x) == 2):
           	l = int(round(float(x[1])/scaling)) # position to cut
            	if (l > renaming[index][1]): # check the length
                    print 'Error: cannot split contig',index,'at position',l,'because it is only',renaming[index][1],'bp long'
                    sys.exit(-1)
                else:
                    location[renaming[index][0]]=[l] # location[contig_name]->position
	    elif (len(x) == 3):
  		l1 = int(round(float(x[1])/scaling)) # position to cut
		l2 = int(round(float(x[2])/scaling)) # position to cut
                if (l1 > renaming[index][1]) or (l2 > renaming[index][1]): # check the length
                    print 'Error: cannot split contig',index,'at position',l1,l2,'because it is only',renaming[index][1],'bp long'
                    sys.exit(-1)
                else:
                    location[renaming[index][0]]=[l1,l2] # location[contig_name]->position
        else: 
           print 'Error: contig',index,'does not exist'
           sys.exit(-1)
    print location
 
    # open the fasta file for reading
    fas = Fasta(prefix+'.fasta')
    # open the new fasta file for writing
    ofa = open(prefix+'_new.fasta','w')
    for x in sorted(fas.keys()): # process all the contigs one by one
        if x in location: # if it needs to be split
	    if len(location[x]) == 1:
		l = location[x][0]
            	print 'Splitting',x,'at location',l
            	ofa.write('>'+x+'|chimeric1\n')
            	ofa.write(fas[x][:l]+'\n') # prefix
            	ofa.write('>'+x+'|chimeric2\n')
            	ofa.write(fas[x][l:]+'\n') # suffix
	    elif len(location[x]) == 2:
                l1 = location[x][0]
                l2 = location[x][1]
		if (l1 > l2):
		    temp = l2
		    l2 = l1
		    l1 = temp
                print 'Splitting',x,'at location',l1,'and',l2
                ofa.write('>'+x+'|chimeric1\n')
                ofa.write(fas[x][:l1]+'\n') # prefix
                ofa.write('>'+x+'|chimeric2\n')
                ofa.write(fas[x][l1:l2]+'\n') # middle
                ofa.write('>'+x+'|chimeric3\n')
                ofa.write(fas[x][l2:]+'\n') # suffix
        else: 
            #print 'Not splitting',x
            ofa.write('>'+x+'\n')
            ofa.write(fas[x][:]+'\n')
    ofa.close()
 
if __name__ == '__main__':
    main()
