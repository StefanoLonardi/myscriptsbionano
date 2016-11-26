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
            print 'Skipping', header, 'lines of header from file', inputStr
	    for i in range(header):
                csvreader.next()
        for row in csvreader:
            rows.append(row)
    return rows

def main():
    """
    select specific contigs from FASTA file
    """
    if len(sys.argv) == 2:
        prefix = sys.argv[1]
    else:
        print "Usage: python select.py <prefix>; assume that <prefix>_BspQI_key.txt <prefix>.fasta and <prefix>_list.txt exist; output will be <prefix>_selected.fasta"
        return 0
     
    ren = ReadTable(prefix+'_BspQI_key.txt', 4, '\t') # 4 lines of header 
    print 'renaming table',ren
    select = ReadTable(prefix+'_list.txt', 0) # no header, text file of contigs numbers, one per line
    print 'select list',select

    # create a dictionary between contig id x[0] and (FASTA id x[1])
    renaming = {}
    for x in ren:
        renaming[int(x[0])]=x[1] # contigs names are converted into integers, as well as length
    print 'renaming dictionary', renaming
  
    # collect the names of the contigs to be cut 
    selected_list = []
    for x in select:
        index = int(x[0]) # name of the contig to select, convert contig name into integer so we can match it
        #print 'index',index
        if index in renaming:
           selected_list.append(renaming[index]) # add the name of the contig
        else: 
           print 'Error: contig',index,'does not exist'
           sys.exit(-1)
    print 'selected_list', selected_list
 
    # open the fasta file for reading
    fas = Fasta(prefix+'.fasta')
    # open the new fasta file for writing
    ofa = open(prefix+'_new.fasta','w')
    print 'writing new fasta'
    for x in sorted(fas.keys()): # process all the contigs one by one
        if x in selected_list: # if it needs to be split
            print 'Selecting',x
            ofa.write('>'+x+'\n')
            ofa.write(fas[x][:]+'\n') # entire contig
        else: 
            print 'Not selecting',x
    ofa.close()
 
if __name__ == '__main__':
    main()
