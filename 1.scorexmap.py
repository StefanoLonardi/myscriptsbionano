#!/usr/bin/python
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
#            print 'Skipping', header, 'lines of header from file', inputStr
	    for i in range(header):
                csvreader.next()
        for row in csvreader:
            rows.append(row)
    return rows

def main():
    """
     Reads xmap/cmap files produced by RefAligner and score it
    """
    if len(sys.argv) == 2:
        myfile = sys.argv[1]
    else:
        print "Usage: python scorexmap.py <file>, where <file> is the prefix, scorexmap will open <file>.xmap and <file>.chimeric"
        return 0

    # process the reference map to see what fragments have not been mapped
    mreflen = {}
    refaligned = {}
    table0 = ReadTable(myfile+'_r.cmap', 11, '\t') # 11 lines of header
    for row in table0:
        refi = int(row[0]) # integer
        refl = float(row[1]) # float
        mreflen[refi] = refl
        refaligned[refi] = False

    # process the query map to see what fragments have not been mapped
    mqrylen = {}
    qryaligned = {}
    table2 = ReadTable(myfile+'_q.cmap', 11, '\t') # 11 lines of header
    for row in table2:
        qryi = int(row[0]) # integer
        qryl = float(row[1]) # float
        mqrylen[qryi] = qryl
        qryaligned[qryi] = False

    # process the xmap     
    table1 = ReadTable(myfile+'.xmap', 10, '\t') # 10 lines of header 
    s = 0.0
    min_confidence = 20.0
    minrefoverhang = 100000
    minqryoverhang = 100000
    chimeric = 0
    for row in table1:
        qry = int(row[1]) # integer
        ref = int(row[2]) # integer
        qrystartpos = float(row[3]) # float
        qryendpos = float(row[4]) # float
        refstartpos = float(row[5]) # float
        refendpos = float(row[6]) # float
        orientation = row[7] # '+' or '-'
        confidence = float(row[8]) # float
        qrylen = float(row[10]) # float
        reflen = float(row[11]) # float
	if (confidence > min_confidence):
            refaligned[ref] = True # we mark this ref moleculed aligned
       	    qryaligned[qry] = True # we mark this qry contig aligned
        ref_left_overlen = refstartpos
        ref_right_overlen = reflen - refendpos
        if (orientation == '+'):
            qry_left_overlen = qrystartpos
            qry_right_overlen = qrylen - qryendpos
        else:
            qry_left_overlen = qrylen - qrystartpos
            qry_right_overlen = qryendpos
        if (qry_left_overlen > minqryoverhang and ref_left_overlen > minrefoverhang): 
            chimeric += 1
            ###print 'ref',ref,'reflen',reflen,'refstartpos',refstartpos,'refendpos',refendpos
            ###print 'qry',qry,'qrylen',qrylen,'qrystartpos',qrystartpos,'qryendpos',qryendpos
            ###print 'orientation', orientation, 'confidence', confidence
            ###print 'ref_left_overlen', ref_left_overlen, 'ref_right_overlen', ref_right_overlen
            ###print 'qry_left_overlen', qry_left_overlen, 'qry_right_overlen', qry_right_overlen
            print 'ref',ref,'qry',qry,'orient',orientation,'left overhang at qry pos', qrystartpos
            ###print '--------------------'
        if (qry_right_overlen > minqryoverhang and ref_right_overlen > minrefoverhang):
            chimeric += 1
            ###print 'ref',ref,'reflen',reflen,'refstartpos',refstartpos,'refendpos',refendpos
            ###print 'qry',qry,'qrylen',qrylen,'qrystartpos',qrystartpos,'qryendpos',qryendpos
            ###print 'orientation', orientation, 'confidence', confidence
            ###print 'ref_left_overlen', ref_left_overlen, 'ref_right_overlen', ref_right_overlen
            ###print 'qry_left_overlen', qry_left_overlen, 'qry_right_overlen', qry_right_overlen
            print 'ref',ref,'qry',qry,'orient',orientation,'right overhang at qry pos', qryendpos
            ###print '--------------------'
        s += confidence 
    print '@ Xmap', myfile, 'has', len(table1),'alignments with an average confidence of',s/len(table1)
    print '@ Xmap', myfile, 'detected', chimeric, 'chimeric unitigs'
    #print reflen
    #print refaligned

    # check how many opt map molecules are not aligned
    count_unaligned = 0
    len_unaligned = 0.0
    min_unaligned = 100000000.0
    max_unaligned = 0.0
    for refi in refaligned.keys():
        if (refaligned[refi] == False):
            count_unaligned += 1
	    length = mreflen[refi]
            if (length > max_unaligned):
		max_unaligned = length
            if (length < min_unaligned):
		min_unaligned = length
            len_unaligned += length
    print '*',count_unaligned,'molecules in the ref are not aligned with contigs with min conf',min_confidence,'for a total length of',len_unaligned,'bps'
    if (count_unaligned > 0):
        print '* average len', len_unaligned/count_unaligned, 'min', min_unaligned,'max',max_unaligned

    # check how many contigs are not aligned
    count_unaligned = 0
    len_unaligned = 0.0
    min_unaligned = 100000000.0
    max_unaligned = 0.0
    print 'list unaligned contigs'
    for qryi in qryaligned.keys():
        if (qryaligned[qryi] == False):
            print qryi
            count_unaligned += 1
	    length = mqrylen[qryi]
            if (length > max_unaligned):
		max_unaligned = length
            if (length < min_unaligned):
		min_unaligned = length
            len_unaligned += length
    print '-----------------'
    print '+',count_unaligned,'contigs not aligned with opt map with min conf',min_confidence,'for a total length of',len_unaligned,'bps'
    if (count_unaligned > 0):
        print '+ average len', len_unaligned/count_unaligned, 'min', min_unaligned,'max',max_unaligned

    # process the indel table
    table3 = ReadTable(myfile+'.indel', 10, '\t') # 10 lines of header
    indels = 0
    indelsthreshold = 30.0
    for row in table3:
        confidence = float(row[9])
        if confidence > indelsthreshold:
            print '^ ref',row[2],'qry',row[1],'start',row[4],'end',row[5],'type',row[10],'confidence',confidence
            indels += 1
    print '^ Indel', myfile, 'has', indels, 'compressed repeats with a confidence score of at least',indelsthreshold

    #process the chimeric table
    #table2 = ReadTable(myfile+'.chimeric', 9, '\t') # 10 lines of header
    #threshold = 0.0
    #chimeric = 0
    #for row in table2:
    #    confidence1 = float(row[9])
    #    confidence2 = float(row[15])
    #    if (confidence1>threshold and confidence2>threshold):
    #        #print 'contig',int(row[1]),'matches',int(row[2]),'and',int(row[3]),'with confidence',confidence1,'and',confidence2
    #        chimeric += 1
    #    #else:
    #    #    print 'Discared', confidence1, confidence2
    # print 'Chimeric', myfile, 'has', chimeric, 'chimeric contigs with a confidence score of at least',threshold

if __name__ == '__main__':
    main()
