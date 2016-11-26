#!/usr/bin/python
import csv
import sys

class Alignment:
    def __init__(self, qry, ref, qrystartpos, qryendpos, refstartpos, refendpos, orientation, confidence, hitenum, qrylen, reflen, channel, alignmentstring): 
        self.qry = qry
        self.ref = ref
        self.qrystartpos = qrystartpos  
        self.qryendpos = qryendpos
        self.refstartpos = refstartpos
        self.refendpos = refendpos
        self.orientation = orientation
        self.confidence = confidence
        self.qrylen = qrylen
        self.reflen = reflen
        self.hitenum = hitenum
        self.channel = channel
        self.alignmentstring = alignmentstring
        self.start=0 # will be filled later
        self.end=0   # will be filled later
        self.qry_left_overlen=0 # will be filled later
        self.qrt_right_overlen=0 # will be filled later
        self.chimeric = False # will be change later
        self.contained = False # will be change later

    def __str__(self):
        return 'refstart %.1f refend %.1f qry %d qrylen %.1f qrystart %.1f qryend %.1f %s conf %.2f start %.2f end %2.f qry_left_overlen %2.f qrt_right_overlen %2.f chimeric %s contained %s' % (self.refstartpos, self.refendpos, self.qry, self.qrylen, self.qrystartpos, self.qryendpos, self.orientation, self.confidence, self.start, self.end, self.qry_left_overlen, self.qrt_right_overlen, self.chimeric, self.contained)
        #return 'ref %d, reflen %.1f, refstartpos %.1f, refendpos %.1f, qry %d, qrylen %.1f, qrystartpos %.1f, qryendpos %.1f, orientation %s, confidence %.2f' % (self.ref, self.reflen, self.refstartpos, self.refendpos, self.qry, self.qrylen, self.qrystartpos, self.qryendpos, self.orientation, self.confidence)

    def unpack(self):
        return [self.qry, self.ref, self.qrystartpos, self.qryendpos, self.refstartpos, self.refendpos, self.orientation, self.confidence, self.hitenum, self.qrylen, self.reflen, self.channel, self.alignmentstring]


def main():
    """
     reads a xmap file and produce another xmap file with the MTP
    """
    if len(sys.argv) == 2:
        myfile = sys.argv[1]
    else:
        print "Usage: python scorexmap.py <file>, where <file> is the prefix, scorexmap will open <file>.xmap and produce <file>_mtp.xmap"
        return 0
 
    # discard alignments below min_confidence
    min_confidence = 20
    header_lines = 10
    header = []
    outalignments = []
    # alignment overhangs above this number of bps are considered chimeric
    minrefoverhang = 100000
    minqryoverhang = 100000

    with open(myfile+'.xmap', 'rb') as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        for i in range(header_lines): # 10 lines of header
            header.append(csvreader.next()) # save them
        # read the first non-header line
        try:
            firstrow = csvreader.next()
        except StopIteration:
            sys.exit()
        while True:
            # save the params of the alignment in group array
            first = Alignment(int(firstrow[1]),int(firstrow[2]),float(firstrow[3]),float(firstrow[4]),float(firstrow[5]),
                              float(firstrow[6]),firstrow[7],float(firstrow[8]),firstrow[9],float(firstrow[10]),
                              float(firstrow[11]),int(firstrow[12]),firstrow[13])
            group = [first]
            # get the next row
            try:
                nextrow = csvreader.next()
            except StopIteration:
                break
            # collect all the alignmenth with the same ref number
            while (first.ref == int(nextrow[2])): # int(nextrow[2]) is ref id
                group.append(Alignment(int(nextrow[1]),int(nextrow[2]),float(nextrow[3]),float(nextrow[4]),
                                       float(nextrow[5]),float(nextrow[6]),nextrow[7],float(nextrow[8]),nextrow[9],
                                       float(nextrow[10]),float(nextrow[11]),int(nextrow[12]),nextrow[13]))
                try:
                    nextrow = csvreader.next()
                except StopIteration:
                    break
            # process the group of all contigs aligned to the same opt map molecule
	    print 'processing contigs aligned optical molecule',group[0].ref
            print 'optical molecule is',group[0].reflen,'bases long'
	    max_qry_left_overlen = 0
	    max_qry_right_overlen = 0
	    # find the reference-based coordinates for each contig
            for x in group:
		if (x.orientation == '+'):
                    x.qry_left_overlen = x.qrystartpos
                    x.qry_right_overlen = x.qrylen - x.qryendpos
		else:
                    x.qry_left_overlen = x.qrylen - x.qrystartpos
                    x.qry_right_overlen = x.qryendpos
                x.start = x.refstartpos - x.qry_left_overlen
                x.end = x.refendpos + x.qry_right_overlen 
            # check for chimeric contigs
            for x in group:
                if (x.confidence > min_confidence):
                    ref_left_overlen = x.refstartpos
                    ref_right_overlen = x.reflen - x.refendpos
                    if (x.qry_left_overlen > minqryoverhang and ref_left_overlen > minrefoverhang): 
                        x.chimeric = True
                    if (x.qry_right_overlen > minqryoverhang and ref_right_overlen > minrefoverhang):
                        x.chimeric = True
            for x in group:
                if (x.chimeric):
                    print 'contig', x.qry, 'is chimeric'
            # check for containments
      	    for i in range(len(group)):
                for j in range(i+1,len(group)):
                    #print 'comparing',group[i],'and',group[j]
                    if (group[i].start > group[j].start) and (group[i].end < group[j].end):
                        #print group[i].qry,'is contained in',group[j].qry
                        group[i].contained = True
                    if (group[j].start > group[i].start) and (group[j].end < group[i].end):
                        #print group[j].qry,'is contained in',group[i].qry
                        group[j].contained = True
            for x in group:
                if (x.contained):
                    print 'contig',x.qry,'is contained'
            # collect all coordinates
            starts,ends = [],[]
            for x in group:
                if (x.confidence > min_confidence) and (x.chimeric == False) and (x.contained == False):
                    starts.append(x.start)
                    ends.append(x.end)
            allc = sorted(starts+ends)
            print 'allcoordinates',allc
            if (allc == []):
                print 'group is EMPTY!'
                print '---------------ENDOFGROUP-------------------'
                firstrow = nextrow
                continue
            # computes coverage
            coverage = [0] * len(allc) # coverage (int) for each interval i
            who = [[] for y in xrange(len(allc))] # list of contigs overlapping interval i
	    for x in group:
                if (x.confidence > min_confidence) and (x.chimeric == False) and (x.contained == False):
		    pos_start = allc.index(x.start)
                    pos_end = allc.index(x.end)
		    for j in range(pos_start, pos_end):
			coverage[j] += 1
			who[j].append(x)
                #print coverage, who
            # compute the binary coverage (covered->True or not->False)
            binary_coverage = [False] * len(allc)
            for i in range(len(allc)):
                if (coverage[i] > 0):
                    binary_coverage[i] = True
	    for i in range(len(allc)-1):
		print 'from',allc[i],'to',allc[i+1],'coverage',coverage[i],'binary',binary_coverage[i],'by [',
                for x in who[i]:
                    print x.qry,
                print ']' 
            # if there are intervals not covered by anything, in that case add fake ones so we can get a solid
	    for i in range(len(allc)-1):
                if (coverage[i] == 0):
                    # dummy contigs have ID 0
                    newx = Alignment(0, 0, 0, 0, 0, 0, '', 2*min_confidence, '', 0, 0, 0, '')
                    newx.start = allc[i]
                    newx.end = allc[i+1]
                    newx.qry_left_overlen = 0 
                    newx.qrt_right_overlen = 0 
                    newx.chimeric = False 
                    newx.contained = False
                    print '+ adding dummy contig ',newx
                    group.append(newx)
            # compute the MTP (greedy)
            mtp = []
            current = allc[0]
            while (current < allc[-1]): # stop when we covered the entire interval
                # among all intervals before current, find the one with the rightmost endpoint
                candidates = []
                for x in group:
                    if (x.confidence > min_confidence) and (x.chimeric == False) and (x.contained == False):
                        # print 'processing contig', x.qry, 'start',x.start,'end',x.end
                        if (x.start <= current): # endpoint of the last element added to mtp
                            candidates.append(x)
                current = candidates[0].end
                rightmost = candidates[0]
                for x in candidates:
                    if (x.end > current):
                        rightmost = x
                        current = x.end
                mtp.append(rightmost)
                print '+ adding contig', rightmost.qry,'coordinate',current
            # compute the coverage from the MTP
            mtp_coverage = [False] * len(allc)
            for x in mtp:
                #print 'processing MTP contig', x.qry
                pos_start = allc.index(x.start)
                pos_end = allc.index(x.end)
                for j in range(pos_start, pos_end):
                    mtp_coverage[j] = True
            # compare whether the coverage of the MTP is the same of all intervals
            if (mtp_coverage == binary_coverage):
                print 'mtp is optimal'
            # save outalignments
	    for x in mtp:
                if (x.qry == 0):
                    print 'contig',x,'is skipped'
                else:
                    outalignments.append(x)
            print '---------------ENDOFGROUP-------------------'
            firstrow = nextrow

    # save the MTP in a new xmap file and count the number of unitigs in each assembly
    # assembly_names =   [(1,878,'low_default'),
    #                    (879,1824,'low_outcov100'),
    #                    (1825,2877,'high_default'),
    #                    (2878,3970,'high_erate.15_outcov100'),
    #                    (3971,5089,'normal_erate.15_outcov100')]
    assembly_names =   [(1,498,'abruijn'),
                        (499,1376,'low_default1'),
                        (1377,3165, 'falcon'),
                        (3166,4043,'low_default2'),
                        (4044,4989,'low_outcov100'),
                        (4990,6042,'high_default'),
                        (6043,7135,'high_erate.15_outcov100'),
                        (7136,8254,'normal_erate.15_outcov100')]
    counts = {}
    for y in assembly_names:
	counts[y[2]] = 0.0
    total_qry = 0.0
    with open(myfile+'_list.txt', 'wb') as listfile:
        with open(myfile+'_mtp.xmap', 'wb') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter='\t')
            # copies the old xmap header
            for x in header: 
                csvwriter.writerow(x)
            i = 1 # progressive number
            for x in outalignments:
                for y in assembly_names:
                    if (x.qry >= y[0] and x.qry <= y[1]):
                        counts[y[2]] += 1.0
                        total_qry += 1.0
                        break
                listfile.write(str(x.qry)+'\n')
                csvwriter.writerow([i]+x.unpack()) 
                i += 1
    csvfile.close()
    listfile.close()
    # prints statistics
    for y in assembly_names:
        name = y[2]
	print name,'assembly',counts[name],'fraction',counts[name]/total_qry

if __name__ == '__main__':
    main()
