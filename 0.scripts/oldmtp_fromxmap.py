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
            print allc
            # find the leftmost contig
            candidates = []
	    for x in group:
                if (x.confidence > min_confidence) and (x.chimeric == False) and (x.contained == False):
                    #print '*processing contig', x.qry, 'start',x.start,'end',x.end
                    if (x.start == allc[0]): # this is the first point 
                        candidates.append(x)
                        print '*appending',x.qry
            # among all contigs that start at allc[0], find the longest one
            rightmost_endcoordinate = candidates[0].end
            rightmost = candidates[0]
            for x in candidates:
                if (x.end > rightmost_endcoordinate):
                    rightmost = x
                    rightmost_endcoordinate = x.end
            newmtp = [rightmost]
            print 'first contig',newmtp[0].qry,'rightmost_coordinate',rightmost_endcoordinate
            while (rightmost_endcoordinate < allc[-1]): # stop when we covered the entire interval
                # among all intervals intersecting with newmtp[-1], find the one with the rightmost endpoint
                candidates = []
                for x in group:
                    if (x.confidence > min_confidence) and (x.chimeric == False) and (x.contained == False):
                        print 'processing contig', x.qry, 'start',x.start,'end',x.end
                        if (x.start <= newmtp[-1].end): # endpoint of the last element added to newmtp
                            candidates.append(x)
                rightmost_endcoordinate = candidates[0].end
                rightmost = candidates[0]
                for x in candidates:
                    if (x.end > rightmost_endcoordinate):
                        rightmost = x
                        rightmost_endcoordinate = x.end
                newmtp.append(rightmost)
                print 'adding contig', newmtp[-1].qry,'rightmost_coordinate',rightmost_endcoordinate

            # computes coverage
            coverage = [0] * len(allc) # coverage (int) for each interval i
            who = [[] for y in xrange(len(allc))] # list of contigs overlapping interval i
	    for x in group:
                if (x.confidence > min_confidence) and (x.chimeric == False) and (x.contained == False):
                    print 'processing contig', x.qry
		    try:
		        pos_start = allc.index(x.start)
                    except ValueError:
                        print 'location',x.start,'not found in',allc
			sys.exit()
		    try:
			pos_end = allc.index(x.end)
                    except ValueError:
                        print 'location',x.end,'not found in',allc
                        sys.exit()
		    for j in range(pos_start, pos_end):
			coverage[j] += 1
			who[j].append(x)
                #print coverage, who
            # compute the binary coverage (covered - True or not - False)
            binary_coverage = [False] * len(allc)
            for i in range(len(allc)):
                if (coverage[i] > 0):
                    binary_coverage[i] = True
	    for i in range(len(allc)-1):
		print 'from',allc[i],'to',allc[i+1],'coverage',coverage[i],'binary',binary_coverage[i],'by ['
                for x in who[i]:
                    print x.qry
                print ']'
            # the initial MTP has obligatory contigs (which are the only one covering an interval)
            mtp = []
            for i in range(len(allc)):
                if (coverage[i] == 1):
                    mtp += who[i]
            # print 'mtp',mtp
            # does the MTP cover the same coverage of the binary coverage?
            while (True):
                mtp_coverage = [False] * len(allc)
                for x in mtp:
                    print 'processing MTP contig', x.qry
                    try:
                        pos_start = allc.index(x.start)
                    except ValueError:
                        print 'location',x.start,'not found in',allc
                        sys.exit()
                    try:
                        pos_end = allc.index(x.end)
                    except ValueError:
                        print 'location',x.end,'not found in',allc
                        sys.exit()
                    for j in range(pos_start, pos_end):
                        mtp_coverage[j] = True

                print 'binary coverage', binary_coverage
                print 'mtp coverage', mtp_coverage
                # print 'mtp',mtp
                if (mtp_coverage == binary_coverage):
                    print 'mtp is optimal'
                    break
                # add the first contig that covers the difference
                for i in range(len(allc)):
                    if (binary_coverage[i]==True and mtp_coverage[i]==False):
                        #print 'who[i]',who[i]
                        mtp.append(who[i][0])

            newmtp_coverage = [False] * len(allc)
            for x in newmtp:
                print 'processing MTP contig', x.qry
                try:
                    pos_start = allc.index(x.start)
                except ValueError:
                    print 'location',x.start,'not found in',allc
                    sys.exit()
                try:
                    pos_end = allc.index(x.end)
                except ValueError:
                    print 'location',x.end,'not found in',allc
                    sys.exit()
                for j in range(pos_start, pos_end):
                    newmtp_coverage[j] = True
            if (newmtp_coverage == binary_coverage):
                print 'newmtp is optimal'
            print '---------------ENDOFGROUP-------------------'
            # save outalignments
	    for x in mtp:
                outalignments.append(x)
            firstrow = nextrow

    # save the MTP in a new xmap file
    with open(myfile+'_mtp.xmap', 'wb') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter='\t')
        for x in header: 
            csvwriter.writerow(x)
        i = 1
        for x in outalignments:
            csvwriter.writerow([i]+x.unpack()) 
            i += 1

if __name__ == '__main__':
    main()
