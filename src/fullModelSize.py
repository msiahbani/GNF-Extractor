## The SCFG rules of Hiero systems are typically learnt for several smaller sets of initial phrase pairs. ##
## These rules are then filtered for the tuning/test set, whose features (probabilities) are computed. ##
## Thus the full (without filtering) model size is neven known. ##
## This program finds the full model size without building it. ##

import os
import sys
import time
import heapq

# Global variables
candLst = []
ruleDict = {}


def read_n_merge(fileLst):
    '''Read entries from the individual files and them counts on the fly'''

    global ruleDict
    global candLst
    candLst = []
    total_rules = 0
    fileTrackLst = [ 1 for file in fileLst ]
    stop_iteration = False

    print "Reading rules and finding the total model size ..."
    fHLst = [ open(file, 'r') for file in fileLst ]
    
    while True:
        stop_iteration = True
        for f_track in fileTrackLst:
            if f_track != 9:
                stop_iteration = False
                break

        if stop_iteration:
            break

        for indx, f_track in enumerate( fileTrackLst ):
            if f_track == 0 or f_track == 9:
                continue

            fileTrackLst[indx] = 0
            line = fHLst[indx].readline()
            line = line.strip()
            if line == '':
                fileTrackLst[indx] = 9          # Set 9 if 'EOF' is reached
                continue
	    
            (src, tgt, _) = line.split(' ||| ', 2)
            rule = src + ' ||| ' + tgt

            if ruleDict.has_key(rule):
                ruleDict[rule] += [indx]
            else:
                ruleDict[rule] = [indx]
                heapq.heappush(candLst, rule)

        if len(candLst) == 0: continue

        rule = heapq.heappop(candLst)
        total_rules += 1
        for indx1 in ruleDict.pop(rule):
            fileTrackLst[indx1] = 1

    for fH in fHLst:
        fH.close()
    sys.stdout.write( "Size of the full model   : %d\n" % (total_rules) )


def main():

    inDir = sys.argv[1]
    if not inDir.endswith('/'): inDir += '/'
    fileLst = []
    for file in os.listdir(inDir):
        if file.startswith("tgt"): is_tgt_file = True
        else: is_tgt_file = False

        file = inDir + file
        if os.path.isfile(file):
            if not is_tgt_file: fileLst.append(file)

    sys.stderr.write( "Total grammar files found: %d\n" % (len(fileLst)) )
    t_beg = time.time()
    read_n_merge(fileLst)
    sys.stderr.write( "Total time taken         : %g\n" % (time.time() - t_beg) )


if __name__ == '__main__':
    main()
