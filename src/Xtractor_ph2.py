# Merge the filtered rule files and consolidate their counts
# Also compute forward & reverse lexical probs

import os, math
import sys
import heapq


# Global variables
min_lprob = -13.8155        # for log (natural)
max_lprob = -0.000100005
null_term = ''
candLst = []
lexProbDict = {}
ruleDict = {}
tgtDict = {}
phrDict = {}


def loadLexProbDistrib(lexDFile):
    '''Loads the lexical probability distribution in dictionary lexProbDict'''

    print 'Loading lexical probability distribution (both forward and reverse) ...'
    global lexProbDict
    # The src & tgt lexemes and forward & reverse lexical probabilities are stored in dictionary as:
    # key: (src_word, tgt_word); value: (fwrd_lex_prob, rvrs_lex_prob)
    lexF = open(lexDFile, 'r')
    for line in lexF:
        line = line.strip()
        src, tgt, f_prob, r_prob = line.split(' ||| ')
        lexProbDict[(src, tgt)] = (float(f_prob), float(r_prob))
    lexF.close()


def loadFwrdLexProbDistrib(lexDFile):
    '''Loads the forward lexical probability distribution in dictionary lexProbDict'''

    print 'Loading forward lexical probability distribution ...'
    global lexProbDict
    # The src & tgt lexemes and forward & reverse lexical probabilities are stored in dictionary as:
    # key: (src_word, tgt_word); value: (fwrd_lex_prob, rvrs_lex_prob)
    # Only forward lex prob is stored in lexProbDict now and reverse lex prob is updated later by loadRvrsLexProbDistrib()
    lexF = open(lexDFile, 'r')
    line_cnt = 0
    for line in lexF:
        line = line.strip()
        try:
            tgt, src, f_prob = line.split()
            line_cnt += 1
        except ValueError:
            sys.stderr.write( "Encountered error in line # %d :: %s\n" % (line_cnt, line) )
            continue
        lexProbDict[(src, tgt)] = float(f_prob)
    lexF.close()


def loadRvrsLexProbDistrib(lexDFile):
    '''Updates lexProbDict with the reverse lexical probabilities'''

    print 'Loading reverse lexical probability distribution ...'
    global lexProbDict
    # Updates the lexProbDict with correct reverse lex prob
    lexF = open(lexDFile, 'r')
    line_cnt = 0
    for line in lexF:
        line = line.strip()
        try:
            src, tgt, r_prob = line.split()
            line_cnt += 1
        except ValueError:
            sys.stderr.write( "Encountered error in line # %d :: %s\n" % (line_cnt, line) )
            continue
        try:
            f_prob = lexProbDict[(src, tgt)]
            lexProbDict[(src, tgt)] = ( f_prob, float(r_prob) )
        except KeyError:
            print "Key (%s, %s) not found in the lexProbDict." % (src, tgt)
            print "Check if the forward lexical probability has been loaded first. Exiting!!"
            sys.exit(1)
    lexF.close()


def mergeTgtCounts(fileLst, tgtFile):
    '''Read target counts from individual files and merge counts on the fly'''

    global tgtDict
    total_rules = 0
    candLst = []
    fileTrackLst = [ 1 for file in fileLst ]
    stop_iteration = False

    print "Reading target rules and merging their counts ..."
    fHLst = [ open(file, 'r') for file in fileLst ]
    oF = open(tgtFile, 'w')
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

            (tgt, r_cnt) = line.split(' ||| ')
            r_count = float( r_cnt )

            if tgtDict.has_key(tgt):
                for indx1, cand in enumerate( candLst ):
                    if cand[0] == tgt:
                        r_count_new = cand[1] + r_count
                        indxLst = cand[2] + [indx]
                        candLst[indx1] = (tgt, r_count_new, indxLst)
                        break
            else:
                tgtDict[tgt] = 1
                candLst.append( (tgt, r_count, [indx]) )

        if len(candLst) == 0: continue

        heapq.heapify(candLst)
        (tgt, r_count, indxLst) = heapq.heappop(candLst)
        oF.write( "%s ||| %g\n" % (tgt, r_count) )
        total_rules += 1
        tgtDict.pop(tgt)
        for indx1 in indxLst:
            fileTrackLst[indx1] = 1

    for fH in fHLst:
        fH.close()
    oF.close()
    print( "Total # of tgt rules : %d" % (total_rules) )


def read_n_merge(fileLst, outFile1, outFile2):
    '''Read entries from the individual files and merge counts on the fly'''

    global ruleDict
    global candLst
    candLst = []
    total_rules = 0
    fileTrackLst = [ 1 for file in fileLst ]
    stop_iteration = False

    print "Reading rules and merging their counts & alignments ..."
    fHLst = [ open(file, 'r') for file in fileLst ]
    oF1 = open(outFile1, 'w')
    oF2 = open(outFile2, 'w')
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

            (src, tgt, r_count, r_align, f_align) = line.split(' ||| ')
            rule = src + ' ||| ' + tgt
            r_count = float( r_count )

            if ruleDict.has_key(rule):
                mergeValues(rule, r_count, f_align, r_align, [indx], 1)
            else:
                ruleDict[rule] = 1
                mergeValues(rule, r_count, f_align, r_align, [indx], 0)

        if len(candLst) == 0: continue

        heapq.heapify(candLst)
        (rule, r_count, f_align, r_align, indxLst) = heapq.heappop(candLst)
        f_lex_prob = aggregateLexProb(rule, f_align, 1)       # 1 (True) denotes forward alignments
        r_lex_prob = aggregateLexProb(rule, r_align, 0)       # 0 (False) denotes reverse alignments
        #oF1.write( "%s ||| %g ||| %s ||| %s\n" % (rule, r_count, r_align, f_align) )
        oF2.write( "%s ||| %g ||| %g ||| %g\n" % (rule, r_count, r_lex_prob, f_lex_prob) )
        total_rules += 1
        ruleDict.pop(rule)
        for indx1 in indxLst:
            fileTrackLst[indx1] = 1

    for fH in fHLst:
        fH.close()
    oF1.close()
    oF2.close()
    print( "Total # of rules    : %d" % (total_rules) )


def mergeValues(rule, r_count, f_align, r_align, indxLst, rule_exists):
    '''Adds/ merges the values in the heap'''

    global ruleDict
    global candLst

    if not rule_exists:
        candLst.append( (rule, r_count, f_align, r_align, indxLst) )
    else:
        for indx, cand in enumerate( candLst ):
            if cand[0] == rule:
                r_count_new = cand[1] + r_count
                f_align_new = unifyAlignment(cand[2], f_align)
                r_align_new = unifyAlignment(cand[3], r_align)
                indxLst += cand[4]
                candLst[indx] = (rule, r_count_new, f_align_new, r_align_new, indxLst)
                break


def unifyAlignment(align_str, curr_align):
    '''Unifies the current alignment with existing alignments'''

    alignment_found = False
    for prev_align in align_str.split(' ## '):
        if prev_align == curr_align:
            alignment_found = True
            break

    if alignment_found:
        return align_str
    else:
        return align_str + ' ## ' + curr_align


def aggregateLexProb(rule, align_str, align_direction):
    '''Aggregates the lexical probabilities of multiple alignments for a given rule and calculate its average/ worst lexical prob'''

    worst_lex_prob = float('inf')
    #align_count = 0
    #lex_prob_sum = 0.0
    ruleT = rule.split(' ||| ')                   # calcLexProb() expects ruleT
    for align_indx in align_str.split(' ## '):
        lex_prob = calcLexProb(ruleT, align_indx, align_direction)
        if lex_prob < worst_lex_prob:
            worst_lex_prob = lex_prob
        #lex_prob_sum += calcLexProb(ruleT, align_indx, align_direction)
        #align_count += 1

    #return lex_prob_sum / float(align_count)
    return worst_lex_prob


def calcLexProb(ruleT, align_indx, align_direction):
    '''Calculates the lexical probability for the given rule and its alignment'''

    global null_term
    global lexProbDict
    lex_prob = 1.0
    srcWrds = []
    tgtWrds = []
    countLst = []
    probLst = []
    alignPosDict = {}

    srcWrds = ruleT[0].split()
    tgtWrds = ruleT[1].split()

    # Initialize the count and prob lists based on the alignment direction
    if align_direction:
        plen = len(srcWrds)
    else:
        plen = len(tgtWrds)
    for i in range(0, plen):
        countLst.append(0)
        probLst.append(0.0)

    for i in range(0, len(align_indx), 2):
        # Get source and target words from the alignment
        if align_indx[i] == 'Z':
            src_wrd = null_term
        else:
            indx = int( align_indx[i] )
            src_wrd = srcWrds[indx]

        if align_indx[i+1] == 'Z':
            tgt_wrd = null_term
        else:
            indx = int( align_indx[i+1] )
            tgt_wrd = tgtWrds[indx]

        lexT = (src_wrd, tgt_wrd)
        # Count the alignments depending on the alignment direction and
        # incrementally compute the lexical probability
        if align_direction:
            pos2count = int( align_indx[i] )
            prob = lexProbDict[lexT][0]
        else:
            pos2count = int( align_indx[i + 1] )
            prob = lexProbDict[lexT][1]

        if alignPosDict.has_key(pos2count):
            countLst[pos2count] += 1
            probLst[pos2count] += prob
        else:
            alignPosDict[pos2count] = 1
            countLst[pos2count] = 1
            probLst[pos2count] = prob

    for indx, sum_prob in enumerate( probLst ):
        if sum_prob > 0.0:
            align_count = countLst[indx]
            lex_prob *= sum_prob / float( align_count )

    return lex_prob

def mergePhrCounts(fileLst, phrFile):
    '''Read LRM counts from individual phrase files and merge counts on the fly'''

    global phrDict
    phrDict = {}
    total_phrs = 0
    candLst = []
    fileTrackLst = [ 1 for file in fileLst ]
    stop_iteration = False

    print "Reading phraes and merging their LRM counts ..."
    fHLst = [ open(file, 'r') for file in fileLst ]
    oF = open(phrFile, 'w')
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

            (src, tgt, l2r_cnt, r2l_cnt) = line.split(' ||| ')
            l2r = [float(i) for i in l2r_cnt.strip().split()]
            r2l = [float(i) for i in r2l_cnt.strip().split()]
            phr = src + ' ||| ' + tgt
            if phr in phrDict:
                for indx1, cand in enumerate( candLst ):
                    if cand[0] == phr:
                        l2r_new = [cand[1][i] + l2r[i] for i in range(3)]
                        r2l_new = [cand[2][i] + r2l[i] for i in range(3)]
                        indxLst = cand[3] + [indx]
                        candLst[indx1] = (phr, l2r_new, r2l_new, indxLst)
                        break
            else:
                phrDict[phr] = 1
                candLst.append( (phr, l2r, r2l, [indx]) )

        if len(candLst) == 0: continue
        heapq.heapify(candLst)
        (phr, l2r, r2l, indxLst) = heapq.heappop(candLst)
        ## double check counts
        assert(sum(l2r)==sum(r2l)), ('l2r and r2l counts should be the same, phr: %s' %phr)
        ## compute probabilty distribution of orientations based on counts using maximum likelihood principle
        #l2r = computeLRM(l2r)
        #r2l = computeLRM(r2l)
        oF.write( "%s ||| %s ||| %s\n" % (phr, " ".join([str(i) for i in l2r]), " ".join([str(i) for i in r2l])) )
        total_phrs += 1
        phrDict.pop(phr)
        for indx1 in indxLst:
            fileTrackLst[indx1] = 1

    for fH in fHLst:
        fH.close()
    oF.close()
    print( "Total # of phrases : %d" % (total_phrs) )

def computeLRM(countLst):
    '''Given the counts, it computes probabilty distribution of orientations (log probability)'''
    tot = sum(countLst)
    probLst = []
    for c in countLst:
        if c == 0:        # use min_lprob if c is zero
            p = min_lprob
        else:
            if c == tot: p = max_lprob    # use max_lprob if prob = 1 
            else: p = math.log( float(c)/ tot )
        probLst.append(p)
    return probLst

def mergeTgtPhrCounts(fileLst, tgtPhrFile):
    '''Read target_phrase counts from individual files and merge counts on the fly'''

    global tgtPhrDict
    tgtPhrDict = {}
    total_rules = 0
    candLst = []
    fileTrackLst = [ 1 for file in fileLst ]
    stop_iteration = False

    print "Reading target phrases and merging their counts ..."
    fHLst = [ open(file, 'r') for file in fileLst ]
    oF = open(tgtPhrFile, 'w')
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

            (tgt, l2r, r2l) = line.split(' ||| ')
            l2r = map(lambda x: float(x), l2r.split())
            r2l = map(lambda x: float(x), r2l.split())

            if tgtPhrDict.has_key(tgt):
                for indx1, cand in enumerate( candLst ):
                    if cand[0] == tgt:
                        l2r_new = [cand[1][i] + l2r[i] for i in range(3)]
                        r2l_new = [cand[2][i] + r2l[i] for i in range(3)]
                        indxLst = cand[-1] + [indx]
                        candLst[indx1] = (tgt, l2r_new, r2l_new, indxLst)
                        break
            else:
                tgtPhrDict[tgt] = 1
                candLst.append( (tgt, l2r, r2l, [indx]) )

        if len(candLst) == 0: continue

        heapq.heapify(candLst)
        (tgt, l2r, r2l, indxLst) = heapq.heappop(candLst)
        l2r = ' '.join( map(lambda x: str(x), l2r))
        r2l = ' '.join( map(lambda x: str(x), r2l))
        oF.write( "%s ||| %s ||| %s\n" % (tgt, l2r, r2l) )
        total_rules += 1
        tgtPhrDict.pop(tgt)
        for indx1 in indxLst:
            fileTrackLst[indx1] = 1

    for fH in fHLst:
        fH.close()
    oF.close()
    print( "Total # of tgt phrases : %d" % (total_rules) )

def mergeCounts(fileLst, cntFile):
    '''Read individual count files and get summation'''

    print "Reading counts ..."
    fHLst = [ open(file, 'r') for file in fileLst ]
    oF = open(cntFile, 'w')
    L2R = [0, 0, 0]
    R2L = [0, 0, 0]
    totCnt = 0
    for inF in fHLst:
        for line in inF:
            line = line.strip()
            (c, l2r, r2l) = line.split(" ||| ")
            totCnt += float(c)
            for i, c in enumerate(l2r.split()):
               L2R[i] += float(c)
            for i, c in enumerate(r2l.split()):
               R2L[i] += float(c)
    totCnt = sum(R2L)
    #assert (sum(L2R) == totCnt and sum(R2L) == totCnt), ("phrase counts are not compatible, %d != %d" %(totCnt, sum(L2R)))

    
    l2r = ' '.join( map(lambda x: str(x), L2R))
    r2l = ' '.join( map(lambda x: str(x), R2L))
    oF.write( "%s ||| %s ||| %s\n" % (str(totCnt), l2r, r2l) )

def main():

    global null_term

    inDir = sys.argv[1]
    outDir = sys.argv[2]
    if len(sys.argv) == 4:
        lexDFile = sys.argv[3]
        null_term = '-NULL-'
    elif len(sys.argv) == 5:
        fwrdLexDist = sys.argv[3]
        rvrsLexDist = sys.argv[4]
        null_term = 'NULL'

    if len(sys.argv) == 4:
        loadLexProbDistrib(lexDFile)
    elif len(sys.argv) == 5:
        pass
        loadFwrdLexProbDistrib(fwrdLexDist)
        loadRvrsLexProbDistrib(rvrsLexDist)

    if not inDir.endswith('/'): inDir += '/'
    if not outDir.endswith('/'): outDir += '/'
    fileLst = []
    tgtFileLst = []
    cntFileLst = []
    tgtPhrFileLst = []
    phrFileLst = []
    lexicalizedReorderingFlag = False
    for file in os.listdir(inDir):
        pathFile = inDir + file
        if os.path.isfile(pathFile):
            if file.startswith("tgt"):        tgtFileLst.append(pathFile)
            elif file.startswith("cnt"):      cntFileLst.append(pathFile)
            elif file.startswith("phr.tgt"):  tgtPhrFileLst.append(pathFile)
            elif file.startswith("phr."):     
                phrFileLst.append(pathFile)
                lexicalizedReorderingFlag = True
            else: fileLst.append(pathFile)
    
    outFile1 = outDir + 'rules_cnt_align.out'
    outFile2 = outDir + 'rules_cnt_lprob.out'
    tgtFile = outDir + 'tgt_rules.all.out'
    if lexicalizedReorderingFlag:
        mergeCounts(cntFileLst, outDir + 'cnt.all.out') 
        tgtPhrFile = outDir + 'tgt_phr.all.out'
        mergeTgtPhrCounts(tgtPhrFileLst, tgtPhrFile)
        phrFile = outDir + 'phr_cnt.all.out'
        mergePhrCounts(phrFileLst, phrFile)

    mergeTgtCounts(tgtFileLst, tgtFile)

    read_n_merge(fileLst, outFile1, outFile2)


if __name__ == '__main__':
    main()
