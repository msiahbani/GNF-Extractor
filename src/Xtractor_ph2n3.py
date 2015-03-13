# Identify the source phrases in the given dev/test set and filter the rule file for the source phrases
# Get the total counts for the target phrases that co-occur with the filtered source phrases and write them in a temp file
# Filter the phrase file which keeps the lexicalized reordering model (if phrase files exist)


import optparse, sys, os, logging, time 
import heapq
import math
from datetime import timedelta
from multiNonTermTrie import SimpleSuffixTree

new_trie = None
srcRuleDict = {}
tgtCntDict = {}
phrDict = {}
tgtPhrDict = {}


def loadTrie(ruleFile):
    '''Loads the Suffix Trie with rules from ruleFile'''

    global new_trie, opts
    prev_src = ''
    iF = open(ruleFile, 'r')

    print 'Loading phrases from data file : %s ...\n' % ruleFile
    for line in iF:
        line = line.strip()

        (src, _) = line.split(' ||| ', 1)
        if prev_src != src:
            prev_src = src
            if new_trie is None:
                new_trie = SimpleSuffixTree(src, opts.tot_src_terms)
            else:
                new_trie.addText(src)

    iF.close()
#    new_trie.printFullTree()
    return None

def convertRule2Phr(rule):
    '''Convert the rule to phrase format (remove the boundary non-terminals and replace the others with a uniqe token)'''

    srcLst = []
    src = rule.split()
    l = len(src)
    for i,s_tok in enumerate(src):
        if s_tok.startswith("X__"):
            if i!=0 and i!=l-1: srcLst.append("NON_TOK")
        else: srcLst.append(s_tok)
    l = len(srcLst)
    i = l-1
    while i > -1:
        if srcLst[i] == "NON_TOK": l -= 1
        else: break
        i -= 1
    return " ".join(srcLst[:l])

def filterRules(dataFile):
    '''Filter the partial rule file for the specified dev/test set before computing p(f|e) and p(e|f)'''

    global new_trie, opts
    global srcRuleDict, phrDict

    rF = open(dataFile, 'r')
    print 'Filtering rules for file : %s ...\n' % dataFile
    try:
        for line in rF:
            line = line.strip()
            words = []
            words = line.split()

            for i in range ( len(words) ):
                for j in range(i, i + opts.max_phr_len):
                    if j >= len(words): break
                    phr = ' '.join( words[i:j+1] )
                    rulesLst = []

                    # @type new_trie SimpleSuffixTree
                    rulesLst = new_trie.matchPattern(phr)
                    for rule in rulesLst:
                        srcRuleDict[rule[0]] = 1
                        if opts.lex_reorder_model:
                            phr = convertRule2Phr(rule[0])
                            phrDict[phr] = 1

    finally:
        rF.close()


def writeRules(ruleFile, tempOutFile):

    global srcRuleDict
    global tgtCntDict
    tgt_rule_cnt = 0
    tgtCntDict = {}

    rF = open(ruleFile, 'r')
    tF = open(tempOutFile, 'w')
    print 'Filtering rules from file : %s ...\n' % ruleFile
    try:
        for line in rF:
            line = line.strip()

            (src, tgt, _) = line.split(' ||| ', 2)
            if srcRuleDict.has_key(src):
                tF.write( "%s\n" % (line) )
                if not tgtCntDict.has_key(tgt):
                    tgt_rule_cnt += 1
                    tgtCntDict[tgt] = 0.0

    finally:
        rF.close()
        tF.close()
        print "Unique # of tgt_rules found in set :%d" % ( tgt_rule_cnt )


def updateTgtCnts(tgtFile, tempTgtFile):

    global tgtCntDict, tgtPhrDict

    rF = open(tgtFile, 'r')
    print 'Updating target counts from file : %s ...\n' % tgtFile
    try:
        for line in rF:
            line = line.strip()

            (tgt, r_cnt) = line.split(' ||| ')
            if tgtCntDict.has_key(tgt):
                tgtCntDict[tgt] += float( r_cnt )
                tgtPhr = convertRule2Phr(tgt)
                tgtPhrDict[tgtPhr] = 1

    finally:
        rF.close()

    print 'Writing target counts to file : %s ...\n' % tempTgtFile
    gF = open(tempTgtFile, 'w')
    tgtRules = []
    tgtRules = tgtCntDict.keys()
    tgtRules.sort()
    for tgt in tgtRules:
        gF.write( "%s ||| %g\n" % (tgt, tgtCntDict[tgt]) )
    gF.close()

def writePhrLRM(phrFile, tempPhrFile, tempCntFile):
    '''write the filtered phrases'''

    global phrDict 
    totOrientation = [0 for i in range(6)]
    phrFile = open(phrFile, "r")
    outFile = open(tempPhrFile, 'w')
    cntFile = open(tempCntFile, 'w')
    
    try:
        for line in phrFile:
            line = line.strip()
            (src, tgt, l2r, r2l) = line.split(' ||| ')
            if src in phrDict:
                outFile.write( "%s\n" % (line) )
            for i, c in enumerate(l2r.split()):
               totOrientation[i] += int(c)
            for i, c in enumerate(r2l.split()):
               totOrientation[i+3] += int(c)
        totPhrCnt = sum(totOrientation[0:3])
        l2r = ' '.join( map(lambda x: str(x), totOrientation[:3]))
        r2l = ' '.join( map(lambda x: str(x), totOrientation[3:]))
        cntFile.write( "%s ||| %s ||| %s\n" % (str(totPhrCnt), l2r, r2l) )
    finally:
        phrFile.close()
        outFile.close()
        cntFile.close()

def writeTgtPhrLRM(tgtPhrFile, tempTgtPhrFile):
    '''write the filtered tgt phrases'''

    global tgtPhrDict
    tgtFile = open(tgtPhrFile, "r")
    outFile = open(tempTgtPhrFile, 'w')
    
    try:
        for line in tgtFile:
            line = line.strip()
            (tgt, _) = line.split(' ||| ', 1)
            if tgt in tgtPhrDict:
                outFile.write( "%s\n" % (line) )
    finally:
        tgtFile.close()
        outFile.close()

def main():
    global opts

    if not opts.ruledir.endswith('/'): opts.ruledir += '/'
    if not opts.outdir.endswith('/'): opts.outdir += '/'

    ruleFile = opts.ruledir + str(opts.file_prefix) + ".out"
    tgtFile = opts.ruledir + "tgt." + str(opts.file_prefix) + ".out"
    phrFile = opts.ruledir + "phr." + str(opts.file_prefix) + ".out"
    tgtPhrFile = opts.ruledir + "phr.tgt." + str(opts.file_prefix) + ".out"

    tempOutFile = opts.outdir + str(opts.file_prefix) + ".out"
    tempTgtFile = opts.outdir + "tgt." + str(opts.file_prefix) + ".out"
    if os.path.isfile(phrFile):
        opts.lex_reorder_model = True
        tempPhrFile = opts.outdir + "phr." + str(opts.file_prefix) + ".out"
        tempTgtPhrFile = opts.outdir + "phr.tgt." + str(opts.file_prefix) + ".out"
        tempCntFile = opts.outdir + "cnt." + str(opts.file_prefix) + ".out"
    print "Using the maximum phrase lenght            :", opts.max_phr_len
    print "Using the source side total terms to be    :", opts.tot_src_terms
    print "Filtering lexicalized reordering model     :", opts.lex_reorder_model

    # load the development set phrases to a suffix trie
    loadTrie(ruleFile)

    # filter the consolidated rules that match phrases in trie into tempFile
    filterRules(opts.devFile)

    # write the filtered rules (identified earlier by filterRules)
    writeRules(ruleFile, tempOutFile)

    # update tgtCntDict with counts and write them to tempTgtFile
    updateTgtCnts(tgtFile, tempTgtFile)

    # write the filtered phrases (identified earlier by filterRules)
    if opts.lex_reorder_model:
        writePhrLRM(phrFile, tempPhrFile, tempCntFile)
        writeTgtPhrLRM(tgtPhrFile, tempTgtPhrFile)

if __name__ == "__main__":
    global opts
    optparser = optparse.OptionParser()
    optparser.add_option("-d", "--devFile", dest="devFile", default="dev", help="data (dev/test) file (default=dev)")
    optparser.add_option("-r", "--ruledir", dest="ruledir", default="rules", help="rule directory (default=rules)")
    optparser.add_option("-o", "--outdir", dest="outdir", default="temp", help="temporary filtered directory (default=temp)")
    optparser.add_option("-p", "--prefix", dest="file_prefix", default="1", help="prefix of parallel data files (default=1)")
    optparser.add_option("-l", "--logfile", dest="log_file", default=None, help="filename for logging output")
    #optparser.add_option("","--lexReorderingModel", dest="lex_reorder_model", default=False, action="store_true", help="compute lexicalized reordering model for phrases (default=False)")
    optparser.add_option("", "--maxPhrLen", dest="max_phr_len", default=10, type="int", help="maximum initial phrase length (default=10)")
    optparser.add_option("", "--totSrcTrms", dest="tot_src_terms", default=7, type="int", help="maximum number of terms in src side of rules (default=7)")
    (opts, _) = optparser.parse_args()
    opts.lex_reorder_model = False

    if opts.log_file:
        logging.basicConfig(filename=opts.log_file, filemode='w', level=logging.INFO)

    main()

