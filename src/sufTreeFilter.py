#! /usr/bin/python

# Filters the rules and phrases for smaller subsets based on the suffix-tree implementation
# This is an improvement over older way as in cache_rules.py

# Revision history

# 25 Feb 2015 Modifies by Maryam Siahbani to filter lexicalized reordering model for LR-Hiero
# Original version by Baskaran Sankaran

import optparse, sys, os, logging, time, math

sys.path.append( os.path.dirname(sys.path[0]) )
from multiNonTermTrie import SimpleSuffixTree

# Global variables
max_lprob = -0.000100005                             # max log prob (for prob = 1)
unk_lprob = -13.8155                                 # unknown log prob (also for zero probs)
phrDict = {}
lrmDict = {}
ruleDict = {}
new_trie = None

def loadTrie(ruleFile):

    global new_trie, opts
    prev_src = ''

    rF = open(ruleFile, 'r')
    print "Loading rules from file %s into Trie" % (ruleFile)
    for line in rF:
        line = line.strip()
        itemsLst = []
        itemsLst = line.split(' ||| ')
        src = itemsLst[0]
        if prev_src != src:
            if new_trie is None:
                new_trie = SimpleSuffixTree(src, opts.tot_src_terms)
            else:
                new_trie.addText(src)
        prev_src = src

    rF.close()
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
    return " ".join(srcLst)

def processPhrases(inFile):
    '''Reads the text file and processes the phrases'''

    global new_trie, opts
    global phrDict
    global ruleDict, lrmDict
    phrDict = {}
    lrmDict = {}

    t_tot = 0.0
    p_tot = 0
    iF = open(inFile, 'r')
    print 'Processing phrases from file : %s' % inFile
    try:
        for line in iF:
            line = line.strip()

            wLst = line.split()
            for i in range( len(wLst) ):
                for j in range(i, i + opts.max_phr_len):
                    if j >= len(wLst): break
                    src = ' '.join(wLst[i:j+1])
                    if phrDict.has_key(src): continue
                    else: phrDict[src] = 1
                    matchLst = []
                    t_beg = time.time()
                    matchLst = new_trie.matchPattern(src)
                    t_end = time.time()
                    t_tot += t_end - t_beg
                    p_tot += 1
                    for match in matchLst:
                        ruleDict[match[0]] = 1
                        if opts.lex_reorder_model:
                            phr = convertRule2Phr(match[0])
                            lrmDict[phr] = 1

    finally:
        iF.close()
        print "Unique phrases processed : %4d" % (p_tot)
        print "Total time taken         : %f" % (t_tot)
        print "Average time taken       : %f" % (t_tot/ p_tot)

    return None

def writePhrLRM(phrFile, filteredPhrFile):
    '''write the filtered phrases'''

    global lrmDict
    phrFile = open(phrFile, "r")
    outFile = open(filteredPhrFile, 'w')
    try:
        for line in phrFile:
            line = line.strip()
            (src, tgt, _) = line.split(' ||| ', 2)
            if src in lrmDict:
                outFile.write( "%s\n" % (line) )
    finally:
        phrFile.close()
        outFile.close()

def writeFilteredRules(ruleFile, filtFile):

    global ruleDict

    rF = open(ruleFile, 'r')
    oF = open(filtFile, 'w')
    print "Filtering rules from file %s into %s" % (ruleFile, filtFile)
    try:
        for line in rF:
            line = line.strip()
            (src, _) = line.split(' ||| ', 1)
            if ruleDict.has_key(src):
#                feat_str = getLogProb(prob_str, tgt)
#                oF.write("%s ||| %s ||| %s\n" % (src, tgt, feat_str))
#                oF.write("%s ||| %s ||| %s\n" % (src, tgt, prob_str))
                oF.write("%s\n" % (line))

    finally:
        rF.close()
        oF.close()

def getLogProb(prob_str, tgt_rule):
    '''Get the log probability from a string of probs, terminal counts and heuristic prob'''

    global featVec
    global max_lprob
    global unk_lprob
    featVec = []

    for prob in prob_str.split(' '):
        prob = float(prob)                               # Type cast to float
        if prob < 0.0: l_prob = prob
        elif prob == 0.0: l_prob = unk_lprob
        elif prob == 1.0: l_prob = max_lprob
        else:
            l_prob = math.log( prob )                    # get negative log-prob
        featVec.append( l_prob )

    # this part is now commented and the idea is to compute them in run time (at least temporarily)
    # now add the phrase and word penalties to the featVec
#    term_count = 0
#    for tgt_term in tgt_rule.split(' '):
#        if not tgt_term.startswith('X__'):
#            term_count += 1
#    featVec.append( math.exp(-1) )                       # add phrase penalty
#    featVec.append( -math.exp(-term_count) )              # add word penalty

    return ' '.join( map( lambda x: '%6f' % x, featVec ) )

def main():
    global opts
    if opts.file_index is None:             # if file_index is None; sent_file is directly specified
        sentfile = opts.sentfile
        filterfile = opts.filterfile
        filterlrm = "phr." + opts.filterfile
    else:                               # else sent_file is obtained from sent_dir and file_inedx
        sent_dir = opts.sentdir
        rule_dir = opts.filterdir
        if not sent_dir.endswith('/'): sent_dir = sent_dir + '/'
        if not rule_dir.endswith('/'): rule_dir = rule_dir + '/'
        sentfile = sent_dir + opts.file_index + '.out'
        filterfile = rule_dir + opts.file_index + '.out'
        filterlrm = rule_dir + "phr." + opts.file_index + '.out'
    
    loadTrie(opts.rulefile)
    processPhrases(sentfile)

    # write the filtered rules of the dev/test set
    writeFilteredRules(opts.rulefile, filterfile)

    # write the filtered phr rules of the dev/test set
    if opts.lex_reorder_model:
        writePhrLRM(opts.lrmfile, filterlrm)

if __name__ == "__main__":
    global opts
    optparser = optparse.OptionParser()
    optparser.add_option("-f", "--sentFile", dest="sentfile", default=None, help="data file for filtering (default=None)")
    optparser.add_option("-s", "--sentdir", dest="sentdir", default="sent", help="directory of filter files (default=sent)")
    optparser.add_option("-r", "--ruleFile", dest="rulefile", default=None, help="rule file (default=None)")
    optparser.add_option("", "--outfile", dest="filterfile", default="rules", help="filtered directory (default=rules)")
    optparser.add_option("-o", "--outdir", dest="filterdir", default="rules", help="filtered directory (default=rules)")
    optparser.add_option("-p", "--index", dest="file_index", default=None, help="index of filter data files (default=None)")
    optparser.add_option("-l", "--logfile", dest="log_file", default=None, help="filename for logging output")
    optparser.add_option("","--lrmFile", dest="lrmfile", default=None, help="lexicalized reordering model file (default=None)")
    optparser.add_option("", "--maxPhrLen", dest="max_phr_len", default=10, type="int", help="maximum initial phrase length (default=10)")
    optparser.add_option("", "--totSrcTrms", dest="tot_src_terms", default=7, type="int", help="maximum number of terms in src side of rules (default=7)")
    (opts, _) = optparser.parse_args()
    opts.lex_reorder_model = False
    if opts.lrmfile != None:
        opts.lex_reorder_model = True

    if opts.log_file:
        logging.basicConfig(filename=opts.log_file, filemode='w', level=logging.INFO)

    main()
