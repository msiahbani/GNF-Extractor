## This program extracts GNF rules (special type of Synchronous CFG rules) from the word alignments of a parallel corpus ##

import optparse, sys, os, logging, time 
from zgc import zgc

# Constants
weight_rules = False          # When distributing the unit-count among the rules, should it be weighted by the # of rule occurrences
max_terms = 7

# Global Variables
rule_indx = 1
tot_rules_derived = 0

srcWrds = []
tgtWrds = []
srcSentlen = 0
tgtSentLen = 0
ruleDict = {}                 # dictionary of rules for each sentence, ruleDict[(src, tgt)] = estimated rule count (1.0 for initial phrase pairs at the begining)
phrDictL2R = {}               # dictionary of phrases for each sentence, to keep LRM info: phrDictL2R[(src, tgt)] = 0, 1 or 2 (M:0, S:1, D:2) 
phrDictR2L = {}               # dictionary of phrases for each sentence, to keep LRM info: phrDictR2L[(src, tgt)] = 0, 1 or 2 (M:0, S:1, D:2)
LRMDictL2R = {}               # dictionary of phrases to keep LRM info: LRMDictL2R[(src, tgt)] = 0, 1 or 2 (M:0, S:1, D:2) 
LRMDictR2L = {}               # dictionary of phrases to keep LRM info: LRMDictR2L[(src, tgt)] = 0, 1 or 2 (M:0, S:1, D:2)
tgtPhrCntDict = {}            # dictionary of tgt phrases (w/o non-terminals) to keep counts of tgt phrases
ruleDoD = {}                  # dictionary of rules for each span, ruleDoD[(i, j)] = {(src, tgt):1,...}
nonTermRuleDoD = {}           # dictionary of rules for each span, which just contains non-terminals on target side
alignDoD = {}                 # Dict of dict to store fwd alignments
revAlignDoD = {}              # Dict of dict to store rev alignments
sentInitDoD = {}              # Dict of dict for storing initial phrase pairs, sentInitDoD[src_len] = {(src,tgt):1,...} (tuples of source and target spans)
strPhrDict = {}               # Dict of dict for storing initial phrase pairs (including loose phrases on tgt), strPhrDict[src_len] = {(src,tgt):1,...} (tuples of source and target spans)
rightTgtPhraseDict = {}
tgtPhraseDict = {}            # tgtPhraseDict[tgt_tuple] = (src_tuple, tgt_tuple)  ## just tight phrases
fullTgtPhraseDoD = {}         # fullTgtPhraseDoD[tgt_tuple] = (src_tuple, tgt_tuple) ## all phrases (loose and tight)
ppairRulesSet = set([])

tgtCntDict = {}
ruleIndxCntDict = {}
fAlignDoD = {}
rAlignDoD = {}


def readSentAlign():
    'Reads the input phrase span file for src & tgt sentences, alignment and initial phrases'

    global opts, nonTermRuleDoD, ruleDoD, LRMDictL2R, LRMDictR2L, tgtPhrCntDict
    global ruleDict, ruleIndxCntDict, tgtCntDict, phrPairLst, tgtPhraseDict
    global srcWrds, tgtWrds, srcSentlen, tgtSentLen, rightTgtPhraseDict, sentInitDoD, strPhrDict

    file_indx = opts.file_prefix
    dDir = opts.datadir
    oDir = opts.outdir
            
    if not dDir.endswith("/"): dDir += "/"
    if not oDir.endswith("/"): oDir += "/"

    spanFile = dDir + file_indx + '.outspan' 
    outFile = oDir + file_indx + '.out'
    outTgtFile = oDir + 'tgt.' + file_indx + '.out'
    outLRMFile = oDir + 'phr.' + file_indx + '.out'
    outLRMTgtFile = oDir + 'phr.tgt.' + file_indx + '.out'

    sent_count = 0
    phrLst = []
    aTupLst = []

    print "Using the maximum phrase lenght            :", opts.max_phr_len
    print "Enforcing tight phrase-pairs constraint    :", opts.tight_phrases_only
    print "Using the source side total terms to be    :", opts.tot_src_terms
    print "Using the mximum no of non-terminals to be :", opts.max_non_term

    print "Reading the span file                      :", spanFile
    inF = open(spanFile, 'r')
    for line in inF:
        line = line.strip()
        if line.startswith('LOG: SRC: '):             # Read the source sentence
            src = line.replace('LOG: SRC: ', '')
            srcWrds = src.split()
            srcSentlen = len(srcWrds)
        elif line.startswith('LOG: TGT: '):           # Read the target sentence
            tgt = line.replace('LOG: TGT: ', '')
            tgtWrds = tgt.split()
            tgtSentLen = len(tgtWrds)
        elif line.startswith('LOG: ALT:'):           # Read the source-target alignment
            align = line.replace('LOG: ALT:', '')
            align = align.strip()
            for align_pos in align.split():
                m = align_pos.split('-')
                e = -1 if m[0] == 'Z' else int(m[0])
                f = -1 if m[1] == 'Z' else int(m[1])
                aTupLst.append((e, f))
                try:                                  # Store forward alignments
                    alignDoD[m[0]][m[1]] = 1
                except KeyError:
                    alignDoD[m[0]] = {}
                    alignDoD[m[0]][m[1]] = 1
                try:                                  # Store reverse alignments
                    revAlignDoD[m[1]][m[0]] = 1
                except KeyError:
                    revAlignDoD[m[1]] = {}
                    revAlignDoD[m[1]][m[0]] = 1
        elif line.startswith('LOG: PHRASES_BEGIN:'): continue
        elif line.startswith('LOG: PHRASES_END:'):    # End of the current sentence; now extract rules from it
            align_tree = zgc(opts.max_phr_len)
            phrPairLst = align_tree.getAlignTree(srcSentlen, tgtSentLen, aTupLst)
            fullPhrPairLst = addLoosePhrases(phrPairLst)
            if not opts.tight_phrases_only: phrPairLst = fullPhrPairLst
            if opts.max_phr_len >= srcSentlen and not ((0, srcSentlen-1),(0, tgtSentLen-1)) in phrPairLst:
                phrPairLst.append(((0, srcSentlen-1),(0, tgtSentLen-1)))
            sentInitDoD = {}
            tgtPhraseDict = {}
            nonTermRuleDoD = {}
            strPhrDict = {}
            for ppair in phrPairLst:
                unaligned_edge = False
                # If the boundary term of source or target phrase has an unaligned word, ignore the phrase-pair
                # Unless the tight-phrase options is set to False
                if not alignDoD.has_key( str(ppair[0][0]) ) or not revAlignDoD.has_key( str(ppair[1][0]) ) or \
                        not alignDoD.has_key( str(ppair[0][1]) ) or not revAlignDoD.has_key( str(ppair[1][1]) ):
                    if opts.tight_phrases_only: continue
                    else: unaligned_edge = True
                
                
                init_phr_pair = (' '.join( [str(x) for x in xrange(ppair[0][0], ppair[0][1]+1) ] ), \
                        ' '.join( [str(x) for x in xrange(ppair[1][0], ppair[1][1]+1)] ) )
                if unaligned_edge:
                    if checkRuleConfigure(init_phr_pair):
                        if init_phr_pair in ruleDict: ruleDict[init_phr_pair] += 1.0
                        else: ruleDict[init_phr_pair] = 1.0
                    continue
                
                # Create a dict of dict for storing initial phrase pairs (tuples of source and target spans)
                tphr_len = ppair[1][1] - ppair[1][0] + 1
                if not sentInitDoD.has_key(tphr_len):
                    strPhrDict[tphr_len] = {}
                    sentInitDoD[tphr_len] = {}
                sentInitDoD[tphr_len][ppair] = init_phr_pair
                strPhrDict[tphr_len][ppair] = init_phr_pair
                tgtPhraseDict[ppair[1]] = ppair
                nonTermRuleDoD[ppair[1]] = {("X__1","X__1"):1}
            computeRightTgtPhr()
            xtractRules()
            #printSentRules()
            
            # For every extracted rule call the function compFeatureCounts() and compLRMFeature() to:
            # compFeatureCounts():
            #   i. convert the word positions in the rules into lexical entries, and
            #   ii. find the alignment for the rule and compute the joint count p(s, t)
            # compLRMFeature():
            #   i. convert the rule to a phrase (remove non-terminals) and compute L2R and R2L reordering model            
            
            for rule in ruleDict.keys():
                compFeatureCounts(rule)
                if opts.lex_reorder_model: compLRMFeature(rule)
            
            # Update global LRM
            if opts.lex_reorder_model:    updateLRMFeat()
            # Clear the variables at the end of current sentence
            resetStructs()
            del aTupLst[:]
            sent_count += 1
            if sent_count % 1000 == 0:
                print "Sentences processed : %6d ..." % sent_count
        else:
            continue  #temporary
           
    inF.close()

    # Write the rule counts, forward/reverse alignments and LRM to files
    with open(outFile, 'w') as oF:
        for rule in sorted( ruleIndxCntDict.iterkeys() ):
            r_indx, rule_count = ruleIndxCntDict[rule]
            f_alignments = ' ## '.join( fAlignDoD[r_indx].keys() )
            r_alignments = ' ## '.join( rAlignDoD[r_indx].keys() )
            oF.write( "%s ||| %g ||| %s ||| %s\n" % (rule, rule_count, r_alignments, f_alignments) )

    with open(outTgtFile, 'w') as tF:
        for tgt in sorted( tgtCntDict.iterkeys() ):
            tF.write( "%s ||| %g\n" % (tgt, tgtCntDict[tgt]) )
            
    if opts.lex_reorder_model:
        with open(outLRMFile, 'w') as lrmF:
            for rule in sorted( LRMDictL2R.iterkeys() ):
                l2r = ' '.join( map( lambda x: str(LRMDictL2R[rule][x]), range(3) ) )
                r2l = ' '.join( map( lambda x: str(LRMDictR2L[rule][x]), range(3) ) )
                lrmF.write( "%s ||| %s ||| %s ||| %s\n" % (rule[0], rule[1], l2r, r2l) )
        with open(outLRMTgtFile, 'w') as lrmTF:
            for tgt in sorted( tgtPhrCntDict.iterkeys() ):
                #tgtCnt = sum( map( lambda x: tgtPhrCntDict[tgt][x], range(3) ) )
                l2r = ' '.join( map( lambda x: str(tgtPhrCntDict[tgt][x]) , range(3) ) )
                r2l = ' '.join( map( lambda x: str(tgtPhrCntDict[tgt][x]) , range(3,6) ) )
                lrmTF.write( "%s ||| %s ||| %s\n" % (tgt, l2r, r2l) )
    return None

def addLoosePhrases(phr_lst):
    '''Add phrase-pairs with unaligned words on src/tgt boundary'''
    
    global alignDoD, revAlignDoD, tgtSentLen, srcSentlen, opts, fullTgtPhraseDoD
    full_phr_lst = set(phr_lst)
    st_length=[srcSentlen, tgtSentLen]
    alignDict = [alignDoD, revAlignDoD]          #dictionary of alignments (src/tgt)
    for tight_ppair in phr_lst:
        curr_lst = set()
        curr_lst.add(tight_ppair)
        for st_ind in [0,1]:                     #0: src,   1:tgt
            for p_ind,step in enumerate([-1, 1]):
                new_lst = set()
                for ppair in curr_lst:
                    j = tight_ppair[st_ind][p_ind] + step
                    while j >= 0 and j < st_length[st_ind]:
                        if alignDict[st_ind].has_key(str(j)): break                 #boundary reaches an aligned words 
                        if p_ind == 0: new_phr = (j, ppair[st_ind][1])
                        elif p_ind == 1: new_phr = (ppair[st_ind][0], j)
                        if abs(new_phr[1] - new_phr[0]) > opts.max_phr_len: break   #length of the new phrase gets larger than max_phr_len
                        if st_ind == 0: new_ppair = (new_phr, ppair[1])
                        elif st_ind == 1: new_ppair = (ppair[0], new_phr)
                        new_lst.add(new_ppair)
                        j += step
                curr_lst.update(new_lst)
        full_phr_lst.update(curr_lst)
    
    ## store the info for LRM computation
    fullTgtPhraseDoD = {}
    fullTgtPhraseDoD[(-1,-1)] = (-1,-1)     ## <s>
    fullTgtPhraseDoD[(tgtSentLen,tgtSentLen)] = (srcSentlen,srcSentlen)     ## <\s>
    for ppair in full_phr_lst:
        if ppair not in fullTgtPhraseDoD:
            fullTgtPhraseDoD[ppair[1]]={}
        fullTgtPhraseDoD[ppair[1]][ppair[1]] = 1
    return list(full_phr_lst)
                        
def resetStructs():
    '''Clean all data structures'''
    
    global alingDoD, nonTermRuleDoD, revAlignDoD, fullTgtPhraseDoD
    global rightTgtPhraseDict, ruleDict, ruleDoD, phrDictL2R, phrDictR2L
    global sentInitDoD, strPhrDict, tgtPhraseDict, basePhrDict
    alignDoD.clear()
    nonTermRuleDoD.clear()
    revAlignDoD.clear()
    rightTgtPhraseDict.clear()
    ruleDict.clear()
    ruleDoD.clear()
    sentInitDoD.clear()
    strPhrDict.clear()
    tgtPhraseDict.clear()
    basePhrDict.clear()
    fullTgtPhraseDoD.clear()
    phrDictL2R.clear()
    phrDictR2L.clear()

def computeRightTgtPhr():
    ''' fill the table rightTgtPhraseDict. For each possible span it keeps the largest 
    subphrase which shares the right boundary on target side. '''
    
    global rightTgtPhraseDict, tgtSentLen, tgtPhraseDict, sentInitDoD, opts
    rightTgtPhraseDict = {}
    for l in xrange(2, tgtSentLen+1):
        for i in xrange(0,tgtSentLen-l+1):
            j = i+l-1
            if l-1 in sentInitDoD and (i+1, j) in tgtPhraseDict:
                rightTgtPhraseDict[(i,j)] = tgtPhraseDict[(i+1, j)]
            elif (i+1, j) in rightTgtPhraseDict:
                rightTgtPhraseDict[(i,j)] = rightTgtPhraseDict[(i+1, j)]
            #if (i, j) in rightTgtPhraseDict:
            #    print (i, j), " : ", rightTgtPhraseDict[(i, j)]
    #if opts.tight_phrases_only: 
    fixUnAlignTgtWords()    
            
def fixUnAlignTgtWords():
    '''This function maps each phrase pair to phrase pairs which share the same source side but target sides just varies on left boundary
       (left boundaries are unaligned tgt words).
       Fills basePhrDict, a dictionary which maps tight tgt spans (with aligned boundaries) -> a set of phrase pairs 
          these phrase pairs share the same src side but the tgt sides are different, 
          tgt sides share the same right boundary (the same as key of dictionary) but left boundary varies (unaligned tgt words)'''
    
    global tgtSentLen, tgtPhraseDict, sentInitDoD, basePhrDict, strPhrDict
    basePhrDict = {}
    for i in xrange(0, tgtSentLen):
        if not revAlignDoD.has_key(str(i)):
            j = i - 1
            while j >= 0:
                if revAlignDoD.has_key(str(j)):
                    for l in xrange(0,j+1):
                        if l+1 in sentInitDoD and (j-l, j) in tgtPhraseDict:
                            ppair = (tgtPhraseDict[(j-l, j)][0],(j-l, i))
                            if (j-l, j) not in basePhrDict: basePhrDict[(j-l, j)] = {}
                            basePhrDict[(j-l, j)][ppair] = 1
                            #tgtPhraseDict[(j-l, i)] = ppair
                            if i-(j-l)+1 not in strPhrDict: strPhrDict[i-(j-l)+1] = {}
                            init_phr_pair = (' '.join( [str(x) for x in xrange(ppair[0][0], ppair[0][1]+1) ] ), \
                                             ' '.join( [str(x) for x in xrange(ppair[1][0], ppair[1][1]+1)] ) )                            
                            strPhrDict[i-(j-l)+1][ppair] = init_phr_pair
                    break
                j -= 1
    return

def printSentRules():
    global sentInitDoD, ruleDoD, tgtSentLen, ruleDict
    for tphr_len in xrange(1, tgtSentLen+1):
        if tphr_len in sentInitDoD:
            print "\nTARGET LENGTH:  ", tphr_len, "\n"
            for phr_pair in sentInitDoD[tphr_len]:
                print "Rules extracted from: ", phr_pair, sentInitDoD[tphr_len][phr_pair]
                for rule in ruleDoD[phr_pair[1]]:
                    try:
                        print " ||| ".join([rule[0], rule[1], str(ruleDict[rule])])
                    except KeyError:
                        print " ||| ".join([rule[0], rule[1], "0"])

def xtractRules():
    ''' Extract the rules for different phrase lengths (from smallest to longest) '''
    
    global srcSentlen, tgtSentLen, sentInitDoD, ruleDoD, rightTgtPhraseDict
    for tphr_len in xrange(1, tgtSentLen+1):
        if sentInitDoD.has_key(tphr_len):
            for phr_pair in sentInitDoD[tphr_len]:
                tgt_tuple = phr_pair[1]
                ruleDoD[tgt_tuple] = {}                 #target side of the phrase pair
                genRule4Phrase(sentInitDoD[tphr_len][phr_pair], phr_pair)

def genRule4Phrase(phrPairStr, (src_tuple, tgt_tuple)):
    global ruleDict, ruleDoD, rightTgtPhraseDict, ppairRulesSet, tgtPhraseDict, basePhrDict, opts
    
    ruleDoD[tgt_tuple] = {}   
    ppairRulesSet = set()
    if checkRuleConfigure(phrPairStr):
        #ruleDoD[tgt_tuple][phrPairStr] = 1         # we do not need to add terminal rules to ruleDoD :-?
        #ruleDict[phrPairStr] = 1.0
        ppairRulesSet.add(phrPairStr)
        
    (startP, endP) = tgt_tuple
    curr_rule_lst = set()
    curr_rule_lst.add(phrPairStr)         # it is a set, because current rule can be more than one rule (if there is unaligned src words)
    
    while endP >= startP:
        if (startP, endP) not in rightTgtPhraseDict:  break
        rtgt_phr_span = rightTgtPhraseDict[(startP, endP)]
        # generating rules by replacing the right most sub-phrase with its rules
        substituteRuleSet(curr_rule_lst, rtgt_phr_span)       
        # updating current rule list by replacing whole sub-phrase with one (or list of) non-terminal
        curr_rule_lst = substituteRuleSet(curr_rule_lst, rtgt_phr_span, "X")        
        
        if curr_rule_lst is None or len(curr_rule_lst) == 0: # it is not possible to generate more rule
            break
        endP = rtgt_phr_span[1][0] - 1

    # distribute unit rule count among rules extracted from this phrase pair
    ruleNo = 0
    for rule in ppairRulesSet:
        max_x = findMaxNonTerm(rule[0])
        if max_x <= opts.max_non_term:
            ruleNo += 1
    if ruleNo == 0: ruleCount = 0
    else:    ruleCount = 1.0/float(ruleNo)    
    for rule in ppairRulesSet:
        if rule[1].startswith("X__"):
            if not opts.full_aligned_rules: nonTermRuleDoD[tgt_tuple][rule] = 1
            continue    #target side does not have any terminal
        if not checkRuleConfigure(rule): continue
        ruleDoD[tgt_tuple][rule] = 1
        if rule in ruleDict: ruleDict[rule] += ruleCount
        else: ruleDict[rule] = ruleCount
    if phrPairStr in ppairRulesSet and ruleCount < 1.0:
        ruleDict[phrPairStr] += 1.0
        
    # updating right most tgt phrase for this tgt
    rightTgtPhraseDict[tgt_tuple] = (src_tuple, tgt_tuple)
    
    # update phrases with unaligned words on right boundry of tgt side
    if tgt_tuple in basePhrDict:
        for ppair in basePhrDict[tgt_tuple]:
            rightTgtPhraseDict[ppair[1]] = ppair
            nonTermRuleDoD[ppair[1]] = nonTermRuleDoD[tgt_tuple]
            ruleDoD[ppair[1]] = ruleDoD[tgt_tuple]
        
def checkRuleConfigure((src_phr, tgt_phr), isNonTerm=False):
    '''Checks if the rule configuration is compatible with the constraints (for both src & tgt sides)'''
    
    global opts, max_terms
    # source lenght
    src_len = len(src_phr.split())
    if src_len > opts.tot_src_terms: return False
    # no of non-terminals
    if isNonTerm and tgt_phr.startswith("X__"): return True
    max_x = findMaxNonTerm(src_phr)
    if max_x > opts.max_non_term:  return False 
    if src_len - max_x > max_terms: return False
    if len(tgt_phr.split()) - max_x > max_terms+3: return False
    pre_x = False
    
    #check adjacent non-terms on src
    if src_phr.find("X__2") >= 0:
        for w in src_phr.split():
            if w.startswith("X__"):
                if pre_x: return False
                pre_x = True
            else: pre_x = False
    return True

def substituteRuleSet(main_rule_lst, sub_ppair_span, isNonTerm = None):
    ''' input: a list of main rules and a subphrase [if isNonTerm == 'X' substituting subphrase with a list of non-terminals]
        task:  substituting subphrase with its rules in the main rules
               updating the final main rules by substituting the subphrase with a non-terminal
               returning updated main rules (or None, if it is not possible)'''
    
    global ppairRulesSet, ruleDoD, nonTermRuleDoD, strPhrDict, opts

    len_tgt_phr = sub_ppair_span[1][1] - sub_ppair_span[1][0] +1
    sub_phr = strPhrDict[len_tgt_phr][sub_ppair_span]                        # retrive subphrase using its span        
    update_rule_lst = set()
    if isNonTerm: local_rule_dict = nonTermRuleDoD[sub_ppair_span[1]]
    else:         local_rule_dict = ruleDoD[sub_ppair_span[1]]
    
    for main_rule in main_rule_lst:
        max_x = findMaxNonTerm(main_rule[0])
        if not isNonTerm and max_x >= opts.max_non_term:       # rule cannot have more non-terminal
            continue
    
        if main_rule[0].startswith(sub_phr[0]): left_bound_src = 0                            # happens at the begining of src side 
        elif main_rule[0].find(" "+sub_phr[0]+" ") > 0: left_bound_src = main_rule[0].find(" "+sub_phr[0]+" ")+1  
        else: left_bound_src = main_rule[0].find(" "+sub_phr[0])+1                              # find boundaries subphrase on src side 
        right_bound_src = left_bound_src + len(sub_phr[0])
        if main_rule[1].startswith(sub_phr[1]): left_bound_tgt = 0
        else:   left_bound_tgt = main_rule[1].find(" "+sub_phr[1])+1                            # find boundaries subphrase on tgt side
        right_bound_tgt = left_bound_tgt + len(sub_phr[1])
    
        # left and right part of main rule (for later rule construction)
        right_src_main = main_rule[0][right_bound_src:].strip()
        right_tgt_main = main_rule[1][right_bound_tgt:].strip()
        left_src_main = main_rule[0][:left_bound_src].strip()
        left_tgt_main = main_rule[1][:left_bound_tgt].strip()
    
        start_x = findMaxNonTerm(main_rule[0], 0, left_bound_src)
        
        # check if the first (last) term in right (left) part of source rule is non-terminal (to check the filtering constraints)
        right_non_term = True if right_src_main.startswith("X__") else False
        left_non_term = True if left_src_main.endswith("X__"+str(start_x)) else False
    
        for rule in local_rule_dict:
            mid_x = findMaxNonTerm(rule[0])
            if mid_x == 0 or (not isNonTerm and (mid_x+max_x) > opts.max_non_term): continue        # (does not create new rule) or (more than maximum non-terminals)
        
            #check adjacent non-terminals on src side
            if not isNonTerm and ( (right_non_term and rule[0].endswith("X__"+str(mid_x))) or \
               (left_non_term and rule[0].startswith("X__1")) ):
                continue
        
            mid_src_phr = updateNonTerms(rule[0], start_x)
            mid_tgt_phr = updateNonTerms(rule[1], start_x)
            
            right_src_phr = updateNonTerms(right_src_main, mid_x)
            right_tgt_phr = updateNonTerms(right_tgt_main, mid_x, start_x)
            new_rule = ((" ".join([left_src_main, mid_src_phr, right_src_phr]).strip()), \
                        " ".join([left_tgt_main, mid_tgt_phr, right_tgt_phr]).strip())
            #if checkRuleConfigure(new_rule, isNonTerm):
            if (right_non_term and rule[0].endswith("X__"+str(mid_x))) or \
                      (left_non_term and rule[0].startswith("X__1")):
                pass
            else:
                ppairRulesSet.add(new_rule)      
                update_rule_lst.add(new_rule)
                
            # if the subphrase is replaced with just non-terminal:
            #    check if the first two nonterminals on the target side can be merge 
            #    (if there are just unaligned words between corresponding non-terminals on the src side)
            if isNonTerm:
                new_rule_lst = iterativeMerge(new_rule, start_x+1)
                for new_rule in new_rule_lst:
                    #if checkRuleConfigure(new_rule, isNonTerm):
                    ppairRulesSet.add(new_rule)
                    update_rule_lst.add(new_rule)
            
    return update_rule_lst


def iterativeMerge(rule, no_x):
    '''Check if begining nonterminals in target side of rule can be merged to create different rules'''
    
    global alignDoD
    curr_rule = rule
    new_rules = set()
    
    first_x = curr_rule[1].find("X__")   # tgt might start with some terminal term
    if first_x < 0 : return new_rules
    
    while True:
        t_words = curr_rule[1][first_x:].split()
        if len(t_words) < 2: break    
        find_phr_flag = False
        max_x = int(t_words[0][3:])
        min_x = max_x
        for i in xrange(1,len(t_words)):
            x_ind = int(t_words[i][3:])
            if x_ind < min_x: min_x = x_ind
            elif x_ind > max_x: max_x = x_ind
            new_r = mergeNonTerms(curr_rule, t_words[:i+1], min_x, max_x)
            if new_r:
                new_rules.add(new_r)
                find_phr_flag = True
                break
        if not find_phr_flag: return new_rules
        curr_rule = new_r
    return new_rules

def mergeNonTerms(rule, t_words, min_x, max_x):
    ''' Replacing a list of non-terminals t-words (with max and min nonterminals given) in the given rule '''
    
    curr_x = "X__"+str(min_x)
    next_x = "X__"+str(max_x)
    s = rule[0].find(curr_x)+5
    e = rule[0].find(next_x)-1    
    # if the next (or prev) non-term is adjacent, it is not acceptable anyway
    if rule[0][min(e+6, len(rule[0])):].startswith("X__") or \
       (min_x > 1 and rule[0][:max(s-6, 0)].endswith("X__"+str(min_x-1))): return None
    aligned = False
    for w in rule[0][s:e].strip().split():
        if w in alignDoD or (w.startswith("X__") and w not in t_words):
            aligned = True
            break
    if aligned: return None
    # find the right and left part of src (including the min non-terminal of the phrase)
    right_src = rule[0][min(e+6, len(rule[0])):].strip()
    right_src = updateNonTerms(right_src, min_x-max_x)
    left_src = rule[0][:s].strip()
    # find the right and left part of tgt     
    s = rule[1].find(t_words[0])
    #e = rule[1].find(t_words[-1])+5
    e = s+len(" ".join(t_words))
    right_tgt = rule[1][min(e, len(rule[1])):].strip()
    right_tgt = updateNonTerms(right_tgt, min_x-max_x, min_x)
    left_tgt = rule[1][:s].strip()
    return (" ".join([left_src,right_src]).strip(), " ".join([left_tgt, curr_x, right_tgt]).strip())
            
def updateNonTerms(phrase, add_x, start_x=0):
    '''Update non-terminals (from X__(start_x)) by adding add_x to the name of non-terminal'''
    
    if add_x == 0 or phrase.find("X__") < 0:
        return phrase
    tmpLst = []
    for w in phrase.split():
        if w.startswith("X__") and int(w[3:]) > start_x:
            tmpLst.append("X__"+str(add_x+int(w[3:])))
        else:
            tmpLst.append(w)    
    return " ".join(tmpLst)

def findMaxNonTerm(phrase, s=None, e=None):
    '''Find the index of bigest non-terminal in phrase[s:e] (the whole phase if s or e is not defined)'''
    
    if s !=None and e != None:
        start_x = phrase.rfind("X__", s, e)
    else: start_x = phrase.rfind("X__")
    
    if start_x > -1: start_x = int(phrase[start_x+3:min(start_x+5, len(phrase))])
    else: start_x = 0
    
    return start_x

def compLRMFeature(rule):
    '''Compute the position of target side phrase in compare to nighbour tgt phrases (before and after).
       The position is computed in terms of source side, Monotone (0), Swape (1) or Discontinuous (2).
       Information is saved in phrDictL2R (or phrDictR2L) for Left-to-right (or R2L).'''

    global fullTgtPhraseDoD, tgtSentLen, srcSentlen, phrDictL2R, phrDictR2L, srcWrds, tgtWrds
    srcPosLst = []
    tgtPosLst = []
    # Convert the word positions in source side of the rule to corresponding lexemes
    for i,s_tok in enumerate(rule[0].split()):
        if s_tok.startswith('X__'):
            if i!=0 and i!=len(rule[0].split())-1:
                srcPosLst.append("NON_TOK")
        else: srcPosLst.append(int(s_tok))
    # Convert the word positions in target side of the rule to corresponding lexemes
    for t_tok in rule[1].split():
        if t_tok.startswith('X__'):
            break
        tgtPosLst.append(int(t_tok))
        
    tgtPhr = " ".join([str(i) for i in tgtPosLst])
    srcPhr = " ".join([str(i) for i in srcPosLst])
    if (srcPhr, tgtPhr) in phrDictR2L:                ## LRM has been computed for this phr in this sentence
        return       
    
    # Compute L2R lexicalized reordering model
    prev_t = tgtPosLst[0] - 1
    l2r = {}
    while prev_t > -1:
        if (prev_t, tgtPosLst[0]-1) in fullTgtPhraseDoD:
            for src in fullTgtPhraseDoD[(prev_t, tgtPosLst[0]-1)]:
                if src[1] == srcPosLst[0]-1:        ## Monotone
                    l2r[0] = 1
                elif src[0] == srcPosLst[-1]+1:     ## Swape
                    l2r[1] = 1
            #break     ## word level LRM (but considering unaligned word)
                      ## remove break -> phrase LRM
        prev_t -= 1
    ## compare to <s>
    if tgtPosLst[0] == 0 and srcPosLst[0] == 0: l2r[0] = 1
    if len(l2r) == 0: l2r[2] = 1    ## Discontinuous
    phrDictL2R[(srcPhr, tgtPhr)] = l2r
        
    # Compute R2L lexicalized reordering model
    next_t = tgtPosLst[-1] + 1
    r2l = {}
    while next_t < tgtSentLen:
        if (tgtPosLst[-1]+1, next_t) in fullTgtPhraseDoD:
            for src in fullTgtPhraseDoD[(tgtPosLst[-1]+1, next_t)]:
                if src[0] == srcPosLst[-1]+1:        ## Monotone
                    r2l[0] = 1
                elif src[1] == srcPosLst[0]-1:       ## Swape
                    r2l[1] = 1
            #break     ## word level LRM (but considering unaligned word)
                       ## remove break -> phrase LRM
        next_t += 1    
    ## compare to <\s>
    if tgtPosLst[-1] == tgtSentLen-1 and srcPosLst[-1] == srcSentlen-1: r2l[0] = 1
    if len(r2l) == 0: r2l[2] = 1    ## Discontinuous
    phrDictR2L[(srcPhr, tgtPhr)] = r2l

def updateLRMFeat():
    ''' Updates LRM feats for all phrases of the current sentence '''
    
    global phrDictL2R, phrDictR2L, LRMDictL2R, LRMDictR2L, srcWrds, tgtWrds, tgtPhrCntDict
    for (src, tgt) in phrDictL2R:
        srcLst = []
        tgtPhr = " ".join([tgtWrds[int(t_tok)] for t_tok in tgt.split()])
        for s_tok in src.split():
            if s_tok != "NON_TOK": srcLst.append(srcWrds[int(s_tok)])
            else: srcLst.append(s_tok)
        srcPhr = " ".join(srcLst)
        if (srcPhr, tgtPhr) not in LRMDictL2R:
            LRMDictL2R[(srcPhr, tgtPhr)] = {0:0, 1:0, 2:0}
            LRMDictR2L[(srcPhr, tgtPhr)] = {0:0, 1:0, 2:0}
        if tgtPhr not in tgtPhrCntDict: tgtPhrCntDict[tgtPhr] = {0:0, 1:0, 2:0, 3:0, 4:0, 5:0}
        for k in LRMDictL2R[(srcPhr, tgtPhr)]:
            LRMDictL2R[(srcPhr, tgtPhr)][k] += phrDictL2R[(src, tgt)].get(k, 0)
            LRMDictR2L[(srcPhr, tgtPhr)][k] += phrDictR2L[(src, tgt)].get(k, 0)
            # update the tgt count
            tgtPhrCntDict[tgtPhr][k] += phrDictL2R[(src, tgt)].get(k, 0)
            tgtPhrCntDict[tgtPhr][k+3] += phrDictR2L[(src, tgt)].get(k, 0)
            
def compFeatureCounts(rule):
    '''Convert to lexical rule and find the alignment for the entries in the rule. Also compute feature counts P(s|t), P(t|s), P_w(s|t) and P_w(t|s)'''

    global srcWrds, tgtWrds
    global fAlignDoD, rAlignDoD
    srcLexLst = []
    tgtLexLst = []
    alignLst = []

    sPosLst = rule[0].split()
    tPosLst = rule[1].split()
    # Convert the word positions in source side of the rule to corresponding lexemes
    item_indx = 0
    for s_tok in sPosLst:
        if s_tok.startswith('X__'):
            srcLexLst.append(s_tok)
        else:
            srcLexLst.append(srcWrds[int(s_tok)])
            # Find the forward alignment for the lexemes in the rule
            alignment = getFwrdAlignment(item_indx, s_tok, tPosLst)
            alignLst.append(alignment)
            #if len(alignment) > 0:
            #    alignLst.append(alignment)
        item_indx += 1
    fAlignment = ' '.join(alignLst)

    # Convert the word positions in target side of the rule to corresponding lexemes
    del alignLst[:]
    item_indx = 0
    for t_tok in tPosLst:
        if t_tok.startswith('X__'):
            tgtLexLst.append(t_tok)
        else:
            tgtLexLst.append(tgtWrds[int(t_tok)])
            # Find the reverse alignment for the lexemes in the rule
            alignment = getRvrsAlignment(item_indx, t_tok, sPosLst)
            alignLst.append(alignment)
            #if len(alignment) > 0:
            #    alignLst.append(alignment)
        item_indx += 1
    rAlignment = ' '.join(alignLst)

    # Get the lexical rule and add its count from the current sentence to total count so far
    curr_rindx = updateRuleCount(' '.join(srcLexLst), ' '.join(tgtLexLst), rule)

    # Update forward and reverse alignment dicts
    f_align_indx = getAlignIndex(fAlignment)
    r_align_indx = getAlignIndex(rAlignment)
    if not fAlignDoD.has_key(curr_rindx):
        fAlignDoD[curr_rindx] = {}
        rAlignDoD[curr_rindx] = {}
    if not fAlignDoD[curr_rindx].has_key(f_align_indx):
        fAlignDoD[curr_rindx][f_align_indx] = 1
    if not rAlignDoD[curr_rindx].has_key(r_align_indx):
        rAlignDoD[curr_rindx][r_align_indx] = 1


def updateRuleCount(mc_src, mc_tgt, rule):
    ''' Updates rule and target counts '''

    global rule_indx, ruleDict, ruleIndxCntDict, tgtCntDict
    if not mc_tgt in tgtCntDict:
        tgtCntDict[mc_tgt] = 0
    tgtCntDict[mc_tgt] += ruleDict[rule]

    mc_key = mc_src + ' ||| ' + mc_tgt              # ' ||| ' is the delimiter separating items in the key/value
    if mc_key in ruleIndxCntDict:
        curr_rindx, curr_cnt = ruleIndxCntDict[mc_key]
        ruleIndxCntDict[mc_key] = ( curr_rindx, curr_cnt + ruleDict[rule] )
    else:
        ruleIndxCntDict[mc_key] = (rule_indx, ruleDict[rule])
        curr_rindx = rule_indx
        rule_indx += 1
    return curr_rindx

def getAlignIndex(align_str):

    tmpLst = align_str.split(' ')
    tmpLst.sort()
    aindx = ''.join(tmpLst)
    return aindx.replace('-', '')


def getFwrdAlignment(item_indx, s_pos, tPosLst):
    '''Computes the alignment and lexical weights in forward direction'''

    alignLst = []
    if alignDoD.has_key(s_pos):
        alignKeyLst = alignDoD[s_pos].keys()
        alignKeyLst.sort()
        for t_pos in alignKeyLst:
            try:
                # Get the alignment and append it to the list
                alignment = str(item_indx) + '-' + str(tPosLst.index(t_pos))
                alignLst.append(alignment)
            except ValueError:
                pass
    else:
        alignLst.append( str(item_indx) + '-Z' )     # 'Z' represents 'NULL' (i.e. word is unaligned)

    return ' '.join(alignLst)


def getRvrsAlignment(item_indx, t_pos, sPosLst):
    '''Computes the alignment and lexical weights in reverse direction'''

    alignLst = []
    if revAlignDoD.has_key(t_pos):
        alignKeyLst = revAlignDoD[t_pos].keys()
        alignKeyLst.sort()
        for s_pos in alignKeyLst:
            try:
                # Get the alignment and append it to the list
                alignment = str(sPosLst.index(s_pos)) + '-' + str(item_indx)
                alignLst.append(alignment)
            except ValueError:
                pass
    else:
        alignLst.append( 'Z-' + str(item_indx) )     # 'Z' represents 'NULL' (i.e. word is unaligned)

    return ' '.join(alignLst)

if __name__ == '__main__':
    global opts
    optparser = optparse.OptionParser()
    optparser.add_option("-d", "--datadir", dest="datadir", default="data", help="data directory (default=data)")
    optparser.add_option("-o", "--outdir", dest="outdir", default="rules", help="data directory (default=rules)")
    optparser.add_option("-p", "--prefix", dest="file_prefix", default="1", help="prefix of parallel data files (default=1)")
    optparser.add_option("-l", "--logfile", dest="log_file", default=None, help="filename for logging output")
    optparser.add_option("","--tightPhrase", dest="tight_phrases_only", default=False, action="store_true", help="extract just tight-phrases (default=False)")
    optparser.add_option("","--fullAlignedRules", dest="full_aligned_rules", default=False, action="store_true", help="extract just rules with aligned subphrases (Phrasal-Hiero in (Nguyen and Vogel 2013)) (default=False)")
    optparser.add_option("","--lexReorderingModel", dest="lex_reorder_model", default=False, action="store_true", help="compute lexicalized reordering model for phrases (default=False)")
    optparser.add_option("", "--maxPhrLen", dest="max_phr_len", default=10, type="int", help="maximum initial phrase length (default=10)")
    optparser.add_option("", "--totSrcTrms", dest="tot_src_terms", default=7, type="int", help="maximum number of terms in src side of rules (default=7)")
    optparser.add_option("", "--maxNonTrms", dest="max_non_term", default=2, type="int", help="maximum number of non-terminals in rules (default=2)")
    (opts, _) = optparser.parse_args()
    
    if opts.log_file:
        logging.basicConfig(filename=opts.log_file, filemode='w', level=logging.INFO)

    readSentAlign()
