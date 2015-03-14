## This program is the phase-III of the SCFG rules extraction ##
## It estimates the forward (P(s|t)) and reverse (P(t|s)) probabilities using relative frequency estimation ##

import os
import sys
import heapq
import math
import time
from datetime import timedelta

#min_lprob = -6.0           # for log10
#max_lprob = -0.000434
min_lprob = -13.8155        # for log (natural)
max_lprob = -0.000100005
tgt_trie = None
rulesLst = []
tgtCntDict = {}

def loadTgtCnts2Dict(consTgtFile):
    '''Loads target counts in a dict'''

    global tgtCntDict
    tot_rules = 0
    tot_time = 0.0
    tgtCntDict = {}

    print "Reading target counts and loading the counts in a dict ..."
    rF = open(consTgtFile, 'r')
    t_beg = time.time()
    while True:
        line = rF.readline()
        line = line.strip()
        if line == '': break

        tot_rules += 1
        (tgt, r_cnt) = line.split(' ||| ')
        tgtCntDict[tgt] = float( r_cnt )

        # track progress (for every ten million rules)
        if (tot_rules % 10000000) == 0:
	    t_taken = time.time() - t_beg
            tot_time += t_taken
            print "Processed %8d rules in time %s" % ( tot_rules, timedelta(seconds=t_taken) )
	    t_beg = time.time()

    rF.close()
    tot_time += time.time() - t_beg
    print "Total # of unique rules processed   : %d" % tot_rules
    print "Total time taken                    : %s" % timedelta(seconds=tot_time)

def writeFeats(ruleFile, outFile):

    global rulesLst
    tot_rules = 0
    tot_time = 0.0
    prev_src = ''
    src_cnt = 0.0
    rulesLst = []

    print "\n\nComputing source cnt and feature values before writing them to file ..."
    rF = open(ruleFile, 'r')
    oF = open(outFile, 'w')
    t_beg = time.time()
    while True:
        line = rF.readline()
        line = line.strip()
        if line == '': break

        tot_rules += 1
        (src_rule, tgt_rule, r_cnt, r_lprob, f_lprob) = line.split(' ||| ')
        rule_cnt = float(r_cnt)
        f_lex_prob = float(f_lprob)
        r_lex_prob = float(r_lprob)

        if prev_src != src_rule and tot_rules > 1:
            # New unique src_rule found; flush the rulesLst into file
            flush2File(src_cnt, oF)

            # Clear the src_cnt and rulesLst for next unique source rule
            src_cnt = 0.0
            rulesLst = []

        src_cnt += rule_cnt
        prev_src = src_rule
        rulesLst.append( (src_rule, tgt_rule, rule_cnt, f_lex_prob, r_lex_prob) )

        # tracking progress (for every million rules)
        if (tot_rules % 1000000) == 0:
	    t_taken = time.time() - t_beg
            tot_time += t_taken
            print "Processed %d million rules in %.4f sec" % (tot_rules / 1000000, t_taken)
	    t_beg = time.time()

    # flush the final rule after the last line is read
    flush2File(src_cnt, oF)

    rF.close()
    oF.close()
    tot_time += time.time() - t_beg
    print "Total # of unique rules processed   : %d" % tot_rules
    print "Total time taken                    : %f" % tot_time

def flush2File(src_cnt, oF):
    '''Flush the rules accumulated so far to the file'''

    global rulesLst
    global tgtCntDict

    for (src, tgt, r_cnt, r_lprob, f_lprob) in rulesLst:
        tgt_cnt = tgtCntDict[tgt]
        if ( tgt_cnt < r_cnt ):
            tgt_cnt = r_cnt
            #print "Log:: Rule: %s ||| %s :: tgt_cnt (%e) is smaller than r_cnt (%e)\n" % (src, tgt, tgt_cnt, r_cnt)

        # Compute the 4 features
        if r_cnt == 0.0:        # use min_lprob if r_cnt is zero
            f_p = min_lprob
            r_p = min_lprob
        else:
            if r_cnt == tgt_cnt: f_p = max_lprob    # use max_lprob if prob = 1 (for fwrd prob)
            else: f_p = math.log( r_cnt/ tgt_cnt )
            if r_cnt == src_cnt: r_p = max_lprob    # use max_lprob if prob = 1 (for rvrs prob)
            else: r_p = math.log( r_cnt/ src_cnt )

        if f_lprob == 0.0: f_lp = min_lprob
        elif f_lprob == 1.0: f_lp = max_lprob
        else: f_lp = math.log( f_lprob )

        if r_lprob == 0.0: r_lp = min_lprob
        elif r_lprob == 1.0: r_lp = max_lprob
        else: r_lp = math.log( r_lprob )

        if (f_p > 0.0 or r_p > 0.0 or f_lp > 0.0 or r_lp > 0.0):         # Check that the log-prob values are in fact negative
            print "** ", src, " ||| ", tgt, " :: ", src_cnt, tgt_cnt, r_cnt, " : ", f_lprob, r_lprob
            assert (f_p < 0.0), "ERROR: *1 - Log-prob value for forward prob is Positive : %g. Exiting!!\n" % (f_p)
            assert (r_p < 0.0), "ERROR: *2 - Log-prob value for reverse prob is Positive : %g. Exiting!!\n" % (r_p)
            assert (f_lp < 0.0), "ERROR: *3 - Log-prob value for forward lexical prob is Positive : %g. Exiting!!\n" % (f_lp)
            assert (r_lp < 0.0), "ERROR: *4 - Log-prob value for reverse lexical prob is Positive : %g. Exiting!!\n" % (r_lp)

        # Write the features to the file
        oF.write( "%s ||| %s ||| %g %g %g %g\n" % (src, tgt, r_p, f_p, r_lp, f_lp) )

def loadTgtPhr(tgtPhrFile):
    '''Loads target phrase counts in a dict'''

    global tgtPhrDict
    tot_rules = 0
    tot_time = 0.0
    tgtPhrDict = {}

    print "Reading target phrase counts and loading the orientation counts in a dict ..."
    rF = open(tgtPhrFile, 'r')
    t_beg = time.time()
    while True:
        line = rF.readline()
        line = line.strip()
        if line == '': break

        tot_rules += 1
        (tgt, l2r, r2l) = line.split(' ||| ')
        l2r = [float(x) for x in l2r.split()]
        r2l = [float(x) for x in r2l.split()]
        tgtPhrDict[tgt] = (sum(l2r), l2r, r2l)

        # track progress (for every ten million rules)
        if (tot_rules % 10000000) == 0:
	    t_taken = time.time() - t_beg
            tot_time += t_taken
            print "Processed %8d rules in time %s" % ( tot_rules, timedelta(seconds=t_taken) )
	    t_beg = time.time()

    rF.close()
    tot_time += time.time() - t_beg
    print "Total # of unique rules processed   : %d" % tot_rules
    print "Total time taken                    : %s" % timedelta(seconds=tot_time)

def computeFinalLRM(cntFile, phrFile, outFile):

    global tgtPhrDict
    tot_phr = 0
    tot_time = 0.0
    prev_src = ''
    src_cnt = 0.0
    phrLst = []
    alpha_u = 10.0
    alpha_g = 10.0
    alpha_s = 10.0
    alpha_t = 10.0

    print "\nReading total orientation counts for smoothing"
    cF = open(cntFile, 'r')
    line = cF.readline().strip()
    (tot_phr_cnt, tot_l2r, tot_r2l) = line.split(' ||| ')
    tot_phr_cnt = float(tot_phr_cnt)
    tot_l2r = [float(x) for x in tot_l2r.split()]
    tot_r2l = [float(x) for x in tot_r2l.split()]
    prob_l2r = [(i + alpha_u/3.0)/(tot_phr_cnt + alpha_u) for i in tot_l2r]
    prob_r2l = [(i + alpha_u/3.0)/(tot_phr_cnt + alpha_u) for i in tot_r2l]
    cF.close()
    print "\n\nComputing source cnt and feature values before writing them to file ..."
    rF = open(phrFile, 'r')
    oF = open(outFile, 'w')
    t_beg = time.time()
    while True:
        line = rF.readline()
        line = line.strip()
        if line == '': break

        tot_phr += 1
        (src_phr, tgt_phr, cur_l2r, cur_r2l) = line.split(' ||| ')
        cur_l2r = [float(x) for x in cur_l2r.split()]
        cur_r2l = [float(x) for x in cur_r2l.split()]
        cur_cnt = sum(cur_l2r)

        if prev_src != src_phr and tot_phr > 1:
            # New unique src_phr found; flush the phrLst into file
            s_l2r = [0.0, 0.0, 0.0]
            s_r2l = [0.0, 0.0, 0.0]
            for (src, tgt, cnt, l2r, r2l) in phrLst:
                for i in range(3): 
                    s_l2r[i] += l2r[i]
                    s_r2l[i] += r2l[i]

            ps_l2r = [ (s_l2r[i] + alpha_g*prob_l2r[i])/(src_cnt+alpha_g) for i in range(3) ]
            ps_r2l = [ (s_r2l[i] + alpha_g*prob_r2l[i])/(src_cnt+alpha_g) for i in range(3) ]
            for (src, tgt, cnt, l2r, r2l) in phrLst:
                (t_cnt, t_l2r, t_r2l) = tgtPhrDict[tgt]
                p_l2r = []
                p_r2l = []
                for i in range(3):
                    pt_l2r = (t_l2r[i] + alpha_g*prob_l2r[i])/(t_cnt+alpha_g)
                    pt_r2l = (t_r2l[i] + alpha_g*prob_r2l[i])/(t_cnt+alpha_g)

                    p_l2r.append( math.log((l2r[i] + alpha_s*ps_l2r[i] + alpha_t*pt_l2r)/(cnt + alpha_s + alpha_t)) )
                    p_r2l.append( math.log((r2l[i] + alpha_s*ps_r2l[i] + alpha_t*pt_r2l)/(cnt + alpha_s + alpha_t)) )
                
                oF.write( "%s ||| %s ||| %g %g %g ||| %g %g %g\n" % (src, tgt, p_l2r[0], p_l2r[1], p_l2r[2], p_r2l[0], p_r2l[1], p_r2l[2]) )
            # Clear the src_cnt and phrLst for next unique source phr
            src_cnt = 0.0
            phrLst = []

        src_cnt += cur_cnt
        prev_src = src_phr
        phrLst.append( (src_phr, tgt_phr, cur_cnt, cur_l2r, cur_r2l) )

        # tracking progress (for every million phrases)
        if (tot_phr % 1000000) == 0:
	    t_taken = time.time() - t_beg
            tot_time += t_taken
            print "Processed %d million phrases in %.4f sec" % (tot_phr / 1000000, t_taken)
	    t_beg = time.time()

    # flush the final phrase after the last line is read
    s_l2r = [0.0, 0.0, 0.0]
    s_r2l = [0.0, 0.0, 0.0]
    for (src, tgt, cnt, l2r, r2l) in phrLst:
        for i in range(3): 
            s_l2r[i] += l2r[i]
            s_r2l[i] += r2l[i]

    ps_l2r = [ (s_l2r[i] + alpha_g*prob_l2r[i])/(src_cnt+alpha_g) for i in range(3) ]
    ps_r2l = [ (s_r2l[i] + alpha_g*prob_r2l[i])/(src_cnt+alpha_g) for i in range(3) ]
    for (src, tgt, cnt, l2r, r2l) in phrLst:
        (t_cnt, t_l2r, t_r2l) = tgtPhrDict[tgt]
        p_l2r = []
        p_r2l = []
        for i in range(3):
            pt_l2r = (t_l2r[i] + alpha_g*prob_l2r[i])/(t_cnt+alpha_g)
            pt_r2l = (t_r2l[i] + alpha_g*prob_r2l[i])/(t_cnt+alpha_g)

            p_l2r.append( math.log((l2r[i] + alpha_s*ps_l2r[i] + alpha_t*pt_l2r)/(cnt + alpha_s + alpha_t)) )
            p_r2l.append( math.log((r2l[i] + alpha_s*ps_r2l[i] + alpha_t*pt_r2l)/(cnt + alpha_s + alpha_t)) )
                
        oF.write( "%s ||| %s ||| %g %g %g ||| %g %g %g\n" % (src, tgt, p_l2r[0], p_l2r[1], p_l2r[2], p_r2l[0], p_r2l[1], p_r2l[2]) )


    rF.close()
    oF.close()
    tot_time += time.time() - t_beg
    print "Total # of unique phrases processed   : %d" % tot_phr
    print "Total time taken                    : %f" % tot_time


def main():
    ruleFile = sys.argv[1]
    ruleDir = os.path.dirname(ruleFile)
    outFile = ruleDir + "/" + sys.argv[2]
    tgtFile = ruleDir + "/tgt_rules.all.out"

    # load target counts into dict
    loadTgtCnts2Dict(tgtFile)
    # compute and write the features in the outFile
    writeFeats(ruleFile, outFile)
    

    cntFile = ruleDir + "/cnt.all.out" 
    if os.path.isfile(cntFile):
        phrFile = ruleDir + "/phr_cnt.all.out"
        tgtPhrFile = ruleDir + "/tgt_phr.all.out"
        lrmFile = ruleDir + "/phr_lprob.all.out"
        loadTgtPhr(tgtPhrFile)
        computeFinalLRM(cntFile, phrFile, lrmFile)

if __name__ == '__main__':
    main()
