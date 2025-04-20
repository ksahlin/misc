
import random
import argparse
import sys
import re
import os
import errno
import shutil
import glob
import math

# import pysam

def mkdir_p(path):
    print("creating", path)
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def reverse_complement(string):
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}
    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)



def get_breakpoint_site(read, bp_site):
    aln_start = read.reference_start
    read_cigar_tuples = []
    result = re.split(r'[=DXSMINH]+', read.cigarstring)
    i = 0
    for length in result[:-1]:
        i += len(length)
        type_ = read.cigarstring[i]
        i += 1
        read_cigar_tuples.append((int(length), type_ ))  

    t = read_cigar_tuples[0][1]
    l = read_cigar_tuples[0][0]
    if t == "=" or t== "M" or t == "X": # start of read is aligned
        read_bp_site = aln_start + l
        flank = 'left'

    elif t == "N" or t == "I" or t == "S" or t == "H": # end of read is aligned
        read_bp_site = aln_start
        flank = 'right'

    else: # reference skip or soft/hardclip "~", or match =
        print("UNEXPECTED!", t)
        sys.exit()

    return read_bp_site, flank


def overlap(q_a, q_b, p_a, p_b):
    # print(q_a, q_b, p_a, p_b)
    assert q_a < q_b and p_a < p_b
    return  (p_a <= q_a <= p_b) or (p_a <= q_b <= p_b) or (q_a <= p_a <= q_b) or (q_a <= p_b <= q_b)


def get_statistics(true_sites, aligned_sites):
    nr_exact = 0
    smallest_exact = 2**32

    nr_approx = 0
    smallest_approx = 2**32

    for t_start, t_stop in true_sites: # ordered smallest to largest
        for a_start, a_stop in aligned_sites:
            if (t_start == a_start) and (t_stop == a_stop):
                nr_exact += 1
                if t_stop - t_start < smallest_exact:
                    smallest_exact = t_stop - t_start
                
            if overlap(a_start, a_stop, t_start, t_stop):
                nr_approx += 1
                if t_stop - t_start < smallest_approx:
                    smallest_approx = t_stop - t_start               
                # break
    if smallest_exact == 2**32:
        smallest_exact = '-'
    if smallest_approx == 2**32:
        smallest_approx = '-'

    print("{0}\t{1}\t{2}\t{3}".format(nr_exact, nr_approx, smallest_exact, smallest_approx))
    # print('Number of exact exon alignments:', nr_exact)
    # print('Smallest exon with an exact exon alignment:', smallest_exact)
    # print('Number of exact OR approximate exon alignments (alnmt overlaps true exon):', nr_approx)
    # print('Smallest exon with an exact OR approximate exon alignment:', smallest_approx)


class Read:
    def __init__(self, info):
        self.read_acc = info[0]
        self.ref_acc = info[2][:-2]
        self.cigarstring = info[5]
        self.reference_start = int(info[3])
        self.m = int(self.read_acc.split("_")[-1]) 

class Statistics:
    def __init__(self, true_bp):
        self.stats = {}
        self.true_bp = true_bp

    def add(self, m, flank, aln_bp_site):
        if m not in self.stats:
            self.stats[m] = {}
            self.stats[m]['left'] = [0,0] # exact, approx
            self.stats[m]['right'] = [0,0] # exact, approx
        

        if self.true_bp == aln_bp_site:
            self.stats[m][flank][0] += 1
        elif max( (self.true_bp - aln_bp_site), (aln_bp_site - self.true_bp) ) < 5:
            self.stats[m][flank][1] += 1


def main(args):
    SAM_bp = args.B + 1 # convert to 1 indexed
    S = Statistics(SAM_bp)
    for line in open(args.samfile, 'r'):
        if line[0] == "@":
            continue
        else:
            info = line.split()
            read = Read(info)
            # print(info)
            if int(info[1]) == 0 or int(info[1])  == 16: #primary
                aln_bp_site, flank = get_breakpoint_site(read, SAM_bp)
                S.add(read.m, flank, aln_bp_site)
            else: #Suppl
                aln_bp_site, flank = get_breakpoint_site(read, SAM_bp)
                S.add(read.m, flank, aln_bp_site)

    # print stats:
    print("m\te_l\te_r\ta_l\ta_r")
    for m in sorted(S.stats.keys()):
        print(m, S.stats[m]['left'][0], S.stats[m]['right'][0],S.stats[m]['left'][1], S.stats[m]['right'][1])





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Simulate references", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('samfile', type=str, help='alignment file')
    parser.add_argument('--B', type=int, default = 500, help='True breakpoint location (0-indexed).')
    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)