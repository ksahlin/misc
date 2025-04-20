
import random
import argparse
import sys
import re
import os
import errno
import shutil
import glob
import math

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")

import seaborn as sns
import pandas as pd
from matplotlib import pyplot


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
            self.stats[m]['left'] = [0, 0, 0] # exact, approx, error
            self.stats[m]['right'] = [0, 0, 0] # exact, approx, error
        

        if self.true_bp == aln_bp_site:
            self.stats[m][flank][0] += 1
        elif max( (self.true_bp - aln_bp_site), (aln_bp_site - self.true_bp) ) < 5:
            self.stats[m][flank][1] += 1
        else:
            self.stats[m][flank][2] += 1


def plot(infile):
    linewidth = 2.5
    matplotlib.rcParams.update({'font.size': 18})
    # sns.set(font_scale=1.5)
    sns.set(font_scale=1.6) 
    sns.set_style("whitegrid")
    indata = pd.read_csv(infile)
    g = sns.relplot(data=indata, x="m", y="y", hue="type", style="flank", linewidth = linewidth, kind="line")

    g.ax.set_title("Splice alignments") #, #dashes = dashes,
        #col="dataset") # hue_order = tools, # hue="datastructure", style="datastructure",
        # col_wrap=3, col_order=["SIM1", "SIM2", "SIM4"], palette=palette)
        # col_order=["SIM3"],  palette=palette)
    # ax = sns.lineplot(data=indata, x="k", y="unique", hue="datastructure", style="chr", palette = sns.color_palette()[:7])
    # axes = g.axes
    # g.set_titles("Genome size 2^16 (repeats)")
    g.set_axis_labels("M (#nt left of breakpoint)", "Fraction per type")
    # g.set(ylim=(94, 99), xticks=[50,75,100,150,200,250,300,500])
    g.set(xticks=[x for x in range(0, 160, 10)]) #[0, 10, 20, 30, 40, 50, 60]
    # g.set_xticks([12, 16, 20, 24, 28])
    g.set_xticklabels(labels=[x for x in range(0, 160, 10)], rotation=60) #rotation=60, 
    # g.tight_layout()
    # g.set(ylim=(95, 100))
    plt.savefig(infile+'.pdf')
    plt.close()


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
    # print("m\te_l\te_r\ta_l\ta_r")
    outfile = open(args.outfile, 'w')
    outfile.write("tool,m,y,flank,type\n")
    N = float(args.N)
    for m in sorted(S.stats.keys()):
        # outfile.write("{0},{1},{2},{3},{4},{5},{6},{7}\n".format(args.tool, m, S.stats[m]['left'][0]/N, S.stats[m]['left'][1]/N, S.stats[m]['left'][2]/N, S.stats[m]['right'][0]/N, S.stats[m]['right'][1]/N ,S.stats[m]['right'][2]/N))
        outfile.write("{0},{1},{2},{3},{4}\n".format(args.tool, m, S.stats[m]['left'][0]/N, 'left', 'exact'))
        outfile.write("{0},{1},{2},{3},{4}\n".format(args.tool, m, S.stats[m]['left'][1]/N, 'left', 'approx (<5)'))
        outfile.write("{0},{1},{2},{3},{4}\n".format(args.tool, m, S.stats[m]['left'][2]/N, 'left', 'error'))
        outfile.write("{0},{1},{2},{3},{4}\n".format(args.tool, m, S.stats[m]['right'][0]/N, 'right', 'exact'))
        outfile.write("{0},{1},{2},{3},{4}\n".format(args.tool, m, S.stats[m]['right'][1]/N, 'right', 'approx (<5)'))
        outfile.write("{0},{1},{2},{3},{4}\n".format(args.tool, m, S.stats[m]['right'][2]/N, 'right', 'error'))

    outfile.close()
    plot(args.outfile)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate split mappings", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('samfile', type=str, help='alignment file')
    parser.add_argument('outfile', type=str, help='output file')
    parser.add_argument('tool', type=str, help='alignment tool')
    parser.add_argument('--B', type=int, default = 500, help='True breakpoint location (0-indexed).')
    parser.add_argument('--N', type=int, default = 100, help='Number of experiments.')
    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)