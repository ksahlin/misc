
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


def get_true_exon_sites(gtf_file):
    true_sites = []
    for line in gtf_file:
        info = line.split()
        feature, start, stop  = info[2], int(info[3]), int(info[4])
        if feature == 'exon':
            true_sites.append((start, stop+1)) # converting back to exclusive end coordinate

    annotated_splice_junctions = []
    gtf_file.seek(0)
    all_lines = gtf_file.readlines()
    for line1, line2 in zip(all_lines[:-1], all_lines[1:]):
        info = line1.split()
        feature1, start1, stop1  = info[2], int(info[3]), int(info[4])
        info = line2.split()
        feature2, start2, stop2  = info[2], int(info[3]), int(info[4])
        if feature1 == 'exon' and feature2 == 'exon':
            annotated_splice_junctions.append((stop1, start2))
    # print(true_sites)
    # print(annotated_splice_junctions)
    return true_sites, annotated_splice_junctions


def get_exon_sites(cigar_tuples, first_exon_start, annotated_splice_junctions):
    ref_pos = first_exon_start
    exon_sites = [ref_pos]
    # print(cigar_tuples)
    for i, (l,t) in enumerate(cigar_tuples):
        if t == "D":
            if (ref_pos, ref_pos + l) in annotated_splice_junctions:
                exon_sites.append( ref_pos )
                exon_sites.append( ref_pos + l )
                # print("HEERE")
            ref_pos += l

        elif t == "=" or t== "M" or t == "X":
            ref_pos += l

        elif t == "N":
            splice_start = ref_pos
            ref_pos += l
            splice_stop = ref_pos
            exon_sites.append( splice_start)
            exon_sites.append( splice_stop )

        elif t == "I" or t == "S" or t == "H": # insertion or softclip
            ref_pos += 0

        else: # reference skip or soft/hardclip "~", or match =
            print("UNEXPECTED!", t)
            sys.exit()

    exon_sites.append(ref_pos)
    exon_sites = [(exon_sites[i],exon_sites[i+1]) for i in range(0, len(exon_sites), 2) ]  
    return exon_sites

def get_aligned_exon_sites(read, annotated_splice_junctions):
    read_cigar_tuples = []
    result = re.split(r'[=DXSMINH]+', read.cigarstring)
    i = 0
    for length in result[:-1]:
        i += len(length)
        type_ = read.cigarstring[i]
        i += 1
        read_cigar_tuples.append((int(length), type_ ))  

    read_exon_sites = get_exon_sites(read_cigar_tuples, read.reference_start, annotated_splice_junctions)

    return read_exon_sites


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
        
        self.cigarstring = info[5]
        self.reference_start = int(info[3])

def main(args):
    true_sites, annotated_splice_junctions = get_true_exon_sites(open(args.gtf, 'r'))

    # get mapping
    # mappings = {}
    for line in open(args.samfile, 'r'):
        if line[0] == "@":
            continue
        else:
            info = line.split()
            # print(info)
            if int(info[1]) == 0 or int(info[1])  == 16:
                read = Read(info)
                aligned_sites = get_aligned_exon_sites(read, annotated_splice_junctions)
                # print(read.reference_start)
                # print(aligned_sites)
                get_statistics(true_sites, aligned_sites)

    # SAM_file = pysam.AlignmentFile(args.samfile, "r", check_sq=False)
    # for read in SAM_file.fetch(until_eof=True): # TODO: remove pysam and parse manually
    #     if read.flag == 0 or read.flag == 16:
    #         aligned_sites = get_aligned_exon_sites(read, annotated_splice_junctions)
    #         # print(aligned_sites)
    #         print(read.reference_start)
    #         print(aligned_sites)
    #         get_statistics(true_sites, aligned_sites)







if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Simulate references", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--samfile', type=str, help='alignment file')
    parser.add_argument('--gtf', type=str, help='True exon annotation.')
    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)