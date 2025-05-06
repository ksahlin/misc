
import random
import argparse
import sys
# import re
# import os
# import errno
# import shutil
# import glob



# def mkdir_p(path):
#     print("creating", path)
#     try:
#         os.makedirs(path)
#     except OSError as exc:  # Python >2.5
#         if exc.errno == errno.EEXIST and os.path.isdir(path):
#             pass
#         else:
#             raise



# def reverse_complement(string):
#     rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}
#     rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
#     return(rev_comp)

def not_identical(candidate_frag, frags):
    """
    Check if candidate_frag is identical to any of the fragments in frags
    """
    for frag in frags:
        if candidate_frag == frag:
            return False
    return True
    

def main(args):
    frags = [''.join([random.choice('ACGT') for i in range(args.fraglen)])]

    while len(frags) < args.nr_frags:
        frag = random.choice(frags)
        muts = set(random.sample(range(len(frag)),int(len(frag)*args.mut_rate))) # sample fraglen*mut_rate mutation sites
        frag_copy = "".join([frag[i] if i not in muts else  random.choice("ACGT") for i in range(len(frag))])
        if not_identical(frag_copy, frags):
            frags.append(frag_copy)

    frags_out = open(args.outfile, 'w')
    for i, frag in enumerate(frags):
        frags_out.write('>{0}\n{1}\n'.format('frag_'+ str(i), frag))
    frags_out.close()




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Simulate references", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outfile', type=str, help='Output file with fragments')
    parser.add_argument('--mut_rate', type=float, default=0.2, help='mutation rate.')
    parser.add_argument('--fraglen', type=int, default=20, help='segment length.')
    parser.add_argument('--nr_frags', type=int, default=20, help='Number of fragments.')

    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)