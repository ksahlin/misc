
import random
import argparse
import sys
import re
import os
import errno
import shutil
import glob



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


def main(args):
    """
        Simulate a contig pair (c1, c2). Each contig is of size 1000nt 
        and has a site at B (B=500, 1-indexed for gtf), where a read is split mapped over. 
        The first m nucleotides of the read maps to c1, 
        and the last R-m nucleotides map on c2. We pick all m in interval [10, R-10], 
        resulting in R-10 reads over the same region with a unique splicing. 
        
        The above described experiment is repeated N times (default = 100),
        to account for randomness in syncmer/minimizer generation around the region.

        While it is a simple experiment, it should serve as a first test dataset to
        investigate splice mapping performance.
        
    """
    mkdir_p(args.outfolder)
    B = 500
    # simulate N "contig" pairs 

    ctg_pairs = []
    for i in range(args.N):
        c1 = ''.join([random.choice('ACGT') for i in range(args.seglen)])
        c2 = ''.join([random.choice('ACGT') for i in range(args.seglen)])
        ctg_pairs.append((c1,c2))


    # output genome
    genome_out = open(os.path.join(args.outfolder,'genome.fa'), 'w')
    for i, (c1, c2) in enumerate(ctg_pairs):
        genome_out.write('>{0}\n{1}\n'.format('ctg_{0}_A'.format(i),c1))
        genome_out.write('>{0}\n{1}\n'.format('ctg_{0}_B'.format(i),c2))
    genome_out.close()

    
    reads_out = open(os.path.join(args.outfolder,'reads.fa'), 'w')
    for i in range(args.N):
        # simulate reads
        for m in range(10, args.R-10+1):
            c1, c2 = ctg_pairs[i]
            read = c1[B - m: B] # up to B-1; 499 (0 indexed here)
            read += c2[B: B + args.R-m] # start at B; 500 (0 indexed here)
            if args.error_rate > 0:
                raise('To be implemented - simulate errors here according to error rate')
            reads_out.write('>{0}\n{1}\n'.format('read_cgt_{0}_m_{1}'.format(i, m), read))
    reads_out.close()



    # # construct gtf
    # gtf_out = open(os.path.join(args.outfolder,'annotation.gtf'), 'w')
    # gtf_out.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format('genome','sim','gene', 1, 1+args.exon_max*args.seglen, '.', '+', '.', 'gene_id "ENSG00000223972";' ) )
    # gtf_out.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format('genome','sim','transcript', 1, 1+args.exon_max*args.seglen, '.', '+', '.', 'gene_id "ENSG00000223972";' ) )
    # for i in range(args.exon_min, args.exon_max):
    #     start = 1 + i * args.seglen
    #     stop  = 1 + i * args.seglen + i - 1 # end coordinate is inclusive 
    #     gtf_out.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format('genome','sim','exon', start, stop, '.', '+', '.', 'gene_id "ENSG00000223972";' ) )
    # gtf_out.close()





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Simulate references", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outfolder', type=str, help='Output file with references')
    parser.add_argument('--error_rate', type=float, default=0.0, help='Read error rate.')
    parser.add_argument('--seglen', type=int, default=1000, help='segment length.')
    parser.add_argument('--R', type=int, default=150, help='Read length.')
    # parser.add_argument('--canonical', action= 'store_true', help='Only simulate canonical spice sites (GT-AG).')
    parser.add_argument('--N', type=int, default=100, help='Number of experiments.')
    # parser.add_argument('--repeat_mutrate', type=float, default=0.025, help='Number of reads.')

    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)