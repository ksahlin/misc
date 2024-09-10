
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
    mkdir_p(args.outfolder)

    genome = [''.join([random.choice('ACGT') for i in range(args.seglen)]) for j in range(args.exon_max)]

    query = []
    for i in range(args.exon_min, args.exon_max):
        query.append(genome[i][:i]) # takes out subsegment of size i from segment i in genome 
        if args.canonical:
            genome[i] = genome[i][:i] + 'GT' + genome[i][i+2 : args.seglen-2] + 'AG'

    if args.repeats > 0:
        original = ''.join([s for s in genome])
        for j in range(args.repeats):
            muts = set(random.sample(range(len(original)),int(len(original)*args.repeat_mutrate)))
            copy = "".join([original[i] if i not in muts else random.choice(['', random.choice("ACGT"), original[i] + random.choice("ACGT")]) for i in range(len(original))])
            genome.append(copy)

    genome_out = open(os.path.join(args.outfolder,'genome.fa'), 'w')
    genome_out.write('>{0}\n{1}'.format('genome',''.join([s for s in genome])))
    genome_out.close()

    if args.error_rate > 0: # introduce errors
        read = ''.join([s for s in query])
        muts = set(random.sample(range(len(read)),int(len(read)*args.error_rate)))
        read_mod = "".join([read[i] if i not in muts else random.choice(['', random.choice("ACGT"), read[i] + random.choice("ACGT")]) for i in range(len(read))])
        read_out = open(os.path.join(args.outfolder,'query.fa'), 'w')
        read_out.write('>{0}\n{1}'.format('read',read_mod))
        read_out.close()
    else:
        gene_out = open(os.path.join(args.outfolder,'query.fa'), 'w')
        gene_out.write('>{0}\n{1}'.format('gene',''.join([s for s in query])))
        gene_out.close()


    # construct gtf
    gtf_out = open(os.path.join(args.outfolder,'annotation.gtf'), 'w')
    gtf_out.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format('genome','sim','gene', 1, 1+args.exon_max*args.seglen, '.', '+', '.', 'gene_id "ENSG00000223972";' ) )
    gtf_out.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format('genome','sim','transcript', 1, 1+args.exon_max*args.seglen, '.', '+', '.', 'gene_id "ENSG00000223972";' ) )
    for i in range(args.exon_min, args.exon_max):
        start = 1 + i * args.seglen
        stop  = 1 + i * args.seglen + i
        gtf_out.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format('genome','sim','exon', start, stop, '.', '+', '.', 'gene_id "ENSG00000223972";' ) )
    gtf_out.close()


    # simulate reads



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Simulate references", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--outfolder', type=str, help='Output file with references')
    parser.add_argument('--mut_rate', type=float, default=0.01, help='error rate.')
    parser.add_argument('--seglen', type=int, default=100, help='segment length.')
    parser.add_argument('--exon_min', type=int, default=1, help='exon_min.')
    parser.add_argument('--exon_max', type=int, default=100, help='exon_max.')
    parser.add_argument('--error_rate', type=float, default=0.0, help='Number of reads.')
    parser.add_argument('--canonical', action= 'store_true', help='Only simulate canonical spice sites (GT-AG).')
    parser.add_argument('--repeats', type=int, default=0, help='Number of repeat copies of gene.')
    parser.add_argument('--repeat_mutrate', type=float, default=0.025, help='Number of reads.')

    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)