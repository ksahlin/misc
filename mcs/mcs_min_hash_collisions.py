import math
import random
import sys
from collections import Counter

def compute_hashes(seq, k, w, v1, v2, BIT_SPACE):
    hashes = []
    main_hashes = []
    mask = 2**BIT_SPACE - 1
    for i in range(len(seq)):
        h1 = hash(seq[i:i+k])
        min_hash = 2**64
        for j in range(i+k, i+k+w):
            h2 = hash(seq[j:j+k])
            if h1 ^ h2 < min_hash:
                min_hash =  h1 ^ h2
                h2_min = h2

        main_h = min(mask & h1, mask & h2_min)
        aux_h = max(mask & h1, mask & h2_min)
        mask_main = 2**v1-1 << (BIT_SPACE - v1)
        mask_aux = 2**v2-1

        final_hash = (main_h & mask_main) ^ (aux_h & mask_aux)
        hashes.append(final_hash)
        main_hashes.append((main_h & mask_main))
    return hashes, main_hashes

def construct_ref(key, ref_size): 
    # key = 0 fully random
    # key = 1 repeats

    if key == 0:
        # fully random
        seq = ''.join([random.choice('ACGT') for i in range(ref_size)])
    else:
        # repeats
        origin = ''.join([random.choice('ACGT') for i in range(2**9)])
        copies = 0
        nr_mut = int(len(origin)*0.05) # 5% mut
        while len(origin)*copies < ref_size:
            seq_tot = [origin]
            while copies < 375:
                muts = set(random.sample(range(len(origin)),nr_mut))
                cpy = "".join([origin[i] if i not in muts else random.choice("ACGT") for i in range(len(origin))])
                copies += 1
                seq_tot.append(cpy)
            seq = ''.join([s for s in seq_tot])

    # print(len(seq))

    return seq

def e_hits(l):
    s = len(l)
    C = Counter(l)
    e_hits = 0
    for seed in C:
        count = C[seed]
        e_hits += count**2

    e_hits = e_hits / s
    return e_hits


def main():
    """
        G - genome size
        w - window size
        v1 - nr bits allocated for strobe1
        v2 - nr bits allocated for strobe2
        BIT_SPACE - total avialable bits (we have 64 in reality)
    """
    G = 2**16
    # w = 2**6
    # k=8
    # BIT_SPACE=32

    #seq = construct_ref(0, G)
    seq = construct_ref(1, G)
    print('G,w,k,B,E,Exp_B,Pred_B')
    for k in [4,6,8]:
        for w in [2**2, 2**4, 2**6]:
            for BIT_SPACE in [12, 16, 20, 24, 28]:
                exp = []
                for v1 in range(BIT_SPACE//2, BIT_SPACE+1): # iterate over several different bit allocations
                    # print()
                    # print('V1 = ', v1)
                    v2 = BIT_SPACE - v1
                    hashes, main_hashes = compute_hashes(seq, k, w, v1, v2, BIT_SPACE)
                    # m = len(hashes) + len(main_hashes)
                    e_h = e_hits(hashes+main_hashes)
                    exp.append( (round(e_h,2), round(v1/BIT_SPACE, 2)) )
                    # print('E hits:', e_h)
                    # print('tot hashes:', m)
                    # n = len(set(hashes) | set(main_hashes) )
                    # print('Distinct hashes:', n)
                    # print('Collisions: ', m - n)
                best = sorted(exp)[0]
                pred_B = round(math.log(G,2)/( math.log(G,2) + math.log(w,2) ), 2)
                print("{0},{1},{2},{3},{4},{5},{6}".format(G,w,k,BIT_SPACE,best[0],best[1],pred_B) )

main()