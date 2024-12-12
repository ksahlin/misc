import math
import random
import sys
from collections import Counter

def current_method(seq, k, w_min, w_max, v1, v2, BIT_SPACE):
    hashes = []
    main_hashes = []
    mask = 2**BIT_SPACE - 1
    for i in range(len(seq)):
        h1 = hash(seq[i:i+k])
        min_hash = 2**64
        for j in range(i+k+w_min, i+k+w_max):
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

def new_method(seq, k, w_min, w_max, v1, v2, BIT_SPACE):
    hashes = []
    main_hashes = []
    mask = 2**BIT_SPACE - 1
    for i in range(len(seq)):
        h1 = hash(seq[i:i+k])
        min_hash = 2**64
        for j in range(i+k+w_min, i+k+w_max):
            h2 = hash(seq[j:j+k])
            if h1 ^ h2 < min_hash:
                min_hash =  h1 ^ h2
                h2_min = h2

        h1_is_odd = h1 % 2 == 1
        h2_is_odd = h2_min % 2 == 1
        if h1_is_odd and h1_is_odd != h2_is_odd:
            main_h = h1
            aux_h = h2
        elif h2_is_odd and h2_is_odd != h1_is_odd:
            main_h = h2
            aux_h = h1
        elif h1_is_odd == h2_is_odd:
            main_h = min(mask & h1, mask & h2_min)
            aux_h = max(mask & h1, mask & h2_min)

        mask_main = 2**v1-1 << (BIT_SPACE - v1)
        mask_aux = 2**v2-1

        final_hash = (main_h & mask_main) ^ (aux_h & mask_aux)
        hashes.append(final_hash)
        main_hashes.append((main_h & mask_main))
    return hashes, main_hashes



def construct_ref(ref_size, key='random'): 
    # key = 0 fully random
    # key = 1 repeats

    if key == 'random':
        # fully random
        seq = ''.join([random.choice('ACGT') for i in range(ref_size)])
    elif key == 'repeats':
        # repeats
        origin = ''.join([random.choice('ACGT') for i in range(2**9)])
        copies = 1
        nr_mut = int(len(origin)*0.05) # 5% mut
        seq_tot = [origin]
        while len(origin)*copies < ref_size:
            muts = set(random.sample(range(len(origin)),nr_mut))
            cpy = "".join([origin[i] if i not in muts else random.choice("ACGT") for i in range(len(origin))])
            copies += 1
            seq_tot.append(cpy)
        seq = ''.join([s for s in seq_tot])
    else:
        print('Error argument')
        sys.exit()

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
    G = 2**20
    # w = 2**6
    # k=8
    # BIT_SPACE=32

    seq = construct_ref(G, 'random')
    # seq = construct_ref(1, G)
    print('G,w,k,B,E,Exp_B,Pred_B')
    for k in [20]:
        for w_max in [2**6]:
            for BIT_SPACE in [64]:
                # exp = []
                for v1 in [32]: # range(BIT_SPACE//2, BIT_SPACE+1): # iterate over several different bit allocations
                    # print()
                    # print('V1 = ', v1)
                    v2 = BIT_SPACE - v1
                    w_min = w_max - 2
                    hashes, main_hashes = current_method(seq, k, w_min, w_max, v1, v2, BIT_SPACE)
                    # m = len(hashes) + len(main_hashes)
                    total_main = len(main_hashes)
                    distinct_main = len(set(main_hashes))
                    e_h = e_hits(main_hashes)
                    # exp.append( (round(e_h,2), round(v1/BIT_SPACE, 2)) )
                    print('CURRENT METHOD: Total: {0}, distinct: {1}, E-hits: {2}'.format(total_main, distinct_main, round(e_h,2)))

                    hashes, main_hashes = new_method(seq, k, w_min, w_max, v1, v2, BIT_SPACE)
                    # m = len(hashes) + len(main_hashes)
                    total_main = len(main_hashes)
                    distinct_main = len(set(main_hashes))
                    e_h = e_hits(main_hashes)
                    # exp.append( (round(e_h,2), round(v1/BIT_SPACE, 2)) )
                    print('NEW METHOD: Total: {0}, distinct: {1}, E-hits: {2}'.format(total_main, distinct_main, round(e_h,2)))

                    # print('E hits:', e_h)
                    # print('tot hashes:', m)
                    # n = len(set(hashes) | set(main_hashes) )
                    # print('Distinct hashes:', n)
                    # print('Collisions: ', m - n)
                # best = sorted(exp)[0]
                # pred_B = round(math.log(G,2)/( math.log(G,2) + math.log(w,2) ), 2)
                # print("{0},{1},{2},{3},{4},{5},{6}".format(len(seq),w,k,BIT_SPACE,best[0],best[1],pred_B) )

main()