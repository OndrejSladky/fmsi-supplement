#! /usr/bin/env python3

import math
import sys

def GetMinimizerLen4SSHash(k, kmerCount): # using a formula from Giulio
    return min(math.ceil(math.log(kmerCount, 4)) + 1, int(k) - 2)


def main():
    k = int(sys.argv[1])
    #print(k)
    kmer_cnt_file = sys.argv[2]
    with open(kmer_cnt_file) as f:
        kmer_cnt = int(f.readline().strip('\n'))
    #print(kmer_cnt)
    print(GetMinimizerLen4SSHash(k, kmer_cnt))


if __name__ == "__main__":
    main()
