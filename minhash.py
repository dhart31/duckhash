import mmh3
import screed
import heapq
import numpy as np

tab = str.maketrans("ACTG", "TGAC")
def reverse_complement(seq):
    return seq.translate(tab)[::-1]

MAX_HASH = 2**64

def hash_kmer(kmer):
    hash_val = mmh3.hash64(kmer, 42, signed=False)[0]
    return hash_val

def sketch_sequence(sequence, kmer_size=21, sampling_factor=1000):
    hash_heap = []
    rc_seq = reverse_complement(sequence)

    seq_len = len(sequence)
    if seq_len < kmer_size:
        return np.array([], dtype=np.uint64)

    for i in range(seq_len - kmer_size + 1):
        hash_val = hash_kmer(
            min(
                sequence[i:i + kmer_size],
                rc_seq[-(i + kmer_size):-i or None]
                )
            )

        if len(hash_heap) < sampling_factor:
            heapq.heappush(hash_heap, -hash_val)
        else:
            if -hash_heap[0] > hash_val:
                heapq.heappushpop(hash_heap, -hash_val)

    return np.array(sorted(-h for h in hash_heap), dtype=np.uint64)

for record in screed.open('rawdata/ecoliMG1655.fa.gz'):
    print(sketch_sequence(record.sequence))