# import mmh3
# import screed
# MAX_HASH = 2**64
# import numpy as np

# def hash_kmer(kmer,rc_kmer):
#     canonical_kmer = min(kmer, rc_kmer)
#     hash = mmh3.hash64(canonical_kmer,42,signed=False)[0]
#     return hash

# def sketch_sequence(sequence,kmer_size=21,sampling_factor=1000):
#     hash_threshold = MAX_HASH/sampling_factor
#     seq_len = len(sequence)
#     count = 0
#     hashes = np.zeros(seq_len - kmer_size + 1, dtype=np.uint64)
#     rc_seq = screed.rc(sequence)
#     for i in range(seq_len-kmer_size+1):
#         kmer = sequence[i:i + kmer_size]         # Forward k-mer
#         rc_kmer = rc_seq[-(i + kmer_size):-i or None]  # Reverse complement k-mer
#         hash = hash_kmer(kmer,rc_kmer)
#         if hash < hash_threshold:
#             hashes[count]= hash
#             count +=1
#     return hashes[:count]


# for record in screed.open('rawdata/ecoliMG1655.fa.gz'):
#     print(len(sketch_sequence(record.sequence)))

import mmh3
import screed
import heapq
import numpy as np
import time

start = time.time()
MAX_HASH = 2**64

def hash_kmer(kmer):
    hash_val = mmh3.hash64(kmer, 42, signed=False)[0]
    return hash_val

def sketch_sequence(sequence, kmer_size=21, sampling_factor=1000):
    hash_heap = []  # MinHeap of size at most `sampling_factor`
    rc_seq = screed.rc(sequence)

    seq_len = len(sequence)
    if seq_len < kmer_size:
        return np.array([], dtype=np.uint64)  # Edge case: Sequence too short

    for i in range(seq_len - kmer_size + 1):
        kmer = sequence[i:i + kmer_size]         # Forward k-mer
        rc_kmer = rc_seq[-(i + kmer_size):-i or None]  # Reverse complement k-mer
        hash_val = hash_kmer(min(kmer, rc_kmer))

        if len(hash_heap) < sampling_factor:
            heapq.heappush(hash_heap, -hash_val)  # Push negative value to simulate max heap
        else:
            # Only keep smaller hashes
            if -hash_heap[0] > hash_val:
                heapq.heappushpop(hash_heap, -hash_val)

    # Convert back to positive values and return as a sorted NumPy array
    return np.array(sorted(-h for h in hash_heap), dtype=np.uint64)

for record in screed.open('rawdata/ecoliMG1655.fa.gz'):
    print(sketch_sequence(record.sequence))

end = time.time()
print("The time of execution of above program is :",
      (end-start), "s")