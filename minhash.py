"""
Author: Derek Hart
Date: 03/10/2025
Creates sequence sketch of sequences
"""
import heapq
import mmh3
import screed
import duckdb
import typing
from duckdb.typing import *
import numpy as np

def reverse_complement(seq:str) -> str:
    """
    Returns reverse complement string of DNA sequence
    """
    tab = str.maketrans("ACTG", "TGAC")
    return seq.translate(tab)[::-1]

MAX_HASH = 2**64

def hash_kmer(kmer:str) -> int:
    """
    Calculates a 64-bit unsigned hash for input kmer
    """
    hash_val = mmh3.hash64(kmer, 42, signed=False)[0]
    return hash_val

def sketch_sequence(sequence:str, kmer_size:int = 21, max_hashes:int =10000) -> np.array:
    """
    Creates a representative sketch of sequence with hash array
    - keeps smallest hash values
    - array is can be no larger than max_hashes
    - returns numpy array
    """
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

        if len(hash_heap) < max_hashes:
            heapq.heappush(hash_heap, -hash_val)
        else:
            if -hash_heap[0] > hash_val:
                heapq.heappushpop(hash_heap, -hash_val)

    return np.array(sorted(-h for h in hash_heap), dtype=np.uint64)

con = duckdb.connect('sketch_results.db')
for record in screed.open('rawdata/ecoliMG1655.fa.gz'):
    con.create_function(
        "sketch_sequence",
        sketch_sequence,
        [str,int,int],
        list[int])
    con.execute("CREATE OR REPLACE TABLE hash_table (sequence VARCHAR, hash_value BIGINT)")
    con.execute(f"INSERT INTO hash_table SELECT '{record.name}',UNNEST(sketch_sequence('{record.sequence}',21,10000))")
con.commit()
con.close()
