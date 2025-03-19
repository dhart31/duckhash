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
import argparse
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

def make_hash_array(sequence:str, kmer_size:int = 21, max_hashes:int =10000) -> np.array:
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

def sketch_sequence(filename,db_name):
    con = duckdb.connect(db_name)
    con.create_function(
        "make_hash_array",
        make_hash_array,
        [str,int,int],
        list[HUGEINT])
    
    sql_create = """
    CREATE OR REPLACE TABLE hash_table (
        sequence VARCHAR, 
        hash_value HUGEINT
    )
    """
    con.execute(sql_create)
    with screed.open(filename) as seqfile:
            for record in seqfile:
                sql_insert = f"""
                INSERT INTO hash_table 
                SELECT 
                    '{record.name}', 
                    UNNEST(make_hash_array('{record.sequence}', 21, 1000))
                """
                con.execute(sql_insert)
            
    con.commit()
    con.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sequence','-s', type=str, help="Full path to a fasta/fastq file.")
    parser.add_argument('--database','-db', type=str, help="Name of output database")

    args=parser.parse_args()
    sketch_sequence(args.sequence,args.database)


