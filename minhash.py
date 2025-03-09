import mmh3
import screed
MAX_HASH = 2**64


def hash_kmer(kmer):
    rc_kmer = screed.rc(kmer)
    canonical_kmer = kmer if kmer<rc_kmer else rc_kmer
    hash = mmh3.hash64(canonical_kmer,42)[0]
    return hash

def sketch_sequence(sequence,kmer_size=21,sampling_factor=1000):
    hash_threshold = MAX_HASH/sampling_factor
    hashes = []
    for i in range(len(sequence)-kmer_size+1):
        hash = hash_kmer(sequence[i:i+kmer_size])
        if hash < 0:
            hash += 2**64
        if hash < hash_threshold:
            hashes.append(hash)
    return hashes


for record in screed.open('rawdata/ecoli_ref-5m.fastq.gz'):
    print(sketch_sequence(record.sequence))