ATTACH 'sequences.db' as seq;
ATTACH 'ref_seq.db' as ref;

with seq_subsample as (
    SELECT DISTINCT hash_value FROM (
        SELECT * FROM seq.hash_table ORDER BY hash_value LIMIT 10000
    )
),
ref_subsample as (
    SELECT DISTINCT hash_value FROM (
        SELECT * FROM ref.hash_table ORDER BY hash_value LIMIT 10000
    )
),
intersection AS (
    SELECT COUNT(*) as intersection_size
    FROM seq_subsample
    JOIN ref_subsample ON ref_subsample.hash_value = seq_subsample.hash_value
),
seq_size AS (
    SELECT COUNT(*) as seq_count
    FROM seq_subsample
),
ref_size AS (
    SELECT COUNT(*) as ref_count
    FROM ref_subsample
),
union_size AS (
    SELECT (seq_count + ref_count - intersection_size) as union_size
    FROM seq_size, ref_size, intersection
)

SELECT 
    intersection_size,
    seq_count,
    ref_count,
    union_size,
    intersection_size / seq_count AS containment_seq_in_ref,
    intersection_size / ref_count AS containment_ref_in_seq,
    intersection_size / union_size AS jaccard_similarity
FROM intersection, seq_size, ref_size, union_size;