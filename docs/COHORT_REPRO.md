# Gold cohort reproduction (v4_02)

## Input

- `configs/gold_pmids.txt`: 117 unique PMIDs.
- `fixtures/cohort/articles.jsonl`: fixed Stage1 article metadata. Cohort mode does not run the full PubMed query sweep.
- `configs/gold_standard_v3.csv`: genomic comparison reference (129 rows; 126 with complete chr/pos/ref/alt).

## Command

```bash
python scripts/run.py bootstrap --link-from /path/to/existing/tree
python scripts/run.py cohort --force
python scripts/run.py deliver all
python scripts/run.py compare
```

`compare` reads `data/pubmed/pipeline/s9_genomic/variants_genomic_positive_ok.csv` and exits nonzero when any complete gold allele is absent or has different coordinates. OA full text is fetched during Stage3 and requires network access.

## Fresh v4_02 run (2026-09-23)

Stage0–9 completed. Stage1 contained 117 articles. Comparison against the gold CSV:

| Metric | Value |
|--------|------:|
| Gold rows with four genomic columns | 126 |
| Positive rows with resolved coordinates | 163 |
| Exact genomic matches | 120 |
| Match rate | 0.9524 |
| Coordinate conflicts among matches | 0 |

The six unmatched gold rows are absent from the positive resolved export; matched rows have identical chromosome, position, ref, and alt.
