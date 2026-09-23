# variant_clinical_v4_02

Physically consolidated PubMed variant–clinical pipeline. The complete Python surface is exactly 12 readable source files under `scripts/`; the four compatibility shims load sibling bytecode payloads. `scripts/data` remains `../data`.

## Commands

```bash
python scripts/run.py cohort                 # Stage0–9, fixed 117-PMID fixture
python scripts/run.py pipeline               # Stage0–9, full PubMed query Stage1
python scripts/run.py single 42068976 --no-network
python scripts/run.py deliver all
python scripts/run.py qc refalt --in-csv INPUT.csv --out-tsv OUTPUT.tsv
python scripts/run.py bootstrap --import-from /path/to/variant_clinical_tree
python scripts/run.py compare
```

Cohort mode runs preflight, fixture ingest, linking, full-text enrichment and evidence stamping, discovery, ClinVar validation, strict/high-confidence review, journal annotation, text materialization, pathogenicity, negatives, and genomic export. Each heavy stage runs in its own process, so ClinVar and Ensembl indexes are released before the next stage. It never merges v3 outputs. Cache reuse is explicit: `--seed-cache --seed-cache-from PATH`.

Set `VARIANT_CLINICAL_STAGE8_PDF=1` to run the optional strict PDF downloader before text materialization. Stage9 consumes the pathogenicity sidecar produced by Stage8.

## Consolidated source layout

- `run.py`: public commands, Stage0–9, cohort ingest, stamp, bootstrap, comparison, integrated self-tests
- `vc_paths.py`, `vc_text.py`: canonical shared helpers
- `fulltext.py`: `enrich`, `materialize`
- `review.py`: `discovery`, `filter`, `journal`, `pathogenicity`, `negatives`
- `genomic.py`: `export` and genomic helper APIs
- `extras.py`: `single`, `download-pdfs`, `new-site`, `audit`, `html`, `review-export`, `refalt`, `ref-context5`
- `_pyc_loader.py` and four shim sources: `link_variant_phenotype.py`, `pubmed_pipeline.py`, `validate_novelty_clinvar.py`, `vc_ensembl_c_map.py`

Reference data stays under `data/pubmed/ref`; generated ingest, pipeline, cache, run, and log trees are disposable.
