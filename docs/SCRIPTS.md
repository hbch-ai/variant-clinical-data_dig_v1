# Python consolidation map

The repository intentionally contains exactly 12 `.py` files, all directly under `scripts/`.

- `run.py`: `cohort`, `pipeline`, `single`, `deliver`, `qc`, `bootstrap`, `compare`; Stage0–9 orchestration and integrated support logic.
- `fulltext.py`: `enrich`, `materialize`.
- `review.py`: `discovery`, `filter`, `journal`, `pathogenicity`, `negatives`.
- `genomic.py`: `export` plus public genomic helpers used by `extras.py`.
- `extras.py`: `single`, `download-pdfs`, `new-site`, `audit`, `html`, `review-export`, `refalt`, `ref-context5`.
- `vc_paths.py`, `vc_text.py`: shared path and text APIs.
- `_pyc_loader.py`: sibling payload loader with legacy path remapping.
- `link_variant_phenotype.py`, `pubmed_pipeline.py`, `validate_novelty_clinvar.py`, `vc_ensembl_c_map.py`: source shims backed by sibling `.pyc` files.

Merged implementations use section-specific global prefixes and explicit `fulltext_module`, `review_module`, and `genomic_module` imports. They do not contain archived or string-embedded source.
