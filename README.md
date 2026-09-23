# variant_clinical_data_dig_v1

The pipeline finds variant–phenotype links in 2026 literature, keeps candidates that are absent from a pinned ClinVar snapshot and reviewable in the text, and resolves them to GRCh38 `chr / pos / ref / alt`. The published run does not query PubMed again. Stage1 reads a fixed list of 117 PMIDs. The runtime is CPython 3.13. There is no language model. Stages 2–9 are rules, regular expressions, and index lookups.

A positive record passes four checks: a variant and a clinical entity co-occur in the abstract or full text; the exact ClinVar key is absent; delivery drops a bare rsID that has no paired c.HGVS in the text; coordinates are checked against the reference FASTA and annotated with a GFF3 region. `screen_grade` and discovery labels are ranking and annotation fields. ClinVar-absent does not mean the variant is the first report in the literature. Negative candidates are extracted on a side track and are not written back into the positive funnel.

## 1. Layout and entry points

The public entry point is `scripts/run.py`. Stage2, Stage5, PubMed retrieval, and the Ensembl c. mapping load sibling `.pyc` files through same-named `.py` shims. `scripts/data` points at `../data`.

### 1.1 Repository tree

```text
variant_clinical_v1/
├── configs/
│   ├── pipeline_gates.json
│   ├── pipeline_gate_decisions.example.json
│   ├── journal_priority.json
│   ├── queries.json                 # Stage1 full-corpus queries
│   ├── repro.json
│   ├── gold_pmids.txt               # 117 PMIDs for cohort
│   ├── gold_standard_v3.csv         # gold table used by compare
│   └── gold_cohort_manifest.json
├── fixtures/cohort/articles.jsonl   # fixed cohort records
├── scripts/
│   ├── run.py                       # cohort / pipeline / single / deliver / qc / compare
│   ├── fulltext.py                  # Stage3 enrichment; Stage8 text materialization
│   ├── review.py                    # discovery, review, journal, pathogenicity, negatives
│   ├── genomic.py                   # Stage9 genomic export
│   ├── extras.py                    # single PMID, PDF, new-site, review HTML, QC
│   ├── vc_text.py                   # variant/clinical spans, gene windows
│   ├── vc_paths.py                  # data/pubmed path constants
│   ├── link_variant_phenotype.py    # Stage2 (loads the sibling .pyc)
│   ├── pubmed_pipeline.py           # full-corpus Stage1 (loads the sibling .pyc)
│   ├── validate_novelty_clinvar.py  # Stage5 (loads the sibling .pyc)
│   └── vc_ensembl_c_map.py          # Stage9 c. mapping (loads the sibling .pyc)
├── data/pubmed/
│   ├── ref/                         # ClinVar snapshot + Ensembl GFF3/FASTA/FAI (not in Git)
│   │   ├── clinvar/
│   │   └── ensembl/
│   ├── ingest/                      # Stage1 merged records
│   ├── pipeline/                    # Stage2–9 outputs
│   │   ├── s2_linking/
│   │   ├── s4_discovery/
│   │   ├── s5_clinvar/{strict,fulltext,enriched}/
│   │   ├── s6_review/
│   │   ├── s8_negatives/negatives_independent/
│   │   └── s9_genomic/{,new_site,manual_review,audit_trails}/
│   ├── cache/                       # full text, text cache, dbSNP, Ensembl indexes
│   └── runs/                        # cohort report / single-PMID sandbox / compare
└── docs/
```

### 1.2 Path prefixes

Common prefixes are `ingest/merged/`, `pipeline/s2_linking/`, `pipeline/s4_discovery/`, `pipeline/s5_clinvar/enriched/` (positive track), `pipeline/s6_review/`, `cache/fulltext/`, `cache/pdfs/`, `cache/dbsnp/`, `cache/ensembl/`, `ref/ensembl/`, `pipeline/s9_genomic/`, `s9_genomic/new_site/`, `s9_genomic/manual_review/`, and `s9_genomic/audit_trails/`. Literature negatives stay in `pipeline/s8_negatives/negatives_independent/`.

A single-PMID run writes under `data/pubmed/runs/single/<PMID>/`. Inside that sandbox the names are still `merged/`, `linked/`, and `fulltext/`. Those files do not replace the production trees.

### 1.3 Output classes and text tracks

Outputs fall into four classes: **primary** (the default input of the next positive step), **side** (routed-away or negative rows that can be replayed), **mirror** (the same content copied to another directory), and **sidecar** (TSV, counts, HTML). `pairs*.jsonl` is one variant–clinical link. `corpus*.jsonl` is one article-level row. The positive funnel follows the enriched corpus.

Stages 2–5 still keep three text tracks. `strict` is abstract-only. `fulltext` is new full-text links only. `enriched` is abstract plus full text, and it is the positive track.

Ensembl files live inside the repository. The ClinVar 2026 benign/likely-benign negative directory from the earlier tree is not part of this repository. The former per-stage scripts are the 12 Python modules listed above.

## 2. How to run a reproduction

After clone, the reference files are not in Git. Place the ClinVar `variant_summary.txt.gz` and the Ensembl GRCh38.115 GFF3, FASTA, and FAI under `data/pubmed/ref/`, then run the fixed article list. `cohort` runs Stage0–9 in order. It does not read a v3 directory and does not merge a historical high-confidence file.

### 2.1 Commands

```bash
python3.13 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
python scripts/run.py bootstrap --import-from /path/to/tree
python scripts/run.py cohort --force
python scripts/run.py deliver all
python scripts/run.py compare
```

`bootstrap` copies or hardlinks the reference data into the project. `cohort` is the published reproduction entry. `deliver` applies the new-site filter and writes the review pages. `compare` aligns positive rows that have complete coordinates with the 126 complete rows in `configs/gold_standard_v3.csv`.

### 2.2 Other commands

`python scripts/run.py pipeline` still runs the full PubMed search. Stage1 then uses the seven 2026 queries in `configs/queries.json`. `python scripts/run.py single PMID` runs one article in a sandbox and leaves the production corpus unchanged.

## 3. Pipeline diagram

The main path follows the enriched corpus into the high-confidence set, pathogenicity labels, genomic positives, and the new-site filter. Journal grade does not delete rows. Independent negatives are drawn on the right and are not part of the gold comparison.

### 3.1 Stage flow

```text
fixtures/cohort/articles.jsonl + configs/gold_pmids.txt
        |
        v
+---------------------------+
| S0  preflight             |  configs and ClinVar / GFF3 / FASTA must exist
+---------------------------+
        |
        v
+---------------------------+
| S1  pubmed_retrieval      |  cohort: 117 fixed articles
+---------------------------+     pipeline: seven 2026 queries
        |
        |  ingest/merged/articles.jsonl
        v
+---------------------------+
| S2  abstract_linking      |  a specific variant and a clinical signal
+---------------------------+
        |
        |  pipeline/s2_linking/articles_pass.jsonl
        v
+---------------------------+
| S3  fulltext_enrichment   |  OA full text when available; otherwise keep the abstract
+---------------------------+
        |
        |  pipeline/s2_linking/corpus_enriched.jsonl
        v
+---------------------------+
| S4  discovery_tagging     |  novel / expansion; known and unclear rows stay aside
+---------------------------+
        |
        |  pipeline/s4_discovery/corpus_enriched_novel.jsonl
        v
+---------------------------+
| S5  clinvar_validation    |  candidates absent from the snapshot continue
+---------------------------+
        |
        |  …_clinvar_absent.jsonl
        v
+---------------------------+
| S6  review_quality        |  strict (tier A) is merged into highconf
+---------------------------+
        |
        v
+---------------------------+
| S7  journal_priority      |  screen_grade is a label; rows are kept
+---------------------------+
        |
        |  …_review_highconf_journal_annotated.jsonl
        v
+---------------------------+          +---------------------------+
| S8  pathogenicity         |          | S8  independent negatives |
| full text → sentence tags |          | not written back          |
+---------------------------+          +---------------------------+
        |
        v
+---------------------------+
| S9  genomic_export        |  chr / pos / ref / alt, FASTA check, GFF3 region
+---------------------------+
        |
        |  variants_genomic_positive_ok.csv
        v
+---------------------------+
| deliver / new_site        |  drop a bare rsID with no paired c.HGVS
+---------------------------+
        |
        v
+---------------------------+
| compare                   |  126 gold rows with complete coordinates
+---------------------------+
```

### 3.2 Published chain

```text
117 articles
  → abstract links
  → full text or abstract enrichment
  → ClinVar-absent
  → highconf
  → pathogenicity labels
  → GRCh38 coordinates
  → new-site
  → comparison with 126 gold rows
```

## 4. Per-stage procedure

`configs/pipeline_gates.json` still stores the conditions and checklist for each stage. The configured order of one stage is:

```text
1. Read this stage's conditions and checklist
2. Run the stage entry point, or the preflight check
3. Write the pass set and keep rejected or routed rows when the stage splits them
4. Record counts, reject reasons, and samples
5. Use the checklist to decide whether the next stage may start
6. Enter the next stage only after this one has passed
```

### 4.1 What the runner actually does

`cohort` and `pipeline` do not pause for a manual APPROVE, HOLD, or BYPASS, and they do not write `gate_runs/.../gate_report.json`. Stage0 stops when a required file is missing. Stages 1–9 are called in order by `scripts/run.py`. The summary is `data/pubmed/runs/cohort_pipeline_report.json`. Rejected rows remain in that stage's own side files.

### 4.2 Single-PMID sandbox

To send one article through Stage0–9 without changing stage rules, use the isolated workspace:

```bash
python scripts/run.py single 41992294
python scripts/run.py single 41992294 --from stage5_clinvar_validation --force
python scripts/run.py single 41992294 --no-network
```

- Workspace: `data/pubmed/runs/single/<PMID>/` (`merged/`, `linked/`, `fulltext/`, `stage_report.json`).
- File names inside the sandbox match the production names.
- `stage_report.json` records status, input and output counts, and the first drop. Logs are `logs/stage*.log`. The article record is `article.json`.
- The sandbox must not overwrite production `pipeline/s2_linking/` or `ingest/merged/articles.jsonl`.
- Each stage is the existing entry point with different input and output paths.
- If the production merged file already contains the PMID, its metadata is reused. Otherwise the runner fetches that article from NCBI unless `--no-network` is set.
- `status=dropped` and `first_drop` in `stage_report.json` show where the article stopped. After a change, rerun with `--from`.
- Stage5 runs the enriched track by default. That track is the input to Stages 6–9.

## 5. What each stage does

Each stage reads the previous primary output and writes its own primary file. Rejected and side rows are kept. Heavy stages run in separate processes so the ClinVar and Ensembl indexes are released before the next stage.

### 5.1 Stage0 preflight

The runner checks that the query config, journal table, gate config, ClinVar snapshot, Ensembl GFF3, and GRCh38 FASTA/FAI are present inside the repository. A missing file stops the run before retrieval or resolution.

### 5.2 Stage1 article intake

`cohort` writes the 117 PMIDs listed in `fixtures/cohort/articles.jsonl` to `ingest/merged/articles.jsonl`. Each record keeps PMID, title, abstract, journal, publication date, DOI, and PMCID. Only `pipeline` mode calls PubMed with the seven queries.

### 5.3 Stage2 abstract linking

The stage finds a specific variant and a clinical entity in the title and abstract, then assigns link tier A, B, or C from sentence or abstract co-occurrence. An article passes only when both signals are present. An rsID is bound to the PubTator or dbSNP catalog gene. An HGVS expression may use the adjacent sentence and a two-sentence naming window.

### 5.4 Stage3 full-text enrichment

Pass articles from Stage2 are fetched as PMC BioC, Europe PMC JATS, or HTML when open access is available. A same-sentence or same-paragraph link can then add evidence. Articles without full text stay on the abstract track. Each record is stamped with whether the evidence came from the abstract or the full text. The primary output is `corpus_enriched.jsonl`, which feeds discovery tagging and ClinVar comparison.

### 5.5 Stage4 discovery tagging

Wording such as novel, previously reported, citation of a prior report, or phenotype expansion sets a discovery label. The default label is unclear. This stage does not delete rows. Novel and phenotype-expansion rows continue. Phrases such as `we_identified` are high confidence. The earliest occurrence inside this corpus is not a worldwide first report.

### 5.6 Stage5 ClinVar comparison

The pinned `variant_summary.txt.gz` is queried by gene plus c., gene plus p., or rsID. Candidates absent from that snapshot continue to review. Known variants and keys that cannot be compared stay on side files. Absent means the key is missing from this snapshot.

### 5.7 Stage6 review subset

From the absent candidates, the stage keeps rows with a resolved gene, high discovery confidence, and a trusted co-occurrence method, then deduplicates by `gene|variant`. Strict keeps link tier A only. Highconf is the union of strict and the extended set, which also allows tier B. Highconf is the positive input to journal annotation and genomic export.

### 5.8 Stage7 journal grade

`configs/journal_priority.json` labels each journal A, B, C, or WATCH and writes `screen_grade`. C and WATCH rows are kept. Journal grade does not decide whether a variant is valid or whether it is a new site.

### 5.9 Stage8 pathogenicity and negatives

For positives, open-access text already fetched in Stage3 is written to `cache/pdfs/pdf_text_cache/PMID*.txt`, then sentence-level pathogenicity labels are extracted. A label must be anchored to the target variant. A miss is recorded as not evaluated. Setting `VARIANT_CLINICAL_STAGE8_PDF=1` downloads PDFs before that text is written.

Independent negatives are drawn from the full-text pool and exclude any `gene|allele` already present in the highconf positives. This side track is not written back into Stages 5–6.

### 5.10 Stage9 genomic coordinates

The stage reads the journal-annotated highconf rows and the Stage8 pathogenicity sidecar, then resolves each literature allele to four GRCh38 columns. The attempt order is a complete coordinate in the text, ClinVar, a genomic literal in the full text, a MANE protein change, a MANE coding change, dbSNP, coding alleles only, then unresolved. Coordinates are cleared when REF disagrees with the FASTA. A resolved coordinate receives `genomic_region` from the GFF3 overlap: coding, UTR, intron, ncRNA, pseudogene, or intergenic. The region label does not replace the gene name from the literature. Positive rows with complete coordinates are written to `variants_genomic_positive_ok.csv`.

## 6. Delivery and gold comparison

### 6.1 New-site filter

After Stage9, `deliver` drops a bare rsID: a row whose variant text is only an `rs` number, with no `c.… (rs…)` or `rs… (c.…)` pair in the source text, does not enter the new-site set. A c.HGVS filled in later from dbSNP does not put that row back. c., p., NM, and genomic literals are kept. The same command writes the review pages.

### 6.2 Gold comparison

The gold table shipped with the repository is `configs/gold_standard_v3.csv`: 129 rows, of which 126 have all four coordinate columns. `compare` aligns rows by PMID, gene, and variant text, then compares the four columns. The 2026-09-23 cohort matched 120 of 126 exactly. Matched rows had no coordinate conflicts. Six gold rows were absent from the positive set with complete coordinates, so `compare` exits with code 1. Stage3 contacts PMC, so another networked run is not guaranteed to be byte-identical.
