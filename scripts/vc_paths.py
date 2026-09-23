"""Canonical data/pubmed layout (5 top-level role folders).

    data/pubmed/
      ref/        queries + ClinVar snapshot + Ensembl GRCh38 pointers
      ingest/     Stage1 retrieval (raw, pmids, records, merged)
      pipeline/   Stage2–9 products (s2_linking … s9_genomic)
      cache/      bulky reuse (fulltext, PDFs, dbSNP, Ensembl indexes)
      runs/       gate_runs, single-PMID sandboxes, stats
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PUBMED = ROOT / "data" / "pubmed"

REF = PUBMED / "ref"
INGEST = PUBMED / "ingest"
PIPELINE = PUBMED / "pipeline"
CACHE = PUBMED / "cache"
RUNS = PUBMED / "runs"

QUERIES_DIR = REF / "queries"
QUERIES_USED = QUERIES_DIR / "queries_used.json"
CLINVAR_DIR = REF / "clinvar"
CLINVAR_SUMMARY = CLINVAR_DIR / "variant_summary.txt.gz"
ENSEMBL_DIR = REF / "ensembl"
ENSEMBL_GFF3 = ENSEMBL_DIR / "Homo_sapiens.GRCh38.115.gff3"
ENSEMBL_FASTA = ENSEMBL_DIR / "Homo_sapiens.GRCh38.dna.toplevel.fa"
ENSEMBL_FAI = ENSEMBL_DIR / "Homo_sapiens.GRCh38.dna.toplevel.fa.fai"
REFERENCE_MANIFEST = REF / "reference_manifest.json"

RAW_DIR = INGEST / "raw"
PMIDS_DIR = INGEST / "pmids"
RECORDS_DIR = INGEST / "records"
MERGED_DIR = INGEST / "merged"
MERGED_ARTICLES = MERGED_DIR / "articles.jsonl"

S2 = PIPELINE / "s2_linking"
S4 = PIPELINE / "s4_discovery"
S5 = PIPELINE / "s5_clinvar"
S5_ENRICHED = S5 / "enriched"
S5_STRICT = S5 / "strict"
S5_FULLTEXT = S5 / "fulltext"
S6 = PIPELINE / "s6_review"
S8 = PIPELINE / "s8_negatives"
S9 = PIPELINE / "s9_genomic"

PUBTATOR = S2 / "pubtator_parsed.jsonl"
ARTICLES_PASS = S2 / "articles_pass.jsonl"
PAIRS = S2 / "pairs.jsonl"
PAIRS_STRICT = S2 / "pairs_strict.jsonl"
CORPUS_ENRICHED = S2 / "corpus_enriched.jsonl"
CORPUS_FULLTEXT = S2 / "corpus_fulltext.jsonl"

DISCOVERY_NOVEL = S4 / "corpus_enriched_novel.jsonl"

CLINVAR_ABSENT = S5_ENRICHED / "corpus_enriched_novel_clinvar_absent.jsonl"

REVIEW_STRICT = S6 / "corpus_enriched_novel_clinvar_absent_review_strict.jsonl"
REVIEW_HIGHCONF = S6 / "corpus_enriched_novel_clinvar_absent_review_highconf.jsonl"
REVIEW_HIGHCONF_ANN = (
    S6 / "corpus_enriched_novel_clinvar_absent_review_highconf_journal_annotated.jsonl"
)
REVIEW_STRICT_ANN = (
    S6 / "corpus_enriched_novel_clinvar_absent_review_strict_journal_annotated.jsonl"
)

NEGATIVES_DIR = S8 / "negatives_independent"
NEGATIVES_FULL_PASS = NEGATIVES_DIR / "negatives_full_pass.jsonl"
NEGATIVES_HARD = NEGATIVES_DIR / "negatives_full_pass_hard.jsonl"
CLINVAR_2026_BLB_DIR = S8 / "clinvar_2026_blb"

FT_ROOT = CACHE / "fulltext"
FT_RAW = FT_ROOT / "raw"
FT_PARSED = FT_ROOT / "parsed"
FT_LINKED = FT_ROOT / "linked"
PDFS = CACHE / "pdfs"
PDF_TEXT_CACHE = PDFS / "pdf_text_cache"
PDF_PATHO = PDFS / "corpus_enriched_novel_clinvar_absent_review_highconf_pdf_pathogenicity.jsonl"
DBSNP_CACHE = CACHE / "dbsnp"
ENSEMBL_CACHE_DIR = CACHE / "ensembl"
ENSEMBL_REGION_INDEX = ENSEMBL_CACHE_DIR / "ensembl_gff3_region_index.pkl"
ENSEMBL_MANE_INDEX = ENSEMBL_CACHE_DIR / "ensembl_mane_transcript_index.json"

GENOMIC_CSV_DIR = S9
GATE_RUNS = RUNS / "gate_runs"
SINGLE_ROOT = RUNS / "single"
STATS_DIR = RUNS / "stats"


def rewrite_legacy_pubmed_path(path: Path | str) -> Path:
    """Map pre-2026-08-17 pubmed / references / data/cache paths onto the 5-role tree."""
    p = Path(path)
    s = p.as_posix()
    mapped = _map_project_legacy(s)
    return Path(mapped) if mapped != s else p


def _map_project_legacy(s: str) -> str:
    marker = "/data/pubmed/"
    idx = s.find(marker)
    if idx >= 0:
        root = s[: idx + len(marker)]
        rest = s[idx + len(marker) :]
        new_rest = _map_pubmed_rest(rest)
        return root + new_rest if new_rest != rest else s
    if s.startswith("data/pubmed/"):
        rest = s[len("data/pubmed/") :]
        new_rest = _map_pubmed_rest(rest)
        return "data/pubmed/" + new_rest if new_rest != rest else s

    for old, new in (
        ("/references/ensembl/", "/data/pubmed/ref/ensembl/"),
        ("/data/cache/dbsnp/", "/data/pubmed/cache/dbsnp/"),
        (
            "/data/cache/ensembl_gff3_region_index.pkl",
            "/data/pubmed/cache/ensembl/ensembl_gff3_region_index.pkl",
        ),
        (
            "/data/cache/ensembl_mane_transcript_index.json",
            "/data/pubmed/cache/ensembl/ensembl_mane_transcript_index.json",
        ),
        ("/data/cache/", "/data/pubmed/cache/"),
    ):
        if old in s:
            return s.replace(old, new, 1)
    for old, new in (
        ("references/ensembl/", "data/pubmed/ref/ensembl/"),
        ("data/cache/dbsnp/", "data/pubmed/cache/dbsnp/"),
        (
            "data/cache/ensembl_gff3_region_index.pkl",
            "data/pubmed/cache/ensembl/ensembl_gff3_region_index.pkl",
        ),
        (
            "data/cache/ensembl_mane_transcript_index.json",
            "data/pubmed/cache/ensembl/ensembl_mane_transcript_index.json",
        ),
    ):
        if s == old.rstrip("/") or s.startswith(old):
            return new + s[len(old) :] if s.startswith(old) else new.rstrip("/")
    return s


def _map_pubmed_rest(rest: str) -> str:
    def swap(old: str, new: str) -> str | None:
        if rest == old.rstrip("/"):
            return new.rstrip("/")
        if rest.startswith(old):
            return new + rest[len(old) :]
        return None

    enriched = "linked/discovery/clinvar_novelty_enriched/"
    if rest == enriched.rstrip("/"):
        return "pipeline/s5_clinvar/enriched"
    if rest.startswith(enriched):
        tail = rest[len(enriched) :]
        for old, new in (
            ("genomic_allele_csv/", "pipeline/s9_genomic/"),
            ("negatives_independent/", "pipeline/s8_negatives/negatives_independent/"),
            ("pdfs_review_strict/", "cache/pdfs/"),
        ):
            if tail == old.rstrip("/"):
                return new.rstrip("/")
            if tail.startswith(old):
                return new + tail[len(old) :]
        fname = tail.split("/", 1)[0]
        if "review" in fname or "journal" in fname:
            return "pipeline/s6_review/" + tail
        return "pipeline/s5_clinvar/enriched/" + tail

    for old, new in (
        ("linked/discovery/clinvar_novelty_fulltext/", "pipeline/s5_clinvar/fulltext/"),
        ("linked/discovery/clinvar_novelty/", "pipeline/s5_clinvar/strict/"),
        ("linked/discovery/", "pipeline/s4_discovery/"),
        ("linked/", "pipeline/s2_linking/"),
        ("fulltext/", "cache/fulltext/"),
        ("pdfs_review_strict/", "cache/pdfs/"),
        ("merged/", "ingest/merged/"),
        ("pmids/", "ingest/pmids/"),
        ("records/", "ingest/records/"),
        ("raw/", "ingest/raw/"),
        ("queries/", "ref/queries/"),
        ("clinvar/", "ref/clinvar/"),
        ("gate_runs/", "runs/gate_runs/"),
        ("stats/", "runs/stats/"),
        ("single/", "runs/single/"),
    ):
        if old == "linked/discovery/clinvar_novelty/" and rest.startswith(
            "linked/discovery/clinvar_novelty_"
        ):
            continue
        hit = swap(old, new)
        if hit is not None:
            return hit
    return rest


class RemapPubmedRoot(type(Path())):
    """``pubmed / 'linked/pairs.jsonl'`` still resolves after the 5-folder move."""

    def __truediv__(self, key):  # type: ignore[override]
        return rewrite_legacy_pubmed_path(Path(str(self)) / key)
