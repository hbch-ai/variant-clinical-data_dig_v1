#!/usr/bin/env python3
"""Physically consolidated pipeline commands. Generated from the v4_01 stage sources.

Each section retains its source module boundary through descriptive global prefixes;
there is no dynamic source loading or embedded executable source text.
"""
from __future__ import annotations
import fulltext as fulltext_module
import vc_paths as vc_paths_module
import vc_text as vc_text_module

# === tag_discovery.py ===
"""
Tag variant–clinical evidence as novel discovery vs previously reported / citation.

Problem: co-occurrence extraction also captures sentences that merely *cite*
already-known variants (e.g. "previously reported c.349G>T ... ACMG").
Users need labels to keep *new* genotype–phenotype findings.

Labels (discovery_label):
  novel                 – wording claims novelty / first report / not in literature
  previously_reported   – wording says previously reported / known / ClinVar etc.
  citation_of_prior     – cites other authors' cases (Author et al. reported ...)
  phenotype_expansion   – known variant BUT claims new/expanded clinical spectrum
  unclear               – no decisive cue (default; keep for manual review)

Also attaches:
  paper_pub_date        – YYYY / YYYY-MM from merged articles
  variant_first_seen    – earliest pub_date of this variant string in *this corpus*
  is_first_in_corpus    – this record's paper date == variant_first_seen

Outputs under data/pubmed/pipeline/s4_discovery/:
  pairs_fulltext_tagged.jsonl
  corpus_fulltext_tagged.jsonl
  pairs_strict_tagged.jsonl
  corpus_strict_tagged.jsonl
  *_novel.jsonl         – discovery_label in {novel, phenotype_expansion}
  *_previously_reported.jsonl
  summary.json
"""
import argparse
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
tag_discovery_ROOT = Path(__file__).resolve().parents[1]
tag_discovery_MERGED = tag_discovery_ROOT / 'data' / 'pubmed' / 'ingest' / 'merged' / 'articles.jsonl'
tag_discovery_LINKED = tag_discovery_ROOT / 'data' / 'pubmed' / 'pipeline' / 's2_linking'
tag_discovery_OUT = tag_discovery_ROOT / 'data' / 'pubmed' / 'pipeline' / 's4_discovery'
tag_discovery_NOVEL_CUES: list[tuple[str, re.Pattern[str]]] = [('novel', re.compile('\\bnovel\\b', re.I)), ('first_report', re.compile('\\bfirst\\s+(?:report|description|case|identification)\\b', re.I)), ('not_reported', re.compile('(?:not|never|hitherto|previously\\s+un)\\s+(?:been\\s+)?(?:previously\\s+)?(?:report|describ|publish|identif)\\w*|not\\s+reported\\s+in\\s+the\\s+literature|unreported\\b|undescribed\\b', re.I)), ('we_identified', re.compile('\\b(?:we|our\\s+study|this\\s+study|the\\s+present\\s+study)\\s+(?:herein\\s+)?(?:identif|detect|found|discover|report|describ)\\w*', re.I)), ('de_novo_claim', re.compile('\\bde\\s+novo\\b', re.I)), ('newly_identified', re.compile('\\bnewly\\s+(?:identif|detect|found|describ)\\w*', re.I))]
tag_discovery_PREV_CUES: list[tuple[str, re.Pattern[str]]] = [('previously_reported', re.compile('\\bprevious(?:ly)?\\s+(?:report|describ|publish|identif|known)\\w*|\\balready\\s+(?:been\\s+)?(?:report|describ|publish)\\w*|\\bknown\\s+(?:pathogenic\\s+)?(?:variant|mutation)\\b', re.I)), ('database_known', re.compile('\\b(?:ClinVar|HGMD|gnomAD|LOVD|dbSNP)\\b|\\breported\\s+in\\s+(?:ClinVar|the\\s+literature|previous\\s+studies)\\b|\\bpathogenic\\s+according\\s+to\\s+(?:the\\s+)?AC[MG]', re.I)), ('same_as_prior_case', re.compile('\\bsame\\s+(?:mutation|variant)\\s+as\\b|\\bpreviously\\s+described\\s+(?:in|by)\\b', re.I))]
tag_discovery_CITATION_CUES: list[tuple[str, re.Pattern[str]]] = [('et_al_reported', re.compile("\\b[A-Z][A-Za-z\\-']+(?:\\s+[A-Z][A-Za-z\\-']+)?\\s+et\\s+al\\.?\\s+(?:\\(\\d{4}\\)\\s+)?(?:report|describ|identif|found|detect)\\w*", re.I)), ('cited_case', re.compile('\\bin\\s+a\\s+(?:previous|prior|earlier)\\s+(?:report|study|case)\\b|\\bas\\s+reported\\s+by\\b|\\bhave\\s+been\\s+reported\\b', re.I))]
tag_discovery_EXPANSION_CUES: list[tuple[str, re.Pattern[str]]] = [('expands_spectrum', re.compile('\\bexpand(?:s|ed|ing)?\\s+(?:the\\s+)?(?:phenotypic\\s+)?spectrum\\b|\\bnew\\s+clinical\\s+(?:feature|finding|manifestation|presentation)s?\\b|\\bpreviously\\s+unreport(?:ed)?\\s+(?:clinical|phenotyp)\\w*|\\bfirst\\s+(?:description|report)\\s+of\\s+.{0,40}(?:with|in)\\s+this\\s+variant\\b', re.I))]
tag_discovery_MONTH_MAP = {'Jan': '01', 'Feb': '02', 'Mar': '03', 'Apr': '04', 'May': '05', 'Jun': '06', 'Jul': '07', 'Aug': '08', 'Sep': '09', 'Oct': '10', 'Nov': '11', 'Dec': '12'}

def tag_discovery_pub_date_key(art: dict[str, Any] | None) -> str | None:
    if not art:
        return None
    pd = art.get('pub_date') or {}
    y = pd.get('year')
    if not y:
        return None
    m = pd.get('month') or ''
    if m in tag_discovery_MONTH_MAP:
        m = tag_discovery_MONTH_MAP[m]
    elif str(m).isdigit():
        m = f'{int(m):02d}'
    else:
        return str(y)
    return f'{y}-{m}'

def tag_discovery_window_around(text: str, variant: str, radius: int=160) -> str:
    if not text:
        return ''
    if not variant:
        return text
    idx = text.find(variant)
    if idx < 0:
        loose = re.sub('\\s+', '', variant)
        compact = re.sub('\\s+', '', text)
        if loose and loose in compact:
            return text
        return text
    start = max(0, idx - radius)
    end = min(len(text), idx + len(variant) + radius)
    return text[start:end]

def tag_discovery_match_cues(text: str, patterns: list[tuple[str, re.Pattern[str]]]) -> list[str]:
    hits = []
    for name, pat in patterns:
        if pat.search(text or ''):
            hits.append(name)
    return hits

def tag_discovery_classify_discovery(evidence: str, variant: str) -> dict[str, Any]:
    """Classify using cues near the variant mention when possible."""
    full = evidence or ''
    local = tag_discovery_window_around(full, variant or '', radius=180)
    novel = tag_discovery_match_cues(local, tag_discovery_NOVEL_CUES) or tag_discovery_match_cues(full, tag_discovery_NOVEL_CUES)
    prev = tag_discovery_match_cues(local, tag_discovery_PREV_CUES) or tag_discovery_match_cues(full, tag_discovery_PREV_CUES)
    cite = tag_discovery_match_cues(local, tag_discovery_CITATION_CUES) or tag_discovery_match_cues(full, tag_discovery_CITATION_CUES)
    expand = tag_discovery_match_cues(local, tag_discovery_EXPANSION_CUES) or tag_discovery_match_cues(full, tag_discovery_EXPANSION_CUES)
    label = 'unclear'
    confidence = 'low'
    if cite and (not novel):
        label = 'citation_of_prior'
        confidence = 'high' if cite else 'medium'
    if prev and (not novel):
        label = 'previously_reported'
        confidence = 'high'
    if prev and expand and (not novel):
        label = 'phenotype_expansion'
        confidence = 'medium'
    if novel and (not prev) and (not cite):
        label = 'novel'
        high_cues = ('novel', 'not_reported', 'first_report', 'we_identified', 'newly_identified')
        confidence = 'high' if any((c in novel for c in high_cues)) else 'medium'
    if novel and prev:
        local_novel = tag_discovery_match_cues(local, tag_discovery_NOVEL_CUES)
        local_prev = tag_discovery_match_cues(local, tag_discovery_PREV_CUES)
        if local_prev and (not local_novel):
            label = 'previously_reported'
            confidence = 'medium'
        elif local_novel and (not local_prev):
            label = 'novel'
            confidence = 'medium'
        else:
            label = 'unclear'
            confidence = 'low'
            novel, prev = (local_novel or novel, local_prev or prev)
    if expand and novel:
        label = 'novel'
        confidence = 'high'
    return {'discovery_label': label, 'discovery_confidence': confidence, 'discovery_cues': {'novel': novel, 'previously_reported': prev, 'citation': cite, 'phenotype_expansion': expand}, 'discovery_window': local if local != full else None}

def tag_discovery_load_articles(merged_path: Path | None=None) -> dict[str, dict[str, Any]]:
    path = merged_path or tag_discovery_MERGED
    out: dict[str, dict[str, Any]] = {}
    with path.open(encoding='utf-8') as f:
        for line in f:
            if not line.strip():
                continue
            r = json.loads(line)
            out[str(r['pmid'])] = r
    return out

def tag_discovery_earliest_variant_dates(records: list[dict[str, Any]], articles: dict[str, dict[str, Any]]) -> dict[str, str]:
    """Map normalized variant string -> earliest paper_pub_date in this batch."""
    first: dict[str, str] = {}
    for r in records:
        v = (r.get('variant') or '').strip()
        if not v:
            continue
        pmid = str(r.get('pmid') or '')
        d = tag_discovery_pub_date_key(articles.get(pmid))
        if not d:
            continue
        if v not in first or d < first[v]:
            first[v] = d
    return first

def tag_discovery_tag_record(r: dict[str, Any], articles: dict[str, dict[str, Any]], first_seen: dict[str, str]) -> dict[str, Any]:
    out = dict(r)
    pmid = str(r.get('pmid') or '')
    art = articles.get(pmid)
    evidence = r.get('evidence_sentence') or r.get('text') or ''
    if not r.get('evidence_sentence') and art:
        evidence = ' '.join((x for x in [art.get('title') or '', art.get('abstract') or ''] if x))
    variant = r.get('variant') or ''
    disc = tag_discovery_classify_discovery(evidence, variant)
    out.update(disc)
    paper_date = tag_discovery_pub_date_key(art)
    out['paper_pub_date'] = paper_date
    out['paper_pub_year'] = (art.get('pub_date') or {}).get('year') if art else None
    vs = first_seen.get(variant.strip()) if variant else None
    out['variant_first_seen_in_corpus'] = vs
    out['is_first_in_corpus'] = bool(paper_date and vs and (paper_date == vs))
    out['keep_as_new_discovery'] = out['discovery_label'] in {'novel', 'phenotype_expansion'}
    return out

def tag_discovery_save_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8') as f:
        for row in rows:
            f.write(json.dumps(row, ensure_ascii=False) + '\n')

def tag_discovery_process_file(src: Path, articles: dict[str, dict[str, Any]], stem: str, out_dir: Path) -> dict[str, Any]:
    if not src.exists():
        return {'source': str(src), 'missing': True}
    rows = []
    with src.open(encoding='utf-8') as f:
        for line in f:
            if line.strip():
                rows.append(json.loads(line))
    first_seen = tag_discovery_earliest_variant_dates(rows, articles)
    tagged = [tag_discovery_tag_record(r, articles, first_seen) for r in rows]
    label_counts = Counter((r['discovery_label'] for r in tagged))
    conf_counts = Counter((r['discovery_confidence'] for r in tagged))
    novel_rows = [r for r in tagged if r.get('keep_as_new_discovery')]
    prev_rows = [r for r in tagged if r['discovery_label'] in {'previously_reported', 'citation_of_prior'}]
    unclear_rows = [r for r in tagged if r['discovery_label'] == 'unclear']
    first_and_novel = [r for r in novel_rows if r.get('is_first_in_corpus')]
    tag_discovery_save_jsonl(out_dir / f'{stem}_tagged.jsonl', tagged)
    tag_discovery_save_jsonl(out_dir / f'{stem}_novel.jsonl', novel_rows)
    tag_discovery_save_jsonl(out_dir / f'{stem}_previously_reported.jsonl', prev_rows)
    tag_discovery_save_jsonl(out_dir / f'{stem}_unclear.jsonl', unclear_rows)
    tag_discovery_save_jsonl(out_dir / f'{stem}_novel_first_in_corpus.jsonl', first_and_novel)
    tsv = out_dir / f'{stem}_novel_compact.tsv'
    with tsv.open('w', encoding='utf-8') as f:
        f.write('pmid\tpaper_pub_date\tdiscovery_label\tconfidence\tvariant\tclinical\tis_first_in_corpus\tcues\tevidence_sentence\n')
        for r in novel_rows:
            cues = ','.join((r.get('discovery_cues') or {}).get('novel') or (r.get('discovery_cues') or {}).get('phenotype_expansion') or [])
            ev = (r.get('evidence_sentence') or r.get('text') or '').replace('\t', ' ').replace('\n', ' ')
            f.write(f"{r.get('pmid')}\t{r.get('paper_pub_date') or ''}\t{r.get('discovery_label')}\t{r.get('discovery_confidence')}\t{(r.get('variant') or '').replace(chr(9), ' ')}\t{(r.get('clinical') or '').replace(chr(9), ' ')}\t{r.get('is_first_in_corpus')}\t{cues}\t{ev}\n")
    return {'source': str(src), 'n_total': len(tagged), 'label_counts': dict(label_counts), 'confidence_counts': dict(conf_counts), 'n_novel_or_expansion': len(novel_rows), 'n_previously_or_citation': len(prev_rows), 'n_unclear': len(unclear_rows), 'n_novel_and_first_in_corpus': len(first_and_novel)}

def tag_discovery_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--out-dir', type=Path, default=tag_discovery_ROOT / 'data' / 'pubmed' / 'pipeline' / 's4_discovery')
    parser.add_argument('--merged', type=Path, default=tag_discovery_MERGED, help='articles.jsonl for pub dates')
    parser.add_argument('--linked-dir', type=Path, default=tag_discovery_LINKED, help='directory containing pairs_*.jsonl / corpus_*.jsonl to tag')
    parser.add_argument('--stems', nargs='*', default=None, help='optional subset of stems, e.g. pairs_strict corpus_strict')
    parser.add_argument('--no-mirror-novel', action='store_true', help='do not copy corpus_fulltext_novel.jsonl back to linked-dir')
    args = parser.parse_args()
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    linked_dir: Path = args.linked_dir
    articles = tag_discovery_load_articles(args.merged)
    default_jobs = [('pairs_fulltext.jsonl', 'pairs_fulltext'), ('corpus_fulltext.jsonl', 'corpus_fulltext'), ('pairs_strict.jsonl', 'pairs_strict'), ('corpus_strict.jsonl', 'corpus_strict'), ('pairs_enriched.jsonl', 'pairs_enriched'), ('corpus_enriched.jsonl', 'corpus_enriched'), ('pairs.jsonl', 'pairs'), ('corpus.jsonl', 'corpus')]
    jobs = []
    for fname, stem in default_jobs:
        if args.stems and stem not in args.stems:
            continue
        src = linked_dir / fname
        if src.exists():
            jobs.append((src, stem))
    if not jobs:
        for src in sorted(linked_dir.glob('*.jsonl')):
            stem = src.stem
            if args.stems and stem not in args.stems:
                continue
            if stem.endswith('_tagged') or 'novel' in stem or 'previously' in stem or ('unclear' in stem):
                continue
            jobs.append((src, stem))
    summary: dict[str, Any] = {'created_at': datetime.now(timezone.utc).isoformat(), 'merged': str(args.merged), 'linked_dir': str(linked_dir), 'label_definitions': {'novel': 'Text claims novel / first report / not reported in literature / we identified', 'previously_reported': 'Text marks variant as previously reported / known / database-known', 'citation_of_prior': "Text cites other authors' prior cases for this variant–phenotype", 'phenotype_expansion': 'Known variant but claims new/expanded clinical spectrum', 'unclear': 'No decisive novelty/prior cue; needs review or external databases'}, 'files': {}}
    for src, stem in jobs:
        summary['files'][stem] = tag_discovery_process_file(src, articles, stem, out_dir)
    (out_dir / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    novel_ft = out_dir / 'corpus_fulltext_novel.jsonl'
    if novel_ft.exists() and (not args.no_mirror_novel):
        (linked_dir / 'corpus_fulltext_novel.jsonl').write_bytes(novel_ft.read_bytes())
    print(json.dumps(summary, ensure_ascii=False, indent=2))

# === filter_review_corpus.py ===
"""Filter gene-confirmed ClinVar-absent rows into high-confidence manual-review subsets.

Input default: corpus_enriched_novel_clinvar_absent.jsonl (Stage5 enriched absent track).

Two deduplicated outputs (one row per gene|variant locus):
  strict   – tier A + regex co-occurrence + lexicon entity clinical + trusted gene source
  extended – tier A or B with the same quality gates (broader manual review pool)

Also writes companion TSV files and filter_stats.json.
"""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any
import sys
filter_review_corpus_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
filter_review_corpus_DEFAULT_INPUT = vc_paths_module.CLINVAR_ABSENT
filter_review_corpus_DEFAULT_OUT_DIR = vc_paths_module.S6
filter_review_corpus_REGEX_METHODS = frozenset({'same_sentence_cooccurrence', 'same_abstract_cooccurrence', 'fulltext_same_sentence_cooccurrence', 'fulltext_same_paragraph_cooccurrence', 'fulltext_known_variant_sentence'})
filter_review_corpus_ENTITY_CLINICAL_SOURCES = frozenset({'lexicon_entity', 'lexicon_sentence', 'lexicon_paragraph', 'lexicon_title'})
filter_review_corpus_TRUSTED_GENE_SOURCES = frozenset({'row_gene', 'pubtator_or_context'})
filter_review_corpus_METHOD_RANK = {'same_sentence_cooccurrence': 0, 'fulltext_same_sentence_cooccurrence': 1, 'fulltext_known_variant_sentence': 2, 'fulltext_same_paragraph_cooccurrence': 3, 'same_abstract_cooccurrence': 4}
filter_review_corpus_TIER_RANK = {'A': 0, 'B': 1, 'C': 2}
filter_review_corpus_TSV_FIELDS = ['review_tier', 'locus_key', 'pmid', 'gene', 'variant', 'clinical', 'tier', 'method', 'discovery_label', 'discovery_confidence', 'evidence_source', 'text_scope', 'section', 'clinical_source', 'gene_source', 'paper_pub_date', 'is_first_in_corpus', 'evidence_text', 'pubmed_url']

def filter_review_corpus_locus_key(row: dict[str, Any]) -> str:
    gene = ((row.get('clinvar_query') or {}).get('genes') or [''])[0]
    variant = (row.get('variant') or '').lower()
    return f'{gene}|{variant}'

def filter_review_corpus_row_score(row: dict[str, Any]) -> tuple:
    cq = row.get('clinvar_query') or {}
    return (filter_review_corpus_TIER_RANK.get(row.get('tier'), 9), filter_review_corpus_METHOD_RANK.get(row.get('method'), 9), 0 if row.get('evidence_source') == 'fulltext' else 1, 0 if row.get('evidence_sentence') else 1, 0 if cq.get('gene_source') == 'row_gene' else 1)

def filter_review_corpus_passes_quality_gates(row: dict[str, Any]) -> bool:
    genes = (row.get('clinvar_query') or {}).get('genes') or []
    if not genes or not genes[0]:
        return False
    if row.get('discovery_confidence') != 'high':
        return False
    if row.get('method') not in filter_review_corpus_REGEX_METHODS:
        return False
    if row.get('clinical_source') not in filter_review_corpus_ENTITY_CLINICAL_SOURCES:
        return False
    if vc_text_module.is_generic_clinical_text(row.get('clinical')):
        return False
    gene_source = (row.get('clinvar_query') or {}).get('gene_source')
    return gene_source in filter_review_corpus_TRUSTED_GENE_SOURCES

def filter_review_corpus_passes_strict(row: dict[str, Any]) -> bool:
    return filter_review_corpus_passes_quality_gates(row) and row.get('tier') == 'A'

def filter_review_corpus_passes_extended(row: dict[str, Any]) -> bool:
    return filter_review_corpus_passes_quality_gates(row) and row.get('tier') in ('A', 'B')

def filter_review_corpus_pick_best_per_locus(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    best: dict[str, dict[str, Any]] = {}
    for row in rows:
        key = filter_review_corpus_locus_key(row)
        if key not in best or filter_review_corpus_row_score(row) < filter_review_corpus_row_score(best[key]):
            best[key] = row
    return list(best.values())

def filter_review_corpus_evidence_text(row: dict[str, Any], max_len: int=500) -> str:
    for field in ('evidence_sentence', 'discovery_window'):
        val = row.get(field)
        if val and isinstance(val, str) and val.strip():
            text = ' '.join(val.split())
            return text[:max_len] + ('…' if len(text) > max_len else '')
    return ''

def filter_review_corpus_compact_row(row: dict[str, Any], review_tier: str) -> dict[str, Any]:
    cq = row.get('clinvar_query') or {}
    gene = (cq.get('genes') or [''])[0]
    out = {'review_tier': review_tier, 'locus_key': filter_review_corpus_locus_key(row), 'pmid': row.get('pmid'), 'gene': gene, 'variant': row.get('variant'), 'clinical': row.get('clinical'), 'tier': row.get('tier'), 'method': row.get('method'), 'discovery_label': row.get('discovery_label'), 'discovery_confidence': row.get('discovery_confidence'), 'evidence_source': row.get('evidence_source'), 'text_scope': row.get('text_scope'), 'clinical_source': row.get('clinical_source'), 'gene_source': cq.get('gene_source'), 'paper_pub_date': row.get('paper_pub_date'), 'is_first_in_corpus': row.get('is_first_in_corpus'), 'evidence_text': filter_review_corpus_evidence_text(row), 'pubmed_url': f"https://pubmed.ncbi.nlm.nih.gov/{row.get('pmid')}/"}
    if row.get('section'):
        out['section'] = row.get('section')
    return out

def filter_review_corpus_write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open('w', encoding='utf-8') as fh:
        for row in rows:
            fh.write(json.dumps(row, ensure_ascii=False) + '\n')

def filter_review_corpus_write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open('w', encoding='utf-8', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=filter_review_corpus_TSV_FIELDS, extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, '') for k in filter_review_corpus_TSV_FIELDS})

def filter_review_corpus_summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {'n_rows': len(rows), 'tier': dict(Counter((r.get('tier') for r in rows))), 'method': dict(Counter((r.get('method') for r in rows)).most_common()), 'evidence_source': dict(Counter((r.get('evidence_source') for r in rows))), 'gene_source': dict(Counter((r.get('gene_source') for r in rows))), 'clinical_source': dict(Counter((r.get('clinical_source') for r in rows)))}

def filter_review_corpus_main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=filter_review_corpus_DEFAULT_INPUT)
    parser.add_argument('--out-dir', type=Path, default=filter_review_corpus_DEFAULT_OUT_DIR)
    args = parser.parse_args()
    all_rows = [json.loads(line) for line in args.input.open(encoding='utf-8')]
    gene_confirmed = [r for r in all_rows if ((r.get('clinvar_query') or {}).get('genes') or [None])[0]]
    strict_pool = [r for r in gene_confirmed if filter_review_corpus_passes_strict(r)]
    extended_pool = [r for r in gene_confirmed if filter_review_corpus_passes_extended(r)]
    strict_rows = sorted([filter_review_corpus_compact_row(r, 'strict') for r in filter_review_corpus_pick_best_per_locus(strict_pool)], key=lambda r: (r['gene'], r['variant']))
    extended_only = [r for r in filter_review_corpus_pick_best_per_locus(extended_pool) if not filter_review_corpus_passes_strict(r)]
    extended_rows = sorted([filter_review_corpus_compact_row(r, 'extended') for r in extended_only], key=lambda r: (r['gene'], r['variant']))
    combined_rows = sorted(strict_rows + extended_rows, key=lambda r: (0 if r['review_tier'] == 'strict' else 1, r['gene'], r['variant']))
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    strict_jsonl = out_dir / 'corpus_enriched_novel_clinvar_absent_review_strict.jsonl'
    extended_jsonl = out_dir / 'corpus_enriched_novel_clinvar_absent_review_extended.jsonl'
    combined_jsonl = out_dir / 'corpus_enriched_novel_clinvar_absent_review_highconf.jsonl'
    combined_tsv = out_dir / 'corpus_enriched_novel_clinvar_absent_review_highconf.tsv'
    stats_path = out_dir / 'review_filter_stats.json'
    filter_review_corpus_write_jsonl(strict_jsonl, strict_rows)
    filter_review_corpus_write_jsonl(extended_jsonl, extended_rows)
    filter_review_corpus_write_jsonl(combined_jsonl, combined_rows)
    filter_review_corpus_write_tsv(combined_tsv, combined_rows)
    stats = {'input': str(args.input), 'input_rows_total': len(all_rows), 'input_rows_gene_confirmed': len(gene_confirmed), 'input_rows_gene_unconfirmed_skipped': len(all_rows) - len(gene_confirmed), 'filters': {'gene_confirmed': 'clinvar_query.genes[0] non-empty', 'discovery_confidence': 'high', 'method': sorted(filter_review_corpus_REGEX_METHODS), 'clinical_source': sorted(filter_review_corpus_ENTITY_CLINICAL_SOURCES), 'gene_source': sorted(filter_review_corpus_TRUSTED_GENE_SOURCES), 'strict_tier': 'A', 'extended_tier': 'A or B', 'dedupe': 'one row per gene|variant; prefer tier A, tighter method, fulltext evidence', 'excluded_methods': ['pubtator_variant_disease', 'pubtator_gene_disease_relation'], 'excluded_clinical_source': ['lexicon_cue', 'narrative_context', 'pubtator_abstract'], 'excluded_clinical_text': 'generic bare heads (syndrome/disorder/dysplasia/…)', 'excluded_gene_source': ['text_guess']}, 'strict': filter_review_corpus_summarize(strict_rows), 'extended_only': filter_review_corpus_summarize(extended_rows), 'combined': filter_review_corpus_summarize(combined_rows), 'outputs': {'strict_jsonl': str(strict_jsonl), 'extended_jsonl': str(extended_jsonl), 'combined_jsonl': str(combined_jsonl), 'combined_tsv': str(combined_tsv)}}
    stats_path.write_text(json.dumps(stats, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps(stats, indent=2, ensure_ascii=False))

# === filter_journal_priority.py ===
"""Annotate review candidates with a curated journal priority tier.

Journal tier is a review/resource-allocation signal only. It must never replace
variant evidence, discovery, ClinVar, or gene/phenotype quality gates.
"""
import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any
filter_journal_priority_ROOT = Path(__file__).resolve().parents[1]
filter_journal_priority_DEFAULT_INPUT = filter_journal_priority_ROOT / 'data/pubmed/pipeline/s6_review' / 'corpus_enriched_novel_clinvar_absent_review_strict.jsonl'
filter_journal_priority_DEFAULT_ARTICLES = filter_journal_priority_ROOT / 'data/pubmed/ingest/merged/articles.jsonl'
filter_journal_priority_DEFAULT_CONFIG = filter_journal_priority_ROOT / 'configs/journal_priority.json'

def filter_journal_priority_normalize_name(value: str | None) -> str:
    text = (value or '').casefold()
    text = text.replace('&', ' and ')
    return re.sub('[^a-z0-9]+', ' ', text).strip()

def filter_journal_priority_normalize_issn(value: str | None) -> str:
    return re.sub('[^0-9Xx]', '', value or '').upper()

def filter_journal_priority_load_jsonl(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding='utf-8') as fh:
        return [json.loads(line) for line in fh if line.strip()]

def filter_journal_priority_write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open('w', encoding='utf-8') as fh:
        for row in rows:
            fh.write(json.dumps(row, ensure_ascii=False) + '\n')

def filter_journal_priority_build_indexes(config: dict[str, Any]) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    by_issn: dict[str, dict[str, Any]] = {}
    by_name: dict[str, dict[str, Any]] = {}
    for journal in config.get('journals') or []:
        for issn in journal.get('issns') or []:
            key = filter_journal_priority_normalize_issn(str(issn))
            if key:
                by_issn[key] = journal
        names = [journal.get('canonical_name'), *(journal.get('aliases') or [])]
        for name in names:
            key = filter_journal_priority_normalize_name(str(name or ''))
            if key:
                by_name[key] = journal
    return (by_issn, by_name)

def filter_journal_priority_match_journal(article: dict[str, Any], by_issn: dict[str, dict[str, Any]], by_name: dict[str, dict[str, Any]]) -> tuple[dict[str, Any] | None, str]:
    issn = filter_journal_priority_normalize_issn(str(article.get('issn') or ''))
    if issn and issn in by_issn:
        return (by_issn[issn], 'issn')
    for field in ('journal_iso', 'journal'):
        name = filter_journal_priority_normalize_name(str(article.get(field) or ''))
        if name and name in by_name:
            return (by_name[name], f'{field}_alias')
    return (None, 'default')

def filter_journal_priority_annotate_row(row: dict[str, Any], article: dict[str, Any] | None, config: dict[str, Any], by_issn: dict[str, dict[str, Any]], by_name: dict[str, dict[str, Any]]) -> dict[str, Any]:
    out = dict(row)
    article = article or {}
    policy = config.get('policy') or {}
    match, method = filter_journal_priority_match_journal(article, by_issn, by_name)
    tier = str((match or {}).get('tier') or policy.get('default_tier') or 'C')
    score = int((policy.get('tier_scores') or {}).get(tier, 0))
    if row.get('evidence_source') == 'fulltext':
        score += int(policy.get('fulltext_bonus') or 0)
    out.update({'journal': article.get('journal'), 'journal_iso': article.get('journal_iso'), 'journal_issn': article.get('issn'), 'journal_canonical_name': (match or {}).get('canonical_name'), 'journal_tier': tier, 'journal_domains': (match or {}).get('domains') or [], 'journal_match_method': method, 'journal_priority_score': score, 'journal_priority_reason': (match or {}).get('reason') or 'not on curated A/B list; retained as default tier', 'journal_list_version': config.get('list_version')})
    return out

def filter_journal_priority_write_unmatched(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = ['journal', 'journal_iso', 'journal_issn', 'n_rows', 'n_pmids']
    grouped: dict[tuple[str, str, str], dict[str, Any]] = {}
    for row in rows:
        key = (str(row.get('journal') or ''), str(row.get('journal_iso') or ''), str(row.get('journal_issn') or ''))
        bucket = grouped.setdefault(key, {'journal': key[0], 'journal_iso': key[1], 'journal_issn': key[2], 'n_rows': 0, 'pmids': set()})
        bucket['n_rows'] += 1
        bucket['pmids'].add(str(row.get('pmid') or ''))
    output = sorted(grouped.values(), key=lambda x: (-x['n_rows'], x['journal_iso']))
    with path.open('w', encoding='utf-8', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for item in output:
            writer.writerow({'journal': item['journal'], 'journal_iso': item['journal_iso'], 'journal_issn': item['journal_issn'], 'n_rows': item['n_rows'], 'n_pmids': len(item['pmids'])})

def filter_journal_priority_main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=filter_journal_priority_DEFAULT_INPUT)
    parser.add_argument('--articles', type=Path, default=filter_journal_priority_DEFAULT_ARTICLES)
    parser.add_argument('--config', type=Path, default=filter_journal_priority_DEFAULT_CONFIG)
    parser.add_argument('--out-dir', type=Path, default=None)
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding='utf-8'))
    rows = filter_journal_priority_load_jsonl(args.input)
    articles = {str(row.get('pmid')): row for row in filter_journal_priority_load_jsonl(args.articles) if row.get('pmid') is not None}
    by_issn, by_name = filter_journal_priority_build_indexes(config)
    annotated = [filter_journal_priority_annotate_row(row, articles.get(str(row.get('pmid'))), config, by_issn, by_name) for row in rows]
    annotated.sort(key=lambda row: (-int(row.get('journal_priority_score') or 0), str(row.get('gene') or ''), str(row.get('variant') or '')))
    tier_a = [row for row in annotated if row.get('journal_tier') == 'A']
    tier_ab = [row for row in annotated if row.get('journal_tier') in {'A', 'B'}]
    other = [row for row in annotated if row.get('journal_tier') not in {'A', 'B'}]
    out_dir = args.out_dir or args.input.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = args.input.stem
    paths = {'annotated': out_dir / f'{stem}_journal_annotated.jsonl', 'tier_a': out_dir / f'{stem}_journal_A.jsonl', 'tier_ab': out_dir / f'{stem}_journal_AB.jsonl', 'other': out_dir / f'{stem}_journal_other.jsonl', 'unmatched': out_dir / f'{stem}_journal_unmatched.tsv', 'stats': out_dir / f'{stem}_journal_filter_stats.json'}
    filter_journal_priority_write_jsonl(paths['annotated'], annotated)
    filter_journal_priority_write_jsonl(paths['tier_a'], tier_a)
    filter_journal_priority_write_jsonl(paths['tier_ab'], tier_ab)
    filter_journal_priority_write_jsonl(paths['other'], other)
    filter_journal_priority_write_unmatched(paths['unmatched'], other)
    stats = {'created_from': str(args.input), 'articles': str(args.articles), 'config': str(args.config), 'journal_list_version': config.get('list_version'), 'policy': config.get('policy'), 'n_input_rows': len(rows), 'n_input_pmids': len({str(row.get('pmid')) for row in rows}), 'missing_article_metadata': sum((1 for row in rows if str(row.get('pmid')) not in articles)), 'tier_counts_rows': dict(Counter((row.get('journal_tier') for row in annotated))), 'tier_counts_pmids': {tier: len({str(row.get('pmid')) for row in annotated if row.get('journal_tier') == tier}) for tier in sorted({str(row.get('journal_tier')) for row in annotated})}, 'match_method_counts': dict(Counter((row.get('journal_match_method') for row in annotated))), 'top_journals': dict(Counter((row.get('journal_iso') or row.get('journal') or 'MISSING' for row in annotated)).most_common(30)), 'outputs': {key: str(value) for key, value in paths.items()}, 'warning': 'Journal priority controls review order only; it does not establish variant validity, pathogenicity, or literature novelty.'}
    paths['stats'].write_text(json.dumps(stats, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps(stats, indent=2, ensure_ascii=False))

# === annotate_pathogenicity_from_pdfs.py ===
"""
Annotate pathogenicity (pathogenic / likely pathogenic / VUS / likely benign / benign)
by searching text extracted from local review-strict PDFs.

This is a heuristic keyword-based matcher intended to speed up manual review.
It aligns each review row (gene+variant) to the corresponding PMID PDF.
"""
import argparse
import json
import logging
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable
annotate_pathogenicity_from_pdfs_ROOT = Path(__file__).resolve().parents[1]
annotate_pathogenicity_from_pdfs_SCRIPTS = Path(__file__).resolve().parent

def annotate_pathogenicity_from_pdfs_run_pdftotext_layout(pdf_path: Path) -> str:
    """Extract PDF text via pdftotext (preferred)."""
    proc = subprocess.run(['pdftotext', '-layout', '-enc', 'UTF-8', str(pdf_path), '-'], check=False, capture_output=True, text=True, timeout=180)
    out = proc.stdout or ''
    if proc.returncode != 0 and (not out.strip()):
        return ''
    return out

def annotate_pathogenicity_from_pdfs_normalize_variant_for_search(variant: str | None) -> list[str]:
    """Generate a few string variants for text matching."""
    v = (variant or '').strip()
    if not v:
        return []
    v_l = v.lower()
    toks: list[str] = []
    toks.append(v_l)
    m = re.search('(rs\\d{4,})', v_l)
    if m:
        toks.append(m.group(1))
    v_relaxed = re.sub('[\\(\\)\\s]', '', v_l)
    toks.append(v_relaxed)
    if v_relaxed.startswith('p.'):
        toks.append(v_relaxed[2:])
    if v_relaxed.startswith('c.'):
        toks.append(v_relaxed[2:])
    m2 = re.search('\\b([a-z]{1,4}\\d+[^;\\s,]*)\\b', v_relaxed, re.I)
    if m2:
        toks.append(m2.group(1).lower())
    seen: set[str] = set()
    out: list[str] = []
    for t in toks:
        if t and t not in seen:
            out.append(t)
            seen.add(t)
    return out

@dataclass(frozen=True)
class annotate_pathogenicity_from_pdfs_MatchResult:
    label: str
    confidence: str
    evidence: str

def annotate_pathogenicity_from_pdfs_classify_pathogenicity(*, text: str, gene: str, variant: str) -> annotate_pathogenicity_from_pdfs_MatchResult:
    """Heuristic classification from text context near variant mention."""
    re_likely_path = re.compile('\\blikely\\s+pathogenic\\b', re.I)
    re_path = re.compile('(?<!likely\\s)\\bpathogenic\\b', re.I)
    re_likely_benign = re.compile('\\blikely\\s+benign\\b', re.I)
    re_benign = re.compile('\\bbenign\\b', re.I)
    re_vus = re.compile('\\bvariant\\s+of\\s+uncertain\\s+significance\\b|\\bVUS\\b|\\buncertain\\s+significance\\b', re.I)
    variant_terms = annotate_pathogenicity_from_pdfs_normalize_variant_for_search(variant)
    gene_l = (gene or '').strip().lower()
    text_l = text.lower()
    anchor_idx = None
    for term in variant_terms:
        idx = text_l.find(term)
        if idx != -1:
            anchor_idx = idx
            break
    if anchor_idx is None and gene_l:
        anchor_idx = text_l.find(gene_l)
    window_before = 900
    window_after = 1400
    if anchor_idx is None:
        context = text
        confidence_anchor = 'low'
    else:
        lo = max(0, anchor_idx - window_before)
        hi = min(len(text), anchor_idx + window_after)
        context = text[lo:hi]
        confidence_anchor = 'high' if context.strip() else 'medium'

    def _find_snip(rx: re.Pattern[str]) -> str | None:
        m = rx.search(context)
        if not m:
            return None
        s0 = max(0, m.start() - 180)
        s1 = min(len(context), m.end() + 220)
        snip = context[s0:s1].strip()
        return snip
    sn_likely_path = _find_snip(re_likely_path)
    sn_path = _find_snip(re_path)
    sn_likely_benign = _find_snip(re_likely_benign)
    sn_benign = _find_snip(re_benign)
    sn_vus = _find_snip(re_vus)
    path_like = bool(sn_likely_path or sn_path)
    benign_like = bool(sn_likely_benign or sn_benign)
    if path_like and benign_like:
        evidence = sn_likely_path or sn_path or sn_likely_benign or sn_benign or ''
        return annotate_pathogenicity_from_pdfs_MatchResult(label='conflicting', confidence='medium', evidence=evidence)
    if sn_likely_path:
        return annotate_pathogenicity_from_pdfs_MatchResult(label='likely_pathogenic', confidence=confidence_anchor, evidence=sn_likely_path)
    if sn_path:
        return annotate_pathogenicity_from_pdfs_MatchResult(label='pathogenic', confidence=confidence_anchor, evidence=sn_path)
    if sn_vus:
        return annotate_pathogenicity_from_pdfs_MatchResult(label='uncertain', confidence=confidence_anchor, evidence=sn_vus)
    if sn_likely_benign:
        return annotate_pathogenicity_from_pdfs_MatchResult(label='likely_benign', confidence=confidence_anchor, evidence=sn_likely_benign)
    if sn_benign:
        return annotate_pathogenicity_from_pdfs_MatchResult(label='benign', confidence=confidence_anchor, evidence=sn_benign)
    return annotate_pathogenicity_from_pdfs_MatchResult(label='not_evaluated', confidence='low' if anchor_idx is not None else 'very_low', evidence='')

def annotate_pathogenicity_from_pdfs_main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=annotate_pathogenicity_from_pdfs_ROOT / 'data/pubmed/pipeline/s6_review/corpus_enriched_novel_clinvar_absent_review_strict.jsonl')
    parser.add_argument('--pdf-dir', type=Path, default=annotate_pathogenicity_from_pdfs_ROOT / 'data/pubmed/cache/pdfs')
    parser.add_argument('--out-dir', type=Path, default=None)
    parser.add_argument('--text-cache-dir', type=Path, default=None)
    parser.add_argument('--limit', type=int, default=0, help='optional limit for debugging')
    args = parser.parse_args()
    out_dir = args.out_dir or args.pdf_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    text_cache_dir = args.text_cache_dir or out_dir / 'pdf_text_cache'
    text_cache_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    logger = logging.getLogger('pdf_pathogenicity')
    rows: list[dict[str, Any]] = []
    for line in args.input.open(encoding='utf-8'):
        if line.strip():
            rows.append(json.loads(line))
    if args.limit:
        rows = rows[:args.limit]
    pmid_to_text: dict[str, str] = {}
    out_jsonl = out_dir / f'{args.input.stem}_pdf_pathogenicity.jsonl'
    out_tsv = out_dir / f'{args.input.stem}_pdf_pathogenicity.tsv'
    tsv_cols = ['pmid', 'gene', 'variant', 'pdf_pathogenicity', 'pdf_pathogenicity_confidence', 'pdf_pathogenicity_evidence', 'pubmed_url']
    written = 0
    with out_jsonl.open('w', encoding='utf-8') as oj, out_tsv.open('w', encoding='utf-8', newline='') as ot:
        ot.write('\t'.join(tsv_cols) + '\n')
        for r in rows:
            pmid = str(r.get('pmid') or '').strip()
            gene = str(r.get('gene') or '').strip()
            variant = str(r.get('variant') or '').strip()
            pdf_path = args.pdf_dir / f'PMID{pmid}.pdf'
            cache_path = text_cache_dir / f'PMID{pmid}.txt'
            if pmid not in pmid_to_text:
                text = ''
                if cache_path.is_file():
                    text = cache_path.read_text(encoding='utf-8', errors='replace')
                elif pdf_path.is_file():
                    text = annotate_pathogenicity_from_pdfs_run_pdftotext_layout(pdf_path)
                    if text.strip():
                        cache_path.write_text(text, encoding='utf-8', errors='replace')
                pmid_to_text[pmid] = text
            text = pmid_to_text.get(pmid, '')
            if not text.strip():
                r['pdf_pathogenicity'] = 'not_evaluated'
                r['pdf_pathogenicity_confidence'] = 'no_fulltext'
                r['pdf_pathogenicity_evidence'] = ''
            else:
                res = annotate_pathogenicity_from_pdfs_classify_pathogenicity(text=text, gene=gene, variant=variant)
                r['pdf_pathogenicity'] = res.label
                r['pdf_pathogenicity_confidence'] = res.confidence
                r['pdf_pathogenicity_evidence'] = res.evidence
            oj.write(json.dumps(r, ensure_ascii=False) + '\n')
            ot.write('\t'.join([str(r.get('pmid') or ''), str(r.get('gene') or ''), str(r.get('variant') or ''), str(r.get('pdf_pathogenicity') or ''), str(r.get('pdf_pathogenicity_confidence') or ''), str(r.get('pdf_pathogenicity_evidence') or '').replace('\t', ' ').replace('\n', ' '), str(r.get('pubmed_url') or '')]) + '\n')
            written += 1
            if written % 20 == 0:
                logger.info('processed %s/%s', written, len(rows))
    logger.info('done. wrote %s rows', written)
    logger.info('outputs: %s , %s', out_jsonl, out_tsv)

# === extract_negative_candidates.py ===
"""
Extract negative alleles from literature (benign / LB / control / VUS / assoc-null).

Two operating modes:

1) **Anchor mode** (`--pmid-source jsonl`, default): PMIDs come from a positive /
   review JSONL; exclude the indexed pathogenic allele within each paper
   (`same_paper_other_site` / control / association_null).

2) **Independent corpus mode** (`--pmid-source fulltext|full_pass|…`):
   screen a literature pool that is **not** required to overlap the positive set.
   Roles become `corpus_benign_mention` (plus control / association_null).
   Optionally collide-filter against a positive JSONL via `--exclude-positives`.

Main JSONL is deduped by negative_locus_key (gene|allele). Hard/soft splits are
written when `--split-hard-soft` is set (default on for independent pools).
"""
import argparse
import csv
import json
import logging
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable
extract_negative_candidates_ROOT = Path(__file__).resolve().parents[1]
extract_negative_candidates_SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(extract_negative_candidates_SCRIPTS))
extract_negative_candidates_DEFAULT_INPUT = extract_negative_candidates_ROOT / 'data/pubmed/pipeline/s6_review' / 'corpus_enriched_novel_clinvar_absent_review_strict.jsonl'
extract_negative_candidates_DEFAULT_HIGHCONF = extract_negative_candidates_ROOT / 'data/pubmed/pipeline/s6_review' / 'corpus_enriched_novel_clinvar_absent_review_highconf.jsonl'
extract_negative_candidates_DEFAULT_MERGED = extract_negative_candidates_ROOT / 'data/pubmed/ingest/merged/articles.jsonl'
extract_negative_candidates_DEFAULT_ARTICLES_PASS = extract_negative_candidates_ROOT / 'data/pubmed/pipeline/s2_linking/articles_pass.jsonl'
extract_negative_candidates_DEFAULT_PARSED_DIR = extract_negative_candidates_ROOT / 'data/pubmed/cache/fulltext/parsed'
extract_negative_candidates_DEFAULT_RAW_DIR = extract_negative_candidates_ROOT / 'data/pubmed/cache/fulltext/raw'
extract_negative_candidates_DEFAULT_PDF_CACHE = extract_negative_candidates_ROOT / 'data/pubmed/cache/pdfs/pdf_text_cache'
extract_negative_candidates_DEFAULT_INDEPENDENT_OUT = extract_negative_candidates_ROOT / 'data/pubmed/pipeline/s8_negatives/negatives_independent'
extract_negative_candidates_DEFAULT_EXCLUDE_POS = extract_negative_candidates_ROOT / 'data/pubmed/pipeline/s6_review' / 'corpus_enriched_novel_clinvar_absent_review_highconf_journal_annotated.jsonl'
extract_negative_candidates_HARD_LABELS = frozenset({'benign', 'likely_benign', 'not_pathogenic'})
extract_negative_candidates_SOFT_LABELS = frozenset({'uncertain', 'association_null'})
extract_negative_candidates_RE_LIKELY_BENIGN = re.compile('\\blikely\\s+benign\\b', re.I)
extract_negative_candidates_RE_BENIGN = re.compile('\\bclassified\\s+as\\s+benign\\b|\\b(?:is|was|are|were)\\s+benign\\b|\\bbenign\\s+(?:variant|variants|mutation|mutations|allele|polymorphism)\\b|\\bnot\\s+pathogenic\\b|\\bnon[- ]pathogenic\\b', re.I)
extract_negative_candidates_RE_VUS = re.compile('\\bclassified\\s+as\\s+(?:a\\s+)?VUS\\b|\\b(?:a\\s+)?variant\\s+of\\s+uncertain(?:\\s+clinical)?\\s+significance\\b|\\buncertain\\s+(?:clinical\\s+)?significance\\b|\\breclassified\\s+as\\s+(?:a\\s+)?VUS\\b', re.I)
extract_negative_candidates_RE_NULL_ASSOC = re.compile('\\bno\\s+significant\\s+associations?\\b|\\bnot\\s+significantly\\s+associated\\b|\\bnot\\s+associated\\s+with\\b|\\bmay\\s+not\\s+be\\s+related\\b|\\bnot\\s+related\\s+to\\b|\\bno\\s+association\\s+(?:was\\s+)?(?:found|observed|detected)\\b|\\bpolymorphism\\s+may\\s+not\\s+be\\s+related\\b', re.I)
extract_negative_candidates_RE_CONTROL_CUE = re.compile('\\bbenign\\s+control\\b|\\bcontrol\\s+(?:variant|allele|mutation)\\b|\\bnegative[- ]control\\s+(?:variant|allele)\\b|\\bhealthy[- ]control\\s+(?:variant|allele|genotype)\\b|\\bas\\s+a\\s+(?:benign|negative)\\s+control\\b|\\bused\\s+as\\s+(?:a\\s+)?(?:benign|negative)\\s+control\\b|\\bcompared\\s+(?:with|to)\\s+(?:the\\s+)?(?:benign|non[- ]pathogenic)\\b', re.I)
extract_negative_candidates_RE_EXCLUDE = re.compile('\\bbenign\\s+(?:insulinoma|lesions?|gingival|fibromatosis|muscular|tumor|tumour|prostatic|hyperplasia|epilepsy|neonatal)\\b|\\bwithout\\s+benign\\s+variations?\\b|\\bB/LB/P/LP\\b|\\bBA1/BS1\\b|\\bMutationTaster\\b|\\bPolyPhen\\b|\\bcolou?red\\s+based\\s+on\\b|\\bempty\\s+vector\\b|\\bPathogenesity\\b|\\b1000G\\b.*\\bExAC\\b|\\bClinVar\\b.*\\bLoVD\\b', re.I)
extract_negative_candidates_RE_OVERTURNED_VUS = re.compile('\\bpreviously\\s+classified\\s+as\\b|\\bdespite\\b.{0,40}\\b(?:VUS|uncertain)\\b|\\bconfirming\\s+(?:its\\s+)?pathogenic\\b|\\bpathogenic\\s+impact\\b', re.I)
extract_negative_candidates_RE_CLASSIFY_PAREN = re.compile('(?P<label>likely\\s+benign|benign|not\\s+pathogenic|non[- ]pathogenic|VUS|variant\\s+of\\s+uncertain(?:\\s+clinical)?\\s+significance|uncertain(?:\\s+clinical)?\\s+significance)\\s*\\((?P<body>[^)]{3,200})\\)', re.I)
extract_negative_candidates_RE_VARIANT_RECLASS_LB = re.compile('(?P<var>(?:c\\.|p\\.|rs)[\\w()*>.\\-+_/ ]+)[^.]{0,60}?reclassified\\s+as\\s+(?P<label>likely\\s+benign|benign|VUS|variant\\s+of\\s+uncertain(?:\\s+clinical)?\\s+significance)', re.I)
extract_negative_candidates_RE_VARIANT_THEN_LABEL = re.compile('(?P<var>(?:c\\.|p\\.|rs)[\\w()*>.\\-+_/]+)(?P<mid>.{0,40}?)(?P<label>likely\\s+benign|classified\\s+as\\s+benign|not\\s+pathogenic|non[- ]pathogenic|benign\\s+(?:variant|polymorphism)|classified\\s+as\\s+(?:a\\s+)?VUS|variant\\s+of\\s+uncertain(?:\\s+clinical)?\\s+significance)', re.I)
extract_negative_candidates_RE_LABEL_THEN_VARIANT = re.compile('(?P<label>likely\\s+benign|classified\\s+as\\s+benign|not\\s+pathogenic|non[- ]pathogenic|benign\\s+(?:variant|variants|polymorphism|mutation|mutations)|classified\\s+as\\s+(?:a\\s+)?VUS|variant\\s+of\\s+uncertain(?:\\s+clinical)?\\s+significance)(?P<mid>.{0,40}?)(?P<var>(?:c\\.|p\\.|rs)[\\w()*>.\\-+_/]+)', re.I)
extract_negative_candidates_BAD_GENE_TOKENS = {'LP', 'LB', 'VUS', 'P', 'B', 'HGVS', 'IVS', 'EXAC', 'GNOMAD', 'CLINVAR', 'LOVD', 'CADD', 'SIFT', 'FAD', 'FIG', 'FIGURE', 'REF', 'WT', 'WILDTYPE', 'ROS', 'HPP', 'HMNX', 'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6', 'PP1', 'PP2', 'PP3', 'PP4', 'PP5', 'PS1', 'PS2', 'PS3', 'PS4', 'PVS1', 'BA1', 'BS1', 'BS2', 'BS3', 'BS4', 'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7', 'ACMG', 'AMP', 'BENIGN', 'LIKELY', 'PATHOGENIC', 'CONTROL', 'CONTROLS', 'NEGATIVE', 'POSITIVE'}
extract_negative_candidates_LABEL_PRIORITY = {'benign': 0, 'likely_benign': 1, 'not_pathogenic': 2, 'uncertain': 3, 'association_null': 4}
extract_negative_candidates_ROLE_PRIORITY = {'control_variant': 0, 'same_paper_other_site': 1, 'corpus_benign_mention': 2, 'association_null_site': 3}
extract_negative_candidates_CONCRETE_OK = set(vc_text_module.CONCRETE_VARIANT_SUBTYPES) | {'rsid', 'nm_np_with_change', 'ivs'}
extract_negative_candidates_RE_RS_ONLY = re.compile('^rs\\d{4,}$', re.I)

@dataclass(frozen=True)
class extract_negative_candidates_NegativeAlleleHit:
    negative_variant: str
    negative_label: str
    negative_role: str
    confidence: str
    rule: str
    evidence: str
    evidence_source: str
    negative_gene: str | None = None

def extract_negative_candidates_sentences(text: str) -> list[str]:
    text = re.sub('\\s+', ' ', text or '').strip()
    if not text:
        return []
    return re.split('(?<=[.!?])\\s+(?=[A-Z0-9(])', text)

def extract_negative_candidates_normalize_label(raw: str) -> str:
    s = re.sub('\\s+', ' ', (raw or '').strip().lower())
    if 'likely benign' in s:
        return 'likely_benign'
    if 'not pathogenic' in s or 'non-pathogenic' in s or 'non pathogenic' in s:
        return 'not_pathogenic'
    if 'vus' in s or 'uncertain' in s:
        return 'uncertain'
    if 'association_null' in s or 'null association' in s:
        return 'association_null'
    return 'benign'

def extract_negative_candidates_locus_key_of(gene: str | None, variant: str | None) -> str:
    """Canonical gene|allele key for cross-PMID dedupe."""
    g = (gene or '').strip().upper() or 'NA'
    tok = extract_negative_candidates_clean_variant_token(variant or '') or (variant or '').strip()
    terms = annotate_pathogenicity_from_pdfs_normalize_variant_for_search(tok)
    core = sorted((t for t in terms if len(t) >= 4), key=len, reverse=True)
    allele = core[0] if core else tok.lower()
    allele = re.sub('[\\s()]', '', allele).lower()
    return f'{g}|{allele}'

def extract_negative_candidates_candidate_score(c: dict[str, Any]) -> tuple:
    conf_rank = {'high': 0, 'medium': 1, 'low': 2}
    return (extract_negative_candidates_LABEL_PRIORITY.get(str(c.get('negative_label') or ''), 9), extract_negative_candidates_ROLE_PRIORITY.get(str(c.get('negative_role') or ''), 9), conf_rank.get(str(c.get('negative_confidence') or ''), 9))

def extract_negative_candidates_dedupe_by_locus(candidates: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, int]]:
    """Keep one row per negative locus (gene|allele), best label/role/confidence."""
    best: dict[str, dict[str, Any]] = {}
    for c in candidates:
        key = extract_negative_candidates_locus_key_of(c.get('negative_gene'), c.get('negative_variant'))
        c = dict(c)
        c['negative_locus_key'] = key
        prev = best.get(key)
        if prev is None or extract_negative_candidates_candidate_score(c) < extract_negative_candidates_candidate_score(prev):
            best[key] = c
    out = list(best.values())
    out.sort(key=lambda r: (r.get('negative_locus_key') or '', r.get('pmid') or ''))
    counts: dict[str, int] = {}
    for c in out:
        lab = str(c.get('negative_label') or '')
        counts[lab] = counts.get(lab, 0) + 1
    return (out, counts)

def extract_negative_candidates_clean_variant_token(raw: str) -> str | None:
    raw = (raw or '').strip()
    if not raw:
        return None
    for subtype, pat in vc_text_module.VARIANT_PATTERNS:
        if subtype not in extract_negative_candidates_CONCRETE_OK:
            continue
        m = pat.search(raw)
        if not m:
            continue
        tok = vc_text_module.normalize_variant_token(m.group(0))
        if tok and len(tok) >= 4:
            return tok
    return None

def extract_negative_candidates_is_same_allele(a: str, b: str) -> bool:
    ta = annotate_pathogenicity_from_pdfs_normalize_variant_for_search(a)
    tb = annotate_pathogenicity_from_pdfs_normalize_variant_for_search(b)
    if not ta or not tb:
        return False
    set_a = {x for x in ta if len(x) >= 4}
    set_b = {x for x in tb if len(x) >= 4}
    return bool(set_a & set_b) if set_a and set_b else False

def extract_negative_candidates_concrete_variants_in(text: str) -> list[tuple[str, int, int]]:
    out: list[tuple[str, int, int]] = []
    seen: set[str] = set()
    for span in vc_text_module.find_spans(text, vc_text_module.VARIANT_PATTERNS, 'variant', 'scan'):
        if span.subtype not in extract_negative_candidates_CONCRETE_OK:
            continue
        tok = extract_negative_candidates_clean_variant_token(span.text)
        if not tok:
            continue
        key = tok.lower()
        if key in seen:
            continue
        seen.add(key)
        out.append((tok, span.start, span.end))
    return out

def extract_negative_candidates_nearest_variant_to(text: str, label_start: int, label_end: int, *, exclude_terms: list[str], window: int=80) -> str | None:
    candidates = extract_negative_candidates_concrete_variants_in(text)
    mid = (label_start + label_end) // 2
    best: tuple[int, str] | None = None
    for tok, start, end in candidates:
        if any((extract_negative_candidates_is_same_allele(tok, ex) for ex in exclude_terms if ex)):
            continue
        if end < label_start:
            dist = label_start - end
        elif start > label_end:
            dist = start - label_end
        else:
            dist = 0
        if dist > window:
            continue
        if best is None or dist < best[0]:
            best = (dist, tok)
    return best[1] if best else None

def extract_negative_candidates_path_intervening(mid: str) -> bool:
    return bool(re.search('\\b(?:likely\\s+)?pathogenic\\b', mid or '', re.I))

def extract_negative_candidates__acceptable_gene(sym: str | None) -> str | None:
    if not sym:
        return None
    g = str(sym).strip()
    if not g:
        return None
    gu = g.upper()
    if gu in extract_negative_candidates_BAD_GENE_TOKENS or not vc_text_module.is_plausible_gene_symbol(g):
        return None
    return gu

def extract_negative_candidates_guess_gene_near(text: str, variant: str, fallback: str | None) -> str | None:
    """Bind gene for a negative allele; reuse vc_text plausibility + adjacency rules.

    Prefer evidence near this allele; only then fall back to the paper's indexed gene.
    Supports left adjacency (GENE c.xxx) and right ``c.xxx in GENE`` (A2).
    """
    resolved = vc_text_module.resolve_gene_for_variant(variant, evidence_sentence=text, link_gene=None)
    ok = extract_negative_candidates__acceptable_gene(resolved)
    if ok:
        return ok
    v = re.escape(variant or '')
    if v:
        for rx in (f'\\b([A-Z][A-Z0-9]{{1,10}})\\b.{{0,40}}{v}', f'{v}(?:\\s*\\([^)]{{0,40}}\\))?\\s+in\\s+(?:the\\s+)?\\b([A-Z][A-Z0-9]{{1,10}})\\b', f'{v}.{{0,40}}\\b([A-Z][A-Z0-9]{{1,10}})\\b'):
            m = re.search(rx, text or '')
            if not m:
                continue
            ok = extract_negative_candidates__acceptable_gene(m.group(1))
            if ok:
                return ok
    return extract_negative_candidates__acceptable_gene(fallback)

def extract_negative_candidates_variants_in_paren_body(body: str) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for tok, _, _ in extract_negative_candidates_concrete_variants_in(body):
        k = tok.lower()
        if k not in seen:
            seen.add(k)
            out.append(tok)
    return out

def extract_negative_candidates_scan_sentence_for_other_benign(sent: str, *, indexed_variant: str, indexed_gene: str, source: str, allow_loose: bool=True) -> list[extract_negative_candidates_NegativeAlleleHit]:
    hits: list[extract_negative_candidates_NegativeAlleleHit] = []
    if extract_negative_candidates_RE_EXCLUDE.search(sent):
        return hits

    def add(variant: str, label_raw: str, role: str, rule: str, conf: str, *, allow_indexed: bool=False) -> None:
        tok = extract_negative_candidates_clean_variant_token(variant)
        if not tok:
            return
        if not allow_indexed and extract_negative_candidates_is_same_allele(tok, indexed_variant):
            return
        label = extract_negative_candidates_normalize_label(label_raw)
        if label == 'uncertain' and extract_negative_candidates_RE_OVERTURNED_VUS.search(sent):
            return
        hits.append(extract_negative_candidates_NegativeAlleleHit(negative_variant=tok, negative_label=label, negative_role=role, confidence=conf, rule=rule, evidence=sent.strip()[:700], evidence_source=source, negative_gene=extract_negative_candidates_guess_gene_near(sent, tok, indexed_gene)))
    role_base = 'control_variant' if extract_negative_candidates_RE_CONTROL_CUE.search(sent) else 'same_paper_other_site'
    n_before = len(hits)
    for m in extract_negative_candidates_RE_CLASSIFY_PAREN.finditer(sent):
        pre = sent[max(0, m.start() - 30):m.start()]
        if re.search('\\b(?:likely\\s+)?pathogenic\\s*$', pre, re.I):
            continue
        for var in extract_negative_candidates_variants_in_paren_body(m.group('body')):
            add(var, m.group('label'), role_base, 'classify_paren', 'high')
    for m in extract_negative_candidates_RE_VARIANT_RECLASS_LB.finditer(sent):
        add(m.group('var'), m.group('label'), role_base, 'reclassified_lb', 'high')
    for m in re.finditer('(?:classified|reclassified)\\s+as\\s+(?P<label>likely\\s+benign|benign|VUS|a\\s+VUS|a\\s+variant\\s+of\\s+uncertain(?:\\s+clinical)?\\s+significance)', sent, re.I):
        for tok, _, _ in extract_negative_candidates_concrete_variants_in(sent[m.end():m.end() + 80]):
            add(tok, m.group('label'), role_base, 'classify_then_nearby_var', 'high')
            break
    for m in re.finditer('(?P<body>(?:c\\.|p\\.|rs)[^.]{0,120}?)\\b(?:was|were)\\s+classified\\b[^.]{0,60}?(?P<label>likely\\s+benign|benign|VUS|variant\\s+of\\s+uncertain(?:\\s+clinical)?\\s+significance)', sent, re.I):
        for var in extract_negative_candidates_variants_in_paren_body(m.group('body')):
            add(var, m.group('label'), role_base, 'were_classified_lb', 'high')
    for m in re.finditer('(?P<var>(?:c\\.|p\\.|rs)[\\w()*>.\\-+_/ ]{2,40})[^.]{0,80}?(?:classified\\s+as\\s+(?:a\\s+)?(?:VUS|variant\\s+of\\s+uncertain)|was\\s+classified\\s+as\\s+VUS)', sent, re.I):
        add(m.group('var'), 'uncertain', role_base, 'variant_then_vus', 'high')
    if not allow_loose:
        pass
    else:
        for m in extract_negative_candidates_RE_VARIANT_THEN_LABEL.finditer(sent):
            if extract_negative_candidates_path_intervening(m.group('mid')):
                continue
            nearest = extract_negative_candidates_nearest_variant_to(sent, m.start('label'), m.end('label'), exclude_terms=[indexed_variant], window=100)
            if nearest and (not extract_negative_candidates_is_same_allele(nearest, m.group('var'))):
                continue
            add(m.group('var'), m.group('label'), role_base, 'variant_then_label', 'medium')
        for m in extract_negative_candidates_RE_LABEL_THEN_VARIANT.finditer(sent):
            if extract_negative_candidates_path_intervening(m.group('mid')):
                continue
            nearest = extract_negative_candidates_nearest_variant_to(sent, m.start('label'), m.end('label'), exclude_terms=[indexed_variant], window=100)
            if nearest and (not extract_negative_candidates_is_same_allele(nearest, m.group('var'))):
                continue
            add(m.group('var'), m.group('label'), role_base, 'label_then_variant', 'medium')
        if len(hits) == n_before:
            for label_rx, default_label in ((extract_negative_candidates_RE_LIKELY_BENIGN, 'likely_benign'), (extract_negative_candidates_RE_BENIGN, 'benign'), (extract_negative_candidates_RE_VUS, 'uncertain')):
                for m in label_rx.finditer(sent):
                    pre = sent[max(0, m.start() - 25):m.start()]
                    if re.search('\\b(?:likely\\s+)?pathogenic\\b', pre, re.I):
                        continue
                    other = extract_negative_candidates_nearest_variant_to(sent, m.start(), m.end(), exclude_terms=[indexed_variant], window=80)
                    if other:
                        add(other, default_label, role_base, 'nearest_other_to_label', 'low')
        if extract_negative_candidates_RE_CONTROL_CUE.search(sent):
            for tok, _, _ in extract_negative_candidates_concrete_variants_in(sent):
                if not extract_negative_candidates_is_same_allele(tok, indexed_variant):
                    add(tok, 'benign', 'control_variant', 'control_cue_with_variant', 'medium')
    if extract_negative_candidates_RE_NULL_ASSOC.search(sent):
        rs_tokens = [tok for tok, _, _ in extract_negative_candidates_concrete_variants_in(sent) if extract_negative_candidates_RE_RS_ONLY.match(tok)]
        if not rs_tokens:
            iv = extract_negative_candidates_clean_variant_token(indexed_variant)
            if iv and extract_negative_candidates_RE_RS_ONLY.match(iv) and (iv.lower() in sent.lower()):
                rs_tokens = [iv]
        for tok in rs_tokens:
            add(tok, 'association_null', 'association_null_site', 'rs_null_association', 'medium', allow_indexed=True)
    return hits

def extract_negative_candidates_scan_texts_for_other_benign(texts: list[tuple[str, str]], *, indexed_variant: str, indexed_gene: str) -> list[extract_negative_candidates_NegativeAlleleHit]:
    all_hits: list[extract_negative_candidates_NegativeAlleleHit] = []
    for source, text in texts:
        if not (text or '').strip():
            continue
        sents = extract_negative_candidates_sentences(text)
        for i, sent in enumerate(sents):
            all_hits.extend(extract_negative_candidates_scan_sentence_for_other_benign(sent, indexed_variant=indexed_variant, indexed_gene=indexed_gene, source=source, allow_loose=True))
            if i + 1 < len(sents):
                window = f'{sent} {sents[i + 1]}'
                for hit in extract_negative_candidates_scan_sentence_for_other_benign(window, indexed_variant=indexed_variant, indexed_gene=indexed_gene, source=source, allow_loose=False):
                    all_hits.append(extract_negative_candidates_NegativeAlleleHit(negative_variant=hit.negative_variant, negative_label=hit.negative_label, negative_role=hit.negative_role, confidence='medium' if hit.confidence == 'high' else hit.confidence, rule=f'window_{hit.rule}', evidence=hit.evidence, evidence_source=hit.evidence_source, negative_gene=hit.negative_gene))
    return all_hits

def extract_negative_candidates_pick_best_per_allele(hits: list[extract_negative_candidates_NegativeAlleleHit]) -> list[extract_negative_candidates_NegativeAlleleHit]:
    conf_rank = {'high': 0, 'medium': 1, 'low': 2}
    best: dict[str, extract_negative_candidates_NegativeAlleleHit] = {}
    for h in hits:
        key = h.negative_variant.lower()
        prev = best.get(key)
        if prev is None:
            best[key] = h
            continue
        score_new = (extract_negative_candidates_ROLE_PRIORITY.get(h.negative_role, 9), extract_negative_candidates_LABEL_PRIORITY.get(h.negative_label, 9), conf_rank.get(h.confidence, 9))
        score_old = (extract_negative_candidates_ROLE_PRIORITY.get(prev.negative_role, 9), extract_negative_candidates_LABEL_PRIORITY.get(prev.negative_label, 9), conf_rank.get(prev.confidence, 9))
        if score_new < score_old:
            best[key] = h
    return list(best.values())

def extract_negative_candidates_load_abstracts(merged: Path, pmids: set[str]) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    if not merged.exists():
        return out
    for line in merged.open(encoding='utf-8'):
        if not line.strip():
            continue
        r = json.loads(line)
        pmid = str(r.get('pmid') or '')
        if pmid in pmids:
            out[pmid] = {'title': r.get('title') or '', 'abstract': r.get('abstract') or ''}
    return out

def extract_negative_candidates_load_parsed_fulltext(pmid: str, parsed_dir: Path, raw_dir: Path) -> tuple[str, str]:
    path = parsed_dir / f'{pmid}.json'
    if not path.exists():
        return ('', '')
    parsed = json.loads(path.read_text(encoding='utf-8'))
    pmcid = str(parsed.get('pmcid') or '').strip()
    if not pmcid:
        return ('', '')
    passages, source = fulltext_module.load_passages_from_raw(pmcid, parsed, raw_dir)
    if not passages:
        return ('', '')
    parts = [(p.get('text') or '').strip() for p in passages if (p.get('text') or '').strip()]
    return ('\n\n'.join(parts), str(source or parsed.get('source') or 'parsed_fulltext'))

def extract_negative_candidates_build_text_sources(row: dict[str, Any] | None, *, abstract: str, parsed_text: str, pdf_cache_text: str) -> list[tuple[str, str]]:
    texts: list[tuple[str, str]] = []
    if parsed_text.strip():
        texts.append(('parsed_fulltext', parsed_text))
    if abstract.strip():
        texts.append(('abstract', abstract))
    if pdf_cache_text.strip():
        texts.append(('pdf_cache', pdf_cache_text))
    if row:
        ev = str(row.get('evidence_text') or '').strip()
        if ev:
            texts.append(('evidence_text', ev))
    return texts

def extract_negative_candidates_load_positive_exclude_keys(path: Path | None) -> set[str]:
    """gene|allele keys (lower) from a positive review JSONL."""
    if path is None or not path.exists():
        return set()
    keys: set[str] = set()
    with path.open(encoding='utf-8') as f:
        for line in f:
            if not line.strip():
                continue
            r = json.loads(line)
            gene = r.get('gene')
            variant = r.get('variant') or r.get('negative_variant')
            lk = r.get('locus_key') or extract_negative_candidates_locus_key_of(gene, variant)
            if lk:
                keys.add(str(lk).lower())
            if variant:
                keys.add(f'|{str(variant).strip().lower()}')
    return keys

def extract_negative_candidates_remap_role_independent(role: str) -> str:
    if role == 'same_paper_other_site':
        return 'corpus_benign_mention'
    return role

def extract_negative_candidates_collect_pmids_fulltext(*, parsed_dir: Path, articles_pass: Path | None, require_pass: bool) -> list[str]:
    parsed_pmids = {p.stem for p in parsed_dir.glob('*.json')} if parsed_dir.exists() else set()
    if require_pass and articles_pass and articles_pass.exists():
        pass_pmids: set[str] = set()
        with articles_pass.open(encoding='utf-8') as f:
            for line in f:
                if not line.strip():
                    continue
                o = json.loads(line)
                pmid = str(o.get('pmid') or '').strip()
                if pmid:
                    pass_pmids.add(pmid)
        return sorted(parsed_pmids & pass_pmids)
    if require_pass:
        return sorted(parsed_pmids)
    return sorted(parsed_pmids)

def extract_negative_candidates_collect_pmids_from_jsonl(path: Path) -> tuple[list[str], dict[str, list[dict[str, Any]]]]:
    by_pmid: dict[str, list[dict[str, Any]]] = {}
    with path.open(encoding='utf-8') as f:
        for line in f:
            if not line.strip():
                continue
            r = json.loads(line)
            pmid = str(r.get('pmid') or '').strip()
            if pmid:
                by_pmid.setdefault(pmid, []).append(r)
    return (sorted(by_pmid), by_pmid)

def extract_negative_candidates_write_jsonl(path: Path, rows: Iterable[dict[str, Any]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with path.open('w', encoding='utf-8') as f:
        for r in rows:
            f.write(json.dumps(r, ensure_ascii=False) + '\n')
            n += 1
    return n

def extract_negative_candidates_main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pmid-source', choices=('jsonl', 'fulltext', 'full_pass', 'highconf', 'strict'), default='jsonl', help='jsonl: --input anchors (legacy same-paper mode); fulltext: all PMIDs with parsed fulltext; full_pass: articles_pass ∩ parsed fulltext (recommended N-ft); highconf/strict: PMIDs from those review JSONLs as anchors')
    parser.add_argument('--input', type=Path, default=extract_negative_candidates_DEFAULT_INPUT)
    parser.add_argument('--merged', type=Path, default=extract_negative_candidates_DEFAULT_MERGED)
    parser.add_argument('--articles-pass', type=Path, default=extract_negative_candidates_DEFAULT_ARTICLES_PASS)
    parser.add_argument('--parsed-dir', type=Path, default=extract_negative_candidates_DEFAULT_PARSED_DIR)
    parser.add_argument('--raw-dir', type=Path, default=extract_negative_candidates_DEFAULT_RAW_DIR)
    parser.add_argument('--pdf-cache-dir', type=Path, default=extract_negative_candidates_DEFAULT_PDF_CACHE)
    parser.add_argument('--exclude-positives', type=Path, default=None, help='positive JSONL whose gene|allele keys are dropped from negatives')
    parser.add_argument('--no-exclude-positives', action='store_true', help='skip collision filter against positives')
    parser.add_argument('--out-dir', type=Path, default=None)
    parser.add_argument('--out-prefix', type=str, default=None)
    parser.add_argument('--limit', type=int, default=0)
    parser.add_argument('--split-hard-soft', action=argparse.BooleanOptionalAction, default=None, help='write *_hard.jsonl / *_soft.jsonl (default: on for independent pools)')
    args = parser.parse_args()
    independent = args.pmid_source in ('fulltext', 'full_pass')
    if args.split_hard_soft is None:
        args.split_hard_soft = independent
    if args.out_dir is None:
        args.out_dir = extract_negative_candidates_DEFAULT_INDEPENDENT_OUT if independent else args.input.parent
    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.out_prefix is None:
        if independent:
            args.out_prefix = f'negatives_{args.pmid_source}'
        else:
            args.out_prefix = f'{args.input.stem}_other_benign_negatives'
    stem = args.out_prefix
    out_jsonl = args.out_dir / f'{stem}.jsonl'
    out_hits = args.out_dir / f'{stem}_hits.jsonl'
    out_tsv = args.out_dir / f'{stem}.tsv'
    hits_csv = args.out_dir / f'{stem}_hits.csv'
    out_hard = args.out_dir / f'{stem}_hard.jsonl'
    out_soft = args.out_dir / f'{stem}_soft.jsonl'
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    log = logging.getLogger('extract_negative_candidates')
    by_pmid: dict[str, list[dict[str, Any]]] = {}
    if args.pmid_source == 'jsonl':
        _, by_pmid = extract_negative_candidates_collect_pmids_from_jsonl(args.input)
    elif args.pmid_source == 'highconf':
        _, by_pmid = extract_negative_candidates_collect_pmids_from_jsonl(extract_negative_candidates_DEFAULT_HIGHCONF)
    elif args.pmid_source == 'strict':
        _, by_pmid = extract_negative_candidates_collect_pmids_from_jsonl(extract_negative_candidates_DEFAULT_INPUT)
    elif args.pmid_source == 'fulltext':
        for pmid in extract_negative_candidates_collect_pmids_fulltext(parsed_dir=args.parsed_dir, articles_pass=None, require_pass=False):
            by_pmid[pmid] = []
    elif args.pmid_source == 'full_pass':
        for pmid in extract_negative_candidates_collect_pmids_fulltext(parsed_dir=args.parsed_dir, articles_pass=args.articles_pass, require_pass=True):
            by_pmid[pmid] = []
    pmid_list = sorted(by_pmid)
    if args.limit:
        pmid_list = pmid_list[:args.limit]
        by_pmid = {p: by_pmid[p] for p in pmid_list}
    exclude_path = None
    if not args.no_exclude_positives:
        exclude_path = args.exclude_positives or extract_negative_candidates_DEFAULT_EXCLUDE_POS
    exclude_keys = extract_negative_candidates_load_positive_exclude_keys(exclude_path)
    log.info('pmid_source=%s pmids=%s independent=%s exclude_positives=%s keys=%s', args.pmid_source, len(pmid_list), independent, exclude_path, len(exclude_keys))
    abstracts = extract_negative_candidates_load_abstracts(args.merged, set(pmid_list))
    log.info('abstracts loaded=%s', len(abstracts))
    hits_fields = ['pmid', 'negative_gene', 'negative_variant', 'negative_label', 'negative_role', 'negative_confidence', 'negative_rule', 'anchor_gene', 'anchor_variant', 'negative_evidence_source', 'negative_evidence']
    candidates: list[dict[str, Any]] = []
    all_hits_rows: list[dict[str, Any]] = []
    seen_allele: set[tuple[str, str]] = set()
    n_collide = 0
    n_no_text = 0
    for i, pmid in enumerate(pmid_list, 1):
        anchors = by_pmid.get(pmid) or []
        art = abstracts.get(pmid) or {}
        title = art.get('title') or ''
        abstract = art.get('abstract') or ''
        abs_block = (title + '. ' + abstract).strip() if title else abstract
        parsed_text, parsed_meta = extract_negative_candidates_load_parsed_fulltext(pmid, args.parsed_dir, args.raw_dir)
        cache = args.pdf_cache_dir / f'PMID{pmid}.txt'
        pdf_cache_text = cache.read_text(encoding='utf-8', errors='replace') if cache.exists() else ''
        seed_row = anchors[0] if anchors else {}
        texts = extract_negative_candidates_build_text_sources(seed_row or None, abstract=abs_block, parsed_text=parsed_text, pdf_cache_text=pdf_cache_text)
        if independent and (not any((t[1].strip() for t in texts))):
            n_no_text += 1
            continue
        if independent and (not parsed_text.strip()) and (args.pmid_source in ('fulltext', 'full_pass')):
            n_no_text += 1
            continue
        indexed_variants = [str(a.get('variant') or '') for a in anchors] or ['']
        indexed_genes = [str(a.get('gene') or '') for a in anchors] or ['']
        paper_hits: list[extract_negative_candidates_NegativeAlleleHit] = []
        for gene, variant in zip(indexed_genes, indexed_variants):
            paper_hits.extend(extract_negative_candidates_scan_texts_for_other_benign(texts, indexed_variant=variant, indexed_gene=gene))
        filtered: list[extract_negative_candidates_NegativeAlleleHit] = []
        for h in paper_hits:
            role = extract_negative_candidates_remap_role_independent(h.negative_role) if independent else h.negative_role
            h = extract_negative_candidates_NegativeAlleleHit(negative_variant=h.negative_variant, negative_label=h.negative_label, negative_role=role, confidence=h.confidence, rule=h.rule, evidence=h.evidence, evidence_source=h.evidence_source, negative_gene=h.negative_gene)
            if not independent and h.negative_label != 'association_null':
                if any((extract_negative_candidates_is_same_allele(h.negative_variant, iv) for iv in indexed_variants if iv)):
                    continue
            lk = extract_negative_candidates_locus_key_of(h.negative_gene, h.negative_variant).lower()
            bare = f'|{h.negative_variant.lower()}'
            if exclude_keys and (lk in exclude_keys or bare in exclude_keys):
                n_collide += 1
                continue
            filtered.append(h)
        best = extract_negative_candidates_pick_best_per_allele(filtered)
        for h in filtered:
            all_hits_rows.append({'pmid': pmid, 'negative_gene': h.negative_gene, 'negative_variant': h.negative_variant, 'negative_label': h.negative_label, 'negative_role': h.negative_role, 'negative_confidence': h.confidence, 'negative_rule': h.rule, 'anchor_gene': extract_negative_candidates__acceptable_gene(indexed_genes[0]) if indexed_genes and indexed_genes[0] else None, 'anchor_variant': indexed_variants[0] if indexed_variants and indexed_variants[0] else None, 'negative_evidence_source': h.evidence_source, 'negative_evidence': h.evidence})
        for h in best:
            key = (pmid, h.negative_variant.lower())
            if key in seen_allele:
                continue
            seen_allele.add(key)
            chosen = anchors[0] if anchors else {}
            for a in anchors:
                if h.negative_gene and str(a.get('gene') or '') == h.negative_gene:
                    chosen = a
                    break
            out = {'pmid': pmid, 'negative_gene': h.negative_gene, 'negative_variant': h.negative_variant, 'negative_label': h.negative_label, 'negative_role': h.negative_role, 'negative_confidence': h.confidence, 'negative_rule': h.rule, 'negative_evidence': h.evidence, 'negative_evidence_source': h.evidence_source, 'anchor_gene': extract_negative_candidates__acceptable_gene(chosen.get('gene')) if chosen else None, 'anchor_variant': chosen.get('variant') if chosen else None, 'anchor_clinical': chosen.get('clinical') if chosen else None, 'anchor_locus_key': chosen.get('locus_key') if chosen else None, 'negative_locus_key': extract_negative_candidates_locus_key_of(h.negative_gene, h.negative_variant), 'pubmed_url': (chosen.get('pubmed_url') if chosen else None) or f'https://pubmed.ncbi.nlm.nih.gov/{pmid}/', 'text_sources_used': [s for s, _ in texts], 'parsed_fulltext_source': parsed_meta or None, 'design': 'independent_corpus_benign_fulltext' if independent else 'same_paper_other_or_control_benign_plus_vus_assoc_null', 'pmid_source': args.pmid_source}
            candidates.append(out)
        if i % 100 == 0:
            log.info('processed pmids %s/%s candidates=%s hits=%s', i, len(pmid_list), len(candidates), len(all_hits_rows))
    pre_n = len(candidates)
    pre_labels = dict(Counter((c['negative_label'] for c in candidates)))
    pre_roles = dict(Counter((c['negative_role'] for c in candidates)))
    out_raw = args.out_dir / f'{stem}_pre_locus_dedupe.jsonl'
    extract_negative_candidates_write_jsonl(out_raw, candidates)
    candidates, label_counts = extract_negative_candidates_dedupe_by_locus(candidates)
    role_counts = dict(Counter((c['negative_role'] for c in candidates)))
    hard_rows = [c for c in candidates if c['negative_label'] in extract_negative_candidates_HARD_LABELS or c['negative_role'] == 'control_variant']
    soft_rows = [c for c in candidates if c['negative_label'] in extract_negative_candidates_SOFT_LABELS]
    tsv_cols = ['pmid', 'negative_gene', 'negative_variant', 'negative_locus_key', 'negative_label', 'negative_role', 'negative_confidence', 'negative_rule', 'anchor_gene', 'anchor_variant', 'anchor_clinical', 'negative_evidence_source', 'negative_evidence', 'pubmed_url', 'pmid_source', 'design']
    extract_negative_candidates_write_jsonl(out_jsonl, candidates)
    if args.split_hard_soft:
        extract_negative_candidates_write_jsonl(out_hard, hard_rows)
        extract_negative_candidates_write_jsonl(out_soft, soft_rows)
    with out_hits.open('w', encoding='utf-8') as oh:
        for h in all_hits_rows:
            oh.write(json.dumps(h, ensure_ascii=False) + '\n')
    with out_tsv.open('w', encoding='utf-8', newline='') as ot:
        w = csv.DictWriter(ot, fieldnames=tsv_cols, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for c in candidates:
            row = dict(c)
            row['negative_evidence'] = str(row.get('negative_evidence') or '').replace('\n', ' ')
            w.writerow(row)
    with hits_csv.open('w', encoding='utf-8', newline='') as hf:
        w = csv.DictWriter(hf, fieldnames=hits_fields, extrasaction='ignore')
        w.writeheader()
        for h in all_hits_rows:
            row = dict(h)
            row['negative_evidence'] = str(row.get('negative_evidence') or '').replace('\n', ' ')
            w.writerow(row)
    pos_pmids: set[str] = set()
    if exclude_path and exclude_path.exists():
        with exclude_path.open(encoding='utf-8') as f:
            for line in f:
                if line.strip():
                    pos_pmids.add(str(json.loads(line).get('pmid') or ''))
    neg_pmids = {str(c.get('pmid') or '') for c in candidates}
    overlap = len(neg_pmids & pos_pmids)
    summary = {'design': 'independent_corpus_benign_fulltext' if independent else 'same_paper_other_or_control_benign_plus_vus_assoc_null', 'pmid_source': args.pmid_source, 'note': 'Independent mode screens fulltext/pass literature without requiring a positive-sample anchor PMID. Collision-filtered against --exclude-positives. Hard = benign/LB/not_pathogenic/control; soft = VUS/association_null.' if independent else 'Anchor mode: same-paper OTHER / control / association_null relative to indexed positive alleles. Deduped by negative_locus_key.', 'input': str(args.input) if args.pmid_source in ('jsonl', 'strict', 'highconf') else None, 'exclude_positives': str(exclude_path) if exclude_path else None, 'n_pmids_scanned': len(pmid_list), 'n_pmids_no_usable_text': n_no_text, 'n_collide_dropped': n_collide, 'n_candidates_pre_locus_dedupe': pre_n, 'n_candidates': len(candidates), 'n_hard': len(hard_rows), 'n_soft': len(soft_rows), 'n_hits': len(all_hits_rows), 'n_neg_pmids': len(neg_pmids), 'n_pos_pmids_exclude_file': len(pos_pmids), 'n_pmid_overlap_with_positives': overlap, 'label_counts_pre_locus_dedupe': pre_labels, 'role_counts_pre_locus_dedupe': pre_roles, 'label_counts': label_counts, 'role_counts': role_counts, 'outputs': {'jsonl': str(out_jsonl), 'hard_jsonl': str(out_hard) if args.split_hard_soft else None, 'soft_jsonl': str(out_soft) if args.split_hard_soft else None, 'tsv': str(out_tsv), 'pre_locus_dedupe_jsonl': str(out_raw), 'hits_jsonl': str(out_hits), 'hits_csv': str(hits_csv)}}
    (args.out_dir / f'{stem}_summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    log.info('done scanned=%s pre=%s post=%s hard=%s soft=%s collide=%s pmid_overlap=%s', len(pmid_list), pre_n, len(candidates), len(hard_rows), len(soft_rows), n_collide, overlap)
    log.info('outputs: %s', out_jsonl)

# Public helper compatibility aliases.
classify_discovery = tag_discovery_classify_discovery
BAD_GENE_TOKENS = extract_negative_candidates_BAD_GENE_TOKENS
_acceptable_gene = extract_negative_candidates__acceptable_gene
guess_gene_near = extract_negative_candidates_guess_gene_near
normalize_variant_for_search = annotate_pathogenicity_from_pdfs_normalize_variant_for_search
annotate_row = filter_journal_priority_annotate_row
build_indexes = filter_journal_priority_build_indexes

def _dispatch_main(commands: dict[str, object]) -> None:
    import argparse
    import sys
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=sorted(commands))
    args, rest = parser.parse_known_args()
    old = sys.argv
    try:
        sys.argv = [f"{old[0]} {args.command}", *rest]
        commands[args.command]()
    finally:
        sys.argv = old

if __name__ == "__main__":
    _dispatch_main({
        'discovery': tag_discovery_main,
        'filter': filter_review_corpus_main,
        'journal': filter_journal_priority_main,
        'pathogenicity': annotate_pathogenicity_from_pdfs_main,
        'negatives': extract_negative_candidates_main,
    })
