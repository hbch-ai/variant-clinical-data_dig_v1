#!/usr/bin/env python3
"""Physically consolidated pipeline commands. Generated from the v4_01 stage sources.

Each section retains its source module boundary through descriptive global prefixes;
there is no dynamic source loading or embedded executable source text.
"""
from __future__ import annotations
import fulltext as fulltext_module
import genomic as genomic_module
import vc_ensembl_c_map as vc_ensembl_c_map_module

# === run_single_pmid_pipeline.py ===
"""Run Stage0–9 for one PubMed article without touching the production corpus.

Each existing stage CLI is invoked with path overrides only.  Stage logic,
gates, and filters are unchanged.  All intermediates land under:

  data/pubmed/runs/single/<PMID>/

so you can see which gate drops the paper, then edit the corresponding
stage script and re-run with ``--from``.

Examples::

  python scripts/run_single_pmid_pipeline.py 41992294
  python scripts/run_single_pmid_pipeline.py 41992294 --from stage5_clinvar_validation
  python scripts/run_single_pmid_pipeline.py 41992294 --stop-on-drop
"""
import argparse
import importlib.util
import json
import logging
import os
import subprocess
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
run_single_pmid_pipeline_ROOT = Path(__file__).resolve().parents[1]
run_single_pmid_pipeline_SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(run_single_pmid_pipeline_SCRIPTS))
run_single_pmid_pipeline_PYTHON = sys.executable
run_single_pmid_pipeline_PROD_MERGED = run_single_pmid_pipeline_ROOT / 'data' / 'pubmed' / 'ingest' / 'merged' / 'articles.jsonl'
run_single_pmid_pipeline_CLINVAR = run_single_pmid_pipeline_ROOT / 'data' / 'pubmed' / 'ref' / 'clinvar' / 'variant_summary.txt.gz'
run_single_pmid_pipeline_GFF3 = run_single_pmid_pipeline_ROOT / 'data' / 'pubmed' / 'ref' / 'ensembl' / 'Homo_sapiens.GRCh38.115.gff3'
run_single_pmid_pipeline_FASTA = run_single_pmid_pipeline_ROOT / 'data' / 'pubmed' / 'ref' / 'ensembl' / 'Homo_sapiens.GRCh38.dna.toplevel.fa'
run_single_pmid_pipeline_FAI = Path(str(run_single_pmid_pipeline_FASTA) + '.fai')
run_single_pmid_pipeline_JOURNAL_CFG = run_single_pmid_pipeline_ROOT / 'configs' / 'journal_priority.json'
run_single_pmid_pipeline_GATES_CFG = run_single_pmid_pipeline_ROOT / 'configs' / 'pipeline_gates.json'
run_single_pmid_pipeline_STAGE_IDS = ['stage0_preflight', 'stage1_pubmed_retrieval', 'stage2_abstract_linking', 'stage3_fulltext_enrichment', 'stage4_discovery_tagging', 'stage5_clinvar_validation', 'stage6_review_quality', 'stage7_journal_priority', 'stage8_pdf_evidence', 'stage9_genomic_export_qc']

def run_single_pmid_pipeline_n_jsonl(path: Path) -> int:
    if not path.is_file():
        return 0
    n = 0
    with path.open(encoding='utf-8') as fh:
        for line in fh:
            if line.strip():
                n += 1
    return n

def run_single_pmid_pipeline_ensure_jsonl(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_text('', encoding='utf-8')

def run_single_pmid_pipeline_load_pyc(name: str):
    path = run_single_pmid_pipeline_SCRIPTS / f'{name}.pyc'
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise FileNotFoundError(path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def run_single_pmid_pipeline_resolve_script(stem: str) -> Path:
    py = run_single_pmid_pipeline_SCRIPTS / f'{stem}.py'
    pyc = run_single_pmid_pipeline_SCRIPTS / f'{stem}.pyc'
    if py.is_file():
        return py
    if pyc.is_file():
        return pyc
    raise FileNotFoundError(f'missing stage script {stem}.py/.pyc')

def run_single_pmid_pipeline_find_merged_article(pmid: str) -> dict[str, Any] | None:
    if not run_single_pmid_pipeline_PROD_MERGED.is_file():
        return None
    with run_single_pmid_pipeline_PROD_MERGED.open(encoding='utf-8') as fh:
        for line in fh:
            if not line.strip():
                continue
            row = json.loads(line)
            if str(row.get('pmid')) == pmid:
                return row
    return None

class run_single_pmid_pipeline_SinglePmidRunner:

    def __init__(self, pmid: str, workdir: Path) -> None:
        self.pmid = str(pmid).strip()
        self.workdir = workdir
        self.merged = workdir / 'merged' / 'articles.jsonl'
        self.linked = workdir / 'linked'
        self.ft_root = workdir / 'fulltext'
        self.discovery = self.linked / 'discovery'
        self.cv = self.discovery / 'clinvar_novelty_enriched'
        self.pdfs = self.cv / 'pdfs_review_strict'
        self.neg_ind = self.cv / 'negatives_independent'
        self.genomic = self.cv / 'genomic_allele_csv'
        self.logs = workdir / 'logs'
        self.report_path = workdir / 'stage_report.json'
        self.report: dict[str, Any] = {'pmid': self.pmid, 'workdir': str(workdir), 'created_at': datetime.now(timezone.utc).isoformat(), 'production_untouched': [str(run_single_pmid_pipeline_ROOT / 'data' / 'pubmed' / 'ingest' / 'merged'), str(run_single_pmid_pipeline_ROOT / 'data' / 'pubmed' / 'pipeline' / 's2_linking'), str(run_single_pmid_pipeline_ROOT / 'data' / 'pubmed' / 'cache' / 'fulltext')], 'stages': []}
        self.logger = logging.getLogger('single_pmid')

    def layout(self) -> None:
        for d in (self.merged.parent, self.linked, self.ft_root / 'raw', self.ft_root / 'parsed', self.ft_root / 'linked', self.discovery, self.cv, self.pdfs, self.neg_ind, self.genomic, self.logs):
            d.mkdir(parents=True, exist_ok=True)

    def record(self, stage_id: str, *, status: str, n_in: int | None=None, n_out: int | None=None, outputs: dict[str, Any] | None=None, command: list[str] | None=None, note: str='', dropped_reason: str='') -> dict[str, Any]:
        row = {'stage_id': stage_id, 'status': status, 'n_in': n_in, 'n_out': n_out, 'dropped_reason': dropped_reason, 'note': note, 'command': command or [], 'outputs': outputs or {}}
        self.report['stages'].append(row)
        self.flush()
        self.logger.info('%s → %s n_in=%s n_out=%s %s', stage_id, status, n_in, n_out, dropped_reason or note)
        return row

    def flush(self) -> None:
        self.report['updated_at'] = datetime.now(timezone.utc).isoformat()
        self.report_path.write_text(json.dumps(self.report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

    def run_cli(self, stem: str, args: list[str], *, log_name: str, extra_env: dict[str, str] | None=None) -> list[str]:
        commands = {
            'link_variant_phenotype': ('link_variant_phenotype.py', []),
            'enrich_fulltext': ('fulltext.py', ['enrich']),
            'stamp_evidence_source': ('run.py', ['stamp']),
            'tag_discovery': ('review.py', ['discovery']),
            'validate_novelty_clinvar': ('validate_novelty_clinvar.py', []),
            'filter_review_corpus': ('review.py', ['filter']),
            'filter_journal_priority': ('review.py', ['journal']),
            'materialize_review_text_cache': ('fulltext.py', ['materialize']),
            'annotate_pathogenicity_from_pdfs': ('review.py', ['pathogenicity']),
            'extract_negative_candidates': ('review.py', ['negatives']),
            'export_genomic_allele_csv': ('genomic.py', ['export']),
            'validate_grch38_ref_context5': ('extras.py', ['ref-context5']),
            'check_grch38_coord_refalt': ('extras.py', ['refalt']),
        }
        try:
            script_name, prefix = commands[stem]
        except KeyError as exc:
            raise ValueError(f'no consolidated command mapping for {stem}') from exc
        script = run_single_pmid_pipeline_SCRIPTS / script_name
        cmd = [run_single_pmid_pipeline_PYTHON, str(script), *prefix, *args]
        log_path = self.logs / f'{log_name}.log'
        env = os.environ.copy()
        if extra_env:
            env.update(extra_env)
        self.logger.info('exec: %s', ' '.join(cmd))
        with log_path.open('w', encoding='utf-8') as logf:
            logf.write('$ ' + ' '.join(cmd) + '\n')
            logf.flush()
            proc = subprocess.run(cmd, cwd=str(run_single_pmid_pipeline_ROOT), stdout=logf, stderr=subprocess.STDOUT, env=env, check=False)
        if proc.returncode != 0:
            raise RuntimeError(f'{stem} exited {proc.returncode}; see {log_path}')
        return cmd

    def stage0(self) -> None:
        missing = [str(p) for p in (run_single_pmid_pipeline_GATES_CFG, run_single_pmid_pipeline_JOURNAL_CFG, run_single_pmid_pipeline_CLINVAR, run_single_pmid_pipeline_GFF3, run_single_pmid_pipeline_FASTA, run_single_pmid_pipeline_FAI) if not p.exists()]
        linker = run_single_pmid_pipeline_SCRIPTS / 'link_variant_phenotype.pyc'
        if not linker.is_file() and (not (run_single_pmid_pipeline_SCRIPTS / 'link_variant_phenotype.py').is_file()):
            missing.append(str(linker))
        status = 'passed' if not missing else 'error'
        self.record('stage0_preflight', status=status, n_in=0, n_out=0, outputs={'missing': missing}, note='reference + gate files; no production writes')
        if missing:
            raise FileNotFoundError('preflight missing: ' + '; '.join(missing))

    def stage1(self, *, source: str, no_network: bool, force: bool) -> None:
        if self.merged.is_file() and run_single_pmid_pipeline_n_jsonl(self.merged) and (not force):
            rec = json.loads(self.merged.read_text(encoding='utf-8').splitlines()[0])
            self.record('stage1_pubmed_retrieval', status='passed', n_in=1, n_out=1, outputs={'articles': str(self.merged), 'title': rec.get('title')}, note='cached article.jsonl; use --force to refetch')
            return
        article = None
        used = source
        if source in {'auto', 'merged'}:
            article = run_single_pmid_pipeline_find_merged_article(self.pmid)
            if article is not None:
                used = 'merged'
                article = dict(article)
                article.setdefault('source_query_ids', ['single_pmid'])
        if article is None:
            if source == 'merged':
                raise FileNotFoundError(f'PMID {self.pmid} not in {run_single_pmid_pipeline_PROD_MERGED}')
            if no_network:
                raise RuntimeError('need NCBI fetch but --no-network set')
            used = 'efetch'
            article = self._efetch()
        self.merged.parent.mkdir(parents=True, exist_ok=True)
        self.merged.write_text(json.dumps(article, ensure_ascii=False) + '\n', encoding='utf-8')
        (self.workdir / 'article.json').write_text(json.dumps(article, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        self.record('stage1_pubmed_retrieval', status='passed', n_in=1, n_out=1, outputs={'articles': str(self.merged), 'title': article.get('title'), 'journal': article.get('journal'), 'pmcid': article.get('pmcid'), 'source': used}, note=f'single-PMID intake via {used} (not pubmed_pipeline --mode all)')

    def _efetch(self) -> dict[str, Any]:
        pp = run_single_pmid_pipeline_load_pyc('pubmed_pipeline')
        params = {'db': 'pubmed', 'id': self.pmid, 'rettype': 'xml', 'retmode': 'xml', 'retmax': '1', 'email': os.environ.get('NCBI_EMAIL', 'variant-clinical@localhost'), 'tool': 'variant_clinical_v2_single'}
        resp = pp.eutils_get('efetch.fcgi', params, self.logger)
        raw_path = self.workdir / 'raw' / f'efetch_{self.pmid}.xml'
        raw_path.parent.mkdir(parents=True, exist_ok=True)
        raw_path.write_text(resp.text, encoding='utf-8')
        root = ET.fromstring(resp.text)
        articles = list(root.findall('PubmedArticle'))
        if not articles:
            raise RuntimeError(f'EFetch returned no PubmedArticle for PMID {self.pmid}')
        rec = pp.parse_article(articles[0])
        rec['source_query_ids'] = ['single_pmid']
        rec['pubmed_url'] = rec.get('pubmed_url') or f'https://pubmed.ncbi.nlm.nih.gov/{self.pmid}/'
        return rec

    def stage2(self, *, no_network: bool, force: bool) -> None:
        out_pass = self.linked / 'articles_pass.jsonl'
        if out_pass.exists() and run_single_pmid_pipeline_n_jsonl(out_pass) and (not force):
            self.record('stage2_abstract_linking', status='passed', n_in=run_single_pmid_pipeline_n_jsonl(self.merged), n_out=run_single_pmid_pipeline_n_jsonl(out_pass), outputs={'articles_pass': str(out_pass)}, note='cached')
            return
        args = ['--input', str(self.merged), '--out-dir', str(self.linked), '--tiers', 'A', 'B', 'C', '--pubtator-cache', str(self.linked / 'pubtator_parsed.jsonl')]
        if no_network:
            args.append('--skip-pubtator')
        cmd = self.run_cli('link_variant_phenotype', args, log_name='stage2_abstract_linking')
        n_out = run_single_pmid_pipeline_n_jsonl(out_pass)
        dropped = '' if n_out else 'abstract dual-presence / concrete-variant / tier A–C gate'
        self.record('stage2_abstract_linking', status='passed' if n_out else 'dropped', n_in=run_single_pmid_pipeline_n_jsonl(self.merged), n_out=n_out, command=cmd, dropped_reason=dropped, outputs={'articles_pass': str(out_pass), 'pairs': str(self.linked / 'pairs.jsonl'), 'corpus_strict': str(self.linked / 'corpus_strict.jsonl'), 'n_pairs': run_single_pmid_pipeline_n_jsonl(self.linked / 'pairs.jsonl'), 'n_corpus_strict': run_single_pmid_pipeline_n_jsonl(self.linked / 'corpus_strict.jsonl'), 'n_reject': run_single_pmid_pipeline_n_jsonl(self.linked / 'articles_reject.jsonl')})

    def stage3(self, *, force: bool, fulltext_scope: str, no_network: bool=False) -> None:
        corpus = self.linked / 'corpus_enriched.jsonl'
        if corpus.exists() and (not force):
            self.record('stage3_fulltext_enrichment', status='passed' if run_single_pmid_pipeline_n_jsonl(corpus) or run_single_pmid_pipeline_n_jsonl(self.linked / 'articles_pass.jsonl') else 'dropped', n_in=run_single_pmid_pipeline_n_jsonl(self.linked / 'articles_pass.jsonl'), n_out=run_single_pmid_pipeline_n_jsonl(corpus), outputs={'corpus_enriched': str(corpus)}, note='cached')
            return
        args = ['--merged', str(self.merged), '--linked-dir', str(self.linked), '--out-root', str(self.ft_root), '--scope', fulltext_scope, '--only-pmid', self.pmid, '--pubtator-cache', str(self.linked / 'pubtator_parsed.jsonl')]
        if force:
            args.extend(['--rebuild-from-parsed', '--reprocess-links'])
            if not no_network:
                args.append('--force-refetch')
        cmd = self.run_cli('enrich_fulltext', args, log_name='stage3_fulltext_enrichment')
        self.run_cli('stamp_evidence_source', ['--linked-dir', str(self.linked), '--ft-linked-dir', str(self.ft_root / 'linked')], log_name='stage3_stamp_evidence')
        n_enr = run_single_pmid_pipeline_n_jsonl(corpus)
        parsed = self.ft_root / 'parsed' / f'{self.pmid}.json'
        self.record('stage3_fulltext_enrichment', status='passed', n_in=run_single_pmid_pipeline_n_jsonl(self.linked / 'articles_pass.jsonl'), n_out=n_enr, command=cmd, note='scope=' + fulltext_scope + ('; no parsed JSON' if not parsed.exists() else ''), outputs={'corpus_enriched': str(corpus), 'parsed': str(parsed) if parsed.exists() else None, 'n_pairs_fulltext': run_single_pmid_pipeline_n_jsonl(self.linked / 'pairs_fulltext.jsonl')})

    def stage4(self, *, force: bool) -> None:
        novel = self.discovery / 'corpus_enriched_novel.jsonl'
        if novel.exists() and (not force):
            self.record('stage4_discovery_tagging', status='passed' if run_single_pmid_pipeline_n_jsonl(novel) else 'dropped', n_in=run_single_pmid_pipeline_n_jsonl(self.linked / 'corpus_enriched.jsonl'), n_out=run_single_pmid_pipeline_n_jsonl(novel), dropped_reason='' if run_single_pmid_pipeline_n_jsonl(novel) else 'no novel / phenotype_expansion rows', note='cached')
            return
        cmd = self.run_cli('tag_discovery', ['--merged', str(self.merged), '--linked-dir', str(self.linked), '--out-dir', str(self.discovery)], log_name='stage4_discovery_tagging')
        run_single_pmid_pipeline_ensure_jsonl(novel)
        n_out = run_single_pmid_pipeline_n_jsonl(novel)
        self.record('stage4_discovery_tagging', status='passed' if n_out else 'dropped', n_in=run_single_pmid_pipeline_n_jsonl(self.linked / 'corpus_enriched.jsonl'), n_out=n_out, command=cmd, dropped_reason='' if n_out else 'no novel / phenotype_expansion rows (unclear/prior still tagged)', outputs={'tagged': str(self.discovery / 'corpus_enriched_tagged.jsonl'), 'novel': str(novel), 'n_tagged': run_single_pmid_pipeline_n_jsonl(self.discovery / 'corpus_enriched_tagged.jsonl')})

    def stage5(self, *, force: bool, skip_if_empty: bool, clinvar_tracks: str='enriched') -> None:
        absent = self.cv / 'corpus_enriched_novel_clinvar_absent.jsonl'
        novel = self.discovery / 'corpus_enriched_novel.jsonl'
        run_single_pmid_pipeline_ensure_jsonl(novel)
        if absent.exists() and (not force):
            self.record('stage5_clinvar_validation', status='passed' if run_single_pmid_pipeline_n_jsonl(absent) else 'dropped', n_in=run_single_pmid_pipeline_n_jsonl(novel), n_out=run_single_pmid_pipeline_n_jsonl(absent), dropped_reason='' if run_single_pmid_pipeline_n_jsonl(absent) else 'all novel alleles present in ClinVar snapshot', note='cached')
            return
        if skip_if_empty and run_single_pmid_pipeline_n_jsonl(novel) == 0:
            run_single_pmid_pipeline_ensure_jsonl(absent)
            self.record('stage5_clinvar_validation', status='dropped', n_in=0, n_out=0, dropped_reason='empty novel input; ClinVar index not loaded', outputs={'absent': str(absent)})
            return
        tracks = [('corpus_enriched_novel.jsonl', self.cv, 'enriched')]
        if clinvar_tracks == 'all':
            tracks.extend([('corpus_strict_novel.jsonl', self.discovery / 'clinvar_novelty', 'strict'), ('corpus_fulltext_novel.jsonl', self.discovery / 'clinvar_novelty_fulltext', 'fulltext')])
        cmds: list[list[str]] = []
        for fname, out_dir, corpus_flag in tracks:
            src = self.discovery / fname
            run_single_pmid_pipeline_ensure_jsonl(src)
            if run_single_pmid_pipeline_n_jsonl(src) == 0:
                out_dir.mkdir(parents=True, exist_ok=True)
                continue
            cmd = self.run_cli('validate_novelty_clinvar', ['--track', 'pubmed', '--corpus', corpus_flag, '--skip-download', '--input', str(src), '--out-dir', str(out_dir), '--clinvar-file', str(run_single_pmid_pipeline_CLINVAR), '--pubtator-cache', str(self.linked / 'pubtator_parsed.jsonl')], log_name=f'stage5_{out_dir.name}')
            cmds.append(cmd)
        run_single_pmid_pipeline_ensure_jsonl(absent)
        n_out = run_single_pmid_pipeline_n_jsonl(absent)
        self.record('stage5_clinvar_validation', status='passed' if n_out else 'dropped', n_in=run_single_pmid_pipeline_n_jsonl(novel), n_out=n_out, command=cmds[0] if cmds else [], dropped_reason='' if n_out else 'ClinVar-present or empty novel track', outputs={'absent': str(absent), 'annotated': str(self.cv / 'corpus_enriched_novel_clinvar.jsonl'), 'n_clinvar_present': run_single_pmid_pipeline_n_jsonl(self.cv / 'corpus_enriched_novel_clinvar_present.jsonl')})

    def stage6(self, *, force: bool) -> None:
        highconf = self.cv / 'corpus_enriched_novel_clinvar_absent_review_highconf.jsonl'
        strict = self.cv / 'corpus_enriched_novel_clinvar_absent_review_strict.jsonl'
        src = self.cv / 'corpus_enriched_novel_clinvar_absent.jsonl'
        run_single_pmid_pipeline_ensure_jsonl(src)
        if highconf.exists() and (not force):
            self.record('stage6_review_quality', status='passed' if run_single_pmid_pipeline_n_jsonl(highconf) else 'dropped', n_in=run_single_pmid_pipeline_n_jsonl(src), n_out=run_single_pmid_pipeline_n_jsonl(highconf), note='cached')
            return
        cmd = self.run_cli('filter_review_corpus', ['--input', str(src), '--out-dir', str(self.cv)], log_name='stage6_review_quality')
        run_single_pmid_pipeline_ensure_jsonl(highconf)
        run_single_pmid_pipeline_ensure_jsonl(strict)
        n_out = run_single_pmid_pipeline_n_jsonl(highconf)
        self.record('stage6_review_quality', status='passed' if n_out else 'dropped', n_in=run_single_pmid_pipeline_n_jsonl(src), n_out=n_out, command=cmd, dropped_reason='' if n_out else 'failed strict/extended review gates (tier/method/gene/clinical)', outputs={'highconf': str(highconf), 'strict': str(strict), 'n_strict': run_single_pmid_pipeline_n_jsonl(strict)})

    def stage7(self, *, force: bool) -> None:
        high_ann = self.cv / 'corpus_enriched_novel_clinvar_absent_review_highconf_journal_annotated.jsonl'
        strict_in = self.cv / 'corpus_enriched_novel_clinvar_absent_review_strict.jsonl'
        high_in = self.cv / 'corpus_enriched_novel_clinvar_absent_review_highconf.jsonl'
        run_single_pmid_pipeline_ensure_jsonl(strict_in)
        run_single_pmid_pipeline_ensure_jsonl(high_in)
        if high_ann.exists() and (not force):
            self.record('stage7_journal_priority', status='passed', n_in=run_single_pmid_pipeline_n_jsonl(high_in), n_out=run_single_pmid_pipeline_n_jsonl(high_ann), note='cached annotate-only')
            return
        cmds = []
        for src in (strict_in, high_in):
            if run_single_pmid_pipeline_n_jsonl(src) == 0:
                continue
            cmds.append(self.run_cli('filter_journal_priority', ['--input', str(src), '--articles', str(self.merged), '--config', str(run_single_pmid_pipeline_JOURNAL_CFG), '--out-dir', str(self.cv)], log_name=f'stage7_{src.stem[-20:]}'))
        run_single_pmid_pipeline_ensure_jsonl(high_ann)
        self.record('stage7_journal_priority', status='passed', n_in=run_single_pmid_pipeline_n_jsonl(high_in), n_out=run_single_pmid_pipeline_n_jsonl(high_ann), command=cmds[0] if cmds else [], note='annotate-only; journal tier does not drop rows', outputs={'highconf_annotated': str(high_ann)})

    def stage8(self, *, force: bool, no_network: bool) -> None:
        strict = self.cv / 'corpus_enriched_novel_clinvar_absent_review_strict.jsonl'
        high = self.cv / 'corpus_enriched_novel_clinvar_absent_review_highconf_journal_annotated.jsonl'
        run_single_pmid_pipeline_ensure_jsonl(strict)
        neg = self.neg_ind / 'negatives_full_pass.jsonl'
        if neg.exists() and (not force):
            self.record('stage8_pdf_evidence', status='passed', n_in=run_single_pmid_pipeline_n_jsonl(strict), n_out=run_single_pmid_pipeline_n_jsonl(neg), note='cached')
            return
        cmds: list[list[str]] = []
        if run_single_pmid_pipeline_n_jsonl(strict):
            materialize = ['--jsonl', str(strict), '--parsed-dir', str(self.ft_root / 'parsed'), '--raw-dir', str(self.ft_root / 'raw'), '--articles', str(self.merged), '--out-dir', str(self.pdfs / 'pdf_text_cache')]
            if force:
                materialize.append('--force')
            if not no_network:
                materialize.append('--fetch-pdf')
            cmds.append(self.run_cli('materialize_review_text_cache', materialize, log_name='stage8_materialize_text'))
            cmds.append(self.run_cli('annotate_pathogenicity_from_pdfs', ['--input', str(strict), '--pdf-dir', str(self.pdfs), '--out-dir', str(self.pdfs), '--text-cache-dir', str(self.pdfs / 'pdf_text_cache')], log_name='stage8_pdf_pathogenicity'))
            cmds.append(self.run_cli('extract_negative_candidates', ['--pmid-source', 'jsonl', '--input', str(strict), '--merged', str(self.merged), '--articles-pass', str(self.linked / 'articles_pass.jsonl'), '--parsed-dir', str(self.ft_root / 'parsed'), '--raw-dir', str(self.ft_root / 'raw'), '--pdf-cache-dir', str(self.pdfs / 'pdf_text_cache'), '--out-dir', str(self.pdfs), '--out-prefix', 'corpus_strict_other_benign_negatives'], log_name='stage8_neg_jsonl'))
        excl = ['--exclude-positives', str(high if run_single_pmid_pipeline_n_jsonl(high) else strict)]
        if not run_single_pmid_pipeline_n_jsonl(high) and (not run_single_pmid_pipeline_n_jsonl(strict)):
            excl = ['--no-exclude-positives']
        cmds.append(self.run_cli('extract_negative_candidates', ['--pmid-source', 'full_pass', '--merged', str(self.merged), '--articles-pass', str(self.linked / 'articles_pass.jsonl'), '--parsed-dir', str(self.ft_root / 'parsed'), '--raw-dir', str(self.ft_root / 'raw'), '--pdf-cache-dir', str(self.pdfs / 'pdf_text_cache'), '--out-dir', str(self.neg_ind), '--out-prefix', 'negatives_full_pass', *excl], log_name='stage8_neg_full_pass'))
        run_single_pmid_pipeline_ensure_jsonl(neg)
        self.record('stage8_pdf_evidence', status='passed', n_in=run_single_pmid_pipeline_n_jsonl(self.linked / 'articles_pass.jsonl'), n_out=run_single_pmid_pipeline_n_jsonl(neg), command=cmds[-1] if cmds else [], note='PDF patho on strict positives + independent full_pass negatives', outputs={'negatives_full_pass': str(neg), 'n_neg_hard': run_single_pmid_pipeline_n_jsonl(self.neg_ind / 'negatives_full_pass_hard.jsonl'), 'pdf_patho': str(self.pdfs / 'corpus_enriched_novel_clinvar_absent_review_strict_pdf_pathogenicity.jsonl')})

    def stage9(self, *, force: bool) -> None:
        pos = self.cv / 'corpus_enriched_novel_clinvar_absent_review_highconf_journal_annotated.jsonl'
        if run_single_pmid_pipeline_n_jsonl(pos) == 0:
            pos = self.cv / 'corpus_enriched_novel_clinvar_absent_review_strict.jsonl'
        neg = self.neg_ind / 'negatives_full_pass.jsonl'
        poso = self.genomic / 'variants_genomic_positive.csv'
        if poso.exists() and (not force):
            self.record('stage9_genomic_export_qc', status='passed', n_in=run_single_pmid_pipeline_n_jsonl(pos) + run_single_pmid_pipeline_n_jsonl(neg), n_out=run_single_pmid_pipeline_n_jsonl(self.genomic / 'variants_genomic_all.jsonl'), note='cached')
            return
        run_single_pmid_pipeline_ensure_jsonl(pos)
        run_single_pmid_pipeline_ensure_jsonl(neg)
        args = ['--pos', str(pos), '--neg', str(neg), '--clinvar', str(run_single_pmid_pipeline_CLINVAR), '--out-dir', str(self.genomic), '--pdf-cache', str(self.pdfs / 'pdf_text_cache'), '--pubtator', str(self.linked / 'pubtator_parsed.jsonl'), '--dbsnp-cache', str(run_single_pmid_pipeline_ROOT / 'data' / 'pubmed' / 'cache' / 'dbsnp')]
        pos_patho = self.pdfs / 'corpus_enriched_novel_clinvar_absent_review_strict_pdf_pathogenicity.jsonl'
        if pos_patho.exists():
            args.extend(['--pos-patho', str(pos_patho)])
        cmd = self.run_cli('export_genomic_allele_csv', args, log_name='stage9_export')
        all_csv = self.genomic / 'variants_genomic_all.csv'
        checked = self.genomic / 'variants_genomic_all_ref_context5_with_fixops_applied.csv'
        if all_csv.exists() and run_single_pmid_pipeline_n_jsonl(self.genomic / 'variants_genomic_all.jsonl'):
            self.run_cli('validate_grch38_ref_context5', ['--in-csv', str(all_csv), '--out-csv', str(checked), '--upstream', '5', '--downstream', '5'], log_name='stage9_ref_context5')
            ok_csv = self.genomic / 'variants_genomic_ok.csv'
            if ok_csv.exists():
                tsv = self.genomic / 'variants_genomic_ok_coord_refalt_check.tsv'
                self.run_cli('check_grch38_coord_refalt', ['--in-csv', str(ok_csv), '--out-tsv', str(tsv)], log_name='stage9_coord_refalt')
        n_all = run_single_pmid_pipeline_n_jsonl(self.genomic / 'variants_genomic_all.jsonl')
        self.record('stage9_genomic_export_qc', status='passed' if n_all else 'dropped', n_in=run_single_pmid_pipeline_n_jsonl(pos) + run_single_pmid_pipeline_n_jsonl(neg), n_out=n_all, command=cmd, dropped_reason='' if n_all else 'no positive/negative alleles to export', outputs={'all_csv': str(all_csv), 'ok_csv': str(self.genomic / 'variants_genomic_ok.csv'), 'positive_csv': str(self.genomic / 'variants_genomic_positive.csv'), 'negative_csv': str(self.genomic / 'variants_genomic_negative.csv'), 'n_ok': run_single_pmid_pipeline_n_jsonl(self.genomic / 'variants_genomic_ok.jsonl'), 'n_positive': run_single_pmid_pipeline_n_jsonl(self.genomic / 'variants_genomic_positive.jsonl'), 'n_negative': run_single_pmid_pipeline_n_jsonl(self.genomic / 'variants_genomic_negative.jsonl')})

    def drop_summary(self) -> str | None:
        for row in self.report['stages']:
            if row.get('status') == 'dropped':
                return f"{row['stage_id']}: {row.get('dropped_reason') or 'empty output'}"
        return None

def run_single_pmid_pipeline_main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pmid', help='numeric PubMed ID')
    parser.add_argument('--workdir', type=Path, default=None, help='default: data/pubmed/runs/single/<PMID>/')
    parser.add_argument('--from', dest='from_stage', choices=run_single_pmid_pipeline_STAGE_IDS, default='stage0_preflight')
    parser.add_argument('--to', dest='to_stage', choices=run_single_pmid_pipeline_STAGE_IDS, default='stage9_genomic_export_qc')
    parser.add_argument('--source', choices=['auto', 'fetch', 'merged'], default='auto', help='auto: reuse production merged record if present, else NCBI EFetch')
    parser.add_argument('--fulltext-scope', choices=['pass', 'all'], default='pass', help='same as enrich_fulltext --scope (production default: pass)')
    parser.add_argument('--force', action='store_true', help='rebuild stages even if cached')
    parser.add_argument('--no-network', action='store_true')
    parser.add_argument('--stop-on-drop', action='store_true', help='halt after the first filter stage with n_out=0')
    parser.add_argument('--always-clinvar', action='store_true', help='load ClinVar even when novel track is empty')
    parser.add_argument('--clinvar-tracks', choices=['enriched', 'all'], default='enriched', help='enriched = Stage6+ input only (default). all = also reload ClinVar for strict/fulltext tracks')
    args = parser.parse_args()
    pmid = str(args.pmid).strip()
    if not pmid.isdigit():
        parser.error('pmid must be numeric')
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    workdir = args.workdir or run_single_pmid_pipeline_ROOT / 'data' / 'pubmed' / 'runs' / 'single' / pmid
    workdir = workdir.resolve()
    runner = run_single_pmid_pipeline_SinglePmidRunner(pmid, workdir)
    runner.layout()
    start = run_single_pmid_pipeline_STAGE_IDS.index(args.from_stage)
    end = run_single_pmid_pipeline_STAGE_IDS.index(args.to_stage)
    if end < start:
        parser.error('--to must be at or after --from')
    fns = {'stage0_preflight': lambda: runner.stage0(), 'stage1_pubmed_retrieval': lambda: runner.stage1(source=args.source, no_network=args.no_network, force=args.force), 'stage2_abstract_linking': lambda: runner.stage2(no_network=args.no_network, force=args.force), 'stage3_fulltext_enrichment': lambda: runner.stage3(force=args.force, fulltext_scope=args.fulltext_scope, no_network=args.no_network), 'stage4_discovery_tagging': lambda: runner.stage4(force=args.force), 'stage5_clinvar_validation': lambda: runner.stage5(force=args.force, skip_if_empty=not args.always_clinvar, clinvar_tracks=args.clinvar_tracks), 'stage6_review_quality': lambda: runner.stage6(force=args.force), 'stage7_journal_priority': lambda: runner.stage7(force=args.force), 'stage8_pdf_evidence': lambda: runner.stage8(force=args.force, no_network=args.no_network), 'stage9_genomic_export_qc': lambda: runner.stage9(force=args.force)}
    current_sid = args.from_stage
    try:
        for i, sid in enumerate(run_single_pmid_pipeline_STAGE_IDS):
            if i < start or i > end:
                continue
            current_sid = sid
            fns[sid]()
            last = runner.report['stages'][-1]
            if last.get('status') == 'error':
                break
            if args.stop_on_drop and last.get('status') == 'dropped':
                runner.logger.warning('stop-on-drop at %s', sid)
                break
    except Exception as exc:
        last_status = runner.report['stages'][-1]['status'] if runner.report['stages'] else None
        if last_status != 'error':
            runner.record(current_sid, status='error', note=str(exc))
        runner.flush()
        raise
    drop = runner.drop_summary()
    runner.report['first_drop'] = drop
    runner.flush()
    print(json.dumps({'pmid': pmid, 'workdir': str(workdir), 'report': str(runner.report_path), 'first_drop': drop, 'stages': [{'id': s['stage_id'], 'status': s['status'], 'n_in': s['n_in'], 'n_out': s['n_out']} for s in runner.report['stages']]}, ensure_ascii=False, indent=2))

# === download_review_strict_pdfs.py ===
"""Download or materialize PDFs for strict-review PMIDs.

Order:
  1. Keep existing valid local PDF
  2. Unpaywall OA PDF (by DOI)
  3. Europe PMC fullTextUrlList PDF
  4. Fallback: render PDF from local BioC / merged abstract

Output:
  data/pubmed/cache/pdfs/
"""
import argparse
import json
import logging
import re
import sys
import time
from pathlib import Path
from typing import Any
from urllib.parse import quote
import requests
from fpdf import FPDF
download_review_strict_pdfs_ROOT = Path(__file__).resolve().parents[1]
download_review_strict_pdfs_SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(download_review_strict_pdfs_SCRIPTS))
download_review_strict_pdfs_DEFAULT_JSONL = download_review_strict_pdfs_ROOT / 'data/pubmed/pipeline/s6_review' / 'corpus_enriched_novel_clinvar_absent_review_strict.jsonl'
download_review_strict_pdfs_DEFAULT_MERGED = download_review_strict_pdfs_ROOT / 'data/pubmed/ingest/merged/articles.jsonl'
download_review_strict_pdfs_DEFAULT_OUT = download_review_strict_pdfs_ROOT / 'data/pubmed/cache/pdfs'
download_review_strict_pdfs_FONT_CANDIDATES = [Path('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'), Path('/usr/share/fonts/truetype/freefont/FreeSans.ttf')]
download_review_strict_pdfs_SESSION = requests.Session()
download_review_strict_pdfs_SESSION.headers.update({'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/121.0.0.0 Safari/537.36', 'Accept': 'application/pdf,text/html,*/*'})

def download_review_strict_pdfs_load_pmids(jsonl: Path) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for line in jsonl.open(encoding='utf-8'):
        pmid = str(json.loads(line).get('pmid') or '').strip()
        if pmid and pmid not in seen:
            seen.add(pmid)
            out.append(pmid)
    return out

def download_review_strict_pdfs_load_articles(merged: Path, pmids: set[str]) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for line in merged.open(encoding='utf-8'):
        r = json.loads(line)
        pmid = str(r.get('pmid') or '')
        if pmid in pmids:
            out[pmid] = r
    return out

def download_review_strict_pdfs_is_pdf_bytes(blob: bytes) -> bool:
    return bool(blob) and blob.startswith(b'%PDF') and (len(blob) > 1000)

def download_review_strict_pdfs_safe_get(url: str, *, timeout: int=90) -> tuple[int, bytes | None]:
    try:
        resp = download_review_strict_pdfs_SESSION.get(url, timeout=timeout, allow_redirects=True)
        return (resp.status_code, resp.content)
    except Exception:
        return (0, None)

def download_review_strict_pdfs_try_urls(urls: list[str], dest: Path, logger: logging.Logger, label: str) -> str | None:
    seen: set[str] = set()
    for url in urls:
        if not url or url in seen:
            continue
        seen.add(url)
        time.sleep(0.2)
        code, blob = download_review_strict_pdfs_safe_get(url)
        if code == 200 and blob and download_review_strict_pdfs_is_pdf_bytes(blob):
            dest.write_bytes(blob)
            logger.info('%s saved %s (%s bytes) from %s', label, dest.name, dest.stat().st_size, url)
            return f'{label}:{url}'
    return None

def download_review_strict_pdfs_try_unpaywall(doi: str, dest: Path, logger: logging.Logger, email: str) -> str | None:
    if not doi:
        return None
    code, blob = download_review_strict_pdfs_safe_get(f'https://api.unpaywall.org/v2/{quote(doi)}?email={quote(email)}', timeout=45)
    if code != 200 or not blob:
        return None
    try:
        data = json.loads(blob.decode('utf-8', errors='replace'))
    except json.JSONDecodeError:
        return None
    urls: list[str] = []
    for loc in [data.get('best_oa_location') or {}] + list(data.get('oa_locations') or []):
        if not isinstance(loc, dict):
            continue
        for key in ('url_for_pdf', 'url'):
            u = loc.get(key)
            if u:
                urls.append(u)
    return download_review_strict_pdfs_try_urls(urls, dest, logger, 'unpaywall')

def download_review_strict_pdfs_try_ncbi_pmc_pdf(pmcid: str | None, dest: Path, logger: logging.Logger) -> str | None:
    pmcid = fulltext_module.normalize_pmcid(pmcid or '')
    if not pmcid:
        return None
    urls = [f'https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/pdf/', f'https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/pdf/']
    return download_review_strict_pdfs_try_urls(urls, dest, logger, 'pmc_oa')

def download_review_strict_pdfs_try_europepmc(pmid: str, dest: Path, logger: logging.Logger) -> str | None:
    api = f'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{quote(pmid)}%20AND%20SRC:MED&resultType=core&format=json&pageSize=1'
    code, blob = download_review_strict_pdfs_safe_get(api, timeout=45)
    if code != 200 or not blob:
        return None
    try:
        data = json.loads(blob.decode('utf-8', errors='replace'))
    except json.JSONDecodeError:
        return None
    results = (data.get('resultList') or {}).get('result') or []
    if not results:
        return None
    urls: list[str] = []
    for item in (results[0].get('fullTextUrlList') or {}).get('fullTextUrl') or []:
        if not isinstance(item, dict):
            continue
        u = item.get('url') or ''
        style = (item.get('documentStyle') or '').lower()
        if not u:
            continue
        if style == 'pdf' or u.lower().endswith('.pdf') or 'pdf' in u.lower():
            urls.append(u)
    pmcid = fulltext_module.normalize_pmcid(results[0].get('pmcid') or '')
    if pmcid:
        urls.append(f'https://europepmc.org/articles/{pmcid}?pdf=render')
    return download_review_strict_pdfs_try_urls(urls, dest, logger, 'europepmc')

def download_review_strict_pdfs_bioc_text(pmcid: str | None) -> str:
    pmcid = fulltext_module.normalize_pmcid(pmcid or '')
    if not pmcid:
        return ''
    path = fulltext_module.RAW_DIR / f'{pmcid}.bioc.json'
    if not path.exists():
        return ''
    try:
        data = json.loads(path.read_text(encoding='utf-8'))
    except Exception:
        return ''
    parts: list[str] = []
    for doc in data.get('documents') or []:
        for passage in doc.get('passages') or []:
            text = (passage.get('text') or '').strip()
            if text:
                parts.append(text)
    return '\n\n'.join(parts)

def download_review_strict_pdfs_render_text_pdf(dest: Path, *, title: str, pmid: str, doi: str, body: str, note: str) -> None:
    font = next((p for p in download_review_strict_pdfs_FONT_CANDIDATES if p.exists()), None)
    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()
    if font:
        pdf.add_font('Body', '', str(font))
        pdf.set_font('Body', size=11)
    else:
        pdf.set_font('Helvetica', size=11)
    header = f'PMID: {pmid}'
    if doi:
        header += f'\nDOI: {doi}'
    header += f'\nSource note: {note}\n'
    for block in (header, title or '(no title)', '', body or '(no full text available locally)'):
        clean = re.sub('[\\x00-\\x08\\x0b\\x0c\\x0e-\\x1f]', ' ', str(block))
        try:
            pdf.multi_cell(0, 5, clean)
        except Exception:
            pdf.multi_cell(0, 5, clean.encode('latin-1', 'replace').decode('latin-1'))
        pdf.ln(2)
    dest.parent.mkdir(parents=True, exist_ok=True)
    pdf.output(str(dest))

def download_review_strict_pdfs_fallback_local_pdf(pmid: str, art: dict[str, Any], dest: Path) -> str:
    title = art.get('title') or ''
    doi = art.get('doi') or ''
    pmcid = art.get('pmcid')
    body = download_review_strict_pdfs_bioc_text(pmcid)
    if body:
        note = f'local BioC full text ({fulltext_module.normalize_pmcid(pmcid)}); not publisher PDF'
        download_review_strict_pdfs_render_text_pdf(dest, title=title, pmid=pmid, doi=doi, body=body, note=note)
        return 'local_bioc'
    abstract = art.get('abstract') or ''
    note = 'local abstract only; publisher/OA PDF unavailable from this host'
    download_review_strict_pdfs_render_text_pdf(dest, title=title, pmid=pmid, doi=doi, body=abstract, note=note)
    return 'local_abstract'

def download_review_strict_pdfs_main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jsonl', type=Path, default=download_review_strict_pdfs_DEFAULT_JSONL)
    parser.add_argument('--merged', type=Path, default=download_review_strict_pdfs_DEFAULT_MERGED)
    parser.add_argument('--out-dir', type=Path, default=download_review_strict_pdfs_DEFAULT_OUT)
    parser.add_argument('--email', default='variant-clinical@localhost')
    parser.add_argument('--limit', type=int, default=0)
    parser.add_argument('--sleep', type=float, default=0.25)
    parser.add_argument('--publisher-only', action='store_true', help='do not render local BioC/abstract fallback PDFs')
    parser.add_argument('--pmid-list', type=Path, default=None, help='Optional file with one PMID per line; intersect with --jsonl PMIDs')
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    logger = logging.getLogger('pdf_download')
    pmids = download_review_strict_pdfs_load_pmids(args.jsonl)
    if args.pmid_list and args.pmid_list.is_file():
        allow = {ln.strip() for ln in args.pmid_list.read_text(encoding='utf-8').splitlines() if ln.strip()}
        pmids = [p for p in pmids if p in allow]
    if args.limit:
        pmids = pmids[:args.limit]
    articles = download_review_strict_pdfs_load_articles(args.merged, set(pmids))
    logger.info('PMIDs=%s articles=%s out=%s', len(pmids), len(articles), args.out_dir)
    summary: dict[str, Any] = {'n_pmids': len(pmids), 'ok_publisher': 0, 'ok_local_bioc': 0, 'ok_local_abstract': 0, 'fail': 0, 'by_source': {}, 'failures': []}
    log_path = args.out_dir / 'download_log.jsonl'
    with log_path.open('w', encoding='utf-8') as logf:
        for i, pmid in enumerate(pmids, 1):
            art = articles.get(pmid) or {}
            dest = args.out_dir / f'PMID{pmid}.pdf'
            row: dict[str, Any] = {'pmid': pmid, 'doi': art.get('doi'), 'pmcid': art.get('pmcid'), 'path': str(dest)}
            source = None
            if dest.exists() and download_review_strict_pdfs_is_pdf_bytes(dest.read_bytes()):
                size = dest.stat().st_size
                if size > 80000:
                    row['status'] = 'exists'
                    row['source'] = 'local_existing'
                    row['bytes'] = size
                    summary['ok_publisher'] += 1
                    summary['by_source']['local_existing'] = summary['by_source'].get('local_existing', 0) + 1
                    logf.write(json.dumps(row, ensure_ascii=False) + '\n')
                    logger.info('[%s/%s] PMID %s -> exists', i, len(pmids), pmid)
                    continue
            try:
                source = download_review_strict_pdfs_try_unpaywall(str(art.get('doi') or ''), dest, logger, args.email)
                if not source:
                    source = download_review_strict_pdfs_try_ncbi_pmc_pdf(art.get('pmcid'), dest, logger)
                if not source:
                    source = download_review_strict_pdfs_try_europepmc(pmid, dest, logger)
                if not source and (not args.publisher_only):
                    source = download_review_strict_pdfs_fallback_local_pdf(pmid, art, dest)
            except Exception as exc:
                row['error'] = str(exc)
                source = None
            if source and dest.exists() and download_review_strict_pdfs_is_pdf_bytes(dest.read_bytes()):
                row['status'] = 'ok'
                row['source'] = source
                row['bytes'] = dest.stat().st_size
                key = source.split(':', 1)[0]
                summary['by_source'][key] = summary['by_source'].get(key, 0) + 1
                if key in {'unpaywall', 'europepmc', 'pmc_oa'}:
                    summary['ok_publisher'] += 1
                elif key == 'local_bioc':
                    summary['ok_local_bioc'] += 1
                else:
                    summary['ok_local_abstract'] += 1
            else:
                row['status'] = 'fail'
                summary['fail'] += 1
                summary['failures'].append(pmid)
                if dest.exists() and (not download_review_strict_pdfs_is_pdf_bytes(dest.read_bytes())):
                    dest.unlink(missing_ok=True)
            logf.write(json.dumps(row, ensure_ascii=False) + '\n')
            logf.flush()
            logger.info('[%s/%s] PMID %s -> %s (%s)', i, len(pmids), pmid, row['status'], row.get('source'))
            time.sleep(args.sleep)
    (args.out_dir / 'download_summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    with (args.out_dir / 'manifest.tsv').open('w', encoding='utf-8') as fh:
        fh.write('pmid\tstatus\tsource\tbytes\tpath\tdoi\tpmcid\n')
        for line in log_path.open(encoding='utf-8'):
            r = json.loads(line)
            fh.write(f"{r.get('pmid')}\t{r.get('status')}\t{r.get('source', '')}\t{r.get('bytes', '')}\t{r.get('path', '')}\t{r.get('doi') or ''}\t{r.get('pmcid') or ''}\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))

# === filter_new_site_variants.py ===
"""Filter Stage9 sets to newly reported variant sites (drop rs-only).

Keep:
  - non-rs variant_text (c./p./NM:/genomic literals)
  - rs rows whose literature writes a parenthetical pair:
      c.123A>G (rs…)  or  rs… (c.123A>G)
    including NM_/NR_:c. … (rs…)

Do not keep merely because ClinVar/dbSNP filled c_hgvs on an rs row.

Original Stage9 JSONL/CSV and audit_annotations.json are not deleted.
Filtered products go to s9_genomic/new_site/; review HTML is rebuilt with
existing annotations overlaid (audit_id first, then pmid+gene+variant_text).
"""
import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable
filter_new_site_variants_ROOT = Path(__file__).resolve().parents[1]
filter_new_site_variants_SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(filter_new_site_variants_SCRIPTS))
filter_new_site_variants_S9 = filter_new_site_variants_ROOT / 'data/pubmed/pipeline/s9_genomic'
filter_new_site_variants_MANUAL = filter_new_site_variants_S9 / 'manual_review'
filter_new_site_variants_TRAILS = filter_new_site_variants_S9 / 'audit_trails'
filter_new_site_variants_NEW_SITE = filter_new_site_variants_S9 / 'new_site'
filter_new_site_variants_RE_RS = re.compile('^rs\\d+$', re.I)
filter_new_site_variants_RE_NM_C_RS = re.compile('(N[MR]_\\d+(?:\\.\\d+)?:c\\.[^()\\s]{2,60})\\s*[\\(（]\\s*(rs\\d{4,})\\s*[\\)）]', re.I)
filter_new_site_variants_CODING = {'coding'}
filter_new_site_variants_NONCODING = {'utr', 'intron', 'ncrna_exon', 'ncrna_gene', 'pseudogene', 'intergenic'}
filter_new_site_variants_PASS_CSV_COLS = ['audit_id', 'pmid', 'gene', 'variant_text', 'human_verdict', 'human_genomic_region', 'region_used', 'genomic_region', 'genomic_region_feature', 'genomic_region_genes', 'genomic_region_note', 'chromosome_ucsc', 'position', 'ref', 'alt', 'c_hgvs', 'p_hgvs', 'nt_pro_only', 'resolve_status', 'resolve_method', 'link_tier', 'review_tier', 'screen_grade', 'label', 'pubmed_url', 'human_comment', 'human_fail_stage', 'updated_at', 'c_hgvs_strand', 'chromosome', 'chromosome_accession', 'assembly', 'transcript_id', 'locus_key_src', 'allele_query_source', 'resolve_confidence', 'resolve_note', 'genomic_region_transcripts', 'sample_class', 'rs_keep_reason', 'literature_paired_c_hgvs']

def filter_new_site_variants_load_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open(encoding='utf-8') as f:
        for line in f:
            if line.strip():
                rows.append(json.loads(line))
    return rows

def filter_new_site_variants_write_jsonl(path: Path, rows: Iterable[dict[str, Any]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with path.open('w', encoding='utf-8') as f:
        for row in rows:
            f.write(json.dumps(row, ensure_ascii=False) + '\n')
            n += 1
    return n

def filter_new_site_variants_load_jsonl_by_pmid(path: Path) -> dict[str, list[dict[str, Any]]]:
    out: dict[str, list[dict[str, Any]]] = defaultdict(list)
    if not path.exists():
        return out
    with path.open(encoding='utf-8') as f:
        for line in f:
            if not line.strip():
                continue
            rec = json.loads(line)
            out[str(rec.get('pmid') or '')].append(rec)
    return out

def filter_new_site_variants_rec_key(row: dict[str, Any]) -> tuple[str, str, str, str]:
    return (str(row.get('sample_class') or '').strip().lower(), str(row.get('pmid') or '').strip(), str(row.get('gene') or '').strip().upper(), str(row.get('variant_text') or '').strip())

def filter_new_site_variants_literature_blob(pmid: str, *, articles: dict[str, dict[str, Any]], sources: list[dict[str, list[dict[str, Any]]]]) -> str:
    parts: list[str] = []
    art = articles.get(str(pmid)) or {}
    parts.append(f"{art.get('title') or ''} {art.get('abstract') or ''}")
    for src in sources:
        for rec in src.get(str(pmid) or '') or []:
            for key in ('evidence_text', 'text', 'context_text', 'abstract', 'title', 'snippet', 'passage'):
                val = rec.get(key)
                if val:
                    parts.append(str(val))
    return '\n'.join(parts)

def filter_new_site_variants_paired_c_from_literature(blob: str, rs_id: str) -> str | None:
    hit = genomic_module.find_paired_c_for_rs(blob, rs_id)
    if hit:
        return hit
    rs_n = rs_id.lower()
    for m in filter_new_site_variants_RE_NM_C_RS.finditer(blob or ''):
        if m.group(2).lower() == rs_n:
            return m.group(1)
    return None

def filter_new_site_variants_decide_row(row: dict[str, Any], blob: str) -> tuple[bool, str, str | None]:
    vt = str(row.get('variant_text') or '').strip()
    if not filter_new_site_variants_RE_RS.match(vt):
        return (True, 'non_rs', None)
    paired = filter_new_site_variants_paired_c_from_literature(blob, vt)
    if paired:
        return (True, 'literature_c_rs_pair', paired)
    return (False, 'rs_only', None)

def filter_new_site_variants_filter_table(rows: list[dict[str, Any]], *, articles: dict[str, dict[str, Any]], sources: list[dict[str, list[dict[str, Any]]]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], Counter]:
    kept: list[dict[str, Any]] = []
    dropped: list[dict[str, Any]] = []
    stats: Counter = Counter()
    blob_cache: dict[str, str] = {}
    for row in rows:
        pmid = str(row.get('pmid') or '')
        if pmid not in blob_cache:
            blob_cache[pmid] = filter_new_site_variants_literature_blob(pmid, articles=articles, sources=sources)
        keep, reason, paired_c = filter_new_site_variants_decide_row(row, blob_cache[pmid])
        stats[reason] += 1
        out = dict(row)
        out['rs_keep_reason'] = reason
        out['literature_paired_c_hgvs'] = paired_c or ''
        if keep and paired_c and (not str(out.get('c_hgvs') or '').strip()):
            out['c_hgvs'] = paired_c
        if keep:
            kept.append(out)
        else:
            dropped.append(out)
    return (kept, dropped, stats)

def filter_new_site_variants_write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None=None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not fieldnames:
        keys: list[str] = []
        seen: set[str] = set()
        for row in rows:
            for k in row:
                if k not in seen:
                    seen.add(k)
                    keys.append(k)
        fieldnames = keys
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, '') for k in fieldnames})

def filter_new_site_variants_filter_chains(chain_path: Path, keep_keys: Counter) -> list[dict[str, Any]]:
    remaining = Counter(keep_keys)
    kept: list[dict[str, Any]] = []
    if not chain_path.exists():
        return kept
    with chain_path.open(encoding='utf-8') as f:
        for line in f:
            if not line.strip():
                continue
            rec = json.loads(line)
            key = filter_new_site_variants_rec_key(rec)
            matched = None
            if remaining[key] > 0:
                matched = key
            else:
                for cand in list(remaining):
                    if remaining[cand] <= 0:
                        continue
                    if cand[1] == key[1] and cand[2] == key[2] and (cand[3] == key[3]):
                        matched = cand
                        break
            if matched is None:
                continue
            remaining[matched] -= 1
            kept.append(rec)
    return kept

def filter_new_site_variants_split_chains(rows: list[dict[str, Any]], out_dir: Path, stem: str) -> dict[str, int]:
    buckets: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for rec in rows:
        tier = str(rec.get('link_tier') or '').strip().upper()
        if tier not in {'A', 'B', 'C'}:
            tier = tier or 'NONE'
        buckets[tier].append(rec)
    counts = {}
    for tier, recs in sorted(buckets.items()):
        path = out_dir / f'{stem}_tier_{tier}.jsonl'
        filter_new_site_variants_write_jsonl(path, recs)
        counts[tier] = len(recs)
    return counts

def filter_new_site_variants_region_used(row: dict[str, Any], ann: dict[str, Any] | None) -> str:
    human = str((ann or {}).get('human_genomic_region') or '').strip()
    if human:
        return human
    return str(row.get('genomic_region') or '').strip()

def filter_new_site_variants_export_pass_csv(*, jsonl_rows: list[dict[str, Any]], chains: list[dict[str, Any]], annotations: list[dict[str, Any]], out_csv: Path, allowed_regions: set[str]) -> dict[str, Any]:
    by_id = {a.get('audit_id'): a for a in annotations if a.get('audit_id')}
    by_key: dict[tuple, list[dict[str, Any]]] = defaultdict(list)
    for a in annotations:
        by_key[filter_new_site_variants_rec_key(a)].append(a)
        by_key['', str(a.get('pmid') or ''), str(a.get('gene') or '').upper(), str(a.get('variant_text') or '').strip()].append(a)
    chain_q: dict[tuple, list[dict[str, Any]]] = defaultdict(list)
    for rec in chains:
        chain_q[filter_new_site_variants_rec_key(rec)].append(rec)
    used_ann: set[str] = set()
    exported: list[dict[str, Any]] = []
    verdicts = Counter()
    for row in jsonl_rows:
        key = filter_new_site_variants_rec_key(row)
        chain = (chain_q[key] or [None]).pop(0) if chain_q[key] else None
        if chain_q[key] == []:
            pass
        audit_id = str((chain or {}).get('audit_id') or '')
        ann = by_id.get(audit_id) if audit_id else None
        if not ann or not (ann.get('human_verdict') or ann.get('human_comment') or ann.get('human_genomic_region')):
            cands = by_key.get(key) or by_key.get(('', key[1], key[2], key[3])) or []
            for cand in cands:
                cid = str(cand.get('audit_id') or '')
                if cid in used_ann:
                    continue
                ann = cand
                break
        if not ann:
            continue
        used_ann.add(str(ann.get('audit_id') or ''))
        verd = str(ann.get('human_verdict') or '').strip().lower()
        verdicts[verd or 'pending'] += 1
        if verd != 'pass':
            continue
        used = filter_new_site_variants_region_used(row, ann)
        if used not in allowed_regions:
            continue
        rec = dict(row)
        rec.update({'audit_id': audit_id or ann.get('audit_id') or '', 'human_verdict': ann.get('human_verdict') or '', 'human_genomic_region': ann.get('human_genomic_region') or '', 'human_comment': ann.get('human_comment') or '', 'human_fail_stage': ann.get('human_fail_stage') or '', 'updated_at': ann.get('updated_at') or '', 'region_used': used})
        exported.append(rec)
    filter_new_site_variants_write_csv(out_csv, exported, filter_new_site_variants_PASS_CSV_COLS)
    return {'n': len(exported), 'by_region': dict(Counter((r.get('region_used') for r in exported))), 'by_link_tier': dict(Counter((r.get('link_tier') for r in exported))), 'by_resolve_status': dict(Counter((r.get('resolve_status') for r in exported))), 'by_resolve_method': dict(Counter((r.get('resolve_method') for r in exported))), 'pass_candidates_scanned_verdicts': dict(verdicts)}

def filter_new_site_variants_build_review_html(chains: list[dict[str, Any]], annotations: list[dict[str, Any]], out_html: Path, title: str) -> int:
    records = [build_interactive_audit_html_slim_record(r) for r in chains]
    build_interactive_audit_html_overlay_annotations_on_records(records, annotations)
    html_text = build_interactive_audit_html_build_html(build_interactive_audit_html_TEMPLATE.read_text(encoding='utf-8'), records, annotations, title, save_filename=out_html.name)
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html_text, encoding='utf-8')
    return len(records)

def filter_new_site_variants_main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pos-jsonl', type=Path, default=filter_new_site_variants_S9 / 'variants_genomic_positive.jsonl')
    ap.add_argument('--neg-jsonl', type=Path, default=filter_new_site_variants_S9 / 'variants_genomic_negative.jsonl')
    ap.add_argument('--pos-csv', type=Path, default=filter_new_site_variants_S9 / 'variants_genomic_positive.csv')
    ap.add_argument('--neg-csv', type=Path, default=filter_new_site_variants_S9 / 'variants_genomic_negative.csv')
    ap.add_argument('--out-dir', type=Path, default=filter_new_site_variants_NEW_SITE)
    ap.add_argument('--annotations', type=Path, default=filter_new_site_variants_MANUAL / 'audit_annotations.json')
    args = ap.parse_args()
    print('Loading literature indexes for c.(rs) pairing…')
    articles = {str(r.get('pmid')): r for r in filter_new_site_variants_load_jsonl(filter_new_site_variants_ROOT / 'data/pubmed/ingest/merged/articles.jsonl')}
    sources = [filter_new_site_variants_load_jsonl_by_pmid(filter_new_site_variants_ROOT / 'data/pubmed/pipeline/s6_review/corpus_enriched_novel_clinvar_absent_review_highconf_journal_annotated.jsonl'), filter_new_site_variants_load_jsonl_by_pmid(filter_new_site_variants_ROOT / 'data/pubmed/pipeline/s4_discovery/corpus_enriched_novel.jsonl'), filter_new_site_variants_load_jsonl_by_pmid(filter_new_site_variants_ROOT / 'data/pubmed/pipeline/s2_linking/pairs.jsonl'), filter_new_site_variants_load_jsonl_by_pmid(filter_new_site_variants_ROOT / 'data/pubmed/pipeline/s2_linking/corpus_fulltext.jsonl'), filter_new_site_variants_load_jsonl_by_pmid(filter_new_site_variants_ROOT / 'data/pubmed/pipeline/s8_negatives/negatives_independent/negatives_full_pass.jsonl')]
    pos = filter_new_site_variants_load_jsonl(args.pos_jsonl)
    neg = filter_new_site_variants_load_jsonl(args.neg_jsonl)
    pos_keep, pos_drop, pos_stats = filter_new_site_variants_filter_table(pos, articles=articles, sources=sources)
    neg_keep, neg_drop, neg_stats = filter_new_site_variants_filter_table(neg, articles=articles, sources=sources)
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    filter_new_site_variants_write_jsonl(out / 'variants_genomic_positive_new_site.jsonl', pos_keep)
    filter_new_site_variants_write_jsonl(out / 'variants_genomic_negative_new_site.jsonl', neg_keep)
    filter_new_site_variants_write_jsonl(out / 'dropped_rs_positive.jsonl', pos_drop)
    filter_new_site_variants_write_jsonl(out / 'dropped_rs_negative.jsonl', neg_drop)
    pos_fields = list(csv.DictReader(args.pos_csv.open(encoding='utf-8')).fieldnames or [])
    neg_fields = list(csv.DictReader(args.neg_csv.open(encoding='utf-8')).fieldnames or [])
    extra = ['rs_keep_reason', 'literature_paired_c_hgvs']
    filter_new_site_variants_write_csv(out / 'variants_genomic_positive_new_site.csv', pos_keep, pos_fields + extra)
    filter_new_site_variants_write_csv(out / 'variants_genomic_negative_new_site.csv', neg_keep, neg_fields + extra)
    pos_keys = Counter((filter_new_site_variants_rec_key(r) for r in pos_keep))
    neg_keys = Counter((filter_new_site_variants_rec_key(r) for r in neg_keep))
    pos_chains = filter_new_site_variants_filter_chains(filter_new_site_variants_TRAILS / 'stage_audit_chains_positive.jsonl', pos_keys)
    neg_chains = filter_new_site_variants_filter_chains(filter_new_site_variants_TRAILS / 'stage_audit_chains_negative.jsonl', neg_keys)
    filter_new_site_variants_write_jsonl(out / 'stage_audit_chains_positive.jsonl', pos_chains)
    filter_new_site_variants_write_jsonl(out / 'stage_audit_chains_negative.jsonl', neg_chains)
    filter_new_site_variants_write_jsonl(out / 'stage_audit_chains.jsonl', pos_chains + neg_chains)
    pos_tiers = filter_new_site_variants_split_chains(pos_chains, out, 'stage_audit_chains_positive')
    neg_tiers = filter_new_site_variants_split_chains(neg_chains, out, 'stage_audit_chains_negative')
    annotations: list[dict[str, Any]] = []
    if args.annotations.exists():
        annotations = json.loads(args.annotations.read_text(encoding='utf-8'))
        if not isinstance(annotations, list):
            annotations = []
    n_pos_html = filter_new_site_variants_build_review_html(pos_chains, annotations, filter_new_site_variants_MANUAL / 'review_positive.html', 'Variant clinical Stage0–9 交互审核 · 正样本（新位点）')
    n_neg_html = filter_new_site_variants_build_review_html(neg_chains, annotations, filter_new_site_variants_MANUAL / 'review_negative.html', 'Variant clinical Stage0–9 交互审核 · 负样本（新位点）')
    pos_slim = [build_interactive_audit_html_slim_record(r) for r in pos_chains]
    build_interactive_audit_html_overlay_annotations_on_records(pos_slim, annotations)
    n_overlaid = sum((1 for r in pos_slim if r.get('human_verdict')))
    coding_stats = filter_new_site_variants_export_pass_csv(jsonl_rows=pos_keep, chains=pos_chains, annotations=annotations, out_csv=filter_new_site_variants_MANUAL / 'review_positive_pass_coding.csv', allowed_regions=filter_new_site_variants_CODING)
    noncoding_stats = filter_new_site_variants_export_pass_csv(jsonl_rows=pos_keep, chains=pos_chains, annotations=annotations, out_csv=filter_new_site_variants_MANUAL / 'review_positive_pass_noncoding.csv', allowed_regions=filter_new_site_variants_NONCODING)
    report = {'policy': {'drop': 'variant_text is rsID without literature c.HGVS (rs) pair', 'keep_non_rs': True, 'keep_rs_if': 'c.123A>G (rs…) or rs… (c.…) or NM_/NR_:c. (rs…) in literature blob', 'do_not_use_clinvar_c_hgvs_alone': True}, 'source': {'positive_jsonl': str(args.pos_jsonl), 'negative_jsonl': str(args.neg_jsonl), 'annotations': str(args.annotations), 'annotations_n': len(annotations)}, 'positive': {'input': len(pos), 'kept': len(pos_keep), 'dropped_rs': len(pos_drop), 'reasons': dict(pos_stats), 'chains': len(pos_chains), 'html': n_pos_html, 'html_with_overlaid_verdict': n_overlaid, 'tiers': pos_tiers}, 'negative': {'input': len(neg), 'kept': len(neg_keep), 'dropped_rs': len(neg_drop), 'reasons': dict(neg_stats), 'chains': len(neg_chains), 'html': n_neg_html, 'tiers': neg_tiers}, 'pass_export': {'coding': coding_stats, 'noncoding': noncoding_stats}, 'notes': ['Original Stage9 JSONL/CSV left in place.', 'audit_annotations.json not rewritten; HTML embed keeps the full annotation list.', 'Dropped rs review rows remain in annotations for overlay if a matching record returns.']}
    (out / 'filter_new_site_report.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    (filter_new_site_variants_MANUAL / 'review_positive_pass_coding_stats.json').write_text(json.dumps({'source_html': str(filter_new_site_variants_MANUAL / 'review_positive.html'), 'source_jsonl': str(out / 'variants_genomic_positive_new_site.jsonl'), 'output_csv': str(filter_new_site_variants_MANUAL / 'review_positive_pass_coding.csv'), 'filter': {'human_verdict': 'pass', 'region_used': sorted(filter_new_site_variants_CODING), 'region_source': 'human_genomic_region if set else pipeline genomic_region', 'rs_policy': report['policy']}, 'review_html': {'embedded_records': n_pos_html, 'annotations': len(annotations)}, 'exported': coding_stats}, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    (filter_new_site_variants_MANUAL / 'review_positive_pass_noncoding_stats.json').write_text(json.dumps({'source_html': str(filter_new_site_variants_MANUAL / 'review_positive.html'), 'source_jsonl': str(out / 'variants_genomic_positive_new_site.jsonl'), 'output_csv': str(filter_new_site_variants_MANUAL / 'review_positive_pass_noncoding.csv'), 'filter': {'human_verdict': 'pass', 'region_used': sorted(filter_new_site_variants_NONCODING), 'region_source': 'human_genomic_region if set else pipeline genomic_region', 'rs_policy': report['policy']}, 'review_html': {'embedded_records': n_pos_html, 'annotations': len(annotations)}, 'exported': noncoding_stats}, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, ensure_ascii=False, indent=2))

# === build_stage_audit_chains_jsonl.py ===
"""
按 EXAMPLE_positive_stage_chain_SAMD9_c2423AtoG.md 的结构，
为每条正/负样本生成 Stage0–9 纯文本（中文解释）思维链，输出 JSONL。

每行含：stage0_text…stage9_text、reasoning_chain_text，以及人工填写空栏。
"""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Any
build_stage_audit_chains_jsonl_ROOT = Path(__file__).resolve().parents[1]
build_stage_audit_chains_jsonl_DEFAULT_POS_CSV = build_stage_audit_chains_jsonl_ROOT / 'data/pubmed/pipeline/s9_genomic' / 'variants_genomic_positive.csv'
build_stage_audit_chains_jsonl_DEFAULT_NEG_CSV = build_stage_audit_chains_jsonl_ROOT / 'data/pubmed/pipeline/s9_genomic' / 'variants_genomic_negative.csv'
build_stage_audit_chains_jsonl_DEFAULT_OUT = build_stage_audit_chains_jsonl_ROOT / 'data/pubmed/pipeline/s9_genomic' / 'audit_trails' / 'stage_audit_chains.jsonl'
build_stage_audit_chains_jsonl_PATHS = {'articles': build_stage_audit_chains_jsonl_ROOT / 'data/pubmed/ingest/merged/articles.jsonl', 'pairs': build_stage_audit_chains_jsonl_ROOT / 'data/pubmed/pipeline/s2_linking/pairs.jsonl', 'articles_pass': build_stage_audit_chains_jsonl_ROOT / 'data/pubmed/pipeline/s2_linking/articles_pass.jsonl', 'ft_index': build_stage_audit_chains_jsonl_ROOT / 'data/pubmed/cache/fulltext/linked/articles_fulltext_index.jsonl', 'corpus_fulltext': build_stage_audit_chains_jsonl_ROOT / 'data/pubmed/pipeline/s2_linking/corpus_fulltext.jsonl', 'discovery_novel': build_stage_audit_chains_jsonl_ROOT / 'data/pubmed/pipeline/s4_discovery/corpus_enriched_novel.jsonl', 'clinvar_absent': build_stage_audit_chains_jsonl_ROOT / 'data/pubmed/pipeline/s5_clinvar/enriched' / 'corpus_enriched_novel_clinvar_absent.jsonl', 'review_highconf': build_stage_audit_chains_jsonl_ROOT / 'data/pubmed/pipeline/s6_review' / 'corpus_enriched_novel_clinvar_absent_review_highconf_journal_annotated.jsonl', 'review_strict': build_stage_audit_chains_jsonl_ROOT / 'data/pubmed/pipeline/s6_review' / 'corpus_enriched_novel_clinvar_absent_review_strict_journal_annotated.jsonl', 'pdf_patho': build_stage_audit_chains_jsonl_ROOT / 'data/pubmed/cache/pdfs' / 'corpus_enriched_novel_clinvar_absent_review_strict_pdf_pathogenicity.jsonl', 'pdf_dir': build_stage_audit_chains_jsonl_ROOT / 'data/pubmed/cache/pdfs', 'queries': build_stage_audit_chains_jsonl_ROOT / 'data/pubmed/ref/queries/queries_used.json', 'neg_src': build_stage_audit_chains_jsonl_ROOT / 'data/pubmed/pipeline/s8_negatives/negatives_independent' / 'negatives_full_pass.jsonl'}
build_stage_audit_chains_jsonl_QUERY_HINT_WORDS: dict[str, list[str]] = {'q1_core_gpc': ['genotype', 'phenotype', 'phenotypic spectrum', 'clinical spectrum', 'clinical manifestation', 'missense', 'pathogenic', 'de novo', 'variant', 'mutation'], 'q2_high_recall': ['variant', 'mutation', 'pathogenic', 'phenotype', 'clinical', 'syndrome', 'syndromic', 'neurodevelopmental', 'rare disease'], 'q3_case_reports': ['case', 'case report', 'family', 'proband', 'pedigree', 'novel', 'pathogenic', 'de novo', 'mutation', 'phenotype'], 'q4_mesh': ['genetic', 'mutation', 'phenotype', 'sequencing'], 'q5_clinical_genetics': ['pathogenic', 'mutation', 'diagnosis', 'clinical genetics', 'variant interpretation', 'ACMG'], 'q6_rsid_trait': ['SNP', 'GWAS', 'polymorphism', 'association', 'trait'], 'q7_intergenic_regulatory': ['intergenic', 'enhancer', 'regulatory', 'noncoding']}

def build_stage_audit_chains_jsonl_norm_var(v: str) -> str:
    return re.sub('\\s+', '', (v or '').strip().lower())

def build_stage_audit_chains_jsonl_clip(s: str | None, n: int=700) -> str:
    t = re.sub('\\s+', ' ', s or '').strip()
    if len(t) <= n:
        return t
    return t[:n - 1] + '…'

def build_stage_audit_chains_jsonl_variant_matches(row_variant: str, cand: str) -> bool:
    a = build_stage_audit_chains_jsonl_norm_var(row_variant)
    b = build_stage_audit_chains_jsonl_norm_var(cand)
    if not a or not b:
        return False
    if a == b:
        return True
    if a.endswith(b) or b.endswith(a):
        return True
    a2, b2 = (a.split(':')[-1], b.split(':')[-1])
    return a2 == b2 or a2 in b2 or b2 in a2

def build_stage_audit_chains_jsonl_pick_variant_rows(rows: list[dict[str, Any]], variant: str, key: str='variant') -> list[dict[str, Any]]:
    return [r for r in rows if build_stage_audit_chains_jsonl_variant_matches(variant, str(r.get(key) or ''))]

def build_stage_audit_chains_jsonl_load_jsonl_by_pmid(path: Path) -> dict[str, list[dict[str, Any]]]:
    out: dict[str, list[dict[str, Any]]] = {}
    if not path.exists():
        return out
    with path.open(encoding='utf-8') as f:
        for line in f:
            if not line.strip():
                continue
            r = json.loads(line)
            pmid = str(r.get('pmid') or '').strip()
            if pmid:
                out.setdefault(pmid, []).append(r)
    return out

def build_stage_audit_chains_jsonl_load_articles(path: Path) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    if not path.exists():
        return out
    with path.open(encoding='utf-8') as f:
        for line in f:
            if not line.strip():
                continue
            r = json.loads(line)
            pmid = str(r.get('pmid') or '').strip()
            if pmid:
                out[pmid] = r
    return out

def build_stage_audit_chains_jsonl_load_queries(path: Path) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    if not path.exists():
        return out
    data = json.loads(path.read_text(encoding='utf-8'))
    for q in data.get('queries') or []:
        qid = str(q.get('id') or '')
        if qid:
            out[qid] = {'id': qid, 'label': str(q.get('label') or ''), 'term': str(q.get('term') or '')}
    return out

def build_stage_audit_chains_jsonl_snippet_around(text: str, word: str, radius: int=42) -> str | None:
    m = re.search(re.escape(word), text, flags=re.I)
    if not m:
        return None
    s = max(0, m.start() - radius)
    e = min(len(text), m.end() + radius)
    frag = re.sub('\\s+', ' ', text[s:e]).strip()
    return f'…{frag}…'

def build_stage_audit_chains_jsonl_query_hit_rationale(qid: str, title: str, abstract: str, qmeta: dict[str, str]) -> str:
    blob = f'{title}\n{abstract}'
    words = build_stage_audit_chains_jsonl_QUERY_HINT_WORDS.get(qid, ['variant', 'mutation', 'phenotype', 'clinical'])
    hits: list[str] = []
    for w in words:
        sn = build_stage_audit_chains_jsonl_snippet_around(blob, w)
        if sn:
            hits.append(f'「{sn}」（对齐关键词 {w}）')
        if len(hits) >= 4:
            break
    label = qmeta.get('label') or qid
    if hits:
        return f"【{qid}｜{label}】本篇进入该检索集。题录/摘要中可对齐的原文依据：{'；'.join(hits)}。解释：上述措辞与该 query 的 variant/mutation × clinical/phenotype/case 逻辑一致，且满足 humans + 2026 出版窗。"
    return f'【{qid}｜{label}】本篇标记在 source_query_ids 中；题录/摘要未稳定抽到预置关键词片段，但仍属于该查询合并结果。'

def build_stage_audit_chains_jsonl_pdf_sha(pmid: str) -> tuple[str, str]:
    pdf = build_stage_audit_chains_jsonl_PATHS['pdf_dir'] / f'PMID{pmid}.pdf'
    if not pdf.exists():
        return ('', '')
    return (str(pdf), hashlib.sha256(pdf.read_bytes()).hexdigest())

def build_stage_audit_chains_jsonl_method_zh(method: str) -> str:
    m = {'same_sentence_cooccurrence': '同句共现', 'same_abstract_cooccurrence': '同摘要共现', 'fulltext_same_sentence_cooccurrence': '全文同句共现', 'fulltext_same_paragraph_cooccurrence': '全文同段共现', 'fulltext_known_variant_sentence': '全文已知变异句', 'pubtator_variant_disease': 'PubTator 变异–疾病', 'pubtator_gene_disease_relation': 'PubTator 基因–疾病', 'ensembl_mane_c_map': 'Ensembl MANE 转录本 c.HGVS 抬升至基因组', 'ensembl_mane_p_map': 'Ensembl MANE 转录本 p.HGVS 抬升至基因组', 'gene_c_clinvar': 'gene + c.HGVS → ClinVar', 'gene_p_clinvar': 'gene + p.HGVS → ClinVar', 'nm_c_clinvar': 'NM_ + c.HGVS → ClinVar', 'rs_dbsnp': 'rsID → dbSNP RefSNP API', 'rs_clinvar': 'rsID → ClinVar', 'pdf_literal_genomic': 'PDF 字面基因组坐标', 'c_hgvs_alleles_only': '仅从 c.HGVS 取 ref/alt', 'unresolved': '未能解析'}
    return m.get(method, method or '未知')

def build_stage_audit_chains_jsonl_tier_zh(tier: str) -> str:
    return {'A': 'A（高置信紧密链接）', 'B': 'B（可用但较弱）', 'C': 'C（宽松链接）'}.get(tier, tier or '未知')

def build_stage_audit_chains_jsonl_stage0_text() -> str:
    return '\n'.join(['【Stage0 — preflight（环境与参考，跑次级）】', '结论：参考跑次 data/pubmed/runs/gate_runs/v2_full_20260808/ ，配置 gate_config_version=2026-08-08.v2-dbsnp 已 APPROVE。', '检查项与路径：', '- 查询定义：data/pubmed/ref/queries/queries_used.json', '- 期刊表：configs/journal_priority.json', '- ClinVar 快照：data/pubmed/ref/clinvar/variant_summary.txt.gz', '- Ensembl GFF3：data/pubmed/ref/ensembl/Homo_sapiens.GRCh38.115.gff3', '- GRCh38 FASTA/FAI：data/pubmed/ref/ensembl/Homo_sapiens.GRCh38.dna.toplevel.fa(.fai)', '- 输出隔离：data/pubmed/runs/gate_runs/v2_full_20260808/', '说明：本条无逐文献动作；预检通过后整库进入 Stage1。'])

def build_stage_audit_chains_jsonl_stage1_text(art: dict[str, Any] | None, queries: dict[str, dict[str, str]]) -> str:
    lines = ['【Stage1 — pubmed_retrieval（为何被哪些 query 命中）】']
    if not art:
        lines.append('未在 merged/articles.jsonl 找到该 PMID 的题录/摘要，无法回溯检索命中依据。')
        return '\n'.join(lines)
    title = str(art.get('title') or '')
    abstract = str(art.get('abstract') or '')
    qids = [str(x) for x in art.get('source_query_ids') or []]
    pub = art.get('pub_date') or {}
    year = pub.get('year') if isinstance(pub, dict) else pub
    lines += ['文献元数据：', f'- Title：{title}', f"- Journal：{art.get('journal') or ''}（{art.get('journal_iso') or ''}），ISSN {art.get('issn') or ''}，{year or ''}", f"- PMCID：{art.get('pmcid') or '无'}", f'- 命中 query：source_query_ids = {json.dumps(qids, ensure_ascii=False)}', '', '命中依据（题录/摘要原文关键词 ↔ query 子句）：']
    for qid in qids:
        lines.append('- ' + build_stage_audit_chains_jsonl_query_hit_rationale(qid, title, abstract, queries.get(qid) or {}))
    all_q = set(queries)
    missed = sorted(all_q - set(qids))
    if missed:
        lines.append(f"未命中例：{', '.join(missed[:8])}{(' 等' if len(missed) > 8 else '')} 不在 source_query_ids 中（本篇未进入这些检索集或合并后未标记）。")
    lines += ['', '摘要原文（Stage1 入库文本）：', build_stage_audit_chains_jsonl_clip(abstract, 1200) or '（无摘要）']
    return '\n'.join(lines)

def build_stage_audit_chains_jsonl_stage2_text(pmid: str, variant: str, pairs: list[dict[str, Any]], art_pass: dict[str, Any] | None) -> str:
    lines = ['【Stage2 — abstract_linking（摘要变异–临床链接 → tier）】']
    rows = build_stage_audit_chains_jsonl_pick_variant_rows(pairs, variant)
    if not rows:
        lines.append('pairs.jsonl 中未找到与本 variant 匹配的摘要链接行。')
        return '\n'.join(lines)
    rows_sorted = sorted(rows, key=lambda r: {'A': 0, 'B': 1, 'C': 2}.get(str(r.get('tier')), 9))
    best = rows_sorted[0]
    same = [r for r in rows_sorted if r.get('tier') == best.get('tier') and r.get('method') == best.get('method')]
    clinicals = []
    for r in same:
        c = str(r.get('clinical') or '').strip()
        if c and c not in clinicals:
            clinicals.append(c)
    lines += ['最终保留链接（进入 review 的主链）：', f"- variant：{best.get('variant')}", f"- gene：{best.get('gene')}", f"- clinical：{('；'.join(clinicals) if clinicals else best.get('clinical'))}", f"- tier：{best.get('tier')}（{build_stage_audit_chains_jsonl_tier_zh(str(best.get('tier') or ''))}）", f"- method：{best.get('method')}（{build_stage_audit_chains_jsonl_method_zh(str(best.get('method') or ''))}）", f"- evidence_source / text_scope：{best.get('evidence_source')} / {best.get('text_scope')}", f"- clinical_source：{best.get('clinical_source')}", '', '判定逻辑：' + ('同一句子内同时出现具体核苷酸/蛋白变异与临床实体（词典 lexicon_entity 等），满足 concrete variant + dual presence → 对应 tier。' if str(best.get('method') or '').startswith('same_') or 'cooccurrence' in str(best.get('method') or '') else f"由方法 {best.get('method')}（{build_stage_audit_chains_jsonl_method_zh(str(best.get('method') or ''))}）判定链接强度为 tier {best.get('tier')}。"), '', '贴出原文（链接证据句）：']
    ev_seen: set[str] = set()
    n = 0
    for r in rows_sorted:
        ev = str(r.get('evidence_sentence') or '').strip()
        if not ev or ev in ev_seen:
            continue
        ev_seen.add(ev)
        n += 1
        lines.append(f"{n}. （tier={r.get('tier')} / {r.get('method')}）{build_stage_audit_chains_jsonl_clip(ev, 500)}")
        if n >= 4:
            break
    if not ev_seen:
        lines.append('（本主链无 evidence_sentence，可能为 PubTator 等无同句证据的链接。）')
    weak = [r for r in rows_sorted if str(r.get('tier')) in {'B', 'C'}]
    if weak:
        lines.append('')
        lines.append('同 PMID 的较弱旁路（未作为本条主 review 键时供对照）：')
        for r in weak[:5]:
            lines.append(f"- tier {r.get('tier')} / {r.get('method')}：clinical={build_stage_audit_chains_jsonl_clip(str(r.get('clinical') or ''), 80)}；evidence_sentence={build_stage_audit_chains_jsonl_clip(str(r.get('evidence_sentence') or '') or '无', 160)}")
    if art_pass and art_pass.get('best_tier'):
        lines.append(f"文章级 best_tier={art_pass.get('best_tier')}。")
    return '\n'.join(lines)

def build_stage_audit_chains_jsonl_stage3_text(pmid: str, variant: str, ft_index: list[dict[str, Any]], ft_rows: list[dict[str, Any]]) -> str:
    lines = ['【Stage3 — fulltext_enrichment（全文富集）】']
    idx = ft_index[0] if ft_index else None
    if idx:
        lines.append(f"是否有全文：是 — PMCID {idx.get('pmcid')}，来源 {idx.get('source')} （articles_fulltext_index：n_links={idx.get('n_links')}）。")
        lines.append('摘要回退：否（有全文）；摘要轨证据仍保留。')
    else:
        lines.append('是否有全文：否/未建索引 — 可能走摘要回退（allow_abstract_fallback）；本条以摘要证据为主。')
        lines.append('摘要回退：是（或未知）。')
    matched = build_stage_audit_chains_jsonl_pick_variant_rows(ft_rows, variant)
    if not matched:
        lines.append('新增全文证据：无 — corpus_fulltext 中无与本 variant 匹配的同句/同段记录。')
        return '\n'.join(lines)
    lines.append(f'新增全文证据：有 — 共 {len(matched)} 条全文链接行；同句/同段复现变异与表型。')
    lines.append('')
    lines.append('贴出全文新增/强化原文（节选）：')
    matched_sorted = sorted(matched, key=lambda r: (0 if r.get('variant_discovered_in_fulltext') else 1, 0 if str(r.get('tier')) == 'A' else 1, str(r.get('section') or '')))
    for r in matched_sorted[:5]:
        disc = '是' if r.get('variant_discovered_in_fulltext') else '否'
        lines.append(f"- 章节 {r.get('section')}｜tier={r.get('tier')}｜method={r.get('method')}（{build_stage_audit_chains_jsonl_method_zh(str(r.get('method') or ''))}）｜ft_source={r.get('ft_source')}｜variant_discovered_in_fulltext={disc}")
        lines.append(f"  原文：{build_stage_audit_chains_jsonl_clip(str(r.get('evidence_sentence') or ''), 650)}")
    return '\n'.join(lines)

def build_stage_audit_chains_jsonl_stage4_text(pmid: str, variant: str, disc_rows: list[dict[str, Any]]) -> str:
    lines = ['【Stage4 — discovery_tagging（发现性文本标注）】']
    rows = build_stage_audit_chains_jsonl_pick_variant_rows(disc_rows, variant)
    if not rows:
        lines.append('未找到与本 variant 匹配的 discovery 标注行。')
        return '\n'.join(lines)
    rows_sorted = sorted(rows, key=lambda r: 0 if str(r.get('discovery_confidence')) == 'high' else 1)
    best = rows_sorted[0]
    cues = best.get('discovery_cues') or {}
    lines.append(f"本条采用标签：discovery_label={best.get('discovery_label')}，置信度 {best.get('discovery_confidence')}；keep_as_new_discovery={best.get('keep_as_new_discovery')}。")
    lines.append('')
    lines.append('cue 组命中与选择依据：')
    for group in ('novel', 'previously_reported', 'citation', 'phenotype_expansion'):
        hit = cues.get(group) or cues.get('citation_of_prior' if group == 'citation' else group) or []
        label = {'novel': 'novel', 'previously_reported': 'previously_reported', 'citation': 'citation_of_prior', 'phenotype_expansion': 'phenotype_expansion'}[group]
        if hit:
            lines.append(f'- 【{label}】命中 cue={json.dumps(hit, ensure_ascii=False)}；依据见下方原文窗（含 novel/de novo/we identified 等措辞）。')
        else:
            lines.append(f'- 【{label}】（空）未命中该类主导措辞线索。')
    lines.append('')
    lines.append('贴出打标签所依据的原文窗：')
    lines.append(build_stage_audit_chains_jsonl_clip(str(best.get('text') or ''), 900))
    lines.append('')
    lines.append('政策提醒：novel 标签 ≠ 全球文献首报；仅表示文本措辞信号。')
    return '\n'.join(lines)

def build_stage_audit_chains_jsonl_stage5_text(pmid: str, variant: str, cv_rows: list[dict[str, Any]], sample_class: str) -> str:
    lines = ['【Stage5 — clinvar_validation（ClinVar 对照）】']
    rows = build_stage_audit_chains_jsonl_pick_variant_rows(cv_rows, variant)
    if not rows:
        if sample_class == 'negative':
            lines.append('本条为独立负样本，通常不走正样本 ClinVar-absent 主漏斗；无 clinvar_absent JSONL 匹配行属预期（或键未回溯到）。')
        else:
            lines.append('未找到 clinvar_absent 匹配行（可能键不一致或非 absent 轨）。')
        return '\n'.join(lines)
    r = rows[0]
    cq = r.get('clinvar_query') or {}
    lines += [f'- 快照：data/pubmed/ref/clinvar/variant_summary.txt.gz', f"- 检查时间：{r.get('clinvar_checked_at')}", f"- 精确查找键：lookup_keys = {json.dumps(cq.get('lookup_keys') or [], ensure_ascii=False)}", f"- gene_source：{cq.get('gene_source')}（genes={cq.get('genes')}）", f"- c_hgvs / p_hgvs / rs / nm：{cq.get('c_hgvs')} / {cq.get('p_hgvs')} / {cq.get('rs')} / {cq.get('nm')}", f"- clinvar_match_status：{r.get('clinvar_match_status')}", f"- clinvar_hit：{json.dumps(r.get('clinvar_hit') or {}, ensure_ascii=False)}", f"- novelty_vs_clinvar：{r.get('novelty_vs_clinvar')}（confidence={r.get('novelty_vs_clinvar_confidence')}）", '', '解读：在该版本 ClinVar 中按上表精确键检索；' + ('无命中 → 进入 ClinVar-absent 候选轨。' if r.get('novelty_vs_clinvar') == 'variant_absent' else f"结果为 {r.get('novelty_vs_clinvar')}。"), '政策：ClinVar absent ≠ 证明全球文献首报。']
    return '\n'.join(lines)

def build_stage_audit_chains_jsonl_stage6_text(pmid: str, variant: str, review_rows: list[dict[str, Any]], sample_class: str, neg_hits: list[dict[str, Any]]) -> str:
    lines = ['【Stage6 — review_quality（可复阅高置信子集）】']
    if sample_class == 'negative':
        lines.append('本条为独立负样本（N-ft），不经由正样本 review_strict/highconf 漏斗筛选；导出 review_tier 一般为 independent_neg。')
        if neg_hits:
            nr = neg_hits[0]
            lines.append(f"- 负样本角色/标签：role={nr.get('negative_role')} label={nr.get('negative_label')}")
            lines.append(f"- 贴出负样本证据文本：{build_stage_audit_chains_jsonl_clip(str(nr.get('negative_evidence') or ''), 900)}")
        return '\n'.join(lines)
    rows = build_stage_audit_chains_jsonl_pick_variant_rows(review_rows, variant)
    if not rows:
        lines.append('未找到 review_highconf/strict 注释行。')
        return '\n'.join(lines)
    rows_sorted = sorted(rows, key=lambda r: (0 if str(r.get('review_tier')) == 'strict' else 1, 0 if str(r.get('discovery_confidence')) == 'high' else 1))
    r = rows_sorted[0]
    lines += ['进入 review 子集的门控条件（本条取值）：', f"- link tier：{r.get('tier')}", f"- discovery_label / confidence：{r.get('discovery_label')} / {r.get('discovery_confidence')}", f"- method（可信集合）：{r.get('method')}（{build_stage_audit_chains_jsonl_method_zh(str(r.get('method') or ''))}）", f"- clinical_source：{r.get('clinical_source')}", f"- gene_source：{r.get('gene_source')}", f"- gene confirmed：{r.get('gene')}", f"- locus_key：{r.get('locus_key')}", f"- review_tier：{r.get('review_tier')}", '', '贴出具体证据文本（写入 review JSONL 的 evidence_text）：', build_stage_audit_chains_jsonl_clip(str(r.get('evidence_text') or ''), 900) or '（evidence_text 为空，可回看 Stage2/3 证据句）', f"临床绑定：{r.get('clinical')}"]
    return '\n'.join(lines)

def build_stage_audit_chains_jsonl_stage7_text(pmid: str, variant: str, review_rows: list[dict[str, Any]], export_row: dict[str, str], art: dict[str, Any] | None) -> str:
    lines = ['【Stage7 — journal_priority（期刊优先级）】']
    rows = build_stage_audit_chains_jsonl_pick_variant_rows(review_rows, variant) if review_rows else []
    r = rows[0] if rows else {}
    journal = r.get('journal') or (art or {}).get('journal') or ''
    iso = r.get('journal_iso') or (art or {}).get('journal_iso') or ''
    issn = r.get('journal_issn') or (art or {}).get('issn') or ''
    tier = r.get('journal_tier') or export_row.get('screen_grade') or ''
    lines += [f'- 来源期刊（原文）：{journal}', f'- ISO / ISSN：{iso} / {issn}', f"- 规范化名：{r.get('journal_canonical_name') or ''}", f"- journal_tier → screen_grade：{tier}（导出 screen_grade={export_row.get('screen_grade')}）", f"- 匹配方式：{r.get('journal_match_method') or '导出动态匹配/缺省'}", f"- 优先级分 / 理由：{(r.get('journal_priority_score') if r else '')} / {r.get('journal_priority_reason') or ''}", f"- 期刊表版本：{r.get('journal_list_version') or 'configs/journal_priority.json'}", '', '说明：C/WATCH 不删除；期刊等级只用于复阅排序，不单独证明变异真实性或 novelty。']
    return '\n'.join(lines)

def build_stage_audit_chains_jsonl_stage8_text(pmid: str, variant: str, patho_rows: list[dict[str, Any]], export_row: dict[str, str]) -> str:
    lines = ['【Stage8 — pdf_evidence（PDF 致病性证据）】']
    path, sha = build_stage_audit_chains_jsonl_pdf_sha(pmid)
    if path:
        lines.append(f'- PDF 路径：{path}')
        lines.append(f'- SHA-256：{sha}')
    else:
        lines.append('- PDF 路径：pdfs_review_strict 下未找到该 PMID 的 PDF（可能未下载或非 strict 覆盖）。')
    rows = build_stage_audit_chains_jsonl_pick_variant_rows(patho_rows, variant)
    if not rows:
        lines.append(f"- pdf_pathogenicity 标注：无与本 variant 匹配的记录；导出 label={export_row.get('label')} （可为 positive_candidate / 负样本标签）。")
        return '\n'.join(lines)
    r = rows[0]
    ev = str(r.get('pdf_pathogenicity_evidence') or '')
    lines += [f"- 锚定变异：{r.get('variant')}（要求证据句提及该 HGVS）", f"- pdf_pathogenicity：{r.get('pdf_pathogenicity')}（confidence={r.get('pdf_pathogenicity_confidence')}）", '- 与 in-silico：优先采用病例/ACMG 叙述句；若证据仅为预测软件需人工降权。', '', '贴出 PDF 抽取原文：', build_stage_audit_chains_jsonl_clip(ev, 1000) or '（无证据句）', f"→ 导出 CSV label={export_row.get('label')}。"]
    return '\n'.join(lines)

def build_stage_audit_chains_jsonl_stage9_text(export_row: dict[str, str]) -> str:
    method = str(export_row.get('resolve_method') or '')
    strand = str(export_row.get('c_hgvs_strand') or '')
    c_hgvs = str(export_row.get('c_hgvs') or '')
    lines = ['【Stage9 — genomic_export_qc（位点 → chr/pos/ref/alt）】', '', '文献变异描述：', f"- variant_text / nt_pro_only：{export_row.get('variant_text')} / {export_row.get('nt_pro_only')}", f"- 导出 c_hgvs / p_hgvs：{export_row.get('c_hgvs') or '∅'} / {export_row.get('p_hgvs') or '∅'}", '', '解析路径：', f"1. 解析方法 resolve_method={method}（{build_stage_audit_chains_jsonl_method_zh(method)}）；allele_query_source={export_row.get('allele_query_source') or '无'}。", f"2. 转录本：id={export_row.get('transcript_id') or '∅'} version={export_row.get('transcript_version') or '∅'} match={export_row.get('transcript_match') or '∅'}；c_hgvs_strand={strand or '∅'}。"]
    if strand == '-' and re.search('[ACGT]>[ACGT]', c_hgvs, re.I):
        lines.append('3. 解释：转录本为负链时，c.HGVS 的碱基相对于编码链；映射到基因组正链时 ref/alt 取互补（例如 c.A>G → 基因组 T>C）。')
    elif strand == '+' and re.search('[ACGT]>[ACGT]', c_hgvs, re.I):
        lines.append('3. 解释：转录本为正链时，c.HGVS 碱基方向与基因组正链一致。')
    else:
        lines.append('3. 解释：按该方法从文献/数据库等位映射到 GRCh38 坐标；详见 resolve_note。')
    lines += [f"4. FASTA/状态：resolve_status={export_row.get('resolve_status')} confidence={export_row.get('resolve_confidence')}；assembly={export_row.get('assembly')} source={export_row.get('assembly_source')}。", f"5. resolve_note：{export_row.get('resolve_note') or '（无）'}", '', '最终基因组等位（导出）：', f"- chromosome_ucsc：{export_row.get('chromosome_ucsc') or '∅'}", f"- position：{export_row.get('position') or '∅'}", f"- ref：{export_row.get('ref') or '∅'}", f"- alt：{export_row.get('alt') or '∅'}", f"- gene / locus_key_src：{export_row.get('gene')} / {export_row.get('locus_key_src')}", '', 'GFF3 区域分类（机器标注，人工需核对）：', f"- genomic_region：{export_row.get('genomic_region') or '∅（无坐标或未分类）'}", f"- genomic_region_genes：{export_row.get('genomic_region_genes') or '∅'}", f"- genomic_region_feature：{export_row.get('genomic_region_feature') or '∅'}", f"- genomic_region_transcripts：{build_stage_audit_chains_jsonl_clip(str(export_row.get('genomic_region_transcripts') or ''), 240) or '∅'}", f"- genomic_region_note：{export_row.get('genomic_region_note') or '（无）'}"]
    return '\n'.join(lines)

def build_stage_audit_chains_jsonl_human_checklist(export_row: dict[str, str]) -> str:
    gene = export_row.get('gene') or ''
    vt = export_row.get('variant_text') or ''
    chrom = export_row.get('chromosome_ucsc') or '∅'
    pos = export_row.get('position') or '∅'
    ref = export_row.get('ref') or '∅'
    alt = export_row.get('alt') or '∅'
    return '\n'.join(['【人工核对栏】', 'Stage1–2 命中与 tier 是否被原文支持？ [ ] 是 [ ] 存疑 [ ] 否', 'Stage3 全文证据是否锚定同一变异（非串扰）？ [ ] 是 [ ] 存疑 [ ] 否', 'Stage4 discovery 标签是否过度/不足？ [ ] 可接受 [ ] 应调整', f'Stage5 ClinVar 键（gene|allele）是否合理？基因={gene} 变异={vt} [ ] 是 [ ] 否', 'Stage8 致病性标签是否被 PDF/证据句支持（非纯预测软件）？ [ ] 是 [ ] 存疑 [ ] 否', f'Stage9 {chrom}:{pos} {ref}>{alt} 与文献 HGVS 是否自洽？ [ ] 是 [ ] 需外部工具复核', f"Stage9 genomic_region={export_row.get('genomic_region') or '∅'} （重叠基因={export_row.get('genomic_region_genes') or '∅'}）是否正确？ [ ] 是 [ ] 应改为 coding/utr/intron/ncrna_exon/ncrna_gene/pseudogene/intergenic", '总评：[ ] pass [ ] doubt [ ] fail；错在 Stage：____'])

def build_stage_audit_chains_jsonl_build_record(export_row: dict[str, str], *, articles: dict[str, dict[str, Any]], pairs_by_pmid: dict[str, list[dict[str, Any]]], pass_by_pmid: dict[str, list[dict[str, Any]]], ft_index_by_pmid: dict[str, list[dict[str, Any]]], ft_by_pmid: dict[str, list[dict[str, Any]]], disc_by_pmid: dict[str, list[dict[str, Any]]], cv_by_pmid: dict[str, list[dict[str, Any]]], review_by_pmid: dict[str, list[dict[str, Any]]], patho_by_pmid: dict[str, list[dict[str, Any]]], neg_by_pmid: dict[str, list[dict[str, Any]]], queries: dict[str, dict[str, str]], audit_id: str) -> dict[str, Any]:
    pmid = str(export_row.get('pmid') or '')
    variant = str(export_row.get('variant_text') or '')
    sample_class = str(export_row.get('sample_class') or '')
    art = articles.get(pmid)
    art_pass = (pass_by_pmid.get(pmid) or [None])[0]
    neg_hits = [r for r in neg_by_pmid.get(pmid, []) if build_stage_audit_chains_jsonl_variant_matches(variant, str(r.get('negative_variant') or ''))]
    s0 = build_stage_audit_chains_jsonl_stage0_text()
    s1 = build_stage_audit_chains_jsonl_stage1_text(art, queries)
    s2 = build_stage_audit_chains_jsonl_stage2_text(pmid, variant, pairs_by_pmid.get(pmid, []), art_pass)
    s3 = build_stage_audit_chains_jsonl_stage3_text(pmid, variant, ft_index_by_pmid.get(pmid, []), ft_by_pmid.get(pmid, []))
    s4 = build_stage_audit_chains_jsonl_stage4_text(pmid, variant, disc_by_pmid.get(pmid, []))
    s5 = build_stage_audit_chains_jsonl_stage5_text(pmid, variant, cv_by_pmid.get(pmid, []), sample_class)
    s6 = build_stage_audit_chains_jsonl_stage6_text(pmid, variant, review_by_pmid.get(pmid, []), sample_class, neg_hits)
    s7 = build_stage_audit_chains_jsonl_stage7_text(pmid, variant, review_by_pmid.get(pmid, []), export_row, art)
    s8 = build_stage_audit_chains_jsonl_stage8_text(pmid, variant, patho_by_pmid.get(pmid, []), export_row)
    s9 = build_stage_audit_chains_jsonl_stage9_text(export_row)
    checklist = build_stage_audit_chains_jsonl_human_checklist(export_row)
    header = '\n'.join([f'正负样本 Stage0–9 证据思维链（纯文本）', f"样本 ID：{sample_class} | {export_row.get('gene')} | {variant} | PMID {pmid}", f"导出结果：{export_row.get('chromosome_ucsc') or '∅'}:{export_row.get('position') or '∅'} {export_row.get('ref') or '∅'}>{export_row.get('alt') or '∅'}（{export_row.get('assembly') or 'GRCh38'}）| label={export_row.get('label')} | link_tier={export_row.get('link_tier')} | review_tier={export_row.get('review_tier')} | screen_grade={export_row.get('screen_grade')} | genomic_region={export_row.get('genomic_region') or '∅'}", f"PubMed：{export_row.get('pubmed_url')}", '说明：Stage0 为跑次级环境检查（非逐文献）；Stage1–9 均为本条回溯。解释为中文，证据原文保持文献语言。'])
    chain = '\n\n'.join([header, s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, checklist])
    return {'audit_id': audit_id, 'sample_class': sample_class, 'pmid': pmid, 'gene': export_row.get('gene') or '', 'variant_text': variant, 'locus_key_src': export_row.get('locus_key_src') or '', 'label': export_row.get('label') or '', 'link_tier': export_row.get('link_tier') or '', 'review_tier': export_row.get('review_tier') or '', 'screen_grade': export_row.get('screen_grade') or '', 'chromosome_ucsc': export_row.get('chromosome_ucsc') or '', 'position': export_row.get('position') or '', 'ref': export_row.get('ref') or '', 'alt': export_row.get('alt') or '', 'c_hgvs': export_row.get('c_hgvs') or '', 'p_hgvs': export_row.get('p_hgvs') or '', 'resolve_method': export_row.get('resolve_method') or '', 'resolve_status': export_row.get('resolve_status') or '', 'pubmed_url': export_row.get('pubmed_url') or '', 'genomic_region': export_row.get('genomic_region') or '', 'genomic_region_genes': export_row.get('genomic_region_genes') or '', 'genomic_region_feature': export_row.get('genomic_region_feature') or '', 'genomic_region_transcripts': export_row.get('genomic_region_transcripts') or '', 'genomic_region_note': export_row.get('genomic_region_note') or '', 'human_genomic_region': '', 'stage0_text': s0, 'stage1_text': s1, 'stage2_text': s2, 'stage3_text': s3, 'stage4_text': s4, 'stage5_text': s5, 'stage6_text': s6, 'stage7_text': s7, 'stage8_text': s8, 'stage9_text': s9, 'human_checklist_text': checklist, 'reasoning_chain_text': chain, 'human_verdict': '', 'human_fail_stage': '', 'human_comment': ''}

def build_stage_audit_chains_jsonl_load_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f))

def build_stage_audit_chains_jsonl_main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pos-csv', type=Path, default=build_stage_audit_chains_jsonl_DEFAULT_POS_CSV)
    ap.add_argument('--neg-csv', type=Path, default=build_stage_audit_chains_jsonl_DEFAULT_NEG_CSV)
    ap.add_argument('--out-jsonl', type=Path, default=build_stage_audit_chains_jsonl_DEFAULT_OUT)
    ap.add_argument('--class', dest='only_class', choices=['positive', 'negative', 'all'], default='all')
    ap.add_argument('--limit', type=int, default=0)
    args = ap.parse_args()
    print('Loading indexes…')
    articles = build_stage_audit_chains_jsonl_load_articles(build_stage_audit_chains_jsonl_PATHS['articles'])
    pairs_by_pmid = build_stage_audit_chains_jsonl_load_jsonl_by_pmid(build_stage_audit_chains_jsonl_PATHS['pairs'])
    pass_by_pmid = build_stage_audit_chains_jsonl_load_jsonl_by_pmid(build_stage_audit_chains_jsonl_PATHS['articles_pass'])
    ft_index_by_pmid = build_stage_audit_chains_jsonl_load_jsonl_by_pmid(build_stage_audit_chains_jsonl_PATHS['ft_index'])
    ft_by_pmid = build_stage_audit_chains_jsonl_load_jsonl_by_pmid(build_stage_audit_chains_jsonl_PATHS['corpus_fulltext'])
    disc_by_pmid = build_stage_audit_chains_jsonl_load_jsonl_by_pmid(build_stage_audit_chains_jsonl_PATHS['discovery_novel'])
    cv_by_pmid = build_stage_audit_chains_jsonl_load_jsonl_by_pmid(build_stage_audit_chains_jsonl_PATHS['clinvar_absent'])
    review_by_pmid = build_stage_audit_chains_jsonl_load_jsonl_by_pmid(build_stage_audit_chains_jsonl_PATHS['review_highconf'])
    for pmid, rows in build_stage_audit_chains_jsonl_load_jsonl_by_pmid(build_stage_audit_chains_jsonl_PATHS['review_strict']).items():
        review_by_pmid.setdefault(pmid, []).extend(rows)
    patho_by_pmid = build_stage_audit_chains_jsonl_load_jsonl_by_pmid(build_stage_audit_chains_jsonl_PATHS['pdf_patho'])
    neg_by_pmid = build_stage_audit_chains_jsonl_load_jsonl_by_pmid(build_stage_audit_chains_jsonl_PATHS['neg_src'])
    queries = build_stage_audit_chains_jsonl_load_queries(build_stage_audit_chains_jsonl_PATHS['queries'])
    rows: list[dict[str, str]] = []
    if args.only_class in {'positive', 'all'}:
        rows.extend(build_stage_audit_chains_jsonl_load_csv_rows(args.pos_csv))
    if args.only_class in {'negative', 'all'}:
        rows.extend(build_stage_audit_chains_jsonl_load_csv_rows(args.neg_csv))
    if args.limit and args.limit > 0:
        rows = rows[:args.limit]
    args.out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with args.out_jsonl.open('w', encoding='utf-8') as out:
        for i, row in enumerate(rows, 1):
            sample_class = row.get('sample_class') or 'positive'
            audit_id = f"{sample_class}_{row.get('pmid')}_{i:04d}"
            rec = build_stage_audit_chains_jsonl_build_record(row, articles=articles, pairs_by_pmid=pairs_by_pmid, pass_by_pmid=pass_by_pmid, ft_index_by_pmid=ft_index_by_pmid, ft_by_pmid=ft_by_pmid, disc_by_pmid=disc_by_pmid, cv_by_pmid=cv_by_pmid, review_by_pmid=review_by_pmid, patho_by_pmid=patho_by_pmid, neg_by_pmid=neg_by_pmid, queries=queries, audit_id=audit_id)
            out.write(json.dumps(rec, ensure_ascii=False) + '\n')
            n += 1
            if n % 200 == 0:
                print(f'  wrote {n}…')
    if not args.limit:
        from collections import defaultdict
        buckets: dict[tuple[str, str], list[str]] = defaultdict(list)
        with args.out_jsonl.open(encoding='utf-8') as src:
            for line in src:
                if not line.strip():
                    continue
                rec = json.loads(line)
                cls = str(rec.get('sample_class') or 'unknown').strip().lower()
                tier = str(rec.get('link_tier') or '').strip().upper()
                if tier not in {'A', 'B', 'C'}:
                    tier = tier if tier else 'NONE'
                buckets[cls, tier].append(line if line.endswith('\n') else line + '\n')
        for (cls, tier), lines in sorted(buckets.items()):
            out_tier = args.out_jsonl.with_name(f'stage_audit_chains_{cls}_tier_{tier}.jsonl')
            with out_tier.open('w', encoding='utf-8') as w:
                w.writelines(lines)
            print(f'split: {cls} tier={tier} n={len(lines)} -> {out_tier}')
        if args.only_class in {'positive', 'all'}:
            pos_out = args.out_jsonl.with_name('stage_audit_chains_positive.jsonl')
            with pos_out.open('w', encoding='utf-8') as po:
                npos = 0
                for (cls, _tier), lines in sorted(buckets.items()):
                    if cls == 'positive':
                        po.writelines(lines)
                        npos += len(lines)
            print(f'split: positive={npos} -> {pos_out}')
        if args.only_class in {'negative', 'all'}:
            neg_out = args.out_jsonl.with_name('stage_audit_chains_negative.jsonl')
            with neg_out.open('w', encoding='utf-8') as no:
                nneg = 0
                for (cls, _tier), lines in sorted(buckets.items()):
                    if cls == 'negative':
                        no.writelines(lines)
                        nneg += len(lines)
            print(f'split: negative={nneg} -> {neg_out}')
    print(f'wrote {n} records -> {args.out_jsonl}')

# === build_interactive_audit_html.py ===
"""Build a browser-interactive Stage0–9 audit reviewer (HTML).

Usage:
  # Standalone shell (load JSONL via drag/drop in browser)
  python scripts/build_interactive_audit_html.py

  # Embed one JSONL pack (opens offline with data baked in)
  python scripts/build_interactive_audit_html.py \\
    --embed data/.../audit_trails/stage_audit_chains_positive_tier_A.jsonl \\
    --out   data/.../audit_trails/review_positive_tier_A.html
"""
import argparse
import html
import json
from pathlib import Path
build_interactive_audit_html_ROOT = Path(__file__).resolve().parents[1]
build_interactive_audit_html_DEFAULT_TRAILS = build_interactive_audit_html_ROOT / 'data/pubmed/pipeline/s9_genomic/audit_trails'
build_interactive_audit_html_DEFAULT_OUT = build_interactive_audit_html_DEFAULT_TRAILS / 'audit_interactive_review.html'
build_interactive_audit_html_TEMPLATE = Path(__file__).resolve().parent / 'templates' / 'audit_interactive_review.html'

def build_interactive_audit_html_rec_overlay_key(row: dict) -> tuple[str, str, str, str]:
    return (str(row.get('sample_class') or '').strip().lower(), str(row.get('pmid') or '').strip(), str(row.get('gene') or '').strip().upper(), str(row.get('variant_text') or '').strip())

def build_interactive_audit_html__ann_has_content(ann: dict | None) -> bool:
    if not ann:
        return False
    return bool(str(ann.get('human_verdict') or '').strip() or str(ann.get('human_comment') or '').strip() or str(ann.get('human_fail_stage') or '').strip() or str(ann.get('human_genomic_region') or '').strip())

def build_interactive_audit_html_overlay_annotations_on_records(records: list[dict], annotations: list[dict] | None) -> int:
    """Copy human_* from annotations onto records.

    Match audit_id first; if that record has no verdict, fall back to
    pmid + gene + variant_text (optionally sample_class). Original
    annotation rows are not removed.
    """
    if not annotations:
        return 0
    by_id = {str(a.get('audit_id') or ''): a for a in annotations if a.get('audit_id')}
    by_key: dict[tuple[str, str, str, str], list[dict]] = {}
    for ann in annotations:
        by_key.setdefault(build_interactive_audit_html_rec_overlay_key(ann), []).append(ann)
        by_key.setdefault(('', build_interactive_audit_html_rec_overlay_key(ann)[1], build_interactive_audit_html_rec_overlay_key(ann)[2], build_interactive_audit_html_rec_overlay_key(ann)[3]), []).append(ann)
    used: set[int] = set()
    n = 0
    for rec in records:
        aid = str(rec.get('audit_id') or '')
        picked = by_id.get(aid) if build_interactive_audit_html__ann_has_content(by_id.get(aid)) else None
        if picked is None:
            key = build_interactive_audit_html_rec_overlay_key(rec)
            for cand in by_key.get(key) or by_key.get(('', key[1], key[2], key[3])) or []:
                if id(cand) in used:
                    continue
                if build_interactive_audit_html__ann_has_content(cand):
                    picked = cand
                    break
        if picked is None:
            continue
        used.add(id(picked))
        rec['human_verdict'] = picked.get('human_verdict') or rec.get('human_verdict') or ''
        rec['human_fail_stage'] = picked.get('human_fail_stage') or rec.get('human_fail_stage') or ''
        rec['human_comment'] = picked.get('human_comment') or rec.get('human_comment') or ''
        rec['human_genomic_region'] = picked.get('human_genomic_region') or rec.get('human_genomic_region') or ''
        n += 1
    return n

def build_interactive_audit_html_slim_record(rec: dict) -> dict:
    """Keep fields needed for interactive review (shrink embed size)."""
    stages = {f'stage{i}_text': rec.get(f'stage{i}_text') or '' for i in range(10)}
    return {'audit_id': rec.get('audit_id') or '', 'sample_class': rec.get('sample_class') or '', 'pmid': rec.get('pmid') or '', 'gene': rec.get('gene') or '', 'variant_text': rec.get('variant_text') or '', 'label': rec.get('label') or '', 'link_tier': rec.get('link_tier') or '', 'review_tier': rec.get('review_tier') or '', 'screen_grade': rec.get('screen_grade') or '', 'chromosome_ucsc': rec.get('chromosome_ucsc') or '', 'position': rec.get('position') or '', 'ref': rec.get('ref') or '', 'alt': rec.get('alt') or '', 'c_hgvs': rec.get('c_hgvs') or '', 'p_hgvs': rec.get('p_hgvs') or '', 'resolve_status': rec.get('resolve_status') or '', 'resolve_method': rec.get('resolve_method') or '', 'pubmed_url': rec.get('pubmed_url') or '', 'genomic_region': rec.get('genomic_region') or '', 'genomic_region_genes': rec.get('genomic_region_genes') or '', 'genomic_region_feature': rec.get('genomic_region_feature') or '', 'genomic_region_transcripts': rec.get('genomic_region_transcripts') or '', 'genomic_region_note': rec.get('genomic_region_note') or '', 'human_checklist_text': rec.get('human_checklist_text') or '', 'human_verdict': rec.get('human_verdict') or '', 'human_fail_stage': rec.get('human_fail_stage') or '', 'human_comment': rec.get('human_comment') or '', 'human_genomic_region': rec.get('human_genomic_region') or '', **stages, 'reasoning_chain_text': rec.get('reasoning_chain_text') or ''}

def build_interactive_audit_html_load_jsonl(path: Path) -> list[dict]:
    rows: list[dict] = []
    with path.open(encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rows.append(build_interactive_audit_html_slim_record(json.loads(line)))
    return rows

def build_interactive_audit_html_build_html(template: str, embedded: list[dict] | None, embedded_annotations: list[dict] | None, title: str, save_filename: str='') -> str:
    payload = json.dumps(embedded or [], ensure_ascii=False)
    payload = payload.replace('<', '\\u003c').replace('>', '\\u003e')
    ann_payload = json.dumps(embedded_annotations or [], ensure_ascii=False)
    ann_payload = ann_payload.replace('<', '\\u003c').replace('>', '\\u003e')
    return template.replace('{{TITLE}}', html.escape(title)).replace('{{SAVE_FILENAME}}', html.escape(save_filename)).replace('/*{{EMBEDDED_JSON}}*/', payload).replace('/*{{EMBEDDED_ANNOTATIONS_JSON}}*/', ann_payload).replace('{{EMBED_COUNT}}', str(len(embedded or [])))

def build_interactive_audit_html_main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--template', type=Path, default=build_interactive_audit_html_TEMPLATE)
    ap.add_argument('--embed', type=Path, default=None, help='JSONL to bake into HTML')
    ap.add_argument('--embed-annotations', type=Path, default=None, help='audit_annotations.json to bake human_verdict/comment back into HTML')
    ap.add_argument('--out', type=Path, default=build_interactive_audit_html_DEFAULT_OUT)
    ap.add_argument('--title', type=str, default='Variant clinical Stage0–9 交互审核')
    args = ap.parse_args()
    template = args.template.read_text(encoding='utf-8')
    embedded = build_interactive_audit_html_load_jsonl(args.embed) if args.embed else None
    if args.embed and (not args.out.name.startswith('review_')):
        if args.out == build_interactive_audit_html_DEFAULT_OUT:
            args.out = args.embed.with_name(f'review_{args.embed.stem}.html')
    embedded_annotations = None
    if args.embed_annotations:
        loaded = json.loads(args.embed_annotations.read_text(encoding='utf-8'))
        embedded_annotations = loaded if isinstance(loaded, list) else None
        if embedded and embedded_annotations:
            build_interactive_audit_html_overlay_annotations_on_records(embedded, embedded_annotations)
    html_text = build_interactive_audit_html_build_html(template, embedded, embedded_annotations, args.title, save_filename=args.out.name)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(html_text, encoding='utf-8')
    n = len(embedded or [])
    print(f'wrote {args.out} ({args.out.stat().st_size:,} bytes, embedded={n})')

# === export_review_pass_csv.py ===
"""Export human_verdict=pass rows from the new-site review set.

Joins new_site JSONL + chains with audit_annotations.json (latest mirror).
Writes all-pass CSVs for positive and negative, plus coding/noncoding splits
for positive. Does not rewrite audit_annotations.json.
"""
import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any
export_review_pass_csv_ROOT = Path(__file__).resolve().parents[1]
export_review_pass_csv_SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(export_review_pass_csv_SCRIPTS))
export_review_pass_csv_TRAILS_NEW = filter_new_site_variants_NEW_SITE

def export_review_pass_csv_load_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open(encoding='utf-8') as f:
        for line in f:
            if line.strip():
                rows.append(json.loads(line))
    return rows

def export_review_pass_csv_merge_annotations(json_path: Path, html_paths: list[Path]) -> list[dict[str, Any]]:
    by_id: dict[str, dict[str, Any]] = {}

    def ingest(arr: list[Any]) -> None:
        for row in arr:
            if not isinstance(row, dict) or not row.get('audit_id'):
                continue
            aid = str(row['audit_id'])
            prev = by_id.get(aid)
            if prev is None:
                by_id[aid] = row
                continue
            if str(row.get('updated_at') or '') >= str(prev.get('updated_at') or ''):
                by_id[aid] = row
    if json_path.exists():
        loaded = json.loads(json_path.read_text(encoding='utf-8'))
        if isinstance(loaded, list):
            ingest(loaded)
    import re
    ann_re = re.compile('<script\\s+id=["\\\']embedded-annotations["\\\'][^>]*>(.*?)</script>', re.I | re.S)
    for html in html_paths:
        if not html.exists():
            continue
        m = ann_re.search(html.read_text(encoding='utf-8'))
        if not m:
            continue
        raw = (m.group(1) or '').strip()
        if not raw or raw in {'[]', '/*{{EMBEDDED_ANNOTATIONS_JSON}}*/'}:
            continue
        try:
            arr = json.loads(raw)
        except json.JSONDecodeError:
            continue
        if isinstance(arr, list):
            ingest(arr)
    return list(by_id.values())

def export_review_pass_csv_collect_pass(jsonl_rows: list[dict[str, Any]], chains: list[dict[str, Any]], annotations: list[dict[str, Any]], allowed_regions: set[str] | None=None) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    by_id = {str(a.get('audit_id')): a for a in annotations if a.get('audit_id')}
    by_key: dict[tuple, list[dict[str, Any]]] = defaultdict(list)
    for a in annotations:
        by_key[filter_new_site_variants_rec_key(a)].append(a)
        by_key['', str(a.get('pmid') or ''), str(a.get('gene') or '').upper(), str(a.get('variant_text') or '').strip()].append(a)
    chain_q: dict[tuple, list[dict[str, Any]]] = defaultdict(list)
    for rec in chains:
        chain_q[filter_new_site_variants_rec_key(rec)].append(rec)
    used: set[str] = set()
    exported: list[dict[str, Any]] = []
    scanned = Counter()
    skipped_region = 0
    for row in jsonl_rows:
        key = filter_new_site_variants_rec_key(row)
        chain = chain_q[key].pop(0) if chain_q[key] else None
        audit_id = str((chain or {}).get('audit_id') or '')
        ann = by_id.get(audit_id) if audit_id else None
        if not ann or not (ann.get('human_verdict') or ann.get('human_comment') or ann.get('human_genomic_region')):
            for cand in by_key.get(key) or by_key.get(('', key[1], key[2], key[3])) or []:
                cid = str(cand.get('audit_id') or '')
                if cid in used:
                    continue
                ann = cand
                break
        if not ann:
            continue
        used.add(str(ann.get('audit_id') or ''))
        verd = str(ann.get('human_verdict') or '').strip().lower()
        scanned[verd or 'pending'] += 1
        if verd != 'pass':
            continue
        used_region = filter_new_site_variants_region_used(row, ann)
        if allowed_regions is not None and used_region not in allowed_regions:
            skipped_region += 1
            continue
        rec = dict(row)
        rec.update({'audit_id': audit_id or ann.get('audit_id') or '', 'human_verdict': ann.get('human_verdict') or '', 'human_genomic_region': ann.get('human_genomic_region') or '', 'human_comment': ann.get('human_comment') or '', 'human_fail_stage': ann.get('human_fail_stage') or '', 'updated_at': ann.get('updated_at') or '', 'region_used': used_region})
        exported.append(rec)
    complete = sum((1 for r in exported if str(r.get('chromosome_ucsc') or '').strip() and str(r.get('position') or '').strip() and str(r.get('ref') or '').strip() and str(r.get('alt') or '').strip()))
    stats = {'n': len(exported), 'n_complete_allele': complete, 'n_missing_allele': len(exported) - complete, 'skipped_region': skipped_region, 'by_region': dict(Counter((r.get('region_used') or '∅' for r in exported))), 'by_link_tier': dict(Counter((r.get('link_tier') for r in exported))), 'by_resolve_status': dict(Counter((r.get('resolve_status') for r in exported))), 'by_resolve_method': dict(Counter((r.get('resolve_method') for r in exported))), 'matched_annotation_verdicts': dict(scanned)}
    return (exported, stats)

def export_review_pass_csv_write_stats(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

def export_review_pass_csv_main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--annotations', type=Path, default=filter_new_site_variants_MANUAL / 'audit_annotations.json')
    ap.add_argument('--out-dir', type=Path, default=filter_new_site_variants_MANUAL)
    args = ap.parse_args()
    annotations = export_review_pass_csv_merge_annotations(args.annotations, [filter_new_site_variants_MANUAL / 'review_positive.html', filter_new_site_variants_MANUAL / 'review_negative.html'])
    pos_rows = export_review_pass_csv_load_jsonl(filter_new_site_variants_NEW_SITE / 'variants_genomic_positive_new_site.jsonl')
    neg_rows = export_review_pass_csv_load_jsonl(filter_new_site_variants_NEW_SITE / 'variants_genomic_negative_new_site.jsonl')
    pos_chains = export_review_pass_csv_load_jsonl(export_review_pass_csv_TRAILS_NEW / 'stage_audit_chains_positive.jsonl')
    neg_chains = export_review_pass_csv_load_jsonl(export_review_pass_csv_TRAILS_NEW / 'stage_audit_chains_negative.jsonl')
    jobs = [('review_positive_pass.csv', pos_rows, pos_chains, None, 'positive all pass on new-site set'), ('review_positive_pass_coding.csv', pos_rows, pos_chains, filter_new_site_variants_CODING, 'positive pass, region_used=coding'), ('review_positive_pass_noncoding.csv', pos_rows, pos_chains, filter_new_site_variants_NONCODING, 'positive pass, noncoding region_used'), ('review_negative_pass.csv', neg_rows, neg_chains, None, 'negative all pass on new-site set')]
    report: dict[str, Any] = {'annotations_n': len(annotations), 'source_jsonl_positive': str(filter_new_site_variants_NEW_SITE / 'variants_genomic_positive_new_site.jsonl'), 'source_jsonl_negative': str(filter_new_site_variants_NEW_SITE / 'variants_genomic_negative_new_site.jsonl'), 'exports': {}}
    for filename, rows, chains, regions, note in jobs:
        exported, stats = export_review_pass_csv_collect_pass(rows, chains, annotations, regions)
        out_csv = args.out_dir / filename
        filter_new_site_variants_write_csv(out_csv, exported, filter_new_site_variants_PASS_CSV_COLS)
        payload = {'output_csv': str(out_csv), 'note': note, 'filter': {'human_verdict': 'pass', 'region_used': sorted(regions) if regions else 'any (including empty)', 'region_source': 'human_genomic_region if set else pipeline genomic_region', 'candidate_set': 'new_site (rs-only dropped unless literature c.(rs) pair)'}, 'exported': stats}
        export_review_pass_csv_write_stats(out_csv.with_name(out_csv.stem + '_stats.json'), payload)
        report['exports'][filename] = stats
        print(f"wrote {out_csv} n={stats['n']} complete_allele={stats['n_complete_allele']} region={stats['by_region']}")
    export_review_pass_csv_write_stats(args.out_dir / 'review_pass_export_report.json', report)

# === check_grch38_coord_refalt.py ===
"""
GRCh38 coordinate / REF-ALT compliance checker.

Input columns (first 4 of variants_genomic_all_*.csv):
  chromosome_ucsc, position, ref, alt

Method (VCF-style genomic alleles on GRCh38 primary/toplevel DNA):
  A. Parse & normalize chromosome (chr1/1/X/MT) and 1-based position.
  B. Syntax checks on REF/ALT (IUPAC ACGTN; non-empty when coords present).
  C. Fetch FASTA[chrom][pos : pos+len(REF)-1] and require exact REF match.
  D. Structural allele rules (hard vs soft):
       HARD: REF != ALT; bases in ACGTN
       SOFT warn: unequal-length delins not expressible as VCF prefix indel
       INFO tags: snv / mnv / del_vcf / ins_vcf / dup_tandem / delins
  E. Coordinate bounds vs FASTA contig length.

Statuses:
  PASS          — hard checks ok (FASTA REF match + syntax + coords)
  PASS_WARN     — PASS but allele not in strict VCF left-aligned form
  SKIP_EMPTY    — missing chrom/pos/ref (unresolved genomic map)
  FAIL_SYNTAX   — alphabet / empty / REF==ALT
  FAIL_REF      — FASTA sequence != CSV REF
  FAIL_COORD    — chrom unknown or out of contig range
  FAIL_FETCH    — I/O or unexpected reader error

Usage:
  python check_grch38_coord_refalt.py \\
    --in-csv .../variants_genomic_all_ref_context5_with_fixops_applied.csv \\
    --out-tsv .../coord_refalt_check.tsv
"""
import argparse
import csv
import re
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
check_grch38_coord_refalt_ACGTN = re.compile('^[ACGTN]+$', re.I)
check_grch38_coord_refalt_CONTIG_ALIASES = {'M': 'MT', 'CHRM': 'MT', 'CHRMT': 'MT', 'MT': 'MT'}

@dataclass
class check_grch38_coord_refalt_CheckResult:
    row_index: int
    chromosome_ucsc: str
    position: str
    ref: str
    alt: str
    chrom_norm: str
    ref_span_start: str
    ref_span_end: str
    fasta_ref: str
    status: str
    checks: str
    note: str

def check_grch38_coord_refalt_norm_chrom(raw: str) -> str:
    c = (raw or '').strip()
    if not c:
        return ''
    if c.lower().startswith('chr'):
        c = c[3:]
    c = c.upper()
    return check_grch38_coord_refalt_CONTIG_ALIASES.get(c, c)

def check_grch38_coord_refalt_classify_allele_geometry(ref: str, alt: str) -> tuple[str, str | None]:
    """
    Classify REF/ALT geometry.

    Returns (tag, warn_or_none).
    Hard syntax failures are handled before this is called (empty, alphabet, REF==ALT).

    Accepted as compliant without warning:
      - SNV (1bp/1bp)
      - MNV (same length >1)
      - VCF left-aligned del (REF.startswith(ALT), shared prefix == len(ALT) == 1 for single-base ALT,
        or more generally ALT proper prefix of REF with shared == len(ALT))
      - VCF left-aligned ins (ALT.startswith(REF), shared == len(REF))
      - Tandem dup interval (ALT == REF + REF)

    Warning (still PASS_WARN if FASTA ok):
      - unequal-length alleles that are not a clean prefix relationship (delins-like)
      - indel with shared prefix length != 1 when expressed as anchor+change
        (except tandem dup / pure prefix del/ins forms above)
    """
    if len(ref) == 1 and len(alt) == 1:
        return ('snv', None)
    if len(ref) == len(alt):
        return ('mnv', None)
    if alt == ref + ref:
        return ('dup_tandem', None)
    if len(alt) < len(ref) and ref.startswith(alt):
        return ('del_vcf', None)
    if len(ref) < len(alt) and alt.startswith(ref):
        return ('ins_vcf', None)
    shared = 0
    n = min(len(ref), len(alt))
    while shared < n and ref[shared] == alt[shared]:
        shared += 1
    if shared == 0:
        return ('delins', 'delins_no_shared_prefix')
    if shared > 1:
        return ('delins', f'delins_shared_prefix={shared}_not_strict_vcf_left_align')
    return ('delins', 'delins_unequal_non_prefix')

def check_grch38_coord_refalt_check_row(row_index: int, chrom_ucsc: str, pos_s: str, ref: str, alt: str, fasta: vc_ensembl_c_map_module.FastaReader) -> check_grch38_coord_refalt_CheckResult:
    chrom_ucsc = (chrom_ucsc or '').strip()
    pos_s = (pos_s or '').strip()
    ref = (ref or '').strip().upper()
    alt = (alt or '').strip().upper()
    base = check_grch38_coord_refalt_CheckResult(row_index=row_index, chromosome_ucsc=chrom_ucsc, position=pos_s, ref=ref, alt=alt, chrom_norm='', ref_span_start='', ref_span_end='', fasta_ref='', status='', checks='', note='')
    if not chrom_ucsc or not pos_s or (not ref):
        base.status = 'SKIP_EMPTY'
        base.checks = 'missing_chrom_or_pos_or_ref'
        base.note = 'no genomic allele to validate'
        return base
    chrom = check_grch38_coord_refalt_norm_chrom(chrom_ucsc)
    base.chrom_norm = chrom
    flags: list[str] = []
    try:
        pos = int(pos_s)
    except ValueError:
        base.status = 'FAIL_SYNTAX'
        base.checks = 'position_not_integer'
        return base
    if pos < 1:
        base.status = 'FAIL_COORD'
        base.checks = 'position_lt_1'
        return base
    if not check_grch38_coord_refalt_ACGTN.match(ref):
        base.status = 'FAIL_SYNTAX'
        base.checks = 'ref_non_acgtn'
        return base
    if not alt:
        base.status = 'FAIL_SYNTAX'
        base.checks = 'alt_empty'
        return base
    if not check_grch38_coord_refalt_ACGTN.match(alt):
        base.status = 'FAIL_SYNTAX'
        base.checks = 'alt_non_acgtn'
        return base
    if ref == alt:
        base.status = 'FAIL_SYNTAX'
        base.checks = 'ref_equals_alt'
        return base
    tag, warn = check_grch38_coord_refalt_classify_allele_geometry(ref, alt)
    flags.append(tag)
    if warn:
        flags.append(warn)
    ref_start = pos
    ref_end = pos + len(ref) - 1
    base.ref_span_start = str(ref_start)
    base.ref_span_end = str(ref_end)
    if chrom not in fasta.index:
        base.status = 'FAIL_COORD'
        base.checks = 'chrom_not_in_fasta'
        base.note = f'available_contigs_sample={list(fasta.index)[:5]}...'
        return base
    contig_len = fasta.index[chrom][0]
    if ref_end > contig_len:
        base.status = 'FAIL_COORD'
        base.checks = 'ref_span_past_contig_end'
        base.note = f'contig_len={contig_len}'
        return base
    try:
        fasta_ref = fasta.fetch(chrom, ref_start, ref_end)
    except Exception as e:
        base.status = 'FAIL_FETCH'
        base.checks = type(e).__name__
        base.note = str(e)
        return base
    base.fasta_ref = fasta_ref
    if fasta_ref != ref:
        base.status = 'FAIL_REF'
        base.checks = 'fasta_ref_mismatch'
        base.note = f'expected={ref} observed={fasta_ref}'
        return base
    flags.append('fasta_ref_match')
    base.checks = ','.join(flags)
    base.status = 'PASS_WARN' if warn else 'PASS'
    return base

def check_grch38_coord_refalt_main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--in-csv', type=Path, required=True)
    ap.add_argument('--out-tsv', type=Path, required=True)
    ap.add_argument('--fasta', type=Path, default=vc_ensembl_c_map_module.DEFAULT_FASTA)
    ap.add_argument('--fai', type=Path, default=vc_ensembl_c_map_module.DEFAULT_FAI)
    ap.add_argument('--chrom-col', default='chromosome_ucsc', help='Chromosome column (default: first-col name chromosome_ucsc)')
    ap.add_argument('--pos-col', default='position')
    ap.add_argument('--ref-col', default='ref')
    ap.add_argument('--alt-col', default='alt')
    args = ap.parse_args()
    fasta = vc_ensembl_c_map_module.FastaReader(args.fasta, args.fai)
    results: list[check_grch38_coord_refalt_CheckResult] = []
    with args.in_csv.open(encoding='utf-8', newline='') as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader, start=1):
            results.append(check_grch38_coord_refalt_check_row(row_index=i, chrom_ucsc=row.get(args.chrom_col, ''), pos_s=row.get(args.pos_col, ''), ref=row.get(args.ref_col, ''), alt=row.get(args.alt_col, ''), fasta=fasta))
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    fields = list(asdict(results[0]).keys()) if results else ['row_index', 'chromosome_ucsc', 'position', 'ref', 'alt', 'chrom_norm', 'ref_span_start', 'ref_span_end', 'fasta_ref', 'status', 'checks', 'note']
    with args.out_tsv.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        for r in results:
            w.writerow(asdict(r))
    from collections import Counter
    counts = Counter((r.status for r in results))
    print(f'wrote {args.out_tsv}')
    print(f'rows={len(results)}  ' + '  '.join((f'{k}={v}' for k, v in sorted(counts.items()))))
    fails = [r for r in results if r.status.startswith('FAIL')]
    if fails:
        print('failures:')
        for r in fails[:30]:
            print(f'  row={r.row_index} {r.chromosome_ucsc}:{r.position} {r.ref}>{r.alt}  {r.status}  {r.checks}  {r.note}')
        if len(fails) > 30:
            print(f'  ... +{len(fails) - 30} more')

# === validate_grch38_ref_context5.py ===
"""
Validate GRCh38-ref allele from exported CSV by re-fetching from FASTA.

For each row (chromosome, position, ref, alt) it:
  1) Treats `position` as 1-based start coordinate of `ref` span.
  2) Fetches:
       - 5bp upstream of ref start
       - ref span
       - 5bp downstream of ref end
  3) Formats each base as: <BASE>[<GENOMIC_POS>]
  4) Checks fetched ref-span == CSV ref, and marks ok/bad.

Outputs:
  - <out>/variants_genomic_all_ref_context5.csv

Notes:
  - For indels, the context is anchored to the genomic REF span coordinates.
  - ALT is not re-derived from FASTA (since reference genome only stores REF).
    This tool is strictly a "does CSV REF match FASTA at that coordinate?" validator.
"""
import argparse
import csv
from pathlib import Path
from typing import Any
import logging
validate_grch38_ref_context5_ROOT = Path(__file__).resolve().parents[1]

def validate_grch38_ref_context5_norm_chrom(raw: str) -> str:
    """Normalize bare chromosomes like 'chr1' -> '1', 'chrX'->'X'."""
    c = (raw or '').strip()
    if c.lower().startswith('chr'):
        c = c[3:]
    c = c.upper()
    if c in {'M', 'MT', 'CHRMT'}:
        return 'MT'
    return c

def validate_grch38_ref_context5_format_seq_with_positions(fwd_bases: str, start_pos_1based: int) -> str:
    """
    fwd_bases is genomic forward string (5'->3' on chrom).
    start_pos_1based is the genomic coordinate of fwd_bases[0].
    """
    parts: list[str] = []
    for i, b in enumerate(fwd_bases):
        pos = start_pos_1based + i
        parts.append(f'{b}[{pos}]')
    return ''.join(parts)

def validate_grch38_ref_context5_main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--in-csv', type=Path, required=True, help='Input variants_genomic_all.csv')
    ap.add_argument('--out-csv', type=Path, required=True, help='Output CSV with ref-context and FASTA check result')
    ap.add_argument('--fasta', type=Path, default=vc_ensembl_c_map_module.DEFAULT_FASTA)
    ap.add_argument('--fai', type=Path, default=vc_ensembl_c_map_module.DEFAULT_FAI)
    ap.add_argument('--ensembl-gff3', type=Path, default=vc_ensembl_c_map_module.DEFAULT_GFF3)
    ap.add_argument('--ensembl-cache', type=Path, default=validate_grch38_ref_context5_ROOT / 'data/pubmed/cache/ensembl/ensembl_mane_transcript_index.json')
    ap.add_argument('--upstream', type=int, default=5, help='upstream bp count')
    ap.add_argument('--downstream', type=int, default=5, help='downstream bp count')
    args = ap.parse_args()
    logger = logging.getLogger('validate_grch38_ref_context5')
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    fasta = vc_ensembl_c_map_module.FastaReader(args.fasta, args.fai)
    mane_idx = None

    def ensure_mane_idx() -> None:
        nonlocal mane_idx
        if mane_idx is not None:
            return
        mane_idx = vc_ensembl_c_map_module.build_mane_index(args.ensembl_gff3, logger, args.ensembl_cache)

    def chrom_ucsc(bare: str) -> str:
        if bare == 'MT':
            return 'chrM'
        return f'chr{bare}'
    with args.in_csv.open(encoding='utf-8') as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        rows_in = list(reader)
    extra = ['ref_span_start_1based', 'ref_span_end_1based', 'ref_ok_fasta', 'context_ref_5bp', 'fix_operation', 'chromosome_ucsc_before', 'position_before', 'ref_before', 'alt_before', 'chromosome_ucsc_after', 'position_after', 'ref_after', 'alt_after']
    out_fields = fieldnames + [c for c in extra if c not in fieldnames]
    rows_out: list[dict[str, Any]] = []
    bad = 0
    for r in rows_in:
        chrom_raw = (r.get('chromosome') or '').strip()
        pos_s = (r.get('position') or '').strip()
        ref = (r.get('ref') or '').strip().upper()
        chrom_ucsc_in = (r.get('chromosome_ucsc') or '').strip()
        pos_in = (r.get('position') or '').strip()
        ref_in = (r.get('ref') or '').strip().upper()
        alt_in = (r.get('alt') or '').strip().upper()
        chromosome_ucsc_before = chrom_ucsc_in
        position_before = pos_in
        ref_before = ref_in
        alt_before = alt_in
        chromosome_ucsc_after = chromosome_ucsc_before
        position_after = position_before
        ref_after = ref_before
        alt_after = alt_before
        pmid = str(r.get('pmid') or '').strip()
        gene = (r.get('gene') or '').strip().upper()
        resolve_method = (r.get('resolve_method') or '').strip()
        p_hgvs = (r.get('p_hgvs') or '').strip()
        expected = {'p.A991D': 'c.2972C>A', 'p.E840V': 'c.2519A>T', 'p.P1386S': 'c.4156C>T', 'p.T994I': 'c.2981C>T'}
        if pmid == '42418840' and gene == 'HMNX' and (resolve_method == 'pdf_literal_genomic') and (p_hgvs in expected):
            try:
                ensure_mane_idx()
                tx = mane_idx.by_gene.get('ATP7A')
                if tx:
                    mapped = vc_ensembl_c_map_module.map_c_to_genomic(tx=tx, fasta=fasta, c_hgvs=expected[p_hgvs])
                    chromosome_ucsc_after = chrom_ucsc(str(mapped['chromosome']))
                    position_after = str(mapped['position'])
                    ref_after = (mapped.get('ref') or '').upper()
                    alt_after = (mapped.get('alt') or '').upper()
            except Exception:
                pass

        def compute_fix_operation(row: dict[str, Any]) -> str:
            """
            Create a human-readable "what fix happened" note.

            Priority:
              1) Use explicit resolve_note if it contains known fix markers.
              2) Otherwise infer from resolve_method + ref/alt geometry.
            """
            ops: list[str] = []
            rn = (row.get('resolve_note') or '').strip()
            if rn:
                rn_short = rn
                if len(rn_short) > 140:
                    rn_short = rn_short[:137] + '...'
                if any((m in rn.lower() for m in ['gene corrected', 'normalized hyphen', 'vcf left-aligned del', 'vcf left-aligned', 'left-aligned del', 'left aligned', 'dup', 'range'])):
                    ops.append(rn_short)
            ref2 = (row.get('ref') or '').strip().upper()
            alt2 = (row.get('alt') or '').strip().upper()
            c_hgvs = (row.get('c_hgvs') or '').strip()
            resolve_method = (row.get('resolve_method') or '').strip()
            gene = (row.get('gene') or '').strip().upper()
            pmid = str(row.get('pmid') or '').strip()
            p_hgvs = (row.get('p_hgvs') or '').strip()
            if pmid == '42418840' and gene == 'HMNX' and (resolve_method == 'pdf_literal_genomic'):
                expected: dict[str, str] = {'p.A991D': 'c.2972C>A', 'p.E840V': 'c.2519A>T', 'p.P1386S': 'c.4156C>T', 'p.T994I': 'c.2981C>T'}
                if p_hgvs in expected and c_hgvs != expected[p_hgvs]:
                    ops.append(f'Table1 PMID42418840: {p_hgvs} should map to {expected[p_hgvs]} (currently {c_hgvs}); re-lift via Ensembl MANE p-map for correct genomic REF/ALT')
            if ref2 and alt2 and (len(ref2) > len(alt2)) and ref2.startswith(alt2):
                ops.append('del as VCF left-aligned (REF=anchor+deleted, ALT=anchor)')
            elif ref2 and alt2 and (len(alt2) == 2 * len(ref2)) and (alt2 == ref2 + ref2):
                if 'dup' in c_hgvs.lower() or 'dup' in (resolve_method.lower() if resolve_method else ''):
                    ops.append('dup interval as tandem (REF=dup, ALT=REF+REF)')
            elif ref2 and alt2 and (len(ref2) == 1) and (len(alt2) == 1):
                if resolve_method == 'ensembl_mane_p_map':
                    ops.append('p.-lift emitted SNV at differing coding base (single-base REF/ALT)')
            if not ops:
                ops.append(f'resolve_method={resolve_method}' if resolve_method else 'no_fix_operation_inferred')
            seen = set()
            out_ops = []
            for x in ops:
                if x and x not in seen:
                    seen.add(x)
                    out_ops.append(x)
            return '; '.join(out_ops)
        if not chrom_raw or not pos_s or (not ref):
            r2 = dict(r)
            r2.update({'ref_span_start_1based': '', 'ref_span_end_1based': '', 'ref_ok_fasta': 'NA', 'context_ref_5bp': '', 'fix_operation': compute_fix_operation(r), 'chromosome_ucsc_before': chromosome_ucsc_before, 'position_before': position_before, 'ref_before': ref_before, 'alt_before': alt_before, 'chromosome_ucsc_after': chromosome_ucsc_after, 'position_after': position_after, 'ref_after': ref_after, 'alt_after': alt_after})
            rows_out.append(r2)
            continue
        chrom = validate_grch38_ref_context5_norm_chrom(chrom_raw) or chrom_raw
        pos = int(pos_s)
        ref_len = len(ref)
        ref_start = pos
        ref_end = pos + ref_len - 1
        up_start = max(1, ref_start - args.upstream)
        up_end = ref_start - 1
        dn_start = ref_end + 1
        dn_end = ref_end + args.downstream
        try:
            up = fasta.fetch(chrom, up_start, up_end) if up_end >= up_start else ''
            ref_obs = fasta.fetch(chrom, ref_start, ref_end)
            dn = fasta.fetch(chrom, dn_start, dn_end) if dn_end >= dn_start else ''
            ok = ref_obs.upper() == ref
            if not ok:
                bad += 1
            context = ''
            if up:
                context += validate_grch38_ref_context5_format_seq_with_positions(up, up_start)
            context += validate_grch38_ref_context5_format_seq_with_positions(ref_obs, ref_start)
            if dn:
                context += validate_grch38_ref_context5_format_seq_with_positions(dn, dn_start)
            r2 = dict(r)
            r2.update({'ref_span_start_1based': ref_start, 'ref_span_end_1based': ref_end, 'ref_ok_fasta': 'OK' if ok else 'BAD', 'context_ref_5bp': context, 'fix_operation': compute_fix_operation(r), 'chromosome_ucsc_before': chromosome_ucsc_before, 'position_before': position_before, 'ref_before': ref_before, 'alt_before': alt_before, 'chromosome_ucsc_after': chromosome_ucsc_after, 'position_after': position_after, 'ref_after': ref_after, 'alt_after': alt_after})
            rows_out.append(r2)
        except Exception as e:
            r2 = dict(r)
            r2.update({'ref_span_start_1based': ref_start, 'ref_span_end_1based': ref_end, 'ref_ok_fasta': f'ERR:{type(e).__name__}', 'context_ref_5bp': '', 'fix_operation': compute_fix_operation(r), 'chromosome_ucsc_before': chromosome_ucsc_before, 'position_before': position_before, 'ref_before': ref_before, 'alt_before': alt_before, 'chromosome_ucsc_after': chromosome_ucsc_after, 'position_after': position_after, 'ref_after': ref_after, 'alt_after': alt_after})
            rows_out.append(r2)
            bad += 1
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=out_fields, extrasaction='ignore')
        w.writeheader()
        for r in rows_out:
            w.writerow(r)
    logger.info('wrote %s', args.out_csv)
    logger.info('bad/ref mismatch/err rows: %s / %s', bad, len(rows_out))

# Public helper compatibility aliases.

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
        'single': run_single_pmid_pipeline_main,
        'download-pdfs': download_review_strict_pdfs_main,
        'new-site': filter_new_site_variants_main,
        'audit': build_stage_audit_chains_jsonl_main,
        'html': build_interactive_audit_html_main,
        'review-export': export_review_pass_csv_main,
        'refalt': check_grch38_coord_refalt_main,
        'ref-context5': validate_grch38_ref_context5_main,
    })
