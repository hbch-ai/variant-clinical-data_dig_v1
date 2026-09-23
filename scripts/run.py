#!/usr/bin/env python3
"""Unified Stage0-9 runner and maintenance entry point."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from vc_paths import (
    CLINVAR_SUMMARY, ENSEMBL_FAI, ENSEMBL_FASTA, ENSEMBL_GFF3,
    MERGED_ARTICLES, NEGATIVES_FULL_PASS, PDF_PATHO, PDF_TEXT_CACHE,
    REVIEW_HIGHCONF, REVIEW_HIGHCONF_ANN, REVIEW_STRICT, S2, S4,
    S5_ENRICHED, S6, S9,
)

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = Path(__file__).resolve().parent
STAGE_ORDER = [
    'stage0_preflight', 'stage1_pubmed_retrieval', 'stage2_abstract_linking',
    'stage3_fulltext_enrichment', 'stage4_discovery_tagging',
    'stage5_clinvar_validation', 'stage6_review_quality',
    'stage7_journal_priority', 'stage8_pdf_evidence', 'stage8_negatives',
    'stage9_genomic_export',
]
QUERIES = ['q1_core_gpc', 'q3_case_reports', 'q2_high_recall', 'q4_mesh',
           'q5_clinical_genetics', 'q6_rsid_trait', 'q7_intergenic_regulatory']


def n_jsonl(path: Path) -> int:
    if not path.is_file():
        return 0
    with path.open(encoding='utf-8') as handle:
        return sum(bool(line.strip()) for line in handle)


def ensure_data_link() -> None:
    """Bytecode expects scripts/data to resolve inside this repository."""
    link = SCRIPTS / 'data'
    if link.is_symlink() or link.exists():
        return
    link.symlink_to(Path('..') / 'data', target_is_directory=True)


def stage_env(cohort: bool) -> dict[str, str]:
    env = os.environ.copy()
    env['VARIANT_CLINICAL_ROOT'] = str(ROOT)
    env['PYTHONPATH'] = os.pathsep.join([str(SCRIPTS), env.get('PYTHONPATH', '')]).rstrip(os.pathsep)
    if cohort:
        env['VARIANT_CLINICAL_COHORT'] = '1'
    return env


def run_script(argv: list[str], log: Path, *, cohort: bool) -> None:
    """Run one stage in a fresh process so large indexes are released afterward."""
    cmd = [sys.executable, *argv]
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open('w', encoding='utf-8') as handle:
        handle.write('$ ' + ' '.join(cmd) + '\n')
        handle.flush()
        proc = subprocess.run(cmd, cwd=ROOT, env=stage_env(cohort), stdout=handle, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise RuntimeError(f'exit {proc.returncode}; see {log}')


# --- Stage1 fixed-cohort ingest (formerly cohort_ingest.py) ---
def load_pmids(path: Path) -> list[str]:
    values = []
    for line in path.read_text(encoding='utf-8').splitlines():
        value = line.strip()
        if value and not value.startswith('#'):
            values.append(value.split()[0])
    return sorted(set(values))


def ingest_cohort(pmid_list: Path, seed_from: Path | None, out: Path, fetch_missing: bool) -> None:
    pmids = set(load_pmids(pmid_list))
    found: dict[str, dict] = {}
    if seed_from and seed_from.is_file():
        for line in seed_from.open(encoding='utf-8'):
            if line.strip():
                row = json.loads(line)
                pmid = str(row.get('pmid', '')).strip()
                if pmid in pmids:
                    found[pmid] = row
    missing = sorted(pmids - set(found))
    if missing and fetch_missing:
        import pubmed_pipeline as pubmed_pipeline_module
        for pmid in missing:
            params = {'db': 'pubmed', 'id': pmid, 'rettype': 'xml', 'retmode': 'xml',
                      'retmax': '1', 'email': 'variant_clinical_v4@local.genos',
                      'tool': 'variant_clinical_v4_pipeline'}
            response = pubmed_pipeline_module.eutils_get('efetch.fcgi', params, None)
            import xml.etree.ElementTree as ET
            articles = list(ET.fromstring(response.text).findall('PubmedArticle'))
            if articles:
                row = pubmed_pipeline_module.parse_article(articles[0])
                row['source_query_ids'] = ['cohort_pmid_list']
                found[pmid] = row
        missing = sorted(pmids - set(found))
    if missing:
        raise RuntimeError(f'cohort fixture missing {len(missing)} PMIDs')
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('w', encoding='utf-8') as handle:
        for pmid in sorted(found):
            handle.write(json.dumps(found[pmid], ensure_ascii=False) + '\n')
    if len(pmids) != 117 or n_jsonl(out) != 117:
        raise RuntimeError(f'fixed cohort must contain 117 rows, got {n_jsonl(out)}')
    print(f'cohort Stage1 done: n=117 path={out}')


# --- Stage3 evidence stamping (formerly stamp_evidence_source.py) ---
def stamp_evidence_source(linked_dir: Path = S2) -> None:
    abstract = ['pairs.jsonl', 'pairs_strict.jsonl', 'corpus.jsonl', 'corpus_strict.jsonl']
    fulltext = ['pairs_fulltext.jsonl', 'corpus_fulltext.jsonl']
    enriched = ['pairs_enriched.jsonl', 'corpus_enriched.jsonl']
    for name, forced in [(n, 'abstract') for n in abstract] + [(n, 'fulltext') for n in fulltext] + [(n, None) for n in enriched]:
        path = linked_dir / name
        if not path.is_file():
            continue
        rows = []
        for line in path.open(encoding='utf-8'):
            if not line.strip():
                continue
            row = json.loads(line)
            source = forced or row.get('evidence_source') or row.get('text_scope')
            if source not in {'abstract', 'fulltext'}:
                source = 'fulltext' if row.get('ft_source') or row.get('section') else 'abstract'
            row['evidence_source'] = row['text_scope'] = source
            text = row.get('text')
            if source == 'abstract' and isinstance(text, str) and '[abstract]' not in text and '[fulltext:' not in text:
                row['text'] = f'[abstract] {text}'
            rows.append(row)
        with path.open('w', encoding='utf-8') as handle:
            for row in rows:
                handle.write(json.dumps(row, ensure_ascii=False) + '\n')


# --- Explicit cache seed (never implicit in cohort mode) ---
def seed_cache(source_root: Path, pmids: set[str]) -> int:
    src = source_root / 'data/pubmed'
    dst = ROOT / 'data/pubmed'
    copied = 0
    for pmid in pmids:
        for rel in (f'cache/fulltext/parsed/{pmid}.json', f'cache/pdfs/pdf_text_cache/PMID{pmid}.txt', f'cache/pdfs/PMID{pmid}.pdf'):
            source, target = src / rel, dst / rel
            if source.is_file():
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
                copied += 1
        raw = src / 'cache/fulltext/raw'
        if raw.is_dir():
            for source in raw.glob(f'*{pmid}*'):
                target = dst / 'cache/fulltext/raw' / source.name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
                copied += 1
    return copied


class Runner:
    def __init__(self, args: argparse.Namespace, cohort: bool) -> None:
        self.args, self.cohort = args, cohort
        self.logs = ROOT / 'logs'
        self.logs.mkdir(parents=True, exist_ok=True)
        name = 'cohort_pipeline_report.json' if cohort else 'full_pipeline_report.json'
        self.report_path = ROOT / 'data/pubmed/runs' / name
        self.report = {'project': ROOT.name, 'mode': 'cohort' if cohort else 'full_query_sweep',
                       'created_at': datetime.now(timezone.utc).isoformat(), 'stages': []}

    def flush(self) -> None:
        self.report_path.parent.mkdir(parents=True, exist_ok=True)
        self.report['updated_at'] = datetime.now(timezone.utc).isoformat()
        self.report_path.write_text(json.dumps(self.report, indent=2) + '\n', encoding='utf-8')

    def call(self, stage: str, argv: list[str], sentinel: Path | None = None) -> None:
        if sentinel and sentinel.is_file() and n_jsonl(sentinel) and not self.args.force:
            self.report['stages'].append({'stage_id': stage, 'status': 'skipped', 'sentinel_n': n_jsonl(sentinel)})
            self.flush(); print(f'[skip] {stage}'); return
        log = self.logs / f'{stage}.log'
        print(f'[run] {stage}', flush=True)
        try:
            run_script(argv, log, cohort=self.cohort)
        except Exception:
            self.report['stages'].append({'stage_id': stage, 'status': 'error', 'log': str(log)})
            self.flush()
            raise
        self.report['stages'].append({'stage_id': stage, 'status': 'ok', 'log': str(log),
                                      'sentinel_n': n_jsonl(sentinel) if sentinel else None})
        self.flush(); print(f'[ok] {stage}')

    def stage0(self) -> None:
        required = {'queries': ROOT/'configs/queries.json', 'journal': ROOT/'configs/journal_priority.json',
                    'gates': ROOT/'configs/pipeline_gates.json', 'clinvar': CLINVAR_SUMMARY,
                    'gff3': ENSEMBL_GFF3, 'fasta': ENSEMBL_FASTA, 'fai': ENSEMBL_FAI}
        missing = [f'{k}: {p}' for k,p in required.items() if not p.exists()]
        self.report['stages'].append({'stage_id':'stage0_preflight','status':'error' if missing else 'ok','missing':missing})
        self.flush()
        if missing: raise FileNotFoundError('; '.join(missing))
        print('[ok] stage0_preflight')

    def stage1(self) -> None:
        if self.cohort:
            ingest_cohort(self.args.pmid_list, self.args.seed_from, MERGED_ARTICLES, self.args.fetch_missing)
            self.report['stages'].append({'stage_id':'stage1_pubmed_retrieval','status':'ok','sentinel_n':117})
            if self.args.seed_cache:
                copied = seed_cache(self.args.seed_cache_from, set(load_pmids(self.args.pmid_list)))
                self.report['stages'].append({'stage_id':'stage1_optional_seed_cache','status':'ok','copied':copied})
            self.flush(); return
        self.call('stage1_pubmed_retrieval',
                  [str(SCRIPTS/'pubmed_pipeline.py'), '--config', str(ROOT/'configs/queries.json'),
                   '--mode','all','--queries',*QUERIES], MERGED_ARTICLES)

    def stage2(self): self.call('stage2_abstract_linking', [str(SCRIPTS/'link_variant_phenotype.py')], S2/'articles_pass.jsonl')
    def stage3(self):
        argv=[str(SCRIPTS/'fulltext.py'),'enrich','--scope','pass','--retry-missing']
        if self.cohort and any((ROOT/'data/pubmed/cache/fulltext/parsed').glob('*.json')):
            argv += ['--rebuild-from-parsed','--reprocess-links']
        self.call('stage3_fulltext_enrichment', argv, S2/'corpus_enriched.jsonl')
        stamp_evidence_source()
    def stage4(self): self.call('stage4_discovery_tagging', [str(SCRIPTS/'review.py'),'discovery'], S4/'corpus_enriched_novel.jsonl')
    def stage5(self): self.call('stage5_clinvar_validation',
        [str(SCRIPTS/'validate_novelty_clinvar.py'),'--track','pubmed','--corpus','all','--skip-download','--clinvar-file',str(CLINVAR_SUMMARY)],
        S5_ENRICHED/'corpus_enriched_novel_clinvar_absent.jsonl')
    def stage6(self): self.call('stage6_review_quality', [str(SCRIPTS/'review.py'),'filter'], REVIEW_HIGHCONF)
    def stage7(self):
        self.call('stage7_journal_strict',
                  [str(SCRIPTS/'review.py'),'journal','--input',str(REVIEW_STRICT),'--out-dir',str(S6)],
                  S6/'corpus_enriched_novel_clinvar_absent_review_strict_journal_annotated.jsonl')
        self.call('stage7_journal_highconf',
                  [str(SCRIPTS/'review.py'),'journal','--input',str(REVIEW_HIGHCONF),'--out-dir',str(S6)], REVIEW_HIGHCONF_ANN)
    def stage8_pdf(self):
        if os.environ.get('VARIANT_CLINICAL_STAGE8_PDF'):
            self.call('stage8_download_pdfs', [str(SCRIPTS/'extras.py'),'download-pdfs','--jsonl',str(REVIEW_HIGHCONF)])
        argv=[str(SCRIPTS/'fulltext.py'),'materialize','--jsonl',str(REVIEW_HIGHCONF)] + (['--force'] if self.args.force else [])
        self.call('stage8_materialize_text', argv)
        self.call('stage8_pdf_pathogenicity',
                  [str(SCRIPTS/'review.py'),'pathogenicity','--input',str(REVIEW_HIGHCONF),'--text-cache-dir',str(PDF_TEXT_CACHE)])
    def stage8_neg(self): self.call('stage8_negatives',
        [str(SCRIPTS/'review.py'),'negatives','--pmid-source','full_pass','--exclude-positives',str(REVIEW_HIGHCONF_ANN)], NEGATIVES_FULL_PASS)
    def stage9(self): self.call('stage9_genomic_export',
        [str(SCRIPTS/'genomic.py'),'export','--pos',str(REVIEW_HIGHCONF_ANN),'--pos-patho',str(PDF_PATHO)], S9/'variants_genomic_all.jsonl')


STAGE_METHODS = {'stage0_preflight':'stage0','stage1_pubmed_retrieval':'stage1','stage2_abstract_linking':'stage2',
 'stage3_fulltext_enrichment':'stage3','stage4_discovery_tagging':'stage4','stage5_clinvar_validation':'stage5',
 'stage6_review_quality':'stage6','stage7_journal_priority':'stage7','stage8_pdf_evidence':'stage8_pdf',
 'stage8_negatives':'stage8_neg','stage9_genomic_export':'stage9'}


def run_pipeline(args: argparse.Namespace, cohort: bool) -> None:
    if cohort and args.seed_cache and not args.seed_cache_from:
        raise SystemExit('--seed-cache requires --seed-cache-from')
    start, end = STAGE_ORDER.index(args.from_stage), STAGE_ORDER.index(args.to_stage)
    if end < start: raise SystemExit('--to must be at or after --from')
    runner=Runner(args,cohort)
    for stage in STAGE_ORDER[start:end+1]: getattr(runner,STAGE_METHODS[stage])()
    runner.report['status']='completed'; runner.flush()


# --- Reference bootstrap (formerly tools/bootstrap_refs.py) ---
def sha256(path: Path) -> str:
    digest=hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda:handle.read(1<<20),b''): digest.update(chunk)
    return digest.hexdigest()


def transfer_reference(src: Path, dst: Path, link: bool) -> None:
    if not src.is_file(): raise FileNotFoundError(src)
    dst.parent.mkdir(parents=True,exist_ok=True)
    if dst.exists(): dst.unlink()
    os.link(src,dst) if link else shutil.copy2(src,dst)


def cmd_bootstrap(args: argparse.Namespace) -> None:
    if args.import_from and args.link_from: raise SystemExit('choose only one source option')
    ref=ROOT/'data/pubmed/ref'; ref.mkdir(parents=True,exist_ok=True)
    source=args.import_from or args.link_from
    if source:
        src=source/'data/pubmed/ref'
        names=['clinvar/variant_summary.txt.gz']
        if not args.skip_ensembl: names += ['ensembl/Homo_sapiens.GRCh38.115.gff3','ensembl/Homo_sapiens.GRCh38.dna.toplevel.fa','ensembl/Homo_sapiens.GRCh38.dna.toplevel.fa.fai']
        for name in names: transfer_reference(src/name,ref/name,bool(args.link_from))
    elif args.download_clinvar:
        import urllib.request
        dest=ref/'clinvar/variant_summary.txt.gz'; dest.parent.mkdir(parents=True,exist_ok=True)
        urllib.request.urlretrieve('https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz',dest)
    clinvar=ref/'clinvar/variant_summary.txt.gz'
    if not clinvar.is_file(): raise FileNotFoundError(clinvar)
    manifest={'project':ROOT.name,'files':{'clinvar_variant_summary':{'path':'clinvar/variant_summary.txt.gz','bytes':clinvar.stat().st_size,'sha256':sha256(clinvar)}}}
    (ref/'reference_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


# --- Gold genomic comparison (formerly tests/compare_gold_genomic.py) ---
def cmd_compare(args: argparse.Namespace) -> None:
    def load(path):
        with path.open(encoding='utf-8-sig') as handle: return list(csv.DictReader(handle))
    norm=lambda value:(value or '').strip()
    locus=lambda row:(norm(row.get('pmid')),norm(row.get('gene')).upper(),norm(row.get('variant_text') or row.get('variant')).lower().replace(' ',''))
    genomic=lambda row:(norm(row.get('chromosome_ucsc') or row.get('chromosome')),norm(row.get('position')),norm(row.get('ref')).upper(),norm(row.get('alt')).upper())
    gold=load(args.gold); complete=[r for r in gold if all(norm(r.get(c)) for c in ('chromosome_ucsc','position','ref','alt'))]
    if not args.candidate.is_file(): raise FileNotFoundError(args.candidate)
    candidate=load(args.candidate); by_locus={locus(r):r for r in candidate}
    exact=sum(locus(g) in by_locus and genomic(g)==genomic(by_locus[locus(g)]) for g in complete)
    report={'gold_rows':len(gold),'gold_with_4col':len(complete),'candidate_rows':len(candidate),'genomic_exact_match':exact,'genomic_match_rate':round(exact/len(complete),4) if complete else None}
    args.out.parent.mkdir(parents=True,exist_ok=True); args.out.write_text(json.dumps(report,indent=2)+'\n'); print(json.dumps(report,indent=2))
    if exact != len(complete): raise SystemExit(1)


def cmd_selftest(_: argparse.Namespace) -> None:
    import genomic as genomic_module
    import review as review_module
    sentence='While studying blood phenotype-associated variation, we identified a genetic variant, rs17437411, associated with reduced blood counts.'
    result=review_module.classify_discovery(sentence,'rs17437411')
    assert result['discovery_label']=='novel' and result['discovery_confidence']=='high'
    for token in ('WT','ROS','HPP','PM2','VUS','ACMG'):
        assert token in review_module.BAD_GENE_TOKENS and review_module._acceptable_gene(token) is None
    assert review_module.guess_gene_near('c.206G>T (p.Arg69Leu) in EDA','c.206G>T',None)=='EDA'
    blob='KDM4C(NM_015061.6):c.2556C>T p.(Asp852=) 9:7103816 (GRCh38) rs3763651'
    comps=genomic_module.merge_comps_from_evidence(genomic_module.parse_variant_components('rs3763651'),variant_text='rs3763651',evidence_blob=blob)
    assert comps.get('c_hgvs_input')=='c.2556C>T' and genomic_module.try_literal_genomic(blob,comps)['resolve_status']=='partial'
    mini='##gff-version 3\n1\tx\tgene\t100\t300\t.\t+\t.\tID=gene:G;Name=GENE1;biotype=protein_coding\n1\tx\tmRNA\t100\t300\t.\t+\t.\tID=transcript:T;Parent=gene:G\n1\tx\texon\t100\t200\t.\t+\t.\tParent=transcript:T\n1\tx\tCDS\t120\t180\t.\t+\t0\tParent=transcript:T\n'
    with tempfile.TemporaryDirectory() as td:
        path=Path(td)/'mini.gff3'; path.write_text(mini)
        index=genomic_module.build_gff3_region_index(path)
        assert index.annotate('1',150)['genomic_region']=='coding'
    from vc_text import gene_from_naming_cues, resolve_gene_for_variant
    named='We found a long non-coding RNA, which we named HOTSCRAMBL.'
    assert gene_from_naming_cues(named)=='HOTSCRAMBL'
    pub={'variants':[{'text':'rs1','gene_id':3204,'identifier':'CorrespondingGene:3204'}],'genes':[{'name':'HOXA7','identifier':'3204'}]}
    assert resolve_gene_for_variant('rs1',evidence_sentence='HOXA7 rs1 was found',pubtator=pub)=='HOXA7'
    print('ok: four integrated regression groups passed')


def add_pipeline_args(parser: argparse.ArgumentParser, cohort: bool) -> None:
    parser.add_argument('--from',dest='from_stage',choices=STAGE_ORDER,default=STAGE_ORDER[0])
    parser.add_argument('--to',dest='to_stage',choices=STAGE_ORDER,default=STAGE_ORDER[-1])
    parser.add_argument('--force',action='store_true')
    if cohort:
        parser.add_argument('--pmid-list',type=Path,default=ROOT/'configs/gold_pmids.txt')
        parser.add_argument('--seed-from',type=Path,default=ROOT/'fixtures/cohort/articles.jsonl')
        parser.add_argument('--fetch-missing',action='store_true')
        parser.add_argument('--seed-cache',action='store_true')
        parser.add_argument('--seed-cache-from',type=Path)


def main() -> None:
    ensure_data_link()
    parser=argparse.ArgumentParser(description=__doc__); sub=parser.add_subparsers(dest='command',required=True)
    cohort=sub.add_parser('cohort'); add_pipeline_args(cohort,True)
    pipeline=sub.add_parser('pipeline'); add_pipeline_args(pipeline,False)
    single=sub.add_parser('single'); single.add_argument('pmid'); single.add_argument('--from',dest='from_stage',default='stage0_preflight'); single.add_argument('--to',dest='to_stage',default='stage9_genomic_export_qc'); single.add_argument('--force',action='store_true'); single.add_argument('--no-network',action='store_true')
    deliver=sub.add_parser('deliver'); deliver.add_argument('step',nargs='?',default='all',choices=['all','new-site','audit','html','review-export'])
    qc=sub.add_parser('qc'); qsub=qc.add_subparsers(dest='qc_command',required=True)
    refalt=qsub.add_parser('refalt'); refalt.add_argument('--in-csv',required=True); refalt.add_argument('--out-tsv',required=True)
    refctx=qsub.add_parser('ref-context5'); refctx.add_argument('--in-csv',required=True); refctx.add_argument('--out-csv')
    bootstrap=sub.add_parser('bootstrap'); bootstrap.add_argument('--import-from',type=Path); bootstrap.add_argument('--link-from',type=Path); bootstrap.add_argument('--skip-ensembl',action='store_true'); bootstrap.add_argument('--download-clinvar',action='store_true')
    compare=sub.add_parser('compare'); compare.add_argument('--gold',type=Path,default=ROOT/'configs/gold_standard_v3.csv'); compare.add_argument('--candidate',type=Path,default=S9/'variants_genomic_positive_ok.csv'); compare.add_argument('--out',type=Path,default=ROOT/'data/pubmed/runs/cohort_gold/repro_compare.json')
    selftest=sub.add_parser('selftest',help=argparse.SUPPRESS)
    stamp=sub.add_parser('stamp',help=argparse.SUPPRESS); stamp.add_argument('--linked-dir',type=Path,default=S2)
    args=parser.parse_args()
    if args.command=='cohort': run_pipeline(args,True)
    elif args.command=='pipeline': run_pipeline(args,False)
    elif args.command=='single':
        argv=[str(SCRIPTS/'extras.py'),'single',args.pmid,'--from',args.from_stage,'--to',args.to_stage]+(['--force'] if args.force else [])+(['--no-network'] if args.no_network else [])
        run_script(argv, ROOT/'logs'/f'single_{args.pmid}.log', cohort=False)
    elif args.command=='deliver':
        steps=['new-site','audit','html','review-export'] if args.step=='all' else [args.step]
        extra={'audit':['--class','positive']}
        for step in steps:
            run_script([str(SCRIPTS/'extras.py'),step,*extra.get(step,[])], ROOT/'logs'/f'deliver_{step}.log', cohort=False)
    elif args.command=='qc':
        argv=[str(SCRIPTS/'extras.py'),args.qc_command,'--in-csv',args.in_csv]
        argv += ['--out-tsv',args.out_tsv] if args.qc_command=='refalt' else (['--out-csv',args.out_csv] if args.out_csv else [])
        run_script(argv, ROOT/'logs'/f'qc_{args.qc_command}.log', cohort=False)
    elif args.command=='bootstrap': cmd_bootstrap(args)
    elif args.command=='compare': cmd_compare(args)
    elif args.command=='selftest': cmd_selftest(args)
    elif args.command=='stamp': stamp_evidence_source(args.linked_dir)

if __name__=='__main__': main()
