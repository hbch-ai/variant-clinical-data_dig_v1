#!/usr/bin/env python3
"""Physically consolidated pipeline commands. Generated from the v4_01 stage sources.

Each section retains its source module boundary through descriptive global prefixes;
there is no dynamic source loading or embedded executable source text.
"""
from __future__ import annotations
import vc_paths as vc_paths_module
import vc_text as vc_text_module

# === enrich_fulltext.py ===
"""
Full-text enrichment for variant ↔ clinical manifestation links.

Fetch order for each PMCID (any freely available PMC full text, not only
BioC Open-Access subset):
  1. NCBI BioC JSON (PMC OA subset)
  2. Europe PMC JATS XML
  3. PMC HTML (OA + Author Manuscript / free PMC viewer)
  4. PMC PDF (when HTML missing or too thin)

Outputs under data/pubmed/cache/fulltext/:
  raw/{pmcid}.bioc.json          cached BioC payloads
  raw/{pmcid}.jats.xml           EuropePMC JATS
  raw/{pmcid}.pmc.html           PMC HTML viewer page
  raw/{pmcid}.pmc.pdf            PMC PDF (optional)
  parsed/{pmid}.json             passage inventory + local spans summary
  linked/pairs_fulltext.jsonl    new/enriched FT links
  linked/pairs_enriched.jsonl    abstract pairs + FT enrichments merged
  linked/corpus_fulltext.jsonl   FT-derived corpus rows (text first)
  linked/corpus_enriched.jsonl   merged abstract-strict + FT corpus
  linked/summary.json
"""
import argparse
import json
import logging
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from html import unescape
from pathlib import Path
from typing import Any
from urllib.parse import quote, urljoin
import requests
enrich_fulltext_ROOT = Path(__file__).resolve().parents[1]
enrich_fulltext_SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(enrich_fulltext_SCRIPTS))
enrich_fulltext_MERGED = enrich_fulltext_ROOT / 'data' / 'pubmed' / 'ingest' / 'merged' / 'articles.jsonl'
enrich_fulltext_LINKED_DIR = enrich_fulltext_ROOT / 'data' / 'pubmed' / 'pipeline' / 's2_linking'
enrich_fulltext_FT_ROOT = enrich_fulltext_ROOT / 'data' / 'pubmed' / 'cache' / 'fulltext'
enrich_fulltext_RAW_DIR = enrich_fulltext_FT_ROOT / 'raw'
enrich_fulltext_PARSED_DIR = enrich_fulltext_FT_ROOT / 'parsed'
enrich_fulltext_FT_LINKED = enrich_fulltext_FT_ROOT / 'linked'
enrich_fulltext_LOG_DIR = enrich_fulltext_ROOT / 'logs'
enrich_fulltext_BIOC_URL = 'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/{pmcid}/unicode'
enrich_fulltext_EUROPEPMC_XML = 'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML'
enrich_fulltext_PMC_ARTICLE_URL = 'https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/'
enrich_fulltext_PMC_OA_STATUS_URL = 'https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id={pmcid}'
enrich_fulltext_SESSION = requests.Session()
enrich_fulltext_SESSION.headers.update({'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/122.0.0.0 Safari/537.36', 'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8', 'Accept-Language': 'en-US,en;q=0.9'})
enrich_fulltext_MIN_PMC_HTML_CHARS = 2500
enrich_fulltext_MIN_PASSAGE_CHARS = 40

def enrich_fulltext_curl_get(url: str, dest: Path | None=None, timeout: int=90, binary: bool=False) -> tuple[int, bytes | str | None]:
    """
    Fetch via system curl. Returns (http_code, body_or_None).

    NCBI PMC often returns 403 to Python requests (TLS fingerprinting) while
    curl succeeds for freely available pages; use curl as primary transport.
    """
    cmd = ['curl', '-sS', '-L', '--max-time', str(timeout), '-A', 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/122.0.0.0 Safari/537.36', '-H', 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8', '-H', 'Accept-Language: en-US,en;q=0.9', '-w', '%{http_code}']
    import tempfile
    body_path: Path
    tmp_owned = False
    if dest is not None:
        body_path = dest
    else:
        tmp = tempfile.NamedTemporaryFile(delete=False)
        tmp.close()
        body_path = Path(tmp.name)
        tmp_owned = True
    cmd.extend(['-o', str(body_path), url])
    try:
        proc = subprocess.run(cmd, check=False, capture_output=True, timeout=timeout + 15)
        code_s = (proc.stdout or b'').decode('ascii', errors='ignore').strip()
        try:
            code = int(code_s[-3:]) if code_s else 0
        except ValueError:
            code = 0
        if proc.returncode not in (0, 22) and code == 0:
            return (0, None)
        if not body_path.exists() or body_path.stat().st_size < 50:
            return (code or 0, None)
        data = body_path.read_bytes()
        if tmp_owned:
            body_path.unlink(missing_ok=True)
        if binary:
            return (code, data)
        return (code, data.decode('utf-8', errors='replace'))
    except Exception:
        if tmp_owned and body_path.exists():
            body_path.unlink(missing_ok=True)
        return (0, None)

def enrich_fulltext_setup_logger(run_id: str) -> logging.Logger:
    enrich_fulltext_LOG_DIR.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger('enrich_fulltext')
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    fh = logging.FileHandler(enrich_fulltext_LOG_DIR / f'fulltext_{run_id}.log', encoding='utf-8')
    fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger

def enrich_fulltext_normalize_pmcid(pmcid: str | None) -> str | None:
    if not pmcid:
        return None
    pmcid = str(pmcid).strip()
    if not pmcid:
        return None
    if pmcid.upper().startswith('PMC'):
        return 'PMC' + re.sub('\\D', '', pmcid)
    if pmcid.isdigit():
        return f'PMC{pmcid}'
    return pmcid

def enrich_fulltext_fetch_bioc(pmcid: str, logger: logging.Logger, retries: int=3) -> dict[str, Any] | None:
    url = enrich_fulltext_BIOC_URL.format(pmcid=quote(pmcid))
    last_err: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            time.sleep(0.34)
            resp = enrich_fulltext_SESSION.get(url, timeout=120)
            if resp.status_code == 404:
                return None
            ctype = (resp.headers.get('content-type') or '').lower()
            body = resp.text.lstrip()
            if resp.status_code != 200:
                if resp.status_code in (429, 500, 502, 503, 504):
                    time.sleep(min(2 ** attempt, 20))
                    continue
                logger.info('BioC %s HTTP %s', pmcid, resp.status_code)
                return None
            if 'json' not in ctype and (not body.startswith('{')) and (not body.startswith('[')):
                logger.info('BioC %s non-JSON payload (likely non-OA); skip', pmcid)
                return None
            data = resp.json()
            if isinstance(data, list):
                return data[0] if data else None
            return data
        except Exception as exc:
            last_err = exc
            if 'Expecting value' in str(exc) or 'JSON' in type(exc).__name__:
                logger.info('BioC %s JSON parse fail (likely non-OA): %s', pmcid, exc)
                return None
            logger.warning('BioC %s failed (%s); retry %s', pmcid, exc, attempt)
            time.sleep(min(2 ** attempt, 20))
    logger.error('BioC permanently failed for %s: %s', pmcid, last_err)
    return None

def enrich_fulltext_fetch_europepmc_xml(pmcid: str, logger: logging.Logger) -> str | None:
    url = enrich_fulltext_EUROPEPMC_XML.format(pmcid=quote(pmcid))
    try:
        time.sleep(0.2)
        resp = enrich_fulltext_SESSION.get(url, timeout=15)
        if resp.status_code != 200:
            logger.info('EuropePMC miss %s HTTP %s', pmcid, resp.status_code)
            return None
        if b'<?xml' not in resp.content[:200] and b'<article' not in resp.content[:500]:
            return None
        return resp.text
    except requests.exceptions.SSLError as exc:
        logger.info('EuropePMC %s SSL error; skip (%s)', pmcid, exc)
        return None
    except requests.exceptions.RequestException as exc:
        logger.info('EuropePMC %s request error; skip (%s)', pmcid, exc)
        return None
    except Exception as exc:
        logger.info('EuropePMC %s failed; skip (%s)', pmcid, exc)
        return None

def enrich_fulltext_check_pmc_oa_subset(pmcid: str, logger: logging.Logger) -> bool | None:
    """True if in PMC OA subset; False if explicitly not; None if unknown."""
    url = enrich_fulltext_PMC_OA_STATUS_URL.format(pmcid=quote(pmcid))
    try:
        time.sleep(0.15)
        resp = enrich_fulltext_SESSION.get(url, timeout=30)
        if resp.status_code != 200:
            return None
        body = resp.text
        if 'idIsNotOpenAccess' in body:
            return False
        if '<link' in body and 'href=' in body:
            return True
        if 'error' in body.lower():
            return False
        return None
    except Exception as exc:
        logger.info('OA status %s failed: %s', pmcid, exc)
        return None

def enrich_fulltext_fetch_pmc_html(pmcid: str, logger: logging.Logger, retries: int=2) -> str | None:
    """Fetch PMC article viewer HTML (OA or Author Manuscript)."""
    url = enrich_fulltext_PMC_ARTICLE_URL.format(pmcid=pmcid)
    for attempt in range(1, retries + 1):
        time.sleep(0.25)
        code, text = enrich_fulltext_curl_get(url, timeout=90, binary=False)
        if code in (429, 500, 502, 503, 504):
            logger.info('PMC HTML %s HTTP %s; retry %s', pmcid, code, attempt)
            time.sleep(min(2 ** attempt, 20))
            continue
        if code == 403:
            logger.info('PMC HTML %s HTTP 403 (not freely readable); skip', pmcid)
            return None
        if code != 200 or not isinstance(text, str) or (not text):
            logger.info('PMC HTML %s HTTP %s; skip', pmcid, code)
            return None
        low = text.lower()
        if '<article' not in low and 'pmc_sec_title' not in low:
            if 'login' in low or 'subscription' in low:
                logger.info('PMC HTML %s looks gated; skip', pmcid)
            else:
                logger.info('PMC HTML %s missing article body; skip', pmcid)
            return None
        return text
    return None

def enrich_fulltext__strip_html_noise(html: str) -> str:
    html = re.sub('(?is)<script[^>]*>.*?</script>', ' ', html)
    html = re.sub('(?is)<style[^>]*>.*?</style>', ' ', html)
    html = re.sub('(?is)<noscript[^>]*>.*?</noscript>', ' ', html)
    return html

def enrich_fulltext__html_to_text(fragment: str) -> str:
    text = re.sub('(?is)<br\\s*/?>', '\n', fragment)
    text = re.sub('(?is)</p>', '\n', text)
    text = re.sub('(?is)</h[1-6]>', '\n', text)
    text = re.sub('(?is)<[^>]+>', ' ', text)
    text = unescape(text)
    text = text.replace('\xa0', ' ')
    text = re.sub('[ \\t]+', ' ', text)
    text = re.sub('\\n{2,}', '\n', text)
    return text.strip()

def enrich_fulltext__guess_section_type(title: str) -> str:
    t = (title or '').lower()
    if not t:
        return 'OTHER'
    if 'abstract' in t or 'summary' in t:
        return 'ABSTRACT'
    if 'intro' in t:
        return 'INTRO'
    if 'case' in t or 'patient' in t:
        return 'CASE'
    if 'result' in t:
        return 'RESULTS'
    if 'discuss' in t:
        return 'DISCUSS'
    if 'method' in t or 'material' in t:
        return 'METHODS'
    if 'concl' in t:
        return 'CONCL'
    if 'reference' in t or 'bibliograph' in t:
        return 'REF'
    if 'figure' in t or t.startswith('fig'):
        return 'FIG'
    if 'table' in t:
        return 'TABLE'
    if 'supplement' in t or 'appendix' in t:
        return 'SUPPL'
    if 'acknowledge' in t or 'funding' in t or 'conflict' in t or ('disclosure' in t):
        return 'OTHER'
    return 'OTHER'

def enrich_fulltext_passages_from_pmc_html(html: str) -> list[dict[str, Any]]:
    """
    Parse PMC viewer HTML into BioC-like passages.

    Handles both commercial OA layouts and NIH Author Manuscript pages.
    """
    if not html:
        return []
    cleaned = enrich_fulltext__strip_html_noise(html)
    m = re.search('(?is)<article\\b[^>]*>(.*?)</article>', cleaned)
    body_html = m.group(1) if m else cleaned
    body_html = re.sub('(?is)<section[^>]*id=\\"[^\\"\']*ref[^\\"\']*\\"[^>]*>.*?</section>', ' ', body_html)
    body_html = re.sub('(?is)<div[^>]*class="[^"]*ref-list[^"]*"[^>]*>.*?</div>', ' ', body_html)
    out: list[dict[str, Any]] = []
    idx = 0
    title_m = re.search('(?is)<h1[^>]*class="[^"]*title[^"]*"[^>]*>(.*?)</h1>|<title[^>]*>(.*?)</title>', cleaned)
    if title_m:
        title = enrich_fulltext__html_to_text(title_m.group(1) or title_m.group(2) or '')
        title = re.sub('\\s*[-|].*PMC.*$', '', title).strip()
        if title and len(title) > 5:
            out.append({'idx': idx, 'section_type': 'TITLE', 'passage_type': 'title', 'text': title, 'offset': None})
            idx += 1
    section_pat = re.compile('(?is)<section\\b[^>]*>\\s*(?:<h[1-4][^>]*class="[^"]*pmc_sec_title[^"]*"[^>]*>(.*?)</h[1-4]>)?(.*?)</section>')
    sections = list(section_pat.finditer(body_html))
    used_sectioned = False
    if sections:
        for sm in sections:
            heading = enrich_fulltext__html_to_text(sm.group(1) or '')
            sec_type = enrich_fulltext__guess_section_type(heading)
            if sec_type in vc_text_module.SKIP_SECTION_TYPES or sec_type == 'REF':
                continue
            inner = sm.group(2) or ''
            paras = re.findall('(?is)<p\\b[^>]*>(.*?)</p>', inner)
            if not paras:
                blob = enrich_fulltext__html_to_text(inner)
                if len(blob) >= enrich_fulltext_MIN_PASSAGE_CHARS:
                    paras = [blob]
            else:
                paras = [enrich_fulltext__html_to_text(p) for p in paras]
            for para in paras:
                para = re.sub('\\s+', ' ', para).strip()
                if len(para) < enrich_fulltext_MIN_PASSAGE_CHARS:
                    continue
                if para.lower().startswith('open in a new tab'):
                    continue
                out.append({'idx': idx, 'section_type': sec_type, 'passage_type': 'paragraph', 'text': para, 'offset': None, 'section_title': heading or None})
                idx += 1
                used_sectioned = True
    if used_sectioned and len(out) >= 3:
        return out
    current_sec = 'OTHER'
    for block in re.finditer('(?is)<(h[1-4]|p)\\b([^>]*)>(.*?)</\\1>', body_html):
        tag = block.group(1).lower()
        attrs = block.group(2) or ''
        inner = block.group(3) or ''
        text = enrich_fulltext__html_to_text(inner)
        text = re.sub('\\s+', ' ', text).strip()
        if tag.startswith('h'):
            if 'pmc_sec_title' in attrs or len(text) < 120:
                current_sec = enrich_fulltext__guess_section_type(text)
            continue
        if current_sec in vc_text_module.SKIP_SECTION_TYPES or current_sec == 'REF':
            continue
        if len(text) < enrich_fulltext_MIN_PASSAGE_CHARS:
            continue
        if text.lower().startswith('open in a new tab'):
            continue
        out.append({'idx': idx, 'section_type': current_sec, 'passage_type': 'paragraph', 'text': text, 'offset': None})
        idx += 1
    return out

def enrich_fulltext_extract_pmc_pdf_urls(html: str, pmcid: str) -> list[str]:
    """Collect candidate PDF URLs from PMC HTML."""
    urls: list[str] = []
    base = enrich_fulltext_PMC_ARTICLE_URL.format(pmcid=pmcid)
    for m in re.finditer('(?is)href="([^"]+\\.pdf[^"]*)"', html):
        href = unescape(m.group(1).strip())
        if href.startswith('//'):
            href = 'https:' + href
        elif href.startswith('/'):
            href = urljoin('https://www.ncbi.nlm.nih.gov', href)
        elif not href.startswith('http'):
            href = urljoin(base, href)
        urls.append(href)
    urls.extend([f'https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/pdf/', f'https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/pdf/main.pdf'])
    seen: set[str] = set()
    out: list[str] = []
    for u in urls:
        if u in seen:
            continue
        seen.add(u)
        out.append(u)
    return out

def enrich_fulltext_fetch_pmc_pdf(pmcid: str, html: str | None, dest: Path, logger: logging.Logger) -> Path | None:
    """Download PMC PDF if available; return local path."""
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    candidates = enrich_fulltext_extract_pmc_pdf_urls(html or '', pmcid)
    for url in candidates:
        try:
            time.sleep(0.25)
            tmp = dest.with_suffix('.pdf.part')
            code, blob = enrich_fulltext_curl_get(url, dest=tmp, timeout=120, binary=True)
            if code == 403:
                logger.info('PMC PDF %s HTTP 403 for %s', pmcid, url)
                if tmp.exists():
                    tmp.unlink(missing_ok=True)
                continue
            if code != 200 or not isinstance(blob, (bytes, bytearray)):
                if tmp.exists():
                    content = tmp.read_bytes()
                else:
                    continue
            else:
                content = bytes(blob) if not tmp.exists() else tmp.read_bytes()
            if not content.startswith(b'%PDF') or len(content) < 1000:
                if tmp.exists():
                    tmp.unlink(missing_ok=True)
                continue
            if tmp.exists():
                tmp.replace(dest)
            else:
                dest.write_bytes(content)
            logger.info('PMC PDF %s saved (%s bytes) from %s', pmcid, dest.stat().st_size, url)
            return dest
        except Exception as exc:
            logger.info('PMC PDF %s candidate fail %s: %s', pmcid, url, exc)
            continue
    logger.info('PMC PDF %s unavailable', pmcid)
    return None

def enrich_fulltext_pdf_to_text(pdf_path: Path, logger: logging.Logger) -> str:
    """Extract text via pdftotext, fallback to pypdf."""
    try:
        proc = subprocess.run(['pdftotext', '-layout', '-enc', 'UTF-8', str(pdf_path), '-'], check=False, capture_output=True, text=True, timeout=120)
        if proc.returncode == 0 and (proc.stdout or '').strip():
            return proc.stdout
    except Exception as exc:
        logger.info('pdftotext failed for %s: %s', pdf_path.name, exc)
    try:
        from pypdf import PdfReader
        reader = PdfReader(str(pdf_path))
        parts: list[str] = []
        for page in reader.pages:
            parts.append(page.extract_text() or '')
        return '\n'.join(parts)
    except Exception as exc:
        logger.info('pypdf failed for %s: %s', pdf_path.name, exc)
        return ''

def enrich_fulltext_passages_from_plain_text(text: str, source_label: str='pdf') -> list[dict[str, Any]]:
    """Split plain full text into paragraph-like passages."""
    if not text:
        return []
    text = text.replace('\r\n', '\n').replace('\r', '\n')
    text = text.replace('\x0c', '\n\n')
    chunks = re.split('\\n{2,}', text)
    out: list[dict[str, Any]] = []
    idx = 0
    current_sec = 'OTHER'
    for chunk in chunks:
        lines = [re.sub('\\s+', ' ', ln).strip() for ln in chunk.split('\n')]
        lines = [ln for ln in lines if ln]
        if not lines:
            continue
        if len(lines[0]) < 100 and (lines[0].isupper() or re.match('^(Abstract|Introduction|Methods?|Results?|Discussion|Conclusions?|Case\\b|References|Acknowledgments?)\\b', lines[0], re.I)):
            current_sec = enrich_fulltext__guess_section_type(lines[0])
            body = ' '.join(lines[1:]) if len(lines) > 1 else ''
        else:
            body = ' '.join(lines)
        body = re.sub('\\s+', ' ', body).strip()
        if current_sec in vc_text_module.SKIP_SECTION_TYPES or current_sec == 'REF':
            continue
        if len(body) < enrich_fulltext_MIN_PASSAGE_CHARS:
            continue
        out.append({'idx': idx, 'section_type': current_sec, 'passage_type': source_label, 'text': body, 'offset': None})
        idx += 1
    return out

def enrich_fulltext_passages_from_bioc(coll: dict[str, Any]) -> list[dict[str, Any]]:
    docs = coll.get('documents') or []
    if not docs:
        return []
    out: list[dict[str, Any]] = []
    for i, p in enumerate(docs[0].get('passages') or []):
        info = p.get('infons') or {}
        section = info.get('section_type') or info.get('type') or 'OTHER'
        ptype = info.get('type') or ''
        text = re.sub('\\s+', ' ', (p.get('text') or '').strip())
        if not text:
            continue
        out.append({'idx': i, 'section_type': section, 'passage_type': ptype, 'text': text, 'offset': p.get('offset')})
    return out

def enrich_fulltext_passages_from_jats_xml(xml_text: str) -> list[dict[str, Any]]:
    """Lightweight JATS paragraph extraction fallback."""
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return []

    def local(tag: str) -> str:
        return tag.split('}')[-1] if '}' in tag else tag
    out: list[dict[str, Any]] = []
    idx = 0
    for t in root.iter():
        if local(t.tag) == 'article-title':
            text = re.sub('\\s+', ' ', ''.join(t.itertext()).strip())
            if text:
                out.append({'idx': idx, 'section_type': 'TITLE', 'passage_type': 'title', 'text': text, 'offset': None})
                idx += 1
            break
    for abs_el in root.iter():
        if local(abs_el.tag) == 'abstract':
            text = re.sub('\\s+', ' ', ''.join(abs_el.itertext()).strip())
            if text:
                out.append({'idx': idx, 'section_type': 'ABSTRACT', 'passage_type': 'abstract', 'text': text, 'offset': None})
                idx += 1
    body = None
    for el in root.iter():
        if local(el.tag) == 'body':
            body = el
            break
    if body is not None:
        current_sec = 'OTHER'
        for el in body.iter():
            name = local(el.tag)
            if name == 'title':
                current_sec = 'OTHER'
                title = re.sub('\\s+', ' ', ''.join(el.itertext()).strip()).lower()
                if 'case' in title:
                    current_sec = 'CASE'
                elif 'discuss' in title:
                    current_sec = 'DISCUSS'
                elif 'result' in title:
                    current_sec = 'RESULTS'
                elif 'intro' in title:
                    current_sec = 'INTRO'
                elif 'method' in title:
                    current_sec = 'METHODS'
                elif 'concl' in title:
                    current_sec = 'CONCL'
            if name == 'p':
                text = re.sub('\\s+', ' ', ''.join(el.itertext()).strip())
                if len(text) < 40:
                    continue
                out.append({'idx': idx, 'section_type': current_sec, 'passage_type': 'paragraph', 'text': text, 'offset': None})
                idx += 1
    return out

def enrich_fulltext_usable_passage(p: dict[str, Any]) -> bool:
    sec = (p.get('section_type') or '').upper()
    if sec in vc_text_module.SKIP_SECTION_TYPES or sec == 'REF':
        return False
    if p.get('passage_type') == 'ref':
        return False
    text = p.get('text') or ''
    if len(text) < 40:
        return False
    return True

def enrich_fulltext_link_passage(p: dict[str, Any]) -> list[dict[str, Any]]:
    sec = p.get('section_type')
    extracted = vc_text_module.extract_sentence_links(p['text'], source='fulltext_regex', method_prefix='fulltext_', section=sec, require_concrete_variant=True)
    links = extracted['same_sentence_links']
    if not links and extracted['flags']['has_strong_variant_text'] and extracted['flags']['has_clinical_text']:
        var_spans = [s for s in extracted['variant_spans'] if s['subtype'] in vc_text_module.CONCRETE_VARIANT_SUBTYPES]
        clin_spans = extracted['clinical_spans']
        pick = vc_text_module.pick_best_clinical_span(clin_spans, allow_cue=True, join_entities=True)
        if var_spans and pick:
            para = extracted['text']
            v0 = var_spans[0]
            if len(para) > 900:
                left = max(0, v0['start'] - 350)
                right = min(len(para), v0['end'] + 450)
                para = ('...' if left > 0 else '') + para[left:right] + ('...' if right < len(extracted['text']) else '')
            links.append({'variant': v0['text'], 'variant_subtype': v0['subtype'], 'clinical': pick['text'], 'clinical_subtype': pick['subtype'], 'evidence_sentence': para, 'method': 'fulltext_same_paragraph_cooccurrence', 'tier': 'B', 'text_scope': 'fulltext', 'section': sec, 'clinical_source': 'lexicon_entity' if pick['subtype'] in vc_text_module.ENTITY_CLINICAL_SUBTYPES else 'lexicon_cue'})
    return links

def enrich_fulltext_enrich_known_variants(passages: list[dict[str, Any]], known_variants: list[str], pmid: str) -> list[dict[str, Any]]:
    """For known HGVS/rsID from abstract linking, find richer full-text sentences."""
    out: list[dict[str, Any]] = []
    variants = [vc_text_module.normalize_variant_token(v) for v in known_variants if vc_text_module.is_concrete_variant_string(v)]
    if not variants:
        return out
    variants = sorted(set(variants), key=len, reverse=True)

    def _variant_in_sentence(var: str, sent: str) -> bool:
        if var in sent:
            return True
        ns = vc_text_module.normalize_variant_token(sent)
        if var and var in ns:
            return True
        m = re.match('^(c\\.\\S*?)([ACGT])>([ACGT])$', var, re.I)
        if m:
            pat = re.escape(m.group(1) + m.group(2)) + '\\s*>\\s*' + re.escape(m.group(3))
            return bool(re.search(pat, sent, re.I))
        return False
    for p in passages:
        if not enrich_fulltext_usable_passage(p):
            continue
        sec = p.get('section_type')
        text = p['text']
        for sent_tuple in vc_text_module.split_sentences(text):
            sent = sent_tuple[2]
            for var in variants:
                if not _variant_in_sentence(var, sent):
                    continue
                clin_hit = vc_text_module.extract_sentence_links(sent, method_prefix='fulltext_enrich_')
                pick = vc_text_module.pick_best_clinical_span(clin_hit['clinical_spans'], allow_cue=True, join_entities=True)
                if pick:
                    clin_text, clin_subtype = (pick['text'], pick['subtype'])
                    clin_source = 'lexicon_entity' if clin_subtype in vc_text_module.ENTITY_CLINICAL_SUBTYPES else 'lexicon_cue'
                elif len(sent) >= 60:
                    clin_text = 'fulltext_context'
                    clin_subtype = 'narrative_context'
                    clin_source = 'narrative_context'
                else:
                    continue
                out.append({'pmid': pmid, 'variant': var, 'clinical': clin_text, 'clinical_subtype': clin_subtype, 'clinical_source': clin_source, 'evidence_sentence': sent, 'method': 'fulltext_known_variant_sentence', 'tier': 'A' if clin_subtype in vc_text_module.ENTITY_CLINICAL_SUBTYPES else 'B', 'text_scope': 'fulltext', 'section': sec})
                break
    return out

def enrich_fulltext_uniq_keep_order(items: list[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for x in items:
        x = vc_text_module.normalize_variant_token(x)
        if not x or x in seen:
            continue
        seen.add(x)
        out.append(x)
    return out

def enrich_fulltext_load_pubtator_cache(path: Path) -> dict[str, dict[str, Any]]:
    by_pmid: dict[str, dict[str, Any]] = {}
    if not path.exists():
        return by_pmid
    with path.open(encoding='utf-8') as f:
        for line in f:
            if not line.strip():
                continue
            row = json.loads(line)
            pmid = str(row.get('pmid') or '')
            if pmid:
                by_pmid[pmid] = row
    return by_pmid

def enrich_fulltext_passage_for_sentence(passages: list[dict[str, Any]], sentence: str) -> str | None:
    if not sentence:
        return None
    for p in passages:
        text = p.get('text') or ''
        if sentence in text:
            return text
    return None

def enrich_fulltext_build_fulltext_links(discovery_passages: list[dict[str, Any]], abstract_variants: list[str], pubtator: dict[str, Any] | None, pmid: str, *, title: str | None=None, abstract: str | None=None) -> tuple[list[dict[str, Any]], dict[str, int]]:
    """
    Combined full-text linking:
      1) scan passages for concrete variants (scheme 2)
      2) link_passage + enrich_known_variants on merged seeds
      3) upgrade clinical fields via lexicon + PubTator (scheme 1 / PubTator)
      4) resolve gene per variant (PubTator / nearby mention; never genes[0])
    """
    ft_variants = vc_text_module.scan_passages_for_concrete_variants(discovery_passages)
    abs_set = {vc_text_module.normalize_variant_token(v) for v in abstract_variants if v}
    ft_only = {v for v in ft_variants if v not in abs_set}
    seed_variants = enrich_fulltext_uniq_keep_order(list(abstract_variants) + ft_variants)
    meta = {'n_ft_variants_scanned': len(ft_variants), 'n_ft_only_variants': len(ft_only), 'n_seed_variants': len(seed_variants)}
    discovered: list[dict[str, Any]] = []
    for p in discovery_passages:
        for link in enrich_fulltext_link_passage(p):
            link['pmid'] = pmid
            v = vc_text_module.normalize_variant_token(str(link.get('variant') or ''))
            link['variant'] = v
            link['variant_discovered_in_fulltext'] = v in ft_only
            discovered.append(link)
    enriched = enrich_fulltext_enrich_known_variants(discovery_passages, seed_variants, pmid)
    for link in enriched:
        v = vc_text_module.normalize_variant_token(str(link.get('variant') or ''))
        link['variant'] = v
        link['variant_discovered_in_fulltext'] = v in ft_only
    merged_links: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for link in discovered + enriched:
        key = (str(link.get('variant')), str(link.get('evidence_sentence'))[:240])
        if key in seen:
            continue
        seen.add(key)
        para = enrich_fulltext_passage_for_sentence(discovery_passages, str(link.get('evidence_sentence') or ''))
        up = vc_text_module.upgrade_clinical_field(link, paragraph_text=para, title=title, pubtator=pubtator)
        gene = vc_text_module.resolve_gene_for_variant(up.get('variant'), evidence_sentence=up.get('evidence_sentence'), title=title, abstract=abstract, context_text=para, pubtator=pubtator, link_gene=up.get('gene'))
        if gene:
            up['gene'] = gene
        merged_links.append(up)
    return (merged_links, meta)

def enrich_fulltext_load_passages_from_raw(pmcid: str, parsed: dict[str, Any], raw_dir: Path) -> tuple[list[dict[str, Any]], str | None]:
    """Reconstruct passages from cached raw full-text without network fetch."""
    source = parsed.get('source') or ''
    if 'pmc_html' in str(source):
        html_path = raw_dir / f'{pmcid}.pmc.html'
        if html_path.exists():
            html = html_path.read_text(encoding='utf-8', errors='replace')
            passages = enrich_fulltext_passages_from_pmc_html(html)
            if passages:
                return (passages, str(source))
    if 'bioc' in str(source):
        bioc_path = raw_dir / f'{pmcid}.bioc.json'
        if bioc_path.exists():
            coll = json.loads(bioc_path.read_text(encoding='utf-8'))
            passages = enrich_fulltext_passages_from_bioc(coll)
            if passages:
                return (passages, str(source))
    if 'europepmc' in str(source) or 'jats' in str(source):
        xml_path = raw_dir / f'{pmcid}.jats.xml'
        if xml_path.exists():
            passages = enrich_fulltext_passages_from_jats_xml(xml_path.read_text(encoding='utf-8'))
            if passages:
                return (passages, str(source))
    if 'pdf' in str(source):
        txt_path = raw_dir / f'{pmcid}.pmc.pdf.txt'
        if txt_path.exists():
            passages = enrich_fulltext_passages_from_plain_text(txt_path.read_text(encoding='utf-8'), source_label='pdf')
            if passages:
                return (passages, str(source))
    for loader, label in ((lambda: enrich_fulltext_passages_from_pmc_html((raw_dir / f'{pmcid}.pmc.html').read_text(encoding='utf-8', errors='replace')), 'pmc_html_cache'), (lambda: enrich_fulltext_passages_from_bioc(json.loads((raw_dir / f'{pmcid}.bioc.json').read_text(encoding='utf-8'))), 'bioc_cache')):
        try:
            passages = loader()
            if passages:
                return (passages, label)
        except Exception:
            continue
    return ([], None)

def enrich_fulltext_corpus_row(pmid: str, title: str | None, link: dict[str, Any]) -> dict[str, Any]:
    sent = link.get('evidence_sentence') or ''
    variant = link.get('variant')
    clinical = link.get('clinical')
    gene = link.get('gene')
    section = link.get('section')
    if sent and variant and (variant in sent) and (len(sent) < 1200):
        text = f"PMID:{pmid}. [fulltext:{section or 'body'}] {sent}"
    else:
        gene_bit = f' in {gene}' if gene else ''
        text = f"PMID:{pmid}. [fulltext:{section or 'body'}] The human genome variant {variant}{gene_bit} is described in full text together with clinical manifestation/phenotype context: {clinical}."
        if title:
            text += f' Source title: {title}'
    return {'text': text, 'pmid': pmid, 'tier': link.get('tier'), 'method': link.get('method'), 'variant': variant, 'clinical': clinical, 'gene': gene, 'section': section, 'evidence_source': 'fulltext', 'text_scope': 'fulltext', 'evidence_sentence': sent or None, 'variant_in_abstract': link.get('variant_in_abstract')}

def enrich_fulltext_save_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8') as f:
        for row in rows:
            f.write(json.dumps(row, ensure_ascii=False) + '\n')

def enrich_fulltext_emit_article_outputs(pmid: str, pmcid: str | None, title: str | None, source: str | None, concrete_links: list[dict[str, Any]], known_variants: dict[str, list[str]], art: dict[str, Any], ft_pairs: list[dict[str, Any]], ft_corpus: list[dict[str, Any]], stats: dict[str, Any], per_article: list[dict[str, Any]]) -> None:
    abs_vars = {vc_text_module.normalize_variant_token(v) for v in known_variants.get(pmid, []) if v}
    if concrete_links:
        stats['articles_with_ft_links'] += 1
    stats['n_ft_links'] += len([l for l in concrete_links if 'known_variant' not in (l.get('method') or '')])
    stats['n_ft_enrich_links'] += len([l for l in concrete_links if 'known_variant' in (l.get('method') or '')])
    stats['n_ft_only_variant_links'] += sum((1 for l in concrete_links if l.get('variant_discovered_in_fulltext')))
    stats['n_clinical_upgraded'] += sum((1 for l in concrete_links if l.get('clinical_source') not in {None, 'narrative_context'}))
    for link in concrete_links:
        v = vc_text_module.normalize_variant_token(str(link.get('variant') or ''))
        link['variant_in_abstract'] = bool(v and v in abs_vars)
        link['evidence_source'] = 'fulltext'
        link['text_scope'] = 'fulltext'
        pair = {'pmid': pmid, 'pmcid': pmcid, 'tier': link.get('tier'), 'method': link.get('method'), 'variant': link.get('variant'), 'clinical': link.get('clinical'), 'clinical_subtype': link.get('clinical_subtype'), 'clinical_source': link.get('clinical_source'), 'gene': link.get('gene'), 'section': link.get('section'), 'evidence_source': 'fulltext', 'text_scope': 'fulltext', 'variant_in_abstract': link.get('variant_in_abstract'), 'variant_discovered_in_fulltext': link.get('variant_discovered_in_fulltext'), 'evidence_sentence': link.get('evidence_sentence'), 'title': title, 'pubmed_url': art.get('pubmed_url') or f'https://pubmed.ncbi.nlm.nih.gov/{pmid}/', 'ft_source': source}
        ft_pairs.append(pair)
        row = enrich_fulltext_corpus_row(pmid, title, link)
        row['ft_source'] = source
        row['clinical_subtype'] = link.get('clinical_subtype')
        row['clinical_source'] = link.get('clinical_source')
        row['variant_discovered_in_fulltext'] = link.get('variant_discovered_in_fulltext')
        ft_corpus.append(row)
    per_article.append({'pmid': pmid, 'pmcid': pmcid, 'source': source, 'n_links': len(concrete_links)})

def enrich_fulltext_load_articles(path: Path) -> list[dict[str, Any]]:
    rows = []
    with path.open(encoding='utf-8') as f:
        for line in f:
            if line.strip():
                rows.append(json.loads(line))
    return rows

def enrich_fulltext_main() -> None:
    parser = argparse.ArgumentParser(description='PMC full-text variant-clinical enrichment')
    parser.add_argument('--merged', type=Path, default=enrich_fulltext_MERGED)
    parser.add_argument('--linked-dir', type=Path, default=enrich_fulltext_LINKED_DIR)
    parser.add_argument('--out-root', type=Path, default=enrich_fulltext_FT_ROOT)
    parser.add_argument('--scope', choices=['pass', 'all'], default='pass', help='pass=only abstract-gated articles with PMCID; all=every merged article with PMCID')
    parser.add_argument('--limit', type=int, default=0)
    parser.add_argument('--force-refetch', action='store_true')
    parser.add_argument('--only-pmid', action='append', default=[], help='Restrict to one or more PMIDs (repeatable)')
    parser.add_argument('--only-pmcid', action='append', default=[], help='Restrict to one or more PMCIDs (repeatable)')
    parser.add_argument('--retry-missing', action='store_true', help='Only process targets that lack parsed/{pmid}.json')
    parser.add_argument('--merge-existing', action='store_true', help='Merge new FT pairs/corpus into existing linked outputs (default with --retry-missing/--only-*)')
    parser.add_argument('--skip-bioc', action='store_true', help='Skip BioC OA API; go straight to EuropePMC / PMC HTML / PDF')
    parser.add_argument('--skip-europepmc', action='store_true', help='Skip EuropePMC JATS (useful when EPMC is slow/down)')
    parser.add_argument('--prefer-pmc-html', action='store_true', help='Try PMC HTML before BioC (Author Manuscript / free PMC first)')
    parser.add_argument('--rebuild-from-parsed', action='store_true', help='Rebuild pairs/corpus_fulltext from data/pubmed/cache/fulltext/parsed/*.json (no refetch)')
    parser.add_argument('--reprocess-links', action='store_true', default=None, help='Re-run linking from cached raw passages (default with --rebuild-from-parsed)')
    parser.add_argument('--no-reprocess-links', dest='reprocess_links', action='store_false', help='Keep links stored in parsed/*.json when rebuilding outputs')
    parser.add_argument('--pubtator-cache', type=Path, default=enrich_fulltext_LINKED_DIR / 'pubtator_parsed.jsonl', help='PubTator cache for clinical field upgrades on full-text links')
    args = parser.parse_args()
    if args.skip_bioc and (not args.skip_europepmc):
        args.skip_europepmc = True
    only_pmids = {str(x).strip() for x in args.only_pmid if str(x).strip()}
    only_pmcids = {enrich_fulltext_normalize_pmcid(x) for x in args.only_pmcid if enrich_fulltext_normalize_pmcid(x)}
    merge_existing = args.merge_existing or args.retry_missing or bool(only_pmids or only_pmcids)
    if args.rebuild_from_parsed:
        merge_existing = False
    if args.reprocess_links is None:
        args.reprocess_links = args.rebuild_from_parsed
    run_id = datetime.now().strftime('%Y%m%d_%H%M%S')
    logger = enrich_fulltext_setup_logger(run_id)
    raw_dir = args.out_root / 'raw'
    parsed_dir = args.out_root / 'parsed'
    linked_out = args.out_root / 'linked'
    raw_dir.mkdir(parents=True, exist_ok=True)
    parsed_dir.mkdir(parents=True, exist_ok=True)
    linked_out.mkdir(parents=True, exist_ok=True)
    articles = enrich_fulltext_load_articles(args.merged)
    by_pmid = {str(a['pmid']): a for a in articles}
    pass_pmids: set[str] = set()
    pass_path = args.linked_dir / 'articles_pass.jsonl'
    known_variants: dict[str, list[str]] = {}
    if pass_path.exists():
        with pass_path.open(encoding='utf-8') as f:
            for line in f:
                r = json.loads(line)
                pmid = str(r['pmid'])
                pass_pmids.add(pmid)
                known_variants[pmid] = list(r.get('variant_mentions') or [])
    strict_path = args.linked_dir / 'pairs_strict.jsonl'
    if strict_path.exists():
        with strict_path.open(encoding='utf-8') as f:
            for line in f:
                r = json.loads(line)
                pmid = str(r['pmid'])
                known_variants.setdefault(pmid, [])
                if r.get('variant'):
                    known_variants[pmid].append(r['variant'])
    pubtator_by_pmid = enrich_fulltext_load_pubtator_cache(args.pubtator_cache)
    logger.info('Loaded PubTator cache %s (%s docs)', args.pubtator_cache, len(pubtator_by_pmid))
    targets: list[dict[str, Any]] = []
    for art in articles:
        pmcid = enrich_fulltext_normalize_pmcid(art.get('pmcid'))
        if not pmcid:
            continue
        pmid = str(art['pmid'])
        if args.scope == 'pass' and pmid not in pass_pmids:
            continue
        if only_pmids and pmid not in only_pmids:
            continue
        if only_pmcids and pmcid not in only_pmcids:
            continue
        if args.retry_missing and (parsed_dir / f'{pmid}.json').exists():
            continue
        art = dict(art)
        art['pmcid_norm'] = pmcid
        targets.append(art)
        if args.limit and len(targets) >= args.limit:
            break
    logger.info('Full-text targets: %s (scope=%s retry_missing=%s prefer_pmc_html=%s merge_existing=%s rebuild=%s)', len(targets), args.scope, args.retry_missing, args.prefer_pmc_html, merge_existing, args.rebuild_from_parsed)
    stats = {'targets': len(targets), 'bioc_ok': 0, 'europepmc_ok': 0, 'pmc_html_ok': 0, 'pmc_pdf_ok': 0, 'fetch_fail': 0, 'articles_with_ft_links': 0, 'n_ft_links': 0, 'n_ft_enrich_links': 0, 'n_ft_corpus': 0, 'n_ft_only_variant_links': 0, 'n_clinical_upgraded': 0}
    ft_pairs: list[dict[str, Any]] = []
    ft_corpus: list[dict[str, Any]] = []
    per_article: list[dict[str, Any]] = []
    if args.rebuild_from_parsed:
        mode = 'reprocess-links from raw cache' if args.reprocess_links else 'parsed links as-is'
        logger.info('Rebuilding FT outputs from parsed/*.json (%s) …', mode)
        for path in sorted(parsed_dir.glob('*.json')):
            try:
                parsed_doc = json.loads(path.read_text(encoding='utf-8'))
            except Exception as exc:
                logger.warning('Skip bad parsed %s: %s', path, exc)
                continue
            pmid = str(parsed_doc.get('pmid') or path.stem)
            if only_pmids and pmid not in only_pmids:
                continue
            pmcid = enrich_fulltext_normalize_pmcid(parsed_doc.get('pmcid'))
            if only_pmcids and pmcid not in only_pmcids:
                continue
            title = parsed_doc.get('title')
            source = parsed_doc.get('source')
            art = by_pmid.get(pmid) or {}
            if source and str(source).startswith('pmc_html'):
                stats['pmc_html_ok'] += 1
            elif source and str(source).startswith('pmc_pdf'):
                stats['pmc_pdf_ok'] += 1
            elif source and 'bioc' in str(source):
                stats['bioc_ok'] += 1
            elif source and 'europepmc' in str(source):
                stats['europepmc_ok'] += 1
            if args.reprocess_links and pmcid:
                passages, re_source = enrich_fulltext_load_passages_from_raw(pmcid, parsed_doc, raw_dir)
                if not passages:
                    logger.warning('Reprocess skip %s: no raw passages for %s', pmid, pmcid)
                    continue
                source = re_source or source
                body_passages = [p for p in passages if enrich_fulltext_usable_passage(p)]
                discovery_passages = [p for p in body_passages if (p.get('section_type') or '').upper() in vc_text_module.BODY_SECTION_TYPES | {'CASE', 'RESULTS', 'DISCUSS', 'CONCL', 'INTRO', 'OTHER', 'METHODS'}]
                merged_links, _meta = enrich_fulltext_build_fulltext_links(discovery_passages, known_variants.get(pmid, []), pubtator_by_pmid.get(pmid), pmid, title=title or parsed_doc.get('title'), abstract=art.get('abstract') or parsed_doc.get('abstract'))
                concrete_links = [l for l in merged_links if vc_text_module.is_concrete_variant_string(str(l.get('variant') or '')) or l.get('method') == 'fulltext_known_variant_sentence']
                parsed_doc = {**parsed_doc, 'n_links': len(concrete_links), 'links': concrete_links, 'reprocessed_at': datetime.now(timezone.utc).isoformat()}
                path.write_text(json.dumps(parsed_doc, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
            else:
                concrete_links = parsed_doc.get('links') or []
            enrich_fulltext_emit_article_outputs(pmid, pmcid, title, source, concrete_links, known_variants, art, ft_pairs, ft_corpus, stats, per_article)
        stats['targets'] = len(per_article)
        targets = []

    def try_bioc(pmcid: str) -> tuple[list[dict[str, Any]], str | None]:
        raw_path = raw_dir / f'{pmcid}.bioc.json'
        if raw_path.exists() and (not args.force_refetch):
            try:
                coll = json.loads(raw_path.read_text(encoding='utf-8'))
                passages = enrich_fulltext_passages_from_bioc(coll)
                if passages:
                    stats['bioc_ok'] += 1
                    return (passages, 'bioc_cache')
            except Exception as exc:
                logger.warning('Bad BioC cache %s: %s', raw_path, exc)
        if args.skip_bioc:
            return ([], None)
        coll = enrich_fulltext_fetch_bioc(pmcid, logger)
        if coll:
            raw_path.write_text(json.dumps(coll, ensure_ascii=False), encoding='utf-8')
            passages = enrich_fulltext_passages_from_bioc(coll)
            if passages:
                stats['bioc_ok'] += 1
                return (passages, 'bioc')
        return ([], None)

    def try_europepmc(pmcid: str) -> tuple[list[dict[str, Any]], str | None]:
        if args.skip_europepmc:
            return ([], None)
        xml_path = raw_dir / f'{pmcid}.jats.xml'
        xml_text = None
        if xml_path.exists() and (not args.force_refetch):
            xml_text = xml_path.read_text(encoding='utf-8')
        else:
            xml_text = enrich_fulltext_fetch_europepmc_xml(pmcid, logger)
            if xml_text:
                xml_path.write_text(xml_text, encoding='utf-8')
        if not xml_text:
            return ([], None)
        passages = enrich_fulltext_passages_from_jats_xml(xml_text)
        if passages:
            stats['europepmc_ok'] += 1
            return (passages, 'europepmc_jats')
        return ([], None)

    def try_pmc_html(pmcid: str) -> tuple[list[dict[str, Any]], str | None, str | None]:
        html_path = raw_dir / f'{pmcid}.pmc.html'
        html = None
        if html_path.exists() and (not args.force_refetch) and (html_path.stat().st_size > 1000):
            html = html_path.read_text(encoding='utf-8', errors='replace')
            source = 'pmc_html_cache'
        else:
            html = enrich_fulltext_fetch_pmc_html(pmcid, logger)
            source = 'pmc_html'
            if html:
                html_path.write_text(html, encoding='utf-8')
        if not html:
            return ([], None, None)
        passages = enrich_fulltext_passages_from_pmc_html(html)
        body_chars = sum((len(p.get('text') or '') for p in passages))
        if len(passages) < 2 or body_chars < enrich_fulltext_MIN_PMC_HTML_CHARS:
            logger.info('PMC HTML %s thin parse (passages=%s chars=%s); keep for PDF fallback', pmcid, len(passages), body_chars)
            return (passages, source, html)
        stats['pmc_html_ok'] += 1
        return (passages, source, html)

    def try_pmc_pdf(pmcid: str, html: str | None) -> tuple[list[dict[str, Any]], str | None]:
        pdf_path = raw_dir / f'{pmcid}.pmc.pdf'
        got = enrich_fulltext_fetch_pmc_pdf(pmcid, html, pdf_path, logger)
        if not got:
            return ([], None)
        text = enrich_fulltext_pdf_to_text(got, logger)
        txt_path = raw_dir / f'{pmcid}.pmc.pdf.txt'
        if text:
            txt_path.write_text(text, encoding='utf-8')
        passages = enrich_fulltext_passages_from_plain_text(text, source_label='pdf')
        if passages:
            stats['pmc_pdf_ok'] += 1
            return (passages, 'pmc_pdf')
        return ([], None)
    for i, art in enumerate(targets, 1):
        pmid = str(art['pmid'])
        pmcid = art['pmcid_norm']
        title = art.get('title')
        passages: list[dict[str, Any]] = []
        source = None
        html_blob: str | None = None
        html_attempted = False

        def passages_ok(ps: list[dict[str, Any]]) -> bool:
            if len(ps) < 2:
                return False
            return sum((len(p.get('text') or '') for p in ps)) >= enrich_fulltext_MIN_PMC_HTML_CHARS
        order: list[str]
        if args.prefer_pmc_html:
            order = ['pmc_html', 'bioc', 'europepmc', 'pmc_pdf']
        else:
            order = ['bioc', 'europepmc', 'pmc_html', 'pmc_pdf']
        for step in order:
            if passages_ok(passages) and step == 'pmc_pdf':
                break
            if passages_ok(passages) and step != 'pmc_pdf':
                break
            if step == 'bioc':
                ps, src = try_bioc(pmcid)
                if ps:
                    passages, source = (ps, src)
            elif step == 'europepmc':
                ps, src = try_europepmc(pmcid)
                if ps:
                    passages, source = (ps, src)
            elif step == 'pmc_html':
                ps, src, html_blob = try_pmc_html(pmcid)
                html_attempted = True
                if ps and (not passages_ok(passages) or passages_ok(ps)):
                    if not passages_ok(passages) or sum((len(p.get('text') or '') for p in ps)) > sum((len(p.get('text') or '') for p in passages)):
                        passages, source = (ps, src)
            elif step == 'pmc_pdf':
                if html_blob is None and (not html_attempted):
                    html_path = raw_dir / f'{pmcid}.pmc.html'
                    if html_path.exists():
                        html_blob = html_path.read_text(encoding='utf-8', errors='replace')
                    else:
                        html_blob = enrich_fulltext_fetch_pmc_html(pmcid, logger)
                        html_attempted = True
                        if html_blob:
                            html_path.write_text(html_blob, encoding='utf-8')
                elif html_blob is None and html_attempted:
                    pass
                pdf_passages, pdf_source = try_pmc_pdf(pmcid, html_blob)
                if pdf_passages and (not passages_ok(passages) or sum((len(p.get('text') or '') for p in pdf_passages)) > sum((len(p.get('text') or '') for p in passages))):
                    passages, source = (pdf_passages, pdf_source)
        if not passages:
            stats['fetch_fail'] += 1
            if i % 50 == 0 or i == len(targets):
                logger.info('Progress %s/%s fail_so_far=%s', i, len(targets), stats['fetch_fail'])
            continue
        body_passages = [p for p in passages if enrich_fulltext_usable_passage(p)]
        discovery_passages = [p for p in body_passages if (p.get('section_type') or '').upper() in vc_text_module.BODY_SECTION_TYPES | {'CASE', 'RESULTS', 'DISCUSS', 'CONCL', 'INTRO', 'OTHER', 'METHODS'}]
        merged_links, link_meta = enrich_fulltext_build_fulltext_links(discovery_passages, known_variants.get(pmid, []), pubtator_by_pmid.get(pmid), pmid, title=title, abstract=art.get('abstract'))
        for link in merged_links:
            link['pmcid'] = pmcid
            link['title'] = title
        concrete_links = [l for l in merged_links if vc_text_module.is_concrete_variant_string(str(l.get('variant') or '')) or l.get('method') == 'fulltext_known_variant_sentence']
        parsed = {'pmid': pmid, 'pmcid': pmcid, 'title': title, 'source': source, 'n_passages': len(passages), 'n_usable_passages': len(body_passages), 'section_counts': {}, 'n_links': len(concrete_links), 'links': concrete_links, 'link_meta': link_meta, 'parsed_at': datetime.now(timezone.utc).isoformat()}
        sec_counts: dict[str, int] = {}
        for p in passages:
            sec_counts[p.get('section_type') or 'NA'] = sec_counts.get(p.get('section_type') or 'NA', 0) + 1
        parsed['section_counts'] = sec_counts
        (parsed_dir / f'{pmid}.json').write_text(json.dumps(parsed, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        enrich_fulltext_emit_article_outputs(pmid, pmcid, title, source, concrete_links, known_variants, art, ft_pairs, ft_corpus, stats, per_article)
        if i % 25 == 0 or i == len(targets):
            logger.info('Progress %s/%s ft_pairs=%s articles_with_links=%s fail=%s html=%s pdf=%s', i, len(targets), len(ft_pairs), stats['articles_with_ft_links'], stats['fetch_fail'], stats['pmc_html_ok'], stats['pmc_pdf_ok'])
    if merge_existing:
        prev_pairs_path = linked_out / 'pairs_fulltext.jsonl'
        prev_corpus_path = linked_out / 'corpus_fulltext.jsonl'
        prev_index_path = linked_out / 'articles_fulltext_index.jsonl'
        touched = {str(x['pmid']) for x in per_article} | {str(p['pmid']) for p in ft_pairs}

        def load_jsonl(path: Path) -> list[dict[str, Any]]:
            if not path.exists():
                return []
            rows = []
            with path.open(encoding='utf-8') as f:
                for line in f:
                    if line.strip():
                        rows.append(json.loads(line))
            return rows
        old_pairs = [r for r in load_jsonl(prev_pairs_path) if str(r.get('pmid')) not in touched]
        old_corpus = [r for r in load_jsonl(prev_corpus_path) if str(r.get('pmid')) not in touched]
        old_index = [r for r in load_jsonl(prev_index_path) if str(r.get('pmid')) not in touched]
        ft_pairs = old_pairs + ft_pairs
        ft_corpus = old_corpus + ft_corpus
        per_article = old_index + per_article
        logger.info('Merged with existing FT outputs: pairs=%s corpus=%s index=%s (touched_pmids=%s)', len(ft_pairs), len(ft_corpus), len(per_article), len(touched))
    stats['n_ft_corpus'] = len(ft_corpus)
    enrich_fulltext_save_jsonl(linked_out / 'pairs_fulltext.jsonl', ft_pairs)
    enrich_fulltext_save_jsonl(linked_out / 'corpus_fulltext.jsonl', ft_corpus)
    enrich_fulltext_save_jsonl(linked_out / 'articles_fulltext_index.jsonl', per_article)
    enriched_pairs: list[dict[str, Any]] = []
    enriched_corpus: list[dict[str, Any]] = []
    seen_pair: set[tuple[str, str, str]] = set()
    abstract_strict = args.linked_dir / 'pairs_strict.jsonl'
    if abstract_strict.exists():
        with abstract_strict.open(encoding='utf-8') as f:
            for line in f:
                r = json.loads(line)
                r = dict(r)
                r.setdefault('text_scope', 'abstract')
                r.setdefault('evidence_source', 'abstract')
                key = (str(r.get('pmid')), str(r.get('variant')), str(r.get('clinical')))
                if key in seen_pair:
                    continue
                seen_pair.add(key)
                enriched_pairs.append(r)
    for r in ft_pairs:
        r = dict(r)
        r.setdefault('evidence_source', 'fulltext')
        r.setdefault('text_scope', 'fulltext')
        enriched_pairs.append(r)
    abstract_corpus = args.linked_dir / 'corpus_strict.jsonl'
    if abstract_corpus.exists():
        with abstract_corpus.open(encoding='utf-8') as f:
            for line in f:
                r = json.loads(line)
                r = dict(r)
                r.setdefault('text_scope', 'abstract')
                r.setdefault('evidence_source', 'abstract')
                text = r.get('text') or ''
                if text and '[abstract]' not in text and ('[fulltext:' not in text):
                    if text.startswith(f"PMID:{r.get('pmid')}."):
                        r['text'] = text.replace(f"PMID:{r.get('pmid')}.", f"PMID:{r.get('pmid')}. [abstract]", 1)
                    else:
                        r['text'] = f'[abstract] {text}'
                enriched_corpus.append(r)
    for r in ft_corpus:
        r = dict(r)
        r.setdefault('evidence_source', 'fulltext')
        r.setdefault('text_scope', 'fulltext')
        enriched_corpus.append(r)
    enrich_fulltext_save_jsonl(linked_out / 'pairs_enriched.jsonl', enriched_pairs)
    enrich_fulltext_save_jsonl(linked_out / 'corpus_enriched.jsonl', enriched_corpus)
    tsv = linked_out / 'pairs_fulltext_compact.tsv'
    with tsv.open('w', encoding='utf-8') as f:
        f.write('pmid\tpmcid\ttier\tmethod\tsection\tvariant\tclinical\tft_source\tevidence_sentence\n')
        for p in ft_pairs:
            ev = (p.get('evidence_sentence') or '').replace('\t', ' ').replace('\n', ' ')
            f.write(f"{p.get('pmid')}\t{p.get('pmcid')}\t{p.get('tier')}\t{p.get('method')}\t{p.get('section') or ''}\t{(p.get('variant') or '').replace(chr(9), ' ')}\t{(p.get('clinical') or '').replace(chr(9), ' ')}\t{p.get('ft_source') or ''}\t{ev}\n")
    summary = {**stats, 'n_pairs_fulltext': len(ft_pairs), 'n_pairs_enriched': len(enriched_pairs), 'n_corpus_enriched': len(enriched_corpus), 'scope': args.scope, 'retry_missing': args.retry_missing, 'prefer_pmc_html': args.prefer_pmc_html, 'merge_existing': merge_existing, 'created_at': datetime.now(timezone.utc).isoformat()}
    (linked_out / 'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    mirror = args.linked_dir
    enrich_fulltext_save_jsonl(mirror / 'pairs_fulltext.jsonl', ft_pairs)
    enrich_fulltext_save_jsonl(mirror / 'corpus_fulltext.jsonl', ft_corpus)
    enrich_fulltext_save_jsonl(mirror / 'pairs_enriched.jsonl', enriched_pairs)
    enrich_fulltext_save_jsonl(mirror / 'corpus_enriched.jsonl', enriched_corpus)
    (mirror / 'fulltext_summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    logger.info('Summary: %s', json.dumps(summary, ensure_ascii=False))
    logger.info('Done. Outputs in %s and mirrored to %s', linked_out, mirror)

# === materialize_review_text_cache.py ===
"""Build pdf_text_cache/PMID*.txt from Stage3 OA/fulltext (no PDF download by default).

Uses cached fulltext raw + parsed metadata under data/pubmed/cache/fulltext/.
Falls back to abstract from ingest/merged when full text is unavailable.
Optional --fetch-pdf runs download_review_strict_pdfs for gaps only.
"""
import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path
from typing import Any
materialize_review_text_cache_ROOT = Path(__file__).resolve().parents[1]
materialize_review_text_cache_SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(materialize_review_text_cache_SCRIPTS.parent))
sys.path.insert(0, str(materialize_review_text_cache_SCRIPTS))
materialize_review_text_cache_SKIP_SECTIONS = {'REF', 'FIG', 'TABLE', 'COMP_INT', 'AUTH_CONT', 'ABBR'}

def materialize_review_text_cache_load_articles_index(path: Path) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    if not path.is_file():
        return out
    for line in path.open(encoding='utf-8'):
        if not line.strip():
            continue
        row = json.loads(line)
        pmid = str(row.get('pmid') or '').strip()
        if pmid:
            out[pmid] = row
    return out

def materialize_review_text_cache_passages_to_plain(passages: list[dict[str, Any]]) -> str:
    parts: list[str] = []
    for p in passages:
        sec = (p.get('section_type') or '').upper()
        if sec in materialize_review_text_cache_SKIP_SECTIONS:
            continue
        if not enrich_fulltext_usable_passage(p):
            continue
        text = (p.get('text') or '').strip()
        if text:
            parts.append(text)
    return '\n\n'.join(parts)

def materialize_review_text_cache_links_blob(parsed: dict[str, Any]) -> str:
    sents: list[str] = []
    seen: set[str] = set()
    for link in parsed.get('links') or []:
        sent = (link.get('evidence_sentence') or '').strip()
        if sent and sent not in seen:
            seen.add(sent)
            sents.append(sent)
    return '\n\n'.join(sents)

def materialize_review_text_cache_abstract_blob(art: dict[str, Any] | None) -> str:
    if not art:
        return ''
    title = (art.get('title') or '').strip()
    abstract = (art.get('abstract') or '').strip()
    bits = [b for b in (title, abstract) if b]
    return '\n\n'.join(bits)

def materialize_review_text_cache_build_text_for_pmid(pmid: str, *, parsed_dir: Path, raw_dir: Path, articles: dict[str, dict[str, Any]]) -> tuple[str, str]:
    parsed_path = parsed_dir / f'{pmid}.json'
    source = 'none'
    body = ''
    if parsed_path.is_file():
        parsed = json.loads(parsed_path.read_text(encoding='utf-8'))
        pmcid = (parsed.get('pmcid') or '').strip()
        if pmcid:
            passages, src = enrich_fulltext_load_passages_from_raw(pmcid, parsed, raw_dir)
            if passages:
                body = materialize_review_text_cache_passages_to_plain(passages)
                source = src or parsed.get('source') or 'fulltext_raw'
        if len(body) < 400:
            extra = materialize_review_text_cache_links_blob(parsed)
            if extra:
                body = (body + '\n\n' + extra).strip() if body else extra
                source = source if source != 'none' else 'parsed_links'
    if len(body) < 200:
        abs_body = materialize_review_text_cache_abstract_blob(articles.get(pmid))
        if abs_body:
            body = abs_body if not body else body + '\n\n' + abs_body
            source = 'abstract_fallback' if source == 'none' else f'{source}+abstract'
    return (body, source)

def materialize_review_text_cache_pmids_from_jsonl(path: Path) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for line in path.open(encoding='utf-8'):
        if not line.strip():
            continue
        row = json.loads(line)
        pmid = str(row.get('pmid') or '').strip()
        if pmid and pmid not in seen:
            seen.add(pmid)
            out.append(pmid)
    return out

def materialize_review_text_cache_main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jsonl', type=Path, default=materialize_review_text_cache_ROOT / 'data/pubmed/pipeline/s6_review/corpus_enriched_novel_clinvar_absent_review_highconf.jsonl')
    parser.add_argument('--parsed-dir', type=Path, default=enrich_fulltext_PARSED_DIR)
    parser.add_argument('--raw-dir', type=Path, default=enrich_fulltext_RAW_DIR)
    parser.add_argument('--articles', type=Path, default=vc_paths_module.MERGED_ARTICLES)
    parser.add_argument('--out-dir', type=Path, default=vc_paths_module.PDF_TEXT_CACHE)
    parser.add_argument('--force', action='store_true', help='rewrite existing PMID*.txt')
    parser.add_argument('--fetch-pdf', action='store_true', help='After materialize, download PDFs only for PMIDs still lacking text cache')
    parser.add_argument('--min-chars', type=int, default=120, help='treat shorter cache as missing')
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    logger = logging.getLogger('materialize_text_cache')
    args.out_dir.mkdir(parents=True, exist_ok=True)
    articles = materialize_review_text_cache_load_articles_index(args.articles)
    pmids = materialize_review_text_cache_pmids_from_jsonl(args.jsonl)
    logger.info('PMIDs from %s: %s', args.jsonl.name, len(pmids))
    stats = {'written': 0, 'skipped': 0, 'empty': 0, 'sources': {}}
    missing_for_pdf: list[str] = []
    for pmid in pmids:
        cache_path = args.out_dir / f'PMID{pmid}.txt'
        if cache_path.is_file() and (not args.force):
            try:
                if len(cache_path.read_text(encoding='utf-8', errors='replace').strip()) >= args.min_chars:
                    stats['skipped'] += 1
                    continue
            except OSError:
                pass
        body, source = materialize_review_text_cache_build_text_for_pmid(pmid, parsed_dir=args.parsed_dir, raw_dir=args.raw_dir, articles=articles)
        if len(body.strip()) < args.min_chars:
            stats['empty'] += 1
            missing_for_pdf.append(pmid)
            if cache_path.is_file() and args.force:
                cache_path.unlink(missing_ok=True)
            continue
        try:
            ft_rel = str(enrich_fulltext_FT_ROOT.relative_to(materialize_review_text_cache_ROOT))
        except ValueError:
            ft_rel = str(enrich_fulltext_FT_ROOT)
        header = f'# PMID {pmid} source={source} chars={len(body)}\n# materialized from {ft_rel}\n\n'
        cache_path.write_text(header + body, encoding='utf-8', errors='replace')
        stats['written'] += 1
        stats['sources'][source] = stats['sources'].get(source, 0) + 1
    logger.info('done written=%s skipped=%s empty=%s out=%s', stats['written'], stats['skipped'], stats['empty'], args.out_dir)
    if args.fetch_pdf and missing_for_pdf:
        stages_dir = materialize_review_text_cache_SCRIPTS
        legacy = stages_dir / 'legacy' / 'download_review_strict_pdfs.py'
        dl_script = legacy if legacy.is_file() else stages_dir / 'download_review_strict_pdfs.py'
        pmid_file = args.out_dir.parent / '_materialize_missing_pmids.txt'
        pmid_file.write_text('\n'.join(missing_for_pdf) + '\n', encoding='utf-8')
        cmd = [sys.executable, str(dl_script), '--jsonl', str(args.jsonl), '--out-dir', str(vc_paths_module.PDFS), '--pmid-list', str(pmid_file)]
        logger.info('fetch-pdf for %s PMIDs: %s', len(missing_for_pdf), ' '.join(cmd))
        subprocess.run(cmd, cwd=str(materialize_review_text_cache_ROOT), check=False)

# Public helper compatibility aliases.
RAW_DIR = enrich_fulltext_RAW_DIR
normalize_pmcid = enrich_fulltext_normalize_pmcid
load_passages_from_raw = enrich_fulltext_load_passages_from_raw

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
        'enrich': enrich_fulltext_main,
        'materialize': materialize_review_text_cache_main,
    })
