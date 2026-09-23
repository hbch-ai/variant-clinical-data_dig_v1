#!/usr/bin/env python3
"""Physically consolidated pipeline commands. Generated from the v4_01 stage sources.

Each section retains its source module boundary through descriptive global prefixes;
there is no dynamic source loading or embedded executable source text.
"""
from __future__ import annotations
import review as review_module
import validate_novelty_clinvar as validate_novelty_clinvar_module
import vc_ensembl_c_map as vc_ensembl_c_map_module
import vc_paths as vc_paths_module
import vc_text as vc_text_module

# === vc_intergenic.py ===
"""Intergenic / regulatory variant patterns (bypass track; do not merge into coding VARIANT_PATTERNS)."""
import re
from dataclasses import dataclass
from typing import Any
vc_intergenic_INTERGENIC_VARIANT_PATTERNS: list[tuple[str, re.Pattern[str]]] = [('hgvs_g_nc', re.compile('\\b(?:NC_\\d+(?:\\.\\d+)?|NG_\\d+(?:\\.\\d+)?):g\\.(?:\\d+(?:_\\d+)?)(?:delins|del|dup|ins|[ACGT]\\s*>\\s*[ACGT]|[ACGT]+)\\b', re.I)), ('hgvs_g', re.compile('\\bg\\.(?:\\d+(?:_\\d+)?)(?:delins|del|dup|ins|[ACGT]\\s*>\\s*[ACGT]|[ACGT]+)\\b', re.I)), ('genomic_chrpos_change', re.compile('\\b(?:(?:GRCh3[78]|hg(?:19|38))\\s*)?(?:chr)?([0-9]{1,2}|[XYM])\\s*[:\\-]\\s*([\\d,]+)\\s*([ACGT])\\s*(?:>|/|→)\\s*([ACGT])\\b', re.I)), ('genomic_chrpos', re.compile('\\b(?:(?:GRCh3[78]|hg(?:19|38))\\s*)?(?:chr)?([0-9]{1,2}|[XYM])\\s*[:\\-]\\s*([\\d,]{4,})\\b', re.I)), ('rsid', re.compile('\\brs\\d{4,}\\b', re.I))]
vc_intergenic_INTERGENIC_CONCRETE = {'hgvs_g_nc', 'hgvs_g', 'genomic_chrpos_change', 'genomic_chrpos', 'rsid'}
vc_intergenic_REGION_CUE_PATTERNS: list[tuple[str, re.Pattern[str]]] = [('intergenic', re.compile('\\bintergenic\\b', re.I)), ('enhancer', re.compile('\\benhancer(?:s)?\\b|\\benhancer\\s+(?:SNP|variant|mutation)s?\\b|\\bregulatory\\s+(?:variant|SNP|element|region)s?\\b|\\bnon[- ]?coding\\s+regulatory\\b|\\binsulator\\b|\\bTAD\\b', re.I)), ('upstream', re.compile("\\bupstream\\s+of\\b|\\b\\d+\\s*k?bp\\s+upstream\\b|\\b5['′]\\s*(?:flanking|region)\\b", re.I)), ('downstream', re.compile("\\bdownstream\\s+of\\b|\\b\\d+\\s*k?bp\\s+downstream\\b|\\b3['′]\\s*(?:flanking|region)\\b", re.I))]
vc_intergenic_INTRON_UTR_ONLY = re.compile("\\b(?:intronic|deep\\s+intronic|intron\\s+\\d+|3['′]\\s*UTR|5['′]\\s*UTR|untranslated)\\b", re.I)
vc_intergenic_ASSEMBLY_PAT = re.compile('\\b(GRCh3[78]|hg19|hg38)\\b', re.I)
vc_intergenic_ANCHOR_GENE_PATS: list[re.Pattern[str]] = [re.compile('\\b(?:upstream|downstream|near|adjacent\\s+to|close\\s+to)\\s+(?:of\\s+|the\\s+)?([A-Z][A-Z0-9]{1,14}(?:-[A-Z0-9]+)*)\\b'), re.compile('\\bintergenic\\s+(?:region\\s+)?(?:between|of)\\s+([A-Z][A-Z0-9]{1,14})\\s+(?:and|/|-)\\s+([A-Z][A-Z0-9]{1,14})\\b'), re.compile('\\b([A-Z][A-Z0-9]{1,14})-([A-Z][A-Z0-9]{1,14})\\s+intergenic\\b')]

@dataclass
class vc_intergenic_GenomicLocus:
    assembly: str | None
    chrom: str | None
    pos: int | None
    ref: str | None
    alt: str | None
    rs: str | None
    raw: str
    subtype: str

    def locus_key(self) -> str:
        if self.rs:
            return self.rs.lower()
        ass = (self.assembly or 'unknown').lower()
        chrom = self.chrom or 'NA'
        pos = self.pos if self.pos is not None else 'NA'
        if self.ref and self.alt:
            return f'{ass}|{chrom}|{pos}|{self.ref}>{self.alt}'
        return f'{ass}|{chrom}|{pos}'

def vc_intergenic_detect_assembly(text: str) -> str | None:
    m = vc_intergenic_ASSEMBLY_PAT.search(text or '')
    if not m:
        return None
    a = m.group(1)
    if a.lower() == 'hg19':
        return 'GRCh37'
    if a.lower() == 'hg38':
        return 'GRCh38'
    return a

def vc_intergenic_classify_region(text: str) -> tuple[str | None, list[str]]:
    """Return (best region_class, cue labels). Prefer intergenic > enhancer > up/down."""
    cues: list[str] = []
    for label, pat in vc_intergenic_REGION_CUE_PATTERNS:
        if pat.search(text or ''):
            cues.append(label)
    if not cues:
        return (None, [])
    priority = ('intergenic', 'enhancer', 'upstream', 'downstream')
    for p in priority:
        if p in cues:
            return (p, cues)
    return (cues[0], cues)

def vc_intergenic_extract_anchor_genes(text: str) -> tuple[str | None, str | None, str | None]:
    """Return (anchor_gene, anchor_gene_b, distance_cue)."""
    t = text or ''
    for pat in vc_intergenic_ANCHOR_GENE_PATS:
        m = pat.search(t)
        if not m:
            continue
        g1 = m.group(1)
        g2 = m.group(2) if m.lastindex and m.lastindex >= 2 else None
        if not vc_text_module.is_plausible_gene_symbol(g1):
            continue
        if g2 and (not vc_text_module.is_plausible_gene_symbol(g2)):
            g2 = None
        cue = None
        low = t[max(0, m.start() - 30):m.end() + 10].lower()
        if 'upstream' in low:
            cue = 'upstream'
        elif 'downstream' in low:
            cue = 'downstream'
        elif 'intergenic' in low:
            cue = 'intergenic_between' if g2 else 'intergenic'
        elif 'near' in low or 'adjacent' in low:
            cue = 'near'
        return (g1.upper(), g2.upper() if g2 else None, cue)
    return (None, None, None)

def vc_intergenic__norm_chrom(c: str) -> str:
    c = c.upper().lstrip('CHR')
    return c

def vc_intergenic__norm_pos(p: str) -> int | None:
    try:
        return int(p.replace(',', ''))
    except ValueError:
        return None

def vc_intergenic_parse_genomic_locus(raw: str, subtype: str, context: str='') -> vc_intergenic_GenomicLocus:
    text = (raw or '').strip()
    assembly = vc_intergenic_detect_assembly(context) or vc_intergenic_detect_assembly(text)
    rs = None
    chrom = pos = ref = alt = None
    if subtype == 'rsid':
        rs = text.lower()
    elif subtype in {'genomic_chrpos_change', 'genomic_chrpos'}:
        m = re.search('(?:chr)?([0-9]{1,2}|[XYM])\\s*[:\\-]\\s*([\\d,]+)(?:\\s*([ACGT])\\s*(?:>|/|→)\\s*([ACGT]))?', text, re.I)
        if m:
            chrom = vc_intergenic__norm_chrom(m.group(1))
            pos = vc_intergenic__norm_pos(m.group(2))
            ref = m.group(3).upper() if m.group(3) else None
            alt = m.group(4).upper() if m.group(4) else None
    elif subtype in {'hgvs_g', 'hgvs_g_nc'}:
        m = re.search('(?:NC_\\d+(?:\\.\\d+)?:)?g\\.(\\d+)(?:_(\\d+))?(?:([ACGT])\\s*>\\s*([ACGT])|(delins|del|dup|ins))', text, re.I)
        if m:
            pos = int(m.group(1))
            if m.group(3) and m.group(4):
                ref, alt = (m.group(3).upper(), m.group(4).upper())
        nc = re.search('NC_0*(\\d+)', text, re.I)
        if nc:
            n = int(nc.group(1))
            if 1 <= n <= 22:
                chrom = str(n)
            elif n == 23:
                chrom = 'X'
            elif n == 24:
                chrom = 'Y'
    return vc_intergenic_GenomicLocus(assembly=assembly, chrom=chrom, pos=pos, ref=ref, alt=alt, rs=rs, raw=text, subtype=subtype)

def vc_intergenic_find_intergenic_variant_spans(text: str) -> list[dict[str, Any]]:
    """Non-overlapping concrete intergenic/genomic spans (prefer longer / more specific)."""
    hits: list[tuple[int, int, str, str]] = []
    for subtype, pat in vc_intergenic_INTERGENIC_VARIANT_PATTERNS:
        for m in pat.finditer(text or ''):
            hits.append((m.start(), m.end(), m.group(0), subtype))
    hits.sort(key=lambda x: (x[0], -(x[1] - x[0])))
    out: list[dict[str, Any]] = []
    occupied: list[tuple[int, int]] = []
    for start, end, raw, subtype in hits:
        if any((not (end <= a or start >= b) for a, b in occupied)):
            continue
        occupied.append((start, end))
        locus = vc_intergenic_parse_genomic_locus(raw, subtype, context=text)
        out.append({'start': start, 'end': end, 'text': raw.strip(), 'subtype': subtype, 'locus': locus})
    return out

def vc_intergenic_sentence_allows_intergenic_variant(sentence: str, *, subtype: str) -> tuple[bool, str | None, list[str]]:
    """
    Gate: region cue required. rsID alone needs intergenic/enhancer/up/down cue;
    intron/UTR-only without those cues is rejected.
    """
    region, cues = vc_intergenic_classify_region(sentence)
    if region:
        return (True, region, cues)
    if subtype in {'hgvs_g', 'hgvs_g_nc', 'genomic_chrpos', 'genomic_chrpos_change'}:
        if vc_intergenic_INTRON_UTR_ONLY.search(sentence) and (not region):
            return (False, None, [])
        return (True, 'unknown_noncoding', [])
    return (False, None, [])

# === vc_pdf_genomic.py ===
"""Extract literal genomic coordinates from cached PDF full text."""
import re
from pathlib import Path
from typing import Any
vc_pdf_genomic_RE_NC_G = re.compile('\\b(NC_\\d+\\.\\d+):g\\.(?:\\d+(?:_\\d+)?)(?:delins|del|dup|ins|[ACGT]\\s*>\\s*[ACGT]|[ACGT]+)\\b', re.I)
vc_pdf_genomic_RE_G_SUB = re.compile('\\bg\\.(\\d+)\\s*([ACGT])\\s*>\\s*([ACGT])\\b', re.I)
vc_pdf_genomic_RE_G_RANGE = re.compile('\\bg\\.(\\d+)(?:_(\\d+))?(del|dup|ins|delins)([ACGT]*)\\b', re.I)
vc_pdf_genomic_RE_VCFISH = re.compile('\\b(?:chr)?([0-9]{1,2}|X|Y|MT)\\s+(\\d{5,9})\\s+([ACGTN]+)\\s+([ACGTN]+)\\b', re.I)

def vc_pdf_genomic_load_pdf_text(pmid: str, cache_dir: Path) -> str | None:
    if not pmid or not cache_dir.exists():
        return None
    p = cache_dir / f'PMID{pmid}.txt'
    if not p.exists():
        hits = list(cache_dir.glob(f'*{pmid}*'))
        if not hits:
            return None
        p = hits[0]
    try:
        return p.read_text(encoding='utf-8', errors='replace')
    except OSError:
        return None

def vc_pdf_genomic__windows(text: str, anchors: list[str], radius: int=220) -> list[str]:
    if not text:
        return []
    low = text.lower()
    out: list[str] = []
    seen: set[tuple[int, int]] = set()
    for a in anchors:
        if not a or len(a) < 2:
            continue
        start = 0
        needle = a.lower()
        while True:
            i = low.find(needle, start)
            if i < 0:
                break
            lo = max(0, i - radius)
            hi = min(len(text), i + len(a) + radius)
            key = (lo, hi)
            if key not in seen:
                seen.add(key)
                out.append(text[lo:hi])
            start = i + len(needle)
    if not out and text:
        out.append(text)
    return out

def vc_pdf_genomic__parse_nc_g(raw: str, assembly_hint: str | None) -> dict[str, Any] | None:
    m = vc_pdf_genomic_RE_NC_G.search(raw)
    if not m:
        return None
    nc = m.group(1)
    chrom = None
    m2 = re.search('NC_0*(\\d{2})', nc, re.I)
    if m2:
        n = int(m2.group(1))
        chrom = str(n) if n <= 22 else 'X' if n == 23 else 'Y' if n == 24 else 'MT' if n == 129 else None
    loc = vc_intergenic_parse_genomic_locus(raw, 'hgvs_g_nc')
    if not loc:
        return None
    ass = assembly_hint or loc.assembly or 'GRCh38'
    ref = (loc.ref or '').upper()
    alt = (loc.alt or '').upper()
    if not chrom:
        chrom = loc.chrom
    if not loc.pos or not ref or (not alt):
        return None
    return {'chromosome': chrom, 'position': int(loc.pos), 'ref': ref, 'alt': alt, 'assembly': ass, 'raw': raw}

def vc_pdf_genomic__parse_span(sp: dict[str, Any], assembly_hint: str | None) -> dict[str, Any] | None:
    loc = sp.get('locus')
    if not loc or not loc.chrom or (not loc.pos):
        return None
    ref = (loc.ref or '').upper()
    alt = (loc.alt or '').upper()
    ass = assembly_hint or loc.assembly or 'GRCh38'
    if sp['subtype'] in {'genomic_chrpos_change', 'hgvs_g', 'hgvs_g_nc'} and ref and alt:
        return {'chromosome': loc.chrom, 'position': int(loc.pos), 'ref': ref, 'alt': alt, 'assembly': ass, 'raw': loc.raw}
    if sp['subtype'] == 'genomic_chrpos':
        return {'chromosome': loc.chrom, 'position': int(loc.pos), 'ref': ref, 'alt': alt, 'assembly': ass, 'raw': loc.raw, 'partial': not (ref and alt)}
    return None

def vc_pdf_genomic_extract_pdf_genomic_hits(*, pdf_text: str, gene: str | None, variant_text: str, c_hgvs: str | None=None, p_hgvs: str | None=None) -> list[dict[str, Any]]:
    """Return candidate genomic loci mentioned near the variant in PDF text."""
    if not pdf_text:
        return []
    anchors = [variant_text, c_hgvs or '', p_hgvs or '', gene or '']
    hits: list[dict[str, Any]] = []
    for win in vc_pdf_genomic__windows(pdf_text, anchors):
        ass_m = vc_intergenic_ASSEMBLY_PAT.search(win)
        ass = ass_m.group(1) if ass_m else None
        for sp in vc_intergenic_find_intergenic_variant_spans(win):
            parsed = vc_pdf_genomic__parse_span(sp, ass)
            if parsed:
                parsed['window'] = win[:120]
                hits.append(parsed)
        nc = vc_pdf_genomic__parse_nc_g(win, ass)
        if nc:
            nc['window'] = win[:120]
            hits.append(nc)
        for m in vc_pdf_genomic_RE_VCFISH.finditer(win):
            hits.append({'chromosome': m.group(1).upper().replace('M', 'MT'), 'position': int(m.group(2)), 'ref': m.group(3).upper(), 'alt': m.group(4).upper(), 'assembly': ass or 'GRCh38', 'raw': m.group(0), 'window': win[:120]})
    seen: set[tuple] = set()
    uniq: list[dict[str, Any]] = []
    for h in hits:
        key = (h.get('chromosome'), h.get('position'), h.get('ref'), h.get('alt'))
        if key in seen:
            continue
        seen.add(key)
        uniq.append(h)
    return uniq

def vc_pdf_genomic_pick_pdf_genomic_hit(hits: list[dict[str, Any]], *, gene: str | None, variant_text: str) -> dict[str, Any] | None:
    if not hits:
        return None
    complete = [h for h in hits if h.get('ref') and h.get('alt') and h.get('chromosome') and h.get('position')]
    pool = complete or hits

    def score(h: dict[str, Any]) -> tuple[int, int]:
        w = (h.get('window') or '').lower()
        s = 0
        if gene and gene.lower() in w:
            s += 2
        if variant_text and variant_text.lower()[:12] in w:
            s += 2
        if h.get('ref') and h.get('alt'):
            s += 1
        return (s, 1 if h.get('partial') else 0)
    pool.sort(key=score, reverse=True)
    return pool[0]

# === vc_gff3_region.py ===
"""Annotate GRCh38 1-based positions against Ensembl/GENCODE GFF3.

After Stage9 resolves chr/pos, classify:
  coding / utr / intron / ncrna_exon / ncrna_gene / pseudogene / intergenic / "" (no coordinates).

Coordinates: GFF3 and VCF/HGVS are both 1-based inclusive.
"""
import gzip
import logging
import pickle
from pathlib import Path
from typing import Any, Iterable
vc_gff3_region_KEEP_SEQ = {str(i) for i in range(1, 23)} | {'X', 'Y', 'MT'}
vc_gff3_region_GENE_TYPES = {'gene', 'ncRNA_gene', 'pseudogene'}
vc_gff3_region_TX_TYPES = {'mRNA', 'transcript', 'lnc_RNA', 'miRNA', 'snRNA', 'snoRNA', 'rRNA', 'ncRNA', 'scRNA', 'sRNA', 'snoRNA_pseudogene', 'tRNA'}
vc_gff3_region_UTR_TYPES = {'five_prime_UTR', 'three_prime_UTR'}
vc_gff3_region_FEATURE_TYPES = vc_gff3_region_GENE_TYPES | vc_gff3_region_TX_TYPES | {'exon', 'CDS'} | vc_gff3_region_UTR_TYPES
vc_gff3_region_PROTEIN_BIOTYPES = {'protein_coding'}
vc_gff3_region_NCRNA_BIOTYPES = {'lncRNA', 'lincRNA', 'miRNA', 'snRNA', 'snoRNA', 'rRNA', 'misc_RNA', 'Mt_rRNA', 'Mt_tRNA', 'scRNA', 'sRNA', 'scaRNA', 'vault_RNA', 'ribozyme', 'processed_transcript', 'non_coding', '3prime_overlapping_ncRNA', 'sense_intronic', 'sense_overlapping', 'antisense', 'bidirectional_promoter_lncRNA', 'macro_lncRNA', 'TEC', 'tec'}
vc_gff3_region_PSEUDO_BIOTYPES = {'pseudogene', 'processed_pseudogene', 'unprocessed_pseudogene', 'transcribed_processed_pseudogene', 'transcribed_unprocessed_pseudogene', 'transcribed_unitary_pseudogene', 'unitary_pseudogene', 'polymorphic_pseudogene', 'rRNA_pseudogene', 'translated_processed_pseudogene'}
vc_gff3_region_REGION_CODING = 'coding'
vc_gff3_region_REGION_UTR = 'utr'
vc_gff3_region_REGION_INTRON = 'intron'
vc_gff3_region_REGION_NCRNA_EXON = 'ncrna_exon'
vc_gff3_region_REGION_NCRNA_GENE = 'ncrna_gene'
vc_gff3_region_REGION_PSEUDOGENE = 'pseudogene'
vc_gff3_region_REGION_INTERGENIC = 'intergenic'
vc_gff3_region_REGION_UNRESOLVED = ''
vc_gff3_region_CACHE_VERSION = 2
vc_gff3_region_BIN_SIZE = 250000

def vc_gff3_region_parse_gff_attrs(raw: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in raw.strip().split(';'):
        if not part or '=' not in part:
            continue
        key, val = part.split('=', 1)
        out[key] = val
    return out

def vc_gff3_region__strip_id(raw: str | None) -> str:
    if not raw:
        return ''
    text = raw.split(',', 1)[0]
    if ':' in text:
        return text.split(':', 1)[1]
    return text

def vc_gff3_region__norm_chrom(raw: str) -> str | None:
    c = raw.strip()
    if c.lower().startswith('chr'):
        c = c[3:]
    if c in {'M', 'chrM'}:
        c = 'MT'
    return c if c in vc_gff3_region_KEEP_SEQ else None

def vc_gff3_region__biotype_class(biotype: str) -> str:
    b = (biotype or '').strip()
    if b in vc_gff3_region_PROTEIN_BIOTYPES:
        return 'protein'
    if b in vc_gff3_region_NCRNA_BIOTYPES:
        return 'ncrna'
    if b in vc_gff3_region_PSEUDO_BIOTYPES or 'pseudogene' in b.lower():
        return 'pseudo'
    return 'other'

class vc_gff3_region_Gff3RegionIndex:
    """Binned gene / exon / CDS intervals for point-or-span overlap."""

    def __init__(self) -> None:
        self.genes: dict[str, dict[int, list[tuple]]] = {}
        self.exons: dict[str, dict[int, list[tuple]]] = {}
        self.cds: dict[str, dict[int, list[tuple]]] = {}
        self.utr5: dict[str, dict[int, list[tuple]]] = {}
        self.utr3: dict[str, dict[int, list[tuple]]] = {}
        self.tx_strand: dict[str, str] = {}
        self.tx_cds_span: dict[str, tuple[int, int]] = {}
        self.source = ''
        self.n_genes = 0
        self.n_exons = 0
        self.n_cds = 0
        self.n_utr5 = 0
        self.n_utr3 = 0

    def _add(self, store: dict[str, dict[int, list[tuple]]], chrom: str, start: int, end: int, rec: tuple) -> None:
        bins = store.setdefault(chrom, {})
        lo = start // vc_gff3_region_BIN_SIZE
        hi = end // vc_gff3_region_BIN_SIZE
        for b in range(lo, hi + 1):
            bins.setdefault(b, []).append(rec)

    def add_gene(self, chrom: str, start: int, end: int, name: str, biotype: str, gene_id: str) -> None:
        self._add(self.genes, chrom, start, end, (start, end, name, biotype, gene_id))
        self.n_genes += 1

    def add_exon(self, chrom: str, start: int, end: int, name: str, biotype: str, gene_id: str, tx_id: str) -> None:
        self._add(self.exons, chrom, start, end, (start, end, name, biotype, gene_id, tx_id))
        self.n_exons += 1

    def add_cds(self, chrom: str, start: int, end: int, name: str, gene_id: str, tx_id: str) -> None:
        self._add(self.cds, chrom, start, end, (start, end, name, gene_id, tx_id))
        self.n_cds += 1
        lo, hi = self.tx_cds_span.get(tx_id, (start, end))
        self.tx_cds_span[tx_id] = (min(lo, start), max(hi, end))

    def add_utr(self, which: str, chrom: str, start: int, end: int, name: str, gene_id: str, tx_id: str) -> None:
        rec = (start, end, name, gene_id, tx_id)
        if which == '5':
            self._add(self.utr5, chrom, start, end, rec)
            self.n_utr5 += 1
        else:
            self._add(self.utr3, chrom, start, end, rec)
            self.n_utr3 += 1

    def _hits(self, store: dict, chrom: str, start: int, end: int) -> list[tuple]:
        bins = store.get(chrom) or {}
        lo = start // vc_gff3_region_BIN_SIZE
        hi = end // vc_gff3_region_BIN_SIZE
        seen: set[tuple] = set()
        out: list[tuple] = []
        for b in range(lo, hi + 1):
            for rec in bins.get(b) or []:
                if rec[1] < start or rec[0] > end:
                    continue
                if rec in seen:
                    continue
                seen.add(rec)
                out.append(rec)
        return out

    def annotate(self, chrom: str | None, position: Any, ref: str | None=None, catalog_gene: str | None=None) -> dict[str, str]:
        empty = {'genomic_region': vc_gff3_region_REGION_UNRESOLVED, 'genomic_region_genes': '', 'genomic_region_feature': '', 'genomic_region_transcripts': '', 'genomic_region_note': ''}
        chrom_n = vc_gff3_region__norm_chrom(str(chrom)) if chrom else None
        try:
            pos = int(position)
        except (TypeError, ValueError):
            return empty
        if not chrom_n or pos <= 0:
            return empty
        span = max(1, len(str(ref).strip()) if ref else 1)
        start, end = (pos, pos + span - 1)
        cds_hits = self._hits(self.cds, chrom_n, start, end)
        utr5_hits = self._hits(self.utr5, chrom_n, start, end)
        utr3_hits = self._hits(self.utr3, chrom_n, start, end)
        exon_hits = self._hits(self.exons, chrom_n, start, end)
        gene_hits = self._hits(self.genes, chrom_n, start, end)
        protein_exons = [h for h in exon_hits if vc_gff3_region__biotype_class(h[3]) == 'protein']
        ncrna_exons = [h for h in exon_hits if vc_gff3_region__biotype_class(h[3]) == 'ncrna']
        protein_genes = [h for h in gene_hits if vc_gff3_region__biotype_class(h[3]) == 'protein']
        ncrna_genes = [h for h in gene_hits if vc_gff3_region__biotype_class(h[3]) == 'ncrna']
        pseudo_genes = [h for h in gene_hits if vc_gff3_region__biotype_class(h[3]) == 'pseudo']
        region = vc_gff3_region_REGION_INTERGENIC
        feature = ''
        tx_ids: list[str] = []
        names: list[str] = []
        if cds_hits:
            region = vc_gff3_region_REGION_CODING
            feature = 'CDS'
            names = vc_gff3_region__uniq_names((h[2] or h[3] for h in cds_hits))
            tx_ids = vc_gff3_region__uniq_names((h[4] for h in cds_hits))
        elif utr5_hits or utr3_hits or protein_exons:
            region = vc_gff3_region_REGION_UTR
            feature = 'UTR' if utr5_hits or utr3_hits else 'exon'
            names = vc_gff3_region__uniq_names((h[2] or h[3] for h in list(utr5_hits) + list(utr3_hits))) or vc_gff3_region__uniq_names((h[2] or h[4] for h in protein_exons))
            tx_ids = vc_gff3_region__uniq_names((h[4] for h in list(utr5_hits) + list(utr3_hits))) or vc_gff3_region__uniq_names((h[5] for h in protein_exons))
        elif protein_genes:
            region = vc_gff3_region_REGION_INTRON
            feature = 'gene'
            names = vc_gff3_region__uniq_names((h[2] or h[4] for h in protein_genes))
        elif ncrna_exons:
            region = vc_gff3_region_REGION_NCRNA_EXON
            feature = 'exon'
            names = vc_gff3_region__uniq_names((h[2] or h[4] for h in ncrna_exons))
            tx_ids = vc_gff3_region__uniq_names((h[5] for h in ncrna_exons))
        elif ncrna_genes:
            region = vc_gff3_region_REGION_NCRNA_GENE
            feature = 'gene'
            names = vc_gff3_region__uniq_names((h[2] or h[4] for h in ncrna_genes))
        elif pseudo_genes:
            region = vc_gff3_region_REGION_PSEUDOGENE
            feature = 'gene'
            names = vc_gff3_region__uniq_names((h[2] or h[4] for h in pseudo_genes))
        elif gene_hits:
            region = vc_gff3_region_REGION_INTRON
            feature = 'gene'
            names = vc_gff3_region__uniq_names((h[2] or h[4] for h in gene_hits))
        else:
            region = vc_gff3_region_REGION_INTERGENIC
            feature = ''
        note = ''
        cat = (catalog_gene or '').strip().upper()
        if cat and names and (cat not in {n.upper() for n in names}):
            note = f"catalog_gene {cat} not overlapping position; GFF3={','.join(names)}"
        elif cat and region == vc_gff3_region_REGION_INTERGENIC:
            note = f'catalog_gene {cat} not overlapping position; intergenic'
        return {'genomic_region': region, 'genomic_region_genes': ';'.join(names), 'genomic_region_feature': feature, 'genomic_region_transcripts': ';'.join(tx_ids), 'genomic_region_note': note}

def vc_gff3_region__uniq_names(vals: Iterable[str | None]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for v in vals:
        s = str(v or '').strip()
        if not s:
            continue
        key = s.upper()
        if key in seen:
            continue
        seen.add(key)
        out.append(s)
    return out

def vc_gff3_region_load_gff3_region_index(gff3: Path, *, cache_path: Path | None=None, logger: logging.Logger | None=None) -> vc_gff3_region_Gff3RegionIndex:
    gff3 = Path(gff3)
    if cache_path is not None:
        cached = vc_gff3_region__try_load_cache(cache_path, gff3)
        if cached is not None:
            if logger:
                logger.info('GFF3 region index cache hit %s (genes=%s exons=%s CDS=%s utr5=%s utr3=%s)', cache_path, cached.n_genes, cached.n_exons, cached.n_cds, cached.n_utr5, cached.n_utr3)
            return cached
    idx = vc_gff3_region_build_gff3_region_index(gff3, logger=logger)
    if cache_path is not None:
        vc_gff3_region__save_cache(cache_path, gff3, idx)
        if logger:
            logger.info('wrote GFF3 region cache %s', cache_path)
    return idx

def vc_gff3_region_build_gff3_region_index(gff3: Path, *, logger: logging.Logger | None=None) -> vc_gff3_region_Gff3RegionIndex:
    idx = vc_gff3_region_Gff3RegionIndex()
    idx.source = str(gff3)
    genes: dict[str, tuple[str, str]] = {}
    tx_to_gene: dict[str, str] = {}
    opener = gzip.open if str(gff3).endswith('.gz') else open
    n_line = 0
    with opener(gff3, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            n_line += 1
            if not line or line[0] == '#':
                continue
            cols = line.split('\t')
            if len(cols) < 9:
                continue
            ftype = cols[2]
            if ftype not in vc_gff3_region_FEATURE_TYPES:
                continue
            chrom = vc_gff3_region__norm_chrom(cols[0])
            if not chrom:
                continue
            try:
                start = int(cols[3])
                end = int(cols[4])
            except ValueError:
                continue
            attrs = vc_gff3_region_parse_gff_attrs(cols[8])
            fid = vc_gff3_region__strip_id(attrs.get('ID'))
            parent = vc_gff3_region__strip_id(attrs.get('Parent'))
            biotype = attrs.get('biotype') or ''
            name = attrs.get('Name') or ''
            strand = (cols[6] or '').strip()
            if ftype in vc_gff3_region_GENE_TYPES:
                gid = fid or parent
                if not gid:
                    continue
                display = name or gid
                genes[gid] = (display, biotype)
                idx.add_gene(chrom, start, end, display, biotype, gid)
                continue
            if ftype in vc_gff3_region_TX_TYPES:
                gid = parent
                if fid and gid:
                    tx_to_gene[fid] = gid
                if fid and strand in {'+', '-'}:
                    idx.tx_strand[fid] = strand
                continue
            gid = tx_to_gene.get(parent) or (parent if parent in genes else '')
            gname, gbio = genes.get(gid, (name or gid, biotype))
            if ftype == 'exon':
                idx.add_exon(chrom, start, end, gname, gbio or biotype, gid, parent)
            elif ftype == 'CDS':
                idx.add_cds(chrom, start, end, gname, gid, parent)
            elif ftype == 'five_prime_UTR':
                idx.add_utr('5', chrom, start, end, gname, gid, parent)
            elif ftype == 'three_prime_UTR':
                idx.add_utr('3', chrom, start, end, gname, gid, parent)
    if logger:
        logger.info('GFF3 region index from %s: genes=%s exons=%s CDS=%s utr5=%s utr3=%s lines=%s', gff3, idx.n_genes, idx.n_exons, idx.n_cds, idx.n_utr5, idx.n_utr3, n_line)
    return idx

def vc_gff3_region_default_region_cache(root: Path | None=None) -> Path:
    base = root or Path(__file__).resolve().parents[1]
    return base / 'data' / 'pubmed' / 'cache' / 'ensembl' / 'ensembl_gff3_region_index.pkl'

def vc_gff3_region__cache_meta(gff3: Path) -> dict[str, Any]:
    st = gff3.stat()
    return {'version': vc_gff3_region_CACHE_VERSION, 'gff3': str(gff3.resolve()), 'size': st.st_size, 'mtime_ns': st.st_mtime_ns}

def vc_gff3_region__try_load_cache(cache_path: Path, gff3: Path) -> vc_gff3_region_Gff3RegionIndex | None:
    path = Path(cache_path)
    if not path.is_file():
        return None
    try:
        with path.open('rb') as fh:
            blob = pickle.load(fh)
    except Exception:
        return None
    if not isinstance(blob, dict) or blob.get('meta') != vc_gff3_region__cache_meta(gff3):
        return None
    idx = blob.get('index')
    return idx if isinstance(idx, vc_gff3_region_Gff3RegionIndex) else None

def vc_gff3_region__save_cache(cache_path: Path, gff3: Path, idx: vc_gff3_region_Gff3RegionIndex) -> None:
    path = Path(cache_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    with tmp.open('wb') as fh:
        pickle.dump({'meta': vc_gff3_region__cache_meta(gff3), 'index': idx}, fh, protocol=4)
    tmp.replace(path)

# === vc_dbsnp.py ===
"""Fetch GRCh38 genomic alleles for rsIDs from the NCBI dbSNP / Variation API.

Primary endpoint:
  https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/<id>

Coordinates use SPDI (0-based). Exported positions are converted to 1-based
VCF/HGVS style to match the rest of the genomic allele CSV.
"""
import json
import logging
import os
import re
import time
import ssl
import subprocess
import http.client
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any
vc_dbsnp_REFSNP_API = 'https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{rs_num}'
vc_dbsnp_DBSNP_PAGE = 'https://www.ncbi.nlm.nih.gov/snp/{rs_id}'
vc_dbsnp_USER_AGENT = 'variant_clinical_v2/1.0 (dbsnp-fallback)'
vc_dbsnp_NC_TO_CHROM_GRCH38 = {'NC_000001.11': '1', 'NC_000002.12': '2', 'NC_000003.12': '3', 'NC_000004.12': '4', 'NC_000005.10': '5', 'NC_000006.12': '6', 'NC_000007.14': '7', 'NC_000008.11': '8', 'NC_000009.12': '9', 'NC_000010.11': '10', 'NC_000011.10': '11', 'NC_000012.12': '12', 'NC_000013.11': '13', 'NC_000014.9': '14', 'NC_000015.10': '15', 'NC_000016.10': '16', 'NC_000017.11': '17', 'NC_000018.10': '18', 'NC_000019.10': '19', 'NC_000020.11': '20', 'NC_000021.9': '21', 'NC_000022.11': '22', 'NC_000023.11': 'X', 'NC_000024.10': 'Y', 'NC_012920.1': 'MT'}
vc_dbsnp_RE_RS = re.compile('^rs(\\d+)$', re.I)

def vc_dbsnp_normalize_rs_id(raw: str | None) -> str | None:
    if not raw:
        return None
    text = str(raw).strip()
    m = vc_dbsnp_RE_RS.match(text) or re.search('\\brs(\\d{4,})\\b', text, re.I)
    if not m:
        return None
    return f'rs{m.group(1)}'

def vc_dbsnp__cache_path(cache_dir: Path, rs_id: str) -> Path:
    return cache_dir / f'{rs_id.lower()}.json'

def vc_dbsnp__ssl_context() -> ssl.SSLContext:
    try:
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except Exception:
        return ssl.create_default_context()

def vc_dbsnp__fetch_refsnp_via_curl(url: str, timeout: float) -> dict[str, Any] | None:
    """Fallback when Python SSL cannot verify NCBI (common on some hosts)."""
    base = ['curl', '-fsS', '-L', '--max-time', str(max(5, int(timeout))), '-A', vc_dbsnp_USER_AGENT]
    for extra in ([], ['-k']):
        try:
            proc = subprocess.run([*base, *extra, url], capture_output=True, text=True, timeout=timeout + 10)
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return None
        if proc.returncode != 0 or not (proc.stdout or '').strip():
            continue
        try:
            data = json.loads(proc.stdout)
        except json.JSONDecodeError:
            continue
        if isinstance(data, dict):
            return data
    return None

def vc_dbsnp_fetch_refsnp_json(rs_id: str, *, cache_dir: Path | None, logger: logging.Logger | None=None, timeout: float=60.0, force: bool=False, retries: int=3) -> dict[str, Any] | None:
    rs_id = vc_dbsnp_normalize_rs_id(rs_id) or ''
    if not rs_id:
        return None
    rs_num = rs_id[2:]
    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        cached = vc_dbsnp__cache_path(cache_dir, rs_id)
        if cached.exists() and (not force):
            try:
                return json.loads(cached.read_text(encoding='utf-8'))
            except json.JSONDecodeError:
                if logger:
                    logger.warning('corrupt dbSNP cache %s; refetching', cached)
    url = vc_dbsnp_REFSNP_API.format(rs_num=rs_num)
    req = urllib.request.Request(url, headers={'User-Agent': vc_dbsnp_USER_AGENT})
    ctx = vc_dbsnp__ssl_context()
    payload: dict[str, Any] | None = None
    last_exc: Exception | None = None
    for attempt in range(1, max(1, retries) + 1):
        try:
            with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
                payload = json.load(resp)
            break
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, json.JSONDecodeError, http.client.IncompleteRead, ConnectionError, ssl.SSLError) as exc:
            last_exc = exc
            if logger:
                logger.warning('dbSNP fetch failed for %s (attempt %s/%s): %s', rs_id, attempt, retries, exc)
            if attempt < retries:
                time.sleep(0.5 * attempt)
    if payload is None:
        payload = vc_dbsnp__fetch_refsnp_via_curl(url, timeout)
        if payload is not None and logger:
            logger.warning('dbSNP fetch for %s succeeded via curl after urllib SSL/HTTP failure', rs_id)
    if payload is None:
        try:
            insecure = ssl._create_unverified_context()
            with urllib.request.urlopen(req, timeout=timeout, context=insecure) as resp:
                payload = json.load(resp)
            if logger:
                logger.warning('dbSNP fetch for %s used unverified SSL (NCBI RefSNP only)', rs_id)
        except Exception as exc:
            last_exc = exc
    if payload is None:
        if logger and last_exc is not None:
            logger.warning('dbSNP fetch gave up for %s: %s', rs_id, last_exc)
        return None
    if cache_dir is not None:
        vc_dbsnp__cache_path(cache_dir, rs_id).write_text(json.dumps(payload, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    return payload

def vc_dbsnp__alt_alleles_from_placement(placement: dict[str, Any]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for allele in placement.get('alleles') or []:
        spdi = allele.get('allele') or {} if isinstance(allele.get('allele'), dict) else {}
        if 'seq_id' not in spdi and isinstance(allele.get('allele'), dict):
            spdi = allele['allele'].get('spdi') or spdi
        if not isinstance(spdi, dict):
            continue
        deleted = str(spdi.get('deleted_sequence') or '')
        inserted = str(spdi.get('inserted_sequence') or '')
        if deleted == inserted:
            continue
        if spdi.get('position') is None:
            continue
        out.append({'seq_id': str(spdi.get('seq_id') or placement.get('seq_id') or ''), 'position_0based': int(spdi['position']), 'ref': deleted.upper(), 'alt': inserted.upper(), 'hgvs': allele.get('hgvs') or ''})
    return out

def vc_dbsnp__allele_frequency_scores(payload: dict[str, Any]) -> dict[tuple[str, str, str], float]:
    """Map (seq_id, ref, alt) → best available allele frequency (prefer GnomAD)."""
    scores: dict[tuple[str, str, str], float] = {}
    preferred = {'gnomad_genomes', 'gnomad_exomes', 'gnomad', '1000genomes', 'alfa'}
    for ann in (payload.get('primary_snapshot_data') or {}).get('allele_annotations') or []:
        for freq in ann.get('frequency') or []:
            obs = freq.get('observation') or {}
            seq_id = str(obs.get('seq_id') or '')
            deleted = str(obs.get('deleted_sequence') or '').upper()
            inserted = str(obs.get('inserted_sequence') or '').upper()
            if not seq_id or deleted == inserted:
                continue
            total = float(freq.get('total_count') or 0)
            count = float(freq.get('allele_count') or 0)
            if total <= 0:
                continue
            af = count / total
            study = str(freq.get('study_name') or '').lower()
            key = (seq_id, deleted, inserted)
            weight = 2.0 if any((name in study for name in preferred)) else 1.0
            scores[key] = max(scores.get(key, -1.0), af * weight)
    return scores

def vc_dbsnp_pick_grch38_allele(payload: dict[str, Any]) -> dict[str, Any] | None:
    """Return one GRCh38 genomic alt allele, or None if absent/ambiguous."""
    psd = payload.get('primary_snapshot_data') or {}
    placements = psd.get('placements_with_allele') or []
    candidates: list[dict[str, Any]] = []
    for placement in placements:
        seq_id = str(placement.get('seq_id') or '')
        chrom = vc_dbsnp_NC_TO_CHROM_GRCH38.get(seq_id)
        if not chrom:
            continue
        alts = vc_dbsnp__alt_alleles_from_placement(placement)
        for alt in alts:
            candidates.append({'chromosome': chrom, 'chromosome_accession': seq_id, 'position': alt['position_0based'] + 1, 'ref': alt['ref'], 'alt': alt['alt'], 'g_hgvs': alt.get('hgvs') or '', 'is_ptlp': bool(placement.get('is_ptlp'))})
    if not candidates:
        return None
    ptlp = [c for c in candidates if c['is_ptlp']]
    pool = ptlp or candidates
    uniq: dict[tuple[str, int, str, str], dict[str, Any]] = {}
    for c in pool:
        key = (c['chromosome'], int(c['position']), c['ref'], c['alt'])
        uniq[key] = c
    if len(uniq) == 1:
        return next(iter(uniq.values()))
    freq_scores = vc_dbsnp__allele_frequency_scores(payload)
    ranked: list[tuple[float, dict[str, Any]]] = []
    for c in uniq.values():
        score = freq_scores.get((c['chromosome_accession'], c['ref'], c['alt']), -1.0)
        ranked.append((score, c))
    ranked.sort(key=lambda item: item[0], reverse=True)
    if ranked and ranked[0][0] >= 0 and (len(ranked) == 1 or ranked[0][0] > ranked[1][0]):
        chosen = dict(ranked[0][1])
        chosen['selection_note'] = 'multi_allelic_chose_highest_frequency'
        return chosen
    clinvar_counts: dict[tuple[str, str, str], int] = {}
    for mov in payload.get('present_obs_movements') or []:
        allele = mov.get('allele_in_cur_release') or mov.get('observation') or {}
        seq_id = str(allele.get('seq_id') or '')
        deleted = str(allele.get('deleted_sequence') or '').upper()
        inserted = str(allele.get('inserted_sequence') or '').upper()
        if not seq_id or not deleted or deleted == inserted:
            continue
        weight = 0
        for cid in mov.get('component_ids') or []:
            ctype = str((cid or {}).get('type') or '').lower()
            if ctype == 'clinvar':
                weight += 3
            elif ctype == 'subsnp':
                weight += 1
            else:
                weight += 1
        if weight <= 0:
            weight = 1
        key = (seq_id, deleted, inserted)
        clinvar_counts[key] = clinvar_counts.get(key, 0) + weight
    clin_ranked: list[tuple[int, dict[str, Any]]] = []
    for c in uniq.values():
        score = clinvar_counts.get((c['chromosome_accession'], c['ref'], c['alt']), 0)
        clin_ranked.append((score, c))
    clin_ranked.sort(key=lambda item: item[0], reverse=True)
    if clin_ranked and clin_ranked[0][0] > 0 and (len(clin_ranked) == 1 or clin_ranked[0][0] > clin_ranked[1][0]):
        chosen = dict(clin_ranked[0][1])
        chosen['selection_note'] = 'multi_allelic_chose_most_clinvar_obs'
        return chosen
    return None

def vc_dbsnp_gene_loci_from_payload(payload: dict[str, Any]) -> list[dict[str, Any]]:
    """Unique gene loci annotated on a RefSNP payload (symbol + GeneID)."""
    seen: dict[str, dict[str, Any]] = {}
    for ann in (payload.get('primary_snapshot_data') or {}).get('allele_annotations') or []:
        for ag in ann.get('assembly_annotation') or []:
            for g in ag.get('genes') or []:
                locus = str(g.get('locus') or '').strip().upper()
                if not locus:
                    continue
                gid = g.get('id')
                seen[locus] = {'locus': locus, 'gene_id': int(gid) if gid is not None else None, 'name': g.get('name') or ''}
    return list(seen.values())

def vc_dbsnp_gene_symbol_for_gene_id(payload: dict[str, Any], gene_id: Any) -> str | None:
    if gene_id is None:
        return None
    want = str(gene_id).strip()
    for g in vc_dbsnp_gene_loci_from_payload(payload):
        if g.get('gene_id') is not None and str(g['gene_id']) == want:
            return g['locus']
    return None

class vc_dbsnp_DbSnpClient:
    """Small cache-backed client for unresolved rsID genomic fill-in."""

    def __init__(self, cache_dir: Path | None=None, *, logger: logging.Logger | None=None, sleep_s: float=0.34) -> None:
        self.cache_dir = cache_dir
        self.logger = logger
        self.sleep_s = sleep_s
        self._last_network_fetch = 0.0

    def _fetch(self, rs_id: str) -> dict[str, Any] | None:
        rs_norm = vc_dbsnp_normalize_rs_id(rs_id)
        if not rs_norm:
            return None
        cached_path = vc_dbsnp__cache_path(self.cache_dir, rs_norm) if self.cache_dir is not None else None
        need_network = cached_path is None or not cached_path.exists()
        if need_network and self.sleep_s > 0:
            elapsed = time.monotonic() - self._last_network_fetch
            if elapsed < self.sleep_s:
                time.sleep(self.sleep_s - elapsed)
        payload = vc_dbsnp_fetch_refsnp_json(rs_norm, cache_dir=self.cache_dir, logger=self.logger)
        if need_network:
            self._last_network_fetch = time.monotonic()
        return payload

    def lookup_gene(self, rs_id: str) -> str | None:
        """Return unique dbSNP gene locus for an rsID, or None if absent/ambiguous."""
        payload = self._fetch(rs_id)
        if not payload:
            return None
        loci = vc_dbsnp_gene_loci_from_payload(payload)
        if len(loci) != 1:
            return None
        return loci[0]['locus']

    def lookup_gene_for_id(self, rs_id: str, gene_id: Any) -> str | None:
        payload = self._fetch(rs_id)
        if not payload:
            return None
        return vc_dbsnp_gene_symbol_for_gene_id(payload, gene_id)

    def resolve_rs(self, rs_id: str) -> dict[str, Any] | None:
        rs_norm = vc_dbsnp_normalize_rs_id(rs_id)
        if not rs_norm:
            return None
        payload = self._fetch(rs_norm)
        if not payload:
            return None
        allele = vc_dbsnp_pick_grch38_allele(payload)
        if not allele:
            return None
        allele['rs_id'] = rs_norm
        allele['assembly'] = 'GRCh38'
        allele['source_url'] = vc_dbsnp_DBSNP_PAGE.format(rs_id=rs_norm)
        allele['allele_query_source'] = 'dbsnp'
        loci = vc_dbsnp_gene_loci_from_payload(payload)
        if len(loci) == 1:
            allele['gene'] = loci[0]['locus']
            allele['gene_id'] = loci[0].get('gene_id')
        return allele

def vc_dbsnp_default_dbsnp_cache(root: Path | None=None) -> Path:
    base = root or Path(os.environ.get('VARIANT_CLINICAL_ROOT', Path(__file__).resolve().parents[1]))
    return Path(base) / 'data/pubmed/cache/dbsnp'

# === export_genomic_allele_csv.py ===
"""
Export strict positive/negative literature alleles to CSV with genomic columns:

  chromosome, position (1-based), ref, alt

Plus assembly / chromosome-notation identifiers, and transcript mapping labels:

  - literal chr/g. / NC_:g. in variant text
  - rsID → ClinVar (prefer GRCh38 VCF alleles)
  - NM_:c. / gene+c. / gene+p. → ClinVar Name genomic rows (transcript fields)
  - PDF full text → literal g./chr:pos near variant mention
  - Ensembl GRCh38 GFF3 (MANE) + DNA FASTA → c. lift for ClinVar-absent alleles
  - rsID → NCBI dbSNP / Variation API (network fallback when upstream lifts fail);
    rows filled this way set allele_query_source=dbsnp

Default inputs: review_strict positives + strict other_benign negatives.

Outputs (under --out-dir):
  variants_genomic_all.csv / variants_genomic_ok.csv
  variants_genomic_positive.csv / variants_genomic_negative.csv
  variants_genomic_positive_ok.csv / variants_genomic_negative_ok.csv
  matching .jsonl sidecars
  genomic_allele_export_contract.json
"""
import argparse
import csv
import gzip
import json
import logging
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any
export_genomic_allele_csv_ROOT = Path(__file__).resolve().parents[1]
export_genomic_allele_csv_SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(export_genomic_allele_csv_SCRIPTS))
export_genomic_allele_csv_DEFAULT_POS = vc_paths_module.REVIEW_HIGHCONF_ANN
export_genomic_allele_csv_DEFAULT_POS_PATHO = vc_paths_module.PDF_PATHO
export_genomic_allele_csv_DEFAULT_NEG = vc_paths_module.NEGATIVES_FULL_PASS
export_genomic_allele_csv_DEFAULT_NEG_HARD = vc_paths_module.NEGATIVES_HARD
export_genomic_allele_csv_DEFAULT_JOURNAL_CONFIG = export_genomic_allele_csv_ROOT / 'configs/journal_priority.json'
export_genomic_allele_csv_DEFAULT_ARTICLES = vc_paths_module.MERGED_ARTICLES
export_genomic_allele_csv_DEFAULT_CLINVAR = vc_paths_module.CLINVAR_SUMMARY
export_genomic_allele_csv_DEFAULT_OUT_DIR = vc_paths_module.S9
export_genomic_allele_csv_DEFAULT_PDF_CACHE = vc_paths_module.PDF_TEXT_CACHE
export_genomic_allele_csv_DEFAULT_ENSEMBL_CACHE = export_genomic_allele_csv_ROOT / 'data/pubmed/cache/ensembl/ensembl_mane_transcript_index.json'
export_genomic_allele_csv_DEFAULT_PUBTATOR = vc_paths_module.PUBTATOR
export_genomic_allele_csv_GENOS_DB = export_genomic_allele_csv_ROOT.parent / 'data/database/ensembl'
export_genomic_allele_csv_CHROMOSOME_NOTATION = 'bare_ncbi'
export_genomic_allele_csv_ASSEMBLY_PREFER = 'GRCh38'
export_genomic_allele_csv_RE_RS = re.compile('\\brs(\\d{4,})\\b', re.I)
export_genomic_allele_csv_RE_NM = re.compile('\\b(NM_\\d+)(?:\\.(\\d+))?\\b', re.I)
export_genomic_allele_csv_RE_C = re.compile('\\b(c\\.(?:-?\\d+)(?:[+-]\\s*\\d+)?(?:_\\d+(?:[+-]\\s*\\d+)?)?(?:delins[ACGT]+|del|dup|ins[ACGT]+|[ACGT]\\s*>\\s*[ACGT]|[ACGT]+|[*=?]))', re.I)
export_genomic_allele_csv_RE_C_SUB = re.compile('\\bc\\.(?:-?\\d+)(?:[+-]\\d+)?([ACGT])\\s*>\\s*([ACGT])\\b', re.I)
export_genomic_allele_csv_RE_C_INCOMPLETE = re.compile('\\bc\\.(?:-?\\d+)(?:[+-]\\s*\\d+)?([ACGT])(?!\\s*>)', re.I)
export_genomic_allele_csv_RE_P = re.compile('\\b(p\\.\\(?[A-Za-z0-9_*?=]+\\)?)', re.I)
export_genomic_allele_csv_RE_ASSEMBLY_IN_TEXT = re.compile('\\b(GRCh3[78]|hg19|hg38)\\b', re.I)
export_genomic_allele_csv_RE_NM_FROM_CLINVAR_NAME = re.compile('\\b(NM_\\d+)(?:\\.(\\d+))?', re.I)
export_genomic_allele_csv_RE_EVIDENCE_NM_C = re.compile('\\b(NM_\\d+(?:\\.\\d+)?):c\\.([^\\s;,)]+)', re.I)
export_genomic_allele_csv_RE_EVIDENCE_C_FULL = re.compile('\\b(c\\.(?:-?\\d+|\\*\\d+)(?:[+-]\\s*\\d+)?(?:_\\d+(?:[+-]\\s*\\d+)?)?(?:delins|del|dup|ins|[ACGT]\\s*>\\s*[ACGT]|[ACGT]+))', re.I)
export_genomic_allele_csv_RE_EVIDENCE_C_SUB = re.compile('\\bc\\.(?:-?\\d+|\\*\\d+)(?:[+-]\\s*\\d+)?([ACGT])\\s*>\\s*([ACGT])', re.I)
export_genomic_allele_csv_RE_INCOMPLETE_NM_C = re.compile('\\b(NM_\\d+(?:\\.\\d+)?):c\\.(\\d+(?:[+-]\\s*\\d+)?)([ACGT]?)\\b', re.I)
export_genomic_allele_csv_RE_C_P_PAIR = re.compile('(c\\.(?:-?\\d+|\\*\\d+)(?:[+-]\\s*\\d+)?(?:_\\d+(?:[+-]\\s*\\d+)?)?(?:delins[ACGT]+|del|dup|ins[ACGT]+|\\s*[ACGT]\\s*>\\s*[ACGT]|\\s*[ACGT]+))\\s*[\\(（]\\s*(p\\.\\(?[A-Za-z0-9_*?=]+\\)?)', re.I)
export_genomic_allele_csv_RE_P_RS_PAIR = re.compile('(p\\.\\(?[A-Za-z0-9_*?=]+\\)?)\\s*[\\(（]\\s*(rs\\d{4,})\\s*[\\)）]', re.I)
export_genomic_allele_csv_RE_RS_P_PAIR = re.compile('(rs\\d{4,})\\s*[\\(（]\\s*(p\\.\\(?[A-Za-z0-9_*?=]+\\)?)', re.I)
export_genomic_allele_csv_RE_C_RS_PAIR = re.compile('(c\\.(?:-?\\d+|\\*\\d+)(?:[+-]\\s*\\d+)?(?:_\\d+(?:[+-]\\s*\\d+)?)?(?:delins[ACGT]+|del|dup|ins[ACGT]+|\\s*[ACGT]\\s*>\\s*[ACGT]|\\s*[ACGT]+))\\s*[\\(（]\\s*(rs\\d{4,})\\s*[\\)）]', re.I)
export_genomic_allele_csv_RE_RS_C_PAIR = re.compile('(rs\\d{4,})\\s*[\\(（]\\s*(c\\.(?:-?\\d+|\\*\\d+)(?:[+-]\\s*\\d+)?(?:_\\d+(?:[+-]\\s*\\d+)?)?(?:delins[ACGT]+|del|dup|ins[ACGT]+|\\s*[ACGT]\\s*>\\s*[ACGT]|\\s*[ACGT]+))', re.I)
export_genomic_allele_csv__AA3_TO1 = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Ter': '*', 'Stop': '*'}
export_genomic_allele_csv_CSV_PRIMARY_FIELDS = ['chromosome_ucsc', 'position', 'ref', 'alt', 'c_hgvs', 'nt_pro_only', 'p_hgvs', 'c_hgvs_strand', 'gene', 'genomic_region', 'genomic_region_genes', 'variant_text', 'chromosome_accession', 'assembly', 'pubmed_url']
export_genomic_allele_csv_CSV_REST_FIELDS = ['chromosome', 'chromosome_notation', 'assembly_source', 'assembly_confirmed', 'transcript_id', 'transcript_id_input', 'transcript_version', 'transcript_source', 'transcript_match', 'genomic_region_feature', 'genomic_region_transcripts', 'genomic_region_note', 'sample_class', 'label', 'screen_grade', 'link_tier', 'review_tier', 'pmid', 'locus_key_src', 'resolve_method', 'allele_query_source', 'resolve_status', 'resolve_confidence', 'resolve_note']
export_genomic_allele_csv_CSV_FIELDS = export_genomic_allele_csv_CSV_PRIMARY_FIELDS + export_genomic_allele_csv_CSV_REST_FIELDS

def export_genomic_allele_csv_setup_logger() -> logging.Logger:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    return logging.getLogger('export_genomic_allele_csv')

def export_genomic_allele_csv_norm_chrom(raw: str | None) -> str | None:
    if not raw:
        return None
    c = str(raw).strip()
    if c.lower().startswith('chr'):
        c = c[3:]
    c = c.upper()
    if c in {'M', 'MT', 'CHRMT'}:
        return 'MT'
    if re.fullmatch('[0-9]{1,2}|X|Y|MT', c):
        if c.isdigit() and (not 1 <= int(c) <= 22):
            return None
        return c
    return None

def export_genomic_allele_csv_chrom_ucsc(bare: str | None) -> str | None:
    if not bare:
        return None
    if bare == 'MT':
        return 'chrM'
    return f'chr{bare}'

def export_genomic_allele_csv_normalize_assembly(raw: str | None) -> str | None:
    if not raw:
        return None
    a = str(raw).strip()
    low = a.lower()
    if low in {'hg38', 'grch38'}:
        return 'GRCh38'
    if low in {'hg19', 'grch37'}:
        return 'GRCh37'
    if a in {'GRCh38', 'GRCh37'}:
        return a
    return a

def export_genomic_allele_csv_load_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not path.exists():
        return rows
    with path.open(encoding='utf-8') as f:
        for line in f:
            if line.strip():
                rows.append(json.loads(line))
    return rows

def export_genomic_allele_csv_normalize_c(text: str | None) -> str | None:
    if not text:
        return None
    t = validate_novelty_clinvar_module.norm_hgvs(text)
    if not t:
        return None
    t = re.sub('\\s+', '', t)
    if not t.lower().startswith('c.'):
        return None
    return 'c.' + t[2:]

def export_genomic_allele_csv_normalize_p_keys(text: str | None) -> set[str]:
    if not text:
        return set()
    raw = text.strip()
    body = re.sub('^p\\.\\((.+)\\)$', 'p.\\1', raw, flags=re.I)
    keys: set[str] = set()
    for cand in {body, raw, vc_text_module.normalize_p_hgvs_for_match(body)}:
        n = validate_novelty_clinvar_module.norm_hgvs(cand)
        if not n or not n.lower().startswith('p.'):
            continue
        keys |= validate_novelty_clinvar_module.expand_protein_keys(n)
        keys.add(vc_text_module.normalize_p_hgvs_for_match(n))
    extra: set[str] = set()
    for k in list(keys):
        m = re.fullmatch('(p\\.)([A-Za-z]{1,3})(\\d+)([A-Za-z]{1,3}|\\*)', k, flags=re.I)
        if m:
            a1 = export_genomic_allele_csv__aa_token_to1(m.group(2))
            b1 = export_genomic_allele_csv__aa_token_to1(m.group(4))
            if a1 and b1 and (a1 == b1) and (a1 != '*'):
                pos = m.group(3)
                extra.add(f'p.{a1}{pos}=')
                extra.add(f'p.{a1}{pos}{a1}')
                aa3 = next((t for t, v in export_genomic_allele_csv__AA3_TO1.items() if v == a1), None)
                if aa3:
                    extra.add(f'p.{aa3}{pos}=')
                    extra.add(f'p.{aa3}{pos}{aa3}')
            continue
        m_eq = re.fullmatch('(p\\.)([A-Za-z]{1,3})(\\d+)=', k, flags=re.I)
        if m_eq:
            a1 = export_genomic_allele_csv__aa_token_to1(m_eq.group(2))
            if a1 and a1 != '*':
                pos = m_eq.group(3)
                extra.add(f'p.{a1}{pos}{a1}')
                extra.add(f'p.{a1}{pos}=')
                aa3 = next((t for t, v in export_genomic_allele_csv__AA3_TO1.items() if v == a1), None)
                if aa3:
                    extra.add(f'p.{aa3}{pos}{aa3}')
                    extra.add(f'p.{aa3}{pos}=')
    keys |= extra
    return {k for k in keys if k}

def export_genomic_allele_csv_find_paired_c_for_protein(blob: str, p_hgvs: str) -> str | None:
    """Return c.HGVS that literature pairs with this p. (e.g. c.224C>A (p.Ala75Asp))."""
    targets = export_genomic_allele_csv_normalize_p_keys(p_hgvs)
    if not blob or not targets:
        return None
    for m in export_genomic_allele_csv_RE_C_P_PAIR.finditer(blob):
        c_full = export_genomic_allele_csv_normalize_c(m.group(1))
        p_keys = export_genomic_allele_csv_normalize_p_keys(m.group(2))
        if c_full and p_keys & targets:
            return c_full
    return None

def export_genomic_allele_csv_find_paired_p_for_cdna(blob: str, c_hgvs: str) -> str | None:
    """Return p.HGVS that literature pairs with this c."""
    c_norm = (export_genomic_allele_csv_normalize_c(c_hgvs) or '').lower().replace(' ', '')
    if not blob or not c_norm:
        return None
    for m in export_genomic_allele_csv_RE_C_P_PAIR.finditer(blob):
        c_full = export_genomic_allele_csv_normalize_c(m.group(1))
        if not c_full:
            continue
        if c_full.lower().replace(' ', '') != c_norm:
            continue
        p_raw = m.group(2)
        p_hgvs = re.sub('^p\\.\\((.+)\\)$', 'p.\\1', p_raw, flags=re.I)
        return validate_novelty_clinvar_module.norm_hgvs(p_hgvs)
    return None

def export_genomic_allele_csv__norm_rs_token(rs: str) -> str:
    m = export_genomic_allele_csv_RE_RS.search(rs or '')
    return f'rs{m.group(1)}' if m else ''

def export_genomic_allele_csv__norm_p_token(p_raw: str) -> str:
    p_hgvs = re.sub('^p\\.\\((.+)\\)$', 'p.\\1', (p_raw or '').strip(), flags=re.I)
    return validate_novelty_clinvar_module.norm_hgvs(p_hgvs)

def export_genomic_allele_csv_find_paired_p_for_rs(blob: str, rs: str) -> str | None:
    """Return p.HGVS literature-paired with this rs (e.g. p.S1437S (rs1131692092))."""
    rs_n = export_genomic_allele_csv__norm_rs_token(rs).lower()
    if not blob or not rs_n:
        return None
    for m in export_genomic_allele_csv_RE_P_RS_PAIR.finditer(blob):
        if export_genomic_allele_csv__norm_rs_token(m.group(2)).lower() == rs_n:
            return export_genomic_allele_csv__norm_p_token(m.group(1))
    for m in export_genomic_allele_csv_RE_RS_P_PAIR.finditer(blob):
        if export_genomic_allele_csv__norm_rs_token(m.group(1)).lower() == rs_n:
            return export_genomic_allele_csv__norm_p_token(m.group(2))
    return None

def export_genomic_allele_csv_find_paired_c_for_rs(blob: str, rs: str) -> str | None:
    """Return c.HGVS literature-paired with this rs."""
    rs_n = export_genomic_allele_csv__norm_rs_token(rs).lower()
    if not blob or not rs_n:
        return None
    for m in export_genomic_allele_csv_RE_C_RS_PAIR.finditer(blob):
        if export_genomic_allele_csv__norm_rs_token(m.group(2)).lower() == rs_n:
            return export_genomic_allele_csv_normalize_c(m.group(1))
    for m in export_genomic_allele_csv_RE_RS_C_PAIR.finditer(blob):
        if export_genomic_allele_csv__norm_rs_token(m.group(1)).lower() == rs_n:
            return export_genomic_allele_csv_normalize_c(m.group(2))
    return None

def export_genomic_allele_csv__aa_token_to1(tok: str) -> str | None:
    t = (tok or '').strip()
    if not t:
        return None
    if len(t) == 1:
        u = t.upper()
        return u if u.isalpha() or u == '*' else None
    if len(t) == 3:
        return export_genomic_allele_csv__AA3_TO1.get(t[0].upper() + t[1:].lower())
    low = t.lower()
    if low in ('ter', 'stop'):
        return '*'
    return None

def export_genomic_allele_csv_is_synonymous_p_hgvs(p_hgvs: str) -> bool:
    """True for silent protein HGVS: p.Pro1192= / p.P1192P / p.Pro1192Pro."""
    raw = (p_hgvs or '').strip()
    if not raw:
        return False
    body = re.sub('^p\\.\\(?', '', validate_novelty_clinvar_module.norm_hgvs(raw), flags=re.I).rstrip(')')
    if '=' in body:
        return True
    m = re.fullmatch('([A-Za-z]{1,3})(\\d+)([A-Za-z]{1,3}|\\*)', body)
    if not m:
        return False
    a1 = export_genomic_allele_csv__aa_token_to1(m.group(1))
    b1 = export_genomic_allele_csv__aa_token_to1(m.group(3))
    return bool(a1 and b1 and (a1 == b1) and (a1 != '*'))

def export_genomic_allele_csv_enrich_comps_from_blob(comps: dict[str, Any], blob: str, variant_text: str, *, allow_orphan_c: bool=True) -> dict[str, Any]:
    """Fill incomplete NM_/c. from evidence or PDF text (literature actual content).

    When allow_orphan_c is False (p.-only primary variants), never adopt an unrelated
    first c. mention from a multi-variant evidence sentence — that caused identical
    chr/pos/ref/alt across different protein alleles (e.g. ACAT1 p.Ala75Asp ← c.439G>T).
    """
    if not blob:
        return comps
    out = dict(comps)
    for m in export_genomic_allele_csv_RE_EVIDENCE_NM_C.finditer(blob):
        nm_raw = m.group(1)
        c_tail = m.group(2)
        nm_m = export_genomic_allele_csv_RE_NM.match(nm_raw)
        if nm_m and (not out.get('transcript_id_input')):
            out['transcript_id_input'] = nm_m.group(1).upper()
            out['transcript_version_input'] = nm_m.group(2)
        c_full = export_genomic_allele_csv_normalize_c('c.' + c_tail.split(')')[0])
        if c_full and allow_orphan_c and (not out.get('c_hgvs_input')):
            out['c_hgvs_input'] = c_full
    inc = export_genomic_allele_csv_RE_INCOMPLETE_NM_C.search(variant_text or '')
    if allow_orphan_c and (not inc):
        inc = export_genomic_allele_csv_RE_INCOMPLETE_NM_C.search(blob)
    if inc:
        nm_m = export_genomic_allele_csv_RE_NM.match(inc.group(1))
        if nm_m:
            out.setdefault('transcript_id_input', nm_m.group(1).upper())
            if nm_m.group(2):
                out.setdefault('transcript_version_input', nm_m.group(2))
        pos = inc.group(2)
        if allow_orphan_c and (not out.get('c_hgvs_input') or '>' not in (out.get('c_hgvs_input') or '')):
            sub_off = re.search(f'\\bc\\.{re.escape(pos)}\\s*([+-])\\s*(\\d+)\\s*([ACGT])\\s*>\\s*([ACGT])', blob, re.I)
            if sub_off:
                sign = sub_off.group(1)
                dist = sub_off.group(2)
                ref = sub_off.group(3)
                alt = sub_off.group(4)
                out['c_hgvs_input'] = export_genomic_allele_csv_normalize_c(f'c.{pos}{sign}{dist}{ref}>{alt}')
            else:
                sub = re.search(f'\\bc\\.{re.escape(pos)}\\s*([ACGT])\\s*>\\s*([ACGT])', blob, re.I)
                if sub:
                    ref = sub.group(1)
                    alt = sub.group(2)
                    out['c_hgvs_input'] = export_genomic_allele_csv_normalize_c(f'c.{pos}{ref}>{alt}')
                elif inc.group(3):
                    out['c_hgvs_input'] = export_genomic_allele_csv_normalize_c(f'c.{pos}{inc.group(3)}')
    if allow_orphan_c and (not out.get('c_hgvs_input') or len(out.get('c_hgvs_input') or '') < 6):
        for m in export_genomic_allele_csv_RE_EVIDENCE_C_FULL.finditer(blob):
            c_full = export_genomic_allele_csv_normalize_c(m.group(1))
            if c_full and (not variant_text or variant_text.lower() in c_full.lower() or c_full in variant_text):
                out['c_hgvs_input'] = c_full
                break
        if not out.get('c_hgvs_input'):
            for m in export_genomic_allele_csv_RE_EVIDENCE_C_FULL.finditer(blob):
                c_full = export_genomic_allele_csv_normalize_c(m.group(1))
                if c_full:
                    out['c_hgvs_input'] = c_full
                    break
    return out

def export_genomic_allele_csv_merge_comps_from_evidence(comps: dict[str, Any], *, variant_text: str, evidence_blob: str | None) -> dict[str, Any]:
    """Merge evidence into comps without cross-contaminating multi-allele sentences.

    Root cause of identical chr/pos/ref/alt with different p_hgvs:
    evidence lists several alleles (c.439G>T (p.Val147Leu), c.193A>T (p.Thr65Ser), …)
    and a p.-only row previously inherited the *first* orphan c. from the blob.

    Root cause of rs↔p cross-wire (e.g. rs1131692092 + p.P1192P):
    evidence lists p.P1192P (rs766447664) and p.S1437S (rs1131692092); an rs-primary
    row previously inherited the *first* orphan p. from the blob.
    """
    out = dict(comps)
    vt_level = export_genomic_allele_csv_infer_nt_pro_only(variant_text)
    blob = evidence_blob or ''
    if blob:
        ecomp = export_genomic_allele_csv_parse_variant_components(blob)
        for k in ('transcript_id_input', 'transcript_version_input'):
            if not out.get(k) and ecomp.get(k):
                out[k] = ecomp[k]
        if not out.get('rs_input') and ecomp.get('rs_input') and (vt_level != 'other'):
            if len(export_genomic_allele_csv_RE_RS.findall(blob)) == 1:
                out['rs_input'] = ecomp['rs_input']
        if vt_level == 'pro' and out.get('p_hgvs_input'):
            paired_c = export_genomic_allele_csv_find_paired_c_for_protein(blob, str(out['p_hgvs_input']))
            if paired_c:
                out['c_hgvs_input'] = paired_c
            out = export_genomic_allele_csv_enrich_comps_from_blob(out, blob, variant_text, allow_orphan_c=False)
        elif vt_level == 'nt' and out.get('c_hgvs_input'):
            if not out.get('p_hgvs_input'):
                paired_p = export_genomic_allele_csv_find_paired_p_for_cdna(blob, str(out['c_hgvs_input']))
                if paired_p:
                    out['p_hgvs_input'] = paired_p
            out = export_genomic_allele_csv_enrich_comps_from_blob(out, blob, variant_text, allow_orphan_c=True)
        elif vt_level == 'both':
            out = export_genomic_allele_csv_enrich_comps_from_blob(out, blob, variant_text, allow_orphan_c=False)
        else:
            rs = out.get('rs_input') or (f'rs{export_genomic_allele_csv_RE_RS.search(variant_text).group(1)}' if export_genomic_allele_csv_RE_RS.search(variant_text or '') else None)
            if rs:
                if not out.get('p_hgvs_input'):
                    paired_p = export_genomic_allele_csv_find_paired_p_for_rs(blob, str(rs))
                    if paired_p:
                        out['p_hgvs_input'] = paired_p
                if not out.get('c_hgvs_input'):
                    paired_c = export_genomic_allele_csv_find_paired_c_for_rs(blob, str(rs))
                    if paired_c:
                        out['c_hgvs_input'] = paired_c
                if out.get('p_hgvs_input') and (not out.get('c_hgvs_input')):
                    paired_c = export_genomic_allele_csv_find_paired_c_for_protein(blob, str(out['p_hgvs_input']))
                    if paired_c:
                        out['c_hgvs_input'] = paired_c
                if not out.get('c_hgvs_input'):
                    c_hits = [export_genomic_allele_csv_normalize_c(m.group(1)) for m in export_genomic_allele_csv_RE_EVIDENCE_C_FULL.finditer(blob)]
                    c_hits = [c for c in c_hits if c and '>' in c]
                    uniq = list(dict.fromkeys(c_hits))
                    if len(uniq) == 1:
                        out['c_hgvs_input'] = uniq[0]
                out = export_genomic_allele_csv_enrich_comps_from_blob(out, blob, variant_text, allow_orphan_c=False)
            else:
                for k in ('c_hgvs_input', 'p_hgvs_input', 'rs_input'):
                    if not out.get(k) and ecomp.get(k):
                        out[k] = ecomp[k]
                out = export_genomic_allele_csv_enrich_comps_from_blob(out, blob, variant_text, allow_orphan_c=True)
    return out

def export_genomic_allele_csv_parse_variant_components(variant_text: str) -> dict[str, Any]:
    v = variant_text or ''
    nm = export_genomic_allele_csv_RE_NM.search(v)
    c_m = export_genomic_allele_csv_RE_C.search(v)
    p_m = export_genomic_allele_csv_RE_P.search(v)
    rs_m = export_genomic_allele_csv_RE_RS.search(v)
    c_hgvs = export_genomic_allele_csv_normalize_c(c_m.group(1) if c_m else None)
    if not c_hgvs:
        inc = export_genomic_allele_csv_RE_C_INCOMPLETE.search(v)
        if inc:
            c_hgvs = export_genomic_allele_csv_normalize_c(inc.group(0))
    p_hgvs = None
    if p_m:
        p_raw = p_m.group(1)
        p_hgvs = re.sub('^p\\.\\((.+)\\)$', 'p.\\1', p_raw, flags=re.I)
        p_hgvs = validate_novelty_clinvar_module.norm_hgvs(p_hgvs)
    return {'transcript_id_input': nm.group(1).upper() if nm else None, 'transcript_version_input': nm.group(2) if nm and nm.group(2) else None, 'c_hgvs_input': c_hgvs, 'p_hgvs_input': p_hgvs, 'rs_input': f'rs{rs_m.group(1)}' if rs_m else None}

def export_genomic_allele_csv_infer_nt_pro_only(variant_text: str) -> str:
    """Classify literature variant_text by HGVS level present in the string.

    Returns:
      nt   — only c./NM_:c. nucleotide HGVS (no p.)
      pro  — only p. protein HGVS (no c.)
      both — both c. and p. in variant_text
      other — neither (rsID, g./chr literal, gene-only, etc.)
    """
    comps = export_genomic_allele_csv_parse_variant_components(variant_text or '')
    has_c = bool(comps.get('c_hgvs_input'))
    has_p = bool(comps.get('p_hgvs_input'))
    if has_c and has_p:
        return 'both'
    if has_c:
        return 'nt'
    if has_p:
        return 'pro'
    return 'other'

class export_genomic_allele_csv_ClinVarGenomicIndex:

    def __init__(self) -> None:
        self.by_rs: dict[str, list[dict[str, Any]]] = defaultdict(list)
        self.by_gene_c: dict[str, list[dict[str, Any]]] = defaultdict(list)
        self.by_gene_p: dict[str, list[dict[str, Any]]] = defaultdict(list)
        self.by_nm_c: dict[str, list[dict[str, Any]]] = defaultdict(list)
        self.n_rows = 0

def export_genomic_allele_csv__allele_ok(ref: str, alt: str) -> bool:
    if not ref or not alt:
        return False
    if ref in {'NA', 'N/A', '-'} or alt in {'NA', 'N/A', '-'}:
        return False
    return bool(re.fullmatch('[ACGTN]+', ref) and re.fullmatch('[ACGTN]+', alt))

def export_genomic_allele_csv_build_clinvar_genomic_index(clinvar_gz: Path, logger: logging.Logger) -> export_genomic_allele_csv_ClinVarGenomicIndex:
    idx = export_genomic_allele_csv_ClinVarGenomicIndex()
    if not clinvar_gz.exists():
        logger.warning('ClinVar missing: %s', clinvar_gz)
        return idx
    logger.info('Indexing ClinVar genomic+transcript from %s …', clinvar_gz)
    opener = gzip.open if str(clinvar_gz).endswith('.gz') else open
    with opener(clinvar_gz, 'rt', encoding='utf-8', errors='replace', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            idx.n_rows += 1
            if idx.n_rows % 1000000 == 0:
                logger.info('  scanned %s; rs=%s gene_c=%s gene_p=%s nm_c=%s', idx.n_rows, len(idx.by_rs), len(idx.by_gene_c), len(idx.by_gene_p), len(idx.by_nm_c))
            assembly = export_genomic_allele_csv_normalize_assembly(row.get('Assembly'))
            if assembly not in {'GRCh38', 'GRCh37'}:
                continue
            chrom = export_genomic_allele_csv_norm_chrom(row.get('Chromosome'))
            pos_raw = (row.get('PositionVCF') or row.get('Start') or '').strip()
            ref = (row.get('ReferenceAlleleVCF') or row.get('ReferenceAllele') or '').strip().upper()
            alt = (row.get('AlternateAlleleVCF') or row.get('AlternateAllele') or '').strip().upper()
            if not chrom or not pos_raw.isdigit() or (not export_genomic_allele_csv__allele_ok(ref, alt)):
                continue
            name = row.get('Name') or ''
            gene = (row.get('GeneSymbol') or '').strip().upper()
            g_from_name, c_hgvs, p_hgvs = validate_novelty_clinvar_module.extract_from_clinvar_name(name)
            if not gene and g_from_name:
                gene = g_from_name.upper()
            nm_m = export_genomic_allele_csv_RE_NM_FROM_CLINVAR_NAME.search(name)
            tx = nm_m.group(1).upper() if nm_m else None
            tx_ver = nm_m.group(2) if nm_m and nm_m.group(2) else None
            rs = (row.get('RS# (dbSNP)') or row.get('RS#') or '').strip()
            if rs in {'', '-1', '-', '.'}:
                rs = ''
            rec = {'assembly': assembly, 'chromosome': chrom, 'chromosome_accession': (row.get('ChromosomeAccession') or '').strip() or None, 'position': int(pos_raw), 'ref': ref, 'alt': alt, 'gene': gene or None, 'transcript_id': tx, 'transcript_version': tx_ver, 'c_hgvs': c_hgvs, 'p_hgvs': p_hgvs, 'rs': f'rs{rs}' if rs else None, 'name': name}
            if rs:
                idx.by_rs[f'rs{rs}'].append(rec)
            if gene and c_hgvs:
                idx.by_gene_c[f'{gene}|{c_hgvs.lower()}'].append(rec)
            if gene and p_hgvs:
                for pk in validate_novelty_clinvar_module.expand_protein_keys(p_hgvs) | {vc_text_module.normalize_p_hgvs_for_match(p_hgvs)}:
                    if pk:
                        idx.by_gene_p[f'{gene}|{pk.lower()}'].append(rec)
            if tx and c_hgvs:
                idx.by_nm_c[f'{tx}|{c_hgvs.lower()}'].append(rec)
                if tx_ver:
                    idx.by_nm_c[f'{tx}.{tx_ver}|{c_hgvs.lower()}'].append(rec)
    for store in (idx.by_rs, idx.by_gene_c, idx.by_gene_p, idx.by_nm_c):
        for k, hits in list(store.items()):
            pref = [h for h in hits if h['assembly'] == export_genomic_allele_csv_ASSEMBLY_PREFER]
            store[k] = pref or hits
    logger.info('ClinVar genomic index ready: rows=%s rs=%s gene_c=%s gene_p=%s nm_c=%s', idx.n_rows, len(idx.by_rs), len(idx.by_gene_c), len(idx.by_gene_p), len(idx.by_nm_c))
    return idx

def export_genomic_allele_csv_unique_genomic(hits: list[dict[str, Any]]) -> dict[str, Any] | None:
    if not hits:
        return None
    usable = [h for h in hits if (h.get('ref') or '') and (h.get('alt') or '') and (h.get('ref') != h.get('alt'))]
    if not usable:
        return None
    keys = {(h['assembly'], h['chromosome'], h['position'], h['ref'], h['alt']) for h in usable}
    if len(keys) == 1:
        return usable[0]
    pref = [h for h in usable if h['assembly'] == export_genomic_allele_csv_ASSEMBLY_PREFER]
    keys38 = {(h['chromosome'], h['position'], h['ref'], h['alt']) for h in pref}
    if len(keys38) == 1:
        return pref[0]
    return None

def export_genomic_allele_csv_pick_hit(hits: list[dict[str, Any]], *, gene: str | None=None, transcript_id: str | None=None, transcript_version: str | None=None) -> tuple[dict[str, Any] | None, str]:
    """Return (hit, transcript_match label)."""
    if not hits:
        return (None, 'none')
    pool = hits
    g = (gene or '').upper()
    if g:
        gene_hits = [h for h in pool if (h.get('gene') or '').upper() == g]
        if gene_hits:
            pool = gene_hits
    tx = (transcript_id or '').upper() if transcript_id else None
    if tx:
        exact = [h for h in pool if (h.get('transcript_id') or '').upper() == tx and (not transcript_version or str(h.get('transcript_version') or '') == str(transcript_version))]
        if exact:
            hit = export_genomic_allele_csv_unique_genomic(exact)
            if hit:
                return (hit, 'exact_nm_version' if transcript_version else 'exact_nm')
            return (None, 'ambiguous_nm')
        same_nm = [h for h in pool if (h.get('transcript_id') or '').upper() == tx]
        if same_nm:
            hit = export_genomic_allele_csv_unique_genomic(same_nm)
            if hit:
                return (hit, 'nm_version_mismatch' if transcript_version else 'exact_nm')
    hit = export_genomic_allele_csv_unique_genomic(pool)
    if hit:
        if tx and (hit.get('transcript_id') or '').upper() != tx:
            return (hit, 'clinvar_other_transcript')
        if tx:
            return (hit, 'exact_nm')
        return (hit, 'clinvar_preferred_transcript' if hit.get('transcript_id') else 'no_transcript_in_clinvar')
    return (None, 'ambiguous')

def export_genomic_allele_csv_attach_tx_fields(out: dict[str, Any], comps: dict[str, Any], hit: dict[str, Any] | None=None, *, transcript_match: str='none', transcript_source: str='') -> dict[str, Any]:
    out['transcript_id_input'] = comps.get('transcript_id_input')
    out['c_hgvs'] = (hit or {}).get('c_hgvs') or comps.get('c_hgvs_input')
    out['p_hgvs'] = (hit or {}).get('p_hgvs') or comps.get('p_hgvs_input')
    if hit and hit.get('transcript_id'):
        out['transcript_id'] = hit.get('transcript_id')
        out['transcript_version'] = hit.get('transcript_version')
        out['transcript_source'] = transcript_source or 'clinvar_name'
    elif comps.get('transcript_id_input'):
        out['transcript_id'] = comps.get('transcript_id_input')
        out['transcript_version'] = comps.get('transcript_version_input')
        out['transcript_source'] = transcript_source or 'literature_variant_text'
    else:
        out['transcript_id'] = None
        out['transcript_version'] = None
        out['transcript_source'] = transcript_source or 'none'
    out['transcript_match'] = transcript_match
    return out

def export_genomic_allele_csv_try_literal_genomic(variant_text: str, comps: dict[str, Any]) -> dict[str, Any] | None:
    text = variant_text or ''
    spans = vc_intergenic_find_intergenic_variant_spans(text)
    for sp in spans:
        if sp['subtype'] in {'genomic_chrpos_change', 'hgvs_g_nc', 'hgvs_g'}:
            loc = sp['locus']
            chrom = export_genomic_allele_csv_norm_chrom(loc.chrom)
            ass = export_genomic_allele_csv_normalize_assembly(loc.assembly) or export_genomic_allele_csv_normalize_assembly(vc_intergenic_detect_assembly(text))
            if chrom and loc.pos and loc.ref and loc.alt:
                out = {'chromosome': chrom, 'position': int(loc.pos), 'ref': loc.ref.upper(), 'alt': loc.alt.upper(), 'assembly': ass or export_genomic_allele_csv_ASSEMBLY_PREFER, 'assembly_source': 'literal_text' if ass else 'pipeline_default_target', 'assembly_confirmed': bool(ass), 'chromosome_accession': None, 'resolve_method': f"literal_{sp['subtype']}", 'resolve_confidence': 'high' if ass else 'medium', 'resolve_note': 'parsed from variant_text genomic notation' + ('' if ass else f'; assembly defaulted to {export_genomic_allele_csv_ASSEMBLY_PREFER}'), 'resolve_status': 'ok'}
                return export_genomic_allele_csv_attach_tx_fields(out, comps, transcript_match='not_applicable', transcript_source='none')
    for sp in spans:
        if sp['subtype'] == 'genomic_chrpos':
            loc = sp['locus']
            chrom = export_genomic_allele_csv_norm_chrom(loc.chrom)
            ass = export_genomic_allele_csv_normalize_assembly(loc.assembly) or export_genomic_allele_csv_normalize_assembly(vc_intergenic_detect_assembly(text))
            if chrom and loc.pos:
                pos_out = int(loc.pos) if loc.ref and loc.alt else None
                out = {'chromosome': chrom, 'position': pos_out, 'ref': loc.ref.upper() if loc.ref else '', 'alt': loc.alt.upper() if loc.alt else '', 'assembly': ass or export_genomic_allele_csv_ASSEMBLY_PREFER, 'assembly_source': 'literal_text' if ass else 'pipeline_default_target', 'assembly_confirmed': bool(ass), 'chromosome_accession': None, 'resolve_method': 'literal_genomic_chrpos', 'resolve_confidence': 'low', 'resolve_note': 'chr:pos without ref/alt', 'resolve_status': 'partial'}
                return export_genomic_allele_csv_attach_tx_fields(out, comps, transcript_match='not_applicable', transcript_source='none')
    return None

def export_genomic_allele_csv_hit_to_resolved(hit: dict[str, Any], comps: dict[str, Any], *, method: str, transcript_match: str) -> dict[str, Any]:
    out = {'chromosome': hit['chromosome'], 'position': hit['position'], 'ref': hit['ref'], 'alt': hit['alt'], 'assembly': hit['assembly'], 'assembly_source': 'clinvar_variant_summary', 'assembly_confirmed': True, 'chromosome_accession': hit.get('chromosome_accession'), 'resolve_method': method, 'resolve_confidence': 'high' if hit['assembly'] == export_genomic_allele_csv_ASSEMBLY_PREFER else 'medium', 'resolve_note': f'{method}; transcript_match={transcript_match}', 'resolve_status': 'ok'}
    return export_genomic_allele_csv_attach_tx_fields(out, comps, hit, transcript_match=transcript_match, transcript_source='clinvar_name')

def export_genomic_allele_csv_try_clinvar_transcript(*, variant_text: str, gene: str | None, comps: dict[str, Any], idx: export_genomic_allele_csv_ClinVarGenomicIndex) -> dict[str, Any] | None:
    g = (gene or '').upper() or None
    c = comps.get('c_hgvs_input')
    p = comps.get('p_hgvs_input')
    tx = comps.get('transcript_id_input')
    tx_ver = comps.get('transcript_version_input')
    rs = comps.get('rs_input')
    rs_hit_pool: list[dict[str, Any]] | None = None

    def _intersect_rs_pool(pool: list[dict[str, Any]]) -> list[dict[str, Any]]:
        if not rs_hit_pool:
            return pool
        keys = {(h.get('chromosome'), h.get('position'), h.get('ref'), h.get('alt')) for h in rs_hit_pool}
        clipped = [h for h in pool if (h.get('chromosome'), h.get('position'), h.get('ref'), h.get('alt')) in keys]
        return clipped or pool
    if rs:
        hits = list(idx.by_rs.get(rs) or [])
        if len(hits) > 1 and c:
            c_l = str(c).lower().replace(' ', '')
            filtered = [h for h in hits if (h.get('c_hgvs') or '').lower().replace(' ', '') == c_l]
            if filtered:
                hits = filtered
        if len(hits) > 1 and p:
            p_keys = export_genomic_allele_csv_normalize_p_keys(str(p))
            filtered = []
            for h in hits:
                hp = h.get('p_hgvs') or ''
                if hp and export_genomic_allele_csv_normalize_p_keys(hp) & p_keys:
                    filtered.append(h)
            if filtered:
                hits = filtered
        hit, tmatch = export_genomic_allele_csv_pick_hit(hits, gene=g, transcript_id=tx, transcript_version=tx_ver)
        if hit:
            return export_genomic_allele_csv_hit_to_resolved(hit, comps, method='rs_clinvar', transcript_match=tmatch)
        if hits:
            rs_hit_pool = hits
    if tx and c:
        keys = [f'{tx}|{c.lower()}']
        if tx_ver:
            keys.insert(0, f'{tx}.{tx_ver}|{c.lower()}')
        hits = []
        for k in keys:
            hits.extend(idx.by_nm_c.get(k) or [])
        hits = _intersect_rs_pool(hits)
        hit, tmatch = export_genomic_allele_csv_pick_hit(hits, gene=g, transcript_id=tx, transcript_version=tx_ver)
        if hit:
            return export_genomic_allele_csv_hit_to_resolved(hit, comps, method='nm_c_clinvar', transcript_match=tmatch)
    if g and c:
        hits = list(idx.by_gene_c.get(f'{g}|{c.lower()}') or [])
        hits = _intersect_rs_pool(hits)
        hit, tmatch = export_genomic_allele_csv_pick_hit(hits, gene=g, transcript_id=tx, transcript_version=tx_ver)
        if hit:
            return export_genomic_allele_csv_hit_to_resolved(hit, comps, method='gene_c_clinvar', transcript_match=tmatch)
    if g and p:
        hits = []
        for pk in export_genomic_allele_csv_normalize_p_keys(p):
            hits.extend(idx.by_gene_p.get(f'{g}|{pk.lower()}') or [])
        hits = _intersect_rs_pool(hits)
        seen: set[tuple] = set()
        uniq_hits: list[dict[str, Any]] = []
        for h in hits:
            key = (h['assembly'], h['chromosome'], h['position'], h['ref'], h['alt'], h.get('transcript_id'))
            if key in seen:
                continue
            seen.add(key)
            uniq_hits.append(h)
        hit, tmatch = export_genomic_allele_csv_pick_hit(uniq_hits, gene=g, transcript_id=tx, transcript_version=tx_ver)
        if hit:
            return export_genomic_allele_csv_hit_to_resolved(hit, comps, method='gene_p_clinvar', transcript_match=tmatch)
        if uniq_hits:
            out = {'resolve_method': 'gene_p_clinvar', 'resolve_status': 'failed', 'resolve_confidence': 'low', 'resolve_note': f'gene+p matched {len(uniq_hits)} ClinVar genomic alleles (ambiguous)', 'assembly': export_genomic_allele_csv_ASSEMBLY_PREFER, 'assembly_source': 'clinvar_preferred_target', 'assembly_confirmed': False, 'chromosome': None, 'position': None, 'ref': '', 'alt': '', 'chromosome_accession': None, 'failed': True}
            txs = {(h.get('transcript_id'), h.get('transcript_version')) for h in uniq_hits if h.get('transcript_id')}
            tip = uniq_hits[0] if len(txs) == 1 else None
            return export_genomic_allele_csv_attach_tx_fields(out, comps, tip, transcript_match='ambiguous', transcript_source='clinvar_name' if tip else 'literature_variant_text' if tx else 'none')
    return None

def export_genomic_allele_csv_try_pdf_genomic(*, pmid: str | None, gene: str | None, variant_text: str, comps: dict[str, Any], pdf_cache: Path) -> dict[str, Any] | None:
    if not pmid:
        return None
    pdf_text = vc_pdf_genomic_load_pdf_text(str(pmid), pdf_cache)
    if not pdf_text:
        return None
    blob = export_genomic_allele_csv_enrich_comps_from_blob(comps, pdf_text, variant_text)
    hits = vc_pdf_genomic_extract_pdf_genomic_hits(pdf_text=pdf_text, gene=gene, variant_text=variant_text, c_hgvs=blob.get('c_hgvs_input') or comps.get('c_hgvs_input'), p_hgvs=blob.get('p_hgvs_input') or comps.get('p_hgvs_input'))
    hit = vc_pdf_genomic_pick_pdf_genomic_hit(hits, gene=gene, variant_text=variant_text)
    if not hit:
        return None
    ass = export_genomic_allele_csv_normalize_assembly(hit.get('assembly')) or export_genomic_allele_csv_ASSEMBLY_PREFER
    ref = (hit.get('ref') or '').upper()
    alt = (hit.get('alt') or '').upper()
    chrom = export_genomic_allele_csv_norm_chrom(hit.get('chromosome'))
    pos = hit.get('position')
    if not chrom or not pos:
        return None
    pos_out = int(pos) if ref or alt else None
    out = {'chromosome': chrom, 'position': pos_out, 'ref': ref, 'alt': alt, 'assembly': ass, 'assembly_source': 'pdf_fulltext', 'assembly_confirmed': bool(export_genomic_allele_csv_normalize_assembly(hit.get('assembly'))), 'chromosome_accession': None, 'resolve_method': 'pdf_literal_genomic', 'resolve_confidence': 'high' if ref and alt else 'medium', 'resolve_note': f"from PDF near variant; raw={hit.get('raw', '')[:80]}", 'resolve_status': 'ok' if ref and alt else 'partial'}
    return export_genomic_allele_csv_attach_tx_fields(out, blob, transcript_match='not_applicable' if not blob.get('transcript_id_input') else 'input_only', transcript_source='literature_pdf' if blob.get('transcript_id_input') else 'none')

def export_genomic_allele_csv_try_ensembl_c_map(*, gene: str | None, comps: dict[str, Any], ensembl_idx: vc_ensembl_c_map_module.EnsemblManeIndex, fasta: vc_ensembl_c_map_module.FastaReader) -> dict[str, Any] | None:
    c_hgvs = comps.get('c_hgvs_input')
    if not gene or not c_hgvs:
        return None
    tx = ensembl_idx.get(gene)
    if not tx:
        return None
    try:
        mapped = vc_ensembl_c_map_module.map_c_to_genomic(tx=tx, fasta=fasta, c_hgvs=c_hgvs)
    except (KeyError, ValueError) as exc:
        return {'resolve_method': 'ensembl_mane_c_map', 'resolve_status': 'failed', 'resolve_confidence': 'low', 'resolve_note': f'Ensembl MANE c. lift failed: {exc}', 'assembly': export_genomic_allele_csv_ASSEMBLY_PREFER, 'assembly_source': 'ensembl_gff3_grch38', 'assembly_confirmed': True, 'chromosome': None, 'position': None, 'ref': '', 'alt': '', 'chromosome_accession': None, 'failed': True}
    out = {'chromosome': mapped['chromosome'], 'position': mapped['position'], 'ref': mapped['ref'], 'alt': mapped['alt'], 'assembly': export_genomic_allele_csv_ASSEMBLY_PREFER, 'assembly_source': 'ensembl_gff3_grch38', 'assembly_confirmed': True, 'chromosome_accession': None, 'resolve_method': 'ensembl_mane_c_map', 'resolve_confidence': 'medium', 'resolve_note': f'MANE {tx.transcript_id}.{tx.version} c. lift from {c_hgvs}', 'resolve_status': 'ok' if mapped.get('ref') is not None else 'partial'}
    comps_out = dict(comps)
    comps_out['transcript_id_input'] = comps_out.get('transcript_id_input') or tx.transcript_id
    comps_out['transcript_version_input'] = comps_out.get('transcript_version_input') or tx.version
    return export_genomic_allele_csv_attach_tx_fields(out, comps_out, {'transcript_id': tx.transcript_id, 'transcript_version': tx.version, 'c_hgvs': c_hgvs, 'p_hgvs': comps.get('p_hgvs_input')}, transcript_match='ensembl_mane_select', transcript_source='ensembl_mane')

def export_genomic_allele_csv_try_ensembl_p_map(*, gene: str | None, comps: dict[str, Any], ensembl_idx: vc_ensembl_c_map_module.EnsemblManeIndex, fasta: vc_ensembl_c_map_module.FastaReader) -> dict[str, Any] | None:
    p_hgvs = comps.get('p_hgvs_input')
    if not gene or not p_hgvs:
        return None
    if export_genomic_allele_csv_is_synonymous_p_hgvs(str(p_hgvs)):
        return {'resolve_method': 'ensembl_mane_p_map', 'resolve_status': 'failed', 'resolve_confidence': 'low', 'resolve_note': f'synonymous p.HGVS skipped for MANE p. lift: {p_hgvs}; prefer c./rs/ClinVar', 'assembly': export_genomic_allele_csv_ASSEMBLY_PREFER, 'assembly_source': 'ensembl_gff3_grch38', 'assembly_confirmed': True, 'chromosome': None, 'position': None, 'ref': '', 'alt': '', 'chromosome_accession': None, 'failed': True}
    tx = ensembl_idx.get(gene)
    if not tx:
        return None
    try:
        mapped = vc_ensembl_c_map_module.map_p_to_genomic(tx=tx, fasta=fasta, p_hgvs=p_hgvs)
    except ValueError as exc:
        return {'resolve_method': 'ensembl_mane_p_map', 'resolve_status': 'failed', 'resolve_confidence': 'low', 'resolve_note': f'Ensembl MANE p. lift failed: {exc}', 'assembly': export_genomic_allele_csv_ASSEMBLY_PREFER, 'assembly_source': 'ensembl_gff3_grch38', 'assembly_confirmed': True, 'chromosome': None, 'position': None, 'ref': '', 'alt': '', 'chromosome_accession': None, 'failed': True}
    ref = (mapped.get('ref') or '').upper()
    alt = (mapped.get('alt') or '').upper()
    if ref and alt and (ref == alt):
        return {'resolve_method': 'ensembl_mane_p_map', 'resolve_status': 'failed', 'resolve_confidence': 'low', 'resolve_note': f'MANE p. lift produced ref==alt ({ref}); rejected for {p_hgvs}', 'assembly': export_genomic_allele_csv_ASSEMBLY_PREFER, 'assembly_source': 'ensembl_gff3_grch38', 'assembly_confirmed': True, 'chromosome': None, 'position': None, 'ref': '', 'alt': '', 'chromosome_accession': None, 'failed': True}
    out = {'chromosome': mapped['chromosome'], 'position': mapped['position'], 'ref': mapped['ref'], 'alt': mapped['alt'], 'assembly': export_genomic_allele_csv_ASSEMBLY_PREFER, 'assembly_source': 'ensembl_gff3_grch38', 'assembly_confirmed': True, 'chromosome_accession': None, 'resolve_method': 'ensembl_mane_p_map', 'resolve_confidence': 'medium', 'resolve_note': f'MANE {tx.transcript_id}.{tx.version} p. lift from {p_hgvs}', 'resolve_status': 'ok'}
    comps_out = dict(comps)
    comps_out['transcript_id_input'] = comps_out.get('transcript_id_input') or tx.transcript_id
    comps_out['transcript_version_input'] = comps_out.get('transcript_version_input') or tx.version
    return export_genomic_allele_csv_attach_tx_fields(out, comps_out, {'transcript_id': tx.transcript_id, 'transcript_version': tx.version, 'c_hgvs': mapped.get('c_hgvs'), 'p_hgvs': p_hgvs}, transcript_match='ensembl_mane_select', transcript_source='ensembl_mane')

def export_genomic_allele_csv_empty_resolve(comps: dict[str, Any], note: str, method: str='unresolved') -> dict[str, Any]:
    out = {'chromosome': None, 'position': None, 'ref': '', 'alt': '', 'assembly': export_genomic_allele_csv_ASSEMBLY_PREFER, 'assembly_source': 'pipeline_default_target', 'assembly_confirmed': False, 'chromosome_accession': None, 'resolve_method': method, 'allele_query_source': '', 'resolve_confidence': 'low', 'resolve_note': note, 'resolve_status': 'failed', 'failed': True}
    tmatch = 'input_only' if comps.get('transcript_id_input') else 'none'
    return export_genomic_allele_csv_attach_tx_fields(out, comps, transcript_match=tmatch, transcript_source='literature_variant_text' if comps.get('transcript_id_input') else 'none')

def export_genomic_allele_csv_try_dbsnp_rs(*, comps: dict[str, Any], dbsnp: vc_dbsnp_DbSnpClient | None) -> dict[str, Any] | None:
    """Network fallback for rsIDs that ClinVar/PDF/Ensembl could not uniquely map."""
    if dbsnp is None:
        return None
    rs = comps.get('rs_input')
    if not rs:
        return None
    hit = dbsnp.resolve_rs(str(rs))
    if not hit:
        return None
    out = {'chromosome': hit['chromosome'], 'position': hit['position'], 'ref': hit['ref'], 'alt': hit['alt'], 'assembly': hit.get('assembly') or export_genomic_allele_csv_ASSEMBLY_PREFER, 'assembly_source': 'dbsnp_refsnp_api', 'assembly_confirmed': True, 'chromosome_accession': hit.get('chromosome_accession') or '', 'resolve_method': 'rs_dbsnp', 'allele_query_source': 'dbsnp', 'resolve_status': 'ok', 'resolve_confidence': 'high', 'resolve_note': f"dbSNP RefSNP API {hit.get('rs_id')}; {hit.get('source_url')}" + (f"; g_hgvs={hit['g_hgvs']}" if hit.get('g_hgvs') else '') + (f"; {hit['selection_note']}" if hit.get('selection_note') else '')}
    return export_genomic_allele_csv_attach_tx_fields(out, comps, transcript_match='not_applicable', transcript_source='none')

def export_genomic_allele_csv_resolve_allele(*, variant_text: str, gene: str | None, evidence_blob: str | None, idx: export_genomic_allele_csv_ClinVarGenomicIndex, pmid: str | None=None, pdf_cache: Path | None=None, ensembl_idx: vc_ensembl_c_map_module.EnsemblManeIndex | None=None, fasta: vc_ensembl_c_map_module.FastaReader | None=None, dbsnp: vc_dbsnp_DbSnpClient | None=None) -> dict[str, Any]:

    def attach_c_hgvs_strand(resolved: dict[str, Any]) -> None:
        resolved['c_hgvs_strand'] = ''
        if not ensembl_idx:
            return
        if not (resolved.get('c_hgvs') or comps.get('c_hgvs_input')):
            return
        if not gene:
            return
        tx = ensembl_idx.get(str(gene).upper())
        if not tx:
            return
        resolved['c_hgvs_strand'] = tx.strand or ''

    def reconfirm_ref_from_fasta(resolved: dict[str, Any]) -> None:
        if not fasta:
            return
        if resolved.get('assembly') != export_genomic_allele_csv_ASSEMBLY_PREFER:
            return
        chrom = resolved.get('chromosome')
        pos = resolved.get('position')
        ref = (resolved.get('ref') or '').upper()
        if not chrom or pos is None or (not ref):
            return
        try:
            end = int(pos) + len(ref) - 1
            fetched = fasta.fetch(str(chrom), int(pos), end)
        except Exception:
            return
        if fetched.upper() != ref:
            resolved['position'] = None
            resolved['ref'] = ''
            resolved['alt'] = ''
            resolved['resolve_status'] = 'partial'
            resolved['resolve_confidence'] = 'low'
            resolved['resolve_note'] = (resolved.get('resolve_note') or '') + '; ref mismatch with FASTA'
    blob = ' '.join((x for x in [variant_text, evidence_blob or ''] if x))
    comps = export_genomic_allele_csv_parse_variant_components(variant_text)
    comps = export_genomic_allele_csv_merge_comps_from_evidence(comps, variant_text=variant_text, evidence_blob=evidence_blob)
    vt_level = export_genomic_allele_csv_infer_nt_pro_only(variant_text)
    lit = export_genomic_allele_csv_try_literal_genomic(variant_text, comps)
    if not lit and vt_level != 'pro':
        lit = export_genomic_allele_csv_try_literal_genomic(blob, comps)
    lit_partial = None
    if lit and lit.get('resolve_status') == 'ok':
        attach_c_hgvs_strand(lit)
        reconfirm_ref_from_fasta(lit)
        return lit
    if lit and lit.get('resolve_status') == 'partial':
        has_better = bool(comps.get('c_hgvs_input') or comps.get('rs_input') or (comps.get('p_hgvs_input') and gene) or export_genomic_allele_csv_RE_RS.search(variant_text or '') or export_genomic_allele_csv_RE_C.search(variant_text or ''))
        if has_better:
            lit_partial = lit
        else:
            attach_c_hgvs_strand(lit)
            reconfirm_ref_from_fasta(lit)
            return lit
    cv = export_genomic_allele_csv_try_clinvar_transcript(variant_text=variant_text, gene=gene, comps=comps, idx=idx)
    if cv and cv.get('resolve_status') == 'ok':
        attach_c_hgvs_strand(cv)
        reconfirm_ref_from_fasta(cv)
        return cv
    if pdf_cache:
        pdf_hit = export_genomic_allele_csv_try_pdf_genomic(pmid=pmid, gene=gene, variant_text=variant_text, comps=comps, pdf_cache=pdf_cache)
        if pdf_hit and pdf_hit.get('resolve_status') == 'ok':
            attach_c_hgvs_strand(pdf_hit)
            reconfirm_ref_from_fasta(pdf_hit)
            return pdf_hit
    if ensembl_idx and fasta and comps.get('p_hgvs_input') and (not comps.get('c_hgvs_input')):
        ens_p = export_genomic_allele_csv_try_ensembl_p_map(gene=gene, comps=comps, ensembl_idx=ensembl_idx, fasta=fasta)
        if ens_p and ens_p.get('resolve_status') == 'ok':
            attach_c_hgvs_strand(ens_p)
            reconfirm_ref_from_fasta(ens_p)
            return ens_p
    if ensembl_idx and fasta and comps.get('c_hgvs_input'):
        ens = export_genomic_allele_csv_try_ensembl_c_map(gene=gene, comps=comps, ensembl_idx=ensembl_idx, fasta=fasta)
        if ens and ens.get('resolve_status') == 'ok':
            attach_c_hgvs_strand(ens)
            reconfirm_ref_from_fasta(ens)
            return ens
    dbsnp_hit = export_genomic_allele_csv_try_dbsnp_rs(comps=comps, dbsnp=dbsnp)
    if dbsnp_hit and dbsnp_hit.get('resolve_status') == 'ok':
        attach_c_hgvs_strand(dbsnp_hit)
        reconfirm_ref_from_fasta(dbsnp_hit)
        if dbsnp_hit.get('resolve_status') == 'ok':
            return dbsnp_hit
    if cv:
        attach_c_hgvs_strand(cv)
        reconfirm_ref_from_fasta(cv)
        cv.setdefault('allele_query_source', '')
        return cv
    if pdf_cache:
        pdf_hit = export_genomic_allele_csv_try_pdf_genomic(pmid=pmid, gene=gene, variant_text=variant_text, comps=comps, pdf_cache=pdf_cache)
        if pdf_hit:
            attach_c_hgvs_strand(pdf_hit)
            reconfirm_ref_from_fasta(pdf_hit)
            pdf_hit.setdefault('allele_query_source', '')
            return pdf_hit
    if ensembl_idx and fasta and comps.get('p_hgvs_input') and (not comps.get('c_hgvs_input')):
        ens_p = export_genomic_allele_csv_try_ensembl_p_map(gene=gene, comps=comps, ensembl_idx=ensembl_idx, fasta=fasta)
        if ens_p:
            if ens_p.get('resolve_status') != 'ok' and comps.get('rs_input') and ('synonymous' in (ens_p.get('resolve_note') or '').lower()):
                return export_genomic_allele_csv_empty_resolve(comps, f"{comps['rs_input']} not uniquely resolved by ClinVar/dbSNP; synonymous {comps.get('p_hgvs_input')} not liftable via MANE p-map", method='rs_unresolved')
            attach_c_hgvs_strand(ens_p)
            reconfirm_ref_from_fasta(ens_p)
            ens_p.setdefault('allele_query_source', '')
            return ens_p
    if ensembl_idx and fasta and comps.get('c_hgvs_input'):
        ens = export_genomic_allele_csv_try_ensembl_c_map(gene=gene, comps=comps, ensembl_idx=ensembl_idx, fasta=fasta)
        if ens:
            attach_c_hgvs_strand(ens)
            reconfirm_ref_from_fasta(ens)
            ens.setdefault('allele_query_source', '')
            return ens
    m = export_genomic_allele_csv_RE_C_SUB.search(variant_text or '')
    if not m and comps.get('c_hgvs_input'):
        m = export_genomic_allele_csv_RE_C_SUB.search(str(comps['c_hgvs_input']))
    if not m and vt_level not in {'pro'}:
        m = export_genomic_allele_csv_RE_C_SUB.search(blob or '') or export_genomic_allele_csv_RE_EVIDENCE_C_SUB.search(blob or '')
    if m:
        out = {'chromosome': None, 'position': None, 'ref': m.group(1).upper(), 'alt': m.group(2).upper(), 'assembly': export_genomic_allele_csv_ASSEMBLY_PREFER, 'assembly_source': 'pipeline_default_target', 'assembly_confirmed': False, 'chromosome_accession': None, 'resolve_method': 'c_hgvs_alleles_only', 'allele_query_source': '', 'resolve_status': 'partial', 'resolve_confidence': 'low', 'resolve_note': 'ref/alt from c.HGVS; no ClinVar genomic hit (often ClinVar-absent positives)'}
        out = export_genomic_allele_csv_attach_tx_fields(out, comps, transcript_match='input_only' if comps.get('transcript_id_input') else 'none', transcript_source='literature_variant_text' if comps.get('transcript_id_input') else 'none')
        attach_c_hgvs_strand(out)
        reconfirm_ref_from_fasta(out)
        return out
    if lit_partial:
        attach_c_hgvs_strand(lit_partial)
        reconfirm_ref_from_fasta(lit_partial)
        lit_partial.setdefault('allele_query_source', '')
        note = lit_partial.get('resolve_note') or ''
        if 'deferred' not in note:
            lit_partial['resolve_note'] = (note + '; deferred until after c./rs/ClinVar/Ensembl/dbSNP').strip('; ')
        return lit_partial
    ass_in_text = export_genomic_allele_csv_normalize_assembly(export_genomic_allele_csv_RE_ASSEMBLY_IN_TEXT.search(blob).group(1) if export_genomic_allele_csv_RE_ASSEMBLY_IN_TEXT.search(blob) else None)
    note = 'no genomic map via ClinVar/PDF/Ensembl c. lift/dbSNP'
    if comps.get('transcript_id_input') or comps.get('c_hgvs_input') or comps.get('p_hgvs_input'):
        note += '; transcript/HGVS labels retained from literature'
    out = export_genomic_allele_csv_empty_resolve(comps, note, method='unresolved')
    if ass_in_text:
        out['assembly'] = ass_in_text
        out['assembly_source'] = 'literal_text'
        out['assembly_confirmed'] = True
        out['resolve_note'] += f'; assembly in text={ass_in_text}'
    attach_c_hgvs_strand(out)
    reconfirm_ref_from_fasta(out)
    return out

def export_genomic_allele_csv_patho_label_map(patho_rows: list[dict[str, Any]]) -> dict[tuple[str, str], str]:
    out: dict[tuple[str, str], str] = {}
    for r in patho_rows:
        key = (str(r.get('pmid') or ''), str(r.get('variant') or '').lower())
        lab = r.get('pdf_pathogenicity')
        if lab:
            out[key] = str(lab)
    return out
export_genomic_allele_csv_RE_GENE_NEAR_C = re.compile('\\b([A-Z][A-Z0-9]{1,14}(?:-[A-Z0-9]+)*)\\s*(?:\\([^)]{0,40}\\))?\\s*(?:\\(?)(?:NM_\\d+(?:\\.\\d+)?:)?c\\.', re.I)
export_genomic_allele_csv__GENE_CORRECT_STOPWORDS = {'VARIANT', 'VARIANTS', 'MUTATION', 'MUTATIONS', 'NOVEL', 'GENE', 'GENES', 'INSERTION', 'DELETION', 'SUBSTITUTION', 'HOMOZYGOUS', 'HETEROZYGOUS', 'HEMIZYGOUS', 'COMPOUND', 'NAMELY', 'WHICH', 'SHOWED', 'UNREPORTED', 'REPORTED', 'SIGNIFICANCE', 'PATHOGENIC', 'BENIGN', 'LIKELY', 'FRAMESHIFT', 'MISSENSE', 'NONSENSE', 'SPLICE', 'INTRONIC', 'EXONIC', 'CODING', 'PROTEIN', 'DOMAIN', 'PATIENT', 'PATIENTS', 'FAMILY', 'PROBAND', 'DNA', 'RNA', 'CDS', 'UTR', 'NN', 'OR', 'AND', 'THE', 'A', 'AN', 'IN', 'OF', 'TO', 'WITH', 'FROM', 'FOR', 'EACH', 'BOTH', 'SAME', 'OTHER', 'HUMAN', 'MOUSE', 'LOCUS', 'REGION', 'AUTOPHAGY', 'TRANSMEMBRANE', 'SUPPRESSOR', 'ADIPOCYTE', 'CONSERVED', 'LINKER', 'GAMMA', 'ALPHA', 'BETA', 'DELTA', 'KINASE', 'RECEPTOR', 'NUCLEAR', 'CANDIDATE', 'RISK', 'LEAD', 'FIGURE', 'TABLE', 'II', 'III', 'IV', 'POLYMORPHIC', 'PREVALENT', 'COMMON', 'ASSOCIATED', 'SIGNIFICANT', 'PUTATIVE', 'SUSCEPTIBILITY', 'WT', 'WILDTYPE', 'ROS', 'HPP', 'HMNX', 'VUS', 'LB', 'LP', 'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6', 'PP1', 'PP2', 'PP3', 'PP4', 'PP5', 'PS1', 'PS2', 'PS3', 'PS4', 'PVS1', 'BA1', 'BS1', 'BS2', 'BS3', 'BS4', 'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7', 'ACMG', 'CONTROL', 'CONTROLS'}

def export_genomic_allele_csv_load_pubtator_by_pmid(path: Path) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    if not path.exists():
        return out
    with path.open(encoding='utf-8') as f:
        for line in f:
            if not line.strip():
                continue
            row = json.loads(line)
            pmid = str(row.get('pmid') or '').strip()
            if pmid:
                out[pmid] = row
    return out

def export_genomic_allele_csv__pubtator_gene_id_for_variant(pubtator: dict[str, Any] | None, variant_text: str) -> str | None:
    if not pubtator or not variant_text:
        return None
    for v in pubtator.get('variants') or []:
        hv = v.get('hgvs') or v.get('name') or v.get('text') or ''
        if not vc_text_module.variant_tokens_match(str(variant_text), str(hv)):
            continue
        gid = v.get('gene_id')
        if gid is None:
            ident = str(v.get('identifier') or '')
            m = re.search('CorrespondingGene:(\\d+)', ident)
            if m:
                gid = m.group(1)
        if gid is not None and str(gid).strip():
            return str(gid).strip()
    return None

def export_genomic_allele_csv_correct_gene_from_text(*, gene: str | None, variant_text: str, evidence_blob: str | None=None, pdf_text: str | None=None, context_text: str | None=None, known_genes: set[str] | None=None, pubtator: dict[str, Any] | None=None, dbsnp: vc_dbsnp_DbSnpClient | None=None) -> tuple[str | None, str | None]:
    """Prefer literature / dbSNP / PubTator gene over a mismatched corpus gene.

    Priority (may override corpus gene):
      1) gene token immediately before/after the variant (evidence/title; not PDF)
      2) "located in the GENE locus" near this variant (evidence for rsIDs)
      3) HGVS only: ±2 sentence naming / SYMBOL variant window
      4) GENE near matching c. HGVS
      5) PubTator variant gene_id → symbol (catalog gene for rsIDs)
      6) dbSNP unique RefSNP gene locus (rsID only)

    Returns (gene, note). note is set when the corpus gene is overridden.
    Locus / dbSNP / PubTator gene_id hits do not require MANE membership.
    """

    def _ok_candidate(cand: str, *, require_known: bool) -> bool:
        cu = cand.upper()
        if cu in export_genomic_allele_csv__GENE_CORRECT_STOPWORDS:
            return False
        if not vc_text_module.is_plausible_gene_symbol(cand):
            return False
        if require_known and known_genes is not None and (cu not in known_genes):
            return False
        return True

    def _maybe_override(cand: str, reason: str) -> tuple[str | None, str | None]:
        cu = cand.upper()
        if gene and cu != str(gene).upper():
            return (cu, f'gene corrected from {reason}: {gene}→{cu}')
        return (cu, None)
    blobs = [x for x in (evidence_blob, context_text, pdf_text, variant_text) if x]
    rs = vc_dbsnp_normalize_rs_id(variant_text)
    evidence_only = [x for x in (evidence_blob, variant_text) if x]
    if not rs:
        for blob in evidence_only:
            adj = vc_text_module.gene_adjacent_to_variant(blob, variant_text)
            if adj and _ok_candidate(adj, require_known=True):
                return _maybe_override(adj, 'literature adjacency')
    if not rs:
        for blob in blobs:
            neigh = vc_text_module.gene_from_neighbor_window(blob, variant_text)
            if neigh and _ok_candidate(neigh, require_known=False):
                return _maybe_override(neigh, 'literature neighbor window')
    for blob in evidence_only if rs else blobs:
        loc = vc_text_module.gene_in_locus_near_variant(blob, variant_text)
        if not loc:
            continue
        if rs and _ok_candidate(loc, require_known=False):
            return _maybe_override(loc, 'literature locus phrase')
        if not rs and _ok_candidate(loc, require_known=True):
            return _maybe_override(loc, 'literature locus phrase')
    c_m = export_genomic_allele_csv_RE_C.search(variant_text or '')
    c_body = re.sub('\\s+', '', c_m.group(1) if c_m else '')
    if c_body:
        c_flex = re.escape(c_body)
        c_flex = c_flex.replace('>', '\\s*>\\s*')
        c_flex = re.sub('(?<=\\d)-(?=\\d)', '\\\\s*-\\\\s*', c_flex)
        c_flex = re.sub('(?<=\\d)\\+(?=\\d)', '\\\\s*\\\\+\\\\s*', c_flex)
        cre = re.compile(c_flex, re.I)
        for blob in blobs:
            for m in cre.finditer(blob):
                win = blob[max(0, m.start() - 80):m.end()]
                gm = export_genomic_allele_csv_RE_GENE_NEAR_C.search(win)
                if gm:
                    cand = gm.group(1).upper()
                    if _ok_candidate(cand, require_known=True):
                        return _maybe_override(cand, 'literature near c.')
    gid = export_genomic_allele_csv__pubtator_gene_id_for_variant(pubtator, variant_text)
    if gid:
        sym = vc_text_module.gene_id_to_symbol(gid, pubtator)
        if not sym and dbsnp is not None and rs:
            sym = dbsnp.lookup_gene_for_id(rs, gid)
        if sym and _ok_candidate(sym, require_known=False):
            return _maybe_override(sym, f'pubtator gene_id={gid}')
    if dbsnp is not None and rs:
        dgene = dbsnp.lookup_gene(rs)
        if dgene and _ok_candidate(dgene, require_known=False):
            return _maybe_override(dgene, 'dbsnp gene locus')
    return (str(gene).upper() if gene else None, None)

def export_genomic_allele_csv_attach_gff3_region(resolved: dict[str, Any], region_idx: Any, catalog_gene: str | None) -> None:
    """Fill genomic_region* from GFF3 overlap of the resolved GRCh38 allele."""
    blank = {'genomic_region': '', 'genomic_region_genes': '', 'genomic_region_feature': '', 'genomic_region_transcripts': '', 'genomic_region_note': ''}
    if region_idx is None or not resolved.get('position') or (not resolved.get('chromosome')):
        resolved.update(blank)
        return
    ann = region_idx.annotate(resolved.get('chromosome'), resolved.get('position'), ref=resolved.get('ref') or '', catalog_gene=catalog_gene)
    resolved.update(ann)

def export_genomic_allele_csv_row_from_resolved(*, sample_class: str, label: str, gene: str | None, variant_text: str, pmid: str, locus_key_src: str, pubmed_url: str, resolved: dict[str, Any], screen_grade: str='', link_tier: str='', review_tier: str='') -> dict[str, Any]:
    chrom = resolved.get('chromosome')
    return {'chromosome_ucsc': export_genomic_allele_csv_chrom_ucsc(chrom) or '', 'position': resolved.get('position') if resolved.get('position') is not None else '', 'ref': resolved.get('ref') or '', 'alt': resolved.get('alt') or '', 'c_hgvs': resolved.get('c_hgvs') or '', 'nt_pro_only': export_genomic_allele_csv_infer_nt_pro_only(variant_text), 'p_hgvs': resolved.get('p_hgvs') or '', 'c_hgvs_strand': resolved.get('c_hgvs_strand') or '', 'gene': gene or '', 'genomic_region': resolved.get('genomic_region') or '', 'genomic_region_genes': resolved.get('genomic_region_genes') or '', 'variant_text': variant_text, 'chromosome_accession': resolved.get('chromosome_accession') or '', 'assembly': resolved.get('assembly') or export_genomic_allele_csv_ASSEMBLY_PREFER, 'pubmed_url': pubmed_url, 'chromosome': chrom or '', 'chromosome_notation': export_genomic_allele_csv_CHROMOSOME_NOTATION, 'assembly_source': resolved.get('assembly_source') or 'pipeline_default_target', 'assembly_confirmed': 'true' if resolved.get('assembly_confirmed') else 'false', 'transcript_id': resolved.get('transcript_id') or '', 'transcript_id_input': resolved.get('transcript_id_input') or '', 'transcript_version': resolved.get('transcript_version') or '', 'transcript_source': resolved.get('transcript_source') or '', 'transcript_match': resolved.get('transcript_match') or '', 'genomic_region_feature': resolved.get('genomic_region_feature') or '', 'genomic_region_transcripts': resolved.get('genomic_region_transcripts') or '', 'genomic_region_note': resolved.get('genomic_region_note') or '', 'sample_class': sample_class, 'label': label, 'screen_grade': (screen_grade or '').upper(), 'link_tier': (link_tier or '').upper(), 'review_tier': review_tier or '', 'pmid': pmid, 'locus_key_src': locus_key_src, 'resolve_method': resolved.get('resolve_method') or '', 'allele_query_source': resolved.get('allele_query_source') or '', 'resolve_status': resolved.get('resolve_status') or 'failed', 'resolve_confidence': resolved.get('resolve_confidence') or 'low', 'resolve_note': resolved.get('resolve_note') or ''}

def export_genomic_allele_csv_genomic_allele_key(row: dict[str, Any]) -> tuple[str, str, str, str] | None:
    """Identity key for complete genomic alleles (chr/pos/ref/alt)."""
    chrom = str(row.get('chromosome_ucsc') or row.get('chromosome') or '').strip()
    pos = row.get('position')
    ref = str(row.get('ref') or '').strip().upper()
    alt = str(row.get('alt') or '').strip().upper()
    if not chrom or pos in (None, '') or (not ref) or (not alt):
        return None
    return (chrom, str(pos), ref, alt)

def export_genomic_allele_csv_genomic_dedup_rank(row: dict[str, Any]) -> tuple:
    """Lower is better when choosing one row per genomic allele."""
    status_rank = {'ok': 0, 'partial': 1, 'failed': 2}.get(str(row.get('resolve_status') or ''), 9)
    class_rank = 0 if row.get('sample_class') == 'positive' else 1
    review_rank = {'strict': 0, 'extended': 1, 'independent_neg': 2}.get(str(row.get('review_tier') or ''), 9)
    link_rank = {'A': 0, 'B': 1}.get(str(row.get('link_tier') or '').upper(), 9)
    screen_rank = {'A': 0, 'B': 1, 'C': 2}.get(str(row.get('screen_grade') or '').upper(), 9)
    nt_rank = {'nt': 0, 'both': 1, 'pro': 2, 'other': 3}.get(str(row.get('nt_pro_only') or ''), 9)
    return (status_rank, class_rank, review_rank, link_rank, screen_rank, nt_rank, str(row.get('pmid') or ''), str(row.get('locus_key_src') or ''))

def export_genomic_allele_csv_deduplicate_genomic_rows(rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], int, int]:
    """Keep one row per complete genomic allele (chr/pos/ref/alt).

    Incomplete rows (missing any of the four alleles) are always kept.
    Returns (deduped_rows, n_dropped, n_groups_collapsed).
    """
    groups: dict[tuple[str, str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for r in rows:
        key = export_genomic_allele_csv_genomic_allele_key(r)
        if key is not None:
            groups[key].append(r)
    out: list[dict[str, Any]] = []
    emitted: set[tuple[str, str, str, str]] = set()
    n_dropped = 0
    n_groups = 0
    for r in rows:
        key = export_genomic_allele_csv_genomic_allele_key(r)
        if key is None:
            out.append(r)
            continue
        if key in emitted:
            continue
        emitted.add(key)
        cands = groups[key]
        best = min(cands, key=export_genomic_allele_csv_genomic_dedup_rank)
        if len(cands) > 1:
            n_groups += 1
            n_dropped += len(cands) - 1
            dropped_keys = [str(c.get('locus_key_src') or c.get('variant_text') or '') for c in cands if c.get('locus_key_src') != best.get('locus_key_src') or c.get('pmid') != best.get('pmid')]
            best = dict(best)
            best['resolve_note'] = ((best.get('resolve_note') or '') + f'; genomic_dedup n={len(cands)} dropped=' + ','.join((x for x in dropped_keys if x))[:500]).strip('; ')
        out.append(best)
    return (out, n_dropped, n_groups)

def export_genomic_allele_csv_load_articles_by_pmid(path: Path) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    if not path.exists():
        return out
    with path.open(encoding='utf-8') as f:
        for line in f:
            if not line.strip():
                continue
            row = json.loads(line)
            pmid = str(row.get('pmid') or '').strip()
            if pmid:
                out[pmid] = row
    return out

def export_genomic_allele_csv_screen_fields_from_row(row: dict[str, Any], *, articles: dict[str, dict[str, Any]] | None=None, journal_cfg: dict[str, Any] | None=None, by_issn: dict[str, Any] | None=None, by_name: dict[str, Any] | None=None) -> tuple[str, str, str]:
    """Return (screen_grade/journal_tier, link_tier, review_tier)."""
    link_tier = str(row.get('tier') or row.get('link_tier') or '').strip().upper()
    review_tier = str(row.get('review_tier') or '').strip()
    grade = str(row.get('journal_tier') or row.get('screen_grade') or '').strip().upper()
    if grade not in {'A', 'B', 'C'} and journal_cfg is not None and (by_issn is not None) and (by_name is not None):
        pmid = str(row.get('pmid') or '')
        art = (articles or {}).get(pmid)
        annotated = review_module.annotate_row(row, art, journal_cfg, by_issn, by_name)
        grade = str(annotated.get('journal_tier') or 'C').upper()
    if grade not in {'A', 'B', 'C'}:
        grade = 'C'
    return (grade, link_tier, review_tier)
export_genomic_allele_csv__NEG_LINK_TIER_A_RULES = {'classify_paren', 'reclassified_lb', 'classify_then_nearby_var', 'were_classified_lb', 'variant_then_vus', 'variant_then_label', 'label_then_variant', 'control_cue_with_variant'}
export_genomic_allele_csv__NEG_LINK_TIER_B_RULES = {'window_rs_null_association', 'window_variant_then_vus', 'window_classify_then_nearby_var', 'nearest_other_to_label', 'window_label_then_variant', 'association_null'}
export_genomic_allele_csv__RE_NEG_VARIANT_TOKEN = re.compile('\\b(?:rs\\d{3,}|(?:c\\.|p\\.)[A-Za-z0-9_()*?>.\\-+/]+)\\b', re.I)
export_genomic_allele_csv__RE_NEG_CUE = re.compile('\\b(?:likely\\s+benign|benign|not\\s+pathogenic|non-?pathogenic|VUS|uncertain(?:\\s+clinical)?\\s+significance|not\\s+associated|no\\s+(?:significant\\s+)?association|non-significant\\s+association|was\\s+not\\s+associated|control(?:\\s+variant)?s?|wild[- ]?type)\\b', re.I)

def export_genomic_allele_csv__variant_mentioned_in_text(variant: str, text: str) -> bool:
    if not variant or not text:
        return False
    v = re.sub('\\s+', '', variant.strip())
    t = re.sub('\\s+', '', text)
    if not v:
        return False
    if v.lower() in t.lower():
        return True
    core = v.split(':')[-1]
    return bool(core) and core.lower() in t.lower()

def export_genomic_allele_csv_assign_negative_link_tier(row: dict[str, Any]) -> tuple[str, str]:
    """Assign link_tier A/B/C for independent negatives (mirrors positive Stage2 idea).

    A — 同句共现：证据句同时含本变异与良性/VUS/明确无关联等负向描述，且抽取规则为
        分类/对照等高精度同句规则。
    B — 较弱同句/同窗：association_null、window_*、nearest_* 等；或仅有一侧信号。
    C — 证据不足：无证据句，或无法确认变异与负向描述共现。

    Returns (link_tier, reason).
    """
    preset = str(row.get('tier') or row.get('link_tier') or '').strip().upper()
    if preset in {'A', 'B', 'C'}:
        return (preset, 'preset')
    variant = str(row.get('negative_variant') or row.get('variant') or '')
    evidence = str(row.get('negative_evidence') or row.get('evidence_text') or '')
    rule = str(row.get('negative_rule') or row.get('rule') or '').strip()
    has_var = export_genomic_allele_csv__variant_mentioned_in_text(variant, evidence)
    has_cue = bool(export_genomic_allele_csv__RE_NEG_CUE.search(evidence)) if evidence else False
    if not evidence.strip():
        return ('C', 'no_negative_evidence')
    if has_var and has_cue and (rule in export_genomic_allele_csv__NEG_LINK_TIER_A_RULES):
        return ('A', f'same_sentence_classify_or_control:{rule}')
    if has_var and has_cue and (rule in export_genomic_allele_csv__NEG_LINK_TIER_B_RULES):
        return ('B', f'same_sentence_or_window_association:{rule}')
    if has_var and has_cue:
        if 'window' in rule or 'nearest' in rule or 'null' in rule:
            return ('B', f"same_sentence_windowish:{rule or 'unspecified'}")
        return ('A', f"same_sentence_variant_and_neg_cue:{rule or 'unspecified'}")
    if has_var or has_cue:
        return ('B', 'partial_variant_or_cue_only')
    if export_genomic_allele_csv__RE_NEG_VARIANT_TOKEN.search(evidence) and has_cue:
        return ('B', 'generic_variant_token_and_cue')
    return ('C', 'unconfirmed_cooccurrence')

def export_genomic_allele_csv_positive_exclude_keys(rows: list[dict[str, Any]]) -> set[str]:
    keys: set[str] = set()
    for r in rows:
        gene = r.get('gene')
        variant = r.get('variant')
        lk = str(r.get('locus_key') or '').strip().lower()
        if lk:
            keys.add(lk)
        if gene and variant:
            keys.add(f'{str(gene).strip().lower()}|{str(variant).strip().lower()}')
        if variant:
            keys.add(f'|{str(variant).strip().lower()}')
    return keys

def export_genomic_allele_csv_export_rows_csv(rows: list[dict[str, Any]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=export_genomic_allele_csv_CSV_FIELDS, extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow(r)

def export_genomic_allele_csv_export_rows_jsonl(rows: list[dict[str, Any]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open('w', encoding='utf-8') as f:
        for r in rows:
            obj = {k: r.get(k, '') for k in export_genomic_allele_csv_CSV_FIELDS}
            f.write(json.dumps(obj, ensure_ascii=False) + '\n')

def export_genomic_allele_csv_main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pos', type=Path, default=export_genomic_allele_csv_DEFAULT_POS)
    ap.add_argument('--pos-patho', type=Path, default=export_genomic_allele_csv_DEFAULT_POS_PATHO)
    ap.add_argument('--neg', type=Path, default=export_genomic_allele_csv_DEFAULT_NEG)
    ap.add_argument('--neg-hard-only', action='store_true', help=f'use hard independent negatives only (default file: {export_genomic_allele_csv_DEFAULT_NEG_HARD.name})')
    ap.add_argument('--neg-absent', type=Path, default=None, help='optional; omit for strict-only negatives (default)')
    ap.add_argument('--clinvar', type=Path, default=export_genomic_allele_csv_DEFAULT_CLINVAR)
    ap.add_argument('--out-dir', type=Path, default=export_genomic_allele_csv_DEFAULT_OUT_DIR)
    ap.add_argument('--pdf-cache', type=Path, default=export_genomic_allele_csv_DEFAULT_PDF_CACHE)
    ap.add_argument('--ensembl-gff3', type=Path, default=vc_ensembl_c_map_module.DEFAULT_GFF3)
    ap.add_argument('--ensembl-fasta', type=Path, default=vc_ensembl_c_map_module.DEFAULT_FASTA)
    ap.add_argument('--ensembl-fai', type=Path, default=vc_ensembl_c_map_module.DEFAULT_FAI)
    ap.add_argument('--ensembl-cache', type=Path, default=export_genomic_allele_csv_DEFAULT_ENSEMBL_CACHE)
    ap.add_argument('--no-pdf', action='store_true')
    ap.add_argument('--no-ensembl', action='store_true')
    ap.add_argument('--gff3-region-cache', type=Path, default=vc_gff3_region_default_region_cache(export_genomic_allele_csv_ROOT), help='pickle cache for Ensembl GFF3 gene/exon/CDS interval index')
    ap.add_argument('--no-gff3-region', action='store_true', help='skip CDS/exon/gene overlap annotation after coordinate resolve')
    ap.add_argument('--dbsnp-cache', type=Path, default=vc_dbsnp_default_dbsnp_cache(export_genomic_allele_csv_ROOT), help='cache directory for NCBI RefSNP JSON payloads')
    ap.add_argument('--no-dbsnp', action='store_true', help='disable NCBI dbSNP network fallback for unresolved rsIDs')
    ap.add_argument('--pubtator', type=Path, default=export_genomic_allele_csv_DEFAULT_PUBTATOR, help='PubTator parsed JSONL used for rsID→gene_id correction')
    ap.add_argument('--no-pubtator-gene', action='store_true', help='disable PubTator gene_id gene correction')
    args = ap.parse_args()
    logger = export_genomic_allele_csv_setup_logger()
    if args.neg_hard_only:
        if args.neg == export_genomic_allele_csv_DEFAULT_NEG or not args.neg:
            args.neg = export_genomic_allele_csv_DEFAULT_NEG_HARD
        logger.info('neg-hard-only: using %s', args.neg)
    idx = export_genomic_allele_csv_build_clinvar_genomic_index(args.clinvar, logger)
    ensembl_idx = None
    fasta = None
    if not args.no_ensembl and args.ensembl_gff3.exists() and args.ensembl_fasta.exists():
        ensembl_idx = vc_ensembl_c_map_module.build_mane_index(args.ensembl_gff3, logger, args.ensembl_cache)
        fasta = vc_ensembl_c_map_module.FastaReader(args.ensembl_fasta, args.ensembl_fai)
        logger.info('Ensembl FASTA ready: %s', args.ensembl_fasta)
    else:
        logger.warning('Ensembl c. lift disabled (missing GFF3/FASTA or --no-ensembl)')
    region_idx = None
    if not args.no_gff3_region and args.ensembl_gff3.exists():
        region_idx = vc_gff3_region_load_gff3_region_index(args.ensembl_gff3, cache_path=args.gff3_region_cache, logger=logger)
    elif not args.no_gff3_region:
        logger.warning('GFF3 region annotation disabled (missing %s)', args.ensembl_gff3)
    pdf_cache = None if args.no_pdf else args.pdf_cache
    dbsnp_client = None
    if not args.no_dbsnp:
        dbsnp_client = vc_dbsnp_DbSnpClient(cache_dir=args.dbsnp_cache, logger=logger)
        logger.info('dbSNP fallback enabled; cache=%s', args.dbsnp_cache)
    pubtator_by_pmid: dict[str, dict[str, Any]] = {}
    if not args.no_pubtator_gene:
        pubtator_by_pmid = export_genomic_allele_csv_load_pubtator_by_pmid(args.pubtator)
        logger.info('PubTator gene map loaded: %s PMIDs from %s', len(pubtator_by_pmid), args.pubtator)
    known_genes = set(ensembl_idx.by_gene.keys()) if ensembl_idx else None
    patho_map = export_genomic_allele_csv_patho_label_map(export_genomic_allele_csv_load_jsonl(args.pos_patho))
    pos_rows_in = export_genomic_allele_csv_load_jsonl(args.pos)
    journal_cfg = json.loads(export_genomic_allele_csv_DEFAULT_JOURNAL_CONFIG.read_text(encoding='utf-8'))
    by_issn, by_name = review_module.build_indexes(journal_cfg)
    articles_by_pmid = export_genomic_allele_csv_load_articles_by_pmid(export_genomic_allele_csv_DEFAULT_ARTICLES)
    logger.info('journal priority loaded (%s journals); articles=%s', len(journal_cfg.get('journals') or []), len(articles_by_pmid))
    pos_exclude = export_genomic_allele_csv_positive_exclude_keys(pos_rows_in)
    neg_paths = [args.neg]
    if args.neg_absent:
        neg_paths.append(args.neg_absent)
    out_rows: list[dict[str, Any]] = []
    for r in pos_rows_in:
        pmid = str(r.get('pmid') or '')
        variant = str(r.get('variant') or '')
        gene = r.get('gene')
        label = patho_map.get((pmid, variant.lower())) or 'positive_candidate'
        screen_grade, link_tier, review_tier = export_genomic_allele_csv_screen_fields_from_row(r, articles=articles_by_pmid, journal_cfg=journal_cfg, by_issn=by_issn, by_name=by_name)
        pdf_text = vc_pdf_genomic_load_pdf_text(pmid, pdf_cache) if pdf_cache else None
        gene_in = str(gene).upper() if gene else None
        gene, gene_note = export_genomic_allele_csv_correct_gene_from_text(gene=gene_in, variant_text=variant, evidence_blob=str(r.get('evidence_text') or ''), pdf_text=pdf_text, context_text=str((articles_by_pmid.get(pmid) or {}).get('abstract') or ''), known_genes=known_genes, pubtator=pubtator_by_pmid.get(pmid), dbsnp=dbsnp_client)
        resolved = export_genomic_allele_csv_resolve_allele(variant_text=variant, gene=gene, evidence_blob=str(r.get('evidence_text') or ''), idx=idx, pmid=pmid, pdf_cache=pdf_cache, ensembl_idx=ensembl_idx, fasta=fasta, dbsnp=dbsnp_client)
        export_genomic_allele_csv_attach_gff3_region(resolved, region_idx, gene)
        if gene_note:
            resolved['resolve_note'] = ((resolved.get('resolve_note') or '') + '; ' + gene_note).strip('; ')
        locus_key_src = f'{gene}|{variant}' if gene else str(r.get('locus_key') or variant)
        if gene_note and r.get('locus_key') and (str(r.get('locus_key')) != locus_key_src):
            resolved['resolve_note'] = ((resolved.get('resolve_note') or '') + f"; locus_key_src {r.get('locus_key')}→{locus_key_src}").strip('; ')
        out_rows.append(export_genomic_allele_csv_row_from_resolved(sample_class='positive', label=label, gene=gene, variant_text=variant, pmid=pmid, locus_key_src=locus_key_src, pubmed_url=str(r.get('pubmed_url') or ''), resolved=resolved, screen_grade=screen_grade, link_tier=link_tier, review_tier=review_tier))
    n_neg_collide = 0
    for npath in neg_paths:
        if not npath or not npath.exists():
            continue
        for r in export_genomic_allele_csv_load_jsonl(npath):
            pmid = str(r.get('pmid') or '')
            variant = str(r.get('negative_variant') or '')
            gene = r.get('negative_gene')
            lk = str(r.get('negative_locus_key') or '').strip().lower()
            bare = f'|{variant.strip().lower()}' if variant else ''
            gene_lk = f'{str(gene).strip().lower()}|{variant.strip().lower()}' if gene and variant else ''
            if pos_exclude and (lk and lk in pos_exclude or (bare and bare in pos_exclude) or (gene_lk and gene_lk in pos_exclude)):
                n_neg_collide += 1
                continue
            screen_grade, link_tier, review_tier = export_genomic_allele_csv_screen_fields_from_row(r, articles=articles_by_pmid, journal_cfg=journal_cfg, by_issn=by_issn, by_name=by_name)
            if not link_tier:
                link_tier, tier_reason = export_genomic_allele_csv_assign_negative_link_tier(r)
                r = dict(r)
                r['link_tier'] = link_tier
                r['link_tier_reason'] = tier_reason
            pdf_text = vc_pdf_genomic_load_pdf_text(pmid, pdf_cache) if pdf_cache else None
            gene_in = str(gene).upper() if gene else None
            gene, gene_note = export_genomic_allele_csv_correct_gene_from_text(gene=gene_in, variant_text=variant, evidence_blob=str(r.get('negative_evidence') or ''), pdf_text=pdf_text, context_text=str((articles_by_pmid.get(pmid) or {}).get('abstract') or ''), known_genes=known_genes, pubtator=pubtator_by_pmid.get(pmid), dbsnp=dbsnp_client)
            resolved = export_genomic_allele_csv_resolve_allele(variant_text=variant, gene=gene, evidence_blob=str(r.get('negative_evidence') or ''), idx=idx, pmid=pmid, pdf_cache=pdf_cache, ensembl_idx=ensembl_idx, fasta=fasta, dbsnp=dbsnp_client)
            export_genomic_allele_csv_attach_gff3_region(resolved, region_idx, gene)
            if gene_note:
                resolved['resolve_note'] = ((resolved.get('resolve_note') or '') + '; ' + gene_note).strip('; ')
            locus_key_src = f'{gene}|{variant}' if gene else str(r.get('negative_locus_key') or variant)
            if gene_note and r.get('negative_locus_key') and (str(r.get('negative_locus_key')) != locus_key_src):
                resolved['resolve_note'] = ((resolved.get('resolve_note') or '') + f"; locus_key_src {r.get('negative_locus_key')}→{locus_key_src}").strip('; ')
            out_rows.append(export_genomic_allele_csv_row_from_resolved(sample_class='negative', label=str(r.get('negative_label') or 'negative'), gene=gene, variant_text=variant, pmid=pmid, locus_key_src=locus_key_src, pubmed_url=str(r.get('pubmed_url') or ''), resolved=resolved, screen_grade=screen_grade, link_tier=link_tier, review_tier=review_tier or 'independent_neg'))
    if n_neg_collide:
        logger.info('dropped %s negatives colliding with positive locus keys', n_neg_collide)
    n_before_dedup = len(out_rows)
    out_rows, n_genomic_dedup, n_genomic_groups = export_genomic_allele_csv_deduplicate_genomic_rows(out_rows)
    if n_genomic_dedup:
        logger.info('genomic allele dedup: %s → %s rows (dropped %s across %s loci)', n_before_dedup, len(out_rows), n_genomic_dedup, n_genomic_groups)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    all_csv = args.out_dir / 'variants_genomic_all.csv'
    ok_csv = args.out_dir / 'variants_genomic_ok.csv'
    pos_csv = args.out_dir / 'variants_genomic_positive.csv'
    neg_csv = args.out_dir / 'variants_genomic_negative.csv'
    pos_ok_csv = args.out_dir / 'variants_genomic_positive_ok.csv'
    neg_ok_csv = args.out_dir / 'variants_genomic_negative_ok.csv'
    all_jsonl = args.out_dir / 'variants_genomic_all.jsonl'
    ok_jsonl = args.out_dir / 'variants_genomic_ok.jsonl'
    pos_jsonl = args.out_dir / 'variants_genomic_positive.jsonl'
    neg_jsonl = args.out_dir / 'variants_genomic_negative.jsonl'
    ok_rows = [r for r in out_rows if r.get('resolve_status') == 'ok' and r.get('chromosome') and (r.get('position') != '') and r.get('ref')]
    pos_rows = [r for r in out_rows if r.get('sample_class') == 'positive']
    neg_rows = [r for r in out_rows if r.get('sample_class') == 'negative']
    pos_ok_rows = [r for r in ok_rows if r.get('sample_class') == 'positive']
    neg_ok_rows = [r for r in ok_rows if r.get('sample_class') == 'negative']
    export_genomic_allele_csv_export_rows_csv(out_rows, all_csv)
    export_genomic_allele_csv_export_rows_csv(ok_rows, ok_csv)
    export_genomic_allele_csv_export_rows_csv(pos_rows, pos_csv)
    export_genomic_allele_csv_export_rows_csv(neg_rows, neg_csv)
    export_genomic_allele_csv_export_rows_csv(pos_ok_rows, pos_ok_csv)
    export_genomic_allele_csv_export_rows_csv(neg_ok_rows, neg_ok_csv)
    export_genomic_allele_csv_export_rows_jsonl(out_rows, all_jsonl)
    export_genomic_allele_csv_export_rows_jsonl(ok_rows, ok_jsonl)
    export_genomic_allele_csv_export_rows_jsonl(pos_rows, pos_jsonl)
    export_genomic_allele_csv_export_rows_jsonl(neg_rows, neg_jsonl)
    contract = {'scope': 'highconf_positive_plus_independent_full_pass_negatives' if 'negatives_independent' in str(args.neg) else 'highconf_positive_plus_negatives' if not args.neg_absent else 'highconf_plus_absent_negatives', 'output_files': {'all_csv': str(all_csv), 'ok_csv': str(ok_csv), 'positive_csv': str(pos_csv), 'negative_csv': str(neg_csv), 'positive_ok_csv': str(pos_ok_csv), 'negative_ok_csv': str(neg_ok_csv), 'all_jsonl': str(all_jsonl), 'ok_jsonl': str(ok_jsonl), 'positive_jsonl': str(pos_jsonl), 'negative_jsonl': str(neg_jsonl)}, 'column_order': export_genomic_allele_csv_CSV_FIELDS, 'primary_columns': export_genomic_allele_csv_CSV_PRIMARY_FIELDS, 'nt_pro_only': {'field': 'nt_pro_only', 'meaning': 'HGVS level present in literature variant_text only', 'values': {'nt': 'only nucleotide/c. (or NM_:c.) in variant_text', 'pro': 'only protein/p. in variant_text', 'both': 'both c. and p. in variant_text', 'other': 'neither (rsID, genomic literal, etc.)'}}, 'chromosome_field': 'chromosome', 'chromosome_notation': export_genomic_allele_csv_CHROMOSOME_NOTATION, 'chromosome_notation_meaning': "Bare NCBI-style contig name without 'chr' prefix: 1-22, X, Y, MT. chromosome_ucsc uses UCSC style; chromosome_accession holds NC_… when from ClinVar.", 'position_coordinate_system': '1-based (VCF PositionVCF / HGVS g. style start)', 'assembly_prefer': export_genomic_allele_csv_ASSEMBLY_PREFER, 'assembly_fields': {'assembly': 'GRCh38 or GRCh37', 'assembly_source': 'clinvar_variant_summary | literal_text | pipeline_default_target | clinvar_preferred_target', 'assembly_confirmed': 'true only when from ClinVar row or explicit text mention'}, 'transcript_fields': {'transcript_id': 'mapped or literature NM_ accession (without version unless only input)', 'transcript_id_input': 'NM_ parsed from literature variant_text/evidence', 'transcript_version': 'version from ClinVar Name or literature', 'transcript_source': 'clinvar_name | literature_variant_text | none', 'transcript_match': 'exact_nm | exact_nm_version | nm_version_mismatch | clinvar_preferred_transcript | clinvar_other_transcript | input_only | ambiguous | none | not_applicable', 'c_hgvs': 'normalized c. from literature or ClinVar Name', 'p_hgvs': 'normalized p. from literature or ClinVar Name'}, 'n_all': len(out_rows), 'n_ok': len(ok_rows), 'n_positive': len(pos_rows), 'n_negative': len(neg_rows), 'n_positive_ok': len(pos_ok_rows), 'n_negative_ok': len(neg_ok_rows), 'n_before_genomic_dedup': n_before_dedup, 'n_genomic_dedup_dropped': n_genomic_dedup, 'n_genomic_dedup_groups': n_genomic_groups, 'genomic_dedup': {'key': 'chromosome_ucsc|position|ref|alt (complete alleles only)', 'policy': 'keep best row by resolve_status, sample_class(positive>negative), review_tier, link_tier, screen_grade, nt_pro_only(nt>both>pro)'}, 'resolve_status_counts': dict(Counter((r['resolve_status'] for r in out_rows))), 'resolve_method_counts': dict(Counter((r['resolve_method'] for r in out_rows))), 'allele_query_source_counts': dict(Counter((r.get('allele_query_source') or '' for r in out_rows))), 'transcript_match_counts': dict(Counter((r['transcript_match'] for r in out_rows))), 'nt_pro_only_counts': dict(Counter((r['nt_pro_only'] for r in out_rows))), 'genomic_region_counts': dict(Counter((r.get('genomic_region') or '' for r in out_rows))), 'genomic_region': {'field': 'genomic_region', 'source': 'Ensembl/GENCODE 115 GFF3 overlap of resolved GRCh38 chr:pos (1-based)', 'values': {'coding': 'overlaps protein_coding CDS', 'utr': "overlaps 5'/3' UTR (GFF3 UTR features or protein-coding exon outside CDS)", 'intron': 'inside a protein-coding gene body, not in exon/UTR/CDS', 'ncrna_exon': 'overlaps lncRNA/miRNA/etc exon', 'ncrna_gene': 'inside an ncRNA gene body, not in an exon', 'pseudogene': 'inside a pseudogene body, not in an exon', 'intergenic': 'no overlapping gene/exon/CDS/UTR', '': 'no resolved coordinates'}, 'related_fields': ['genomic_region_genes', 'genomic_region_feature', 'genomic_region_transcripts', 'genomic_region_note']}, 'assembly_counts_ok': dict(Counter((r['assembly'] for r in ok_rows))), 'sample_class_counts': dict(Counter((r['sample_class'] for r in out_rows))), 'screen_grade_counts': dict(Counter((r.get('screen_grade') or '' for r in out_rows))), 'screen_grade_by_class': {cls: dict(Counter((r.get('screen_grade') or '' for r in out_rows if r['sample_class'] == cls))) for cls in sorted({r['sample_class'] for r in out_rows})}, 'link_tier_counts': dict(Counter((r.get('link_tier') or '' for r in out_rows))), 'review_tier_counts': dict(Counter((r.get('review_tier') or '' for r in out_rows))), 'screen_grade': {'field': 'screen_grade', 'meaning': 'Journal priority tier A/B/C (review ordering only; not pathogenicity)', 'source': 'configs/journal_priority.json via journal_tier'}, 'link_tier': {'field': 'link_tier', 'meaning': 'For positives: Stage2/6 abstract/fulltext link tier A/B. For independent negatives: assigned at export — A=same-sentence variant+benign/VUS/null-association (classify/control rules); B=window/association_null or partial cue; C=unconfirmed.'}, 'review_tier': {'field': 'review_tier', 'meaning': 'strict|extended for positives; independent_neg for independent negatives'}, 'resolver_order': ['literal_text', 'clinvar', 'pdf_literal_genomic', 'ensembl_mane_p_map', 'ensembl_mane_c_map', 'rs_dbsnp', 'c_hgvs_alleles_only', 'unresolved'], 'allele_query_source': {'field': 'allele_query_source', 'meaning': 'External network source used to fill genomic alleles when local lifts fail', 'values': {'dbsnp': 'NCBI Variation RefSNP API (https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/)', '': 'filled without dbSNP network query (ClinVar/PDF/Ensembl/literal/unresolved)'}}, 'note': 'Strict positives are often ClinVar-absent; Ensembl MANE c. lift, PDF literal coords, and dbSNP rsID lookup fill remaining rows. Pure p. without paired c. in literature may stay unresolved.'}
    (args.out_dir / 'genomic_allele_export_contract.json').write_text(json.dumps(contract, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    logger.info('wrote %s (%s rows)', all_csv, len(out_rows))
    logger.info('wrote %s (%s rows ok)', ok_csv, len(ok_rows))
    logger.info('wrote %s (%s positive)', pos_csv, len(pos_rows))
    logger.info('wrote %s (%s negative)', neg_csv, len(neg_rows))
    logger.info('wrote %s (%s positive ok)', pos_ok_csv, len(pos_ok_rows))
    logger.info('wrote %s (%s negative ok)', neg_ok_csv, len(neg_ok_rows))
    logger.info('wrote %s (%s rows)', all_jsonl, len(out_rows))
    logger.info('wrote %s (%s rows ok)', ok_jsonl, len(ok_rows))
    logger.info('contract: %s', json.dumps(contract, ensure_ascii=False))

# Public helper compatibility aliases.
build_gff3_region_index = vc_gff3_region_build_gff3_region_index
merge_comps_from_evidence = export_genomic_allele_csv_merge_comps_from_evidence
parse_variant_components = export_genomic_allele_csv_parse_variant_components
try_literal_genomic = export_genomic_allele_csv_try_literal_genomic
find_paired_c_for_rs = export_genomic_allele_csv_find_paired_c_for_rs

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
        'export': export_genomic_allele_csv_main,
    })
