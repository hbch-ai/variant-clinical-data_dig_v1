#!/usr/bin/env python3
"""Shared variant / clinical text patterns and span utilities."""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from typing import Any

VARIANT_PATTERNS: list[tuple[str, re.Pattern[str]]] = [
    (
        "hgvs_c",
        # Prefer spaced or tight substitutions (c.2663A > G / c.2663A>G) before bare allele tails.
        re.compile(
            r"\bc\.(?:-?\d+)(?:[+-]\d+)?(?:_\d+(?:[+-]\d+)?)?"
            r"(?:delins|del|dup|ins|[ACGT]\s*>\s*[ACGT]|[ACGT]+|[*=?])\b",
            re.I,
        ),
    ),
    (
        "hgvs_p_3letter",
        re.compile(r"\bp\.(?:[A-Z][a-z]{2})(?:\d+)(?:[A-Z][a-z]{2}|\*|Ter|fs(?:Ter?\d+)?|=|\?)\b"),
    ),
    ("hgvs_p_1letter", re.compile(r"\bp\.[A-Z]\d+(?:[A-Z]|\*|fs)\b")),
    ("hgvs_p_paren", re.compile(r"\bp\.\([A-Za-z0-9_*?=]+\)")),
    ("rsid", re.compile(r"\brs\d{4,}\b", re.I)),
    ("nm_np_with_change", re.compile(r"\b(?:NM_|NP_|NC_|NG_)\d+(?:\.\d+)?:[cp]\.[^\s,;]+")),
    ("ivs", re.compile(r"\bIVS\d+(?:[+-]\d+)?[ACGT]\s*>\s*[ACGT]\b", re.I)),
    (
        "pathogenic_phrase",
        re.compile(
            r"\b(?:pathogenic|likely pathogenic)\s+(?:variant|mutation|variants|mutations)\b"
            r"|\b(?:novel|de novo|biallelic|homozygous|heterozygous|compound heterozygous)\s+"
            r"(?:variant|mutation|variants|mutations)\b"
            r"|\b(?:loss[- ]of[- ]function|frameshift|nonsense|missense|splice[- ]site)\s+"
            r"(?:variant|mutation|variants|mutations)\b",
            re.I,
        ),
    ),
]

SOFT_VARIANT = re.compile(
    r"\b(?:genetic variant|sequence variant|genomic variant|germline mutation|"
    r"somatic mutation|gene mutation|mutational|variant carriers?)\b",
    re.I,
)

# Discourse cues only — not disease/phenotype entities.
CLINICAL_CUE_PATTERN = re.compile(
    r"\b(?:clinical\s+(?:manifestations?|features?|presentations?|phenotype|findings?|symptoms?|"
    r"spectrum|course|picture)|phenotypic\s+(?:spectrum|heterogeneity|features?|presentation)|"
    r"presenting\s+with|presented\s+with|characterized\s+by|symptoms?\s+include|"
    r"manifest(?:s|ed|ing)?\s+as)\b",
    re.I,
)

# Bare heads that must not stand alone as clinical entities (need a named prefix).
GENERIC_CLINICAL_TOKENS = frozenset(
    {
        "syndrome",
        "syndromic",
        "disorder",
        "disease",
        "dystrophy",
        "dysplasia",
        "anomaly",
        "anomalies",
        "malformation",
        "hypoplasia",
        "hyperplasia",
        "phenotype",
        "neuropathy",
        "myopathy",
        "neoplasm",
        "neoplasms",
    }
)

# Leading filler stripped when normalizing "... syndrome" matches.
# Keep disease-name adjectives (hereditary/congenital/severe/…) — they are part of many eponyms.
_CLINICAL_PREFIX_STOPWORDS = frozenset(
    {
        "a",
        "an",
        "the",
        "this",
        "that",
        "these",
        "those",
        "such",
        "any",
        "some",
        "other",
        "another",
        "rare",
        "novel",
        "new",
        "known",
        "common",
        "patient",
        "patients",
        "case",
        "cases",
        "child",
        "children",
        "family",
        "families",
        "with",
        "and",
        "or",
        "of",
        "in",
        "on",
        "to",
        "for",
        "by",
        "from",
        "as",
        "is",
        "are",
        "was",
        "were",
        "including",
        "due",
        "caused",
        "causing",
        "presenting",
        "described",
        "reported",
        "description",
        "suspicion",
        "suspected",
        "led",
        "first",
        "one",
        "two",
        "both",
        "several",
        "multiple",
        "various",
        "different",
        "human",
        "clinical",
        "molecular",
        "pathogenic",
        "pathological",
        "associated",
        "related",  # alone; "CAMRQ-related" stays one token via hyphen
        "complex",
        "form",
        "type",  # handled via type-suffix on head; bare "type" mid-phrase is filler
    }
)

# Inheritance-only adjectives: not enough to make "... disorder" a named disease.
_INHERITANCE_ONLY_TOKENS = frozenset(
    {
        "autosomal",
        "recessive",
        "dominant",
        "x-linked",
        "xlinked",
        "inherited",
        "genetic",
        "hereditary",
        "congenital",
        "familial",
        "sporadic",
        "maternal",
        "paternal",
        "biallelic",
        "monoallelic",
        "homozygous",
        "heterozygous",
        "severe",
        "mild",
        "moderate",
        "typical",
        "atypical",
        "primary",
        "secondary",
        "acquired",
        "progressive",
        "early-onset",
        "late-onset",
    }
)

# Token before disease head (allows Alström, Dystonia-deafness, CAMRQ-related).
_CLIN_PREFIX_TOKEN = r"[A-Za-z][\w'*-]*"
_CLIN_SEP = r"[\s/\-‐‑–—]+"
# Heads that require ≥1 content prefix (no bare "syndrome" / "dysplasia").
_NAMED_DISEASE_HEADS = (
    r"syndrome|disorder|dystrophy|dysplasia|anomaly|anomalies|malformation|"
    r"hypoplasia|hyperplasia|neuropathy|myopathy"
)
_NAMED_DISEASE_HEADS_SET = frozenset(
    {
        "syndrome",
        "disorder",
        "dystrophy",
        "dysplasia",
        "anomaly",
        "anomalies",
        "malformation",
        "hypoplasia",
        "hyperplasia",
        "neuropathy",
        "myopathy",
    }
)
_NAMED_DISEASE_TYPE_SUFFIX = r"(?:[\s/\-‐‑–—]+(?:type[\s/\-‐‑–—]+)?[IVXLC\d]+)?"

NAMED_DISEASE_PATTERN = re.compile(
    rf"\b(?:{_CLIN_PREFIX_TOKEN}(?:{_CLIN_SEP}{_CLIN_PREFIX_TOKEN}){{0,5}}){_CLIN_SEP}"
    rf"(?:{_NAMED_DISEASE_HEADS}){_NAMED_DISEASE_TYPE_SUFFIX}\b",
    re.I,
)
# syndromic + up to 3 following content tokens (stopwords not matched via (?!)).
SYNDROMIC_PHRASE_PATTERN = re.compile(
    r"\bsyndromic"
    r"(?:\s+(?!associated\b|with\b|in\b|of\b|and\b|or\b|to\b|for\b|by\b|from\b|as\b|due\b)"
    r"[A-Za-z][\w'*-]*){1,3}\b",
    re.I,
)

CLINICAL_PATTERNS: list[tuple[str, re.Pattern[str]]] = [
    ("clinical_cue", CLINICAL_CUE_PATTERN),
    ("hpo_like", re.compile(r"\bHP:\d{7}\b|\bHuman Phenotype Ontology\b|\bHPO\b")),
    # Multi-word named diseases first (… syndrome / … dysplasia / syndromic …).
    ("syndrome_disease", NAMED_DISEASE_PATTERN),
    ("syndrome_disease", SYNDROMIC_PHRASE_PATTERN),
    (
        "syndrome_disease",
        re.compile(
            r"\b(?:(?:hereditary\s+)?(?:severe\s+)?insulin[- ]resistance\s+syndrome|H-SIRS|"
            r"acanthosis\s+nigricans|"
            r"dilated\s+cardiomyopathy|hypertrophic\s+cardiomyopathy|"
            r"restrictive\s+cardiomyopathy|arrhythmogenic\s+(?:right\s+ventricular\s+)?"
            r"cardiomyopathy|"
            r"cardiomyopathy|epilepsy|seizures?|"
            r"intellectual disability|developmental delay|short stature|"
            r"hearing loss|visual impairment|ataxia)\b",
            re.I,
        ),
    ),
    (
        "phenotype_term",
        re.compile(
            r"\b(?:hyperinsulinemia|hypoglycemia|hyperglycemia|insulin\s+resistance|"
            r"diabetes\s+mellitus|fatty\s+liver(?:\s+disease)?|hyperandrogenism|hirsutism|"
            r"mental\s+retardation|overweight|obesity|fasting\s+hypoglycemia|"
            r"postprandial\s+hyperglycemia|impaired\s+glucose\s+regulation|"
            r"glucose\s+intolerance|polycystic\s+ovary\s+syndrome|PCOS|"
            r"hypertension|hepatomegaly|lipodystrophy)\b",
            re.I,
        ),
    ),
    (
        "hematology_oncology",
        re.compile(
            r"\b(?:myeloproliferative neoplasm(?:s)?|MPN|myelofibrosis|polycythemia vera|"
            r"essential thrombocythemia|myeloid leukemia|myeloid neoplasm(?:s)?|"
            r"acute myeloid leukemia|AML|clonal hematopoiesis(?: of indeterminate potential)?|"
            r"CHIP|hematopoiet(?:ic|ies)|hematologic(?: neoplasm(?:s)?| malignanc(?:y|ies)| "
            r"disease(?:s)?)?|blood cancer(?:s)?|leukemia(?:s)?|lymphoma(?:s)?|"
            r"loss of Y(?: chromosome)?|LOY|myelodysplastic(?:\s+syndrome)?|MDS)\b",
            re.I,
        ),
    ),
]


# "syndromic genes/forms/features" is not a disease name.
_SYNDROMIC_WEAK_FOLLOWERS = frozenset(
    {
        "gene",
        "genes",
        "form",
        "forms",
        "feature",
        "features",
        "patient",
        "patients",
        "case",
        "cases",
        "spectrum",
        "presentation",
        "presentations",
        "phenotype",
        "phenotypes",
        "disease",
        "diseases",
        "disorder",
        "disorders",
    }
)


def is_generic_clinical_text(text: str | None) -> bool:
    """True if clinical value is only bare heads like 'syndrome' / 'disorder'."""
    if not text or not str(text).strip():
        return True
    parts = [p.strip() for p in re.split(r"\s*;\s*", str(text)) if p.strip()]
    if not parts:
        return True

    def _one(part: str) -> bool:
        toks = [t for t in re.split(r"[\s/\-‐‑–—]+", part.strip()) if t]
        if not toks:
            return True
        if len(toks) == 1:
            return toks[0].lower() in GENERIC_CLINICAL_TOKENS
        if toks[0].lower() == "syndromic":
            followers = [t.lower() for t in toks[1:]]
            if not followers or all(f in _SYNDROMIC_WEAK_FOLLOWERS for f in followers):
                return True
        content = [
            t
            for t in toks
            if t.lower() not in _CLINICAL_PREFIX_STOPWORDS
            and t.lower() not in GENERIC_CLINICAL_TOKENS
        ]
        if not content:
            return True
        # inheritance adjectives alone do not name a disease
        if all(t.lower() in _INHERITANCE_ONLY_TOKENS for t in content):
            return True
        return False

    return all(_one(p) for p in parts)


def normalize_named_disease_span(text: str) -> str | None:
    """
    Trim filler around a named-disease regex match.
    Walks left from the disease head and keeps content tokens until a stopword.
    Returns None for bare / inheritance-only heads (e.g. 'autosomal recessive disorder').
    """
    raw = re.sub(r"\s+", " ", (text or "").strip())
    if not raw:
        return None
    toks = [t for t in raw.split(" ") if t]
    if not toks:
        return None

    # Syndromic phrases: "syndromic" + following content until stopword
    if toks[0].lower() == "syndromic":
        kept = ["syndromic"]
        for t in toks[1:]:
            tl = t.lower().strip(",.;:")
            if tl in _CLINICAL_PREFIX_STOPWORDS or tl in {"associated", "with"}:
                break
            kept.append(t)
        cleaned = " ".join(kept)
        return cleaned if not is_generic_clinical_text(cleaned) and len(kept) >= 2 else None

    # Locate disease head (and optional "type N" / trailing roman/digit).
    head_idx = None
    end_idx = len(toks) - 1
    for i in range(len(toks) - 1, -1, -1):
        tl = toks[i].lower().strip(",.;:")
        if tl in _NAMED_DISEASE_HEADS_SET:
            head_idx = i
            end_idx = i
            # Optional type suffix already included after head in toks
            if i + 1 < len(toks) and toks[i + 1].lower() == "type":
                end_idx = i + 1
                if i + 2 < len(toks) and re.fullmatch(r"[IVXLC\d]+", toks[i + 2], re.I):
                    end_idx = i + 2
            elif i + 1 < len(toks) and re.fullmatch(r"[IVXLC\d]+", toks[i + 1], re.I):
                end_idx = i + 1
            break
    if head_idx is None:
        # Not a head-based phrase; keep if already specific (e.g. Acanthosis Nigricans via other pat)
        cleaned = " ".join(toks)
        return cleaned if not is_generic_clinical_text(cleaned) else None

    content: list[str] = []
    j = head_idx - 1
    while j >= 0:
        tl = toks[j].lower().strip(",.;:")
        if tl in _CLINICAL_PREFIX_STOPWORDS:
            break
        content.append(toks[j])
        j -= 1
        if len(content) >= 5:
            break
    content.reverse()
    if not content:
        return None
    if all(t.lower().strip(",.;:") in _INHERITANCE_ONLY_TOKENS for t in content):
        return None
    phrase = content + toks[head_idx : end_idx + 1]
    cleaned = " ".join(phrase)
    if is_generic_clinical_text(cleaned):
        return None
    return cleaned

# Entities usable as clinical field values (not discourse cues).
ENTITY_CLINICAL_SUBTYPES = {
    "hpo_like",
    "syndrome_disease",
    "phenotype_term",
    "hematology_oncology",
    "pubtator_disease",
}

CUE_CLINICAL_SUBTYPES = {"clinical_cue", "clinical_phrase"}

# Strong enough to keep a link / raise tier when paired with a variant.
STRONG_CLINICAL_SUBTYPES = ENTITY_CLINICAL_SUBTYPES | CUE_CLINICAL_SUBTYPES

CONCRETE_VARIANT_SUBTYPES = {
    "hgvs_c",
    "hgvs_p_3letter",
    "hgvs_p_1letter",
    "hgvs_p_paren",
    "rsid",
    "nm_np_with_change",
    "ivs",
}

SENTENCE_SPLIT = re.compile(r"(?<=[.!?])\s+(?=[A-Z0-9\"'(\[])")

# BioC section_type values treated as body narrative (exclude refs)
BODY_SECTION_TYPES = {
    "TITLE",
    "ABSTRACT",
    "INTRO",
    "METHODS",
    "RESULTS",
    "DISCUSS",
    "CONCL",
    "CASE",
    "SUPPL",
    "OTHER",
}
SKIP_SECTION_TYPES = {"REF", "AUTH_CONT", "COMP_INT"}

# Nearby gene mention patterns (symbol-like).
# NOTE: do not use re.I on the gene-capture group path that matches "variant in X",
# otherwise place names like "Japan" become false genes.
# NOTE: do not treat the alt allele in "c.179G>A p.(…)" as a gene (A before p.).
_GENE_NEAR = [
    # GENE gene
    re.compile(r"\b([A-Z][A-Z0-9]{1,14}(?:-[A-Z0-9]+)*)\s+gene\b"),
    # GENE variant / GENE mutation (e.g. "MAX variant, c.179G>A")
    re.compile(r"\b([A-Z][A-Z0-9]{1,14}(?:-[A-Z0-9]+)*)\s+(?:variants?|mutations?)\b"),
    # GENE immediately before HGVS / rs (min 2 chars; not a base after '>')
    re.compile(
        r"(?<![>])\b([A-Z][A-Z0-9]{1,14}(?:-[A-Z0-9]+)*)\s+"
        r"(?:p\.\(?|[cp]\.|rs\d{4,}|NM_\d+)",
    ),
    # mutation/variant in GENE  (gene token must be ALL-CAPS / digit style, not Title Case)
    re.compile(
        r"\b(?:mutations?|variants?|pathogenic\s+variants?)\s+in\s+(?:the\s+)?"
        r"([A-Z][A-Z0-9]{1,14}(?:-[A-Z0-9]+)*)\b"
    ),
    re.compile(r"\bin\s+(?:the\s+)?([A-Z][A-Z0-9]{1,14}(?:-[A-Z0-9]+)*)\s+gene\b"),
    # HGVS / rs … in GENE  (post-variant; e.g. "c.206G>T (p.Arg69Leu) in EDA")
    re.compile(
        r"(?:[cp]\.\(?[A-Za-z0-9_*>+=.-]+\)?|rs\d{4,})"
        r"(?:\s*\([^)]{0,40}\))?"
        r"\s+in\s+(?:the\s+)?"
        r"([A-Z][A-Z0-9]{1,14}(?:-[A-Z0-9]+)*)\b"
    ),
]

# Generic words / places / study nouns that look like gene symbols but are not.
_GENE_STOP = {
    "THE",
    "AND",
    "FOR",
    "WITH",
    "FROM",
    "THIS",
    "THAT",
    "GENETIC",
    "NOVEL",
    "PATIENT",
    "PATIENTS",
    "DNA",
    "RNA",
    "PCR",
    "NGS",
    "WES",
    "WGS",
    "ACMG",
    "OMIM",
    "HPO",
    "SNP",
    "SNPS",
    "CNV",
    "CNVS",
    "LOH",
    "VUS",
    "SIRS",
    "AN",
    "PCOS",
    "AML",
    "MDS",
    "MPN",
    "CHIP",
    "LOY",
    # study / anatomy nouns often after "variant in …"
    "EXON",
    "EXONS",
    "INTRON",
    "INTRONS",
    "CODON",
    "DOMAIN",
    "REGION",
    "FAMILY",
    "FAMILIES",
    "POPULATION",
    "COHORT",
    "STUDY",
    "CHILDREN",
    "ADULTS",
    "MALES",
    "FEMALES",
    "SIBLINGS",
    "PROBAND",
    "PROBANDS",
    "CONTROLS",
    "CASES",
    "TABLE",
    "FIGURE",
    "SUPPLEMENT",
    "EUROPE",
    "ASIA",
    "AFRICA",
    # BioC / paper section labels mistaken for genes
    "ABSTRACT",
    "INTRO",
    "METHODS",
    "RESULTS",
    "DISCUSS",
    "DISCUSSION",
    "CONCL",
    "CONCLUSION",
    "CONCLUSIONS",
    "CASE",
    "SUPPL",
    "TITLE",
    "AUTH",
    "COMP",
    # single-letter / base tokens (not HGNC coding symbols)
    "A",
    "T",
    "G",
    "C",
    "N",
    "U",
    # geographic / demonym false positives ("variant in Japan")
    "JAPAN",
    "JAPANESE",
    "CHINA",
    "CHINESE",
    "KOREA",
    "KOREAN",
    "INDIA",
    "INDIAN",
    "BRAZIL",
    "BRAZILIAN",
    "FRANCE",
    "FRENCH",
    "GERMANY",
    "GERMAN",
    "ITALY",
    "ITALIAN",
    "SPAIN",
    "SPANISH",
    "CANADA",
    "CANADIAN",
    "AUSTRALIA",
    "AUSTRALIAN",
    "TAIWAN",
    "TAIWANESE",
    "TURKEY",
    "TURKISH",
    "IRAN",
    "IRANIAN",
    "IRAQ",
    "EGYPT",
    "EGYPTIAN",
    "MEXICO",
    "MEXICAN",
    "RUSSIA",
    "RUSSIAN",
    "SWEDEN",
    "SWEDISH",
    "NORWAY",
    "NORWEGIAN",
    "POLAND",
    "POLISH",
    "GREECE",
    "GREEK",
    "ISRAEL",
    "ISRAELI",
    "USA",
    "UK",
    "ENGLAND",
    "BRITAIN",
    "BRITISH",
    "AMERICA",
    "AMERICAN",
    "EUROPEAN",
    "ASIAN",
    "AFRICAN",
    # ACMG / clinical shorthand mistaken for genes (negatives: WT/ROS/HPP/PM2…)
    "WT",
    "WILDTYPE",
    "WILD",
    "ROS",
    "HPP",
    "HMNX",
    "VUS",
    "LB",
    "LP",
    "P",
    "B",
    "PM1",
    "PM2",
    "PM3",
    "PM4",
    "PM5",
    "PM6",
    "PP1",
    "PP2",
    "PP3",
    "PP4",
    "PP5",
    "PS1",
    "PS2",
    "PS3",
    "PS4",
    "PVS1",
    "BA1",
    "BS1",
    "BS2",
    "BS3",
    "BS4",
    "BP1",
    "BP2",
    "BP3",
    "BP4",
    "BP5",
    "BP6",
    "BP7",
    "ACMG",
    "AMP",
    "CLINSIG",
    "BENIGN",
    "LIKELY",
    "PATHOGENIC",
    "CONTROL",
    "CONTROLS",
    "NEGATIVE",
    "POSITIVE",
}

# 3-letter → 1-letter amino acids for p.HGVS matching
_AA3_TO1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "TER": "*",
    "STOP": "*",
}


def normalize_p_hgvs_for_match(text: str) -> str:
    """Normalize p.HGVS for comparison: strip parens, 3-letter→1-letter AA codes."""
    t = re.sub(r"\s+", "", (text or "").strip())
    if not t.lower().startswith("p."):
        return t.lower()
    body = t[2:]
    body = body.strip("()")
    # Digits are word chars, so do not use \\b; bound AA codes by non-letters.
    def repl_aa(m: re.Match[str]) -> str:
        aa = m.group(1).upper()
        return _AA3_TO1.get(aa, m.group(1))

    body = re.sub(r"(?<![A-Za-z])([A-Za-z]{3})(?![A-Za-z])", repl_aa, body)
    return ("p." + body).lower()


def is_plausible_gene_symbol(sym: str | None) -> bool:
    """Reject disease phrases, protein-change-like tokens, places, and stopwords."""
    if not sym:
        return False
    g = str(sym).strip()
    gu = g.upper()
    if gu in _GENE_STOP or gu.isdigit():
        return False
    # Require at least 2 chars (reject lone A/T/G/C and other 1-letter tokens)
    if len(gu) < 2:
        return False
    # Title Case place/person style (Japan, China) — not HGNC-style symbols
    if len(g) >= 4 and g[0].isupper() and g[1:].islower() and g.isalpha():
        return False
    # protein change mistaken as gene: P22H, E676K, G641E
    if re.fullmatch(r"[A-Z]\d+[A-Z*]", gu):
        return False
    if re.fullmatch(r"[A-Z]{3}\d+[A-Z]{3}", gu):
        return False
    # disease compounds: WEILL-MARCHESANI, MULTIDRUG-RESISTANT
    if g.count("-") >= 1 and len(g) > 12:
        return False
    # Allow short symbols like C3 / F2 (letter + digits) and standard A-Z symbols
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9]{0,14}(?:-[A-Za-z0-9]{1,6})?", g):
        return False
    return True


def gene_adjacent_to_variant(text: str | None, variant: str | None) -> str | None:
    """Return gene symbol immediately before/after a concrete variant mention, if any.

    Left:  ``GENE c.179G>A`` / ``GENE p.Arg60Gln``
    Right: ``c.206G>T (p.Arg69Leu) in EDA`` / ``p.E840V in ATP7A``
    """
    if not text or not variant:
        return None
    v = str(variant).strip()
    if not v:
        return None
    # Prefer exact / relaxed forms of the variant in text
    forms = [v]
    v_rel = re.sub(r"[\(\)\s]", "", v)
    if v_rel and v_rel not in forms:
        forms.append(v_rel)
    for form in forms:
        # look behind up to ~40 chars for ALL-CAPS gene token
        for m in re.finditer(re.escape(form), text, flags=re.I):
            pre = text[max(0, m.start() - 40) : m.start()]
            # Skip alt allele of c.HGVS glued before p.HGVS: "...c.179G>A "
            if re.search(r"[cC]\.\d+[ACGT]>[ACGT]\s*$", pre):
                continue
            if re.search(r">[ACGT]\s*$", pre):
                continue
            gm = re.search(
                r"\b([A-Z][A-Z0-9]{1,14}(?:-[A-Z0-9]+)*)\s*$",
                pre,
            )
            if not gm:
                continue
            sym = gm.group(1)
            # Reject if the char before the symbol is '>' (allele, not gene)
            abs_start = m.start() - len(pre) + gm.start(1)
            if abs_start > 0 and text[abs_start - 1] == ">":
                continue
            if is_plausible_gene_symbol(sym):
                return sym.upper()
    # Right side: variant … in GENE (optional parenthetical protein between)
    for form in forms:
        for m in re.finditer(re.escape(form), text, flags=re.I):
            post = text[m.end() : m.end() + 90]
            gm = re.search(
                r"^(?:\s*\([^)]{0,50}\))?"
                r"\s*(?:,|\band\b)?"
                r"\s*\bin\s+(?:the\s+)?"
                r"([A-Z][A-Z0-9]{1,14}(?:-[A-Z0-9]+)*)\b",
                post,
            )
            if gm and is_plausible_gene_symbol(gm.group(1)):
                return gm.group(1).upper()
    return None


# English / figure words that look like gene symbols in "… GENE locus/gene" prose.
_LOCUS_PHRASE_STOP = frozenset(
    {
        "THE",
        "AND",
        "FOR",
        "WITH",
        "FROM",
        "THIS",
        "THAT",
        "EACH",
        "BOTH",
        "SAME",
        "OTHER",
        "HUMAN",
        "MOUSE",
        "GENE",
        "GENES",
        "LOCUS",
        "REGION",
        "CHROMOSOME",
        "AUTOPHAGY",
        "TRANSMEMBRANE",
        "SUPPRESSOR",
        "ADIPOCYTE",
        "CONSERVED",
        "LINKER",
        "GAMMA",
        "ALPHA",
        "BETA",
        "DELTA",
        "KINASE",
        "RECEPTOR",
        "PROTEIN",
        "DOMAIN",
        "NUCLEAR",
        "CODING",
        "NONCODING",
        "INTERGENIC",
        "PROMOTER",
        "INTRONIC",
        "EXONIC",
        "CANDIDATE",
        "NOVEL",
        "RISK",
        "LEAD",
        "FIGURE",
        "TABLE",
        "SUPPLEMENTARY",
        "POLYMORPHIC",
        "PREVALENT",
        "COMMON",
        "ASSOCIATED",
        "SIGNIFICANT",
        "PUTATIVE",
        "SUSCEPTIBILITY",
    }
)


NEIGHBOR_SENTENCE_RADIUS = 2

_NAMED_PAREN_GENE = re.compile(
    r"\b(?:named|termed|called|designated|denoted)\b.{0,240}?\(\s*([A-Z][A-Z0-9]{2,20})\s*\)",
    re.I | re.S,
)
_LNCRNA_PAREN_GENE = re.compile(
    r"\b(?:lncRNA|lincRNA|long\s+non[-\s]?coding\s+RNA)\b.{0,80}?\(\s*([A-Z][A-Z0-9]{2,20})\s*\)",
    re.I | re.S,
)
_NAMED_TOKEN_GENE = re.compile(
    r"\b(?:named|termed|called|designated|denoted)\s+"
    r"(?:(?:this|the)\s+(?:lncRNA|gene|transcript|variant)\s+)?"
    r"(?:as\s+)?"
    r"([A-Z][A-Z0-9]{2,20})\b",
    re.I,
)
_GENE_THEN_VARIANT_WORD = re.compile(
    r"\b([A-Z][A-Z0-9]{1,14}(?:-[A-Z0-9]+)*)\s+(?i:variant|mutation|snp)\b"
)


def neighbor_text_for_variant(
    text: str | None,
    variant: str | None,
    radius: int = NEIGHBOR_SENTENCE_RADIUS,
) -> str:
    """Return the ±radius sentence window around the first mention of variant.

    If the variant is not found, return empty so distant article genes are not used.
    """
    if not text:
        return ""
    blob = str(text)
    v = str(variant or "").strip()
    sents = split_sentences(blob)
    if not sents:
        return blob if v and re.search(re.escape(v), blob, flags=re.I) else ""
    if not v:
        return ""
    idx = None
    for i, (_a, _b, sent) in enumerate(sents):
        if re.search(re.escape(v), sent, flags=re.I):
            idx = i
            break
    if idx is None:
        return ""
    lo = max(0, idx - radius)
    hi = min(len(sents), idx + radius + 1)
    return " ".join(s[2] for s in sents[lo:hi])


def gene_from_naming_cues(text: str | None) -> str | None:
    """Gene from 'named/termed … (SYMBOL)' or a bare named SYMBOL (not a phrase)."""
    if not text:
        return None
    for pat in (_NAMED_PAREN_GENE, _LNCRNA_PAREN_GENE):
        m = pat.search(text)
        if m and is_plausible_gene_symbol(m.group(1)):
            return m.group(1).upper()
    for m in _NAMED_TOKEN_GENE.finditer(text):
        rest = text[m.end() : m.end() + 24]
        if re.match(r"\s+[a-z]", rest):
            continue
        if is_plausible_gene_symbol(m.group(1)):
            return m.group(1).upper()
    return None


def gene_before_variant_word(text: str | None) -> str | None:
    """Left-adjacent SYMBOL before 'variant/mutation/SNP' (neighbor-sentence prose)."""
    if not text:
        return None
    for m in _GENE_THEN_VARIANT_WORD.finditer(text):
        sym = m.group(1)
        if is_plausible_gene_symbol(sym) and str(sym).upper() not in _GENE_STOP:
            return str(sym).upper()
    return None


def gene_from_neighbor_window(text: str | None, variant: str | None) -> str | None:
    """Bind gene from a ±2 sentence rule window: naming cues, then SYMBOL variant."""
    win = neighbor_text_for_variant(text, variant)
    if not win:
        return None
    named = gene_from_naming_cues(win)
    if named:
        return named
    before = gene_before_variant_word(win)
    if before:
        return before
    loc = gene_in_locus_near_variant(win, variant)
    if loc:
        return loc
    return gene_adjacent_to_variant(win, variant)


def gene_in_locus_near_variant(text: str | None, variant: str | None) -> str | None:
    """Return gene from phrases like 'rsX, a variant located in the GENE locus'.

    Only the variant→locus direction is used. The reverse 'GENE locus … rs'
    pattern is skipped because prose like 'polymorphic locus … rsX' false-hits.
    """
    if not text or not variant:
        return None
    v = str(variant).strip()
    if not v:
        return None
    v_esc = re.escape(v)
    gene_tok = r"([A-Z][A-Z0-9]{1,14}(?:-[A-Z0-9]+)*)"
    # rs112720315, a variant located in the TCP10L2 locus on Chromosome 6
    after = re.compile(
        rf"(?:{v_esc})(.{{0,160}}?)\b(?:located\s+in|located\s+at|within|mapped\s+to)"
        rf"\s+(?:the\s+)?{gene_tok}\s+locus\b",
        re.I | re.S,
    )
    for m in after.finditer(text):
        mid = m.group(1) or ""
        if mid.count(".") > 1:
            continue
        sym = (m.group(2) or "").upper()
        if sym in _LOCUS_PHRASE_STOP or sym in _GENE_STOP:
            continue
        # Prefer symbols that look like gene IDs (digit / hyphen), else allow.
        if is_plausible_gene_symbol(sym):
            return sym
    return None


def extract_gene_mentions(text: str | None) -> list[str]:
    if not text:
        return []
    out: list[str] = []
    seen: set[str] = set()
    for pat in _GENE_NEAR:
        for m in pat.finditer(text):
            sym = (m.group(1) or "").upper()
            if not is_plausible_gene_symbol(sym) or sym in seen:
                continue
            seen.add(sym)
            out.append(sym)
    return out


def is_rsid_token(variant: str | None) -> bool:
    if not variant:
        return False
    tok = normalize_variant_token(str(variant))
    return bool(re.match(r"^rs\d{4,}$", tok, re.I))


def gene_from_pubtator_variant(
    variant: str | None, pubtator: dict[str, Any] | None
) -> str | None:
    """RefSNP/HGVS → PubTator CorrespondingGene symbol (catalog locus gene)."""
    pt = pubtator or {}
    if not variant:
        return None
    for v in pt.get("variants") or []:
        hv = v.get("hgvs") or v.get("name") or v.get("text") or ""
        if not variant_tokens_match(str(variant), str(hv)):
            continue
        gid = v.get("gene_id")
        if gid is None:
            ident = str(v.get("identifier") or "")
            m = re.search(r"CorrespondingGene:(\d+)", ident)
            if m:
                gid = m.group(1)
        sym = gene_id_to_symbol(gid, pt)
        if sym:
            return sym
    return None


def resolve_gene_for_variant(
    variant: str | None,
    *,
    evidence_sentence: str | None = None,
    title: str | None = None,
    abstract: str | None = None,
    context_text: str | None = None,
    pubtator: dict[str, Any] | None = None,
    link_gene: str | None = None,
) -> str | None:
    """
    Bind a gene symbol to a variant without falling back to article genes[0].

    Priority:
      1) HGVS: gene token immediately adjacent in evidence/title
      2) locus phrase (variant → located in GENE locus)
      3) rsID: PubTator CorrespondingGene / dbSNP catalog gene. Skip adjacency
         and ±2 sentence naming ("we named HOTSCRAMBL"; "HOTSCRAMBL rs…").
      4) HGVS: ±2 sentence naming
      5) adjacency / locus on abstract or long full-text context (non-rs)
      6) explicit link_gene
      7) PubTator variant → gene (non-rs leftover)
      8) nearby gene mentions
    """
    blobs = (evidence_sentence, title, abstract, context_text)
    rs = is_rsid_token(variant)

    # rsIDs: papers often write "HOTSCRAMBL rs17437411" as an experimental allele
    # nickname. That is not the RefSNP locus gene (dbSNP/PubTator).
    if not rs:
        for blob in (evidence_sentence, title):
            adj = gene_adjacent_to_variant(blob, variant)
            if adj:
                return adj

    for blob in (evidence_sentence, title, abstract):
        loc = gene_in_locus_near_variant(blob, variant)
        if loc:
            return loc

    if rs:
        pt_gene = gene_from_pubtator_variant(variant, pubtator)
        if pt_gene:
            return pt_gene
        if link_gene and is_plausible_gene_symbol(str(link_gene).strip()):
            g = str(link_gene).strip()
            return g.upper() if g.isupper() or len(g) <= 8 else g
        return None

    # Prefer paper naming in any ±2 window before weaker neighbor heuristics
    # (evidence-only sentences must not bind "genetic variant").
    for blob in blobs:
        named = gene_from_naming_cues(neighbor_text_for_variant(blob, variant))
        if named:
            return named
    for blob in blobs:
        neigh = gene_from_neighbor_window(blob, variant)
        if neigh:
            return neigh

    for blob in (abstract, context_text):
        adj = gene_adjacent_to_variant(blob, variant)
        if adj:
            return adj
    for blob in (context_text,):
        loc = gene_in_locus_near_variant(blob, variant)
        if loc:
            return loc

    if link_gene and is_plausible_gene_symbol(str(link_gene).strip()):
        g = str(link_gene).strip()
        return g.upper() if g.isupper() or len(g) <= 8 else g

    pt_gene = gene_from_pubtator_variant(variant, pubtator)
    if pt_gene:
        return pt_gene

    pt = pubtator or {}
    for blob in blobs:
        mentions = extract_gene_mentions(blob)
        if not mentions:
            continue
        pt_syms = {
            (g.get("name") or "").upper()
            for g in (pt.get("genes") or [])
            if is_plausible_gene_symbol(g.get("name"))
        }
        for m in mentions:
            if not pt_syms or m in pt_syms:
                return m
        return mentions[0]
    return None

def variant_tokens_match(a: str, b: str) -> bool:
    """Loose match for HGVS/rsID strings (ignore spaces; p. 1↔3 letter; allow prefix)."""
    na = re.sub(r"\s+", "", normalize_variant_token(a)).lower()
    nb = re.sub(r"\s+", "", normalize_variant_token(b)).lower()
    if not na or not nb:
        return False
    if na == nb:
        return True
    if na.startswith("rs") and nb.startswith("rs"):
        return na == nb
    # c.2663A vs c.2663A>G
    if na.startswith("c.") and nb.startswith("c."):
        return na in nb or nb in na
    if na.startswith("p.") and nb.startswith("p."):
        pa, pb = normalize_p_hgvs_for_match(na), normalize_p_hgvs_for_match(nb)
        if pa == pb:
            return True
        return pa in pb or pb in pa
    return False



@dataclass
class Span:
    kind: str
    subtype: str
    text: str
    start: int
    end: int
    source: str
    normalized: dict[str, Any] = field(default_factory=dict)


def split_sentences(text: str) -> list[tuple[int, int, str]]:
    parts: list[tuple[int, int, str]] = []
    start = 0
    for m in SENTENCE_SPLIT.finditer(text):
        end = m.start()
        sent = text[start:end].strip()
        if sent:
            s0 = text.find(sent, start, end)
            if s0 < 0:
                s0 = start
            parts.append((s0, s0 + len(sent), sent))
        start = m.end()
    tail = text[start:].strip()
    if tail:
        s0 = text.find(tail, start)
        if s0 < 0:
            s0 = start
        parts.append((s0, s0 + len(tail), tail))
    return parts


def normalize_hgvs_text(text: str) -> str:
    """Collapse whitespace around HGVS allele substitutions (c.2663A > G → c.2663A>G)."""
    t = (text or "").strip()
    t = re.sub(r"([ACGT])\s*>\s*([ACGT])", r"\1>\2", t, flags=re.I)
    return t


def find_spans(text: str, patterns: list[tuple[str, re.Pattern[str]]], kind: str, source: str) -> list[Span]:
    spans: list[Span] = []
    seen: set[tuple[int, int, str]] = set()
    for subtype, pat in patterns:
        for m in pat.finditer(text):
            raw = m.group(0)
            start, end = m.start(), m.end()
            if kind == "variant":
                norm = normalize_hgvs_text(raw)
            elif kind == "clinical" and subtype == "syndrome_disease":
                cleaned = normalize_named_disease_span(raw)
                if cleaned is None:
                    continue
                # Re-anchor start if leading stopwords were stripped from the match text.
                if cleaned.lower() != raw.lower():
                    idx = raw.lower().rfind(cleaned.lower())
                    if idx >= 0:
                        start = start + idx
                        end = start + len(cleaned)
                        norm = cleaned
                    else:
                        norm = cleaned
                else:
                    norm = cleaned
            else:
                norm = raw
                if kind == "clinical" and subtype in ENTITY_CLINICAL_SUBTYPES and is_generic_clinical_text(norm):
                    continue

            key = (start, end, norm)
            if key in seen:
                continue
            # Prefer longer overlapping entity spans of the same kind
            if kind == "clinical" and any(
                s.start <= start and s.end >= end and s.subtype in ENTITY_CLINICAL_SUBTYPES
                for s in spans
            ):
                continue
            # Drop shorter entity spans fully covered by this longer match
            if kind == "clinical" and subtype in ENTITY_CLINICAL_SUBTYPES:
                spans = [
                    s
                    for s in spans
                    if not (
                        s.subtype in ENTITY_CLINICAL_SUBTYPES
                        and start <= s.start
                        and end >= s.end
                        and (end - start) > (s.end - s.start)
                    )
                ]
            seen.add(key)
            spans.append(
                Span(
                    kind=kind,
                    subtype=subtype,
                    text=norm,
                    start=start,
                    end=end,
                    source=source,
                )
            )
    cleaned: list[Span] = []
    for s in sorted(spans, key=lambda x: (x.start, -(x.end - x.start), x.subtype)):
        if s.kind == "clinical" and s.subtype in CUE_CLINICAL_SUBTYPES:
            if any(
                e.subtype in ENTITY_CLINICAL_SUBTYPES and e.start <= s.start and e.end >= s.end
                for e in cleaned
            ):
                continue
        cleaned.append(s)
    cleaned.sort(key=lambda s: (s.start, s.end))
    return cleaned


def _span_as_dict(span: Span | dict[str, Any] | Any) -> dict[str, Any]:
    if isinstance(span, dict):
        return dict(span)
    if isinstance(span, Span):
        return asdict(span)
    # Foreign Span-like objects (e.g. link_variant_phenotype.Span)
    return {
        "kind": getattr(span, "kind", None),
        "subtype": getattr(span, "subtype", None),
        "text": getattr(span, "text", None),
        "start": getattr(span, "start", None),
        "end": getattr(span, "end", None),
        "source": getattr(span, "source", None),
    }


def pick_best_clinical_span(
    spans: list[Span] | list[dict[str, Any]],
    *,
    allow_cue: bool = True,
    join_entities: bool = True,
    max_entities: int = 5,
) -> dict[str, Any] | None:
    """Prefer disease/phenotype entities over discourse cues; optionally join entities."""
    if not spans:
        return None
    dicts = [_span_as_dict(s) for s in spans]
    entities = [
        s
        for s in dicts
        if s.get("subtype") in ENTITY_CLINICAL_SUBTYPES and not is_generic_clinical_text(s.get("text"))
    ]
    if entities:
        # Prefer longer / more specific first, then keep unique texts in order of appearance
        entities_sorted = sorted(entities, key=lambda s: (-(s["end"] - s["start"]), s["start"]))
        if join_entities:
            seen: set[str] = set()
            texts: list[str] = []
            # Restore appearance order for the joined string
            for s in sorted(entities, key=lambda x: x["start"]):
                t = (s.get("text") or "").strip()
                key = t.lower()
                if not t or key in seen or is_generic_clinical_text(t):
                    continue
                # Skip a span fully covered by a longer kept span text
                if any(key != k and key in k for k in seen):
                    continue
                seen.add(key)
                texts.append(t)
                if len(texts) >= max_entities:
                    break
            if texts:
                primary = max(entities, key=lambda s: (s["end"] - s["start"], -s["start"]))
                # Prefer primary among kept texts
                for s in entities_sorted:
                    if (s.get("text") or "").strip() in texts:
                        primary = s
                        break
                return {
                    "text": "; ".join(texts),
                    "subtype": primary.get("subtype") or "phenotype_term",
                    "start": primary.get("start"),
                    "end": primary.get("end"),
                }
        return entities_sorted[0]
    if allow_cue:
        cues = [s for s in dicts if s.get("subtype") in CUE_CLINICAL_SUBTYPES]
        if cues:
            return cues[0]
        # Do not fall back to bare generic entity tokens
        return None
    return None


def extract_sentence_links(
    text: str,
    *,
    source: str = "regex",
    method_prefix: str = "",
    section: str | None = None,
    require_concrete_variant: bool = True,
) -> dict[str, Any]:
    """Find variant+clinical same-sentence and same-paragraph style links in text."""
    text = re.sub(r"\s+", " ", (text or "").strip())
    if not text:
        return {
            "text": text,
            "variant_spans": [],
            "clinical_spans": [],
            "same_sentence_links": [],
            "flags": {
                "has_strong_variant_text": False,
                "has_clinical_text": False,
                "n_same_sentence_links": 0,
            },
        }

    var_spans = find_spans(text, VARIANT_PATTERNS, "variant", source)
    for m in SOFT_VARIANT.finditer(text):
        var_spans.append(Span("variant", "soft_variant_phrase", m.group(0), m.start(), m.end(), source))
    clin_spans = find_spans(text, CLINICAL_PATTERNS, "clinical", source)

    sentences = split_sentences(text)
    same_sentence_links: list[dict[str, Any]] = []
    for s0, s1, sent in sentences:
        vs_all = [v for v in var_spans if v.start >= s0 and v.end <= s1 and v.subtype != "soft_variant_phrase"]
        cs_all = [c for c in clin_spans if c.start >= s0 and c.end <= s1]
        if not vs_all or not cs_all:
            continue
        vs_concrete = [v for v in vs_all if v.subtype in CONCRETE_VARIANT_SUBTYPES]
        cs_entity = [c for c in cs_all if c.subtype in ENTITY_CLINICAL_SUBTYPES]
        cs_strong = [c for c in cs_all if c.subtype in STRONG_CLINICAL_SUBTYPES]
        vs_path = [v for v in vs_all if v.subtype == "pathogenic_phrase"]
        pick = pick_best_clinical_span(cs_all, allow_cue=True, join_entities=True)
        if not pick:
            continue
        if vs_concrete and cs_all:
            v = vs_concrete[0]
            c_text, c_sub = pick["text"], pick["subtype"]
        elif (not require_concrete_variant) and vs_path and (cs_entity or cs_strong):
            v = vs_path[0]
            c_text, c_sub = pick["text"], pick["subtype"]
        elif vs_path and (cs_entity or cs_strong):
            v = vs_path[0]
            c_text, c_sub = pick["text"], pick["subtype"]
        else:
            continue
        link = {
            "variant": v.text,
            "variant_subtype": v.subtype,
            "clinical": c_text,
            "clinical_subtype": c_sub,
            "evidence_sentence": sent,
            "method": f"{method_prefix}same_sentence_cooccurrence" if method_prefix else "same_sentence_cooccurrence",
            "tier": "A" if c_sub in ENTITY_CLINICAL_SUBTYPES else "B",
            "text_scope": "fulltext" if method_prefix.startswith("fulltext") else "abstract",
            "clinical_source": "lexicon_entity" if c_sub in ENTITY_CLINICAL_SUBTYPES else "lexicon_cue",
        }
        if section:
            link["section"] = section
        same_sentence_links.append(link)

    strong = any(s.subtype in CONCRETE_VARIANT_SUBTYPES or s.subtype == "pathogenic_phrase" for s in var_spans)
    return {
        "text": text,
        "variant_spans": [asdict(s) for s in var_spans],
        "clinical_spans": [asdict(s) for s in clin_spans],
        "same_sentence_links": same_sentence_links,
        "flags": {
            "has_strong_variant_text": strong,
            "has_clinical_text": len(clin_spans) > 0,
            "n_same_sentence_links": len(same_sentence_links),
        },
    }


def normalize_variant_token(text: str) -> str:
    t = normalize_hgvs_text(text or "")
    if re.match(r"^rs\d+$", t, re.I):
        return t.lower()
    return t


def gene_id_to_symbol(gene_id: Any, pubtator: dict[str, Any] | None) -> str | None:
    if gene_id is None:
        return None
    gid = str(gene_id).strip()
    pt = pubtator or {}
    for g in pt.get("genes") or []:
        ids = {
            str(g.get("identifier") or "").strip(),
            str(g.get("normalized_id") or "").strip(),
            str(g.get("gene_id") or "").strip(),
        }
        if gid not in ids:
            continue
        for cand in ((g.get("name") or "").strip(), (g.get("text") or "").strip()):
            if is_plausible_gene_symbol(cand):
                return cand.upper() if cand.isupper() or len(cand) <= 8 else cand
    return None


def scan_concrete_variants(text: str) -> list[str]:
    """Collect unique concrete variant strings (rsID/HGVS) from free text."""
    if not text:
        return []
    extracted = extract_sentence_links(text, require_concrete_variant=True)
    seen: set[str] = set()
    out: list[str] = []
    for span in extracted["variant_spans"]:
        if span["subtype"] not in CONCRETE_VARIANT_SUBTYPES:
            continue
        tok = normalize_variant_token(span["text"])
        if not tok or not is_concrete_variant_string(tok) or tok in seen:
            continue
        seen.add(tok)
        out.append(tok)
    return out


def scan_passages_for_concrete_variants(passages: list[dict[str, Any]]) -> list[str]:
    """Scan full-text passages for concrete variants not limited to abstract seeds."""
    seen: set[str] = set()
    out: list[str] = []
    for p in passages:
        text = p.get("text") or ""
        for tok in scan_concrete_variants(text):
            if tok in seen:
                continue
            seen.add(tok)
            out.append(tok)
    return out


def _is_weak_clinical(clin: str, clin_sub: str) -> bool:
    if not clin or clin == "fulltext_context":
        return True
    if clin_sub in CUE_CLINICAL_SUBTYPES:
        return True
    if clin_sub == "narrative_context":
        return True
    if is_generic_clinical_text(clin):
        return True
    return False


def upgrade_clinical_field(
    link: dict[str, Any],
    *,
    paragraph_text: str | None = None,
    title: str | None = None,
    pubtator: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Upgrade weak clinical placeholders / discourse cues using lexicon entities,
    same-paragraph context, title, or PubTator disease annotations.
    """
    out = dict(link)
    clin_sub = out.get("clinical_subtype") or ""
    clin = out.get("clinical") or ""

    if (
        clin_sub in ENTITY_CLINICAL_SUBTYPES
        and clin
        and clin != "fulltext_context"
        and not is_generic_clinical_text(clin)
    ):
        out.setdefault("clinical_source", "lexicon_entity")
        return out

    def _apply_from_spans(spans: list[dict[str, Any]], source: str) -> bool:
        pick = pick_best_clinical_span(spans, allow_cue=False, join_entities=True)
        if not pick:
            return False
        out["clinical"] = pick["text"]
        out["clinical_subtype"] = pick["subtype"]
        out["clinical_source"] = source
        if pick["subtype"] in ENTITY_CLINICAL_SUBTYPES:
            out["tier"] = "A" if out.get("tier") == "A" or source.endswith("sentence") else out.get("tier") or "B"
            if source in {"lexicon_sentence", "lexicon_title"}:
                out["tier"] = "A"
        return True

    sent = out.get("evidence_sentence") or ""
    if sent:
        sent_hit = extract_sentence_links(sent, method_prefix="fulltext_enrich_")
        if _apply_from_spans(sent_hit["clinical_spans"], "lexicon_sentence"):
            return out

    para = paragraph_text or ""
    if para and para != sent:
        para_hit = extract_sentence_links(para, method_prefix="fulltext_enrich_")
        if _apply_from_spans(para_hit["clinical_spans"], "lexicon_paragraph"):
            if out.get("clinical_subtype") in ENTITY_CLINICAL_SUBTYPES:
                out["tier"] = out.get("tier") or "B"
            return out

    if title:
        title_hit = extract_sentence_links(title, method_prefix="fulltext_enrich_")
        if _apply_from_spans(title_hit["clinical_spans"], "lexicon_title"):
            return out

    pt = pubtator or {}
    variant = normalize_variant_token(str(out.get("variant") or "")).lower()
    for rel in pt.get("relations") or []:
        for role in (rel.get("role1") or {}, rel.get("role2") or {}):
            role_name = str(role.get("name") or "")
            if role.get("type") == "Variant" and role_name.lower() == variant:
                for other in (rel.get("role1") or {}, rel.get("role2") or {}):
                    if other.get("type") == "Disease" and other.get("name"):
                        out["clinical"] = other["name"]
                        out["clinical_subtype"] = "pubtator_disease"
                        out["clinical_source"] = "pubtator_relation"
                        out["tier"] = "A"
                        return out

    diseases = pt.get("diseases") or []
    if diseases and _is_weak_clinical(clin, clin_sub):
        body_dis = [
            d
            for d in diseases
            if (d.get("section") or "").lower() not in {"title", "abstract"}
        ]
        pick = body_dis[0] if body_dis else diseases[0]
        label = pick.get("name") or pick.get("text")
        if label:
            out["clinical"] = label
            out["clinical_subtype"] = "pubtator_disease"
            out["clinical_source"] = "pubtator_abstract"
            out["tier"] = "B"
            return out

    if _is_weak_clinical(clin, clin_sub):
        out.setdefault("clinical_source", "narrative_context" if clin == "fulltext_context" else "lexicon_cue")
    else:
        out.setdefault("clinical_source", "lexicon_entity")
    return out


def is_concrete_variant_string(text: str) -> bool:
    if not text:
        return False
    return bool(
        re.search(
            r"(?:^c\.|^p\.|^rs\d|^(?:NM_|NP_|NC_|NG_)|IVS\d|c\.\d|p\.[A-Z])",
            text,
            re.I,
        )
    )
