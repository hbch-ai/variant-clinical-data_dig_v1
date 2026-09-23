#!/usr/bin/env python3
"""Compatibility shim backed by the sibling pubmed_pipeline.pyc payload."""
from _pyc_loader import load_pyc
_mod = load_pyc("pubmed_pipeline")
globals().update({k: getattr(_mod, k) for k in dir(_mod) if not k.startswith("_")})
if __name__ == "__main__":
    _mod.main()
