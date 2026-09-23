#!/usr/bin/env python3
"""Compatibility shim backed by the sibling vc_ensembl_c_map.pyc payload."""
from _pyc_loader import load_pyc
_mod = load_pyc("vc_ensembl_c_map")
globals().update({k: getattr(_mod, k) for k in dir(_mod) if not k.startswith("_")})
if __name__ == "__main__":
    if hasattr(_mod, "main"):
        _mod.main()
    else:
        import argparse
        argparse.ArgumentParser(description=__doc__).parse_args()
