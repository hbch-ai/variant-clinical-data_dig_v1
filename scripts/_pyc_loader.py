"""Load required sibling bytecode payloads and remap legacy PubMed paths."""
from __future__ import annotations
import importlib.util
import sys
from pathlib import Path
from types import ModuleType
from vc_paths import PUBMED, RemapPubmedRoot, rewrite_legacy_pubmed_path

def _rewrite_value(value):
    if isinstance(value, Path): return rewrite_legacy_pubmed_path(value)
    if isinstance(value, dict): return {k: _rewrite_value(v) for k, v in value.items()}
    if isinstance(value, list): return [_rewrite_value(v) for v in value]
    if isinstance(value, tuple): return tuple(_rewrite_value(v) for v in value)
    return value

def _patch_legacy_pubmed_paths(mod: ModuleType) -> None:
    for name, value in list(vars(mod).items()):
        if name.startswith("_") or name in {"ROOT", "SCRIPTS", "Path"}: continue
        rewritten = _rewrite_value(value)
        if rewritten is not value: setattr(mod, name, rewritten)
    if hasattr(mod, "DATA") and Path(str(mod.DATA)).as_posix().endswith("/data/pubmed"):
        mod.DATA = RemapPubmedRoot(str(PUBMED))

def load_pyc(stem: str) -> ModuleType:
    pyc = Path(__file__).resolve().parent / f"{stem}.pyc"
    if not pyc.is_file(): raise FileNotFoundError(pyc)
    payload_name = f"_vc_payload_{stem}"
    spec = importlib.util.spec_from_file_location(payload_name, pyc)
    if spec is None or spec.loader is None: raise ImportError(f"cannot load {pyc}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[payload_name] = mod
    spec.loader.exec_module(mod)
    _patch_legacy_pubmed_paths(mod)
    return mod
