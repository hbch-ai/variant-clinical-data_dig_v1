# 固定输入（GitHub 与本地一致）

| 文件 | 用途 |
|------|------|
| `cohort/articles.jsonl` | 117 个 gold PMID 的 Stage1 题录（与 `configs/gold_pmids.txt` 一致） |

ClinVar / Ensembl 体积过大，不在 git 中；通过 `scripts/run.py bootstrap` 下载或拷贝到 `data/pubmed/ref/`，并以 `ref/reference_manifest.json` 的 SHA256 校验。

可选：Release 资产 `cohort-fulltext-cache.tar.zst` 解压到 `data/pubmed/cache/fulltext/` 可与某一验收跑次逐字节对齐。
