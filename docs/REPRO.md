# 可复现性约定（本地 = GitHub clone）

本 pipeline **无大模型**；Stage2–9 为规则/正则/索引查找。在**相同输入快照**下，每一步 JSONL/CSV 应 **byte-level 一致**（除时间戳字段）。

## 1. 复现三要素

| 要素 | 位置 | 说明 |
|------|------|------|
| **代码版本** | git tag / `configs/repro.json` | `pipeline_version` 与 `gate_config_version` 一起记录 |
| **参考数据** | `data/pubmed/ref/` | ClinVar `variant_summary.txt.gz`（SHA256 写入 `ref/reference_manifest.json`）；Ensembl GFF3/FASTA/FAI **拷入仓库内路径**，禁止指向项目外软链 |
| **运行输入** | cohort：`fixtures/cohort/articles.jsonl` + `configs/gold_pmids.txt` | 117 PMID 题录固定进库；不依赖本机 v3 路径 |

网络阶段（PubMed EFetch、PMC、dbSNP）若现场拉取，**无法保证与历史跑次字节相同**；gold cohort 应用 **fixtures + bootstrap_refs**，Stage2 起不依赖全库检索。

## 2. 标准命令（clone 后与本地相同）

```bash
cd variant_clinical_v4_02
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt

python scripts/run.py bootstrap --link-from /path/to/existing/tree
# Gold cohort（与文档计数一致的可复现输入）
python scripts/run.py cohort \
  --seed-from fixtures/cohort/articles.jsonl \
  --force

python scripts/run.py deliver all
python scripts/run.py compare
```

默认不读取 v3 或其它项目。只有显式同时传入
`--seed-cache --seed-cache-from /path/to/tree` 才会复制开发用全文/PDF cache；
不会隐式合并历史 highconf 行。Stage3 在空 cache 上按 PMC/BioC/JATS 拉取
OA 全文（需网络）。

## 3. 每步产物契约（cohort 主路径）

| Stage | 主 sentinel 文件 |
|-------|------------------|
| 1 | `data/pubmed/ingest/merged/articles.jsonl`（117 行） |
| 2 | `data/pubmed/pipeline/s2_linking/articles_pass.jsonl` |
| 3 | `data/pubmed/pipeline/s2_linking/corpus_enriched.jsonl` |
| 4 | `data/pubmed/pipeline/s4_discovery/corpus_enriched_novel.jsonl` |
| 5 | `data/pubmed/pipeline/s5_clinvar/enriched/…_clinvar_absent.jsonl` |
| 6 | `data/pubmed/pipeline/s6_review/…_review_highconf.jsonl` |
| 7 | `…_journal_annotated.jsonl` |
| 8 | `cache/pdfs/*_pdf_pathogenicity.jsonl` |
| 9 | `data/pubmed/pipeline/s9_genomic/variants_genomic_positive_ok.csv` |

每步 `data/pubmed/runs/` 下可写 `stage_manifest.json`（计划：`run.py` 汇总 sentinel SHA256）。

## 4. 确定性实现要点

- JSONL 写出前对 PMID / locus 键 **排序**（逐步对齐；legacy pyc 段仍按原序）。
- 配置只读 `configs/*.json`，版本号变更必须 bump `configs/repro.json`。
- 禁止代码中写死 `/home/user/...`；仅允许 `VARIANT_CLINICAL_ROOT` 或仓库相对路径。
- 随机性：仅 `build_allele_audit_trail.py` 使用 `random`（非 cohort 主路径；默认不跑）。

## 5. 与 v3 的关系

v3 是历史全库跑次；**v4 GitHub 复现不读取 v3 目录**。开发机可用 `--seed-cache-from /path/to/v3` 加速对齐，**不计入**官方复现路径。
