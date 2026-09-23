# variant_clinical_v3 策略说明书

> 文档版本：**2026-08-31**（以正样本主漏斗为解释主体；金标准 = 新位点人工 pass）  
> 配置版本：`configs/pipeline_gates.json` → `gate_config_version = 2026-08-08.v2-dbsnp`  
> 参考全量门控跑次：`data/pubmed/runs/gate_runs/v2_full_20260808/`  
> **当前正样本金标准**：人工 pass 完整等位 **126**（coding 93 + intron 33；115 PMID / 113 基因）  
> 上游漏斗：highconf 800 → Stage9 genomic 正 → new_site 去掉纯 rs（700）→ 人工审核  
> 负样本是旁路（文献 24 + ClinVar 2026 B/LB），**不混入正样本漏斗**，见 §4.5。

本文以 **正样本** 为主线记录：目标定位、Stage0–9 门控、基因组解析、新位点过滤与人工审核。负样本、解析边角与已知缺口放在对应小节，不作为漏斗主语。

**阅读约定**：打开一份产物时，先问「这条是不是正样本主漏斗上的行」。Stage8 独立负、ClinVar B/LB 补集、Stage9 正负合并表都是旁路或中间件。

---

## 1. 项目定位

### 1.1 目标

从 **2026 年正式发表的 PubMed 文献**中，抽取「变异 ↔ 临床表型」的 **正样本**：相对当前 ClinVar 快照 **未见（ClinVar-absent）**、文献声称新发现、能解析到 GRCh38 `chr / pos / ref / alt` 的等位基因，供下游评测 / 挖掘。

正样本要同时满足四件事：

1. **文献证据**：2026 PubMed 正式文，摘要或全文里变异与临床共现；
2. **ClinVar-absent**：相对版本化 `variant_summary.txt.gz` 用 gene+c./p./rs 对不上已知记录（**≠** 全球文献首报）；
3. **新位点**：交付金标准去掉「仅 rsID、正文没有配对 c.HGVS」的旧位点新表型（§4.3）；
4. **基因组等位**：尽量补全四列，FASTA 复核 REF，GFF3 标 `genomic_region`。

负样本（文献独立库、ClinVar 2026 B/LB）只用于配比与 held-out，**另开轨道**，见 §4.5。

### 1.2 相对 `variant_clinical_2026` 的关系

| 项 | 说明 |
|----|------|
| 形态 | PubMed 轨 **clean-room 拷贝**（`variant_clinical_v3`），不混入历史文献轨交付 |
| 控制方式 | Stage0–9 **门控确认**（APPROVE / HOLD / BYPASS），每阶段写 `gate_report.json` |
| 范围 | 仅 PubMed 正式文献；预印本支路不作为主交付 |
| 正样本主漏斗 | Stage1 检索 → 2 摘要链接 → 3 全文富集 → 4 discovery → 5 ClinVar-absent → 6 highconf → 7 期刊标注 → 8 PDF 致病性 → 9 基因组 → **new_site** → **人工 pass** |
| 核心增强 | ① rsID → **dbSNP RefSNP API**；② Stage9 **基因纠正**；③ highconf 扩至 800 再筛新位点；④ **screen_grade A/B/C** 只排序不删；⑤ **基因组去重** + 证据 **c./p. 配对**；⑥ **基因绑定分流**（rs = 目录基因；HGVS = ±2 句命名窗，§5.8.3）；⑦ Stage9 **GFF3 区域分类**（§5.4.1）；⑧ **new_site** 丢掉纯 rs（§4.3） |

### 1.3 硬性政策声明

- **期刊优先级**（`screen_grade`）、**discovery 文本标签**只是排序 / 注释信号；
- **不得单独**作为变异真实性、致病性或「全球文献首次」的证据；
- **ClinVar absent ≠ 文献首次报道**（仅表示相对该 ClinVar 快照未见）；
- **rsID 出现在 2026 论文 ≠ 新位点**（多为旧位点新表型，金标准默认剔除，除非正文写成 `c.… (rs…)` 配对）；
- **ClinVar B/LB 负样本 ≠ PubMed 正样本**，不得写回 Stage5–6 正样本漏斗。

---

## 2. 目录与入口

```
variant_clinical_v3/
├── configs/
│   ├── pipeline_gates.json
│   ├── pipeline_gate_decisions.example.json
│   ├── journal_priority.json
│   └── queries.json                      # Stage1 检索式（已从 queries_used.json 恢复）
├── scripts/
│   ├── run_controlled_pipeline.pyc
│   ├── run_single_pmid_pipeline.py       # 单篇 PMID 隔离全流程（不改写生产语料）
│   ├── export_genomic_allele_csv.py      # Stage9 基因组导出（当前主交付）
│   ├── extract_negative_candidates.py    # 独立 / 锚定负样本（旁路）
│   ├── extract_clinvar_2026_blb.py       # ClinVar 2026 B/LB 负样本（旁路）
│   ├── filter_review_corpus.py           # review_strict / highconf（正样本）
│   ├── filter_journal_priority.py        # journal_tier → screen_grade
│   ├── filter_new_site_variants.py       # 丢掉纯 rs，保留新位点
│   ├── export_review_pass_csv.py         # 人工 pass → 金标准 CSV
│   ├── check_grch38_coord_refalt.py      # 前四列 vs GRCh38 FASTA
│   ├── validate_grch38_ref_context5.py
│   ├── build_stage_audit_chains_jsonl.py
│   ├── build_interactive_audit_html.py
│   ├── vc_gff3_region.py                 # Stage9：chr:pos ∩ GFF3 CDS/exon/gene
│   ├── test_gff3_region.py
│   ├── vc_text.py                        # 变异/临床 span、基因邻接与 ±2 句命名窗
│   ├── test_a3_neighbor_gene.py          # rs 用目录基因；HGVS 可用邻句命名
│   ├── test_neg_gene_a1a2_b1b2.py
│   ├── vc_paths.py                       # data/pubmed 五顶层路径常量
│   ├── vc_dbsnp.py / …
│   └── …
├── data/pubmed/
│   ├── ref/                              # queries + ClinVar 快照 + Ensembl GRCh38 指针
│   │   ├── queries/
│   │   ├── clinvar/
│   │   └── ensembl/                      # GFF3/FASTA/FAI → genos 共享库符号链接
│   ├── ingest/                           # Stage1：raw / pmids / records / merged
│   ├── pipeline/                         # Stage2–9 产物（按阶段分子目录）
│   │   ├── s2_linking/
│   │   ├── s4_discovery/
│   │   ├── s5_clinvar/{strict,fulltext,enriched}/
│   │   ├── s6_review/
│   │   ├── s8_negatives/{negatives_independent,clinvar_2026_blb}/
│   │   └── s9_genomic/{,new_site,manual_review,audit_trails}/
│   ├── cache/                            # fulltext、review PDF、dbSNP、Ensembl 索引
│   └── runs/                             # gate_runs / single / stats
└── docs/STRATEGY_v3.md
```

路径前缀（下文表格常用）：`ingest/merged/`、`pipeline/s2_linking/`、`pipeline/s4_discovery/`、`pipeline/s5_clinvar/enriched/`（**正样本主轨 CV**）、`pipeline/s6_review/`、`cache/fulltext/`、`cache/pdfs/`、`cache/dbsnp/`、`cache/ensembl/`、`ref/ensembl/`、`pipeline/s9_genomic/`、`s9_genomic/new_site/`、`s9_genomic/manual_review/`、`s9_genomic/audit_trails/`。旁路：`pipeline/s8_negatives/negatives_independent/`、`pipeline/s8_negatives/clinvar_2026_blb/`。  
单篇隔离根：`data/pubmed/runs/single/<PMID>/`（沙箱**内部**仍用历史同名树 `merged/` `linked/` `fulltext/`，不改写生产五顶层）。

产物分四类：**主**（正样本漏斗下游默认输入）/ **旁**（门控分走或负样本轨、可回放）/ **镜**（同一内容拷到另一目录）/ **侧**（TSV、统计、hits、HTML）。  
`pairs*.jsonl` = 一条变异–临床链接；`corpus*.jsonl` = 一篇文献聚合行（**正样本主漏斗走 corpus enriched**）。  
Stage2–5 三条文本轨：**strict**（仅摘要）/ **fulltext**（仅全文新链）/ **enriched**（摘要+全文，**正样本主轨**）。

---

## 3. 门控总策略（Stage0–9）

### 3.1 设计原理

v2 不用「一键跑完再事后挑错」，而是把流水线切成 **Stage0–9**，每个阶段都是一次可审计的质量门：

| 概念 | 含义 |
|------|------|
| **门控（gate）** | 本阶段的输入、条件、动作、输出与人工确认打包成一次决策 |
| **条件状态** | `required`（不过则停/拒）/ `enabled`（生效过滤）/ `disabled`（明确关闭）/ `audit_only`（只记日志不删行） |
| **动作** | `stop`（预检失败即停）/ `filter`（过门保留）/ `route`（分流但保留旁路）/ `annotate`（只加标签不删） |
| **确认决策** | `APPROVE` 进入下一阶段；`HOLD` 暂停；`BYPASS` 记录风险后跳过本门（须书面理由） |
| **失败策略** | `on_gate_error=stop`；`preserve_rejected_rows=true`（拒绝行进 `routes/reject_or_other.jsonl`，不静默丢弃） |

配置源：`configs/pipeline_gates.json`（`gate_config_version = 2026-08-08.v2-dbsnp`）。  
每次门控写入：`data/pubmed/runs/gate_runs/<run_id>/stageN_*/gate_report.json`（含决策人、时间、配置 SHA-256、pass/reject 抽样）。

**Profile**

| Profile | 用途 |
|---------|------|
| `strict`（默认） | 正式交付：required/enabled 条件全部执行 |
| `audit` | 全量评条件但不删行，用于摸底 |
| `custom` | 按单条件声明覆盖 |
| `resume` | 从已 APPROVE 阶段续跑，要求上游与配置 hash 一致 |

**硬性原则（贯穿各 Stage）**

1. 期刊优先级、discovery 文本标签 **只用于排序/注释**，不得单独证明变异真实、致病或「全球首次」。  
2. ClinVar absent **≠** 文献首次报道。  
3. 旁路（reject/unclear/partial）必须可回放；正样本金标准从 highconf → genomic → **new_site** → 人工 pass 取出，不是 highconf 800 本身。

### 3.2 端到端漏斗（正样本主线）

参考全量门控：`data/pubmed/runs/gate_runs/v2_full_20260808/`（各 Stage 均 `APPROVE`）。下表前半为该跑次计数；后半为在此基础上扩 highconf、基因组解析、**新位点过滤与人工审核**后的正样本口径（§4）。独立负 / ClinVar 负不画进这条竖线。

```text
Stage0  preflight            环境/参考就绪
   |
   v
Stage1  pubmed_retrieval     11,795 篇 2026 PubMed 合并文章
   |
   v
Stage2  abstract_linking     4,993 pass / 6,802 reject
   |                         （须同时有变异 + 临床信号）
   v
Stage3  fulltext_enrichment  PMC 富集（无 PMC 则摘要回退，不砍摘要轨）
   |
   v
Stage4  discovery_tagging    65,081 链接行 → novel/expansion 17,800
   |                         （previously_reported / unclear 旁路保留）
   v
Stage5  clinvar_validation   9,883 ClinVar-absent 候选 / 7,917 已知或不可解析
   |
   v
Stage6  review_quality       review_strict 159 → highconf 800（正样本输入）
   |
   v
Stage7  journal_priority     标注 screen_grade A/B/C（C/WATCH 仍保留，不删正样本）
   |
   v
Stage8  pdf_evidence         正样本句级致病性（覆盖 review_strict；其余多为 positive_candidate）
   |
   v
Stage9  genomic_export_qc    chr/pos/ref/alt + GFF3 区域；正样本 JSONL/CSV
   |
   v
new_site                     丢掉「仅 rs、正文无 c.(rs) 配对」→ 正 700
   |
   v
人工审核                     human_verdict=pass 且四列齐全 → **126**（coding 93 / intron 33）
```

**文件从哪一步写出**：见 **§9 ASCII 总览大图**；每步文件释义在 **§3.4–3.13** 与 **§10.1**。金标准 CSV：`s9_genomic/manual_review/review_positive_pass.csv`。

### 3.3 单阶段通用流程

```text
1. 读取本 Stage 的 conditions + checklist
2. 执行 entrypoint 脚本（或预检）
3. 写 routes/pass.jsonl 与 routes/reject_or_other.jsonl（若适用）
4. 写 gate_report.json（计数、拒绝原因、10 条 pass/reject 样例）
5. 人工按 checklist 抽检 → APPROVE / HOLD / BYPASS
6. 仅 APPROVE（或合规 BYPASS）后进入下一 Stage
```

### 3.3.1 单篇 PMID 隔离全流程（检视漏斗 / 改 pipeline）

批量 Stage1 没有「只注入一篇」的插座；`resume` 是全量续跑。要用**一篇文献作为整条 Stage0–9 的初始输入**、且**不改各阶段规则**，走隔离工作区：

```bash
python scripts/run_single_pmid_pipeline.py 41992294
python scripts/run_single_pmid_pipeline.py 41992294 --stop-on-drop
python scripts/run_single_pmid_pipeline.py 41992294 --from stage5_clinvar_validation --force
```

- 工作区：`data/pubmed/runs/single/<PMID>/`（`merged/` `linked/` `fulltext/` `stage_report.json`）
- 目录内文件与生产 **同名同义**（见各 Stage「产物文件」表）
- 额外：`stage_report.json`（每步 status / n_in / n_out / first_drop）、`logs/stage*.log`、`article.json`
- **不要**用单篇目录覆盖生产 `pipeline/s2_linking/`。
- PMID **41992294**：S6 highconf 17 → S9 正样本行；是否进金标准还要过 new_site（纯 rs 丢掉）。独立负是旁路。
- PMID **42068976**（rs 基因回归）：`rs17437411` 绑定 **HOXA7**（dbSNP/PubTator 目录基因）；邻句 `named … (HOTSCRAMBL)` **不得**当 RefSNP 基因。`we_identified` → discovery **high**，可过 Stage6。chr/pos/ref/alt 走 dbSNP。见 §5.8.3。
- **不写**生产 `data/pubmed/ingest/merged/articles.jsonl`、全量 `pipeline/s2_linking/`、`runs/gate_runs/v2_full_*`
- 各 Stage 仍调用现有 CLI，只改 `--input` / `--merged` / `--out-dir` 等路径
- Stage1：生产合并库已有该 PMID 则复用元数据，否则 NCBI EFetch（`pubmed_pipeline.parse_article`）
- `stage_report.json` 的 `status=dropped` + `first_drop` 标明卡在哪一扇门，再改对应脚本后 `--from` 重跑
- Stage3 `--force`：用本地 `fulltext/raw` + `--rebuild-from-parsed --reprocess-links` 重绑基因；**仅在未加 `--no-network` 时**才 `--force-refetch`。否则沙箱/无网会把 `--force-refetch` 打成 fetch_fail，旧 pairs 原样合并回来。
- Stage5 默认只跑 **enriched** 轨（Stage6–9 的输入）；`--clinvar-tracks all` 才与生产 `--corpus all` 一样再扫 strict/fulltext（ClinVar 索引会加载三次）

`--fulltext-scope all` 可在摘要门失败后仍抓 PMC，仅用于诊断，**不是**生产默认（生产 `--scope pass`）。

---

### 3.4 Stage0 — preflight（预检）

**原理**：在烧计算/网络之前，确认可复现基线：查询配置、ClinVar 快照、期刊表、Ensembl GFF3、GRCh38 FASTA/FAI、输出目录可写、禁止无版本覆盖。

**条件要点**：`queries_config_exists`、`clinvar_snapshot_exists`、`ensembl_gff3_exists`、`grch38_fasta_exists` 等为 `required`；记录 ClinVar / PubMed 检索日 / 期刊表 / 词典 / GFF3 / FASTA 版本。

**失败**：`action=stop`，整条流水线不启动。

**实例**：跑次 `v2_full_20260808` Stage0 `decision=APPROVE`，确认输出落在隔离 `gate_runs/` 目录且配置 SHA-256 入库。

**产物文件**（无语料 JSONL）

| 文件 | 类 | 说明 |
|------|----|------|
| `data/pubmed/runs/gate_runs/<run_id>/stage0_preflight/gate_report.json` | 侧 | 生产门控记录 |
| `data/pubmed/runs/single/<PMID>/stage_report.json` | 侧 | 单篇 runner 逐步记 preflight（不写 gate_report） |

只读输入：`configs/pipeline_gates.json`、`configs/queries.json`、`journal_priority.json`、`data/pubmed/ref/clinvar/variant_summary.txt.gz`、`data/pubmed/ref/ensembl/Homo_sapiens.GRCh38.*`。

---

### 3.5 Stage1 — pubmed_retrieval（文献检索）

**原理**：按 `configs/queries.json` 多查询（如 `q1_core_gpc`、`q3_case_reports`、`q6_rsid_trait` 等）拉取 **2026-01-01 起** 正式文献，合并去重 PMID，规范化题录与摘要。

该文件在 2026-08-07 门控预检时存在（9622 bytes），之后与 `pubmed_pipeline.py` 等源码一并从 `configs/` / `scripts/` 消失；检索式完整副本一直在 `data/pubmed/ref/queries/queries_used.json`（Stage1 运行时写出）。**2026-08-14 已从该快照恢复** `configs/queries.json`（7 条 query 原文不变）。NCBI `email`/`tool` 为恢复时补全的 E-utilities 字段。API 密钥走环境变量 `NCBI_API_KEY` 或 `ENTREZ_API_KEY`。

批量全流程（会写生产 `ingest/merged/articles.jsonl`，先确认不要覆盖当前快照）：

```bash
python scripts/pubmed_pipeline.py --config configs/queries.json --mode all \
  --queries q1_core_gpc q3_case_reports q2_high_recall q4_mesh \
            q5_clinical_genetics q6_rsid_trait q7_intergenic_regulatory
python scripts/run_controlled_pipeline.py --profile strict --from stage2_abstract_linking --yes
```

单篇不依赖 query：`python scripts/run_single_pmid_pipeline.py <PMID>`。

**条件要点**：日期窗 `required`；`require_pmid` / `require_title` `required`；Humans MeSH、查询 ID 列表 `enabled`；缺摘要率等 `audit_only`。

**产物文件**

| 文件 | 类 | 说明 |
|------|----|------|
| `data/pubmed/ingest/merged/articles.jsonl` | 主 | 每行一篇：pmid、title、abstract、journal、pub_date、doi、pmcid、source_query_ids |
| `data/pubmed/runs/single/<PMID>/article.json` | 侧 | 同一条的 pretty JSON |
| `data/pubmed/runs/single/<PMID>/raw/efetch_<pmid>.xml` | 侧 | NCBI EFetch 原文（仅走 fetch 时） |

生产 `pubmed_pipeline` 另写 `data/pubmed/ingest/pmids/`、`data/pubmed/ingest/records/` 等中间件；单篇 runner **不写**这些，以免污染全量库。

**实例（门控）**：输入/通过 **11,795** 篇，reject 0。样例 PMID 进入后续链接（如代谢/药物基因组类摘要含 c./rs 与表型）。

---

### 3.6 Stage2 — abstract_linking（摘要变异–临床链接）

**原理**：在题录/摘要中检测**具体变异**与**临床实体**，并用同句/同摘要共现、PubTator 等打上链接 tier。  
**基因绑定**（`vc_text.resolve_gene_for_variant`）：**rsID** 用 PubTator `CorrespondingGene` / dbSNP 目录基因（可加 `located in the GENE locus`）。**不得**用邻句 `we named (SYMBOL)` 或同句 `SYMBOL rs…` 给 rs 绑基因（被扰动的 lncRNA / 实验简写，不是位点所属基因）。**HGVS** 仍可用邻接与 ±2 句命名窗（§5.8.3）。

| Tier | 直觉 |
|------|------|
| **A** | 高置信具体变异 + 紧密临床共现 |
| **B** | 可用但仍弱于 A（距离/方法较弱） |
| **C** | 更松的链接，后续 review 默认不进 strict |

**条件要点**：`variant_and_clinical_dual_presence` `required`；允许 A/B/C；具体变异子集含 c./p./rs/NM_/NC_ 等；排除「仅软变异词」；启用基因解析。

**失败**：进 `stage2_rejected`，不进入主链接语料。

**实例（门控）**：11,795 → **4,993 pass / 6,802 reject**（`variant_clinical_gate_failed`）。  
- Pass 直觉：摘要同时出现 `CYP2D6*10` 类位点与心率/用药结局。  
- Reject 直觉：仅有味觉受体多态性综述、缺少可解析「变异↔疾病」双信号。

**入口**：`scripts/link_variant_phenotype.py`

**产物文件**（`data/pubmed/pipeline/s2_linking/`）

| 文件 | 类 | 说明 |
|------|----|------|
| `articles_pass.jsonl` | 主 | 过双信号门的文章 + `variant_mentions` |
| `articles_reject.jsonl` | 旁 | 未过门文章 |
| `pairs.jsonl` | 主 | 全部链接（含未进 strict 的） |
| `pairs_strict.jsonl` | 主 | 具体变异子集（c./p./rs/NM_…） |
| `corpus.jsonl` | 主 | 文章级语料（全链接） |
| `corpus_strict.jsonl` | 主 | 文章级、具体变异 |
| `pairs_compact.tsv` | 侧 | pairs 瘦表 |
| `summary.json` | 侧 | 计数 |
| `pubtator_parsed.jsonl` | 侧 | PubTator 基因/疾病缓存，Stage3/5/9 会再用 |

---

### 3.7 Stage3 — fulltext_enrichment（全文富集）

**原理**：对 Stage2 pass 集尝试 PMC BioC / EuropePMC JATS / HTML / PDF，在**同句或同段**补变异–临床证据；排除参考文献区；无全文时 **摘要轨继续**（`allow_abstract_fallback=required`）。  
全文每条 link 再走同一套 `resolve_gene_for_variant`（传入 **abstract** + 段落 `context_text`）。rsID 不把后文 `named (HOTSCRAMBL)` 或机制句里的 lncRNA 当成位点基因；HGVS 仍可走邻句命名。

**动作**：`route`（不是简单删文）。门控里 `n_pass=0 / n_reject=4993` 表示「无 PMC、走摘要回退」记在 `reject_or_other`，**注释写明不是删除**。

**实例**：PMID 有摘要链接但无 PMCID → `no_pmcid_abstract_fallback`；同时仍写入 `corpus_enriched.jsonl` 供 Stage4 使用。  
后续独立负样本库（N-ft）也依赖「articles_pass ∩ 已解析全文」池（见 §4.5），**不是**正样本漏斗的下一步。

**入口**：`scripts/fulltext.py enrich --scope pass`

**产物文件**

| 文件 | 类 | 说明 |
|------|----|------|
| `cache/fulltext/raw/*` | 缓存 | BioC / HTML / PDF 原片 |
| `cache/fulltext/parsed/{pmid}.json` | 缓存 | 切段全文 + 已抽链接 |
| `cache/fulltext/linked/articles_fulltext_index.jsonl` | 侧 | 哪些 PMID 有全文、来源 |
| `cache/fulltext/linked/pairs_fulltext.jsonl` 等 | 镜 | 与 `pipeline/s2_linking/` 下同名文件镜像 |
| `cache/fulltext/linked/pairs_fulltext_compact.tsv` | 侧 | 全文 pairs 瘦表 |
| `cache/fulltext/linked/summary.json` | 侧 | 抓取/链接计数 |
| `pipeline/s2_linking/pairs_fulltext.jsonl` | 主 | 全文新 pairs，给 Stage4 |
| `pipeline/s2_linking/corpus_fulltext.jsonl` | 主 | 全文 corpus |
| `pipeline/s2_linking/pairs_enriched.jsonl` | 主 | 摘要 + 全文 pairs |
| `pipeline/s2_linking/corpus_enriched.jsonl` | 主 | **Stage4 默认入口** |
| `pipeline/s2_linking/fulltext_summary.json` | 侧 | 镜像 summary |

`stamp_evidence_source.py` 不新增文件，只给已有 JSONL 打 `evidence_source` / `text_scope`（abstract \| fulltext）。

---

### 3.8 Stage4 — discovery_tagging（发现性文本标注）

**原理**：根据 novel / previously_reported / citation_of_prior / phenotype_expansion 等**措辞线索**打标签；默认 `unclear`；**不删除记录**。`keep_as_new_discovery_labels = novel | phenotype_expansion` 进入「新发现」候选轨。  
`first_in_corpus` 仅 `audit_only`，不能当全球首报。  
`discovery_confidence`：`novel` / `not_reported` / `first_report` / **`we_identified`**（如 “we identified a genetic variant”）/ `newly_identified` 为 **high**。`de_novo_claim` 单独出现仍为 medium。

**实例（门控）**：65,081 链接行 → pass（novel/expansion 轨）**17,800**；旁路 unclear 44,529、previously_reported 2,331、citation_of_prior 421。  
- Pass 例：`IGSF1 c.3467T>A` + intellectual disability，文题含 “Novel … Variant”。  
- 旁路例：`IL13 rs20541` 与过敏性鼻炎，文本像已知关联而非新发现。

**入口**：`scripts/review.py discovery`

**产物文件**（`pipeline/s4_discovery/`）：对每个已有茎  
`{pairs, pairs_strict, pairs_fulltext, pairs_enriched, corpus, corpus_strict, corpus_fulltext, corpus_enriched}` 写一套后缀。

| 后缀 | 类 | 说明 |
|------|----|------|
| `_tagged.jsonl` | 主 | 全量 + `discovery_label` |
| `_novel.jsonl` | 主 | `novel` 或 `phenotype_expansion`；**Stage5 输入** |
| `_previously_reported.jsonl` | 旁 | 文中写已知/已报道 |
| `_unclear.jsonl` | 旁 | 无明确新颖/已知线索 |
| `_novel_first_in_corpus.jsonl` | 侧 | 本语料内该变异字符串最早出现 |
| `_novel_compact.tsv` | 侧 | novel 瘦表 |
| `summary.json` | 侧 | 各茎计数 |

`corpus_fulltext_novel.jsonl` 会再镜一份到 `pipeline/s2_linking/corpus_fulltext_novel.jsonl`（除非 `--no-mirror-novel`）。  
**下游真正要用的是** `pipeline/s4_discovery/corpus_enriched_novel.jsonl`；其余茎是审计/对照。

---

### 3.9 Stage5 — clinvar_validation（ClinVar 对照）

**原理**：用**版本化** `variant_summary.txt.gz` 判断候选是否已在 ClinVar；支持 gene+c. / gene+p. / rsID 精确键；模糊匹配仅审计。  
显式 **`clinvar_absent_proves_literature_first = disabled`**。

**分流直觉**

| 路由 | 含义 |
|------|------|
| variant_absent（pass 主轨） | 快照中未见该解析键 → ClinVar-absent 候选 |
| variant_known + pheno_* | 变异已在 ClinVar（表型已知/新/不清） |
| unresolvable | 键不足以可靠对照 |

**实例（门控）**：17,800 → **9,883 pass**；reject：known+pheno_known 4,279、known+pheno_new 2,238、pheno_unclear 779、unresolvable 621。  
Pass 例：`INSR c.2663A>G` 与遗传性严重胰岛素抵抗等表型，相对该快照为 absent 候选。

**入口**：`scripts/validate_novelty_clinvar.py`

**产物文件**：每条轨一个目录——`pipeline/s5_clinvar/strict/`、`pipeline/s5_clinvar/fulltext/`、`pipeline/s5_clinvar/enriched/`（**主**）。  
文件名把 `corpus_{strict|fulltext|enriched}_novel` 换成带 `_clinvar*` 的茎：

| 文件模式 | 类 | 说明 |
|----------|----|------|
| `{stem}_clinvar.jsonl` | 主 | 每条带 ClinVar 查询结果 |
| `{stem}_clinvar_absent.jsonl` | 主 | 快照中未见；**Stage6 只吃 enriched 这份** |
| `{stem}_clinvar_known_same_pheno.jsonl` | 旁 | ClinVar 已有且表型同类 |
| `{stem}_clinvar_known_new_pheno.jsonl` | 旁 | ClinVar 已有但文称新表型 |
| `{stem}_clinvar_known_pheno_unclear.jsonl` | 旁 | 已知变异、表型关系不清 |
| `{stem}_clinvar_unresolvable.jsonl` | 旁 | 无法用 rs/HGVS 对齐 ClinVar |
| `{stem}_clinvar_absent_compact.tsv` | 侧 | absent 瘦表 |
| `summary.json` | 侧 | 本轨计数 |

---

### 3.10 Stage6 — review_quality（复阅质量子集）

**原理**：从 ClinVar-absent 中筛「可复阅」高置信行：基因已确认、discovery_confidence=high、可信 method（同句/同摘要/全文同句或同段等）、实体临床来源与可信基因来源、排除泛化临床套话；按 `gene|variant` 去重。

| 子集 | 门控口径 |
|------|----------|
| **strict** | link tier **A** only → 门控跑次 **159** 行 |
| **extended** | tier A+B（后续启用） |
| **highconf** | strict ∪ extended → **当前正样本输入 800**（见 §4.1） |

**实例（门控）**：9,387 → **159 pass / 9,228 reject**（`review_quality_gate_failed`）。  
Pass 例：`ACTA1 NM_001100.4:c.1001C` + adult-onset scapuloperoneal myopathy，tier A，fulltext 同句共现，discovery=novel/high。  
PMID **42068976**：基因 **HOXA7**（rs 目录基因，不是邻句 HOTSCRAMBL）、ClinVar-absent、同句实体临床为 tier A；“we identified a genetic variant” 现计为 `we_identified` → **high**，可过 Stage6。

**入口**：`scripts/review.py filter`

**产物文件**（`pipeline/s6_review/`；输入 `pipeline/s5_clinvar/enriched/corpus_enriched_novel_clinvar_absent.jsonl`）

| 文件 | 类 | 说明 |
|------|----|------|
| `…_review_strict.jsonl` | 主 | tier A + 共现正则 + 词表临床 + 可信基因 |
| `…_review_extended.jsonl` | 主 | 放宽到 A 或 B，其余质量门同 |
| `…_review_highconf.jsonl` | 主 | strict ∪ extended，按 locus 去重（当前约 800） |
| `…_review_highconf.tsv` | 侧 | highconf 表 |
| `review_filter_stats.json` | 侧 | 过滤计数与规则摘要 |

---

### 3.11 Stage7 — journal_priority（期刊优先级）

**原理**：用版本化 `configs/journal_priority.json`（ISSN → 别名 → 默认）标注 **A/B/C/WATCH**，导出为 `screen_grade`。  
**`retain_tier_c` / `retain_watch` = required**：C/WATCH **不删**，只影响复阅与全文获取顺序。  
明确禁用：期刊 tier 决定变异有效性或 novelty。

**实例（门控）**：159 → 记为 journal A/B「优先包」36 行，其余 123 为 tier C 旁路（仍保留在 annotated 语料中）。  
当前交付对 highconf 800 打 `screen_grade`（金标准 126 的分布见 §6）。独立负样本另标，见 §4.5。

**入口**：`scripts/filter_journal_priority.py`

对 **strict** 和 **highconf** 各跑一次。输出文件名 = 输入茎 + 后缀：

| 后缀 | 类 | 说明 |
|------|----|------|
| `_journal_annotated.jsonl` | 主 | 全量 + `journal_tier` / `screen_grade`；**Stage9 正样本默认 highconf 这份** |
| `_journal_A.jsonl` | 侧 | 仅期刊档 A |
| `_journal_AB.jsonl` | 侧 | A+B |
| `_journal_other.jsonl` | 侧 | 未进 A/B（C/WATCH 仍保留在 annotated 里） |
| `_journal_unmatched.tsv` | 侧 | 期刊名对不上表 |
| `_journal_filter_stats.json` | 侧 | 计数 |

---

### 3.12 Stage8 — pdf_evidence（PDF 致病性证据，正样本）

**原理（正样本）**：对 review 子集 PDF 抽句级致病性，要求来源可追溯（路径 + SHA-256）、短语锚定到目标变异、区分 ACMG 分类与 in-silico 预测；失败则 `not_evaluated` / 导出时 `positive_candidate`。这是正样本漏斗在 Stage7 之后的证据层，**不是**负样本入口。

**实例（门控）**：159 → **81** 抽到可评标签，**78** `pdf_not_evaluated`。highconf 800 中多数仍为 `positive_candidate`（PDF 只覆盖 strict）。

同文其它位点、独立 N-ft 负样本是旁路（8b/8c），计数与策略见 **§4.5**，不进入正样本金标准。

**产物文件**

8a 下载 + 致病性（输入 review_strict）→ `cache/pdfs/`

| 文件 | 类 | 说明 |
|------|----|------|
| `{pmid}.pdf` | 缓存 | 全文 PDF（或摘要渲染） |
| `pdf_text_cache/PMID{pmid}.txt` | 缓存 | 抽文本，Stage9 字面 g./chr 用 |
| `download_log.jsonl` / `download_summary.json` / `manifest.tsv` | 侧 | 下载过程 |
| `…_review_strict_pdf_pathogenicity.jsonl` | 主 | 句级致病/良性/未评；Stage9 `--pos-patho` |
| 同上 `.tsv` | 侧 | 人读 |

8b 同文负（`--pmid-source jsonl`）→ 前缀 `corpus_strict_other_benign_negatives`

| 后缀 | 类 | 说明 |
|------|----|------|
| `.jsonl` | 主 | 同文其它位点 benign/LB/VUS/对照；locus 去重后 |
| `_pre_locus_dedupe.jsonl` | 侧 | 去重前 |
| `_hits.jsonl` / `_hits.csv` | 侧 | 规则命中明细 |
| `.tsv` / `_summary.json` | 侧 | 人读表与扫描统计 |

8b / 8c 为负样本旁路，产物在 `pipeline/s8_negatives/`，见 §4.5。独立轨 `--exclude-positives` 丢掉与 highconf 同 `gene|allele` 的键，避免把正样本位点写成负。

---

### 3.13 Stage9 — genomic_export_qc（基因组导出与 QC）

**原理（正样本）**：把文献等位解析为 GRCh38 `chr/pos/ref/alt`，按固定 resolver 序尝试，并用 FASTA 复核 REF；不全则 `partial`/`failed`，宁可空坐标也不写错碱基。随后 GFF3 写 `genomic_region`。

**Resolver 序（当前配置）**：literal → ClinVar → PDF 字面 → Ensembl MANE **p.** → Ensembl MANE **c.** → **dbSNP rsID** → c. 仅等位 → unresolved。

**额外现行规则**：证据 c./p. 配对防串扰；基因纠正；完整位点去重；rs↔p 配对；同义 p. 拒 MANE 抬升；裸 `chr:pos` 不得抢先完整 c./rs（§5.8）。

**正样本出口（不要停在正负合并表）**

| 层 | 路径 | 口径 |
|----|------|------|
| 漏斗导出 | `s9_genomic/variants_genomic_positive.csv` | highconf 解析后的正样本中间件 |
| 新位点 | `s9_genomic/new_site/variants_genomic_positive_new_site.jsonl` | 去掉纯 rs |
| **金标准** | `s9_genomic/manual_review/review_positive_pass.csv` | 人工 pass + 四列齐全 **126** |

人工核对走 §10（`review_positive.html`），而不是只看四列数字。负样本 CSV 是旁路，见 §4.5。

**入口**：`scripts/export_genomic_allele_csv.py`；`filter_new_site_variants.py`；`export_review_pass_csv.py`；`validate_grch38_ref_context5.py` / `check_grch38_coord_refalt.py`（可选审计）

**产物文件**（`pipeline/s9_genomic/`；列见 §7，正样本计数见 §4）

| 文件 | 类 | 说明 |
|------|----|------|
| `variants_genomic_positive.csv` / `.jsonl` | 主 | 正样本中间件 |
| `variants_genomic_positive_ok.csv` | 主 | 正 ∩ ok |
| `new_site/variants_genomic_positive_new_site.*` | 主 | 新位点正样本 |
| `manual_review/review_positive_pass.csv` | **金** | 人工 pass 全量 |
| `manual_review/review_positive_pass_coding.csv` | **金** | coding 93 |
| `manual_review/review_positive_pass_noncoding.csv` | **金** | intron 33 |
| `genomic_allele_export_contract.json` | 侧 | 列契约 + 当时合并计数 |
| `variants_genomic_all.csv` / `_ok.csv` / `*_negative.*` | 旁 | 正负合并或负样本，非金标准 |

审核叙事包见 §10.1。

---

### 3.14 如何读一份 `gate_report.json`

必看字段：`stage_id`、`decision`、`approved_by`/`approved_at`、`profile`、`n_input`/`n_pass`/`n_reject`、`reject_reason_counts`、`pass_sample`/`reject_sample`、`gate_config_sha256`、`notes`（若写明 “reject_or_other is a preserved route”，则旁路≠删除）。

确认清单模板见各 Stage 的 `checklist`（配置内英文条目）；中文流程以上述小节为准。

---

## 4. 正样本漏斗与当前交付

主语是 **ClinVar-absent 正样本**。800 / Stage9 CSV / 正负合并表都是中间层；对外口径是人工 pass 的 **126** 个新位点。

### 4.1 Stage6–7：可复阅正样本（highconf）

从 Stage5 `corpus_enriched_novel_clinvar_absent.jsonl` 筛「可复阅」行：基因已确认、discovery_confidence=high、可信 method、实体临床与可信基因来源、排除泛化套话；按 `gene|variant` 去重。期刊只标注 `screen_grade`，不删行。

| 子集 | 行数 | link tier | review_tier | 角色 |
|------|------|-----------|-------------|------|
| **review_strict** | 159 | A | strict | 最严；PDF 致病性主要覆盖这一层 |
| **extended_only** | 641 | B | extended | 放宽到 A 或 B |
| **review_highconf** | **800** | A:159 / B:641 | strict+extended | **Stage9 正样本输入**，不是金标准 |

输入：`pipeline/s6_review/corpus_enriched_novel_clinvar_absent_review_highconf_journal_annotated.jsonl`

### 4.2 Stage9：解析为 GRCh38 等位

`export_genomic_allele_csv.py` 读 highconf，写出 `variants_genomic_positive.*`（resolver / 去重 / GFF3 见 §5）。这是 **正样本中间件**：仍含大量纯 rs 行（旧位点新表型），不能当「新发」交付。

new_site 过滤前正样本输入 **776** 行（`filter_new_site_report.json`）。

### 4.3 new_site：只留新报道位点

脚本：`filter_new_site_variants.py`。政策：

- **丢掉** `variant_text` 仅为 `rs数字`、且正文没有 `c.… (rs…)` / `rs… (c.…)` / `NM_:c. (rs…)` 配对的行；
- **不**因为 ClinVar/dbSNP 给 rs 填了 `c_hgvs` 就留下；
- 非 rs 的 c./p./NM/基因组字面一律保留。

正样本：776 → **kept 700**（丢掉纯 rs 76）。产物：`s9_genomic/new_site/variants_genomic_positive_new_site.{jsonl,csv}`，审核页 `manual_review/review_positive.html`。

### 4.4 人工审核金标准（当前交付）

在 new_site 上点选 `human_verdict=pass`，导出 `export_review_pass_csv.py`。

| 文件 | n | 完整四列 | 说明 |
|------|--:|--------:|------|
| `review_positive_pass.csv` | 129 | **126** | 金标准全集（3 条缺等位，不进四列交付） |
| `review_positive_pass_coding.csv` | 93 | 93 | `region_used=coding` |
| `review_positive_pass_noncoding.csv` | 33 | 33 | 当前均为 intron |

**126** 条：coding 93 / intron 33；link_tier A 79 / B 47；review_tier strict 79 / extended 47；screen_grade B 29 / C 97（无期刊 A）；115 PMID、113 基因；`resolve_status=ok` 全部。label 以文献/PDF 为准：pathogenic 34、likely_pathogenic 15、positive_candidate 47、not_evaluated 28，另有 uncertain/conflicting 各 1——**金标准是「新位点 + 文献临床关联」**，不是 ACMG 金标签。

相对 2.1 周报：比 highconf/rs 混入的旧口径更严，所以条数远小于 800。

### 4.5 负样本旁路（不进入正样本漏斗）

文献渠道抽不到与正样本同口径的「2026 新发负样本」：Stage8 独立库多为旧良性对照 / rs 关联阴性（full_pass 约 1021 locus，软负占多数），new_site 后仅 **24** 条人工 pass（coding 22 / intron 2），相对正样本 126 严重偏少。因此：

| 轨 | 口径 | 交付 | 用途 |
|----|------|------|------|
| 文献独立负 | N-ft 全文池，`--exclude-positives` | `s8_negatives/negatives_independent/`；人工 pass 24 | 文献对照，软负为主 |
| ClinVar A 档 248 | LastEvaluated=2026 + expert panel | `clinvar_2026_blb/variants_genomic_negative_grade_A_ok.*` | 高置信 B/LB；**仅 5 条** First in ClinVar=2026，其余为再评估 |
| ClinVar B 档 419 | VariationID≥4550000 近似 First in ClinVar=2026 + 2 星无冲突 | `…/variants_genomic_negative_first2026_B_ok.*`（列 `first_in_clinvar`） | 与正样本同「2026 才进公共注释」的 held-out 负样本 |

三条负轨都 **不写回** Stage5–6。与文献 24 条坐标：A 档 0 重叠、first2026 B 档 0 重叠。LastEvaluated≠First in ClinVar，A 档 7 条 P/LP→B/LB 只作再分类评估，勿与 419 混成同一 held-out。

正负解耦：文献集合不必同一批 PMID；抽取时按 `gene|allele` 丢掉与正样本同键。

---

## 5. Stage9 基因组解析（正样本中间件）

本节讲 **如何把正样本从文献记号变成四列**。合并表里的负样本行、以及 1821→1594 去重计数，是当时导出快照，**不是** §4.4 金标准。

### 5.1 默认输入 / 输出

**输入（正样本主线）**

| 角色 | 路径 |
|------|------|
| 正样本 | `…/corpus_enriched_novel_clinvar_absent_review_highconf_journal_annotated.jsonl`（800） |
| 致病性（部分覆盖） | `cache/pdfs/…_pdf_pathogenicity.jsonl`（159 行 / 131 PMID） |

Stage9 脚本仍可同时读独立负 JSONL（旁路，§4.5），避免与正样本同键碰撞。

**输出** — `pipeline/s9_genomic/`（正样本看 positive / new_site / manual_review）

| 文件 | 说明 |
|------|------|
| `variants_genomic_positive.csv` / `_ok.csv` | 正样本中间件（历史快照约 655 / ok 530；new_site 源 JSONL 为 776） |
| `new_site/variants_genomic_positive_new_site.*` | 去掉纯 rs 后 **700** |
| `manual_review/review_positive_pass*.csv` | **金标准 126**（§4.4） |
| `variants_genomic_all.csv` / `_negative.*` | 正负合并或负样本旁路 |
| `genomic_allele_export_contract.json` | 列契约 + 当时合并计数 |

> REF±5 QC（`validate_grch38_ref_context5.py`）为可选审计产物，**不作为主交付**；Stage9 导出已含 `reconfirm_ref_from_fasta`。需要时再按需生成，不必常驻目录。

**坐标系**：GRCh38；`position` 1-based；`chromosome` 裸 NCBI；`chromosome_ucsc` 带 `chr` 前缀。

### 5.2 管道快照计数（非金标准）

当时默认一次导出把 800 正 + 1021 独立负写进同一 Stage9 跑次，仅作解析/去重审计：

| 指标 | 值 |
|------|-----|
| scope | `highconf_positive_plus_independent_full_pass_negatives` |
| n_before_genomic_dedup | **1821** |
| n_genomic_dedup_dropped / groups | **227 / 194** |
| n_all | **1594** |
| n_ok | **1355** |
| failed / partial | **180 / 59** |
| 正 / 负（该快照） | **655 / 939** |
| 正 ok / 负 ok | **530 / 825** |
| `allele_query_source=dbsnp` | **465**（多为 rs 行，金标准会丢掉其中纯 rs） |
| `genomic_region`（合并表） | coding 815 / intron 284 / utr 107 / … |

**resolve_method（节选）**：rs_dbsnp 465；ensembl_mane_c_map 380；rs_clinvar 239；ensembl_mane_p_map 180；gene_p_clinvar 110；gene_c_clinvar 73；…

**金标准 126** 的 resolve_method 以 `ensembl_mane_c_map`（78）为主，其次 p-map / gene_p_clinvar；`rs_dbsnp` 仅 1（正文有 c.(rs) 配对才留下）。

### 5.3 Resolver 顺序

1. `literal_text`（**完整** chr:pos/ref/alt 或 g.HGVS；裸 `chr:pos` 仅在无 c./rs 可抬时才作为 `partial` 返回，见 §5.8 B2）  
2. `clinvar`  
3. `pdf_literal_genomic`  
4. `ensembl_mane_p_map`（同义 p. 跳过，见 §5.8.1）  
5. `ensembl_mane_c_map`  
6. **`rs_dbsnp`**（NCBI RefSNP API；`allele_query_source=dbsnp`）  
7. `c_hgvs_alleles_only`  
8. `unresolved`  / 推迟的裸 `chr:pos` `partial`  

成功后 **`reconfirm_ref_from_fasta`**：REF 与 GRCh38 FASTA 不符则清空坐标 → `partial`。

有坐标后立刻做 **GFF3 区域分类**（§5.4.1）：与 Ensembl/GENCODE 115 的 gene/exon/CDS 求交，写入 `genomic_region`（不改 `gene` 目录绑定）。

### 5.4 dbSNP 回退

| 项 | 说明 |
|----|------|
| API | `https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{id}` |
| 缓存 | `data/pubmed/cache/dbsnp/{rs}.json` |
| 多等位 | 优先按频率选最高频 alt；频率打平时再按 **ClinVar/obs 计数**消歧（例：`rs1131692092` A>G vs A>C） |
| 开关 | `--no-dbsnp` / `--dbsnp-cache` |

### 5.4.1 GFF3 区域分类（编码 / 非编码）

**不是**用 rsID 或 dbSNP Gene 栏判定编码区。Stage9 在 `resolve_allele` 得到 GRCh38 `chr:pos` 之后，与 `Homo_sapiens.GRCh38.115.gff3` 求交。

| `genomic_region` | 含义 |
|------------------|------|
| `coding` | 蛋白编码转录本 **CDS** |
| `utr` | **UTR**（GFF3 `five_prime_UTR` / `three_prime_UTR`，或蛋白编码 exon 但不在 CDS） |
| `intron` | 蛋白编码基因体内、不在 exon/UTR/CDS |
| `ncrna_exon` | lncRNA/miRNA 等 **exon** |
| `ncrna_gene` | ncRNA 基因体内、不在 exon |
| `pseudogene` | 假基因体内、不在 exon |
| `intergenic` | 不与任何 gene/exon/CDS/UTR 重叠 |
| （空） | 无解析坐标 |

附加列写入同一份 Stage9 CSV/JSONL：`genomic_region_genes`（重叠符号或 ENSG）、`genomic_region_feature`、`genomic_region_transcripts`、`genomic_region_note`（目录基因与重叠基因不一致时，如 HOXA7 vs `ENSG00000270182`）。5′/3′ 不拆成独立 `genomic_region` 值；需要时可看 `genomic_region_feature`。

例：`rs17437411` @ `7:27158467` → `ncrna_exon` / `ENST00000602610`，不是 HOXA7 CDS。  
实现：`scripts/vc_gff3_region.py`；缓存 `data/pubmed/cache/ensembl/ensembl_gff3_region_index.pkl`；`--no-gff3-region` 可关。回归：`python scripts/test_gff3_region.py`。

`gene` 仍是文献/PubTator/dbSNP 目录绑定；**区域分类是坐标的独立注释**。

### 5.5 基因纠正（Stage9）

摘要链接可能误挂篇级基因（例：`APOL1|rs112720315` → 应为 **TCP10L2**）。`correct_gene_from_text` 与 Stage2/3 的 `resolve_gene_for_variant` 对齐，优先级：

1. **HGVS**：证据句/题名上的文献邻接。**rs 跳过**同句 `SYMBOL rs`（避免 `HOTSCRAMBL rs17437411`）。  
2. 证据上的 locus 短语（`located in the GENE locus`，仅 variant→locus 方向）  
3. **HGVS only**：±2 句命名窗（`named/termed … (SYMBOL)`、`SYMBOL variant`）  
4. 近邻 c.HGVS 基因  
5. PubTator `gene_id` → 符号（可经 dbSNP GeneID 映射）——**rsID 的主目录基因**  
6. dbSNP 唯一 `genes[].locus`（rsID；并提供 chr/pos/ref/alt）

纠正后更新 `gene`、`locus_key_src`，并写入 `resolve_note`。候选若落在停用词（`WT/ROS/HPP/HMNX/PM2/ACMG/…`）则拒绝。

### 5.6 证据 c./p. 串扰问题（已修复）

**现象**：不同 `p_hgvs` 行却共享同一 `c_hgvs` 与前四列（`chr/pos/ref/alt`）。例如旧导出 CSV 行 380–383（ACAT1）与 786–787（UBA1）。

**根因**：同一 PMID 的 `evidence_text` 常并列多个等位基因，例如：

```text
c.439G > T (p.Val147Leu), c.193 A > T (p.Thr65Ser), and c.224 C > A (p.Ala75Asp)
```

旧逻辑对 **p.-only** 行会从 evidence 中取**第一个 orphan c.**（此处为 `c.439G>T`）填入 `c_hgvs_input`，再走 `ensembl_mane_c_map` / ClinVar gene+c.，导致：

| variant_text | 错误继承的 c. | 正确配对 |
|--------------|---------------|----------|
| p.Ala75Asp | c.439G>T | c.224C>A |
| p.Thr65Ser | c.439G>T | c.193A>T |
| p.Val147Leu | c.439G>T | c.439G>T（碰巧正确） |
| p.S56P / p.S621C（UBA1） | c.118-1G>C（并列另一变异） | 无配对 → 应走 p. lift |

**修复**：

1. **配对优先**：仅当 evidence 出现 `c.XXX (p.YYY)` 且 `p.YYY` 与当前 `variant_text` 匹配时，才为 p.-only 行补 `c_hgvs_input`（允许 `c.193 A > T` 空格写法）。
2. **禁止 orphan c.**：p.-only 行不再从多等位句子里抓第一个 c.。
3. **无配对时走 p. lift**：如 UBA1 `p.S56P` → `ensembl_mane_p_map` → `c.166T>C` / chrX:47199096。

### 5.7 基因组位点去重

对 **完整** 等位基因键 `chromosome_ucsc|position|ref|alt` 只保留 1 行（坐标不全的行不去重）。

**保留优先级**（低者优先）：`resolve_status` ok > partial > failed → `positive` > `negative` → `review_tier` strict > extended > independent_neg → `link_tier` A > B → `screen_grade` A > B > C → `nt_pro_only` **nt > both > pro**（同一蛋白突变的 c. 与 p. 行保留核苷酸表述）。

实测：1821 → **1594**（丢弃 227 行 / 194 个位点组）。例如 ACAT1 的 `c.439G>T` 与 `p.Val147Leu` 解析到同一基因组位点后只留 `c.439G>T`。

### 5.8 解析优化（抽检发现，正负共用）

人工审核（含负样本 tier A）时发现的系统错误，已写入抽取/导出层。规则对 **正样本解析同样生效**（c./p. 配对、rs 基因绑定、裸坐标让位）。独立负的 `link_tier` 仍由 Stage9 `assign_negative_link_tier()` 赋值（§6.3）。

#### 5.8.1 rs↔p 串线 + 同义 p. 乱抬升

**现象（PMID 41992294）**：`variant_text=rs1131692092` 却挂上 `p.P1192P`，MANE p-map 造出 `AGG>AGG` / 假 `c.3574_3576delinsCCT`。  
证据句实为 `p.P1192P (rs766447664) and p.S1437S (rs1131692092)` 两个独立等位。

**修复（`export_genomic_allele_csv.py` / `vc_dbsnp.py`）**：

1. rs 行只吸收 **`p.XXX (rsN)` / `c.YYY (rsN)` 配对**，禁止 orphan 首个 p./c.。
2. 同义 p.（`p.P1192P` / `p.Pro1192=`）**跳过** `ensembl_mane_p_map`；`ref==alt` 直接拒绝。
3. `p.S1437S` ↔ ClinVar `p.Ser1437=` 蛋白键互通。
4. rs 在 ClinVar 多等位时不再提前 hard-fail；后续 gene+c/p 限制在该 rs 命中集。
5. dbSNP 多等位：频率打平后按 ClinVar/obs 计数选一（本例 `17:43082450 A>G`）。

全量负表该位点与同文正样本 `p.S1437S` 去重后只留正样本行。

#### 5.8.2 A1+A2+B1+B2（先做、改动小）

针对 11 条人工修正行中的「疾病缩写当基因 / 变异后才写基因 / 裸坐标抢先」：

| 编号 | 层 | 改动 | 修掉的典型错 |
|------|----|------|----------------|
| **A1** | 抽取 `extract_negative_candidates` + `vc_text._GENE_STOP` | 基因停用词：`WT/ROS/HPP/HMNX/PM2/PM3/VUS/ACMG/…` | HPP、PM2、WT、ROS 当基因 |
| **A2** | 抽取 `guess_gene_near` + `_GENE_NEAR` | 支持 `c.xxx in GENE`、变异后基因 | `c.206G>T … in EDA` |
| **B1** | Stage9 `gene_adjacent_to_variant` | 右侧 `in GENE` | 与 A2 同构，导出层纠正 |
| **B2** | Stage9 `resolve_allele` | 裸 `chr:pos`（无 ref/alt）**不得**抢先完整 c./rs；证据仅一个完整 `c.` 时可吸入 rs 行 | `rs3763651` 旁 `9:7103816` 截成 partial |

回归：`scripts/test_neg_gene_a1a2_b1b2.py`。

#### 5.8.3 基因绑定分流（rs 目录基因 vs HGVS 邻句命名）

**错误 A（过度信任 PubTator）**：把 `CorrespondingGene` 当成「论文在讲哪个基因」的全文语义。对 **c./p. HGVS** 且论文当场命名新基因时，目录基因可能落后于正文。

**错误 B（过度信任邻句命名，PMID 42068976）**：把 ±2 句窗里的 `we named … (HOTSCRAMBL)` 绑到 **rs17437411**。该句只说明变异**扰动**一条新命名的反义 lncRNA；语义上**推不出** RefSNP 属于 HOTSCRAMBL。dbSNP / PubTator 对该 rs 的基因是 **HOXA7**。同句也没有「rs17437411 是 HOTSCRAMBL 的变异」。

| 项 | 值 |
|----|-----|
| 变异 | `rs17437411` |
| dbSNP / PubTator | **HOXA7**（GeneID 3204）；chr/pos/ref/alt 一律 RefSNP API |
| 摘要同句 | `… genetic variant, rs17437411, associated with … blood cancers …`（无基因符号） |
| 下一句 | 该变异 disrupt 一条位于 HOXA7 与 HOXA9 之间的 lncRNA，`named … (HOTSCRAMBL)` |

**现行规则**

| 变异类型 | 基因 | 坐标 |
|----------|------|------|
| **rsID** | **PubTator CorrespondingGene** / **dbSNP `genes[].locus`**；证据句 `located in the GENE locus` 可作补充。**跳过**邻句 `named/termed`、`SYMBOL variant`，以及同句 `SYMBOL rs…`（机制简写，如 `HOTSCRAMBL rs17437411`）。 | Stage9：`vc_dbsnp` RefSNP（本环境 NCBI 可能被 SSL/NAC 拦截） |
| **HGVS** | 仍可用 ±2 句 `named … (SYMBOL)` / `SYMBOL variant`，再邻接、PubTator | ClinVar / PDF / Ensembl MANE 等原顺序 |

| 层 | 改动 |
|----|------|
| Stage2/3 | `resolve_gene_for_variant`：rs **跳过**邻接与命名窗，locus 短语后走 PubTator；HGVS 才走邻接 / `named` / neighbor window。 |
| Stage9 | `correct_gene_from_text`：rs 不跑邻句窗；PubTator/dbSNP 可覆盖误绑的 lncRNA 名。 |
| 单篇 runner | `--force` 且 `--no-network` 时 Stage3 用 `--rebuild-from-parsed --reprocess-links`。坐标需要网络才能打 dbSNP。 |

回归测试：`python scripts/test_a3_neighbor_gene.py`（rs → HOXA7；HGVS 抽象替换仍 → HOTSCRAMBL）。  
复跑：`python scripts/run_single_pmid_pipeline.py 42068976 --from stage2_abstract_linking --force`。

**仍未做**：A4 题名/摘要唯一基因回退；C 负样本 PDF 表解析。生产全量 `pipeline/s2_linking/` 未按本分流重跑。

#### 5.8.4 全量重跑口径（本快照）

```bash
python scripts/extract_negative_candidates.py --pmid-source full_pass \
  --exclude-positives …/corpus_enriched_novel_clinvar_absent_review_highconf_journal_annotated.jsonl
python scripts/genomic.py export
python scripts/build_stage_audit_chains_jsonl.py
```

实测：抽取 2557 PMID → **1021** locus（hard 95 / soft 926）；Stage9 1821→**1591**；负 **938**（ok **824**）；坏基因 token 在负导出中为 0。

#### 5.8.5 前四列 FASTA 复核

脚本：`scripts/check_grch38_coord_refalt.py`（比对 GRCh38 FASTA：坐标、REF 跨度、SNV/indel 形态）。

人工 pass 子集：`…/audit_trails/variants_genomic_negative_tier_A_pass.csv` → **28/28 PASS**（2026-08-14；明细 `*_coord_refalt_check.tsv`）。

```bash
python scripts/check_grch38_coord_refalt.py \
  --in-csv …/audit_trails/variants_genomic_negative_tier_A_pass.csv \
  --out-tsv …/audit_trails/variants_genomic_negative_tier_A_pass_coord_refalt_check.tsv
```

正样本金标准同样复核：`manual_review/review_positive_pass_coord_refalt_check.tsv`（及 coding / noncoding 拆表）。

### 5.9 new_site 过滤（正样本金标准前置）

详见 §4.3。Stage9 原表保留不动；过滤产物在 `s9_genomic/new_site/`。正样本 776→700（丢纯 rs 76）。保留的唯一 rs 必须在文献 blob 里与 c.HGVS 括号配对。

---

## 6. 筛选等级列（screen_grade A/B/C）

### 6.1 字段定义

| 列 | 含义 | 正样本金标准 126 | 管道正样本 / 负样本旁路 |
|----|------|------------------|------------------------|
| **`screen_grade`** | 期刊优先级 **A/B/C** | B 29 / C 97（无 A） | highconf 与独立负另计 |
| **`link_tier`** | 证据链接质量 **A/B** | A 79 / B 47 | 管道正：Stage6 `tier` |
| **`review_tier`** | review 档位 | strict 79 / extended 47 | 负旁路：`independent_neg` |
| `label` | 文献/PDF 临床标签 | pathogenic 34、LP 15、positive_candidate 47、not_evaluated 28 | 负：benign / VUS / association_null 等 |

**注意**：`screen_grade` 仅用于复阅排序与分层统计，**不代表**致病性或变异有效性。金标准 126 绝大多数是期刊 C，不能按期刊档砍正样本。

### 6.2 当前分布

**金标准正样本 126**

| | A | B | C |
|--|--:|--:|--:|
| screen_grade | 0 | 29 | 97 |
| link_tier | 79 | 47 | — |
| review_tier | strict 79 | extended 47 | — |
| region | coding 93 | intron 33 | — |

**管道快照（Stage9 正负合并表，非金标准）**

screen_grade：A 13 / B 264 / C 1314；positive (653) A5/B151/C497；negative (938) A8/B113/C817。  
link_tier 正：A 146 / B 507；负：A 159 / B 656 / C 123。  
负 label：association_null 650；uncertain 197；benign 60；likely_benign 28；not_pathogenic 3。

### 6.3 等级来源

- **screen_grade**：`configs/journal_priority.json` → `filter_journal_priority.annotate_row()` → `journal_tier`  
- **link_tier（正样本）**：Stage6 review 输入字段 `tier`（A/B）  
- **link_tier（独立负样本）**：Stage9 导出时 `assign_negative_link_tier()`，与正样本同构的「证据共现强度」：
  - **A**：证据句同时出现本变异 + 良性/VUS/明确无关联等负向描述，且抽取规则为分类/对照类高精度同句规则（如 `classify_paren`、`variant_then_label`、`variant_then_vus`、`control_cue_with_variant`、`reclassified_lb` 等）
  - **B**：同句/同窗较弱证据——`window_rs_null_association`、`nearest_other_to_label`、`window_*` 等；或仅能确认一侧信号
  - **C**：无证据句，或无法确认变异与负向描述共现
- **review_tier**：Stage6 `review_tier`（strict/extended）；独立负固定 `independent_neg`

正样本 highconf 输入已含 journal 注释；负样本在 Stage9 导出时按 PMID 动态匹配期刊表，并写入 `link_tier`。

---

## 7. 完整列契约（CSV / JSONL）

### 7.1 主列（`primary_columns`）

`chromosome_ucsc, position, ref, alt, c_hgvs, nt_pro_only, p_hgvs, c_hgvs_strand, gene, genomic_region, genomic_region_genes, variant_text, chromosome_accession, assembly, pubmed_url`

### 7.2 样本与筛选列

`sample_class, label, screen_grade, link_tier, review_tier, pmid, locus_key_src`

### 7.3 解析溯源列

`resolve_method, allele_query_source, resolve_status, resolve_confidence, resolve_note, chromosome, chromosome_notation, assembly_source, assembly_confirmed, transcript_*, genomic_region_feature, genomic_region_transcripts, genomic_region_note`

---

## 8. 复跑命令

```bash
# rs vs HGVS 基因绑定（不跑全库）
python scripts/test_a3_neighbor_gene.py
python scripts/test_gff3_region.py
python scripts/run_single_pmid_pipeline.py 42068976 --from stage2_abstract_linking --force --no-network

# ── 正样本主线：Stage9 基因组 ──
python scripts/genomic.py export \
  --pos data/pubmed/pipeline/s6_review/corpus_enriched_novel_clinvar_absent_review_highconf_journal_annotated.jsonl

# 新位点过滤 + 重建 review_positive.html
python scripts/filter_new_site_variants.py

# 人工 pass → 金标准 CSV
python scripts/export_review_pass_csv.py

python scripts/check_grch38_coord_refalt.py \
  --in-csv data/pubmed/pipeline/s9_genomic/manual_review/review_positive_pass.csv \
  --out-tsv data/pubmed/pipeline/s9_genomic/manual_review/review_positive_pass_coord_refalt_check.tsv

# 思维链（正样本）
python scripts/build_stage_audit_chains_jsonl.py --class positive

# ── 负样本旁路（勿与正样本漏斗混跑成金标准）──
python scripts/extract_negative_candidates.py \
  --pmid-source full_pass \
  --exclude-positives data/pubmed/pipeline/s6_review/corpus_enriched_novel_clinvar_absent_review_highconf_journal_annotated.jsonl
python scripts/extract_clinvar_2026_blb.py
python scripts/audit_grade_a_first_in_clinvar.py
```

---

## 9. 数据流总览（ASCII）

文件多是因为 **strict / fulltext / enriched 三轨**、**pairs vs corpus**、以及 TSV/旁路/审核侧车，不是每一步都有新规则。  
**主路径只看 enriched → highconf → PDF patho → genomic positive → new_site → 人工 pass。** 负样本在图右侧虚线，见 §4.5。每步释义在 §3.4–3.13。

```text
 configs/queries.json
        |
        v
 +----------------------+
 | S1 pubmed_retrieval  |
 +----------------------+
        |
        |  主  ingest/merged/articles.jsonl
        |  侧  runs/single/<PMID>/article.json
        v
 +---------------------------+
 | S2 abstract_linking       |  基因：邻接/locus；rs→PubTator/dbSNP；HGVS→邻句命名
 +---------------------------+
        |
        |  主  pipeline/s2_linking/articles_pass.jsonl
        |  主  pipeline/s2_linking/pairs.jsonl , pairs_strict.jsonl
        |  主  pipeline/s2_linking/corpus.jsonl , corpus_strict.jsonl
        |  旁  pipeline/s2_linking/articles_reject.jsonl
        |  侧  pairs_compact.tsv , summary.json , pubtator_parsed.jsonl
        v
 +---------------------------+
 | S3 fulltext_enrichment    |   仅 pass ∩ PMC；无 PMC → 摘要回退
 +---------------------------+
        |
        |  缓存  cache/fulltext/raw/* , parsed/{pmid}.json
        |  主    pipeline/s2_linking/pairs_fulltext.jsonl , corpus_fulltext.jsonl
        |  主    pipeline/s2_linking/pairs_enriched.jsonl , corpus_enriched.jsonl  <--+
        |  镜    cache/fulltext/linked/{同上} + articles_fulltext_index.jsonl         |
        |  侧    pairs_fulltext_compact.tsv , fulltext_summary.json                   |
        v                                                                            |
 +---------------------------+                                                       |
 | S4 tag_discovery          |  每个 pairs/corpus 茎写一套后缀                        |
 +---------------------------+                                                       |
        |                                                                            |
        |  主  pipeline/s4_discovery/{stem}_tagged.jsonl                             |
        |  主  pipeline/s4_discovery/{stem}_novel.jsonl  -----> S5                    |
        |  旁  {stem}_previously_reported.jsonl , {stem}_unclear.jsonl               |
        |  侧  {stem}_novel_first_in_corpus.jsonl , *_compact.tsv                    |
        |                                                                            |
        |  **主漏斗入口**  corpus_enriched_novel.jsonl                               |
        v                                                                            |
 +---------------------------+                                                       |
 | S5 clinvar_validation     |  s5_clinvar/{strict,fulltext,enriched}                |
 +---------------------------+                                           |
        |                                                                |
        |  主  {stem}_clinvar.jsonl                                      |
        |  主  {stem}_clinvar_absent.jsonl  -----> S6（只吃 enriched）    |
        |  旁  *_known_same_pheno / *_known_new_pheno / *_pheno_unclear   |
        |  旁  *_clinvar_unresolvable.jsonl                              |
        v                                                                |
 +---------------------------+                                           |
 | S6 review_quality         |                                           |
 +---------------------------+                                           |
        |                                                                |
        |  主  ..._review_strict.jsonl          (tier A, ~159)           |
        |  主  ..._review_extended.jsonl                                 |
        |  主  ..._review_highconf.jsonl        (strict∪extended, ~800)  |
        |  侧  review_filter_stats.json , highconf.tsv                   |
        v                                                                |
 +---------------------------+                                           |
 | S7 journal_priority       |  只标注不删行                              |
 +---------------------------+                                           |
        |                                                                |
        |  主  ..._journal_annotated.jsonl   (+ screen_grade)            |
        |  侧  ..._journal_{A,AB,other}.jsonl , unmatched.tsv            |
        v
 +---------------------------+
 | S8 PDF patho（正样本）     |  覆盖 review_strict；其余 positive_candidate
 +---------------------------+
        |
        |  缓存  cache/pdfs/*.pdf , pdf_text_cache
        |  主    *_pdf_pathogenicity.jsonl
        v
 +---------------------------+
 | S9 genomic_export_qc      |  正样本 chr/pos/ref/alt + GFF3
 +---------------------------+
        |
        |  主  variants_genomic_positive.csv/.jsonl
        v
 +---------------------------+
 | new_site                  |  丢掉纯 rs（正 776→700）
 +---------------------------+
        |
        |  主  new_site/variants_genomic_positive_new_site.*
        |  主  manual_review/review_positive.html
        v
 +---------------------------+
 | 人工审核 pass             |  §10
 +---------------------------+
        |
        |  金  review_positive_pass.csv          **126** 完整等位
        |  金  review_positive_pass_coding.csv    93
        |  金  review_positive_pass_noncoding.csv 33
        |  侧  *_coord_refalt_check.tsv , *_ref_context5.csv

        （旁路，不进金标准）
        S8c 独立负 → negatives_full_pass.jsonl → 人工 pass 24
        ClinVar B/LB → clinvar_2026_blb/  grade A 248 / first2026 B 419
```

**正样本生产计数（当前）**

```text
S1–S5 ClinVar-absent
    --> S6 highconf 800  (strict 159 + extended 641)
    --> S7 screen_grade（只标注）
    --> S8 PDF patho（strict 部分覆盖）
    --> S9 genomic positive（中间件；new_site 源 776）
    --> new_site 700（丢纯 rs 76）
    --> 人工 pass 完整等位 **126**（coding 93 / intron 33；115 PMID）
```

**嫌文件多时的打开顺序（正样本）**

1. 金标准：`pipeline/s9_genomic/manual_review/review_positive_pass.csv`
2. 审核页：`manual_review/review_positive.html`
3. 新位点中间件：`s9_genomic/new_site/variants_genomic_positive_new_site.jsonl`
4. 漏斗输入：`s6_review/…_review_highconf_journal_annotated.jsonl`
5. 单篇：`data/pubmed/runs/single/<PMID>/stage_report.json`
6. 过门原因：`new_site` / `audit_trails` 思维链，或 `runs/gate_runs/.../gate_report.json`

可整目录当缓存忽略：`cache/fulltext/raw|parsed/`、`cache/pdfs/*.pdf` 与 `pdf_text_cache/`、`cache/dbsnp/`、`cache/ensembl/`、`data/pubmed/ref/clinvar/`（输入快照）、`data/pubmed/ref/ensembl/`（基因组指针）、`runs/gate_runs/`（门控审计不是变异行）。

---

## 10. 人工审核：逐步思维链（正样本）

目标：不重跑 resolver，让审核员读 **每条正样本如何从文献走到四列**，判断逻辑是否成立。现行金标准在 `manual_review/review_positive.html` 点选 pass 后导出 §4.4 CSV。

### 10.1 产物

脚本：`scripts/build_stage_audit_chains_jsonl.py`、`build_interactive_audit_html.py`、`filter_new_site_variants.py`、`export_review_pass_csv.py`

| 文件 | 用途 |
|------|------|
| `manual_review/review_positive.html` | **现行审核页**（new_site 正样本） |
| `manual_review/review_positive_pass.csv` | **金标准 126** |
| `new_site/variants_genomic_positive_new_site.jsonl` | 审核候选（700） |
| `new_site/stage_audit_chains*.jsonl` | 正样本思维链 |
| `audit_trails/stage_audit_chains_positive.jsonl` | 管道正样本链（含后来被丢掉的 rs） |
| `audit_annotations.json` | 交互审核标注 |
| `EXAMPLE_positive_stage_chain_SAMD9_c2423AtoG.md` | 格式示例 |

负样本审核页 `review_negative.html` / pass 24 为旁路，见 §4.5。

每条思维链固定步骤：

1. **S1 输入身份**（PMID / gene / variant_text / 正负标签）
2. **S2 HGVS 层级**（nt/pro/both/other）
3. **S3 证据片段**（回挂源 JSONL evidence，防 c./p. 串扰）
4. **S4 基因纠正**
5. **S5 Resolver 路径**
6. **S6 转录本**
7. **S7 基因组等位结果**
8. **S8 机器备注 / QC**
9. **S9 建议核对点**（按 resolve_method 定制问题）

### 10.2 推荐审核节奏（正样本）

| 阶段 | 动作 |
|------|------|
| 范围 | **只审 new_site 正样本**（已去掉纯 rs），打开 `review_positive.html` |
| 分层 | 先 link_tier A / coding，再 intron 与 tier B |
| 填表 | `human_verdict`=`pass`/`doubt`/`fail`，`human_genomic_region`，`human_comment` |
| 导出 | `export_review_pass_csv.py` → `review_positive_pass.csv` |
| 金标准抽查 | 对 pass 再抽对照原文 PDF / ClinVar / 四列 FASTA（`review_positive_pass_coord_refalt_check.tsv`） |

### 10.3 生成命令

```bash
python scripts/filter_new_site_variants.py
python scripts/build_stage_audit_chains_jsonl.py --class positive
python scripts/build_interactive_audit_html.py
python scripts/export_review_pass_csv.py
```

审核时打开 `review_positive.html`：点选 **pass / doubt / fail**，填写理由与问题 Stage。「保存进度到 HTML」把标注写回文件（勿用 VSCode HTML 预览，外链会被改写）。

### 10.4 边界

思维链是**对已导出字段 + 源证据的可解释重述**，不是重新推断坐标；若导出本身错了，链会忠实地“讲错”。人工要核对的是：**叙事步骤是否自洽、证据是否支持最终四列**。

---

## 11. 前几列不全行（策略说明）


核心列 `chromosome_ucsc / position / ref / alt` 不全时记 `partial` 或 `failed`，宁可留空也不写入 FASTA 不一致坐标。

常见原因：PDF 坐标残缺；Ensembl c. 抬失败（coding REF 不符 / 超 CDS）；ClinVar 多等位歧义；仅 p. 无 c.；RefSNP API 404（历史 rs）。

去重后 failed / partial 见当时 Stage9 合并快照（§5.2）；金标准 126 条 `resolve_status=ok`。更宽漏斗里的 rs 行多数不会进入 new_site。

---

## 12. 已知缺口

1. **金标准正样本仍少（126）**：highconf 800 含大量纯 rs；new_site 丢掉后人工 pass 才 126。扩容应继续挖 2026 文献里的 **c./p. 新位点**，而不是把 rs 放回金标准。  
2. **PDF 致病性**仅覆盖 strict 159；126 里仍有大量 `positive_candidate` / `not_evaluated`。  
3. **生产 s2_linking 未按 rs/HGVS 分流全量回放**：邻句窗绑错 lncRNA 的 rs 需 Stage2–9 重跑。金标准已基本不含纯 rs，影响小于管道中间件。  
4. 文献负样本 **不够**（人工 pass 仅 24，软负为主）；2026 held-out 负样本改走 ClinVar First-in-ClinVar B 档 419（§4.5）。A 档 248 多数不是 2026 首报。  
5. 部分历史 rs RefSNP 404；部分脚本仍为 `.pyc`。  
6. `screen_grade` / discovery 标签 **不得**当作金标准致病性（126 里期刊 C 占 97）。  
7. **A4/C 未做**：题名/摘要唯一基因回退、独立负 PDF 表解析。  
8. **`we_identified` 已升为 high**（§3.8）：全量生产语料未按此重标；重跑 Stage4–6 可能增加 highconf，但仍须过 new_site + 人工。  
9. Stage9 须带 `genomic_region*`；无坐标的行这些字段为空。

---

## 13. 关键路径速查

| 用途 | 路径 |
|------|------|
| 策略本文 | `docs/STRATEGY_v3.md` |
| **金标准正样本** | `pipeline/s9_genomic/manual_review/review_positive_pass.csv`（coding / noncoding 拆表同目录） |
| 正样本审核页 | `pipeline/s9_genomic/manual_review/review_positive.html` |
| 新位点 JSONL | `pipeline/s9_genomic/new_site/variants_genomic_positive_new_site.jsonl` |
| 正样本 highconf（漏斗输入） | `pipeline/s6_review/corpus_enriched_novel_clinvar_absent_review_highconf_journal_annotated.jsonl` |
| 正样本 strict | `pipeline/s6_review/corpus_enriched_novel_clinvar_absent_review_strict.jsonl` |
| 门控配置 / 期刊表 | `configs/pipeline_gates.json` / `journal_priority.json` |
| 文献负（旁路） | `pipeline/s8_negatives/negatives_independent/negatives_full_pass*.jsonl`；pass 24：`manual_review/review_negative_pass.csv` |
| ClinVar 负（旁路） | `pipeline/s8_negatives/clinvar_2026_blb/`（A 248 / first2026 B 419） |
| GFF3 区域回归 | `scripts/test_gff3_region.py` |
| 单篇 PMID 全流程 | `scripts/run_single_pmid_pipeline.py` → `data/pubmed/runs/single/<PMID>/stage_report.json` |

---

*文档版本 2026-08-31。以正样本主漏斗为解释主体。修改门控、new_site 政策或人工 pass 后，请同步 §3.2、§4、§6.2、§9 计数；负样本旁路变更同步 §4.5。基因绑定规则变更时同步 §3.6–3.7、§5.5、§5.8.3。*
