# PDF 下载逻辑（Stage8）

脚本：`scripts/extras.py download-pdfs`（仅作可选 fallback）  
输出：`data/pubmed/cache/pdfs/PMID{pmid}.pdf` 与 `pdf_text_cache/PMID{pmid}.txt`（后者由致病性脚本生成）。

Stage8 首先由 `materialize_review_text_cache.py` 从 Stage3 OA/fulltext 生成
highconf 文本并标注。只有设置 `VARIANT_CLINICAL_STAGE8_PDF=1` 时，下面
流程才用于仍缺文本的 PMID。

## PDF fallback 顺序

1. **本地已有**：若 `PMID{pmid}.pdf` 存在且为有效 PDF（`%PDF` 且 >1KB），大于 ~80KB 则视为已下载并跳过。
2. **Unpaywall OA**：有 DOI 时请求 `api.unpaywall.org`，依次尝试 `url_for_pdf` / `url`。
3. **PMC OA PDF**（v4 新增）：有 PMCID 时尝试  
   `https://www.ncbi.nlm.nih.gov/pmc/articles/{PMCID}/pdf/` 与 `pmc.ncbi.nlm.nih.gov` 镜像。
4. **Europe PMC**：检索 `fullTextUrlList` 中 PDF 链接，或 `europepmc.org/articles/{PMCID}?pdf=render`。
5. **回退（非 `--publisher-only`）**：本地 BioC 全文（`cache/fulltext/raw/{pmcid}.bioc.json`）或摘要渲染为 PDF（fpdf）。

Gold cohort 模式下 Stage8 对 **highconf** 全部 PMID 执行上述流程（不仅 strict 159 行），以便 extended 层金标准也有 PDF/OA 证据。

## 依赖

- 网络：Unpaywall / Europe PMC / PMC PDF  
- OA 文章：通常由步骤 2–4 命中；仅有 PMCID 的 OA 应优先命中步骤 3。  
- 无 OA 时仍可用 BioC/摘要回退供 `annotate_pathogenicity_from_pdfs.py` 抽句（非出版社 PDF）。
