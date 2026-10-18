#!/usr/bin/env python3
"""
Boilerplate QC + differential expression analysis on nf-core/rnaseq outputs.

Covers:
  1. PCA (log2-CPM of the most variable genes)
  2. Number of detected genes per sample
  3. Number of counts (library size) per sample
  4. DEG analysis with PyDESeq2

Expected inputs
---------------
- Gene count matrix from nf-core/rnaseq, e.g.
    results/star_salmon/salmon.merged.gene_counts.tsv   (default aligner)
    results/star_rsem/rsem.merged.gene_counts.tsv       (--aligner star_rsem)
    results/salmon/salmon.merged.gene_counts.tsv        (--pseudo_aligner salmon)
  Columns: gene_id, gene_name, <sample1>, <sample2>, ...
  (Use the raw *gene_counts.tsv for DESeq2, NOT the TPM or length-scaled files.)

- A metadata CSV you create yourself (the nf-core samplesheet has no conditions):
    sample,condition
    ctrl_1,control
    ctrl_2,control
    trt_1,treated
    trt_2,treated
  Sample names must match the column names in the count matrix.

Install
-------
    pip install pandas numpy matplotlib seaborn scikit-learn "pydeseq2>=0.5"

Usage
-----
    python rnaseq_qc_deg.py \
        --counts results/star_salmon/salmon.merged.gene_counts.tsv \
        --metadata metadata.csv \
        --condition-col condition \
        --reference control --test treated \
        --outdir analysis_out
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference


# ----------------------------------------------------------------------------
# Loading
# ----------------------------------------------------------------------------
def load_counts(path: Path) -> tuple[pd.DataFrame, pd.Series]:
    """Load an nf-core merged gene count table -> (genes x samples int matrix, gene_id->gene_name)."""
    df = pd.read_csv(path, sep="\t")
    id_col = "gene_id" if "gene_id" in df.columns else df.columns[0]
    df = df.set_index(id_col)

    gene_names = df["gene_name"] if "gene_name" in df.columns else pd.Series(df.index, index=df.index)
    # Drop any non-sample annotation columns
    counts = df.drop(columns=[c for c in ("gene_name", "transcript_id(s)") if c in df.columns])

    # Salmon/RSEM give estimated (non-integer) counts; DESeq2 needs integers
    counts = counts.round().astype(int)
    return counts, gene_names


def load_metadata(path: Path, samples: list[str], condition_col: str) -> pd.DataFrame:
    meta = pd.read_csv(path).set_index("sample")
    missing = set(samples) - set(meta.index)
    if missing:
        raise ValueError(f"Samples in count matrix missing from metadata: {sorted(missing)}")
    meta = meta.loc[samples]  # same order as count columns
    meta[condition_col] = meta[condition_col].astype(str)
    return meta


# ----------------------------------------------------------------------------
# QC
# ----------------------------------------------------------------------------
def plot_library_size(counts: pd.DataFrame, meta: pd.DataFrame, cond: str, outdir: Path):
    lib = counts.sum(axis=0).rename("total_counts").to_frame().join(meta[cond])
    lib.to_csv(outdir / "counts_per_sample.csv")

    fig, ax = plt.subplots(figsize=(max(6, 0.4 * len(lib)), 4))
    sns.barplot(data=lib.rename_axis("sample").reset_index(), x="sample",
                y="total_counts", hue=cond, dodge=False, ax=ax)
    ax.set_ylabel("Total assigned counts")
    ax.set_xlabel("")
    ax.set_title("Counts per sample (library size)")
    ax.tick_params(axis="x", rotation=90)
    fig.tight_layout()
    fig.savefig(outdir / "counts_per_sample.png", dpi=150)
    plt.close(fig)


def plot_genes_detected(counts: pd.DataFrame, meta: pd.DataFrame, cond: str, outdir: Path,
                        min_count: int = 1):
    det = (counts >= min_count).sum(axis=0).rename("genes_detected").to_frame().join(meta[cond])
    det.to_csv(outdir / "genes_per_sample.csv")

    fig, ax = plt.subplots(figsize=(max(6, 0.4 * len(det)), 4))
    sns.barplot(data=det.rename_axis("sample").reset_index(), x="sample",
                y="genes_detected", hue=cond, dodge=False, ax=ax)
    ax.set_ylabel(f"Genes with >= {min_count} count(s)")
    ax.set_xlabel("")
    ax.set_title("Detected genes per sample")
    ax.tick_params(axis="x", rotation=90)
    fig.tight_layout()
    fig.savefig(outdir / "genes_per_sample.png", dpi=150)
    plt.close(fig)


def plot_pca(counts: pd.DataFrame, meta: pd.DataFrame, cond: str, outdir: Path, n_top: int = 500):
    # log2-CPM normalisation
    cpm = counts / counts.sum(axis=0) * 1e6
    logcpm = np.log2(cpm + 1)

    # Top-N most variable genes (as DESeq2's plotPCA does)
    top = logcpm.var(axis=1).sort_values(ascending=False).index[:n_top]
    X = logcpm.loc[top].T  # samples x genes
    X = X - X.mean(axis=0)

    pca = PCA(n_components=min(5, X.shape[0]))
    pcs = pca.fit_transform(X)
    var = pca.explained_variance_ratio_ * 100

    pc_df = pd.DataFrame(pcs[:, :2], index=X.index, columns=["PC1", "PC2"]).join(meta[cond])
    pc_df.to_csv(outdir / "pca_coordinates.csv")

    fig, ax = plt.subplots(figsize=(6, 5))
    sns.scatterplot(data=pc_df, x="PC1", y="PC2", hue=cond, s=80, ax=ax)
    for s, row in pc_df.iterrows():
        ax.annotate(s, (row.PC1, row.PC2), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
    ax.set_title(f"PCA (top {n_top} variable genes, log2 CPM)")
    fig.tight_layout()
    fig.savefig(outdir / "pca.png", dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------------
# Differential expression
# ----------------------------------------------------------------------------
def run_deseq2(counts: pd.DataFrame, meta: pd.DataFrame, gene_names: pd.Series, cond: str,
               reference: str, test: str, outdir: Path, min_total: int = 10,
               padj_thr: float = 0.05, lfc_thr: float = 1.0, n_cpus: int = 4) -> pd.DataFrame:
    # Pre-filter lowly expressed genes
    keep = counts.sum(axis=1) >= min_total
    counts_f = counts.loc[keep]
    print(f"[DEG] Keeping {keep.sum()} / {len(keep)} genes with >= {min_total} total counts")

    inference = DefaultInference(n_cpus=n_cpus)
    dds = DeseqDataSet(
        counts=counts_f.T,           # PyDESeq2 expects samples x genes
        metadata=meta,
        design=f"~{cond}",           # add covariates here, e.g. "~batch + condition"
        refit_cooks=True,
        inference=inference,
    )
    dds.deseq2()

    ds = DeseqStats(dds, contrast=[cond, test, reference], inference=inference)
    ds.summary()
    # Optional LFC shrinkage (coefficient name depends on design/levels):
    # ds.lfc_shrink(coeff=f"{cond}[T.{test}]")

    res = ds.results_df.copy()
    res.insert(0, "gene_name", gene_names.reindex(res.index))
    res = res.sort_values("padj")
    res.to_csv(outdir / f"deseq2_{test}_vs_{reference}.csv")

    sig = res[(res.padj < padj_thr) & (res.log2FoldChange.abs() >= lfc_thr)]
    sig.to_csv(outdir / f"deseq2_{test}_vs_{reference}_significant.csv")
    print(f"[DEG] {len(sig)} significant genes (padj < {padj_thr}, |log2FC| >= {lfc_thr}): "
          f"{(sig.log2FoldChange > 0).sum()} up, {(sig.log2FoldChange < 0).sum()} down")

    plot_volcano(res, test, reference, outdir, padj_thr, lfc_thr)
    plot_ma(res, test, reference, outdir, padj_thr)
    return res


def plot_volcano(res, test, reference, outdir, padj_thr, lfc_thr, n_label=15):
    df = res.dropna(subset=["padj"]).copy()
    df["-log10padj"] = -np.log10(df.padj.clip(lower=1e-300))
    df["status"] = "ns"
    df.loc[(df.padj < padj_thr) & (df.log2FoldChange >= lfc_thr), "status"] = "up"
    df.loc[(df.padj < padj_thr) & (df.log2FoldChange <= -lfc_thr), "status"] = "down"

    fig, ax = plt.subplots(figsize=(6, 5))
    sns.scatterplot(data=df, x="log2FoldChange", y="-log10padj", hue="status",
                    palette={"ns": "lightgrey", "up": "firebrick", "down": "steelblue"},
                    s=8, linewidth=0, ax=ax)
    ax.axhline(-np.log10(padj_thr), ls="--", c="k", lw=0.7)
    ax.axvline(lfc_thr, ls="--", c="k", lw=0.7)
    ax.axvline(-lfc_thr, ls="--", c="k", lw=0.7)
    for g, row in df[df.status != "ns"].nsmallest(n_label, "padj").iterrows():
        ax.annotate(row.gene_name if pd.notna(row.gene_name) else g,
                    (row.log2FoldChange, row["-log10padj"]), fontsize=6)
    ax.set_title(f"Volcano: {test} vs {reference}")
    fig.tight_layout()
    fig.savefig(outdir / f"volcano_{test}_vs_{reference}.png", dpi=150)
    plt.close(fig)


def plot_ma(res, test, reference, outdir, padj_thr):
    df = res.dropna(subset=["padj"])
    sig = df.padj < padj_thr
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(np.log10(df.baseMean[~sig] + 1), df.log2FoldChange[~sig], s=4, c="lightgrey")
    ax.scatter(np.log10(df.baseMean[sig] + 1), df.log2FoldChange[sig], s=4, c="firebrick")
    ax.axhline(0, c="k", lw=0.7)
    ax.set_xlabel("log10(mean normalised count + 1)")
    ax.set_ylabel("log2 fold change")
    ax.set_title(f"MA plot: {test} vs {reference}")
    fig.tight_layout()
    fig.savefig(outdir / f"ma_{test}_vs_{reference}.png", dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts", required=True, type=Path, help="nf-core *.merged.gene_counts.tsv")
    p.add_argument("--metadata", required=True, type=Path, help="CSV with 'sample' + condition column")
    p.add_argument("--condition-col", default="condition")
    p.add_argument("--reference", required=True, help="Reference level (e.g. control)")
    p.add_argument("--test", required=True, help="Test level (e.g. treated)")
    p.add_argument("--outdir", default=Path("analysis_out"), type=Path)
    p.add_argument("--n-top-pca", type=int, default=500)
    p.add_argument("--min-count-detected", type=int, default=1,
                   help="Min count for a gene to count as 'detected' in a sample")
    p.add_argument("--padj", type=float, default=0.05)
    p.add_argument("--lfc", type=float, default=1.0)
    p.add_argument("--cpus", type=int, default=4)
    args = p.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid")

    counts, gene_names = load_counts(args.counts)
    meta = load_metadata(args.metadata, list(counts.columns), args.condition_col)
    counts.columns.name = "sample"
    meta.index.name = "sample"
    print(f"Loaded {counts.shape[0]} genes x {counts.shape[1]} samples")

    plot_library_size(counts, meta, args.condition_col, args.outdir)                          # (3)
    plot_genes_detected(counts, meta, args.condition_col, args.outdir, args.min_count_detected)  # (2)
    plot_pca(counts, meta, args.condition_col, args.outdir, args.n_top_pca)                   # (1)
    run_deseq2(counts, meta, gene_names, args.condition_col, args.reference, args.test,       # (4)
               args.outdir, padj_thr=args.padj, lfc_thr=args.lfc, n_cpus=args.cpus)

    print(f"Done. Results in {args.outdir.resolve()}")


if __name__ == "__main__":
    main()