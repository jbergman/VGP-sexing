#!/usr/bin/env python3
"""
Visualize results from the VGP sex-chromosome workflow.

Reads all GCA*/contigs.tsv and GCA*/biosample_metadata.tsv files under
the given output directory and produces a multi-panel figure saved as a PDF.

Usage:
    python visualize_results.py <output_dir> [<figure.pdf>]

Default output: <output_dir>/sex_assignment_summary.pdf
"""

import sys
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import numpy as np


# ── helpers ───────────────────────────────────────────────────────────────────

def load_all(output_dir: Path):
    """Load contigs.tsv and biosample_metadata.tsv from every GCA* folder."""
    contigs_frames, meta_frames = [], []

    for acc_dir in sorted(output_dir.glob("GCA*")):
        acc = acc_dir.name

        ctg = acc_dir / "contigs.tsv"
        if ctg.exists():
            df = pd.read_csv(ctg, sep="\t", low_memory=False)
            df["accession"] = acc
            contigs_frames.append(df)

        meta = acc_dir / "biosample_metadata.tsv"
        if meta.exists():
            df = pd.read_csv(meta, sep="\t", low_memory=False)
            df["accession"] = acc
            meta_frames.append(df)

    contigs = pd.concat(contigs_frames, ignore_index=True) if contigs_frames else pd.DataFrame()
    meta    = pd.concat(meta_frames,   ignore_index=True) if meta_frames   else pd.DataFrame()
    return contigs, meta


def sex_system_per_assembly(contigs: pd.DataFrame) -> pd.Series:
    """Map each accession to 'XY', 'ZW', or 'none annotated'."""
    def classify(type_series):
        types = set(type_series.dropna())
        if "chrX" in types:
            return "XY"
        if "chrZ" in types:
            return "ZW"
        return "none annotated"
    return contigs.groupby("accession")["TYPE"].apply(classify)


def normalise_reported_sex(s: pd.Series) -> pd.Series:
    """Lower-case and strip reported sex; treat actual NaN / blank as missing."""
    return (
        s.where(s.notna(), other="not provided")
         .astype(str)
         .str.strip()
         .replace("", "not provided")
         .str.lower()
    )


MISSING_LABELS = {
    "not collected", "missing", "not provided", "not applicable",
    "na", "nan", "", "none",
}

PALETTE = {
    "XY":             "#4C72B0",
    "ZW":             "#DD8452",
    "none annotated": "#C0C0C0",
    "MALE":           "#4C72B0",
    "FEMALE":         "#DD8452",
    "UNKNOWN":        "#8DA0CB",
}


def _chrom_ids_for_types(contigs: pd.DataFrame, acc: str,
                         types: set[str]) -> set[str]:
    """
    Return the set of chromosome identifiers (genbank + refseq accessions)
    for contigs of the given TYPEs in one assembly.  Using both columns means
    assemblies whose mosdepth files carry NC (RefSeq) names are handled the
    same as those using CM (GenBank) names.
    """
    grp = contigs[(contigs["accession"] == acc) & (contigs["TYPE"].isin(types))]
    ids: set[str] = set()
    for col in ("genbank_accession", "refseq_accession"):
        if col in grp.columns:
            ids |= set(grp[col].dropna())
    return ids


def load_coverage(output_dir: Path, meta: pd.DataFrame,
                  contigs: pd.DataFrame) -> pd.DataFrame:
    """
    For each sample, compute bases-weighted mean coverage for:
      • AUTO contigs         → auto_mean_cov
      • chrX / chrZ contigs  → sex_chr_mean_cov

    Returns a DataFrame with columns
    [sample_id, accession, auto_mean_cov, sex_chr_mean_cov].
    """
    SEX_CHR_TYPES = {"chrX", "chrZ"}
    FALLBACK_EXCL  = {"chrX", "chrY", "chrZ", "chrW", "chrM", "total"}

    # Pre-build lookup dicts per assembly
    auto_sets:    dict[str, set[str]] = {}
    sexchr_sets:  dict[str, set[str]] = {}
    if "TYPE" in contigs.columns:
        for acc in contigs["accession"].unique():
            auto_sets[acc]   = _chrom_ids_for_types(contigs, acc, {"AUTO"})
            sexchr_sets[acc] = _chrom_ids_for_types(contigs, acc, SEX_CHR_TYPES)

    rows = []
    for _, row in meta[["sample_id", "accession"]].drop_duplicates().iterrows():
        sid = row["sample_id"]
        acc = row["accession"]
        path = (output_dir / acc / "callable_sites" / "depths"
                / f"{sid}.mosdepth.summary.txt")
        if not path.exists():
            continue
        try:
            df = pd.read_csv(path, sep="\t")

            def weighted_mean(chrom_set, fallback=False):
                if chrom_set:
                    sub = df[df["chrom"].isin(chrom_set)]
                elif fallback:
                    sub = df[~df["chrom"].isin(FALLBACK_EXCL)]
                else:
                    return float("nan")
                if sub.empty or sub["length"].sum() == 0:
                    return float("nan")
                return sub["bases"].sum() / sub["length"].sum()

            auto_cov   = weighted_mean(auto_sets.get(acc, set()), fallback=True)
            sexchr_cov = weighted_mean(sexchr_sets.get(acc, set()), fallback=False)

            if not (auto_cov != auto_cov):  # skip if auto_cov is NaN
                rows.append({"sample_id": sid, "accession": acc,
                             "auto_mean_cov": auto_cov,
                             "sex_chr_mean_cov": sexchr_cov})
        except Exception:
            pass

    cols = ["sample_id", "accession", "auto_mean_cov", "sex_chr_mean_cov"]
    return pd.DataFrame(rows, columns=cols) if rows else pd.DataFrame(columns=cols)


# ── plots ─────────────────────────────────────────────────────────────────────

def plot_assembly_sex_chr(ax, contigs: pd.DataFrame):
    """Bar: assemblies by sex-chromosome system."""
    sys_map  = sex_system_per_assembly(contigs)
    counts   = sys_map.value_counts().reindex(["XY", "ZW", "none annotated"], fill_value=0)
    colors   = [PALETTE[k] for k in counts.index]

    bars = ax.bar(counts.index, counts.values, color=colors, edgecolor="white", width=0.5)
    for bar, val in zip(bars, counts.values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.4,
                str(val), ha="center", va="bottom", fontsize=9)

    ax.set_title("Reference assemblies by sex-chromosome type", fontsize=10, fontweight="bold")
    ax.set_ylabel("Number of assemblies")
    ax.set_ylim(0, counts.max() * 1.15)
    ax.spines[["top", "right"]].set_visible(False)


def plot_reported_sex(ax, meta: pd.DataFrame):
    """Horizontal bar: all unique reported-sex values."""
    if "sex" not in meta.columns:
        ax.text(0.5, 0.5, "no 'sex' column", ha="center", va="center",
                transform=ax.transAxes)
        return

    norm   = normalise_reported_sex(meta["sex"])
    counts = norm.value_counts()
    counts = counts[counts.index != ""]

    def cat_color(label):
        if label in MISSING_LABELS:
            return "#C0C0C0"
        if label in {"male", "1"}:
            return PALETTE["MALE"]
        if label in {"female", "0"}:
            return PALETTE["FEMALE"]
        return "#8DA0CB"

    colors = [cat_color(lbl) for lbl in counts.index]
    y = range(len(counts))
    ax.barh(list(y), counts.values, color=colors, edgecolor="white")
    ax.set_yticks(list(y))
    ax.set_yticklabels(counts.index, fontsize=8)
    for i, val in enumerate(counts.values):
        ax.text(val + counts.max() * 0.01, i, str(val), va="center", fontsize=7)

    ax.set_title("Reported sex (NCBI BioSample)", fontsize=10, fontweight="bold")
    ax.set_xlabel("Number of samples")
    ax.invert_yaxis()
    ax.spines[["top", "right"]].set_visible(False)

    legend_patches = [
        mpatches.Patch(color=PALETTE["MALE"],   label="male"),
        mpatches.Patch(color=PALETTE["FEMALE"], label="female"),
        mpatches.Patch(color="#8DA0CB",          label="other / ambiguous"),
        mpatches.Patch(color="#C0C0C0",          label="missing / not reported"),
    ]
    ax.legend(handles=legend_patches, fontsize=7, loc="lower right")


def plot_genetic_sex_by_system(ax, meta: pd.DataFrame, contigs: pd.DataFrame):
    """Grouped bar: GENETIC_SEX counts split by XY vs ZW sex system."""
    if "GENETIC_SEX" not in meta.columns:
        ax.text(0.5, 0.5, "no GENETIC_SEX column", ha="center", va="center",
                transform=ax.transAxes)
        return

    sys_map = sex_system_per_assembly(contigs)
    df = meta.copy()
    df["sex_system"] = df["accession"].map(sys_map).fillna("none annotated")

    systems   = ["XY", "ZW"]
    sex_cats  = ["MALE", "FEMALE", "UNKNOWN"]
    # Drop UNKNOWN if zero across both systems
    has_unknown = any(
        len(df[(df["sex_system"] == s) & (df["GENETIC_SEX"] == "UNKNOWN")]) > 0
        for s in systems
    )
    if not has_unknown:
        sex_cats = ["MALE", "FEMALE"]

    x      = np.arange(len(systems))
    width  = 0.8 / len(sex_cats)
    offsets = np.linspace(-(0.8 - width) / 2, (0.8 - width) / 2, len(sex_cats))

    for offset, cat in zip(offsets, sex_cats):
        vals = [
            len(df[(df["sex_system"] == s) & (df["GENETIC_SEX"] == cat)])
            for s in systems
        ]
        bars = ax.bar(x + offset, vals, width=width,
                      color=PALETTE.get(cat, "#C0C0C0"),
                      label=cat, edgecolor="white")
        for bar, val in zip(bars, vals):
            if val > 0:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + ax.get_ylim()[1] * 0.01,
                        str(val), ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(systems)
    ax.set_title("Genetic sex assignment by sex system", fontsize=10, fontweight="bold")
    ax.set_ylabel("Number of samples")
    ax.legend(fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)


def plot_ratio_histogram(ax, meta: pd.DataFrame):
    """Histogram of XA_RATIO and ZA_RATIO values."""
    plotted = False
    for col, lbl, color in [("XA_RATIO", "XY (XA)", PALETTE["XY"]),
                              ("ZA_RATIO", "ZW (ZA)", PALETTE["ZW"])]:
        if col in meta.columns:
            vals = pd.to_numeric(meta[col], errors="coerce").dropna()
            if len(vals):
                ax.hist(vals, bins=40, color=color, alpha=0.7,
                        label=f"{lbl} (n={len(vals)})", edgecolor="white")
                plotted = True

    if plotted:
        ax.axvline(0.75, color="black", linestyle="--", linewidth=1,
                   label="threshold (0.75)")
        ax.set_xlabel("Sex-chr / autosome coverage ratio")
        ax.set_ylabel("Number of samples")
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "no ratio columns found", ha="center", va="center",
                transform=ax.transAxes)

    ax.set_title("Coverage ratio distribution", fontsize=10, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)


def _scatter_cov_vs_ratio(axes, meta: pd.DataFrame, cov_df: pd.DataFrame,
                          contigs: pd.DataFrame, x_col: str, x_label: str):
    """
    Core scatter helper: x_col (coverage column) vs XA/ZA ratio, split by
    XY / ZW.  Discordant samples (reported ≠ genetic sex) are circled.
    """
    sys_map = sex_system_per_assembly(contigs)

    if cov_df.empty or x_col not in cov_df.columns:
        for ax in axes:
            ax.text(0.5, 0.5, "no coverage data", ha="center", va="center",
                    transform=ax.transAxes)
        return

    df = meta.merge(cov_df[["sample_id", "accession", x_col]],
                    on=["sample_id", "accession"], how="inner")
    df["sex_system"] = df["accession"].map(sys_map)
    df["reported_sex_norm"] = (normalise_reported_sex(df["sex"])
                               if "sex" in df.columns else "not provided")

    disc_mask = (
        ((df["reported_sex_norm"] == "male")   & (df.get("GENETIC_SEX") == "FEMALE")) |
        ((df["reported_sex_norm"] == "female") & (df.get("GENETIC_SEX") == "MALE"))
    )

    pairs      = [("XY", "XA_RATIO", axes[0]), ("ZW", "ZA_RATIO", axes[1])]
    sex_colors = {"MALE": PALETTE["MALE"], "FEMALE": PALETTE["FEMALE"],
                  "UNKNOWN": PALETTE["UNKNOWN"]}

    for system, ratio_col, ax in pairs:
        sub = df[df["sex_system"] == system].copy()
        if ratio_col in sub.columns:
            sub[ratio_col] = pd.to_numeric(sub[ratio_col], errors="coerce")
        sub = sub.dropna(subset=[x_col])
        if ratio_col in sub.columns:
            sub = sub.dropna(subset=[ratio_col])

        if sub.empty:
            ax.text(0.5, 0.5, f"no data for {system}", ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_title(f"{system} system: {x_label} vs ratio",
                         fontsize=10, fontweight="bold")
            continue

        gs_col = "GENETIC_SEX" if "GENETIC_SEX" in sub.columns else None
        for sex_val, color in sex_colors.items():
            if gs_col and ratio_col in sub.columns:
                grp = sub[sub[gs_col] == sex_val]
                if len(grp):
                    ax.scatter(grp[x_col], grp[ratio_col],
                               c=color, s=8, alpha=0.5, label=sex_val,
                               linewidths=0)
            elif ratio_col in sub.columns:
                ax.scatter(sub[x_col], sub[ratio_col], c="#888", s=8, alpha=0.5)

        # Hollow circles for discordant samples
        if ratio_col in sub.columns and "GENETIC_SEX" in sub.columns:
            disc         = sub[disc_mask.reindex(sub.index, fill_value=False)]
            full_sys     = df[df["sex_system"] == system]
            n_total_disc = int(disc_mask.reindex(full_sys.index, fill_value=False).sum())
            n_plotted    = len(disc)
            if n_plotted > 0 or n_total_disc > 0:
                suffix = (f", {n_total_disc - n_plotted} no cov. data"
                          if n_total_disc > n_plotted else "")
                ax.scatter(disc[x_col], disc[ratio_col],
                           facecolors="none", edgecolors="black",
                           s=40, linewidths=0.8, alpha=0.9,
                           label=f"sex mismatch (n={n_total_disc}{suffix})",
                           zorder=3)

        if ratio_col in sub.columns:
            ax.axhline(0.75, color="black", linestyle="--", linewidth=0.8,
                       label="threshold (0.75)")

        ax.set_xlabel(x_label)
        ax.set_ylabel(ratio_col.replace("_", " ").replace("RATIO", "ratio"))
        ax.set_title(f"{system} system: {x_label.lower()} vs ratio",
                     fontsize=10, fontweight="bold")
        ax.legend(fontsize=8, markerscale=2)
        ax.spines[["top", "right"]].set_visible(False)


def plot_autosome_cov_vs_ratio(axes, meta, cov_df, contigs):
    _scatter_cov_vs_ratio(axes, meta, cov_df, contigs,
                          x_col="auto_mean_cov",
                          x_label="Autosome mean coverage depth")


def plot_sex_chr_cov_vs_ratio(axes, meta, cov_df, contigs):
    _scatter_cov_vs_ratio(axes, meta, cov_df, contigs,
                          x_col="sex_chr_mean_cov",
                          x_label="Sex chromosome mean coverage depth")


def _draw_concordance_heatmap(ax, df: pd.DataFrame, title: str,
                              row_order: list[str], col_order: list[str]):
    """Draw a single reported-sex × GENETIC_SEX heatmap on ax."""
    if df.empty:
        ax.text(0.5, 0.5, "no data", ha="center", va="center",
                transform=ax.transAxes, fontsize=9)
        ax.set_title(title, fontsize=10, fontweight="bold")
        return

    df = df.copy()
    df["reported_sex_norm"] = normalise_reported_sex(df["sex"])
    ct = pd.crosstab(df["reported_sex_norm"], df["GENETIC_SEX"])

    # Apply shared ordering, filling zeros for categories absent in this panel
    ct = ct.reindex(index=row_order, columns=col_order, fill_value=0)

    im   = ax.imshow(ct.values, aspect="auto", cmap="Blues")
    vmax = max(ct.values.max(), 1)
    ax.set_xticks(range(len(ct.columns)))
    ax.set_xticklabels(ct.columns, fontsize=9)
    ax.set_yticks(range(len(ct.index)))
    ax.set_yticklabels(ct.index, fontsize=8)

    for i in range(len(ct.index)):
        for j in range(len(ct.columns)):
            val = ct.values[i, j]
            if val > 0:
                color = "white" if val > vmax * 0.6 else "black"
                ax.text(j, i, str(val), ha="center", va="center",
                        fontsize=8, color=color)

    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.set_xlabel("Genetic sex (assigned by workflow)")
    ax.set_ylabel("Reported sex (NCBI BioSample)")
    plt.colorbar(im, ax=ax, shrink=0.8, label="sample count")


def plot_concordance(axes, meta: pd.DataFrame, contigs: pd.DataFrame):
    """Two heatmaps side by side: reported sex × GENETIC_SEX for XY and ZW."""
    ax_xy, ax_zw = axes

    if "sex" not in meta.columns or "GENETIC_SEX" not in meta.columns:
        for ax in axes:
            ax.text(0.5, 0.5, "data not available", ha="center", va="center",
                    transform=ax.transAxes)
        return

    sys_map = sex_system_per_assembly(contigs)
    df = meta[meta["GENETIC_SEX"].notna() & (meta["GENETIC_SEX"] != "")].copy()
    df["sex_system"] = df["accession"].map(sys_map)
    df["reported_sex_norm"] = normalise_reported_sex(df["sex"])

    # Build a shared row order: female → male → rest sorted by combined count
    combined_counts = df["reported_sex_norm"].value_counts()
    priority = ["female", "male"]
    rest = [lbl for lbl in combined_counts.index
            if lbl not in priority]  # already sorted by count descending
    row_order = [lbl for lbl in priority if lbl in combined_counts.index] + rest

    # Shared column order for GENETIC_SEX
    all_gs = df["GENETIC_SEX"].dropna().unique().tolist()
    col_priority = ["FEMALE", "MALE", "UNKNOWN"]
    col_order = [c for c in col_priority if c in all_gs] + \
                [c for c in all_gs if c not in col_priority]

    _draw_concordance_heatmap(ax_xy, df[df["sex_system"] == "XY"],
                               "Reported vs genetic sex — XY system",
                               row_order, col_order)
    _draw_concordance_heatmap(ax_zw, df[df["sex_system"] == "ZW"],
                               "Reported vs genetic sex — ZW system",
                               row_order, col_order)


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: visualize_results.py <output_dir> [<figure.pdf>]")

    output_dir  = Path(sys.argv[1])
    figure_path = Path(sys.argv[2]) if len(sys.argv) > 2 \
                  else output_dir / "sex_assignment_summary.pdf"

    print(f"Loading contigs and metadata from {output_dir} …", flush=True)
    contigs, meta = load_all(output_dir)
    n_assemblies  = contigs["accession"].nunique() if not contigs.empty else 0
    n_samples     = len(meta)
    print(f"  {n_assemblies} assemblies, {n_samples} samples loaded", flush=True)

    print("Loading coverage from mosdepth …", flush=True)
    cov_df = load_coverage(output_dir, meta, contigs)
    print(f"  coverage computed for {len(cov_df)} samples", flush=True)

    # ── layout: 5 rows × 2 cols ───────────────────────────────────────────────
    fig = plt.figure(figsize=(14, 25))
    gs  = GridSpec(5, 2, figure=fig, hspace=0.50, wspace=0.35)

    ax_assm  = fig.add_subplot(gs[0, 0])   # assembly sex-chr types
    ax_rep   = fig.add_subplot(gs[0, 1])   # reported sex
    ax_gen   = fig.add_subplot(gs[1, 0])   # genetic sex by system
    ax_hist  = fig.add_subplot(gs[1, 1])   # ratio histogram
    ax_cxy   = fig.add_subplot(gs[2, 0])   # concordance XY
    ax_czw   = fig.add_subplot(gs[2, 1])   # concordance ZW
    ax_scax  = fig.add_subplot(gs[3, 0])   # scatter: autosome cov XY
    ax_scaz  = fig.add_subplot(gs[3, 1])   # scatter: autosome cov ZW
    ax_scsx  = fig.add_subplot(gs[4, 0])   # scatter: sex-chr cov XY
    ax_scsz  = fig.add_subplot(gs[4, 1])   # scatter: sex-chr cov ZW

    plot_assembly_sex_chr(ax_assm, contigs)
    plot_reported_sex(ax_rep, meta)
    plot_genetic_sex_by_system(ax_gen, meta, contigs)
    plot_ratio_histogram(ax_hist, meta)
    plot_concordance([ax_cxy, ax_czw], meta, contigs)
    plot_autosome_cov_vs_ratio([ax_scax, ax_scaz], meta, cov_df, contigs)
    plot_sex_chr_cov_vs_ratio([ax_scsx, ax_scsz], meta, cov_df, contigs)

    fig.suptitle("VGP sex-chromosome workflow — results summary",
                 fontsize=13, fontweight="bold", y=0.99)

    fig.savefig(figure_path, bbox_inches="tight")
    print(f"Figure saved → {figure_path}", flush=True)


if __name__ == "__main__":
    main()
