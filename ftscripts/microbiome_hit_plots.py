#!/usr/bin/env python3
"""
Post-run visualizations for FastTarget MGnify microbiome off-target hits.
"""

import argparse
import html
import math
import os
import re
import sys

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import pyarrow.parquet as pq


HIT_COLUMNS = [
    "gene",
    "representative_genome_id",
    "subject_protein_id",
    "pident",
    "qcovhsp",
    "evalue",
    "bitscore",
]
TAXONOMY_COLUMNS = [
    "representative_genome_id",
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]
TAXONOMY_RANKS = ["domain", "phylum", "class", "order", "family", "genus"]
UNCLASSIFIED = "Unclassified"
OTHER = "Other"


def fail(message):
    """
    Stop script execution with a formatted error message.

    :param message: Error message to report.
    """
    raise SystemExit(f"ERROR: {message}")


def warn(warnings, message):
    """
    Store and print a warning message.

    :param warnings: Mutable list of warning messages.
    :param message: Warning message to append and print.
    """
    warnings.append(message)
    print(f"WARNING: {message}", file=sys.stderr)


def output_default(args):
    """
    Build the default output directory for microbiome hit plots.

    :param args: Parsed command line arguments.
    :return: Default plot output directory.
    """
    return os.path.join(
        args.output_path,
        args.organism_name,
        "offtarget",
        "microbiomes",
        args.catalogue,
        "plots",
    )


def hits_path(args):
    """
    Build the expected consolidated microbiome hits Parquet path.

    :param args: Parsed command line arguments.
    :return: Path to the consolidated hits Parquet file.
    """
    return os.path.join(
        args.output_path,
        args.organism_name,
        "offtarget",
        "microbiomes",
        args.catalogue,
        "species_blast_results",
        f"{args.catalogue}_offtarget_hits.parquet",
    )


def taxonomy_path(args):
    """
    Build the expected representative taxonomy Parquet path.

    :param args: Parsed command line arguments.
    :return: Path to the representative taxonomy Parquet file.
    """
    return os.path.join(
        args.databases_path,
        "microbiomes",
        args.catalogue,
        "species_catalogue",
        "representative_taxonomy.parquet",
    )


def safe_filename(value):
    """
    Convert a user-provided value into a filesystem-safe filename component.

    :param value: Value to sanitize.
    :return: Sanitized filename component.
    """
    safe_value = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("._")
    return safe_value or "locus_tag"


def validate_parquet_columns(file_path, expected_columns, label):
    """
    Validate that a Parquet file contains the expected columns.

    :param file_path: Path to the Parquet file.
    :param expected_columns: Columns required by the script.
    :param label: Human-readable label used in error messages.
    """
    schema_names = set(pq.ParquetFile(file_path).schema_arrow.names)
    missing = [column for column in expected_columns if column not in schema_names]
    if missing:
        fail(f"{label} is missing expected columns {missing}: {file_path}")


def load_hits_for_loci(hits_parquet, locus_tags):
    """
    Load consolidated microbiome hits for the requested locus tags.

    :param hits_parquet: Path to the consolidated hits Parquet file.
    :param locus_tags: Locus tags to load.
    :return: DataFrame with one hit row per gene and representative genome.
    """
    validate_parquet_columns(hits_parquet, HIT_COLUMNS, "Hits parquet")
    dataset = ds.dataset(hits_parquet, format="parquet")
    table = dataset.to_table(
        columns=HIT_COLUMNS,
        filter=ds.field("gene").isin(list(locus_tags)),
    )
    hits = table.to_pandas()
    if hits.empty:
        return pd.DataFrame(columns=HIT_COLUMNS)
    hits = hits.drop_duplicates(["gene", "representative_genome_id"], keep="first")
    for column in ["pident", "qcovhsp", "evalue", "bitscore"]:
        hits[column] = pd.to_numeric(hits[column], errors="coerce")
    return hits


def load_taxonomy(taxonomy_parquet):
    """
    Load representative taxonomy from Parquet.

    :param taxonomy_parquet: Path to the representative taxonomy Parquet file.
    :return: Normalized taxonomy DataFrame.
    """
    validate_parquet_columns(taxonomy_parquet, TAXONOMY_COLUMNS, "Taxonomy parquet")
    table = ds.dataset(taxonomy_parquet, format="parquet").to_table(
        columns=TAXONOMY_COLUMNS
    )
    taxonomy = table.to_pandas()
    return normalize_taxonomy(taxonomy)


def normalize_taxonomy(taxonomy):
    """
    Normalize missing and unclassified taxonomy values.

    :param taxonomy: Raw taxonomy DataFrame.
    :return: Taxonomy DataFrame with consistent unclassified labels.
    """
    taxonomy = taxonomy.copy()
    for column in TAXONOMY_COLUMNS:
        taxonomy[column] = taxonomy[column].fillna("unclassified").astype(str)
    for column in TAXONOMY_RANKS + ["species"]:
        values = taxonomy[column].str.strip()
        taxonomy[column] = values.mask(
            values.eq("") | values.str.lower().eq("unclassified"),
            UNCLASSIFIED,
        )
    return taxonomy


def build_species_hits_table(hits, taxonomy, locus_tags, warnings):
    """
    Join hit rows with taxonomy and collect validation warnings.

    :param hits: Filtered microbiome hits DataFrame.
    :param taxonomy: Representative taxonomy DataFrame.
    :param locus_tags: Requested locus tags.
    :param warnings: Mutable list of warning messages.
    :return: Hit table enriched with taxonomy columns.
    """
    missing_loci = [locus for locus in locus_tags if locus not in set(hits["gene"])]
    for locus in missing_loci:
        warn(warnings, f"Locus tag {locus} has no microbiome hits.")

    taxonomy_ids = set(taxonomy["representative_genome_id"])
    hit_ids = set(hits["representative_genome_id"]) if not hits.empty else set()
    missing_taxonomy_ids = sorted(hit_ids - taxonomy_ids)
    if missing_taxonomy_ids:
        warn(
            warnings,
            f"{len(missing_taxonomy_ids)} representative_genome_id values from hits "
            "were not found in representative_taxonomy.parquet.",
        )

    species_hits = hits.merge(taxonomy, on="representative_genome_id", how="left")
    for column in TAXONOMY_RANKS + ["species"]:
        species_hits[column] = species_hits[column].fillna(UNCLASSIFIED)
    return species_hits


def summarize_by_gene(species_hits, locus_tags):
    """
    Summarize microbiome hit counts and alignment metrics by gene.

    :param species_hits: Hit table enriched with taxonomy columns.
    :param locus_tags: Requested locus tags.
    :return: DataFrame with one summary row per requested locus tag.
    """
    rows = []
    for gene in locus_tags:
        data = species_hits[species_hits["gene"] == gene]
        rows.append(
            {
                "gene": gene,
                "total_representative_hits": data["representative_genome_id"].nunique(),
                "families_with_hits": count_classified(data, "family"),
                "genera_with_hits": count_classified(data, "genus"),
                "mean_pident": data["pident"].mean(),
                "sd_pident": data["pident"].std(),
                "mean_qcovhsp": data["qcovhsp"].mean(),
                "sd_qcovhsp": data["qcovhsp"].std(),
            }
        )
    return pd.DataFrame(rows)


def count_classified(data, rank):
    """
    Count classified taxa for a taxonomy rank.

    :param data: Hit table subset.
    :param rank: Taxonomy rank column to count.
    :return: Number of non-unclassified taxa.
    """
    if data.empty:
        return 0
    return data.loc[data[rank] != UNCLASSIFIED, rank].nunique()


def summarize_rank(species_hits, taxonomy, locus_tags, rank):
    """
    Summarize hits and penetrance for a taxonomy rank.

    :param species_hits: Hit table enriched with taxonomy columns.
    :param taxonomy: Representative taxonomy DataFrame.
    :param locus_tags: Requested locus tags.
    :param rank: Taxonomy rank to summarize.
    :return: Summary DataFrame for the selected taxonomy rank.
    """
    # Use the displayed taxon name for readable plots and TSVs. The main pipeline
    # uses hierarchical paths when strict taxonomy group identity is required.
    total_column = f"total_{rank}_representatives"
    percent_column = f"percent_{rank}_representatives_hit"
    columns = [
        "gene",
        rank,
        "hit_representatives",
        "percent_of_gene_hits",
        total_column,
        percent_column,
        "mean_pident",
        "sd_pident",
        "mean_qcovhsp",
        "sd_qcovhsp",
    ]
    totals = taxonomy.groupby(rank)["representative_genome_id"].nunique()
    rows = []
    for gene in locus_tags:
        gene_hits = species_hits[species_hits["gene"] == gene]
        gene_total = gene_hits["representative_genome_id"].nunique()
        if gene_hits.empty:
            continue
        grouped = gene_hits.groupby(rank, dropna=False)
        for taxon, group in grouped:
            hit_representatives = group["representative_genome_id"].nunique()
            total_representatives = int(totals.get(taxon, 0))
            rows.append(
                {
                    "gene": gene,
                    rank: taxon,
                    "hit_representatives": hit_representatives,
                    "percent_of_gene_hits": percentage(hit_representatives, gene_total),
                    total_column: total_representatives,
                    percent_column: percentage(hit_representatives, total_representatives),
                    "mean_pident": group["pident"].mean(),
                    "sd_pident": group["pident"].std(),
                    "mean_qcovhsp": group["qcovhsp"].mean(),
                    "sd_qcovhsp": group["qcovhsp"].std(),
                }
            )
    return pd.DataFrame(rows, columns=columns)


def percentage(numerator, denominator):
    """
    Calculate a percentage while handling empty denominators.

    :param numerator: Numerator value.
    :param denominator: Denominator value.
    :return: Percentage value.
    """
    if denominator in (0, None) or pd.isna(denominator):
        return 0.0
    return float(numerator) / float(denominator) * 100.0


def prepare_matplotlib(outdir):
    """
    Configure matplotlib for headless plot generation.

    :param outdir: Output directory used for matplotlib cache if needed.
    :return: matplotlib.pyplot module.
    """
    os.environ.setdefault("MPLCONFIGDIR", os.path.join(outdir, ".matplotlib"))
    os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def save_figure(plt, fig, outdir, basename):
    """
    Save a matplotlib figure as PNG and SVG.

    :param plt: matplotlib.pyplot module.
    :param fig: matplotlib figure object.
    :param outdir: Output directory.
    :param basename: Output file basename without extension.
    """
    for extension in ["png", "svg"]:
        fig.savefig(
            os.path.join(outdir, f"{basename}.{extension}"),
            dpi=220,
            bbox_inches="tight",
        )
    plt.close(fig)


def plot_total_hits(summary_by_gene, outdir):
    """
    Plot total representative hits by gene.

    :param summary_by_gene: Gene-level summary DataFrame.
    :param outdir: Output directory.
    """
    plt = prepare_matplotlib(outdir)
    fig_width = max(6, len(summary_by_gene) * 0.8)
    fig, ax = plt.subplots(figsize=(fig_width, 4.5))
    bars = ax.bar(
        summary_by_gene["gene"],
        summary_by_gene["total_representative_hits"],
        color="#4c78a8",
    )
    ax.set_xlabel("locus_tag")
    ax.set_ylabel("Representatives with hit")
    ax.set_title("Total microbiome representative hits by gene")
    ax.tick_params(axis="x", rotation=45)
    for bar in bars:
        height = bar.get_height()
        ax.annotate(
            f"{int(height)}",
            (bar.get_x() + bar.get_width() / 2, height),
            ha="center",
            va="bottom",
            fontsize=9,
        )
    ax.margins(y=0.15)
    save_figure(plt, fig, outdir, "total_species_hits_barplot")


def plot_composition_stacked_bar(summary_rank_df, rank, locus_tags, threshold, outdir):
    """
    Plot a 100 percent stacked composition barplot for a taxonomy rank.

    :param summary_rank_df: Rank-level summary DataFrame.
    :param rank: Taxonomy rank column to plot.
    :param locus_tags: Requested locus tags.
    :param threshold: Minimum percentage for taxa to remain in the legend.
    :param outdir: Output directory.
    """
    if summary_rank_df.empty:
        return
    plt = prepare_matplotlib(outdir)
    data = summary_rank_df.copy()
    keep_taxa = set(
        data.loc[data["percent_of_gene_hits"] >= threshold, rank].astype(str)
    )
    data["plot_taxon"] = data[rank].where(
        data[rank].isin(keep_taxa) | data[rank].eq(UNCLASSIFIED),
        OTHER,
    )
    plot_data = (
        data.groupby(["gene", "plot_taxon"], as_index=False)["hit_representatives"]
        .sum()
        .rename(columns={"plot_taxon": rank})
    )
    totals = plot_data.groupby("gene")["hit_representatives"].transform("sum")
    plot_data["percent"] = np.where(
        totals > 0,
        plot_data["hit_representatives"] / totals * 100.0,
        0.0,
    )
    pivot = (
        plot_data.pivot(index="gene", columns=rank, values="percent")
        .reindex(locus_tags)
        .fillna(0.0)
    )
    ordered_columns = order_composition_columns(pivot.columns)

    fig_width = max(7, len(locus_tags) * 0.8)
    fig, ax = plt.subplots(figsize=(fig_width, 5))
    bottom = np.zeros(len(pivot))
    colors = color_map_for_categories(ordered_columns)
    x = np.arange(len(pivot))
    for taxon in ordered_columns:
        values = pivot[taxon].to_numpy()
        ax.bar(x, values, bottom=bottom, label=taxon, color=colors[taxon], width=0.7)
        bottom += values
    ax.set_xticks(x)
    ax.set_xticklabels(pivot.index, rotation=45, ha="right")
    ax.set_ylim(0, 100)
    ax.set_xlabel("locus_tag")
    ax.set_ylabel("% of representatives with hit")
    ax.set_title(f"{rank.capitalize()} hit composition")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    save_figure(plt, fig, outdir, f"{rank}_hit_composition_stacked_barplot")


def order_composition_columns(columns):
    """
    Order composition categories with special groups at the end.

    :param columns: Taxon names present in the composition table.
    :return: Ordered list of taxon names.
    """
    columns = list(columns)
    tail = [taxon for taxon in [UNCLASSIFIED, OTHER] if taxon in columns]
    head = sorted([taxon for taxon in columns if taxon not in tail])
    return head + tail


def color_map_for_categories(categories):
    """
    Build a stable color map for composition categories.

    :param categories: Category names to color.
    :return: Dictionary mapping category names to colors.
    """
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt

    base = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("Set3").colors)
    colors = {}
    for index, category in enumerate(categories):
        if category == OTHER:
            colors[category] = "#b8b8b8"
        elif category == UNCLASSIFIED:
            colors[category] = "#5f6368"
        else:
            colors[category] = mcolors.to_hex(base[index % len(base)])
    return colors


def plot_family_penetrance_heatmap(summary_family, locus_tags, min_hits, min_coverage, outdir):
    """
    Plot family-level hit penetrance across requested genes.

    :param summary_family: Family-level summary DataFrame.
    :param locus_tags: Requested locus tags.
    :param min_hits: Minimum hit representatives for family inclusion.
    :param min_coverage: Minimum family coverage percentage for inclusion.
    :param outdir: Output directory.
    """
    if summary_family.empty:
        return
    plt = prepare_matplotlib(outdir)
    selected = summary_family[
        (summary_family["hit_representatives"] >= min_hits)
        | (summary_family["percent_family_representatives_hit"] >= min_coverage)
    ]
    if selected.empty:
        selected = summary_family.nlargest(
            min(25, len(summary_family)),
            "hit_representatives",
        )
    families = (
        selected.groupby("family")
        .agg(
            max_penetrance=("percent_family_representatives_hit", "max"),
            total_hits=("hit_representatives", "sum"),
        )
        .sort_values(["max_penetrance", "total_hits"], ascending=False)
        .index
    )
    matrix = (
        selected.pivot_table(
            index="family",
            columns="gene",
            values="percent_family_representatives_hit",
            aggfunc="max",
            fill_value=0.0,
        )
        .reindex(index=families, columns=locus_tags)
        .fillna(0.0)
    )
    fig_height = min(18, max(4, 0.35 * len(matrix) + 1.5))
    fig_width = max(7, 0.75 * len(locus_tags) + 3)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    im = ax.imshow(matrix.to_numpy(), aspect="auto", cmap="Reds", vmin=0)
    ax.set_xticks(np.arange(len(matrix.columns)))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(matrix.index)))
    ax.set_yticklabels(matrix.index, fontsize=8)
    ax.set_xlabel("locus_tag")
    ax.set_ylabel("family")
    ax.set_title("Family penetrance")
    if matrix.shape[0] <= 35 and matrix.shape[1] <= 20:
        for row_index, family in enumerate(matrix.index):
            for col_index, gene in enumerate(matrix.columns):
                ax.text(
                    col_index,
                    row_index,
                    f"{matrix.loc[family, gene]:.1f}",
                    ha="center",
                    va="center",
                    fontsize=7,
                    color="black",
                )
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("% family representatives hit")
    save_figure(plt, fig, outdir, "family_penetrance_heatmap")


def plot_identity_boxplot(species_hits, locus_tags, outdir):
    """
    Plot identity distributions by gene.

    :param species_hits: Hit table enriched with taxonomy columns.
    :param locus_tags: Requested locus tags.
    :param outdir: Output directory.
    """
    plt = prepare_matplotlib(outdir)
    data = [species_hits.loc[species_hits["gene"] == gene, "pident"].dropna() for gene in locus_tags]
    fig_width = max(7, len(locus_tags) * 0.8)
    fig, ax = plt.subplots(figsize=(fig_width, 5))
    try:
        ax.boxplot(data, tick_labels=locus_tags, showmeans=True, patch_artist=True)
    except TypeError:
        ax.boxplot(data, labels=locus_tags, showmeans=True, patch_artist=True)
    rng = np.random.default_rng(7)
    for index, values in enumerate(data, start=1):
        if values.empty:
            continue
        jitter = rng.normal(0, 0.04, len(values))
        ax.scatter(
            np.full(len(values), index) + jitter,
            values,
            s=14,
            alpha=0.45,
            color="#333333",
            linewidths=0,
        )
        mean = values.mean()
        sd = values.std()
        if not pd.isna(mean):
            label = f"{mean:.1f}"
            if not pd.isna(sd):
                label = f"{mean:.1f} +/- {sd:.1f}"
            ax.annotate(label, (index, mean), xytext=(0, 10), textcoords="offset points",
                        ha="center", fontsize=8, color="#111111")
    ax.set_xlabel("locus_tag")
    ax.set_ylabel("% identity")
    ax.set_title("Identity distribution by gene")
    ax.tick_params(axis="x", rotation=45)
    save_figure(plt, fig, outdir, "identity_boxplot_by_gene")


def plot_radial_tree(tree_path, species_hits, locus_tags, outdir, warnings):
    """
    Plot radial trees highlighting representatives with hits.

    :param tree_path: Path to a Newick tree.
    :param species_hits: Hit table enriched with taxonomy columns.
    :param locus_tags: Requested locus tags.
    :param outdir: Output directory.
    :param warnings: Mutable list of warning messages.
    """
    if not tree_path:
        return
    try:
        from ete3 import NodeStyle, Tree, TreeStyle
    except Exception as exc:
        warn(
            warnings,
            f"Tree plotting skipped because ete3 is not available ({exc}). "
            "Install ete3 to generate radial tree outputs.",
        )
        return

    try:
        tree = Tree(tree_path, format=1)
    except Exception as exc:
        warn(warnings, f"Tree plotting skipped because Newick could not be loaded: {exc}")
        return

    tree_tip_ids = set(tree.get_leaf_names())
    grey_style = NodeStyle()
    grey_style["fgcolor"] = "#555555"
    grey_style["hz_line_color"] = "#555555"
    grey_style["vt_line_color"] = "#555555"
    grey_style["size"] = 0
    red_style = NodeStyle()
    red_style["fgcolor"] = "#d62728"
    red_style["hz_line_color"] = "#d62728"
    red_style["vt_line_color"] = "#d62728"
    red_style["size"] = 4

    for gene in locus_tags:
        safe_gene = safe_filename(gene)
        hit_ids = set(
            species_hits.loc[
                species_hits["gene"] == gene,
                "representative_genome_id",
            ].astype(str)
        )
        found = hit_ids & tree_tip_ids
        missing = hit_ids - tree_tip_ids
        if missing:
            warn(
                warnings,
                f"{gene}: {len(found)}/{len(hit_ids)} hit representatives were found "
                f"in the tree; {len(missing)} were not found.",
            )

        for node in tree.traverse("postorder"):
            if node.is_leaf():
                node_has_hit = node.name in found
            else:
                node_has_hit = any(
                    getattr(child, "has_hit", False)
                    for child in node.children
                )
            node.add_feature("has_hit", node_has_hit)
            node.set_style(red_style if node_has_hit else grey_style)

        tree_style = TreeStyle()
        tree_style.mode = "c"
        tree_style.show_leaf_name = False
        tree_style.show_scale = False
        tree_style.force_topology = False
        tree.render(
            os.path.join(outdir, f"{safe_gene}_tree_hits_radial.svg"),
            tree_style=tree_style,
        )


def write_tsvs(outdir, species_hits, summary_gene, summary_family, summary_genus, warnings):
    """
    Write reproducible auxiliary TSV outputs.

    :param outdir: Output directory.
    :param species_hits: Hit table enriched with taxonomy columns.
    :param summary_gene: Gene-level summary DataFrame.
    :param summary_family: Family-level summary DataFrame.
    :param summary_genus: Genus-level summary DataFrame.
    :param warnings: Warning messages collected during execution.
    """
    species_hits.to_csv(
        os.path.join(outdir, "species_hits.tsv"),
        sep="\t",
        index=False,
        columns=HIT_COLUMNS + TAXONOMY_COLUMNS[1:],
    )
    summary_gene.to_csv(os.path.join(outdir, "summary_by_gene.tsv"), sep="\t", index=False)
    summary_family.to_csv(os.path.join(outdir, "summary_family.tsv"), sep="\t", index=False)
    summary_genus.to_csv(os.path.join(outdir, "summary_genus.tsv"), sep="\t", index=False)
    pd.DataFrame({"warning": warnings}).to_csv(
        os.path.join(outdir, "warnings.tsv"),
        sep="\t",
        index=False,
    )


def write_interactive_html(
    outdir,
    species_hits,
    summary_gene,
    summary_family,
    locus_tags,
    warnings,
):
    """
    Write the interactive Plotly HTML microbiome hit report.

    :param outdir: Output directory.
    :param species_hits: Hit table enriched with taxonomy columns.
    :param summary_gene: Gene-level summary DataFrame.
    :param summary_family: Family-level summary DataFrame.
    :param locus_tags: Requested locus tags.
    :param warnings: Warning messages collected during execution.
    """
    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.io import to_html

    sections = []
    plotly_loaded = False
    for index, gene in enumerate(locus_tags):
        gene_hits = species_hits[species_hits["gene"] == gene].copy()
        summary = summary_gene[summary_gene["gene"] == gene].iloc[0].to_dict()
        cards = render_cards(summary)
        top_family = summary_family[summary_family["gene"] == gene].sort_values(
            ["hit_representatives", "percent_family_representatives_hit"],
            ascending=False,
        )
        table_html = top_family.head(50).to_html(
            index=False,
            classes="summary-table",
            float_format=lambda value: f"{value:.2f}",
            border=0,
        )

        if gene_hits.empty:
            treemap_html = "<p class='empty'>No microbiome hits were found for this locus tag.</p>"
            box_html = ""
        else:
            treemap_data = build_treemap_data(gene_hits)
            fig_tree = px.treemap(
                treemap_data,
                ids="id",
                names="label",
                parents="parent",
                values="hits",
                color="mean_pident",
                color_continuous_scale="RdYlBu_r",
                hover_data={
                    "hits": True,
                    "percent_of_gene_hits": ":.2f",
                    "mean_pident": ":.2f",
                    "mean_qcovhsp": ":.2f",
                    "id": False,
                    "parent": False,
                },
            )
            fig_tree.update_layout(margin=dict(t=10, l=10, r=10, b=10))
            treemap_html = to_html(
                fig_tree,
                include_plotlyjs=True if not plotly_loaded else False,
                full_html=False,
            )
            plotly_loaded = True

            family_count = gene_hits["family"].nunique()
            if family_count <= 12:
                fig_box = px.violin(
                    gene_hits,
                    x="family",
                    y="pident",
                    color="family",
                    box=True,
                    points="all",
                    labels={"pident": "% identity", "family": "family"},
                )
            else:
                fig_box = go.Figure()
                fig_box.add_trace(
                    go.Box(
                        y=gene_hits["pident"],
                        name=gene,
                        boxmean="sd",
                        boxpoints="all",
                        jitter=0.35,
                        marker=dict(color="#4c78a8"),
                    )
                )
            fig_box.update_layout(
                showlegend=False,
                margin=dict(t=20, l=45, r=10, b=80),
                yaxis_title="% identity",
            )
            box_html = to_html(fig_box, include_plotlyjs=False, full_html=False)

        sections.append(
            f"""
            <section class="gene-section" id="section-{html.escape(gene)}"
                     style="display: {'block' if index == 0 else 'none'}">
                <h2>{html.escape(gene)}</h2>
                {cards}
                <div class="panel"><h3>Taxonomic treemap</h3>{treemap_html}</div>
                <div class="panel"><h3>Top families</h3>{table_html}</div>
                <div class="panel"><h3>Identity distribution</h3>{box_html}</div>
            </section>
            """
        )

    warning_items = "".join(f"<li>{html.escape(message)}</li>" for message in warnings)
    if not warning_items:
        warning_items = "<li>No warnings.</li>"
    options = "".join(
        f"<option value='{html.escape(gene)}'>{html.escape(gene)}</option>"
        for gene in locus_tags
    )
    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>FastTarget microbiome hit report</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 24px; color: #1f2933; }}
header {{ display: flex; justify-content: space-between; gap: 16px; align-items: center; }}
select {{ font-size: 15px; padding: 6px 10px; }}
h1 {{ font-size: 24px; margin: 0 0 16px; }}
h2 {{ margin-top: 24px; }}
.cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 10px; }}
.card {{ border: 1px solid #d8dee4; border-radius: 6px; padding: 10px; background: #f8fafc; }}
.card .label {{ color: #52606d; font-size: 12px; }}
.card .value {{ font-size: 20px; font-weight: 700; margin-top: 4px; }}
.panel {{ margin-top: 18px; border-top: 1px solid #e5e7eb; padding-top: 12px; }}
.summary-table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
.summary-table th, .summary-table td {{ border-bottom: 1px solid #e5e7eb; padding: 6px 8px; text-align: left; }}
.summary-table th {{ background: #f3f4f6; }}
.warnings {{ margin-top: 24px; background: #fff7ed; border: 1px solid #fed7aa; border-radius: 6px; padding: 10px 14px; }}
.empty {{ color: #6b7280; }}
</style>
</head>
<body>
<header>
  <h1>FastTarget microbiome hit report</h1>
  <label>locus_tag <select id="gene-selector">{options}</select></label>
</header>
{''.join(sections)}
<div class="warnings"><strong>Warnings</strong><ul>{warning_items}</ul></div>
<script>
document.getElementById('gene-selector').addEventListener('change', function(event) {{
  document.querySelectorAll('.gene-section').forEach(function(section) {{
    section.style.display = 'none';
  }});
  document.getElementById('section-' + event.target.value).style.display = 'block';
  window.dispatchEvent(new Event('resize'));
}});
</script>
</body>
</html>
"""
    with open(
        os.path.join(outdir, "microbiome_hits_interactive_report.html"),
        "w",
        encoding="utf-8",
    ) as output_file:
        output_file.write(document)


def render_cards(summary):
    """
    Render HTML summary cards for a single gene.

    :param summary: Gene-level summary values.
    :return: HTML string containing summary cards.
    """
    labels = [
        ("total_representative_hits", "representative hits"),
        ("families_with_hits", "families"),
        ("genera_with_hits", "genera"),
        ("mean_pident", "mean identity"),
        ("mean_qcovhsp", "mean qcovhsp"),
    ]
    cards = []
    for key, label in labels:
        value = summary.get(key)
        if isinstance(value, float) and not math.isnan(value):
            text = f"{value:.2f}"
        elif pd.isna(value):
            text = "NA"
        else:
            text = str(int(value)) if isinstance(value, (int, np.integer)) else str(value)
        cards.append(
            f"<div class='card'><div class='label'>{label}</div>"
            f"<div class='value'>{html.escape(text)}</div></div>"
        )
    return "<div class='cards'>" + "".join(cards) + "</div>"


def build_treemap_data(gene_hits):
    """
    Build hierarchical taxonomic data for a Plotly treemap.

    :param gene_hits: Hit table subset for a single gene.
    :return: DataFrame with treemap nodes.
    """
    rows = []
    total = gene_hits["representative_genome_id"].nunique()
    root_id = "hits"
    rows.append(
        {
            "id": root_id,
            "label": "hits",
            "parent": "",
            "hits": total,
            "percent_of_gene_hits": 100.0,
            "mean_pident": gene_hits["pident"].mean(),
            "mean_qcovhsp": gene_hits["qcovhsp"].mean(),
        }
    )
    for depth, rank in enumerate(TAXONOMY_RANKS):
        group_cols = TAXONOMY_RANKS[: depth + 1]
        grouped = gene_hits.groupby(group_cols, dropna=False)
        for keys, group in grouped:
            if not isinstance(keys, tuple):
                keys = (keys,)
            node_id = "|".join(keys)
            parent = root_id if depth == 0 else "|".join(keys[:-1])
            hits = group["representative_genome_id"].nunique()
            rows.append(
                {
                    "id": node_id,
                    "label": keys[-1],
                    "parent": parent,
                    "hits": hits,
                    "percent_of_gene_hits": percentage(hits, total),
                    "mean_pident": group["pident"].mean(),
                    "mean_qcovhsp": group["qcovhsp"].mean(),
                }
            )
    return pd.DataFrame(rows)


def parse_args():
    """
    Parse command line arguments.

    :return: Parsed command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Generate post-run FastTarget MGnify microbiome hit plots."
    )
    parser.add_argument("--output-path", required=True)
    parser.add_argument("--organism-name", required=True)
    parser.add_argument("--databases-path", required=True)
    parser.add_argument("--catalogue", required=True)
    parser.add_argument("--locus-tags", required=True, nargs="+")
    parser.add_argument("--outdir")
    parser.add_argument("--tree")
    parser.add_argument("--top-threshold", type=float, default=10.0)
    parser.add_argument("--min-hits", type=int, default=5)
    parser.add_argument("--min-taxon-coverage", type=float, default=10.0)
    return parser.parse_args()


def main():
    """
    Run the microbiome hit plotting command line workflow.
    """
    args = parse_args()
    outdir = args.outdir or output_default(args)
    os.makedirs(outdir, exist_ok=True)
    warnings = []

    hits_parquet = hits_path(args)
    tax_parquet = taxonomy_path(args)
    if not os.path.isfile(hits_parquet):
        fail(f"Consolidated microbiome hits parquet not found: {hits_parquet}")
    if not os.path.isfile(tax_parquet):
        fail(f"Representative taxonomy parquet not found: {tax_parquet}")
    if args.tree and not os.path.isfile(args.tree):
        warn(warnings, f"Tree file was provided but does not exist: {args.tree}")
        args.tree = None

    locus_tags = list(dict.fromkeys(args.locus_tags))
    hits = load_hits_for_loci(hits_parquet, locus_tags)
    taxonomy = load_taxonomy(tax_parquet)
    species_hits = build_species_hits_table(hits, taxonomy, locus_tags, warnings)
    summary_gene = summarize_by_gene(species_hits, locus_tags)
    summary_family = summarize_rank(species_hits, taxonomy, locus_tags, "family")
    summary_genus = summarize_rank(species_hits, taxonomy, locus_tags, "genus")

    write_tsvs(outdir, species_hits, summary_gene, summary_family, summary_genus, warnings)
    plot_total_hits(summary_gene, outdir)
    plot_composition_stacked_bar(
        summary_family,
        "family",
        locus_tags,
        args.top_threshold,
        outdir,
    )
    plot_composition_stacked_bar(
        summary_genus,
        "genus",
        locus_tags,
        args.top_threshold,
        outdir,
    )
    plot_family_penetrance_heatmap(
        summary_family,
        locus_tags,
        args.min_hits,
        args.min_taxon_coverage,
        outdir,
    )
    plot_identity_boxplot(species_hits, locus_tags, outdir)
    plot_radial_tree(args.tree, species_hits, locus_tags, outdir, warnings)
    if warnings:
        pd.DataFrame({"warning": warnings}).to_csv(
            os.path.join(outdir, "warnings.tsv"),
            sep="\t",
            index=False,
        )
    write_interactive_html(
        outdir,
        species_hits,
        summary_gene,
        summary_family,
        locus_tags,
        warnings,
    )
    print(f"Microbiome hit plots written to: {outdir}")


if __name__ == "__main__":
    main()
