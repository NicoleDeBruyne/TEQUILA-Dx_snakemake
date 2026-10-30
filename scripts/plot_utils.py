# plotting helpers shared by the group merge scripts

import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np


# outlier genes per sample, one box per outlier type
def plot_outlier_counts(samples, dfs, legends, suptitle, outfile):

    counts_list = []
    for df in dfs:
        counts_list.append(df.groupby('sample')['gene'].nunique().sort_index())

    sample_order = sorted(samples)
    for i, counts in enumerate(counts_list):
        counts_list[i] = counts.reindex(sample_order, fill_value=0)

    barplot_file = os.path.splitext(outfile)[0] + "_barplot.pdf"
    fig, axes = plt.subplots(len(dfs), 1, figsize=(max(16, len(sample_order)*0.2), max(6, len(dfs)*4)), sharex=True)
    if len(dfs) == 1:
        axes = [axes]

    for ax, counts, (color, label) in zip(axes, counts_list, legends):
        ax.bar(counts.index, counts.values, color=color)
        ax.set_title(f'Average={counts.mean():.2f}, Median={counts.median()}, Range={counts.min()}–{counts.max()}')
        ax.set_ylabel('Count')
        ax.tick_params(axis='x', rotation=90)
        ax.legend(handles=[mpatches.Patch(color=color, label=label)],
                  loc='center left', bbox_to_anchor=(1.0, 0.5))

    axes[-1].set_xlabel('Sample')
    fig.suptitle(suptitle, fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(barplot_file)
    plt.close()
    print(f"Bar plot saved to {barplot_file}")

    boxplot_file = os.path.splitext(outfile)[0] + "_boxplot.pdf"
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    data = [counts.values for counts in counts_list]
    colors = [c for c, _ in legends]

    ax.boxplot(
        data,
        patch_artist=False,
        medianprops=dict(color='black'),
        boxprops=dict(color='black'),
        whiskerprops=dict(color='black'),
        capprops=dict(color='black'),
        showfliers=False
    )

    jitter_strength = 0.08
    for i, counts in enumerate(data):
        x = np.full(len(counts), i + 1) + np.random.uniform(-jitter_strength, jitter_strength, len(counts))
        ax.scatter(x, counts, color=colors[i], alpha=0.7, s=20)

    ax.set_xticks(range(1, len(dfs) + 1))
    ax.set_xticklabels([])

    ax.set_ylabel('Number of genes')
    ax.set_title(suptitle)

    legend_handles = []
    for counts, color, (_, label) in zip(counts_list, colors, legends):
        stats_text = f"{label}\nAvg={counts.mean():.2f}, Med={counts.median()}, Range={counts.min()}–{counts.max()}"
        legend_handles.append(mpatches.Patch(color=color, label=stats_text))
    ax.legend(handles=legend_handles, loc='center left', bbox_to_anchor=(1.05, 0.5))

    plt.tight_layout()
    plt.savefig(boxplot_file, bbox_inches='tight')
    plt.close()
    print(f"Box plot saved to {boxplot_file}")
