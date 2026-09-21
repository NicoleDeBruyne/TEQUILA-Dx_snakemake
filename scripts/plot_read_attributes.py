
import argparse

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42
import seaborn as sns

from sample_alias import add_alias_map_arg, parse_alias_map, resolve, resolve_all


def parse_args():
    parser = argparse.ArgumentParser(description="Merge per-sample read-length summaries and plot a boxplot across a cohort")
    parser.add_argument("--infiles", nargs="+", required=True,
        help="Every sample's {sample}_read_attributes.tsv, in the desired plot order")
    parser.add_argument("--outprefix", required=True,
        help="Prefix for output files: <outprefix>_read_attributes.tsv, _read_lengths.pdf")
    parser.add_argument('--title')
    add_alias_map_arg(parser)
    return parser.parse_args()


def plot_boxplot(df, sample_order, target_order, palette_dict, outfile, title):
    fig, ax = plt.subplots(figsize=(max(6, len(sample_order) * 0.5), 6))

    n_types = max(len(target_order), 1)
    box_width = 0.7 / n_types
    stats_list, positions, colors = [], [], []
    for i, sample in enumerate(sample_order):
        sub = df[df['sample'] == sample]
        for j, tt in enumerate(target_order):
            row = sub[sub['targetType'] == tt]
            if row.empty:
                continue
            row = row.iloc[0]
            stats_list.append(dict(
                med=row['median'], q1=row['q1'], q3=row['q3'],
                whislo=row['min'], whishi=row['max'], fliers=[],
            ))
            offset = (j - (n_types - 1) / 2) * box_width
            positions.append(i + offset)
            colors.append(palette_dict[tt])

    if stats_list:
        bxp_artists = ax.bxp(stats_list, positions=positions, widths=box_width * 0.9,
                              patch_artist=True, showfliers=False)
        for patch, color in zip(bxp_artists['boxes'], colors):
            patch.set_facecolor(color)

    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=palette_dict[tt]) for tt in target_order]
    ax.legend(handles, target_order, title="TargetType")

    ax.set_xticks(range(len(sample_order)))
    ax.set_xticklabels(sample_order, rotation=45, ha="right")
    ax.set_ylabel("ReadLength")
    ax.set_xlabel("Sample")
    ymax = max((s['whishi'] for s in stats_list), default=1)
    ax.set_ylim(0, ymax * 1.05)
    ax.margins(y=0.05)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(outfile)
    print(f"\u2705 Read length box plot saved to {outfile}")


def main():
    args = parse_args()

    per_sample_dfs = [pd.read_csv(f, sep="\t") for f in args.infiles]
    df = pd.concat(per_sample_dfs, ignore_index=True)

    sample_order = [d['sample'].iloc[0] for d in per_sample_dfs if not d.empty]

    default_target_order = ["on_target", "off_target", "mapped", "unmapped"]
    target_order = default_target_order + [
        t for t in df['targetType'].drop_duplicates() if t not in default_target_order
    ]
    target_order = [t for t in target_order if t in df['targetType'].unique()]

    out_tsv = f"{args.outprefix}_read_attributes.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\u2705 Attributes saved to {out_tsv}")

    palette_dict = dict(zip(target_order, sns.color_palette("pastel")))

    print("\nPlotting read length boxplot...")
    plot_boxplot(df, sample_order, target_order, palette_dict, f"{args.outprefix}_read_lengths.pdf", args.title)

    alias_map = parse_alias_map(args.alias_map)
    alias_df = df.copy()
    alias_df['sample'] = alias_df['sample'].apply(lambda s: resolve(s, alias_map))
    alias_sample_order = resolve_all(sample_order, alias_map)
    plot_boxplot(alias_df, alias_sample_order, target_order, palette_dict,
                 f"{args.outprefix}_read_lengths_alias.pdf", args.title)


if __name__ == "__main__":
    main()
