"""
Shared plot colours for the pipeline (Nature journal palette).

SELECT_ORDER is the order colours are handed out in: categorical plots take colours from the
start of the list, and sample types not listed in config's sample_type_colors take the next
unused ones (see sample_type_color in the Snakefile, which mirrors default_type_colors below).
Continuous scales run red (#C23637) - white - blue (#0068A9).
"""
from matplotlib.colors import LinearSegmentedColormap, to_hex, to_rgb

SELECT_ORDER = [
    "#c23637", "#009099", "#9d4588", "#c69528", "#3d892e", "#ea6302", "#0068a9",
    "#d95d5b", "#42b4b5", "#b271ab", "#e8c048", "#57aa3e", "#f48f3e", "#4c8fca",
    "#ea9a9c", "#8cccce", "#cca0ca", "#f4d77e", "#95c36e", "#fbb875", "#91c4e9",
    "#f8c9c6", "#c5e3ea", "#e8d1e6", "#ffefbb", "#d3e4be", "#fddbb4", "#c1e3fa",
    "#912321", "#005d6e", "#6e2769", "#926d17", "#1e662a", "#ae450b", "#002b58",
]

# colour families (strongest first), for graded categories
RED_FAMILY = ["#c23637", "#d95d5b", "#ea9a9c", "#f8c9c6"]
TEAL_FAMILY = ["#009099", "#42b4b5", "#8cccce", "#c5e3ea"]
PURPLE_FAMILY = ["#9d4588", "#b271ab", "#cca0ca", "#e8d1e6"]
GOLD_FAMILY = ["#c69528", "#e8c048", "#f4d77e", "#ffefbb"]

GREY = "#a0a0a0"        # "other" / unmapped categories
LIGHT_GREY = "#c4c4c4"  # off-target / unmapped bar segments
RED, BLUE = "#c23637", "#0068a9"

# red (low) - white - blue (high); use with vmin/vmax (or a TwoSlopeNorm) so white is the centre
RED_WHITE_BLUE = LinearSegmentedColormap.from_list("red_white_blue", [RED, "#ffffff", BLUE])
# blue (low) - white - red (high)
BLUE_WHITE_RED = LinearSegmentedColormap.from_list("blue_white_red", [BLUE, "#ffffff", RED])


def select(n):
    """The first n colours of SELECT_ORDER (cycling)."""
    return [SELECT_ORDER[i % len(SELECT_ORDER)] for i in range(n)]


def default_type_colors(names, configured=None):
    """name -> colour: configured colours first, every other name the next unused SELECT_ORDER colour."""
    configured = {k: v for k, v in (configured or {}).items() if v}
    used = {c.lower() for c in configured.values()}
    free = [c for c in SELECT_ORDER if c.lower() not in used] or list(SELECT_ORDER)
    out, k = {}, 0
    for n in names:
        if n in configured:
            out[n] = configured[n]
        elif n not in out:
            out[n] = free[k % len(free)]
            k += 1
    return out


def darker(color, f=0.5):
    """`color` mixed with black (f = fraction of black)."""
    r, g, b = to_rgb(color)
    return to_hex((r * (1 - f), g * (1 - f), b * (1 - f)))
