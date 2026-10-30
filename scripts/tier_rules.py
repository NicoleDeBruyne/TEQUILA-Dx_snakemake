# rules for evidence strength and gene tiers (tier 1 = strongest)

STRONG, MODERATE, WEAK, NONE = "Strong", "Moderate", "Weak", "None"


# Strong/Moderate/Weak/None from the biggest delta-PSI across bulk, hap1 and hap2
def junction_strength(max_bulk, max_hap1, max_hap2):
    def ge(x, t):
        return x == x and x >= t

    if ge(max_bulk, 0.5) or (ge(max_hap1, 0.5) and ge(max_hap2, 0.5)):
        return STRONG
    if ge(max_bulk, 0.2) or ge(max_hap1, 0.5) or ge(max_hap2, 0.5):
        return MODERATE
    if ge(max_bulk, 1e-300) or ge(max_hap1, 1e-300) or ge(max_hap2, 1e-300):
        return WEAK
    return NONE


def is_dominant(inheritance_patterns):
    s = str(inheritance_patterns or "")
    return ("AD" in s) or ("XLD" in s)


# the main tiering logic. isoform evidence only ever counts as weak on purpose
def assign_tier(vc, ase, haploinsufficient, strength, dominant, isoform=False):
    if isoform:
        strength = strongest(strength, WEAK)
    strong, moderate, weak = strength == STRONG, strength == MODERATE, strength == WEAK
    if dominant:
        if vc["has_pathogenic"]:
            return 1
        if (ase and haploinsufficient) or strong or moderate:
            return 2
        if weak:
            return 3
        if ase:
            return 5
        if vc["has_vus"]:
            return 6
        return None
    if vc["has_2_pathogenic_trans"] or vc["has_homozygous_pathogenic"]:
        return 1
    if vc["has_1_pathogenic_1_vus_trans"] or vc["has_2_pathogenic_unclear"] or strong:
        return 2
    if vc["has_pathogenic"] or vc["has_2_vus_trans"] or vc["has_homozygous_vus"]:
        if moderate:
            return 2
        if ase or weak:
            return 3
        return 4
    if vc["has_vus"]:
        if moderate:
            return 3
        if ase or weak:
            return 4
        return 6
    if moderate:
        return 4
    if ase or weak:
        return 5
    return None


def strongest(*labels):
    order = {STRONG: 3, MODERATE: 2, WEAK: 1}
    best = max(labels, key=lambda l: order.get(l, 0))
    return best if best in order else NONE


# walk every evidence combo so the tier heatmaps can show cells that never got any genes
_EMPTY = dict(has_pathogenic=False, has_vus=False, has_homozygous_pathogenic=False, has_homozygous_vus=False,
              has_2_pathogenic_trans=False, has_2_pathogenic_unclear=False,
              has_1_pathogenic_1_vus_trans=False, has_2_vus_trans=False)
_GENOTYPES = [
    {},
    dict(has_vus=True),
    dict(has_vus=True, has_homozygous_vus=True),
    dict(has_vus=True, has_2_vus_trans=True),
    dict(has_pathogenic=True),
    dict(has_pathogenic=True, has_homozygous_pathogenic=True),
    dict(has_pathogenic=True, has_vus=True),
    dict(has_pathogenic=True, has_vus=True, has_1_pathogenic_1_vus_trans=True),
    dict(has_pathogenic=True, has_2_pathogenic_trans=True),
    dict(has_pathogenic=True, has_2_pathogenic_unclear=True),
]


def possible_cells():
    cells = set()
    for g in _GENOTYPES:
        vc = {**_EMPTY, **g}
        has_variant = vc["has_pathogenic"] or vc["has_vus"]
        for ase in (False, True):
            for hapi in (False, True):
                for strength in (STRONG, MODERATE, WEAK, NONE):
                    for dominant in (False, True):
                        for isoform in (False, True):
                            tier = assign_tier(vc, ase, hapi, strength, dominant, isoform)
                            cats = frozenset(c for c, on in (("variant", has_variant), ("ASE", ase),
                                                             ("junction", strength != NONE),
                                                             ("isoform", isoform)) if on)
                            if tier is not None and cats:
                                cells.add((cats, tier))
    return cells
