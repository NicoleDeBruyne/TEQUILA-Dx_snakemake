"""
tier_rules.py

Outlier-junction strength and candidate-gene tiers, shared by merge_hits.py (which assigns them)
and make_sample_hit_report.py (which draws the tier x hit-category grid and greys out the cells the
rules can never produce, via possible_cells()).

Junction strength depends only on the size / allelic pattern of the splicing change:
  Strong   : bulk |dPSI| >= 0.5, or both haplotypes >= 0.5
  Moderate : bulk |dPSI| >= 0.2, or one haplotype >= 0.5
  Weak     : any other outlier change

Tiers (1 = almost certainly a diagnosis ... 6 = VUS without RNA dysregulation). Genes with any
dominant pattern (AD or XLD, including mixed AD/AR genes) use the dominant rules; AR, XLR and
genes without a listed inheritance use the recessive rules. A hemizygous variant (male chrX) is
called 1/1 and so counts as homozygous.

  AD / XLD
    any P/LP                                                   -> 1
    ASE in a haploinsufficient gene, or Strong / Moderate      -> 2
    Weak                                                       -> 3
    ASE in a non-haploinsufficient gene                        -> 4
    any VUS, no RNA dysregulation                              -> 6

  AR / XLR / unknown
    2 P/LP in trans, or a homozygous P/LP                      -> 1
    P/LP + VUS in trans, 2 P/LP phase unclear, or Strong       -> 2
    1 P/LP, 2 VUS in trans, or a homozygous VUS, and:
        ASE or Moderate -> 3;  Weak -> 4;  no RNA dysregulation -> 5
    otherwise (another VUS, or no variant), and:
        ASE or Moderate -> 4;  Weak -> 5;  VUS without RNA dysregulation -> 6
"""

STRONG, MODERATE, WEAK, NONE = "Strong", "Moderate", "Weak", "None"


def junction_strength(max_bulk, max_hap1, max_hap2):
    """Strength label from the largest |dPSI| in bulk / hap1 / hap2 (NaN = no outlier there)."""
    def ge(x, t):
        return x == x and x >= t        # x == x is False for NaN

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


def assign_tier(vc, ase, haploinsufficient, strength, dominant):
    """vc: variant classification (dict with has_pathogenic, has_vus, has_homozygous_pathogenic,
    has_homozygous_vus, has_2_pathogenic_trans, has_2_pathogenic_unclear,
    has_1_pathogenic_1_vus_trans, has_2_vus_trans). strength: Strong / Moderate / Weak / None (the
    stronger of the GTEx and cohort calls). Returns 1-6, or None (not a candidate)."""
    strong, moderate, weak = strength == STRONG, strength == MODERATE, strength == WEAK
    if dominant:
        if vc["has_pathogenic"]:
            return 1
        if (ase and haploinsufficient) or strong or moderate:
            return 2
        if weak:
            return 3
        if ase:
            return 4
        if vc["has_vus"]:
            return 6
        return None
    if vc["has_2_pathogenic_trans"] or vc["has_homozygous_pathogenic"]:
        return 1
    if vc["has_1_pathogenic_1_vus_trans"] or vc["has_2_pathogenic_unclear"] or strong:
        return 2
    if vc["has_pathogenic"] or vc["has_2_vus_trans"] or vc["has_homozygous_vus"]:
        if ase or moderate:
            return 3
        if weak:
            return 4
        return 5
    if ase or moderate:
        return 4
    if weak:
        return 5
    if vc["has_vus"]:
        return 6
    return None


def strongest(*labels):
    order = {STRONG: 3, MODERATE: 2, WEAK: 1}
    best = max(labels, key=lambda l: order.get(l, 0))
    return best if best in order else NONE


# variant genotypes a gene can have, as variant classifications
_EMPTY = dict(has_pathogenic=False, has_vus=False, has_homozygous_pathogenic=False, has_homozygous_vus=False,
              has_2_pathogenic_trans=False, has_2_pathogenic_unclear=False,
              has_1_pathogenic_1_vus_trans=False, has_2_vus_trans=False)
_GENOTYPES = [
    {},                                                                            # no variant
    dict(has_vus=True),                                                            # 1 VUS
    dict(has_vus=True, has_homozygous_vus=True),                                   # homozygous VUS
    dict(has_vus=True, has_2_vus_trans=True),                                      # 2 VUS in trans
    dict(has_pathogenic=True),                                                     # 1 P/LP
    dict(has_pathogenic=True, has_homozygous_pathogenic=True),                     # homozygous P/LP
    dict(has_pathogenic=True, has_vus=True),                                       # P/LP + VUS (cis / unclear)
    dict(has_pathogenic=True, has_vus=True, has_1_pathogenic_1_vus_trans=True),    # P/LP + VUS in trans
    dict(has_pathogenic=True, has_2_pathogenic_trans=True),                        # 2 P/LP in trans
    dict(has_pathogenic=True, has_2_pathogenic_unclear=True),                      # 2 P/LP, phase unclear
]


def possible_cells():
    """Set of (hit categories as a frozenset of 'variant' / 'ASE' / 'junction', tier) that the rules
    can produce, over every genotype, ASE status, haploinsufficiency, junction strength and
    inheritance group."""
    cells = set()
    for g in _GENOTYPES:
        vc = {**_EMPTY, **g}
        has_variant = vc["has_pathogenic"] or vc["has_vus"]
        for ase in (False, True):
            for hapi in (False, True):
                for strength in (STRONG, MODERATE, WEAK, NONE):
                    for dominant in (False, True):
                        tier = assign_tier(vc, ase, hapi, strength, dominant)
                        cats = frozenset(c for c, on in (("variant", has_variant), ("ASE", ase),
                                                         ("junction", strength != NONE)) if on)
                        if tier is not None and cats:
                            cells.add((cats, tier))
    return cells
