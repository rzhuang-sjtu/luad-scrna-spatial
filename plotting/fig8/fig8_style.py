"""Shared plotting style for Fig 8 + S11 — translates the project R style guide
to matplotlib/seaborn (sans 8pt, no titles, ANP-like palette, DPI 300, white
heatmap stroke, viridis/magma spatial, RdBu diverging).

Usage:
    from fig8_style import setup, COL, save_panel
    setup()
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    ...
    save_panel(fig, OUT/"fig8x_panel")   # writes both .pdf and .png at 300 dpi
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

# ST gene-expression colormap: blue → red sequential, no white middle
# (white middle would blend into H&E pink background and lose mid-range spots).
CMAP_ST = LinearSegmentedColormap.from_list(
    "st_blue_red",
    ["#4DBBD5", "#7C95B8", "#A07499", "#C7547A", "#E64B35"],
    N=256,
)

# ---- palette (ANP / Lancet Oncology family, matched to R guide) ----
COL = {
    # Macrophage subtypes
    "Macro_C1QC": "#4DBBD5", "Macro_FCN1": "#E64B35", "Macro_FOLR2": "#00A087",
    "Macro_MARCO": "#3C5488", "Macro_SPP1": "#F39B7F",
    "Macro_general": "#8491B4", "Macro_prolif": "#91D1C2",
    # Malignant programs
    "MP1": "#E64B35", "MP2": "#4DBBD5", "MP3": "#00A087", "MP4": "#3C5488",
    # Neutrophil subtypes
    "Neu_Inflammatory": "#E64B35", "Neu_Angiogenic": "#F39B7F",
    "Neu_Metastatic": "#3C5488", "Neu_ECM_remodeling": "#4DBBD5",
    "Neu_OSM_priming": "#00A087", "Neu_OSM_low": "#8491B4",
    "Neu_IFN_response": "#91D1C2", "Neu_unclassified": "#D9D9D9",
    # Generic semantic pairs (consistent across all Fig 8 panels)
    "tumor": "#E64B35",   "normal": "#4DBBD5",
    "high":  "#E64B35",   "low":    "#4DBBD5",
    "NR":    "#E64B35",   "R":      "#4DBBD5",
    "MPR":   "#4DBBD5",   "NMPR":   "#E64B35",
    "ROI":   "#E64B35",   "nonROI": "#4DBBD5",
    "LUAD":  "#E64B35",   "other":  "#8491B4",
    # Hazard direction (Cox forest)
    "hr_up": "#E64B35", "hr_down": "#4DBBD5",
    # Reference threshold
    "ref_red": "#cb181d",
    # Annotation grays
    "grid":  "#cccccc", "rule": "black",
    # Transitions for Venn
    "venn_macro": "#4DBBD5", "venn_mal": "#E64B35", "venn_neu": "#00A087",
}

def _ensure_arial():
    """Make sure matplotlib actually loaded Arial (otherwise picks fallback)."""
    import matplotlib.font_manager as fm
    arial_paths = [p for p in fm.findSystemFonts() if "arial" in p.lower()]
    for p in arial_paths:
        try: fm.fontManager.addfont(p)
        except Exception: pass


def setup(base_size=8):
    """Apply the style. Call once at top of every plotting script."""
    _ensure_arial()
    sns.set_theme(style="ticks", context="paper")
    mpl.rcParams.update({
        # font — Arial first per project guide
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans", "Helvetica"],
        "font.size": base_size,
        "axes.titlesize": base_size,            # we'll set titles to "" manually anyway
        "axes.labelsize": base_size,
        "xtick.labelsize": base_size - 1,
        "ytick.labelsize": base_size - 1,
        "legend.fontsize": base_size - 1,
        "legend.title_fontsize": base_size - 1,
        # axis lines
        "axes.linewidth": 0.4,
        "axes.edgecolor": "black",
        "axes.labelcolor": "black",
        "xtick.color": "black", "ytick.color": "black",
        "xtick.major.width": 0.3, "ytick.major.width": 0.3,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        # remove top/right spines
        "axes.spines.top": False, "axes.spines.right": False,
        # legend & save
        "legend.frameon": False,
        "savefig.dpi": 300, "savefig.bbox": "tight",
        "pdf.fonttype": 42, "ps.fonttype": 42,    # editable text in vector outputs
        "figure.dpi": 120,
    })


def save_panel(fig, stem, also_png=True):
    """Save fig to <stem>.pdf and <stem>.png at 300 dpi (no `device=cairo_pdf`)."""
    fig.savefig(f"{stem}.pdf")
    if also_png:
        fig.savefig(f"{stem}.png", dpi=300)


def sig_stars(p):
    """Return *** / ** / * / ns string."""
    if p is None or p != p:   # NaN
        return ""
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"


def add_subtitle(ax, text, y=1.04, **kw):
    """Add a small italic subtitle above the axes (no plot title).
    Use larger y to lift it off any cohort-level header text."""
    ax.text(0.0, y, text, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=7, style="italic", color="black", **kw)


def style_axes(ax, ylabel=None, xlabel=None):
    """Apply consistent axis styling. Removes top/right spine, sets labels."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if ylabel is not None: ax.set_ylabel(ylabel)
    if xlabel is not None: ax.set_xlabel(xlabel)
    ax.set_title("")        # always blank — panel labels added in layout software
