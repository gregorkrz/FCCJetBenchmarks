"""Stage 2 of the event displays: draw one multi-page PDF in the eta-phi plane.

Reads the payload written by src/event_displays.py. Pure numpy + matplotlib, so
it runs on the login node and can be iterated on without re-clustering.

One page per event:
  - stable gen particles, coloured by PID group, marker *area* proportional to pT
  - the gen jets: axis marker, a convex hull around each jet's constituents, and
    a label; the jets matched to a Higgs parton are marked
  - the hard Higgs partons, with a connector to the jet each matched to
  - an annotation box with the mH values, the per-jet table and the invisible energy

Usage:
    python src/plotting/event_display_plots.py --payload <payload.pkl>
"""
import argparse
import os
import pickle
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Polygon

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from src.process_config import HUMAN_READABLE_PROCESS_NAMES  # noqa: E402

# stable_particles() cuts at |eta| < 2.56 (functions.h), so nothing can appear
# outside it; the limit is drawn so the acceptance edge is visible.
ETA_ACCEPTANCE = 2.56

# matplotlib's scatter `s` *is* area in pt^2, so area proportional to pT is
# simply s = AREA_PER_GEV * pt. Clipped so sub-GeV particles stay visible and a
# 50 GeV particle does not swamp the panel. Shared across all pages.
AREA_PER_GEV = 25.0
AREA_MIN, AREA_MAX = 4.0, 1200.0
PT_LEGEND_VALUES = (0.5, 2.0, 10.0, 40.0)

# Durham clustering has no radius: it assigns *every* particle to one of the N
# jets, so a jet's full constituent set can stretch across the whole eta-phi
# plane and its convex hull then says nothing useful. The hull is therefore drawn
# around the jet's energy core - the smallest set of constituents, taken in
# descending pT, carrying this fraction of the jet's constituent pT.
HULL_CORE_FRACTION = 0.90

# Colour = PID group (as requested), marker shape = the same grouping again, so
# the encoding survives colour-vision deficiency and greyscale printing.
# Neutrinos never appear: particle_filter.py calls stable_particles(Particle,
# true), which drops them. Their energy goes in the annotation box instead.
PID_GROUPS_FIVE = [
    ("Charged hadrons", "#2a78d6", "o"),
    ("Photons", "#eda100", "*"),
    ("Neutral hadrons", "#1baf7a", "s"),
    ("Electrons", "#e34948", "^"),
    ("Muons", "#4a3aa7", "v"),
]
PID_GROUPS_FOUR = [
    ("Charged hadrons", "#2a78d6", "o"),
    ("Photons", "#eb6834", "*"),
    ("Neutral hadrons", "#1baf7a", "s"),
    ("Leptons (e/mu)", "#4a3aa7", "^"),
]

MARKER_EDGE = "#3f3f3f"


def pid_group_index(pdg, groups):
    """Index into `groups` for one PDG code."""
    a = abs(int(pdg))
    five = len(groups) == 5
    if a == 22:
        return 1
    if a == 11:
        return 3
    if a == 13:
        return 4 if five else 3
    # Everything else is a hadron; split by whether it carries charge.
    neutral = a in (111, 130, 310, 311, 2112, 3122, 3212, 3322)
    return 2 if neutral else 0


def wrap_pi(x):
    """Wrap an angle difference into (-pi, pi]."""
    return (np.asarray(x) + np.pi) % (2 * np.pi) - np.pi


def convex_hull(points):
    """Andrew monotone chain. Returns the hull vertices, or None if degenerate."""
    pts = sorted(map(tuple, points))
    if len(pts) < 3:
        return None

    def half(seq):
        out = []
        for p in seq:
            while len(out) >= 2:
                (ox, oy), (ax, ay) = out[-2], out[-1]
                if (ax - ox) * (p[1] - oy) - (ay - oy) * (p[0] - ox) > 0:
                    break
                out.pop()
            out.append(p)
        return out

    lower, upper = half(pts), half(reversed(pts))
    hull = lower[:-1] + upper[:-1]
    return np.array(hull) if len(hull) >= 3 else None


def marker_area(pt):
    return np.clip(AREA_PER_GEV * np.asarray(pt, dtype=float), AREA_MIN, AREA_MAX)


def core_mask(pt, fraction=HULL_CORE_FRACTION):
    """Boolean mask of the hardest constituents carrying `fraction` of the pT."""
    pt = np.asarray(pt, dtype=float)
    if pt.size == 0 or pt.sum() <= 0:
        return np.ones(pt.size, dtype=bool)
    order = np.argsort(-pt)
    cumulative = np.cumsum(pt[order]) / pt.sum()
    n_keep = int(np.searchsorted(cumulative, fraction) + 1)
    mask = np.zeros(pt.size, dtype=bool)
    mask[order[:min(n_keep, pt.size)]] = True
    return mask


def draw_event(event, meta, window, page_index, groups):
    fig = plt.figure(figsize=(11.7, 8.3))
    grid = fig.add_gridspec(1, 2, width_ratios=[2.15, 1], wspace=0.04,
                            left=0.06, right=0.99, top=0.93, bottom=0.08)
    ax = fig.add_subplot(grid[0, 0])
    side = fig.add_subplot(grid[0, 1])
    side.axis("off")

    eta = np.asarray(event["part_eta"], dtype=float)
    phi = np.asarray(event["part_phi"], dtype=float)
    pt = np.asarray(event["part_pt"], dtype=float)
    pdg = np.asarray(event["part_pdg"], dtype=int)
    assign = np.asarray(event["part_jet_index"], dtype=int)

    jet_eta = np.asarray(event["jet_eta"], dtype=float)
    jet_phi = np.asarray(event["jet_phi"], dtype=float)
    jet_e = np.asarray(event["jet_energy"], dtype=float)
    jet_pt = np.asarray(event["jet_pt"], dtype=float)
    higgs_jets = [int(j) for j in event["parton_to_jet"] if int(j) >= 0]
    n_jets = min(meta["n_jets"], len(jet_eta))

    # --- acceptance edges
    for edge in (-ETA_ACCEPTANCE, ETA_ACCEPTANCE):
        ax.axvline(edge, color="0.6", lw=0.7, ls=(0, (3, 3)), zorder=0)
    ax.text(ETA_ACCEPTANCE, np.pi * 0.99, r" $|\eta|=2.56$", fontsize=7,
            color="0.45", ha="left", va="top", rotation=90)

    # --- jet hulls, in a phi frame unwrapped around each jet axis so a jet
    # straddling the +-pi seam is still enclosed. The +-2pi copies are drawn too
    # and clipped by the axes, which is what makes the seam case look right.
    for jet in range(n_jets):
        mask = assign == jet
        if not np.any(mask):
            continue
        # Restrict to the jet's energy core, else an unbounded Durham jet's hull
        # sprawls across the plane (see HULL_CORE_FRACTION).
        core = np.zeros_like(mask)
        core[np.flatnonzero(mask)[core_mask(pt[mask])] ] = True
        unwrapped = jet_phi[jet] + wrap_pi(phi[core] - jet_phi[jet])
        hull = convex_hull(np.column_stack([eta[core], unwrapped]))
        is_higgs = jet in higgs_jets
        style = dict(
            facecolor="0.55", alpha=0.13 if is_higgs else 0.05,
            edgecolor="0.30", lw=1.8 if is_higgs else 0.9,
            ls="-" if is_higgs else (0, (4, 3)), zorder=1,
        )
        for shift in (-2 * np.pi, 0.0, 2 * np.pi):
            if hull is not None:
                ax.add_patch(Polygon(hull + [0.0, shift], closed=True, **style))
            else:
                ax.add_patch(Circle((eta[core].mean(), unwrapped.mean() + shift),
                                    0.15, **style))

    # --- particles: face colour = PID, edge colour = jet, area = pT.
    # Hues are reserved for PID, so jet membership uses the edge and the hull.
    jet_edge_colors = plt.cm.tab10(np.linspace(0, 0.9, max(1, n_jets)))
    order = np.argsort(-pt)  # hard particles drawn last, so never buried
    for gi, (name, color, marker) in enumerate(groups):
        sel = np.array([pid_group_index(p, groups) == gi for p in pdg[order]])
        if not np.any(sel):
            continue
        idx = order[sel]
        edges = np.array([
            jet_edge_colors[assign[i] % len(jet_edge_colors)] if assign[i] >= 0
            else (0.25, 0.25, 0.25, 1.0)
            for i in idx
        ])
        ax.scatter(eta[idx], phi[idx], s=marker_area(pt[idx]), marker=marker,
                   facecolor=color, edgecolors=edges, linewidths=0.8,
                   alpha=0.85, zorder=3)
    # Particles in no kept jet (should be none) flagged explicitly.
    lost = assign < 0
    if np.any(lost):
        ax.scatter(eta[lost], phi[lost], s=marker_area(pt[lost]), marker="x",
                   color="k", linewidths=0.9, zorder=4)

    # --- jet axes and labels
    for jet in range(n_jets):
        is_higgs = jet in higgs_jets
        ax.scatter([jet_eta[jet]], [jet_phi[jet]], s=260, marker="P",
                   facecolor="none", edgecolors="0.15", linewidths=1.8, zorder=5)
        if is_higgs:
            ax.scatter([jet_eta[jet]], [jet_phi[jet]], s=560, marker="o",
                       facecolor="none", edgecolors="0.15", linewidths=1.4, zorder=5)
        label = f"j{jet}  {jet_e[jet]:.1f} GeV"
        ax.annotate(
            label + ("  [H]" if is_higgs else ""),
            (jet_eta[jet], jet_phi[jet]), textcoords="offset points",
            xytext=(9, 9), fontsize=7.5, zorder=6,
            bbox=dict(boxstyle="round,pad=0.2", fc="0.92" if is_higgs else "white",
                      ec="0.6", lw=0.5, alpha=0.9),
        )

    # --- hard Higgs partons, with a connector to the jet each matched to
    parton_eta = np.asarray(event["parton_eta"], dtype=float)
    parton_phi = np.asarray(event["parton_phi"], dtype=float)
    parton_pdg = np.asarray(event["parton_pdg"], dtype=int)
    parton_jet = np.asarray(event["parton_to_jet"], dtype=int)
    if len(parton_eta):
        ax.scatter(parton_eta, parton_phi, s=200, marker=(5, 2),
                   color="0.05", linewidths=1.1, zorder=7)
    for k in range(len(parton_eta)):
        ax.annotate(pdg_name(parton_pdg[k]), (parton_eta[k], parton_phi[k]),
                    textcoords="offset points", xytext=(6, -11), fontsize=7,
                    color="0.05", zorder=7)
        j = int(parton_jet[k])
        if 0 <= j < n_jets:
            target_phi = parton_phi[k] + wrap_pi(jet_phi[j] - parton_phi[k])
            ax.plot([parton_eta[k], jet_eta[j]], [parton_phi[k], target_phi],
                    color="0.35", lw=0.7, ls=(0, (1, 2)), zorder=2)

    ax.set_xlim(-2.9, 2.9)
    ax.set_ylim(-np.pi, np.pi)
    ax.set_xlabel(r"$\eta$")
    ax.set_ylabel(r"$\phi$")
    ax.set_yticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
    ax.set_yticklabels([r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    ax.grid(color="0.9", lw=0.5)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color("0.35")
        spine.set_linewidth(0.6)
    ax.set_title(
        f"{window['label']}  [{window['low']:.0f}, {window['high']:.0f}] GeV"
        f"   -   event {page_index + 1},   "
        rf"$m_H^{{\rm gen}} = {event['mH_gen']:.2f}$ GeV",
        fontsize=11, pad=8,
    )

    # --- legends and the annotation box, stacked in the right-hand column
    pid_handles = [
        Line2D([], [], marker=m, color="none", markerfacecolor=c,
               markeredgecolor=MARKER_EDGE, markeredgewidth=0.5, markersize=8,
               label=name)
        for name, c, m in groups
    ]
    leg_pid = side.legend(handles=pid_handles, loc="upper left",
                          bbox_to_anchor=(0.0, 1.0), frameon=True, fontsize=8,
                          labelspacing=0.35, borderpad=0.5, handletextpad=0.5,
                          title="Particle type (colour)", title_fontsize=8.5)
    leg_pid._legend_box.align = "left"
    side.add_artist(leg_pid)

    # Two columns, so the 40 GeV reference marker does not push the box tall
    # enough to collide with the annotation text below it.
    pt_handles = [
        Line2D([], [], marker="o", color="none", markerfacecolor="0.75",
               markeredgecolor=MARKER_EDGE, markeredgewidth=0.5,
               markersize=np.sqrt(marker_area(v)), label=f"{v:g}")
        for v in PT_LEGEND_VALUES
    ]
    leg_pt = side.legend(handles=pt_handles, loc="upper left",
                         bbox_to_anchor=(0.0, 0.795), frameon=True, fontsize=8,
                         ncol=2, columnspacing=1.6, labelspacing=1.1,
                         borderpad=0.6, handletextpad=0.7,
                         title=r"$p_T$ [GeV]  (marker area $\propto p_T$)",
                         title_fontsize=8.5)
    leg_pt._legend_box.align = "left"
    side.add_artist(leg_pt)

    side.text(0.0, 0.615, annotation_text(event, meta, window, page_index),
              transform=side.transAxes, fontsize=5.9, family="monospace",
              va="top", ha="left", linespacing=1.25,
              bbox=dict(boxstyle="round,pad=0.35", fc="#f7f9fc", ec="#9cc0ea", lw=0.7))
    return fig


PDG_NAMES = {1: "d", 2: "u", 3: "s", 4: "c", 5: "b", 6: "t", 21: "g",
             11: "e", 13: "mu", 15: "tau", 22: "gamma"}


def pdg_name(pdg):
    a = abs(int(pdg))
    name = PDG_NAMES.get(a, str(a))
    if a in (1, 2, 3, 4, 5, 6) and int(pdg) < 0:
        return name + "bar"
    return name


def annotation_text(event, meta, window, page_index):
    process = meta["process"]
    lines = [
        f"{process}",
        f"{HUMAN_READABLE_PROCESS_NAMES.get(process, '')}",
        "gen jets from stable gen particles (|eta|<2.56, no nu)",
        f"{meta['jet_algorithm']} clustering, N={meta['n_jets']};"
        f" H jets by dR<{0.3:.2f} parton match",
        f"window {window['label']} [{window['low']:.0f}, {window['high']:.0f}] GeV"
        f"   event {page_index + 1}",
        "",
        f"mH (Higgs gen jets)     = {event['mH_gen']:8.2f} GeV   <- h_mH_gen",
        f"mH (all {meta['n_jets']} gen jets)     = {event['mH_gen_all_jets']:8.2f} GeV",
        f"mH (visible truth part) = {event['mH_visible_truth']:8.2f} GeV",
        f"mH (hard partons)       = {event['mH_hard_partons']:8.2f} GeV",
        "",
        "jet   E[GeV]    pT     eta     phi  nconst   H",
    ]
    assign = np.asarray(event["part_jet_index"], dtype=int)
    n_jets = min(meta["n_jets"], len(event["jet_eta"]))
    higgs_jets = [int(j) for j in event["parton_to_jet"] if int(j) >= 0]
    for jet in range(n_jets):
        lines.append(
            f" j{jet}  {event['jet_energy'][jet]:7.1f} {event['jet_pt'][jet]:6.1f} "
            f"{event['jet_eta'][jet]:+7.2f} {event['jet_phi'][jet]:+7.2f} "
            f"{int(np.sum(assign == jet)):6d}   {'H' if jet in higgs_jets else '-'}"
        )
    lines += ["", "H parton -> gen jet   (dR)"]
    parton_eta = np.asarray(event["parton_eta"], dtype=float)
    parton_phi = np.asarray(event["parton_phi"], dtype=float)
    for k in range(len(parton_eta)):
        j = int(event["parton_to_jet"][k])
        if 0 <= j < n_jets:
            d_eta = parton_eta[k] - event["jet_eta"][j]
            d_phi = wrap_pi(parton_phi[k] - event["jet_phi"][j])
            lines.append(f"  {pdg_name(event['parton_pdg'][k]):>6s} -> j{j}  "
                         f"({np.hypot(d_eta, d_phi):.3f})")
        else:
            lines.append(f"  {pdg_name(event['parton_pdg'][k]):>6s} -> UNMATCHED")

    e_vis = float(event["E_vis_clustered"])
    e_nu = float(event["E_nu_total"])
    lost = 240.0 - e_vis - e_nu
    n_lost = int(np.sum(assign < 0))
    lines += [
        "",
        f"E_vis (clustered)   = {e_vis:7.1f} GeV",
        f"E_nu  (all)         = {e_nu:7.1f} GeV  from H: {float(event['E_nu_from_H']):.1f}",
        f"E_lost (|eta|>2.56) = {lost:7.1f} GeV  (240 - vis - nu)",
        f"particles clustered = {int(np.sum(assign >= 0)):5d}   unassigned = {n_lost}",
        f"gen jets found      = {int(event['n_gen_jets_total']):5d}",
    ]
    return "\n".join(lines)


def draw_cover(payload):
    """Overview page: which windows were selected and how much they hold."""
    fig = plt.figure(figsize=(11.7, 8.3))
    ax = fig.add_subplot(111)
    ax.axis("off")
    process = payload["process"]
    text = [
        "Event displays for the mH 'Physics' curve",
        "",
        f"process   {process}   {HUMAN_READABLE_PROCESS_NAMES.get(process, '')}",
        f"quantity  {payload['hist']} - Higgs mass from GEN jets built out of stable",
        "          gen particles, i.e. a perfect detector. Its width is therefore a",
        "          pure jet-definition / jet-assignment effect, which is what these",
        "          pages are meant to explain.",
        f"jets      {payload['jet_algorithm']}, exclusive to N={payload['n_jets']};"
        f" {payload['n_higgs_jets']} Higgs partons matched by dR",
        f"selection fully-matched-jets filter: "
        f"{'ON (as in the published figures)' if payload['matched_filter'] else 'OFF'}",
        f"read      {payload['n_input_events']} events",
        "",
        "windows",
    ]
    for window in payload["windows"]:
        text.append(
            f"  {window['label']:10s} [{window['low']:6.1f}, {window['high']:6.1f}] GeV"
            f"   {window['n_candidates']:6d} candidates"
            f" ({window.get('yield_percent', 0.0):5.2f}% of events read)"
            f"   {len(window['events']):3d} drawn"
        )
    text += [
        "",
        "reading a page",
        "  colour  = particle type      area = pT (linear, shared across pages)",
        f"  hull    = the {HULL_CORE_FRACTION:.0%} pT core of one gen jet; solid + [H] badge =",
        "            matched to a Higgs parton, dashed = not in the Higgs candidate.",
        "            Durham has no radius and assigns every particle to some jet,",
        "            so the full constituent set would sprawl across the plane;",
        "            the hull is drawn around the core to stay readable.",
        "  edge    = which jet a particle belongs to (all constituents, not just core)",
        "  star    = hard Higgs parton, dotted line to the jet it matched",
        "  neutrinos and |eta| > 2.56 particles are absent by construction;",
        "  their energy is given in each page's annotation box",
    ]
    ax.text(0.02, 0.98, "\n".join(text), transform=ax.transAxes, fontsize=10,
            family="monospace", va="top", ha="left")
    return fig


def draw_divider(window):
    fig = plt.figure(figsize=(11.7, 8.3))
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.text(0.5, 0.5,
            f"{window['label']}\n[{window['low']:.0f}, {window['high']:.0f}] GeV\n"
            f"{len(window['events'])} events",
            transform=ax.transAxes, fontsize=26, family="monospace",
            va="center", ha="center")
    return fig


def draw_payload(payload, output, pid_groups="five", cover=True):
    groups = PID_GROUPS_FIVE if pid_groups == "five" else PID_GROUPS_FOUR
    meta = {k: payload[k] for k in
            ("process", "jet_algorithm", "n_jets", "n_higgs_jets")}
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    n_pages = 0
    with PdfPages(output) as pdf:
        if cover:
            fig = draw_cover(payload)
            pdf.savefig(fig)
            plt.close(fig)
            n_pages += 1
        for window in payload["windows"]:
            if not window["events"]:
                continue
            fig = draw_divider(window)
            pdf.savefig(fig)
            plt.close(fig)
            n_pages += 1
            for i, event in enumerate(
                sorted(window["events"], key=lambda e: float(e["mH_gen"]))
            ):
                fig = draw_event(event, meta, window, i, groups)
                pdf.savefig(fig)
                plt.close(fig)
                n_pages += 1
    print(f"Wrote {n_pages} page(s) to {output}")
    return output


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--payload", required=True)
    p.add_argument("--output", default=None)
    p.add_argument("--pid-groups", choices=["five", "four"], default="five")
    p.add_argument("--no-cover-page", action="store_true")
    args = p.parse_args(argv)

    with open(args.payload, "rb") as handle:
        payload = pickle.load(handle)
    output = args.output or os.path.splitext(args.payload)[0] + ".pdf"
    draw_payload(payload, output, args.pid_groups, not args.no_cover_page)
    return 0


if __name__ == "__main__":
    sys.exit(main())
