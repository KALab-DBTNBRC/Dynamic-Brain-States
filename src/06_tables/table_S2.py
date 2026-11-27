#!/usr/bin/env python3
"""
Supplementary Table S2 - within-pair differences in brain-state temporal metrics, MZ vs DZ pairs.

Input : inputs/2026-11-13-master-sheet.xlsx, sheet "Master"
        (columns Subject, TwinPairID, Zygosity, FO_State1-4, DT_State1-4, AR_State1-4)
Method: for every twin pair and each of the 12 state x metric values,
            within-pair difference = |twin 1 - twin 2|
        MZ pairs vs DZ pairs: two-sided Mann-Whitney U test (pairs are the units).
        Effect sizes: mean difference DZ - MZ and rank-biserial r = P(DZ > MZ) - P(MZ > DZ),
        each with a 95% percentile bootstrap CI (10,000 resamples of pairs, stratified by
        zygosity, seed 20261112). No correction for multiple testing.
Output: outputs/SuppTableX.docx, outputs/SuppTableX.csv, outputs/SuppTableX_per_pair.csv
Run   : python 2026-11-13-SuppTableX.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt
from docx.table import _Cell
from scipy import stats

# ============================================================================ SETTINGS
MASTER = Path("inputs/2026-11-13-master-sheet.xlsx")
OUT_DIR = Path("outputs")
STATES = {1: "State 1: Visual Network-A", 2: "State 2: DMN",
          3: "State 3: Visual Network-B", 4: "State 4: Hippocampal–Temporoparietal"}
METRICS = [("FO", "Fractional occupancy", "", 4),        # code, name, unit, decimals shown
           ("DT", "Dwell time", "s", 3),
           ("AR", "Appearance rate", "per min", 3)]
# Alcohol and odor variables (group, column, label, decimals); pairs with a missing value are left out
EXTRA = [("Alcohol consumption", "SSAGA_Alc_12_Drinks_Per_Day", "Drinks per drinking occasion", 2),
         ("Alcohol consumption", "Total_Drinks_7days", "Drinks per week", 2),
         ("Odor identification", "Odor_AgeAdj", "Age-adjusted score", 2),
         ("Odor identification", "Odor_Unadj", "Unadjusted score", 2)]
N_BOOT, SEED = 10_000, 20261112
MINUS = "−"


# ============================================================================ analysis
def load_master():
    d = pd.read_excel(MASTER, sheet_name="Master")
    d["Zygosity"] = d["Zygosity"].str.strip().str.upper()
    assert set(d.Zygosity) == {"MZ", "DZ"} and (d.TwinPairID.value_counts() == 2).all()
    return d


def within_pair_differences(d):
    """One row per pair x state: |twin 1 - twin 2| for FO, DT and AR."""
    rows = []
    for (pair, zyg), g in d.groupby(["TwinPairID", "Zygosity"]):
        g = g.sort_values("Subject")
        for s in STATES:
            rec = {"TwinPairID": pair, "Zygosity": zyg, "State": s}
            for m, *_ in METRICS:
                rec[m] = abs(g[f"{m}_State{s}"].iloc[0] - g[f"{m}_State{s}"].iloc[1])
            rows.append(rec)
    return pd.DataFrame(rows)


def extra_differences(d):
    """One row per complete pair x alcohol/odor variable: |twin 1 - twin 2|."""
    rows = []
    for (pair, zyg), g in d.groupby(["TwinPairID", "Zygosity"]):
        g = g.sort_values("Subject")
        for grp, col, lab, _ in EXTRA:
            a, b = g[col].iloc[0], g[col].iloc[1]
            if pd.notna(a) and pd.notna(b):
                rows.append({"TwinPairID": pair, "Zygosity": zyg, "Variable": col, "value": abs(a - b)})
    return pd.DataFrame(rows)


def rank_biserial(dz, mz):
    """r = P(DZ > MZ) - P(MZ > DZ) over all DZ x MZ pair combinations."""
    return np.sign(dz[..., :, None] - mz[..., None, :]).mean(axis=(-2, -1))


def compare(wp, ex):
    """Brain-state rows first (random-number stream unchanged), then the alcohol and odor rows."""
    rng = np.random.default_rng(SEED)
    rows = []
    jobs = [(sname, s, mname + (f" ({unit})" if unit else ""), nd,
             wp.loc[(wp.State == s) & (wp.Zygosity == "MZ"), m].to_numpy(),
             wp.loc[(wp.State == s) & (wp.Zygosity == "DZ"), m].to_numpy())
            for s, sname in STATES.items() for m, mname, unit, nd in METRICS]
    jobs += [(grp, 5 + i // 2, lab, nd,
              ex.loc[(ex.Variable == col) & (ex.Zygosity == "MZ"), "value"].to_numpy(),
              ex.loc[(ex.Variable == col) & (ex.Zygosity == "DZ"), "value"].to_numpy())
             for i, (grp, col, lab, nd) in enumerate(EXTRA)]
    for sname, s, metric, nd, mz, dz in jobs:
        p = stats.mannwhitneyu(dz, mz, alternative="two-sided").pvalue
        bmz = mz[rng.integers(0, len(mz), (N_BOOT, len(mz)))]       # bootstrap of pairs
        bdz = dz[rng.integers(0, len(dz), (N_BOOT, len(dz)))]
        bdiff = bdz.mean(1) - bmz.mean(1)
        br = np.concatenate([rank_biserial(bdz[i:i + 1000], bmz[i:i + 1000])
                             for i in range(0, N_BOOT, 1000)])
        rows.append({"State": sname, "state": s, "Metric": metric,
                     "decimals": nd, "MZ pairs": len(mz), "DZ pairs": len(dz),
                     "MZ mean": mz.mean(), "MZ SD": mz.std(ddof=1),
                     "DZ mean": dz.mean(), "DZ SD": dz.std(ddof=1),
                     "Mean diff (DZ-MZ)": dz.mean() - mz.mean(),
                     "Mean diff CI low": np.percentile(bdiff, 2.5),
                     "Mean diff CI high": np.percentile(bdiff, 97.5),
                     "Rank-biserial r": rank_biserial(dz, mz),
                     "r CI low": np.percentile(br, 2.5), "r CI high": np.percentile(br, 97.5),
                     "p": p})
    return pd.DataFrame(rows)


# ============================================================================ Word table
def f(x, nd):
    return f"{x:.{nd}f}".replace("-", MINUS)


def set_borders(cell, top=0, bottom=0, grey_top=False):
    tcPr = cell._tc.get_or_add_tcPr()
    b = tcPr.find(qn("w:tcBorders"))
    if b is None:
        b = OxmlElement("w:tcBorders"); tcPr.append(b)
    for edge, sz in (("top", top), ("left", 0), ("bottom", bottom), ("right", 0)):
        el = OxmlElement(f"w:{edge}")
        if edge == "top" and grey_top:
            for k, v in (("val", "single"), ("sz", "4"), ("color", "A6A6A6"), ("space", "0")):
                el.set(qn(f"w:{k}"), v)
        else:
            el.set(qn("w:val"), "single" if sz else "nil")
            if sz:
                el.set(qn("w:sz"), str(sz)); el.set(qn("w:color"), "000000"); el.set(qn("w:space"), "0")
        b.append(el)


def write_docx(res, path):
    doc = Document()
    sec = doc.sections[0]
    sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = sec.page_height, sec.page_width
    for side in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(sec, side, Cm(1.8))
    st = doc.styles["Normal"]; st.font.name = "Arial"; st.font.size = Pt(9)
    st.paragraph_format.space_after = Pt(0)
    t = doc.add_paragraph()
    t.add_run("Supplementary Table X. ").bold = True
    t.add_run("Within-pair differences in brain-state temporal metrics in monozygotic (MZ) and "
              "dizygotic (DZ) twin pairs.")
    t.paragraph_format.space_after = Pt(6)
    nmz, ndz = int(res["MZ pairs"].iloc[0]), int(res["DZ pairs"].iloc[0])
    heads = ["Brain state", "Metric", f"MZ pairs\n(n = {nmz})", f"DZ pairs\n(n = {ndz})",
             f"Mean difference\nDZ {MINUS} MZ [95% CI]", "Rank-biserial r\n[95% CI]", "p"]
    table = doc.add_table(rows=1, cols=len(heads))
    for c, h in zip(table.rows[0].cells, heads):
        c.text = h
    group_rows, any_sig = set(), False
    for s, g in res.groupby("state", sort=True):
        if s != 1:
            group_rows.add(len(table.rows))
        first = None
        for k, (_, r) in enumerate(g.iterrows()):
            nd = int(r.decimals)
            cells = table.add_row().cells
            if k == 0:
                first = cells[0]; first.text = r.State
            else:
                first = first.merge(cells[0])
            cells[1].text = r.Metric
            cells[2].text = f"{f(r['MZ mean'], nd)} ± {f(r['MZ SD'], nd)}"
            cells[3].text = f"{f(r['DZ mean'], nd)} ± {f(r['DZ SD'], nd)}"
            cells[4].text = (f"{f(r['Mean diff (DZ-MZ)'], nd)} "
                             f"[{f(r['Mean diff CI low'], nd)}, {f(r['Mean diff CI high'], nd)}]")
            cells[5].text = f"{f(r['Rank-biserial r'], 2)} [{f(r['r CI low'], 2)}, {f(r['r CI high'], 2)}]"
            cells[6].text = ("< 0.001" if r.p < 0.001 else f"{r.p:.3f}") + ("*" if r.p < 0.05 else "")
            any_sig |= r.p < 0.05
    for row in table.rows[1:]:                       # merged state cell: keep the label once
        c = row.cells[0]
        txt = c.text.split("\n")[0]
        for p in c.paragraphs[1:]:
            p._element.getparent().remove(p._element)
        c.paragraphs[0].text = txt

    widths = [4.4, 4.1, 2.8, 2.8, 4.4, 3.4, 1.8]
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False
    for gc, w in zip(table._tbl.tblGrid.findall(qn("w:gridCol")), widths):
        gc.set(qn("w:w"), str(int(Cm(w).twips)))
    last = len(table.rows) - 1
    for ri, row in enumerate(table.rows):
        col = 0
        for tc in row._tr.tc_lst:
            cell, span, ci = _Cell(tc, table), tc.grid_span, col
            col += span
            cell.width = Cm(sum(widths[ci:ci + span]))
            set_borders(cell, top=12 if ri == 0 else 0, bottom=8 if ri == 0 else (12 if ri == last else 0),
                        grey_top=(ri in group_rows and ci > 0))
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            for p in cell.paragraphs:
                p.alignment = WD_ALIGN_PARAGRAPH.LEFT if ci < 2 else WD_ALIGN_PARAGRAPH.CENTER
                p.paragraph_format.space_before = Pt(1.5); p.paragraph_format.space_after = Pt(1.5)
                for run in p.runs:
                    run.font.size = Pt(8.5); run.bold = ri == 0

    note = doc.add_paragraph()
    note.paragraph_format.space_before = Pt(4)
    note.add_run("Note. ").italic = True
    note.add_run(
        f"For each twin pair and each state–metric combination, the within-pair difference was the "
        f"absolute difference between co-twins (|twin 1 {MINUS} twin 2|); smaller values indicate greater "
        "within-pair similarity. MZ and DZ columns show mean ± SD of within-pair differences across pairs. "
        f"Mean difference, DZ minus MZ (positive values: MZ co-twins more similar). Rank-biserial r = "
        f"P(DZ > MZ) {MINUS} P(MZ > DZ) across all MZ–DZ pair combinations (positive values: MZ co-twins more "
        "similar). 95% CIs: percentile bootstrap, 10,000 resamples of twin pairs stratified by zygosity. "
        "p, two-sided Mann–Whitney U test (MZ vs DZ pairs). FO, fractional occupancy (proportion); "
        "DT, dwell time (s); AR, appearance rate (visits per min)." + (" *p < 0.05." if any_sig else ""))
    doc.save(path)


# ============================================================================ main
def main():
    OUT_DIR.mkdir(exist_ok=True)
    master = load_master()
    wp = within_pair_differences(master)
    ex = extra_differences(master)
    res = compare(wp, ex)
    ex.to_csv(OUT_DIR / "SuppTableX_per_pair_alcohol_odor.csv", index=False)
    wp.to_csv(OUT_DIR / "SuppTableX_per_pair.csv", index=False)
    res.drop(columns=["state", "decimals"]).to_csv(OUT_DIR / "SuppTableX.csv", index=False)
    write_docx(res, OUT_DIR / "SuppTableX.docx")
    cols = ["State", "Metric", "MZ mean", "DZ mean", "Mean diff (DZ-MZ)", "Rank-biserial r", "p"]
    print(res[cols].round(4).to_string(index=False))
    print("Saved", OUT_DIR / "SuppTableX.docx")


if __name__ == "__main__":
    main()
