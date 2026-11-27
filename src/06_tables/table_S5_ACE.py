#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supplementary Table S4 - ACE twin-model variance components, fit statistics and model comparisons
for the brain-state dynamic metrics (FO, DT and AR, States 1-4).

Builds the table directly from the ACE results file (ACE_Master_Results.csv): for each metric and state,
the ACE, AE, CE and E models with a2, c2, e2 (95% CIs for ACE), -2LL, number of parameters, chi-square
fit against the saturated model, AIC, and likelihood-ratio tests against the ACE model. The best-fitting
model (lowest AIC) is marked. Alcohol and odor variables in the results file are not included.

HOW TO RUN
  1. Edit the two lines under "EDIT THESE TWO LINES" (full path to the ACE results CSV and an output folder).
  2. Press Run (Spyder, VS Code, IDLE, Jupyter) or type:  python "2026-11-14-Table_X3_ACE_models.py"

OUTPUT  2026-11-14-Table-X3-ACE-models.docx   journal-format Word table (landscape A4)
        2026-11-14-Table-X3-ACE-models.xlsx   full-precision numbers (sheet Table_X3) and the
                                              table text exactly as printed (sheet Table_X3_formatted)

This script is self-contained. Requires Python 3.10+ with pandas, openpyxl and python-docx.
It only reads and reformats the results file; it does not refit any model.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt, Twips

# =============================================================================
# EDIT THESE TWO LINES
# =============================================================================
ACE_RESULTS_FILE = r"inputs/ACE_Master_Results.csv"   # full path, e.g. r"C:\Users\me\Documents\twin\ACE_Master_Results.csv"
OUTPUT_FOLDER = r"outputs"                      # folder for the results (created if it does not exist)
# =============================================================================

FILE_PREFIX = "2026-11-14"
METRICS = [("FO", "Fractional occupancy (FO)"), ("DT", "Dwell time (DT)"), ("AR", "Appearance rate (AR)")]
STATES = ["State 1", "State 2", "State 3", "State 4"]
STATE_NAMES = {"State 1": "Visual Network-A", "State 2": "Default Mode Network", "State 3": "Visual Network-B",
               "State 4": "Hippocampal–Temporoparietal"}
# Alcohol and odor variables: (group heading, [(variable in the ACE file, row label), ...])
EXTRA_GROUPS = [("Alcohol consumption", [("SSAGA_Alc_12_Drinks_Per_Day", "Drinks per drinking occasion"),
                                         ("Total_Drinks_7days", "Drinks per week")]),
                ("Odor identification", [("Odor_AgeAdj", "Age-adjusted"),
                                         ("Odor_Unadj", "Unadjusted")])]
MODELS = [("ACE (Base)", "ACE"), ("AE (Drop C)", "AE"), ("CE (Drop A)", "CE"), ("E (Drop A & C)", "E")]
FIXED_ZERO = {"ACE": set(), "AE": {"C"}, "CE": {"A"}, "E": {"A", "C"}}
ALPHA = 0.05
FONT, TABLE_PT, NOTE_PT, TITLE_PT = "Arial", 7.5, 8, 10
MINUS = "\u2212"

TITLE = ("Supplementary Table X3. ACE twin-model variance components, fit statistics and model comparisons "
         "for the brain-state dynamic metrics, alcohol consumption and odor identification")
NOTES = [
    ("Supplementary Table X3 : ", "Univariate twin models fitted to 80 MZ and 40 DZ complete twin pairs for each "
     "brain-state metric and odor identification score; drinks per drinking occasion was available for 79 MZ and 38 DZ pairs "
     "(231 participants) and drinks per week for 80 MZ and 39 DZ pairs (235 participants). a\u00b2, c\u00b2 "
     "and e\u00b2 are the proportions of phenotypic variance attributed to additive genetic (A), shared environmental (C) "
     "and unique environmental (E, including measurement error) influences; 95% confidence intervals (CIs) are given "
     "for the full ACE model of the brain-state metrics (CIs were not estimated for the alcohol and odor variables). \u22122LL, minus twice the log-likelihood; k, number of free parameters; \u03c7\u00b2 (df) "
     "and p compare each model with the saturated model; AIC = \u22122LL + 2k (lower values indicate better fit). "
     "\u0394\u22122LL, \u0394df and p are likelihood-ratio tests against the ACE model; because a variance component "
     "fixed at zero lies on the boundary of the parameter space, these p-values are conservative. "
     "\u2014, parameter fixed at zero or not applicable; n.e., CI bound not estimable (estimate at or near zero). "
     "\u2020Best-fitting model (lowest AIC). Bold p-values are < 0.05 (uncorrected)."),
]


# ----------------------------------------------------------------------------- data
def load_results(path):
    path = Path(path)
    if not path.exists():
        sys.exit(f"ACE results file not found: {path}\nEdit ACE_RESULTS_FILE at the top of this script.")
    d = pd.read_csv(path, float_precision="round_trip")
    need = ["Variable", "Model", "Parameters", "minus2LL", "df", "chi_square", "chi_square_p", "AIC", "diffLL",
            "diffdf", "p_value", "prop_A", "prop_C", "prop_E", "A2_SE", "C2_SE", "E2_SE", "A2_CI_lower",
            "A2_CI_upper", "C2_CI_lower", "C2_CI_upper", "E2_CI_lower", "E2_CI_upper", "MZ_pairs", "DZ_pairs"]
    missing = [c for c in need if c not in d.columns]
    if missing:
        sys.exit(f"Columns missing from the ACE results file: {missing}")
    rows = []
    for m, _ in METRICS:
        for s in STATES:
            v = f"{m}_{s}"
            block = d[d.Variable == v]
            for model, short in MODELS:
                r = block[block.Model == model]
                if len(r) != 1:
                    sys.exit(f"Expected exactly one '{model}' row for {v}; found {len(r)}.")
                rows.append(r.iloc[0].to_dict() | {"Metric": m, "State": s, "Model_short": short})
    for gname, items in EXTRA_GROUPS:
        for v, lab in items:
            block = d[d.Variable == v]
            for model, short in MODELS:
                r = block[block.Model == model]
                if len(r) != 1:
                    sys.exit(f"Expected exactly one '{model}' row for {v}; found {len(r)}.")
                rows.append(r.iloc[0].to_dict() | {"Metric": gname, "State": lab, "Model_short": short})
    t = pd.DataFrame(rows)
    t["Best_by_AIC"] = t.groupby("Variable")["AIC"].transform("min").eq(t["AIC"])
    return t


# ----------------------------------------------------------------------------- formatting
def num(x, nd=2):
    if pd.isna(x):
        return "n.e."
    if round(x, nd) == 0:
        x = 0.0
    return f"{x:.{nd}f}".replace("-", MINUS)


def fmt_p(p):
    """Returns (text, superscript exponent, bold)."""
    if pd.isna(p):
        return "\u2014", "", False
    if p < 0.001:
        m, e = f"{p:.1e}".split("e")
        return f"{m} \u00d7 10", str(int(e)).replace("-", MINUS), True
    txt = f"{p:.3f}"
    if p < ALPHA <= float(txt):
        txt = f"{p:.4f}"
    return txt, "", p < ALPHA


def component(r, comp):
    """Estimate (with 95% CI for the ACE model) or a dash when the component is fixed at zero."""
    if comp in FIXED_ZERO[r["Model_short"]]:
        return "\u2014"
    est = num(r[f"prop_{comp}"])
    if r["Model_short"] == "ACE" and pd.notna(r[f"{comp}2_CI_upper"]):
        return f"{est} [{num(r[f'{comp}2_CI_lower'])}, {num(r[f'{comp}2_CI_upper'])}]"
    return est


def table_rows(t):
    """List of (group label or None, [cell values]) in print order."""
    out = []
    groups = [(m, label, STATES) for m, label in METRICS] + \
             [(g, g, [lab for _, lab in items]) for g, items in EXTRA_GROUPS]
    for m, label, states in groups:
        out.append((label, None))
        for s in states:
            for i, (_, r) in enumerate(t[(t.Metric == m) & (t.State == s)].iterrows()):
                ace = r["Model_short"] == "ACE"
                out.append((None, [
                    STATE_NAMES.get(s, s) if i == 0 else "",
                    r["Model_short"] + ("\u2020" if r["Best_by_AIC"] else ""),
                    component(r, "A"), component(r, "C"), component(r, "E"),
                    num(r["minus2LL"]), str(int(r["Parameters"])),
                    f"{num(r['chi_square'])} ({int(r['df'])})", fmt_p(r["chi_square_p"]),
                    num(r["AIC"]),
                    "\u2014" if ace else num(r["diffLL"]),
                    "\u2014" if ace else str(int(r["diffdf"])),
                    ("\u2014", "", False) if ace else fmt_p(r["p_value"]),
                ]))
    return out


# ----------------------------------------------------------------------------- Word
def _run(par, text, size, bold=False, italic=False, sup=False):
    r = par.add_run(text)
    r.font.name = FONT
    r._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    r.font.size = Pt(size)
    r.bold, r.italic, r.font.superscript = bold, italic, sup
    return r


def _borders(cell, top=None, bottom=None):
    tcPr = cell._tc.get_or_add_tcPr()
    old = tcPr.find(qn("w:tcBorders"))
    if old is not None:
        tcPr.remove(old)
    b = OxmlElement("w:tcBorders")
    for edge, w in (("top", top), ("left", None), ("bottom", bottom), ("right", None)):
        el = OxmlElement(f"w:{edge}")
        if w:
            el.set(qn("w:val"), "single"); el.set(qn("w:sz"), str(w)); el.set(qn("w:color"), "000000")
        else:
            el.set(qn("w:val"), "nil")
        b.append(el)
    tcPr.append(b)


def _margins(cell, top=25, bottom=25, left=50, right=50):
    tcPr = cell._tc.get_or_add_tcPr()
    mar = OxmlElement("w:tcMar")
    for k, v in (("top", top), ("left", left), ("bottom", bottom), ("right", right)):
        el = OxmlElement(f"w:{k}"); el.set(qn("w:w"), str(v)); el.set(qn("w:type"), "dxa"); mar.append(el)
    tcPr.append(mar)


def _write(cell, value, align="center", bold=False, italic=False, top=None, bottom=None, keep_next=False, pad_top=25):
    for extra in cell.paragraphs[1:]:
        extra._element.getparent().remove(extra._element)
    par = cell.paragraphs[0]
    par.alignment = WD_ALIGN_PARAGRAPH.LEFT if align == "left" else WD_ALIGN_PARAGRAPH.CENTER
    pf = par.paragraph_format
    pf.space_before = Pt(0); pf.space_after = Pt(0); pf.keep_with_next = keep_next
    if isinstance(value, tuple):                       # p-value: (text, exponent, bold)
        txt, sup, b = value
        _run(par, txt, TABLE_PT, b)
        if sup:
            _run(par, sup, TABLE_PT, b, sup=True)
    else:
        _run(par, value, TABLE_PT, bold, italic)
    cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
    _borders(cell, top, bottom); _margins(cell, top=pad_top)


def write_word(t, path):
    doc = Document()
    sec = doc.sections[0]
    sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = Twips(16838), Twips(11906)
    for side in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(sec, side, Twips(900))
    st = doc.styles["Normal"]; st.font.name = FONT; st.font.size = Pt(10)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), FONT)

    p = doc.add_paragraph(); p.paragraph_format.space_after = Pt(6)
    _run(p, TITLE.split(".", 1)[0] + ".", TITLE_PT, bold=True)
    _run(p, TITLE.split(".", 1)[1], TITLE_PT)

    widths = [800, 650, 1650, 1650, 1650, 1050, 400, 1150, 1050, 1050, 900, 550, 1050]
    head1 = [("", 1, True), ("Model", 1, True), ("Variance components", 3, False), ("Model fit", 5, False),
             ("Comparison with ACE", 3, False)]
    head2 = ["a\u00b2 [95% CI]", "c\u00b2 [95% CI]", "e\u00b2 [95% CI]", f"{MINUS}2LL", "k", "\u03c7\u00b2 (df)", "p",
             "AIC", f"\u0394{MINUS}2LL", "\u0394df", "p"]
    ncol = len(widths)
    tbl = doc.add_table(rows=2, cols=ncol)
    tbl.alignment = WD_TABLE_ALIGNMENT.CENTER
    tbl.autofit = False
    lay = OxmlElement("w:tblLayout"); lay.set(qn("w:type"), "fixed"); tbl._tbl.tblPr.append(lay)
    for i, w in enumerate(widths):
        tbl.columns[i].width = Twips(w)
        for rr in range(2):
            tbl.cell(rr, i).width = Twips(w)
    col, rowspan_cols = 0, set()
    for text, n, rowspan in head1:
        c = tbl.cell(0, col)
        if n > 1:
            c = c.merge(tbl.cell(0, col + n - 1))
        if rowspan:
            c = c.merge(tbl.cell(1, col)); rowspan_cols.add(col)
        _write(c, text, align="left" if col == 0 else "center", top=12, bottom=4)
        col += n
    it = iter(head2)
    for c in range(ncol):
        if c not in rowspan_cols:
            _write(tbl.cell(1, c), next(it), bottom=4)
    from docx.table import _Cell
    grid = 0
    for tc in tbl.rows[1]._tr.tc_lst:                 # rule under the row-spanning header cells
        if grid in rowspan_cols:
            _borders(_Cell(tc, tbl), bottom=4)
        grid += tc.grid_span
    for r in tbl.rows[:2]:
        trPr = r._tr.get_or_add_trPr(); h = OxmlElement("w:tblHeader"); trPr.append(h)

    rows = table_rows(t)
    block_start = None
    for k, (group, vals) in enumerate(rows):
        r = tbl.add_row()
        for c, w in enumerate(widths):
            r.cells[c].width = Twips(w)
        last = k == len(rows) - 1
        if group:
            g = r.cells[0].merge(r.cells[ncol - 1])
            _write(g, group, align="left", italic=True, keep_next=True, pad_top=70)
            continue
        first_of_block = vals[0] != ""
        block_pos = [i for i, (gg, vv) in enumerate(rows[:k + 1]) if gg is None and vv[0] != ""][-1]
        keep = (k - block_pos) < 3                     # keep the 4 model rows of a state together
        for c, v in enumerate(vals):
            _write(r.cells[c], v, align="left" if c < 2 else "center", bottom=12 if last else None,
                   keep_next=keep, pad_top=55 if first_of_block else 15)
        trPr = r._tr.get_or_add_trPr(); cs = OxmlElement("w:cantSplit"); trPr.append(cs)
        if first_of_block:
            block_start = r
        elif vals[0] == "" and block_start is not None:     # one label per state, merged over its model rows
            top = block_start.cells[0]
            merged = top.merge(r.cells[0])
            for extra in merged.paragraphs[1:]:
                extra._element.getparent().remove(extra._element)
    doc.add_paragraph().paragraph_format.space_after = Pt(2)
    for lead, body in NOTES:
        p = doc.add_paragraph(); p.paragraph_format.space_after = Pt(3)
        _run(p, lead, NOTE_PT, italic=True); _run(p, body, NOTE_PT)
    doc.save(path)


# ----------------------------------------------------------------------------- Excel
def write_excel(t, path):
    from openpyxl import Workbook
    from openpyxl.styles import Font, Alignment
    from openpyxl.utils import get_column_letter
    numeric = pd.DataFrame({
        "Metric": t.Metric, "State": t.State, "Model": t.Model, "Best by AIC": np.where(t.Best_by_AIC, "Yes", ""),
        "a2": t.prop_A, "a2 95% CI low": t.A2_CI_lower, "a2 95% CI high": t.A2_CI_upper, "a2 SE": t.A2_SE,
        "c2": t.prop_C, "c2 95% CI low": t.C2_CI_lower, "c2 95% CI high": t.C2_CI_upper, "c2 SE": t.C2_SE,
        "e2": t.prop_E, "e2 95% CI low": t.E2_CI_lower, "e2 95% CI high": t.E2_CI_upper, "e2 SE": t.E2_SE,
        "-2LL": t.minus2LL, "k (parameters)": t.Parameters, "chi2 vs saturated": t.chi_square, "df": t.df,
        "p (vs saturated)": t.chi_square_p, "AIC": t.AIC, "delta -2LL vs ACE": t.diffLL, "delta df": t.diffdf,
        "p (LRT vs ACE)": t.p_value, "MZ pairs": t.MZ_pairs, "DZ pairs": t.DZ_pairs})
    shown = []
    for group, vals in table_rows(t):
        if group:
            shown.append([group] + [""] * 12)
        else:
            shown.append([v if not isinstance(v, tuple) else (v[0] + (f"^{v[1]}" if v[1] else "")) for v in vals])
    formatted = pd.DataFrame(shown, columns=["State", "Model", "a2 [95% CI]", "c2 [95% CI]", "e2 [95% CI]", "-2LL", "k",
                                             "chi2 (df)", "p", "AIC", "delta -2LL", "delta df", "p (vs ACE)"])
    wb = Workbook(); wb.remove(wb.active)
    for name, df in (("Table_X3", numeric), ("Table_X3_formatted", formatted)):
        ws = wb.create_sheet(name)
        for j, col in enumerate(df.columns, 1):
            c = ws.cell(1, j, col); c.font = Font(name="Arial", size=10, bold=True)
            c.alignment = Alignment(wrap_text=True, vertical="top")
            ws.column_dimensions[get_column_letter(j)].width = 18 if name.endswith("formatted") else 13
        for i, row in enumerate(df.itertuples(index=False), 2):
            for j, v in enumerate(row, 1):
                if isinstance(v, (float, np.floating)):
                    v = None if not np.isfinite(v) else float(v)
                elif isinstance(v, np.integer):
                    v = int(v)
                ws.cell(i, j, v).font = Font(name="Arial", size=10)
        ws.freeze_panes = "A2"
    wb.save(path)
    return numeric, formatted


def main():
    t = load_results(ACE_RESULTS_FILE)
    out = Path(OUTPUT_FOLDER); out.mkdir(parents=True, exist_ok=True)
    docx_path = out / f"{FILE_PREFIX}-Table-X3-ACE-models.docx"
    xlsx_path = out / f"{FILE_PREFIX}-Table-X3-ACE-models.xlsx"
    write_word(t, docx_path)
    numeric, formatted = write_excel(t, xlsx_path)
    print(f"Table X3: {t.Variable.nunique()} metrics x {len(MODELS)} models")
    print("Saved:", str(docx_path), str(xlsx_path), sep="\n  ")
    return t, numeric, formatted


if __name__ == "__main__":
    main()
