#!/usr/bin/env python3
"""
Table 2 - temporal metrics of the four brain states in the full cohort, MZ twins and DZ twins.

Input : inputs/2026-11-13-master-sheet.xlsx, sheet "Master"
        (columns Subject, TwinPairID, Zygosity, FO_State1-4, DT_State1-4, AR_State1-4)
Output: outputs/Table2.docx  (Word table, landscape)
        outputs/Table2.csv   (same values)
Run   : python 2026-11-13-Table2.py
"""
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt
from docx.table import _Cell

# ============================================================================ SETTINGS
MASTER = Path("inputs/2026-11-13-master-sheet.xlsx")
OUT_DIR = Path("outputs")
STATES = {1: "State 1: Visual Network-A", 2: "State 2: DMN",
          3: "State 3: Visual Network-B", 4: "State 4: Hippocampal–Temporoparietal"}
METRICS = [("FO", "Fractional occupancy", "", 3),      # code, name, unit, decimals
           ("DT", "Dwell time", "s", 2),
           ("AR", "Appearance rate", "per min", 2)]
COHORTS = [("Full", "Full cohort"), ("MZ", "MZ twins"), ("DZ", "DZ twins")]


# ============================================================================ analysis
def load_master():
    d = pd.read_excel(MASTER, sheet_name="Master")
    d["Zygosity"] = d["Zygosity"].str.strip().str.upper()
    assert set(d.Zygosity) == {"MZ", "DZ"} and d.Subject.is_unique
    return d


def make_table(d):
    n = {"Full": len(d), "MZ": int((d.Zygosity == "MZ").sum()), "DZ": int((d.Zygosity == "DZ").sum())}
    rows = []
    for s, sname in STATES.items():
        row = {"State": sname}
        for m, _, _, nd in METRICS:
            col = f"{m}_State{s}"
            for key, cname in COHORTS:
                x = d[col] if key == "Full" else d.loc[d.Zygosity == key, col]
                row[f"{m} | {cname}"] = f"{x.mean():.{nd}f} ± {x.std(ddof=1):.{nd}f}"
        rows.append(row)
    return pd.DataFrame(rows), n


# ============================================================================ Word table
def set_borders(cell, top=0, bottom=0):
    tcPr = cell._tc.get_or_add_tcPr()
    b = tcPr.find(qn("w:tcBorders"))
    if b is None:
        b = OxmlElement("w:tcBorders"); tcPr.append(b)
    for edge, sz in (("top", top), ("left", 0), ("bottom", bottom), ("right", 0)):
        el = OxmlElement(f"w:{edge}")
        el.set(qn("w:val"), "single" if sz else "nil")
        if sz:
            el.set(qn("w:sz"), str(sz)); el.set(qn("w:color"), "000000"); el.set(qn("w:space"), "0")
        b.append(el)


def write_docx(wide, n, path):
    doc = Document()
    sec = doc.sections[0]
    sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = sec.page_height, sec.page_width
    for side in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(sec, side, Cm(1.8))
    st = doc.styles["Normal"]; st.font.name = "Arial"; st.font.size = Pt(9)
    st.paragraph_format.space_after = Pt(0)
    t = doc.add_paragraph()
    t.add_run("Table 2. ").bold = True
    t.add_run("Temporal metrics of the four brain states in the full cohort, MZ twins and DZ twins.")
    t.paragraph_format.space_after = Pt(6)

    table = doc.add_table(rows=2, cols=1 + 3 * len(METRICS))
    h0, h1 = table.rows[0].cells, table.rows[1].cells
    h0[0].merge(h1[0]).text = "Brain state"
    for i, (m, mname, unit, _) in enumerate(METRICS):
        h0[1 + 3 * i].merge(h0[3 + 3 * i]).text = f"{mname} ({m})" + (f", {unit}" if unit else "")
        for j, (k, cname) in enumerate(COHORTS):
            h1[1 + 3 * i + j].text = f"{cname}\n(n = {n[k]})"
    for _, r in wide.iterrows():
        cells = table.add_row().cells
        cells[0].text = r["State"]
        for i, (m, *_) in enumerate(METRICS):
            for j, (k, cname) in enumerate(COHORTS):
                cells[1 + 3 * i + j].text = r[f"{m} | {cname}"]

    widths = [4.2] + [2.23] * 9
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
            set_borders(cell, top=12 if ri == 0 else 0, bottom=8 if ri == 1 else (12 if ri == last else 0))
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            for p in cell.paragraphs:
                p.alignment = WD_ALIGN_PARAGRAPH.LEFT if ci == 0 else WD_ALIGN_PARAGRAPH.CENTER
                p.paragraph_format.space_before = Pt(1.5); p.paragraph_format.space_after = Pt(1.5)
                for run in p.runs:
                    run.font.size = Pt(8.5); run.bold = ri < 2

    note = doc.add_paragraph()
    note.paragraph_format.space_before = Pt(4)
    note.add_run("Note. ").italic = True
    note.add_run("Values are mean ± SD across participants. FO, fractional occupancy (proportion of "
                 "volumes assigned to the state); DT, mean dwell time (mean duration of uninterrupted "
                 "visits to the state); AR, appearance rate (number of visits to the state per minute). "
                 "States were identified once in the full cohort and the same four states were applied "
                 f"to every participant. Full cohort = {n['Full']} participants "
                 f"({n['MZ'] // 2} MZ and {n['DZ'] // 2} DZ twin pairs).")
    doc.save(path)


# ============================================================================ main
def main():
    OUT_DIR.mkdir(exist_ok=True)
    d = load_master()
    wide, n = make_table(d)
    wide.to_csv(OUT_DIR / "Table2.csv", index=False)
    write_docx(wide, n, OUT_DIR / "Table2.docx")
    print(f"n: full = {n['Full']}, MZ = {n['MZ']}, DZ = {n['DZ']}")
    print(wide.to_string(index=False))
    print("Saved", OUT_DIR / "Table2.docx")


if __name__ == "__main__":
    main()
