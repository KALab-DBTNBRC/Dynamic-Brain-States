#!/usr/bin/env python3
"""
Table 1 - demographic characteristics of participants (MZ vs DZ twins).

Input : inputs/2026-11-13-master-sheet.xlsx, sheet "Master"
        (Zygosity, Gender [0 = female], Age_in_Yrs, SSAGA_Alc_Age_1st_Use, Odor_AgeAdj,
         SSAGA_TB_Still_Smoking [1 = yes], Race, Ethnicity)
Tests : age and olfaction score - Mann-Whitney U (two-sided)
        age of onset of alcohol use - t-test (participants with a value)
        sex and smoking - chi-square with Yates' continuity correction
        race and ethnicity - chi-square
Output: outputs/Table1.docx, outputs/Table1.csv
Run   : python 2026-11-13-Table1.py
"""
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.shared import Cm, Pt
from scipy import stats

# ============================================================================ SETTINGS
MASTER = Path(r"C:\Users\Tanya\Desktop\Figures and Tables\Excels\Master-sheet.xlsx")
OUT_DIR = Path(r"C:\Users\Tanya\Desktop\Figures and Tables")
RACE = [("White", "White"), ("Black or African Am.", "Black"),
        ("Asian/Nat. Hawaiian/Othr Pacific Is.", "Asian"), ("Unknown or Not Reported", "Unknown")]
ETHNICITY = [("Hispanic/Latino", "Hispanic"), ("Not Hispanic/Latino", "Non-Hispanic"),
             ("Unknown or Not Reported", "Unknown")]


# ============================================================================ analysis
def make_table(d):
    g = {z: d[d.Zygosity == z] for z in ("MZ", "DZ")}
    n = {z: len(g[z]) for z in g}
    pct = lambda k, z: f"{k} ({round(100 * k / n[z], 1):.1f}%)"
    p_txt = lambda p: "p < 0.001" if p < 0.001 else f"p = {p:.3f}"
    rows = []

    def continuous(label, col, test):
        a, b = g["MZ"][col].dropna(), g["DZ"][col].dropna()
        p = (stats.mannwhitneyu(a, b, alternative="two-sided").pvalue if test == "mwu"
             else stats.ttest_ind(a, b).pvalue)
        rows.append([label, f"{a.mean():.3f} ± {a.std():.3f}", f"{b.mean():.3f} ± {b.std():.3f}",
                     f"P = {p:.3f}"])

    def chi2(counts, yates):
        c = stats.chi2_contingency(counts, correction=yates)
        return f"χ2 = {c.statistic:.3f}; df = {c.dof}; {p_txt(c.pvalue)}"

    continuous("Age (years) (mean ± SD)", "Age_in_Yrs", "mwu")
    fem = {z: int((g[z].Gender == 0).sum()) for z in g}
    rows.append(["Sex (Females) (n, %)", pct(fem["MZ"], "MZ"), pct(fem["DZ"], "DZ"),
                 chi2([[fem["MZ"], n["MZ"] - fem["MZ"]], [fem["DZ"], n["DZ"] - fem["DZ"]]], True)])
    continuous("Age of onset (Alcohol use) (mean ± SD)ᵃ", "SSAGA_Alc_Age_1st_Use", "t")
    continuous("Olfaction Score (mean ± SD)ᵇ", "Odor_AgeAdj", "mwu")
    smk = {z: int((g[z].SSAGA_TB_Still_Smoking == 1).sum()) for z in g}
    rows.append(["Current Smoking Status, n (%)", "", "", ""])
    rows.append(["Yes, n (%)", pct(smk["MZ"], "MZ"), pct(smk["DZ"], "DZ"),
                 chi2([[smk["MZ"], n["MZ"] - smk["MZ"]], [smk["DZ"], n["DZ"] - smk["DZ"]]], True)])
    rows.append(["No, n (%)", pct(n["MZ"] - smk["MZ"], "MZ"), pct(n["DZ"] - smk["DZ"], "DZ"), ""])
    for title, col, cats in (("Race, n (%)", "Race", RACE), ("Ethnicity, n (%)", "Ethnicity", ETHNICITY)):
        counts = [[int((g[z][col] == raw).sum()) for z in ("MZ", "DZ")] for raw, _ in cats]
        assert sum(map(sum, counts)) == len(d), f"{col}: category not listed in SETTINGS"
        rows.append([title, "", "", ""])
        for i, ((_, lab), (k1, k2)) in enumerate(zip(cats, counts)):
            rows.append([lab, pct(k1, "MZ"), pct(k2, "DZ"), chi2(counts, False) if i == 0 else ""])
    cols = ["Characteristics", f"Monozygotic Twins (n = {n['MZ']})", f"Dizygotic Twins (n = {n['DZ']})",
            "Test Statistics"]
    n_onset = {z: int(g[z].SSAGA_Alc_Age_1st_Use.notna().sum()) for z in g}
    return pd.DataFrame(rows, columns=cols), n, n_onset


# ============================================================================ Word table
def write_docx(df, n, n_onset, path):
    doc = Document()
    for side in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(doc.sections[0], side, Cm(2.0))
    st = doc.styles["Normal"]; st.font.name = "Arial"; st.font.size = Pt(9.5)
    st.paragraph_format.space_after = Pt(0)
    t = doc.add_paragraph()
    t.add_run("Table 1 | ").bold = True
    t.add_run(f"Demographic characteristics of participants (n = {n['MZ'] + n['DZ']})")
    t.paragraph_format.space_after = Pt(6)
    table = doc.add_table(rows=1, cols=4)
    table.style = "Table Grid"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False
    widths = [5.8, 3.4, 3.4, 4.6]
    for gc, w in zip(table._tbl.tblGrid.findall(qn("w:gridCol")), widths):
        gc.set(qn("w:w"), str(int(Cm(w).twips)))
    header = list(df.columns)
    header[3] = ("Test Statistics\nt-test or Mann–Whitney for continuous / "
                 "Chi-square for categorical variables")
    for c, h in zip(table.rows[0].cells, header):
        c.text = h
    for _, r in df.iterrows():
        for c, v in zip(table.add_row().cells, r):
            c.text = str(v)
    for ri, row in enumerate(table.rows):
        group = ri > 0 and row.cells[1].text == "" and row.cells[3].text == ""
        for ci, cell in enumerate(row.cells):
            cell.width = Cm(widths[ci])
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            for p in cell.paragraphs:
                p.alignment = WD_ALIGN_PARAGRAPH.LEFT if ci == 0 else WD_ALIGN_PARAGRAPH.CENTER
                p.paragraph_format.space_before = Pt(1.5); p.paragraph_format.space_after = Pt(1.5)
                for run in p.runs:
                    run.font.size = Pt(9.5); run.bold = ri == 0 or group
    note = doc.add_paragraph()
    note.paragraph_format.space_before = Pt(4)
    note.add_run(
        f"ᵃHCP SSAGA ordinal code for age at first alcohol use (not years); participants with a value: "
        f"MZ n = {n_onset['MZ']}, DZ n = {n_onset['DZ']}. ᵇNIH Toolbox Odor Identification Test, "
        "age-adjusted score. Age and olfaction score: Mann–Whitney U test; age of onset: t-test; sex and "
        "smoking: χ² with Yates' correction; race and ethnicity: χ².")
    doc.save(path)


# ============================================================================ main
def main():
    OUT_DIR.mkdir(exist_ok=True)
    d = pd.read_excel(MASTER, sheet_name="Master")
    d["Zygosity"] = d["Zygosity"].str.strip().str.upper()
    assert d.Subject.is_unique and len(d) == 240
    df, n, n_onset = make_table(d)
    df.to_csv(OUT_DIR / "Table1.csv", index=False)
    write_docx(df, n, n_onset, OUT_DIR / "Table1.docx")
    print(df.to_string(index=False))
    print("Saved", OUT_DIR / "Table1.docx")


if __name__ == "__main__":
    main()
