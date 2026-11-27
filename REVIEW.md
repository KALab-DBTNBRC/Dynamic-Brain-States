# Final review — unresolved items

## Confirmed
- Latest manuscript and supplementary titles match. Main figure legends are embedded in the manuscript. The centroid matching direction in section 3.2 is corrected. S5 age-at-first-use rows now include degrees of freedom.
- Figure 1 corrected source uses DMN for CAP 2 and retains atlas DMN-A/B/C and LIM-B labels.
- All packaged Python sources were checked with AST parsing; this does not test runtime behavior. R/MATLAB runtimes and required participant inputs were not available.

## Before release
1. Step 03 skips header rows and then treats the first data row as a header, losing an IC and potentially shifting alignment. Its own comment warns about this. S1 CSV contains 298 rows (73 Signal), S3 contains 99 rows (38 Signal); manuscript reports 64 and 43 retained components. Visual selection may explain some differences, but verify original indices and provide explicit retained lists. Step 04 hard-codes a 64-column selection rather than reading classification results. Step 01 uses signed sums and step 02 returns signed maximum absolute motion correlation, whereas step 03 thresholds positive values only; verify intended feature definitions. No scientific changes made.
2. Supplementary File S2 contains 17 x 4 correlations (four targets), not 17 x 64 IC correlations. Manuscript section describing IC-level correspondence should point to the correct file or describe CAP-level correlations.
3. Figure 3 source explicitly reads ordinal mixed-effects results, while manuscript describes continuous-coded linear mixed models. The updated PNG is preserved, but the exact source/workbook used for the final PNG was not supplied. Verify panel C estimates and standardized weekly alcohol units. Both continuous and ordinal R scripts included.
4. Step 16 uses Pearson r and one-sided Mann–Whitney, not final ICC and two-sided tests; moved to legacy. Table S3 script is the more relevant final analysis. Step 09 uses 125 timepoints and a toolbox call rather than documented 200 prototypes; moved to legacy. Missing explicit subject-level 200-prototype generation, final representative k-means partition, and FO/DT/AR extraction scripts prevent end-to-end reproduction.
5. Step 08 originally imported nonexistent consensus_core; changed import to step07_consensus_core. Step 06 and 08 use different input MAT schemas (X_prototypes versus P). Step 08 documents the Python run used for final consensus results.
6. Head-motion table script originally imported itself under its supplied filename. Redirected to table_S3.py (ICC helper); runtime verification needed. ACE table generator lacks age-at-first-use in EXTRA_GROUPS, and older output captions/numbers remain in several scripts. The final Word tables are preserved, not regenerated.
7. S5 notes call 79 MZ + 38 DZ and 80 MZ + 39 DZ complete pairs despite 231 and 235 observations; clarify families with singletons. Age-at-first-use counts need reconciliation with original ACE output. Manuscript still states mean dwell time approximately 3.2 s; Table 2 means imply approximately 3.1 s.
8. Figure 1 and Table 1 contain hard-coded Windows paths; Figure 2 has Windows defaults. Other scripts require workbooks, NIfTI maps, atlas files and prototype matrices not in current upload. No path refactoring or model edits made.
9. The cover letter and manuscript remain private submission materials. Verify editor designation, declarations, and data permissions before submission/release.

## Changes made
Folder/file naming, source selection (latest individual files and corrected archive), State 2/4 display names in regression scripts, and two local import references only. No numerical or statistical method changes. See SOURCE_MANIFEST.txt.
