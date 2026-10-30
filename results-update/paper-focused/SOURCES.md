# Sources and scientific limits

All packaged observations are existing public outputs from [sanithuheen/ProteinMPNN-Ablang-antibody-design at bc3a1070a4e63489ef785763ae2193c4ffd5aa19](https://github.com/sanithuheen/ProteinMPNN-Ablang-antibody-design/tree/bc3a1070a4e63489ef785763ae2193c4ffd5aa19). [inputs/provenance.json](inputs/provenance.json) maps each selected local filename to the exact public URL, source revision and SHA-256. These are retained observations, not newly generated predictions. No private manuscript, messages, coordination records or PDB collection is included.

## Main cohort

Four BioPhi workbooks supply Heavy OASis Identity from the `Overview` sheet; the `VH` sheet reconstructs the scored domain. Four corresponding final-filled FASTAs establish the selected index linkage. Baseline sample minus one and unchanged ensemble seed map to zero-based structural design IDs. Equal index values across methods do **not** establish paired random draws. The main analysis therefore reports marginal summaries only. The original framework/CDR CSV also contains other loops and comparator arms; this reader selects baseline/ensemble H3 for Table 1 and all nine case/loop combinations for the supporting table, not comparator arms.

The current selected main inputs are no longer all unbound. The preceding audit linked all selected VH sequences to FASTAs and all selected full predicted sequences to those FASTAs, and replayed coordinate measurements. [Prior coordinate receipt](engineering/prior_coordinate_replay.json) records the measured scope and numeric comparison; it is prior evidence, not a run performed by this notebook. The same-author portable recalculation and tests here are computational checks, not independent scientific acceptance.

Heavy OASis Identity describes a whole numbered heavy domain at the reported relaxed threshold. Do not substitute pooled heavy/light Identity, OASis percentile, Jacob `promb` human-oas peptide-content values or clinical immunogenicity. No scoring database was run or authenticated here.

Structural RMSD is the implemented operational H3 backbone metric after a **per-domain framework fit of N/CA/C atoms to deposited experimental references**. It is not an original RFdiffusion-design-model distance. The motif/sequence-alignment mapping is not independently validated ANARCI/AHo/Chothia numbering. Predicted coordinates and confidence do not establish the true fold, binding or contact geometry. Historical generator settings, weights and execution identity remain unauthenticated. Separate domain fits remove relative VH/VL orientation from this endpoint and cannot establish whole-scFv or interface preservation. The public compute_RSMD.py:37–110,218–240,324–405 defines reference chains E/C, motif spans and domain framework N/CA/C fitting with CDRs excluded. Mean loop pLDDT below 70, not any low-confidence residue, defines the confidence flag. All draws remain primary; low confidence does not explain away every structural penalty.

## Separate ablation cohort

The public `Ablation_results/Oasis_Ablation_results/oasis_scores.csv` supplies the separate loop-group scores, joined by explicit case/condition/reported seed. The public [run_Ablation.py](https://github.com/sanithuheen/ProteinMPNN-Ablang-antibody-design/blob/bc3a1070a4e63489ef785763ae2193c4ffd5aa19/run_Ablation.py) defines a same-seed zero-logit control through the same decoder. The named `get_ensemble_sequences` helper/import resolution and exact executed revision/weights remain unresolved. This is a source-defined matched control and descriptive seed correspondence, not authenticated causal execution. It does not pair the main all-loop cohorts retroactively. No ablation structural outputs are used.

The earlier scoring-input audit compared the filled ablation sequences and retained workbook values; this portable reader recalculates the score CSV only. Designed-position normalization is a heuristic because overlapping 9-mer windows are not additive residue effects. Other CDRs means H1/H2 in 9NH7 and H1/H2 + L1/L2/L3 in 9NFU (run_Ablation.py:58–78). The 9NFU heavy score is a joint heavy/light redesign outcome, not an isolated heavy-loop intervention. Light-domain gain and the H3-only unchanged-light control are reported separately in the reader and supporting table; no cross-chain efficiency inference is made.

## Coordinate reproduction, optional and separate

To repeat the measurement rather than only summarize retained scores, acquire these public paths at the pinned revision into a separate workspace:

- `compute_RSMD.py` (SHA-256 recorded in the prior receipt).
- `esmfold_results/{baseline,ensemble}/{9NH7,9NFU}/` predictions.
- `RFDiffusion_structures/raw/{9NH7,9NFU}.pdb` references.
- `RSMD_RESULTS/cdr_rmsd_results.csv` and `RSMD_RESULTS/cdr_rmsd_framework.csv` for comparison.

Inspect the complete source before execution. The prior replay extracted only original function definitions and literal constants via Python AST, injected NumPy, `pathlib.Path`, `csv` and Biopython `PairwiseAligner`, and changed only `REFERENCE_PDB`, `PREDICTIONS`, `OUT_CSV`, `OUT_FRAMEWORK_CSV`. Call the extracted `main()` with baseline and ensemble only. Compare CDR rows by `(structure,method,design_id,cdr)` on RMSD, residue/atom counts, pLDDT and framework RMSD; compare framework rows by `(structure,method,design_id,domain)` on RMSD, residue count, identity and pLDDT, tolerance 1e-8. The receipt records reference hashes and replay dependency versions. This procedure remeasures coordinates, not model inference; do not import uninspected top-level team scripts.

No inferential method is chosen automatically. Main-cohort index pairing is not justified; a paired test cannot be inherited from ablation. A non-inferiority claim requires a meaningful pre-agreed margin and justified uncertainty analysis. Two cases do not establish a general format effect. These qualifications bound interpretation without preventing descriptive reporting of the selected linked inputs.
