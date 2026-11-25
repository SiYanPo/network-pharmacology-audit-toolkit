# Ultimate v5 Methods — Publication-Ready English Text

Six Methods paragraphs for the final冲刺 capabilities. Paste into Methods; replace placeholders.

---

## 1. Environments and Reproducibility

> All computational environments were containerised (Dockerfile, Ubuntu 24.04) and an idempotent
> provisioning script (`setup_env.sh`) installed Scanpy/AnnData, AutoDock Vina, OpenBabel, GROMACS,
> NCBI BLAST+, isPcr and acpype. A pre-flight audit (`np_envcheck.py`) was run before each analysis
> to verify tool availability and versions; tools installed outside the system PATH (e.g., Vina
> v1.2.5 under `D:/MDtools`) were registered explicitly rather than assumed missing. Ligand
> topologies for molecular dynamics were generated automatically with acpype (GAFF2, AM1-BCC
> charges) from the Vina-docked SDF/PDB, producing GROMACS-compatible `.itp` and `.gro` files; no
> ligand force-field parameters were guessed manually.

---

## 2. Physiological Barrier Filtering (BBB)

> For central nervous system indications, constituents predicted to be BBB-negative were excluded
> from the docking queue because they are physically unable to reach intracranial targets. The
> filter used the ADMET BBB classification, with a literature-based whitelist exemption for
> constituents known to be transported across the BBB. Excluded constituents were recorded as
> `Physically Unreachable` and are listed in [Table Sx]. For non-CNS target tissues, membrane
> permeability was reported as an advisory note rather than an absolute filter. This spatial
> bridge between exposure and target localisation prevents mechanistically impossible
> constituent–target pairings from being presented as core findings.

---

## 3. Multi-Scoring Consensus and Uncertainty Quantification

> Binding affinity was not inferred from a single scoring function. For each ligand, Vina and
> Vinardo scores (and, where available, GNINA) were computed with identical receptor and search
> parameters. The consensus score was reported as the mean, and the uncertainty as the sample
> standard deviation across scorers. A ligand was flagged `Unstable_Binding` when the standard
> deviation exceeded [1.5] kcal/mol, indicating that the predicted affinity is scoring-function
> dependent; such ligands were excluded from quantitative structure–activity claims. Stable
> ligands additionally required a mean score ≤ [−7] kcal/mol. All per-scorer values are provided
> in [Table Sy].

---

## 4. Cell–Cell Communication via LIANA

> To determine whether core targets mediate intercellular signalling, we projected them onto a
> single-cell atlas and applied LIANA, a Python framework consolidating CellPhoneDB-, CellChat-
> and NATMI-style ligand–receptor inference. Cell-type annotations were retained from the source
> dataset and ligand–receptor interactions were ranked by the consensus method. Only interactions
> involving a core target as ligand or receptor were retained, and the communication strength was
> visualised as a source–target heat map [Figure Sx]. Targets without ligand/receptor annotation were
> reported as such rather than force-fit into the analysis.

---

## 5. Tissue-Specificity Validation Using the Human Protein Atlas

> To exclude cell-line artefacts, the tissue-level expression of each core target was checked
> against the Human Protein Atlas consensus RNA-seq data (nTPM across 62 tissues). For the disease
> tissue [tissue], an enrichment ratio was computed as nTPM(disease tissue) divided by the median
> nTPM across all tissues; targets with a ratio ≥ [2] were considered tissue-enriched. Targets that
> failed this criterion were flagged for discussion as potential false positives. Because HPA
> represents healthy tissue baseline, disease-associated upregulation was additionally cross-checked
> against GEO disease cohorts where available [Table Sz].

---

## 6. Negative-Control Decoys for Non-Canonical Targets

> For targets not covered by the DUD-E benchmark, decoys were generated with DeepCoy, a
> deep-learning model that produces physicochemical-property-matched but structurally distinct
> negatives conditioned on the binding pocket. Where DeepCoy was unavailable, DUD-E subsets or
> property-matched local decoys were used and explicitly labelled with their provenance tier.
> The specificity Z-score was computed as for the canonical targets; decoys of weaker provenance
> were not used to support definitive claims of selectivity.
