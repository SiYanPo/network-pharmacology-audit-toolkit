# Null Models — Ready-to-Use Methodology Text (for SCI Manuscripts)

This file provides publication-ready English text explaining **why four null models** are used in
`np_null.py`, what each controls for, and how they jointly exclude a "random collision" explanation.
Paste the relevant paragraphs into the Methods/Results section and cite the software version recorded
in `run_manifest.json`. Replace bracketed placeholders with your actual numbers.

---

## 1. Why null models are required (rationale paragraph)

> Network-pharmacology intersections are prone to two trivial explanations that must be ruled out
> before any biological claim is made. First, hub bias: disease gene sets are enriched in highly
> connected proteins, so an arbitrary set of drug targets will tend to overlap them simply because
> hubs are numerous and well-connected. Second, size bias: larger compound-target sets intersect
> disease sets by chance more often than small ones. We therefore evaluated the observed overlap
> against four complementary null models, each preserving a different feature of the observed data
> while randomizing the feature of interest. An enrichment was considered robust only when it
> exceeded the expectation under **all** models that control for the corresponding confounder.

## 2. The four models and what each controls

### Model A — Uniform random gene sampling (control for set size only)

> In Model A, the same number of genes as the observed compound-target set was sampled uniformly at
> random from the background proteome (n = [N_background] genes), and the overlap with the disease
> gene set was recorded. This model controls for **set size alone** and tests whether the observed
> overlap exceeds chance expectation given the number of targets. Because it ignores the connectivity
> of the sampled genes, a significant result under Model A is necessary but not sufficient: it does
> not exclude hub bias.

### Model B — Degree-matched sampling (control for network connectivity / hub bias)

> Model B addressed hub bias by drawing random gene sets whose members were matched to the observed
> targets on their interaction degree in the protein-protein interaction network. For each observed
> gene with degree d, a random background gene with the same (or nearest) degree was selected, and
> the overlap with the disease set was recomputed. This preserves the **degree distribution** of the
> query set while randomizing gene identity, thereby asking whether the overlap is explained merely
> by the fact that the compound targets are highly connected. If the observed overlap is not greater
> than the degree-matched expectation, any apparent enrichment should be attributed to hub bias and
> the phrase "core disease module" must be withheld.

### Model C — Component/module-preserving rewiring (control for network topology)

> Model C preserved the modular architecture of the interaction network. Genes were randomly
> reassigned among network communities (Louvain modules) while keeping the number of genes per
> module fixed, so that the query set retained its **module membership distribution**. This tests
> whether the overlap arises from a genuine concentration in disease-relevant network neighborhoods
> rather than from a non-random distribution across topological modules.

### Model D — Compound-count-matched permuted targets (control for promiscuity / number of compounds)

> Model D controlled for target promiscuity and the number of candidate compounds. Targets were
> permuted within the set of proteins targeted by the same number of compounds, so the **distribution
> of compound counts per target** was preserved. This addresses the concern that multi-target
> compounds acting on promiscuous proteins inflate overlap independently of disease relevance.

## 3. Joint decision rule (the "B > C > A" hierarchy)

> The four models are nested in stringency with respect to the confounder they remove. We applied a
> conservative hierarchical decision rule: an overlap was interpreted as disease-relevant only if it
> was significant under the degree-matched model (Model B), which is the most direct control for the
> dominant confounder in PPI-based analyses. Significance under Model C further supports a
> module-level (rather than single-hub) signal, and significance under Model D excludes a
> promiscuity-driven artifact. Model A serves as the baseline. In formal terms, the evidentiary
> strength follows B > C > A for connectivity-controlled claims, with Model D as an independent
> promiscuity check. Empirical P values were computed as
> P_emp = (1 + number of null replicates ≥ observed) / (1 + number of replicates), with
> [R = 1000] permutations per model and random seed 42 for reproducibility.

## 4. Results sentence template (fill in your numbers)

> The observed intersection of [k] genes was significantly greater than expected under the uniform
> model (mean [mA], P_emp = [pA]), the degree-matched model (mean [mB], P_emp = [pB]), the
> module-preserving model (mean [mC], P_emp = [pC]), and the compound-count-matched model
> (mean [mD], P_emp = [pD]). Because the overlap remained significant under the degree-matched null,
> the enrichment cannot be explained by hub bias alone, supporting a disease-module-associated
> mechanism rather than random collision.

## 5. Caveats to state explicitly (reviewer-proofing)

> Null models assess statistical non-randomness, not biological causality. A significant enrichment
> indicates that the overlap is unlikely under the specified randomization schemes; it does not by
> itself establish that the identified targets mediate the therapeutic effect. These results should
> be read together with the dose-exposure feasibility analysis (`np_dose.py`) and the external
> gold-standard blind test (`np_goldtest.py`), and confirmed experimentally.
