# Discussion Defense Paragraphs — Publication-Ready English Text

Four defensive paragraphs addressing the four most damaging critiques of network-pharmacology
studies of TCM formulae. Each is self-contained, publication-ready, and explicitly states the
limitation rather than concealing it. Replace bracketed placeholders with study-specific values.

---

## 1. Monomer ADMET vs. Decoction Mixture (Soup Effect)

> Several limitations of the ADMET assessment warrant explicit acknowledgement. Our pharmacokinetic
> predictions were derived for individual compounds using monomer-based models (SwissADME and
> pkCSM), whereas the intervention evaluated experimentally is a multi-component aqueous decoction.
> In such a mixture, co-dissolution may enhance the apparent solubility of poorly soluble
> constituents, and simultaneous exposure to multiple phytochemicals may produce competitive
> inhibition of CYP450 and UGT enzymes, thereby altering systemic exposure relative to the
> single-compound prediction. Furthermore, intestinal microbiota can biotransform glycosides and
> polysaccharides into aglycones and short-chain metabolites that were not represented in our
> monomer-based exposure model. Although a simplified linear exposure proxy (daily dose × content ×
> extraction recovery) was used to weight the compound–target network, this proxy does not capture
> these physicochemical and microbial transformations. No open-source substitute currently reproduces
> the mixture-level PBPK capabilities of commercial platforms, and we therefore deliberately refrain
> from quantitative claims about the bioavailability of the whole decoction. The predicted exposure
> rankings should be regarded as hypothesis-generating and require confirmation by plasma
> pharmacokinetic profiling of the administered formulation.

---

## 2. Temporal Bias / Data Leakage in Gold-Standard Validation

> To guard against circular validation, we acknowledged the risk of temporal information leakage in
> the external gold-standard comparison. Disease–target repositories such as DrugBank and CTD
> accumulate annotations over time; a target annotated as disease-associated in the current release
> may have been unknown at the time the constituent targets were first compiled. If the evidence
> underpinning our predictions and the gold standard derive from the same literature corpus, the
> apparent accuracy may be inflated. We therefore performed a temporal-split validation in which
> only associations established before [cutoff year, e.g., 2020] were used to define the reference
> set for method development, and associations reported after that date were reserved as an
> independent blind test set. [n] targets approved or annotated after the cutoff were detected in
> our prediction set and were flagged as potential leakage; these were excluded from the primary
> recall estimate. This temporal separation provides a more conservative and unbiased assessment
> than random cross-validation and addresses the reviewer's concern that the validation may have
> been contaminated by future knowledge.

---

## 3. Apo/Holo Conformational Selection in Molecular Docking

> The molecular docking results should be interpreted in light of conformational selection. The
> protein structures retrieved from the PDB are frequently in the apo (ligand-free) state, and
> alpha-helical or loop rearrangements of side chains occur upon binding of different ligands.
> Docking against a single apo conformation, or against an AlphaFold-predicted model that represents
> one static state, can therefore produce false-negative scores by failing to sample the
> binding-competent conformation. For each target, we annotated the conformational state of the
> structure used (apo, holo, or AlphaFold model) and found that [n] of [N] targets were represented
> only by apo or predicted structures. For these targets, the reported binding energies should be
> regarded as indicative rather than definitive. Where a homologous holo structure was available,
> we recommend ensemble docking across multiple conformational states or homology-model
> replacement of the binding site; in the absence of such a structure, experimental validation
> (e.g., surface plasmon resonance or isothermal titration calorimetry) is required before any
> structure–activity claim is made. We have deliberately avoided interpreting small differences in
> docking score as evidence of selectivity.

---

## 4. Dose–Response Opacity of the Linear Exposure Proxy

> The dose-weighting scheme represents a deliberate simplification and should not be misconstrued
> as a quantitative pharmacokinetic model. The weight assigned to each constituent was computed as
> a linear proxy based on daily herbal dose, reported content, and extraction recovery, and was
> intended only to differentiate constituents on a relative exposure scale. This formulation does
> not account for thermal degradation of labile constituents during decoction, hydrolysis of
> glycosides, or the formation of insoluble complexes through chelation and Maillard-type reactions
> between components. In particular, macromolecular constituents such as polysaccharides and
> peptides—which may contribute to immunomodulatory and gut-microbiota-mediated effects—were not
> incorporated because their molecular targets and exposure are not amenable to the same
> compound-level weighting. Consequently, the dose-weighted topology should be interpreted as a
> relative prioritization of candidate constituents under the specified posology rather than a
> prediction of absolute tissue concentrations. We acknowledge that the true dose–response
> relationship of a multi-component decoction remains a black box, and that resolving it will
> require quantitative chemical profiling of the decoction coupled with exposure–response
> modelling—an effort beyond the scope of the present network analysis.
