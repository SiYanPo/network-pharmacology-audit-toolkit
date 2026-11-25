# Ultimate-Form Methods — Publication-Ready English Text

Five additional Methods paragraphs for the final-form capabilities. Paste into the Methods section;
replace bracketed placeholders.

---

## 1. Disease-Context Exemption from the Promiscuity Penalty

> The promiscuity penalty was applied with a disease-context exemption to avoid discarding
> genuinely disease-relevant high-degree targets. For each target t with multiplicity m_t and
> baseline penalty w_base = 1/(1+log2(m_t/m_median)), a normalised disease-association score d_t
> (0–1) was obtained from CTD gene–disease associations for [disease]. The exemption coefficient
> was defined as exempt = max(0, d_t − c)/(1 − c) with c = [0.5], and the final weight as
> w_t = w_base + (1 − w_base) × exempt. Thus, targets strongly associated with the disease (high
> d_t) were restored towards unit weight and labelled `disease_exempt`, whereas promiscuous
> targets lacking disease evidence retained the full penalty and were labelled `nonspecific_hub`.
> For example, [PTGS2], despite connecting [110] compounds, was exempted in an inflammation
> context because of its high CTD association, preventing the erroneous exclusion of a genuine
> anti-inflammatory target. Both raw multiplicity and final weights are reported in [Table Sx].

---

## 2. Decoy Negative Control and Z-score for Docking Specificity

> To distinguish specific binding from non-specific hydrophobic adsorption, a decoy negative
> control was applied. For each target, [ten] decoy molecules with comparable physicochemical
> properties (molecular weight, LogP, hydrogen-bond donors/acceptors) but no annotated activity
> were docked using identical Vina parameters. The specificity Z-score was computed as
> Z = (E_core − mean(E_decoy)) / SD(E_decoy), where E denotes the best binding energy. A target
> was considered specifically bound only when Z < −2; −2 ≤ Z < −1 was treated as borderline, and
> Z ≥ −1 was flagged as `Red Flag: High False Positive Risk`. Where DUD-E curated decoys were
> unavailable, property-matched decoys were generated locally and explicitly labelled as weak
> controls; formal claims of specificity therefore rely on the DUD-E subset. All decoy sources,
> energies, and Z-scores are provided in [Table Sy].

---

## 3. Single-Cell Projection of Core Targets

> To resolve the cell-type specificity of the predicted mechanism, core targets were projected
> onto a single-cell transcriptomic atlas. A publicly available GEO dataset ([GSE accession]) was
> processed with Scanpy; low-quality cells and genes were filtered, data were normalised and
> log-transformed, and cell subtypes were annotated using [marker genes/reference]. For each core
> target, the mean expression and the fraction of expressing cells were computed per subtype and
> visualised as a dot plot, with violin and UMAP expression plots as complementary views. Targets
> not detected in the atlas were reported separately rather than imputed. This analysis localises
> the predicted pharmacological action to specific cell populations (e.g., Kupffer cells versus
> endothelial cells) and is presented in [Figure Sx]. When a processed h5ad object was unavailable,
> the pipeline produced only a synthetic demonstration output, which was excluded from the
> manuscript.

---

## 4. Local Primer Specificity and MD Bridging

> Primer specificity was verified in silico before synthesis. Candidate primer pairs were aligned
> against the mouse reference genome (GRCm39) using local BLAST+/in-silico PCR; pairs producing
> more than one amplicon of the expected size range were discarded. When local alignment tools
> were unavailable, an executable verification script was generated and the result was reported as
> `TOOL_MISSING` rather than assumed specific; such primers were not used for experimental
> ordering. For structural validation, the top-ranked Vina complexes were subjected to short
> molecular-dynamics simulations. GROMACS input files for energy minimisation, 100 ps NVT and
> 100 ps NPT equilibration, and a 10 ns production run were generated automatically from the
> docked complex, with ligand topologies prepared by [acpype/antechamber] rather than guessed.
> Root-mean-square deviation and fluctuation were calculated to assess complex stability and
> flexible regions; binding free energies were estimated by MM/PBSA where available.

---

## 5. Subcellular Localisation Layered Rendering

> Subcellular localisation was incorporated to avoid biologically implausible compound–target
> pairings. Localisation terms for each target were retrieved from UniProt (or the Human Protein
> Atlas) and grouped into membrane, cytoplasm, nucleus, mitochondrion, and organelle layers.
> The compound–target network was rendered as a layered graph according to these categories.
> A red flag was raised when a compound predicted to have poor membrane/brain penetration was
> paired with a nuclear target, prompting explicit discussion of nuclear uptake. [n] such
> potential mismatches were identified and are listed in [Table Sz]; these were retained but
> annotated rather than silently removed.
