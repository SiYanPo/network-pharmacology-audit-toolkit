# Safety, Promiscuity, Pocket and Wet-Lab Methods — Publication-Ready Text

Methodology paragraphs for the five new capabilities added in the final engineering round.
Each is self-contained and can be pasted into the Methods section; replace bracketed placeholders.

---

## 1. Compatibility Taboo and Toxicity Screening (合规红线)

> Before network construction, the prescription was screened against traditional Chinese medicine
> compatibility taboos and toxicity classifications. The eighteen incompatible pairs (Shi Ba Fan)
> and nineteen antagonistic pairs (Shi Jiu Wei) were encoded from the classical formulae, and herb
> names were normalised to account for processed forms (e.g., processed Pinellia → Pinellia). Toxic
> materials were graded according to the Chinese Pharmacopoeia (2020 edition) as highly toxic,
> toxic, or mildly toxic. Any formula triggering an incompatibility pair was flagged as a hard
> error and analysis was terminated; toxic herbs were retained but annotated and discussed with
> respect to dosage and safety monitoring. The screening was implemented as a deterministic rule
> engine (`np_contraindication.py`) and is reported in [Table Sx]. We note that the classical
> taboos remain scientifically debated; the screen is presented as a regulatory and safety
> declaration rather than as a criterion for component–target inclusion. This check does not
> substitute for review by a licensed pharmacist.

---

## 2. Promiscuity Penalty for Hub Prioritisation (多效性惩罚)

> To avoid the well-known degree bias in network pharmacology, where promiscuous compounds such as
> quercetin and kaempferol inflate the connectivity of non-specific targets, we applied a
> promiscuity penalty when prioritising docking targets. For each target t, we defined
> multiplicity m_t as the number of distinct compounds connected to t. Targets exceeding an
> absolute threshold (m_t > [200]) or the [95th] percentile of the multiplicity distribution were
> flagged as `nonspecific_hub`. A multiplicative penalty weight was assigned as
> w_t = 1 / (1 + log2(m_t / m_median)) for m_t above the median, and 1.0 otherwise; thus a target
> connected to twice the median number of compounds received approximately half weight, four
> times received one-third, and so on. Penalised targets were retained but deprioritised in the
> docking queue, and both raw multiplicity and penalty weight are reported in [Table Sy]. This
> correction prevents promiscuous, low-specificity proteins from being mistaken for disease
> hub genes solely on the basis of degree.

---

## 3. Pocket-Centred Docking Box Definition (精准对接口袋)

> Docking was not performed as blind docking over the entire protein surface. For each target, the
> binding pocket was first localised using a hierarchical strategy. When a co-crystallised ligand
> was present in the PDB entry, the docking box centre was set to the centroid of the ligand's
> heavy atoms, and the box dimensions were derived from the ligand's spatial span with a margin
> [of 8 Å, minimum 20 Å]. Ligand identity was determined either from the annotated ligand code
> or automatically by selecting the largest non-solvent organic hetero-residue. When no
> co-crystal ligand was available, fpocket [or P2Rank] was used to predict candidate pockets and
> the top-ranked pocket(s) were used to define the centre. If neither a co-crystal ligand nor a
> pocket-prediction tool was available, the target was flagged as `pocket_unresolved` and manual
> specification of the centre was required rather than resorting to blind docking. The pocket
> source and confidence level for each target are recorded in [Table Sz], and redocking of the
> co-crystal ligand yielded an RMSD of ≤2 Å for all included structures.

---

## 4. PCR Primer Design for Animal Validation (湿实验闭环)

> To bridge computational predictions with in-vivo verification, orthologues of the core human
> targets were mapped to mouse using a reviewed one-to-one orthologue table with MyGene.info
> cross-checking; one-to-many mappings were retained in full and flagged as ambiguous for separate
> experimental validation. For each mouse target, the corresponding RefSeq mRNA sequence was
> retrieved through NCBI E-utilities, and PCR primer pairs were designed with constraints of
> 18–24 nt length, 40–60% GC content, a melting temperature near 60 °C, and an amplicon of
> 60–250 bp. Melting temperatures were calculated using a modified Wallace approximation and are
> reported as approximate values; a Primer3-compatible batch file was additionally generated for
> exact design in an environment running Primer3. NCBI Primer-BLAST does not provide a public REST
> API; therefore all primer pairs must be subjected to in-silico PCR or Primer-BLAST specificity
> checking, and the optimal annealing temperature must be determined experimentally by gradient
> PCR. Primer sequences and predicted Tm values are provided in [Table Sa].

---

## 5. Pre-Submission Integrity Gate (投稿门禁)

> Before generating the submission-ready report, an automated integrity gate was applied. Three
> hard criteria were required: (i) a completed temporal-split blind validation, (ii) explicit
> annotation of the conformational state (apo/holo/AlphaFold) of every docking structure, and
> (iii) a dose-weighted topology analysis. Four additional soft checks (compatibility screening,
> external gold-standard blind test, figure QC with zero failures, and evidence-anchor validation
> of the Discussion) were reported but did not block generation. Only when all hard criteria were
> satisfied was a `Ready_for_Submission` watermark applied; otherwise the report was stamped
> `DRAFT` with the specific unmet criteria listed. This watermark attests solely to completion of
> the internal checklist and does not imply peer acceptance; all conclusions remain subject to
> experimental validation and responsible interpretation by the authors.
