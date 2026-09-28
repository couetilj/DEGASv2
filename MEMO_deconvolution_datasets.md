# Memo — Datasets for Composition-Aware DEGAS in Cutaneous Melanoma

**To:** Justin (PI)
**From:** Claude, reviewing the Composition-Aware DEGAS plan
**Date:** 2026-05-07
**Branch:** `claude/review-deconvolution-datasets-t90ve`

---

## TL;DR

The plan is methodologically careful, but **the cohort it stands on cannot answer the question it asks**. None of the in-project melanoma assets (Gide, MIA/Garg, Jerby-Arnon, IUSM Xenium) are paired same-patient bulk + scRNA-seq, and none has matched normal skin. Every cross-modal step in Phases 1–4 is therefore between independent cohorts — i.e., the deconvolution will never be directly validated against patient-resolved ground truth in melanoma. Two correctable issues compound this:

1. **MIA/Garg 2021 is targeted RNA capture (rnaaccess), not whole-transcriptome.** That is documented in DATASETS.md but its consequence — that BayesPrism's gene intersection with a scRNA reference will collapse to ~the panel size — is not flagged in the plan.
2. **The plan refers to "IUSM bulk RNA-seq + outcomes," but DATASETS.md only lists IUSM *Xenium*.** Either an IUSM bulk cohort exists and isn't documented, or it doesn't exist and the plan needs adjusting.

Below: what we have, what's missing, what the literature offers, and a proposed path that adds the smallest amount of new data needed to make the deconvolution claim defensible.

---

## 1. What's actually in the project (melanoma)

| Asset | Type | N | WT or targeted | Outcome | Same patients as anything else? |
|---|---|---|---|---|---|
| Gide / Wilmott 2019 | Bulk RNA-seq | 91 samples | Whole transcriptome | RECIST + PFS + OS | **No** |
| Garg / MIA 2021 | Bulk **rnaaccess targeted capture** | 374 annotated (199 primary CM + 175 LN) | **Targeted (~2.5K probes class)** | Mets Y/N, 6-yr follow-up | **No** |
| Jerby-Arnon 2018 | scRNA-seq | 32 tumors / 7,186 cells / 10 cell types | n/a | ICI response (study-level) | **No** |
| IUSM Xenium | Spatial (Xenium) | 4 samples / 347,843 cells | Targeted Xenium panel | ICI TME | **No** (not paired with any bulk listed) |

**Implications**:
- The plan's "MIA bulk + outcomes" + "MIA scRNA-seq" framing is not supported — there is **no** MIA scRNA-seq in the project, and the MIA bulk is targeted capture. The plan acknowledges this conditionally ("if MIA scRNA-seq is unavailable") but it should be treated as a hard constraint, not a contingency.
- The Jerby-Arnon reference (n=7,186 cells) will hit BayesPrism's "≥50 cells per type" floor for minor cell types (CAFs, NK, plasmacytoid DC, mast cells). Phase 1's confidence-weighting plan handles this, but expect the reference itself to be the dominant source of attribution noise for stromal/rare populations.
- The Spatial example in the repo (`Spatial_example/H2_1/`) is from Erickson **prostate** ST, not melanoma — there is currently no melanoma spatial walkthrough in the codebase.

---

## 2. The data gap your plan needs to close

The composition-aware DEGAS framework will produce per-cell risk scores. The validation question — "are those scores reflecting cell-state biology or just composition?" — has two distinct flavors:

- **In silico (Phase 5a synthetic recovery)**: pseudobulks from sc data with known composition. The plan does this. It's necessary but not sufficient — it tests *the model*, not *the deconvolution on real bulk*.
- **In vivo (real patient)**: real bulk RNA-seq from a patient where we *also* have sc from that patient's tumor. This is what your "paired bulk + sc same patient" instinct is asking for. It is currently absent for melanoma in the project.

Without the in vivo arm, your final claims will be of the form: *"On Gide bulk, BayesPrism estimates X% T-cells; we cannot independently verify this for any individual patient."* That's the standard caveat in the deconvolution literature, and reviewers know it.

**For matched normal skin paired with tumor**: this is essentially absent in the public cutaneous melanoma literature. TCGA-SKCM has 1 normal sample. The acral melanoma single-cell paper (Li et al. 2022, eLife) has 4 paired tumor/paracancerous. Otherwise the field uses either (a) GTEx skin as a population-level normal, or (b) low-purity tumor samples as quasi-normals. Generating a small N+T paired cutaneous melanoma cohort would itself be a publishable contribution.

---

## 3. What the public literature offers

### Tier 1 — Truly paired same-patient tumor bulk + scRNA-seq (cutaneous or extracranial melanoma)

**This tier is much thinner than I initially reported.** On closer reading of the data-availability sections, the two cohorts I would have nominated (Pozniak 2024, Biermann 2022) deposit bulk RNA-seq only from **derivative cell lines**, not from the same physical tumor that was scRNA'd. They are NOT paired in the sense you want.

| Dataset | What's actually paired | N (truly paired tumor sc + tumor bulk) | Outcome | Access | Use as deconvolution validation? |
|---|---|---|---|---|---|
| **"Treatment Resistance" cohort (Sci Rep 2024, s41598-024-72255-9)** | sc + bulk + WES on the **same physical specimens** | ~15 patients / 18 tumors | Mixed (neoadjuvant ICI / BRAF-MEKi) | **Public, GEO** | **Yes — currently the only cutaneous melanoma cohort that fits the strict same-patient same-tumor criterion at non-trivial N.** Verify exact citation and accession before relying on it; the n is small and treatments are mixed. |
| **HTAN Melanoma Atlas** (`phs002371`) | Multimodal by design — bulk + scRNA + spatial intended on same patients | Per-patient pairing varies; needs Synapse manifest pull | Pre-malig → mets trajectories | Mixed: public + dbGaP-credentialed | **Maybe — depends on current manifest.** Pull the manifest before committing. |
| Tirosh 2016 (`GSE72056` sc + `GSE77940` bulk) | sc + bulk deposited separately | Verify per-patient overlap from sample IDs | Pre-treatment metastatic | **Public** | **Maybe — only if sample IDs match patient-by-patient.** Often cited as quasi-paired; treat as unproven until verified. |
| Pozniak 2024 (`EGAS00001006488` sc; cell-line bulk on KU Leuven RDR) | sc + Visium **on tumor**, bulk on **derived cell lines only** | **0 patient-level tumor bulk + tumor sc** | ICI response | **EGA credentialed** | **No — bulk is from cell lines, not the index tumor.** Still useful as an sc + spatial reference. |
| Biermann 2022 (`GSE200218` sc; bulk on 2 derived cell lines) | sc/snRNA + Slide-seq on tumor, bulk on cell lines | **2 patients (and bulk is from cell lines)** | Treatment-naive mets | **Public GEO** | **No — same problem.** Still a strong sc + spatial reference for brain mets. |

### Tier 2 — Larger sc references (no paired bulk) for BayesPrism reference building

| Dataset | N | Cells | Notes |
|---|---|---|---|
| Tirosh et al. 2016 (`GSE72056`) | 19 patients | ~4,600 cells | Smart-seq2; canonical reference; complements Jerby-Arnon. |
| Sade-Feldman et al. 2018 (`GSE120575`) | 48 patients | ~16,000 T-cells | ICI-treated, **immune-only** (CD45+ sort). Useful for immune sub-typing but not as a full TME reference. |
| Jerby-Arnon 2018 (`GSE115978`) | 32 patients | 7,186 | Already in project. |
| Pozniak 2024 sc fraction | several dozen | tens of thousands | If we acquire it, immediately upgrades the reference. |
| Li et al. 2022 acral melanoma (`HRA001262` / GSA) | 9 acral | ~140K cells | Acral subtype — useful as a contrast cohort, less for cutaneous reference. |

### Tier 3 — Bulk-only outcome cohorts (for training, not deconvolution validation)

Confirmed bulk-only with no same-patient sc: Gide 2019, Garg/MIA 2021, Hugo 2016, Riaz 2017, Liu 2019, Van Allen 2015, **TCGA-SKCM** (n=472, public; only 1 normal).

### Matched normal skin with tumor

Functionally absent in cutaneous melanoma. TCGA-SKCM (1), acral subset (4). Path forward is either GTEx skin as population-level normal, or institutional generation.

---

## 4. The simplest experiments that actually answer the question

The plan currently has 6 phases and ~5 architectural variants. As your PI, my pushback is that Phase 2's variant zoo can wait until Phase 1's deconvolution is validated *on real paired data*. The cleanest, smallest experiment that closes the central gap:

### Experiment A — Patient-level deconvolution validation (the missing experiment)

The path here is harder than I originally claimed because the obvious candidates (Pozniak, Biermann) turn out to be sc-on-tumor + bulk-on-cell-lines, not what we need. Realistic options, in order of preference:

1. **Sci Rep 2024 paired bulk+sc cohort (n≈15)** — verify the citation and accession, then use it as a small but truly same-patient validation set. n=15 is borderline for a per-patient correlation analysis but enough to compute aggregate Pearson across patients × cell types.
2. **Tirosh 2016 sample-ID cross-check (`GSE72056` sc vs `GSE77940` bulk)** — pull both metadata tables and check whether sample IDs overlap at the patient level. If yes, this is a free public paired cohort. If no, drop it.
3. **HTAN Melanoma Atlas (`phs002371`) manifest pull** — confirm whether any HTAN-Mel sub-cohort has bulk + sc deposited per patient. This is plausible but needs verification.
4. **In silico fallback** — if 1–3 all fail, the Phase 1 plan as written (synthetic mixtures from Jerby-Arnon) is the only validation we get. Document this honestly as the project's central methodological limitation in the manuscript.

The procedure for whichever cohort we land on:
- From the sc data, count cells per type per patient → "ground-truth" patient × cell-type proportion matrix `Π_true`.
- Run BayesPrism on the *same patients'* bulk using a **different** sc reference (Jerby-Arnon + Tirosh, NOT the cohort's own sc) — this is the deployment-realistic scenario.
- Compute per-patient agreement (Pearson r, RMSE) between `Π_BayesPrism` and `Π_true`.
- Pre-register a threshold: r > 0.6 for major cell types (tumor, T-cell, myeloid) on real paired patients. (Lower than 0.7 synthetic threshold — real data routinely loses 0.1 vs. synthetic.)

If this works on n=15 but not on synthetic, the synthetic generator is wrong. If it works on synthetic but not on n=15, BayesPrism is wrong for melanoma at our reference depth. Either outcome is informative.

### Experiment B — Demonstrate composition confounding exists in our outcome cohort

Already in Phase 0, but elevate it: on Gide bulk, plot per-patient mean DEGAS-baseline score against (a) ESTIMATE tumor purity *and* (b) BayesPrism tumor-cell fraction. If the correlation is weak, your motivation for the entire pipeline is weakened — flag this honestly.

### Experiment C — The IUSM scope question

If IUSM holds tissue from melanoma resections and you have block-level access, a small **n=15–20 paired cohort** (one block → bulk RNA-seq, adjacent block → scRNA-seq, third punch → matched normal skin from the same surgical field) would be:
- the only patient-paired T+N+sc cutaneous melanoma cohort in the field, to my knowledge;
- enough for a methods-paper validation arm;
- small enough to run alongside the in silico work without slipping the modeling timeline.

This is the elegant scientific-method answer to your gap. Cost is real (3 libraries × 20 = 60 libraries plus prep), but tractable; reviewer impact is high.

### Experiment D — Drop or downscope MIA/Garg

The rnaaccess targeted-capture issue is not something the plan can engineer around. Either:
- **Restrict** BayesPrism on MIA to the gene intersection between the panel and the sc reference (likely <2,000 genes), and **flag deconvolved proportions on MIA as low-confidence** by design;
- **Keep MIA only as an outcome cohort** (which is its strength — n=199 primary CM with binary mets) and don't deconvolve it;
- **Replace MIA with TCGA-SKCM** as the second outcome cohort: 472 patients, whole-transcriptome, public, includes survival.

I'd recommend the second: use MIA for survival/binary-mets training, do not deconvolve it, and bring TCGA-SKCM in for the cross-cohort attribution work in Phase 5d.

---

## 5. Recommended dataset acquisitions (in priority order, revised)

1. **TCGA-SKCM bulk** — public via GDC; whole-transcriptome; 472 patients with PFI/OS/distant mets via Liu 2018 CDR. Replaces or complements MIA for outcome modeling.
2. **Tirosh 2016 (`GSE72056` sc + `GSE77940` bulk)** — public; doubles sc reference cell count and patient diversity. Cross-check sample IDs to see if it actually constitutes a same-patient paired set (worth ~30 min of work for the answer).
3. **HTAN Melanoma Atlas manifest (`phs002371`)** — pull index; confirm whether any sub-cohort is per-patient paired bulk+sc. Could be the surprise win.
4. **Sci Rep 2024 "Treatment Resistance" paired cohort (n≈15)** — verify citation, then acquire as the validation set for Experiment A. Small but currently the cleanest public option.
5. **Pozniak 2024 sc + Visium (`EGAS00001006488`)** — EGA credentialed. Even though the bulk is cell-line-only, the sc + Visium is the strongest melanoma sc reference and aligns with the project's ICI focus. Start access application early.
6. **Biermann 2022 sc + Slide-seq (`GSE200218`, `GSE200278`)** — public; strong reference for brain-met sc but does *not* solve the paired-bulk gap.
7. **(Stretch — high impact) IUSM paired T+N+sc cohort generation** — if institutionally feasible, the highest-leverage addition we could make. The public field genuinely lacks this. Justifies a separate aim and is publishable on its own as a resource paper.

---

## 6. Specific edits I'd make to the existing plan

- **Phase 1 validation**: add Experiment A (real-data per-patient deconvolution check on Biermann or Pozniak) as a **prerequisite** for proceeding to Phase 2. Don't sweep architectures until the deconvolution is validated.
- **Phase 1 thresholds**: the current "r > 0.7 for dominant types on synthetic mixtures" is fine for in silico, but add a separate pre-registered threshold for in vivo (real paired patients) — e.g., r > 0.6 for tumor + T-cell + myeloid. Synthetic recovery routinely beats real-data recovery; setting two thresholds prevents inflated confidence.
- **Phase 2 variant pruning**: 5 variants × 5 folds × 5 seeds is a lot. Drop 2c (residualization) if 2d (adversarial) clearly beats 2a (baseline) on the leakage probe — the plan already has 2c flagged as "likely to lose information."
- **Phase 2e (composition-aware OT)**: keep as stretch unless you want to claim algorithmic novelty — the standard variant for a methods paper is 2d.
- **Outcome definition**: pick Cox time-to-event as primary (Phase 4c perturbation attribution is much cleaner with a continuous risk output, as the plan notes). Use binary mets-at-3-yr only as a sensitivity check on MIA.
- **Documentation**: update DATASETS.md to (a) mark MIA as "not for deconvolution use" and (b) clarify whether an IUSM bulk melanoma cohort exists or is planned.

---

## 7. Open questions for you

1. Does an IUSM **bulk RNA-seq melanoma cohort with outcomes** actually exist? The composition_degas plan cites it but DATASETS.md only documents IUSM Xenium.
2. Is there institutional appetite for generating a small (n≈20) paired T+N+sc cutaneous melanoma cohort? Given that the public field genuinely has nothing at meaningful N, this is even higher-leverage than I first thought.
3. EGA credentialing for Pozniak 2024 — do we already have UZ Leuven-VIB DAC access for related work, or do we need to apply fresh?
4. Cox vs. binary outcome — if we go Cox-primary, MIA's 6-yr-no-censoring binary becomes the sensitivity arm rather than a primary trainer. OK?
5. The "Sci Rep 2024 n=15 paired cohort" the literature scout surfaced — I haven't independently verified the citation/accession yet. Want me to track that down before we treat it as the Experiment A target?

---

## 8. Bottom line

The plan's modeling instincts are right. The data are not yet adequate to back the claims at the level you'd want for a methods paper, and **the public-data gap is more severe than I initially thought** — most "paired" melanoma cohorts in the literature pair sc-on-tumor with bulk-on-derived-cell-lines, not what we need. The two cheap moves are:

1. **Cross-check Tirosh 2016 + GSE77940 sample IDs** to see if there's a free public paired set hiding in plain sight.
2. **Acquire the Sci Rep 2024 n=15 paired cohort** (verify citation first), use it as Experiment A's validation set, and proceed to Phase 2 only after BayesPrism clears a real-data threshold there.

The high-leverage move is **IUSM-generated paired T+N+sc cutaneous melanoma**. The public field has nothing at meaningful N, and producing it would convert the project from "another adversarial-DA paper" to "the paper that introduced and validated composition-aware DEGAS on the only patient-paired melanoma cohort." That's the elegant scientific-method answer to your gap.
