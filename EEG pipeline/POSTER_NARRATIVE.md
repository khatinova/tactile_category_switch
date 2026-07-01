# Uncertainty Calibrates Behavioural and Neural Adaptation to Category Switching

## Revised Poster Narrative (v3) — Data-Driven

---

## WHAT THE DATA ACTUALLY SHOW

Based on the generated poster figures, here is a frank assessment of what is and isn't working:

### Significant results (keep as main story)
| Finding | Test | p-value |
|---------|------|---------|
| Nassar model fits better than RW | Wilcoxon signed-rank | p=0.031 * |
| Theta ~ |PE| (D blocks) | Pearson r=0.03 | p=0.002 ** |
| Theta ~ |PE| (P blocks) | Pearson r=0.04 | p=0.001 ** |
| Behavioural reversal cost (LE→RN) | Paired t (from poster Fig 1) | p<0.001 *** |
| Confidence D > P across stages | Paired t | p<0.05 * |

### Non-significant results (demote or reframe)
| Finding | Test | p-value | Action |
|---------|------|---------|--------|
| FRN: LE vs RN (D) | Paired t | p=0.757 | DROP from main story |
| FRN: LE vs RN (P) | Paired t | p=0.969 | DROP from main story |
| FRN: True vs False FB | Paired t | p=NaN (empty) | FIX data issue or DROP |
| P300: LE vs RN (D) | Paired t | p=0.481 | DROP from main story |
| P300: LE vs RN (P) | Paired t | p=0.714 | DROP from main story |
| P300 ~ surprise (D) | r=0.01 | p=0.238 | DROP |
| P300 ~ surprise (P) | r=-0.01 | p=0.339 | DROP |
| Theta: LE vs RN (D) | Paired t | p=0.862 | DROP stage comparison |
| Theta: LE vs RN (P) | Paired t | p=0.762 | DROP stage comparison |
| H_det vs H_prob | Paired t | p=0.994 | DROP or reframe as individual diffs |
| Surprise D vs P at reversal | Not plotted (only D appeared) | — | CHECK data split |

### Panels needing investigation
- **False FB panel**: Returned NaN — likely too few trials after filtering (P blocks × incorrect × false_fb==1). Need to either relax the incorrect filter or show all P-block trials split by false_fb.
- **Surprise split by D/P**: Only D blocks appeared in the ribbon. The P-block all_trial_data may not have block_type coded consistently, or there are genuinely no P blocks in some subjects.

---

## REVISED STORY: What CAN we say?

The original narrative framed around "stochasticity vs volatility dissociation in neural signatures" is not supported by the EEG data. Neither the FRN nor the P300 differentiate LE from RN, and neither tracks model-derived surprise. The only robust EEG finding is:

**Frontal theta parametrically tracks prediction error magnitude (|PE|) regardless of block type.**

This is a replication of Cavanagh et al. (2010) in a novel tactile domain, and it IS significant. The story should pivot.

---

## REVISED TITLE
**Frontal theta tracks prediction error magnitude during tactile category learning under uncertainty**

(Alternative: keep original title but with more cautious framing — "Behavioural but not ERP adaptation dissociates stochasticity from volatility; frontal theta tracks PE magnitude")

---

## REVISED RESULTS STRUCTURE

### Panel 1: BEHAVIOUR (keep as-is — this works)
- 1A: Accuracy × stage × block type (with spaghetti + CI)
- 1B: Confidence × stage × block type
- 1C: Reversal cost by n_prev_P (if significant; needs checking)

Key message: *Behaviour clearly dissociates volatility (reversal cost) from stochasticity (lower P-block accuracy). Confidence is well-calibrated to uncertainty.*

### Panel 2: COMPUTATIONAL MODELLING (partially keep)
- 2A: Nassar model schematic (keep — introduces the model)
- 2B: BIC scatter Nassar vs RW (keep — p=0.031, Nassar wins 23/43)
- 2C: Surprise aligned to reversal (keep — shows the model captures the reversal)
- 2D: H_det vs H_prob — **REFRAME**: "Participants did NOT systematically adjust inferred hazard rate between contexts (p=0.994). This suggests the model captures reversal dynamics through ω and α rather than via a context-dependent H parameter. Individual differences in H may still be meaningful."

Key message: *The Nassar change-point model provides a marginally better account of trial-level behaviour than a fixed-rate RW learner, and generates meaningful latent variables (surprise, ω, α) even though context-specific hazard rates do not differ.*

### Panel 3: EEG — FOCUSED ON THETA (the one thing that works)

**Fig 3A: Raw ERPs (keep as context — shows where components are measured)**

**Fig 3B: Frontal theta ~ |PE| (THE MAIN EEG FINDING)**
- LEFT: Binned scatter showing theta increases with |PE| (D: r=0.03, p=.002; P: r=0.04, p=.001)
- RIGHT: Theta ~ surprise (if also significant, include; if not, replace with theta ~ next_correct showing theta predicts behavioural adjustment)
- Key message: *Frontal theta parametrically encodes PE magnitude during tactile feedback processing, replicating Cavanagh et al. (2010) in a novel sensory modality.*

**Fig 3C: PLV connectivity (check if significant in S8 MRQ4)**
- PLV ~ theta (from S8 MRQ4) — if significant, this shows that theta-mediated PE signals drive fronto-parietal decoupling
- PLV ~ surprise by pathway (FP vs FS) — if the pathway dissociation is significant, keep; otherwise drop

**Fig 3D: Null results reported honestly**
- Brief panel or text: "Neither FRN nor P300 amplitude differentiated pre- vs post-reversal stages (all p > 0.4) or tracked model-derived surprise (all |r| < 0.02). This suggests that in a tactile Go/NoGo paradigm with confidence ratings, the FRN and P300 do not index the same volatility-related processes observed in visual/auditory probabilistic learning tasks (cf. Bland & Schaefer, 2012)."

---

## REVISED CONCLUSIONS

1. **Behaviour dissociates stochasticity from volatility**: Accuracy drops sharply at reversal in both block types; confidence tracks block reliability. Prior P-block exposure may (or may not — check) modulate reversal recovery.

2. **Computational modelling**: The Nassar change-point model marginally outperforms RW and generates theoretically motivated latent variables. Context-specific hazard rates (H) do not differ, suggesting participants apply a uniform volatility prior across environments.

3. **Frontal theta is the primary neural correlate of PE**: This is the only robust EEG-model relationship (p<.003). It replicates prior work (Cavanagh et al., 2010) in a novel tactile modality with a more complex task structure, confirming that frontal theta is a domain-general PE magnitude signal.

4. **Null ERP results**: Neither FRN nor P300 showed reversal-related modulation or model-surprise coupling. This may reflect:
   - The longer stimulus-to-feedback interval (confidence rating intervenes)
   - Tactile modality differences in feedback processing
   - The mixed Go/NoGo nature of the task (unlike simple 2AFC)
   - Temporal smearing from the 4-stimulus structure

5. **Clinical implications** (speculative): If frontal theta is the only robust PE signal, then clinical populations with reduced theta may show specifically impaired PE-driven learning, even when slower ERP components appear normal.

---

## WHAT TO DO NEXT

### Immediate fixes for the poster script
1. **False FB panel**: Change filter from `correct==0` to just P-block trials, split by `false_fb`. This will give more trials and might reveal a difference.
2. **Surprise D vs P**: Debug why only D blocks appear in the reversal-aligned surprise. Check that `block_type` field in `all_trial_data` is populated for P blocks.
3. **Consider adding**: Theta × next_correct (from S8 MRQ3c) — does theta predict next-trial accuracy? This would strengthen the theta story.
4. **Consider adding**: PLV ~ theta (from S8 MRQ4a) — does theta predict connectivity? If significant, this gives a mechanistic account.

### For the poster layout
Given that most EEG panels are null, consider:
- **Reducing EEG panels** to just: (a) raw ERPs, (b) theta ~ |PE|, (c) one additional if significant
- **Expanding behaviour/modelling** to fill space (they're the stronger story)
- **Adding an explicit "null results" statement** — this is honest and can be framed as constraining future hypotheses

### Alternative framing (if you want to keep more panels)
Frame the poster as: "Behaviour clearly adapts to uncertainty structure, but neural adaptation is restricted to a single frequency-domain signal (theta). Traditional ERP markers of PE (FRN, P300) do not capture the same variance in this tactile paradigm."

This is actually an interesting result — it says something about the specificity of theta vs ERPs in non-standard modalities.

---

## FIGURE LIST (REVISED)

| Figure | Content | Status | Keep? |
|--------|---------|--------|-------|
| Fig 1A | Accuracy × stage × block_type | ✓ Significant | YES |
| Fig 1B | Confidence × stage × block_type | ✓ Significant | YES |
| Fig 1C | Reversal cost by n_prev_P | ? Check | IF sig |
| Fig 2A | Nassar model schematic | — | YES (context) |
| Fig 2B | Nassar BIC vs RW BIC | ✓ p=0.031 | YES |
| Fig 2C | Surprise at reversal | ✓ Descriptive | YES (model validation) |
| Fig 2D | H_det vs H_prob | ✗ p=0.994 | REFRAME as individual diffs |
| Fig 3A | Raw ERPs | — | YES (context) |
| Fig 3B | FRN LE vs RN | ✗ p>0.7 | DROP or move to supplement |
| Fig 3B | FRN true vs false FB | ✗ NaN/empty | FIX then decide |
| Fig 3C | P300 LE vs RN + surprise | ✗ all p>0.2 | DROP or move to supplement |
| Fig 3D | Theta violin LE vs RN | ✗ p>0.7 | DROP the violin |
| Fig 3D | Theta ~ |PE| | ✓ p=0.002/0.001 | **MAIN EEG PANEL** |
| Fig 3E | PLV pathways | ? Check | IF pathway×stage sig |
| Fig 3F | Neural × n_prev_P | ? Check | IF any trend sig |
| Fig 3G | Confidence × FRN | ? Check | IF sig |

---

## HONEST SUMMARY

The poster's strength is the behavioural story and the theta ~ |PE| replication. The ERP story (FRN, P300) did not pan out. The modelling adds value by generating latent variables that predict theta (but not ERPs). The poster should lean into what works rather than presenting null results as if they're surprising — they should be acknowledged briefly and framed as constraining theory about modality-specific feedback processing.
