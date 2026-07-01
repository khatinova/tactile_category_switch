# Tactile Category-Switch — Analysis Pipeline

Clean, ordered pipeline for the tactile category-switch reversal-learning study
(deterministic **D** vs probabilistic **P** feedback blocks; KH and RR cohorts).

Run the scripts **in numbered order**. Each `S*` script declares its inputs and
outputs in its header. Shared helpers live in `pipeline/utils/` and must be on
the MATLAB path (`addpath(genpath('pipeline'))`).

> The original scripts are left untouched in the repo root so you can diff.
> The cleaned, renamed, variable-unified versions live in `pipeline/`.

---

## Order of execution

| Step | Script | Input | Output |
|------|--------|-------|--------|
| **S1** | `S1_extract_behaviour.m` | raw `trial_data.mat` per subject (KH + RR) | `all_trial_data.mat`, `behav_table.mat` (`group_T`) |
| **S2** | `S2_preprocess_epoch_eeg.m` | raw EEG (Curry/MFF) + `all_trial_data.mat` | cleaned/ICA `.set`, outcome/theta/phase epochs, `*_trial2epoch.mat`, **per-trial spine table** (`group_stage_table_<cohort>.mat`) |
| **S3** | `S3_extract_eeg_features.m` | epochs + spine table | per-trial features (FCz neg-peak/mean, P300, Theta, PLV) **+ per-stage FRN/RewP difference waves** → `group_stage_table_features_<cohort>.mat` |
| **S4** | `S4_merge_feature_tables.m` | KH + RR feature tables | `group_feature_table_combined.mat` |
| **S5** | `S5_behavioural_modelling.m` | `all_trial_data.mat`, `behav_table.mat` | Nassar / RW / HGF fits, `nassar_results.mat`, `model_comparison.mat`, `behav_table_RL.mat` |
| **S6** | `S6_behaviour_plots_stats.m` | `behav_table.mat` (+ RL latents) | behaviour figures + LME stats (transitions, `n_prev_P`, reversal) |
| **S7** | `S7_eeg_rq_analysis.m` | `group_feature_table_combined.mat` | 5-RQ EEG figures + manuscript-ready stats |

KH preprocessing (S2) is already complete; rerun S2 only for RR (`COHORT = 'RR'`).
S3 then produces an RR feature table that S4 merges with the existing KH table.

---

## Unified variable conventions (applied across all scripts)

The original scripts used many different names for the same quantity. The
cleaned pipeline standardises on the following. Use `pipeline/utils/kh_subject_id.m`
and `kh_to_numeric.m` rather than re-deriving these.

| Concept | **Canonical name** | Legacy names replaced |
|---------|--------------------|------------------------|
| Subject label (string) | `subj_id` = `"Ox03"` / `"Nc07"` (zero-padded) | `subjID`, `subject`, `subj_id` |
| Subject number (double) | `subj` | `subj`, `subj_num` |
| Cohort | `cohort` = `"KH"`/`"RR"` | `researcher`, `source` |
| Block index | `block` | `blocknum`, `block_number` |
| Within-block trial | `trial` | `trialnum`, `trial_in_block` |
| Continuous trial index | `trial_continuous` | — |
| Block type | `block_type` ∈ {`D`,`P`} | `block_structure`, `structCode`, `curr_block_type`; legacy `V` → `P` |
| Reversal stage | `stage` ∈ {`LN`,`LE`,`RN`,`RE`} | — |
| Outcome | `correct` (0/1 or Incorrect/Correct) | `perceivedCorrect` kept separate |
| Feedback validity | `false_fb` (logical; true = flipped) | `trueFB` (inverse), `pTrueFB` |
| Confidence | `confidence` (1–10) | `conf`; `conf_z` = within-subject z |
| Reaction time | `RT` (seconds) | `rt` |
| Reversal trial | `revTrial` | — |
| Previous block type | `prev_block_type` | — |
| Transition | `transition` = `prev→curr` (e.g. `D→P`) | — |
| Cumulative prior P blocks | `n_prev_P` | — |

`kh_subject_id('standardise', T)` repairs any table to guarantee
`subj_id` / `subj` / `cohort` are present and consistent (and keeps a `subjID`
alias for backward compatibility).

---

## EEG feature definitions (S3)

### Per-**trial** measures (single-trial; in the per-trial table)
- `FCz_neg_peak_amp` / `_norm` / `_lat` — minimum of the baseline-corrected FCz
  waveform in the FRN window (diagnostic single-trial negativity). `_norm` =
  divided by the subject's pre-stimulus baseline RMS.
- `FCzCz_mean_amp` / `_norm` — mean of the FCz+Cz waveform in the FRN window.
- `P300_amp` / `_norm`, `Theta_amp`, `PLV_fp` / `PLV_fs` (+ `_pairwise`).

> FRN and RewP are **not** stored per trial — they are difference waves and only
> make sense per condition (see below).

### Per-**stage** measures (difference waves; in the stage table)
Computed by `pipeline/utils/kh_compute_frn_rewp_by_stage.m` for each
`subj_id × block_type × stage` (true-feedback trials only):

```
FRN_amp  = mean over FRN window  of (ERP_incorrect − ERP_correct)   → NEGATIVE
RewP_amp = mean over RewP window of (ERP_correct − ERP_incorrect)   → POSITIVE
```

- **FRN** = loss − win (incorrect − correct): relative negativity to
  worse-than-expected outcomes (~250–350 ms, FCz/Cz).
- **RewP** = win − loss (correct − incorrect): the same component viewed as the
  reward-related positivity; reduced/absent for losses.
- Default windows: `FRN_win = [250 300]`, `RewP_win = [250 350]` ms. The classic
  literature window is ~250–350 ms for both — inspect grand averages and adjust.
- References: Holroyd & Coles (2002, *Psych. Review*); Proudfit (2015,
  *Psychophysiology*); Sambrook & Goslin (2015, *Psych. Bulletin*, meta-analysis).

---

## Figure styling

All plotting calls `pipeline/utils/apply_fig_style.m` (and `save_fig`):
ticks **outside**, no top/right box, Arial, consistent font/line sizes, white
background, vector-PDF + 300-dpi-PNG export — suitable for posters/thesis/talks.

---

## Cohort-specific settings (reference)

| | KH (Curry/Neuroscan + ANT) | RR (EGI/MFF) |
|--|--|--|
| Import | `loadcurry` | `pop_mffimport` |
| FCz / Cz | `FCz` / `Cz` | `E11` / `E7` |
| Parietal (P300) | `Pz,P1,P2` | `E62,E67,E72` |
| Frontal/ACC (theta) | `FCz,Fz,AFz` | `E11,E6,E16` |
| Somatosensory | `C3,C4,CP3,CP1,C5,CP5` | `E36,E104,E41,E103` |
| Outcome codes | cohort1 `10/11`; cohort2 `31/32/33/34` | event string `'rewa'` etc. |
| Reference | average ref (protect EOG/mastoid) | Cz (ch 61) → CAR |

---

## Modelling (S5)

- **Nassar (2010)** adaptive learning-rate change-point model — fitted `H`
  (hazard) and `β` (softmax). Latents: `omega`, `alpha_nassar`, `PE_nassar`,
  `certainty` (|θ−0.5|, prospective), `surprise` (ω·|δ|), `theta_nassar`.
- **Rescorla–Wagner** — dual learning rate (`α_pos`/`α_neg`) variant added.
- **HGF** — 3-level Hierarchical Gaussian Filter added as a competing model.
- Pipeline for each model: grid-init → `fmincon` → parameter recovery → posterior
  predictive checks → BIC/AIC model comparison (matched trial-level likelihood).

---

## Status / TODO

See the PR description for the per-step completion status and any items that
still require running in MATLAB to verify (no MATLAB runtime is available in the
build environment, so all scripts are static-reviewed, not executed).
