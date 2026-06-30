% =============================================================================
% RUN_PIPELINE.m  —  Master runner: analysis pipeline steps S5 → S10
%
% Tactile category-switch reversal-learning study
% Deterministic (D) vs probabilistic (P) feedback blocks; KH and RR cohorts.
%
% HOW TO USE
% ----------
%   1. Edit the CONFIGURATION block below to point at your data folder.
%   2. Comment out (%) any steps you want to skip.
%   3. Run this file: >> run('RUN_PIPELINE.m')
%
% PREREQUISITES
% -------------
%   Steps S1–S4 must have been completed first (behaviour extraction,
%   EEG preprocessing, feature extraction, table merge).
%   See PIPELINE.md for the full pipeline description.
%
% PIPELINE ORDER (S5 onwards)
% ---------------------------
%   S5   — Behavioural computational modelling (Nassar / RW)
%   S5b  — Additional models (RW-dual, HGF) + 3-way comparison   [optional]
%   S6   — Behavioural figures, statistics, sequential-block analyses
%   S7b  — EEG feature-table validation / QC          [run BEFORE S7]
%   S7   — EEG research-question analysis (RQ1–RQ5)
%   S7e  — Hierarchical confidence models (RQ2/RQ3 standalone)    [optional]
%   S9   — RewP raw-ERP waveforms + poster figures                 [optional]
%   S10  — Morlet wavelet time-frequency analysis (ERSP/ITPC)       [optional]
%
% OUTPUTS
% -------
%   Each script saves figures and .mat files to sub-folders under
%   <base_path>/Results/.  Exact paths are printed by each script on run.
%
% NOTE: S7d_rq2_confidence.m is a helper FUNCTION called internally by S7.
%       It is not run directly from here.
% =============================================================================

%% ── 0. CONFIGURATION ────────────────────────────────────────────────────────
% Set base_path to your local project root. This value is printed to the
% console at startup and passed to each sub-script via the shared workspace.
% Each S-script also defines base_path near its own top — keep them in sync.

base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH';

% Add pipeline utilities (subject-ID helpers, figure style, etc.)
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'pipeline')));

fprintf('\n');
fprintf('================================================================\n');
fprintf('  TACTILE CATEGORY-SWITCH — ANALYSIS PIPELINE  (S5 → S10)\n');
fprintf('  Base path : %s\n', base_path);
fprintf('  Started   : %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('================================================================\n\n');

%% ── S5 | Behavioural computational modelling ────────────────────────────────
%
%  Models fitted:  Nassar (2010) change-point model, Rescorla-Wagner (single α)
%  Outputs:  nassar_results.mat, behav_table_June2026_RL.mat,
%            sim_data.mat, figures (sim_fig*, fig1–fig8, fig_blocktype_H)
%
%  INPUT  : <DATA>/all_trial_data.mat
%           <DATA>/behav_table_June2026.mat  (or behav_table.mat / group_T.mat)
%  OUTPUT : <Results>/Simulation results/Figures/

fprintf('─── S5: Behavioural computational modelling ────────────────────\n');
run(fullfile(fileparts(mfilename('fullpath')), 'S5_behavioural_modelling.m'));
fprintf('S5 complete.\n\n');

%% ── S5b | Additional models (RW-dual, HGF) — OPTIONAL ──────────────────────
%
%  Adds two field-standard alternatives for BIC model comparison:
%    (1) Dual-α Rescorla-Wagner (α_pos / α_neg)
%    (2) 3-level binary Hierarchical Gaussian Filter (HGF; Mathys 2011/2014)
%        Uses TAPAS tapas_hgf_binary if on path, else built-in implementation.
%
%  INPUT  : all_trial_data.mat, nassar_results.mat (from S5)
%  OUTPUT : model_comparison_RW_HGF.mat, model_comparison_RW_HGF.pdf

% fprintf('─── S5b: Additional models (RW-dual, HGF) ──────────────────────\n');
% run(fullfile(fileparts(mfilename('fullpath')), 'S5b_models_RW_HGF.m'));
% fprintf('S5b complete.\n\n');

%% ── S6 | Behavioural figures, statistics, sequential-block analyses ─────────
%
%  Produces stage-based performance, confidence, RT, stay-behaviour and LME
%  statistics (accuracy ~ stage × transition, ~ stage × n_prev_P).
%  Includes optional E1–E11 sequential-block figures (requires all_trial_data).
%
%  Key settings inside S6 (edit near top of file):
%    STAGE_WIN              — window width in trials (default: 15)
%    STAGE_ASSIGNMENT_POLICY — 'protect_edges' (recommended)
%
%  INPUT  : <DATA>/behav_table_June2026_RL.mat  (from S5; falls back to behav_table.mat)
%           <DATA>/all_trial_data.mat            (for E1–E11 sequential figs, optional)
%  OUTPUT : <Results>/Behav results/Stage Figures_<WIN>trial_protect_edges/

fprintf('─── S6: Behavioural figures + statistics ───────────────────────\n');
run(fullfile(fileparts(mfilename('fullpath')), 'S6_behaviour_plots_stats.m'));
fprintf('S6 complete.\n\n');

%% ── S7b | EEG feature-table validation / QC  ────────────────────────────────
%
%  Run THIS BEFORE S7.  Produces trial-count tables, amplitude distributions,
%  split-half reliability, outlier flags, grand-average waveforms, and
%  cohort (KH vs RR) comparison plots.
%  Based on Clayson et al. (2020) Psychophysiology best-practice recommendations.
%
%  INPUT  : group_feature_table_combined.mat  (from S4)
%  OUTPUT : <Figures>/Validation_S7b/*, S7b_validation_workspace.mat, CSV tables

fprintf('─── S7b: EEG validation / QC ───────────────────────────────────\n');
run(fullfile(fileparts(mfilename('fullpath')), 'S7b_EEG_validation_checks.m'));
fprintf('S7b complete.\n\n');

%% ── S7 | EEG research-question analysis (RQ1–RQ5) ──────────────────────────
%
%  Five research questions:
%    RQ1 — FRN / P300 under uncertainty (D vs P blocks × stage)
%    RQ2 — Confidence × prefrontal negativity (FCz mean, true-FB trials)
%    RQ3 — Frontal theta × stage × block type
%    RQ4 — Fronto-parietal PLV (4-8 Hz) × reversal stage
%    RQ5 — Fronto-parietal vs fronto-somatosensory PLV pathway comparison
%
%  Key setting inside S7 (edit near top of file):
%    PLOT_SCALE  — 'norm' (recommended) | 'z' | 'raw'
%
%  INPUT  : group_feature_table_combined.mat  (from S4)
%           frn_rewp_by_stage_combined.mat     (optional; for grand-avg waveforms)
%  OUTPUT : <Figures>/RQ_analysis_combined_extended_electrodes/S7_plot_<SCALE>/
%           manuscript_stats.txt

fprintf('─── S7: EEG RQ analysis (RQ1–RQ5) ─────────────────────────────\n');
run(fullfile(fileparts(mfilename('fullpath')), 'S7_eeg_rq_analysis.m'));
fprintf('S7 complete.\n\n');

%% ── S7e | Hierarchical confidence models (RQ2/RQ3) — OPTIONAL ──────────────
%
%  Standalone nested LME model sequence for RQ2 (prefrontal mean) and RQ3
%  (frontal theta), comparing M0–M4 with likelihood-ratio tests.
%  Useful for pre-registration or as an extended methods section.
%
%  INPUT  : group_feature_table_combined.mat  (from S4; or gt already in workspace)
%  OUTPUT : printed model-comparison tables (no saved figures)

% fprintf('─── S7e: Hierarchical confidence models ────────────────────────\n');
% run(fullfile(fileparts(mfilename('fullpath')), 'S7_RQ2_RQ3_hierarchical_confidence_models.m'));
% fprintf('S7e complete.\n\n');

%% ── S9 | RewP raw-ERP waveforms + poster figures — OPTIONAL ─────────────────
%
%  Produces:
%    Fig 1 — RewP difference waves per stage × block type
%    Fig 2 — Grand-average ERP (correct / incorrect / false-correct /
%             false-incorrect) with RewP difference wave below each panel
%    Fig 3 — First-10-trials accuracy learning curves (D vs P)
%    Fig 4 — Mock raw EEG traces for poster illustration
%    Fig 5 — Parietal P300 grand-average by correct/incorrect × D/P
%    Fig 6 — RewP modulated by transition type and n_prev_P history
%
%  FRN/RewP window mode (edit inside S9):
%    WINDOW_MODE = 'broad'   → 200–400 ms  (recommended; covers observed peaks)
%    WINDOW_MODE = 'narrow'  → 250–350 ms  (original pipeline window)
%
%  INPUT  : group_feature_table_combined.mat, frn_rewp_by_stage_combined.mat,
%           all_trial_data.mat, behav_table.mat
%  OUTPUT : <Figures>/Poster_P9/

% fprintf('─── S9: RewP raw-ERP and poster figures ────────────────────────\n');
% run(fullfile(fileparts(mfilename('fullpath')), 'S9_RewP_rawERP.m'));
% fprintf('S9 complete.\n\n');

%% ── S10 | Morlet wavelet time-frequency analysis — OPTIONAL ─────────────────
%
%  Computes ERSP + ITPC (Morlet wavelet) on the outcome-locked epochs and
%  appends single-trial TF features (theta_ersp, alpha_ersp, beta_ersp,
%  theta_itpc + _z versions) to the SAME combined table used by S7.
%  Also runs a self-contained 2-D cluster permutation test (no FieldTrip).
%
%  REQUIRES: EEGLAB on path + the epoched .set files in
%            Epoched_data_noisefiltering/ (not just the feature table).
%            KH cohort only (valid_participants edited inside S10).
%
%  INPUT  : group_feature_table_combined.mat (from S4) + per-subject
%           Ox##_outcome_trimmed.set epoch files
%  OUTPUT : group_feature_table_combined_wavelet.mat, grand_tf.mat,
%           tf_perm_stats.mat, and PDFs in <Figures>/TF_analysis/

% fprintf('─── S10: Wavelet time-frequency analysis ───────────────────────\n');
% run(fullfile(fileparts(mfilename('fullpath')), 'S10_wavelet_TF_analysis.m'));
% fprintf('S10 complete.\n\n');

%% ── DONE ────────────────────────────────────────────────────────────────────
fprintf('================================================================\n');
fprintf('  PIPELINE RUN COMPLETE\n');
fprintf('  Finished  : %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('================================================================\n');
