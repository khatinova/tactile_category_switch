% =============================================================================
% S3_extract_eeg_features.m
%
% PIPELINE STEP 3 of 7 — EXTEND the per-trial EEG feature table built by S2.
%
% S3 does NOT re-extract features and does NOT reload raw/epoched .set files.
% S2 already loaded each subject's outcome epochs ONCE and produced the full
% per-trial feature table (N2, FCzCz mean, P300, Theta, PLV) together with the
% per-trial waveforms, the per-stage summary, and the per-stage FRN/RewP
% difference-wave table. S3 simply LOADS that S2 output and ADDS the two things
% S2 leaves out:
%
%   (1) WITHIN-SUBJECT Z-SCORES of the canonical single-trial features
%       (the *_z columns the LME models in S7 expect).
%   (2) GROUP-LEVEL GRAND-AVERAGE FIGURES:
%         - correct vs incorrect frontocentral ERP, raw and baseline-RMS
%           normalised (built from the stored FCzCz_waveform + baseline_rms);
%         - FRN/RewP difference waves by block type (from frn_rewp_stage_table).
%
% This removes the old duplication: S3 previously reloaded a DIFFERENT spine
% (group_stage_table_mICA.mat) and recomputed the same features with a
% mismatched input. It now consumes S2's output directly.
%
% INPUT  (per cohort, from S2):
%   group_feature_table_combined_<COHORT>.mat
%       -> all_trials_table, stage_feature_table, frn_rewp_stage_table, t_ax
% OUTPUT (per cohort, for S4):
%   group_feature_table_<COHORT>_final.mat
%       -> all_trials_table (now with *_z columns), stage_feature_table,
%          frn_rewp_stage_table, t_ax
%   + grand-average figures (PDF/PNG) under Figures/RQ_analysis.
%
% Run once per cohort: set COHORT = 'KH' (this file) or use S3_RR for 'RR'.
% =============================================================================

clear; close all; clc;
addpath(genpath(fileparts(mfilename('fullpath'))));

% -------------------------------------------------------------------------
%% COHORT + PATHS
% -------------------------------------------------------------------------
COHORT = 'KH';        % 'KH' here; S3_RR sets 'RR'. Everything else is identical.
tag    = upper(COHORT);


remote = 0;
if remote == 1
    base_path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch';
else
    base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch';
end

if strcmpi(COHORT, 'RR')
    cohort_results = fullfile(base_path, 'Salient mod switch RR', 'Results', 'EEG analysis');
else
    cohort_results = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis');
end

feature_folder = fullfile(cohort_results, 'Outcome_feature_tables_v4_merged');
figure_folder  = fullfile(cohort_results, 'Figures', 'RQ_analysis');
if ~exist(figure_folder, 'dir'), mkdir(figure_folder); end

% -------------------------------------------------------------------------
%% LOAD S2 OUTPUT (the single source of truth for features)
% -------------------------------------------------------------------------
combined_file = fullfile(feature_folder, sprintf('group_feature_table_combined_%s.mat', tag));
if ~exist(combined_file, 'file')
    error(['S3: S2 output not found:\n  %s\n' ...
           'Run S2 (COHORT=%s) first so the per-trial feature table exists.'], ...
           combined_file, tag);
end

S2out = load(combined_file);
all_trials_table    = S2out.all_trials_table;
t_ax                = S2out.t_ax;
if isfield(S2out, 'stage_feature_table');  stage_feature_table  = S2out.stage_feature_table;  else, stage_feature_table  = table(); end
if isfield(S2out, 'frn_rewp_stage_table'); frn_rewp_stage_table = S2out.frn_rewp_stage_table; else, frn_rewp_stage_table = table(); end

fprintf('S3 (%s): loaded S2 table — %d trials, %d subjects.\n', tag, ...
    height(all_trials_table), numel(unique(string(all_trials_table.subj_id))));

% -------------------------------------------------------------------------
%% SETTINGS
% -------------------------------------------------------------------------
ERP_plot_window = [-200 800];
FRN_win  = [250 300];   % shaded on the difference-wave figure (matches S2)
RewP_win = [250 350];

% Canonical single-trial features to z-score within subject (only those that
% actually exist in the S2 table are used).
features_to_zscore = {'N2_amp','N2_norm','FCzCz_mean_amp','FCzCz_mean_norm', ...
    'P300_amp','P300_norm','Theta_amp', ...
    'PLV_fp','PLV_fs','PLV_fp_pairwise','PLV_fs_pairwise'};


% -------------------------------------------------------------------------
%% STEP 1: WITHIN-SUBJECT Z-SCORES  (the *_z columns S7's LMEs use)
% -------------------------------------------------------------------------
% Unify subject naming first so z-scoring groups by the canonical subj_id.
all_trials_table = kh_subject_id('standardise', all_trials_table);

subj_list = unique(string(all_trials_table.subj_id));
n_rows = height(all_trials_table);

for f = 1:numel(features_to_zscore)
    fn = features_to_zscore{f};
    if ~ismember(fn, all_trials_table.Properties.VariableNames); continue; end

    fn_z = [fn '_z'];
    all_trials_table.(fn_z) = nan(n_rows, 1);

    for si = 1:numel(subj_list)
        mask = string(all_trials_table.subj_id) == subj_list(si);
        vals = all_trials_table.(fn)(mask);
        sd   = std(vals, 'omitnan');
        if sd > 0
            all_trials_table.(fn_z)(mask) = (vals - mean(vals, 'omitnan')) / sd;
        end
    end
end
fprintf('  Added within-subject z-scores for %d features.\n', ...
    sum(ismember(features_to_zscore, all_trials_table.Properties.VariableNames)));

% -------------------------------------------------------------------------
%% SAVE the extended table for S4 (per-trial table now carries *_z columns)
% -------------------------------------------------------------------------
final_file = fullfile(feature_folder, sprintf('group_feature_table_%s_final.mat', tag));
save(final_file, 'all_trials_table', 'stage_feature_table', ...
     'frn_rewp_stage_table', 't_ax', '-v7.3');
fprintf('  Saved extended table: %s\n', final_file);


% -------------------------------------------------------------------------
%% STEP 2a: GRAND-AVERAGE ERP FIGURE (raw vs baseline-RMS normalised)
% Built from the per-trial FCzCz_waveform stored by S2 — no .set reload.
% -------------------------------------------------------------------------
T = all_trials_table;
have_wave = ismember('FCzCz_waveform', T.Properties.VariableNames);

raw_correct = []; raw_incorrect = [];
norm_correct = []; norm_incorrect = [];

if have_wave
    for si = 1:numel(subj_list)
        sm = string(T.subj_id) == subj_list(si) & T.has_eeg_epoch & ...
             ~T.epoch_artifact_flag & ~T.false_fb;
        if ~any(sm); continue; end

        % Subject baseline RMS (constant per subject; stored per trial by S2).
        if ismember('baseline_rms', T.Properties.VariableNames)
            brms = mean(T.baseline_rms(sm), 'omitnan');
        else
            brms = NaN;
        end

        erp_c = mean_waveform(T.FCzCz_waveform(sm & T.correct_num == 1));
        erp_i = mean_waveform(T.FCzCz_waveform(sm & T.correct_num == 0));

        if ~isempty(erp_c)
            raw_correct(end+1,:)  = erp_c;                       %#ok<AGROW>
            norm_correct(end+1,:) = erp_c / max(brms, eps);      %#ok<AGROW>
        end
        if ~isempty(erp_i)
            raw_incorrect(end+1,:)  = erp_i;                     %#ok<AGROW>
            norm_incorrect(end+1,:) = erp_i / max(brms, eps);    %#ok<AGROW>
        end
    end
end

if ~isempty(raw_correct) || ~isempty(raw_incorrect)
    in_w = t_ax >= ERP_plot_window(1) & t_ax <= ERP_plot_window(2);
    fig_norm = figure('Position', [50 50 1200 500]);
    sgtitle(sprintf('%s: frontocentral ERP (raw vs baseline-RMS normalised)', tag));

    subplot(1,2,1); hold on; title('Raw (\muV)');
    shaded_erp(t_ax(in_w), pick(raw_correct, in_w),   [0.10 0.60 0.10], 'Correct');
    shaded_erp(t_ax(in_w), pick(raw_incorrect, in_w), [0.70 0.10 0.10], 'Incorrect');
    xline(0,'k--'); yline(0,'k:'); set(gca,'YDir','reverse');
    xlabel('Time (ms)'); ylabel('Amplitude (\muV)'); legend('Box','off');

    subplot(1,2,2); hold on; title('Normalised (\muV / baseline RMS)');
    shaded_erp(t_ax(in_w), pick(norm_correct, in_w),   [0.10 0.60 0.10], 'Correct');
    shaded_erp(t_ax(in_w), pick(norm_incorrect, in_w), [0.70 0.10 0.10], 'Incorrect');
    xline(0,'k--'); yline(0,'k:'); set(gca,'YDir','reverse');
    xlabel('Time (ms)'); ylabel('Amplitude (normalised)'); legend('Box','off');

    apply_fig_style(fig_norm);
    exportgraphics(fig_norm, fullfile(figure_folder, ...
        sprintf('Normalisation_comparison_%s.pdf', tag)), 'ContentType', 'vector');
end


% -------------------------------------------------------------------------
%% STEP 2b: FRN/RewP DIFFERENCE-WAVE GRAND AVERAGES (by block type)
% Built from frn_rewp_stage_table.diff_wave (incorrect - correct) made by S2.
% -------------------------------------------------------------------------
if ~isempty(frn_rewp_stage_table) && ismember('diff_wave', frn_rewp_stage_table.Properties.VariableNames)
    bts = {'D','P'};
    fig_fr = figure('Position', [60 60 1100 460]);
    sgtitle(sprintf('%s: FRN/RewP difference waves (incorrect - correct), by block type', tag));

    for bi = 1:numel(bts)
        ax = subplot(1, numel(bts), bi); hold(ax, 'on'); title(ax, sprintf('%s blocks', bts{bi}));

        sel = string(frn_rewp_stage_table.block_type) == bts{bi};
        dw  = frn_rewp_stage_table.diff_wave(sel);
        dw  = dw(~cellfun(@isempty, dw));
        if isempty(dw); continue; end

        M  = cell2mat(cellfun(@(x) x(:)', dw, 'UniformOutput', false));
        mn = mean(M, 1, 'omitnan');
        se = std(M, 0, 1, 'omitnan') / sqrt(size(M,1));
        yl = [min(mn)-0.5, max(mn)+0.5];

        patch(ax, [FRN_win fliplr(FRN_win)],  [yl(1) yl(1) yl(2) yl(2)], [0.85 0.90 1.00], 'EdgeColor','none','FaceAlpha',0.4);
        patch(ax, [RewP_win fliplr(RewP_win)],[yl(1) yl(1) yl(2) yl(2)], [1.00 0.90 0.85], 'EdgeColor','none','FaceAlpha',0.3);
        fill(ax, [t_ax fliplr(t_ax)], [mn+se fliplr(mn-se)], [0.2 0.2 0.2], 'FaceAlpha',0.15, 'EdgeColor','none');
        plot(ax, t_ax, mn, 'k', 'LineWidth', 2);
        yline(ax,0,'k:'); xline(ax,0,'k:');
        xlabel(ax,'Time (ms)'); ylabel(ax,'incorrect - correct (\muV)');
        xlim(ax, ERP_plot_window); set(ax,'YDir','reverse');   % EEG negative-up
    end

    apply_fig_style(fig_fr);
    exportgraphics(fig_fr, fullfile(figure_folder, ...
        sprintf('FRN_RewP_difference_waves_%s.pdf', tag)), 'ContentType', 'vector');
end

fprintf('S3 (%s) complete: extended table + grand-average figures saved.\n', tag);

% =========================================================================
%% LOCAL FUNCTIONS
% =========================================================================
function m = mean_waveform(wave_cells)
% Grand-average a cell column of equal-length per-trial waveforms.
wave_cells = wave_cells(~cellfun(@isempty, wave_cells));
if isempty(wave_cells); m = []; return; end
M = cell2mat(cellfun(@(x) x(:)', wave_cells, 'UniformOutput', false));
m = mean(M, 1, 'omitnan');
end

function out = pick(M, cols)
% Safely column-subset a (possibly empty) grand-average matrix.
if isempty(M); out = M; else, out = M(:, cols); end
end

function shaded_erp(t, data, clr, lbl)
% Mean +/- SEM ribbon for a [subjects x time] matrix.
if isempty(data); return; end
mn = mean(data, 1, 'omitnan');
se = std(data, 0, 1, 'omitnan') / sqrt(size(data,1));
fill([t, fliplr(t)], [mn+se, fliplr(mn-se)], clr, ...
    'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(t, mn, 'Color', clr, 'LineWidth', 2.5, ...
    'DisplayName', sprintf('%s (n=%d)', lbl, size(data,1)));
end
