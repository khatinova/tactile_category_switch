% =============================================================================
% OUTCOME ERP ANALYSIS v9 — COMBINED KH + RR FEATURE EXTRACTION TEMPLATE
%
% PURPOSE
% -------
% One shared extraction pipeline for KH (Ox) and RR (Nc) cohorts.
% Cohort differences are handled only through a cohort configuration table:
%   - paths
%   - subject prefix
%   - participant IDs
%   - channel labels
%   - alignment mode
%   - output filenames
%
% The actual ERP feature extraction is identical across cohorts.
%
% IMPORTANT
% ---------
% This is designed to replace separate KH/RR feature extraction scripts.
% Run this script after checking all helper functions are on the MATLAB path:
%   define_trial_stages_v3.m
%   validate_stage_table.m
%   KH_align_epochs_with_offset.m
%   permutation_test_cluster_correction.m
%
% RR-specific assumptions retained:
%   - Direct sequential alignment: epoch k = trial k
%   - Subject prefix Nc%02d
%   - EGI channel labels, e.g. E11 for FCz, E7 for Cz
%   - false_fb derived from beh.trueFB
%   - fb_shown_correct derived from EEG triggers rewa/puni
%
% KH-specific assumptions retained:
%   - Offset alignment via KH_align_epochs_with_offset
%   - Subject prefix Ox%02d
%   - 10-20 channel labels, e.g. FCz/Cz
%   - trimmed file preferred when present
%
% Primary RewP/FRN scalar features are MEAN amplitudes in a priori windows.
% Robust negative peak measures are diagnostic/sensitivity features only.
% =============================================================================

clear; close all;

% -------------------------------------------------------------------------
%% PATHS
% -------------------------------------------------------------------------
remote = 0;
if remote == 1
    base_path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch';
else
    base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch';
end

fieldtrip_path = 'C:\Users\khatinova\OneDrive - Nexus365\Pre_2026_Folders\Documents\MATLAB\fieldtrip-20240110';
eeglab_path    = 'C:\Users\khatinova\OneDrive - Nexus365\Pre_2026_Folders\Documents\MATLAB\eeglab2025.1.0';

addpath(fieldtrip_path); ft_defaults;
addpath(eeglab_path);    eeglab nogui;

KH_data_path = fullfile(base_path, 'Salient mod switch KH', 'Data');
RR_data_path = fullfile(base_path, 'Salient mod switch RR');
addpath(genpath(KH_data_path));
addpath(genpath(RR_data_path));

% This appears to contain both Ox and Nc trial_data structs in your current RR script.
load(fullfile(KH_data_path, 'all_trial_data.mat'));

% KH behav_table is used as the metadata spine if available.
KH_behav_file = fullfile(KH_data_path, 'behav_table.mat');
if exist(KH_behav_file, 'file')
    S = load(KH_behav_file, 'group_T');
    KH_behav_table = S.group_T;
else
    KH_behav_table = table();
end

% -------------------------------------------------------------------------
%% GLOBAL PARAMETERS
% -------------------------------------------------------------------------
save_tables = true;
ERP_plot_window = [-200 1000];
rm_baseline     = [-200 0];

% A priori windows
N2_win       = [120 350];
FRN_win      = [250 350];   % primary robust feedback window
RewP_win     = [250 350];   % same window; interpreted via outcome contrast
P300_win     = [300 600];
Theta_win    = [200 500];
PLV_win      = [200 400];
PLV_baseline = [-200 0];

% Diagnostic robust peak detection
NEG_PEAK_SMOOTH_SAMPLES = 5;      % 10 ms at 500 Hz
NEG_PEAK_MIN_PROM_RMS   = 0.50;   % prominence threshold relative to baseline RMS
NEG_PEAK_EDGE_MARGIN_MS = 20;

% PLV settings
MIN_TRIALS_PLV = 5;
PLV_WINDOW_HALF = 7;
MIN_TRIALS_PLV_WINDOW = 5;

stage_names  = {'LN','LE','RN','RE'};
BTYPE_LABELS = {'D','P'};
STAGE_COLORS = [0.12 0.62 0.47; 0.85 0.65 0.00; 0.80 0.27 0.13; 0.40 0.25 0.65];
LINE_STYLES  = {'-','--'};

% -------------------------------------------------------------------------
%% COHORT CONFIGURATION
% -------------------------------------------------------------------------
cohorts = struct([]);

cohorts(1).name = 'KH';
cohorts(1).prefix = 'Ox';
cohorts(1).valid_participants = [3:8, 10:12, 14:23, 27];
cohorts(1).data_path = KH_data_path;
cohorts(1).epoch_file_folder = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis', 'Epoched_data');
cohorts(1).figure_output_folder = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis', 'Figures', 'outcome_v9_KH');
cohorts(1).alignment_mode = 'offset';
cohorts(1).trimmed_preferred = true;
cohorts(1).behav_table = KH_behav_table;
cohorts(1).fcz_label = 'FCz';
cohorts(1).cz_label = 'Cz';
cohorts(1).frontocentral_channels = {'FCz','Cz','F1','F2'};
cohorts(1).acc_channels = {'FCz','Fz','AFz','F1','F2'};
cohorts(1).par_channels = {'Pz','P1','P2'};
cohorts(1).som_channels = {'C3','C4','CP3','CP1','C5','CP5'};
cohorts(1).handoff_file = 'group_stage_table_features_KH_v9.mat';
cohorts(1).full_file = 'group_table_all_trials_KH_v9.mat';
cohorts(1).grand_file = 'grand_KH_v9.mat';

cohorts(2).name = 'RR';
cohorts(2).prefix = 'Nc';
cohorts(2).valid_participants = 1:15;
cohorts(2).data_path = RR_data_path;
cohorts(2).epoch_file_folder = fullfile(base_path, 'Salient mod switch RR', 'Results', 'EEG_analysis', 'Epoched_data');
cohorts(2).figure_output_folder = fullfile(base_path, 'Salient mod switch RR', 'Results', 'EEG_analysis', 'Outcome_ERP_analysis_figures_v9');
cohorts(2).alignment_mode = 'direct';
cohorts(2).trimmed_preferred = false;
cohorts(2).behav_table = table();
cohorts(2).fcz_label = 'E11';
cohorts(2).cz_label = 'E7';           % verify montage
cohorts(2).frontocentral_channels = {'E11','E7','E6'};
cohorts(2).acc_channels = {'E11','E6','E16'};
cohorts(2).par_channels = {'E62','E67','E72'};
cohorts(2).som_channels = {'E36','E104','E41','E103'};
cohorts(2).handoff_file = 'group_stage_table_features_RR_v9.mat';
cohorts(2).full_file = 'group_table_all_trials_RR_v9.mat';
cohorts(2).grand_file = 'grand_RR_v9.mat';

% -------------------------------------------------------------------------
%% RUN EACH COHORT THROUGH SAME EXTRACTION ENGINE
% -------------------------------------------------------------------------
all_handoff_tables = {};
all_debug_tables = {};

for ci = 1:numel(cohorts)
    cfg = cohorts(ci);
    fprintf('\n\n############################################################\n');
    fprintf('RUNNING COHORT: %s\n', cfg.name);
    fprintf('############################################################\n');

    if ~exist(cfg.figure_output_folder, 'dir'), mkdir(cfg.figure_output_folder); end
    if ~exist(cfg.epoch_file_folder, 'dir')
        error('Epoch folder not found for %s: %s', cfg.name, cfg.epoch_file_folder);
    end

    [group_table, all_trials_table, grand, t_ax, debug_table] = run_one_cohort_extraction( ...
        cfg, all_trial_data, rm_baseline, N2_win, FRN_win, RewP_win, P300_win, Theta_win, ...
        PLV_win, PLV_baseline, MIN_TRIALS_PLV, PLV_WINDOW_HALF, MIN_TRIALS_PLV_WINDOW, ...
        NEG_PEAK_SMOOTH_SAMPLES, NEG_PEAK_MIN_PROM_RMS, NEG_PEAK_EDGE_MARGIN_MS, ...
        ERP_plot_window, stage_names, BTYPE_LABELS, STAGE_COLORS, LINE_STYLES);

    if save_tables
        save(fullfile(cfg.epoch_file_folder, cfg.handoff_file), 'group_table');
        save(fullfile(cfg.epoch_file_folder, cfg.full_file), 'all_trials_table', 't_ax');
        save(fullfile(cfg.epoch_file_folder, cfg.grand_file), 'grand', 't_ax');
        writetable(debug_table, fullfile(cfg.epoch_file_folder, sprintf('debug_%s_v9.csv', cfg.name)));
    end

    all_handoff_tables{end+1} = group_table;
    all_debug_tables{end+1} = debug_table;
end

% -------------------------------------------------------------------------
%% COMBINE KH + RR HANDOFF TABLES
% -------------------------------------------------------------------------
fprintf('\nCombining KH + RR handoff tables...\n');
[group_table_combined, merge_debug] = combine_handoff_tables_v9(all_handoff_tables);

combined_out_folder = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis', 'Epoched_data');
save(fullfile(combined_out_folder, 'group_feature_table_combined_v9.mat'), 'group_table_combined');
writetable(merge_debug, fullfile(combined_out_folder, 'merge_debug_combined_v9.csv'));

fprintf('Saved combined table: %s\n', fullfile(combined_out_folder, 'group_feature_table_combined_v9.mat'));
fprintf('Combined rows: %d, columns: %d\n', height(group_table_combined), width(group_table_combined));

% =============================================================================
% LOCAL FUNCTIONS
% =============================================================================

function [group_table, all_trials_table, grand, t_ax, debug_table] = run_one_cohort_extraction( ...
    cfg, all_trial_data, rm_baseline, N2_win, FRN_win, RewP_win, P300_win, Theta_win, ...
    PLV_win, PLV_baseline, MIN_TRIALS_PLV, PLV_WINDOW_HALF, MIN_TRIALS_PLV_WINDOW, ...
    NEG_PEAK_SMOOTH_SAMPLES, NEG_PEAK_MIN_PROM_RMS, NEG_PEAK_EDGE_MARGIN_MS, ...
    ERP_plot_window, stage_names, BTYPE_LABELS, STAGE_COLORS, LINE_STYLES)

empty_container = struct('data', [], 'subj', []);
for s = 1:4
    for bt = 1:2
        bt_s = BTYPE_LABELS{bt};
        grand.FCz.(stage_names{s}).(bt_s).correct   = empty_container;
        grand.FCz.(stage_names{s}).(bt_s).incorrect = empty_container;
        grand.FCz.(stage_names{s}).(bt_s).false_cor = empty_container;
        grand.FCz.(stage_names{s}).(bt_s).false_inc = empty_container;
        grand.Par.(stage_names{s}).(bt_s).correct   = empty_container;
        grand.Par.(stage_names{s}).(bt_s).incorrect = empty_container;
        grand.Theta.(stage_names{s}).(bt_s).correct = empty_container;
        grand.Theta.(stage_names{s}).(bt_s).incorrect = empty_container;
        grand.PLV_fp.(stage_names{s}).(bt_s) = empty_container;
        grand.PLV_fs.(stage_names{s}).(bt_s) = empty_container;
    end
end

all_trials_table = table();
t_ax = [];
debug_rows = {};

for participant = cfg.valid_participants
    subj = sprintf('%s%02d', cfg.prefix, participant);
    fprintf('\n============ %s ============\n', subj);

    if ~isfield(all_trial_data, subj)
        warning('%s missing from all_trial_data. Skipping.', subj);
        continue;
    end

    beh = all_trial_data.(subj).trial_data;
    if strcmp(cfg.name, 'RR') && isfield(beh, 'structCode')
        beh.block_structure = beh.structCode;
    end

    % Flatten behavioural vectors, with KH 6-block practice handling.
    [beh, beh_correct, beh_conf, beh_trueFB, total_trials] = flatten_behaviour_for_cohort(beh, cfg);

    % Load EEG files.
    EEGp = load_eeg_for_subject(cfg, subj, 'broadband');
    if isempty(EEGp), continue; end
    if isempty(t_ax), t_ax = EEGp.times; end

    EEGp_theta = load_eeg_for_subject(cfg, subj, 'theta');
    EEGp_phase = load_eeg_for_subject(cfg, subj, 'phase');

    total_trials = min(total_trials, EEGp.trials);

    % Alignment.
    if strcmp(cfg.alignment_mode, 'direct')
        trial2epoch = (1:total_trials)';
        fprintf('  Alignment: direct sequential (%d trials)\n', total_trials);
    else
        trial2epoch = align_kh_epochs(cfg, subj, EEGp, beh_correct, total_trials);
    end

    % Stage table.
    stage_table = define_trial_stages_v3(beh, trial2epoch, EEGp, participant, EEGp_theta, EEGp_phase);

    if strcmp(cfg.alignment_mode, 'direct')
        [stage_table, fb_shown_correct_vec, false_fb_vec] = add_rr_feedback_columns(stage_table, EEGp, beh_trueFB, total_trials);
    else
        stage_table = ensure_kh_feedback_columns(stage_table);
    end

    try
        validate_stage_table(stage_table, EEGp, beh, participant);
    catch ME
        warning('%s validate_stage_table warning/failure: %s', subj, ME.message);
    end

    stage_table.confidence = nan(height(stage_table), 1);
    conf_vec = beh_conf(:);
    for r = 1:height(stage_table)
        tc = stage_table.trial_continuous(r);
        if tc >= 1 && tc <= numel(conf_vec)
            stage_table.confidence(r) = conf_vec(tc);
        end
    end
    stage_table.subj = repmat(participant, height(stage_table), 1);
    stage_table.subj_id = repmat(string(subj), height(stage_table), 1);
    stage_table.cohort = repmat(string(cfg.name), height(stage_table), 1);
    stage_table.is_cohort1 = repmat(strcmp(cfg.name, 'KH') && participant <= 8, height(stage_table), 1);

    % Channel indices and masks.
    fcz_idx = safe_chan(EEGp, {cfg.fcz_label});
    cz_idx  = safe_chan(EEGp, {cfg.cz_label});
    par_idx = safe_chan(EEGp, cfg.par_channels);
    acc_idx = safe_chan(EEGp, cfg.acc_channels);
    som_idx = safe_chan(EEGp, cfg.som_channels);

    bl_mask   = EEGp.times >= rm_baseline(1) & EEGp.times <= rm_baseline(2);
    n2_mask   = EEGp.times >= N2_win(1)      & EEGp.times <= N2_win(2);
    frn_mask  = EEGp.times >= FRN_win(1)     & EEGp.times <= FRN_win(2);
    rewp_mask = EEGp.times >= RewP_win(1)    & EEGp.times <= RewP_win(2);
    p300_mask = EEGp.times >= P300_win(1)    & EEGp.times <= P300_win(2);
    th_mask   = EEGp.times >= Theta_win(1)   & EEGp.times <= Theta_win(2);
    plv_mask  = EEGp.times >= PLV_win(1)     & EEGp.times <= PLV_win(2);
    plv_bl    = EEGp.times >= PLV_baseline(1) & EEGp.times <= PLV_baseline(2);

    if ~isempty(fcz_idx)
        bl_data = squeeze(double(EEGp.data(fcz_idx, bl_mask, :)));
        bline_rms = rms(bl_data(:), 'omitnan');
    else
        bline_rms = NaN;
    end

    % Build unified feature table for this subject.
    subj_features = build_subject_feature_spine(cfg, subj, participant, beh, stage_table, trial2epoch, total_trials);

    subj_features = initialise_feature_columns(subj_features);

    % Extract ERP features identically for both cohorts.
    for ti = 1:height(subj_features)
        ep = subj_features.epoch(ti);
        if isnan(ep) || ep < 1 || ep > EEGp.trials, continue; end
        ep = round(ep);

        if ~isempty(fcz_idx) && ~isempty(cz_idx) && ~isnan(bline_rms)
            sig = mean(double(EEGp.data([fcz_idx cz_idx], :, ep)), 1, 'omitnan');
            sig = sig - mean(sig(bl_mask), 'omitnan');
            subj_features.FCzCz_waveform{ti} = sig;
            subj_features.FCzCz_signal{ti} = sig; % compatibility alias

            % Primary mean-amplitude features.
            subj_features.FCzCz_mean_250_350(ti) = mean(sig(frn_mask), 'omitnan');
            subj_features.FCzCz_mean_250_350_norm(ti) = subj_features.FCzCz_mean_250_350(ti) / bline_rms;
            subj_features.FRN_mean_amp(ti) = subj_features.FCzCz_mean_250_350(ti);
            subj_features.FRN_mean_norm(ti) = subj_features.FCzCz_mean_250_350_norm(ti);
            subj_features.RewP_mean_amp(ti) = mean(sig(rewp_mask), 'omitnan');
            subj_features.RewP_mean_norm(ti) = subj_features.RewP_mean_amp(ti) / bline_rms;

            % Backward-compatible aliases now use mean amplitude, not peaks.
            subj_features.FRN_amp(ti) = subj_features.FRN_mean_amp(ti);
            subj_features.FRN_norm(ti) = subj_features.FRN_mean_norm(ti);
            subj_features.RewP_amp(ti) = subj_features.RewP_mean_amp(ti);
            subj_features.RewP_norm(ti) = subj_features.RewP_mean_norm(ti);

            % Legacy/diagnostic N2 minimum.
            win_vals = sig(n2_mask); win_t = EEGp.times(n2_mask);
            if any(~isnan(win_vals))
                [pk, ix] = min(win_vals, [], 'omitnan');
                subj_features.N2_amp(ti) = pk;
                subj_features.N2_lat(ti) = win_t(ix);
                subj_features.N2_norm(ti) = pk / bline_rms;
            end
        end

        % Diagnostic robust negative peak at FCz.
        if ~isempty(fcz_idx) && ~isnan(bline_rms)
            sig_fcz = double(EEGp.data(fcz_idx, :, ep));
            sig_fcz = sig_fcz - mean(sig_fcz(bl_mask), 'omitnan');
            [valid_pk, pk_amp, pk_lat, pk_prom] = robust_negative_peak( ...
                sig_fcz, EEGp.times, frn_mask, bline_rms, NEG_PEAK_SMOOTH_SAMPLES, ...
                NEG_PEAK_MIN_PROM_RMS, NEG_PEAK_EDGE_MARGIN_MS, FRN_win);
            subj_features.FCz_neg_peak_valid(ti) = valid_pk;
            subj_features.FCz_neg_peak_amp_robust(ti) = pk_amp;
            subj_features.FCz_neg_peak_lat_robust(ti) = pk_lat;
            subj_features.FCz_neg_peak_prom(ti) = pk_prom;

            % Backward-compatible diagnostic aliases.
            subj_features.FCz_neg_peak_amp(ti) = pk_amp;
            subj_features.FCz_neg_peak_lat(ti) = pk_lat;
            if ~isnan(pk_amp), subj_features.FCz_neg_peak_norm(ti) = pk_amp / bline_rms; end
        end

        % P300 peak amplitude.
        if ~isempty(par_idx) && ~isnan(bline_rms)
            sig_p = mean(double(EEGp.data(par_idx, :, ep)), 1, 'omitnan');
            sig_p = sig_p - mean(sig_p(bl_mask), 'omitnan');
            subj_features.P300_waveform{ti} = sig_p;
            subj_features.P300_signal{ti} = sig_p;
            win_vals = sig_p(p300_mask); win_t = EEGp.times(p300_mask);
            if any(~isnan(win_vals))
                [pk, ix] = max(win_vals, [], 'omitnan');
                subj_features.P300_amp(ti) = pk;
                subj_features.P300_norm(ti) = pk / bline_rms;
                subj_features.P300_peak_lat(ti) = win_t(ix);
            end
        end

        % Theta envelope.
        if ~isempty(EEGp_theta) && ~isempty(acc_idx) && ep <= EEGp_theta.trials
            sig_th = mean(double(EEGp_theta.data(acc_idx, :, ep)), 1, 'omitnan');
            env = abs(hilbert(sig_th));
            env = env - mean(env(bl_mask), 'omitnan');
            subj_features.Theta_amp(ti) = mean(env(th_mask), 'omitnan');
            subj_features.Theta_waveform{ti} = env;
            subj_features.Theta_signal{ti} = env;
        end
    end

    % Sliding-window local PLV scalar columns.
    subj_features = assign_sliding_window_plv(subj_features, trial2epoch, EEGp_phase, acc_idx, par_idx, som_idx, ...
        plv_mask, plv_bl, PLV_WINDOW_HALF, MIN_TRIALS_PLV_WINDOW);

    all_trials_table = [all_trials_table; subj_features]; %#ok<AGROW>

    % Grand averages using stage_table, identical logic across cohorts.
    grand = accumulate_grand_averages(grand, cfg, participant, stage_table, EEGp, EEGp_theta, EEGp_phase, ...
        fcz_idx, par_idx, acc_idx, som_idx, bl_mask, plv_mask, plv_bl, MIN_TRIALS_PLV, ...
        stage_names, BTYPE_LABELS);

    has_P = any(stage_table.block_type == 'P');
    plot_subject_erps_unified(cfg, subj, participant, has_P, EEGp, EEGp_theta, stage_table, ...
        fcz_idx, par_idx, acc_idx, bl_mask, stage_names, STAGE_COLORS, LINE_STYLES, BTYPE_LABELS, ERP_plot_window);

    debug_rows(end+1,:) = {string(cfg.name), string(subj), participant, height(subj_features), ...
        sum(subj_features.has_eeg_epoch), sum(~cellfun(@isempty, subj_features.FCzCz_waveform)), ...
        sum(~isnan(subj_features.FRN_amp)), sum(subj_features.FCz_neg_peak_valid), ...
        mean(subj_features.FCz_neg_peak_valid, 'omitnan'), bline_rms}; %#ok<AGROW>

    clear EEGp EEGp_theta EEGp_phase
end

group_table = all_trials_table;
group_table = finalise_group_table_types_and_zscore(group_table);

debug_table = cell2table(debug_rows, 'VariableNames', {'cohort','subj_id','subj','n_rows', ...
    'n_has_epoch','n_FCzCz_waveform','n_FRN_amp','n_valid_neg_peak','prop_valid_neg_peak','baseline_rms'});
end

function [beh, beh_correct, beh_conf, beh_trueFB, total_trials] = flatten_behaviour_for_cohort(beh, cfg)
num_blocks = height(beh.correct);
beh_correct = [];
beh_conf = [];
beh_trueFB = [];

if strcmp(cfg.name, 'KH') && num_blocks >= 6
    start_block = 2;
    beh.correct = beh.correct(2:end,:);
    beh.confidence = beh.confidence(2:end,:);
    if isfield(beh, 'trueFB'), beh.trueFB = beh.trueFB(2:end,:); end
    if isfield(beh, 'revTrial'), beh.revTrial = beh.revTrial(2:end); end
else
    start_block = 1;
end

num_blocks = height(beh.correct);
for b = 1:num_blocks
    beh_correct = [beh_correct, beh.correct(b,:)]; %#ok<AGROW>
    if isfield(beh, 'confidence'), beh_conf = [beh_conf, beh.confidence(b,:)]; end %#ok<AGROW>
    if isfield(beh, 'trueFB'), beh_trueFB = [beh_trueFB, beh.trueFB(b,:)]; end %#ok<AGROW>
end
if isempty(beh_trueFB), beh_trueFB = ones(size(beh_correct)); end
if isempty(beh_conf), beh_conf = nan(size(beh_correct)); end

total_trials = numel(beh_correct);
end

function EEGp = load_eeg_for_subject(cfg, subj, kind)
EEGp = [];
switch kind
    case 'broadband'
        candidates = {sprintf('%s_outcome_trimmed.set', subj), sprintf('%s_outcome.set', subj)};
    case 'theta'
        candidates = {sprintf('%s_outcome_theta_trimmed.set', subj), sprintf('%s_outcome_theta.set', subj)};
    case 'phase'
        candidates = {sprintf('%s_outcome_phase_trimmed.set', subj), sprintf('%s_outcome_phase.set', subj)};
end
if ~cfg.trimmed_preferred
    candidates = fliplr(candidates);
end
for i = 1:numel(candidates)
    f = fullfile(cfg.epoch_file_folder, candidates{i});
    if exist(f, 'file')
        EEGp = pop_loadset(candidates{i}, cfg.epoch_file_folder);
        fprintf('  %s: %s (%d epochs)\n', kind, candidates{i}, EEGp.trials);
        return;
    end
end
if strcmp(kind, 'broadband')
    warning('%s file missing for %s.', kind, subj);
else
    warning('%s file missing for %s; features depending on it will be NaN.', kind, subj);
end
end

function trial2epoch = align_kh_epochs(cfg, subj, EEGp, beh_correct, total_trials)
trial2epoch_file = fullfile(cfg.epoch_file_folder, sprintf('%s_trial2epoch.mat', subj));
if exist(trial2epoch_file, 'file')
    S = load(trial2epoch_file);
    if isfield(S, 'trial2epoch_out')
        trial2epoch = S.trial2epoch_out;
    elseif isfield(S, 'trial2epoch')
        trial2epoch = S.trial2epoch;
    else
        error('%s has unexpected field names.', trial2epoch_file);
    end
else
    beh_cv = beh_correct(:);
    beh_cv = beh_cv(~isnan(beh_cv));
    [trial2epoch, diag] = KH_align_epochs_with_offset(EEGp, beh_cv);
    fprintf('  trial2epoch: %d/%d matched (%.1f%%), offset=%d\n', ...
        diag.n_matched, diag.n_trials, diag.match_rate*100, diag.best_offset);
    trial2epoch_out = trial2epoch;
    save(trial2epoch_file, 'trial2epoch_out');
end
trial2epoch = trial2epoch(:);
if numel(trial2epoch) < total_trials
    trial2epoch(end+1:total_trials) = NaN;
elseif numel(trial2epoch) > total_trials
    trial2epoch = trial2epoch(1:total_trials);
end
end

function [stage_table, fb_shown_correct_vec, false_fb_vec] = add_rr_feedback_columns(stage_table, EEGp, beh_trueFB, total_trials)
fb_shown_correct_vec = false(total_trials, 1);
for k = 1:total_trials
    ep_events = {EEGp.event([EEGp.event.epoch] == k).type};
    if any(strcmp(ep_events, 'rewa'))
        fb_shown_correct_vec(k) = true;
    end
end
false_fb_vec = ~logical(beh_trueFB(1:total_trials)');
stage_table.fb_shown_correct = false(height(stage_table), 1);
stage_table.false_fb = false(height(stage_table), 1);
for r = 1:height(stage_table)
    tc = stage_table.trial_continuous(r);
    if tc >= 1 && tc <= total_trials
        stage_table.fb_shown_correct(r) = fb_shown_correct_vec(tc);
        stage_table.false_fb(r) = false_fb_vec(tc);
    end
end
end

function stage_table = ensure_kh_feedback_columns(stage_table)
if ~ismember('false_fb', stage_table.Properties.VariableNames)
    stage_table.false_fb = false(height(stage_table), 1);
end
if ~ismember('fb_shown_correct', stage_table.Properties.VariableNames)
    stage_table.fb_shown_correct = nan(height(stage_table), 1);
end
end

function subj_features = build_subject_feature_spine(cfg, subj, participant, beh, stage_table, trial2epoch, total_trials)

% Start from stage_table when possible
if isempty(stage_table)
    error('stage_table is empty.');
end

subj_features = stage_table;

% If there is a cohort behav_table row subset, use it only for extra columns
if ~isempty(cfg.behav_table) && ismember('subjID', cfg.behav_table.Properties.VariableNames)
    rows = string(cfg.behav_table.subjID) == string(subj);
    if any(rows)
        beh_rows = cfg.behav_table(rows, :);
        n = min(height(subj_features), height(beh_rows));
    else
        beh_rows = table();
        n = height(subj_features);
    end
else
    beh_rows = table();
    n = height(subj_features);
end

subj_features = subj_features(1:n,:);
trial2epoch   = trial2epoch(1:n);

% Overwrite canonical columns from stage_table only
core_vars = {'stage','block_type','correct','false_fb', ...
             'fb_shown_correct','confidence','block_number','trial_in_block'};

for vi = 1:numel(core_vars)
    vn = core_vars{vi};
    if ismember(vn, stage_table.Properties.VariableNames)
        subj_features.(vn) = stage_table.(vn)(1:n);
    end
end

% Add any extra behav_table columns that are missing from stage_table
if ~isempty(beh_rows)
    extra_vars = setdiff(beh_rows.Properties.VariableNames, subj_features.Properties.VariableNames);
    for vi = 1:numel(extra_vars)
        vn = extra_vars{vi};
        subj_features.(vn) = beh_rows.(vn)(1:n);
    end
end

% Set epoch explicitly
subj_features.epoch = trial2epoch(1:n);

end

function T = initialise_feature_columns(T)
n = height(T);
num_cols = {'baseline_rms','N2_amp','N2_lat','N2_norm', ...
    'FCz_neg_peak_amp','FCz_neg_peak_norm','FCz_neg_peak_lat', ...
    'FCz_neg_peak_amp_robust','FCz_neg_peak_lat_robust','FCz_neg_peak_prom', ...
    'FCzCz_mean_250_350','FCzCz_mean_250_350_norm', ...
    'FRN_amp','FRN_norm','FRN_mean_amp','FRN_mean_norm', ...
    'RewP_amp','RewP_norm','RewP_mean_amp','RewP_mean_norm', ...
    'P300_amp','P300_norm','P300_peak_lat','Theta_amp', ...
    'PLV_fp','PLV_fs','PLV_fp_pairwise','PLV_fs_pairwise'};
for i = 1:numel(num_cols)
    if ~ismember(num_cols{i}, T.Properties.VariableNames)
        T.(num_cols{i}) = nan(n,1);
    end
end
if ~ismember('FCz_neg_peak_valid', T.Properties.VariableNames), T.FCz_neg_peak_valid = false(n,1); end
cell_cols = {'FCzCz_waveform','FCzCz_signal','P300_waveform','P300_signal','Theta_waveform','Theta_signal'};
for i = 1:numel(cell_cols)
    if ~ismember(cell_cols{i}, T.Properties.VariableNames)
        T.(cell_cols{i}) = repmat({[]}, n, 1);
    end
end
end

function T = assign_sliding_window_plv(T, trial2epoch, EEGp_phase, acc_idx, par_idx, som_idx, plv_mask, plv_bl, halfwin, min_trials)
if isempty(EEGp_phase) || isempty(acc_idx), return; end
row_idx = (1:height(T))';
valid_ep = ~isnan(trial2epoch(1:height(T))) & trial2epoch(1:height(T)) >= 1 & trial2epoch(1:height(T)) <= EEGp_phase.trials;
unique_buckets = unique([string(T.stage(valid_ep)), string(T.block_type(valid_ep)), string(T.correct(valid_ep))], 'rows');
for ub = 1:size(unique_buckets,1)
    bucket_mask = valid_ep & string(T.stage)==unique_buckets(ub,1) & string(T.block_type)==unique_buckets(ub,2) & string(T.correct)==unique_buckets(ub,3);
    bucket_rows = row_idx(bucket_mask);
    if numel(bucket_rows) < min_trials, continue; end
    if ismember('trial_continuous', T.Properties.VariableNames)
        [~, ord] = sort(T.trial_continuous(bucket_rows));
        bucket_rows = bucket_rows(ord);
    end
    for bi = 1:numel(bucket_rows)
        wrows = bucket_rows(max(1,bi-halfwin):min(numel(bucket_rows),bi+halfwin));
        if numel(wrows) < min_trials, continue; end
        eps = clean_epochs_checked(trial2epoch(wrows), EEGp_phase);
        if numel(eps) < min_trials, continue; end
        cr = bucket_rows(bi);
        if ~isempty(par_idx)
            [plv_ts, ~] = compute_cross_trial_plv(EEGp_phase, acc_idx, par_idx, eps);
            plv_ts = plv_ts - mean(plv_ts(plv_bl), 'omitnan');
            T.PLV_fp(cr) = mean(plv_ts(plv_mask), 'omitnan');
            T.PLV_fp_pairwise(cr) = T.PLV_fp(cr);
        end
        if ~isempty(som_idx)
            [plv_ts, ~] = compute_cross_trial_plv(EEGp_phase, acc_idx, som_idx, eps);
            plv_ts = plv_ts - mean(plv_ts(plv_bl), 'omitnan');
            T.PLV_fs(cr) = mean(plv_ts(plv_mask), 'omitnan');
            T.PLV_fs_pairwise(cr) = T.PLV_fs(cr);
        end
    end
end
end

function grand = accumulate_grand_averages(grand, cfg, participant, stage_table, EEGp, EEGp_theta, EEGp_phase, fcz_idx, par_idx, acc_idx, som_idx, bl_mask, plv_mask, plv_bl, MIN_TRIALS_PLV, stage_names, BTYPE_LABELS)
has_P = any(stage_table.block_type == 'P');
for s_i = 1:4
    for bt_i = 1:2
        bt = BTYPE_LABELS{bt_i};
        if strcmp(bt,'P') && ~has_P, continue; end
        base_mask = stage_table.block_type==bt & stage_table.stage==stage_names{s_i};
        specs = {'correct', base_mask & stage_table.correct==1 & ~stage_table.false_fb, fcz_idx, 'FCz';
                 'incorrect', base_mask & stage_table.correct==0 & ~stage_table.false_fb, fcz_idx, 'FCz';
                 'correct', base_mask & stage_table.correct==1 & ~stage_table.false_fb, par_idx, 'Par';
                 'incorrect', base_mask & stage_table.correct==0 & ~stage_table.false_fb, par_idx, 'Par'};
        for ci = 1:size(specs,1)
            cname = specs{ci,1}; cmask = specs{ci,2}; ch_idx = specs{ci,3}; gfield = specs{ci,4};
            if isempty(ch_idx), continue; end
            dat = extract_epochs(EEGp, ch_idx, stage_table.epoch(cmask), bl_mask);
            if ~isempty(dat)
                grand.(gfield).(stage_names{s_i}).(bt).(cname).data(end+1,:) = mean(dat,1,'omitnan');
                grand.(gfield).(stage_names{s_i}).(bt).(cname).subj(end+1,1) = participant;
            end
        end
        if strcmp(bt,'P') && ~isempty(fcz_idx)
            ff_specs = {'false_cor', base_mask & stage_table.false_fb & stage_table.fb_shown_correct==1;
                        'false_inc', base_mask & stage_table.false_fb & stage_table.fb_shown_correct==0};
            for fi = 1:2
                dat = extract_epochs(EEGp, fcz_idx, stage_table.epoch(ff_specs{fi,2}), bl_mask);
                if ~isempty(dat)
                    grand.FCz.(stage_names{s_i}).(bt).(ff_specs{fi,1}).data(end+1,:) = mean(dat,1,'omitnan');
                    grand.FCz.(stage_names{s_i}).(bt).(ff_specs{fi,1}).subj(end+1,1) = participant;
                end
            end
        end
        if ~isempty(EEGp_theta) && ~isempty(acc_idx)
            for oc = {'correct','incorrect'}
                oc_val = strcmp(oc{1}, 'correct');
                eps_th = clean_epochs_checked(stage_table.epoch(base_mask & stage_table.correct==oc_val & ~stage_table.false_fb), EEGp_theta);
                th_mat = compute_theta_envelope(EEGp_theta, acc_idx, eps_th, bl_mask);
                if ~isempty(th_mat)
                    grand.Theta.(stage_names{s_i}).(bt).(oc{1}).data(end+1,:) = mean(th_mat,1,'omitnan');
                    grand.Theta.(stage_names{s_i}).(bt).(oc{1}).subj(end+1,1) = participant;
                end
            end
        end
        if ~isempty(EEGp_phase) && ~isempty(acc_idx)
            eps_plv = clean_epochs_checked(stage_table.epoch(base_mask & ~stage_table.false_fb), EEGp_phase);
            if numel(eps_plv) >= MIN_TRIALS_PLV
                if ~isempty(par_idx)
                    [plv_ts, ~] = compute_cross_trial_plv(EEGp_phase, acc_idx, par_idx, eps_plv);
                    plv_ts = plv_ts - mean(plv_ts(plv_bl), 'omitnan');
                    grand.PLV_fp.(stage_names{s_i}).(bt).data(end+1,:) = plv_ts;
                    grand.PLV_fp.(stage_names{s_i}).(bt).subj(end+1,1) = participant;
                end
                if ~isempty(som_idx)
                    [plv_ts, ~] = compute_cross_trial_plv(EEGp_phase, acc_idx, som_idx, eps_plv);
                    plv_ts = plv_ts - mean(plv_ts(plv_bl), 'omitnan');
                    grand.PLV_fs.(stage_names{s_i}).(bt).data(end+1,:) = plv_ts;
                    grand.PLV_fs.(stage_names{s_i}).(bt).subj(end+1,1) = participant;
                end
            end
        end
    end
end
end

function T = finalise_group_table_types_and_zscore(T)
% Convert key vars.
T.subj_id = categorical(string(T.subj_id));
T.cohort = categorical(string(T.cohort));
if ismember('stage', T.Properties.VariableNames), T.stage = categorical(string(T.stage), {'LN','LE','RN','RE'}, 'Ordinal', false); end
if ismember('block_type', T.Properties.VariableNames), T.block_type = categorical(string(T.block_type), {'D','P'}); end
if ismember('subj', T.Properties.VariableNames), T.subj = categorical(T.subj); end

features_to_zscore = {'FCz_neg_peak_amp','FCz_neg_peak_amp_robust','FCzCz_mean_250_350', ...
    'FRN_amp','FRN_mean_amp','RewP_amp','RewP_mean_amp', ...
    'P300_amp','Theta_amp','PLV_fp','PLV_fs','PLV_fp_pairwise','PLV_fs_pairwise'};
for f = 1:numel(features_to_zscore)
    fn = features_to_zscore{f};
    if ~ismember(fn, T.Properties.VariableNames), continue; end
    fn_z = [fn '_z'];
    T.(fn_z) = nan(height(T),1);
    subjs = categories(T.subj_id);
    for si = 1:numel(subjs)
        mask = T.subj_id == subjs{si};
        vals = T.(fn)(mask);
        mn = mean(vals, 'omitnan');
        sd = std(vals, 'omitnan');
        if sd > 0
            T.(fn_z)(mask) = (vals - mn) ./ sd;
        end
    end
end
end

function [combined, debug] = combine_handoff_tables_v9(tables)
% Union columns, adding missing columns as NaN/missing to each cohort.
all_vars = {};
for i = 1:numel(tables)
    all_vars = union(all_vars, tables{i}.Properties.VariableNames, 'stable');
end
for i = 1:numel(tables)
    T = tables{i};
    missing = setdiff(all_vars, T.Properties.VariableNames, 'stable');
    for m = 1:numel(missing)
        T.(missing{m}) = make_missing_column(height(T));
    end
    tables{i} = T(:, all_vars);
end
combined = vertcat(tables{:});
combined.subj_id = categorical(string(combined.subj_id));
debug = table(string(all_vars(:)), 'VariableNames', {'variable'});
for i = 1:numel(tables)
    cohort_name = char(string(tables{i}.cohort(1)));
    present = ismember(all_vars, tables{i}.Properties.VariableNames)';
    debug.([cohort_name '_present']) = present;
end
end

function col = make_missing_column(n)
col = nan(n,1);
end

function [valid_pk, pk_amp, pk_lat, pk_prom] = robust_negative_peak(sig, times, win_mask, bline_rms, smooth_n, min_prom_rms, edge_margin_ms, win)
valid_pk = false; pk_amp = NaN; pk_lat = NaN; pk_prom = NaN;
win_vals = sig(win_mask); win_t = times(win_mask);
if isempty(win_vals) || all(isnan(win_vals)) || isnan(bline_rms) || bline_rms <= 0, return; end
smooth_vals = movmean(win_vals, smooth_n, 'omitnan');
try
    [pks, locs, ~, proms] = findpeaks(-smooth_vals, win_t, 'MinPeakProminence', min_prom_rms*bline_rms, 'MinPeakDistance', 40);
catch
    [pks, locs, ~, proms] = findpeaks(-smooth_vals, 'MinPeakProminence', min_prom_rms*bline_rms, 'MinPeakDistance', round(40/mean(diff(win_t))));
    locs = win_t(locs);
end
if isempty(pks), return; end
[~, best] = max(proms);
pk_amp = -pks(best); pk_lat = locs(best); pk_prom = proms(best);
valid_pk = pk_lat > win(1)+edge_margin_ms && pk_lat < win(2)-edge_margin_ms;
if ~valid_pk, pk_amp = NaN; pk_lat = NaN; pk_prom = NaN; end
end

function data_bl = extract_epochs(EEGp, ch_idx, epoch_col, bl_mask)
ep = clean_epochs_checked(epoch_col, EEGp);
if isempty(ep), data_bl = zeros(0, size(EEGp.data,2)); return; end
if isscalar(ch_idx)
    raw = squeeze(double(EEGp.data(ch_idx,:,ep)))';
else
    raw = squeeze(mean(double(EEGp.data(ch_idx,:,ep)), 1, 'omitnan'))';
end
if isvector(raw) && numel(ep)==1, raw = raw(:)'; end
data_bl = raw - mean(raw(:,bl_mask),2,'omitnan');
end

function theta_mat = compute_theta_envelope(EEGp_theta, ch_idx, ep, bl_mask)
if isempty(ep), theta_mat = zeros(0,size(EEGp_theta.data,2)); return; end
if isscalar(ch_idx)
    raw = squeeze(double(EEGp_theta.data(ch_idx,:,ep)))';
else
    raw = squeeze(mean(double(EEGp_theta.data(ch_idx,:,ep)),1,'omitnan'))';
end
if isvector(raw) && numel(ep)==1, raw = raw(:)'; end
theta_mat = nan(size(raw));
for k = 1:size(raw,1)
    env = abs(hilbert(double(raw(k,:))));
    theta_mat(k,:) = env - mean(env(bl_mask),'omitnan');
end
end

function [plv_ts, n_used] = compute_cross_trial_plv(EEGp_phase, ref_idx, tgt_idx, ep)
n_t = size(EEGp_phase.data,2); n_ep = numel(ep);
if n_ep == 0, plv_ts = zeros(1,n_t); n_used = 0; return; end
if isscalar(ref_idx)
    phi_ref = squeeze(double(EEGp_phase.data(ref_idx,:,ep)))';
else
    phi_ref = squeeze(angle(mean(exp(1i*double(EEGp_phase.data(ref_idx,:,ep))),1,'omitnan')))';
end
if isscalar(tgt_idx)
    phi_tgt = squeeze(double(EEGp_phase.data(tgt_idx,:,ep)))';
else
    phi_tgt = squeeze(angle(mean(exp(1i*double(EEGp_phase.data(tgt_idx,:,ep))),1,'omitnan')))';
end
if isvector(phi_ref)&&n_ep==1, phi_ref=phi_ref(:)'; end
if isvector(phi_tgt)&&n_ep==1, phi_tgt=phi_tgt(:)'; end
plv_ts = abs(mean(exp(1i*(phi_ref-phi_tgt)),1,'omitnan'));
n_used = n_ep;
end

function ep = clean_epochs_checked(epoch_col, EEGp)
ep = epoch_col(~isnan(epoch_col) & epoch_col>=1 & epoch_col<=EEGp.trials);
ep = round(ep(:));
end

function idx = safe_chan(EEGp, labels)
idx = find(ismember(lower({EEGp.chanlocs.labels}), lower(labels)));
end

function plot_subject_erps_unified(cfg, subj, participant, has_P, EEGp, EEGp_theta, stage_table, fcz_idx, par_idx, acc_idx, bl_mask, stage_names, STAGE_COLORS, LINE_STYLES, BTYPE_LABELS, ERP_plot_window)
fig = figure('Position',[50 50 1400 560],'Visible','off');
sgtitle(sprintf('%s — FCz ERP by stage', subj), 'Interpreter','none');
for s_i = 1:4
    for oc_i = 1:2
        oc = {'correct','incorrect'};
        ax = subplot(2,4,(oc_i-1)*4+s_i); hold(ax,'on');
        title(ax, sprintf('%s | %s', stage_names{s_i}, oc{oc_i}));
        xline(ax,0,'k:','HandleVisibility','off'); yline(ax,0,'k:','HandleVisibility','off');
        for bt_i = 1:2
            bt = BTYPE_LABELS{bt_i};
            if strcmp(bt,'P') && ~has_P, continue; end
            oc_val = strcmp(oc{oc_i}, 'correct');
            m = stage_table.block_type==bt & stage_table.stage==stage_names{s_i} & stage_table.correct==oc_val & ~stage_table.false_fb;
            dat = extract_epochs(EEGp, fcz_idx, stage_table.epoch(m), bl_mask);
            if isempty(dat), continue; end
            plot_erp_trace(ax, EEGp.times, dat, STAGE_COLORS(s_i,:), sprintf('%s-%s', oc{oc_i}, bt), ERP_plot_window, LINE_STYLES{bt_i});
        end
        set(ax,'YDir','reverse'); xlabel(ax,'Time (ms)'); ylabel(ax,'\muV'); xlim(ax,ERP_plot_window); legend(ax,'Box','off','FontSize',7);
    end
end
saveas(fig, fullfile(cfg.figure_output_folder, sprintf('%s_FCz_stage_v9.pdf', subj)));
close(fig);
end

function h = plot_erp_trace(ax, times, data_mat, clr, lbl, xlims, ls)
if nargin < 7, ls = '-'; end
if isempty(data_mat), h = []; return; end
in_win = times >= xlims(1) & times <= xlims(2);
mn = mean(data_mat(:,in_win),1,'omitnan');
se = std(data_mat(:,in_win),0,1,'omitnan') ./ sqrt(size(data_mat,1));
t = times(in_win);
fill(ax,[t fliplr(t)],[mn+se fliplr(mn-se)], clr, 'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
h = plot(ax,t,mn,'Color',clr,'LineWidth',1.5,'LineStyle',ls,'DisplayName',sprintf('%s (n=%d)',lbl,size(data_mat,1)));
end
