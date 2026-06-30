% =============================================================================
% OUTCOME ERP ANALYSIS v8
%
% INTEGRATION of v7 (grand average / PLV / per-subject ERPs) with
% feature extraction script (single-trial scalar features + waveform storage).
%
% What was redundant in the feature extraction script:
%   - EEG loading:             IDENTICAL to v7 — merged
%   - beh struct handling:     IDENTICAL to v7 — merged
%   - define_trial_stages_v3:  IDENTICAL to v7 — merged (one call)
%   - extract_epochs():        IDENTICAL logic — one shared function
%   - compute_theta_envelope:  IDENTICAL — one shared function
%   - compute_cross_trial_plv: IDENTICAL — one shared function
%   - Channel index lookups:   IDENTICAL — computed once
%   - bl_mask:                 IDENTICAL — computed once
%
% What the feature extraction script added that v7 lacked:
%   [A] Single-trial scalar features: N2_amp, N2_lat, N2_norm,
%       P300_amp, P300_peak_lat, P300_norm, Theta_amp — now added
%   [B] Per-trial waveform storage in group_table:
%       FCzCz_waveform, P300_waveform, Theta_waveform — now added
%   [C] trial2epoch .mat caching (load if exists, create+save if not)
%   [D] Part 2: flexible ERP aggregation by grouping schemes — kept as-is
%   [E] Behav table (behav_table / group_T) as the metadata spine —
%       v7 used stage_table only; v8 uses behav_table and appends EEG cols
%
% What was in the feature extraction script that is NOT merged:
%   - Part 2 aggregation loop: kept at the bottom as its own section.
%     It needs group_table_all_trials + time_vector. Both are now saved
%     during the participant loop.
%
% PLV baseline: changed to [-500 0] ms as requested.
% NOTE on PLV baseline: -500 to 0 is the full pre-stimulus window and is
% the correct choice IF your epochs start at -500 ms. If there is
% substantial edge artefact in the first ~100 ms of the epoch (common with
% Hilbert on short epochs), you may see a non-flat baseline. Watch for this
% in the PLV plots. The previous [-200 -20] choice was conservative for
% exactly that reason. If PLV looks non-flat before t=0, revert to [-200 0].
%
% BUG FIXES from v7 bug report applied here:
%   M1: save_tables defined at top (default true)
%   M2: t_ax assigned explicitly after first EEG load
%   M3: epoch count checked after fallback file loads
%   M4: D/P pooling averaging fixed (count conditions actually found)
%   M5: subgroup_plot_grand calls validated against valid_participants
%   v6-M5: false_fb column guard added after stage_table is built
%
% REQUIRES:
%   FieldTrip, EEGLAB, define_trial_stages_v3, validate_stage_table,
%   KH_align_epochs_with_offset, permutation_test_cluster_correction
% =============================================================================

clear; close all;

% -------------------------------------------------------------------------
%% PATHS
% -------------------------------------------------------------------------
fieldtrip_path       = 'C:\Users\khatinova\OneDrive - Nexus365\Pre_2026_Folders\Documents\MATLAB\fieldtrip-20240110';
eeglab_path          = 'C:\Users\khatinova\OneDrive - Nexus365\Pre_2026_Folders\Documents\MATLAB\eeglab2025.1.0';
data_path            = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Data';
epoch_file_folder    = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\EEG analysis\Epoched_data';
figure_output_folder = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\EEG analysis\Figures';

addpath(fieldtrip_path); ft_defaults;
addpath(eeglab_path);    eeglab nogui;
addpath(genpath(data_path));

% -------------------------------------------------------------------------
%% PARAMETERS
% -------------------------------------------------------------------------
valid_participants = [3:8, 10:12, 14:23, 27];
save_tables        = true;          % BUG FIX v7-M1

ERP_plot_window    = [-200 1000];
rm_baseline        = [-200 0];      % ERP baseline window

% Single-trial feature windows
N2_win    = [120 350];
P300_win  = [300 600];
Theta_win = [200 500];
PLV_win   = [200 400];

% Group-level windows (condition-average only — NOT for single-trial)
FRN_win_group  = [120 350];
RewP_win_group = [200 400];

% PLV baseline: full pre-stimulus window.
% NOTE: if PLV plots show non-flat baselines, narrow this to [-200 0].
PLV_baseline   = [-500 0];

MIN_TRIALS_PLV = 5;

% Single-trial scalar PLV (windowed): true cross-trial PLV requires
% multiple trials, so single-trial PLV_fp/PLV_fs/PLV_fp_pairwise/
% PLV_fs_pairwise are defined here as a SLIDING-WINDOW cross-trial PLV:
% for each trial, take the +/- PLV_WINDOW_HALF surrounding trials
% (same subject, same stage x block_type x outcome bucket), compute
% cross-trial PLV across that local window in PLV_win, and assign the
% resulting scalar to the center trial. This means PLV_fp/PLV_fs are
% NOT independent across nearby trials (they share most of their
% window) — treat as a smoothed local estimate, not a true single-trial
% measure. PLV_fp_pairwise/PLV_fs_pairwise use the same method but are
% kept as separate columns in case you later want a genuinely different
% pairwise definition (e.g. specific channel pairs rather than averaged
% channel groups) — currently they duplicate PLV_fp/PLV_fs.
PLV_WINDOW_HALF = 7;   % +/- 7 trials = up to 15-trial window
MIN_TRIALS_PLV_WINDOW = 5;  % minimum trials in window to compute a value

% Channel groups
fcz_label              = 'FCz';
frontocentral_channels = {'FCz','Cz','F1','F2'};
acc_channels           = {'FCz','Fz','AFz','F1','F2'};
par_channels           = {'Pz','P1','P2'};
som_channels           = {'C3','C4','CP3','CP1','C5','CP5'};

stage_names  = {'LN','LE','RN','RE'};
STAGE_COLORS = [0.12 0.62 0.47; 0.85 0.65 0.00; 0.80 0.27 0.13; 0.40 0.25 0.65];
LINE_STYLES  = {'-','--'};
BTYPE_LABELS = {'D','P'};

% Sub-groups for re-plotting (must be subsets of valid_participants)
SUBGROUPS = struct( ...
    'cohort1_audio',  {{3,4,5,6,7,8}}, ...
    'cohort2_pilot',  {{10}}, ...            % 9 removed — not in valid_participants
    'cohort2_main',   {{11,12,14:20}});      % BUG FIX v7-M5: removed 21-23,27

% -------------------------------------------------------------------------
%% LOAD DATA
% -------------------------------------------------------------------------
load(fullfile(data_path, 'all_trial_data.mat'));
load(fullfile(data_path, 'behav_table.mat'), 'group_T');
behav_table = group_T; clear group_T;

for i = valid_participants
    ps = sprintf('Ox%02d', i);
    data_by_participant{i} = all_trial_data.(ps).trial_data;
end

% -------------------------------------------------------------------------
%% INITIALISE GRAND AVERAGE CONTAINERS
% Each leaf: struct with .data (n_subj × n_time) and .subj (n_subj × 1)
% -------------------------------------------------------------------------
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
        grand.Theta.(stage_names{s}).(bt_s).correct  = empty_container;
        grand.Theta.(stage_names{s}).(bt_s).incorrect= empty_container;
        grand.PLV_fp.(stage_names{s}).(bt_s)         = empty_container;
        grand.PLV_fs.(stage_names{s}).(bt_s)         = empty_container;
    end
end

all_stage_tables  = {};
all_trials_table  = table();     % [NEW] single-trial feature table
perm_stats        = struct();
t_ax              = [];          % BUG FIX v7-M2: explicit initialisation

% =========================================================================
%% PARTICIPANT LOOP
% =========================================================================
for participant = valid_participants

    subj = sprintf('Ox%02d', participant);
    fprintf('\n============ %s ============\n', subj);

    is_cohort1 = (participant <= 8);
    beh        = data_by_participant{participant};

    % ------------------------------------------------------------------
    % Behavioural vectors (handle 5 vs 6 block case)
    % ------------------------------------------------------------------
    num_blocks  = height(beh.correct);
    beh_correct = [];
    beh_conf    = [];
    if num_blocks < 6
        for b = 1:num_blocks
            beh_correct = [beh_correct, beh.correct(b,:)];
            beh_conf    = [beh_conf,    beh.confidence(b,:)];
        end
    else
        for b = 2:num_blocks
            beh_correct = [beh_correct, beh.correct(b,:)];
            beh_conf    = [beh_conf,    beh.confidence(b,:)];
        end
        num_blocks     = num_blocks - 1;
        beh.correct    = beh.correct(2:end,:);
        beh.confidence = beh.confidence(2:end,:);
        beh.trueFB     = beh.trueFB(2:end,:);
        beh.revTrial   = beh.revTrial(2:end);
    end

    % ------------------------------------------------------------------
    % Load EEG files
    % ------------------------------------------------------------------
    broadband_fname = sprintf('%s_outcome_trimmed.set', subj);
    if ~exist(fullfile(epoch_file_folder, broadband_fname),'file')
        broadband_fname = sprintf('%s_outcome.set', subj);
    end
    EEGp = pop_loadset(broadband_fname, epoch_file_folder);
    fprintf('  Broadband: %d epochs\n', EEGp.trials);

    % BUG FIX v7-M2: assign t_ax from first successful load
    if isempty(t_ax)
        t_ax = EEGp.times;
        fprintf('  t_ax assigned: %d timepoints [%.0f %.0f] ms\n', ...
            numel(t_ax), t_ax(1), t_ax(end));
    end

    theta_fname = sprintf('%s_outcome_theta_trimmed.set', subj);
    if ~exist(fullfile(epoch_file_folder, theta_fname),'file')
        theta_fname = sprintf('%s_outcome_theta.set', subj);
    end
    EEGp_theta = pop_loadset(theta_fname, epoch_file_folder);
    fprintf('  Theta:     %d epochs\n', EEGp_theta.trials);

    % BUG FIX v7-M3: warn if epoch count mismatches broadband
    if EEGp_theta.trials ~= EEGp.trials
        warning('%s: theta epochs (%d) ≠ broadband (%d) — fallback file loaded', ...
            subj, EEGp_theta.trials, EEGp.trials);
    end

    phase_fname = sprintf('%s_outcome_phase_trimmed.set', subj);
    if ~exist(fullfile(epoch_file_folder, phase_fname),'file')
        phase_fname = sprintf('%s_outcome_phase.set', subj);
    end
    EEGp_phase = pop_loadset(phase_fname, epoch_file_folder);
    fprintf('  Phase:     %d epochs\n', EEGp_phase.trials);

    if EEGp_phase.trials ~= EEGp.trials
        warning('%s: phase epochs (%d) ≠ broadband (%d) — fallback file loaded', ...
            subj, EEGp_phase.trials, EEGp.trials);
    end

    % ------------------------------------------------------------------
    % [NEW from feature script] Trial2epoch: load cache or create+save
    % ------------------------------------------------------------------
    trial2epoch_file = fullfile(epoch_file_folder, sprintf('%s_trial2epoch.mat', subj));
    if exist(trial2epoch_file, 'file')
        S = load(trial2epoch_file);
        if isfield(S,'trial2epoch_out'), trial2epoch = S.trial2epoch_out;
        elseif isfield(S,'trial2epoch'), trial2epoch = S.trial2epoch;
        else, error('%s: trial2epoch file has unexpected field names.', trial2epoch_file);
        end
        fprintf('  trial2epoch: loaded from cache\n');
    else
        beh_cv = beh.correct'; beh_cv = beh_cv(:);
        beh_cv = beh_cv(~isnan(beh_cv));
        [trial2epoch, diag] = KH_align_epochs_with_offset(EEGp, beh_cv);
        fprintf('  trial2epoch: %d/%d matched (%.1f%%), offset=%d\n', ...
            diag.n_matched, diag.n_trials, diag.match_rate*100, diag.best_offset);
        trial2epoch_out = trial2epoch;
        save(trial2epoch_file, 'trial2epoch_out');
    end
    trial2epoch = trial2epoch(:);

    % ------------------------------------------------------------------
    % Stage table
    % ------------------------------------------------------------------
    stage_table = define_trial_stages_v3(beh, trial2epoch, EEGp, participant, ...
        EEGp_theta, EEGp_phase);
    validate_stage_table(stage_table, EEGp, beh, participant);

    % BUG FIX v6-M5: ensure false_fb column exists
    if ~ismember('false_fb', stage_table.Properties.VariableNames)
        stage_table.false_fb = false(height(stage_table), 1);
        warning('%s: false_fb column missing from stage_table — filled with false', subj);
    end
    if ~ismember('fb_shown_correct', stage_table.Properties.VariableNames)
        stage_table.fb_shown_correct = nan(height(stage_table), 1);
    end

    % Add confidence
    conf_vec = beh_conf(:);
    stage_table.confidence = nan(height(stage_table), 1);
    for r = 1:height(stage_table)
        tc = stage_table.trial_continuous(r);
        if tc >= 1 && tc <= numel(conf_vec)
            stage_table.confidence(r) = conf_vec(tc);
        end
    end
    stage_table.subj       = repmat(participant, height(stage_table), 1);
    stage_table.is_cohort1 = repmat(is_cohort1, height(stage_table), 1);
    all_stage_tables{end+1} = stage_table;

    % ------------------------------------------------------------------
    % Channel indices and time masks (computed once, shared everywhere)
    % ------------------------------------------------------------------
    fcz_idx = safe_chan(EEGp, {fcz_label});
    cz_idx  = safe_chan(EEGp, {'Cz'});
    fc_idx  = safe_chan(EEGp, frontocentral_channels);
    par_idx = safe_chan(EEGp, par_channels);
    acc_idx = safe_chan(EEGp, acc_channels);
    som_idx = safe_chan(EEGp, som_channels);

    bl_mask   = EEGp.times >= rm_baseline(1) & EEGp.times <= rm_baseline(2);
    n2_mask   = EEGp.times >= N2_win(1)      & EEGp.times <= N2_win(2);
    p300_mask = EEGp.times >= P300_win(1)    & EEGp.times <= P300_win(2);
    th_mask   = EEGp.times >= Theta_win(1)   & EEGp.times <= Theta_win(2);
    plv_bl    = EEGp.times >= PLV_baseline(1) & EEGp.times <= PLV_baseline(2);
    plv_mask  = EEGp.times >= PLV_win(1)      & EEGp.times <= PLV_win(2);
    % [NEW] single-trial FRN/RewP windows (same windows as group-level,
    % applied per trial — see header note on what "single-trial FRN" means)
    frn_mask_st  = EEGp.times >= FRN_win_group(1)  & EEGp.times <= FRN_win_group(2);
    rewp_mask_st = EEGp.times >= RewP_win_group(1) & EEGp.times <= RewP_win_group(2);

    has_P = any(stage_table.block_type == 'P');

    % Baseline RMS (single scalar per subject, used for normalisation)
    if ~isempty(fcz_idx)
        bl_data  = squeeze(EEGp.data(fcz_idx, bl_mask, :));
        bline_rms = rms(bl_data(:), 'omitnan');
    else
        bline_rms = NaN;
    end

    % ==================================================================
    %% [NEW] BUILD SINGLE-TRIAL FEATURE TABLE for this subject
    % Spine = rows from behav_table; EEG columns appended trial-by-trial.
    % ==================================================================
    subj_beh_rows = string(behav_table.subjID) == string(subj);
    subj_beh      = behav_table(subj_beh_rows, :);
    n_beh_trials  = height(subj_beh);

    if sum(subj_beh_rows) == 0
        warning('No behavioral data found for %s, skipping single-trial table build', subj);
        continue;
    end

    % Ensure subj_num exists (numeric subject ID) — needed by merge script
    if ~ismember('subj_num', subj_beh.Properties.VariableNames)
        subj_beh.subj_num = repmat(participant, n_beh_trials, 1);
    end
    % Ensure 'subj' (numeric) exists — merge script reads double(group_table_KH.subj)
    if ~ismember('subj', subj_beh.Properties.VariableNames)
        subj_beh.subj = repmat(participant, n_beh_trials, 1);
    end

    % Ensure trial2epoch matches behav_table height
    if numel(trial2epoch) < n_beh_trials
        trial2epoch(end+1:n_beh_trials) = NaN;
    elseif numel(trial2epoch) > n_beh_trials
        trial2epoch = trial2epoch(1:n_beh_trials);
    end

    subj_features = subj_beh;
    subj_features.epoch         = trial2epoch;
    subj_features.has_eeg_epoch = ~isnan(trial2epoch);

    % Stage labels from stage_table
    subj_features.stage           = strings(n_beh_trials, 1); subj_features.stage(:) = missing;
    subj_features.trigger_source  = strings(n_beh_trials, 1); subj_features.trigger_source(:) = missing;
    subj_features.had_overlap     = nan(n_beh_trials, 1);
    subj_features.false_fb        = nan(n_beh_trials, 1);
    subj_features.fb_shown_correct= nan(n_beh_trials, 1);

    % [NEW] block_type / block_number / trial_in_block — required by the
    % merge script (block_number == 1 mask) and RQ script (block_number
    % for next_correct, trial_in_block for the block-transition analysis).
    % Mirrors the fallback logic from the feature-extraction script.
    if ~ismember('block_type', subj_features.Properties.VariableNames)
        subj_features.block_type = strings(n_beh_trials, 1);
        subj_features.block_type(:) = missing;
    else
        subj_features.block_type = string(subj_features.block_type);
    end

    if ~ismember('block_number', subj_features.Properties.VariableNames)
        if ismember('block', subj_features.Properties.VariableNames)
            subj_features.block_number = subj_features.block;
        else
            subj_features.block_number = nan(n_beh_trials, 1);
        end
    end

    if ~ismember('trial_in_block', subj_features.Properties.VariableNames)
        subj_features.trial_in_block = nan(n_beh_trials, 1);
        if ismember('block', subj_features.Properties.VariableNames)
            blocks_here = unique(subj_features.block);
            blocks_here = blocks_here(~isnan(blocks_here));
            for bb = 1:numel(blocks_here)
                idx = find(subj_features.block == blocks_here(bb));
                subj_features.trial_in_block(idx) = (1:numel(idx))';
            end
        end
    end

    % false_fb fallback from trueFB if behav_table carries it
    if all(isnan(subj_features.false_fb)) && ismember('trueFB', subj_features.Properties.VariableNames)
        subj_features.false_fb = double(subj_features.trueFB == 0);
    end

    % EEG scalar features
    subj_features.baseline_rms    = repmat(bline_rms, n_beh_trials, 1);
    subj_features.N2_amp          = nan(n_beh_trials, 1);
    subj_features.N2_lat          = nan(n_beh_trials, 1);
    subj_features.N2_norm         = nan(n_beh_trials, 1);
    subj_features.FCz_neg_peak_amp = nan(n_beh_trials, 1);  % backward-compat alias
    subj_features.FCz_neg_peak_lat = nan(n_beh_trials, 1);
    subj_features.FCz_neg_peak_norm= nan(n_beh_trials, 1);
    subj_features.P300_amp        = nan(n_beh_trials, 1);
    subj_features.P300_peak_lat   = nan(n_beh_trials, 1);
    subj_features.P300_norm       = nan(n_beh_trials, 1);
    subj_features.Theta_amp       = nan(n_beh_trials, 1);

    % [NEW] Single-trial windowed FRN/RewP — measured on the same
    % FCzCz_waveform as N2, just a different time window and direction.
    % These exist so group_table downstream has FRN_amp/RewP_amp columns
    % matching what the merging + RQ scripts expect. They are single-trial
    % analogs of the group-level FRN/RewP windows, NOT the same as the
    % group-level condition-averaged FRN_amp computed in aggregate_erp_data.
    subj_features.FRN_amp         = nan(n_beh_trials, 1);
    subj_features.FRN_norm        = nan(n_beh_trials, 1);
    subj_features.RewP_amp        = nan(n_beh_trials, 1);
    subj_features.RewP_norm       = nan(n_beh_trials, 1);

    % [NEW] Single-trial PLV (sliding-window cross-trial estimate, see
    % PLV_WINDOW_HALF note above). Filled in a second pass after the
    % main trial loop below, once all FCzCz/phase data is available.
    subj_features.PLV_fp           = nan(n_beh_trials, 1);
    subj_features.PLV_fs           = nan(n_beh_trials, 1);
    subj_features.PLV_fp_pairwise  = nan(n_beh_trials, 1);
    subj_features.PLV_fs_pairwise  = nan(n_beh_trials, 1);

    % Waveform storage (cell arrays)
    subj_features.FCzCz_waveform  = repmat({[]}, n_beh_trials, 1);
    subj_features.P300_waveform   = repmat({[]}, n_beh_trials, 1);
    subj_features.Theta_waveform  = repmat({[]}, n_beh_trials, 1);

    % Fill stage labels from stage_table
    for r = 1:height(stage_table)
        tc = stage_table.trial_continuous(r);
        if tc < 1 || tc > n_beh_trials, continue; end
        if ismember('stage',           stage_table.Properties.VariableNames), subj_features.stage(tc)           = string(stage_table.stage(r)); end
        if ismember('trigger_source',  stage_table.Properties.VariableNames), subj_features.trigger_source(tc)  = string(stage_table.trigger_source(r)); end
        if ismember('had_overlap',     stage_table.Properties.VariableNames), subj_features.had_overlap(tc)     = stage_table.had_overlap(r); end
        if ismember('false_fb',        stage_table.Properties.VariableNames), subj_features.false_fb(tc)        = double(stage_table.false_fb(r)); end
        if ismember('fb_shown_correct',stage_table.Properties.VariableNames), subj_features.fb_shown_correct(tc)= stage_table.fb_shown_correct(r); end
        % [NEW] block_type/block_number from stage_table take priority
        % over behav_table fallbacks, since stage_table is built from the
        % same alignment pipeline used for ERP figures.
        if ismember('block_type',   stage_table.Properties.VariableNames), subj_features.block_type(tc)   = string(stage_table.block_type(r)); end
        if ismember('block_number', stage_table.Properties.VariableNames), subj_features.block_number(tc) = stage_table.block_number(r); end
    end

    % Trial-by-trial EEG feature extraction
    for ti = 1:n_beh_trials
        ep = trial2epoch(ti);
        if isnan(ep) || ep < 1 || ep > EEGp.trials, continue; end
        ep = round(ep);

        % FCz+Cz waveform and N2 peak
        if ~isempty(fcz_idx) && ~isempty(cz_idx)
            sig = mean(double(EEGp.data([fcz_idx cz_idx], :, ep)), 1);
            sig = sig - mean(sig(bl_mask), 'omitnan');
            subj_features.FCzCz_waveform{ti} = sig;
            win_vals = sig(n2_mask); win_t = EEGp.times(n2_mask);
            if any(~isnan(win_vals))
                [pk, ix] = min(win_vals, [], 'omitnan');
                subj_features.N2_amp(ti) = pk;
                subj_features.N2_lat(ti) = win_t(ix);
                subj_features.FCz_neg_peak_amp(ti) = pk;
                subj_features.FCz_neg_peak_lat(ti) = win_t(ix);
                if bline_rms > 0
                    subj_features.N2_norm(ti) = pk / bline_rms;
                    subj_features.FCz_neg_peak_norm(ti) = pk / bline_rms;
                end
            end

            % [NEW] Single-trial FRN (min in FRN_win_group) and RewP
            % (max in RewP_win_group), measured on the same waveform.
            % These are genuinely single-trial (no window-sharing issue
            % like PLV below), just different windows from N2_win.
            frn_vals = sig(frn_mask_st); frn_t = EEGp.times(frn_mask_st);
            if any(~isnan(frn_vals))
                [pk_frn, ix_frn] = min(frn_vals, [], 'omitnan');
                subj_features.FRN_amp(ti) = pk_frn;
                if bline_rms > 0, subj_features.FRN_norm(ti) = pk_frn / bline_rms; end
            end
            rewp_vals = sig(rewp_mask_st); rewp_t = EEGp.times(rewp_mask_st);
            if any(~isnan(rewp_vals))
                [pk_rewp, ix_rewp] = max(rewp_vals, [], 'omitnan');
                subj_features.RewP_amp(ti) = pk_rewp;
                if bline_rms > 0, subj_features.RewP_norm(ti) = pk_rewp / bline_rms; end
            end
        end

        % P300
        if ~isempty(par_idx)
            sig = mean(double(EEGp.data(par_idx, :, ep)), 1);
            sig = sig - mean(sig(bl_mask), 'omitnan');
            subj_features.P300_waveform{ti} = sig;
            win_vals = sig(p300_mask); win_t = EEGp.times(p300_mask);
            if any(~isnan(win_vals))
                [pk, ix] = max(win_vals, [], 'omitnan');
                subj_features.P300_amp(ti) = pk;
                subj_features.P300_peak_lat(ti) = win_t(ix);
                if bline_rms > 0, subj_features.P300_norm(ti) = pk / bline_rms; end
            end
        end

        % Theta amplitude (from theta-filtered file)
        if ~isempty(acc_idx) && ep <= EEGp_theta.trials
            sig_th = mean(double(EEGp_theta.data(acc_idx, :, ep)), 1);
            env    = abs(hilbert(sig_th));
            env    = env - mean(env(bl_mask), 'omitnan');
            subj_features.Theta_amp(ti) = mean(env(th_mask), 'omitnan');
            subj_features.Theta_waveform{ti} = env;
        end
    end

    all_trials_table = [all_trials_table; subj_features]; %#ok<AGROW>

    % ==================================================================
    %% [NEW] SLIDING-WINDOW SINGLE-TRIAL PLV
    %
    % True PLV needs multiple trials (it's cross-trial phase consistency),
    % so a single trial cannot have its own independent PLV value. Here
    % we compute a LOCAL windowed PLV: for each trial, take a window of
    % +/- PLV_WINDOW_HALF neighboring trials within the SAME stage x
    % block_type x outcome bucket (so the window doesn't mix conditions),
    % compute cross-trial PLV across that window in PLV_win, and assign
    % the resulting scalar (mean PLV across PLV_win, baseline-corrected)
    % to the center trial. Neighboring trials share most of their window,
    % so these values are smoothed/non-independent — treat PLV_fp/PLV_fs
    % as a local estimate, not a true single-trial measure.
    %
    % PLV_fp_pairwise / PLV_fs_pairwise currently duplicate PLV_fp/PLV_fs
    % (same method). Kept as separate columns per the merge script's
    % expected schema in case a genuinely different pairwise definition
    % is added later.
    % ==================================================================
    if ~isempty(acc_idx) && EEGp_phase.trials > 0
        % Work in the same row-index space as subj_features
        row_idx = (1:n_beh_trials)';
        valid_ep = ~isnan(trial2epoch) & trial2epoch >= 1 & trial2epoch <= EEGp_phase.trials;

        bucket_stage = subj_features.stage;
        bucket_bt    = subj_features.block_type;
        bucket_corr  = subj_features.correct;

        unique_buckets = unique([string(bucket_stage(valid_ep)), ...
                                  string(bucket_bt(valid_ep)), ...
                                  string(bucket_corr(valid_ep))], 'rows');

        for ub = 1:size(unique_buckets,1)
            s_lbl  = unique_buckets(ub,1);
            bt_lbl = unique_buckets(ub,2);
            c_lbl  = unique_buckets(ub,3);

            bucket_mask = valid_ep & ...
                string(bucket_stage)==s_lbl & string(bucket_bt)==bt_lbl & ...
                string(bucket_corr)==c_lbl;

            bucket_rows = row_idx(bucket_mask);
            if numel(bucket_rows) < MIN_TRIALS_PLV_WINDOW, continue; end

            % Order by trial_continuous so "neighboring trials" means
            % temporally adjacent within this bucket, not table order
            if ismember('trial_continuous', subj_features.Properties.VariableNames)
                [~, ord] = sort(subj_features.trial_continuous(bucket_rows));
                bucket_rows = bucket_rows(ord);
            end

            n_bucket = numel(bucket_rows);
            for bi = 1:n_bucket
                win_lo = max(1, bi - PLV_WINDOW_HALF);
                win_hi = min(n_bucket, bi + PLV_WINDOW_HALF);
                window_rows = bucket_rows(win_lo:win_hi);
                if numel(window_rows) < MIN_TRIALS_PLV_WINDOW, continue; end

                eps_window = clean_epochs_checked(trial2epoch(window_rows), EEGp_phase);
                if numel(eps_window) < MIN_TRIALS_PLV_WINDOW, continue; end

                center_row = bucket_rows(bi);

                if ~isempty(par_idx)
                    [plv_ts, ~] = compute_cross_trial_plv(EEGp_phase, acc_idx, par_idx, eps_window);
                    plv_ts = plv_ts - mean(plv_ts(plv_bl), 'omitnan');
                    plv_scalar = mean(plv_ts(plv_mask), 'omitnan');
                    subj_features.PLV_fp(center_row) = plv_scalar;
                    subj_features.PLV_fp_pairwise(center_row) = plv_scalar;
                end
                if ~isempty(som_idx)
                    [plv_ts, ~] = compute_cross_trial_plv(EEGp_phase, acc_idx, som_idx, eps_window);
                    plv_ts = plv_ts - mean(plv_ts(plv_bl), 'omitnan');
                    plv_scalar = mean(plv_ts(plv_mask), 'omitnan');
                    subj_features.PLV_fs(center_row) = plv_scalar;
                    subj_features.PLV_fs_pairwise(center_row) = plv_scalar;
                end
            end
        end

        % Write back into all_trials_table (subj_features was appended
        % by value above, so re-write the just-appended block directly)
        rows_this_subj = (height(all_trials_table)-n_beh_trials+1):height(all_trials_table);
        all_trials_table.PLV_fp(rows_this_subj)          = subj_features.PLV_fp;
        all_trials_table.PLV_fs(rows_this_subj)          = subj_features.PLV_fs;
        all_trials_table.PLV_fp_pairwise(rows_this_subj) = subj_features.PLV_fp_pairwise;
        all_trials_table.PLV_fs_pairwise(rows_this_subj) = subj_features.PLV_fs_pairwise;

        fprintf('  Single-trial PLV: %d/%d trials assigned (windowed, +/-%d trials)\n', ...
            sum(~isnan(subj_features.PLV_fp)), n_beh_trials, PLV_WINDOW_HALF);
    end

    % ==================================================================
    %% GRAND AVERAGE ACCUMULATION (unchanged from v7, uses stage_table)
    % ==================================================================
    for s_i = 1:4
        for bt_i = 1:2
            bt = BTYPE_LABELS{bt_i};
            if strcmp(bt,'P') && ~has_P, continue; end

            base_mask = stage_table.block_type==bt & stage_table.stage==stage_names{s_i};

            % FCz and Par (both from EEGp)
            cond_specs = { ...
                'correct',   base_mask & stage_table.correct==1 & ~stage_table.false_fb, fcz_idx, 'FCz'; ...
                'incorrect', base_mask & stage_table.correct==0 & ~stage_table.false_fb, fcz_idx, 'FCz'; ...
                'correct',   base_mask & stage_table.correct==1 & ~stage_table.false_fb, par_idx, 'Par'; ...
                'incorrect', base_mask & stage_table.correct==0 & ~stage_table.false_fb, par_idx, 'Par'; ...
            };
            for ci = 1:size(cond_specs,1)
                cname  = cond_specs{ci,1};
                cmask  = cond_specs{ci,2};
                ch_idx = cond_specs{ci,3};
                gfield = cond_specs{ci,4};
                if isempty(ch_idx), continue; end
                dat = extract_epochs(EEGp, ch_idx, stage_table.epoch(cmask), bl_mask);
                if ~isempty(dat)
                    grand.(gfield).(stage_names{s_i}).(bt).(cname).data(end+1,:) = mean(dat,1,'omitnan');
                    grand.(gfield).(stage_names{s_i}).(bt).(cname).subj(end+1,1) = participant;
                end
            end

            % False feedback (cohort 2, P blocks)
            if strcmp(bt,'P') && ~is_cohort1 && ~isempty(fcz_idx)
                for ff_spec = {{'false_cor', base_mask & stage_table.false_fb & stage_table.fb_shown_correct==1}, ...
                               {'false_inc', base_mask & stage_table.false_fb & stage_table.fb_shown_correct==0}}
                    ff_name = ff_spec{1}{1};
                    ff_mask = ff_spec{1}{2};
                    dat = extract_epochs(EEGp, fcz_idx, stage_table.epoch(ff_mask), bl_mask);
                    if ~isempty(dat)
                        grand.FCz.(stage_names{s_i}).(bt).(ff_name).data(end+1,:) = mean(dat,1,'omitnan');
                        grand.FCz.(stage_names{s_i}).(bt).(ff_name).subj(end+1,1) = participant;
                    end
                end
            end

            % Theta
            if ~isempty(acc_idx) && EEGp_theta.trials > 0
                for oc = {'correct','incorrect'}
                    oc_val = strcmp(oc{1},'correct');
                    m_th   = base_mask & stage_table.correct==oc_val & ~stage_table.false_fb;
                    eps_th = clean_epochs_checked(stage_table.epoch(m_th), EEGp_theta);
                    th_mat = compute_theta_envelope(EEGp_theta, acc_idx, eps_th, bl_mask);
                    if ~isempty(th_mat)
                        grand.Theta.(stage_names{s_i}).(bt).(oc{1}).data(end+1,:) = mean(th_mat,1,'omitnan');
                        grand.Theta.(stage_names{s_i}).(bt).(oc{1}).subj(end+1,1) = participant;
                    end
                end
            end

            % PLV
            if ~isempty(acc_idx) && EEGp_phase.trials > 0
                eps_plv = clean_epochs_checked( ...
                    stage_table.epoch(base_mask & ~stage_table.false_fb), EEGp_phase);
                if numel(eps_plv) >= MIN_TRIALS_PLV
                    for pp_spec = {{'PLV_fp', par_idx}, {'PLV_fs', som_idx}}
                        pname   = pp_spec{1}{1};
                        tgt_idx = pp_spec{1}{2};
                        if isempty(tgt_idx), continue; end
                        [plv_ts, ~] = compute_cross_trial_plv(EEGp_phase, acc_idx, tgt_idx, eps_plv);
                        plv_ts = plv_ts - mean(plv_ts(plv_bl), 'omitnan');
                        grand.(pname).(stage_names{s_i}).(bt).data(end+1,:) = plv_ts;
                        grand.(pname).(stage_names{s_i}).(bt).subj(end+1,1) = participant;
                    end
                end
            end
        end
    end

    % Per-subject ERP figures
    plot_subject_erps(subj, participant, is_cohort1, has_P, EEGp, EEGp_theta, ...
        stage_table, fcz_idx, par_idx, acc_idx, bl_mask, ...
        stage_names, STAGE_COLORS, LINE_STYLES, BTYPE_LABELS, ...
        ERP_plot_window, figure_output_folder);

    fprintf('  Mapped %d/%d trials to EEG epochs\n', sum(~isnan(trial2epoch)), n_beh_trials);
    clear EEGp EEGp_theta EEGp_phase

end % participant loop

% =========================================================================
%% SAVE OUTPUTS
% =========================================================================
assert(~isempty(t_ax), 't_ax is empty — participant loop produced no data');

group_table_ERPfigs = vertcat(all_stage_tables{:});
group_table_ERPfigs.subj       = categorical(group_table_ERPfigs.subj);
group_table_ERPfigs.block_type = categorical(group_table_ERPfigs.block_type);
group_table_ERPfigs.stage      = categorical(group_table_ERPfigs.stage, {'LN','LE','RN','RE'}, 'Ordinal',true);

% Stabilize types in all_trials_table for CSV export
if ismember('stage', all_trials_table.Properties.VariableNames)
    all_trials_table.stage = string(all_trials_table.stage);
    all_trials_table.stage(all_trials_table.stage == "" | ...
        all_trials_table.stage == "NaN") = missing;
end

if save_tables
    % ERP/grand-figure table (from stage_table) — used only by THIS script
    % for plotting. NOT the table consumed by the merge/RQ pipeline.
    group_table = group_table_ERPfigs; %#ok<NASGU> saved under its own name below
    save(fullfile(epoch_file_folder,'group_table_v8_ERPfigs.mat'), 'group_table');
    clear group_table
    save(fullfile(epoch_file_folder,'grand_v8.mat'), 'grand', 't_ax', 'SUBGROUPS');

    % [PIPELINE HANDOFF] Single-trial feature table is the one the
    % merge script and RQ script actually need: it has FCz_neg_peak_amp,
    % FRN_amp, RewP_amp, P300_amp, Theta_amp, PLV_fp, PLV_fs,
    % PLV_fp_pairwise, PLV_fs_pairwise, block_type, stage, correct,
    % false_fb, confidence, block_number, trial_in_block, subj — all as
    % per-trial scalar columns. Renamed to 'group_table' on disk so the
    % merge script's `load(..., 'group_table')` calls work unmodified.
    % Filename matches the merge script's KH branch:
    %   load(fullfile(epoch_file_folder,'group_stage_table_features.mat'), 'group_table')
    %
    % IMPORTANT: the merge script does
    %   vars = group_table_KH.Properties.VariableNames;
    %   group_table_RR = group_table_RR(:, vars);
    % which means EVERY column saved here must also exist in the RR
    % cohort's table. The RR cohort's pipeline does not store waveform
    % cell arrays, so waveform columns are dropped from this handoff
    % table (they remain in group_table_all_trials_v8.mat below if you
    % need them for KH-only analyses).
    waveform_vars = {'FCzCz_waveform','P300_waveform','Theta_waveform'};
    handoff_vars  = setdiff(all_trials_table.Properties.VariableNames, waveform_vars);
    group_table   = all_trials_table(:, handoff_vars); %#ok<NASGU> intentional rename for handoff
    save(fullfile(epoch_file_folder,'group_stage_table_features.mat'), 'group_table');
    clear group_table

    % Keep the full single-trial table (with waveforms) under its own name
    csv_vars = setdiff(all_trials_table.Properties.VariableNames, waveform_vars);
    save(fullfile(epoch_file_folder,'group_table_all_trials_v8.mat'), 'all_trials_table', 't_ax');
    writetable(all_trials_table(:, csv_vars), ...
        fullfile(epoch_file_folder, 'all_trials_features_v8.csv'));

    clear group_table  % avoid confusing the ERP-figures group_table below

    fprintf('\nSaved: group_table_v8_ERPfigs.mat, grand_v8.mat\n');
    fprintf('       group_stage_table_features.mat   <- feeds merge script (KH cohort)\n');
    fprintf('       group_table_all_trials_v8.mat, all_trials_features_v8.csv\n');
end

fprintf('\nSingle-trial summary:\n');
fprintf('  Total behavioral trials: %d\n', height(all_trials_table));
fprintf('  With EEG epoch:          %d\n', sum(all_trials_table.has_eeg_epoch));
fprintf('  Without EEG epoch:       %d\n', sum(~all_trials_table.has_eeg_epoch));
valid_stage = all_trials_table.stage ~= "" & ~ismissing(all_trials_table.stage);
fprintf('  With stage assignment:   %d\n', sum(valid_stage));

% =========================================================================
%% GRAND AVERAGE FIGURES
% =========================================================================
in_win = t_ax >= ERP_plot_window(1) & t_ax <= ERP_plot_window(2);

fig_ga = figure('Position',[50 50 1400 560]);
sgtitle('Grand average — FCz ERP by stage and block type (v8)');
for s_i = 1:4
    for oc_i = 1:2
        oc = {'correct','incorrect'};
        ax = subplot(2,4,(oc_i-1)*4+s_i);
        hold(ax,'on');
        title(ax, sprintf('%s | %s', stage_names{s_i}, oc{oc_i}));
        xline(ax,0,'k:','HandleVisibility','off');
        yline(ax,0,'k:','HandleVisibility','off');
        for bt_i = 1:2
            bt = BTYPE_LABELS{bt_i};
            d  = grand.FCz.(stage_names{s_i}).(bt).(oc{oc_i}).data;
            if isempty(d), continue; end
            plot_grand_ribbon(ax, t_ax, d, in_win, STAGE_COLORS(s_i,:), LINE_STYLES{bt_i}, bt);
        end
        set(ax,'YDir','reverse');
        legend(ax,'FontSize',8,'Box','off');
        xlabel(ax,'Time (ms)'); ylabel(ax,'\muV'); xlim(ax,ERP_plot_window);
    end
end
saveas(fig_ga, fullfile(figure_output_folder,'v8_Grand_FCz_stage.pdf'));
saveas(fig_ga, fullfile(figure_output_folder,'v8_Grand_FCz_stage.png'));

% =========================================================================
%% CLUSTER-BASED PERMUTATION TESTS (unchanged from v7 except pooling fix)
% =========================================================================
MIN_N_SUBJ = 8;
N_PERM     = 10000;
fprintf('\n\n=== CLUSTER PERMUTATION TESTS ===\n');

% Contrast A: Correct vs Incorrect per stage × block_type
for s_i = 1:4
    for bt_i = 1:2
        bt = BTYPE_LABELS{bt_i};
        d_cor = grand.FCz.(stage_names{s_i}).(bt).correct;
        d_inc = grand.FCz.(stage_names{s_i}).(bt).incorrect;
        shared = intersect(d_cor.subj, d_inc.subj);
        if numel(shared) < MIN_N_SUBJ
            fprintf('  SKIP %s %s CorVsInc — %d subjects\n', stage_names{s_i}, bt, numel(shared));
            continue
        end
        [~,ia] = ismember(shared, d_cor.subj);
        [~,ib] = ismember(shared, d_inc.subj);
        dat1 = d_cor.data(ia, in_win)';
        dat2 = d_inc.data(ib, in_win)';
        cname = sprintf('%s_%s_CorVsInc', stage_names{s_i}, bt);
        fprintf('  Running: %s (n=%d)\n', cname, numel(shared));
        fig_p = figure('Name',cname,'Position',[50 50 900 400]); hold on;
        title(sprintf('%s | %s — Correct vs Incorrect', stage_names{s_i}, bt));
        try
            [p,~,stat] = permutation_test_cluster_correction( ...
                t_ax(in_win), dat1, dat2, 'Correct','Incorrect','paired', ...
                [min([dat1(:);dat2(:)])-0.5, max([dat1(:);dat2(:)])+0.5], ...
                cname, [[0.10 0.60 0.10];[0.70 0.10 0.10]], [1 2]);
            perm_stats.(cname).p = p; perm_stats.(cname).stat = stat;
            perm_stats.(cname).n = numel(shared);
        catch ME
            warning('  %s FAILED: %s', cname, ME.message);
        end
        saveas(fig_p, fullfile(figure_output_folder, sprintf('v8_perm_%s.pdf', cname)));
        saveas(fig_p, fullfile(figure_output_folder, sprintf('v8_perm_%s.png', cname)));
    end
end

% Contrast B: D vs P blocks — BUG FIX v7-M4: correct averaging
for s_i = 1:4
    d_D_c  = grand.FCz.(stage_names{s_i}).D.correct;
    d_D_ic = grand.FCz.(stage_names{s_i}).D.incorrect;
    d_P_c  = grand.FCz.(stage_names{s_i}).P.correct;
    d_P_ic = grand.FCz.(stage_names{s_i}).P.incorrect;

    shared_all = intersect(intersect(d_D_c.subj, d_D_ic.subj), ...
                           intersect(d_P_c.subj,  d_P_ic.subj));
    if numel(shared_all) < MIN_N_SUBJ
        shared_all = intersect(d_D_c.subj, d_P_c.subj);  % fallback: correct only
    end
    if numel(shared_all) < MIN_N_SUBJ
        fprintf('  SKIP D vs P %s — %d subjects\n', stage_names{s_i}, numel(shared_all));
        continue
    end

    n_sh  = numel(shared_all);
    dat_D = nan(n_sh, sum(in_win));
    dat_P = nan(n_sh, sum(in_win));

    for si = 1:n_sh
        s_id = shared_all(si);
        % BUG FIX v7-M4: count found conditions and divide correctly
        for [dat_out, d_c, d_ic] = deal_pairs({dat_D, d_D_c, d_D_ic; dat_P, d_P_c, d_P_ic})
            % (helper below handles the accumulation)
        end
        % Inline the fix directly:
        n_D = 0; tmp_D = zeros(1,sum(in_win));
        row = find(d_D_c.subj==s_id,1);
        if ~isempty(row), tmp_D = tmp_D + d_D_c.data(row,in_win); n_D = n_D+1; end
        row = find(d_D_ic.subj==s_id,1);
        if ~isempty(row), tmp_D = tmp_D + d_D_ic.data(row,in_win); n_D = n_D+1; end
        if n_D > 0, dat_D(si,:) = tmp_D / n_D; end

        n_P = 0; tmp_P = zeros(1,sum(in_win));
        row = find(d_P_c.subj==s_id,1);
        if ~isempty(row), tmp_P = tmp_P + d_P_c.data(row,in_win); n_P = n_P+1; end
        row = find(d_P_ic.subj==s_id,1);
        if ~isempty(row), tmp_P = tmp_P + d_P_ic.data(row,in_win); n_P = n_P+1; end
        if n_P > 0, dat_P(si,:) = tmp_P / n_P; end
    end

    % Remove subjects where we got no data for either condition
    valid_rows = ~all(isnan(dat_D),2) & ~all(isnan(dat_P),2);
    dat_D = dat_D(valid_rows,:); dat_P = dat_P(valid_rows,:);
    if sum(valid_rows) < MIN_N_SUBJ
        fprintf('  SKIP D vs P %s after cleaning — %d subjects\n', stage_names{s_i}, sum(valid_rows));
        continue
    end

    cname = sprintf('%s_DvsP', stage_names{s_i});
    fprintf('  Running: %s (n=%d)\n', cname, sum(valid_rows));
    fig_p = figure('Name',cname,'Position',[50 50 900 400]); hold on;
    title(sprintf('%s — D vs P block type', stage_names{s_i}));
    try
        [p,~,stat] = permutation_test_cluster_correction( ...
            t_ax(in_win), dat_D', dat_P', ...
            'Deterministic','Probabilistic','paired', ...
            [min([dat_D(:);dat_P(:)])-0.5, max([dat_D(:);dat_P(:)])+0.5], ...
            cname, [[0.15 0.45 0.70];[0.80 0.30 0.10]], [1 2]);
        perm_stats.(cname).p = p; perm_stats.(cname).stat = stat;
        perm_stats.(cname).n = sum(valid_rows);
    catch ME
        warning('  %s FAILED: %s', cname, ME.message);
    end
    saveas(fig_p, fullfile(figure_output_folder, sprintf('v8_perm_%s.pdf', cname)));
    saveas(fig_p, fullfile(figure_output_folder, sprintf('v8_perm_%s.png', cname)));
end

% Contrast C: False incorrect vs True incorrect (P blocks)
for s_i = 1:4
    d_fi = grand.FCz.(stage_names{s_i}).P.false_inc;
    d_ti = grand.FCz.(stage_names{s_i}).P.incorrect;
    shared = intersect(d_fi.subj, d_ti.subj);
    if numel(shared) < MIN_N_SUBJ
        fprintf('  SKIP false vs true inc P %s — %d subjects\n', stage_names{s_i}, numel(shared));
        continue
    end
    [~,ia] = ismember(shared, d_fi.subj);
    [~,ib] = ismember(shared, d_ti.subj);
    dat1 = d_fi.data(ia, in_win)';
    dat2 = d_ti.data(ib, in_win)';
    cname = sprintf('%s_P_FalseInc_vs_TrueInc', stage_names{s_i});
    fprintf('  Running: %s (n=%d)\n', cname, numel(shared));
    fig_p = figure('Name',cname,'Position',[50 50 900 400]); hold on;
    title(sprintf('%s P — False Inc vs True Inc', stage_names{s_i}));
    try
        [p,~,stat] = permutation_test_cluster_correction( ...
            t_ax(in_win), dat1, dat2, ...
            'False incorrect','True incorrect','paired', ...
            [min([dat1(:);dat2(:)])-0.5, max([dat1(:);dat2(:)])+0.5], ...
            cname, [[0.90 0.50 0.00];[0.70 0.10 0.10]], [1 2]);
        perm_stats.(cname).p = p; perm_stats.(cname).stat = stat;
        perm_stats.(cname).n = numel(shared);
    catch ME
        warning('  %s FAILED: %s', cname, ME.message);
    end
    saveas(fig_p, fullfile(figure_output_folder, sprintf('v8_perm_%s.pdf', cname)));
    saveas(fig_p, fullfile(figure_output_folder, sprintf('v8_perm_%s.png', cname)));
end

save(fullfile(epoch_file_folder,'perm_stats_v8.mat'), 'perm_stats');

% =========================================================================
%% PERMUTATION RESULTS SUMMARY
% =========================================================================
fprintf('\n=== PERMUTATION RESULTS SUMMARY ===\n');
for cn = fieldnames(perm_stats)'
    ps_r = perm_stats.(cn{1});
    min_pos = NaN; min_neg = NaN;
    if isfield(ps_r,'stat') && isfield(ps_r.stat,'posclusters') && ~isempty(ps_r.stat.posclusters)
        min_pos = min([ps_r.stat.posclusters.prob]);
    end
    if isfield(ps_r,'stat') && isfield(ps_r.stat,'negclusters') && ~isempty(ps_r.stat.negclusters)
        min_neg = min([ps_r.stat.negclusters.prob]);
    end
    sig = ''; if (~isnan(min_pos)&&min_pos<0.025)||(~isnan(min_neg)&&min_neg<0.025), sig='  ***'; end
    fprintf('  %-45s  n=%2d  pos_p=%.3f  neg_p=%.3f%s\n', cn{1}, ps_r.n, min_pos, min_neg, sig);
end

% =========================================================================
%% PART 2: FLEXIBLE ERP AGGREGATION (from feature extraction script)
% Runs on all_trials_table which was built during the participant loop.
% =========================================================================
fprintf('\n\n=== PART 2: FLEXIBLE ERP AGGREGATION ===\n');

if ~exist('t_ax','var') || isempty(t_ax)
    error('t_ax missing — Part 2 cannot run standalone without loading grand_v8.mat first.');
end
time_vector = t_ax;

GROUPING_SCHEMES = {
    {'correct'}, ...
    {'block_type'}, ...
    {'block_type','correct'}, ...
    {'block_type','confidence','correct'}, ...
    {'stage'}, ...
    {'block_type','stage'}, ...
    {'block_type','stage','correct'}, ...
    {'block_type','false_fb'}, ...
    {'block_type','false_fb','stage'}
};

aggregated_erps_cell = {};

% Stabilize types
if ismember('subjID', all_trials_table.Properties.VariableNames)
    all_trials_table.subjID = string(all_trials_table.subjID);
end
if ismember('block_type', all_trials_table.Properties.VariableNames)
    all_trials_table.block_type = string(all_trials_table.block_type);
end

for scheme_idx = 1:numel(GROUPING_SCHEMES)
    grouping_vars = GROUPING_SCHEMES{scheme_idx};
    fprintf('  Grouping: %s\n', strjoin(grouping_vars,' x '));

    missing_vars = setdiff(grouping_vars, all_trials_table.Properties.VariableNames);
    if ~isempty(missing_vars)
        warning('  Skipping — missing: %s', strjoin(missing_vars,', '));
        continue;
    end

    valid = all_trials_table.has_eeg_epoch & ...
            ~cellfun(@isempty, all_trials_table.FCzCz_waveform);
    if ismember('stage',      grouping_vars), valid = valid & ~ismissing(all_trials_table.stage); end
    if ismember('block_type', grouping_vars), valid = valid & all_trials_table.block_type ~= "" & all_trials_table.block_type ~= "unknown"; end
    if ismember('confidence', grouping_vars), valid = valid & ~isnan(all_trials_table.confidence); end
    if ismember('correct',    grouping_vars), valid = valid & ~isnan(all_trials_table.correct); end
    if ismember('false_fb',   grouping_vars), valid = valid & ~isnan(all_trials_table.false_fb); end

    if sum(valid) == 0, continue; end
    combos = unique(all_trials_table(valid, grouping_vars), 'rows');

    for comb_idx = 1:height(combos)
        trial_mask = valid;
        label_parts = strings(1, numel(grouping_vars));
        for vi = 1:numel(grouping_vars)
            vn  = grouping_vars{vi};
            vval= combos.(vn)(comb_idx);
            trial_mask = trial_mask & compare_table_value(all_trials_table.(vn), vval);
            label_parts(vi) = sprintf('%s=%s', vn, value_to_label(vval));
        end
        if sum(trial_mask) < 3, continue; end
        cond_label = strjoin(label_parts, '_');
        agg = aggregate_erp_data(all_trials_table(trial_mask,:), time_vector, FRN_win_group, RewP_win_group);
        agg.condition_label = char(cond_label);
        agg.grouping_scheme = scheme_idx;
        agg.grouping_vars   = grouping_vars;
        agg.grouping_values = combos(comb_idx,:);
        agg.time_vector     = time_vector;
        if ismember('false_fb', grouping_vars)
            ff = combos.false_fb(comb_idx);
            if ff==1, agg.fb_type='falseFB'; elseif ff==0, agg.fb_type='trueFB'; else, agg.fb_type='unknownFB'; end
        else
            agg.fb_type = 'standard';
        end
        aggregated_erps_cell{end+1,1} = agg;
        fprintf('    %s: n=%d trials\n', cond_label, sum(trial_mask));
    end
end

if isempty(aggregated_erps_cell)
    aggregated_erps = struct([]);
    warning('No aggregated ERP conditions created.');
else
    aggregated_erps = standardize_struct_array(aggregated_erps_cell);
end

save(fullfile(epoch_file_folder,'aggregated_erps_v8.mat'), 'aggregated_erps');
if ~isempty(aggregated_erps)
    export_aggregated_summary(aggregated_erps, epoch_file_folder);
end
fprintf('\n✓ Saved aggregated_erps_v8.mat with %d conditions\n', numel(aggregated_erps));

% =========================================================================
%% SUBGROUP RE-PLOTTING
% =========================================================================
subgroup_plot_grand(grand, t_ax, 'cohort1 auditory FB',          [3:8],    stage_names, STAGE_COLORS, ERP_plot_window);
subgroup_plot_grand(grand, t_ax, 'cohort2a timing slightly off', [10],     stage_names, STAGE_COLORS, ERP_plot_window);
subgroup_plot_grand(grand, t_ax, 'cohort2b feedback trigger',    [11:18],  stage_names, STAGE_COLORS, ERP_plot_window);
subgroup_plot_grand(grand, t_ax, 'cohort2bv2 better impedance',  [17:20],  stage_names, STAGE_COLORS, ERP_plot_window);
subgroup_plot_grand(grand, t_ax, 'cohort LaMB all visual',       [10:20],  stage_names, STAGE_COLORS, ERP_plot_window);

fprintf('\nAll done.\n');


% =========================================================================
%% LOCAL FUNCTIONS
% =========================================================================

function plot_grand_ribbon(ax, t_ax, data_mat, in_win, clr, ls, lbl)
    if isempty(data_mat); return; end
    nS = size(data_mat,1);
    mn = mean(data_mat(:,in_win),1,'omitnan');
    se = std(data_mat(:,in_win),0,1,'omitnan')/sqrt(nS);
    fill(ax,[t_ax(in_win),fliplr(t_ax(in_win))],[mn+se,fliplr(mn-se)], ...
        clr,'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
    plot(ax,t_ax(in_win),mn,'Color',clr,'LineWidth',2,'LineStyle',ls, ...
        'DisplayName',sprintf('%s (n=%d subj)',lbl,nS));
end

function h = plot_erp_trace(ax, times, data_mat, clr, lbl, xlims, ls)
    if nargin < 7; ls = '-'; end
    n = size(data_mat,1);
    if n == 0; h=[]; return; end
    in_win = times >= xlims(1) & times <= xlims(2);
    mn = mean(data_mat(:,in_win),1,'omitnan');
    se = std(data_mat(:,in_win),0,1,'omitnan')/sqrt(n);
    t  = times(in_win);
    fill(ax,[t,fliplr(t)],[mn+se,fliplr(mn-se)],clr,'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
    h = plot(ax,t,mn,'Color',clr,'LineWidth',1.5,'LineStyle',ls, ...
        'DisplayName',sprintf('%s (%d trials)',lbl,n));  % "trials" explicit
end

function data_bl = extract_epochs(EEGp, ch_idx, epoch_col, bl_mask)
    ep = clean_epochs_checked(epoch_col, EEGp);
    if isempty(ep); data_bl = zeros(0,size(EEGp.data,2)); return; end
    if isscalar(ch_idx)
        raw = squeeze(EEGp.data(ch_idx,:,ep))';
    else
        raw = squeeze(mean(EEGp.data(ch_idx,:,ep),1,'omitnan'))';
    end
    if isvector(raw) && numel(ep)==1; raw=raw(:)'; end
    bl      = mean(raw(:,bl_mask),2,'omitnan');
    data_bl = raw - bl;
end

function theta_mat = compute_theta_envelope(EEGp_theta, ch_idx, ep, bl_mask)
    if isempty(ep); theta_mat=zeros(0,size(EEGp_theta.data,2)); return; end
    if isscalar(ch_idx)
        raw = squeeze(EEGp_theta.data(ch_idx,:,ep))';
    else
        raw = squeeze(mean(EEGp_theta.data(ch_idx,:,ep),1,'omitnan'))';
    end
    if isvector(raw)&&numel(ep)==1; raw=raw(:)'; end
    n_ep=size(raw,1); theta_mat=nan(n_ep,size(raw,2));
    for k=1:n_ep
        env=abs(hilbert(double(raw(k,:))));
        theta_mat(k,:)=env-mean(env(bl_mask),'omitnan');
    end
end

function [plv_ts,n_used] = compute_cross_trial_plv(EEGp_phase, ref_idx, tgt_idx, ep)
    n_t=size(EEGp_phase.data,2); n_ep=numel(ep);
    if n_ep==0; plv_ts=zeros(1,n_t); n_used=0; return; end
    if isscalar(ref_idx)
        phi_ref=squeeze(double(EEGp_phase.data(ref_idx,:,ep)))';
    else
        phi_ref=squeeze(angle(mean(exp(1i*double(EEGp_phase.data(ref_idx,:,ep))),1,'omitnan')))';
    end
    if isscalar(tgt_idx)
        phi_tgt=squeeze(double(EEGp_phase.data(tgt_idx,:,ep)))';
    else
        phi_tgt=squeeze(angle(mean(exp(1i*double(EEGp_phase.data(tgt_idx,:,ep))),1,'omitnan')))';
    end
    if isvector(phi_ref)&&n_ep==1; phi_ref=phi_ref(:)'; end
    if isvector(phi_tgt)&&n_ep==1; phi_tgt=phi_tgt(:)'; end
    plv_ts=abs(mean(exp(1i*(phi_ref-phi_tgt)),1,'omitnan'));
    n_used=n_ep;
end

function ep=clean_epochs_checked(epoch_col, EEGp)
    ep=epoch_col(~isnan(epoch_col) & epoch_col>=1 & epoch_col<=EEGp.trials);
end

function idx=safe_chan(EEGp, labels)
    idx=find(ismember(lower({EEGp.chanlocs.labels}),lower(labels)));
end

function plot_subject_erps(subj, participant, is_cohort1, has_P, EEGp, EEGp_theta, ...
    stage_table, fcz_idx, par_idx, acc_idx, bl_mask, ...
    stage_names, STAGE_COLORS, LINE_STYLES, BTYPE_LABELS, ERP_plot_window, outdir)

    nFigRows = 1 + double(has_P);
    fig1 = figure('Position',[50 50 1400 300*nFigRows],'Visible','off');
    sgtitle(sprintf('%s — FCz ERP by stage', subj),'Interpreter','none');
    for bt_i=1:numel(BTYPE_LABELS)
        bt=BTYPE_LABELS{bt_i};
        if strcmp(bt,'P')&&~has_P; continue; end
        for s_i=1:4
            ax=subplot(nFigRows,4,(bt_i-1)*4+s_i); hold(ax,'on');
            title(ax,sprintf('%s | %s',stage_names{s_i},bt));
            xline(ax,0,'k:','HandleVisibility','off');
            yline(ax,0,'k:','HandleVisibility','off');
            if ~isempty(fcz_idx)
                bm=stage_table.block_type==bt & stage_table.stage==stage_names{s_i};
                dat=extract_epochs(EEGp,fcz_idx,stage_table.epoch(bm&stage_table.correct==1&~stage_table.false_fb),bl_mask);
                plot_erp_trace(ax,EEGp.times,dat,[0.1 0.6 0.1],'Correct',ERP_plot_window);
                dat=extract_epochs(EEGp,fcz_idx,stage_table.epoch(bm&stage_table.correct==0&~stage_table.false_fb),bl_mask);
                plot_erp_trace(ax,EEGp.times,dat,[0.7 0.1 0.1],'Incorrect',ERP_plot_window);
                if strcmp(bt,'P')&&~is_cohort1
                    dat=extract_epochs(EEGp,fcz_idx,stage_table.epoch(bm&stage_table.false_fb&stage_table.fb_shown_correct==1),bl_mask);
                    plot_erp_trace(ax,EEGp.times,dat,[0.0 0.4 0.9],'False cor',ERP_plot_window);
                    dat=extract_epochs(EEGp,fcz_idx,stage_table.epoch(bm&stage_table.false_fb&stage_table.fb_shown_correct==0),bl_mask);
                    plot_erp_trace(ax,EEGp.times,dat,[0.9 0.5 0.0],'False inc',ERP_plot_window);
                end
            end
            set(ax,'YDir','reverse');
            legend(ax,'FontSize',7,'Box','off'); xlabel(ax,'Time (ms)'); ylabel(ax,'\muV'); xlim(ax,ERP_plot_window);
        end
    end
    saveas(fig1,fullfile(outdir,sprintf('%s_FCz_stage_v8.pdf',subj)));
    saveas(fig1,fullfile(outdir,sprintf('%s_FCz_stage_v8.png',subj)));
    close(fig1);
end

function subgroup_plot_grand(grand, t_ax, label, subj_ids, stage_names, STAGE_COLORS, ERP_plot_window)
    in_win=t_ax>=ERP_plot_window(1)&t_ax<=ERP_plot_window(2);
    BTYPE_LABELS={'D','P'}; LINE_STYLES={'-','--'};
    fig=figure('Position',[50 50 1400 560]);
    sgtitle(sprintf('Grand average FCz — %s',label));
    for s_i=1:4
        for oc_i=1:2
            oc={'correct','incorrect'};
            ax=subplot(2,4,(oc_i-1)*4+s_i); hold(ax,'on');
            title(ax,sprintf('%s | %s',stage_names{s_i},oc{oc_i}));
            xline(ax,0,'k:','HandleVisibility','off'); yline(ax,0,'k:','HandleVisibility','off');
            for bt_i=1:2
                bt=BTYPE_LABELS{bt_i};
                g=grand.FCz.(stage_names{s_i}).(bt).(oc{oc_i});
                if isempty(g.data); continue; end
                mask=ismember(g.subj,subj_ids);
                dat=g.data(mask,:);
                if isempty(dat); continue; end
                nS=size(dat,1);
                mn=mean(dat(:,in_win),1,'omitnan');
                se=std(dat(:,in_win),0,1,'omitnan')/sqrt(nS);
                fill(ax,[t_ax(in_win),fliplr(t_ax(in_win))],[mn+se,fliplr(mn-se)], ...
                    STAGE_COLORS(s_i,:),'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
                plot(ax,t_ax(in_win),mn,'Color',STAGE_COLORS(s_i,:),'LineWidth',2, ...
                    'LineStyle',LINE_STYLES{bt_i},'DisplayName',sprintf('%s (n=%d subj)',bt,nS));
            end
            set(ax,'YDir','reverse');
            legend(ax,'FontSize',8,'Box','off'); xlabel(ax,'Time (ms)'); ylabel(ax,'\muV'); xlim(ax,ERP_plot_window);
        end
    end
end

% ---- Part 2 helpers (unchanged from feature extraction script) ----------

function agg = aggregate_erp_data(trial_subset, time_vector, FRN_win, RewP_win)
    agg.n_trials   = height(trial_subset);
    agg.n_subjects = numel(unique(string(trial_subset.subjID)));
    agg.subjects   = unique(string(trial_subset.subjID))';
    frn_mask  = time_vector>=FRN_win(1)  & time_vector<=FRN_win(2);
    rewp_mask = time_vector>=RewP_win(1) & time_vector<=RewP_win(2);
    waves = trial_subset.FCzCz_waveform(~cellfun(@isempty,trial_subset.FCzCz_waveform));
    if ~isempty(waves)
        wf = vertcat(waves{:});
        agg.FCzCz_mean        = mean(wf,1,'omitnan');
        agg.FCzCz_sem         = std(wf,0,1,'omitnan')./sqrt(size(wf,1));
        agg.FCzCz_all_trials  = wf;
        fv=agg.FCzCz_mean(frn_mask); ft=time_vector(frn_mask);
        if any(~isnan(fv)); [agg.FRN_amp,ix]=min(fv,[],'omitnan'); agg.FRN_lat=ft(ix); agg.FRN_mean_amp=mean(fv,'omitnan');
        else; agg.FRN_amp=NaN; agg.FRN_lat=NaN; agg.FRN_mean_amp=NaN; end
        rv=agg.FCzCz_mean(rewp_mask); rt=time_vector(rewp_mask);
        if any(~isnan(rv)); [agg.RewP_peak_amp,ix]=max(rv,[],'omitnan'); agg.RewP_peak_lat=rt(ix); agg.RewP_mean_amp=mean(rv,'omitnan');
        else; agg.RewP_peak_amp=NaN; agg.RewP_peak_lat=NaN; agg.RewP_mean_amp=NaN; end
    else
        agg.FCzCz_mean=[];  agg.FCzCz_sem=[];  agg.FCzCz_all_trials=[];
        agg.FRN_amp=NaN; agg.FRN_lat=NaN; agg.FRN_mean_amp=NaN;
        agg.RewP_peak_amp=NaN; agg.RewP_peak_lat=NaN; agg.RewP_mean_amp=NaN;
    end
    for wtype = {'P300','Theta'}
        wf_col = [wtype{1},'_waveform'];
        waves  = trial_subset.(wf_col)(~cellfun(@isempty,trial_subset.(wf_col)));
        if ~isempty(waves)
            wf = vertcat(waves{:});
            agg.([wtype{1},'_mean_waveform']) = mean(wf,1,'omitnan');
            agg.([wtype{1},'_sem_waveform'])  = std(wf,0,1,'omitnan')./sqrt(size(wf,1));
        else
            agg.([wtype{1},'_mean_waveform']) = [];
            agg.([wtype{1},'_sem_waveform'])  = [];
        end
    end
    for fn = {'N2_amp','N2_norm','N2_lat','FCz_neg_peak_amp','FCz_neg_peak_norm', ...
              'FCz_neg_peak_lat','P300_amp','P300_norm','P300_peak_lat','Theta_amp'}
        f = fn{1};
        agg.([f,'_mean'])=NaN; agg.([f,'_std'])=NaN; agg.([f,'_sem'])=NaN;
        if ismember(f,trial_subset.Properties.VariableNames)
            vals=double(trial_subset.(f)); nv=sum(~isnan(vals));
            if nv>0
                agg.([f,'_mean'])=mean(vals,'omitnan');
                agg.([f,'_std'])=std(vals,'omitnan');
                if nv>1, agg.([f,'_sem'])=std(vals,'omitnan')/sqrt(nv); end
            end
        end
    end
end

function export_aggregated_summary(aggregated_erps, outdir)
    n=numel(aggregated_erps);
    condition_labels=cell(n,1); n_trials=zeros(n,1); n_subjects=zeros(n,1);
    N2_amp_mean=nan(n,1); N2_amp_sem=nan(n,1);
    FCz_neg_peak_amp_mean=nan(n,1); FCz_neg_peak_amp_sem=nan(n,1);
    P300_amp_mean=nan(n,1); P300_amp_sem=nan(n,1);
    Theta_amp_mean=nan(n,1); Theta_amp_sem=nan(n,1);
    FRN_amp=nan(n,1); FRN_lat=nan(n,1); FRN_mean_amp=nan(n,1);
    RewP_peak_amp=nan(n,1); RewP_peak_lat=nan(n,1); RewP_mean_amp=nan(n,1);
    fb_type=cell(n,1);
    for i=1:n
        condition_labels{i}=aggregated_erps(i).condition_label;
        n_trials(i)=aggregated_erps(i).n_trials;
        n_subjects(i)=aggregated_erps(i).n_subjects;
        for [f,v]=deal({'N2_amp_mean','N2_amp_sem','FCz_neg_peak_amp_mean','FCz_neg_peak_amp_sem', ...
                        'P300_amp_mean','P300_amp_sem','Theta_amp_mean','Theta_amp_sem', ...
                        'FRN_amp','FRN_lat','FRN_mean_amp','RewP_peak_amp','RewP_peak_lat','RewP_mean_amp'})
            if isfield(aggregated_erps(i),f{1}), eval([f{1},'(i)=aggregated_erps(i).',f{1},';']); end
        end
        fb_type{i}=ternary(isfield(aggregated_erps(i),'fb_type'), aggregated_erps(i).fb_type, 'standard');
    end
    T=table(condition_labels,n_trials,n_subjects, ...
        N2_amp_mean,N2_amp_sem,FCz_neg_peak_amp_mean,FCz_neg_peak_amp_sem, ...
        P300_amp_mean,P300_amp_sem,Theta_amp_mean,Theta_amp_sem, ...
        FRN_amp,FRN_lat,FRN_mean_amp,RewP_peak_amp,RewP_peak_lat,RewP_mean_amp,fb_type);
    writetable(T,fullfile(outdir,'aggregated_erps_summary_v8.csv'));
    fprintf('  Exported aggregated_erps_summary_v8.csv\n');
end

function mask=compare_table_value(col,val)
    if iscell(col)||iscategorical(col), col=string(col); end
    if iscell(val)||iscategorical(val), val=string(val); end
    if isstring(col)||ischar(col), mask=string(col)==string(val);
    elseif isnumeric(col)||islogical(col), mask=col==double(val);
    else, try; mask=col==val; catch; mask=string(col)==string(val); end; end
    mask=mask(:);
end

function label=value_to_label(val)
    if iscell(val), val=val{1}; end
    if iscategorical(val), val=string(val); end
    if isstring(val)||ischar(val), label=char(string(val));
    elseif isnumeric(val)||islogical(val)
        if isnan(double(val)), label='NaN'; else, label=num2str(double(val)); end
    else, label=char(string(val)); end
    label=strrep(strrep(label,' ',''),'/','–');
end

function S=standardize_struct_array(C)
    if isempty(C); S=struct([]); return; end
    all_fields={};
    for i=1:numel(C), all_fields=union(all_fields,fieldnames(C{i})); end
    for i=1:numel(C)
        for f=setdiff(all_fields,fieldnames(C{i}))', C{i}.(f{1})=[]; end
        C{i}=orderfields(C{i},all_fields);
    end
    S=vertcat(C{:});
end

function out=ternary(cond,a,b)
    if cond, out=a; else, out=b; end
end