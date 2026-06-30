% =============================================================================
% OUTCOME ERP ANALYSIS v9c
%
% CHANGES vs v9b (uploaded document version):
%
%  BUG FIX A — Stage assignment column-name mismatch
%     assign_stages_preserve_LE_RN required 'block_number' and
%     'trial_in_block', but behav_table uses 'block' and 'trial'.
%     The function now accepts the actual column names as arguments and
%     maps them internally. The call site passes the correct names.
%
%  BUG FIX B — trial2epoch shape not guaranteed to be a column vector
%     For KH the cached file or KH_align_epochs_with_offset may return
%     a row vector. A single (:) reshape is applied immediately after
%     loading/computing so downstream code can always assume a column.
%
%  BUG FIX C — Double assignment of subj_features.epoch
%     The spine-build block set epoch correctly, then a stale line
%     'subj_features.epoch = trial2epoch(1:n_rows)' later overwrote it
%     unconditionally — and could crash if trial2epoch was shorter than
%     n_rows (for KH the mapping vector length is not guaranteed to equal
%     the behavioural row count). The stale line is removed; the
%     spine-block assignment is the only one.
%
%  BUG FIX D — Participant loop skipped first two subjects
%     valid_participants(3:end) was used — almost certainly a debugging
%     leftover. Restored to valid_participants.
%
%  BUG FIX E — beh.revTrial indexing in assign_stages_preserve_LE_RN
%     revTrial is indexed by sequential block position (1 = first real
%     block). The behav_table 'block' column may start at 1 or 2
%     depending on whether practice was included. We now pass the
%     per-subject revTrial vector directly to the function and index it
%     by sequential block position within that subject's data, not by
%     the raw block number from the table.
%
%  BUG FIX F — Grand average and ERP figure used stage_table but
%     stage_table.epoch was from define_trial_stages_v3, while
%     subj_features.epoch is the authoritative alignment. The grand
%     average now reads stage information from subj_features, not
%     from a separately-built stage_table. define_trial_stages_v3 is
%     removed from the pipeline. (The grand average loop is refactored
%     to use subj_features directly.)
%
%  BUG FIX G — has_P flag still read stage_table which no longer exists.
%     Now derived from subj_features.block_type directly.
%
%  BUG FIX H — RR_data_path had '/Data' appended in the document version
%     but the RR epoch/figure folders are set explicitly and don't depend
%     on RR_data_path; the path is corrected to match the v9b original.
%
%  NOTE on trial2epoch semantics (for reference):
%     trial2epoch is a column vector of length n_beh_trials.
%     trial2epoch(t) = e means behavioural trial t corresponds to EEG
%     epoch e. trial2epoch(t) = NaN means trial t has no matched epoch
%     (e.g. rejected by ICA, or outside the alignment window for KH).
%     For RR: direct sequential — trial t maps to epoch t (up to the
%     shorter of the two).
%     For KH: offset alignment — the EEG trigger stream and the
%     behavioural log may have a constant offset due to practice-block
%     triggers being included in the EEG file. KH_align_epochs_with_offset
%     finds the best integer offset by maximising correct/incorrect
%     agreement between the behavioural log and EEG event markers.
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

addpath(fieldtrip_path); 
addpath(eeglab_path);    eeglab nogui;

KH_data_path = fullfile(base_path, 'Salient mod switch KH', 'Data');
RR_data_path = fullfile(base_path, 'Salient mod switch RR');   % BUG H: no /Data subfolder

addpath(genpath(KH_data_path));
addpath(genpath(RR_data_path));

load(fullfile(KH_data_path, 'all_trial_data.mat'));

% KH_behav_table is the single joint table (all subjects, both cohorts).
% It has a 'researcher' column ('KH' or 'RR') used to split cohorts.
KH_behav_file = fullfile(KH_data_path, 'behav_table.mat');
if exist(KH_behav_file, 'file')
    S = load(KH_behav_file, 'group_T');
    KH_behav_table = S.group_T;
else
    error('behav_table.mat not found at %s', KH_behav_file);
end

KH_epoch_file_folder    = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis', 'Epoched_data');
KH_figure_output_folder = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis', 'Figures', 'outcome_v9c_KH');

RR_epoch_file_folder    = fullfile(base_path, 'Salient mod switch RR', 'Results', 'EEG analysis', 'Epoched_data');
RR_figure_output_folder = fullfile(base_path, 'Salient mod switch RR', 'Results', 'EEG analysis', 'Outcome_ERP_analysis_figures_v9c');

if ~exist(KH_figure_output_folder, 'dir'), mkdir(KH_figure_output_folder); end
if ~exist(RR_figure_output_folder, 'dir'), mkdir(RR_figure_output_folder); end
if ~exist(KH_epoch_file_folder, 'dir'), error('KH epoch folder not found: %s', KH_epoch_file_folder); end
if ~exist(RR_epoch_file_folder, 'dir'), error('RR epoch folder not found: %s', RR_epoch_file_folder); end

% -------------------------------------------------------------------------
%% GLOBAL PARAMETERS
% -------------------------------------------------------------------------
save_tables = true;
ERP_plot_window = [-200 1000];
rm_baseline     = [-200 0];

N2_win    = [120 350];
FRN_win   = [250 350];
RewP_win  = [250 350];
P300_win  = [300 600];
Theta_win = [200 500];
PLV_win   = [200 400];
PLV_baseline = [-200 0];

MIN_TRIALS_PLV       = 5;
PLV_WINDOW_HALF      = 7;
MIN_TRIALS_PLV_WINDOW = 5;

stage_names  = {'LN','LE','RN','RE'};
BTYPE_LABELS = {'D','P'};
STAGE_COLORS = [0.12 0.62 0.47; 0.85 0.65 0.00; 0.80 0.27 0.13; 0.40 0.25 0.65];
LINE_STYLES  = {'-','--'};

% -------------------------------------------------------------------------
%% COHORT LIST
% -------------------------------------------------------------------------
cohort_names = {'KH', 'RR'};

KH_valid_participants = [3:8, 10:12, 14:23, 27];
RR_valid_participants = 1:15;

all_handoff_tables = {};
all_debug_rows     = {};

% =========================================================================
%% MAIN LOOP
% =========================================================================
for cohort_idx = 1:numel(cohort_names)

    cohort_name = cohort_names{cohort_idx};
    fprintf('\n\n############################################################\n');
    fprintf('RUNNING COHORT: %s\n', cohort_name);
    fprintf('############################################################\n');

    % ---------------------------------------------------------------
    % Cohort-specific settings
    % ---------------------------------------------------------------
    if strcmp(cohort_name, 'KH')
        valid_participants   = KH_valid_participants;   % BUG D: no (3:end)
        subj_prefix          = 'Ox';
        epoch_file_folder    = KH_epoch_file_folder;
        figure_output_folder = KH_figure_output_folder;
        trimmed_preferred    = true;
        mask = strcmp(KH_behav_table.researcher, 'KH');
        behav_table_this = KH_behav_table(mask, :);
        fcz_label    = 'FCz'; cz_label = 'Cz';
        par_channels = {'Pz','P1','P2'};
        acc_channels = {'FCz','Fz','AFz','F1','F2'};
        som_channels = {'C3','C4','CP3','CP1','C5','CP5'};
        handoff_file = 'group_stage_table_features_KH_v9c.mat';
        full_file    = 'group_table_all_trials_KH_v9c.mat';
        grand_file   = 'grand_KH_v9c.mat';

        % Column names in behav_table for block and trial position
        block_col       = 'block';
        trial_col       = 'trial';

    elseif strcmp(cohort_name, 'RR')
        valid_participants   = RR_valid_participants;
        subj_prefix          = 'Nc';
        epoch_file_folder    = RR_epoch_file_folder;
        figure_output_folder = RR_figure_output_folder;
        trimmed_preferred    = false;
        mask = strcmp(KH_behav_table.researcher, 'RR');
        behav_table_this = KH_behav_table(mask, :);
        fcz_label    = 'E11'; cz_label = 'E7';
        par_channels = {'E62','E67','E72'};
        acc_channels = {'E11','E6','E16'};
        som_channels = {'E36','E104','E41','E103'};
        handoff_file = 'group_stage_table_features_RR_v9c.mat';
        full_file    = 'group_table_all_trials_RR_v9c.mat';
        grand_file   = 'grand_RR_v9c.mat';

        block_col = 'block';
        trial_col = 'trial';
    end

    % ---------------------------------------------------------------
    % Grand-average containers for this cohort
    % BUG F: grand average is now accumulated from subj_features rows,
    % not from a separately-built stage_table. Containers are identical.
    % ---------------------------------------------------------------
    empty_container = struct('data', [], 'subj', []);
    grand = struct();
    for s = 1:4
        for bt = 1:2
            bt_s = BTYPE_LABELS{bt};
            grand.FCz.(stage_names{s}).(bt_s).correct    = empty_container;
            grand.FCz.(stage_names{s}).(bt_s).incorrect  = empty_container;
            grand.FCz.(stage_names{s}).(bt_s).false_cor  = empty_container;
            grand.FCz.(stage_names{s}).(bt_s).false_inc  = empty_container;
            grand.Par.(stage_names{s}).(bt_s).correct    = empty_container;
            grand.Par.(stage_names{s}).(bt_s).incorrect  = empty_container;
            grand.Theta.(stage_names{s}).(bt_s).correct  = empty_container;
            grand.Theta.(stage_names{s}).(bt_s).incorrect = empty_container;
            grand.PLV_fp.(stage_names{s}).(bt_s)         = empty_container;
            grand.PLV_fs.(stage_names{s}).(bt_s)         = empty_container;
        end
    end

    all_trials_table = table();
    t_ax = [];

    % ===================================================================
    %% PARTICIPANT LOOP
    % ===================================================================
    for participant = valid_participants   % BUG D fixed: full list

        subj = sprintf('%s%02d', subj_prefix, participant);
        fprintf('\n============ %s ============\n', subj);

        if ~isfield(all_trial_data, subj)
            warning('%s missing from all_trial_data. Skipping.', subj);
            continue;
        end

        beh = all_trial_data.(subj).trial_data;
        if strcmp(cohort_name, 'RR') && isfield(beh, 'structCode')
            beh.block_structure = beh.structCode;
        end

        % -----------------------------------------------------------
        % Flatten behavioural vectors from all_trial_data.
        % This is used ONLY for:
        %   (a) computing trial2epoch alignment for KH
        %   (b) providing beh_trueFB for RR (event-based)
        %   (c) providing revTrial to assign_stages_preserve_LE_RN
        % The actual behavioural data used in subj_features comes from
        % behav_table_this, not from these flattened vectors.
        % -----------------------------------------------------------
        num_blocks  = height(beh.correct);
        beh_correct = [];
        beh_trueFB  = [];

        if strcmp(cohort_name, 'KH') && num_blocks >= 6
            % Remove practice block from beh struct before flattening.
            % revTrial is a per-block vector; trimming removes block 1.
            beh.correct = beh.correct(2:end, :);
            if isfield(beh, 'confidence'), beh.confidence = beh.confidence(2:end, :); end
            if isfield(beh, 'trueFB'),     beh.trueFB     = beh.trueFB(2:end, :);     end
            if isfield(beh, 'revTrial'),   beh.revTrial   = beh.revTrial(2:end);      end
        end

        num_blocks = height(beh.correct);
        for b = 1:num_blocks
            beh_correct = [beh_correct, beh.correct(b, :)]; %#ok<AGROW>
            if isfield(beh, 'trueFB')
                beh_trueFB = [beh_trueFB, beh.trueFB(b, :)]; %#ok<AGROW>
            end
        end
        if isempty(beh_trueFB), beh_trueFB = ones(size(beh_correct)); end
        total_trials = numel(beh_correct);

        % -----------------------------------------------------------
        % Load EEG files
        % -----------------------------------------------------------
        if strcmp(cohort_name, 'KH')
            broadband_candidates = {sprintf('%s_outcome_trimmed.set', subj), sprintf('%s_outcome.set', subj)};
            theta_candidates     = {sprintf('%s_outcome_theta_trimmed.set', subj), sprintf('%s_outcome_theta.set', subj)};
            phase_candidates     = {sprintf('%s_outcome_phase_trimmed.set', subj), sprintf('%s_outcome_phase.set', subj)};
        elseif strcmp(cohort_name, 'RR')
            broadband_candidates = {sprintf('%s_outcome.set', subj), sprintf('%s_outcome_trimmed.set', subj)};
            theta_candidates     = {sprintf('%s_outcome_theta.set', subj), sprintf('%s_outcome_theta_trimmed.set', subj)};
            phase_candidates     = {sprintf('%s_outcome_phase.set', subj), sprintf('%s_outcome_phase_trimmed.set', subj)};
        end

        EEGp = [];
        for ci2 = 1:numel(broadband_candidates)
            f = fullfile(epoch_file_folder, broadband_candidates{ci2});
            if exist(f, 'file')
                EEGp = pop_loadset(broadband_candidates{ci2}, epoch_file_folder);
                fprintf('  broadband: %s (%d epochs)\n', broadband_candidates{ci2}, EEGp.trials);
                break;
            end
        end
        if isempty(EEGp)
            warning('broadband file missing for %s. Skipping.', subj);
            continue;
        end
        if isempty(t_ax), t_ax = EEGp.times; end

        EEGp_theta = [];
        for ci2 = 1:numel(theta_candidates)
            f = fullfile(epoch_file_folder, theta_candidates{ci2});
            if exist(f, 'file')
                EEGp_theta = pop_loadset(theta_candidates{ci2}, epoch_file_folder);
                fprintf('  theta: %s (%d epochs)\n', theta_candidates{ci2}, EEGp_theta.trials);
                break;
            end
        end
        if isempty(EEGp_theta)
            warning('theta file missing for %s; Theta_amp will be NaN.', subj);
        end

        EEGp_phase = [];
        for ci2 = 1:numel(phase_candidates)
            f = fullfile(epoch_file_folder, phase_candidates{ci2});
            if exist(f, 'file')
                EEGp_phase = pop_loadset(phase_candidates{ci2}, epoch_file_folder);
                fprintf('  phase: %s (%d epochs)\n', phase_candidates{ci2}, EEGp_phase.trials);
                break;
            end
        end
        if isempty(EEGp_phase)
            warning('phase file missing for %s; PLV columns will be NaN.', subj);
        end

        fprintf('  Behaviour trials: %d | EEG epochs: %d\n', total_trials, EEGp.trials);

        % -----------------------------------------------------------
        % ALIGNMENT — build trial2epoch
        %
        % trial2epoch is a column vector of length n_beh_trials.
        % trial2epoch(t) = e   : behavioural trial t → EEG epoch e
        % trial2epoch(t) = NaN : no matched epoch for trial t
        %
        % RR (direct sequential):
        %   EEG epochs are in the same order as behavioural trials.
        %   If there are fewer EEG epochs than behavioural trials
        %   (e.g. some epochs rejected before epoching), the remaining
        %   entries are NaN. If there are fewer behavioural trials than
        %   EEG epochs (should not happen), we take only the first
        %   n_beh_trials epochs.
        %
        % KH (offset alignment):
        %   The EEG file includes the practice block triggers, so the
        %   EEG epoch index for the first real behavioural trial is not
        %   necessarily 1. KH_align_epochs_with_offset searches for the
        %   constant integer offset that best aligns the correct/incorrect
        %   sequence in the behavioural log to the event markers in the
        %   EEG file, then returns a trial2epoch vector for all
        %   n_beh_trials real trials (NaN where no match found).
        %   The result is cached to disk as trial2epoch_out.
        % -----------------------------------------------------------
        if strcmp(cohort_name, 'RR')
            % Direct sequential: trial t → epoch t
            trial2epoch = nan(total_trials, 1);
            n_map = min(total_trials, EEGp.trials);
            trial2epoch(1:n_map) = (1:n_map)';
            fprintf('  Alignment: direct sequential (%d/%d mapped)\n', n_map, total_trials);

        elseif strcmp(cohort_name, 'KH')
            trial2epoch_file = fullfile(epoch_file_folder, sprintf('%s_trial2epoch.mat', subj));
            if exist(trial2epoch_file, 'file')
                S = load(trial2epoch_file);
                if isfield(S, 'trial2epoch_out')
                    trial2epoch = S.trial2epoch_out;
                elseif isfield(S, 'trial2epoch')
                    trial2epoch = S.trial2epoch;
                else
                    error('%s: unexpected field names in cached trial2epoch file.', trial2epoch_file);
                end
                fprintf('  trial2epoch: loaded from cache\n');
            else
                beh_cv = beh_correct(:);
                beh_cv = beh_cv(~isnan(beh_cv));
                [trial2epoch, diag_out] = KH_align_epochs_with_offset(EEGp, beh_cv);
                fprintf('  trial2epoch: %d/%d matched (%.1f%%), offset=%d\n', ...
                    diag_out.n_matched, diag_out.n_trials, diag_out.match_rate*100, diag_out.best_offset);
                trial2epoch_out = trial2epoch; %#ok<NASGU>
                save(trial2epoch_file, 'trial2epoch_out');
            end

            % BUG B FIX: guarantee column vector regardless of cache format
            trial2epoch = trial2epoch(:);

            % Pad or truncate to match total_trials so indexing is safe
            if numel(trial2epoch) < total_trials
                trial2epoch(end+1:total_trials) = NaN;
            elseif numel(trial2epoch) > total_trials
                trial2epoch = trial2epoch(1:total_trials);
            end
        end

        % -----------------------------------------------------------
        % Build behavioural spine from joint behav_table.
        % Behaviour is authoritative — keep ALL rows for this subject.
        % -----------------------------------------------------------
        if ~ismember('subjID', behav_table_this.Properties.VariableNames)
            error('behav_table_this lacks subjID column. Cannot build spine for %s.', subj);
        end

        subj_rows = string(behav_table_this.subjID) == string(subj);
        if ~any(subj_rows)
            warning('%s has no rows in behav_table_this. Skipping.', subj);
            continue;
        end

        subj_features = behav_table_this(subj_rows, :);
        n_rows = height(subj_features);
        fprintf('  Behavioural spine: %d rows\n', n_rows);

        % Core identifiers
        subj_features.subj_id = repmat(string(subj), n_rows, 1);
        subj_features.cohort  = repmat(string(cohort_name), n_rows, 1);
        subj_features.subj    = repmat(participant, n_rows, 1);

        % trial_continuous: sequential index 1..n_rows within this subject
        if ~ismember('trial_continuous', subj_features.Properties.VariableNames)
            subj_features.trial_continuous = (1:n_rows)';
        end

        % -----------------------------------------------------------
        % Assign epoch indices — SINGLE assignment (BUG C fixed).
        % trial2epoch has length total_trials (from beh flattening).
        % behav_table rows should equal total_trials for this subject.
        % If they differ (e.g. practice block included in one but not
        % the other), map only the overlapping range and NaN the rest.
        % -----------------------------------------------------------
        subj_features.epoch = nan(n_rows, 1);
        n_map = min(numel(trial2epoch), n_rows);
        subj_features.epoch(1:n_map) = trial2epoch(1:n_map);

        subj_features.has_eeg_epoch = ~isnan(subj_features.epoch) & ...
                                       subj_features.epoch >= 1 & ...
                                       subj_features.epoch <= EEGp.trials;

        % -----------------------------------------------------------
        % Ensure required columns exist before stage assignment
        % -----------------------------------------------------------
        if ~ismember('block_type', subj_features.Properties.VariableNames)
            subj_features.block_type = categorical(repmat(missing, n_rows, 1), {'D','P'});
        else
            subj_features.block_type = categorical(string(subj_features.block_type), {'D','P'});
        end

        if ~ismember('false_fb', subj_features.Properties.VariableNames)
            subj_features.false_fb = false(n_rows, 1);
        end
        if ~ismember('fb_shown_correct', subj_features.Properties.VariableNames)
            subj_features.fb_shown_correct = nan(n_rows, 1);
        end
        if ~ismember('confidence', subj_features.Properties.VariableNames)
            subj_features.confidence = nan(n_rows, 1);
        end

        % -----------------------------------------------------------
        % RR feedback-validity columns
        % (KH false_fb is already in behav_table)
        % -----------------------------------------------------------
        if strcmp(cohort_name, 'RR')
            fb_shown_correct_vec = false(total_trials, 1);
            for k = 1:min(total_trials, EEGp.trials)
                ep_events = {EEGp.event([EEGp.event.epoch] == k).type};
                if any(strcmp(ep_events, 'rewa'))
                    fb_shown_correct_vec(k) = true;
                end
            end
            false_fb_vec = ~logical(beh_trueFB(1:total_trials)');

            subj_features.fb_shown_correct = nan(n_rows, 1);
            subj_features.false_fb         = false(n_rows, 1);
            for r = 1:n_rows
                tc = subj_features.trial_continuous(r);
                if tc >= 1 && tc <= total_trials
                    subj_features.fb_shown_correct(r) = fb_shown_correct_vec(tc);
                    subj_features.false_fb(r)         = false_fb_vec(tc);
                end
            end
        end

        % -----------------------------------------------------------
        % STAGE ASSIGNMENT (BUG A, E fixed)
        %
        % Stage labels (LN, LE, RN, RE) are derived from position
        % within each block relative to the reversal trial.
        %
        % assign_stages_preserve_LE_RN now takes:
        %   block_col  — name of the column in subj_features that
        %                holds the block number (e.g. 'block')
        %   trial_col  — name of the column that holds trial position
        %                within the block (e.g. 'trial')
        %   rev_trials — vector where rev_trials(b) is the reversal
        %                trial number (within-block) for block b.
        %                This is beh.revTrial after practice trimming,
        %                so index 1 = first real block.
        %
        % Block numbers in subj_features.(block_col) may not start at 1
        % if practice was block 0 or 1 in the original log. We resolve
        % this by mapping block numbers to sequential positions (1-based)
        % before indexing into rev_trials. See function for details.
        % -----------------------------------------------------------
        rev_trials_vec = [];
        if isfield(beh, 'revTrial') && ~isempty(beh.revTrial)
            rev_trials_vec = beh.revTrial(:);  % column, already trimmed
        end

        subj_features = assign_stages_preserve_LE_RN(subj_features, ...
            block_col, trial_col, rev_trials_vec, stage_names);

        subj_features.subj_id = repmat(string(subj), n_rows, 1);
        subj_features.cohort  = repmat(string(cohort_name), n_rows, 1);
        subj_features.is_cohort1 = repmat(strcmp(cohort_name,'KH') && participant <= 8, n_rows, 1);

        % -----------------------------------------------------------
        % Channel indices and time masks
        % -----------------------------------------------------------
        fcz_idx = find(ismember(lower({EEGp.chanlocs.labels}), lower({fcz_label})));
        cz_idx  = find(ismember(lower({EEGp.chanlocs.labels}), lower({cz_label})));
        par_idx = find(ismember(lower({EEGp.chanlocs.labels}), lower(par_channels)));
        acc_idx = find(ismember(lower({EEGp.chanlocs.labels}), lower(acc_channels)));
        som_idx = find(ismember(lower({EEGp.chanlocs.labels}), lower(som_channels)));

        bl_mask   = EEGp.times >= rm_baseline(1)  & EEGp.times <= rm_baseline(2);
        n2_mask   = EEGp.times >= N2_win(1)       & EEGp.times <= N2_win(2);
        frn_mask  = EEGp.times >= FRN_win(1)      & EEGp.times <= FRN_win(2);
        rewp_mask = EEGp.times >= RewP_win(1)     & EEGp.times <= RewP_win(2);
        p300_mask = EEGp.times >= P300_win(1)     & EEGp.times <= P300_win(2);
        th_mask   = EEGp.times >= Theta_win(1)    & EEGp.times <= Theta_win(2);
        plv_mask  = EEGp.times >= PLV_win(1)      & EEGp.times <= PLV_win(2);
        plv_bl    = EEGp.times >= PLV_baseline(1) & EEGp.times <= PLV_baseline(2);

        if ~isempty(fcz_idx)
            bl_data   = squeeze(double(EEGp.data(fcz_idx, bl_mask, :)));
            bline_rms = rms(bl_data(:), 'omitnan');
        else
            bline_rms = NaN;
        end

        % -----------------------------------------------------------
        % Initialise EEG feature columns
        % -----------------------------------------------------------
        subj_features.baseline_rms = repmat(bline_rms, n_rows, 1);

        subj_features.N2_amp   = nan(n_rows, 1);
        subj_features.N2_lat   = nan(n_rows, 1);
        subj_features.N2_norm  = nan(n_rows, 1);

        subj_features.FRN_mean_amp  = nan(n_rows, 1);
        subj_features.FRN_mean_norm = nan(n_rows, 1);
        subj_features.FRN_peak_amp  = nan(n_rows, 1);
        subj_features.FRN_peak_lat  = nan(n_rows, 1);
        subj_features.FRN_peak_norm = nan(n_rows, 1);
        subj_features.FRN_excluded  = false(n_rows, 1);

        subj_features.RewP_mean_amp  = nan(n_rows, 1);
        subj_features.RewP_mean_norm = nan(n_rows, 1);
        subj_features.RewP_peak_amp  = nan(n_rows, 1);
        subj_features.RewP_peak_lat  = nan(n_rows, 1);
        subj_features.RewP_peak_norm = nan(n_rows, 1);
        subj_features.RewP_excluded  = false(n_rows, 1);

        subj_features.P300_amp      = nan(n_rows, 1);
        subj_features.P300_peak_lat = nan(n_rows, 1);
        subj_features.P300_norm     = nan(n_rows, 1);
        subj_features.Theta_amp     = nan(n_rows, 1);

        subj_features.PLV_fp           = nan(n_rows, 1);
        subj_features.PLV_fs           = nan(n_rows, 1);
        subj_features.PLV_fp_pairwise  = nan(n_rows, 1);
        subj_features.PLV_fs_pairwise  = nan(n_rows, 1);

        subj_features.FCzCz_waveform = repmat({[]}, n_rows, 1);
        subj_features.P300_waveform  = repmat({[]}, n_rows, 1);
        subj_features.Theta_waveform = repmat({[]}, n_rows, 1);

        % -----------------------------------------------------------
        % Per-trial EEG feature extraction
        % -----------------------------------------------------------
        for ti = 1:n_rows
            ep = subj_features.epoch(ti);
            if isnan(ep) || ep < 1 || ep > EEGp.trials, continue; end
            ep = round(ep);

            % ---- FCz+Cz: N2, FRN, RewP ----
            if ~isempty(fcz_idx) && ~isempty(cz_idx) && ~isnan(bline_rms)
                sig = mean(double(EEGp.data([fcz_idx cz_idx], :, ep)), 1, 'omitnan');
                sig = sig - mean(sig(bl_mask), 'omitnan');
                subj_features.FCzCz_waveform{ti} = sig;

                % N2 diagnostic (global min, not FRN measure)
                win_vals = sig(n2_mask); win_t = EEGp.times(n2_mask);
                if any(~isnan(win_vals))
                    [pk, ix] = min(win_vals, [], 'omitnan');
                    subj_features.N2_amp(ti) = pk;
                    subj_features.N2_lat(ti) = win_t(ix);
                    if bline_rms > 0, subj_features.N2_norm(ti) = pk / bline_rms; end
                end

                % FRN: mean always + genuine local minimum peak
                frn_vals = sig(frn_mask);
                frn_t    = EEGp.times(frn_mask);
                if any(~isnan(frn_vals))
                    subj_features.FRN_mean_amp(ti) = mean(frn_vals, 'omitnan');
                    if bline_rms > 0
                        subj_features.FRN_mean_norm(ti) = subj_features.FRN_mean_amp(ti) / bline_rms;
                    end
                    is_min = islocalmin(frn_vals);
                    if any(is_min)
                        cand_vals = frn_vals(is_min); cand_t = frn_t(is_min);
                        [pk_amp, ix] = min(cand_vals);
                        subj_features.FRN_peak_amp(ti) = pk_amp;
                        subj_features.FRN_peak_lat(ti) = cand_t(ix);
                        if bline_rms > 0
                            subj_features.FRN_peak_norm(ti) = pk_amp / bline_rms;
                        end
                    else
                        subj_features.FRN_excluded(ti) = true;
                    end
                else
                    subj_features.FRN_excluded(ti) = true;
                end

                % RewP: mean always + genuine local maximum peak
                rewp_vals = sig(rewp_mask);
                rewp_t    = EEGp.times(rewp_mask);
                if any(~isnan(rewp_vals))
                    subj_features.RewP_mean_amp(ti) = mean(rewp_vals, 'omitnan');
                    if bline_rms > 0
                        subj_features.RewP_mean_norm(ti) = subj_features.RewP_mean_amp(ti) / bline_rms;
                    end
                    is_max = islocalmax(rewp_vals);
                    if any(is_max)
                        cand_vals = rewp_vals(is_max); cand_t = rewp_t(is_max);
                        [pk_amp, ix] = max(cand_vals);
                        subj_features.RewP_peak_amp(ti) = pk_amp;
                        subj_features.RewP_peak_lat(ti) = cand_t(ix);
                        if bline_rms > 0
                            subj_features.RewP_peak_norm(ti) = pk_amp / bline_rms;
                        end
                    else
                        subj_features.RewP_excluded(ti) = true;
                    end
                else
                    subj_features.RewP_excluded(ti) = true;
                end
            end

            % ---- P300 ----
            if ~isempty(par_idx) && ~isnan(bline_rms)
                sig_p = mean(double(EEGp.data(par_idx, :, ep)), 1, 'omitnan');
                sig_p = sig_p - mean(sig_p(bl_mask), 'omitnan');
                subj_features.P300_waveform{ti} = sig_p;
                win_vals = sig_p(p300_mask); win_t = EEGp.times(p300_mask);
                if any(~isnan(win_vals))
                    [pk, ix] = max(win_vals, [], 'omitnan');
                    subj_features.P300_amp(ti)      = pk;
                    subj_features.P300_peak_lat(ti) = win_t(ix);
                    if bline_rms > 0, subj_features.P300_norm(ti) = pk / bline_rms; end
                end
            end

            % ---- Theta ----
            if ~isempty(EEGp_theta) && ~isempty(acc_idx) && ep <= EEGp_theta.trials
                sig_th = mean(double(EEGp_theta.data(acc_idx, :, ep)), 1, 'omitnan');
                env = abs(hilbert(sig_th));
                env = env - mean(env(bl_mask), 'omitnan');
                subj_features.Theta_amp(ti)      = mean(env(th_mask), 'omitnan');
                subj_features.Theta_waveform{ti} = env;
            end
        end

        % -----------------------------------------------------------
        % Sliding-window single-trial PLV
        % -----------------------------------------------------------
        if ~isempty(acc_idx) && ~isempty(EEGp_phase) && EEGp_phase.trials > 0
            row_idx  = (1:n_rows)';
            valid_ep = subj_features.has_eeg_epoch & ...
                       subj_features.epoch <= EEGp_phase.trials;

            % PLV bucketed by stage × block_type × correct.
            % Rows with missing stage are excluded (not in a stage window).
            stage_str = string(subj_features.stage);
            btype_str = string(subj_features.block_type);
            corr_str  = string(subj_features.correct);

            valid_bucket = valid_ep & ~ismissing(subj_features.stage);
            unique_buckets = unique([stage_str(valid_bucket), ...
                                     btype_str(valid_bucket), ...
                                     corr_str(valid_bucket)], 'rows');

            for ub = 1:size(unique_buckets, 1)
                bucket_mask = valid_bucket & ...
                    stage_str == unique_buckets(ub,1) & ...
                    btype_str == unique_buckets(ub,2) & ...
                    corr_str  == unique_buckets(ub,3);

                bucket_rows = row_idx(bucket_mask);
                if numel(bucket_rows) < MIN_TRIALS_PLV_WINDOW, continue; end

                [~, ord] = sort(subj_features.trial_continuous(bucket_rows));
                bucket_rows = bucket_rows(ord);
                n_bucket = numel(bucket_rows);

                for bi = 1:n_bucket
                    win_lo = max(1, bi - PLV_WINDOW_HALF);
                    win_hi = min(n_bucket, bi + PLV_WINDOW_HALF);
                    window_rows = bucket_rows(win_lo:win_hi);
                    if numel(window_rows) < MIN_TRIALS_PLV_WINDOW, continue; end

                    eps_window = subj_features.epoch(window_rows);
                    eps_window = eps_window(~isnan(eps_window) & eps_window>=1 & eps_window<=EEGp_phase.trials);
                    eps_window = round(eps_window);
                    if numel(eps_window) < MIN_TRIALS_PLV_WINDOW, continue; end

                    center_row = bucket_rows(bi);

                    if ~isempty(par_idx)
                        phi_ref = squeeze(angle(mean(exp(1i*double(EEGp_phase.data(acc_idx,:,eps_window))),1,'omitnan')))';
                        phi_tgt = squeeze(angle(mean(exp(1i*double(EEGp_phase.data(par_idx,:,eps_window))),1,'omitnan')))';
                        if isvector(phi_ref) && numel(eps_window)==1, phi_ref = phi_ref(:)'; end
                        if isvector(phi_tgt) && numel(eps_window)==1, phi_tgt = phi_tgt(:)'; end
                        plv_ts = abs(mean(exp(1i*(phi_ref-phi_tgt)),1,'omitnan'));
                        plv_ts = plv_ts - mean(plv_ts(plv_bl),'omitnan');
                        subj_features.PLV_fp(center_row)          = mean(plv_ts(plv_mask),'omitnan');
                        subj_features.PLV_fp_pairwise(center_row) = subj_features.PLV_fp(center_row);
                    end
                    if ~isempty(som_idx)
                        phi_ref = squeeze(angle(mean(exp(1i*double(EEGp_phase.data(acc_idx,:,eps_window))),1,'omitnan')))';
                        phi_tgt = squeeze(angle(mean(exp(1i*double(EEGp_phase.data(som_idx,:,eps_window))),1,'omitnan')))';
                        if isvector(phi_ref) && numel(eps_window)==1, phi_ref = phi_ref(:)'; end
                        if isvector(phi_tgt) && numel(eps_window)==1, phi_tgt = phi_tgt(:)'; end
                        plv_ts = abs(mean(exp(1i*(phi_ref-phi_tgt)),1,'omitnan'));
                        plv_ts = plv_ts - mean(plv_ts(plv_bl),'omitnan');
                        subj_features.PLV_fs(center_row)          = mean(plv_ts(plv_mask),'omitnan');
                        subj_features.PLV_fs_pairwise(center_row) = subj_features.PLV_fs(center_row);
                    end
                end
            end
        end

        all_trials_table = [all_trials_table; subj_features]; %#ok<AGROW>

        % -----------------------------------------------------------
        % Grand average accumulation (BUG F, G fixed)
        % Reads directly from subj_features; no separate stage_table.
        % has_P derived from subj_features.block_type (BUG G fixed).
        % -----------------------------------------------------------
        has_P = any(string(subj_features.block_type) == 'P');

        for s_i = 1:4
            for bt_i = 1:2
                bt = BTYPE_LABELS{bt_i};
                if strcmp(bt,'P') && ~has_P, continue; end

                base_mask = string(subj_features.block_type) == bt & ...
                            string(subj_features.stage) == stage_names{s_i};

                % Helper: extract epochs from subj_features mask
                get_eps = @(m) round(subj_features.epoch( ...
                    m & subj_features.has_eeg_epoch));
                get_eps_clamp = @(m, maxep) intersect(get_eps(m), 1:maxep);

                % FCz correct / incorrect
                for oc_i = 1:2
                    oc_val  = (oc_i == 2);   % 1=incorrect (oc_i==1), 2=correct
                    oc_name = {'incorrect','correct'};
                    corr_mask = base_mask & ...
                        subj_features.correct == oc_val & ...
                        ~subj_features.false_fb;

                    for gf = {'FCz','Par'}
                        if strcmp(gf{1},'FCz'), ch = fcz_idx; else, ch = par_idx; end
                        if isempty(ch), continue; end
                        eps_c = get_eps_clamp(corr_mask, EEGp.trials);
                        if isempty(eps_c), continue; end
                        if isscalar(ch)
                            raw = squeeze(double(EEGp.data(ch,:,eps_c)))';
                        else
                            raw = squeeze(mean(double(EEGp.data(ch,:,eps_c)),1,'omitnan'))';
                        end
                        if isvector(raw) && numel(eps_c)==1, raw = raw(:)'; end
                        dat = raw - mean(raw(:,bl_mask),2,'omitnan');
                        grand.(gf{1}).(stage_names{s_i}).(bt).(oc_name{oc_i}).data(end+1,:) = mean(dat,1,'omitnan');
                        grand.(gf{1}).(stage_names{s_i}).(bt).(oc_name{oc_i}).subj(end+1,1) = participant;
                    end
                end

                % False feedback (P blocks only)
                if strcmp(bt,'P') && ~isempty(fcz_idx)
                    ff_specs = { ...
                        'false_cor', base_mask & subj_features.false_fb & subj_features.fb_shown_correct==1; ...
                        'false_inc', base_mask & subj_features.false_fb & subj_features.fb_shown_correct==0};
                    for fi = 1:2
                        eps_f = get_eps_clamp(ff_specs{fi,2}, EEGp.trials);
                        if isempty(eps_f), continue; end
                        raw = squeeze(double(EEGp.data(fcz_idx,:,eps_f)))';
                        if isvector(raw) && numel(eps_f)==1, raw = raw(:)'; end
                        dat = raw - mean(raw(:,bl_mask),2,'omitnan');
                        grand.FCz.(stage_names{s_i}).(bt).(ff_specs{fi,1}).data(end+1,:) = mean(dat,1,'omitnan');
                        grand.FCz.(stage_names{s_i}).(bt).(ff_specs{fi,1}).subj(end+1,1) = participant;
                    end
                end

                % Theta
                if ~isempty(EEGp_theta) && ~isempty(acc_idx)
                    for oc_i = 1:2
                        oc_val  = (oc_i == 2);
                        oc_name = {'incorrect','correct'};
                        tmask = base_mask & subj_features.correct==oc_val & ~subj_features.false_fb;
                        eps_th = get_eps_clamp(tmask, EEGp_theta.trials);
                        if isempty(eps_th), continue; end
                        raw = squeeze(mean(double(EEGp_theta.data(acc_idx,:,eps_th)),1,'omitnan'))';
                        if isvector(raw) && numel(eps_th)==1, raw = raw(:)'; end
                        th_mat = nan(size(raw));
                        for kk = 1:size(raw,1)
                            env = abs(hilbert(double(raw(kk,:))));
                            th_mat(kk,:) = env - mean(env(bl_mask),'omitnan');
                        end
                        grand.Theta.(stage_names{s_i}).(bt).(oc_name{oc_i}).data(end+1,:) = mean(th_mat,1,'omitnan');
                        grand.Theta.(stage_names{s_i}).(bt).(oc_name{oc_i}).subj(end+1,1) = participant;
                    end
                end

                % PLV
                if ~isempty(EEGp_phase) && ~isempty(acc_idx)
                    plv_base_mask = base_mask & ~subj_features.false_fb;
                    eps_plv = get_eps_clamp(plv_base_mask, EEGp_phase.trials);
                    if numel(eps_plv) >= MIN_TRIALS_PLV
                        if ~isempty(par_idx)
                            phi_ref = squeeze(angle(mean(exp(1i*double(EEGp_phase.data(acc_idx,:,eps_plv))),1,'omitnan')))';
                            phi_tgt = squeeze(angle(mean(exp(1i*double(EEGp_phase.data(par_idx,:,eps_plv))),1,'omitnan')))';
                            plv_ts = abs(mean(exp(1i*(phi_ref-phi_tgt)),1,'omitnan'));
                            plv_ts = plv_ts - mean(plv_ts(plv_bl),'omitnan');
                            grand.PLV_fp.(stage_names{s_i}).(bt).data(end+1,:) = plv_ts;
                            grand.PLV_fp.(stage_names{s_i}).(bt).subj(end+1,1) = participant;
                        end
                        if ~isempty(som_idx)
                            phi_ref = squeeze(angle(mean(exp(1i*double(EEGp_phase.data(acc_idx,:,eps_plv))),1,'omitnan')))';
                            phi_tgt = squeeze(angle(mean(exp(1i*double(EEGp_phase.data(som_idx,:,eps_plv))),1,'omitnan')))';
                            plv_ts = abs(mean(exp(1i*(phi_ref-phi_tgt)),1,'omitnan'));
                            plv_ts = plv_ts - mean(plv_ts(plv_bl),'omitnan');
                            grand.PLV_fs.(stage_names{s_i}).(bt).data(end+1,:) = plv_ts;
                            grand.PLV_fs.(stage_names{s_i}).(bt).subj(end+1,1) = participant;
                        end
                    end
                end
            end
        end

        % -----------------------------------------------------------
        % Per-subject ERP figure (reads from subj_features)
        % -----------------------------------------------------------
        fig = figure('Position',[50 50 1400 560],'Visible','off');
        sgtitle(sprintf('%s — FCz ERP by stage', subj), 'Interpreter','none');
        for s_i = 1:4
            for oc_i = 1:2
                oc_val  = (oc_i == 2);
                oc_name = {'incorrect','correct'};
                ax = subplot(2,4,(oc_i-1)*4+s_i); hold(ax,'on');
                title(ax, sprintf('%s | %s', stage_names{s_i}, oc_name{oc_i}));
                xline(ax,0,'k:','HandleVisibility','off');
                yline(ax,0,'k:','HandleVisibility','off');
                for bt_i = 1:2
                    bt = BTYPE_LABELS{bt_i};
                    if strcmp(bt,'P') && ~has_P, continue; end
                    m = string(subj_features.block_type)==bt & ...
                        string(subj_features.stage)==stage_names{s_i} & ...
                        subj_features.correct==oc_val & ~subj_features.false_fb & ...
                        subj_features.has_eeg_epoch;
                    eps_plot = round(subj_features.epoch(m));
                    eps_plot = eps_plot(eps_plot>=1 & eps_plot<=EEGp.trials);
                    if isempty(eps_plot) || isempty(fcz_idx), continue; end
                    raw = squeeze(double(EEGp.data(fcz_idx,:,eps_plot)))';
                    if isvector(raw) && numel(eps_plot)==1, raw = raw(:)'; end
                    dat = raw - mean(raw(:,bl_mask),2,'omitnan');
                    in_win = EEGp.times >= ERP_plot_window(1) & EEGp.times <= ERP_plot_window(2);
                    mn = mean(dat(:,in_win),1,'omitnan');
                    se = std(dat(:,in_win),0,1,'omitnan') ./ sqrt(size(dat,1));
                    tt = EEGp.times(in_win);
                    fill(ax,[tt fliplr(tt)],[mn+se fliplr(mn-se)],STAGE_COLORS(s_i,:), ...
                        'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
                    plot(ax,tt,mn,'Color',STAGE_COLORS(s_i,:),'LineWidth',1.5, ...
                        'LineStyle',LINE_STYLES{bt_i}, ...
                        'DisplayName',sprintf('%s-%s (n=%d)',oc_name{oc_i},bt,size(dat,1)));
                end
                set(ax,'YDir','reverse'); xlabel(ax,'Time (ms)'); ylabel(ax,'\muV');
                xlim(ax,ERP_plot_window); legend(ax,'Box','off','FontSize',7);
            end
        end
        saveas(fig, fullfile(figure_output_folder, sprintf('%s_FCz_stage_v9c.pdf', subj)));
        % close(fig);

        % -----------------------------------------------------------
        % Debug row
        % -----------------------------------------------------------
        all_debug_rows(end+1,:) = {string(cohort_name), string(subj), participant, n_rows, ...
            sum(subj_features.has_eeg_epoch), ...
            sum(~cellfun(@isempty, subj_features.FCzCz_waveform)), ...
            sum(~isnan(subj_features.FRN_mean_amp)), ...
            sum(~subj_features.FRN_excluded), ...
            mean(~subj_features.FRN_excluded, 'omitnan'), bline_rms}; %#ok<AGROW>

        clear EEGp EEGp_theta EEGp_phase
    end % participant loop

    % -------------------------------------------------------------------
    %% Finalise cohort table: types + z-scoring
    % -------------------------------------------------------------------
    group_table = all_trials_table;
    group_table.subj_id    = categorical(string(group_table.subj_id));
    group_table.cohort     = categorical(string(group_table.cohort));
    if ismember('stage', group_table.Properties.VariableNames)
        group_table.stage = categorical(string(group_table.stage), {'LN','LE','RN','RE'}, 'Ordinal', false);
    end
    if ismember('block_type', group_table.Properties.VariableNames)
        group_table.block_type = categorical(string(group_table.block_type), {'D','P'});
    end

    features_to_zscore = {'N2_amp','FRN_mean_amp','FRN_peak_amp','RewP_mean_amp','RewP_peak_amp', ...
        'P300_amp','Theta_amp','PLV_fp','PLV_fs','PLV_fp_pairwise','PLV_fs_pairwise'};
    for f = 1:numel(features_to_zscore)
        fn = features_to_zscore{f};
        if ~ismember(fn, group_table.Properties.VariableNames), continue; end
        fn_z = [fn '_z'];
        group_table.(fn_z) = nan(height(group_table), 1);
        subjs = categories(group_table.subj_id);
        for si = 1:numel(subjs)
            mask2 = group_table.subj_id == subjs{si};
            vals = group_table.(fn)(mask2);
            mn = mean(vals,'omitnan'); sd = std(vals,'omitnan');
            if sd > 0
                group_table.(fn_z)(mask2) = (vals - mn) ./ sd;
            end
        end
    end

    if save_tables
        waveform_vars = {'FCzCz_waveform','P300_waveform','Theta_waveform'};
        handoff_vars  = setdiff(group_table.Properties.VariableNames, waveform_vars);
        group_table_save = group_table(:, handoff_vars); %#ok<NASGU>
        save(fullfile(epoch_file_folder, handoff_file), 'group_table_save');
        save(fullfile(epoch_file_folder, full_file), 'all_trials_table', 't_ax');
        save(fullfile(epoch_file_folder, grand_file), 'grand', 't_ax');
        fprintf('Saved %s, %s, %s\n', handoff_file, full_file, grand_file);
    end

    all_handoff_tables{end+1} = group_table; %#ok<AGROW>

    fprintf('\nCohort %s complete: %d trials, %d subjects\n', ...
        cohort_name, height(group_table), numel(unique(group_table.subj_id)));

end % cohort loop

% =========================================================================
%% COMBINE KH + RR HANDOFF TABLES
% =========================================================================
fprintf('\nCombining KH + RR handoff tables...\n');

waveform_vars = {'FCzCz_waveform','P300_waveform','Theta_waveform'};
all_vars = {};
for i = 1:numel(all_handoff_tables)
    these_vars = setdiff(all_handoff_tables{i}.Properties.VariableNames, waveform_vars);
    all_vars = union(all_vars, these_vars, 'stable');
end

tables_for_combine = cell(numel(all_handoff_tables), 1);
for i = 1:numel(all_handoff_tables)
    T = all_handoff_tables{i}(:, setdiff(all_handoff_tables{i}.Properties.VariableNames, waveform_vars));
    missing_cols = setdiff(all_vars, T.Properties.VariableNames, 'stable');
    for m = 1:numel(missing_cols)
        T.(missing_cols{m}) = nan(height(T), 1);
    end
    tables_for_combine{i} = T(:, all_vars);
end

group_table_combined = vertcat(tables_for_combine{:});
group_table_combined.subj_id = categorical(string(group_table_combined.subj_id));

combined_out_folder = KH_epoch_file_folder;
save(fullfile(combined_out_folder, 'group_feature_table_combined_v9c.mat'), 'group_table_combined');

fprintf('Saved combined table: %s\n', fullfile(combined_out_folder, 'group_feature_table_combined_v9c.mat'));
fprintf('Combined rows: %d, columns: %d\n', height(group_table_combined), width(group_table_combined));
fprintf('KH subjects: %d, RR subjects: %d\n', ...
    numel(unique(all_handoff_tables{1}.subj_id)), ...
    numel(unique(all_handoff_tables{2}.subj_id)));

% -------------------------------------------------------------------------
%% FRN exclusion summary
% -------------------------------------------------------------------------
fprintf('\n=== FRN PEAK EXCLUSION SUMMARY ===\n');
for i = 1:numel(all_handoff_tables)
    T = all_handoff_tables{i};
    cname = char(string(T.cohort(1)));
    fprintf('  %s: %d/%d trials excluded from FRN_peak (%.1f%%) — FRN_mean_amp available for all\n', ...
        cname, sum(T.FRN_excluded), height(T), 100*sum(T.FRN_excluded)/height(T));
end

fprintf('\nAll done.\n');


function T = assign_stages_preserve_LE_RN(T, block_col, trial_col, rev_trials_vec, stage_names)
%ASSIGN_STAGES_PRESERVE_LE_RN
%
% Assign LN, LE, RN, RE to a behavioural trial table.
%
% One row = one behavioural trial.
% Stage is an annotation.
% Trials outside the four stage windows remain missing.
%
% Critical rule:
%   LE and RN are protected.
%   If stage windows overlap, overlap is removed from LN/RE, not LE/RN.
%
% Required columns:
%   block_col: block number column, e.g. 'block'
%   trial_col: within-block trial number column, e.g. 'trial'
%
% rev_trials_vec:
%   Vector indexed by sequential real block position.
%   Example: if table block values are [2 3 4 5], then:
%       block 2 -> rev_trials_vec(1)
%       block 3 -> rev_trials_vec(2)
%       block 4 -> rev_trials_vec(3)
%       block 5 -> rev_trials_vec(4)

STAGE_LEN = 20;

if ~ismember(block_col, T.Properties.VariableNames)
    error('assign_stages_preserve_LE_RN: column "%s" not found.', block_col);
end

if ~ismember(trial_col, T.Properties.VariableNames)
    error('assign_stages_preserve_LE_RN: column "%s" not found.', trial_col);
end

% Initialise output columns
T.stage = categorical(repmat(missing, height(T), 1), stage_names, 'Ordinal', false);
T.in_stage_window = false(height(T), 1);
T.stage_overlap_resolved = false(height(T), 1);

% Convert block and trial columns robustly to numeric
block_nums = local_to_numeric(T.(block_col));
trial_nums = local_to_numeric(T.(trial_col));

if all(isnan(block_nums))
    error('assign_stages_preserve_LE_RN: block column "%s" could not be converted to numeric.', block_col);
end

if all(isnan(trial_nums))
    error('assign_stages_preserve_LE_RN: trial column "%s" could not be converted to numeric.', trial_col);
end

unique_blocks = unique(block_nums(~isnan(block_nums)));
unique_blocks = sort(unique_blocks(:)');

for bi = 1:numel(unique_blocks)

    raw_block = unique_blocks(bi);

    % Full-table row indices for this block
    block_rows = find(block_nums == raw_block);

    if isempty(block_rows)
        continue
    end

    % Trial numbers within this block, same length as block_rows
    tib = trial_nums(block_rows);

    valid_tib = ~isnan(tib);
    if ~any(valid_tib)
        continue
    end

    max_trial = max(tib(valid_tib), [], 'omitnan');

    if isnan(max_trial) || max_trial < 1
        continue
    end

    % Reversal trial for this sequential block
    rev_trial = NaN;
    if ~isempty(rev_trials_vec) && bi <= numel(rev_trials_vec)
        rev_trial = rev_trials_vec(bi);
    end

    % Candidate windows
    raw_ln = 1:min(STAGE_LEN, max_trial);
    raw_re = max(1, max_trial - STAGE_LEN + 1):max_trial;

    if isnan(rev_trial)
        % No reversal info: only assign LN and RE
        ln_trials = raw_ln;
        le_trials = [];
        rn_trials = [];
        re_trials = raw_re;

    else
        % Protected windows
        le_start = max(1, rev_trial - STAGE_LEN);
        le_end   = max(1, rev_trial - 1);

        rn_start = rev_trial;
        rn_end   = min(max_trial, rev_trial + STAGE_LEN - 1);

        le_trials = le_start:le_end;
        rn_trials = rn_start:rn_end;

        protected = unique([le_trials, rn_trials]);

        % LN and RE lose any overlap with LE/RN
        ln_trials = setdiff(raw_ln, protected);
        re_trials = setdiff(raw_re, protected);

        % Diagnostic: rows that would have been LN/RE but were protected
        removed_from_ln = intersect(raw_ln, protected);
        removed_from_re = intersect(raw_re, protected);
        overlap_trials = unique([removed_from_ln, removed_from_re]);

        if ~isempty(overlap_trials)
            overlap_rows_local = ismember(tib, overlap_trials);
            overlap_rows_global = block_rows(overlap_rows_local);
            T.stage_overlap_resolved(overlap_rows_global) = true;
        end
    end

    % Assign stages.
    % This uses block_rows to map within-block logicals back to the full table.
    T = set_stage_for_block(T, block_rows, tib, ln_trials, 'LN', stage_names);
    T = set_stage_for_block(T, block_rows, tib, le_trials, 'LE', stage_names);
    T = set_stage_for_block(T, block_rows, tib, rn_trials, 'RN', stage_names);
    T = set_stage_for_block(T, block_rows, tib, re_trials, 'RE', stage_names);

    % Debug print for this block
    fprintf(['    block raw=%g seq=%d | rev=%g | LN=%d LE=%d RN=%d RE=%d | ' ...
             'missing=%d | overlap_resolved=%d\n'], ...
        raw_block, bi, rev_trial, ...
        sum(T.stage(block_rows) == 'LN'), ...
        sum(T.stage(block_rows) == 'LE'), ...
        sum(T.stage(block_rows) == 'RN'), ...
        sum(T.stage(block_rows) == 'RE'), ...
        sum(ismissing(T.stage(block_rows))), ...
        sum(T.stage_overlap_resolved(block_rows)));
end

end


function T = set_stage_for_block(T, block_rows, tib, trial_list, label, stage_names)

if isempty(trial_list)
    return
end

% Local mask within this block
local_mask = ismember(tib, trial_list);

% Convert local block rows back to full-table rows
global_rows = block_rows(local_mask);

if isempty(global_rows)
    return
end

T.stage(global_rows) = categorical({label}, stage_names);
T.in_stage_window(global_rows) = true;

end


function x = local_to_numeric(v)

if isnumeric(v)
    x = double(v);

elseif islogical(v)
    x = double(v);

elseif iscategorical(v)
    x = str2double(string(v));

elseif isstring(v)
    x = str2double(v);

elseif iscell(v)
    x = str2double(string(v));

else
    error('Unsupported column class: %s', class(v));
end

x = x(:);

end