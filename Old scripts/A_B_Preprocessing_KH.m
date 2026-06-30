% =============================================================================
% MERGED_KH_PREPROCESS_AND_OUTCOME_FEATURES_v4_debugged.m
%
% PURPOSE
% -------
% One-script KH pipeline:
%   1. Load raw KH EEG.
%   2. Preprocess using v4-style noise harmonisation / burst interpolation.
%   3. Automatically reject ICLabel components with Brain probability < 0.50.
%   4. Export OUTCOME broadband, theta-amplitude, and theta-phase epochs.
%   5. Trim practice/misaligned trials while preserving trial alignment.
%   6. Save trial2epoch + hw_filter_label + epoch_artifact_flag_trimmed.
%   7. Build subsequent OUTCOME ERP/theta/PLV feature tables.
%
% TO RUN ALL VALID KH SUBJECTS
% ----------------------------
%   Set RUN_ALL_VALID_KH = true.
%
% IMPORTANT DESIGN CHOICE
% -----------------------
% This merged script exports only outcome epochs, because the subsequent
% feature-table script only needs outcome broadband/theta/phase epochs.
% Stimulus/confidence/response epoch exports from the original preprocessing
% script are intentionally omitted here to keep the merged pipeline coherent.
%
% ICLabel CHANGE REQUESTED
% ------------------------
% Manual review is commented out/removed. Components are rejected if:
%   Brain probability < 0.50
%
% =============================================================================

clear; close all; clc;

% -------------------------------------------------------------------------
%% USER SWITCHES
% -------------------------------------------------------------------------

RUN_PREPROCESSING       = true;
RUN_FEATURE_EXTRACTION  = true;

RUN_ALL_VALID_KH = true;

KH_VALID_FULL =  [3:12, 14:23, 27,28];
KH_PILOT_ONLY = [12, 14];

if RUN_ALL_VALID_KH
    valid_participants = KH_VALID_FULL;
else
    valid_participants = KH_PILOT_ONLY;
end

save_tables  = true;
save_figures = true;

% -------------------------------------------------------------------------
%% PATHS
% -------------------------------------------------------------------------

remote = 0;

if remote == 1
    base_path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch';
else
    base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch';
end

eeglab_path = 'C:\Users\khatinova\OneDrive - Nexus365\Pre_2026_Folders\Documents\MATLAB\eeglab2025.1.0';

KH_data_path      = fullfile(base_path, 'Salient mod switch KH', 'Data');
study_filepath    = fullfile(base_path, 'Salient mod switch KH', 'Results', ...
    'EEG analysis', 'Winter 2026');
epoch_outpath     = fullfile(base_path, 'Salient mod switch KH', 'Results', ...
    'EEG analysis', 'Epoched_data_noisefiltering');
feature_outpath   = fullfile(base_path, 'Salient mod switch KH', 'Results', ...
    'EEG analysis', 'Outcome_feature_tables_v4_merged');
figure_outpath    = fullfile(base_path, 'Salient mod switch KH', 'Results', ...
    'EEG analysis', 'Figures', 'outcome_v4_merged_QC');

if ~exist(study_filepath, 'dir'),  mkdir(study_filepath);  end
if ~exist(epoch_outpath, 'dir'),   mkdir(epoch_outpath);   end
if ~exist(feature_outpath, 'dir'), mkdir(feature_outpath); end
if ~exist(figure_outpath, 'dir'),  mkdir(figure_outpath);  end

addpath(eeglab_path);
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab('nogui'); %#ok<ASGLU>

load(fullfile(KH_data_path, 'all_trial_data_June2026.mat'));

KH_behav_file = fullfile(KH_data_path, 'behav_table_June2026.mat');
if exist(KH_behav_file, 'file')
    S_beh = load(KH_behav_file, 'group_T');
    KH_behav_table = S_beh.group_T;
else
    error('behav_table.mat not found at %s', KH_behav_file);
end

if ismember('researcher', KH_behav_table.Properties.VariableNames)
    KH_behav_table = KH_behav_table(strcmp(string(KH_behav_table.researcher), 'KH'), :);
end

% -------------------------------------------------------------------------
%% PREPROCESSING PARAMETERS
% -------------------------------------------------------------------------

do_theta_epochs = true;
do_trimming     = true;

% Outcome epoch window, in seconds.
outcome_window = [-0.5 2.0];

% Baseline in ms.
baseline_win_ms = [-500 0];

% Cohort-2 outcome trigger correction.
% Ox09 onward in KH had outcome trigger after visual onset.
OUTCOME_DELAY_C2_S = 0.06;
misaligned_trigger_IDs = {'Ox09', 'Ox10', 'Ox11', 'Ox12', 'Ox14', ...
                          'Ox15', 'Ox16', 'Ox17'};

% Hardware-filter registry.
% EDIT THESE LISTS once the acquisition log is known.
subject_hw_filter_30 = [];
subject_hw_filter_70 = [];

HARMONISE_LOCUTOFF_HZ     = 30;
APPLY_HARMONISING_FILTER  = true;

% Theta filter parameters.
THETA_LO  = 4;
THETA_HI  = 8;
THETA_ORD = 3;

% Continuous burst interpolation parameters.
BURST_SD_THRESH       = 6;
BURST_GRAD_SD_THRESH  = 6;
BURST_PAD_MS          = 20;
BURST_MAX_RUN_MS      = 300;

% Per-epoch artifact flagging parameters.
EPOCH_PTP_SD_THRESH  = 5;
EPOCH_GRAD_SD_THRESH = 5;

% Subjects with LM/RM mastoids.
LMRM_SUBJECT_IDS = [17:23, 27,28];

protect_labels = {'VEOG','HEOG','EOG','TRIGGER','LM','RM'};

% ICLabel auto-rejection threshold.
ICLABEL_MIN_BRAIN_PROB = 0.50;

% -------------------------------------------------------------------------
%% FEATURE EXTRACTION PARAMETERS
% -------------------------------------------------------------------------

ERP_plot_window = [-200 1000];
rm_baseline     = [-200 0];

N2_win    = [120 350];
FRN_win   = [250 350];
RewP_win  = [250 350];
P300_win  = [300 600];
Theta_win = [200 500];
PLV_win   = [200 400];
PLV_baseline = [-200 0];

MIN_TRIALS_PLV_WINDOW = 5;
PLV_WINDOW_HALF       = 7;

stage_names  = {'LN','LE','RN','RE'};
BTYPE_LABELS = {'D','P'};

fcz_label    = 'FCz';
cz_label     = 'Cz';
par_channels = {'Pz','P1','P2'};
acc_channels = {'FCz','Fz','AFz','F1','F2'};
som_channels = {'C3','C4','CP3','CP1','C5','CP5'};

block_col = 'block';
trial_col = 'trial';

% =============================================================================
%% PART 1: PREPROCESS AND EXPORT OUTCOME EPOCHS
% =============================================================================

long_burst_log = table();

if RUN_PREPROCESSING

    data_dirs = dir(fullfile(KH_data_path, 'Ox*'));
    data_dirs = data_dirs([data_dirs.isdir]);
    all_subj_names = {data_dirs.name};

    pop_editoptions('option_storedisk', 0);

    for i = valid_participants

        subjID = sprintf('Ox%02d', i);
        is_cohort1 = i <= 8;
        has_LMRM = ismember(i, LMRM_SUBJECT_IDS);

        fprintf('\n======================================================\n');
        fprintf('PREPROCESSING %s\n', subjID);
        fprintf('======================================================\n');

        % -----------------------------------------------------------------
        % Hardware-filter bookkeeping.
        % -----------------------------------------------------------------
        if ismember(i, subject_hw_filter_30)
            hw_filter_label = '0_30Hz';
        elseif ismember(i, subject_hw_filter_70)
            hw_filter_label = '0_70Hz';
        else
            hw_filter_label = 'unknown';
            fprintf(['  [hardware filter] %s not listed in subject_hw_filter_30 ' ...
                     'or subject_hw_filter_70. Conservative harmonisation will run.\n'], subjID);
        end

        % -----------------------------------------------------------------
        % Locate raw subject directory.
        % -----------------------------------------------------------------
        subj_dir_idx = find(strcmp(all_subj_names, subjID), 1);
        if isempty(subj_dir_idx)
            warning('%s: directory not found, skipping.', subjID);
            continue;
        end

        subjPath = fullfile(KH_data_path, data_dirs(subj_dir_idx).name);

        curryFile = dir(fullfile(subjPath, 'Acquisition*.dat'));
        if isempty(curryFile)
            warning('%s: no Curry file found, skipping.', subjID);
            continue;
        end

        % -----------------------------------------------------------------
        % 1. Load raw Curry data.
        % -----------------------------------------------------------------
        EEG = loadcurry(fullfile(subjPath, curryFile(1).name), ...
            'KeepTriggerChannel', 'True', ...
            'CurryLocations', 'False');

        EEG.subject = subjID;
        EEG.etc.hw_filter_label = hw_filter_label;

        fprintf('  Loaded raw: %d channels, %.1f seconds\n', ...
            EEG.nbchan, EEG.pnts / EEG.srate);

        % -----------------------------------------------------------------
        % 2. Trim trailing data and resample.
        % -----------------------------------------------------------------
        if ~isempty(EEG.event)
            last_latency = max([EEG.event.latency]);
            trim_sample = min(last_latency + 10 * EEG.srate, EEG.pnts);
            EEG = pop_select(EEG, 'point', [1 trim_sample]);
        end

        EEG = pop_resample(EEG, 500);
        EEG = eeg_checkset(EEG);

        % -----------------------------------------------------------------
        % 3a. Harmonising low-pass.
        % -----------------------------------------------------------------
        if APPLY_HARMONISING_FILTER
            fprintf('  Harmonising low-pass at %g Hz, hw_filter=%s\n', ...
                HARMONISE_LOCUTOFF_HZ, hw_filter_label);
            EEG = pop_eegfiltnew(EEG, 'hicutoff', HARMONISE_LOCUTOFF_HZ);
            EEG = eeg_checkset(EEG);
        end

        % -----------------------------------------------------------------
        % 3b. Standard analysis filters.
        % -----------------------------------------------------------------
        EEG = pop_eegfiltnew(EEG, 'locutoff', 0.5);
        EEG = pop_eegfiltnew(EEG, 'hicutoff', 40);
        EEG = eeg_checkset(EEG);

        % -----------------------------------------------------------------
        % 4. Protect EOG/trigger/mastoid channels.
        % -----------------------------------------------------------------
        orig_chanlocs = EEG.chanlocs;
        orig_labels   = {orig_chanlocs.labels};
        all_labels    = {EEG.chanlocs.labels};

        protect_idx = find(ismember(lower(all_labels), lower(protect_labels)));

        if ~isempty(protect_idx)
            EEG_protect = pop_select(EEG, 'channel', protect_idx);
            EEG_scalp   = pop_select(EEG, 'nochannel', protect_idx);
        else
            EEG_protect = [];
            EEG_scalp   = EEG;
        end
        EEG_scalp = eeg_checkset(EEG_scalp);

        % -----------------------------------------------------------------
        % 5. LM/RM motion regression.
        % -----------------------------------------------------------------
        if has_LMRM
            lmrm_labels = {'LM','RM'};
            lmrm_idx_full = find(ismember(lower(orig_labels), lower(lmrm_labels)));

            if numel(lmrm_idx_full) >= 1
                fprintf('  LM/RM motion regression...\n');

                motion_data = EEG.data(lmrm_idx_full, :);
                motion_ref  = mean(double(motion_data), 1, 'omitnan');
                motion_ref  = motion_ref - mean(motion_ref, 'omitnan');

                scalp_in_full = setdiff(1:EEG.nbchan, protect_idx);

                denom = motion_ref * motion_ref';
                if denom > 0
                    for ch = scalp_in_full
                        sig = double(EEG.data(ch, :));
                        sig0 = sig - mean(sig, 'omitnan');
                        beta = sig0 * motion_ref' / denom;
                        EEG.data(ch, :) = sig - beta * motion_ref;
                    end
                else
                    warning('%s: LM/RM motion reference variance was zero; skipping regression.', subjID);
                end

                if ~isempty(protect_idx)
                    EEG_scalp = pop_select(EEG, 'nochannel', protect_idx);
                else
                    EEG_scalp = EEG;
                end
                EEG_scalp = eeg_checkset(EEG_scalp);
            else
                warning('%s listed as LMRM subject, but LM/RM channels were not found.', subjID);
            end
        end

        % -----------------------------------------------------------------
        % 5b. Alignment-preserving continuous burst interpolation.
        % -----------------------------------------------------------------
        [EEG_scalp, burst_report] = interpolate_amplitude_bursts( ...
            EEG_scalp, ...
            BURST_SD_THRESH, ...
            BURST_GRAD_SD_THRESH, ...
            BURST_PAD_MS, ...
            BURST_MAX_RUN_MS);

        fprintf('  Burst interpolation: %d short bursts repaired, %d long bursts flagged.\n', ...
            burst_report.n_short_bursts, burst_report.n_long_bursts);

        if burst_report.n_long_bursts > 0
            new_rows = table( ...
                repmat({subjID}, burst_report.n_long_bursts, 1), ...
                burst_report.long_run_start_s, ...
                burst_report.long_run_end_s, ...
                burst_report.long_run_channels, ...
                'VariableNames', {'subjID','start_s','end_s','channels'});
            long_burst_log = [long_burst_log; new_rows]; %#ok<AGROW>
        end

        % -----------------------------------------------------------------
        % 6. ASR channel-level cleaning only; no segment rejection.
        % -----------------------------------------------------------------
        EEG_ASR = pop_clean_rawdata(EEG_scalp, ...
            'FlatlineCriterion',  5, ...
            'ChannelCriterion',   0.8, ...
            'LineNoiseCriterion', 4, ...
            'BurstCriterion',     'off', ...
            'WindowCriterion',    'off', ...
            'BurstRejection',     'off');
        EEG_ASR = eeg_checkset(EEG_ASR);

        scalp_before = {EEG_scalp.chanlocs.labels};
        scalp_after  = {EEG_ASR.chanlocs.labels};
        removed = setdiff(scalp_before, scalp_after, 'stable');

        if ~isempty(removed)
            fprintf('  Interpolating ASR-removed channels: %s\n', strjoin(removed, ', '));
            EEG_ASR = pop_interp(EEG_ASR, EEG_scalp.chanlocs, 'spherical');
            EEG_ASR = eeg_checkset(EEG_ASR);
        else
            fprintf('  No channels removed by ASR.\n');
        end

        % -----------------------------------------------------------------
        % 7. Reinsert protected channels.
        % -----------------------------------------------------------------
        if ~isempty(protect_idx)
            nFullChan = numel(orig_chanlocs);
            nSamp     = size(EEG_ASR.data, 2);
            newdata   = zeros(nFullChan, nSamp);

            scalp_labels = {EEG_ASR.chanlocs.labels};

            for c = 1:numel(scalp_labels)
                orig_idx = find(strcmp(orig_labels, scalp_labels{c}), 1);
                if ~isempty(orig_idx)
                    newdata(orig_idx, :) = EEG_ASR.data(c, :);
                end
            end

            for pp = 1:numel(protect_idx)
                newdata(protect_idx(pp), :) = EEG_protect.data(pp, :);
            end

            EEG.data     = newdata;
            EEG.nbchan   = nFullChan;
            EEG.chanlocs = orig_chanlocs;
        else
            EEG = EEG_ASR;
        end
        EEG = eeg_checkset(EEG);

        % -----------------------------------------------------------------
        % 8. Average reference, excluding protected channels.
        % -----------------------------------------------------------------
        if ~isempty(protect_idx)
            EEG = pop_reref(EEG, [], 'exclude', protect_idx);
        else
            EEG = pop_reref(EEG, []);
        end
        EEG = eeg_checkset(EEG);

        pop_saveset(EEG, 'filename', [subjID '_cleaned_v4_merged.set'], ...
            'filepath', study_filepath);

        % -----------------------------------------------------------------
        % 9. ICA on 1 Hz copy, transfer weights to 0.5 Hz data.
        % -----------------------------------------------------------------
        EEG_ica = pop_eegfiltnew(EEG, 'locutoff', 1.0);
        EEG_ica = pop_runica(EEG_ica, 'extended', 1);

        EEG.icaweights  = EEG_ica.icaweights;
        EEG.icasphere   = EEG_ica.icasphere;
        EEG.icawinv     = EEG_ica.icawinv;
        EEG.icachansind = EEG_ica.icachansind;
        EEG = eeg_checkset(EEG);

        EEG = pop_iclabel(EEG, 'default');

        % -----------------------------------------------------------------
        % ICLabel auto-rejection requested:
        % Commented-out/manual route removed. Reject all components with
        % Brain probability < 0.50.
        % -----------------------------------------------------------------

        % OLD MANUAL ROUTE:
        EEG = pop_icflag(EEG, ...
            [0   0.2;   % Brain
             0.5 1;     % Muscle
             0.5 1;     % Eye
             0.5 1;     % Heart
             0.5 1;     % Line noise
             0.5 1;     % Channel noise
             0.5 1]);   % Other
        % EEG = review_iclabel_interactive(EEG, 'nCols', 10);
        EEG = pop_subcomp(EEG, []);

        if ~isfield(EEG.etc, 'ic_classification') || ...
           ~isfield(EEG.etc.ic_classification, 'ICLabel') || ...
           ~isfield(EEG.etc.ic_classification.ICLabel, 'classifications')
            error('%s: ICLabel classifications not found after pop_iclabel.', subjID);
        end

        ic_probs = EEG.etc.ic_classification.ICLabel.classifications;
        brain_prob = ic_probs(:, 1);

        reject_comps = find(brain_prob < ICLABEL_MIN_BRAIN_PROB);

        fprintf('  ICLabel auto-rejection: rejecting %d/%d ICs with Brain probability < %.2f.\n', ...
            numel(reject_comps), numel(brain_prob), ICLABEL_MIN_BRAIN_PROB);

        EEG.etc.iclabel_auto_reject.brain_prob = 0.5;
        EEG.etc.iclabel_auto_reject.threshold = ICLABEL_MIN_BRAIN_PROB;
        EEG.etc.iclabel_auto_reject.rejected_components = reject_comps(:);

        if ~isempty(reject_comps)
            EEG = pop_subcomp(EEG, reject_comps, 0);
        end
        EEG = eeg_checkset(EEG);

        pop_saveset(EEG, 'filename', [subjID '_ICA_pruned_June26.set'], ...
            'filepath', study_filepath, 'savemode', 'onefile');

        % -----------------------------------------------------------------
        % 10. Theta-filter continuous data and Hilbert phase.
        % -----------------------------------------------------------------
        if do_theta_epochs
            fprintf('  Theta-filtering continuous data %g-%g Hz...\n', THETA_LO, THETA_HI);

            fs_now = EEG.srate;
            [bf, af] = butter(THETA_ORD, [THETA_LO THETA_HI] / (fs_now / 2), 'bandpass');

            all_labels = {EEG.chanlocs.labels};
            scalp_all = find(~ismember(lower(all_labels), lower(protect_labels)));

            EEG_theta = EEG;
            for ch = scalp_all
                EEG_theta.data(ch, :) = filtfilt(bf, af, double(EEG.data(ch, :)));
            end
            EEG_theta = eeg_checkset(EEG_theta);

            EEG_phase = EEG_theta;
            for ch = scalp_all
                analytic = hilbert(double(EEG_theta.data(ch, :)));
                EEG_phase.data(ch, :) = angle(analytic);
            end
            EEG_phase = eeg_checkset(EEG_phase);
        else
            EEG_theta = [];
            EEG_phase = [];
        end

        % -----------------------------------------------------------------
        % 11. Outcome event codes and timing correction.
        % -----------------------------------------------------------------
        if is_cohort1
            outcome_codes = {'10','11', 10, 11};
        else
            outcome_codes = {'31','32','33','34', 31, 32, 33, 34};
        end

        should_shift_outcome = ismember(subjID, misaligned_trigger_IDs);

        % Epoch window relative to the recorded trigger.
        % For cohort-2 delayed trigger, include pre-trigger offset so that
        % after relabelling, t = 0 means visual outcome onset.
        if should_shift_outcome
            outcome_window_epoch = outcome_window - OUTCOME_DELAY_C2_S;
        else
            outcome_window_epoch = outcome_window;
        end

        % -----------------------------------------------------------------
        % 12. Epoch outcome broadband/theta/phase.
        % -----------------------------------------------------------------
        fprintf('  Epoching outcome broadband...\n');
        EEG_out = pop_epoch(EEG, outcome_codes, outcome_window_epoch, ...
            'epochinfo', 'yes');
        EEG_out = pop_rmbase(EEG_out, baseline_win_ms);
        EEG_out = eeg_checkset(EEG_out);

        if should_shift_outcome
            EEG_out.times = EEG_out.times + OUTCOME_DELAY_C2_S * EEGp.srate;
            EEG_out.xmin  = EEG_out.xmin  + OUTCOME_DELAY_C2_S;
            EEG_out.xmax  = EEG_out.xmax  + OUTCOME_DELAY_C2_S;
        end

        pop_saveset(EEG_out, 'filename', [subjID '_outcome.set'], ...
            'filepath', epoch_outpath, 'savemode', 'onefile');

        if do_theta_epochs
            fprintf('  Epoching outcome theta amplitude source...\n');
            EEG_theta_out = pop_epoch(EEG_theta, outcome_codes, outcome_window_epoch, ...
                'epochinfo', 'yes');
            EEG_theta_out = pop_rmbase(EEG_theta_out, baseline_win_ms);
            EEG_theta_out = eeg_checkset(EEG_theta_out);

            fprintf('  Epoching outcome theta phase source...\n');
            EEG_phase_out = pop_epoch(EEG_phase, outcome_codes, outcome_window_epoch, ...
                'epochinfo', 'yes');
            EEG_phase_out = eeg_checkset(EEG_phase_out);

            if should_shift_outcome
                EEG_theta_out.times = EEG_theta_out.times + OUTCOME_DELAY_C2_S * EEGp.srate;
                EEG_theta_out.xmin  = EEG_theta_out.xmin  + OUTCOME_DELAY_C2_S;
                EEG_theta_out.xmax  = EEG_theta_out.xmax  + OUTCOME_DELAY_C2_S;

                EEG_phase_out.times = EEG_phase_out.times + OUTCOME_DELAY_C2_S * EEGp.srate;
                EEG_phase_out.xmin  = EEG_phase_out.xmin  + OUTCOME_DELAY_C2_S;
                EEG_phase_out.xmax  = EEG_phase_out.xmax  + OUTCOME_DELAY_C2_S;
            end

            pop_saveset(EEG_theta_out, 'filename', [subjID '_outcome_theta.set'], ...
                'filepath', epoch_outpath, 'savemode', 'onefile');

            pop_saveset(EEG_phase_out, 'filename', [subjID '_outcome_phase.set'], ...
                'filepath', epoch_outpath, 'savemode', 'onefile');
        end

        % -----------------------------------------------------------------
        % 13. Build trial2epoch and trim practice/misaligned epochs.
        % -----------------------------------------------------------------
        if do_trimming
            if ~isfield(all_trial_data, subjID)
                warning('%s: missing all_trial_data; cannot build behavioural alignment. Saving untrimmed only.', subjID);
                continue;
            end

            beh = all_trial_data.(subjID).trial_data;
            beh_correct = [];

            for b = 1:height(beh.correct)
                beh_correct = [beh_correct, beh.correct(b, :)]; %#ok<AGROW>
            end

            beh_correct = beh_correct(:);
            beh_correct = beh_correct(~isnan(beh_correct));

            [trial2epoch_original, diag_out] = KH_align_epochs_with_offset(EEG_out, beh_correct);

            fprintf('  trial2epoch: %d/%d matched, offset=%d, match rate %.1f%%\n', ...
                diag_out.n_matched, diag_out.n_trials, diag_out.best_offset, ...
                100 * diag_out.match_rate);

            % valid_ep are epoch indices in the untrimmed outcome file.
            valid_ep = unique(trial2epoch_original(~isnan(trial2epoch_original)));
            valid_ep = valid_ep(:);
            valid_ep = valid_ep(valid_ep >= 1 & valid_ep <= EEG_out.trials);

            if isempty(valid_ep)
                warning('%s: no valid epochs after alignment; not saving trimmed files.', subjID);
                continue;
            end

            EEG_out_trim = pop_select(EEG_out, 'trial', valid_ep);
            EEG_out_trim = eeg_checkset(EEG_out_trim);

            pop_saveset(EEG_out_trim, 'filename', [subjID '_outcome_trimmed.set'], ...
                'filepath', epoch_outpath, 'savemode', 'onefile');

            if do_theta_epochs
                valid_ep_theta = valid_ep(valid_ep <= EEG_theta_out.trials);
                EEG_theta_trim = pop_select(EEG_theta_out, 'trial', valid_ep_theta);
                EEG_theta_trim = eeg_checkset(EEG_theta_trim);

                EEG_phase_trim = pop_select(EEG_phase_out, 'trial', valid_ep_theta);
                EEG_phase_trim = eeg_checkset(EEG_phase_trim);

                pop_saveset(EEG_theta_trim, 'filename', [subjID '_outcome_theta_trimmed.set'], ...
                    'filepath', epoch_outpath, 'savemode', 'onefile');

                pop_saveset(EEG_phase_trim, 'filename', [subjID '_outcome_phase_trimmed.set'], ...
                    'filepath', epoch_outpath, 'savemode', 'onefile');
            end

            % Per-epoch artifact flagging on the trimmed broadband outcome file.
            epoch_artifact_flag_trimmed = flag_outcome_epochs( ...
                EEG_out_trim, ...
                EPOCH_PTP_SD_THRESH, ...
                EPOCH_GRAD_SD_THRESH, ...
                protect_labels);

            % Critical: trial2epoch_original maps behaviour rows to original
            % untrimmed epoch indices. Downstream code remaps using valid_ep.
            trial2epoch = trial2epoch_original(:);
            trial2epoch_out = trial2epoch; %#ok<NASGU>

            save(fullfile(epoch_outpath, [subjID '_trial2epoch.mat']), ...
                'trial2epoch', ...
                'trial2epoch_out', ...
                'valid_ep', ...
                'hw_filter_label', ...
                'epoch_artifact_flag_trimmed', ...
                'diag_out');

            fprintf('  Saved trimmed outcome files and trial2epoch cache.\n');
        end

        clear EEG EEG_ica EEG_theta EEG_phase EEG_out EEG_theta_out EEG_phase_out
    end

    if ~isempty(long_burst_log)
        writetable(long_burst_log, fullfile(study_filepath, 'long_burst_log_v4_merged.csv'));
        save(fullfile(study_filepath, 'long_burst_log_v4_merged.mat'), 'long_burst_log');
    end
end

% =============================================================================
%% PART 2: OUTCOME FEATURE TABLE CONSTRUCTION
% =============================================================================

all_trials_table = table();
stage_feature_table = table();

if RUN_FEATURE_EXTRACTION

    fprintf('\n======================================================\n');
    fprintf('BUILDING OUTCOME FEATURE TABLES\n');
    fprintf('======================================================\n');

    for participant = valid_participants

        subj = sprintf('Ox%02d', participant);
        fprintf('\n============ FEATURE EXTRACTION %s ============\n', subj);

        if ~isfield(all_trial_data, subj)
            warning('%s missing from all_trial_data. Skipping.', subj);
            continue;
        end

        % -----------------------------------------------------------------
        % Load broadband/theta/phase outcome epochs.
        % -----------------------------------------------------------------
        EEGp = load_first_existing_set(epoch_outpath, { ...
            sprintf('%s_outcome_trimmed.set', subj), ...
            sprintf('%s_outcome.set', subj)});

        if isempty(EEGp)
            warning('%s: no broadband outcome file found. Skipping.', subj);
            continue;
        end

        EEGp_theta = load_first_existing_set(epoch_outpath, { ...
            sprintf('%s_outcome_theta_trimmed.set', subj), ...
            sprintf('%s_outcome_theta.set', subj)});

        EEGp_phase = load_first_existing_set(epoch_outpath, { ...
            sprintf('%s_outcome_phase_trimmed.set', subj), ...
            sprintf('%s_outcome_phase.set', subj)});

        % -----------------------------------------------------------------
        % Behavioural spine.
        % -----------------------------------------------------------------
        subj_rows = string(KH_behav_table.subjID) == string(subj);
        if ~any(subj_rows)
            warning('%s has no rows in KH_behav_table. Skipping.', subj);
            continue;
        end

        subj_features = KH_behav_table(subj_rows, :);
        n_rows = height(subj_features);

        subj_features.subj_id = repmat(string(subj), n_rows, 1);
        subj_features.cohort  = repmat("KH", n_rows, 1);
        subj_features.subj    = repmat(participant, n_rows, 1);
        subj_features.is_cohort1 = repmat(participant <= 8, n_rows, 1);

        if ~ismember('trial_continuous', subj_features.Properties.VariableNames)
            subj_features.trial_continuous = (1:n_rows)';
        end

        if ~ismember('block_type', subj_features.Properties.VariableNames)
            subj_features.block_type = categorical(repmat(missing, n_rows, 1), BTYPE_LABELS);
        else
            subj_features.block_type = categorical(string(subj_features.block_type), BTYPE_LABELS);
        end

        if ~ismember('false_fb', subj_features.Properties.VariableNames)
            subj_features.false_fb = false(n_rows, 1);
        else
            subj_features.false_fb = logical(local_to_numeric(subj_features.false_fb));
        end

        if ~ismember('fb_shown_correct', subj_features.Properties.VariableNames)
            subj_features.fb_shown_correct = nan(n_rows, 1);
        end

        if ~ismember('confidence', subj_features.Properties.VariableNames)
            subj_features.confidence = nan(n_rows, 1);
        end

        if ~ismember('correct', subj_features.Properties.VariableNames)
            error('%s: subj_features lacks correct column.', subj);
        end

        subj_features.correct_num = local_correct_to_numeric(subj_features.correct);

        % -----------------------------------------------------------------
        % Load and remap trial2epoch.
        % -----------------------------------------------------------------
        trial2epoch_file = fullfile(epoch_outpath, sprintf('%s_trial2epoch.mat', subj));

        hw_filter_label_subj = "unknown";
        epoch_artifact_flag_eeg = false(EEGp.trials, 1);

        if exist(trial2epoch_file, 'file')
            S2 = load(trial2epoch_file);

            if isfield(S2, 'trial2epoch_out')
                trial2epoch_raw = S2.trial2epoch_out(:);
            elseif isfield(S2, 'trial2epoch')
                trial2epoch_raw = S2.trial2epoch(:);
            else
                error('%s lacks trial2epoch/trial2epoch_out.', trial2epoch_file);
            end

            trial2epoch_raw = pad_or_truncate(trial2epoch_raw, n_rows);

            if isfield(S2, 'valid_ep')
                valid_ep = S2.valid_ep(:);
            else
                valid_ep = [];
            end

            using_trimmed_file = EEGp.trials <= numel(valid_ep) || ...
                exist(fullfile(epoch_outpath, sprintf('%s_outcome_trimmed.set', subj)), 'file');

            trial2epoch = remap_trial2epoch_to_loaded_file( ...
                trial2epoch_raw, valid_ep, EEGp.trials, using_trimmed_file);

            if isfield(S2, 'hw_filter_label') && ~isempty(S2.hw_filter_label)
                hw_filter_label_subj = normalize_hw_filter_label(S2.hw_filter_label);
            end

            epoch_artifact_flag_eeg = read_epoch_artifact_flags(S2, valid_ep, EEGp.trials);

        else
            warning('%s: trial2epoch cache missing. Falling back to sequential mapping.', subj);
            trial2epoch = nan(n_rows, 1);
            n_map = min(n_rows, EEGp.trials);
            trial2epoch(1:n_map) = (1:n_map)';
        end

        subj_features.epoch = nan(n_rows, 1);
        n_map = min(numel(trial2epoch), n_rows);
        subj_features.epoch(1:n_map) = trial2epoch(1:n_map);

        subj_features.has_eeg_epoch = ~isnan(subj_features.epoch) & ...
            subj_features.epoch >= 1 & subj_features.epoch <= EEGp.trials;

        subj_features.hw_filter_label = repmat(categorical( ...
            hw_filter_label_subj, {'0_30Hz','0_70Hz','unknown'}), n_rows, 1);

        subj_features.epoch_artifact_flag = false(n_rows, 1);
        valid_ep_mask = subj_features.has_eeg_epoch;
        ep_indices = round(subj_features.epoch(valid_ep_mask));
        ep_indices = min(ep_indices, numel(epoch_artifact_flag_eeg));
        subj_features.epoch_artifact_flag(valid_ep_mask) = epoch_artifact_flag_eeg(ep_indices);

        % -----------------------------------------------------------------
        % Stage assignment.
        % -----------------------------------------------------------------
        beh = all_trial_data.(subj).trial_data;
        rev_trials_vec = [];
        if isfield(beh, 'revTrial') && ~isempty(beh.revTrial)
            rev_trials_vec = beh.revTrial(:);
        end

        subj_features = assign_stages_preserve_LE_RN( ...
            subj_features, block_col, trial_col, rev_trials_vec, stage_names);

        % -----------------------------------------------------------------
        % Channel indices and time windows.
        % -----------------------------------------------------------------
        chan_labels_lower = lower(string({EEGp.chanlocs.labels}));

        fcz_idx = find(chan_labels_lower == lower(string(fcz_label)));
        cz_idx  = find(chan_labels_lower == lower(string(cz_label)));
        par_idx = find(ismember(chan_labels_lower, lower(string(par_channels))));
        acc_idx = find(ismember(chan_labels_lower, lower(string(acc_channels))));
        som_idx = find(ismember(chan_labels_lower, lower(string(som_channels))));

        if isempty(fcz_idx)
            warning('%s: FCz missing. Skipping feature extraction.', subj);
            continue;
        end

        if isempty(cz_idx)
            warning('%s: Cz missing. FCz/Cz features will use FCz only.', subj);
            cz_idx = fcz_idx;
        end

        bl_mask   = EEGp.times >= rm_baseline(1)  & EEGp.times <= rm_baseline(2);
        n2_mask   = EEGp.times >= N2_win(1)       & EEGp.times <= N2_win(2);
        frn_mask  = EEGp.times >= FRN_win(1)      & EEGp.times <= FRN_win(2);
        rewp_mask = EEGp.times >= RewP_win(1)     & EEGp.times <= RewP_win(2);
        p300_mask = EEGp.times >= P300_win(1)     & EEGp.times <= P300_win(2);
        th_mask   = EEGp.times >= Theta_win(1)    & EEGp.times <= Theta_win(2);
        plv_mask  = EEGp.times >= PLV_win(1)      & EEGp.times <= PLV_win(2);
        plv_bl    = EEGp.times >= PLV_baseline(1) & EEGp.times <= PLV_baseline(2);

        bl_data = squeeze(double(EEGp.data(fcz_idx, bl_mask, :)));
        bline_rms = rms(bl_data(:), 'omitnan');

        % -----------------------------------------------------------------
        % Initialise feature columns.
        % -----------------------------------------------------------------
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

        subj_features.PLV_fp          = nan(n_rows, 1);
        subj_features.PLV_fs          = nan(n_rows, 1);
        subj_features.PLV_fp_pairwise = nan(n_rows, 1);
        subj_features.PLV_fs_pairwise = nan(n_rows, 1);

        subj_features.FCzCz_waveform = repmat({[]}, n_rows, 1);
        subj_features.P300_waveform  = repmat({[]}, n_rows, 1);
        subj_features.Theta_waveform = repmat({[]}, n_rows, 1);

        % -----------------------------------------------------------------
        % Per-trial features.
        % -----------------------------------------------------------------
        for ti = 1:n_rows

            ep = subj_features.epoch(ti);
            if isnan(ep) || ep < 1 || ep > EEGp.trials
                continue;
            end
            ep = round(ep);

            % FCz/Cz N2, FRN, RewP.
            if ~isempty(fcz_idx) && ~isempty(cz_idx) && ~isnan(bline_rms)
                sig = mean(double(EEGp.data([fcz_idx(:)' cz_idx(:)'], :, ep)), 1, 'omitnan');
                sig = sig - mean(sig(bl_mask), 'omitnan');
                subj_features.FCzCz_waveform{ti} = sig;

                win_vals = sig(n2_mask);
                win_t = EEGp.times(n2_mask);
                if any(~isnan(win_vals))
                    [pk, ix] = min(win_vals, [], 'omitnan');
                    subj_features.N2_amp(ti) = pk;
                    subj_features.N2_lat(ti) = win_t(ix);
                    if bline_rms > 0
                        subj_features.N2_norm(ti) = pk / bline_rms;
                    end
                end

                frn_vals = sig(frn_mask);
                frn_t = EEGp.times(frn_mask);
                if any(~isnan(frn_vals))
                    subj_features.FRN_mean_amp(ti) = mean(frn_vals, 'omitnan');
                    if bline_rms > 0
                        subj_features.FRN_mean_norm(ti) = ...
                            subj_features.FRN_mean_amp(ti) / bline_rms;
                    end

                    is_min = islocalmin(frn_vals);
                    if any(is_min)
                        cand_vals = frn_vals(is_min);
                        cand_t = frn_t(is_min);
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

                rewp_vals = sig(rewp_mask);
                rewp_t = EEGp.times(rewp_mask);
                if any(~isnan(rewp_vals))
                    subj_features.RewP_mean_amp(ti) = mean(rewp_vals, 'omitnan');
                    if bline_rms > 0
                        subj_features.RewP_mean_norm(ti) = ...
                            subj_features.RewP_mean_amp(ti) / bline_rms;
                    end

                    is_max = islocalmax(rewp_vals);
                    if any(is_max)
                        cand_vals = rewp_vals(is_max);
                        cand_t = rewp_t(is_max);
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

            % P300.
            if ~isempty(par_idx) && ~isnan(bline_rms)
                sig_p = mean(double(EEGp.data(par_idx, :, ep)), 1, 'omitnan');
                sig_p = sig_p - mean(sig_p(bl_mask), 'omitnan');
                subj_features.P300_waveform{ti} = sig_p;

                win_vals = sig_p(p300_mask);
                win_t = EEGp.times(p300_mask);
                if any(~isnan(win_vals))
                    [pk, ix] = max(win_vals, [], 'omitnan');
                    subj_features.P300_amp(ti) = pk;
                    subj_features.P300_peak_lat(ti) = win_t(ix);
                    if bline_rms > 0
                        subj_features.P300_norm(ti) = pk / bline_rms;
                    end
                end
            end

            % Theta amplitude envelope.
            if ~isempty(EEGp_theta) && ~isempty(acc_idx) && ep <= EEGp_theta.trials
                sig_th = mean(double(EEGp_theta.data(acc_idx, :, ep)), 1, 'omitnan');
                env = abs(hilbert(sig_th));
                env = env - mean(env(bl_mask), 'omitnan');
                subj_features.Theta_amp(ti) = mean(env(th_mask), 'omitnan');
                subj_features.Theta_waveform{ti} = env;
            end
        end

        % -----------------------------------------------------------------
        % Sliding-window PLV.
        % -----------------------------------------------------------------
        if ~isempty(EEGp_phase) && ~isempty(acc_idx) && EEGp_phase.trials > 0

            row_idx = (1:n_rows)';
            valid_ep = subj_features.has_eeg_epoch & ...
                subj_features.epoch <= EEGp_phase.trials;

            stage_str = string(subj_features.stage);
            btype_str = string(subj_features.block_type);
            corr_str  = string(subj_features.correct_num);

            valid_bucket = valid_ep & ~ismissing(subj_features.stage) & ...
                ~ismissing(subj_features.block_type) & ~isnan(subj_features.correct_num);

            unique_buckets = unique([stage_str(valid_bucket), ...
                                     btype_str(valid_bucket), ...
                                     corr_str(valid_bucket)], 'rows');

            for ub = 1:size(unique_buckets, 1)

                bucket_mask = valid_bucket & ...
                    stage_str == unique_buckets(ub, 1) & ...
                    btype_str == unique_buckets(ub, 2) & ...
                    corr_str  == unique_buckets(ub, 3);

                bucket_rows = row_idx(bucket_mask);
                if numel(bucket_rows) < MIN_TRIALS_PLV_WINDOW
                    continue;
                end

                [~, ord] = sort(subj_features.trial_continuous(bucket_rows));
                bucket_rows = bucket_rows(ord);
                n_bucket = numel(bucket_rows);

                for bi2 = 1:n_bucket

                    win_lo = max(1, bi2 - PLV_WINDOW_HALF);
                    win_hi = min(n_bucket, bi2 + PLV_WINDOW_HALF);
                    window_rows = bucket_rows(win_lo:win_hi);

                    if numel(window_rows) < MIN_TRIALS_PLV_WINDOW
                        continue;
                    end

                    eps_window = subj_features.epoch(window_rows);
                    eps_window = eps_window(~isnan(eps_window) & ...
                        eps_window >= 1 & eps_window <= EEGp_phase.trials);
                    eps_window = round(eps_window);

                    if numel(eps_window) < MIN_TRIALS_PLV_WINDOW
                        continue;
                    end

                    center_row = bucket_rows(bi2);

                    phi_ref = squeeze(angle(mean(exp(1i * double( ...
                        EEGp_phase.data(acc_idx, :, eps_window))), 1, 'omitnan')))';

                    if ~isempty(par_idx)
                        phi_tgt = squeeze(angle(mean(exp(1i * double( ...
                            EEGp_phase.data(par_idx, :, eps_window))), 1, 'omitnan')))';

                        plv_ts = abs(mean(exp(1i * (phi_ref - phi_tgt)), 1, 'omitnan'));
                        plv_ts = plv_ts - mean(plv_ts(plv_bl), 'omitnan');

                        subj_features.PLV_fp(center_row) = mean(plv_ts(plv_mask), 'omitnan');
                        subj_features.PLV_fp_pairwise(center_row) = subj_features.PLV_fp(center_row);
                    end

                    if ~isempty(som_idx)
                        phi_tgt = squeeze(angle(mean(exp(1i * double( ...
                            EEGp_phase.data(som_idx, :, eps_window))), 1, 'omitnan')))';

                        plv_ts = abs(mean(exp(1i * (phi_ref - phi_tgt)), 1, 'omitnan'));
                        plv_ts = plv_ts - mean(plv_ts(plv_bl), 'omitnan');

                        subj_features.PLV_fs(center_row) = mean(plv_ts(plv_mask), 'omitnan');
                        subj_features.PLV_fs_pairwise(center_row) = subj_features.PLV_fs(center_row);
                    end
                end
            end
        end

        % -----------------------------------------------------------------
        % Append all-trial table.
        % -----------------------------------------------------------------
        all_trials_table = safe_vertcat_tables(all_trials_table, subj_features);

        % -----------------------------------------------------------------
        % Build stage-level summary table.
        % -----------------------------------------------------------------
        subj_stage_table = build_stage_summary_table(subj_features, stage_names, BTYPE_LABELS);
        stage_feature_table = safe_vertcat_tables(stage_feature_table, subj_stage_table);

        % -----------------------------------------------------------------
        % Optional subject QC figure.
        % -----------------------------------------------------------------
        if save_figures
            make_subject_qc_figure(subj_features, EEGp.times, ERP_plot_window, ...
                figure_outpath, subj);
        end

        clear EEGp EEGp_theta EEGp_phase
    end

    % ---------------------------------------------------------------------
    % Save feature tables.
    % ---------------------------------------------------------------------
    if save_tables
        all_trials_file = fullfile(feature_outpath, 'group_table_all_trials_KH_v4_merged.mat');
        stage_file      = fullfile(feature_outpath, 'group_stage_table_features_KH_v4_merged.mat');
        combined_file   = fullfile(feature_outpath, 'group_feature_table_combined_KH_v4_merged.mat');

        save(all_trials_file, 'all_trials_table', '-v7.3');
        save(stage_file, 'stage_feature_table', '-v7.3');
        save(combined_file, 'all_trials_table', 'stage_feature_table', '-v7.3');

        writetable(remove_cell_waveforms_for_csv(all_trials_table), ...
            fullfile(feature_outpath, 'group_table_all_trials_KH_v4_merged.csv'));

        writetable(stage_feature_table, ...
            fullfile(feature_outpath, 'group_stage_table_features_KH_v4_merged.csv'));

        fprintf('\nSaved:\n  %s\n  %s\n  %s\n', ...
            all_trials_file, stage_file, combined_file);
    end
end

fprintf('\nMerged KH preprocessing + outcome feature construction complete.\n');

% =============================================================================
%% LOCAL FUNCTIONS
% =============================================================================

function EEG_loaded = load_first_existing_set(folder, candidates)

EEG_loaded = [];

for ci = 1:numel(candidates)
    f = fullfile(folder, candidates{ci});
    if exist(f, 'file')
        EEG_loaded = pop_loadset(candidates{ci}, folder);
        fprintf('  Loaded %s\n', candidates{ci});
        return;
    end
end

end

function x = pad_or_truncate(x, n)

x = x(:);

if numel(x) < n
    x(end+1:n, 1) = NaN;
elseif numel(x) > n
    x = x(1:n);
end

end

function hw = normalize_hw_filter_label(raw)

hw = string(raw);
hw = strtrim(hw);

if hw == "0-30Hz" || hw == "0_30Hz"
    hw = "0_30Hz";
elseif hw == "0-70Hz" || hw == "0_70Hz"
    hw = "0_70Hz";
elseif hw == "" || ismissing(hw)
    hw = "unknown";
else
    hw = "unknown";
end

end

function trial2epoch_mapped = remap_trial2epoch_to_loaded_file(trial2epoch_raw, valid_ep, n_loaded_trials, using_trimmed_file)

trial2epoch_raw = trial2epoch_raw(:);
trial2epoch_mapped = nan(size(trial2epoch_raw));

if using_trimmed_file && ~isempty(valid_ep)
    valid_ep = valid_ep(:);

    max_old = max(valid_ep);
    old_to_new = nan(max_old, 1);
    old_to_new(valid_ep) = (1:numel(valid_ep))';

    ok = ~isnan(trial2epoch_raw) & ...
        trial2epoch_raw >= 1 & ...
        trial2epoch_raw <= numel(old_to_new) & ...
        ~isnan(old_to_new(round(trial2epoch_raw)));

    trial2epoch_mapped(ok) = old_to_new(round(trial2epoch_raw(ok)));
else
    trial2epoch_mapped = trial2epoch_raw;
end

bad = ~isnan(trial2epoch_mapped) & ...
    (trial2epoch_mapped < 1 | trial2epoch_mapped > n_loaded_trials);
trial2epoch_mapped(bad) = NaN;

end

function flags = read_epoch_artifact_flags(S2, valid_ep, n_loaded_trials)

flags = false(n_loaded_trials, 1);

if isfield(S2, 'epoch_artifact_flag_trimmed') && ~isempty(S2.epoch_artifact_flag_trimmed)

    raw = logical(S2.epoch_artifact_flag_trimmed(:));

    if numel(raw) >= n_loaded_trials
        flags = raw(1:n_loaded_trials);
    else
        flags(1:numel(raw)) = raw;
        warning('epoch_artifact_flag_trimmed shorter than loaded EEG trials. Padding false.');
    end

elseif isfield(S2, 'epoch_artifact_flag') && ~isempty(S2.epoch_artifact_flag)

    raw = logical(S2.epoch_artifact_flag(:));

    if ~isempty(valid_ep) && max(valid_ep) <= numel(raw)
        raw = raw(valid_ep);
    end

    if numel(raw) >= n_loaded_trials
        flags = raw(1:n_loaded_trials);
    else
        flags(1:numel(raw)) = raw;
        warning('epoch_artifact_flag shorter than loaded EEG trials. Padding false.');
    end
else
    fprintf('  No epoch artifact flag found. Defaulting to false.\n');
end

end

function y = local_correct_to_numeric(v)

if isnumeric(v) || islogical(v)
    y = double(v);
    y = y(:);
    return;
end

sv = lower(strtrim(string(v)));
y = nan(numel(sv), 1);

y(sv == "1" | sv == "true" | sv == "correct" | sv == "corr") = 1;
y(sv == "0" | sv == "false" | sv == "incorrect" | ...
  sv == "incorr" | sv == "wrong") = 0;

tmp = str2double(sv);
fillable = isnan(y) & ~isnan(tmp);
y(fillable) = tmp(fillable);

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

function [trial2epoch, diag_out] = KH_align_epochs_with_offset(EEG_out, beh_correct)

% Align behaviour rows to outcome epochs by matching reward/loss event codes.
% This is a robust fallback implementation for KH outcome epochs.
%
% For cohort 1:
%   10/11 typically map to incorrect/correct feedback.
% For cohort 2:
%   31/32/33/34 include true/false feedback variants.
%
% The function tests offsets and chooses the one with best event/behaviour
% agreement. If event codes cannot be decoded, it falls back to sequential
% mapping.

beh_correct = beh_correct(:);
n_trials = numel(beh_correct);
n_epochs = EEG_out.trials;

epoch_correct = nan(n_epochs, 1);

for e = 1:n_epochs
    evtypes = {};
    if isfield(EEG_out, 'epoch') && numel(EEG_out.epoch) >= e && ...
            isfield(EEG_out.epoch(e), 'eventtype')
        tmp = EEG_out.epoch(e).eventtype;
        if iscell(tmp)
            evtypes = tmp;
        else
            evtypes = {tmp};
        end
    end

    evstr = string(evtypes);
    evnum = str2double(evstr);

    % Conservative mapping:
    % 11 and 32/34 are often reward/correct-style outcomes.
    % 10 and 31/33 are often loss/incorrect-style outcomes.
    if any(ismember(evnum, [11 32 34]))
        epoch_correct(e) = 1;
    elseif any(ismember(evnum, [10 31 33]))
        epoch_correct(e) = 0;
    end
end

if all(isnan(epoch_correct))
    n_map = min(n_trials, n_epochs);
    trial2epoch = nan(n_trials, 1);
    trial2epoch(1:n_map) = (1:n_map)';

    diag_out = struct();
    diag_out.n_matched = n_map;
    diag_out.n_trials = n_trials;
    diag_out.n_epochs = n_epochs;
    diag_out.best_offset = 0;
    diag_out.match_rate = n_map / max(1, n_trials);
    diag_out.note = 'No decodable event codes; used sequential mapping.';
    return;
end

offset_range = -20:20;
best_score = -Inf;
best_offset = 0;

for off = offset_range
    score = 0;
    denom = 0;

    for t = 1:n_trials
        e = t + off;
        if e >= 1 && e <= n_epochs && ~isnan(epoch_correct(e)) && ~isnan(beh_correct(t))
            denom = denom + 1;
            score = score + double(epoch_correct(e) == beh_correct(t));
        end
    end

    if denom > 0
        score_rate = score / denom;
    else
        score_rate = -Inf;
    end

    if score_rate > best_score
        best_score = score_rate;
        best_offset = off;
    end
end

trial2epoch = nan(n_trials, 1);
for t = 1:n_trials
    e = t + best_offset;
    if e >= 1 && e <= n_epochs
        trial2epoch(t) = e;
    end
end

diag_out = struct();
diag_out.n_matched = sum(~isnan(trial2epoch));
diag_out.n_trials = n_trials;
diag_out.n_epochs = n_epochs;
diag_out.best_offset = best_offset;
diag_out.match_rate = best_score;
diag_out.note = 'Offset selected by outcome-code agreement.';

end

function [EEG, report] = interpolate_amplitude_bursts(EEG, amp_sd_thresh, grad_sd_thresh, pad_ms, max_run_ms)

% Alignment-preserving in-place interpolation of short high-amplitude bursts.
% No samples are deleted.

data = double(EEG.data);
[n_ch, n_samp] = size(data);

pad_samp = round((pad_ms / 1000) * EEG.srate);
max_run_samp = round((max_run_ms / 1000) * EEG.srate);

n_short = 0;
n_long = 0;
total_repaired_samples = 0;

long_start_s = [];
long_end_s = [];
long_channels = {};

for ch = 1:n_ch

    x = data(ch, :);

    med_x = median(x, 'omitnan');
    robust_sd_x = 1.4826 * median(abs(x - med_x), 'omitnan');
    if robust_sd_x <= 0 || isnan(robust_sd_x)
        robust_sd_x = std(x, 0, 'omitnan');
    end
    if robust_sd_x <= 0 || isnan(robust_sd_x)
        continue;
    end

    dx = [0 diff(x)];
    med_dx = median(dx, 'omitnan');
    robust_sd_dx = 1.4826 * median(abs(dx - med_dx), 'omitnan');
    if robust_sd_dx <= 0 || isnan(robust_sd_dx)
        robust_sd_dx = std(dx, 0, 'omitnan');
    end
    if robust_sd_dx <= 0 || isnan(robust_sd_dx)
        robust_sd_dx = Inf;
    end

    bad = abs(x - med_x) > amp_sd_thresh * robust_sd_x | ...
          abs(dx - med_dx) > grad_sd_thresh * robust_sd_dx;

    if ~any(bad)
        continue;
    end

    bad = pad_binary_mask(bad, pad_samp);
    runs = mask_to_runs(bad);

    for rr = 1:size(runs, 1)
        s = runs(rr, 1);
        e = runs(rr, 2);
        run_len = e - s + 1;

        if run_len > max_run_samp
            n_long = n_long + 1;
            long_start_s(end+1, 1) = s / EEG.srate; %#ok<AGROW>
            long_end_s(end+1, 1) = e / EEG.srate; %#ok<AGROW>
            long_channels{end+1, 1} = EEG.chanlocs(ch).labels; %#ok<AGROW>
            continue;
        end

        left = s - 1;
        right = e + 1;

        if left < 1 && right > n_samp
            continue;
        elseif left < 1
            x(s:e) = x(right);
        elseif right > n_samp
            x(s:e) = x(left);
        else
            x(s:e) = interp1([left right], [x(left) x(right)], s:e, 'linear');
        end

        n_short = n_short + 1;
        total_repaired_samples = total_repaired_samples + run_len;
    end

    data(ch, :) = x;
end

EEG.data = cast(data, class(EEG.data));
EEG = eeg_checkset(EEG);

report = struct();
report.n_short_bursts = n_short;
report.n_long_bursts = n_long;
report.total_repaired_s = total_repaired_samples / EEG.srate;
report.long_run_start_s = long_start_s;
report.long_run_end_s = long_end_s;
report.long_run_channels = long_channels;

end

function padded = pad_binary_mask(mask, pad_samp)

mask = logical(mask(:)');
if pad_samp <= 0
    padded = mask;
    return;
end

kernel = true(1, 2 * pad_samp + 1);
padded = conv(double(mask), double(kernel), 'same') > 0;

end

function runs = mask_to_runs(mask)

mask = logical(mask(:));
d = diff([false; mask; false]);
starts = find(d == 1);
ends = find(d == -1) - 1;
runs = [starts ends];

end


function T = assign_stages_preserve_LE_RN(T, block_col, trial_col, rev_trials_vec, stage_names)

STAGE_LEN = 20;

T.stage = categorical(repmat(missing, height(T), 1), stage_names, 'Ordinal', false);
T.in_stage_window = false(height(T), 1);
T.stage_overlap_resolved = false(height(T), 1);

block_nums = local_to_numeric(T.(block_col));
trial_nums = local_to_numeric(T.(trial_col));

unique_blocks = unique(block_nums(~isnan(block_nums)));
unique_blocks = sort(unique_blocks(:)');

for bi = 1:numel(unique_blocks)

    raw_block  = unique_blocks(bi);
    block_rows = find(block_nums == raw_block);
    if isempty(block_rows), continue; end

    tib = trial_nums(block_rows);
    valid_tib = ~isnan(tib);
    if ~any(valid_tib), continue; end

    max_trial = max(tib(valid_tib), [], 'omitnan');
    if isnan(max_trial) || max_trial < 1, continue; end

    rev_trial = NaN;
    if ~isempty(rev_trials_vec) && bi <= numel(rev_trials_vec)
        rev_trial = rev_trials_vec(bi);
    end

    raw_ln = 1:min(STAGE_LEN, max_trial);
    raw_re = max(1, max_trial - STAGE_LEN + 1):max_trial;

    if isnan(rev_trial)
        ln_trials = raw_ln;
        le_trials = [];
        rn_trials = [];
        re_trials = raw_re;
    else
        le_start = max(1, rev_trial - STAGE_LEN);
        le_end   = max(1, rev_trial - 1);
        rn_start = rev_trial;
        rn_end   = min(max_trial, rev_trial + STAGE_LEN - 1);

        le_trials = le_start:le_end;
        rn_trials = rn_start:rn_end;

        protected = unique([le_trials, rn_trials]);

        ln_trials = setdiff(raw_ln, protected);
        re_trials = setdiff(raw_re, protected);

        overlap_trials = unique([intersect(raw_ln, protected), ...
                                  intersect(raw_re, protected)]);

        if ~isempty(overlap_trials)
            overlap_rows_global = block_rows(ismember(tib, overlap_trials));
            T.stage_overlap_resolved(overlap_rows_global) = true;
        end
    end

    T = set_stage_for_block(T, block_rows, tib, ln_trials, 'LN', stage_names);
    T = set_stage_for_block(T, block_rows, tib, le_trials, 'LE', stage_names);
    T = set_stage_for_block(T, block_rows, tib, rn_trials, 'RN', stage_names);
    T = set_stage_for_block(T, block_rows, tib, re_trials, 'RE', stage_names);
end

end

function T = set_stage_for_block(T, block_rows, tib, trial_list, label, stage_names)

if isempty(trial_list), return; end

local_mask = ismember(tib, trial_list);
global_rows = block_rows(local_mask);

if isempty(global_rows), return; end

T.stage(global_rows) = categorical({label}, stage_names);
T.in_stage_window(global_rows) = true;

end

function Tout = build_stage_summary_table(T, stage_names, btype_labels)

rows = {};

for s = 1:numel(stage_names)
    st = stage_names{s};

    for b = 1:numel(btype_labels)
        bt = btype_labels{b};

        for corr_val = 0:1

            m = string(T.stage) == st & ...
                string(T.block_type) == bt & ...
                T.correct_num == corr_val & ...
                ~T.false_fb;

            if ~any(m)
                continue;
            end

            row = table();

            row.subj_id = T.subj_id(find(m, 1));
            row.subj = T.subj(find(m, 1));
            row.cohort = T.cohort(find(m, 1));
            row.stage = categorical({st}, stage_names);
            row.block_type = categorical({bt}, btype_labels);
            row.correct_num = corr_val;

            row.n_trials_total = sum(m);
            row.n_trials_eeg = sum(m & T.has_eeg_epoch);
            row.n_epoch_artifact_flag = sum(m & T.epoch_artifact_flag);
            row.epoch_artifact_flag_rate = mean(T.epoch_artifact_flag(m), 'omitnan');

            row.hw_filter_label = T.hw_filter_label(find(m, 1));

            row.N2_amp_mean = mean(T.N2_amp(m), 'omitnan');
            row.N2_lat_mean = mean(T.N2_lat(m), 'omitnan');
            row.N2_norm_mean = mean(T.N2_norm(m), 'omitnan');

            row.FRN_mean_amp_mean = mean(T.FRN_mean_amp(m), 'omitnan');
            row.FRN_mean_norm_mean = mean(T.FRN_mean_norm(m), 'omitnan');
            row.FRN_peak_amp_mean = mean(T.FRN_peak_amp(m), 'omitnan');
            row.FRN_peak_lat_mean = mean(T.FRN_peak_lat(m), 'omitnan');
            row.FRN_peak_norm_mean = mean(T.FRN_peak_norm(m), 'omitnan');
            row.FRN_excluded_rate = mean(T.FRN_excluded(m), 'omitnan');

            row.RewP_mean_amp_mean = mean(T.RewP_mean_amp(m), 'omitnan');
            row.RewP_mean_norm_mean = mean(T.RewP_mean_norm(m), 'omitnan');
            row.RewP_peak_amp_mean = mean(T.RewP_peak_amp(m), 'omitnan');
            row.RewP_peak_lat_mean = mean(T.RewP_peak_lat(m), 'omitnan');
            row.RewP_peak_norm_mean = mean(T.RewP_peak_norm(m), 'omitnan');
            row.RewP_excluded_rate = mean(T.RewP_excluded(m), 'omitnan');

            row.P300_amp_mean = mean(T.P300_amp(m), 'omitnan');
            row.P300_peak_lat_mean = mean(T.P300_peak_lat(m), 'omitnan');
            row.P300_norm_mean = mean(T.P300_norm(m), 'omitnan');

            row.Theta_amp_mean = mean(T.Theta_amp(m), 'omitnan');

            row.PLV_fp_mean = mean(T.PLV_fp(m), 'omitnan');
            row.PLV_fs_mean = mean(T.PLV_fs(m), 'omitnan');
            row.PLV_fp_pairwise_mean = mean(T.PLV_fp_pairwise(m), 'omitnan');
            row.PLV_fs_pairwise_mean = mean(T.PLV_fs_pairwise(m), 'omitnan');

            rows{end+1, 1} = row; %#ok<AGROW>
        end
    end
end

if isempty(rows)
    Tout = table();
else
    Tout = rows{1};
    for r = 2:numel(rows)
        Tout = safe_vertcat_tables(Tout, rows{r});
    end
end

end

function T = safe_vertcat_tables(A, B)

if isempty(A)
    T = B;
    return;
end

if isempty(B)
    T = A;
    return;
end

varsA = A.Properties.VariableNames;
varsB = B.Properties.VariableNames;

missingInA = setdiff(varsB, varsA);
missingInB = setdiff(varsA, varsB);

for i = 1:numel(missingInA)
    A.(missingInA{i}) = missing_column_like(B.(missingInA{i}), height(A));
end

for i = 1:numel(missingInB)
    B.(missingInB{i}) = missing_column_like(A.(missingInB{i}), height(B));
end

B = B(:, A.Properties.VariableNames);
T = [A; B];

end

function col = missing_column_like(example, n)

if isnumeric(example)
    col = nan(n, size(example, 2));
elseif islogical(example)
    col = false(n, size(example, 2));
elseif isstring(example)
    col = strings(n, size(example, 2));
elseif iscategorical(example)
    col = categorical(repmat(missing, n, 1), categories(example));
elseif iscell(example)
    col = cell(n, size(example, 2));
else
    col = repmat(missing, n, 1);
end

end

function Tcsv = remove_cell_waveforms_for_csv(T)

Tcsv = T;
vars = Tcsv.Properties.VariableNames;

for i = 1:numel(vars)
    if iscell(Tcsv.(vars{i}))
        Tcsv.(vars{i}) = [];
    end
end

end

function make_subject_qc_figure(T, times, plot_window, fig_folder, subj)

if ~ismember('FCzCz_waveform', T.Properties.VariableNames)
    return;
end

in_win = times >= plot_window(1) & times <= plot_window(2);
tt = times(in_win);

fig = figure('Position', [100 100 1200 700], 'Visible', 'off');
sgtitle(sprintf('%s merged v4 outcome QC', subj), 'Interpreter', 'none');

stage_names = categories(T.stage);
if isempty(stage_names)
    stage_names = {'LN','LE','RN','RE'};
end

for s = 1:min(4, numel(stage_names))
    st = stage_names{s};

    ax = subplot(2, 2, s);
    hold(ax, 'on');
    title(ax, st);

    xline(ax, 0, 'k:');
    yline(ax, 0, 'k:');

    for corr_val = 0:1
        m = string(T.stage) == st & ...
            T.correct_num == corr_val & ...
            T.has_eeg_epoch & ...
            ~T.epoch_artifact_flag & ...
            ~T.false_fb;

        waves = T.FCzCz_waveform(m);
        waves = waves(~cellfun(@isempty, waves));

        if isempty(waves), continue; end

        M = cell2mat(cellfun(@(x) x(:)', waves, 'UniformOutput', false));
        M = M(:, in_win);

        mn = mean(M, 1, 'omitnan');

        if corr_val == 1
            plot(ax, tt, mn, 'LineWidth', 2, 'DisplayName', 'correct');
        else
            plot(ax, tt, mn, '--', 'LineWidth', 2, 'DisplayName', 'incorrect');
        end
    end

    set(ax, 'YDir', 'reverse');
    xlabel(ax, 'Time (ms)');
    ylabel(ax, 'FCz/Cz \muV');
    xlim(ax, plot_window);
    legend(ax, 'Box', 'off');
end

if ~exist(fig_folder, 'dir')
    mkdir(fig_folder);
end

exportgraphics(fig, fullfile(fig_folder, sprintf('%s_merged_v4_outcome_QC.pdf', subj)), ...
    'ContentType', 'vector');

close(fig);

end