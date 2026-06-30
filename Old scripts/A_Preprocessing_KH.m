% =============================================================================
% EEG PREPROCESSING v3
%
% Full pipeline: load → resample → filter → LM/RM motion regression (Ox17-18)
% → ASR channel rejection → interpolate → rereference → ICA (1Hz copy) →
% ICA pruning → theta filter continuous → Hilbert phase continuous →
% epoch all event types with PER-COHORT timing corrections →
% practice trial trimming
%
% KEY COHORT DIFFERENCES CORRECTED HERE:
%
%  Ox03-08  (neuroscan_v2, cohort 1, Neuroscan system):
%    - Outcome trigger string-based ('correct_trueFB' etc.) → numeric 10/11
%    - Trigger fires ~16 ms BEFORE visual onset (correct direction)
%    - Confidence trigger 200 fires AT confidence screen onset (valid)
%    - Response trigger fires BEFORE getResponse() — NOT valid for response ERPs
%    - No false feedback, no stimID field, no isVisual field
%    - Epoch correction: NONE needed for outcome; confidence valid
%
%  Ox09-18  (category_switch v2/v4/v6, cohort 2, ANT EEG):
%    - Outcome trigger numeric 31-34
%    - showOutcome() calls Screen('Flip') THEN trigger fires → trigger is
%      ~1-2 screen frames (~33 ms) AFTER visual onset
%    - Correction: epoch starting OUTCOME_DELAY_C2 earlier, then shift
%      EEG.times forward so t=0 genuinely = visual feedback onset
%    - Confidence trigger 40 fires AFTER confidence_scale_function() returns
%      (i.e. after subject has already made their rating) → NOT valid for
%      confidence-locked EEG. Confidence epochs SKIPPED for cohort 2.
%    - Response triggers 51-55 fire correctly at button press → valid
%    - Ox17 and Ox18 have LM/RM mastoid channels → motion regression applied
%
%  REVERSAL ALIGNMENT:
%    For cohort 2 the task-stored revTrial (first trial under the new rule)
%    is used as ground truth. detect_reversal_KH() validates this and warns
%    if detection differs by >5 trials. For cohort 1 (no stored revTrial)
%    detect_reversal_KH() is used exclusively.
% =============================================================================
% 
clear; close all; clc;

% -------------------------------------------------------------------------
%% PATHS
% -------------------------------------------------------------------------
eeglab_path    = 'C:\Users\khatinova\OneDrive - Nexus365\Pre_2026_Folders\Documents\MATLAB\eeglab2025.1.0';
addpath(eeglab_path);
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

data_path      = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Data';
study_filepath = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\EEG analysis\Winter 2026';
epoch_outpath  = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\EEG analysis\Epoched_data';

load("\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Data\all_trial_data.mat")

% -------------------------------------------------------------------------
%% PROCESSING FLAGS
% -------------------------------------------------------------------------
do_preprocess   = 1;
do_theta_epochs = 1;
do_trimming     = 1;
count_triggers  = 0;

% -------------------------------------------------------------------------
%% EPOCH WINDOWS (relative to trigger, before timing correction)
% -------------------------------------------------------------------------
outcome_window  = [-0.5  2.0];
stim_window     = [-0.5  2.0];
conf_window     = [-2.0  4.0];   % cohort 1 only — see notes above
response_window = [-1.5  1.5];   % cohort 2 only — see notes above
baseline_win    = [-500   0];

% -------------------------------------------------------------------------
%% COHORT-SPECIFIC TIMING CORRECTIONS
%
% OUTCOME (cohort 2 only):
%   The task calls showOutcome() (which flips the screen) then sends the
%   trigger. The delay between visual onset and trigger is approximately
%   one screen frame at the lab monitor refresh rate. ANT EEG system at
%   60Hz display → 16.67ms per frame. Empirical estimate from similar
%   setups is 1-2 frames → 17-33ms. We use 33ms (2 frames) as a
%   conservative correction. The epoch window is shifted back by this
%   amount so that the actually-recorded trigger falls at t = +33ms
%   in the epoch, and we then relabel times so t=0 = visual onset.
%
% CONFIDENCE (cohort 2):
%   Trigger 40 fires after the rating window closes (up to 4s jitter).
%   No valid EEG marker exists for confidence onset in cohort 2.
%   Confidence epoching is therefore SKIPPED for cohort 2.
%
% RESPONSE (cohort 1):
%   In neuroscan_v2 the response trigger fires before getResponse() is
%   called. Response-locked epochs from cohort 1 are therefore invalid
%   and response epoching is SKIPPED for cohort 1.
% -------------------------------------------------------------------------
OUTCOME_DELAY_C2_S = 0.0;   % seconds; visual onset to trigger delay for cohort 2 (note!

% -------------------------------------------------------------------------
%% THETA FILTER PARAMETERS
% -------------------------------------------------------------------------
THETA_LO  = 4;
THETA_HI  = 8;
THETA_ORD = 3;

% -------------------------------------------------------------------------
%% SUBJECT SETUP
% -------------------------------------------------------------------------
all_subjects   = dir(fullfile(data_path, 'Ox*'));
all_subjects   = all_subjects([all_subjects.isdir]);
all_subj_names = {all_subjects.name};   % {'Ox03','Ox04',...,'Ox18'}

valid_participants = [3:12, 14:23,27];

protect_labels = {'VEOG','HEOG','EOG','TRIGGER','LM','RM'};
% NOTE: LM and RM are initially protected so they are not affected by
% ASR or average reference. They are then used for motion regression
% (Ox17, Ox18 only) and subsequently removed before ICA.

pop_editoptions('option_storedisk', 0);
cd(study_filepath)

for i = valid_participants(20:end)

    subjID    = sprintf('Ox%02d', i);
    is_cohort1 = (i <= 8);
    has_LMRM  = ismember(i, [17 18]);   % only these two have LM/RM channels

    fprintf('\n========== %s  (cohort %d) ==========\n', subjID, 1 + ~is_cohort1);

    % Find this subject's directory safely (avoids indexing bug from v1/v2)
    subj_dir_idx = find(strcmp(all_subj_names, subjID));
    if isempty(subj_dir_idx)
        warning('%s: directory not found, skipping.', subjID);
        continue
    end
    subjPath = fullfile(data_path, all_subjects(subj_dir_idx).name);

    if do_preprocess

        % ------------------------------------------------------------------
        % 1. LOAD RAW DATA
        % ------------------------------------------------------------------
        curryFile = dir(fullfile(subjPath, 'Acquisition*.dat'));
        if isempty(curryFile)
            warning('%s: no Curry file found, skipping.', subjID);
            continue
        end
        EEG = loadcurry(fullfile(subjPath, curryFile(1).name), ...
            'KeepTriggerChannel','True', 'CurryLocations','False');
        EEG.subject = subjID;
        fprintf('  Loaded raw: %d channels, %.0f s\n', EEG.nbchan, EEG.pnts/EEG.srate);

        % ------------------------------------------------------------------
        % 2. TRIM TRAILING DATA AND RESAMPLE
        % ------------------------------------------------------------------
        last_latency = max([EEG.event.latency]);
        trim_sample  = min(last_latency + 10*EEG.srate, EEG.pnts);
        EEG = pop_select(EEG, 'point', [1 trim_sample]);
        EEG = pop_resample(EEG, 500);

        % ------------------------------------------------------------------
        % 3. BANDPASS FILTER (0.5–40 Hz)
        % 0.5 Hz HP preserves slow cortical potentials.
        % ICA will be run on a 1 Hz HP copy — see step 7.
        % ------------------------------------------------------------------
        EEG = pop_eegfiltnew(EEG, 'locutoff', 0.5);
        EEG = pop_eegfiltnew(EEG, 'hicutoff', 40);
        EEG = eeg_checkset(EEG);

        % ------------------------------------------------------------------
        % 4. SEPARATE AND PROTECT EOG / TRIGGER / REFERENCE CHANNELS
        % LM and RM are kept in protect_idx so they survive unmodified
        % through ASR and rereferencing, ready for motion regression.
        % ------------------------------------------------------------------
        orig_chanlocs = EEG.chanlocs;
        orig_labels   = {orig_chanlocs.labels};
        all_labels    = {EEG.chanlocs.labels};

        protect_idx = find(ismember(lower(all_labels), lower(protect_labels)));

        if ~isempty(protect_idx)
            EEG_protect = pop_select(EEG, 'channel', protect_idx);
            EEG_scalp   = pop_select(EEG, 'nochannel', protect_idx);
        else
            EEG_scalp   = EEG;
            EEG_protect = [];
        end
        EEG_scalp = eeg_checkset(EEG_scalp);

        % ------------------------------------------------------------------
        % 5. LM/RM MOTION REGRESSION (Ox17 and Ox18 only)
        %
        % Background: LM (left mastoid) and RM (right mastoid) sit on bony
        % prominences close to the response hand (via conducted body motion
        % through the braille device). They pick up low-frequency motion
        % artefacts correlated with button presses and device vibration.
        % Regressing them out of scalp channels removes this common-mode
        % motion noise without distorting EEG signal, because mastoid
        % channels have minimal neural signal of their own.
        %
        % Method: for each scalp channel, compute OLS regression coefficient
        % of (LM+RM)/2 onto the channel signal, subtract the fit.
        % This is equivalent to the standard EOG regression approach
        % (Gratton et al. 1983) applied to motion rather than eye channels.
        % ------------------------------------------------------------------
        if has_LMRM
            lmrm_labels = {'LM','RM'};
            lmrm_idx_full = find(ismember(lower(orig_labels), lower(lmrm_labels)));

            if numel(lmrm_idx_full) >= 1
                fprintf('  LM/RM motion regression for %s...\n', subjID);

                % Build motion reference from full EEG (before scalp split)
                motion_data = EEG.data(lmrm_idx_full, :);
                motion_ref  = mean(motion_data, 1, 'omitnan');   % 1 x time
                motion_ref  = motion_ref - mean(motion_ref);     % zero-mean

                % Regress from all scalp channels (not protect channels)
                scalp_in_full = setdiff(1:EEG.nbchan, protect_idx);
                for ch = scalp_in_full
                    sig  = double(EEG.data(ch, :));
                    % OLS: beta = cov(sig, ref) / var(ref)
                    beta = (sig - mean(sig)) * motion_ref' / (motion_ref * motion_ref');
                    EEG.data(ch,:) = sig - beta * motion_ref;
                end
                EEG_scalp = pop_select(EEG, 'nochannel', protect_idx);
                EEG_scalp = eeg_checkset(EEG_scalp);
                fprintf('  LM/RM regression complete.\n');
            else
                fprintf('  LM/RM not found in channel list for %s — skipping regression.\n', subjID);
            end
        end

        % ------------------------------------------------------------------
        % 6. ASR — BAD CHANNEL REMOVAL ONLY
        % BurstCriterion off: we do not delete segments (preserves alignment)
        % ------------------------------------------------------------------
        EEG_ASR = pop_clean_rawdata(EEG_scalp, ...
            'FlatlineCriterion',  5,     ...
            'ChannelCriterion',   0.8,   ...
            'LineNoiseCriterion', 4,     ...
            'BurstCriterion',     'off', ...
            'WindowCriterion',    'off', ...
            'BurstRejection',     'off');
        EEG_ASR = eeg_checkset(EEG_ASR);

        scalp_before = {EEG_scalp.chanlocs.labels};
        scalp_after  = {EEG_ASR.chanlocs.labels};
        removed      = setdiff(scalp_before, scalp_after, 'stable');
        if ~isempty(removed)
            fprintf('  Interpolating removed channels: %s\n', strjoin(removed, ', '));
            EEG_ASR = pop_interp(EEG_ASR, EEG_scalp.chanlocs, 'spherical');
            EEG_ASR = eeg_checkset(EEG_ASR);
        else
            fprintf('  No channels removed by ASR.\n');
        end

        % ------------------------------------------------------------------
        % 7. REINSERT PROTECTED CHANNELS
        % ------------------------------------------------------------------
        if ~isempty(protect_idx)
            nFullChan = length(orig_chanlocs);
            nSamp     = size(EEG_ASR.data, 2);
            newdata   = zeros(nFullChan, nSamp);
            scalp_labels = {EEG_ASR.chanlocs.labels};
            for c = 1:length(scalp_labels)
                orig_idx = find(strcmp(orig_labels, scalp_labels{c}));
                if ~isempty(orig_idx)
                    newdata(orig_idx, :) = EEG_ASR.data(c, :);
                end
            end
            for p = 1:length(protect_idx)
                newdata(protect_idx(p), :) = EEG_protect.data(p, :);
            end
            EEG.data     = newdata;
            EEG.nbchan   = nFullChan;
            EEG.chanlocs = orig_chanlocs;
        else
            EEG = EEG_ASR;
        end
        EEG = eeg_checkset(EEG);

        % ------------------------------------------------------------------
        % 8. AVERAGE REFERENCE (scalp channels only, once)
        % LM/RM are still in protect_idx so they are excluded.
        % ------------------------------------------------------------------
        if ~isempty(protect_idx)
            EEG = pop_reref(EEG, [], 'exclude', protect_idx);
        else
            EEG = pop_reref(EEG, []);
        end
        EEG = eeg_checkset(EEG);

        pop_saveset(EEG, 'filename', [subjID '_cleaned_v3.set'], ...
            'filepath', study_filepath, 'savemode', 'onefile');
        fprintf('  Saved cleaned continuous data.\n');

        % ------------------------------------------------------------------
        % 9. ICA — run on 1 Hz HP copy, apply weights to 0.5 Hz data
        % ------------------------------------------------------------------
        EEG_ica = pop_eegfiltnew(EEG, 'locutoff', 1.0);
        EEG_ica = pop_runica(EEG_ica, 'extended', 1);

        EEG.icaweights  = EEG_ica.icaweights;
        EEG.icasphere   = EEG_ica.icasphere;
        EEG.icawinv     = EEG_ica.icawinv;
        EEG.icachansind = EEG_ica.icachansind;
        EEG = eeg_checkset(EEG);


        EEG = pop_iclabel(EEG, 'default');
    end
        if ~do_preprocess 
            EEG = pop_loadset([subjID '_ICA_laxASR.set']);
        end
        EEG = pop_icflag(EEG, ...
            [0   0.2;   % Brain: keep
             0.5 1;     % Muscle
             0.5 1;     % Eye
             0.5 1;     % Heart
             0.5 1;     % Line noise
             0.5 1;     % Channel noise
             0.5 1]);   % Other

        flagged = find(EEG.reject.gcompreject);
        fprintf('  ICA: flagging %d/%d components.\n', numel(flagged), size(EEG.icaweights,1));

        EEG = review_iclabel_interactive(EEG, 'nCols', 10);          % blocks until you close it
        % Optional: pass 'nCols', 10 for wider screens

        EEG = pop_subcomp(EEG, []);
        EEG = eeg_checkset(EEG);
        pop_saveset(EEG, 'filename', [subjID '_ICA_pruned_v3.set'], ...
            'filepath', study_filepath, 'savemode', 'onefile');
        fprintf('  ICA pruning complete.\n');

        % ------------------------------------------------------------------
        % 10. THETA FILTER CONTINUOUS DATA
        % Filtering the ENTIRE continuous recording means no filtfilt edge
        % artefacts exist within any epoch window.
        % ------------------------------------------------------------------
        if do_theta_epochs
            fprintf('  Theta-filtering continuous data (%d-%d Hz)...\n', THETA_LO, THETA_HI);
            fs_now = EEG.srate;
            [bf, af] = butter(THETA_ORD, [THETA_LO THETA_HI]/(fs_now/2), 'bandpass');

            scalp_all = find(~ismember(lower({EEG.chanlocs.labels}), lower(protect_labels)));

            EEG_theta = EEG;
            for ch = scalp_all
                EEG_theta.data(ch,:) = filtfilt(bf, af, double(EEG.data(ch,:)));
            end
            EEG_theta = eeg_checkset(EEG_theta);

            % Hilbert phase from the theta-filtered signal
            fprintf('  Extracting instantaneous theta phase...\n');
            EEG_phase = EEG_theta;
            for ch = scalp_all
                analytic = hilbert(double(EEG_theta.data(ch,:)));
                EEG_phase.data(ch,:) = angle(analytic);   % radians [-pi, pi]
            end
            EEG_phase = eeg_checkset(EEG_phase);
            fprintf('  Phase extraction done.\n');
        end

        % ------------------------------------------------------------------
        % 11. DEFINE EVENT CODES PER COHORT
        % ------------------------------------------------------------------
        if is_cohort1
            outcome_codes  = [10 11];                          % no false FB
            stim_codes     = [1 2 3 4 5 6 7 8];               % stim triggers
            conf_codes     = 200;                              % confidence onset fires AT screen
            response_codes = [];                               % SKIP: trigger fires before response
        else
            outcome_codes  = [31 32 33 34];                   % true/false correct/incorrect
            stim_codes     = [11 12 13 14 21 22 23 24];
            conf_codes     = [];                               % SKIP: trigger fires after rating
            response_codes = [51 52 53 54];                   % fires correctly
        end

        % Convert to string if EEG stores types as strings
        evt0 = EEG.event(1).type;
        
        if isnumeric(evt0)
            outcome_arg = outcome_codes;
            stim_arg    = stim_codes;
            conf_arg    = conf_codes;
            resp_arg    = response_codes;
        else
            to_str = @(v) arrayfun(@num2str, v, 'UniformOutput', false);
        
            outcome_arg = to_str(outcome_codes);
            stim_arg    = to_str(stim_codes);
        
            if isempty(conf_codes)
                conf_arg = {};
            else
                conf_arg = to_str(conf_codes);
            end
        
            if isempty(response_codes)
                resp_arg = {};
            else
                resp_arg = to_str(response_codes);
            end
        end

        % ------------------------------------------------------------------
        % 12. EPOCH ALL EVENT TYPES WITH PER-COHORT TIMING CORRECTIONS
        %
        % OUTCOME CORRECTION (cohort 2):
        %   The outcome trigger fires after showOutcome() → trigger is
        %   ~OUTCOME_DELAY_C2_S seconds AFTER visual onset.
        %   Solution: epoch using a window shifted back by this amount,
        %   then relabel EEG.times so that t=0 = estimated visual onset.
        %
        %   Concretely: if normal window is [-0.5, 2.0]s and delay is 0.033s,
        %   we epoch at [-0.533, 1.967]s. The trigger fires at t=0.033s in
        %   this epoch. We then shift times by -0.033s so the trigger falls
        %   at t=0, and visual onset (which preceded the trigger by 0.033s)
        %   would be at t=-0.033s. BUT we want visual onset = t=0, so we
        %   shift times by +0.033s → visual onset goes from t=-0.033s to
        %   t=0.
        %
        %   Net result: epoch window [-0.533, 1.967]s, times relabelled so
        %   the trigger is at t=-0.033s and visual onset is at t=0.
        %
        % CONFIDENCE CORRECTION (cohort 1 only):
        %   Trigger fires at screen onset — valid, no correction needed.
        %
        % RESPONSE CORRECTION (cohort 2 only):
        %   Triggers 51-54 fire at button press — valid, no correction needed.
        % ------------------------------------------------------------------

        if is_cohort1
            outcome_win_ep = outcome_window;          % no correction
        else
            % Shift epoch start back by delay, end back by delay
            outcome_win_ep = outcome_window - OUTCOME_DELAY_C2_S;
        end

        % Broadband outcome
        EEG_out = [];
        if ~isempty(outcome_codes)
            try
                EEG_out = pop_epoch(EEG, outcome_arg, outcome_win_ep, ...
                    'newname', [subjID '_outcome']);
                EEG_out = pop_rmbase(EEG_out, baseline_win);

                % For cohort 2: relabel times so visual onset = t=0
                if ~is_cohort1
                    EEG_out.times  = EEG_out.times  + OUTCOME_DELAY_C2_S*1000;
                    EEG_out.xmin   = EEG_out.times(1)/1000;
                    EEG_out.xmax   = EEG_out.times(end)/1000;
                    fprintf('  Cohort 2: outcome times shifted +%.0f ms (visual onset = t=0)\n', ...
                        OUTCOME_DELAY_C2_S*1000);
                end

                pop_saveset(EEG_out, [subjID '_outcome.set'], epoch_outpath);
                fprintf('  Outcome epochs: %d\n', EEG_out.trials);
            catch ME
                warning('  %s outcome epoching failed: %s', subjID, ME.message);
            end
        end

        % Broadband stimulus
        EEG_stim = epoch_and_save(EEG, stim_arg, stim_window, baseline_win, ...
            [subjID '_stimulus'], epoch_outpath, [subjID '_stimulus']);

        % Broadband confidence (cohort 1 only)
        EEG_conf = [];
        if ~isempty(conf_codes)
            EEG_conf = epoch_and_save(EEG, conf_arg, conf_window, baseline_win, ...
                [subjID '_confidence'], epoch_outpath, [subjID '_confidence']);
        else
            fprintf('  Confidence epoching skipped for cohort 2 (trigger fires after rating).\n');
        end

        % Broadband response (cohort 2 only)
        EEG_resp = [];
        if ~isempty(response_codes)
            EEG_resp = epoch_and_save(EEG, resp_arg, response_window, [], ...
                [subjID '_response'], epoch_outpath, [subjID '_response']);
        else
            fprintf('  Response epoching skipped for cohort 1 (trigger fires before response).\n');
        end

        % Theta amplitude epochs
        if do_theta_epochs
            % Outcome theta — same timing correction as broadband outcome
            if ~isempty(outcome_codes)
                try
                    EEGt = pop_epoch(EEG_theta, outcome_arg, outcome_win_ep, ...
                        'newname', [subjID '_outcome_theta']);
                    EEGt = pop_rmbase(EEGt, baseline_win);
                    if ~is_cohort1
                        EEGt.times = EEGt.times + OUTCOME_DELAY_C2_S*1000;
                        EEGt.xmin  = EEGt.times(1)/1000;
                        EEGt.xmax  = EEGt.times(end)/1000;
                    end
                    pop_saveset(EEGt, [subjID '_outcome_theta.set'], epoch_outpath);
                    fprintf('  Outcome theta epochs: %d\n', EEGt.trials);
                catch ME
                    warning('  %s outcome_theta failed: %s', subjID, ME.message);
                end
            end

            epoch_and_save(EEG_theta, stim_arg, stim_window, baseline_win, ...
                [subjID '_stimulus_theta'], epoch_outpath, [subjID '_stimulus_theta']);

            if ~isempty(conf_codes)
                epoch_and_save(EEG_theta, conf_arg, conf_window, baseline_win, ...
                    [subjID '_confidence_theta'], epoch_outpath, [subjID '_confidence_theta']);
            end

            if ~isempty(response_codes)
                epoch_and_save(EEG_theta, resp_arg, response_window, [], ...
                    [subjID '_response_theta'], epoch_outpath, [subjID '_response_theta']);
            end

            % Phase epochs (no baseline correction — phase is circular)
            if ~isempty(outcome_codes)
                try
                    EEGph = pop_epoch(EEG_phase, outcome_arg, outcome_win_ep, ...
                        'newname', [subjID '_outcome_phase']);
                    if ~is_cohort1
                        EEGph.times = EEGph.times + OUTCOME_DELAY_C2_S*1000;
                        EEGph.xmin  = EEGph.times(1)/1000;
                        EEGph.xmax  = EEGph.times(end)/1000;
                    end
                    pop_saveset(EEGph, [subjID '_outcome_phase.set'], epoch_outpath);
                    fprintf('  Outcome phase epochs: %d\n', EEGph.trials);
                catch ME
                    warning('  %s outcome_phase failed: %s', subjID, ME.message);
                end
            end

            if ~isempty(conf_codes)
                try
                    EEGcph = pop_epoch(EEG_phase, conf_arg, conf_window);
                    pop_saveset(EEGcph, [subjID '_confidence_phase.set'], epoch_outpath);
                catch; end
            end

            if ~isempty(response_codes)
                try
                    EEGrph = pop_epoch(EEG_phase, resp_arg, response_window);
                    pop_saveset(EEGrph, [subjID '_response_phase.set'], epoch_outpath);
                catch; end
            end
        end

   

    % ======================================================================
    %% PRACTICE TRIAL TRIMMING
    %
    % For each epoch type, apply the SAME valid_ep index vector derived from
    % aligning the outcome file to behaviour. This guarantees consistent
    % trial counts across all epoch types.
    %
    % For cohort 2 the beh struct may have 5 or 6 blocks. If 6, the first
    % is a practice block that must be excluded from beh_correct.
    % ======================================================================
    if do_trimming

        fprintf('\n  --- Practice trimming for %s ---\n', subjID);

        beh = all_trial_data.(subjID).trial_data;
        num_blocks = height(beh.correct);
        beh_correct = [];

        if num_blocks < 6
            for b = 1:num_blocks
                beh_correct = [beh_correct, beh.correct(b,:)];
            end
        else
            for b = 2:num_blocks
                beh_correct = [beh_correct, beh.correct(b,:)];
            end
        end

        out_fname = [subjID '_outcome.set'];
        if ~exist('EEG_out','var') || isempty(EEG_out)
            if exist(fullfile(epoch_outpath, out_fname),'file')
                EEG_out = pop_loadset(out_fname, epoch_outpath);
            else
                warning('  %s: outcome file not found, skipping trimming.', subjID);
                continue
            end
        end

        [trial2epoch, diagnostics] = KH_align_epochs_with_offset(EEG_out, beh_correct);
        fprintf('  Alignment: %d/%d (%.1f%%), offset=%d\n', ...
            diagnostics.n_matched, diagnostics.n_trials, ...
            100*diagnostics.match_rate, diagnostics.best_offset);

        valid_ep   = sort(unique(trial2epoch(~isnan(trial2epoch))));
        n_practice = EEG_out.trials - numel(valid_ep);
        fprintf('  Removing %d practice epochs.\n', n_practice);

        if n_practice == 0
            fprintf('  No practice epochs — files unchanged.\n');
        else
            % Collect all epoch types to trim
            trim_types = {'outcome','stimulus'};
            if ~isempty(conf_codes);     trim_types{end+1} = 'confidence'; end
            if ~isempty(response_codes); trim_types{end+1} = 'response';   end
            if do_theta_epochs
                theta_phase_types = {'outcome_theta','outcome_phase'};
                if ~isempty(conf_codes)
                    theta_phase_types = [theta_phase_types, {'confidence_theta','confidence_phase'}];
                end
                if ~isempty(response_codes)
                    theta_phase_types = [theta_phase_types, {'response_theta','response_phase'}];
                end
                trim_types = [trim_types, theta_phase_types, {'stimulus_theta'}];
            end

            for tt = 1:numel(trim_types)
                etype = trim_types{tt};
                fname = [subjID '_' etype '.set'];
                fpath = fullfile(epoch_outpath, fname);
                if ~exist(fpath,'file')
                    fprintf('    %s: not found, skipping.\n', fname);
                    continue
                end
                EEGt = pop_loadset(fname, epoch_outpath);
                if max(valid_ep) > EEGt.trials
                    warning('    %s: index exceeds trials (%d > %d)', fname, max(valid_ep), EEGt.trials);
                    continue
                end
                EEGt = pop_select(EEGt, 'trial', valid_ep);
                trimmed_name = [subjID '_' etype '_trimmed.set'];
                pop_saveset(EEGt, trimmed_name, epoch_outpath);
                fprintf('    Saved %s (%d epochs)\n', trimmed_name, EEGt.trials);
            end

            % Save alignment vectors for analysis script
            save(fullfile(epoch_outpath, [subjID '_trial2epoch.mat']), ...
                'trial2epoch', 'diagnostics', 'valid_ep');
            fprintf('  Saved trial2epoch.mat\n');
        end

    end % do_trimming

    clear EEG EEG_ica EEG_theta EEG_phase EEG_out EEG_stim EEG_conf EEG_resp ...
          EEG_ASR EEG_scalp EEG_protect EEGt EEGph EEGcph EEGrph

end % participant loop

fprintf('\n\nAll participants processed.\n');


% =========================================================================
%% HELPER: epoch_and_save
% =========================================================================
function EEG_ep = epoch_and_save(EEGin, codes, window, bwin, name, outpath, tag)
EEG_ep = [];
if isempty(codes); return; end
try
    EEG_ep = pop_epoch(EEGin, codes, window, 'newname', name);
    if ~isempty(bwin)
        EEG_ep = pop_rmbase(EEG_ep, bwin);
    end
    pop_saveset(EEG_ep, [tag '.set'], outpath);
    fprintf('    %s: %d epochs\n', tag, EEG_ep.trials);
catch ME
    warning('    %s failed: %s', tag, ME.message);
end
end