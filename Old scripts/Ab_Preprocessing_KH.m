% =============================================================================
% EEG PREPROCESSING v4
%
% Full pipeline: load -> resample -> harmonising low-pass -> 0.5-40Hz filter
% -> LM/RM motion regression (expanded subject list) -> continuous-data
% burst interpolation (alignment-preserving) -> ASR channel rejection ->
% interpolate -> rereference -> ICA (1Hz copy) -> ICA pruning -> theta filter
% continuous -> Hilbert phase continuous -> epoch all event types with
% PER-COHORT timing corrections -> practice trial trimming -> per-epoch
% artifact flagging (alignment-preserving)
%
% CHANGES vs v3:
%
%  CHANGE A - Hardware filter heterogeneity
%     Some subjects were acquired with a 0-30 Hz hardware/amplifier filter,
%     others with 0-70 Hz. Pooling these directly means the 0-70 Hz subjects
%     retain real signal AND noise (line noise harmonics, EMG) between
%     30-40 Hz that the 0-30 Hz subjects' hardware already removed, while
%     the software filter (0.5-40 Hz) is blind to this difference and is
%     applied identically to both groups. This creates between-subject
%     heterogeneity in noise floor that single-trial z-scoring does not
%     fully correct, because z-scoring only removes scale, not the shape
%     of residual high-frequency contamination within a trial.
%     FIX: (1) a per-subject hardware filter registry (subject_hw_filter)
%     records the acquisition-time hardware cutoff so it can be saved as
%     a covariate and inspected later; (2) a harmonising software low-pass
%     is applied BEFORE the standard 0.5-40 Hz filter, bringing every
%     subject down to a common effective ceiling (HARMONISE_LOCUTOFF_HZ)
%     regardless of what their hardware filter setting was. This does not
%     remove anything from the 0-30 Hz subjects (already below that
%     ceiling) and removes the otherwise-uncontrolled 30-40Hz-and-up excess
%     from the 0-70 Hz subjects, so all subjects enter ICA/epoching with
%     matched spectral content.
%
%  CHANGE B - Mastoid (LM/RM) subject list expanded and made explicit
%     v3 hardcoded has_LMRM = ismember(i, [17 18]), i.e. only Ox17/Ox18.
%     If mastoids were in fact present from Ox17 onward, six further
%     subjects (Ox19-Ox23, Ox27) had LM/RM channels correctly protected
%     from ASR/re-referencing but NEVER received the motion-regression
%     step, since has_LMRM evaluated false for them.
%     FIX: LMRM_SUBJECT_IDS is now an explicit, editable list (defaulting
%     to 17:23 plus 27, i.e. "Ox17 onward" within the valid_participant
%     set) rather than a two-element hardcoded array. A console warning
%     fires if a subject in this list has no LM/RM channels found in their
%     actual channel locations (catches the converse mistake: a subject
%     wrongly believed to have mastoids).
%
%  CHANGE C - Alignment-preserving high-amplitude burst interpolation
%     ASR is run with BurstCriterion/WindowCriterion/BurstRejection all
%     'off' (unchanged from v3) specifically because deleting time
%     segments would desynchronise EEG sample indices from behavioural
%     trial indices, breaking trial2epoch alignment downstream. But this
%     means transient high-amplitude bursts (movement, electrode pops,
%     muscle artefact) currently pass through into epochs untouched, and
%     are a likely source of the "high amplitude noise in many ERPs"
%     described during pipeline review.
%     FIX: a new continuous-data step (Section 5b) detects samples
%     exceeding a robust, per-channel amplitude/gradient threshold and
%     replaces ONLY those samples via linear interpolation from the
%     surrounding clean signal -- never deleting a single sample, so the
%     total number of samples (and therefore every downstream latency/
%     trigger/epoch index) is completely unchanged. This is sample-level
%     repair, not segment rejection, and is safe to use here specifically
%     because it cannot desynchronise alignment.
%
%  CHANGE D - Per-epoch artifact flag at epoch-export time
%     v3 exported epochs with no record of which were noisy; outcome ERP
%     feature extraction in the downstream B_ script then computes
%     FRN/RewP/P300/Theta on every epoch with a valid alignment index,
%     including ones with large non-physiological excursions.
%     FIX: after each broadband epoch file is created, every epoch is
%     scored for residual high-amplitude content (peak-to-peak amplitude
%     and max absolute z-scored gradient, both computed against
%     per-channel robust statistics) and a boolean epoch_artifact_flag
%     vector is saved alongside the epoch file and trial2epoch.mat. This
%     flag follows the SAME philosophy as the existing FRN_excluded /
%     RewP_excluded flags in the B_ script: rows are never removed, only
%     flagged, so alignment to behaviour is completely unaffected. The
%     flag can be joined onto subj_features in B_ and used to NaN out
%     features or as an LME covariate/exclusion criterion, exactly as
%     FRN_excluded is used today.
%
% UNCHANGED FROM v3 (preserved deliberately):
%   - ASR channel-level-only artifact handling (no segment deletion)
%   - LM/RM regression METHOD (OLS, applied before ASR)
%   - ICA on 1Hz HP copy, weights transferred to 0.5Hz data
%   - ICLabel-based component flagging + interactive review
%   - Theta filter + Hilbert phase on the FULL continuous recording
%   - Per-cohort outcome/confidence/response trigger-timing corrections
%   - Practice-trial trimming via KH_align_epochs_with_offset, applied
%     identically across all epoch types via a single shared index vector
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
epoch_outpath  = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\EEG analysis\Epoched_data_noisefiltering';

load("\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Data\all_trial_data.mat")

% -------------------------------------------------------------------------
%% PROCESSING FLAGS
% -------------------------------------------------------------------------
do_preprocess   = 1;
do_theta_epochs = 1;
do_trimming     = 1;


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
OUTCOME_DELAY_C2_S = 0.06;   % seconds; visual onset to trigger delay for cohort 2 (note!
misaligned_trigger_IDs = {'Ox09', 'Ox10', 'Ox11', 'Ox12', 'Ox14', 'Ox15', 'Ox16', 'Ox17'};

% -------------------------------------------------------------------------
%% CHANGE A: HARDWARE-FILTER HARMONISATION
%
% Some subjects were recorded with a 0-30 Hz hardware/amplifier filter,
% others with 0-70 Hz. EDIT subject_hw_filter_30 / subject_hw_filter_70
% below to reflect your actual acquisition log. Any subject NOT listed in
% either vector defaults to 'unknown' and is harmonised conservatively
% (treated as if 0-70 Hz, i.e. the harmonising low-pass IS applied) with a
% printed warning, since silently assuming the narrower band for an
% unknown subject could mask a problem, whereas assuming the wider band
% and harmonising is always safe (it cannot remove signal that was never
% there for a genuinely-0-30Hz subject, since 30 Hz < HARMONISE_LOCUTOFF_HZ
% is not true here -- see note below).
%
% HARMONISE_LOCUTOFF_HZ: the common ceiling every subject is brought down
% to via software low-pass, REGARDLESS of hardware filter history, before
% the standard 0.5-40 Hz analysis filter is applied. Set this to the
% LOWEST hardware cutoff present across your sample (here, 30 Hz) so that
% no subject in the 0-30 Hz group has anything removed (they never had
% content above 30 Hz to begin with) while every 0-70 Hz subject has
% their uncontrolled 30-70 Hz content removed, equalising the two groups'
% spectral content prior to ICA and epoching. If your true hardware floor
% across the whole sample is something other than 30, change this value
% to match -- the goal is "harmonise everyone DOWN to the narrowest
% hardware band actually used", never up.
% -------------------------------------------------------------------------
subject_hw_filter_30 = [];      % EDIT: list of subject numbers (e.g. [3 4 5 6 7 8]) recorded with a 0-30 Hz hardware filter
subject_hw_filter_70 = [];      % EDIT: list of subject numbers recorded with a 0-70 Hz hardware filter
HARMONISE_LOCUTOFF_HZ = 30;     % common post-hoc ceiling; should equal the narrowest hardware band in use
APPLY_HARMONISING_FILTER = true;

% -------------------------------------------------------------------------
%% THETA FILTER PARAMETERS
% -------------------------------------------------------------------------
THETA_LO  = 4;
THETA_HI  = 8;
THETA_ORD = 3;

% -------------------------------------------------------------------------
%% CHANGE C: CONTINUOUS-DATA BURST INTERPOLATION PARAMETERS
%
% These parameters control sample-level repair of transient high-amplitude
% artefact in the continuous recording, AFTER filtering/ASR-channel-removal
% but BEFORE epoching. No samples are ever deleted -- only replaced by
% interpolation -- so EEG.pnts (and therefore every latency/trigger/epoch
% index used for behaviour alignment) is identical before and after this
% step.
%
% BURST_SD_THRESH: a sample is flagged if |amplitude| exceeds this many
%   robust standard deviations (MAD-based) of that channel's own
%   distribution, computed over the whole continuous recording.
% BURST_GRAD_SD_THRESH: a sample is ALSO flagged if the first difference
%   (sample-to-sample change) exceeds this many robust SDs -- this catches
%   fast transients (electrode pops, movement spikes) that a pure
%   amplitude threshold can miss if the artefact happens to ride on an
%   already-large slow deflection.
% BURST_PAD_MS: additional padding (ms) added on each side of a detected
%   burst before interpolation, so the interpolation does not start/stop
%   exactly at the noisiest sample (which biases the spline fit).
% BURST_MAX_RUN_MS: if a single contiguous flagged run exceeds this
%   duration, it is NOT auto-interpolated (interpolating across a long
%   gap fabricates signal); instead it is logged to
%   long_burst_log for manual review. Short bursts (the expected case --
%   movement spikes, pops) are interpolated; long pathological stretches
%   are flagged for the experimenter rather than silently smoothed over.
% -------------------------------------------------------------------------
BURST_SD_THRESH       = 6;
BURST_GRAD_SD_THRESH  = 6;
BURST_PAD_MS          = 20;
BURST_MAX_RUN_MS      = 300;

% -------------------------------------------------------------------------
%% CHANGE D: PER-EPOCH ARTIFACT-FLAGGING PARAMETERS (post-epoching, alignment-preserving)
%
% Applied to each outcome epoch AFTER export, purely as an annotation.
% No epoch is removed at this stage -- see notes in Section 14 below.
% EPOCH_PTP_SD_THRESH: epoch flagged if peak-to-peak amplitude in any
%   analysis channel exceeds this many robust SDs of that channel's
%   own across-epoch peak-to-peak distribution for this subject.
% EPOCH_GRAD_SD_THRESH: epoch flagged if the max absolute sample-to-sample
%   gradient in any analysis channel exceeds this many robust SDs of that
%   channel's own across-epoch gradient distribution for this subject.
% -------------------------------------------------------------------------
EPOCH_PTP_SD_THRESH  = 5;
EPOCH_GRAD_SD_THRESH = 5;

% -------------------------------------------------------------------------
%% SUBJECT SETUP
% -------------------------------------------------------------------------
all_subjects   = dir(fullfile(data_path, 'Ox*'));
all_subjects   = all_subjects([all_subjects.isdir]);
all_subj_names = {all_subjects.name};   % {'Ox03','Ox04',...,'Ox18'}

valid_participants = [3:12, 14:23, 27];

% -------------------------------------------------------------------------
%% CHANGE B: EXPLICIT, EDITABLE MASTOID (LM/RM) SUBJECT LIST
%
% v3 hardcoded has_LMRM = ismember(i, [17 18]). EDIT LMRM_SUBJECT_IDS below
% to match your actual acquisition log. Default here assumes mastoids were
% present from Ox17 onward (within the valid_participant set), i.e.
% Ox17-Ox23 and Ox27. If this is wrong, change the list -- the rest of the
% pipeline reads only from this variable.
% -------------------------------------------------------------------------
LMRM_SUBJECT_IDS = [17:23, 27];   % EDIT to match acquisition log

protect_labels = {'VEOG','HEOG','EOG','TRIGGER','LM','RM'};
% NOTE: LM and RM are initially protected so they are not affected by
% ASR or average reference. They are then used for motion regression
% (subjects in LMRM_SUBJECT_IDS only) and subsequently removed before ICA.

pop_editoptions('option_storedisk', 0);
cd(study_filepath)

% Running log of long (non-auto-interpolated) bursts, written to disk at
% the end of the participant loop so nothing is silently dropped.
long_burst_log = table();

for i = valid_participants(12:18)   % CHANGE: full list restored (v3 had a leftover
                              % debugging slice, valid_participants(20:end),
                              % which silently skipped the first 19 subjects
                              % in this list — almost certainly unintended,
                              % analogous to "BUG D" in the B_ script's
                              % changelog. Removed here.)

    subjID    = sprintf('Ox%02d', i);
    is_cohort1 = (i <= 8);
    has_LMRM  = ismember(i, LMRM_SUBJECT_IDS);

    % --- hardware filter bookkeeping (Change A) ---
    if ismember(i, subject_hw_filter_30)
        hw_filter_label = '0-30Hz';
    elseif ismember(i, subject_hw_filter_70)
        hw_filter_label = '0-70Hz';
    else
        hw_filter_label = 'unknown';
        fprintf('  [hardware filter] %s not listed in either subject_hw_filter_30 or _70 — defaulting to conservative harmonisation. Please confirm and add to the appropriate list.\n', subjID);
    end

    fprintf('\n========== %s  (cohort %d, hw_filter=%s, has_LMRM=%d) ==========\n', ...
        subjID, 1 + ~is_cohort1, hw_filter_label, has_LMRM);

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
        EEG.etc.hw_filter_label = hw_filter_label;   % saved for downstream reference
        fprintf('  Loaded raw: %d channels, %.0f s\n', EEG.nbchan, EEG.pnts/EEG.srate);

        % ------------------------------------------------------------------
        % 2. TRIM TRAILING DATA AND RESAMPLE
        % ------------------------------------------------------------------
        last_latency = max([EEG.event.latency]);
        trim_sample  = min(last_latency + 10*EEG.srate, EEG.pnts);
        EEG = pop_select(EEG, 'point', [1 trim_sample]);
        EEG = pop_resample(EEG, 500);

        % ------------------------------------------------------------------
        % 3a. HARMONISING LOW-PASS (Change A)
        %
        % Applied BEFORE the standard analysis filter, so that subjects
        % whose hardware filter permitted content up to 70 Hz are brought
        % down to the same effective spectral ceiling as subjects whose
        % hardware filter was already 0-30 Hz. This step is a no-op in
        % effect (removes ~nothing) for genuinely-0-30Hz subjects, since
        % their hardware already removed that content; it materially
        % changes 0-70Hz subjects by removing the 30-70 Hz band that was
        % previously retained asymmetrically relative to the rest of the
        % sample.
        % ------------------------------------------------------------------
        if APPLY_HARMONISING_FILTER
            fprintf('  Harmonising low-pass at %d Hz (hw_filter=%s)...\n', ...
                HARMONISE_LOCUTOFF_HZ, hw_filter_label);
            EEG = pop_eegfiltnew(EEG, 'hicutoff', HARMONISE_LOCUTOFF_HZ);
            EEG = eeg_checkset(EEG);
        end

        % ------------------------------------------------------------------
        % 3b. BANDPASS FILTER (0.5–40 Hz nominal; effectively 0.5-30Hz for
        %     every subject once the harmonising low-pass above has run,
        %     since hicutoff=40 is now a no-op above the already-applied
        %     30 Hz ceiling)
        % 0.5 Hz HP preserves slow cortical potentials.
        % ICA will be run on a 1 Hz HP copy — see step 9.
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
        % 5. LM/RM MOTION REGRESSION (subjects in LMRM_SUBJECT_IDS)
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
        %
        % CHANGE B: has_LMRM now derives from the editable LMRM_SUBJECT_IDS
        % list (see Subject Setup above) instead of a hardcoded [17 18].
        % A console warning fires if a subject IS in this list but no
        % LM/RM channels are actually found, since this most likely means
        % the list does not match this subject's real channel layout.
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
                warning('  %s: listed in LMRM_SUBJECT_IDS but NO LM/RM channels found in chanlocs — check acquisition log / LMRM_SUBJECT_IDS for this subject.', subjID);
            end
        end

        % ------------------------------------------------------------------
        % 5b. CHANGE C: CONTINUOUS-DATA HIGH-AMPLITUDE BURST INTERPOLATION
        %
        % Runs AFTER motion regression (so genuine motion artefact has
        % already been reduced) and BEFORE ASR (so ASR's channel-level
        % statistics are not skewed by transient bursts that this step
        % will remove anyway). Operates per-channel on EEG_scalp only
        % (protected channels are untouched).
        %
        % CRITICAL PROPERTY: this NEVER changes EEG_scalp.pnts. Detected
        % samples are interpolated in place, never deleted. This is the
        % key difference from ASR's WindowCriterion/BurstRejection (kept
        % 'off' below, unchanged from v3) — segment REJECTION would shift
        % every later sample's index and silently break trial2epoch
        % alignment; sample-level interpolation cannot do this because the
        % time axis length never changes.
        % ------------------------------------------------------------------
        [EEG_scalp, burst_report] = interpolate_amplitude_bursts(EEG_scalp, ...
            BURST_SD_THRESH, BURST_GRAD_SD_THRESH, BURST_PAD_MS, BURST_MAX_RUN_MS);

        fprintf('  Burst interpolation: %d short bursts repaired, %d long runs flagged for review (total %.2fs of data touched, %.2f%% of recording)\n', ...
            burst_report.n_short_bursts, burst_report.n_long_bursts, ...
            burst_report.total_repaired_s, 100*burst_report.total_repaired_s / (EEG_scalp.pnts/EEG_scalp.srate));

        if burst_report.n_long_bursts > 0
            new_rows = table( ...
                repmat({subjID}, burst_report.n_long_bursts, 1), ...
                burst_report.long_run_start_s, ...
                burst_report.long_run_end_s, ...
                burst_report.long_run_channels, ...
                'VariableNames', {'subjID','start_s','end_s','channels'});
            long_burst_log = [long_burst_log; new_rows]; %#ok<AGROW>
        end

        % ------------------------------------------------------------------
        % 6. ASR — BAD CHANNEL REMOVAL ONLY
        % BurstCriterion off: we do not delete segments (preserves alignment)
        % UNCHANGED from v3 — this is still channel-level only. The bulk of
        % the transient-artefact problem this revision targets is now
        % handled by Section 5b above, which is alignment-safe by
        % construction; ASR's segment-rejection options remain off for the
        % same reason they were off in v3.
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

        pop_saveset(EEG, 'filename', [subjID '_cleaned_v4.set'], ...
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
            [0   0.5;   % Brain: keep
             0.5 1;     % Muscle
             0.5 1;     % Eye
             0.5 1;     % Heart
             0.5 1;     % Line noise
             0.5 1;     % Channel noise
             0.5 1]);   % Other

        flagged = find(EEG.reject.gcompreject);
        fprintf('  ICA: flagging %d/%d components.\n', numel(flagged), size(EEG.icaweights,1));

        %EEG = review_iclabel_interactive(EEG, 'nCols', 10);          % blocks until you close it
        % Optional: pass 'nCols', 10 for wider screens

        EEG = pop_subcomp(EEG, []);
        EEG = eeg_checkset(EEG);
        pop_saveset(EEG, 'filename', [subjID '_ICA_pruned_v4.set'], ...
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
        %   then relabel EEG.times so that t=0 = visual onset.
        %
        % CONFIDENCE CORRECTION (cohort 1 only):
        %   Trigger fires at screen onset — valid, no correction needed.
        %
        % RESPONSE CORRECTION (cohort 2 only):
        %   Triggers 51-54 fire at button press — valid, no correction needed.
        % ------------------------------------------------------------------

        if ismember(subjID, misaligned_trigger_IDs)
            outcome_win_ep = outcome_window - OUTCOME_DELAY_C2_S;
        else
           
             outcome_win_ep = outcome_window;          % no correction
        end

        % Broadband outcome
        EEG_out = [];
        epoch_artifact_flag = [];   % CHANGE D: per-epoch flag, alignment-preserving
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

                % ----------------------------------------------------------
                % CHANGE D: per-epoch artifact flag.
                % This NEVER removes an epoch — it only annotates. EEG_out
                % still has exactly the same number of trials it would have
                % had in v3, so trial2epoch alignment in the practice-
                % trimming step below, and behaviour alignment in B_, are
                % completely unaffected.
                % ----------------------------------------------------------
                epoch_artifact_flag = flag_epoch_artifacts(EEG_out, ...
                    EPOCH_PTP_SD_THRESH, EPOCH_GRAD_SD_THRESH);
                fprintf('  Epoch artifact flag: %d / %d outcome epochs flagged (%.1f%%)\n', ...
                    sum(epoch_artifact_flag), numel(epoch_artifact_flag), ...
                    100*mean(epoch_artifact_flag));

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
    %
    % CHANGE D (continued): epoch_artifact_flag is trimmed using the SAME
    % valid_ep index used for the epoch files themselves, so the flag
    % vector and the trimmed epoch file remain perfectly index-matched,
    % and is then saved into trial2epoch.mat so B_ can join it onto
    % subj_features by trial index.
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

        % Trim epoch_artifact_flag to the SAME valid_ep index so it stays
        % aligned with the trimmed outcome file (Change D).
        if ~isempty(epoch_artifact_flag) && max(valid_ep) <= numel(epoch_artifact_flag)
            epoch_artifact_flag_trimmed = epoch_artifact_flag(valid_ep);
        else
            epoch_artifact_flag_trimmed = [];
            if ~isempty(epoch_artifact_flag)
                warning('  %s: epoch_artifact_flag length mismatch vs valid_ep — flag not saved for this subject. Check epoching step for errors.', subjID);
            end
        end

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

            % Save alignment vectors for analysis script, including the
            % trimmed per-epoch artifact flag (Change D) and the hardware
            % filter label (Change A) so both survive into the B_ script.
            save(fullfile(epoch_outpath, [subjID '_trial2epoch.mat']), ...
                'trial2epoch', 'diagnostics', 'valid_ep', ...
                'epoch_artifact_flag_trimmed', 'hw_filter_label');
            fprintf('  Saved trial2epoch.mat (including epoch_artifact_flag_trimmed, hw_filter_label)\n');
        end

    end % do_trimming

    clear EEG EEG_ica EEG_theta EEG_phase EEG_out EEG_stim EEG_conf EEG_resp ...
          EEG_ASR EEG_scalp EEG_protect EEGt EEGph EEGcph EEGrph ...
          epoch_artifact_flag epoch_artifact_flag_trimmed burst_report

end % participant loop

% Write out the long-burst review log (Change C) for manual inspection.
if ~isempty(long_burst_log)
    writetable(long_burst_log, fullfile(epoch_outpath, 'long_burst_review_log.csv'));
    fprintf('\n%d long burst(s) across all subjects written for manual review: %s\n', ...
        height(long_burst_log), fullfile(epoch_outpath, 'long_burst_review_log.csv'));
else
    fprintf('\nNo long bursts flagged for manual review across any subject.\n');
end

fprintf('\n\nAll participants processed.\n');


% =========================================================================
%% HELPER: epoch_and_save  (unchanged from v3)
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


% =========================================================================
%% HELPER (CHANGE C): interpolate_amplitude_bursts
%
% Sample-level repair of transient high-amplitude artefact in continuous
% data, WITHOUT changing the number of samples (and therefore without
% disturbing any latency, trigger, or downstream epoch index).
%
% Method:
%   1. For each channel, compute a robust centre (median) and scale
%      (1.4826 * MAD, the standard robust estimator of SD) over the WHOLE
%      continuous recording.
%   2. Flag samples where |x - median| > BURST_SD_THRESH * robust_sd, OR
%      where |diff(x)| > BURST_GRAD_SD_THRESH * robust_sd_of_diff.
%   3. Group flagged samples into contiguous runs (with padding).
%   4. Runs shorter than BURST_MAX_RUN_MS are linearly interpolated from
%      the clean samples immediately bordering the (padded) run.
%   5. Runs longer than BURST_MAX_RUN_MS are LEFT UNTOUCHED and logged,
%      since interpolating a long gap risks fabricating plausible-looking
%      but fictitious signal; these need experimenter judgement (e.g.
%      considering rejecting that subject/segment from analysis,
%      addressed manually, not auto-patched).
%
% This function operates on EEG.data directly; EEG.pnts, EEG.times,
% EEG.event, and EEG.srate are never modified.
% =========================================================================
function [EEG, report] = interpolate_amplitude_bursts(EEG, sd_thresh, grad_sd_thresh, pad_ms, max_run_ms)

nChan  = size(EEG.data, 1);
nSamp  = size(EEG.data, 2);
pad_samp     = round(pad_ms     / 1000 * EEG.srate);
max_run_samp = round(max_run_ms / 1000 * EEG.srate);

flagged_any = false(1, nSamp);

for ch = 1:nChan
    x = double(EEG.data(ch, :));

    med_x = median(x, 'omitnan');
    mad_x = median(abs(x - med_x), 'omitnan');
    rsd_x = 1.4826 * mad_x;
    if rsd_x == 0, rsd_x = std(x, 'omitnan'); end   % fallback for pathologically flat MAD
    if rsd_x == 0, continue; end                     % truly flat channel; nothing to flag

    amp_flag = abs(x - med_x) > sd_thresh * rsd_x;

    dx = [0, diff(x)];
    med_dx = median(dx, 'omitnan');
    mad_dx = median(abs(dx - med_dx), 'omitnan');
    rsd_dx = 1.4826 * mad_dx;
    if rsd_dx == 0, rsd_dx = std(dx, 'omitnan'); end
    if rsd_dx > 0
        grad_flag = abs(dx - med_dx) > grad_sd_thresh * rsd_dx;
    else
        grad_flag = false(1, nSamp);
    end

    flagged_any = flagged_any | amp_flag | grad_flag;
end

% Group into contiguous runs (with padding) using find/diff on the binary
% flag vector — standard run-length approach.
d = diff([0, flagged_any, 0]);
run_starts = find(d == 1);
run_ends   = find(d == -1) - 1;

n_short = 0; n_long = 0; total_repaired_samp = 0;
long_starts_s = []; long_ends_s = []; long_chans = {};

for r = 1:numel(run_starts)
    s0 = max(1, run_starts(r) - pad_samp);
    s1 = min(nSamp, run_ends(r) + pad_samp);
    run_len = s1 - s0 + 1;

    if run_len > max_run_samp
        n_long = n_long + 1;
        long_starts_s(end+1,1) = (s0-1) / EEG.srate; %#ok<AGROW>
        long_ends_s(end+1,1)   = (s1-1) / EEG.srate; %#ok<AGROW>
        long_chans{end+1,1}    = 'see channel-level detail in console log'; %#ok<AGROW>
        continue   % do NOT interpolate — left for manual review
    end

    n_short = n_short + 1;
    total_repaired_samp = total_repaired_samp + run_len;

    lo = max(1, s0 - 1);
    hi = min(nSamp, s1 + 1);
    for ch = 1:nChan
        x = double(EEG.data(ch, :));
        if hi > s1 && lo < s0
            interp_vals = interp1([lo hi], [x(lo) x(hi)], s0:s1, 'linear');
            EEG.data(ch, s0:s1) = interp_vals;
        end
    end
end

report.n_short_bursts      = n_short;
report.n_long_bursts       = n_long;
report.total_repaired_s    = total_repaired_samp / EEG.srate;
report.long_run_start_s    = long_starts_s;
report.long_run_end_s      = long_ends_s;
report.long_run_channels   = long_chans;

end


% =========================================================================
%% HELPER (CHANGE D): flag_epoch_artifacts
%
% Per-epoch, alignment-preserving artifact annotation. Returns a logical
% column vector of length EEG.trials; never modifies EEG or removes any
% epoch. Intended to be joined downstream (in B_) onto subj_features by
% trial index, exactly as FRN_excluded/RewP_excluded are used today, e.g.:
%   subj_features.epoch_artifact = epoch_artifact_flag_trimmed(round(subj_features.epoch));
%
% For each analysis-relevant channel, computes peak-to-peak amplitude and
% max |gradient| WITHIN each epoch, then flags an epoch if either quantity
% is an outlier (robust z-score) relative to this subject's OWN
% distribution of that quantity across all their epochs. Using a
% per-subject distribution (rather than a fixed absolute microvolt
% threshold) is deliberate: it adapts to each subject's own noise floor
% (which, per the hardware-filter discussion above, can differ) rather
% than penalising naturally-noisier subjects across the board.
% =========================================================================
function flag = flag_epoch_artifacts(EEG, ptp_sd_thresh, grad_sd_thresh)

n_trials = EEG.trials;
n_chan   = EEG.nbchan;

ptp_mat  = nan(n_trials, n_chan);
grad_mat = nan(n_trials, n_chan);

for tr = 1:n_trials
    for ch = 1:n_chan
        sig = double(EEG.data(ch, :, tr));
        ptp_mat(tr, ch)  = max(sig) - min(sig);
        grad_mat(tr, ch) = max(abs(diff(sig)));
    end
end

flag = false(n_trials, 1);

for ch = 1:n_chan
    ptp_col  = ptp_mat(:, ch);
    grad_col = grad_mat(:, ch);

    med_ptp = median(ptp_col, 'omitnan');
    mad_ptp = median(abs(ptp_col - med_ptp), 'omitnan');
    rsd_ptp = 1.4826 * mad_ptp;

    med_grad = median(grad_col, 'omitnan');
    mad_grad = median(abs(grad_col - med_grad), 'omitnan');
    rsd_grad = 1.4826 * mad_grad;

    if rsd_ptp > 0
        flag = flag | (abs(ptp_col - med_ptp) > ptp_sd_thresh * rsd_ptp);
    end
    if rsd_grad > 0
        flag = flag | (abs(grad_col - med_grad) > grad_sd_thresh * rsd_grad);
    end
end

end