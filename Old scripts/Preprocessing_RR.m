% =============================================================================
% EEG PREPROCESSING — RR COHORT (Nc subjects, EGI NetStation .mff files)
%
% Follows the same structure as Preprocessing_v3.m for KH cohort.
%
% KEY DIFFERENCES FROM KH COHORT:
%
%  System:       EGI 128-channel + EOG, 500 Hz, NetStation .mff format
%  Trigger type: 4-character NetStation event codes (strings in MFF)
%
%  TRIGGER CODES AND TIMING:
%
%    Stimulus:   'T1','T2','T3','T4' — sent AFTER stimulus loop ends with
%                onset = Wsot (start of loop). NetStation records the event
%                with the true onset time passed explicitly as Wsot, so
%                the marker in the MFF IS at the correct stimulus onset.
%                Epoch correction: NONE needed for stimulus.
%
%    Response:   'Go' — sent AT button press (gotime = GetSecs inside KbCheck
%                loop). This fires correctly. No correction needed.
%                NoGo responses have NO trigger (rtvector=5 but no event).
%                Go-only response epoching is valid; NoGo = absence of 'Go'.
%
%    Outcome:    'reward' / 'punish' — startoutcome = GetSecs is set BEFORE
%                DrawFormattedText+Screen('Flip'), so the trigger fires
%                BEFORE the visual. This is the correct direction (trigger
%                precedes visual onset). No delay correction needed.
%                Unlike KH cohort 2, outcome triggers are semantically
%                labelled as reward/punish rather than correctness × trueFB.
%                See mapping table below.
%
%    Confidence: NO EEG trigger sent. confidence_matrix is stored in the
%                behavioural file only. Confidence epoching SKIPPED.
%
%  OUTCOME TRIGGER MAPPING for this cohort:
%    'reward' = participant received positive feedback (Outcome==1)
%    'punish' = participant received negative feedback (Outcome==-1)
%    These correspond to perceivedCorrect=1 and perceivedCorrect=0 respectively.
%    Whether it was true or false feedback requires cross-referencing with
%    trial_data.trueFB (from Convert_RR_to_trial_data.m).
%
%  EPOCH COUNTS EXPECTED:
%    4 blocks × 100 trials = 400 of each event type per subject.
%    These are printed for verification. No practice blocks for RR subjects.
%
%  CHANNEL COUNT:
%    EGI 128-channel cap. Channel layout: EGI GSN-HydroCel-128.
%    EOG channels: identified by label (e.g. E125, E126, E127, E128 or
%    explicit VEOG/HEOG in the cap file — check per subject).
% =============================================================================

clear; close all; clc;

% -------------------------------------------------------------------------
%% PATHS
% -------------------------------------------------------------------------
eeglab_path = 'C:\Users\khatinova\OneDrive - Nexus365\Pre_2026_Folders\Documents\MATLAB\eeglab2025.1.0';
addpath(eeglab_path);
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% change rr_data_path to the indivdidual det/prob/det to prob subfolders,
% otherwise everything stays the same
rr_data_path   = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch RR\Data\probabilistic';
study_filepath = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch RR\Results\EEG_analysis';
epoch_outpath  = fullfile(study_filepath, 'Epoched_data');

if ~exist(epoch_outpath,'dir'), mkdir(epoch_outpath); end

% EGI channel location file — 128-channel HydroCel
% Adjust path to your EEGLAB/plugins location
egi_chanlocs_file = fullfile(eeglab_path, 'plugins', 'dipfit', 'standard_BEM', ...
    'elec', 'GSN-HydroCel-128.sfp');

% -------------------------------------------------------------------------
%% FLAGS AND PARAMETERS
% -------------------------------------------------------------------------
do_preprocess   = 1;
do_theta_epochs = 1;

THETA_LO  = 4;
THETA_HI  = 8;
THETA_ORD = 3;

% Epoch windows
stim_window    = [-0.5  2.0];
outcome_window = [-0.5  2.0];
response_window = [-1.5  1.5];
baseline_win   = [-500   0];

% -------------------------------------------------------------------------
%% TRIGGER CODES (NetStation 4-char string codes from MFF)
% -------------------------------------------------------------------------
stim_codes    = {'T1  ','T2  ','T3  ','T4  '};   % note: NetStation pads to 4 chars
outcome_codes = {'rewa','puni'};                  % 'reward' and 'punish' truncated to 4
% If MFF stores full strings use: {'reward','punish'}
% Run trigger_count block first to verify exact strings in your files

response_codes = {'Go  '};   % only Go responses trigger; NoGo = absent

% -------------------------------------------------------------------------
%% EOG / REFERENCE CHANNELS TO PROTECT
% EGI 128-channel: peripheral channels used as reference/EOG
% Adjust labels to match your specific cap file
% -------------------------------------------------------------------------
protect_labels = {'E125','E126','E127','E128','VREF','Cz'};
% E125-E128 are typically peri-ocular in 128-channel GSN HydroCel

% -------------------------------------------------------------------------
%% FIND ALL SUBJECT FOLDERS
% -------------------------------------------------------------------------
subj_dirs = dir(fullfile(rr_data_path, 'Nc*'));
subj_dirs = subj_dirs([subj_dirs.isdir]);

fprintf('Found %d RR subjects in %s\n', numel(subj_dirs), rr_data_path);

pop_editoptions('option_storedisk', 0);

% -------------------------------------------------------------------------
%% MAIN LOOP
% -------------------------------------------------------------------------
for s = 1:numel(subj_dirs)

    subjID   = subj_dirs(s).name;   % e.g. 'Nc01_p1_DPDP'
    subjPath = fullfile(rr_data_path, subjID);

    % Parse clean Nc## label for file naming
    nc_match = regexp(subjID, '^Nc\d+', 'match', 'once');
    if isempty(nc_match)
        warning('Could not parse Nc## from folder %s — skipping', subjID);
        continue
    end
    nc_label = nc_match;   % e.g. 'Nc01'

    fprintf('\n========== %s (%s) ==========\n', nc_label, subjID);

    % ------------------------------------------------------------------
    % 1. FIND .MFF FILE
    % The .mff is a folder (NetStation format) inside the subject folder.
    % ------------------------------------------------------------------
    mff_dirs = dir(fullfile(subjPath, '*.mff'));
    if isempty(mff_dirs)
        warning('%s: no .mff folder found — skipping', nc_label);
        continue
    end
    mff_path = fullfile(subjPath, mff_dirs(1).name);
    fprintf('  MFF: %s\n', mff_dirs(1).name);

    if do_preprocess

        % ------------------------------------------------------------------
        % 2. IMPORT MFF
        % pop_mffimport reads the EGI MFF format. The 'code' parameter
        % controls whether event codes are read as strings (0) or numeric (1).
        % Use 0 to preserve the 4-char string codes.
        % ------------------------------------------------------------------
        % try
        %     EEG = pop_mffimport(mff_path, 'typefield', 'code');
        % catch
            % Fallback: try without options if plugin version differs
            EEG = pop_mffimport(mff_path);
        % end
        EEG.subject = nc_label;
        fprintf('  Imported: %d channels, %.0f s, %d Hz\n', ...
            EEG.nbchan, EEG.pnts/EEG.srate, EEG.srate);

        % Print trigger counts to verify 400 per type
        trigger_count_table = count_triggers(EEG);
        fprintf('  --- Trigger counts ---\n');
        disp(trigger_count_table);
        expected_n = 400;
        for tc = 1:height(trigger_count_table)
            code = string(trigger_count_table.EventType(tc));
            cnt  = trigger_count_table.Count(tc);
            if any(contains([stim_codes, outcome_codes, response_codes], strtrim(code)))
                if cnt ~= expected_n
                    warning('  %s: %s has %d epochs', nc_label, code, cnt);
                end
            end
        end

        % ------------------------------------------------------------------
        % 3. LOAD CHANNEL LOCATIONS
        % ------------------------------------------------------------------
        if exist(egi_chanlocs_file, 'file')
            EEG = pop_chanedit(EEG, 'lookup', egi_chanlocs_file);
            EEG = eeg_checkset(EEG);
            fprintf('  Channel locations loaded.\n');
        else
            warning('  EGI chanlocs file not found at %s', egi_chanlocs_file);
            fprintf('  Proceeding without explicit channel locations.\n');
        end

        % ------------------------------------------------------------------
        % 4. TRIM AND RESAMPLE
        % If already at 500 Hz, resample is a no-op.
        % ------------------------------------------------------------------

        % ------------------------------------------------------------------
        % 5. BANDPASS FILTER (0.5–40 Hz)
        % ------------------------------------------------------------------
        EEG = pop_eegfiltnew(EEG, 'locutoff', 0.5);
        EEG = pop_eegfiltnew(EEG, 'hicutoff', 40);
        EEG = eeg_checkset(EEG);

        % ------------------------------------------------------------------
        % 6. SEPARATE EOG / REFERENCE CHANNELS
        % ------------------------------------------------------------------
        orig_chanlocs = EEG.chanlocs;
        orig_labels   = {orig_chanlocs.labels};
        all_labels    = {EEG.chanlocs.labels};

        protect_idx = find(ismember(lower(all_labels), lower(protect_labels)));
        if ~isempty(protect_idx)
            EEG_protect = pop_select(EEG, 'channel', protect_idx);
            EEG_scalp   = pop_select(EEG, 'nochannel', protect_idx);
            fprintf('  Protected %d channels: %s\n', numel(protect_idx), ...
                strjoin({orig_chanlocs(protect_idx).labels}, ', '));
        else
            EEG_scalp   = EEG;
            EEG_protect = [];
            fprintf('  No EOG/reference channels identified — check protect_labels.\n');
        end
        EEG_scalp = eeg_checkset(EEG_scalp);

        % ------------------------------------------------------------------
        % 7. ASR — BAD CHANNEL REMOVAL ONLY
        % No burst rejection (preserves alignment to behaviour).
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
            fprintf('  Interpolating: %s\n', strjoin(removed, ', '));
            EEG_ASR = pop_interp(EEG_ASR, EEG_scalp.chanlocs, 'spherical');
            EEG_ASR = eeg_checkset(EEG_ASR);
        else
            fprintf('  No channels removed by ASR.\n');
        end

        % ------------------------------------------------------------------
        % 8. REINSERT PROTECTED CHANNELS
        % ------------------------------------------------------------------
        if ~isempty(protect_idx)
            nFullChan = length(orig_chanlocs);
            nSamp     = size(EEG_ASR.data, 2);
            newdata   = zeros(nFullChan, nSamp);
            for c = 1:length({EEG_ASR.chanlocs.labels})
                oi = find(strcmp(orig_labels, EEG_ASR.chanlocs(c).labels));
                if ~isempty(oi), newdata(oi,:) = EEG_ASR.data(c,:); end
            end
            for p = 1:length(protect_idx)
                newdata(protect_idx(p),:) = EEG_protect.data(p,:);
            end
            EEG.data     = newdata;
            EEG.nbchan   = nFullChan;
            EEG.chanlocs = orig_chanlocs;
        else
            EEG = EEG_ASR;
        end
        EEG = eeg_checkset(EEG);

        % ------------------------------------------------------------------
        % 9. AVERAGE REFERENCE (scalp channels only)
        % ------------------------------------------------------------------
        if ~isempty(protect_idx)
            EEG = pop_reref(EEG, [], 'exclude', protect_idx);
        else
            EEG = pop_reref(EEG, []);
        end
        EEG = eeg_checkset(EEG);

        pop_saveset(EEG, 'filename', [nc_label '_cleaned.set'], ...
            'filepath', study_filepath, 'savemode', 'onefile');
        fprintf('  Saved cleaned.\n');

        % ------------------------------------------------------------------
        % 10. ICA — run on 1 Hz HP copy, weights back to 0.5 Hz data
        % ------------------------------------------------------------------
        EEG_ica = pop_eegfiltnew(EEG, 'locutoff', 1.0);
        EEG_ica = pop_runica(EEG_ica, 'extended', 1);

        EEG.icaweights  = EEG_ica.icaweights;
        EEG.icasphere   = EEG_ica.icasphere;
        EEG.icawinv     = EEG_ica.icawinv;
        EEG.icachansind = EEG_ica.icachansind;
        EEG = eeg_checkset(EEG);

        EEG = pop_iclabel(EEG, 'default');
        EEG = pop_icflag(EEG, ...
            [0   0.2;   % Brain
             0.8 1;     % Muscle
             0.8 1;     % Eye
             0.8 1;     % Heart
             0.8 1;     % Line noise
             0.8 1;     % Channel noise
             0.8 1]);   % Other

        flagged = find(EEG.reject.gcompreject);
        fprintf('  ICA: flagging %d/%d components.\n', numel(flagged), size(EEG.icaweights,1));

        %% ---------- OPTIONAL : MANUAL ICA rejection -----------------
        % EEG = review_iclabel_interactive(EEG);          % blocks until you close it
        % Optional: pass 'nCols', 10 for wider screens


        pop_saveset(EEG, 'filename', [nc_label '_ICA.set'], ...
            'filepath', study_filepath, 'savemode', 'onefile');
        EEG = pop_subcomp(EEG, []);
        EEG = eeg_checkset(EEG);
        pop_saveset(EEG, 'filename', [nc_label '_ICA_pruned.set'], ...
            'filepath', study_filepath, 'savemode', 'onefile');
        fprintf('  ICA pruning done.\n');

        % ------------------------------------------------------------------
        % 11. THETA FILTER + HILBERT PHASE ON CONTINUOUS DATA
        % ------------------------------------------------------------------
        if do_theta_epochs
            fprintf('  Theta filtering (%d-%d Hz)...\n', THETA_LO, THETA_HI);
            fs_now = EEG.srate;
            [bf, af] = butter(THETA_ORD, [THETA_LO THETA_HI]/(fs_now/2), 'bandpass');

            scalp_all = find(~ismember(lower({EEG.chanlocs.labels}), lower(protect_labels)));

            EEG_theta = EEG;
            for ch = scalp_all
                EEG_theta.data(ch,:) = filtfilt(bf, af, double(EEG.data(ch,:)));
            end
            EEG_theta = eeg_checkset(EEG_theta);

            EEG_phase = EEG_theta;
            for ch = scalp_all
                EEG_phase.data(ch,:) = angle(hilbert(double(EEG_theta.data(ch,:))));
            end
            EEG_phase = eeg_checkset(EEG_phase);
            fprintf('  Theta filter and phase extraction done.\n');
        end

        % ------------------------------------------------------------------
        % 12. EPOCH — STIMULUS
        % Trigger: 'T1','T2','T3','T4' sent with onset = Wsot (correct).
        % No timing correction needed.
        % ------------------------------------------------------------------
        EEG_stim = epoch_and_save_rr(EEG, stim_codes, stim_window, ...
            baseline_win, nc_label, epoch_outpath, 'stimulus', 400);

        % ------------------------------------------------------------------
        % 13. EPOCH — OUTCOME
        % 'reward'/'punish' sent at startoutcome = GetSecs BEFORE Screen flip.
        % Trigger is therefore BEFORE visual onset — correct direction.
        % No timing correction needed.
        % ------------------------------------------------------------------
        EEG_out = epoch_and_save_rr(EEG, outcome_codes, outcome_window, ...
            baseline_win, nc_label, epoch_outpath, 'outcome', 400);

        % ------------------------------------------------------------------
        % 14. EPOCH — RESPONSE (Go only)
        % 'Go' sent at button press time — correct.
        % Only Go trials trigger; NoGo trials have no marker.
        % Expected count < 400 (depends on accuracy rate).
        % NO BASELINE CORRECTION because I'll look at pre-response ramp in
        % activity
        % ------------------------------------------------------------------
        epoch_and_save_rr(EEG, response_codes, response_window, ...
            [], nc_label, epoch_outpath, 'response', NaN);

        % No confidence epoching — no EEG trigger in RR task.
        fprintf('  Confidence epoching skipped (no EEG trigger in RR task).\n');

        % ------------------------------------------------------------------
        % 15. THETA AND PHASE EPOCHS
        % ------------------------------------------------------------------
        if do_theta_epochs
            epoch_and_save_rr(EEG_theta, stim_codes, stim_window, ...
                baseline_win, nc_label, epoch_outpath, 'stimulus_theta', 400);
            epoch_and_save_rr(EEG_theta, outcome_codes, outcome_window, ...
                baseline_win, nc_label, epoch_outpath, 'outcome_theta', 400);

            % Phase: no baseline (circular data)
            epoch_and_save_rr(EEG_phase, stim_codes, stim_window, ...
                [], nc_label, epoch_outpath, 'stimulus_phase', 400);
            epoch_and_save_rr(EEG_phase, outcome_codes, outcome_window, ...
                [], nc_label, epoch_outpath, 'outcome_phase', 400);
        end

    end % do_preprocess

    clear EEG EEG_ica EEG_theta EEG_phase EEG_out EEG_stim ...
          EEG_ASR EEG_scalp EEG_protect

end % subject loop

fprintf('\n\nAll RR subjects processed.\n');




% =========================================================================
%% HELPER: count_triggers
% =========================================================================
function T = count_triggers(EEG)
types = {EEG.event.type};
types_str = cellfun(@(x) string(x), types, 'UniformOutput', true);
[u, ~, ic] = unique(types_str);
counts = accumarray(ic, 1);
T = table(u(:), counts(:), 'VariableNames', {'EventType','Count'});
T = sortrows(T, 'Count', 'descend');
end