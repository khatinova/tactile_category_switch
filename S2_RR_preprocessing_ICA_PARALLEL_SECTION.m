% =============================================================================
% REFACTORED: ICA SECTION WITH PARFOR PARALLELIZATION
% Replace the existing "if RUN_ICA" block (starting around line 280) with this
% =============================================================================

if RUN_ICA

    % -------------------------------------------------------------------------
    % Setup Parallel Pool
    % -------------------------------------------------------------------------
    % Check if parallel pool exists; if not, create one with default workers
    if isempty(gcp('nocreate'))
        pool = parpool('local');  % Default: uses all available cores
    end
    
    % For explicit control, use:
    % delete(gcp('nocreate'));  % Close any existing pool
    % parpool('local', 4);       % Use 4 workers
    
    % -------------------------------------------------------------------------
    % Parallel ICA Processing Loop
    % -------------------------------------------------------------------------
    % Each subject is processed independently by a separate worker.
    % IMPORTANT: parfor has restrictions on variable usage (see comments below).
    %
    % Broadcast variables (read-only on all workers):
    %   - study_filepath: Directory where cleaned sets are loaded/saved
    %   - ICLABEL_MIN_BRAIN_PROB: Threshold for ICLabel classification
    %
    % Sliced variable (distributed across workers):
    %   - s_idx: Loop index (parfor automatically distributes)
    %
    % No fprintf/disp statements (parfor restriction).
    % All file I/O uses subject-specific identifiers (no collisions).
    % -------------------------------------------------------------------------
    
    parfor s_idx = 4:numel(rr_subjects)

        subjID   = rr_subjects(s_idx).nc_label;   % 'Nc##'
        subjPath = rr_subjects(s_idx).subjPath;   % (currently unused in this loop)

        % -----------------------------------------------------------------
        % Load cleaned preprocessed data (from RUN_PREPROCESSING stage).
        % -----------------------------------------------------------------
        EEG = pop_loadset(fullfile(study_filepath, [subjID '_cleaned_v4_merged.set']));

        % -----------------------------------------------------------------
        % 9. ICA on 1 Hz copy, transfer weights to 0.5 Hz data.
        %    (ICA is memory-intensive, so we filter to 1 Hz to reduce rank).
        % -----------------------------------------------------------------
        EEG_ica = pop_eegfiltnew(EEG, 'locutoff', 1.0);
        EEG_ica = pop_runica(EEG_ica, 'extended', 1);

        % Transfer ICA decomposition from 1 Hz copy back to original 0.5 Hz data
        EEG.icaweights  = EEG_ica.icaweights;
        EEG.icasphere   = EEG_ica.icasphere;
        EEG.icawinv     = EEG_ica.icawinv;
        EEG.icachansind = EEG_ica.icachansind;
        EEG             = eeg_checkset(EEG);

        % -----------------------------------------------------------------
        % Classify ICA components using ICLabel.
        % -----------------------------------------------------------------
        EEG = pop_iclabel(EEG, 'default');

        % -----------------------------------------------------------------
        % ICLabel auto-rejection.
        % Apply pop_icflag to mark components for rejection, then remove them.
        % Threshold: Brain probability < ICLABEL_MIN_BRAIN_PROB (default 0.50)
        % -----------------------------------------------------------------
        
        % Apply ICLabel flag thresholds
        EEG = pop_icflag(EEG, ...
            [0   0.2;   % Brain
             0.5 1;     % Muscle
             0.5 1;     % Eye
             0.5 1;     % Heart
             0.5 1;     % Line noise
             0.5 1;     % Channel noise
             0.5 1]);   % Other
        
        % Remove components flagged as non-brain
        EEG = pop_subcomp(EEG, []);

        % Extract brain probabilities and identify low-confidence components
        if ~isfield(EEG.etc, 'ic_classification') || ...
           ~isfield(EEG.etc.ic_classification, 'ICLabel') || ...
           ~isfield(EEG.etc.ic_classification.ICLabel, 'classifications')
            error('%s: ICLabel classifications not found after pop_iclabel.', subjID);
        end

        ic_probs = EEG.etc.ic_classification.ICLabel.classifications;
        brain_prob = ic_probs(:, 1);

        % Identify components with brain probability below threshold
        reject_comps = find(brain_prob < ICLABEL_MIN_BRAIN_PROB);

        % Store rejection metadata for QC/reporting
        EEG.etc.iclabel_auto_reject.brain_prob = 0.5;
        EEG.etc.iclabel_auto_reject.threshold = ICLABEL_MIN_BRAIN_PROB;
        EEG.etc.iclabel_auto_reject.rejected_components = reject_comps(:);

        % Remove low-confidence components
        if ~isempty(reject_comps)
            EEG = pop_subcomp(EEG, reject_comps, 0);
        end
        EEG = eeg_checkset(EEG);

        % -----------------------------------------------------------------
        % Save ICA-pruned dataset.
        % Subject-specific filename prevents collisions across parallel workers.
        % -----------------------------------------------------------------
        pop_saveset(EEG, 'filename', [subjID '_ICA_pruned_June26.set'], ...
            'filepath', study_filepath, 'savemode', 'onefile');

        % -----------------------------------------------------------------
        % 10. Theta-filter continuous data and Hilbert phase.
        % -----------------------------------------------------------------
        if do_theta_epochs
            % (Same as sequential version - processed here because it depends
            %  on the ICA-cleaned EEG and is part of the same data pipeline)

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
        outcome_codes = {'rewa','puni'};
        should_shift_outcome = ismember(subjID, misaligned_trigger_IDs);

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
            EEG_out.times = EEG_out.times + OUTCOME_DELAY_C2_S * EEG_out.srate;
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
                EEG_theta_out.times = EEG_theta_out.times + OUTCOME_DELAY_C2_S * EEG_theta_out.srate;
                EEG_theta_out.xmin  = EEG_theta_out.xmin  + OUTCOME_DELAY_C2_S;
                EEG_theta_out.xmax  = EEG_theta_out.xmax  + OUTCOME_DELAY_C2_S;

                EEG_phase_out.times = EEG_phase_out.times + OUTCOME_DELAY_C2_S * EEG_phase_out.srate;
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

            trial2epoch_original = [1:400]';

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

            epoch_artifact_flag_trimmed = flag_outcome_epochs( ...
                EEG_out_trim, ...
                EPOCH_PTP_SD_THRESH, ...
                EPOCH_GRAD_SD_THRESH, ...
                protect_labels);

            trial2epoch = trial2epoch_original(:);
            trial2epoch_out = trial2epoch; %#ok<NASGU>

            save(fullfile(epoch_outpath, [subjID '_trial2epoch.mat']), ...
                'trial2epoch', ...
                'trial2epoch_out', ...
                'valid_ep', ...
                'epoch_artifact_flag_trimmed');

            fprintf('  Saved trimmed outcome files and trial2epoch cache.\n');
        end

        clear EEG EEG_ica EEG_theta EEG_phase EEG_out EEG_theta_out EEG_phase_out
    end

    % -------------------------------------------------------------------------
    % (End of parfor loop)
    % -------------------------------------------------------------------------
    % Note: fprintf statements above still appear because they can be used
    % after parfor (in the main thread) but NOT inside the parfor loop itself.
    % If you need progress updates inside the parfor loop, remove all fprintf
    % statements or use dataQueue (advanced feature).
    % -------------------------------------------------------------------------

    if ~isempty(long_burst_log)
        writetable(long_burst_log, fullfile(study_filepath, 'long_burst_log_v4_merged.csv'));
        save(fullfile(study_filepath, 'long_burst_log_v4_merged.mat'), 'long_burst_log');
    end
end
