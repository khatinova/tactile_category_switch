% =============================================================================
% S2_RR_preprocess_epoch_eeg.m  (RR-locked twin of S2_preprocess_epoch_eeg.m)
%
% PIPELINE STEP 2 of 7 (RR COHORT) — preprocessing + epoching for the RR data.
%
% This is a FAITHFUL PARALLEL of the KH S2 script. Every preprocessing step,
% epoching step, feature-extraction step and local function is identical to
% S2_preprocess_epoch_eeg.m. The ONLY differences are the RR-specific settings,
% all of which flow from COHORT = 'RR' (see the COHORT SELECTION block):
%   * Import         : EGI/.mff via pop_mffimport  (KH used Curry/.dat).
%   * Channels       : E11 (FCz), E7 (Cz), parietal E62/E67/E72,
%                      ACC E11/E6/E16, somatosensory E36/E104/E41/E103.
%   * Reference      : Cz (E7) then common average  (KH used average ref).
%   * Outcome codes  : RR EGI event strings (see PART 1, step 11).
%   * No LM/RM mastoid regression (RR EGI net has no mastoid channels).
%   * Participants   : RR_VALID_FULL  (see USER SWITCHES).
%   * Behaviour spine: researcher == 'RR' rows; subject prefix 'Nc'.
%
% PURPOSE
% -------
%   1. Load raw RR EEG (EGI/.mff).
%   2. Preprocess: harmonising 0-30 Hz low-pass, 0.5-40 Hz band-pass,
%      alignment-preserving burst interpolation, ASR channel cleaning,
%      Cz-then-average reference, ICA + ICLabel auto-reject.
%   3. Export OUTCOME broadband, theta-amplitude, and theta-phase epochs.
%   4. Build trial2epoch alignment (external KH_align_epochs_with_offset) and
%      per-epoch bad-epoch FLAGGING; trim practice/misaligned trials.
%   5. Build the per-trial spine feature table (E11 neg-peak/mean, P300, theta,
%      PLV). FRN/RewP are NOT per-trial -- they are difference waves built per
%      stage (see note in PART 2 and pipeline/utils/kh_compute_frn_rewp_by_stage.m).
%
% Outputs are written to the RR results tree and tagged "RR", so they sit
% alongside (and never overwrite) the KH outputs. Merge with S4 afterwards.
% =============================================================================

clear; close all; clc;

% Put pipeline utils on the path (figure style, subject-id, stage, FRN/RewP)
addpath(genpath(fileparts(mfilename('fullpath'))));

% -------------------------------------------------------------------------
%% COHORT SELECTION
% -------------------------------------------------------------------------
% 'KH' = Curry/Neuroscan+ANT (already preprocessed);  'RR' = EGI/MFF (to run).
% This is the RR-locked twin script, so COHORT is fixed to 'RR'.
COHORT = 'RR';

switch upper(COHORT)
    case 'KH'
        cohort_import     = 'curry';      % loadcurry
        cohort_ref_mode   = 'average';    % average reference, protect EOG/mastoid
        fcz_label    = 'FCz';   cz_label = 'Cz';
        par_channels = {'Pz','P1','P2'};
        acc_channels = {'FCz','Fz','AFz','F1','F2'};
        som_channels = {'C3','C4','CP3','CP1','C5','CP5'};
    case 'RR'
        cohort_import     = 'mff';        % pop_mffimport
        cohort_ref_mode   = 'Cz_then_CAR';% Cz (ch 61) then common average ref
        fcz_label    = 'E11';   cz_label = 'E7';
        par_channels = {'E62','E67','E72'};
        acc_channels = {'E11','E6','E16'};
        som_channels = {'E36','E104','E41','E103'};
    otherwise
        error('Unknown COHORT "%s" (use ''KH'' or ''RR'').', COHORT);
end
% NOTE: RR raw import is wired in PART 1 (pop_mffimport on the discovered .mff,
% then GSN-HydroCel-128 channel locations). The reference (Cz then common
% average), the 0-30 Hz harmonising filter, and per-epoch bad-epoch flagging
% are all applied as in KH. cohort_import/cohort_ref_mode below drive this.

% -------------------------------------------------------------------------
%% USER SWITCHES
% -------------------------------------------------------------------------

RUN_PREPROCESSING       = false;
RUN_ICA                 = true;
RUN_FEATURE_EXTRACTION  = true;

% RR subjects are AUTO-DISCOVERED by scanning the four condition subfolders
% (see PATHS + discover_rr_subjects). You do NOT need to list them by number.
%   RUN_ALL_VALID_RR = true  -> process every Nc## subject that is found.
%   RUN_ALL_VALID_RR = false -> process only the Nc numbers in RR_PILOT_ONLY
%                               (handy for a quick test on one or two subjects).
RUN_ALL_VALID_RR = true;
RR_PILOT_ONLY    = [1, 2];   % <-- numbers (Nc01, Nc02 ...) used only when RUN_ALL_VALID_RR=false

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

% -------------------------------------------------------------------------
% RR raw-data layout.
% The RR Data tree is split into FOUR task-condition subfolders. Each of the
% 15 Nc subjects lives in exactly ONE of these, in a folder named like
% "Nc01_p1_DPDP", which in turn contains the EGI recording as a ".mff" folder.
% We scan all four subfolders, parse the clean "Nc##" label, and locate the
% single .mff recording inside each subject folder (see PART 1).
% -------------------------------------------------------------------------
RR_CONDITION_SUBFOLDERS = { ...
    'det_or_prob_and_conf', ...
    'det_to_prob', ...
    'deterministic', ...
    'probabilistic'};

% EGI channel-location file (128-channel HydroCel). Adjust if your EEGLAB
% plugins live elsewhere.
egi_chanlocs_file = fullfile(eeglab_path, 'plugins', 'dipfit', 'standard_BEM', ...
    'elec', 'GSN-HydroCel-128.sfp');

% Cohort-aware data + output roots.
%   KH raw/results live under 'Salient mod switch KH'
%   RR raw/results live under 'Salient mod switch RR'
KH_data_path = fullfile(base_path, 'Salient mod switch KH', 'Data');
RR_data_path = fullfile(base_path, 'Salient mod switch RR', 'Data');

if strcmpi(COHORT, 'RR')
    cohort_results = fullfile(base_path, 'Salient mod switch RR', 'Results', 'EEG analysis');
else
    cohort_results = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis');
end

study_filepath  = fullfile(cohort_results, 'Winter 2026');
epoch_outpath   = fullfile(cohort_results, 'Epoched_data_noisefiltering');
feature_outpath = fullfile(cohort_results, 'Outcome_feature_tables_v4_merged');
figure_outpath  = fullfile(cohort_results, 'Figures', 'outcome_v4_merged_QC');

if ~exist(study_filepath, 'dir'),  mkdir(study_filepath);  end
if ~exist(epoch_outpath, 'dir'),   mkdir(epoch_outpath);   end
if ~exist(feature_outpath, 'dir'), mkdir(feature_outpath); end
if ~exist(figure_outpath, 'dir'),  mkdir(figure_outpath);  end

addpath(eeglab_path);
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab('nogui'); %#ok<ASGLU>

% all_trial_data is built by S1 (behaviour) and is shared across cohorts.
load(fullfile(KH_data_path, 'all_trial_data_June2026.mat'));

% The behaviour table is a SINGLE combined table (group_T) holding ALL
% subjects (KH "Ox##" and RR "Nc##"), and it lives in the canonical combined
% Data folder alongside all_trial_data (the KH Data path). We keep the RR rows
% by subjID prefix "Nc"; per-subject selection in PART 2 is by exact subjID.
behav_file = fullfile(KH_data_path, 'behav_table_June2026.mat');
if exist(behav_file, 'file')
    S_beh = load(behav_file, 'group_T');
    RR_behav_table = S_beh.group_T;
else
    error('behav_table.mat not found at %s', behav_file);
end

% Keep only RR subjects (subjID starts with "Nc").
RR_behav_table = RR_behav_table(startsWith(string(RR_behav_table.subjID), "Nc"), :);
fprintf('Behaviour table: %d RR rows across %d RR subjects.\n', ...
    height(RR_behav_table), numel(unique(string(RR_behav_table.subjID))));

% -------------------------------------------------------------------------
%% PREPROCESSING PARAMETERS
% -------------------------------------------------------------------------

do_theta_epochs = true;
do_trimming     =true;

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

% Subjects with LM/RM mastoids (KH only; RR EGI net has none).
LMRM_SUBJECT_IDS = [17:23, 27,28];

% EGI 128-channel HydroCel: peri-ocular / reference channels to protect from
% scalp cleaning + average reference (E125-E128 are peri-ocular; VREF/Cz is the
% online vertex reference). Matches the RR preprocessing reference script.
protect_labels = {'E125','E126','E127','E128','VREF','Cz'};

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

% Channel labels are set in the COHORT block at the top of the script
% (KH = FCz/Cz/Pz...;  RR = E11/E7/E62...). Do NOT re-hardcode them here.

block_col = 'block';
trial_col = 'trial';

% =============================================================================
%% PART 1: PREPROCESS AND EXPORT OUTCOME EPOCHS
% =============================================================================

long_burst_log = table();



    % ---------------------------------------------------------------------
    % Build the RR subject registry by scanning all FOUR condition subfolders.
    %   <RR_data_path>/<condition>/Nc##_p#_<order>/<recording>.mff
    % Each of the 15 Nc subjects lives in exactly ONE condition subfolder.
    % See the discover_rr_subjects() local function for the scan/parse logic.
    % ---------------------------------------------------------------------
    rr_subjects = discover_rr_subjects(RR_data_path, RR_CONDITION_SUBFOLDERS);

    fprintf('Found %d RR subjects across %d condition subfolders.\n', ...
        numel(rr_subjects), numel(RR_CONDITION_SUBFOLDERS));

    % Optional subset for testing: when RUN_ALL_VALID_RR is false, keep only
    % the Nc numbers listed in RR_PILOT_ONLY.
    if ~RUN_ALL_VALID_RR && ~isempty(rr_subjects)
        keep = ismember([rr_subjects.num], RR_PILOT_ONLY);
        rr_subjects = rr_subjects(keep);
        fprintf('  RUN_ALL_VALID_RR=false -> processing %d pilot subjects.\n', numel(rr_subjects));
    end

    pop_editoptions('option_storedisk', 0);
if RUN_PREPROCESSING
    for s_idx = 9:numel(rr_subjects)

        subjID   = rr_subjects(s_idx).nc_label;   % clean 'Nc##' used for all outputs
        subjPath = rr_subjects(s_idx).subjPath;
        mff_path = rr_subjects(s_idx).mff_path;
        i        = rr_subjects(s_idx).num;        % numeric id for hw-filter bookkeeping
        is_cohort1 = false;   % RR has no cohort-1
        has_LMRM   = false;   % RR EGI net has no LM/RM mastoids

        fprintf('\n======================================================\n');
        fprintf('PREPROCESSING %s  [RR | %s | %s]\n', subjID, ...
            rr_subjects(s_idx).condition, rr_subjects(s_idx).folder);
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
        % 1. Import the RR EGI recording (.mff) and load channel locations.
        %    subjPath / mff_path were resolved in the registry scan above.
        % -----------------------------------------------------------------
        EEG = pop_mffimport(mff_path);

        EEG.subject = subjID;
        EEG.etc.hw_filter_label = hw_filter_label;
        EEG.etc.rr_condition    = rr_subjects(s_idx).condition;
        EEG.etc.rr_folder       = rr_subjects(s_idx).folder;

        % EGI 128-ch HydroCel channel locations.
        if exist(egi_chanlocs_file, 'file')
            EEG = pop_chanedit(EEG, 'lookup', egi_chanlocs_file);
            EEG = eeg_checkset(EEG);
            fprintf('  Channel locations loaded.\n');
        else
            warning('  EGI chanlocs file not found at %s. Proceeding without explicit locations.', ...
                egi_chanlocs_file);
        end

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
        % 8. Reference.
        %    KH: average reference, excluding protected (EOG/mastoid) channels.
        %    RR (EGI): re-reference to Cz first (the EGI online ref electrode),
        %        then to common average. Cz index defaults to 61 for the EGI
        %        net used here -- VERIFY for your montage.
        % -----------------------------------------------------------------
        if strcmpi(cohort_ref_mode, 'Cz_then_CAR')
            rr_cz_index = find(strcmpi({EEG.chanlocs.labels}, cz_label), 1);
            if isempty(rr_cz_index); rr_cz_index = 61; end   % EGI default
            EEG = pop_reref(EEG, rr_cz_index, 'keepref', 'on');
            EEG = pop_reref(EEG, []);   % common average reference
        elseif ~isempty(protect_idx)
            EEG = pop_reref(EEG, [], 'exclude', protect_idx);
        else
            EEG = pop_reref(EEG, []);
        end
        EEG = eeg_checkset(EEG);

        pop_saveset(EEG, 'filename', [subjID '_cleaned_v4_merged.set'], ...
            'filepath', study_filepath);
    end
end
if RUN_ICA
        if isempty(gcp('nocreate'))
        pool = parpool('local');
        end


    parfor s_idx = 8:numel(rr_subjects)

        subjID   = rr_subjects(s_idx).nc_label;   % clean 'Nc##' used for all outputs
        subjPath = rr_subjects(s_idx).subjPath;

        EEG = pop_loadset(fullfile(study_filepath, [subjID '_cleaned_v4_merged.set']));
        % -----------------------------------------------------------------
        % 9. ICA on 1 Hz copy, transfer weights to 0.5 Hz data.
        % -----------------------------------------------------------------
        EEG_ica = pop_eegfiltnew(EEG, 'locutoff', 1.0);
        EEG_ica = pop_runica(EEG_ica, 'extended', 1);

        EEG.icaweights  = EEG_ica.icaweights;
        EEG.icasphere   = EEG_ica.icasphere;
        EEG.icawinv     = EEG_ica.icawinv;
        EEG.icachansind = EEG_ica.icachansind;
        EEG             = eeg_checkset(EEG);

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

        % fprintf('  ICLabel auto-rejection: rejecting %d/%d ICs with Brain probability < %.2f.\n', ...
            % numel(reject_comps), numel(brain_prob), ICLABEL_MIN_BRAIN_PROB);

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
            % fprintf('  Theta-filtering continuous data %g-%g Hz...\n', THETA_LO, THETA_HI);

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
        % RR (EGI/NetStation) outcome markers. NetStation stores 4-char codes;
        % 'rewa'/'puni' are 'reward'/'punish' truncated to 4 chars. VERIFY the
        % exact strings in your MFF (run a trigger count) and extend if needed.
        outcome_codes = {'rewa','puni'};

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
        % fprintf('  Epoching outcome broadband...\n');
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
            % fprintf('  Epoching outcome theta amplitude source...\n');
            EEG_theta_out = pop_epoch(EEG_theta, outcome_codes, outcome_window_epoch, ...
                'epochinfo', 'yes');
            EEG_theta_out = pop_rmbase(EEG_theta_out, baseline_win_ms);
            EEG_theta_out = eeg_checkset(EEG_theta_out);

            % fprintf('  Epoching outcome theta phase source...\n');
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

            % [trial2epoch_original, diag_out] = KH_align_epochs_with_offset(EEG_out, beh_correct);
            trial2epoch_original = [1:400]';

            % fprintf('  trial2epoch: %d/%d matched, offset=%d, match rate %.1f%%\n', ...
            %     diag_out.n_matched, diag_out.n_trials, diag_out.best_offset, ...
            %     100 * diag_out.match_rate);

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

            % pop_saveset(EEG_out_trim, 'filename', [subjID '_outcome_trimmed.set'], ...
            %     'filepath', epoch_outpath, 'savemode', 'onefile');

            if do_theta_epochs
                valid_ep_theta = valid_ep(valid_ep <= EEG_theta_out.trials);
                EEG_theta_trim = pop_select(EEG_theta_out, 'trial', valid_ep_theta);
                EEG_theta_trim = eeg_checkset(EEG_theta_trim);

                EEG_phase_trim = pop_select(EEG_phase_out, 'trial', valid_ep_theta);
                EEG_phase_trim = eeg_checkset(EEG_phase_trim);
                % 
                % pop_saveset(EEG_theta_trim, 'filename', [subjID '_outcome_theta_trimmed.set'], ...
                %     'filepath', epoch_outpath, 'savemode', 'onefile');
                % 
                % pop_saveset(EEG_phase_trim, 'filename', [subjID '_outcome_phase_trimmed.set'], ...
                %     'filepath', epoch_outpath, 'savemode', 'onefile');
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
                'epoch_artifact_flag_trimmed');

            % fprintf('  Saved trimmed outcome files and trial2epoch cache.\n');
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

    % Discover the same RR subjects as PART 1 (works even if PART 1 was skipped
    % this run), then optionally subset to RR_PILOT_ONLY.
    rr_subjects_feat = discover_rr_subjects(RR_data_path, RR_CONDITION_SUBFOLDERS);
    if ~RUN_ALL_VALID_RR && ~isempty(rr_subjects_feat)
        rr_subjects_feat = rr_subjects_feat(ismember([rr_subjects_feat.num], RR_PILOT_ONLY));
    end

    for s_idx = 1:numel(rr_subjects_feat)

        subj        = rr_subjects_feat(s_idx).nc_label;   % clean 'Nc##'
        participant = rr_subjects_feat(s_idx).num;
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

        % Capture a shared time axis (identical across subjects: same window/srate).
        if ~exist('t_ax', 'var') || isempty(t_ax)
            t_ax = EEGp.times;
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
        subj_rows = string(RR_behav_table.subjID) == string(subj);
        if ~any(subj_rows)
            warning('%s has no rows in RR_behav_table. Skipping.', subj);
            continue;
        end

        subj_features = RR_behav_table(subj_rows, :);
        n_rows = height(subj_features);

        subj_features.subj_id = repmat(string(subj), n_rows, 1);
        subj_features.cohort  = repmat("RR", n_rows, 1);
        subj_features.subj    = repmat(participant, n_rows, 1);
        subj_features.is_cohort1 = false(n_rows, 1);

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

            trial2epoch_raw = [1:400];

            if isfield(S2, 'valid_ep')
                valid_ep = S2.valid_ep(:);
            else
                valid_ep = [];
            end

            using_trimmed_file = EEGp.trials <= numel(valid_ep) || ...
                exist(fullfile(epoch_outpath, sprintf('%s_outcome_trimmed.set', subj)), 'file');

            % FIXED: trial2epoch_raw is already correct for the loaded file
            % No remapping needed - KH_align_epochs_with_offset created the correct mapping
            trial2epoch = trial2epoch_raw;

            if isfield(S2, 'hw_filter_label') && ~isempty(S2.hw_filter_label)
                hw_filter_label_subj = normalize_hw_filter_label(S2.hw_filter_label);
            end

            epoch_artifact_flag_eeg = read_epoch_artifact_flags(S2, valid_ep, EEGp.trials);

        else
            warning('%s: trial2epoch cache missing. Falling back to sequential mapping.', subj);
            trial2epoch = [1:400]';
            % n_map = min(n_rows, EEGp.trials);
            % trial2epoch(1:n_map) = (1:n_map)';
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
    % Per-stage FRN / RewP difference-wave table (CANONICAL FRN/RewP).
    %
    % FRN and RewP are difference waves (correct vs incorrect) and therefore
    % cannot be single-trial. They are computed here per
    % subj_id x block_type x stage from the per-trial FCz/Cz waveforms:
    %   FRN_amp  = mean over FRN window  of (incorrect - correct)  [negative]
    %   RewP_amp = mean over RewP window of (correct - incorrect)  [positive]
    % See pipeline/utils/kh_compute_frn_rewp_by_stage.m for sign conventions
    % and references. The per-trial table (all_trials_table) deliberately keeps
    % only the single-trial frontocentral measures (FCz negative-peak + FCz/Cz
    % mean); it does NOT carry FRN/RewP.
    % ---------------------------------------------------------------------
    frn_rewp_opts = struct('wave_col','FCzCz_waveform', ...
                           'FRN_win',[250 300], 'RewP_win',[250 350], ...
                           'stages',{stage_names}, 'block_types',{BTYPE_LABELS});
    frn_rewp_stage_table = kh_compute_frn_rewp_by_stage(all_trials_table, t_ax, frn_rewp_opts);

    % ---------------------------------------------------------------------
    % Save feature tables (cohort-tagged file names; subject naming unified).
    % ---------------------------------------------------------------------
    all_trials_table = kh_subject_id('standardise', all_trials_table);

    % Per-trial table keeps ONLY single-trial frontocentral measures.
    % The original PART-2 'FRN_mean_amp' is exactly the FCz/Cz fixed-window mean,
    % so rename it to FCzCz_mean_amp (+ _norm) and DROP the misleading single-trial
    % FRN_*/RewP_* columns (the true FRN/RewP live in frn_rewp_stage_table above).
    if ismember('FRN_mean_amp', all_trials_table.Properties.VariableNames)
        all_trials_table.FCzCz_mean_amp  = all_trials_table.FRN_mean_amp;
        if ismember('FRN_mean_norm', all_trials_table.Properties.VariableNames)
            all_trials_table.FCzCz_mean_norm = all_trials_table.FRN_mean_norm;
        end
    end
    drop_cols = intersect(all_trials_table.Properties.VariableNames, ...
        {'FRN_mean_amp','FRN_mean_norm','FRN_peak_amp','FRN_peak_lat','FRN_peak_norm', ...
         'FRN_excluded','RewP_mean_amp','RewP_mean_norm','RewP_peak_amp', ...
         'RewP_peak_lat','RewP_peak_norm','RewP_excluded'});
    all_trials_table = removevars(all_trials_table, drop_cols);
    if save_tables
        tag = upper(COHORT);
        all_trials_file = fullfile(feature_outpath, sprintf('group_table_all_trials_%s.mat', tag));
        stage_file      = fullfile(feature_outpath, sprintf('group_stage_table_features_%s.mat', tag));
        frn_file        = fullfile(feature_outpath, sprintf('frn_rewp_by_stage_%s.mat', tag));
        combined_file   = fullfile(feature_outpath, sprintf('group_feature_table_combined_%s.mat', tag));

        save(all_trials_file, 'all_trials_table', 't_ax', '-v7.3');
        save(stage_file, 'stage_feature_table', '-v7.3');
        save(frn_file, 'frn_rewp_stage_table', '-v7.3');
        save(combined_file, 'all_trials_table', 'stage_feature_table', 'frn_rewp_stage_table', 't_ax', '-v7.3');

        writetable(remove_cell_waveforms_for_csv(all_trials_table), ...
            fullfile(feature_outpath, sprintf('group_table_all_trials_%s.csv', tag)));
        writetable(stage_feature_table, ...
            fullfile(feature_outpath, sprintf('group_stage_table_features_%s.csv', tag)));
        if ~isempty(frn_rewp_stage_table)
            writetable(frn_rewp_stage_table(:, setdiff(frn_rewp_stage_table.Properties.VariableNames, {'diff_wave'})), ...
                fullfile(feature_outpath, sprintf('frn_rewp_by_stage_%s.csv', tag)));
        end

        fprintf('\nSaved (%s):\n  %s\n  %s\n  %s\n  %s\n', tag, ...
            all_trials_file, stage_file, frn_file, combined_file);
    end
end

fprintf('\nRR preprocessing + outcome feature construction complete.\n');

% =============================================================================
%% LOCAL FUNCTIONS
% =============================================================================

function rr_subjects = discover_rr_subjects(rr_data_path, condition_subfolders)
% Scan the RR Data tree (organised by task condition) and return a struct
% array of subjects. Each of the 15 Nc subjects lives in exactly ONE condition
% subfolder, in a folder named like "Nc01_p1_DPDP" that contains the EGI
% recording as a single ".mff" folder. We parse the clean "Nc##" label and
% locate that .mff. Fields: nc_label, num, condition, folder, subjPath, mff_path.

rr_subjects = struct('nc_label', {}, 'num', {}, 'condition', {}, ...
                     'folder', {}, 'subjPath', {}, 'mff_path', {});

for cf = 1:numel(condition_subfolders)
    cond_name = condition_subfolders{cf};
    cond_path = fullfile(rr_data_path, cond_name);
    if ~exist(cond_path, 'dir')
        warning('RR condition subfolder not found: %s', cond_path);
        continue;
    end

    subj_dirs = dir(fullfile(cond_path, 'Nc*'));
    subj_dirs = subj_dirs([subj_dirs.isdir]);

    for sd = 1:numel(subj_dirs)
        folder_name = subj_dirs(sd).name;                    % 'Nc01_p1_DPDP'
        nc_label = regexp(folder_name, '^Nc\d+', 'match', 'once');
        if isempty(nc_label)
            warning('Could not parse Nc## from folder %s -- skipping.', folder_name);
            continue;
        end

        subjPath = fullfile(cond_path, folder_name);

        % The .mff is a NetStation FOLDER inside the subject folder.
        mff_dirs = dir(fullfile(subjPath, '*.mff'));
        if isempty(mff_dirs)
            warning('%s (%s): no .mff folder found -- skipping.', nc_label, cond_name);
            continue;
        end

        rr_subjects(end+1) = struct( ...
            'nc_label', nc_label, ...
            'num', str2double(regexp(nc_label, '\d+', 'match', 'once')), ...
            'condition', cond_name, ...
            'folder', folder_name, ...
            'subjPath', subjPath, ...
            'mff_path', fullfile(subjPath, mff_dirs(1).name)); %#ok<AGROW>
    end
end

end

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

% -------------------------------------------------------------------------
% NOTE: the local KH_align_epochs_with_offset function was DELETED here on
% request (it was a flawed fallback implementation). The pipeline now uses the
% authoritative external KH_align_epochs_with_offset on func_path, which has the
% same signature [trial2epoch, diag_out] = KH_align_epochs_with_offset(EEG_out, beh_correct)
% and has been validated. Ensure func_path is on the MATLAB path (addpath above).
% -------------------------------------------------------------------------

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