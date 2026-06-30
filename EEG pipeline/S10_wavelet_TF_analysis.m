% =============================================================================
% S10_wavelet_TF_analysis.m  —  PIPELINE STEP 10 (optional)
%
% MORLET WAVELET TIME-FREQUENCY ANALYSIS  (BOTH cohorts: KH Ox## + RR Nc##)
%
% Computes event-related spectral perturbation (ERSP) and inter-trial
% phase coherence (ITPC) for outcome-locked epochs using complex Morlet
% wavelet convolution, then appends single-trial TF features to the SAME
% combined feature table used by S7 (group_feature_table_combined.mat).
%
% COHORTS
% -------
%   Both KH (Ox##) and RR (Nc##) subjects are processed in one loop. Each
%   cohort's epoched outcome .set files live under its own results tree
%   (KH_epoch_folder / RR_epoch_folder), and the two EEG nets use different
%   channel labels for the frontal midline / parietal ROIs:
%       KH (Curry/ANT)        : FCz ; Pz/P1/P2
%       RR (EGI 128 HydroCel) : E11 ; E62/E67/E72
%   (RR labels match S2_RR / S3_RR.) subj_id alone selects cohort, folder,
%   and channel set. Grand-average TF maps pool only subjects sharing the
%   reference epoch time grid; single-trial features are extracted for every
%   subject regardless.
%
% ALIGNED WITH S7
% ---------------
%   - Loads group_table from group_feature_table_combined.mat (from S4),
%     using the same workspace-aware load guard as S7 (reuses gt/group_table
%     if already in the workspace, otherwise loads from saved_tables_folder).
%   - Uses the canonical column names: subj_id (e.g. "Ox03"/"Nc07"),
%     block_type, stage (ordinal LN/LE/RN/RE), correct (0/1),
%     false_fb (logical), epoch.
%   - Within-subject z-scored features carry the _z suffix, as in S7.
%
% WHAT THIS DOES:
%   For every epoch, at every time point and every frequency between
%   MIN_FREQ and MAX_FREQ, the signal is convolved with a complex Morlet
%   wavelet. The result is a complex number whose:
%       - magnitude squared  = instantaneous power
%       - angle              = instantaneous phase
%   Power is then baseline-corrected to dB (10*log10(power/baseline)).
%   ITPC (also called PLF/ITC) is computed as |mean(exp(i*phase))| across
%   trials, ranging 0 (random phases) to 1 (perfectly locked phases).
%
% OUTPUTS (saved into saved_tables_folder, alongside the S7 table):
%   1. TF maps (ERSP and ITPC) per subject / stage / block type / outcome,
%      and grand-average maps — saved as PDFs to figure_output_folder.
%   2. group_feature_table_combined_wavelet.mat — group_table with appended
%      single-trial feature columns:
%         .theta_ersp    — mean dB power in theta band x FRN window
%         .alpha_ersp    — mean dB power in alpha band x P300 window
%         .beta_ersp     — mean dB power in beta band x post-feedback window
%         .theta_itpc    — mean ITPC in theta band x FRN window
%         .theta_ersp_z / .theta_itpc_z / ... — within-subject z-scores
%   3. grand_tf.mat — grand-average TF containers.
%   4. tf_perm_stats.mat — 2-D cluster-based permutation tests (correct vs
%      incorrect, D vs P) on the time-frequency plane.
%
% WHAT TO EXPECT:
%   - Theta (4-8 Hz) increase at FCz after incorrect outcomes, ~250-400 ms
%     (the oscillatory signature of the FRN).
%   - Alpha (8-12 Hz) desynchronisation at parietal channels for correct.
%   - Beta (15-25 Hz) suppression post-feedback at sensorimotor channels.
%   - ITPC peak in theta for incorrect: the FRN is phase-locked.
%
% LITERATURE:
%   Cohen (2014) Analyzing Neural Time Series Data, MIT Press.
%   Cohen & Donner (2013) NeuroImage 73:174-180. [wavelet definition]
%   Cavanagh & Frank (2014) Trends Cogn Sci 18:414-421. [frontal theta]
%   Trujillo & Allen (2007) Clin Neurophysiol 118:645-668. [FRN as theta]
%   Maris & Oostenveld (2007) J Neurosci Methods 164:177-190. [cluster test]
%
% REQUIRES: EEGLAB on path (pop_loadset) and the epoched outcome .set files
%           in Epoched_data_noisefiltering/. The cluster permutation test is
%           self-contained (no FieldTrip dependency) and uses bwconncomp
%           (Image Processing Toolbox) + tinv/tcdf (Statistics Toolbox).
% =============================================================================

close all;
addpath(genpath(fileparts(mfilename('fullpath'))));   % pipeline utils on path

% -------------------------------------------------------------------------
%% PATHS  (same base_path / saved_tables_folder layout as S7)
% -------------------------------------------------------------------------
base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch';
saved_tables_folder  = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis', 'Outcome_feature_tables_v4_merged');
KH_epoch_folder      = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis', 'Epoched_data_noisefiltering');
RR_epoch_folder      = fullfile(base_path, 'Salient mod switch RR', 'Results', 'EEG analysis', 'Epoched_data_noisefiltering');
figure_output_folder = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis', 'Figures', 'TF_analysis');
if ~exist(figure_output_folder,'dir'), mkdir(figure_output_folder); end

% Toolboxes (machine-specific; edit to match your environment).
eeglab_path = 'C:\Users\khatinova\OneDrive - Nexus365\Pre_2026_Folders\Documents\MATLAB\eeglab2025.1.0';
addpath(eeglab_path);    eeglab nogui;

% -------------------------------------------------------------------------
%% LOAD THE COMBINED FEATURE TABLE (same source/guard as S7)
% -------------------------------------------------------------------------
% S7b loads the table as 'gt'; this guard handles standalone use, chaining
% after S7/S7b, or having group_table already in the workspace.
if exist('gt', 'var') && ~exist('group_table', 'var')
    group_table = gt;
elseif ~exist('group_table', 'var')
    load(fullfile(saved_tables_folder, 'group_feature_table_combined.mat'), 'group_table');
end

% -------------------------------------------------------------------------
%% ENSURE CANONICAL TYPES / CODING (matches S7 + S7b conventions)
% -------------------------------------------------------------------------
group_table.subj_id = categorical(string(group_table.subj_id));

bt_str = string(group_table.block_type);
bt_str(bt_str == "V") = "P";                       % V -> P (probabilistic)
group_table.block_type = categorical(bt_str, {'D','P'});

group_table.stage = categorical(string(group_table.stage), ...
    {'LN','LE','RN','RE'}, 'Ordinal', true);

% correct -> numeric 0/1
if iscategorical(group_table.correct) || isstring(group_table.correct) || iscellstr(group_table.correct)
    cstr = lower(string(group_table.correct));
    cnum = nan(height(group_table),1);
    cnum(cstr == "correct"   | cstr == "1" | cstr == "true")  = 1;
    cnum(cstr == "incorrect" | cstr == "0" | cstr == "false") = 0;
    group_table.correct = cnum;
else
    group_table.correct = double(group_table.correct);
end
if all(ismember(unique(group_table.correct(~isnan(group_table.correct))), [1 2]))
    group_table.correct = group_table.correct - 1;  % 1/2 -> 0/1
end

% false_fb -> logical
if ismember('false_fb', group_table.Properties.VariableNames)
    group_table.false_fb = logical(group_table.false_fb);
else
    group_table.false_fb = false(height(group_table),1);
end

% epoch column is required to map table rows to .set trials
if ~ismember('epoch', group_table.Properties.VariableNames)
    error('S10: group_table has no ''epoch'' column — cannot map rows to epoched .set trials.');
end

fprintf('Loaded group_table: %d rows, %d subjects\n', ...
    height(group_table), numel(unique(group_table.subj_id)));

% -------------------------------------------------------------------------
%% WAVELET PARAMETERS
% -------------------------------------------------------------------------
MIN_FREQ  = 2;    % Hz
MAX_FREQ  = 30;   % Hz
N_FREQS   = 40;   % number of frequencies (log-spaced)
MIN_CYCL  = 3;    % number of cycles at MIN_FREQ
MAX_CYCL  = 7;    % number of cycles at MAX_FREQ (scales linearly with freq)

%   Frequency vector (log-spaced for equal relative bandwidth)
frex = logspace(log10(MIN_FREQ), log10(MAX_FREQ), N_FREQS);

%   Number of cycles scales linearly from MIN_CYCL to MAX_CYCL.
%   This gives better time resolution at low frequencies (where you need it)
%   and better frequency resolution at high frequencies.
n_cycles = logspace(log10(MIN_CYCL), log10(MAX_CYCL), N_FREQS);

%   Bands of interest (for single-trial feature extraction)
theta_band = frex >= 4  & frex <= 8;
alpha_band = frex >= 8  & frex <= 12;
beta_band  = frex >= 15 & frex <= 25;

%   Time windows of interest (ms)
BL_WIN    = [-200  0];   % baseline window for dB correction
FRN_WIN   = [200 400];   % FRN/theta feature window (Cavanagh & Frank 2014)
P300_WIN  = [300 600];   % P300/alpha window
BETA_WIN  = [400 700];   % beta rebound window
PLOT_WIN  = [-200 800];  % display window

%   Minimum trials per condition to compute a meaningful TF map
MIN_TRIALS = 10;

%   Stage and block info
stage_names  = {'LN','LE','RN','RE'};
BTYPE_LABELS = {'D','P'};
STAGE_COLORS = [0.12 0.62 0.47; 0.85 0.65 0.00; 0.80 0.27 0.13; 0.40 0.25 0.65];

% -------------------------------------------------------------------------
%% COHORT-AWARE PARTICIPANTS + CHANNEL LABELS
%
% Both cohorts are processed. They use DIFFERENT EEG nets, so the frontal
% midline / parietal channel labels differ:
%   KH (Curry/ANT)        : frontal midline = FCz, parietal = Pz/P1/P2
%   RR (EGI 128 HydroCel) : frontal midline = E11, parietal = E62/E67/E72
% (RR labels match S2_RR_preprocess_epoch_eeg.m / S3_RR_extract_eeg_features.m.)
%
% Each cohort's epoched outcome .set files live under its own results tree
% (KH_epoch_folder / RR_epoch_folder, defined in PATHS above). The combined
% feature table already carries both Ox## and Nc## rows with cohort-relative
% epoch indices, so subj_id alone determines where to look.
% -------------------------------------------------------------------------
KH_PARTICIPANTS = [3:12, 14:23, 27:28];   % Ox## (KH)
RR_PARTICIPANTS = 1:15;                    % Nc## (RR)

KH_FCZ_LABELS = {'FCz'};            KH_PZ_LABELS = {'Pz','P1','P2'};
RR_FCZ_LABELS = {'E11'};            RR_PZ_LABELS = {'E62','E67','E72'};

% Combined, cohort-prefixed subject list driving the participant loop.
ALL_SUBJECTS = [compose("Ox%02d", KH_PARTICIPANTS(:)); ...
                compose("Nc%02d", RR_PARTICIPANTS(:))];

% -------------------------------------------------------------------------
%% INITIALISE GRAND-AVERAGE TF CONTAINERS
%
% Each leaf stores: .ersp (n_subj x n_freq x n_time) and .subj (string ids)
% for the FCz channel. .subj holds the canonical subj_id string ("Ox03" /
% "Nc07") — NOT the numeric participant number — so KH and RR subjects with
% the same number never collide in the cross-subject grand averages / stats.
% -------------------------------------------------------------------------
for s = 1:4
    for bt = 1:2
        bt_s = BTYPE_LABELS{bt};
        for oc = {'correct','incorrect'}
            grand_tf.FCz.(stage_names{s}).(bt_s).(oc{1}).ersp = [];
            grand_tf.FCz.(stage_names{s}).(bt_s).(oc{1}).itpc = [];
            grand_tf.FCz.(stage_names{s}).(bt_s).(oc{1}).subj = strings(0,1);
        end
    end
end

% -------------------------------------------------------------------------
%% ADD WAVELET FEATURE COLUMNS TO GROUP TABLE
% -------------------------------------------------------------------------
n_rows = height(group_table);
group_table.theta_ersp = nan(n_rows,1);
group_table.alpha_ersp = nan(n_rows,1);
group_table.beta_ersp  = nan(n_rows,1);
group_table.theta_itpc = nan(n_rows,1);

t_global = [];   % time axis, filled from first loaded subject

% =========================================================================
%% PARTICIPANT LOOP  (KH Ox## + RR Nc##)
% =========================================================================
for sidx = 1:numel(ALL_SUBJECTS)

    subj = char(ALL_SUBJECTS(sidx));            % canonical subj_id ("Ox03"/"Nc07")
    if startsWith(subj, 'Ox')
        cohort = 'KH'; epoch_folder = KH_epoch_folder;
        fcz_labels = KH_FCZ_LABELS; pz_labels = KH_PZ_LABELS;
    else
        cohort = 'RR'; epoch_folder = RR_epoch_folder;
        fcz_labels = RR_FCZ_LABELS; pz_labels = RR_PZ_LABELS;
    end
    participant = str2double(regexp(subj, '\d+', 'match', 'once'));
    subj_rows   = group_table.subj_id == subj;
    fprintf('\n=== %s (%s) ===\n', subj, cohort);

    if ~any(subj_rows)
        fprintf('  no rows in group_table for %s, skipping.\n', subj); continue
    end

    % Locate outcome .set: prefer the trimmed file, fall back to untrimmed.
    fname = '';
    if     exist(fullfile(epoch_folder, sprintf('%s_outcome_trimmed.set', subj)),'file')
        fname = sprintf('%s_outcome_trimmed.set', subj);
    elseif exist(fullfile(epoch_folder, sprintf('%s_outcome.set', subj)),'file')
        fname = sprintf('%s_outcome.set', subj);
    end
    if isempty(fname)
        warning('%s: no outcome .set found in %s, skipping.', subj, epoch_folder); continue
    end
    EEGp = pop_loadset(fname, epoch_folder);
    fprintf('  Loaded %d epochs, %d channels, Fs=%d Hz\n', ...
        EEGp.trials, EEGp.nbchan, EEGp.srate);

    if isempty(t_global), t_global = EEGp.times; end

    % Subjects whose epoch time grid matches the reference can contribute to
    % the cross-subject grand-average TF maps; mismatched grids still get
    % single-trial features (STEP 1b) but are excluded from grand averages.
    grid_ok = (numel(EEGp.times) == numel(t_global));
    if ~grid_ok
        warning(['%s: epoch time-axis length %d ~= reference %d; excluded from ' ...
            'grand-average TF maps (single-trial features still computed).'], ...
            subj, numel(EEGp.times), numel(t_global));
    end

    % Time masks
    bl_mask   = EEGp.times >= BL_WIN(1)   & EEGp.times <= BL_WIN(2);
    frn_mask  = EEGp.times >= FRN_WIN(1)  & EEGp.times <= FRN_WIN(2);
    p300_mask = EEGp.times >= P300_WIN(1) & EEGp.times <= P300_WIN(2);
    beta_mask = EEGp.times >= BETA_WIN(1) & EEGp.times <= BETA_WIN(2);
    plot_mask = EEGp.times >= PLOT_WIN(1) & EEGp.times <= PLOT_WIN(2);

    % Channel indices (cohort-specific labels)
    fcz_idx = safe_chan(EEGp, fcz_labels);
    pz_idx  = safe_chan(EEGp, pz_labels);

    % ------------------------------------------------------------------
    % STEP 1: COMPUTE MORLET WAVELET TRANSFORM FOR ALL EPOCHS
    %
    % Strategy: precompute the full TF decomposition once for the ROI
    % channels, then extract condition-specific epochs afterwards.
    % This avoids redundant computation across conditions.
    %
    % tf_power : n_freq x n_time x n_epoch  (raw power, before baseline)
    % tf_phase : n_freq x n_time x n_epoch  (phase angles, radians)
    % ------------------------------------------------------------------
    n_times  = size(EEGp.data, 2);
    n_epochs = EEGp.trials;

    if ~isempty(fcz_idx)
        fprintf('  Computing wavelet TF at FCz (%d freqs x %d times x %d epochs)...\n', ...
            N_FREQS, n_times, n_epochs);

        % Average channel data across ROI (for multi-channel idx)
        if isscalar(fcz_idx)
            chan_data = squeeze(double(EEGp.data(fcz_idx,:,:)));
        else
            chan_data = squeeze(mean(double(EEGp.data(fcz_idx,:,:)), 1, 'omitnan'));
        end
        % chan_data: n_times x n_epochs

        [tf_power_fcz, tf_phase_fcz] = morlet_wavelet_tf( ...
            chan_data, EEGp.srate, frex, n_cycles);
        % tf_power_fcz: N_FREQS x n_times x n_epochs

        fprintf('  Done.\n');
    else
        tf_power_fcz = [];
        tf_phase_fcz = [];
        warning('%s: FCz not found.', subj);
    end

    % Parietal channel TF (for P300/alpha)
    if ~isempty(pz_idx)
        if isscalar(pz_idx)
            chan_data_pz = squeeze(double(EEGp.data(pz_idx,:,:)));
        else
            chan_data_pz = squeeze(mean(double(EEGp.data(pz_idx,:,:)), 1, 'omitnan'));
        end
        [tf_power_pz, ~] = morlet_wavelet_tf(chan_data_pz, EEGp.srate, frex, n_cycles);
    else
        tf_power_pz = [];
    end

    % ------------------------------------------------------------------
    % STEP 1b: SINGLE-TRIAL ERSP POWER FEATURES FOR EVERY EPOCH
    %
    % theta/alpha/beta ERSP are per-trial quantities: a trial's band power
    % does not depend on how many other trials share its condition cell.
    % They are therefore computed here for ALL epochs and are deliberately
    % NOT gated by MIN_TRIALS (that gate only governs the grand-average TF
    % maps and the cell-level ITPC below, which need enough trials to be
    % stable). Decoupling means sparse cells — e.g. the few incorrect trials
    % in deterministic blocks — still contribute, which is exactly what the
    % downstream difference figure / LME n depends on.
    %
    % Single-trial ERSP uses the subject-level mean baseline (Cohen 2014,
    % Ch.18) so trial-wise variance is preserved.
    % ------------------------------------------------------------------
    if ~isempty(tf_power_fcz)
        bl_mean_subj = mean(mean(tf_power_fcz(:,bl_mask,:), 2, 'omitnan'), 3, 'omitnan'); % N_FREQS x 1
        for ep = 1:n_epochs
            rows_ep = find(subj_rows & group_table.epoch == ep);
            if isempty(rows_ep), continue; end
            trial_ersp = 10 * log10(double(tf_power_fcz(:,:,ep)) ./ bl_mean_subj);
            group_table.theta_ersp(rows_ep) = mean(trial_ersp(theta_band, frn_mask),  'all', 'omitnan');
            group_table.alpha_ersp(rows_ep) = mean(trial_ersp(alpha_band, p300_mask), 'all', 'omitnan');
            group_table.beta_ersp(rows_ep)  = mean(trial_ersp(beta_band,  beta_mask), 'all', 'omitnan');
        end
    end

    % ------------------------------------------------------------------
    % STEP 2: ACCUMULATE GRAND AVERAGES + CELL-LEVEL ITPC
    % ------------------------------------------------------------------
    has_P = any(group_table.block_type(subj_rows) == 'P');

    for s_i = 1:4
        for bt_i = 1:2
            bt = BTYPE_LABELS{bt_i};
            if strcmp(bt,'P') && ~has_P; continue; end

            base_mask_gt = subj_rows & ...
                group_table.block_type == bt & ...
                group_table.stage      == stage_names{s_i};

            for oc_i = 1:2
                if oc_i == 1
                    oc_str = 'correct';   oc_val = 1;
                else
                    oc_str = 'incorrect'; oc_val = 0;
                end

                cond_mask_gt = base_mask_gt & ...
                    group_table.correct  == oc_val & ...
                    ~group_table.false_fb;

                ep_vec = group_table.epoch(cond_mask_gt);
                ep_vec = ep_vec(~isnan(ep_vec) & ep_vec>=1 & ep_vec<=n_epochs);

                if numel(ep_vec) < MIN_TRIALS
                    fprintf('  SKIP %s %s %s — only %d trials\n', ...
                        stage_names{s_i}, bt, oc_str, numel(ep_vec));
                    continue
                end

                %% --- ERSP (dB baseline-corrected power) ---
                if ~isempty(tf_power_fcz)
                    cond_power = tf_power_fcz(:,:,ep_vec);
                    % N_FREQS x n_times x n_trials

                    ersp = compute_ersp(cond_power, bl_mask);
                    % ersp: N_FREQS x n_times (mean across trials, baseline-corrected dB)

                    %% --- ITPC (inter-trial phase coherence) ---
                    cond_phase = tf_phase_fcz(:,:,ep_vec);
                    itpc = abs(mean(exp(1i * cond_phase), 3, 'omitnan'));
                    % itpc: N_FREQS x n_times

                    %% --- Accumulate grand average (common time grid only) ---
                    % .subj stores the canonical subj_id STRING so KH/RR
                    % subjects with the same number never collide.
                    if grid_ok
                        grand_tf.FCz.(stage_names{s_i}).(bt).(oc_str).ersp(end+1,:,:) = ersp;
                        grand_tf.FCz.(stage_names{s_i}).(bt).(oc_str).itpc(end+1,:,:) = itpc;
                        grand_tf.FCz.(stage_names{s_i}).(bt).(oc_str).subj(end+1)     = string(subj);
                    end

                    %% --- Cell-level ITPC feature ---
                    % ITPC is computed ACROSS the trials of this cell, so it
                    % is genuinely a cell-level (not single-trial) measure and
                    % stays gated by MIN_TRIALS. Every trial in the cell is
                    % tagged with the cell's theta ITPC. (The per-trial ERSP
                    % power features were computed once per epoch in STEP 1b,
                    % decoupled from this gate.)
                    itpc_rows = find(cond_mask_gt & ~isnan(group_table.epoch) & ...
                        group_table.epoch >= 1 & group_table.epoch <= n_epochs);
                    group_table.theta_itpc(itpc_rows) = mean(itpc(theta_band, frn_mask), 'all', 'omitnan');

                end % ~isempty fcz
            end % outcome loop
        end % block type loop
    end % stage loop

    % ------------------------------------------------------------------
    % STEP 3: PER-SUBJECT TF FIGURES
    %
    % One figure per subject: 4 stages x 2 outcomes (correct/incorrect),
    % correct and incorrect overlaid on same colorscale for comparison.
    % Bottom row: ITPC.
    % ------------------------------------------------------------------
    if ~isempty(tf_power_fcz) && grid_ok
        plot_subject_tf(subj, participant, has_P, group_table, subj_rows, ...
            tf_power_fcz, tf_phase_fcz, tf_power_pz, ...
            frex, t_global, bl_mask, plot_mask, ...
            stage_names, BTYPE_LABELS, STAGE_COLORS, ...
            theta_band, FRN_WIN, N_FREQS, MIN_TRIALS, ...
            figure_output_folder, n_epochs);
    end

    clear EEGp tf_power_fcz tf_phase_fcz tf_power_pz chan_data chan_data_pz

end % participant loop

% =========================================================================
%% WITHIN-SUBJECT Z-SCORE OF WAVELET FEATURES  (same _z convention as S7)
% =========================================================================
fprintf('\nZ-scoring wavelet features within subject...\n');
tf_features = {'theta_ersp','alpha_ersp','beta_ersp','theta_itpc'};
subj_list   = unique(group_table.subj_id);

for f = 1:numel(tf_features)
    fn   = tf_features{f};
    fn_z = [fn '_z'];
    group_table.(fn_z) = nan(n_rows, 1);
    for si = 1:numel(subj_list)
        mask = group_table.subj_id == subj_list(si);
        vals = group_table.(fn)(mask);
        mn   = mean(vals, 'omitnan');
        sd   = std(vals,  'omitnan');
        if sd > 0
            group_table.(fn_z)(mask) = (vals - mn) / sd;
        end
    end
end

save(fullfile(saved_tables_folder,'group_feature_table_combined_wavelet.mat'), 'group_table');
fprintf('Saved group_feature_table_combined_wavelet.mat\n');

save(fullfile(saved_tables_folder,'grand_tf.mat'), 'grand_tf', 't_global', 'frex');
fprintf('Saved grand_tf.mat\n');

% =========================================================================
%% GRAND-AVERAGE TF FIGURES
% =========================================================================
fprintf('\nPlotting grand-average TF maps...\n');

t_plot = t_global(t_global >= PLOT_WIN(1) & t_global <= PLOT_WIN(2));
clim_ersp = [-3 3];   % dB colour scale
clim_itpc = [0 0.5];

for s_i = 1:4
    fig = figure('Position',[50 50 1200 700]);
    sgtitle(sprintf('Grand average TF — FCz — %s', stage_names{s_i}), 'FontSize',11);

    for bt_i = 1:2
        bt = BTYPE_LABELS{bt_i};

        for oc_i = 1:2
            oc_str = {'correct','incorrect'};
            oc     = oc_str{oc_i};

            % ERSP subplot
            ax_ersp = subplot(2, 4, (bt_i-1)*4 + oc_i);
            g = grand_tf.FCz.(stage_names{s_i}).(bt).(oc);
            if ~isempty(g.ersp)
                plot_mask_g = t_global >= PLOT_WIN(1) & t_global <= PLOT_WIN(2);
                grand_ersp  = squeeze(mean(g.ersp(:,:,plot_mask_g), 1, 'omitnan'));
                % n_subj x n_freq x n_time -> freq x time after mean
                imagesc(ax_ersp, t_plot, frex, grand_ersp);
                axis(ax_ersp,'xy');
                set(ax_ersp,'YTick',[4 8 12 20 30],'YScale','log');
                clim(ax_ersp, clim_ersp);
                colormap(ax_ersp, tf_colormap());
                colorbar(ax_ersp);
                xline(ax_ersp, 0, 'w:', 'LineWidth',1);
                yline(ax_ersp, 4, 'w--', 'LineWidth',0.8);
                yline(ax_ersp, 8, 'w--', 'LineWidth',0.8);
                title(ax_ersp, sprintf('%s | %s | %s | n=%d', ...
                    stage_names{s_i}, bt, oc, size(g.ersp,1)), 'FontSize',8);
            else
                title(ax_ersp, sprintf('%s | %s | %s | no data', stage_names{s_i}, bt, oc));
            end
            xlabel(ax_ersp,'Time (ms)');
            ylabel(ax_ersp,'Frequency (Hz)');

            % ITPC subplot
            ax_itpc = subplot(2, 4, (bt_i-1)*4 + oc_i + 2);
            if ~isempty(g.itpc)
                grand_itpc = squeeze(mean(g.itpc(:,:,t_global>=PLOT_WIN(1)&t_global<=PLOT_WIN(2)), 1, 'omitnan'));
                imagesc(ax_itpc, t_plot, frex, grand_itpc);
                axis(ax_itpc,'xy');
                set(ax_itpc,'YTick',[4 8 12 20 30],'YScale','log');
                clim(ax_itpc, clim_itpc);
                colormap(ax_itpc, parula);
                colorbar(ax_itpc);
                xline(ax_itpc,0,'w:','LineWidth',1);
                yline(ax_itpc,4,'w--','LineWidth',0.8);
                yline(ax_itpc,8,'w--','LineWidth',0.8);
                title(ax_itpc, sprintf('ITPC — %s | %s | %s', stage_names{s_i}, bt, oc), 'FontSize',8);
            end
            xlabel(ax_itpc,'Time (ms)');
            ylabel(ax_itpc,'Frequency (Hz)');
        end
    end

    saveas(fig, fullfile(figure_output_folder, ...
        sprintf('S10_TF_grand_%s.pdf', stage_names{s_i})));
    % close(fig);
end

% =========================================================================
%% GRAND-AVERAGE: CORRECT vs INCORRECT DIFFERENCE TF MAPS
%
% The difference map (incorrect minus correct) isolates the error-related
% theta signal. This is analogous to the difference-wave approach used
% to isolate the FRN in ERPs, but in time-frequency space.
% Expected result: positive cluster in theta band, 200-500 ms, at FCz.
% =========================================================================
fprintf('Plotting correct vs incorrect difference TF maps...\n');

fig_diff = figure('Position',[50 50 1400 600]);
sgtitle('Grand average: Incorrect − Correct ERSP (FCz)', 'FontSize',11);

for s_i = 1:4
    for bt_i = 1:2
        bt  = BTYPE_LABELS{bt_i};
        ax  = subplot(2, 4, (bt_i-1)*4 + s_i);

        gc = grand_tf.FCz.(stage_names{s_i}).(bt).correct;
        gi = grand_tf.FCz.(stage_names{s_i}).(bt).incorrect;

        % Use only subjects present in both conditions
        shared = intersect(gc.subj, gi.subj);
        if numel(shared) < 3
            title(ax, sprintf('%s %s — insufficient n', stage_names{s_i}, bt));
            continue
        end

        [~,ia] = ismember(shared, gc.subj);
        [~,ib] = ismember(shared, gi.subj);

        plot_t_mask = t_global >= PLOT_WIN(1) & t_global <= PLOT_WIN(2);

        ersp_c = gc.ersp(ia,:,plot_t_mask);   % n x freq x time
        ersp_i = gi.ersp(ib,:,plot_t_mask);

        diff_ersp = mean(ersp_i, 1, 'omitnan') - mean(ersp_c, 1, 'omitnan');
        diff_ersp = squeeze(diff_ersp);   % freq x time

        imagesc(ax, t_plot, frex, diff_ersp);
        axis(ax,'xy');
        set(ax,'YTick',[4 8 12 20 30],'YScale','log');
        clim(ax, [-2 2]);
        colormap(ax, tf_colormap());
        colorbar(ax);
        xline(ax,0,'k:','LineWidth',1);
        yline(ax,4,'k--','LineWidth',0.8);
        yline(ax,8,'k--','LineWidth',0.8);
        title(ax, sprintf('%s | %s | n=%d', stage_names{s_i}, bt, numel(shared)), 'FontSize',9);
        xlabel(ax,'Time (ms)');
        ylabel(ax,'Frequency (Hz)');
    end
end

saveas(fig_diff, fullfile(figure_output_folder,'S10_TF_diff_IncVsCor.pdf'));
close(fig_diff);

% =========================================================================
%% TF CLUSTER PERMUTATION TESTS
%
% Two-dimensional cluster permutation test on the freq x time plane,
% implemented natively (Maris & Oostenveld 2007) so it does NOT require
% FieldTrip / ft_freqstatistics (which needs a spatial channel dimension).
%
%   1. Paired t-statistic at every (freq, time) point.
%   2. Threshold at |t| > t_crit (alpha = 0.025 per tail, df = n-1).
%   3. 8-connected components in the thresholded map = clusters.
%   4. Cluster statistic = sum of t-values within each cluster.
%   5. Sign-flip permutation (paired design) builds the null distribution.
%   6. p = proportion of null max-cluster statistics >= observed.
%
% Output tf_perm_stats.(contrast).{stat,mask,posclusters,negclusters,...}
% mirrors the fields plot_tf_cluster_result expects.
% =========================================================================
MIN_N_SUBJ   = 8;
N_PERM       = 10000;
ALPHA_CLUST  = 0.025;  % two-tailed, so each tail uses 0.025

tf_perm_stats = struct();

plot_t_mask = t_global >= PLOT_WIN(1) & t_global <= PLOT_WIN(2);
t_plot      = t_global(plot_t_mask);

fprintf('\n\n=== TF CLUSTER PERMUTATION TESTS ===\n');

% ── Contrast A: Correct vs Incorrect at each stage x block_type ─────────────
for s_i = 1:4
    for bt_i = 1:2
        bt = BTYPE_LABELS{bt_i};

        gc = grand_tf.FCz.(stage_names{s_i}).(bt).correct;
        gi = grand_tf.FCz.(stage_names{s_i}).(bt).incorrect;

        shared = intersect(gc.subj, gi.subj);
        if numel(shared) < MIN_N_SUBJ
            fprintf('  SKIP %s %s — only %d shared subjects\n', ...
                stage_names{s_i}, bt, numel(shared));
            continue
        end

        [~, ia] = ismember(shared, gc.subj);
        [~, ib] = ismember(shared, gi.subj);

        % diff_maps: n_subj x n_freq x n_time  (incorrect − correct)
        diff_maps = gi.ersp(ib,:,plot_t_mask) - gc.ersp(ia,:,plot_t_mask);

        contrast_name = sprintf('%s_%s_CorVsInc_TF', stage_names{s_i}, bt);
        fprintf('  Running: %s (n=%d)...\n', contrast_name, numel(shared));

        stat_struct = run_tf_cluster_permtest(diff_maps, N_PERM, ALPHA_CLUST);

        tf_perm_stats.(contrast_name)   = stat_struct;
        tf_perm_stats.(contrast_name).n = numel(shared);

        % Report and plot significant clusters
        report_tf_clusters(stat_struct, contrast_name);
        if any(stat_struct.mask(:))
            plot_tf_cluster_result(stat_struct, frex, t_plot, contrast_name, ...
                figure_output_folder, PLOT_WIN);
        end
    end
end

% ── Contrast B: D vs P blocks ────────────────────────────────────────────────
for s_i = 1:4
    d_D_c  = grand_tf.FCz.(stage_names{s_i}).D.correct;
    d_P_c  = grand_tf.FCz.(stage_names{s_i}).P.correct;
    d_D_ic = grand_tf.FCz.(stage_names{s_i}).D.incorrect;
    d_P_ic = grand_tf.FCz.(stage_names{s_i}).P.incorrect;

    shared_all = intersect(intersect(d_D_c.subj, d_D_ic.subj), ...
                           intersect(d_P_c.subj, d_P_ic.subj));
    if numel(shared_all) < MIN_N_SUBJ
        shared_all = intersect(d_D_c.subj, d_P_c.subj);
    end
    if numel(shared_all) < MIN_N_SUBJ
        fprintf('  SKIP D vs P %s — only %d subjects\n', ...
            stage_names{s_i}, numel(shared_all));
        continue
    end

    n_sh  = numel(shared_all);
    n_tf  = sum(plot_t_mask);
    n_frq = numel(frex);
    dat_D = zeros(n_sh, n_frq, n_tf);
    dat_P = zeros(n_sh, n_frq, n_tf);

    for si = 1:n_sh
        s_id = shared_all(si);
        % Pool correct + incorrect per block type per subject
        row_Dc = find(d_D_c.subj  == s_id, 1);
        row_Di = find(d_D_ic.subj == s_id, 1);
        row_Pc = find(d_P_c.subj  == s_id, 1);
        row_Pi = find(d_P_ic.subj == s_id, 1);

        n_D = 0;
        if ~isempty(row_Dc), dat_D(si,:,:) = dat_D(si,:,:) + d_D_c.ersp(row_Dc,:,plot_t_mask);  n_D=n_D+1; end
        if ~isempty(row_Di), dat_D(si,:,:) = dat_D(si,:,:) + d_D_ic.ersp(row_Di,:,plot_t_mask); n_D=n_D+1; end
        if n_D > 1, dat_D(si,:,:) = dat_D(si,:,:) / n_D; end

        n_P = 0;
        if ~isempty(row_Pc), dat_P(si,:,:) = dat_P(si,:,:) + d_P_c.ersp(row_Pc,:,plot_t_mask);  n_P=n_P+1; end
        if ~isempty(row_Pi), dat_P(si,:,:) = dat_P(si,:,:) + d_P_ic.ersp(row_Pi,:,plot_t_mask); n_P=n_P+1; end
        if n_P > 1, dat_P(si,:,:) = dat_P(si,:,:) / n_P; end
    end

    diff_maps = dat_D - dat_P;   % positive = D > P

    contrast_name = sprintf('%s_DvsP_TF', stage_names{s_i});
    fprintf('  Running: %s (n=%d)...\n', contrast_name, n_sh);

    stat_struct = run_tf_cluster_permtest(diff_maps, N_PERM, ALPHA_CLUST);
    tf_perm_stats.(contrast_name)   = stat_struct;
    tf_perm_stats.(contrast_name).n = n_sh;

    report_tf_clusters(stat_struct, contrast_name);
    if any(stat_struct.mask(:))
        plot_tf_cluster_result(stat_struct, frex, t_plot, contrast_name, ...
            figure_output_folder, PLOT_WIN);
    end
end

% ── Save ─────────────────────────────────────────────────────────────────────
save(fullfile(saved_tables_folder,'tf_perm_stats.mat'), 'tf_perm_stats');
fprintf('\nTF permutation stats saved to tf_perm_stats.mat\n');

% =========================================================================
%% THETA ERSP OVER STAGES: LINE PLOT SUMMARY
%
% After extracting single-trial features, plot mean theta ERSP (±SEM)
% across stages for correct vs incorrect, separately for D and P blocks.
% This gives a compact summary of the learning trajectory in theta power.
% =========================================================================
fprintf('\nPlotting theta ERSP stage summary...\n');

fig_summary = figure('Position',[50 50 900 400]);
sgtitle('Mean theta ERSP (4-8 Hz, 200-400 ms) by stage', 'FontSize',11);

for bt_i = 1:2
    bt  = BTYPE_LABELS{bt_i};
    ax  = subplot(1,2,bt_i);
    hold(ax,'on');
    title(ax, sprintf('Block type: %s', bt));
    xlabel(ax,'Stage'); ylabel(ax,'Theta ERSP (dB)');
    xline(ax,2.5,'k:','HandleVisibility','off');
    clrs = [0.1 0.6 0.1; 0.7 0.1 0.1];

    for oc_i = 1:2
        oc_val  = 2 - oc_i;   % 1=correct, 0=incorrect
        oc_name = {'Correct','Incorrect'};
        means   = nan(1,4); sems = nan(1,4);
        for s_i = 1:4
            m = group_table.block_type==bt & ...
                group_table.stage==stage_names{s_i} & ...
                group_table.correct==oc_val & ...
                ~group_table.false_fb & ...
                ~isnan(group_table.theta_ersp);
            vals = group_table.theta_ersp(m);
            if ~isempty(vals)
                means(s_i) = mean(vals,'omitnan');
                sems(s_i)  = std(vals,'omitnan') / sqrt(sum(~isnan(vals)));
            end
        end
        errorbar(ax, 1:4, means, sems, 'o-', ...
            'Color', clrs(oc_i,:), 'LineWidth',2, ...
            'MarkerFaceColor', clrs(oc_i,:), ...
            'MarkerSize',7, 'DisplayName', oc_name{oc_i});
    end
    set(ax,'XTick',1:4,'XTickLabel',stage_names);
    legend(ax,'Location','best','Box','off');
    yline(ax,0,'k:');
end

saveas(fig_summary, fullfile(figure_output_folder,'S10_TF_theta_stage_summary.pdf'));
close(fig_summary);

% =========================================================================
%% THETA POWER DIFFERENCE FIGURE
%
% Four requested contrasts on theta ERSP (4-8 Hz, 200-400 ms, FCz):
%   1. D blocks : Incorrect − Correct
%   2. P blocks : Incorrect − Correct
%   3. Correct  : D − P
%   4. Incorrect: D − P
%
% Each panel shows the two conditions being contrasted (group mean bar +
% per-subject paired points/lines), with the mean difference (Delta, dB)
% and a paired t-test annotated. The difference IS the paired contrast, so
% the connecting lines and Delta text make it directly readable.
% Trials are restricted to true-feedback trials (~false_fb) and collapsed
% across stages (per-subject means).
% =========================================================================
fprintf('\nPlotting theta-power difference figure...\n');

% Minimum trials per subject per cell to contribute a (stable) subject mean.
% Lower this toward 1 to maximise n at the cost of noisier subject means.
MIN_TRIALS_DIFF = 5;

subjs = unique(group_table.subj_id);
nS    = numel(subjs);
[thD_c, thD_i, thP_c, thP_i] = deal(nan(nS,1));

for si = 1:nS
    sm = group_table.subj_id == subjs(si) & ...
         ~group_table.false_fb & ~isnan(group_table.theta_ersp);
    thD_c(si) = mean_if_enough(group_table.theta_ersp(sm & group_table.block_type=='D' & group_table.correct==1), MIN_TRIALS_DIFF);
    thD_i(si) = mean_if_enough(group_table.theta_ersp(sm & group_table.block_type=='D' & group_table.correct==0), MIN_TRIALS_DIFF);
    thP_c(si) = mean_if_enough(group_table.theta_ersp(sm & group_table.block_type=='P' & group_table.correct==1), MIN_TRIALS_DIFF);
    thP_i(si) = mean_if_enough(group_table.theta_ersp(sm & group_table.block_type=='P' & group_table.correct==0), MIN_TRIALS_DIFF);
end

% Transparency: report how many subjects contribute (and why this can be
% fewer than the number of EEG datasets available).
fprintf('  Difference-figure n (subjects with >= %d true-FB trials in the cell):\n', MIN_TRIALS_DIFF);
fprintf('    D: Correct=%d, Incorrect=%d | P: Correct=%d, Incorrect=%d\n', ...
    sum(~isnan(thD_c)), sum(~isnan(thD_i)), sum(~isnan(thP_c)), sum(~isnan(thP_i)));
fprintf('    NOTE: wavelet loop attempted %d KH + %d RR participants;\n', ...
    numel(KH_PARTICIPANTS), numel(RR_PARTICIPANTS));
fprintf('          combined table has %d subjects. Each contrast is PAIRED, so n is\n', nS);
fprintf('          the subjects with data in BOTH of its conditions; sparse incorrect\n');
fprintf('          cells (esp. deterministic blocks) are the main limiter.\n');

CLR_COR = [0.10 0.60 0.10];   % correct (green)
CLR_INC = [0.70 0.10 0.10];   % incorrect (red)
CLR_D   = [0.15 0.45 0.70];   % deterministic (blue)
CLR_P   = [0.80 0.30 0.10];   % probabilistic (orange)

fig_thdiff = figure('Position',[50 50 1100 800]);
sgtitle({'Theta ERSP differences (4-8 Hz, 200-400 ms, FCz)', ...
    'Lines = within-subject pairs; \Delta = mean difference (dB)'}, 'FontSize',12);

% 1. D blocks: Incorrect − Correct
ax1 = subplot(2,2,1);
paired_theta_panel(ax1, thD_c, thD_i, 'Correct','Incorrect', CLR_COR, CLR_INC, ...
    'D blocks: Incorrect − Correct');

% 2. P blocks: Incorrect − Correct
ax2 = subplot(2,2,2);
paired_theta_panel(ax2, thP_c, thP_i, 'Correct','Incorrect', CLR_COR, CLR_INC, ...
    'P blocks: Incorrect − Correct');

% 3. Correct trials: D − P
ax3 = subplot(2,2,3);
paired_theta_panel(ax3, thP_c, thD_c, 'P','D', CLR_P, CLR_D, ...
    'Correct trials: D − P');

% 4. Incorrect trials: D − P
ax4 = subplot(2,2,4);
paired_theta_panel(ax4, thP_i, thD_i, 'P','D', CLR_P, CLR_D, ...
    'Incorrect trials: D − P');

annotation(fig_thdiff,'textbox',[0.01 0.005 0.98 0.05],'String', ...
    sprintf(['n = subjects with >= %d true-feedback trials in BOTH conditions of the contrast ' ...
    '(per-subject means, collapsed across stages). The wavelet loop processed %d KH + %d RR datasets; ' ...
    'paired incorrect-trial cells (especially in deterministic blocks) are the main reason n is ' ...
    'smaller than the total number of datasets. Significance bar: * p<.05, ** p<.01, *** p<.001 (paired t-test).'], ...
    MIN_TRIALS_DIFF, numel(KH_PARTICIPANTS), numel(RR_PARTICIPANTS)), ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig_thdiff, fullfile(figure_output_folder,'S10_TF_theta_difference.pdf'));
fprintf('Saved S10_TF_theta_difference.pdf\n');

fprintf('\nAll done.\n');


% =========================================================================
%% LOCAL FUNCTIONS  (all definitions live at the END of the file)
% =========================================================================

% --------------------------------------------------------------------------
function paired_theta_panel(ax, va, vb, lbl_a, lbl_b, clr_a, clr_b, ttl)
%PAIRED_THETA_PANEL  Two-condition paired comparison of theta ERSP.
%   Left bar  = va (lbl_a); right bar = vb (lbl_b).
%   Delta and the paired t-test are reported as (vb - va), matching the
%   "<lbl_b> - <lbl_a>" convention used in the panel title. A significance
%   bracket (* p<.05, ** p<.01, *** p<.001) is drawn when the paired
%   difference is significant.
ok = ~isnan(va) & ~isnan(vb);
va = va(ok);  vb = vb(ok);
n  = numel(va);
hold(ax,'on');

mA = mean(va);  seA = std(va,'omitnan')/sqrt(max(n,1));
mB = mean(vb);  seB = std(vb,'omitnan')/sqrt(max(n,1));

bar(ax, 1, mA, 0.6, 'FaceColor', clr_a, 'FaceAlpha', 0.65, 'EdgeColor','none');
bar(ax, 2, mB, 0.6, 'FaceColor', clr_b, 'FaceAlpha', 0.65, 'EdgeColor','none');

% Within-subject connecting lines
for i = 1:n
    hl = plot(ax, [1 2], [va(i) vb(i)], '-', 'Color', [0.5 0.5 0.5], ...
        'LineWidth', 0.6, 'HandleVisibility','off');
    try, hl.Color(4) = 0.30; catch, end   % transparency (R2019b+)
end
scatter(ax, ones(n,1),   va, 20, clr_a, 'filled', 'MarkerFaceAlpha', 0.55, 'HandleVisibility','off');
scatter(ax, 2*ones(n,1), vb, 20, clr_b, 'filled', 'MarkerFaceAlpha', 0.55, 'HandleVisibility','off');

errorbar(ax, [1 2], [mA mB], [seA seB], ...
    'k.', 'LineWidth', 1.3, 'CapSize', 7, 'HandleVisibility','off');

yline(ax, 0, 'k:', 'HandleVisibility','off');
set(ax, 'XTick', [1 2], 'XTickLabel', {lbl_a, lbl_b}, 'TickDir','out');
xlim(ax, [0.4 2.6]);
ylabel(ax, 'Theta ERSP (dB)');
title(ax, sprintf('%s  (n=%d)', ttl, n), 'FontSize', 9);

if n < 2, return; end

d_mean        = mean(vb - va, 'omitnan');
[~, p, ~, st] = ttest(vb, va);
if     p < 0.001, star = '***';
elseif p < 0.01,  star = '**';
elseif p < 0.05,  star = '*';
else,             star = 'ns';
end

% y-extent of all plotted data, with headroom for the significance bracket
y_hi = max([va; vb; mA+seA; mB+seB]);
y_lo = min([va; vb; mA-seA; mB-seB; 0]);
rng  = max(y_hi - y_lo, eps);
ylim(ax, [y_lo - 0.08*rng, y_hi + 0.30*rng]);

% Significance bracket (only when significant)
if p < 0.05
    yb   = y_hi + 0.12*rng;
    tick = 0.04*rng;
    line(ax, [1 1 2 2], [yb-tick yb yb yb-tick], ...
        'Color','k', 'LineWidth',1.2, 'HandleVisibility','off');
    text(ax, 1.5, yb + 0.01*rng, star, ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'FontSize',13, 'FontWeight','bold');
end

% Stats line under the title
subtitle(ax, sprintf('\\Delta = %.2f dB,  t(%d) = %.2f, p = %.3f %s', ...
    d_mean, st.df, st.tstat, p, star), 'FontSize',8, 'Color',[0.3 0.3 0.3]);
end


% --------------------------------------------------------------------------
function m = mean_if_enough(v, min_n)
%MEAN_IF_ENOUGH  Mean of v (NaNs ignored), or NaN if fewer than min_n values.
v = v(~isnan(v));
if numel(v) >= min_n
    m = mean(v);
else
    m = NaN;
end
end


% --------------------------------------------------------------------------
function stat = run_tf_cluster_permtest(diff_maps, n_perm, alpha_clust)
%RUN_TF_CLUSTER_PERMTEST  Two-dimensional paired cluster permutation test.
%
%   diff_maps : n_subj x n_freq x n_time  (condition A minus condition B,
%               one difference map per subject — already the paired contrast)
%   n_perm    : number of sign-flip permutations (default 10000)
%   alpha_clust: p-value threshold for a cluster to be reported (e.g. 0.025)
%
%   SIGN-FLIP PERMUTATION (paired design):
%   Under H0 (no condition effect), A-B and B-A are equally likely.
%   Randomly multiplying each subject's difference map by ±1 simulates this.

if nargin < 2 || isempty(n_perm),       n_perm      = 10000; end
if nargin < 3 || isempty(alpha_clust),  alpha_clust = 0.025; end

[n_subj, n_freq, n_time] = size(diff_maps);

% ── Step 1: observed t-statistic map ─────────────────────────────────────
t_obs = paired_tmap(diff_maps);          % n_freq x n_time

% Critical t-value (two-tailed, df = n_subj-1)
t_crit = tinv(1 - alpha_clust, n_subj - 1);

% ── Step 2: observed clusters ─────────────────────────────────────────────
[pos_label, pos_stats] = find_clusters(t_obs,  t_crit);   % t > +t_crit
[neg_label, neg_stats] = find_clusters(-t_obs, t_crit);   % t < -t_crit

% ── Step 3: null distribution via sign-flip permutations ─────────────────
null_max_pos = zeros(1, n_perm);
null_max_neg = zeros(1, n_perm);

for p = 1:n_perm
    signs     = sign(randn(n_subj, 1, 1));   % each ±1 with p=0.5
    perm_maps = diff_maps .* signs;          % broadcast across freq x time

    t_perm = paired_tmap(perm_maps);

    [~, perm_pos_stats] = find_clusters(t_perm,  t_crit);
    [~, perm_neg_stats] = find_clusters(-t_perm, t_crit);

    null_max_pos(p) = max([0, perm_pos_stats]);   % max cluster t-sum
    null_max_neg(p) = max([0, perm_neg_stats]);
end

% ── Step 4: cluster p-values (proportion of null >= observed) ─────────────
n_pos_obs = numel(pos_stats);
n_neg_obs = numel(neg_stats);

posclusters = struct('prob', {}, 'clusterstat', {});
negclusters = struct('prob', {}, 'clusterstat', {});

for k = 1:n_pos_obs
    posclusters(k).clusterstat = pos_stats(k);
    posclusters(k).prob        = mean(null_max_pos >= pos_stats(k));
end
for k = 1:n_neg_obs
    negclusters(k).clusterstat = neg_stats(k);
    negclusters(k).prob        = mean(null_max_neg >= neg_stats(k));
end

% Sort by p-value (ascending) — most significant first
if ~isempty(posclusters)
    [~, ord] = sort([posclusters.prob]);
    posclusters = posclusters(ord);
    [~, ord_lab] = sort(pos_stats, 'descend');  %#ok<NASGU>
end
if ~isempty(negclusters)
    [~, ord] = sort([negclusters.prob]);
    negclusters = negclusters(ord);
end

% ── Step 5: significance mask ─────────────────────────────────────────────
% Note: posclusters has been re-sorted, so test against clusterstat-matched
% labels rather than the (now stale) ordinal index k.
sig_mask = false(n_freq, n_time);
for k = 1:numel(posclusters)
    if posclusters(k).prob < alpha_clust
        lab = find(pos_stats == posclusters(k).clusterstat, 1);
        if ~isempty(lab), sig_mask = sig_mask | (pos_label == lab); end
    end
end
for k = 1:numel(negclusters)
    if negclusters(k).prob < alpha_clust
        lab = find(neg_stats == negclusters(k).clusterstat, 1);
        if ~isempty(lab), sig_mask = sig_mask | (neg_label == lab); end
    end
end

% ── Assemble output struct (mirrors FieldTrip stat struct fields used by
%    plot_tf_cluster_result) ────────────────────────────────────────────────
stat.stat                = t_obs;
stat.mask                = sig_mask;
stat.posclusters         = posclusters;
stat.negclusters         = negclusters;
stat.posclusterslabelmat = pos_label;
stat.negclusterslabelmat = neg_label;
stat.null_max_pos        = null_max_pos;
stat.null_max_neg        = null_max_neg;
stat.t_crit              = t_crit;
stat.df                  = n_subj - 1;
stat.n_perm              = n_perm;
stat.alpha_clust         = alpha_clust;

end


% --------------------------------------------------------------------------
function [t_map, p_map] = paired_tmap(diff_maps)
%PAIRED_TMAP  Paired t-statistic at each freq x time point.
%   diff_maps : n_subj x n_freq x n_time
%   t = mean(diff) / (std(diff) / sqrt(n))

n_subj = size(diff_maps, 1);
mn     = mean(diff_maps, 1, 'omitnan');   % 1 x n_freq x n_time
sd     = std(diff_maps,  0, 1, 'omitnan');

se     = sd / sqrt(n_subj);
t_map  = squeeze(mn ./ max(se, 1e-10));   % n_freq x n_time

if nargout > 1
    p_map = 2 * tcdf(-abs(t_map), n_subj - 1);
end

end


% --------------------------------------------------------------------------
function [label_map, cluster_stats] = find_clusters(t_map, t_crit)
%FIND_CLUSTERS  Connected components above t_crit in a 2-D t-map.
%   8-connectivity in the freq x time plane (FieldTrip / SPM default).

thresh_map = t_map > t_crit;   % binary map of suprathreshold points

if ~any(thresh_map(:))
    label_map     = zeros(size(t_map));
    cluster_stats = [];
    return
end

CC = bwconncomp(thresh_map, 8);

label_map     = zeros(size(t_map));
cluster_stats = zeros(1, CC.NumObjects);

for k = 1:CC.NumObjects
    idx              = CC.PixelIdxList{k};
    label_map(idx)   = k;
    cluster_stats(k) = sum(t_map(idx));   % cluster mass statistic
end

end


% --------------------------------------------------------------------------
function report_tf_clusters(stat, contrast_name)
%REPORT_TF_CLUSTERS  Print a one-line summary of significant clusters.

fprintf('  %s:\n', contrast_name);

if ~isempty(stat.posclusters)
    for k = 1:numel(stat.posclusters)
        sig = '';
        if stat.posclusters(k).prob < stat.alpha_clust, sig = '  ***'; end
        fprintf('    pos cluster %d: p=%.4f  mass=%.1f%s\n', k, ...
            stat.posclusters(k).prob, stat.posclusters(k).clusterstat, sig);
    end
else
    fprintf('    no positive clusters above threshold\n');
end

if ~isempty(stat.negclusters)
    for k = 1:numel(stat.negclusters)
        sig = '';
        if stat.negclusters(k).prob < stat.alpha_clust, sig = '  ***'; end
        fprintf('    neg cluster %d: p=%.4f  mass=%.1f%s\n', k, ...
            stat.negclusters(k).prob, stat.negclusters(k).clusterstat, sig);
    end
else
    fprintf('    no negative clusters above threshold\n');
end

end


% =========================================================================
% CORE WAVELET FUNCTION
%
% morlet_wavelet_tf: full time-frequency decomposition.
%   sig       — n_times x n_epochs  (double, mean-centred recommended)
%   Fs        — sampling rate (Hz)
%   frex      — 1 x N_FREQS  frequency vector (Hz)
%   n_cycles  — 1 x N_FREQS  number of cycles per frequency
% Returns:
%   tf_power  — N_FREQS x n_times x n_epochs  (instantaneous power)
%   tf_phase  — N_FREQS x n_times x n_epochs  (instantaneous phase, rad)
% =========================================================================
function [tf_power, tf_phase] = morlet_wavelet_tf(sig, Fs, frex, n_cycles)

[n_times, n_epochs] = size(sig);
n_freqs = numel(frex);

tf_power = zeros(n_freqs, n_times, n_epochs, 'single');
tf_phase = zeros(n_freqs, n_times, n_epochs, 'single');

% FFT of signal (all epochs simultaneously); zero-pad to power of 2.
nConv_base = n_times * 2 - 1;
nConv      = 2^nextpow2(nConv_base);

sig_fft = fft(sig, nConv, 1);   % nConv x n_epochs

% Wavelet time axis: centred on zero, long enough for the lowest frequency.
sigma_t_max = n_cycles(1) / (2 * pi * frex(1));
wt_len      = ceil(3 * sigma_t_max * Fs);   % 3 standard deviations
wt          = (-wt_len : wt_len) / Fs;
n_wt        = length(wt);
half_wt     = floor(n_wt / 2);

for fi = 1:n_freqs
    f       = frex(fi);
    sigma_t = n_cycles(fi) / (2 * pi * f);

    % Complex Morlet wavelet
    gauss_env = exp(-wt.^2 / (2 * sigma_t^2));
    sine_wave = exp(2 * pi * 1i * f * wt);
    wavelet   = sine_wave .* gauss_env;
    wavelet   = wavelet / sum(abs(wavelet));   % normalise

    wavelet_fft = fft(wavelet(:), nConv);

    conv_result = ifft(wavelet_fft .* sig_fft, nConv, 1);   % nConv x n_epochs

    % Trim to valid region (remove edge artefacts)
    trim_start   = half_wt + 1;
    trim_end     = trim_start + n_times - 1;
    conv_trimmed = conv_result(trim_start:trim_end, :);     % n_times x n_epochs

    tf_power(fi,:,:) = single(abs(conv_trimmed).^2);
    tf_phase(fi,:,:) = single(angle(conv_trimmed));
end

end


% =========================================================================
% ERSP COMPUTATION (baseline-corrected dB power)
%   cond_power — N_FREQS x n_times x n_trials
%   bl_mask    — logical, n_times x 1
%   ersp       — N_FREQS x n_times  (mean across trials, in dB)
% =========================================================================
function ersp = compute_ersp(cond_power, bl_mask)

mean_power = mean(cond_power, 3, 'omitnan');             % N_FREQS x n_times
bl_power   = mean(mean_power(:, bl_mask), 2, 'omitnan'); % N_FREQS x 1
ersp       = 10 * log10(mean_power ./ bl_power);

end


% =========================================================================
% PER-SUBJECT TF FIGURE
% =========================================================================
function plot_subject_tf(subj, participant, has_P, group_table, subj_rows, ...
    tf_power, tf_phase, tf_power_pz, ...
    frex, t_global, bl_mask, plot_mask, ...
    stage_names, BTYPE_LABELS, STAGE_COLORS, ...
    theta_band, FRN_WIN, N_FREQS, MIN_TRIALS, ...
    figure_output_folder, n_epochs) %#ok<INUSD>

t_plot  = t_global(plot_mask);

fig = figure('Position',[50 50 1400 600], 'Visible','off');
sgtitle(sprintf('%s — TF maps: FCz, D blocks', subj), 'Interpreter','none');

clim_ersp = [-3 3];
clim_itpc = [0 0.4];

bt = 'D';
for s_i = 1:4
    for oc_i = 1:2
        oc_val  = 2 - oc_i;
        oc_str  = {'correct','incorrect'};

        cond_mask = subj_rows & ...
            group_table.block_type==bt & ...
            group_table.stage==stage_names{s_i} & ...
            group_table.correct==oc_val & ...
            ~group_table.false_fb;

        ep_vec = group_table.epoch(cond_mask);
        ep_vec = ep_vec(~isnan(ep_vec) & ep_vec>=1 & ep_vec<=n_epochs);

        ax = subplot(2, 8, (oc_i-1)*8 + s_i);

        if numel(ep_vec) >= MIN_TRIALS
            cond_power = tf_power(:,:,ep_vec);
            ersp = compute_ersp(cond_power, bl_mask);
            imagesc(ax, t_plot, frex, ersp(:,plot_mask));
            axis(ax,'xy');
            set(ax,'YTick',[4 8 12 20 30],'YScale','log');
            clim(ax, clim_ersp);
            colormap(ax, tf_colormap());
            xline(ax,0,'w:'); yline(ax,4,'w--'); yline(ax,8,'w--');
        else
            text(ax,0.5,0.5,sprintf('n=%d\n(too few)',numel(ep_vec)), ...
                'HorizontalAlignment','center','Units','normalized','FontSize',8);
            set(ax,'XTick',[],'YTick',[]);
        end
        title(ax, sprintf('%s|%s|%s', stage_names{s_i}, bt, oc_str{oc_i}), 'FontSize',7);
        xlabel(ax,'ms'); ylabel(ax,'Hz');

        % ITPC
        ax_it = subplot(2, 8, (oc_i-1)*8 + s_i + 4);
        if numel(ep_vec) >= MIN_TRIALS
            cond_phase = tf_phase(:,:,ep_vec);
            itpc = abs(mean(exp(1i*cond_phase), 3, 'omitnan'));
            imagesc(ax_it, t_plot, frex, itpc(:,plot_mask));
            axis(ax_it,'xy');
            set(ax_it,'YTick',[4 8 12 20 30],'YScale','log');
            clim(ax_it, clim_itpc);
            colormap(ax_it, parula);
            xline(ax_it,0,'w:'); yline(ax_it,4,'w--'); yline(ax_it,8,'w--');
            title(ax_it, sprintf('ITPC|%s|%s|%s', stage_names{s_i}, bt, oc_str{oc_i}), 'FontSize',7);
        end
        xlabel(ax_it,'ms'); ylabel(ax_it,'Hz');
    end
end

saveas(fig, fullfile(figure_output_folder, sprintf('%s_TF_S10.pdf', subj)));
close(fig);
end


% =========================================================================
% PLOT TF CLUSTER PERMUTATION RESULT
% =========================================================================
function plot_tf_cluster_result(stat, frex, t_plot, name, outfolder, PLOT_WIN) %#ok<INUSD>

fig = figure('Position',[50 50 900 400]);
sgtitle(sprintf('TF cluster result: %s', name), 'Interpreter','none');

% T-statistic map
ax1 = subplot(1,2,1);
imagesc(ax1, t_plot, frex, stat.stat);
axis(ax1,'xy');
set(ax1,'YTick',[4 8 12 20 30],'YScale','log');
colormap(ax1, tf_colormap());
colorbar(ax1);
title(ax1,'T-statistic map','FontSize',9);
xlabel(ax1,'Time (ms)'); ylabel(ax1,'Frequency (Hz)');
xline(ax1,0,'w:'); yline(ax1,4,'w--'); yline(ax1,8,'w--');

% Significant cluster mask
ax2 = subplot(1,2,2);
sig_mask = zeros(size(stat.stat));
if isfield(stat,'posclusters')
    for c = 1:numel(stat.posclusters)
        if stat.posclusters(c).prob < 0.025
            lab = find([stat.posclusters.clusterstat] == stat.posclusters(c).clusterstat, 1);
            sig_mask = sig_mask + double(stat.posclusterslabelmat == lab);
        end
    end
end
if isfield(stat,'negclusters')
    for c = 1:numel(stat.negclusters)
        if stat.negclusters(c).prob < 0.025
            lab = find([stat.negclusters.clusterstat] == stat.negclusters(c).clusterstat, 1);
            sig_mask = sig_mask - double(stat.negclusterslabelmat == lab);
        end
    end
end
imagesc(ax2, t_plot, frex, sig_mask);
axis(ax2,'xy');
set(ax2,'YTick',[4 8 12 20 30],'YScale','log');
colormap(ax2, tf_colormap());
clim(ax2,[-1 1]);
title(ax2,'Significant clusters (p<0.025)','FontSize',9);
xlabel(ax2,'Time (ms)'); ylabel(ax2,'Frequency (Hz)');
xline(ax2,0,'w:'); yline(ax2,4,'w--'); yline(ax2,8,'w--');

saveas(fig, fullfile(outfolder, sprintf('S10_TF_cluster_%s.pdf', name)));
close(fig);
end


% =========================================================================
% CUSTOM TF COLORMAP (blue-white-red)
% Blue = power decrease, white = no change, red = power increase.
% =========================================================================
function cmap = tf_colormap(n)
if nargin < 1; n = 256; end
half = floor(n/2);
r1 = linspace(0.20, 1.0, half)';
g1 = linspace(0.40, 1.0, half)';
b1 = linspace(0.80, 1.0, half)';
r2 = linspace(1.0, 0.80, n-half)';
g2 = linspace(1.0, 0.20, n-half)';
b2 = linspace(1.0, 0.15, n-half)';
cmap = [r1,g1,b1; r2,g2,b2];
end


% =========================================================================
% SAFE CHANNEL INDEX HELPER
% =========================================================================
function idx = safe_chan(EEGp, labels)
idx = find(ismember(lower({EEGp.chanlocs.labels}), lower(labels)));
end
