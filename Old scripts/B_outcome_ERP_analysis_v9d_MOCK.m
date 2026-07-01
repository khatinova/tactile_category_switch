% =============================================================================
% PLOT_OUTCOME_ERPS_KH_v4_QC_debugged.m
%
% PURPOSE
% -------
% Dedicated QC plotting script for KH outcome ERPs produced by
% A_Preprocessing_KH_v4.m in Epoched_data_noisefiltering/.
%
% Current setting:
%   - plots only KH Ox12 and Ox14
%
% Ready for full KH cohort:
%   - set PLOT_ALL_VALID_KH = true
%
% What it checks:
%   1. Outcome ERPs at FCz and FCz/Cz by block type, accuracy, and stage.
%   2. Artifact-flagged epochs from v4 preprocessing.
%   3. Whether residual high-frequency activity is still visible after the
%      harmonising low-pass/noise-filtering step.
%
% Major debug fixes vs previous plotting/analysis script:
%   FIX 1: v4 preprocessing saves epoch_artifact_flag_trimmed, not always
%          epoch_artifact_flag.
%   FIX 2: trial2epoch saved by preprocessing is indexed to the original
%          untrimmed outcome file. If plotting *_trimmed.set, remap original
%          epoch indices to trimmed epoch indices using valid_ep.
%   FIX 3: hw_filter_label may be saved as '0-30Hz'/'0-70Hz', while the
%          analysis script expected '0_30Hz'/'0_70Hz'. Normalise both.
%   FIX 4: robust participant switch: all valid KH participants are enabled,
%          but currently disabled in favour of Ox12/Ox14 pilot.
%   FIX 5: avoid RR code paths and feature-extraction dependencies when the
%          goal is only plotting/QC.
% =============================================================================

clear; close all; clc;

% -------------------------------------------------------------------------
%% USER SWITCHES
% -------------------------------------------------------------------------

PLOT_ALL_VALID_KH = false;          % true = plot all valid KH subjects
KH_VALID_FULL     = [3:12, 14:23, 27];
KH_PILOT_ONLY     = [15];       % current setting requested

if PLOT_ALL_VALID_KH
    participants_to_plot = KH_VALID_FULL;
else
    participants_to_plot = KH_PILOT_ONLY;
end

save_figures  = true;
show_figures  = false;              % set true for interactive inspection
make_group_qc = true;

% ERP/QC parameters
ERP_plot_window = [-200 1000];       % ms
baseline_win    = [-200 0];          % ms
FRN_win         = [250 350];         % ms
HF_BAND         = [30 45];           % Hz; should be strongly reduced after harmonising LP
LF_REF_BAND     = [1 30];            % Hz
PSD_MAX_HZ      = 80;

stage_names  = {'LN','LE','RN','RE'};
btype_labels = {'D','P'};

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

KH_data_path = fullfile(base_path, 'Salient mod switch KH', 'Data');

epoch_file_folder = fullfile(base_path, 'Salient mod switch KH', 'Results', ...
    'EEG analysis', 'Epoched_data_noisefiltering');

figure_output_folder = fullfile(base_path, 'Salient mod switch KH', 'Results', ...
    'EEG analysis', 'Figures', 'outcome_v4_noise_QC_debugged');

if ~exist(figure_output_folder, 'dir')
    mkdir(figure_output_folder);
end

if ~exist(epoch_file_folder, 'dir')
    error('Epoch folder not found: %s', epoch_file_folder);
end

addpath(eeglab_path);
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab('nogui'); %#ok<ASGLU>

load(fullfile(KH_data_path, 'all_trial_data.mat'));

behav_file = fullfile(KH_data_path, 'behav_table.mat');
if exist(behav_file, 'file')
    S_beh = load(behav_file, 'group_T');
    KH_behav_table = S_beh.group_T;
else
    error('behav_table.mat not found at %s', behav_file);
end

if ismember('researcher', KH_behav_table.Properties.VariableNames)
    KH_behav_table = KH_behav_table(strcmp(string(KH_behav_table.researcher), 'KH'), :);
end

% -------------------------------------------------------------------------
%% OUTPUT CONTAINERS
% -------------------------------------------------------------------------

qc_summary = table();

grand_clean_FCz = struct();
grand_all_FCz   = struct();

for si = 1:numel(stage_names)
    st = stage_names{si};
    for bi = 1:numel(btype_labels)
        bt = btype_labels{bi};
        grand_clean_FCz.(st).(bt).correct   = [];
        grand_clean_FCz.(st).(bt).incorrect = [];
        grand_all_FCz.(st).(bt).correct     = [];
        grand_all_FCz.(st).(bt).incorrect   = [];
    end
end

% -------------------------------------------------------------------------
%% MAIN PARTICIPANT LOOP
% -------------------------------------------------------------------------

for p = participants_to_plot

    subj = sprintf('Ox%02d', p);
    fprintf('\n==================== %s ====================\n', subj);

    if ~isfield(all_trial_data, subj)
        warning('%s missing from all_trial_data. Skipping.', subj);
        continue;
    end

    % ---------------------------------------------------------------------
    % Load broadband outcome epochs. Prefer trimmed v4 output.
    % ---------------------------------------------------------------------
    epoch_candidates = { ...
        sprintf('%s_outcome_trimmed.set', subj), ...
        sprintf('%s_outcome.set', subj)};

    EEGp = [];
    epoch_file_used = '';
    used_trimmed_file = false;

    for ci = 1:numel(epoch_candidates)
        fpath = fullfile(epoch_file_folder, epoch_candidates{ci});
        if exist(fpath, 'file')
            EEGp = pop_loadset(epoch_candidates{ci}, epoch_file_folder);
            epoch_file_used = epoch_candidates{ci};
            used_trimmed_file = contains(epoch_file_used, '_trimmed');
            break;
        end
    end

    if isempty(EEGp)
        warning('%s: no outcome epoch file found. Skipping.', subj);
        continue;
    end

    fprintf('Loaded %s: %d channels x %d timepoints x %d epochs\n', ...
        epoch_file_used, EEGp.nbchan, EEGp.pnts, EEGp.trials);

    % ---------------------------------------------------------------------
    % Channel lookup
    % ---------------------------------------------------------------------
    labels_lower = lower(string({EEGp.chanlocs.labels}));

    fcz_idx = find(labels_lower == "fcz", 1);
    cz_idx  = find(labels_lower == "cz", 1);

    if isempty(fcz_idx)
        warning('%s: FCz not found. Skipping.', subj);
        continue;
    end

    if isempty(cz_idx)
        warning('%s: Cz not found. FCz/Cz plots will use FCz only.', subj);
        fczcz_idx = fcz_idx;
    else
        fczcz_idx = [fcz_idx cz_idx];
    end

    % ---------------------------------------------------------------------
    % Time masks
    % ---------------------------------------------------------------------
    in_win   = EEGp.times >= ERP_plot_window(1) & EEGp.times <= ERP_plot_window(2);
    bl_mask  = EEGp.times >= baseline_win(1)    & EEGp.times <= baseline_win(2);
    frn_mask = EEGp.times >= FRN_win(1)         & EEGp.times <= FRN_win(2);

    tt = EEGp.times(in_win);

    if ~any(bl_mask)
        error('%s: baseline window [%g %g] ms not present in EEGp.times.', ...
            subj, baseline_win(1), baseline_win(2));
    end

    % ---------------------------------------------------------------------
    % Behavioural spine for this subject
    % ---------------------------------------------------------------------
    if ~ismember('subjID', KH_behav_table.Properties.VariableNames)
        error('KH_behav_table lacks subjID column.');
    end

    subj_rows = string(KH_behav_table.subjID) == string(subj);
    if ~any(subj_rows)
        warning('%s has no rows in KH_behav_table. Skipping.', subj);
        continue;
    end

    sf = KH_behav_table(subj_rows, :);
    n_rows = height(sf);

    sf.subj_id = repmat(string(subj), n_rows, 1);
    sf.subj    = repmat(p, n_rows, 1);

    if ~ismember('trial_continuous', sf.Properties.VariableNames)
        sf.trial_continuous = (1:n_rows)';
    end

    if ~ismember('block_type', sf.Properties.VariableNames)
        warning('%s: block_type missing. Assigning all rows to D for plotting.', subj);
        sf.block_type = categorical(repmat("D", n_rows, 1), btype_labels);
    else
        sf.block_type = categorical(string(sf.block_type), btype_labels);
    end

    if ~ismember('false_fb', sf.Properties.VariableNames)
        sf.false_fb = false(n_rows, 1);
    else
        sf.false_fb = logical(local_to_numeric(sf.false_fb));
    end

    if ~ismember('correct', sf.Properties.VariableNames)
        error('%s: behavioural table lacks correct column.', subj);
    end
    sf.correct_num = local_correct_to_numeric(sf.correct);

    % Assign/reassign stages using the same logic as the analysis script.
    block_col = 'block';
    trial_col = 'trial';

    if ~ismember(block_col, sf.Properties.VariableNames)
        error('%s: behavioural table lacks block column.', subj);
    end
    if ~ismember(trial_col, sf.Properties.VariableNames)
        error('%s: behavioural table lacks trial column.', subj);
    end

    rev_trials_vec = [];
    beh = all_trial_data.(subj).trial_data;
    if isfield(beh, 'revTrial') && ~isempty(beh.revTrial)
        rev_trials_vec = beh.revTrial(:);
    end

    sf = assign_stages_preserve_LE_RN(sf, block_col, trial_col, rev_trials_vec, stage_names);

    % ---------------------------------------------------------------------
    % Load trial2epoch cache and remap to the loaded epoch file
    % ---------------------------------------------------------------------
    cache_file = fullfile(epoch_file_folder, sprintf('%s_trial2epoch.mat', subj));

    if ~exist(cache_file, 'file')
        warning('%s: trial2epoch cache not found. Falling back to sequential mapping.', subj);
        trial2epoch_mapped = nan(n_rows, 1);
        n_map = min(n_rows, EEGp.trials);
        trial2epoch_mapped(1:n_map) = (1:n_map)';
        valid_ep = [];
        hw_filter_label = "unknown";
        epoch_artifact_flag_eeg = false(EEGp.trials, 1);
    else
        S2 = load(cache_file);

        if isfield(S2, 'trial2epoch_out')
            trial2epoch_raw = S2.trial2epoch_out(:);
        elseif isfield(S2, 'trial2epoch')
            trial2epoch_raw = S2.trial2epoch(:);
        else
            error('%s: cache lacks trial2epoch/trial2epoch_out.', cache_file);
        end

        if isfield(S2, 'valid_ep')
            valid_ep = S2.valid_ep(:);
        else
            valid_ep = [];
        end

        if isfield(S2, 'hw_filter_label') && ~isempty(S2.hw_filter_label)
            hw_filter_label = normalize_hw_filter_label(S2.hw_filter_label);
        else
            hw_filter_label = "unknown";
        end

        trial2epoch_raw = pad_or_truncate(trial2epoch_raw, n_rows);

        % Critical fix:
        % If the loaded file is *_trimmed.set and valid_ep exists, the
        % saved trial2epoch values refer to original untrimmed epoch IDs.
        % Remap original epoch IDs to the current trimmed epoch order.
        trial2epoch_mapped = remap_trial2epoch_to_loaded_file( ...
            trial2epoch_raw, valid_ep, EEGp.trials, used_trimmed_file);

        % Critical fix:
        % v4 preprocessing saves epoch_artifact_flag_trimmed after trimming.
        epoch_artifact_flag_eeg = read_epoch_artifact_flags(S2, valid_ep, EEGp.trials);
    end

    sf.epoch = trial2epoch_mapped;
    sf.has_eeg_epoch = ~isnan(sf.epoch) & sf.epoch >= 1 & sf.epoch <= EEGp.trials;

    sf.epoch_artifact_flag = false(n_rows, 1);
    idx_valid = find(sf.has_eeg_epoch);
    ep_idx = round(sf.epoch(idx_valid));
    sf.epoch_artifact_flag(idx_valid) = epoch_artifact_flag_eeg(ep_idx);

    fprintf('Mapped behavioural rows to EEG epochs: %d/%d\n', ...
        sum(sf.has_eeg_epoch), n_rows);
    fprintf('Artifact-flagged mapped trials: %d/%d (%.1f%%)\n', ...
        sum(sf.epoch_artifact_flag), sum(sf.has_eeg_epoch), ...
        100 * sum(sf.epoch_artifact_flag) / max(1, sum(sf.has_eeg_epoch)));
    fprintf('hw_filter_label: %s\n', hw_filter_label);

    % ---------------------------------------------------------------------
    % Precompute single-trial FCz and FCz/Cz waveforms, baseline corrected
    % ---------------------------------------------------------------------
    fczt = nan(n_rows, EEGp.pnts);
    fczczt = nan(n_rows, EEGp.pnts);
    frn_mean = nan(n_rows, 1);

    for r = 1:n_rows
        if ~sf.has_eeg_epoch(r), continue; end
        ep = round(sf.epoch(r));

        sig_fcz = double(EEGp.data(fcz_idx, :, ep));
        sig_fcz = sig_fcz - mean(sig_fcz(bl_mask), 'omitnan');
        fczt(r, :) = sig_fcz;

        sig_fczcz = squeeze(mean(double(EEGp.data(fczcz_idx, :, ep)), 1, 'omitnan'));
        sig_fczcz = sig_fczcz(:)';
        sig_fczcz = sig_fczcz - mean(sig_fczcz(bl_mask), 'omitnan');
        fczczt(r, :) = sig_fczcz;

        frn_mean(r) = mean(sig_fczcz(frn_mask), 'omitnan');
    end

    sf.FRN_mean_amp = frn_mean;

    % ---------------------------------------------------------------------
    % High-frequency diagnostic
    % ---------------------------------------------------------------------
    eps_all = round(sf.epoch(sf.has_eeg_epoch));
    eps_all = eps_all(eps_all >= 1 & eps_all <= EEGp.trials);

    hf_ratio_FCz = NaN;
    hf_abs_FCz   = NaN;
    psd_f        = [];
    psd_mean     = [];

    if ~isempty(eps_all)
        [psd_f, psd_mean, hf_ratio_FCz, hf_abs_FCz] = compute_subject_psd_hf_ratio( ...
            EEGp, fcz_idx, eps_all, LF_REF_BAND, HF_BAND, PSD_MAX_HZ);
    end

    fprintf('HF diagnostic FCz: mean PSD %g-%g Hz = %.4g; HF/LF ratio = %.4g\n', ...
        HF_BAND(1), HF_BAND(2), hf_abs_FCz, hf_ratio_FCz);

    % ---------------------------------------------------------------------
    % Per-subject ERP plot: 4 stages x 2 accuracy panels
    % ---------------------------------------------------------------------
    fig1 = figure('Position', [50 50 1500 760], 'Visible', onoff(show_figures));
    sgtitle(sprintf('%s outcome ERP QC | %s | %s', subj, epoch_file_used, hw_filter_label), ...
        'Interpreter', 'none');

    for s_i = 1:numel(stage_names)
        st = stage_names{s_i};

        for oc_i = 1:2
            oc_val = oc_i - 1;  % 0 incorrect, 1 correct
            oc_label = "incorrect";
            if oc_val == 1, oc_label = "correct"; end

            ax = subplot(2, 4, (oc_i - 1) * 4 + s_i);
            hold(ax, 'on');

            title(ax, sprintf('%s | %s', st, oc_label), 'Interpreter', 'none');
            xline(ax, 0, 'k:', 'HandleVisibility', 'off');
            yline(ax, 0, 'k:', 'HandleVisibility', 'off');

            for bi = 1:numel(btype_labels)
                bt = btype_labels{bi};

                m_clean = string(sf.stage) == st & ...
                    string(sf.block_type) == bt & ...
                    sf.correct_num == oc_val & ...
                    ~sf.false_fb & ...
                    sf.has_eeg_epoch & ...
                    ~sf.epoch_artifact_flag;

                dat = fczczt(m_clean, :);
                dat = dat(:, in_win);
                dat = dat(any(~isnan(dat), 2), :);

                if isempty(dat), continue; end

                mn = mean(dat, 1, 'omitnan');
                se = std(dat, 0, 1, 'omitnan') ./ sqrt(size(dat, 1));

                c = line_color_for(bt, oc_val);

                fill(ax, [tt fliplr(tt)], [mn + se fliplr(mn - se)], c, ...
                    'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
                    'HandleVisibility', 'off');

                ls = '-';
                if strcmp(bt, 'P'), ls = '--'; end

                plot(ax, tt, mn, 'Color', c, 'LineStyle', ls, 'LineWidth', 1.8, ...
                    'DisplayName', sprintf('%s n=%d', bt, size(dat, 1)));

                % store group clean FCz, not FCz/Cz
                m_group = m_clean;
                dat_group = fczt(m_group, :);
                dat_group = dat_group(:, in_win);
                dat_group = dat_group(any(~isnan(dat_group), 2), :);
                if ~isempty(dat_group)
                    subj_mean = mean(dat_group, 1, 'omitnan');
                    if oc_val == 1
                        grand_clean_FCz.(st).(bt).correct(end+1, :) = subj_mean; %#ok<SAGROW>
                    else
                        grand_clean_FCz.(st).(bt).incorrect(end+1, :) = subj_mean; %#ok<SAGROW>
                    end
                end
            end

            set(ax, 'YDir', 'reverse');
            xlabel(ax, 'Time (ms)');
            ylabel(ax, 'FCz/Cz \muV');
            xlim(ax, ERP_plot_window);
            legend(ax, 'Box', 'off', 'FontSize', 7, 'Location', 'best');
        end
    end

    if save_figures
        save_figure(fig1, fullfile(figure_output_folder, sprintf('%s_ERP_by_stage_debugged.pdf', subj)));
    end
    if ~show_figures, close(fig1); end

    % ---------------------------------------------------------------------
    % Per-subject preprocessing/noise QC figure
    % ---------------------------------------------------------------------
    fig2 = figure('Position', [80 80 1600 700], 'Visible', onoff(show_figures));
    sgtitle(sprintf('%s v4 noise-filtering QC | artifact flags + spectrum', subj), ...
        'Interpreter', 'none');

    % Panel 1: clean mean + artifact flagged single trials
    ax1 = subplot(2, 2, 1);
    hold(ax1, 'on');
    title(ax1, 'FCz all clean means; grey = artifact-flagged single trials');

    xline(ax1, 0, 'k:', 'HandleVisibility', 'off');
    yline(ax1, 0, 'k:', 'HandleVisibility', 'off');

    m_flag = sf.has_eeg_epoch & sf.epoch_artifact_flag & ~sf.false_fb;
    flag_dat = fczt(m_flag, in_win);
    flag_dat = flag_dat(any(~isnan(flag_dat), 2), :);

    for k = 1:size(flag_dat, 1)
        plot(ax1, tt, flag_dat(k, :), 'Color', [0.78 0.78 0.78], ...
            'LineWidth', 0.4, 'HandleVisibility', 'off');
    end

    for bi = 1:numel(btype_labels)
        bt = btype_labels{bi};
        for oc_val = 0:1
            m = string(sf.block_type) == bt & ...
                sf.correct_num == oc_val & ...
                ~sf.false_fb & ...
                sf.has_eeg_epoch & ...
                ~sf.epoch_artifact_flag;

            dat = fczt(m, in_win);
            dat = dat(any(~isnan(dat), 2), :);
            if isempty(dat), continue; end

            mn = mean(dat, 1, 'omitnan');
            c = line_color_for(bt, oc_val);
            ls = '-';
            if strcmp(bt, 'P'), ls = '--'; end

            oc_name = 'incorrect';
            if oc_val == 1, oc_name = 'correct'; end

            plot(ax1, tt, mn, 'Color', c, 'LineStyle', ls, 'LineWidth', 2, ...
                'DisplayName', sprintf('%s %s n=%d', bt, oc_name, size(dat, 1)));

            % all-trials group store
            m_all = string(sf.block_type) == bt & ...
                sf.correct_num == oc_val & ...
                ~sf.false_fb & ...
                sf.has_eeg_epoch;

            dat_all = fczt(m_all, in_win);
            dat_all = dat_all(any(~isnan(dat_all), 2), :);
            if ~isempty(dat_all)
                subj_mean_all = mean(dat_all, 1, 'omitnan');
                if oc_val == 1
                    grand_all_FCz.LN.(bt).correct(end+1, :) = nan(1, numel(tt)); %#ok<SAGROW>
                end
                % Stored by stage in the earlier loop only for clean; all-trial
                % group summary is kept subject-level in qc_summary.
            end
        end
    end

    set(ax1, 'YDir', 'reverse');
    xlabel(ax1, 'Time (ms)');
    ylabel(ax1, 'FCz \muV');
    xlim(ax1, ERP_plot_window);
    legend(ax1, 'Box', 'off', 'FontSize', 8, 'Location', 'best');

    % Panel 2: FRN scatter with artifact flags
    ax2 = subplot(2, 2, 2);
    hold(ax2, 'on');
    title(ax2, 'FRN mean amp by stage; grey = artifact flagged');

    for s_i = 1:numel(stage_names)
        st = stage_names{s_i};

        m_flag_stage = string(sf.stage) == st & sf.epoch_artifact_flag & ~isnan(sf.FRN_mean_amp);
        if any(m_flag_stage)
            scatter(ax2, s_i * ones(sum(m_flag_stage), 1), sf.FRN_mean_amp(m_flag_stage), ...
                42, [0.72 0.72 0.72], 'o', 'LineWidth', 1.1, ...
                'DisplayName', 'flagged');
        end

        for bi = 1:numel(btype_labels)
            bt = btype_labels{bi};
            for oc_val = 0:1
                m = string(sf.stage) == st & ...
                    string(sf.block_type) == bt & ...
                    sf.correct_num == oc_val & ...
                    ~sf.epoch_artifact_flag & ...
                    ~isnan(sf.FRN_mean_amp);

                if ~any(m), continue; end

                xj = s_i + (bi - 1.5) * 0.16 + (oc_val - 0.5) * 0.06;
                c = line_color_for(bt, oc_val);
                mk = 'o';
                if oc_val == 1, mk = '^'; end

                scatter(ax2, xj * ones(sum(m), 1), sf.FRN_mean_amp(m), ...
                    28, c, mk, 'filled', 'MarkerFaceAlpha', 0.45, ...
                    'HandleVisibility', 'off');
            end
        end
    end

    yline(ax2, 0, 'k--', 'HandleVisibility', 'off');
    set(ax2, 'XTick', 1:numel(stage_names), 'XTickLabel', stage_names, 'YDir', 'reverse');
    xlabel(ax2, 'Stage');
    ylabel(ax2, sprintf('FRN mean %d-%d ms, FCz/Cz \\muV', FRN_win(1), FRN_win(2)));

    % Panel 3: artifact flag rate
    ax3 = subplot(2, 2, 3);
    hold(ax3, 'on');
    title(ax3, 'Artifact flag rate by stage and block type');

    bw = 0.35;
    for bi = 1:numel(btype_labels)
        bt = btype_labels{bi};
        for s_i = 1:numel(stage_names)
            st = stage_names{s_i};
            m = string(sf.stage) == st & string(sf.block_type) == bt & sf.has_eeg_epoch;
            if ~any(m), continue; end

            rate = mean(sf.epoch_artifact_flag(m), 'omitnan');
            xpos = s_i + (bi - 1.5) * bw;

            bar(ax3, xpos, rate, bw * 0.9, ...
                'FaceColor', line_color_for(bt, 1), ...
                'EdgeColor', 'none');
        end
    end

    set(ax3, 'XTick', 1:numel(stage_names), 'XTickLabel', stage_names);
    xlabel(ax3, 'Stage');
    ylabel(ax3, 'Flag rate');
    ylim(ax3, [0 1]);

    % Panel 4: PSD high-frequency diagnostic
    ax4 = subplot(2, 2, 4);
    hold(ax4, 'on');
    title(ax4, sprintf('FCz PSD; HF/LF ratio = %.4g', hf_ratio_FCz));

    if ~isempty(psd_f)
        plot(ax4, psd_f, 10 * log10(psd_mean), 'k', 'LineWidth', 1.6);
        xline(ax4, HF_BAND(1), 'r--', 'HandleVisibility', 'off');
        xline(ax4, HF_BAND(2), 'r--', 'HandleVisibility', 'off');
    else
        text(ax4, 0.5, 0.5, 'No valid epochs for PSD', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center');
    end

    xlabel(ax4, 'Frequency (Hz)');
    ylabel(ax4, 'Power/Frequency (dB/Hz)');
    xlim(ax4, [0 PSD_MAX_HZ]);
    grid(ax4, 'on');

    if save_figures
        save_figure(fig2, fullfile(figure_output_folder, sprintf('%s_v4_noise_QC_debugged.pdf', subj)));
    end
    if ~show_figures, close(fig2); end

    % ---------------------------------------------------------------------
    % Subject summary row
    % ---------------------------------------------------------------------
    this_summary = table( ...
        string(subj), ...
        p, ...
        string(epoch_file_used), ...
        string(hw_filter_label), ...
        EEGp.trials, ...
        n_rows, ...
        sum(sf.has_eeg_epoch), ...
        sum(sf.epoch_artifact_flag), ...
        mean(sf.epoch_artifact_flag(sf.has_eeg_epoch), 'omitnan'), ...
        hf_abs_FCz, ...
        hf_ratio_FCz, ...
        'VariableNames', { ...
            'subjID', ...
            'subj_num', ...
            'epoch_file_used', ...
            'hw_filter_label', ...
            'n_epochs_file', ...
            'n_behav_rows', ...
            'n_mapped_rows', ...
            'n_artifact_flagged_rows', ...
            'artifact_flag_rate', ...
            'FCz_HF_mean_PSD', ...
            'FCz_HF_to_LF_ratio'});

    qc_summary = [qc_summary; this_summary]; %#ok<AGROW>

    clear EEGp
end

% -------------------------------------------------------------------------
%% Save summary table
% -------------------------------------------------------------------------

summary_file = fullfile(figure_output_folder, 'KH_v4_ERP_noise_QC_summary_debugged.csv');
writetable(qc_summary, summary_file);
fprintf('\nSaved QC summary: %s\n', summary_file);

mat_summary_file = fullfile(figure_output_folder, 'KH_v4_ERP_noise_QC_summary_debugged.mat');
save(mat_summary_file, 'qc_summary');
fprintf('Saved QC summary MAT: %s\n', mat_summary_file);

% -------------------------------------------------------------------------
%% Optional group-level clean ERP plot
% -------------------------------------------------------------------------

if make_group_qc && ~isempty(qc_summary)

    figG = figure('Position', [100 100 1500 760], 'Visible', onoff(show_figures));
    sgtitle('KH v4 group QC: clean FCz subject means only', 'Interpreter', 'none');

    for s_i = 1:numel(stage_names)
        st = stage_names{s_i};

        for oc_i = 1:2
            oc_label = "incorrect";
            field_label = 'incorrect';
            if oc_i == 2
                oc_label = "correct";
                field_label = 'correct';
            end

            ax = subplot(2, 4, (oc_i - 1) * 4 + s_i);
            hold(ax, 'on');

            title(ax, sprintf('%s | %s', st, oc_label), 'Interpreter', 'none');
            xline(ax, 0, 'k:', 'HandleVisibility', 'off');
            yline(ax, 0, 'k:', 'HandleVisibility', 'off');

            for bi = 1:numel(btype_labels)
                bt = btype_labels{bi};
                dat = grand_clean_FCz.(st).(bt).(field_label);

                if isempty(dat), continue; end

                mn = mean(dat, 1, 'omitnan');
                se = std(dat, 0, 1, 'omitnan') ./ sqrt(size(dat, 1));

                c = line_color_for(bt, oc_i - 1);
                ls = '-';
                if strcmp(bt, 'P'), ls = '--'; end

                fill(ax, [tt fliplr(tt)], [mn + se fliplr(mn - se)], c, ...
                    'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
                    'HandleVisibility', 'off');

                plot(ax, tt, mn, 'Color', c, 'LineStyle', ls, 'LineWidth', 2, ...
                    'DisplayName', sprintf('%s N=%d', bt, size(dat, 1)));
            end

            set(ax, 'YDir', 'reverse');
            xlabel(ax, 'Time (ms)');
            ylabel(ax, 'FCz \muV');
            xlim(ax, ERP_plot_window);
            legend(ax, 'Box', 'off', 'FontSize', 8, 'Location', 'best');
        end
    end

    if save_figures
        save_figure(figG, fullfile(figure_output_folder, 'KH_group_clean_FCz_ERP_debugged.pdf'));
    end
    if ~show_figures, close(figG); end
end

fprintf('\nAll requested ERP/noise QC plots complete.\n');

% =============================================================================
%% LOCAL FUNCTIONS
% =============================================================================

function y = onoff(tf)
if tf
    y = 'on';
else
    y = 'off';
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

function trial2epoch_mapped = remap_trial2epoch_to_loaded_file(trial2epoch_raw, valid_ep, n_loaded_trials, used_trimmed_file)

trial2epoch_raw = trial2epoch_raw(:);
trial2epoch_mapped = nan(size(trial2epoch_raw));

if used_trimmed_file && ~isempty(valid_ep)
    valid_ep = valid_ep(:);

    % valid_ep(k) in the original untrimmed file is now epoch k in the
    % trimmed file. This is the important remapping step.
    max_old = max(valid_ep);
    old_to_new = nan(max_old, 1);
    old_to_new(valid_ep) = (1:numel(valid_ep))';

    ok = ~isnan(trial2epoch_raw) & ...
         trial2epoch_raw >= 1 & ...
         trial2epoch_raw <= numel(old_to_new) & ...
         ~isnan(old_to_new(round(trial2epoch_raw)));

    trial2epoch_mapped(ok) = old_to_new(round(trial2epoch_raw(ok)));

else
    % Untrimmed file, or no valid_ep available.
    trial2epoch_mapped = trial2epoch_raw;
end

% Final safety clamp against the actually loaded file.
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
        warning('epoch_artifact_flag_trimmed shorter than loaded EEG trials. Padding with false.');
    end

elseif isfield(S2, 'epoch_artifact_flag') && ~isempty(S2.epoch_artifact_flag)
    raw = logical(S2.epoch_artifact_flag(:));

    if ~isempty(valid_ep) && max(valid_ep) <= numel(raw)
        raw_trimmed = raw(valid_ep);
        if numel(raw_trimmed) >= n_loaded_trials
            flags = raw_trimmed(1:n_loaded_trials);
        else
            flags(1:numel(raw_trimmed)) = raw_trimmed;
        end
    elseif numel(raw) >= n_loaded_trials
        flags = raw(1:n_loaded_trials);
    else
        flags(1:numel(raw)) = raw;
        warning('epoch_artifact_flag shorter than loaded EEG trials. Padding with false.');
    end
else
    fprintf('No epoch artifact flag found in cache. Using all false.\n');
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
y(sv == "0" | sv == "false" | sv == "incorrect" | sv == "incorr" | sv == "wrong") = 0;

if any(isnan(y))
    tmp = str2double(sv);
    fillable = isnan(y) & ~isnan(tmp);
    y(fillable) = tmp(fillable);
end

end

function c = line_color_for(bt, correct_val)

% Fixed colours for visual interpretability:
%   D incorrect, D correct, P incorrect, P correct
if strcmp(char(bt), 'D') && correct_val == 0
    c = [0.20 0.35 0.80];
elseif strcmp(char(bt), 'D') && correct_val == 1
    c = [0.05 0.55 0.85];
elseif strcmp(char(bt), 'P') && correct_val == 0
    c = [0.80 0.30 0.20];
else
    c = [0.85 0.55 0.10];
end

end

function [f_keep, pxx_mean_keep, hf_ratio, hf_abs] = compute_subject_psd_hf_ratio(EEGp, ch_idx, eps_use, lf_band, hf_band, max_hz)

fs = EEGp.srate;

x = squeeze(double(EEGp.data(ch_idx, :, eps_use)))';

if isvector(x)
    x = x(:)';
end

x = x - mean(x, 2, 'omitnan');

% pwelch expects observations by column, so transpose to time x trials.
win = min(round(fs), size(x, 2));
if win < 8
    f_keep = [];
    pxx_mean_keep = [];
    hf_ratio = NaN;
    hf_abs = NaN;
    return;
end

noverlap = floor(win / 2);
nfft = max(512, 2^nextpow2(win));

[pxx, f] = pwelch(x', win, noverlap, nfft, fs);
pxx_mean = mean(pxx, 2, 'omitnan');

keep = f <= max_hz;
f_keep = f(keep);
pxx_mean_keep = pxx_mean(keep);

lf_mask = f >= lf_band(1) & f <= lf_band(2);
hf_mask = f >= hf_band(1) & f <= hf_band(2);

lf_abs = mean(pxx_mean(lf_mask), 'omitnan');
hf_abs = mean(pxx_mean(hf_mask), 'omitnan');

hf_ratio = hf_abs / lf_abs;

end

function save_figure(fig, fname)

[folder, ~, ~] = fileparts(fname);
if ~exist(folder, 'dir')
    mkdir(folder);
end

try
    exportgraphics(fig, fname, 'ContentType', 'vector');
catch
    warning('exportgraphics failed for %s. Falling back to saveas.', fname);
    saveas(fig, fname);
end

fprintf('Saved figure: %s\n', fname);

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