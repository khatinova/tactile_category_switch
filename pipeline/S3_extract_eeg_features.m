% =============================================================================
% S3_extract_eeg_features.m
%
% PIPELINE STEP 3 of 7 — (re)build the per-trial EEG feature table so that its
% schema MATCHES the S2_RR output exactly, allowing KH + RR to merge in S4.
%
% WHY THIS SCRIPT EXISTS
% ----------------------
% S2 / S2_RR export the per-trial spine (behaviour + epoch + stage) and the
% outcome/theta/phase epochs. S3 reloads each subject's epochs ONCE and
% recomputes the EEG features over the PREFRONTAL CHANNEL CLUSTER (not single
% FCz), writing the SAME columns S2_RR produces:
%
%   PER-TRIAL (stored in all_trials_table):
%     N2_amp/_lat/_norm            frontal negativity (min in N2 window)
%     prefrontal_mean_amp/_norm    prefrontal fixed-window (FRN-win) mean
%     prefrontal_waveform          per-trial prefrontal-cluster waveform (cell)
%     P300_amp/_norm/_peak_lat, P300_waveform
%     Theta_amp, Theta_waveform
%     PLV_fp/_fs(_pairwise)        sliding-window fronto-parietal / -sensory PLV
%   PER-STAGE difference waves (frn_rewp_stage_table):
%     FRN_amp  = mean over FRN window  of (incorrect - correct)  [negative]
%     RewP_amp = mean over RewP window of (correct - incorrect)  [positive]
%   built from prefrontal_waveform via kh_compute_frn_rewp_by_stage.m.
%
% Cohort-aware twin of S2_RR PART 2: set COHORT = 'KH' (default) to recompute
% the KH features with the prefrontal cluster; 'RR' is also supported.
%
% INPUT  : group_feature_table_combined_<COHORT>.mat (all_trials_table spine
%          with epoch + stage) + outcome/theta/phase .set (from S2 / S2_RR).
% OUTPUT : group_feature_table_combined_<COHORT>.mat (features now on the
%          prefrontal cluster) + group_table_all_trials/stage/frn _<TAG>
%          + grand-average figures.
% =============================================================================

clear; close all; clc;
addpath(genpath(fileparts(mfilename('fullpath'))));   % pipeline utils on path

% -------------------------------------------------------------------------
%% COHORT SELECTION  (channel clusters identical to S2_RR)
% -------------------------------------------------------------------------
COHORT = 'KH';                 % 'KH' (default) or 'RR'
tag    = upper(COHORT);

switch tag
    case 'KH'
        cz_label            = 'Cz';
        prefrontal_channels = {'FCz','Cz','Fz','FC1','FC2'};
        par_channels        = {'Pz','P1','P2'};
        acc_channels        = {'FCz','Fz','AFz','F1','F2'};
        som_channels        = {'C1','C5'};   % e.g. {'C3','C4','CP3','CP1','C5','CP5'}
    case 'RR'
        cz_label            = 'E129';
        prefrontal_channels = {'E7','E6','E5','E12','E106','E129'};
        par_channels        = {'E61','E62','E78'};
        acc_channels        = {'E11','E19','E4','E5','E6','E12'};
        som_channels        = {'E30','E36','E37','E41'};
    otherwise
        error('Unknown COHORT "%s" (use ''KH'' or ''RR'').', COHORT);
end


% -------------------------------------------------------------------------
%% PATHS
% -------------------------------------------------------------------------
remote = 0;
if remote == 1
    base_path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch';
else
    base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch';
end

if strcmpi(COHORT,'RR')
    cohort_results = fullfile(base_path, 'Salient mod switch RR', 'Results', 'EEG analysis');
else
    cohort_results = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis');
end

epoch_file_folder    = fullfile(cohort_results, 'Epoched_data_noisefiltering');
feature_folder       = fullfile(cohort_results, 'Outcome_feature_tables_v4_merged');
figure_output_folder = fullfile(cohort_results, 'Figures', 'RQ_analysis');
if ~exist(figure_output_folder,'dir'), mkdir(figure_output_folder); end
if ~exist(feature_folder,'dir'),       mkdir(feature_folder);       end

% -------------------------------------------------------------------------
%% SETTINGS  (windows identical to S2_RR)
% -------------------------------------------------------------------------
ERP_plot_window = [-200 800];
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

% -------------------------------------------------------------------------
%% LOAD THE PER-TRIAL SPINE (from S2 / S2_RR)
% -------------------------------------------------------------------------
combined_file = fullfile(feature_folder, sprintf('group_feature_table_combined_%s.mat', tag));
if ~exist(combined_file,'file')
    error(['S3: spine not found:\n  %s\nRun S2 (COHORT=%s) first so the ' ...
           'per-trial table (with epoch + stage) exists.'], combined_file, tag);
end
S = load(combined_file);
if isfield(S,'all_trials_table')
    group_table = S.all_trials_table;
elseif isfield(S,'group_table')
    group_table = S.group_table;
else
    error('S3: %s lacks all_trials_table/group_table.', combined_file);
end
group_table = kh_subject_id('standardise', group_table);   % canonical subj_id/subj/cohort
group_table.subj_id_s = string(group_table.subj_id);

% Robust accessors the spine must provide.
if ~ismember('correct_num', group_table.Properties.VariableNames)
    group_table.correct_num = local_correct_to_numeric(group_table.correct);
end
if ~ismember('false_fb', group_table.Properties.VariableNames)
    group_table.false_fb = false(height(group_table),1);
end
if ~ismember('epoch', group_table.Properties.VariableNames)
    error('S3: spine lacks an "epoch" column; cannot align EEG to trials.');
end
if ~ismember('trial_continuous', group_table.Properties.VariableNames)
    % needed only to order trials inside the sliding-window PLV buckets
    group_table.trial_continuous = (1:height(group_table))';
end


% -------------------------------------------------------------------------
%% (RE)INITIALISE FEATURE COLUMNS  (same names/shape as S2_RR)
% -------------------------------------------------------------------------
n_rows = height(group_table);
group_table.baseline_rms         = nan(n_rows,1);
group_table.N2_amp               = nan(n_rows,1);
group_table.N2_lat               = nan(n_rows,1);
group_table.N2_norm              = nan(n_rows,1);
group_table.prefrontal_mean_amp  = nan(n_rows,1);
group_table.prefrontal_mean_norm = nan(n_rows,1);
group_table.P300_amp             = nan(n_rows,1);
group_table.P300_peak_lat        = nan(n_rows,1);
group_table.P300_norm            = nan(n_rows,1);
group_table.Theta_amp            = nan(n_rows,1);
group_table.PLV_fp               = nan(n_rows,1);
group_table.PLV_fs               = nan(n_rows,1);
group_table.PLV_fp_pairwise      = nan(n_rows,1);
group_table.PLV_fs_pairwise      = nan(n_rows,1);
group_table.prefrontal_waveform  = repmat({[]}, n_rows, 1);
group_table.P300_waveform        = repmat({[]}, n_rows, 1);
group_table.Theta_waveform       = repmat({[]}, n_rows, 1);

t_ax = [];
% grand-average storage for the normalisation figure (no .set reload later)
raw_correct = []; raw_incorrect = []; norm_correct = []; norm_incorrect = [];

% -------------------------------------------------------------------------
%% MASTER LOOP — load each subject's epochs ONCE, recompute features
% -------------------------------------------------------------------------
subj_list = unique(group_table.subj_id_s);
fprintf('S3 (%s): recomputing prefrontal-cluster features for %d subjects.\n', ...
    tag, numel(subj_list));

for si = 1:numel(subj_list)

    subj      = subj_list(si);                     % e.g. "Ox03" / "Nc07"
    subj_rows = group_table.subj_id_s == subj;
    fprintf('  %s (%d trials)...\n', subj, sum(subj_rows));

    % --- locate + load the broadband outcome epochs (trimmed preferred) ---
    EEGp = load_first_existing_set(epoch_file_folder, { ...
        sprintf('%s_outcome_trimmed.set', subj), ...
        sprintf('%s_outcome.set', subj)});
    if isempty(EEGp)
        warning('%s: no outcome .set found. Skipping.', subj); continue;
    end
    if isempty(t_ax), t_ax = EEGp.times; end

    EEGp_theta = load_first_existing_set(epoch_file_folder, { ...
        sprintf('%s_outcome_theta_trimmed.set', subj), ...
        sprintf('%s_outcome_theta.set', subj)});
    EEGp_phase = load_first_existing_set(epoch_file_folder, { ...
        sprintf('%s_outcome_phase_trimmed.set', subj), ...
        sprintf('%s_outcome_phase.set', subj)});

    % --- channel indices (prefrontal CLUSTER, parietal, ACC, somatosensory) ---
    chan_labels_lower = lower(string({EEGp.chanlocs.labels}));
    prefrontal_idx = find(ismember(chan_labels_lower, lower(string(prefrontal_channels))));
    cz_idx         = find(chan_labels_lower == lower(string(cz_label)), 1);
    par_idx        = find(ismember(chan_labels_lower, lower(string(par_channels))));
    acc_idx        = find(ismember(chan_labels_lower, lower(string(acc_channels))));
    som_idx        = find(ismember(chan_labels_lower, lower(string(som_channels))));

    if isempty(prefrontal_idx)
        warning('%s: no prefrontal-cluster channels found. Skipping.', subj); continue;
    end

    % --- time masks ---
    bl_mask   = EEGp.times >= rm_baseline(1)  & EEGp.times <= rm_baseline(2);
    n2_mask   = EEGp.times >= N2_win(1)       & EEGp.times <= N2_win(2);
    frn_mask  = EEGp.times >= FRN_win(1)      & EEGp.times <= FRN_win(2);
    p300_mask = EEGp.times >= P300_win(1)     & EEGp.times <= P300_win(2);
    th_mask   = EEGp.times >= Theta_win(1)    & EEGp.times <= Theta_win(2);
    plv_mask  = EEGp.times >= PLV_win(1)      & EEGp.times <= PLV_win(2);
    plv_bl    = EEGp.times >= PLV_baseline(1) & EEGp.times <= PLV_baseline(2);

    % --- baseline RMS (single Cz channel, as in S2_RR) for normalisation ---
    if ~isempty(cz_idx)
        bl_data   = squeeze(double(EEGp.data(cz_idx, bl_mask, :)));
        bline_rms = rms(bl_data(:), 'omitnan');
    else
        bl_data   = squeeze(mean(double(EEGp.data(prefrontal_idx, bl_mask, :)),1,'omitnan'));
        bline_rms = rms(bl_data(:), 'omitnan');
    end
    group_table.baseline_rms(subj_rows) = bline_rms;


    % --- per-trial feature extraction ---
    row_indices = find(subj_rows);
    for ri = 1:numel(row_indices)
        r  = row_indices(ri);
        ep = group_table.epoch(r);
        if isnan(ep) || ep < 1 || ep > EEGp.trials, continue; end
        ep = round(ep);

        % Prefrontal-cluster waveform (baseline-corrected), N2 + FRN-win mean.
        sig = mean(double(EEGp.data(prefrontal_idx, :, ep)), 1, 'omitnan');
        sig = sig - mean(sig(bl_mask), 'omitnan');
        group_table.prefrontal_waveform{r} = sig;

        win_vals = sig(n2_mask);  win_t = EEGp.times(n2_mask);
        if any(~isnan(win_vals))
            [pk, ix] = min(win_vals, [], 'omitnan');
            group_table.N2_amp(r) = pk;
            group_table.N2_lat(r) = win_t(ix);
            if bline_rms > 0, group_table.N2_norm(r) = pk / bline_rms; end
        end

        frn_vals = sig(frn_mask);
        if any(~isnan(frn_vals))
            group_table.prefrontal_mean_amp(r) = mean(frn_vals, 'omitnan');
            if bline_rms > 0
                group_table.prefrontal_mean_norm(r) = group_table.prefrontal_mean_amp(r) / bline_rms;
            end
        end

        % Parietal P300.
        if ~isempty(par_idx)
            sig_p = mean(double(EEGp.data(par_idx, :, ep)), 1, 'omitnan');
            sig_p = sig_p - mean(sig_p(bl_mask), 'omitnan');
            group_table.P300_waveform{r} = sig_p;

            win_vals = sig_p(p300_mask);  win_t = EEGp.times(p300_mask);
            if any(~isnan(win_vals))
                [pk, ix] = max(win_vals, [], 'omitnan');
                group_table.P300_amp(r)      = pk;
                group_table.P300_peak_lat(r) = win_t(ix);
                if bline_rms > 0, group_table.P300_norm(r) = pk / bline_rms; end
            end
        end

        % Frontal theta amplitude envelope (ACC cluster).
        if ~isempty(EEGp_theta) && ~isempty(acc_idx) && ep <= EEGp_theta.trials
            sig_th = mean(double(EEGp_theta.data(acc_idx, :, ep)), 1, 'omitnan');
            env    = abs(hilbert(sig_th));
            env    = env - mean(env(bl_mask), 'omitnan');
            group_table.Theta_amp(r)      = mean(env(th_mask), 'omitnan');
            group_table.Theta_waveform{r} = env;
        end
    end


    % --- sliding-window PLV (fronto-parietal / fronto-sensory), as in S2_RR ---
    if ~isempty(EEGp_phase) && ~isempty(acc_idx) && EEGp_phase.trials > 0
        eprows   = row_indices;
        valid_ep = group_table.epoch(eprows) >= 1 & ...
                   group_table.epoch(eprows) <= EEGp_phase.trials & ...
                   ~isnan(group_table.epoch(eprows));

        stg = string(group_table.stage(eprows));
        bty = string(group_table.block_type(eprows));
        cor = string(group_table.correct_num(eprows));
        ok  = valid_ep & ~ismissing(group_table.stage(eprows)) & ...
              ~ismissing(group_table.block_type(eprows)) & ~isnan(group_table.correct_num(eprows));

        buckets = unique([stg(ok), bty(ok), cor(ok)], 'rows');
        for ub = 1:size(buckets,1)
            bmask = ok & stg==buckets(ub,1) & bty==buckets(ub,2) & cor==buckets(ub,3);
            brows = eprows(bmask);
            if numel(brows) < MIN_TRIALS_PLV_WINDOW, continue; end

            [~, ord] = sort(group_table.trial_continuous(brows));
            brows = brows(ord);
            nb = numel(brows);

            for bi2 = 1:nb
                lo = max(1, bi2-PLV_WINDOW_HALF);  hi = min(nb, bi2+PLV_WINDOW_HALF);
                wrows = brows(lo:hi);
                eps_w = round(group_table.epoch(wrows));
                eps_w = eps_w(~isnan(eps_w) & eps_w>=1 & eps_w<=EEGp_phase.trials);
                if numel(eps_w) < MIN_TRIALS_PLV_WINDOW, continue; end
                cr = brows(bi2);

                phi_ref = squeeze(angle(mean(exp(1i*double( ...
                    EEGp_phase.data(acc_idx,:,eps_w))),1,'omitnan')))';

                if ~isempty(par_idx)
                    phi_t  = squeeze(angle(mean(exp(1i*double( ...
                        EEGp_phase.data(par_idx,:,eps_w))),1,'omitnan')))';
                    pts    = abs(mean(exp(1i*(phi_ref-phi_t)),1,'omitnan'));
                    pts    = pts - mean(pts(plv_bl),'omitnan');
                    group_table.PLV_fp(cr)          = mean(pts(plv_mask),'omitnan');
                    group_table.PLV_fp_pairwise(cr) = group_table.PLV_fp(cr);
                end
                if ~isempty(som_idx)
                    phi_t  = squeeze(angle(mean(exp(1i*double( ...
                        EEGp_phase.data(som_idx,:,eps_w))),1,'omitnan')))';
                    pts    = abs(mean(exp(1i*(phi_ref-phi_t)),1,'omitnan'));
                    pts    = pts - mean(pts(plv_bl),'omitnan');
                    group_table.PLV_fs(cr)          = mean(pts(plv_mask),'omitnan');
                    group_table.PLV_fs_pairwise(cr) = group_table.PLV_fs(cr);
                end
            end
        end
    end

    % --- store subject grand-average ERP (correct/incorrect) for the figure ---
    if bline_rms > 0
        pf = group_table.prefrontal_waveform(subj_rows & ~group_table.false_fb & ...
             group_table.correct_num==1);
        erp_c = mean_waveform(pf);
        pf = group_table.prefrontal_waveform(subj_rows & ~group_table.false_fb & ...
             group_table.correct_num==0);
        erp_i = mean_waveform(pf);
        if ~isempty(erp_c), raw_correct(end+1,:)=erp_c; norm_correct(end+1,:)=erp_c/bline_rms; end %#ok<AGROW>
        if ~isempty(erp_i), raw_incorrect(end+1,:)=erp_i; norm_incorrect(end+1,:)=erp_i/bline_rms; end %#ok<AGROW>
    end

    clear EEGp EEGp_theta EEGp_phase
end


% -------------------------------------------------------------------------
%% STEP 2: WITHIN-SUBJECT Z-SCORES (for the LME models in S6/S7)
% -------------------------------------------------------------------------
features_to_zscore = {'N2_amp','N2_norm','prefrontal_mean_amp','prefrontal_mean_norm', ...
    'P300_amp','P300_norm','Theta_amp','PLV_fp','PLV_fs','PLV_fp_pairwise','PLV_fs_pairwise'};
for f = 1:numel(features_to_zscore)
    fn = features_to_zscore{f};
    if ~ismember(fn, group_table.Properties.VariableNames), continue; end
    fn_z = [fn '_z'];
    group_table.(fn_z) = nan(height(group_table),1);
    for si = 1:numel(subj_list)
        m  = group_table.subj_id_s == subj_list(si);
        v  = group_table.(fn)(m);
        sd = std(v,'omitnan');
        if sd > 0, group_table.(fn_z)(m) = (v - mean(v,'omitnan')) / sd; end
    end
end

% -------------------------------------------------------------------------
%% PER-STAGE FRN / RewP difference waves (from the prefrontal cluster)
%   FRN  = mean over FRN window  of (incorrect - correct)  [negative]
%   RewP = mean over RewP window of (correct - incorrect)  [positive]
% -------------------------------------------------------------------------
frn_rewp_opts = struct('wave_col','prefrontal_waveform', ...
                       'FRN_win',[250 300], 'RewP_win',[250 350], ...
                       'stages',{stage_names}, 'block_types',{BTYPE_LABELS});
frn_rewp_stage_table = kh_compute_frn_rewp_by_stage(group_table, t_ax, frn_rewp_opts);

% Per-stage summary of the single-trial measures (matches S2_RR output set).
stage_feature_table = build_stage_summary(group_table, stage_names, BTYPE_LABELS);

% -------------------------------------------------------------------------
%% SAVE (cohort-tagged; overwrites the combined file with prefrontal features)
% -------------------------------------------------------------------------
all_trials_table = group_table;
save(fullfile(feature_folder, sprintf('group_table_all_trials_%s.mat', tag)), ...
     'all_trials_table','t_ax','-v7.3');
save(fullfile(feature_folder, sprintf('group_stage_table_features_%s.mat', tag)), ...
     'stage_feature_table','-v7.3');
save(fullfile(feature_folder, sprintf('frn_rewp_by_stage_%s.mat', tag)), ...
     'frn_rewp_stage_table','-v7.3');
save(fullfile(feature_folder, sprintf('group_feature_table_combined_%s.mat', tag)), ...
     'all_trials_table','stage_feature_table','frn_rewp_stage_table','t_ax','-v7.3');
fprintf('S3 (%s): saved feature tables to\n  %s\n', tag, feature_folder);


% -------------------------------------------------------------------------
%% FIGURE 1: raw vs baseline-RMS normalised grand-average ERP
% -------------------------------------------------------------------------
if ~isempty(raw_correct) || ~isempty(raw_incorrect)
    in_w = t_ax >= ERP_plot_window(1) & t_ax <= ERP_plot_window(2);
    fig = figure('Position',[50 50 1200 500]);
    sgtitle(sprintf('%s: prefrontal ERP (raw vs baseline-RMS normalised)', tag));

    subplot(1,2,1); hold on; title('Raw (\muV)');
    shaded_erp(t_ax(in_w), pick(raw_correct,in_w),   [0.20 0.60 0.20], 'Correct');
    shaded_erp(t_ax(in_w), pick(raw_incorrect,in_w), [0.80 0.20 0.20], 'Incorrect');
    xline(0,'k--'); yline(0,'k:'); set(gca,'YDir','reverse');
    xlabel('Time (ms)'); ylabel('Amplitude (\muV)'); legend('Box','off');

    subplot(1,2,2); hold on; title('Normalised (\muV / baseline RMS)');
    shaded_erp(t_ax(in_w), pick(norm_correct,in_w),   [0.20 0.60 0.20], 'Correct');
    shaded_erp(t_ax(in_w), pick(norm_incorrect,in_w), [0.80 0.20 0.20], 'Incorrect');
    xline(0,'k--'); yline(0,'k:'); set(gca,'YDir','reverse');
    xlabel('Time (ms)'); ylabel('Amplitude (normalised)'); legend('Box','off');

    apply_fig_style(fig);
    exportgraphics(fig, fullfile(figure_output_folder, ...
        sprintf('Normalisation_comparison_%s.pdf', tag)), 'ContentType','vector');
end

% -------------------------------------------------------------------------
%% FIGURE 2: FRN / RewP difference-wave grand averages, by block type
% -------------------------------------------------------------------------
if ~isempty(frn_rewp_stage_table) && ismember('diff_wave', frn_rewp_stage_table.Properties.VariableNames)
    fig = figure('Position',[60 60 1100 460]);
    sgtitle(sprintf('%s: FRN / RewP difference waves (incorrect - correct), by block type', tag));
    FRNw = [250 300]; RewPw = [250 350]; bts = {'D','P'};
    for bi = 1:2
        ax = subplot(1,2,bi); hold(ax,'on'); title(ax, sprintf('%s blocks', bts{bi}));
        sel = string(frn_rewp_stage_table.block_type) == bts{bi};
        dw  = frn_rewp_stage_table.diff_wave(sel); dw = dw(~cellfun(@isempty,dw));
        if isempty(dw); continue; end
        M  = cell2mat(cellfun(@(x) x(:)', dw, 'UniformOutput', false));
        mn = mean(M,1,'omitnan'); se = std(M,0,1,'omitnan')/sqrt(size(M,1));
        yl = [min(mn)-0.5, max(mn)+0.5];
        patch(ax,[FRNw fliplr(FRNw)],[yl(1) yl(1) yl(2) yl(2)],[0.85 0.9 1],'EdgeColor','none','FaceAlpha',0.4);
        patch(ax,[RewPw fliplr(RewPw)],[yl(1) yl(1) yl(2) yl(2)],[1 0.9 0.85],'EdgeColor','none','FaceAlpha',0.3);
        fill(ax,[t_ax fliplr(t_ax)],[mn+se fliplr(mn-se)],[0.2 0.2 0.2],'FaceAlpha',0.15,'EdgeColor','none');
        plot(ax, t_ax, mn, 'k', 'LineWidth', 2);
        yline(ax,0,'k:'); xline(ax,0,'k:');
        xlabel(ax,'Time (ms)'); ylabel(ax,'incorrect - correct (\muV)');
        xlim(ax, ERP_plot_window); set(ax,'YDir','reverse');
    end
    apply_fig_style(fig);
    exportgraphics(fig, fullfile(figure_output_folder, ...
        sprintf('FRN_RewP_difference_waves_%s.pdf', tag)), 'ContentType','vector');
end

fprintf('S3 (%s) complete.\n', tag);


% =========================================================================
%% LOCAL FUNCTIONS
% =========================================================================
function EEG_loaded = load_first_existing_set(folder, candidates)
% Load the first .set in CANDIDATES that exists in FOLDER ([] if none).
EEG_loaded = [];
for ci = 1:numel(candidates)
    if exist(fullfile(folder, candidates{ci}), 'file')
        EEG_loaded = pop_loadset(candidates{ci}, folder);
        return;
    end
end
end

function m = mean_waveform(wave_cells)
% Grand-average a cell column of equal-length per-trial waveforms ([] if none).
wave_cells = wave_cells(~cellfun(@isempty, wave_cells));
if isempty(wave_cells), m = []; return; end
M = cell2mat(cellfun(@(x) x(:)', wave_cells, 'UniformOutput', false));
m = mean(M, 1, 'omitnan');
end

function out = pick(M, cols)
% Safely column-subset a (possibly empty) [subjects x time] matrix.
if isempty(M), out = M; else, out = M(:, cols); end
end

function shaded_erp(t, data, clr, lbl)
% Mean +/- SEM ribbon for a [subjects x time] matrix.
if isempty(data), return; end
mn = mean(data,1,'omitnan');
se = std(data,0,1,'omitnan')/sqrt(size(data,1));
fill([t,fliplr(t)],[mn+se,fliplr(mn-se)],clr,'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
plot(t,mn,'Color',clr,'LineWidth',2.5,'DisplayName',sprintf('%s (n=%d)',lbl,size(data,1)));
end

function y = local_correct_to_numeric(v)
% Coerce a 'correct' column (numeric/logical/categorical/string) to 0/1.
if isnumeric(v) || islogical(v), y = double(v(:)); return; end
sv = lower(strtrim(string(v)));  y = nan(numel(sv),1);
y(sv=="1"|sv=="true"|sv=="correct"|sv=="corr") = 1;
y(sv=="0"|sv=="false"|sv=="incorrect"|sv=="incorr"|sv=="wrong") = 0;
tmp = str2double(sv); fillable = isnan(y) & ~isnan(tmp); y(fillable) = tmp(fillable);
end

function Tout = build_stage_summary(T, stage_names, btype_labels)
% Per subj_id x block_type x stage x outcome mean of the single-trial measures
% that S3 computes (true-feedback trials only). Matches the S2_RR stage table
% for the measures that exist here.
meas = {'N2_amp','N2_lat','N2_norm','prefrontal_mean_amp','prefrontal_mean_norm', ...
        'P300_amp','P300_peak_lat','P300_norm','Theta_amp', ...
        'PLV_fp','PLV_fs','PLV_fp_pairwise','PLV_fs_pairwise'};
meas = meas(ismember(meas, T.Properties.VariableNames));
rows = {};
for s = 1:numel(stage_names)
    for b = 1:numel(btype_labels)
        for cv = 0:1
            m = string(T.stage)==stage_names{s} & string(T.block_type)==btype_labels{b} & ...
                T.correct_num==cv & ~T.false_fb;
            if ~any(m), continue; end
            row = table();
            row.subj_id     = T.subj_id(find(m,1));
            if ismember('subj',T.Properties.VariableNames),   row.subj   = T.subj(find(m,1));   end
            if ismember('cohort',T.Properties.VariableNames), row.cohort = T.cohort(find(m,1)); end
            row.stage       = categorical(stage_names(s), stage_names);
            row.block_type  = categorical(btype_labels(b), btype_labels);
            row.correct_num = cv;
            row.n_trials    = sum(m);
            for mm = 1:numel(meas)
                row.([meas{mm} '_mean']) = mean(T.(meas{mm})(m),'omitnan');
            end
            rows{end+1,1} = row; %#ok<AGROW>
        end
    end
end
if isempty(rows), Tout = table(); return; end
Tout = rows{1};
for r = 2:numel(rows)
    Tout = [Tout; rows{r}]; %#ok<AGROW>
end
end
