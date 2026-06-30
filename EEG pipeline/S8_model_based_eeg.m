% =============================================================================
% S8_model_based_eeg.m  — Model-based single-trial EEG analysis
%
% PIPELINE STEP 8 (optional, after S7)
%
% Uses computational model latents (Nassar & RW, from S5) as trial-level
% regressors predicting single-trial EEG features. Addresses four questions
% motivated by the peer-reviewed literature:
%
%   MRQ1 — Does the P300 track surprise (ω) and predict learning rate
%           adjustment? Does this differ for true vs false feedback?
%           (Nassar et al., 2019 eLife; Nassar et al., 2022 J Neurosci)
%
%   MRQ2 — Does the FRN (prefrontal negativity) scale with signed or
%           unsigned prediction error? Is it attenuated for false feedback
%           (feedback discounting)?
%           (Chase et al., 2011 JOCN; Worthy et al., 2018 CBB;
%            Talmi et al., 2013 J Neurosci)
%
%   MRQ3 — Does frontal theta power scale with unsigned PE (surprise)
%           and predict next-trial behavioural adjustment?
%           (Cavanagh et al., 2010 J Neurosci; Reteig et al., 2020)
%
%   MRQ4 — Does fronto-parietal PLV mediate the PE → behaviour link?
%           (theta-band functional coupling as the mechanism for
%            communicating need-for-adjustment signals)
%
% INPUT:  group_feature_table_combined.mat (from S4, contains EEG features +
%         model latents merged in via S5 → S4 pipeline)
% OUTPUT: Poster-quality figures + manuscript_stats_S8.txt in figure_output_folder
%
% REQUIREMENTS: Statistics & Machine Learning Toolbox (fitlme, fitglme)
%
% REFERENCES:
%   Nassar et al. (2019) eLife 8:e46975 — P300 as context-dependent surprise
%   Nassar et al. (2022) J Neurosci 42:2524 — bidirectional P300 learning signal
%   Chase et al. (2011) JOCN 23:936 — FRN = PE, P300 = behavioural adjustment
%   Cavanagh et al. (2010) J Neurosci 30:3051 — theta tracks PE magnitude
%   Worthy et al. (2018) Comput Brain Behav 1:1 — feedback discounting
%   Kirsch et al. (2021) Commun Biol 4:910 — separating overlapping PEs
% =============================================================================
close all; clc;

addpath(genpath(fileparts(mfilename('fullpath'))));

%% ── 0. PATHS & CONFIGURATION ───────────────────────────────────────────────
base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch';
saved_tables_folder = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis', 'Outcome_feature_tables_v4_merged');
figure_output_folder = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis', 'Figures', 'S8_model_based_EEG');

if ~exist(figure_output_folder, 'dir'), mkdir(figure_output_folder); end

% ── Load data ────────────────────────────────────────────────────────────────
if exist('gt', 'var') && istable(gt)
    fprintf('Using gt already in workspace (%d rows).\n', height(gt));
elseif exist('group_table', 'var') && istable(group_table)
    gt = group_table;
else
    load(fullfile(saved_tables_folder, 'group_feature_table_combined.mat'), 'group_table');
    gt = group_table;
end

% ── Poster colour palette ────────────────────────────────────────────────────
CLR_D    = [0.15 0.45 0.70];   % deterministic blocks — blue
CLR_P    = [0.80 0.30 0.10];   % probabilistic blocks — orange
CLR_TRUE = [0.20 0.60 0.30];   % true feedback — green
CLR_FALSE= [0.75 0.20 0.55];   % false feedback — magenta
CLR_GREY = [0.55 0.55 0.55];

fprintf('\n');
fprintf('================================================================\n');
fprintf('  S8: MODEL-BASED SINGLE-TRIAL EEG ANALYSIS\n');
fprintf('  Output: %s\n', figure_output_folder);
fprintf('================================================================\n\n');

%% ── 1. ENSURE REQUIRED VARIABLES ────────────────────────────────────────────
% Nassar latents (from S5 → S4 pipeline)
required_nassar = {'PE_nassar','PE_unsigned','omega','alpha_nassar', ...
                   'certainty','surprise','theta_nassar'};
required_rw     = {'RW_PE','RW_PE_unsigned','RW_value'};
required_eeg    = {'prefrontal_neg_peak_norm','P300_norm','Theta_amp', ...
                   'PLV_fp'};

missing = {};
for c = [required_nassar, required_rw, required_eeg]
    if ~ismember(c{1}, gt.Properties.VariableNames)
        missing{end+1} = c{1}; %#ok<SAGROW>
    end
end
if ~isempty(missing)
    warning('Missing columns: %s\nSome analyses will be skipped.', strjoin(missing, ', '));
end

% Ensure categoricals
gt.subj_id    = categorical(gt.subj_id);
gt.block_type = categorical(gt.block_type);
gt.stage      = categorical(gt.stage, {'LN','LE','RN','RE'}, 'Ordinal', true);

% Within-subject z-scoring of model and EEG variables
vars_to_z = {'PE_nassar','PE_unsigned','omega','alpha_nassar','certainty', ...
             'surprise','theta_nassar','RW_PE','RW_PE_unsigned','RW_value', ...
             'prefrontal_neg_peak_norm','P300_norm','Theta_amp','PLV_fp'};
subj_list = unique(gt.subj_id);

for f = 1:numel(vars_to_z)
    fn   = vars_to_z{f};
    fn_z = [fn '_z'];
    if ~ismember(fn, gt.Properties.VariableNames), continue; end
    if ismember(fn_z, gt.Properties.VariableNames), continue; end  % already done
    gt.(fn_z) = nan(height(gt), 1);
    for si = 1:numel(subj_list)
        mask = gt.subj_id == subj_list(si);
        vals = gt.(fn)(mask);
        mn = mean(vals, 'omitnan');
        sd = std(vals, 'omitnan');
        if sd > 0
            gt.(fn_z)(mask) = (vals - mn) / sd;
        end
    end
end

% Compute next-trial accuracy (within block boundaries)
if ~ismember('next_correct', gt.Properties.VariableNames)
    gt.next_correct = nan(height(gt), 1);
    blk_col = 'block';
    if ~ismember(blk_col, gt.Properties.VariableNames)
        if ismember('block_number', gt.Properties.VariableNames)
            blk_col = 'block_number';
        end
    end
    for si = 1:numel(subj_list)
        for b = 1:10
            mask = gt.subj_id == subj_list(si) & gt.(blk_col) == b;
            idx  = find(mask);
            if numel(idx) < 2, continue; end
            gt.next_correct(idx(1:end-1)) = gt.correct(idx(2:end));
        end
    end
end

% false_fb as categorical for LME
if ismember('false_fb', gt.Properties.VariableNames)
    gt.false_fb_cat = categorical(double(gt.false_fb), [0 1], {'True','False'});
end

fprintf('Data ready: %d trials, %d subjects.\n', height(gt), numel(subj_list));
fprintf('Model vars present: PE_nassar=%d, RW_PE=%d, omega=%d, surprise=%d\n', ...
    ismember('PE_nassar', gt.Properties.VariableNames), ...
    ismember('RW_PE', gt.Properties.VariableNames), ...
    ismember('omega', gt.Properties.VariableNames), ...
    ismember('surprise', gt.Properties.VariableNames));


%% =========================================================================
%  MRQ1: P300 AS A SURPRISE SIGNAL — MODULATED BY FEEDBACK VALIDITY
%  =========================================================================
%  Hypothesis (Nassar 2019, 2022): P300 amplitude tracks change-point
%  probability (ω). In change-point contexts (genuine reversals), P300
%  predicts learning rate increase. In noise contexts (false feedback in
%  P blocks), P300 should NOT predict increased learning.
%
%  Key test: P300 × false_fb interaction on alpha_nassar (next-trial LR)
% =========================================================================
fprintf('\n=== MRQ1: P300 ~ surprise × feedback validity ===\n');

if all(ismember({'P300_norm_z','omega_z','surprise_z','false_fb'}, gt.Properties.VariableNames))

    % --- Model 1a: P300 predicted by surprise (all trials) ---
    gt_m1 = gt(~isnan(gt.P300_norm_z) & ~isnan(gt.surprise_z), :);

    mdl_p3_surprise = fitlme(gt_m1, ...
        'P300_norm_z ~ surprise_z * block_type + stage + (1 + surprise_z | subj_id)', ...
        'FitMethod','REML');
    fprintf('  Model 1a: P300 ~ surprise × block_type\n');
    disp(mdl_p3_surprise.Coefficients);

    % --- Model 1b: P300 predicted by omega (change-point probability) ---
    gt_m1b = gt(~isnan(gt.P300_norm_z) & ~isnan(gt.omega_z), :);

    mdl_p3_omega = fitlme(gt_m1b, ...
        'P300_norm_z ~ omega_z * block_type + stage + (1 | subj_id)', ...
        'FitMethod','REML');
    fprintf('  Model 1b: P300 ~ omega × block_type\n');
    disp(mdl_p3_omega.Coefficients);

    % --- Model 1c: Does P300 predict learning differently for true vs false FB? ---
    % Only P blocks have false feedback
    gt_p = gt(gt.block_type == 'P' & ~isnan(gt.P300_norm_z) & ~isnan(gt.alpha_nassar_z), :);

    if height(gt_p) > 50
        mdl_p3_learn = fitlme(gt_p, ...
            'alpha_nassar_z ~ P300_norm_z * false_fb_cat + stage + (1 | subj_id)', ...
            'FitMethod','REML');
        fprintf('  Model 1c: alpha_nassar ~ P300 × false_fb (P blocks)\n');
        disp(mdl_p3_learn.Coefficients);
    end

    % ── FIGURE MRQ1: P300 vs Surprise, split by block type ───────────────────
    fig1 = figure('Position', [50 50 1200 500], 'Color', 'w');
    sgtitle('MRQ1: P300 tracks surprise — modulated by feedback context', ...
        'FontSize', 14, 'FontWeight', 'bold');

    % Panel A: P300 ~ surprise, D vs P
    ax1 = subplot(1,3,1); hold(ax1, 'on');
    plot_binned_relationship(ax1, gt_m1, 'surprise_z', 'P300_norm_z', ...
        'block_type', {'D','P'}, {CLR_D, CLR_P}, {'Deterministic','Probabilistic'});
    xlabel(ax1, 'Surprise (\omega \times |\delta|) [z]');
    ylabel(ax1, 'P300 amplitude [norm z]');
    title(ax1, 'A. P300 ~ Surprise');
    legend(ax1, 'Box','off','Location','northwest','FontSize',9);

    % Panel B: P300 ~ omega, D vs P
    ax2 = subplot(1,3,2); hold(ax2, 'on');
    plot_binned_relationship(ax2, gt_m1b, 'omega_z', 'P300_norm_z', ...
        'block_type', {'D','P'}, {CLR_D, CLR_P}, {'Deterministic','Probabilistic'});
    xlabel(ax2, 'Change-point probability \omega [z]');
    ylabel(ax2, 'P300 amplitude [norm z]');
    title(ax2, 'B. P300 ~ \omega');

    % Panel C: P300 → learning rate, true vs false FB (P blocks)
    ax3 = subplot(1,3,3); hold(ax3, 'on');
    if height(gt_p) > 50
        plot_binned_relationship(ax3, gt_p, 'P300_norm_z', 'alpha_nassar_z', ...
            'false_fb_cat', {'True','False'}, {CLR_TRUE, CLR_FALSE}, ...
            {'True feedback','False feedback'});
    end
    xlabel(ax3, 'P300 amplitude [norm z]');
    ylabel(ax3, 'Learning rate \alpha [z]');
    title(ax3, 'C. P300 \rightarrow \alpha (P blocks)');
    legend(ax3, 'Box','off','Location','northwest','FontSize',9);

    apply_fig_style(fig1);
    save_fig(fig1, fullfile(figure_output_folder, 'MRQ1_P300_surprise'));
    fprintf('  Figure saved: MRQ1_P300_surprise\n');

else
    warning('MRQ1 skipped: missing required columns.');
    mdl_p3_surprise = []; mdl_p3_omega = []; mdl_p3_learn = [];
end

