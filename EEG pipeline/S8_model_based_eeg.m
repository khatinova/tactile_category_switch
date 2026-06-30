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


%% =========================================================================
%  MRQ2: FRN AND PREDICTION ERROR — FEEDBACK DISCOUNTING
%  =========================================================================
%  Hypothesis (Chase et al., 2011; Worthy et al., 2018):
%    - FRN (prefrontal negativity) scales with signed PE (more negative for
%      worse-than-expected outcomes).
%    - In P blocks, participants who learn to discount unreliable feedback
%      show attenuated FRN over time (Worthy 2018).
%    - False feedback should produce a PE signal but adaptive learners
%      should discount it (FRN attenuated for false vs true negative FB).
%
%  Key test: FRN × false_fb × trial_within_block (discounting dynamics)
% =========================================================================
fprintf('\n=== MRQ2: FRN ~ prediction error × feedback validity ===\n');

if all(ismember({'prefrontal_neg_peak_norm_z','PE_nassar_z','RW_PE_z','false_fb'}, gt.Properties.VariableNames))

    % --- Model 2a: FRN ~ Nassar PE (signed) × block type ---
    gt_m2 = gt(~isnan(gt.prefrontal_neg_peak_norm_z) & ~isnan(gt.PE_nassar_z), :);

    mdl_frn_pe = fitlme(gt_m2, ...
        'prefrontal_neg_peak_norm_z ~ PE_nassar_z * block_type + correct + stage + (1 + PE_nassar_z | subj_id)', ...
        'FitMethod','REML');
    fprintf('  Model 2a: FRN ~ Nassar PE × block_type\n');
    disp(mdl_frn_pe.Coefficients);

    % --- Model 2b: FRN ~ unsigned PE (salience) × block type ---
    gt_m2b = gt(~isnan(gt.prefrontal_neg_peak_norm_z) & ~isnan(gt.PE_unsigned_z), :);

    mdl_frn_upe = fitlme(gt_m2b, ...
        'prefrontal_neg_peak_norm_z ~ PE_unsigned_z * block_type + correct + stage + (1 | subj_id)', ...
        'FitMethod','REML');
    fprintf('  Model 2b: FRN ~ |PE| × block_type\n');
    disp(mdl_frn_upe.Coefficients);

    % --- Model 2c: False feedback dissociation (P blocks only) ---
    % Does FRN differentiate true from false negative feedback?
    gt_p_inc = gt(gt.block_type == 'P' & gt.correct == 0 & ...
                  ~isnan(gt.prefrontal_neg_peak_norm_z), :);

    if height(gt_p_inc) > 30
        mdl_frn_false = fitlme(gt_p_inc, ...
            'prefrontal_neg_peak_norm_z ~ false_fb_cat * stage + (1 | subj_id)', ...
            'FitMethod','REML');
        fprintf('  Model 2c: FRN ~ false_fb × stage (P-block incorrect trials)\n');
        disp(mdl_frn_false.Coefficients);
    else
        mdl_frn_false = [];
    end

    % --- Model 2d: Feedback discounting over time in P blocks ---
    % Does FRN attenuate within P blocks as participants learn to discount?
    gt_p_all = gt(gt.block_type == 'P' & ~isnan(gt.prefrontal_neg_peak_norm_z), :);
    if ismember('trial', gt_p_all.Properties.VariableNames)
        gt_p_all.trial_z = nan(height(gt_p_all), 1);
        for si = 1:numel(subj_list)
            mask = gt_p_all.subj_id == subj_list(si);
            tv = double(gt_p_all.trial(mask));
            mn = mean(tv,'omitnan'); sd = std(tv,'omitnan');
            if sd > 0, gt_p_all.trial_z(mask) = (tv - mn) / sd; end
        end

        mdl_frn_discount = fitlme(gt_p_all, ...
            'prefrontal_neg_peak_norm_z ~ trial_z * correct + stage + (1 + trial_z | subj_id)', ...
            'FitMethod','REML');
        fprintf('  Model 2d: FRN ~ trial (within-block) × correct (P blocks) — discounting\n');
        disp(mdl_frn_discount.Coefficients);
    else
        mdl_frn_discount = [];
    end

    % ── FIGURE MRQ2: FRN and PE ──────────────────────────────────────────────
    fig2 = figure('Position', [50 50 1400 500], 'Color', 'w');
    sgtitle('MRQ2: Prefrontal negativity (FRN) tracks prediction error — feedback discounting', ...
        'FontSize', 14, 'FontWeight', 'bold');

    % Panel A: FRN ~ signed PE, D vs P
    ax1 = subplot(1,4,1); hold(ax1, 'on');
    plot_binned_relationship(ax1, gt_m2, 'PE_nassar_z', 'prefrontal_neg_peak_norm_z', ...
        'block_type', {'D','P'}, {CLR_D, CLR_P}, {'Deterministic','Probabilistic'});
    xlabel(ax1, 'Signed PE (\delta) [z]');
    ylabel(ax1, 'Prefrontal negativity [norm z]');
    title(ax1, 'A. FRN ~ signed PE');
    set(ax1, 'YDir', 'reverse');

    % Panel B: FRN ~ unsigned PE (salience), D vs P
    ax2 = subplot(1,4,2); hold(ax2, 'on');
    plot_binned_relationship(ax2, gt_m2b, 'PE_unsigned_z', 'prefrontal_neg_peak_norm_z', ...
        'block_type', {'D','P'}, {CLR_D, CLR_P}, {'Deterministic','Probabilistic'});
    xlabel(ax2, '|PE| (salience) [z]');
    ylabel(ax2, 'Prefrontal negativity [norm z]');
    title(ax2, 'B. FRN ~ |PE| (salience)');
    set(ax2, 'YDir', 'reverse');

    % Panel C: FRN for true vs false negative FB (P blocks, incorrect)
    ax3 = subplot(1,4,3); hold(ax3, 'on');
    if ~isempty(mdl_frn_false) && height(gt_p_inc) > 30
        stages = {'LN','LE','RN','RE'};
        for fi = 0:1
            fb_mask = gt_p_inc.false_fb == fi;
            means = nan(1,4); sems = nan(1,4);
            for si2 = 1:4
                sm = fb_mask & gt_p_inc.stage == stages{si2};
                vals = gt_p_inc.prefrontal_neg_peak_norm_z(sm);
                means(si2) = mean(vals,'omitnan');
                sems(si2) = std(vals,'omitnan') / sqrt(max(sum(~isnan(vals)),1));
            end
            clr = ternary_s8(fi==0, CLR_TRUE, CLR_FALSE);
            lbl = ternary_s8(fi==0, 'True negative FB', 'False negative FB');
            x_off = (fi-0.5)*0.1;
            errorbar(ax3, (1:4)+x_off, means, sems, 'o-', ...
                'Color', clr, 'LineWidth', 1.8, 'MarkerFaceColor', clr, ...
                'MarkerSize', 7, 'DisplayName', lbl);
        end
        set(ax3, 'XTick', 1:4, 'XTickLabel', stages);
        xlabel(ax3, 'Stage');
        ylabel(ax3, 'Prefrontal negativity [norm z]');
        title(ax3, 'C. True vs False negative FB');
        set(ax3, 'YDir', 'reverse');
        legend(ax3, 'Box','off','Location','best','FontSize',9);
    else
        text(ax3, 0.5, 0.5, 'Insufficient data', 'HorizontalAlignment','center');
    end

    % Panel D: Discounting — FRN across trial position (P blocks)
    ax4 = subplot(1,4,4); hold(ax4, 'on');
    if ~isempty(mdl_frn_discount) && ismember('trial', gt_p_all.Properties.VariableNames)
        % Bin trials into thirds: early, middle, late
        for ci = 0:1
            c_mask = gt_p_all.correct == ci;
            trial_vals = double(gt_p_all.trial(c_mask));
            frn_vals   = gt_p_all.prefrontal_neg_peak_norm_z(c_mask);
            edges = quantile(trial_vals, [0 1/3 2/3 1]);
            bin_means = nan(1,3); bin_sems = nan(1,3);
            for bi = 1:3
                bm = trial_vals >= edges(bi) & trial_vals < edges(bi+1) + (bi==3);
                v = frn_vals(bm);
                bin_means(bi) = mean(v,'omitnan');
                bin_sems(bi) = std(v,'omitnan')/sqrt(max(sum(~isnan(v)),1));
            end
            clr = ternary_s8(ci==1, CLR_TRUE, CLR_FALSE);
            lbl = ternary_s8(ci==1, 'Correct', 'Incorrect');
            x_off = (ci-0.5)*0.12;
            errorbar(ax4, (1:3)+x_off, bin_means, bin_sems, 'o-', ...
                'Color', clr, 'LineWidth', 1.8, 'MarkerFaceColor', clr, ...
                'MarkerSize', 7, 'DisplayName', lbl);
        end
        set(ax4, 'XTick', 1:3, 'XTickLabel', {'Early','Middle','Late'});
        xlabel(ax4, 'Trial position (within block)');
        ylabel(ax4, 'Prefrontal negativity [norm z]');
        title(ax4, 'D. Feedback discounting (P blocks)');
        set(ax4, 'YDir', 'reverse');
        legend(ax4, 'Box','off','Location','best','FontSize',9);
    end

    apply_fig_style(fig2);
    save_fig(fig2, fullfile(figure_output_folder, 'MRQ2_FRN_prediction_error'));
    fprintf('  Figure saved: MRQ2_FRN_prediction_error\n');

else
    warning('MRQ2 skipped: missing required columns.');
    mdl_frn_pe = []; mdl_frn_upe = []; mdl_frn_false = []; mdl_frn_discount = [];
end


%% =========================================================================
%  MRQ3: FRONTAL THETA ~ UNSIGNED PE → BEHAVIOURAL ADJUSTMENT
%  =========================================================================
%  Hypothesis (Cavanagh et al., 2010; Reteig et al., 2020):
%    - Mid-frontal theta power scales parametrically with unsigned PE
%      (the "need for control" signal).
%    - Theta predicts trial-to-trial behavioural adjustment (stay/switch).
%    - This relationship should be stronger in D blocks (where feedback
%      is always valid) than in P blocks (where learned discounting should
%      attenuate the theta → behaviour link).
%
%  Key test: Theta ~ |PE| and Theta → next_correct
% =========================================================================
fprintf('\n=== MRQ3: Frontal theta ~ |PE| → behavioural adjustment ===\n');

if all(ismember({'Theta_amp_z','PE_unsigned_z','surprise_z'}, gt.Properties.VariableNames))

    % --- Model 3a: Theta ~ unsigned PE × block type ---
    gt_m3 = gt(~isnan(gt.Theta_amp_z) & ~isnan(gt.PE_unsigned_z), :);

    mdl_theta_pe = fitlme(gt_m3, ...
        'Theta_amp_z ~ PE_unsigned_z * block_type + correct + stage + (1 + PE_unsigned_z | subj_id)', ...
        'FitMethod','REML');
    fprintf('  Model 3a: Theta ~ |PE| × block_type\n');
    disp(mdl_theta_pe.Coefficients);

    % --- Model 3b: Theta ~ surprise (Nassar) × block type ---
    gt_m3b = gt(~isnan(gt.Theta_amp_z) & ~isnan(gt.surprise_z), :);

    mdl_theta_surp = fitlme(gt_m3b, ...
        'Theta_amp_z ~ surprise_z * block_type + stage + (1 | subj_id)', ...
        'FitMethod','REML');
    fprintf('  Model 3b: Theta ~ surprise × block_type\n');
    disp(mdl_theta_surp.Coefficients);

    % --- Model 3c: Does theta predict next-trial accuracy? ---
    gt_m3c = gt(~isnan(gt.Theta_amp_z) & ~isnan(gt.next_correct), :);

    if height(gt_m3c) > 50
        mdl_theta_next = fitglme(gt_m3c, ...
            'next_correct ~ Theta_amp_z * block_type + correct + stage + (1 | subj_id)', ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        fprintf('  Model 3c: next_correct ~ Theta × block_type (logistic)\n');
        disp(mdl_theta_next.Coefficients);
    else
        mdl_theta_next = [];
    end

    % ── FIGURE MRQ3: Theta and PE ────────────────────────────────────────────
    fig3 = figure('Position', [50 50 1200 500], 'Color', 'w');
    sgtitle('MRQ3: Frontal theta tracks |PE| and predicts behavioural adjustment', ...
        'FontSize', 14, 'FontWeight', 'bold');

    % Panel A: Theta ~ |PE|, D vs P
    ax1 = subplot(1,3,1); hold(ax1, 'on');
    plot_binned_relationship(ax1, gt_m3, 'PE_unsigned_z', 'Theta_amp_z', ...
        'block_type', {'D','P'}, {CLR_D, CLR_P}, {'Deterministic','Probabilistic'});
    xlabel(ax1, '|PE| (unsigned) [z]');
    ylabel(ax1, 'Frontal theta power [z]');
    title(ax1, 'A. Theta ~ |PE|');
    legend(ax1, 'Box','off','Location','northwest','FontSize',9);

    % Panel B: Theta ~ surprise (Nassar), D vs P
    ax2 = subplot(1,3,2); hold(ax2, 'on');
    plot_binned_relationship(ax2, gt_m3b, 'surprise_z', 'Theta_amp_z', ...
        'block_type', {'D','P'}, {CLR_D, CLR_P}, {'Deterministic','Probabilistic'});
    xlabel(ax2, 'Surprise (\omega \times |\delta|) [z]');
    ylabel(ax2, 'Frontal theta power [z]');
    title(ax2, 'B. Theta ~ Surprise');

    % Panel C: Theta predicts next-trial accuracy (binned)
    ax3 = subplot(1,3,3); hold(ax3, 'on');
    if ~isempty(mdl_theta_next)
        for bt_i = 1:2
            bt_lbl = ternary_s8(bt_i==1, 'D', 'P');
            bt_mask = gt_m3c.block_type == bt_lbl;
            theta_vals = gt_m3c.Theta_amp_z(bt_mask);
            next_vals  = gt_m3c.next_correct(bt_mask);
            edges_t = quantile(theta_vals(~isnan(theta_vals)), linspace(0,1,6));
            bin_x = nan(1,5); bin_y = nan(1,5); bin_se = nan(1,5);
            for bi = 1:5
                bm = theta_vals >= edges_t(bi) & theta_vals < edges_t(bi+1) + (bi==5)*0.01;
                nv = next_vals(bm);
                bin_x(bi) = mean(theta_vals(bm),'omitnan');
                bin_y(bi) = mean(nv,'omitnan');
                bin_se(bi) = std(nv,'omitnan')/sqrt(max(sum(~isnan(nv)),1));
            end
            clr = ternary_s8(bt_i==1, CLR_D, CLR_P);
            errorbar(ax3, bin_x, bin_y, bin_se, 'o-', ...
                'Color', clr, 'LineWidth', 1.8, 'MarkerFaceColor', clr, ...
                'MarkerSize', 7, 'DisplayName', [bt_lbl ' blocks']);
        end
        xlabel(ax3, 'Frontal theta [z]');
        ylabel(ax3, 'P(correct next trial)');
        title(ax3, 'C. Theta \rightarrow next-trial accuracy');
        legend(ax3, 'Box','off','Location','best','FontSize',9);
    end

    apply_fig_style(fig3);
    save_fig(fig3, fullfile(figure_output_folder, 'MRQ3_theta_PE_adjustment'));
    fprintf('  Figure saved: MRQ3_theta_PE_adjustment\n');

else
    warning('MRQ3 skipped: missing required columns.');
    mdl_theta_pe = []; mdl_theta_surp = []; mdl_theta_next = [];
end


%% =========================================================================
%  MRQ4: FRONTO-PARIETAL PLV — THETA COUPLING AS MECHANISM
%  =========================================================================
%  Hypothesis: frontal theta is the mechanism through which mPFC
%  communicates "need for cognitive control" to parietal areas. If so:
%    - Trial-level theta should predict trial-level FP PLV
%    - FP PLV should mediate the relationship between PE and next-trial
%      behavioural adjustment
%    - This coupling should be upregulated after surprising events
%      (high omega / high surprise) especially during genuine change-points
%
%  Key test: PLV ~ Theta and PLV mediates PE → behaviour
% =========================================================================
fprintf('\n=== MRQ4: Fronto-parietal PLV ~ theta coupling ===\n');

if all(ismember({'PLV_fp_z','Theta_amp_z','surprise_z','omega_z'}, gt.Properties.VariableNames))

    % --- Model 4a: PLV ~ theta × block type ---
    gt_m4 = gt(~isnan(gt.PLV_fp_z) & ~isnan(gt.Theta_amp_z), :);

    mdl_plv_theta = fitlme(gt_m4, ...
        'PLV_fp_z ~ Theta_amp_z * block_type + stage + (1 + Theta_amp_z | subj_id)', ...
        'FitMethod','REML');
    fprintf('  Model 4a: PLV ~ Theta × block_type\n');
    disp(mdl_plv_theta.Coefficients);

    % --- Model 4b: PLV ~ surprise × block type ---
    gt_m4b = gt(~isnan(gt.PLV_fp_z) & ~isnan(gt.surprise_z), :);

    mdl_plv_surp = fitlme(gt_m4b, ...
        'PLV_fp_z ~ surprise_z * block_type + stage + (1 | subj_id)', ...
        'FitMethod','REML');
    fprintf('  Model 4b: PLV ~ surprise × block_type\n');
    disp(mdl_plv_surp.Coefficients);

    % --- Model 4c: Mediation test (simplified) ---
    % Step 1: PE → PLV (path a)
    % Step 2: PLV → next_correct controlling for PE (path b)
    gt_m4c = gt(~isnan(gt.PLV_fp_z) & ~isnan(gt.PE_unsigned_z) & ~isnan(gt.next_correct), :);

    if height(gt_m4c) > 50
        % Path a: PE → PLV
        mdl_path_a = fitlme(gt_m4c, ...
            'PLV_fp_z ~ PE_unsigned_z + block_type + stage + (1 | subj_id)', ...
            'FitMethod','REML');
        fprintf('  Mediation path a (PE → PLV):\n');
        disp(mdl_path_a.Coefficients(2,:));  % PE_unsigned_z coefficient

        % Path b + c': PLV → next_correct controlling for PE
        mdl_path_bc = fitglme(gt_m4c, ...
            'next_correct ~ PLV_fp_z + PE_unsigned_z + block_type + stage + (1 | subj_id)', ...
            'Distribution','Binomial','Link','logit','FitMethod','Laplace');
        fprintf('  Mediation path b+c'' (PLV + PE → next_correct):\n');
        disp(mdl_path_bc.Coefficients(2:3,:));  % PLV and PE coefficients
    else
        mdl_path_a = []; mdl_path_bc = [];
    end

    % ── FIGURE MRQ4: PLV and theta coupling ──────────────────────────────────
    fig4 = figure('Position', [50 50 1200 500], 'Color', 'w');
    sgtitle('MRQ4: Fronto-parietal PLV tracks theta and mediates PE \rightarrow behaviour', ...
        'FontSize', 14, 'FontWeight', 'bold');

    % Panel A: PLV ~ Theta
    ax1 = subplot(1,3,1); hold(ax1, 'on');
    plot_binned_relationship(ax1, gt_m4, 'Theta_amp_z', 'PLV_fp_z', ...
        'block_type', {'D','P'}, {CLR_D, CLR_P}, {'Deterministic','Probabilistic'});
    xlabel(ax1, 'Frontal theta power [z]');
    ylabel(ax1, 'Fronto-parietal PLV [z]');
    title(ax1, 'A. PLV ~ Theta');
    legend(ax1, 'Box','off','Location','northwest','FontSize',9);

    % Panel B: PLV ~ surprise
    ax2 = subplot(1,3,2); hold(ax2, 'on');
    plot_binned_relationship(ax2, gt_m4b, 'surprise_z', 'PLV_fp_z', ...
        'block_type', {'D','P'}, {CLR_D, CLR_P}, {'Deterministic','Probabilistic'});
    xlabel(ax2, 'Surprise (\omega \times |\delta|) [z]');
    ylabel(ax2, 'Fronto-parietal PLV [z]');
    title(ax2, 'B. PLV ~ Surprise');

    % Panel C: Mediation path diagram (schematic + stats)
    ax3 = subplot(1,3,3); hold(ax3, 'on'); axis(ax3, 'off');
    title(ax3, 'C. Mediation summary', 'FontSize', 12);
    if ~isempty(mdl_path_a) && ~isempty(mdl_path_bc)
        a_coef = mdl_path_a.Coefficients.Estimate(2);
        a_p    = mdl_path_a.Coefficients.pValue(2);
        b_coef = mdl_path_bc.Coefficients.Estimate(2);
        b_p    = mdl_path_bc.Coefficients.pValue(2);
        cp_coef = mdl_path_bc.Coefficients.Estimate(3);
        cp_p   = mdl_path_bc.Coefficients.pValue(3);

        med_txt = { ...
            'PE → PLV (path a):', ...
            sprintf('  b = %.3f, p = %s', a_coef, format_p_s8(a_p)), ...
            '', ...
            'PLV → next correct (path b):', ...
            sprintf('  OR = %.2f, p = %s', exp(b_coef), format_p_s8(b_p)), ...
            '', ...
            'PE → next correct (path c''):', ...
            sprintf('  OR = %.2f, p = %s', exp(cp_coef), format_p_s8(cp_p)), ...
            '', ...
            sprintf('Indirect effect (a×b) = %.4f', a_coef * b_coef)};
        text(ax3, 0.1, 0.9, med_txt, 'VerticalAlignment', 'top', ...
            'FontSize', 10, 'FontName', 'Arial');
    else
        text(ax3, 0.5, 0.5, 'Insufficient data for mediation', ...
            'HorizontalAlignment','center','FontSize',11);
    end

    apply_fig_style(fig4);
    save_fig(fig4, fullfile(figure_output_folder, 'MRQ4_PLV_theta_mediation'));
    fprintf('  Figure saved: MRQ4_PLV_theta_mediation\n');

else
    warning('MRQ4 skipped: missing required columns.');
    mdl_plv_theta = []; mdl_plv_surp = []; mdl_path_a = []; mdl_path_bc = [];
end

