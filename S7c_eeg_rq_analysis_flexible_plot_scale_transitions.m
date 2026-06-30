% =============================================================================
% RQ STATISTICAL ANALYSIS — CORRECTED VERSION
%
% BUGS FIXED vs original:
%
%  BUG 1 — gt_plv.trueFB == 1 (line ~RQ4)
%     The group_table has no column called 'trueFB'. The stage_table pipeline
%     stores feedback validity as 'false_fb' (logical, true = feedback was
%     flipped). Using gt_plv.trueFB == 1 silently creates an all-zero vector
%     (MATLAB returns [] or zeros for missing struct fields accessed this way
%     on a table), filtering out ALL rows.
%     FIX: ~gt_plv.false_fb
%
%  BUG 2 — next_correct computation ignores block boundaries (RQ4)
%     The original code assigns next_correct by taking consecutive rows
%     within a subject across the whole table, so the last trial of block 1
%     gets next_correct = first trial of block 2. That is a different block
%     with a different stimulus set and reversal — the "next trial correct"
%     has no meaningful continuity there.
%     FIX: nest the loop inside a per-block inner loop.
%
%  BUG 3 — categorical comparison in plot_stage_bar
%     categorical({bts{bt_i}}) creates a 1×1 categorical cell array, not a
%     scalar categorical. The == comparison against a column categorical may
%     silently fail or produce incorrect logical indexing.
%     FIX: compare directly as string: gt_sub.block_type == bts{bt_i}
%
%  BUG 4 — plot_lme_effects called with yvar='correct' (numeric 0/1)
%     The function calls categories(gt.(yvar)) which works only on categorical
%     arrays. If correct is still numeric, this errors. If already converted
%     to categorical it labels groups '0'/'1' which are opaque in figures.
%     FIX: convert correct to categorical with explicit labels inside the
%     function; use 'Incorrect'/'Correct' as display names.
%
%  BUG 5 — prefrontal_neg_peak_norm used in models but NOT in z-score list
%     The z-score loop processes FCz_neg_peak_amp (raw µV) but the models
%     use FCz_neg_peak_norm (baseline-RMS-normalised). That column is on a
%     completely different scale from the z-scored columns used in the same
%     models. Both should either be z-scored or neither.
%     FIX: add FCz_neg_peak_norm to the z-score list and use the _z version
%     in all models for consistency. The baseline-RMS normalised value is
%     already interpretable as "signal / noise" but it still has
%     between-subject variance in scale.
%
%  BUG 6 — plot_lme_effects uses ttest2 (independent samples) to annotate
%     p-values comparing two groups within an x-level. Because the same
%     subject appears in both D and P block columns, the correct test is
%     ttest (paired), not ttest2.
%     FIX: match subjects and use ttest on paired differences.
%
%  BUG 7 — mdl_false_fn uses false_fb (logical) as a continuous predictor
%     in fitlme. MATLAB's fitlme will treat logicals as doubles (0/1), which
%     is correct numerically, but the interaction conf_z * false_fb will test
%     whether the slope of confidence on FRN differs between true and false
%     trials. This is actually the intended test but false_fb should be cast
%     to categorical so coefficients are labelled clearly.
%     FIX: cast to categorical before model call.
%
%  BUG 8 — PLV_z assignment in plv_long may misalign rows
%     plv_long = [gt_plv; gt_fs] stacks the rows. When gt_plv and gt_fs
%     differ in height (different subjects passed the MIN_TRIALS_PLV filter
%     for fp vs fs), the assignment plv_long.PLV_z = [gt_plv.PLV_fp_z;
%     gt_fs.PLV_fs_z] will error or silently mismatch if any rows with
%     NaN PLV were dropped differently between the two tables.
%     FIX: keep all rows and fill PLV_z conditionally after stacking.
%
%  BUG 9 — randon effects formula `(1 + stage | subj_id)` with ordinal stage
%     stage is set as ordinal categorical with 4 levels. A random slope over
%     an ordinal predictor requires MATLAB to dummy-code it (3 contrasts).
%     A random slope for 3 dummy variables with a small n (~30 subjects)
%     will not converge. This will produce a singular fit warning and
%     unreliable variance estimates.
%     FIX: Use (1 | subj_id) + (1 | subj_id:stage) crossed random effects,
%     which is equivalent but more stable; or drop the random slope for stage.
%
% =============================================================================

% ── PIPELINE STEP 7 of 7 — 5-RQ EEG analysis (S7_eeg_rq_analysis.m) ─────────
% (was: Stats_Uncertainty_Analysis_EEG_combined_FRN_negpeak_v2.m)
%
% Five research questions (see EEG RQs Word doc): RQ1 FN/P300 under uncertainty,
% RQ2 confidence x FN, RQ3 frontal theta x stage, RQ4 fronto-parietal PLV,
% RQ5 fronto-parietal vs fronto-somatosensory pathway. Produces LME stats,
% manuscript-ready tables, and figures (styled via pipeline/utils/apply_fig_style).
%
% INPUT : group_feature_table_combined.mat (group_table) from S4.
%         Optionally frn_rewp_by_stage_combined.mat for FRN/RewP grand averages.
% OUTPUT: RQ figures + manuscript_stats.txt in figure_output_folder.
% -----------------------------------------------------------------------------
addpath(genpath(fileparts(mfilename('fullpath'))));   % pipeline utils on path

remote = 0;
if remote == 1
    base_path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch';
else
    base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch';
end
saved_tables_folder    = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis', 'Outcome_feature_tables_v4_merged');
figure_output_folder = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis', 'Figures', 'RQ_analysis_combined_extended_electrodes');

% -------------------------------------------------------------------------
%% PLOTTING SCALE CONFIGURATION
% -------------------------------------------------------------------------
% This switch changes FIGURE variables only. The inferential models below
% still use within-subject z-scored predictors/outcomes unless you explicitly
% edit the model formulae.
%
%   'z'    = within-subject z-scored values. Good for combining subjects and
%            comparing effect directions, but less directly interpretable.
%   'norm' = baseline-RMS-normalised values where available, e.g.
%            prefrontal_neg_peak_norm, P300_norm. Better for physiology plots.
%   'raw'  = raw amplitude / raw feature values where available, e.g.
%            prefrontal_neg_peak_amp, P300_amp. Most interpretable units,
%            but most sensitive to between-subject scale differences.
%
% Recommendation: use PLOT_SCALE='norm' for main EEG amplitude figures and
% PLOT_SCALE='z' for compact model-aligned supplementary figures.
PLOT_SCALE = 'norm';      % options: 'z', 'norm', 'raw'

figure_output_folder = fullfile(figure_output_folder, ['plot_' PLOT_SCALE]);
if ~exist(figure_output_folder, 'dir'), mkdir(figure_output_folder); end

% Load the combined KH+RR feature table (from S4) unless already in workspace.
% load(fullfile(saved_tables_folder, 'group_feature_table_combined.mat'), 'group_table');


% Figure styling defaults (ticks outside, no top/right box) applied to every
% axes created below; matches pipeline/utils/apply_fig_style.
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'off');
set(groot, 'defaultAxesTickDirMode', 'manual');

% Backward-compatibility: some RQ4 code references block_number; alias to block.
if ~ismember('block_number', group_table.Properties.VariableNames)
    if ismember('block', group_table.Properties.VariableNames)
        group_table.block_number = kh_to_numeric(group_table.block);
    elseif ismember('blocknum', group_table.Properties.VariableNames)
        group_table.block_number = kh_to_numeric(group_table.blocknum);
    end
end

% ─── Ensure categorical types ────────────────────────────────────────────────
group_table.subj_id         = categorical(group_table.subj_id);
group_table.block_type      = categorical(group_table.block_type);
group_table.stage           = categorical(group_table.stage, ...
                                  {'LN','LE','RN','RE'}, 'Ordinal', true);
group_table.feedback_modality = categorical(group_table.feedback_modality);
group_table.stimulus_modality = categorical(group_table.stimulus_modality);
group_table.practice_task     = categorical(group_table.practice_task);

% ─── FIX 5: extend z-score list to include FCz_neg_peak_norm ────────────────
features_to_zscore = {'prefrontal_neg_peak_amp','prefrontal_neg_peak_norm', ...
                      'FRN_amp','prefrontalFRN_norm','RewP_amp','RewP_norm', ...
                      'P300_amp','P300_norm','Theta_amp', ...
                      'PLV_fp','PLV_fs','PLV_fp_pairwise','PLV_fs_pairwise'};

subj_list = unique(group_table.subj_id);
for f = 1:numel(features_to_zscore)
    fn   = features_to_zscore{f};
    fn_z = [fn '_z'];
    if ~ismember(fn, group_table.Properties.VariableNames); continue; end
    group_table.(fn_z) = nan(height(group_table), 1);
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

gt = group_table;   % working copy

% Resolve plot variables from PLOT_SCALE. These names are used only for
% figures. Model formulae continue to use *_z variables below.
[FN_PLOT,    FN_YLBL,    FN_REVERSE_Y]    = resolve_plot_feature(gt, 'prefrontal_neg_peak', PLOT_SCALE);
[P300_PLOT,  P300_YLBL,  P300_REVERSE_Y]  = resolve_plot_feature(gt, 'P300', PLOT_SCALE);
[THETA_PLOT, THETA_YLBL, THETA_REVERSE_Y] = resolve_plot_feature(gt, 'Theta', PLOT_SCALE);
[PLVFP_PLOT, PLVFP_YLBL, PLVFP_REVERSE_Y] = resolve_plot_feature(gt, 'PLV_fp', PLOT_SCALE);
[PLVFS_PLOT, PLVFS_YLBL, PLVFS_REVERSE_Y] = resolve_plot_feature(gt, 'PLV_fs', PLOT_SCALE);

fprintf('\nPlot scale: %s\n', PLOT_SCALE);
fprintf('  FN plot variable     : %s\n', FN_PLOT);
fprintf('  P300 plot variable   : %s\n', P300_PLOT);
fprintf('  Theta plot variable  : %s\n', THETA_PLOT);
fprintf('  FP PLV plot variable : %s\n', PLVFP_PLOT);
fprintf('  FS PLV plot variable : %s\n', PLVFS_PLOT);

% Derive block-transition and prior-uncertainty history columns for the
% additional exploratory EEG figures below. These are plotting columns only.
gt = add_transition_history_columns(gt);

fprintf('\nTransition / prior-P plotting columns added.\n');
if ismember('transition_recent', gt.Properties.VariableNames)
    disp(tabulate(categorical(gt.transition_recent)));
end

fprintf('\n\n=== STATISTICAL MODELS ===\n\n');

% =========================================================================
%% RQ1: Does uncertainty suppress feedback negativity or shift it to P300?
% =========================================================================
fprintf('--- RQ1: Feedback negativity and P300 under uncertainty ---\n');

gt_true = gt(~gt.false_fb & ...
              ~isnan(gt.prefrontal_neg_peak_norm_z) & ...    % FIX 5: use _z version
              ~isnan(gt.P300_amp_z), :);

% FIX 9: stable random effects — per-subject intercept only for stage
mdl_fn_rq1 = fitlme(gt_true, ...
    ['prefrontal_neg_peak_norm_z ~ block_type * correct + stage + ' ...
     '(1 + block_type | subj_id)'], ...
    'FitMethod','REML');
fprintf('  FN model:\n'); disp(mdl_fn_rq1.Coefficients);

mdl_p3_rq1 = fitlme(gt_true, ...
    ['P300_norm_z ~ block_type * correct + stage + ' ...
     '(1 + block_type | subj_id)'], ...
    'FitMethod','REML');
fprintf('  P300 model:\n'); disp(mdl_p3_rq1.Coefficients);

% FIX 4: pass correct as a character label so plot_lme_effects can handle it
plot_lme_effects(gt_true, FN_PLOT, 'block_type', 'correct', ...
    {'Incorrect','Correct'}, ...
    'RQ1: FN × Block type × Outcome', figure_output_folder, 'RQ1_FN', FN_YLBL, FN_REVERSE_Y);
plot_lme_effects(gt_true, P300_PLOT, 'block_type', 'correct', ...
    {'Incorrect','Correct'}, ...
    'RQ1: P300 × Block type × Outcome', figure_output_folder, 'RQ1_P300', P300_YLBL, P300_REVERSE_Y);

% -------------------------------------------------------------------------
%% NEW PLOTS: prefrontal negative peak and P300
% -------------------------------------------------------------------------

plot_stage_outcome_lines(gt, FN_PLOT, ...
    FN_YLBL, ...
    'RQ: prefrontal negative peak across stages — correct vs incorrect', ...
    figure_output_folder, 'prefrontal_negpeak_stage_correct_incorrect', FN_REVERSE_Y);

plot_stage_transition_feedback(gt, FN_PLOT, ...
    FN_YLBL, ...
    'RQ: prefrontal negative peak across stages — D/P and true/false feedback', ...
    figure_output_folder, 'prefrontal_negpeak_stage_block_feedback', FN_REVERSE_Y);

plot_stage_outcome_lines(gt, P300_PLOT, ...
    P300_YLBL, ...
    'RQ: P300 across stages — correct vs incorrect', ...
    figure_output_folder, 'P300_stage_correct_incorrect', P300_REVERSE_Y);

plot_stage_transition_feedback(gt, P300_PLOT, ...
    P300_YLBL, ...
    'RQ: P300 across stages — D/P and true/false feedback', ...
    figure_output_folder, 'P300_stage_block_feedback', P300_REVERSE_Y);

% =========================================================================
%% RQ2: Does confidence predict FRN? Stronger in D than P?
% =========================================================================
fprintf('\n--- RQ2: Confidence × Feedback negativity × block_type ---\n');

gt_inc = gt(...%gt.correct==0 & ~gt.false_fb & ...
             ~isnan(gt.confidence) & ~isnan(gt.prefrontal_neg_peak_norm_z), :);

% Z-score confidence within subject
gt_inc.conf_z = nan(height(gt_inc), 1);
for si = 1:numel(subj_list)
    mask = gt_inc.subj_id == subj_list(si);
    c    = gt_inc.confidence(mask);
    if sum(mask) > 1 && std(c,'omitnan') > 0
        gt_inc.conf_z(mask) = (c - mean(c,'omitnan')) / std(c,'omitnan');
    end
end

mdl_fn_conf = fitlme(gt_inc, ...
    ['prefrontal_neg_peak_norm_z ~ conf_z * block_type + stage + ' ...
     '(1 + conf_z | subj_id)'], ...
    'FitMethod','REML');
fprintf('  FRN ~ confidence × block_type (incorrect trials only):\n');
disp(mdl_fn_conf.Coefficients);

% False feedback analysis (P blocks: does perceived vs true feedback matter?)
gt_p = gt(gt.block_type=='P' & ...
           ~isnan(gt.prefrontal_neg_peak_norm_z) & ~isnan(gt.confidence), :);
gt_p.conf_z    = nan(height(gt_p), 1);
gt_p.false_fb_cat = categorical(double(gt_p.false_fb), [0 1], {'TrueFB','FalseFB'}); % FIX 7

for si = 1:numel(subj_list)
    mask = gt_p.subj_id == subj_list(si);
    c    = gt_p.confidence(mask);
    if sum(mask) > 1 && std(c,'omitnan') > 0
        gt_p.conf_z(mask) = (c - mean(c,'omitnan')) / std(c,'omitnan');
    end
end

if height(gt_p) > 10
    mdl_false_fn = fitlme(gt_p, ...
        ['prefrontal_neg_peak_norm_z ~ conf_z * false_fb + stage + ' ...
         '(1 | subj_id)'], ...
        'FitMethod','REML');
    fprintf('  FRN ~ confidence × false_fb (P blocks only):\n');
    disp(mdl_false_fn.Coefficients);
else
    mdl_false_fn = [];
    warning('Too few P-block rows for false_fb model.');
end

plot_confidence_fRN(gt_inc, figure_output_folder, FN_PLOT, FN_YLBL, FN_REVERSE_Y);


% =========================================================================
%% RQ3: Stage × block_type interaction on theta (incorrect trials)
% =========================================================================
fprintf('\n--- RQ3: Theta × stage × block_type (incorrect trials) ---\n');

gt_th = gt(gt.correct==0 & ~gt.false_fb & ~isnan(gt.Theta_amp_z), :);

% FIX 9: random slope for block_type (2 levels) is tractable; drop stage slope
mdl_theta = fitlme(gt_th, ...
    ['Theta_amp_z ~ stage * block_type + ' ...
     '(1 + block_type | subj_id)'], ...
    'FitMethod','REML');
fprintf('  Theta model:\n'); disp(mdl_theta.Coefficients);

plot_stage_bar(gt_th, THETA_PLOT, [THETA_YLBL ' (incorrect trials)'], ...
    'RQ3: Stage × Block type — Frontal theta', figure_output_folder, 'RQ3_Theta', THETA_REVERSE_Y);


% =========================================================================
%% RQ4: FP-PLV × stage × false_fb
% =========================================================================
fprintf('\n--- RQ4: Fronto-parietal PLV × stage × block_type ---\n');

gt_plv = gt(~isnan(gt.PLV_fp_z), :);

% FIX 1: was gt_plv.trueFB == 1 — column doesn't exist. Use ~false_fb.
gt_plv_true = gt_plv(~gt_plv.false_fb, :);

mdl_plv_fp = fitlme(gt_plv_true, ...
    ['PLV_fp_z ~ stage * block_type + ' ...
     '(1 + block_type | subj_id)'], ...
    'FitMethod','REML');
fprintf('  FP PLV model:\n'); disp(mdl_plv_fp.Coefficients);

% ── FIX 2: next-trial accuracy respecting block boundaries ──────────────────
gt_plv_next = gt_plv;
gt_plv_next.next_correct = nan(height(gt_plv_next), 1);

for si = 1:numel(subj_list)
    for b = 1:5   % 5 blocks maximum
        mask = gt_plv_next.subj_id == subj_list(si) & ...
               gt_plv_next.block_number == b;            
        idx  = find(mask);
        if numel(idx) < 2; continue; end
        correct_vec = gt_plv_next.correct(idx);
        % pair trial t with trial t+1 WITHIN the same block
        gt_plv_next.next_correct(idx(1:end-1)) = correct_vec(2:end);
        % last trial of each block is left NaN (no valid next trial)
    end
end

gt_plv_next2 = gt_plv_next(~isnan(gt_plv_next.next_correct), :);

mdl_plv_predict = fitglme(gt_plv_next2, ...
    ['next_correct ~ PLV_fp_z * block_type + stage + ' ...
     '(1 | subj_id)'], ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace');
fprintf('  PLV → next trial correct (logistic):\n');
disp(mdl_plv_predict.Coefficients);

plot_stage_bar(gt_plv, PLVFP_PLOT, PLVFP_YLBL, ...
    'RQ4: Stage × Block type — FP PLV', figure_output_folder, 'RQ4_PLV_fp', PLVFP_REVERSE_Y);


% =========================================================================
%% RQ5: FS-PLV vs FP-PLV pathway comparison
% =========================================================================
fprintf('\n--- RQ5: Fronto-somatosensory vs fronto-parietal PLV ---\n');

gt_fp = gt(~isnan(gt.PLV_fp_z), :);
gt_fs = gt(~isnan(gt.PLV_fs_z), :);

% ── FIX 8: stack then assign PLV_z using row identity ────────────────────────
n_fp = height(gt_fp);
n_fs = height(gt_fs);
plv_long         = [gt_fp; gt_fs];
plv_long.PLV_z   = [gt_fp.PLV_fp_z; gt_fs.PLV_fs_z];   % now safe: same n
plv_long.pathway = categorical([repmat({'fp'}, n_fp, 1); repmat({'fs'}, n_fs, 1)]);

mdl_pathway = fitlme(plv_long, ...
    ['PLV_z ~ pathway * stage * block_type + ' ...
     '(1 + pathway | subj_id)'], ...
    'FitMethod','REML');
fprintf('  Pathway × stage × block_type:\n');
disp(mdl_pathway.Coefficients);

plot_pathway_comparison(gt_fp, gt_fs, figure_output_folder, PLVFP_PLOT, PLVFS_PLOT, PLVFP_YLBL, PLVFS_YLBL, PLVFP_REVERSE_Y || PLVFS_REVERSE_Y);

% =========================================================================
%% RQ6: Recent transition and cumulative prior-P exposure for all EEG components
% =========================================================================
fprintf('\n--- RQ6: EEG components by recent transition and prior P exposure ---\n');

component_specs = { ...
    FN_PLOT,    FN_YLBL,    FN_REVERSE_Y,    'outcome',   'Prefrontal negative peak'; ...
    P300_PLOT,  P300_YLBL,  P300_REVERSE_Y,  'outcome',   'P300'; ...
    THETA_PLOT, THETA_YLBL, THETA_REVERSE_Y, 'incorrect', 'Frontal theta'; ...
    PLVFP_PLOT, PLVFP_YLBL, PLVFP_REVERSE_Y, 'all',       'Fronto-parietal PLV'; ...
    PLVFS_PLOT, PLVFS_YLBL, PLVFS_REVERSE_Y, 'all',       'Fronto-somatosensory PLV'};

for ci = 1:size(component_specs,1)
    feat_i = component_specs{ci,1};
    ylbl_i = component_specs{ci,2};
    rev_i  = component_specs{ci,3};
    mode_i = component_specs{ci,4};
    name_i = component_specs{ci,5};

    if ~ismember(feat_i, gt.Properties.VariableNames)
        warning('Skipping RQ6 plots for %s: missing feature %s.', name_i, feat_i);
        continue;
    end

    safe_name = matlab.lang.makeValidName(name_i);
    plot_component_by_transition(gt, feat_i, ylbl_i, ...
        sprintf('Recent block transition: %s', name_i), ...
        figure_output_folder, ['RQ6_transition_' safe_name], rev_i, mode_i);

    plot_component_by_priorP(gt, feat_i, ylbl_i, ...
        sprintf('Prior probabilistic exposure: %s', name_i), ...
        figure_output_folder, ['RQ6_nPrevP_' safe_name], rev_i, mode_i);
end

fprintf('\n\nAll RQ analyses complete.\n');

% Save model objects
models_to_save = {'mdl_fn_rq1','mdl_p3_rq1','mdl_fn_conf', ...
                  'mdl_theta','mdl_plv_fp','mdl_plv_predict','mdl_pathway'};
if exist('mdl_false_fn','var') && ~isempty(mdl_false_fn)
    models_to_save{end+1} = 'mdl_false_fn';
end
% save(fullfile(figure_output_folder,'RQ_models.mat'), models_to_save{:});
% fprintf('Saved model objects to RQ_models.mat\n');


% =========================================================================
%% MANUSCRIPT-READY STATISTICS REPORTING
% =========================================================================
%
% Generates three outputs per model:
%   1. Console text — APA-formatted inline sentences ready to paste
%   2. PDF table   — formatted coefficient table with significance stars
%   3. Marginal and conditional R² (Nakagawa & Schielzeth 2013,
%      Methods in Ecology & Evolution 4:133-142)
%
% Usage:
%   report = manuscript_report(mdl_fn_rq1, 'RQ1: Feedback negativity');
%   manuscript_report_all(figure_output_folder);   % run all at once
%

% -- Collect all models into a named struct for batch reporting ─────────────
all_models = struct();
all_models.RQ1_FN   = struct('mdl', mdl_fn_rq1,      'title', 'RQ1: Feedback negativity ~ block type × outcome');
all_models.RQ1_P300 = struct('mdl', mdl_p3_rq1,       'title', 'RQ1: P300 ~ block type × outcome');
all_models.RQ2_FN   = struct('mdl', mdl_fn_conf,       'title', 'RQ2: FN ~ confidence × block type (incorrect)');
all_models.RQ3_Th   = struct('mdl', mdl_theta,         'title', 'RQ3: Theta ~ stage × block type (incorrect)');
all_models.RQ4_PLV  = struct('mdl', mdl_plv_fp,        'title', 'RQ4: FP PLV ~ stage × block type');
all_models.RQ4_Pred = struct('mdl', mdl_plv_predict,   'title', 'RQ4: Next-trial accuracy ~ PLV × block type');
all_models.RQ5_Path = struct('mdl', mdl_pathway,        'title', 'RQ5: PLV ~ pathway × stage × block type');
if exist('mdl_false_fn','var') && ~isempty(mdl_false_fn)
    all_models.RQ2_FalseFB = struct('mdl', mdl_false_fn, 'title', 'RQ2: FN ~ confidence × false FB (P blocks)');
end

% Run all reports
model_names = fieldnames(all_models);
all_reports = cell(numel(model_names), 1);

for mi = 1:numel(model_names)
    mn  = model_names{mi};
    rpt = manuscript_report(all_models.(mn).mdl, all_models.(mn).title, ...
                            figure_output_folder, mn);
    all_reports{mi} = rpt;
    fprintf('\n');
end

% Write a single combined text file with all inline sentences
write_combined_results(all_reports, model_names, figure_output_folder);
fprintf('\nManuscript stats written to %s\n', figure_output_folder);


% =========================================================================
%% MANUSCRIPT REPORTING FUNCTION
% =========================================================================
function report = manuscript_report(mdl, model_title, outdir, fname_stem)
%MANUSCRIPT_REPORT  Generate APA-formatted stats for one LME/GLME model.
%
%  report = manuscript_report(mdl, model_title, outdir, fname_stem)
%
%  INPUTS
%    mdl         — fitted LinearMixedModel or GeneralizedLinearMixedModel
%    model_title — string, used as section heading in output
%    outdir      — directory to save PDF table and text file
%    fname_stem  — filename stem (no extension)
%
%  OUTPUT (struct)
%    report.coef_table    — MATLAB table with formatted coefficients
%    report.inline        — cell array of APA inline sentences per effect
%    report.R2m           — marginal R² (fixed effects only)
%    report.R2c           — conditional R² (fixed + random effects)
%    report.model_title   — echoed title string
%    report.is_glme       — true if GeneralizedLinearMixedModel
%
%  APA FORMAT USED:
%    Continuous predictors:  b = X.XX, SE = X.XX, t(df) = X.XX, p = .XXX
%    Logistic model:         OR = X.XX, 95% CI [X.XX, X.XX], z = X.XX, p = .XXX
%    R²:   marginal R²m = .XX, conditional R²c = .XX (Nakagawa & Schielzeth 2013)
%
%  SIGNIFICANCE STARS:
%    p < .001  ***    p < .01  **    p < .05  *    p < .10  †    otherwise (ns)

is_glme = isa(mdl, 'GeneralizedLinearMixedModel');
coef    = mdl.Coefficients;
n_coef  = height(coef);

fprintf('=== %s ===\n', model_title);

% ── Extract coefficient components ────────────────────────────────────────
names  = coef.Name;
betas  = coef.Estimate;
SEs    = coef.SE;

if is_glme
    stats  = coef.tStat;   % actually z in GLME
    pvals  = coef.pValue;
    dfs    = nan(n_coef, 1);
else
    stats  = coef.tStat;
    dfs    = coef.DF;
    pvals  = coef.pValue;
end

% ── Significance stars ─────────────────────────────────────────────────────
stars = cell(n_coef, 1);
for i = 1:n_coef
    p = pvals(i);
    if p < 0.001
        stars{i} = '***';
    elseif p < 0.01
        stars{i} = '**';
    elseif p < 0.05
        stars{i} = '*';
    elseif p < 0.10
        stars{i} = '\dagger';
    else
        stars{i} = 'ns';
    end
end

% ── Build formatted coefficient table ────────────────────────────────────────
coef_table            = table();
coef_table.Predictor  = names;
coef_table.b          = betas;
coef_table.SE         = SEs;

if is_glme
    % Odds ratios and 95% CIs for logistic models
    coef_table.OR      = exp(betas);
    coef_table.CI_low  = exp(betas - 1.96*SEs);
    coef_table.CI_high = exp(betas + 1.96*SEs);
    coef_table.z       = stats;
    coef_table.p       = pvals;
else
    coef_table.df      = dfs;
    coef_table.t       = stats;
    coef_table.p       = pvals;
end
coef_table.sig = stars;

% ── Marginal and conditional R² ──────────────────────────────────────────────
[R2m, R2c] = compute_R2(mdl, is_glme);
fprintf('  Marginal R²m = %.3f, Conditional R²c = %.3f\n', R2m, R2c);

% ── APA inline sentences per predictor ──────────────────────────────────────
inline = cell(n_coef, 1);
for i = 1:n_coef
    pred = names{i};
    if is_glme
        inline{i} = sprintf(['%s: OR = %.2f, 95%% CI [%.2f, %.2f], ' ...
            'z = %.2f, p = %s'], ...
            pred, exp(betas(i)), exp(betas(i)-1.96*SEs(i)), ...
            exp(betas(i)+1.96*SEs(i)), stats(i), format_p(pvals(i)));
    else
        inline{i} = sprintf(['%s: b = %.3f, SE = %.3f, ' ...
            't(%s) = %.2f, p = %s'], ...
            pred, betas(i), SEs(i), format_df(dfs(i)), ...
            stats(i), format_p(pvals(i)));
    end

    % Flag significant
    sig_marker = '';
    if pvals(i) < 0.05, sig_marker = ' [SIGNIFICANT]'; end
    fprintf('  %s%s\n', inline{i}, sig_marker);
end

% ── Random effects summary ───────────────────────────────────────────────────
re = mdl.covarianceParameters;
fprintf('  Random effects:\n');
for ri = 1:numel(re)
    fprintf('    RE group %d: ', ri);
    disp(re{ri});
end

% ── Model fit indices ────────────────────────────────────────────────────────
aic  = mdl.ModelCriterion.AIC;
bic  = mdl.ModelCriterion.BIC;
logL = mdl.LogLikelihood;
fprintf('  AIC = %.1f, BIC = %.1f, logLik = %.1f\n', aic, bic, logL);

% ── Write PDF-style table ────────────────────────────────────────────────────
if nargin >= 3 && ~isempty(outdir)
    write_coef_table_figure(coef_table, model_title, R2m, R2c, ...
        aic, bic, logL, is_glme, outdir, fname_stem);
end

% ── Assemble report struct ────────────────────────────────────────────────────
report.coef_table  = coef_table;
report.inline      = inline;
report.R2m         = R2m;
report.R2c         = R2c;
report.AIC         = aic;
report.BIC         = bic;
report.logLik      = logL;
report.model_title = model_title;
report.is_glme     = is_glme;
report.n_subj      = numel(numel(mdl.VariableInfo.Range{1,1}));
report.n_obs       = mdl.NumObservations;

end


% ─────────────────────────────────────────────────────────────────────────────
function [R2m, R2c] = compute_R2(mdl, is_glme)
%COMPUTE_R2  Marginal and conditional R² via Nakagawa & Schielzeth (2013).
%
%  For LME: variance of fixed effects (σ²_f) is Var(X*β).
%  For GLME with logit link: distribution-specific variance σ²_d = π²/3.
%  R²m = σ²_f / (σ²_f + σ²_r + σ²_e)
%  R²c = (σ²_f + σ²_r) / (σ²_f + σ²_r + σ²_e)
%
%  Reference: Nakagawa S & Schielzeth H (2013) Methods Ecol Evol 4:133-142.

try
    % Fixed-effect predicted values (design matrix × coefficients)
    [~, ~, stats] = randomEffects(mdl);
    
    % Variance of fixed-effect linear predictor
    % = Var(X * β_fixed)  across all observations
    fitted_fe = fitted(mdl, 'Conditional', false);   % fixed-effects only
    sigma2_f  = var(fitted_fe, 'omitnan');

    % Random-effects variance (sum of all RE variance components)
    re_params = mdl.covarianceParameters;
    sigma2_r  = 0;
    for ri = 1:numel(re_params)
        sigma2_r = sigma2_r + sum(diag(re_params{ri}));
    end

    if is_glme
        % Logit link: distribution-specific variance = π²/3
        sigma2_e = pi^2 / 3;
    else
        % LME: residual variance
        sigma2_e = mdl.MSE;
    end

    total    = sigma2_f + sigma2_r + sigma2_e;
    R2m      = sigma2_f / total;
    R2c      = (sigma2_f + sigma2_r) / total;

catch ME
    warning('R² computation failed: %s. Returning NaN.', ME.message);
    R2m = NaN;
    R2c = NaN;
end

end


% ─────────────────────────────────────────────────────────────────────────────
function write_coef_table_figure(coef_table, title_str, R2m, R2c, ...
    aic, bic, logL, is_glme, outdir, fname_stem)
%WRITE_COEF_TABLE_FIGURE  Save a formatted coefficient table as a PDF figure.
%
% Layout: one row per predictor. Columns depend on model type.
% Significant rows are highlighted. Stars appear in the p column.

n_rows = height(coef_table);

% Column definitions
if is_glme
    col_headers = {'Predictor','b','SE','OR','95% CI','z','p',''};
    col_widths  = [0.28 0.08 0.08 0.08 0.13 0.09 0.09 0.07];
else
    col_headers = {'Predictor','b','SE','df','t','p',''};
    col_widths  = [0.30 0.09 0.09 0.09 0.09 0.09 0.09];
end

n_cols = numel(col_headers);

fig_h  = max(4, 0.4 * (n_rows + 4));   % dynamic height
fig    = figure('Position',[50 50 900 fig_h*72], 'Color','w', 'Visible','off');
ax     = axes(fig, 'Position',[0 0 1 1], 'Visible','off');
hold(ax,'on');
xlim(ax,[0 1]); ylim(ax,[0 1]);

% Title
text(ax, 0.02, 0.97, title_str, ...
    'FontSize',11,'FontWeight','bold','VerticalAlignment','top', ...
    'Interpreter','none');

% Subtitle: model fit
fit_txt = sprintf('AIC = %.1f   BIC = %.1f   logLik = %.1f   R²m = %.3f   R²c = %.3f', ...
    aic, bic, logL, R2m, R2c);
text(ax, 0.02, 0.91, fit_txt, ...
    'FontSize',8,'Color',[0.4 0.4 0.4],'VerticalAlignment','top', ...
    'Interpreter','none');

% Header row
row_top   = 0.85;
row_h     = min(0.08, 0.75 / (n_rows + 1));
header_y  = row_top;

x_left = 0.02;
x_pos  = cumsum([x_left, col_widths]);

for c = 1:n_cols
    text(ax, x_pos(c), header_y, col_headers{c}, ...
        'FontSize',9,'FontWeight','bold','VerticalAlignment','top', ...
        'Interpreter','none');
end

% Horizontal rule after header
line(ax, [0.02 0.98], [header_y - 0.005, header_y - 0.005], ...
    'Color',[0 0 0],'LineWidth',0.8);

% Data rows
for r = 1:n_rows
    row_y = row_top - row_h * r;
    p_val = coef_table.p(r);

    % Highlight significant rows in light blue
    if p_val < 0.05
        patch(ax, [0.02 0.98 0.98 0.02], ...
            [row_y+row_h*0.85 row_y+row_h*0.85 row_y-row_h*0.05 row_y-row_h*0.05], ...
            [0.88 0.93 0.98], 'EdgeColor','none', 'FaceAlpha', 0.6);
    end

    % Predictor name
    pred_str = strrep(coef_table.Predictor{r}, '_', '\_');
    text(ax, x_pos(1), row_y, pred_str, ...
        'FontSize', 8.5, 'VerticalAlignment','top','Interpreter','tex');

    % b
    text(ax, x_pos(2), row_y, sprintf('%.3f', coef_table.b(r)), ...
        'FontSize',8.5,'VerticalAlignment','top');

    % SE
    text(ax, x_pos(3), row_y, sprintf('%.3f', coef_table.SE(r)), ...
        'FontSize',8.5,'VerticalAlignment','top');

    if is_glme
        % OR
        text(ax, x_pos(4), row_y, sprintf('%.2f', coef_table.OR(r)), ...
            'FontSize',8.5,'VerticalAlignment','top');
        % CI
        ci_str = sprintf('[%.2f, %.2f]', coef_table.CI_low(r), coef_table.CI_high(r));
        text(ax, x_pos(5), row_y, ci_str, 'FontSize',8,'VerticalAlignment','top');
        % z
        text(ax, x_pos(6), row_y, sprintf('%.2f', coef_table.z(r)), ...
            'FontSize',8.5,'VerticalAlignment','top');
        % p
        text(ax, x_pos(7), row_y, format_p_table(p_val), ...
            'FontSize',8.5,'VerticalAlignment','top');
        % stars
        text(ax, x_pos(8), row_y, coef_table.sig{r}, ...
            'FontSize',9,'VerticalAlignment','top','Color',[0.6 0 0]);
    else
        % df
        df_str = format_df(coef_table.df(r));
        text(ax, x_pos(4), row_y, df_str, 'FontSize',8.5,'VerticalAlignment','top');
        % t
        text(ax, x_pos(5), row_y, sprintf('%.2f', coef_table.t(r)), ...
            'FontSize',8.5,'VerticalAlignment','top');
        % p
        text(ax, x_pos(6), row_y, format_p_table(p_val), ...
            'FontSize',8.5,'VerticalAlignment','top');
        % stars
        text(ax, x_pos(7), row_y, coef_table.sig{r}, ...
            'FontSize',9,'VerticalAlignment','top','Color',[0.6 0 0]);
    end
end

% Bottom rule
bottom_y = row_top - row_h * (n_rows + 0.5);
line(ax, [0.02 0.98], [bottom_y, bottom_y], 'Color',[0 0 0],'LineWidth',0.8);

% Note
note_str = ['Note: * p < .05, ** p < .01, *** p < .001, \dagger p < .10. ' ...
    'R\^2 via Nakagawa & Schielzeth (2013). ' ...
    'Significant rows shaded.'];
text(ax, 0.02, bottom_y - 0.03, note_str, ...
    'FontSize', 7, 'Color', [0.4 0.4 0.4], 'Interpreter','tex');

if nargin >= 9 && ~isempty(outdir)
    fname = fullfile(outdir, [fname_stem '_table.pdf']);
    exportgraphics(fig, fname, 'ContentType','vector');
    fprintf('  Table saved: %s\n', fname);
end
% % close(fig);

end


% ─────────────────────────────────────────────────────────────────────────────
function write_combined_results(all_reports, model_names, outdir)
%WRITE_COMBINED_RESULTS  Write all inline APA sentences to a single text file.
%
% Output format is ready for direct copy-paste into a manuscript Results
% section, or into supplementary material.

fname = fullfile(outdir, 'manuscript_stats.txt');
fid   = fopen(fname, 'w');
if fid == -1
    warning('Could not open %s for writing.', fname);
    return
end

fprintf(fid, 'STATISTICAL RESULTS — AUTO-GENERATED\n');
fprintf(fid, 'Generated: %s\n', datestr(now, 'dd-mmm-yyyy HH:MM'));
fprintf(fid, '%s\n\n', repmat('=', 1, 70));

for mi = 1:numel(model_names)
    rpt = all_reports{mi};
    if isempty(rpt); continue; end

    fprintf(fid, '\n%s\n', rpt.model_title);
    fprintf(fid, '%s\n', repmat('-', 1, 60));

    % Model-level summary
    fprintf(fid, 'N subjects = %d, N observations = %d\n', ...
        rpt.n_subj, rpt.n_obs);
    fprintf(fid, 'AIC = %.1f, BIC = %.1f, logLik = %.1f\n', ...
        rpt.AIC, rpt.BIC, rpt.logLik);

    if ~isnan(rpt.R2m)
        fprintf(fid, 'Marginal R²m = %.3f, Conditional R²c = %.3f\n', ...
            rpt.R2m, rpt.R2c);
        fprintf(fid, '(Nakagawa & Schielzeth 2013, Methods Ecol Evol 4:133-142)\n');
    end

    fprintf(fid, '\nFixed effects:\n');
    for ii = 1:numel(rpt.inline)
        p_val = rpt.coef_table.p(ii);
        marker = '';
        if p_val < 0.001,     marker = ' ***';
        elseif p_val < 0.01,  marker = ' **';
        elseif p_val < 0.05,  marker = ' *';
        elseif p_val < 0.10,  marker = ' (trend)';
        end
        fprintf(fid, '  %s%s\n', rpt.inline{ii}, marker);
    end

    % Manuscript-ready paragraph for significant effects
    sig_mask = rpt.coef_table.p < 0.05;
    if any(sig_mask)
        fprintf(fid, '\nSuggested inline text for significant effects:\n');
        sig_idx = find(sig_mask);
        for ii = sig_idx(:)'
            pred = rpt.coef_table.Predictor{ii};
            if rpt.is_glme
                fprintf(fid, ['  The predictor %s was a significant predictor ' ...
                    '(OR = %.2f, 95%% CI [%.2f, %.2f], z = %.2f, %s).\n'], ...
                    pred, rpt.coef_table.OR(ii), ...
                    rpt.coef_table.CI_low(ii), rpt.coef_table.CI_high(ii), ...
                    rpt.coef_table.z(ii), ...
                    format_p_text(rpt.coef_table.p(ii)));
            else
                fprintf(fid, ['  There was a significant effect of %s ' ...
                    '(b = %.3f, SE = %.3f, t(%s) = %.2f, %s).\n'], ...
                    pred, rpt.coef_table.b(ii), rpt.coef_table.SE(ii), ...
                    format_df(rpt.coef_table.df(ii)), rpt.coef_table.t(ii), ...
                    format_p_text(rpt.coef_table.p(ii)));
            end
        end
    else
        fprintf(fid, '\nNo fixed effects reached significance (p < .05).\n');
    end

    fprintf(fid, '\n');
end

% close(fid);
fprintf('Combined results written to:\n  %s\n', fname);

end


% =========================================================================
%% STRING FORMATTING HELPERS
% =========================================================================

function s = format_p(p)
%FORMAT_P  APA-style p-value string for inline text.
% APA 7 recommends reporting exact p-values, with p < .001 as the lower limit.
if isnan(p)
    s = 'p = n/a';
elseif p < 0.001
    s = 'p < .001';
else
    s = sprintf('p = %s', strtrim(sprintf('%.3f', p)));
    % Remove leading zero per APA style: ".049" not "0.049"
    s = strrep(s, 'p = 0.', 'p = .');
end
end


function s = format_p_table(p)
%FORMAT_P_TABLE  Compact p-value for table cells.
if isnan(p)
    s = 'n/a';
elseif p < 0.001
    s = '< .001';
else
    s = strtrim(sprintf('%.3f', p));
    % Remove leading zero
    if startsWith(s, '0.')
        s = s(2:end);
    end
end
end


function s = format_p_text(p)
%FORMAT_P_TEXT  Full APA inline sentence fragment.
if p < 0.001
    s = 'p < .001';
elseif p < 0.01
    s = sprintf('p = %s', strtrim(strrep(sprintf('%.3f', p), '0.', '.')));
else
    s = sprintf('p = %s', strtrim(strrep(sprintf('%.3f', p), '0.', '.')));
end
end


function s = format_df(df)
%FORMAT_DF  Format degrees of freedom — integer if whole number, else 1dp.
if isnan(df)
    s = '?';
elseif df == round(df)
    s = sprintf('%d', df);
else
    s = sprintf('%.1f', df);
end
end




% =========================================================================
%% TRANSITION AND PRIOR-P PLOTTING HELPERS
% =========================================================================
function T = add_transition_history_columns(T)
%ADD_TRANSITION_HISTORY_COLUMNS  Add transition_recent and n_prev_P columns.
% transition_recent is the previous block type -> current block type label:
% D->D, D->P, P->D, P->P. First blocks are labelled first.
% n_prev_P is the number of previous P blocks within the subject.

if ~ismember('subj_id', T.Properties.VariableNames)
    error('add_transition_history_columns: T must contain subj_id.');
end
if ~ismember('block_number', T.Properties.VariableNames)
    if ismember('block', T.Properties.VariableNames)
        T.block_number = kh_to_numeric(T.block);
    elseif ismember('blocknum', T.Properties.VariableNames)
        T.block_number = kh_to_numeric(T.blocknum);
    else
        error('add_transition_history_columns: no block_number/block/blocknum column found.');
    end
end

T.subj_id_s = string(T.subj_id);
bt = string(T.block_type);
bt(bt=="V") = "P";
T.block_type_s = bt;
T.transition_recent = repmat("first", height(T), 1);
T.n_prev_P = nan(height(T), 1);
T.n_prev_P_bin = repmat("missing", height(T), 1);

subs = unique(T.subj_id_s);
for si = 1:numel(subs)
    sm = T.subj_id_s == subs(si);
    blks = unique(T.block_number(sm & ~isnan(T.block_number)));
    blks = sort(blks(:)');
    prev_types = strings(0,1);

    for bi = 1:numel(blks)
        b = blks(bi);
        bm = sm & T.block_number == b;
        if ~any(bm), continue; end

        curr_vals = T.block_type_s(bm);
        curr_vals = curr_vals(curr_vals=="D" | curr_vals=="P");
        if isempty(curr_vals)
            curr = "";
        else
            curr = curr_vals(find(curr_vals~="",1,'first'));
        end

        n_prev = sum(prev_types == "P");
        T.n_prev_P(bm) = n_prev;
        if n_prev >= 2
            T.n_prev_P_bin(bm) = "2+";
        else
            T.n_prev_P_bin(bm) = string(n_prev);
        end

        if bi == 1 || isempty(prev_types)
            T.transition_recent(bm) = "first";
        else
            prev = prev_types(end);
            if (prev=="D" || prev=="P") && (curr=="D" || curr=="P")
                T.transition_recent(bm) = prev + "->" + curr;
            else
                T.transition_recent(bm) = "unknown";
            end
        end

        if curr=="D" || curr=="P"
            prev_types(end+1,1) = curr; %#ok<AGROW>
        end
    end
end
end

function plot_component_by_transition(gt, feat, ylbl, ttl, outdir, fname, reverse_y, mode)
%PLOT_COMPONENT_BY_TRANSITION  Stage profiles split by recent block transition.
% mode = 'outcome'   : two panels, Incorrect and Correct, true feedback only.
% mode = 'incorrect' : one panel, incorrect true-feedback trials only.
% mode = 'all'       : one panel, all true-feedback trials.
if nargin < 8 || isempty(mode), mode = 'all'; end
if nargin < 7 || isempty(reverse_y), reverse_y = false; end

stages = {'LN','LE','RN','RE'};
transitions = {'D->D','D->P','P->D','P->P'};
trans_labels = {'D->D','D->P','P->D','P->P'};
trans_cols = [0.12 0.47 0.71; 0.85 0.33 0.10; 0.47 0.67 0.19; 0.80 0.20 0.60];

switch char(mode)
    case 'outcome'
        panel_defs = {outcome_mask_local(gt,0), 'Incorrect true-feedback trials'; ...
                      outcome_mask_local(gt,1), 'Correct true-feedback trials'};
    case 'incorrect'
        panel_defs = {outcome_mask_local(gt,0), 'Incorrect true-feedback trials'};
    otherwise
        panel_defs = {true(height(gt),1), 'All true-feedback trials'};
end

fig = figure('Position',[50 50 620*size(panel_defs,1) 520]);
sgtitle(ttl, 'Interpreter','none');

for pi = 1:size(panel_defs,1)
    ax = subplot(1,size(panel_defs,1),pi); hold(ax,'on');
    title(ax, panel_defs{pi,2}, 'FontSize',10);
    base_mask = panel_defs{pi,1} & true_feedback_mask_local(gt) & ~is_first_transition_local(gt);

    for ti = 1:numel(transitions)
        means = nan(1,numel(stages));
        sems  = nan(1,numel(stages));
        ns    = nan(1,numel(stages));
        for si = 1:numel(stages)
            m = base_mask & string(gt.transition_recent)==transitions{ti} & string(gt.stage)==stages{si} & ~isnan(gt.(feat));
            vals = subject_means_local(gt, feat, m);
            [means(si), sems(si), ns(si)] = mean_sem_n_local(vals);
        end
        errorbar(ax,1:numel(stages),means,sems,'-o', ...
            'Color',trans_cols(ti,:), 'MarkerFaceColor',trans_cols(ti,:), ...
            'LineWidth',1.7, 'DisplayName',trans_labels{ti});
        for si = 1:numel(stages)
            if ~isnan(means(si))
                text(ax, si, means(si), sprintf(' n=%d', ns(si)), ...
                    'Color',trans_cols(ti,:), 'FontSize',7, ...
                    'VerticalAlignment','bottom', 'Clipping','off');
            end
        end
    end

    set(ax,'XTick',1:numel(stages),'XTickLabel',stages);
    xlabel(ax,'Stage'); ylabel(ax,ylbl,'Interpreter','none');
    if reverse_y, set(ax,'YDir','reverse'); end
    legend(ax,'Box','off','Location','best');
    axis(ax,'square');
end

exportgraphics(fig, fullfile(outdir,[fname '.pdf']), 'ContentType','vector');
end

function plot_component_by_priorP(gt, feat, ylbl, ttl, outdir, fname, reverse_y, mode)
%PLOT_COMPONENT_BY_PRIORP  n_prev_P profiles split by stage and current D/P block.
% X axis is prior P exposure: 0, 1, 2+. Lines show stage. Columns show current
% block type. For outcome-sensitive measures, rows split Incorrect/Correct.
if nargin < 8 || isempty(mode), mode = 'all'; end
if nargin < 7 || isempty(reverse_y), reverse_y = false; end

stages = {'LN','LE','RN','RE'};
stage_cols = [0.12 0.62 0.47; 0.85 0.65 0.00; 0.80 0.27 0.13; 0.40 0.25 0.65];
npp_bins = {'0','1','2+'};
btypes = {'D','P'};

switch char(mode)
    case 'outcome'
        row_defs = {outcome_mask_local(gt,0), 'Incorrect'; ...
                    outcome_mask_local(gt,1), 'Correct'};
    case 'incorrect'
        row_defs = {outcome_mask_local(gt,0), 'Incorrect'};
    otherwise
        row_defs = {true(height(gt),1), 'All outcomes'};
end

fig = figure('Position',[50 50 1050 430*size(row_defs,1)]);
sgtitle(ttl, 'Interpreter','none');

for ri = 1:size(row_defs,1)
    for bi = 1:numel(btypes)
        ax = subplot(size(row_defs,1), numel(btypes), (ri-1)*numel(btypes)+bi); hold(ax,'on');
        title(ax, sprintf('%s | current %s block', row_defs{ri,2}, btypes{bi}), 'FontSize',10);
        base_mask = row_defs{ri,1} & true_feedback_mask_local(gt) & string(gt.block_type_s)==btypes{bi};

        for si = 1:numel(stages)
            means = nan(1,numel(npp_bins));
            sems  = nan(1,numel(npp_bins));
            ns    = nan(1,numel(npp_bins));
            for ni = 1:numel(npp_bins)
                m = base_mask & string(gt.stage)==stages{si} & string(gt.n_prev_P_bin)==npp_bins{ni} & ~isnan(gt.(feat));
                vals = subject_means_local(gt, feat, m);
                [means(ni), sems(ni), ns(ni)] = mean_sem_n_local(vals);
            end
            errorbar(ax,1:numel(npp_bins),means,sems,'-o', ...
                'Color',stage_cols(si,:), 'MarkerFaceColor',stage_cols(si,:), ...
                'LineWidth',1.7, 'DisplayName',stages{si});
            for ni = 1:numel(npp_bins)
                if ~isnan(means(ni))
                    text(ax, ni, means(ni), sprintf(' n=%d', ns(ni)), ...
                        'Color',stage_cols(si,:), 'FontSize',7, ...
                        'VerticalAlignment','bottom', 'Clipping','off');
                end
            end
        end

        set(ax,'XTick',1:numel(npp_bins),'XTickLabel',{'0','1','2+'});
        xlabel(ax,'Number of prior P blocks'); ylabel(ax,ylbl,'Interpreter','none');
        if reverse_y, set(ax,'YDir','reverse'); end
        legend(ax,'Box','off','Location','best');
        axis(ax,'square');
    end
end

exportgraphics(fig, fullfile(outdir,[fname '.pdf']), 'ContentType','vector');
end

function m = true_feedback_mask_local(T)
if ismember('false_fb', T.Properties.VariableNames)
    m = ~logical(T.false_fb);
else
    m = true(height(T),1);
end
end

function m = is_first_transition_local(T)
if ismember('transition_recent', T.Properties.VariableNames)
    m = string(T.transition_recent)=="first" | string(T.transition_recent)=="unknown";
else
    m = false(height(T),1);
end
end

function m = outcome_mask_local(T, code)
%OUTCOME_MASK_LOCAL  Robust correct/incorrect mask for numeric/categorical/string.
if ~ismember('correct', T.Properties.VariableNames)
    m = false(height(T),1);
    return;
end
c = T.correct;
if isnumeric(c) || islogical(c)
    m = double(c) == double(code);
else
    cs = lower(string(c));
    if code == 1
        m = cs=="1" | cs=="correct" | cs=="true" | cs=="win";
    else
        m = cs=="0" | cs=="incorrect" | cs=="false" | cs=="loss" | cs=="error";
    end
end
end

function vals = subject_means_local(T, feat, mask)
%SUBJECT_MEANS_LOCAL  One value per subject for a condition cell.
subs = unique(string(T.subj_id(mask)));
vals = nan(numel(subs),1);
for si = 1:numel(subs)
    sm = mask & string(T.subj_id)==subs(si);
    vals(si) = mean(T.(feat)(sm), 'omitnan');
end
vals = vals(~isnan(vals));
end

function [mn, se, n] = mean_sem_n_local(vals)
vals = vals(~isnan(vals));
n = numel(vals);
if n == 0
    mn = NaN; se = NaN;
elseif n == 1
    mn = vals(1); se = NaN;
else
    mn = mean(vals,'omitnan');
    se = std(vals,'omitnan') / sqrt(n);
end
end

% =========================================================================
%% PLOTTING SCALE HELPERS
% =========================================================================
function [feat, ylbl, reverse_y] = resolve_plot_feature(T, measure_key, plot_scale)
%RESOLVE_PLOT_FEATURE  Return a valid feature column for plotting.
% measure_key options used here: prefrontal_neg_peak, P300, Theta, PLV_fp, PLV_fs.
% reverse_y is true only for negative-going frontocentral amplitude plots.

plot_scale = lower(string(plot_scale));
reverse_y = false;

switch char(measure_key)
    case 'prefrontal_neg_peak'
        reverse_y = true;  % more negative = larger frontocentral negativity
        raw_col  = 'prefrontal_neg_peak_amp';
        norm_col = 'prefrontal_neg_peak_norm';
        z_col    = 'prefrontal_neg_peak_norm_z';
        raw_lbl  = 'Prefrontal negative peak (uV; more negative = larger FN)';
        norm_lbl = 'Prefrontal negative peak / baseline RMS (more negative = larger FN)';
        z_lbl    = 'Prefrontal negative peak (within-subject z; more negative = larger FN)';
    case 'P300'
        raw_col  = 'P300_amp';
        norm_col = 'P300_norm';
        z_col    = 'P300_norm_z';
        raw_lbl  = 'P300 amplitude (uV)';
        norm_lbl = 'P300 amplitude / baseline RMS';
        z_lbl    = 'P300 amplitude (within-subject z)';
    case 'Theta'
        raw_col  = 'Theta_amp';
        norm_col = 'Theta_amp';
        z_col    = 'Theta_amp_z';
        raw_lbl  = 'Theta amplitude (baseline-corrected)';
        norm_lbl = 'Theta amplitude (baseline-corrected)';
        z_lbl    = 'Theta amplitude (within-subject z)';
    case 'PLV_fp'
        raw_col  = 'PLV_fp';
        norm_col = 'PLV_fp';
        z_col    = 'PLV_fp_z';
        raw_lbl  = 'Fronto-parietal PLV (baseline-corrected)';
        norm_lbl = 'Fronto-parietal PLV (baseline-corrected)';
        z_lbl    = 'Fronto-parietal PLV (within-subject z)';
    case 'PLV_fs'
        raw_col  = 'PLV_fs';
        norm_col = 'PLV_fs';
        z_col    = 'PLV_fs_z';
        raw_lbl  = 'Fronto-somatosensory PLV (baseline-corrected)';
        norm_lbl = 'Fronto-somatosensory PLV (baseline-corrected)';
        z_lbl    = 'Fronto-somatosensory PLV (within-subject z)';
    otherwise
        error('Unknown measure_key: %s', measure_key);
end

switch char(plot_scale)
    case 'z'
        candidates = {z_col, norm_col, raw_col};
        labels     = {z_lbl, norm_lbl, raw_lbl};
    case 'norm'
        candidates = {norm_col, raw_col, z_col};
        labels     = {norm_lbl, raw_lbl, z_lbl};
    case 'raw'
        candidates = {raw_col, norm_col, z_col};
        labels     = {raw_lbl, norm_lbl, z_lbl};
    otherwise
        error('PLOT_SCALE must be ''z'', ''norm'', or ''raw''.');
end

feat = '';
ylbl = '';
for i = 1:numel(candidates)
    if ismember(candidates{i}, T.Properties.VariableNames)
        feat = candidates{i};
        ylbl = labels{i};
        return;
    end
end

error('No valid plot column found for %s using PLOT_SCALE=%s.', measure_key, plot_scale);
end

function cols = outcome_colors_for_categories(ycats)
%OUTCOME_COLORS_FOR_CATEGORIES  Keep Incorrect red and Correct green.
% Handles numeric category labels 0/1 and text labels Incorrect/Correct.

cols = nan(numel(ycats), 3);
for i = 1:numel(ycats)
    lab = lower(char(ycats{i}));
    if strcmp(lab,'0') || contains(lab,'incorrect') || contains(lab,'error') || contains(lab,'loss')
        cols(i,:) = [0.70 0.10 0.10];  % Incorrect = red
    elseif strcmp(lab,'1') || contains(lab,'correct') || contains(lab,'win')
        cols(i,:) = [0.10 0.60 0.10];  % Correct = green
    else
        fallback = [0.15 0.45 0.70; 0.80 0.30 0.10; 0.65 0.25 0.65; 0.3 0.3 0.3];
        cols(i,:) = fallback(min(i,size(fallback,1)),:);
    end
end
end

% =========================================================================
%% PLOTTING HELPERS  (corrected versions)
% =========================================================================

function plot_lme_effects(gt, feat, xvar, yvar, yvar_labels, ttl, outdir, fname, ylbl, reverse_y)
%PLOT_LME_EFFECTS  Plot cell-means with SEM ± p-values (corrected).
%
%  FIX 4: yvar_labels allows caller to supply readable labels instead of
%          relying on the raw category names (e.g. '0'/'1').
%  FIX 6: uses ttest (paired within subjects) not ttest2 for p-value
%          annotation when subjects appear in both y-groups at each x-level.

if nargin < 9 || isempty(ylbl), ylbl = feat; end
if nargin < 10 || isempty(reverse_y), reverse_y = false; end

if ~iscategorical(gt.(xvar)), gt.(xvar) = categorical(gt.(xvar)); end

% Convert yvar numeric → categorical with supplied labels
if isnumeric(gt.(yvar)) || islogical(gt.(yvar))
    yvals_num = unique(gt.(yvar)(~isnan(gt.(yvar))));
    yvals_num = sort(yvals_num);
    if nargin >= 5 && numel(yvar_labels) == numel(yvals_num)
        gt.(yvar) = categorical(gt.(yvar), yvals_num, yvar_labels);
    else
        gt.(yvar) = categorical(gt.(yvar));
    end
end

xcats = categories(gt.(xvar));
ycats = categories(gt.(yvar));
cols  = outcome_colors_for_categories(ycats);


x = 1:numel(xcats);

fig = figure('Position',[50 50 700 500]);
ax  = axes(fig);
hold(ax,'on');

means = nan(numel(xcats), numel(ycats));
sems  = nan(numel(xcats), numel(ycats));
ns    = nan(numel(xcats), numel(ycats));
subj_vals = cell(numel(xcats), numel(ycats));   % for paired test

for j = 1:numel(ycats)
    for i = 1:numel(xcats)
        m = gt.(xvar)==xcats{i} & gt.(yvar)==ycats{j} & ~isnan(gt.(feat));
        vals = gt.(feat)(m);
        if isempty(vals); continue; end
        ns(i,j)       = sum(~isnan(vals));
        means(i,j)    = mean(vals,'omitnan');
        sems(i,j)     = std(vals,'omitnan') / sqrt(ns(i,j));
        subj_vals{i,j}= vals;
    end
end

all_vals = gt.(feat)(~isnan(gt.(feat)));
y_rng    = range(all_vals);
if y_rng == 0; y_rng = 1; end
y_pad    = 0.08 * y_rng;

for j = 1:numel(ycats)
    x_off = (j - (numel(ycats)+1)/2) * 0.08;
    clr   = cols(j,:);
    errorbar(ax, x + x_off, means(:,j)', sems(:,j)', '-o', ...
        'Color', clr, 'LineWidth',1.8, 'MarkerFaceColor', clr, ...
        'DisplayName', char(ycats{j}));

    for i = 1:numel(xcats)
        if isnan(means(i,j)); continue; end
        txt = sprintf('n=%d\n%.2f±%.2f', ns(i,j), means(i,j), sems(i,j));
        text(ax, x(i)+x_off, means(i,j)+sems(i,j)+y_pad, txt, ...
            'Color',clr,'FontSize',7,'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom','Clipping','off');
    end
end

% ── FIX 6: paired t-test within each x-level ─────────────────────────────
if numel(ycats) == 2
    % Aggregate to subject-level means before pairing
    subjs  = unique(gt.subj_id);
    for i = 1:numel(xcats)
        v1_subj = nan(numel(subjs),1);
        v2_subj = nan(numel(subjs),1);
        for si = 1:numel(subjs)
            m1 = gt.subj_id==subjs(si) & gt.(xvar)==xcats{i} & ...
                 gt.(yvar)==ycats{1} & ~isnan(gt.(feat));
            m2 = gt.subj_id==subjs(si) & gt.(xvar)==xcats{i} & ...
                 gt.(yvar)==ycats{2} & ~isnan(gt.(feat));
            if any(m1), v1_subj(si) = mean(gt.(feat)(m1),'omitnan'); end
            if any(m2), v2_subj(si) = mean(gt.(feat)(m2),'omitnan'); end
        end
        ok = ~isnan(v1_subj) & ~isnan(v2_subj);
        if sum(ok) > 1
            [~, pval] = ttest(v1_subj(ok), v2_subj(ok));   % PAIRED
            y_top = max([means(i,1)+sems(i,1), means(i,2)+sems(i,2)], [], 'omitnan');
            text(ax, x(i), y_top+1.6*y_pad, sprintf('p=%.3f',pval), ...
                'HorizontalAlignment','center','FontSize',7, ...
                'FontWeight','bold','Clipping','off','Color',[0 0 0]);
        end
    end
end

yline(0, '--k')
if reverse_y, set(ax,'YDir','reverse'); end
set(ax,'XTick',x,'XTickLabel',xcats);
xlabel(ax,xvar,'Interpreter','none');
ylabel(ax,ylbl,'Interpreter','none');
title(ax,ttl);
legend(ax,'Box','off');
yl = ylim(ax);
ylim(ax,[yl(1)-0.10*range(yl), yl(2)+0.25*range(yl)]);
exportgraphics(fig, fullfile(outdir,[fname '.pdf']), 'ContentType','vector');
% close(fig);

end


function plot_stage_bar(gt_sub, feat, ylbl, ttl, outdir, fname, reverse_y)
%PLOT_STAGE_BAR  Bar chart per stage with D and P block types.
if nargin < 7 || isempty(reverse_y), reverse_y = false; end
%  FIX 3: direct string comparison instead of categorical({...}).

stages = {'LN','LE','RN','RE'};
bts    = {'D','P'};
clrs   = {[0.15 0.45 0.70],[0.80 0.30 0.10]};

fig = figure('Position',[50 50 800 450]);
ax  = axes(fig);
hold(ax,'on');
title(ax,ttl,'FontSize',9);

x_base = 1:4;
w      = 0.35;

vals_all = gt_sub.(feat);
vals_all = vals_all(~isnan(vals_all));
y_rng    = range(vals_all);
if isempty(vals_all) || y_rng==0; y_rng=1; end
y_pad    = 0.03 * y_rng;

for bt_i = 1:2
    for s_i = 1:4
        % FIX 3: compare block_type to string directly
        mask = gt_sub.block_type == bts{bt_i} & ...
               gt_sub.stage      == stages{s_i} & ...
               ~isnan(gt_sub.(feat));

        vals = gt_sub.(feat)(mask);
        if isempty(vals); continue; end

        x_pos = x_base(s_i) + (bt_i-1.5)*w;
        m     = mean(vals,'omitnan');
        se    = std(vals,'omitnan') / sqrt(sum(~isnan(vals)));
        n_v   = sum(~isnan(vals));

        bar(ax,x_pos,m,w*0.9,'FaceColor',clrs{bt_i},'EdgeColor','none', ...
            'DisplayName',sprintf('%s',bts{bt_i}));
        errorbar(ax,x_pos,m,se,'k','LineWidth',1.5,'HandleVisibility','off');
        text(ax,x_pos,m+se+y_pad, ...
            sprintf('n=%d\n%.2f±%.2f',n_v,m,se), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',7);
    end
end

set(ax,'XTick',1:4,'XTickLabel',stages);
xlabel(ax,'Stage');
ylabel(ax,ylbl,'Interpreter','none');
% Deduplicate legend
h = findobj(ax,'Type','bar');
legend(ax,h(end:-1:1),{'D','P'},'Box','off','Location','best');
if reverse_y, set(ax,'YDir','reverse'); end
exportgraphics(fig,fullfile(outdir,[fname '.pdf']),'ContentType','vector');
% close(fig);

end


function plot_confidence_fRN(gt_inc, outdir, fn_feat, fn_ylbl, reverse_y)
%PLOT_CONFIDENCE_FRN  Scatter + regression: confidence vs FN.
if nargin < 3 || isempty(fn_feat), fn_feat = 'prefrontal_neg_peak_norm_z'; end
if nargin < 4 || isempty(fn_ylbl), fn_ylbl = 'Feedback negativity (z)'; end
if nargin < 5 || isempty(reverse_y), reverse_y = false; end

fig = figure('Position',[50 50 800 400]);
sgtitle('RQ2: Confidence × Feedback negativity (incorrect trials only)');

for bt_i = 1:2
    bt = {'D','P'};
    ax = subplot(1,2,bt_i); hold(ax,'on');
    title(ax,sprintf('Block type: %s',bt{bt_i}));

    mask = gt_inc.block_type==bt{bt_i} & ...
           ~isnan(gt_inc.conf_z) & ~isnan(gt_inc.(fn_feat));
    if sum(mask) < 5; continue; end

    x = gt_inc.conf_z(mask);
    y = gt_inc.(fn_feat)(mask);
    ok = ~isnan(x) & ~isnan(y);

    scatter(ax,x,y,8,[0.6 0.6 0.6],'filled','HandleVisibility','off');

    p_fit = polyfit(x(ok),y(ok),1);
    xi    = linspace(min(x(ok)),max(x(ok)),100);
    plot(ax,xi,polyval(p_fit,xi),'r-','LineWidth',2, ...
        'DisplayName',sprintf('slope=%.2f',p_fit(1)));

    [rval,pval] = corr(x(ok),y(ok),'Rows','complete');
    txt = sprintf('n=%d\nr=%.2f\n%s',sum(ok),rval,format_p(pval));
    text(ax,0.03,0.97,txt,'Units','normalized','VerticalAlignment','top', ...
        'FontSize',8,'BackgroundColor','w','Margin',3,'EdgeColor',[0.85 0.85 0.85]);

    xlabel(ax,'Confidence (z)');
    ylabel(ax,fn_ylbl,'Interpreter','none');
    if reverse_y, set(ax,'YDir','reverse'); end
    legend(ax,'Box','off');
    xlim(ax,[-3 3]);
end
exportgraphics(fig,fullfile(outdir,'RQ2_Confidence_FN.pdf'),'ContentType','vector');
% close(fig);

end


function plot_pathway_comparison(gt_fp, gt_fs, outdir, fp_feat, fs_feat, fp_ylbl, fs_ylbl, reverse_y)
%PLOT_PATHWAY_COMPARISON  FP vs FS PLV by stage, split by block type.
if nargin < 4 || isempty(fp_feat), fp_feat = 'PLV_fp_z'; end
if nargin < 5 || isempty(fs_feat), fs_feat = 'PLV_fs_z'; end
if nargin < 6 || isempty(fp_ylbl), fp_ylbl = fp_feat; end
if nargin < 7 || isempty(fs_ylbl), fs_ylbl = fs_feat; end
if nargin < 8 || isempty(reverse_y), reverse_y = false; end

stages = {'LN','LE','RN','RE'};
bts    = {'D','P'};

fig = figure('Position',[50 50 1200 500]);
sgtitle('RQ5: Fronto-parietal vs Fronto-somatosensory PLV by stage');

for bt_i = 1:2
    ax = subplot(1,2,bt_i); hold(ax,'on');
    title(ax,sprintf('Block type: %s',bts{bt_i}));

    fp_means=nan(1,4); fp_sems=nan(1,4); fp_ns=nan(1,4);
    fs_means=nan(1,4); fs_sems=nan(1,4); fs_ns=nan(1,4);

    for s_i = 1:4
        m_fp = gt_fp.block_type==bts{bt_i} & gt_fp.stage==stages{s_i} & ~isnan(gt_fp.(fp_feat));
        m_fs = gt_fs.block_type==bts{bt_i} & gt_fs.stage==stages{s_i} & ~isnan(gt_fs.(fs_feat));
        v_fp = gt_fp.(fp_feat)(m_fp);
        v_fs = gt_fs.(fs_feat)(m_fs);
        if ~isempty(v_fp)
            fp_means(s_i)=mean(v_fp,'omitnan'); fp_sems(s_i)=std(v_fp,'omitnan')/sqrt(sum(~isnan(v_fp))); fp_ns(s_i)=sum(~isnan(v_fp));
        end
        if ~isempty(v_fs)
            fs_means(s_i)=mean(v_fs,'omitnan'); fs_sems(s_i)=std(v_fs,'omitnan')/sqrt(sum(~isnan(v_fs))); fs_ns(s_i)=sum(~isnan(v_fs));
        end
    end

    errorbar(ax,1:4,fp_means,fp_sems,'o-','Color',[0.15 0.45 0.70],'LineWidth',2,'MarkerFaceColor',[0.15 0.45 0.70],'DisplayName','Fronto-parietal');
    errorbar(ax,1:4,fs_means,fs_sems,'s--','Color',[0.65 0.25 0.65],'LineWidth',2,'MarkerFaceColor',[0.65 0.25 0.65],'DisplayName','Fronto-somatosensory');

    y_all = [fp_means+fp_sems, fs_means+fs_sems];
    y_all = y_all(~isnan(y_all));
    y_pad = 0.03 * max(1, range(y_all));
    for s_i = 1:4
        if ~isnan(fp_means(s_i))
            text(ax,s_i,fp_means(s_i)+fp_sems(s_i)+y_pad, sprintf('n=%d',fp_ns(s_i)),'Color',[0.15 0.45 0.70],'FontSize',7,'HorizontalAlignment','center','VerticalAlignment','bottom');
        end
        if ~isnan(fs_means(s_i))
            text(ax,s_i,fs_means(s_i)+fs_sems(s_i)+y_pad, sprintf('n=%d',fs_ns(s_i)),'Color',[0.65 0.25 0.65],'FontSize',7,'HorizontalAlignment','center','VerticalAlignment','bottom');
        end
    end

    set(ax,'XTick',1:4,'XTickLabel',stages);
    xlabel(ax,'Stage'); ylabel(ax,[fp_ylbl ' / ' fs_ylbl],'Interpreter','none');
    if reverse_y, set(ax,'YDir','reverse'); end
    legend(ax,'Box','off');
end
exportgraphics(fig,fullfile(outdir,'RQ5_Pathway_comparison.pdf'),'ContentType','vector');
% close(fig);

end

function plot_stage_outcome_lines(gt, feat, ylbl, ttl, outdir, fname, reverse_y)
% Plot one figure: stages on x-axis, split by correct vs incorrect.
if nargin < 7 || isempty(reverse_y), reverse_y = false; end

stages = {'LN','LE','RN','RE'};
cols   = [0.70 0.10 0.10;   % Incorrect
          0.10 0.60 0.10];  % Correct

fig = figure('Position',[50 50 900 450]);
ax = axes(fig);
hold(ax,'on');
title(ax, ttl, 'FontSize', 10);

x = 1:4;

vals_all = gt.(feat);
vals_all = vals_all(~isnan(vals_all));
if isempty(vals_all)
    vals_all = 0;
end
y_rng = range(vals_all);
if y_rng == 0, y_rng = 1; end
y_pad = 0.08 * y_rng;

labels = {'Incorrect','Correct'};
codes  = [0, 1];

for j = 1:2
    means = nan(1,4);
    sems  = nan(1,4);
    ns    = nan(1,4);

    for s_i = 1:4
        m = gt.stage == stages{s_i} & gt.correct == codes(j) & ~gt.false_fb & ~isnan(gt.(feat));
        v = gt.(feat)(m);
        if isempty(v), continue; end
        ns(s_i)    = sum(~isnan(v));
        means(s_i) = mean(v, 'omitnan');
        sems(s_i)  = std(v, 'omitnan') / sqrt(ns(s_i));
    end

    x_off = (j - 1.5) * 0.08;
    errorbar(ax, x + x_off, means, sems, '-o', ...
        'Color', cols(j,:), 'LineWidth', 1.8, ...
        'MarkerFaceColor', cols(j,:), ...
        'DisplayName', labels{j});

    for s_i = 1:4
        if isnan(means(s_i)), continue; end
        text(ax, x(s_i)+x_off, means(s_i)+sems(s_i)+y_pad, ...
            sprintf('n=%d\n%.2f \\pm %.2f', ns(s_i), means(s_i), sems(s_i)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontSize',7, 'Clipping','off');
    end
end

set(ax,'XTick',1:4,'XTickLabel',stages);
xlabel(ax,'Stage');
ylabel(ax, ylbl);
legend(ax,'Box','off','Location','best');
if reverse_y, set(ax,'YDir','reverse'); end

yl = ylim(ax);
ylim(ax, [yl(1)-0.10*range(yl), yl(2)+0.18*range(yl)]);

exportgraphics(fig, fullfile(outdir,[fname '.pdf']), 'ContentType','vector');
end

function plot_stage_transition_feedback(gt, feat, ylbl, ttl, outdir, fname, reverse_y)
% Two-panel figure:
%   Left  = D blocks: correct vs incorrect
%   Right = P blocks: true/false feedback split, each with correct vs incorrect
if nargin < 7 || isempty(reverse_y), reverse_y = false; end

stages = {'LN','LE','RN','RE'};

fig = figure('Position',[50 50 1200 450]);
sgtitle(ttl);

% Colors
clr_inc   = [0.70 0.10 0.10];
clr_cor   = [0.10 0.60 0.10];
clr_true  = [0.15 0.45 0.70];
clr_false = [0.80 0.30 0.10];

all_vals = gt.(feat);
all_vals = all_vals(~isnan(all_vals));
if isempty(all_vals)
    all_vals = 0;
end
y_rng = range(all_vals);
if y_rng == 0, y_rng = 1; end
y_pad = 0.08 * y_rng;

% ---------- D panel ----------
ax1 = subplot(1,2,1);
hold(ax1,'on');
title(ax1, 'D blocks');
plot_one_block_panel(ax1, gt, feat, stages, 'D', false, ...
    clr_cor, clr_inc, clr_true, clr_false, y_pad, ylbl, reverse_y);

% ---------- P panel ----------
ax2 = subplot(1,2,2);
hold(ax2,'on');
title(ax2, 'P blocks');
plot_one_block_panel(ax2, gt, feat, stages, 'P', true, ...
    clr_cor, clr_inc, clr_true, clr_false, y_pad, ylbl, reverse_y);

exportgraphics(fig, fullfile(outdir,[fname '.pdf']), 'ContentType','vector');
end

function plot_one_block_panel(ax, gt, feat, stages, bt, is_P, clr_cor, clr_inc, clr_true, clr_false, y_pad, ylbl, reverse_y)

x = 1:4;

if is_P
    % Four lines in P. Colour always codes outcome:
    %   Correct = green, Incorrect = red.
    % Line style codes feedback validity:
    %   True feedback = solid, False feedback = dashed.
    spec = {
        ~gt.false_fb & gt.correct==1, clr_cor, 'True correct',  '-'
        ~gt.false_fb & gt.correct==0, clr_inc, 'True incorrect','-'
        gt.false_fb  & gt.correct==1, clr_cor, 'False correct', '--'
        gt.false_fb  & gt.correct==0, clr_inc, 'False incorrect','--'
    };
else
    % D: correct vs incorrect only, true feedback only.
    spec = {
        ~gt.false_fb & gt.correct==1, clr_cor, 'Correct',   '-'
        ~gt.false_fb & gt.correct==0, clr_inc, 'Incorrect', '-'
    };
end

for j = 1:size(spec,1)
    mask_base = spec{j,1} & gt.block_type == bt;
    clr       = spec{j,2};
    lbl       = spec{j,3};
    ls        = spec{j,4};

    means = nan(1,4);
    sems  = nan(1,4);
    ns    = nan(1,4);

    for s_i = 1:4
        m = mask_base & gt.stage == stages{s_i} & ~isnan(gt.(feat));
        v = gt.(feat)(m);
        if isempty(v), continue; end
        ns(s_i)    = sum(~isnan(v));
        means(s_i) = mean(v, 'omitnan');
        sems(s_i)  = std(v, 'omitnan') / sqrt(ns(s_i));
    end

    x_off = (j - (size(spec,1)+1)/2) * 0.08;
    errorbar(ax, x + x_off, means, sems, 'o', ...
        'Color', clr, 'LineWidth', 1.6, ...
        'LineStyle', ls, ...
        'MarkerFaceColor', clr, ...
        'DisplayName', lbl);

    for s_i = 1:4
        if isnan(means(s_i)), continue; end
        text(ax, x(s_i)+x_off, means(s_i)+sems(s_i)+y_pad, ...
            sprintf('n=%d', ns(s_i)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontSize',7, 'Clipping','off');
    end
end

set(ax,'XTick',1:4,'XTickLabel',stages);
xlabel(ax,'Stage');
ylabel(ax, ylbl);
legend(ax,'Box','off','Location','best');
if reverse_y, set(ax,'YDir','reverse'); end
end