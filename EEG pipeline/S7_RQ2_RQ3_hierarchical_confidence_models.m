% =========================================================================
% S7_RQ2_RQ3_hierarchical_confidence_models.m
%
% Clean standalone RQ2/RQ3 script.
%
% Aim:
%   Test whether confidence predicts EEG components using transparent,
%   hierarchical linear mixed-effects models.
%
% Components:
%   RQ2: prefrontal mean / feedback negativity proxy
%   RQ3: frontal theta power
%
% Model sequence, fitted with ML for valid nested comparisons:
%   M0 confidence only:
%      y_z ~ conf_z + (1|subj_id)
%
%   M1 add outcome:
%      y_z ~ conf_z + outcome + (1|subj_id)
%
%   M2 add block type:
%      y_z ~ conf_z + outcome + block_type + (1|subj_id)
%
%   M3 add block_type x outcome interaction:
%      y_z ~ conf_z + outcome*block_type + (1|subj_id)
%
%   M4 add false-feedback/block-type model, if estimable:
%      y_z ~ conf_z + outcome + block_type*false_fb + (1|subj_id)
%
% Notes:
%   - outcome = Correct vs Incorrect.
%   - false_fb uses the numeric/logical gt.false_fb column directly.
%     No false_fb_cat variable is created.
%   - Model comparison is printed using compare(previous,next).
%   - Simple plots show the main confidence effect only; they are descriptive
%     and correspond to M0. Model-comparison tables are the inferential part.
% ======================================================fal===================

clearvars -except group_table gt
close all
clc

%% ------------------------------------------------------------------------
% PATHS
% -------------------------------------------------------------------------
base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch';

saved_tables_folder = fullfile(base_path, 'Salient mod switch KH', ...
    'Results', 'EEG analysis', 'Outcome_feature_tables_v4_merged');

figure_output_folder = fullfile(base_path, 'Salient mod switch KH', ...
    'Results', 'EEG analysis', 'Figures', 'S7_RQ2_RQ3_hierarchical_confidence_models');

if ~exist(figure_output_folder, 'dir'), mkdir(figure_output_folder); end

% Load table if needed.
if ~exist('group_table','var') && ~exist('gt','var')
    in_file = fullfile(saved_tables_folder, 'group_feature_table_combined.mat');
    fprintf('Loading %s\n', in_file);
    load(in_file, 'group_table');
end
if exist('gt','var') && ~exist('group_table','var')
    group_table = gt;
end
gt = group_table;

%% ------------------------------------------------------------------------
% BASIC CLEANUP / REQUIRED VARIABLES
% -------------------------------------------------------------------------
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'off');
set(groot, 'defaultAxesTickDirMode', 'manual');

% Subject ID.
if ~ismember('subj_id', gt.Properties.VariableNames)
    error('Missing required column subj_id.');
end
gt.subj_id = categorical(string(gt.subj_id));

% Confidence.
if ~ismember('confidence', gt.Properties.VariableNames)
    error('Missing required column confidence.');
end
gt.confidence = double(gt.confidence);

% Outcome.
if ismember('correct_num', gt.Properties.VariableNames)
    gt.correct_num = double(gt.correct_num);
elseif ismember('correct', gt.Properties.VariableNames)
    if iscategorical(gt.correct) || isstring(gt.correct) || iscellstr(gt.correct)
        cstr = lower(string(gt.correct));
        gt.correct_num = nan(height(gt),1);
        gt.correct_num(cstr == "correct" | cstr == "1" | cstr == "true") = 1;
        gt.correct_num(cstr == "incorrect" | cstr == "0" | cstr == "false") = 0;
    else
        gt.correct_num = double(gt.correct);
    end
else
    error('Missing required column correct or correct_num.');
end
% Repair categorical code 1/2 if needed.
u = unique(gt.correct_num(~isnan(gt.correct_num)));
if all(ismember(u, [1 2]))
    gt.correct_num = gt.correct_num - 1;
end
gt.outcome = categorical(gt.correct_num, [0 1], {'Incorrect','Correct'});

% Block type.
if ~ismember('block_type', gt.Properties.VariableNames)
    error('Missing required column block_type.');
end
bt = string(gt.block_type);
bt(bt == "V") = "P";
gt.block_type = categorical(bt, {'D','P'});

% False feedback. Use numeric/logical false_fb directly. Do not create false_fb_cat.
if ismember('trueFB', gt.Properties.VariableNames)
    gt.false_fb = double(~logical(gt.trueFB));
else
    warning('No false_fb or trueFB column found. M4 block_type*false_fb model will be skipped.');
    gt.false_fb = nan(height(gt),1);
end

% Stage is not used in this script, by design. This keeps the first pass clear.

% Z-score confidence within subject.
gt.conf_z = zscore_within_subject(gt.confidence, gt.subj_id);

%% ------------------------------------------------------------------------
% CHOOSE COMPONENT COLUMNS
% -------------------------------------------------------------------------
% RQ2: prefrontal mean / FN proxy. Prefer normalised/z-scored prefrontal mean.
prefrontal_candidates = {'prefrontal_mean_amp', 'prefrontal_mean', ...
    'prefrontal_mean_norm', 'prefrontal_neg_peak_norm', 'prefrontal_neg_peak_amp'};
prefrontal_col = first_existing_col(gt, prefrontal_candidates);
if isempty(prefrontal_col)
    error('Could not find a prefrontal mean/FN column. Tried: %s', strjoin(prefrontal_candidates, ', '));
end

% RQ3: theta.
theta_candidates = {'Theta_amp', 'theta_amp', 'Theta', 'theta', 'frontal_theta', 'prefrontal_theta'};
theta_col = first_existing_col(gt, theta_candidates);
if isempty(theta_col)
    error('Could not find a theta column. Tried: %s', strjoin(theta_candidates, ', '));
end

% Ensure model z columns exist.
[gt, prefrontal_z] = ensure_z_col_in_table(gt, prefrontal_col);
[gt, theta_z]      = ensure_z_col_in_table(gt, theta_col);

fprintf('\n============================================================\n');
fprintf('S7 RQ2/RQ3 hierarchical confidence models\n');
fprintf('============================================================\n');
fprintf('RQ2 prefrontal column: %s | model column: %s\n', prefrontal_col, prefrontal_z);
fprintf('RQ3 theta column     : %s | model column: %s\n', theta_col, theta_z);
fprintf('Model confidence predictor: conf_z = within-subject z-scored confidence\n');
fprintf('Random effect: random intercept for subj_id in all models\n');
fprintf('Output folder: %s\n', figure_output_folder);

%% ------------------------------------------------------------------------
% RUN COMPONENTS
% -------------------------------------------------------------------------
results_RQ2 = run_component_hierarchy(gt, prefrontal_z, prefrontal_col, ...
    'RQ2_PrefrontalMean_FN', 'Prefrontal mean / FN proxy', figure_output_folder, true);

results_RQ3 = run_component_hierarchy(gt, theta_z, theta_col, ...
    'RQ3_Theta', 'Frontal theta power', figure_output_folder, false);

save(fullfile(figure_output_folder, 'S7_RQ2_RQ3_hierarchical_model_results.mat'), ...
    'results_RQ2', 'results_RQ3', 'prefrontal_col', 'prefrontal_z', 'theta_col', 'theta_z');

fprintf('\n============================================================\n');
fprintf('Complete. Results saved to:\n  %s\n', figure_output_folder);
fprintf('============================================================\n');

%% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function out = first_existing_col(T, candidates)
out = '';
for ii = 1:numel(candidates)
    if ismember(candidates{ii}, T.Properties.VariableNames)
        out = candidates{ii};
        return;
    end
end
end

function z = zscore_within_subject(x, subj)
x = double(x);
z = nan(size(x));
subs = categories(categorical(subj));
for si = 1:numel(subs)
    m = categorical(subj) == subs{si};
    vals = x(m);
    sd = std(vals, 'omitnan');
    if sum(m) > 1 && sd > 0
        z(m) = (vals - mean(vals, 'omitnan')) ./ sd;
    end
end
end

function [T, z_col] = ensure_z_col_in_table(T, raw_col)
z_col = [raw_col '_z'];
if ismember(z_col, T.Properties.VariableNames)
    return;
end
T.(z_col) = zscore_within_subject(T.(raw_col), T.subj_id);
end

function results = run_component_hierarchy(gt, y_model_col, y_plot_col, tag, label, outdir, reverse_y)
% Fit and compare hierarchical confidence models for one component.
if nargin < 7, reverse_y = false; end
comp_dir = fullfile(outdir, tag);
if ~exist(comp_dir, 'dir'), mkdir(comp_dir); end

% Analysis table.
keep = ~isnan(gt.(y_model_col)) & ~isnan(gt.(y_plot_col)) & ...
       ~isnan(gt.confidence) & ~isnan(gt.conf_z) & ...
       ~isnan(gt.correct_num) & ~isundefined(gt.block_type) & ...
       ~isundefined(gt.outcome) & ~isundefined(gt.subj_id);
T = gt(keep, :);
T.y_model = double(T.(y_model_col));
T.y_plot  = double(T.(y_plot_col));
T.conf_z  = double(T.conf_z);
T.false_fb = double(T.false_fb);

fprintf('\n\n############################################################\n');
fprintf('%s\n', label);
fprintf('############################################################\n');
fprintf('Rows available: %d | Subjects: %d\n', height(T), numel(unique(T.subj_id)));
fprintf('Outcome counts: Incorrect=%d, Correct=%d\n', sum(T.correct_num==0), sum(T.correct_num==1));
fprintf('Block counts  : D=%d, P=%d\n', sum(T.block_type=='D'), sum(T.block_type=='P'));
if any(~isnan(T.false_fb))
    fprintf('False feedback rows: false_fb=0: %d | false_fb=1: %d\n', sum(T.false_fb==0), sum(T.false_fb==1));
end

% Plot the simplest confidence effect first: corresponds to M0.
plot_main_confidence_effect(T, label, tag, comp_dir, reverse_y);

% Model formulas. FitMethod = ML for nested comparisons.
formulas = struct();
formulas.M0 = 'y_model ~ conf_z + (1|subj_id)';
formulas.M1 = 'y_model ~ conf_z + outcome + (1|subj_id)';
formulas.M2 = 'y_model ~ conf_z + outcome + block_type + (1|subj_id)';
formulas.M3 = 'y_model ~ conf_z + outcome*block_type + (1|subj_id)';
formulas.M4 = 'y_model ~ conf_z + outcome + block_type*false_fb + (1|subj_id)';

fprintf('\nMODEL SEQUENCE FOR %s\n', label);
fprintf('  M0 confidence only       : %s\n', formulas.M0);
fprintf('  M1 + outcome             : %s\n', formulas.M1);
fprintf('  M2 + block_type          : %s\n', formulas.M2);
fprintf('  M3 + outcome*block_type  : %s\n', formulas.M3);
fprintf('  M4 + block_type*false_fb : %s\n', formulas.M4);
fprintf('     Note: M4 is skipped if false_fb is missing or not estimable.\n');

models = struct();
models.M0 = fit_lme_safe(T, formulas.M0, 'M0');
models.M1 = fit_lme_safe(T, formulas.M1, 'M1');
models.M2 = fit_lme_safe(T, formulas.M2, 'M2');
models.M3 = fit_lme_safe(T, formulas.M3, 'M3');

% M4 may be non-estimable if false feedback has no variation or is completely
% confounded with block_type. Try, but do not let it crash the script.
if any(T.false_fb==1) && numel(unique(T.false_fb(~isnan(T.false_fb)))) > 1
    models.M4 = fit_lme_safe(T(~isnan(T.false_fb),:), formulas.M4, 'M4');
else
    models.M4 = [];
    fprintf('\nM4 skipped: false_fb has no usable variation.\n');
end

% Print coefficient summaries.
fprintf('\n--- Fixed-effect coefficients ---\n');
print_coef_block(models.M0, 'M0 confidence only');
print_coef_block(models.M1, 'M1 + outcome');
print_coef_block(models.M2, 'M2 + block_type');
print_coef_block(models.M3, 'M3 + outcome*block_type');
if ~isempty(models.M4), print_coef_block(models.M4, 'M4 + block_type*false_fb'); end

% Model comparison.
fprintf('\n--- Hierarchical model comparisons ---\n');
compare_and_print(models.M0, models.M1, 'M0 vs M1: does adding outcome improve fit?');
compare_and_print(models.M1, models.M2, 'M1 vs M2: does adding block_type improve fit?');
compare_and_print(models.M2, models.M3, 'M2 vs M3: does outcome*block_type improve fit?');
if ~isempty(models.M4)
    compare_and_print(models.M2, models.M4, 'M2 vs M4: does block_type*false_fb improve fit?');
end

% Save tables.
write_model_tables(models, comp_dir, tag);

% Additional simple diagnostic plots for the added terms.
plot_outcome_and_block_summary(T, label, tag, comp_dir, reverse_y);
if ~isempty(models.M4)
    plot_false_feedback_summary(T, label, tag, comp_dir, reverse_y);
end

results = struct();
results.tag = tag;
results.label = label;
results.formulas = formulas;
results.models = models;
results.n_rows = height(T);
results.n_subjects = numel(unique(T.subj_id));
end

function mdl = fit_lme_safe(T, formula, name)
try
    mdl = fitlme(T, formula, 'FitMethod', 'ML');
    fprintf('\n%s fit OK. AIC=%.2f | BIC=%.2f | LogLik=%.2f\n', ...
        name, mdl.ModelCriterion.AIC, mdl.ModelCriterion.BIC, mdl.LogLikelihood);
catch ME
    warning('%s failed: %s', name, ME.message);
    mdl = [];
end
end

function print_coef_block(mdl, title_str)
if isempty(mdl)
    fprintf('\n%s: model unavailable.\n', title_str);
    return;
end
fprintf('\n%s\n', title_str);
C = mdl.Coefficients;
for ii = 1:height(C)
    nm = string(C.Name{ii});
    fprintf('  %-35s beta=% .4f  SE=%.4f  t=% .3f  p=%s  %s\n', ...
        nm, C.Estimate(ii), C.SE(ii), C.tStat(ii), pfmt(C.pValue(ii)), stars(C.pValue(ii)));
end
end

function compare_and_print(m0, m1, title_str)
fprintf('\n%s\n', title_str);
if isempty(m0) || isempty(m1)
    fprintf('  Comparison unavailable because one model failed.\n');
    return;
end
try
    cmp = compare(m0, m1);
    disp(cmp);
catch ME
    warning('Model comparison failed: %s', ME.message);
end
end

function write_model_tables(models, comp_dir, tag)
fnames = fieldnames(models);
for ii = 1:numel(fnames)
    nm = fnames{ii};
    mdl = models.(nm);
    if isempty(mdl), continue; end
    try
        writetable(mdl.Coefficients, fullfile(comp_dir, sprintf('%s_%s_coefficients.csv', tag, nm)));
        try
            writetable(anova(mdl), fullfile(comp_dir, sprintf('%s_%s_anova.csv', tag, nm)));
        catch
        end
    catch ME
        warning('Could not write model table for %s: %s', nm, ME.message);
    end
end
end

function plot_main_confidence_effect(T, label, tag, comp_dir, reverse_y)
fig = figure('Position',[80 80 1050 420]);
sgtitle({sprintf('%s: main confidence effect', label), ...
    'Descriptive plot corresponding to M0: y_z ~ conf_z + (1|subj_id)'}, 'FontSize',12);

ax1 = subplot(1,2,1); hold(ax1,'on');
simple_binned_conf_plot(ax1, T.conf_z, T.y_plot, [0.2 0.2 0.2]);
xlabel(ax1, 'Confidence within subject (z)');
ylabel(ax1, label);
title(ax1, 'All trials: raw plotted scale', 'FontSize',10);
if reverse_y, set(ax1,'YDir','reverse'); end
axis(ax1,'square');

ax2 = subplot(1,2,2); hold(ax2,'on');
simple_binned_conf_plot(ax2, T.conf_z, T.y_model, [0.2 0.2 0.2]);
xlabel(ax2, 'Confidence within subject (z)');
ylabel(ax2, sprintf('%s within-subject z', label));
title(ax2, 'All trials: model scale', 'FontSize',10);
if reverse_y, set(ax2,'YDir','reverse'); end
axis(ax2,'square');

% Add M0 slope summary.
try
    mdl0 = fitlme(T, 'y_model ~ conf_z + (1|subj_id)', 'FitMethod','ML');
    C = mdl0.Coefficients;
    row = strcmp(C.Name, 'conf_z');
    txt = sprintf('M0: y_z ~ conf_z + (1|subj)\nconf_z beta=%.3f, p=%s %s', ...
        C.Estimate(row), pfmt(C.pValue(row)), stars(C.pValue(row)));
    text(ax2, 0.04, 0.96, txt, 'Units','normalized', 'VerticalAlignment','top', ...
        'FontSize',8, 'BackgroundColor','w', 'EdgeColor',[0.8 0.8 0.8]);
catch
end

exportgraphics(fig, fullfile(comp_dir, sprintf('%s_M0_main_confidence_effect.pdf', tag)), 'ContentType','vector');
exportgraphics(fig, fullfile(comp_dir, sprintf('%s_M0_main_confidence_effect.png', tag)), 'Resolution',300);
end

function plot_outcome_and_block_summary(T, label, tag, comp_dir, reverse_y)
fig = figure('Position',[80 80 1200 420]);
sgtitle({sprintf('%s: descriptive checks for added fixed effects', label), ...
    'These plots correspond to the terms added in M1-M3; inference is from the command-line model comparisons'}, 'FontSize',12);

ax1 = subplot(1,3,1); hold(ax1,'on');
plot_group_points(ax1, T.outcome, T.y_model, {'Incorrect','Correct'}, [0.75 0.25 0.20; 0.20 0.55 0.25]);
ylabel(ax1, sprintf('%s z', label));
title(ax1, 'M1 added outcome', 'FontSize',10);
if reverse_y, set(ax1,'YDir','reverse'); end
axis(ax1,'square');

ax2 = subplot(1,3,2); hold(ax2,'on');
plot_group_points(ax2, T.block_type, T.y_model, {'D','P'}, [0.15 0.45 0.70; 0.80 0.30 0.10]);
ylabel(ax2, sprintf('%s z', label));
title(ax2, 'M2 added block type', 'FontSize',10);
if reverse_y, set(ax2,'YDir','reverse'); end
axis(ax2,'square');

ax3 = subplot(1,3,3); hold(ax3,'on');
plot_outcome_by_block(ax3, T, label);
title(ax3, 'M3 outcome × block type', 'FontSize',10);
if reverse_y, set(ax3,'YDir','reverse'); end
axis(ax3,'square');

exportgraphics(fig, fullfile(comp_dir, sprintf('%s_M1_M3_added_terms_summary.pdf', tag)), 'ContentType','vector');
exportgraphics(fig, fullfile(comp_dir, sprintf('%s_M1_M3_added_terms_summary.png', tag)), 'Resolution',300);
end

function plot_false_feedback_summary(T, label, tag, comp_dir, reverse_y)
Tp = T(~isnan(T.false_fb), :);
if height(Tp) < 10 || numel(unique(Tp.false_fb)) < 2, return; end
fig = figure('Position',[80 80 920 420]);
sgtitle({sprintf('%s: false feedback descriptive plot', label), ...
    'Corresponds to M4: y_z ~ conf_z + outcome + block_type*false_fb + (1|subj_id)'}, 'FontSize',12);

ax1 = subplot(1,2,1); hold(ax1,'on');
plot_group_points(ax1, categorical(Tp.false_fb, [0 1], {'True feedback','False feedback'}), ...
    Tp.y_model, {'True feedback','False feedback'}, [0.35 0.35 0.35; 0.75 0.35 0.15]);
ylabel(ax1, sprintf('%s z', label));
title(ax1, 'Main false-feedback contrast', 'FontSize',10);
if reverse_y, set(ax1,'YDir','reverse'); end
axis(ax1,'square');

ax2 = subplot(1,2,2); hold(ax2,'on');
% P block only is the meaningful place for false feedback.
Tp2 = Tp(Tp.block_type=='P', :);
if height(Tp2) > 5
    cats = categorical(Tp2.false_fb + 2*Tp2.correct_num, [0 1 2 3], ...
        {'Inc TrueFB','Inc FalseFB','Cor TrueFB','Cor FalseFB'});
    plot_group_points(ax2, cats, Tp2.y_model, ...
        {'Inc TrueFB','Inc FalseFB','Cor TrueFB','Cor FalseFB'}, ...
        [0.70 0.20 0.18; 0.90 0.45 0.25; 0.20 0.55 0.25; 0.45 0.75 0.35]);
    title(ax2, 'P blocks: false_fb by outcome', 'FontSize',10);
else
    text(ax2,0.5,0.5,'Too few P-block rows','Units','normalized','HorizontalAlignment','center');
end
ylabel(ax2, sprintf('%s z', label));
if reverse_y, set(ax2,'YDir','reverse'); end
axis(ax2,'square');

exportgraphics(fig, fullfile(comp_dir, sprintf('%s_M4_false_feedback_summary.pdf', tag)), 'ContentType','vector');
exportgraphics(fig, fullfile(comp_dir, sprintf('%s_M4_false_feedback_summary.png', tag)), 'Resolution',300);
end

function simple_binned_conf_plot(ax, x, y, clr)
ok = ~isnan(x) & ~isnan(y);
x = x(ok); y = y(ok);
if isempty(x), return; end
scatter(ax, x, y, 8, [0.55 0.55 0.55], 'filled', 'MarkerFaceAlpha',0.12, 'HandleVisibility','off');
% Binned means.
edges = quantile_omitnan(x, [0 .25 .5 .75 1]);
if numel(unique(edges)) < 3
    edges = linspace(min(x), max(x), 5);
end
xb = nan(4,1); yb = nan(4,1); ci = nan(4,1);
for bi = 1:4
    if bi < 4
        m = x >= edges(bi) & x < edges(bi+1);
    else
        m = x >= edges(bi) & x <= edges(bi+1);
    end
    vals = y(m);
    xb(bi) = mean(x(m), 'omitnan');
    yb(bi) = mean(vals, 'omitnan');
    ci(bi) = 1.96 * std(vals, 'omitnan') / sqrt(max(sum(~isnan(vals)),1));
end
errorbar(ax, xb, yb, ci, 'o', 'Color', clr, 'MarkerFaceColor', clr, ...
    'MarkerSize',7, 'LineWidth',1.5, 'CapSize',6, 'HandleVisibility','off');
% Simple fitted line for visual guidance only.
if numel(x) > 5
    p = polyfit(x, y, 1);
    xx = linspace(min(x), max(x), 100);
    plot(ax, xx, polyval(p,xx), '-', 'Color', clr, 'LineWidth',1.2, 'HandleVisibility','off');
end
end

function plot_group_points(ax, group_var, y, labels, clrs)
g = string(group_var);
y = double(y);
for ii = 1:numel(labels)
    m = g == string(labels{ii});
    vals = y(m & ~isnan(y));
    if isempty(vals), continue; end
    x = ii + 0.18*(rand(numel(vals),1)-0.5);
    scatter(ax, x, vals, 14, [0.45 0.45 0.45], 'filled', 'MarkerFaceAlpha',0.18, 'HandleVisibility','off');
    ci = 1.96 * std(vals,'omitnan') / sqrt(max(numel(vals),1));
    errorbar(ax, ii, mean(vals,'omitnan'), ci, ...
        'o', 'Color', clrs(ii,:), 'MarkerFaceColor', clrs(ii,:), 'MarkerSize',9, ...
        'LineWidth',1.8, 'CapSize',8, 'HandleVisibility','off');
end
set(ax, 'XTick', 1:numel(labels), 'XTickLabel', labels, 'TickDir','out');
xlim(ax, [0.5 numel(labels)+0.5]);
yline(ax,0,'k:','HandleVisibility','off');
end

function plot_outcome_by_block(ax, T, label)
outcomes = {'Incorrect','Correct'};
blocks = {'D','P'};
clrs = [0.75 0.25 0.20; 0.20 0.55 0.25];
for oi = 1:2
    for bi = 1:2
        x = bi + (oi-1.5)*0.18;
        m = string(T.outcome)==outcomes{oi} & string(T.block_type)==blocks{bi};
        vals = T.y_model(m & ~isnan(T.y_model));
        if isempty(vals), continue; end
        scatter(ax, x + 0.08*(rand(numel(vals),1)-0.5), vals, 12, [0.45 0.45 0.45], ...
            'filled', 'MarkerFaceAlpha',0.12, 'HandleVisibility','off');
        errorbar(ax, x, mean(vals,'omitnan'), 1.96*std(vals,'omitnan')/sqrt(numel(vals)), ...
            'o', 'Color', clrs(oi,:), 'MarkerFaceColor', clrs(oi,:), ...
            'MarkerSize',8, 'LineWidth',1.6, 'CapSize',6, 'HandleVisibility','off');
    end
end
plot(ax, NaN, NaN, 'o', 'Color', clrs(1,:), 'MarkerFaceColor', clrs(1,:), 'DisplayName','Incorrect');
plot(ax, NaN, NaN, 'o', 'Color', clrs(2,:), 'MarkerFaceColor', clrs(2,:), 'DisplayName','Correct');
set(ax, 'XTick',1:2, 'XTickLabel',blocks, 'TickDir','out');
xlim(ax,[0.5 2.5]);
ylabel(ax, sprintf('%s z', label));
yline(ax,0,'k:','HandleVisibility','off');
legend(ax,'Box','off','Location','best');
end

function q = quantile_omitnan(x, p)
x = sort(x(~isnan(x)));
q = nan(size(p));
if isempty(x), return; end
for ii = 1:numel(p)
    idx = 1 + (numel(x)-1)*p(ii);
    lo = floor(idx); hi = ceil(idx);
    if lo == hi
        q(ii) = x(lo);
    else
        q(ii) = x(lo) + (idx-lo)*(x(hi)-x(lo));
    end
end
end

function s = stars(p)
if isnan(p), s = '';
elseif p < 0.001, s = '***';
elseif p < 0.01, s = '**';
elseif p < 0.05, s = '*';
elseif p < 0.10, s = '†';
else, s = 'ns';
end
end

function s = pfmt(p)
if isnan(p), s = 'NaN';
elseif p < 0.001, s = '<.001';
else, s = sprintf('%.3f', p);
end
end
