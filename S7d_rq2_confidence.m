% =============================================================================
% S7_RQ2_confidence_uncertainty.m
%
% SIMPLIFIED REPLACEMENT for plot_confidence_fRN_extended and its helpers.
%
% PURPOSE
% -------
% Three focused figures + LME tables for RQ2 and the uncertainty history
% sub-analyses. Everything is generated from one call:
%
%   plot_rq2_confidence_uncertainty(gt, outdir, fn_feat, fn_ylbl, reverse_y)
%
% FIGURES PRODUCED
% ----------------
%   RQ2_Fig1_confidence_FN_by_blocktype_transition.pdf
%     Panel 1 — FN vs confidence (z), incorrect trials, split D / P block.
%              Points binned into 5 quantile bins; group mean ± SEM per bin.
%              LME annotation (conf_z slope significance per block type).
%     Panel 2 — Same but split by RECENT TRANSITION (D→D, D→P, P→D, P→P).
%              Four colours, two columns (incorrect / correct), no lines
%              connecting categorical X-axis points.
%
%   RQ2_Fig2_nPrevP_FN_by_stage.pdf
%     How cumulative prior P-block exposure (n_prev_P = 0, 1, 2+) modulates
%     FN within each stage × current block type.
%     Layout: 2 columns (D, P) × 4 rows (LN, LE, RN, RE).
%     Each cell: dot + 95% CI per n_prev_P level; subject jitter overlaid.
%     LME: FN ~ n_prev_P_bin * stage + (1|subj_id) within each block type.
%
%   RQ2_Fig3_nassar_uncertainty_FN.pdf
%     How Nassar model latents (certainty |θ−0.5|, omega ω, surprise ω×|δ|,
%     alpha α) relate to FN per trial.
%     Layout: 4 sub-panels (one per latent), each showing:
%       - Binned scatter: FN on Y, latent (z-scored within subject) on X.
%       - Two colours: correct (green) vs incorrect (red) trials.
%       - Horizontal null line; LME slope annotation.
%
% LME TABLES
% ----------
% Each figure saves a matching *_LME.csv with all fixed-effect coefficients
% (b, SE, t, df, p) and model fit indices (AIC, BIC, R2m, R2c).
%
% INPUTS
% ------
%   gt         — combined per-trial feature table (group_table_combined from S4).
%                Required columns: subj_id, block_type, stage, correct,
%                false_fb, confidence, prefrontal_mean_norm (or alt name
%                passed via fn_feat), transition_recent, n_prev_P (or
%                n_prev_P_bin). Nassar columns optional: certainty,
%                omega_trial, surprise, alpha_nassar (or alpha_trial).
%   outdir     — save folder (created if absent).
%   fn_feat    — FN feature column to plot (default: 'prefrontal_mean_norm').
%   fn_ylbl    — y-axis label string.
%   reverse_y  — logical; true flips y-axis so more negative = up (default true).
%
% DESIGN NOTES
% ------------
% * All LME models use within-subject z-scored predictors so betas are
%   comparable across measures. The plot variable (fn_feat) may be norm or
%   raw; the MODEL always targets prefrontal_mean_norm_z for consistency
%   with the parent script.
% * Subject-level means are computed before group means to avoid
%   trial-count imbalance inflating precision.
% * No categorical lines are drawn between X-axis positions when X is
%   nominal (transitions, n_prev_P bins) — only dots + error bars.
% * Nassar latent columns are searched by multiple plausible names and
%   skipped gracefully if absent, with a text notice in the figure.
% =============================================================================

function plot_rq2_confidence_uncertainty(gt, outdir, fn_feat, fn_ylbl, reverse_y)
%PLOT_RQ2_CONFIDENCE_UNCERTAINTY  Top-level entry point.

% ── Defaults ────────────────────────────────────────────────────────────
if nargin < 3 || isempty(fn_feat),   fn_feat   = rq2_best_fn_col(gt);        end
if nargin < 4 || isempty(fn_ylbl),   fn_ylbl   = 'Prefrontal neg. peak / baseline RMS'; end
if nargin < 5 || isempty(reverse_y), reverse_y = true;                        end
if ~exist(outdir,'dir'), mkdir(outdir); end

% Model always uses z-scored version for coefficient comparability
model_y = rq2_model_col(gt, fn_feat);

% ── Prepare shared columns ───────────────────────────────────────────────
gt = rq2_add_working_cols(gt);

% ── Run three figures ────────────────────────────────────────────────────
rq2_fig1_confidence_blocktype_transition(gt, outdir, fn_feat, fn_ylbl, model_y, reverse_y);
rq2_fig2_nPrevP_by_stage(gt, outdir, fn_feat, fn_ylbl, model_y, reverse_y);
rq2_fig3_nassar_latents(gt, outdir, fn_feat, fn_ylbl, model_y, reverse_y);

fprintf('\nRQ2 confidence/uncertainty figures saved to:\n  %s\n', outdir);
end


% =============================================================================
%% FIGURE 1 — Confidence × FN, split by block type and recent transition
% =============================================================================
function rq2_fig1_confidence_blocktype_transition(gt, outdir, fn_feat, fn_ylbl, model_y, reverse_y)

CLR_D = [0.15 0.45 0.70];
CLR_P = [0.80 0.30 0.10];
TRANS_COLS = [0.12 0.47 0.71;   % D→D
              0.85 0.33 0.10;   % D→P
              0.47 0.67 0.19;   % P→D
              0.80 0.20 0.60];  % P→P
TRANS_LABS = {'D→D','D→P','P→D','P→P'};
OUTCOME_LABS = {'Incorrect','Correct'};
OUTCOME_CODES = [0, 1];

fig = figure('Position',[40 40 1400 700], 'Color','w');
sgtitle({'RQ2  —  Confidence × Prefrontal Negativity', ...
    'Left: by block type (incorrect trials).   Right: by recent block transition, split by outcome.'}, ...
    'FontSize',12,'FontWeight','bold');

% Layout: 2 rows (incorrect, correct) × 3 cols.
%   Col 1 = D-block confidence slope.
%   Col 2 = P-block confidence slope.
%   Col 3 = transition × FN (four coloured dots, one per transition).
%
% Grid: subplot(2, 3, panel_index)
%   Row 1 (incorrect): panels 1, 2, 3
%   Row 2 (correct):   panels 4, 5, 6

bts  = {'D','P'};
clrs = {CLR_D, CLR_P};
trans_list = {'D→D','D→P','P→D','P→P'};
n_tr = numel(trans_list);

% Precompute masks for both outcomes so they can be reused
outcome_masks = {
    gt.correct_num == 0 & ~gt.false_fb & ~isnan(gt.conf_z) & ~isnan(gt.(fn_feat));  % incorrect
    gt.correct_num == 1 & ~gt.false_fb & ~isnan(gt.conf_z) & ~isnan(gt.(fn_feat));  % correct
};
inc_mask = outcome_masks{1};   % used for LME table at end

for oi = 1:2   % row: outcome (1=incorrect, 2=correct)
    base_oi = (oi-1)*3;   % panel offset for this row

    % ── Cols 1-2: D and P confidence slopes ──────────────────────────────
    for bi = 1:2
        ax = subplot(2, 3, base_oi + bi);
        hold(ax,'on');
        title(ax, sprintf('%s trials — %s block', OUTCOME_LABS{oi}, bts{bi}), 'FontSize',10);

        m = outcome_masks{oi} & string(gt.block_type) == bts{bi};
        [sig, n_obs] = rq2_conf_lme(gt(m,:), model_y, bts{bi});

        rq2_binned_dot_ci(ax, gt.conf_z(m), gt.(fn_feat)(m), clrs{bi}, 5);
        add_sig_text(ax, sig, n_obs);
        rq2_axis_labels(ax, 'Confidence (z)', fn_ylbl, reverse_y);
        xline(ax,0,'k:','LineWidth',0.8,'HandleVisibility','off');
        yline(ax,0,'k--','LineWidth',0.8,'HandleVisibility','off');
        axis(ax,'square');
    end

    % ── Col 3: transition × FN (dots, no connecting lines) ───────────────
    ax = subplot(2, 3, base_oi + 3);
    hold(ax,'on');
    title(ax, sprintf('Recent transition — %s trials', OUTCOME_LABS{oi}), 'FontSize',10);

    m_base = outcome_masks{oi} & ~isnan(gt.conf_z) & ~isnan(gt.(fn_feat));

    for ti = 1:n_tr
        m = m_base & string(gt.transition_recent) == trans_list{ti};
        vals = rq2_subj_means_conf_fn(gt(m,:), fn_feat);
        if isempty(vals), continue; end
        n_s  = numel(vals);
        mn   = mean(vals,'omitnan');
        se   = std(vals,'omitnan') / sqrt(n_s);
        ci   = 1.96 * se;
        jx   = ti + (rand(n_s,1)-0.5)*0.15;
        scatter(ax, jx, vals, 18, TRANS_COLS(ti,:), 'filled', ...
            'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0,'HandleVisibility','off');
        errorbar(ax, ti, mn, ci, 'o','LineStyle','none', ...
            'Color',TRANS_COLS(ti,:),'MarkerFaceColor',TRANS_COLS(ti,:), ...
            'MarkerEdgeColor','k','MarkerSize',8,'LineWidth',1.5,'CapSize',6, ...
            'DisplayName',sprintf('%s (n=%d)', trans_list{ti}, n_s));
    end

    yline(ax,0,'k--','LineWidth',0.8,'HandleVisibility','off');
    set(ax,'XTick',1:n_tr,'XTickLabel',trans_list,'XTickLabelRotation',25);
    xlabel(ax,'Recent block transition');
    ylabel(ax, fn_ylbl,'Interpreter','none');
    if reverse_y, set(ax,'YDir','reverse'); end
    legend(ax,'Box','off','Location','best','FontSize',8);
    axis(ax,'square');
end   % oi loop

% ── Save ─────────────────────────────────────────────────────────────────
fname = 'RQ2_Fig1_confidence_FN_blocktype_transition';
save_both(fig, outdir, fname);

% ── LME table for block-type confidence slopes ───────────────────────────
% Model: FN_z ~ conf_z * block_type + stage + (1 + conf_z | subj_id)
% Run on incorrect, true-feedback trials only
try
    Tm = rq2_prep_lme_table(gt, model_y, inc_mask);
    if height(Tm) >= 20
        mdl = fitlme(Tm, [model_y ' ~ conf_z * block_type + stage + (1 | subj_id)'], ...
            'FitMethod','REML');
        save_lme_csv(mdl, fullfile(outdir, [fname '_LME.csv']), ...
            'FN ~ conf_z * block_type + stage');
    end
catch ME
    fprintf('  Fig1 LME failed: %s\n', ME.message);
end
end


% =============================================================================
%% FIGURE 2 — n_prev_P × FN, split by stage and current block type
% =============================================================================
function rq2_fig2_nPrevP_by_stage(gt, outdir, fn_feat, fn_ylbl, model_y, reverse_y)

CLR_D = [0.15 0.45 0.70];
CLR_P = [0.80 0.30 0.10];
BT_CLRS  = {CLR_D, CLR_P};
BT_LABS  = {'D','P'};
STAGES   = {'LN','LE','RN','RE'};
NPPBINS  = {'0','1','2+'};
NPPX     = [1 2 3];           % x positions for the three bins
BIN_COLS = [0.30 0.30 0.30;   % 0 prior P blocks
            0.60 0.40 0.10;   % 1
            0.80 0.15 0.15];  % 2+

fig = figure('Position',[50 50 1100 900], 'Color','w');
sgtitle({'RQ2  —  Cumulative prior P-block exposure × Prefrontal Negativity', ...
    'Each cell: FN by n\_prev\_P (0 / 1 / 2+).  Rows = stage.  Cols = current block type.'}, ...
    'FontSize',12,'FontWeight','bold');

n_stage = numel(STAGES);
n_bt    = numel(BT_LABS);
n_bins  = numel(NPPBINS);

all_mdl_rows = {};   % collect LME output

for bi = 1:n_bt
    for si = 1:n_stage
        ax = subplot(n_stage, n_bt, (si-1)*n_bt + bi);
        hold(ax,'on');
        title(ax, sprintf('%s  |  %s block', STAGES{si}, BT_LABS{bi}), 'FontSize',9);

        m_cell = ~gt.false_fb & ~isnan(gt.(fn_feat)) & ...
                 string(gt.block_type) == BT_LABS{bi} & ...
                 string(gt.stage)      == STAGES{si};

        all_subj_vals = cell(n_bins,1);

        for ni = 1:n_bins
            m = m_cell & string(gt.n_prev_P_bin) == NPPBINS{ni};
            subs = unique(string(gt.subj_id(m)));
            sv   = nan(numel(subs),1);
            for ki = 1:numel(subs)
                sm = m & string(gt.subj_id) == subs(ki);
                sv(ki) = mean(gt.(fn_feat)(sm),'omitnan');
            end
            sv = sv(~isnan(sv));
            all_subj_vals{ni} = sv;

            if isempty(sv), continue; end
            n_s = numel(sv);
            mn  = mean(sv,'omitnan');
            se  = std(sv,'omitnan') / sqrt(n_s);
            ci  = 1.96 * se;

            % Subject jitter
            jx = NPPX(ni) + (rand(n_s,1)-0.5)*0.12;
            scatter(ax, jx, sv, 15, BIN_COLS(ni,:), 'filled', ...
                'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0,'HandleVisibility','off');

            % Group CI dot
            errorbar(ax, NPPX(ni), mn, ci, 'o','LineStyle','none', ...
                'Color',BIN_COLS(ni,:),'MarkerFaceColor',BIN_COLS(ni,:), ...
                'MarkerEdgeColor','k','MarkerSize',9,'LineWidth',1.6,'CapSize',7, ...
                'DisplayName',sprintf('n_prevP=%s  (n=%d)', NPPBINS{ni}, n_s));
        end

        yline(ax,0,'k--','LineWidth',0.8,'HandleVisibility','off');

        % Compact LME for this cell: FN ~ n_prev_P (linear) + (1|subj)
        try
            Tc = gt(m_cell & ~isnan(gt.(fn_feat)), :);
            Tc.n_prev_P_bin = categorical(string(Tc.n_prev_P_bin), {'0','1','2+'}, 'Ordinal', false);
            Tc.subj_id      = categorical(string(Tc.subj_id));
            if height(Tc) >= 15 && numel(unique(string(Tc.subj_id))) >= 4 && ...
               numel(categories(removecats(Tc.n_prev_P_bin))) >= 2
                fm  = fitlme(Tc, [fn_feat ' ~ n_prev_P_bin + (1|subj_id)'], 'FitMethod','REML');
                C   = fm.Coefficients;
                % Find the n_prev_P_bin linear/ordinal term
                p_main = min(C.pValue(2:end));   % smallest p among bin levels
                sig_str = p_to_stars(p_main);
                text(ax,0.03,0.97, sprintf('LME: n_prevP %s', sig_str), ...
                    'Units','normalized','VerticalAlignment','top','FontSize',7, ...
                    'BackgroundColor','w','Margin',1,'Interpreter','none');
                all_mdl_rows(end+1,:) = {BT_LABS{bi}, STAGES{si}, ...
                    C.Estimate(1), C.Estimate(end), C.SE(end), ...
                    C.tStat(end), C.DF(end), C.pValue(end)};
            end
        catch
        end

        set(ax,'XTick',NPPX,'XTickLabel',{'0','1','2+'});
        if si == n_stage, xlabel(ax,'n prior P blocks'); end
        if bi == 1, ylabel(ax, fn_ylbl,'Interpreter','none'); end
        if reverse_y, set(ax,'YDir','reverse'); end
        axis(ax,'square');
    end
end

% Legend (once, for the block-type column headers)
annotation(fig,'textbox',[0.01 0.01 0.98 0.04], ...
    'String',['Grey dots = individual subjects (averaged within subject×cell).  ' ...
              'Error bars = mean ± 95% CI.  Stars from LME: FN ~ n\_prev\_P\_bin + (1|subj\_id) per cell.'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

fname = 'RQ2_Fig2_nPrevP_FN_by_stage';
save_both(fig, outdir, fname);

% ── LME table: omnibus model across all cells ─────────────────────────────
% FN ~ n_prev_P_bin * stage * block_type + (1|subj_id)
try
    m_all = ~gt.false_fb & ~isnan(gt.(fn_feat)) & ~isnan(gt.conf_z);
    Tall  = rq2_prep_lme_table(gt, model_y, m_all);
    Tall.n_prev_P_bin = categorical(string(Tall.n_prev_P_bin), {'0','1','2+'}, 'Ordinal', false);
    if height(Tall) >= 30
        mdl_omni = fitlme(Tall, ...
            [model_y ' ~ n_prev_P_bin * stage * block_type + (1 | subj_id)'], ...
            'FitMethod','REML');
        save_lme_csv(mdl_omni, fullfile(outdir, [fname '_omnibus_LME.csv']), ...
            'FN ~ n_prev_P_bin * stage * block_type');
    end
catch ME
    fprintf('  Fig2 omnibus LME failed: %s\n', ME.message);
end

% Per-cell summary
if ~isempty(all_mdl_rows)
    try
        T_cell = cell2table(all_mdl_rows, 'VariableNames', ...
            {'block_type','stage','intercept_b','nPrevP_b','nPrevP_SE','nPrevP_t','nPrevP_df','nPrevP_p'});
        writetable(T_cell, fullfile(outdir, [fname '_cell_LME.csv']));
    catch
    end
end
end


% =============================================================================
%% FIGURE 3 — Nassar model latents × FN
% =============================================================================
function rq2_fig3_nassar_latents(gt, outdir, fn_feat, fn_ylbl, model_y, reverse_y)

% ── Discover available Nassar columns ────────────────────────────────────
latent_specs = {
    'certainty',          'certainty_trial',    'certainty',       '|\theta - 0.5|  (prospective certainty)';
    'omega',              'omega_trial',         'omega',           '\omega_t  (change-point probability)';
    'surprise',           'surprise',            'surprise',        '\omega \times |\delta|  (retrospective surprise)';
    'alpha_nassar',       'alpha_trial',         'alpha',           '\alpha_t  (effective learning rate)';
};

% resolved_specs: {col_name, label, unit_label} for columns that exist
resolved = {};
for ri = 1:size(latent_specs,1)
    candidates = {latent_specs{ri,1}, latent_specs{ri,2}, latent_specs{ri,3}};
    for ci = 1:numel(candidates)
        if ismember(candidates{ci}, gt.Properties.VariableNames)
            resolved(end+1,:) = {candidates{ci}, latent_specs{ri,4}}; %#ok<AGROW>
            break;
        end
    end
end

n_lat = size(resolved,1);

if n_lat == 0
    fig = figure('Position',[50 50 700 300],'Color','w');
    text(0.5,0.5, {'No Nassar latent columns found in gt.', ...
        'Expected: certainty / certainty_trial, omega / omega_trial,', ...
        'surprise, alpha_nassar / alpha_trial.'}, ...
        'Units','normalized','HorizontalAlignment','center','FontSize',10);
    axis off;
    save_both(fig, outdir, 'RQ2_Fig3_nassar_latents_MISSING');
    return;
end

CLR_INC = [0.75 0.10 0.10];   % incorrect
CLR_COR = [0.10 0.60 0.10];   % correct
OUTCOME_LABS  = {'Incorrect','Correct'};
OUTCOME_CODES = [0, 1];
OUTCOME_CLRS  = {CLR_INC, CLR_COR};

% Figure layout: n_lat columns × 2 rows (incorrect / correct)
fig = figure('Position',[40 40 350*n_lat 720], 'Color','w');
sgtitle({'RQ2  —  Nassar model latents × Prefrontal Negativity', ...
    'X = latent variable (z-scored within subject).  Rows = outcome.  Stars from LME slope.'}, ...
    'FontSize',12,'FontWeight','bold');

lme_rows = {};

for li = 1:n_lat
    col = resolved{li,1};
    lbl = resolved{li,2};

    % Z-score within subject
    gt_z_col = ['z_' matlab.lang.makeValidName(col)];
    gt.(gt_z_col) = nan(height(gt),1);
    subs = unique(string(gt.subj_id));
    for si = 1:numel(subs)
        sm = string(gt.subj_id) == subs(si) & ~isnan(gt.(col));
        v  = gt.(col)(sm);
        if sum(~isnan(v)) > 1 && std(v,'omitnan') > 0
            gt.(gt_z_col)(sm) = (v - mean(v,'omitnan')) ./ std(v,'omitnan');
        end
    end

    % Fit one LME for the slope (both outcomes, block types pooled)
    % FN_z ~ latent_z * correct_cat + block_type + stage + (1|subj_id)
    m_lme   = ~gt.false_fb & ~isnan(gt.(gt_z_col)) & ~isnan(gt.(model_y));
    % Apply the same double-filter that rq2_prep_lme_table applies internally
    m_full  = m_lme & ~isnan(gt.(model_y));
    lme_sig = 'n/a';
    try
        Tlm = rq2_prep_lme_table(gt, model_y, m_full);
        % Add the latent z-column — use the SAME double-mask so lengths agree
        Tlm.(gt_z_col) = gt.(gt_z_col)(m_full);
        % Build formula string using the actual column name (gt_z_col is a string)
        formula = sprintf('%s ~ %s * correct_cat + block_type + stage + (1|subj_id)', ...
                          model_y, gt_z_col);
        if height(Tlm) >= 20 && numel(unique(string(Tlm.subj_id))) >= 4
            fm = fitlme(Tlm, formula, 'FitMethod','REML');
            C   = fm.Coefficients;
            p_slope = rq2_coef_p(C, gt_z_col, false);
            p_inter = rq2_coef_p(C, gt_z_col, true);
            lme_sig = sprintf('slope %s  ×outcome %s', p_to_stars(p_slope), p_to_stars(p_inter));
            save_lme_csv(fm, fullfile(outdir, sprintf('RQ2_Fig3_%s_LME.csv', matlab.lang.makeValidName(col))), ...
                formula);
            lme_rows(end+1,:) = {col, p_slope, p_inter};
        end
    catch ME
        lme_sig = sprintf('LME failed: %s', ME.message(1:min(60,numel(ME.message))));
    end

    % Plot: one row per outcome
    for oi = 1:2
        ax = subplot(2, n_lat, (oi-1)*n_lat + li);
        hold(ax,'on');
        title(ax, sprintf('%s\n%s', lbl, OUTCOME_LABS{oi}), ...
            'FontSize',8,'Interpreter','tex');

        m = ~gt.false_fb & gt.correct_num == OUTCOME_CODES(oi) & ...
            ~isnan(gt.(gt_z_col)) & ~isnan(gt.(fn_feat));
        clr = OUTCOME_CLRS{oi};

        rq2_binned_dot_ci(ax, gt.(gt_z_col)(m), gt.(fn_feat)(m), clr, 5);

        if oi == 1
            % LME annotation only on top row to avoid duplication
            add_sig_text(ax, lme_sig, sum(m));
        end

        xline(ax,0,'k:','LineWidth',0.8,'HandleVisibility','off');
        yline(ax,0,'k--','LineWidth',0.8,'HandleVisibility','off');
        xlabel(ax,[strrep(col,'_',' ') ' (z)'],'Interpreter','none');
        if li == 1, ylabel(ax, fn_ylbl,'Interpreter','none'); end
        if reverse_y, set(ax,'YDir','reverse'); end
        axis(ax,'square');
    end
end

annotation(fig,'textbox',[0.01 0.01 0.98 0.04], ...
    'String',['Binned dots = quartile means ± 95% CI of per-trial FN.  ' ...
              'Top row = incorrect; bottom = correct.  ' ...
              'LME stars: slope = main effect of latent; ×outcome = interaction.  ' ...
              'Latents z-scored within subject before fitting.'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

fname = 'RQ2_Fig3_nassar_latents_FN';
save_both(fig, outdir, fname);
end


% =============================================================================
%% SHARED HELPERS — all private to this file
% =============================================================================

function gt = rq2_add_working_cols(gt)
%RQ2_ADD_WORKING_COLS  Ensure derived columns exist in gt.

% correct_num
if ~ismember('correct_num', gt.Properties.VariableNames)
    if ismember('correct', gt.Properties.VariableNames)
        c = gt.correct;
        if isnumeric(c) || islogical(c)
            gt.correct_num = double(c);
        else
            cs = lower(string(c));
            gt.correct_num = double(cs=="1" | cs=="correct" | cs=="true");
        end
    else
        gt.correct_num = nan(height(gt),1);
    end
end

% false_fb as logical
if ~ismember('false_fb', gt.Properties.VariableNames)
    gt.false_fb = false(height(gt),1);
else
    gt.false_fb = logical(gt.false_fb);
end

% conf_z (within-subject z of confidence)
if ~ismember('conf_z', gt.Properties.VariableNames)
    gt.conf_z = nan(height(gt),1);
    if ismember('confidence', gt.Properties.VariableNames)
        subs = unique(string(gt.subj_id));
        for si = 1:numel(subs)
            sm = string(gt.subj_id) == subs(si);
            c  = double(gt.confidence(sm));
            if sum(~isnan(c)) > 1 && std(c,'omitnan') > 0
                gt.conf_z(sm) = (c - mean(c,'omitnan')) ./ std(c,'omitnan');
            end
        end
    end
end

% block_type as string for comparison
if ~ismember('block_type_s', gt.Properties.VariableNames)
    bt = string(gt.block_type);
    bt(bt=="V") = "P";
    gt.block_type_s = bt;
end

% transition_recent — add if missing
if ~ismember('transition_recent', gt.Properties.VariableNames)
    % Inline the transition/n_prev_P computation rather than calling the
    % parent script's add_transition_history_columns (which is a local
    % function only visible within that file).
    gt = rq2_add_transition_cols(gt);
end

% n_prev_P_bin — binned version of n_prev_P
if ~ismember('n_prev_P_bin', gt.Properties.VariableNames)
    if ismember('n_prev_P', gt.Properties.VariableNames)
        npp = double(gt.n_prev_P);
        bins = repmat("0", height(gt),1);
        bins(npp == 1) = "1";
        bins(npp >= 2) = "2+";
        bins(isnan(npp)) = "0";
        gt.n_prev_P_bin = bins;
    else
        gt.n_prev_P_bin = repmat("0", height(gt), 1);
        warning('rq2_add_working_cols: n_prev_P not found; n_prev_P_bin set to 0 for all rows.');
    end
end

% correct_cat for LME interaction
if ~ismember('correct_cat', gt.Properties.VariableNames)
    gt.correct_cat = categorical(gt.correct_num, [0 1], {'Incorrect','Correct'});
end

% block_type as categorical
if ~iscategorical(gt.block_type)
    gt.block_type = categorical(string(gt.block_type), {'D','P'});
end

% stage as categorical
if ~iscategorical(gt.stage)
    gt.stage = categorical(string(gt.stage), {'LN','LE','RN','RE'}, 'Ordinal', false);
end
end

% ─────────────────────────────────────────────────────────────────────────────
function col = rq2_best_fn_col(gt)
%RQ2_BEST_FN_COL  Pick the most informative FN column available.
candidates = {'prefrontal_mean_norm','prefrontal_mean_amp', ...
              'FRN_mean_norm','FRN_mean_amp','prefrontal_mean_norm_z'};
for c = candidates
    if ismember(c{1}, gt.Properties.VariableNames)
        col = c{1}; return;
    end
end
error('rq2_best_fn_col: no FN amplitude column found in gt.');
end

% ─────────────────────────────────────────────────────────────────────────────
function col = rq2_model_col(gt, fn_feat)
%RQ2_MODEL_COL  Prefer the z-scored version of fn_feat for LME.
z_col = [fn_feat '_z'];
if ismember(z_col, gt.Properties.VariableNames)
    col = z_col;
else
    col = fn_feat;   % fall back to the raw/norm column itself
end
end

% ─────────────────────────────────────────────────────────────────────────────
function Tm = rq2_prep_lme_table(gt, model_y, mask)
%RQ2_PREP_LME_TABLE  Subset and type-coerce for fitlme.
Tm = gt(mask & ~isnan(gt.(model_y)), :);
Tm.subj_id    = categorical(string(Tm.subj_id));
Tm.block_type = categorical(string(Tm.block_type), {'D','P'});
Tm.stage      = categorical(string(Tm.stage), {'LN','LE','RN','RE'}, 'Ordinal', false);
if ~ismember('correct_cat', Tm.Properties.VariableNames)
    Tm.correct_cat = categorical(Tm.correct_num, [0 1], {'Incorrect','Correct'});
end
end

% ─────────────────────────────────────────────────────────────────────────────
function rq2_binned_dot_ci(ax, x, y, clr, n_bins)
%RQ2_BINNED_DOT_CI  Quantile-binned group means ± 95% CI.
% Trial-level scatter in background; no regression line drawn.
x = double(x); y = double(y);
ok = ~isnan(x) & ~isnan(y);
x  = x(ok); y = y(ok);
if numel(x) < 6, return; end

% Light scatter
scatter(ax, x, y, 6, [0.75 0.75 0.75], 'filled', ...
    'MarkerFaceAlpha',0.12,'MarkerEdgeAlpha',0,'HandleVisibility','off');

% Quantile bin edges
edges = unique(quantile(x, linspace(0,1,n_bins+1)));
if numel(edges) < 2, return; end

for bi = 1:numel(edges)-1
    if bi == numel(edges)-1
        m = x >= edges(bi) & x <= edges(bi+1);
    else
        m = x >= edges(bi) & x < edges(bi+1);
    end
    xv = x(m); yv = y(m);
    xv = xv(~isnan(yv)); yv = yv(~isnan(yv));
    if numel(yv) < 2, continue; end
    xm  = mean(xv,'omitnan');
    ym  = mean(yv,'omitnan');
    ci  = 1.96 * std(yv,'omitnan') / sqrt(numel(yv));
    errorbar(ax, xm, ym, ci, 'o','LineStyle','none', ...
        'Color',clr,'MarkerFaceColor',clr,'MarkerEdgeColor','k', ...
        'MarkerSize',7,'LineWidth',1.3,'CapSize',5,'HandleVisibility','off');
end
end

% ─────────────────────────────────────────────────────────────────────────────
function vals = rq2_subj_means_conf_fn(T, fn_feat)
%RQ2_SUBJ_MEANS_CONF_FN  One per-subject mean for a condition subset.
subs = unique(string(T.subj_id));
vals = nan(numel(subs),1);
for si = 1:numel(subs)
    sm = string(T.subj_id) == subs(si) & ~isnan(T.(fn_feat));
    if any(sm), vals(si) = mean(T.(fn_feat)(sm),'omitnan'); end
end
vals = vals(~isnan(vals));
end

% ─────────────────────────────────────────────────────────────────────────────
function [sig_str, n_obs] = rq2_conf_lme(Tm, model_y, bt_label)
%RQ2_CONF_LME  Fit conf_z slope in one block type; return annotation string.
sig_str = 'LME n/a';
n_obs   = 0;
try
    Tm = Tm(~isnan(Tm.conf_z) & ~isnan(Tm.(model_y)), :);
    n_obs = height(Tm);
    if n_obs < 20 || numel(unique(string(Tm.subj_id))) < 4
        sig_str = 'too few rows'; return;
    end
    Tm.subj_id = categorical(string(Tm.subj_id));
    Tm.stage   = categorical(string(Tm.stage), {'LN','LE','RN','RE'}, 'Ordinal', false);
    mdl = fitlme(Tm, [model_y ' ~ conf_z + stage + (1|subj_id)'], 'FitMethod','REML');
    C   = mdl.Coefficients;
    idx = find(strcmp(C.Name,'conf_z'),1);
    if isempty(idx), sig_str = 'term not found'; return; end
    b   = C.Estimate(idx);
    p   = C.pValue(idx);
    sig_str = sprintf('conf_z: b=%.3f, %s  [%s block]', b, p_to_stars(p), bt_label);
catch ME
    sig_str = sprintf('LME: %s', ME.message(1:min(50,end)));
end
end

% ─────────────────────────────────────────────────────────────────────────────
function p = rq2_coef_p(C, term_name, want_interaction)
%RQ2_COEF_P  Extract p-value for a named term; interaction if requested.
p = NaN;
names = string(C.Name);
if want_interaction
    idx = find(contains(names, term_name) & contains(names, ':'), 1);
else
    idx = find(strcmp(names, term_name), 1);
    if isempty(idx)
        idx = find(contains(names, term_name) & ~contains(names, ':'), 1);
    end
end
if ~isempty(idx), p = C.pValue(idx); end
end

% ─────────────────────────────────────────────────────────────────────────────
function add_sig_text(ax, sig_str, n_obs)
%ADD_SIG_TEXT  Compact annotation in top-left corner.
if nargin >= 3 && ~isempty(n_obs) && n_obs > 0
    txt = sprintf('n=%d\n%s', n_obs, sig_str);
else
    txt = sig_str;
end
text(ax, 0.03, 0.97, txt, 'Units','normalized','VerticalAlignment','top', ...
    'HorizontalAlignment','left','FontSize',7,'BackgroundColor','w', ...
    'Margin',2,'Interpreter','none','Clipping','off');
end

% ─────────────────────────────────────────────────────────────────────────────
function rq2_axis_labels(ax, xlbl, ylbl, reverse_y)
xlabel(ax, xlbl,'Interpreter','none');
ylabel(ax, ylbl,'Interpreter','none');
if reverse_y, set(ax,'YDir','reverse'); end
end

% ─────────────────────────────────────────────────────────────────────────────
function s = p_to_stars(p)
if isnan(p),      s = 'n/a';
elseif p < 0.001, s = '***';
elseif p < 0.01,  s = '**';
elseif p < 0.05,  s = '*';
elseif p < 0.10,  s = '†';
else,             s = 'ns';
end
end

% ─────────────────────────────────────────────────────────────────────────────
function save_lme_csv(mdl, fpath, model_desc)
%SAVE_LME_CSV  Write coefficient table + fit indices to CSV.
try
    C = mdl.Coefficients;
    C.model = repmat(string(model_desc), height(C), 1);
    C.AIC   = repmat(mdl.ModelCriterion.AIC, height(C), 1);
    C.BIC   = repmat(mdl.ModelCriterion.BIC, height(C), 1);
    writetable(C, fpath);
    fprintf('  LME table saved: %s\n', fpath);
catch ME
    fprintf('  LME CSV save failed: %s\n', ME.message);
end
end

% ─────────────────────────────────────────────────────────────────────────────
% ─────────────────────────────────────────────────────────────────────────────
function T = rq2_add_transition_cols(T)
%RQ2_ADD_TRANSITION_COLS  Self-contained version of add_transition_history_columns.
%
% Adds three columns to T:
%   transition_recent — "first", "D→D", "D→P", "P→D", "P→P", or "unknown"
%   n_prev_P          — count of prior P blocks within the subject sequence
%   n_prev_P_bin      — "0", "1", "2+"

if ~ismember('block_number', T.Properties.VariableNames)
    if ismember('block', T.Properties.VariableNames)
        T.block_number = double(T.block);
    elseif ismember('blocknum', T.Properties.VariableNames)
        T.block_number = double(T.blocknum);
    else
        error('rq2_add_transition_cols: no block column found (tried block_number, block, blocknum).');
    end
end

% Normalise block_type to string; recode V→P
bt_s = string(T.block_type);
bt_s(bt_s == "V") = "P";

T.transition_recent = repmat("first",   height(T), 1);
T.n_prev_P          = nan(height(T), 1);
T.n_prev_P_bin      = repmat("0",       height(T), 1);

subs = unique(string(T.subj_id));

for si = 1:numel(subs)
    sm   = string(T.subj_id) == subs(si);
    blks = unique(T.block_number(sm & ~isnan(T.block_number)));
    blks = sort(blks(:)');

    prev_types = strings(0,1);

    for bi = 1:numel(blks)
        bm = sm & T.block_number == blks(bi);
        if ~any(bm), continue; end

        % Determine current block type from first non-empty entry in this block
        curr_vals = bt_s(bm);
        curr_vals = curr_vals(curr_vals == "D" | curr_vals == "P");
        if isempty(curr_vals)
            curr = "";
        else
            curr = curr_vals(1);
        end

        % n_prev_P
        n_prev = sum(prev_types == "P");
        T.n_prev_P(bm)    = n_prev;
        if n_prev >= 2
            T.n_prev_P_bin(bm) = "2+";
        else
            T.n_prev_P_bin(bm) = string(n_prev);
        end

        % Transition label
        if bi == 1 || isempty(prev_types)
            T.transition_recent(bm) = "first";
        else
            prev = prev_types(end);
            if (prev == "D" || prev == "P") && (curr == "D" || curr == "P")
                T.transition_recent(bm) = prev + "→" + curr;
            else
                T.transition_recent(bm) = "unknown";
            end
        end

        % Update history
        if curr == "D" || curr == "P"
            prev_types(end+1, 1) = curr; %#ok<AGROW>
        end
    end
end

try
    exportgraphics(fig, fullfile(outdir,[fname '.pdf']), 'ContentType','vector');
    exportgraphics(fig, fullfile(outdir,[fname '.png']), 'Resolution',300);
catch
    saveas(fig, fullfile(outdir,[fname '.pdf']));
    saveas(fig, fullfile(outdir,[fname '.png']));
end
fprintf('  Saved: %s\n', fname);
end