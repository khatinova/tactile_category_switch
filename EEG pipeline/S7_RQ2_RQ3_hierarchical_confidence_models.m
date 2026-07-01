% =========================================================================
% S7_RQ2_RQ3_hierarchical_confidence_models.m
%
% Hierarchical analysis of confidence → EEG (FRN proxy, theta).
%
% PLOTTING PHILOSOPHY — hierarchical, mirroring the model sequence:
%   Each figure corresponds to exactly one model step. The sequence goes
%   from global to specific, adding one term at a time:
%
%   Fig 1  M0  — Overall scatter: confidence → component (all trials, no splits)
%   Fig 2  M1  — Split by OUTCOME (correct vs incorrect)
%   Fig 3  M2  — Split by BLOCK TYPE (D vs P)
%   Fig 4  M3  — 2×2: outcome × block type interaction
%   Fig 5  M4  — P-block only: true vs false feedback × outcome
%   Fig 6  M5  — Split by STAGE (LN/LE/RN/RE) — where in learning does
%                the confidence relationship change?
%   Fig 7  M6  — Split by N_PREV_P (0 / 1 / 2+) — does prior uncertainty
%                exposure modulate the relationship?
%   Fig 8  M7  — Split by RECENT TRANSITION (D→D / D→P / P→D / P→P)
%
%   Each figure panel shows:
%     - Grey background scatter (all individual trials)
%     - Quantile-binned means ± 95% CI
%     - Fitted OLS regression line with equation (b, SE, t, p) from the LME
%     - Significance annotation from the corresponding LME term
%
% MODELS (ML for comparison, REML for final inference):
%   M0  y_z ~ conf_z + (1|subj_id)
%   M1  y_z ~ conf_z + outcome + (1|subj_id)
%   M2  y_z ~ conf_z + outcome + block_type + (1|subj_id)
%   M3  y_z ~ conf_z + outcome*block_type + (1|subj_id)
%   M4  y_z ~ conf_z + outcome + block_type*false_fb + (1|subj_id)  [P-block only]
%   M5  y_z ~ conf_z * stage + (1|subj_id)
%   M6  y_z ~ conf_z * n_prev_P_bin + (1|subj_id)
%   M7  y_z ~ conf_z * transition_recent + (1|subj_id)
%
% COMPONENTS: prefrontal mean amp (FRN proxy) and frontal theta power.
%
% OUTPUT: one subfolder per component, one PDF per model step.
% =========================================================================

clearvars -except group_table gt
close all; clc

%% -----------------------------------------------------------------------
%  PATHS
%  -----------------------------------------------------------------------
base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch';
saved_tables_folder = fullfile(base_path, 'Salient mod switch KH', ...
    'Results', 'EEG analysis', 'Outcome_feature_tables_v4_merged');
figure_output_folder = fullfile(base_path, 'Salient mod switch KH', ...
    'Results', 'EEG analysis', 'Figures', 'S7_RQ2_RQ3_hierarchical');
if ~exist(figure_output_folder,'dir'), mkdir(figure_output_folder); end

% Load table
if ~exist('gt','var') && ~exist('group_table','var')
    load(fullfile(saved_tables_folder,'group_feature_table_combined.mat'),'group_table');
end
if exist('gt','var') && ~exist('group_table','var'), group_table = gt; end
gt = group_table;
%% EXCLUDE AUDITORY-FEEDBACK SUBJECTS (Ox01–Ox08 had auditory FB, distorts ERP results)
% Filter by feedback_modality if available, otherwise by subject ID pattern
if ismember('feedback_modality',gt.Properties.VariableNames)
    aud_mask = gt.feedback_modality=='auditory' | gt.feedback_modality=='Auditory';
    if any(aud_mask)
        fprintf('Excluding %d trials from %d auditory-feedback subjects.\n',...
            sum(aud_mask), numel(unique(gt.subj_id(aud_mask))));
        gt = gt(~aud_mask,:);
    end
else
    % Fallback: exclude subjects with ID number < 10 (Ox01-Ox08 pattern)
    sid_str = string(gt.subj_id);
    % Extract numeric part from IDs like 'Ox01', 'Ox10', 'sub01', etc.
    nums = regexp(sid_str, '\d+', 'match', 'once');
    nums_d = str2double(nums);
    exclude_mask = nums_d < 10 & ~isnan(nums_d);
    if any(exclude_mask)
        excl_subs = unique(gt.subj_id(exclude_mask));
        fprintf('Excluding %d early subjects (auditory FB): %s\n',...
            numel(excl_subs), strjoin(string(excl_subs),', '));
        gt = gt(~exclude_mask,:);
    end
end
% Rebuild subject list after exclusion
subj_list=unique(gt.subj_id); N_subj=numel(subj_list);
fprintf('Analysing %d subjects (visual/tactile feedback only).\n', N_subj);

%% -----------------------------------------------------------------------
%  FIGURE STYLE DEFAULTS
%  -----------------------------------------------------------------------
set(groot,'defaultAxesTickDir','out');
set(groot,'defaultAxesBox','off');
set(groot,'defaultAxesTickDirMode','manual');
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultAxesFontSize',11);

CLR_D    = [0.15 0.45 0.70];
CLR_P    = [0.80 0.30 0.10];
CLR_INC  = [0.75 0.20 0.18];
CLR_COR  = [0.20 0.58 0.25];
CLR_TFB  = [0.20 0.58 0.25];
CLR_FFB  = [0.75 0.35 0.15];
CLR_BASE = [0.30 0.30 0.30];
STAGE_COLS = [0.55 0.75 0.90; 0.15 0.45 0.70; 0.90 0.65 0.40; 0.80 0.30 0.10];
NPP_COLS   = [0.55 0.55 0.55; 0.65 0.42 0.15; 0.80 0.15 0.15];
TRANS_COLS = [0.12 0.47 0.71; 0.85 0.33 0.10; 0.47 0.67 0.19; 0.80 0.20 0.60];

%% -----------------------------------------------------------------------
%  PREPARE TABLE
%  -----------------------------------------------------------------------
gt.subj_id   = categorical(string(gt.subj_id));
gt.confidence = double(gt.confidence);

% Outcome
if ismember('correct',gt.Properties.VariableNames)
    gt.correct_num = double(gt.correct);
end
u = unique(gt.correct_num(~isnan(gt.correct_num)));
if all(ismember(u,[1 2])), gt.correct_num = gt.correct_num - 1; end
gt.outcome = categorical(gt.correct_num,[0 1],{'Incorrect','Correct'});

% Block type
bt = string(gt.block_type); bt(bt=="V") = "P";
gt.block_type = categorical(bt,{'D','P'});

% False feedback
if ismember('trueFB',gt.Properties.VariableNames)
    gt.false_fb = double(~logical(gt.trueFB));
elseif ~ismember('false_fb',gt.Properties.VariableNames)
    gt.false_fb = nan(height(gt),1);
end

% Stage
if ismember('stage',gt.Properties.VariableNames)
    gt.stage = categorical(string(gt.stage),{'LN','LE','RN','RE'},'Ordinal',false);
end

% n_prev_P_bin
if ismember('n_prev_P',gt.Properties.VariableNames) && ~ismember('n_prev_P_bin',gt.Properties.VariableNames)
    npp = double(gt.n_prev_P);
    bins = repmat("0",height(gt),1);
    bins(npp==1) = "1"; bins(npp>=2) = "2+"; bins(isnan(npp)) = "0";
    gt.n_prev_P_bin = categorical(bins,{'0','1','2+'});
elseif ismember('n_prev_P_bin',gt.Properties.VariableNames)
    gt.n_prev_P_bin = categorical(string(gt.n_prev_P_bin),{'0','1','2+'});
end

% transition_recent
if ismember('transition_recent',gt.Properties.VariableNames)
    tr_str = string(gt.transition_recent);
    tr_str(~ismember(tr_str,["D→D","D→P","P→D","P→P"])) = "first";
    gt.transition_recent = categorical(tr_str);
end

% Within-subject z-scores
gt.conf_z = zscore_ws(gt.confidence, gt.subj_id);

% Auto-select EEG columns
prefrontal_candidates = {'prefrontal_mean_amp','prefrontal_mean_norm',...
    'prefrontal_neg_peak_norm','prefrontal_neg_peak_amp','prefrontal_mean'};
theta_candidates = {'Theta_amp','theta_amp','Theta','theta','frontal_theta'};
fn_col   = first_col(gt, prefrontal_candidates);
th_col   = first_col(gt, theta_candidates);
if isempty(fn_col),  error('No prefrontal/FRN column found.'); end
if isempty(th_col),  error('No theta column found.'); end
[gt, fn_z]  = ensure_z(gt, fn_col);
[gt, th_z]  = ensure_z(gt, th_col);

fprintf('\nRQ2 column: %s (model: %s)\n', fn_col, fn_z);
fprintf('RQ3 column: %s (model: %s)\n', th_col, th_z);
fprintf('Output:     %s\n\n', figure_output_folder);

%% -----------------------------------------------------------------------
%  RUN HIERARCHICAL ANALYSIS FOR EACH COMPONENT
%  -----------------------------------------------------------------------
comps = {
    fn_col, fn_z, 'RQ2_PrefrontalMean', 'Prefrontal mean / FRN proxy', true;
    th_col, th_z, 'RQ3_Theta',          'Frontal theta power',         false;
};

all_results = struct();

for ci = 1:size(comps,1)
    y_plot   = comps{ci,1};
    y_model  = comps{ci,2};
    tag      = comps{ci,3};
    label    = comps{ci,4};
    rev_y    = comps{ci,5};

    comp_dir = fullfile(figure_output_folder, tag);
    if ~exist(comp_dir,'dir'), mkdir(comp_dir); end

    % Build analysis table (drop rows with any key column missing)
    keep = ~isnan(gt.(y_model)) & ~isnan(gt.(y_plot)) & ...
           ~isnan(gt.conf_z) & ~isnan(gt.correct_num) & ...
           ~isundefined(gt.block_type) & ~isundefined(gt.subj_id);
    T = gt(keep,:);
    T.y_model = double(T.(y_model));
    T.y_plot  = double(T.(y_plot));
    fprintf('\n=== %s ===\n', label);
    fprintf('  Rows: %d  |  Subjects: %d\n', height(T), numel(unique(T.subj_id)));

    % Fit all models
    M = fit_model_hierarchy(T, y_model, label);

    % Print comparisons
    print_hierarchy(M, label);

    % Save coefficient CSVs
    write_coef_tables(M, comp_dir, tag);

    % ─── FIG 1: M0 — overall effect, all trials ────────────────────────
    fig1 = figure('Position',[60 60 600 540],'Color','w');
    ax = axes(fig1); hold(ax,'on');
    hier_scatter(ax, T.conf_z, T.y_plot, CLR_BASE, 8, true);
    add_lme_line_and_stats(ax, M.M0, 'conf_z', T.conf_z, T.y_model, CLR_BASE);
    decorate(ax,'Confidence (z)','y_z','All trials','M0',label,rev_y);
    save_fig(fig1, comp_dir, sprintf('%s_Fig1_M0_overall',tag));

    % ─── FIG 2: M1 — split by outcome ──────────────────────────────────
    fig2 = figure('Position',[60 60 860 480],'Color','w');
    sgtitle_str(fig2, sprintf('%s: Fig 2 — M1: + Outcome',label));
    outcomes = {'Incorrect','Correct'}; oclrs = {CLR_INC, CLR_COR};
    for oi = 1:2
        ax = subplot(1,2,oi); hold(ax,'on');
        m = string(T.outcome)==outcomes{oi};
        Tsub = T(m,:); Tsub.y_model = double(Tsub.(y_model));
        hier_scatter(ax, Tsub.conf_z, Tsub.y_plot, oclrs{oi}, 8, true);
        add_lme_line_and_stats(ax, M.M1, 'conf_z', Tsub.conf_z, Tsub.y_model, oclrs{oi});
        decorate(ax,'Confidence (z)',label,outcomes{oi},'M1',label,rev_y);
    end
    % Annotation: overall outcome effect from M1
    add_model_annotation(fig2, M.M1, 'outcome_Correct', sprintf('M1: + outcome\n%s',label));
    save_fig(fig2, comp_dir, sprintf('%s_Fig2_M1_outcome',tag));

    % ─── FIG 3: M2 — split by block type ───────────────────────────────
    fig3 = figure('Position',[60 60 860 480],'Color','w');
    sgtitle_str(fig3, sprintf('%s: Fig 3 — M2: + Block type',label));
    bts = {'D','P'}; bclrs = {CLR_D, CLR_P};
    for bi = 1:2
        ax = subplot(1,2,bi); hold(ax,'on');
        m = string(T.block_type)==bts{bi};
        Tsub = T(m,:); Tsub.y_model = double(Tsub.(y_model));
        hier_scatter(ax, Tsub.conf_z, Tsub.y_plot, bclrs{bi}, 8, true);
        add_lme_line_and_stats(ax, M.M2, 'conf_z', Tsub.conf_z, Tsub.y_model, bclrs{bi});
        decorate(ax,'Confidence (z)',label,sprintf('%s blocks',bts{bi}),'M2',label,rev_y);
    end
    add_model_annotation(fig3, M.M2, 'block_typeP', sprintf('M2: + block_type\n%s',label));
    save_fig(fig3, comp_dir, sprintf('%s_Fig3_M2_blocktype',tag));

    % ─── FIG 4: M3 — 2×2 outcome × block type ──────────────────────────
    fig4 = figure('Position',[60 60 1100 820],'Color','w');
    sgtitle_str(fig4, sprintf('%s: Fig 4 — M3: outcome × block type interaction',label));
    panel = 0;
    for oi = 1:2
        for bi = 1:2
            panel = panel + 1;
            ax = subplot(2,2,panel); hold(ax,'on');
            m = string(T.outcome)==outcomes{oi} & string(T.block_type)==bts{bi};
            Tsub = T(m,:); Tsub.y_model = double(Tsub.(y_model));
            clr_mix = 0.5*(oclrs{oi} + bclrs{bi});
            hier_scatter(ax, Tsub.conf_z, Tsub.y_plot, clr_mix, 7, true);
            add_lme_line_and_stats(ax, M.M3, 'conf_z', Tsub.conf_z, Tsub.y_model, clr_mix);
            decorate(ax,'Confidence (z)',label,...
                sprintf('%s | %s blocks',outcomes{oi},bts{bi}),'M3',label,rev_y);
        end
    end
    % Interaction term annotation — extract conf_z:block_typeP or outcome:block_type
    add_interaction_annotation(fig4, M.M3, ...
        {'conf_z:block_typeP','outcome_Correct:block_typeP','outcome_CorrectDblock_typeP'}, ...
        sprintf('M3: conf_z + outcome*block_type\n%s',label));
    save_fig(fig4, comp_dir, sprintf('%s_Fig4_M3_outcome_x_blocktype',tag));

    % ─── FIG 5: M4 — P-block only: true vs false feedback × outcome ────
    has_ffb = any(T.false_fb==1,'all') && numel(unique(T.false_fb(~isnan(T.false_fb))))>1;
    if has_ffb
        Tp = T(string(T.block_type)=='P' & ~isnan(T.false_fb),:);
        Tp.y_model = double(Tp.(y_model));
        fig5 = figure('Position',[60 60 1100 820],'Color','w');
        sgtitle_str(fig5, sprintf('%s: Fig 5 — M4: P-block false feedback × outcome',label));
        fblabs = {'True FB','False FB'}; fbcols = {CLR_TFB, CLR_FFB};
        panel = 0;
        for oi = 1:2
            for fi = 0:1
                panel = panel + 1;
                ax = subplot(2,2,panel); hold(ax,'on');
                m = string(Tp.outcome)==outcomes{oi} & Tp.false_fb==fi;
                Tsub2 = Tp(m,:);
                clr = fbcols{fi+1};
                hier_scatter(ax, Tsub2.conf_z, Tsub2.y_plot, clr, 7, true);
                add_lme_line_and_stats(ax, M.M4, 'conf_z', Tsub2.conf_z, Tsub2.y_model, clr);
                decorate(ax,'Confidence (z)',label,...
                    sprintf('%s | %s',outcomes{oi},fblabs{fi+1}),'M4',label,rev_y);
                % n trials annotation
                text(ax,.97,.04,sprintf('n=%d',height(Tsub2)),'Units','normalized',...
                    'HorizontalAlignment','right','FontSize',8,'Color',[.5 .5 .5]);
            end
        end
        add_interaction_annotation(fig5, M.M4, ...
            {'false_fb','block_typeP:false_fb','conf_z:false_fb'}, ...
            sprintf('M4: conf_z + outcome + block_type*false_fb (P only)\n%s',label));
        save_fig(fig5, comp_dir, sprintf('%s_Fig5_M4_falsefb_outcome',tag));
    else
        fprintf('  Fig5 M4 skipped: no false feedback variation.\n');
    end

    % ─── FIG 6: M5 — split by STAGE ────────────────────────────────────
    if ismember('stage',T.Properties.VariableNames) && ...
            numel(unique(string(T.stage(~isundefined(T.stage))))) > 1
        stages = {'LN','LE','RN','RE'};
        fig6 = figure('Position',[60 60 1400 680],'Color','w');
        sgtitle_str(fig6, sprintf('%s: Fig 6 — M5: confidence × stage',label));
        for si_s = 1:4
            ax = subplot(2,4,si_s); hold(ax,'on');
            m = string(T.stage)==stages{si_s} & string(T.block_type)=='D';
            Tsub = T(m,:); Tsub.y_model = double(Tsub.(y_model));
            hier_scatter(ax, Tsub.conf_z, Tsub.y_plot, STAGE_COLS(si_s,:), 7, false);
            add_lme_line_and_stats(ax, M.M5, 'conf_z', Tsub.conf_z, Tsub.y_model, STAGE_COLS(si_s,:));
            decorate(ax,'Conf (z)',label,sprintf('%s | D',stages{si_s}),'M5',label,rev_y);

            ax = subplot(2,4,si_s+4); hold(ax,'on');
            m = string(T.stage)==stages{si_s} & string(T.block_type)=='P';
            Tsub = T(m,:); Tsub.y_model = double(Tsub.(y_model));
            hier_scatter(ax, Tsub.conf_z, Tsub.y_plot, STAGE_COLS(si_s,:), 7, false);
            add_lme_line_and_stats(ax, M.M5, 'conf_z', Tsub.conf_z, Tsub.y_model, STAGE_COLS(si_s,:));
            decorate(ax,'Conf (z)',label,sprintf('%s | P',stages{si_s}),'M5',label,rev_y);
        end
        add_interaction_annotation(fig6, M.M5, ...
            {'conf_z:stage_LE','conf_z:stage_RN','conf_z:stage_RE'}, ...
            sprintf('M5: conf_z * stage\n%s',label));
        save_fig(fig6, comp_dir, sprintf('%s_Fig6_M5_stage',tag));
    else
        fprintf('  Fig6 M5 skipped: stage column missing or single-valued.\n');
    end

    % ─── FIG 7: M6 — split by N_PREV_P ─────────────────────────────────
    if ismember('n_prev_P_bin',T.Properties.VariableNames) && ...
            numel(unique(string(T.n_prev_P_bin(~isundefined(T.n_prev_P_bin))))) > 1
        npp_labs = {'0','1','2+'};
        fig7 = figure('Position',[60 60 1300 480],'Color','w');
        sgtitle_str(fig7, sprintf('%s: Fig 7 — M6: confidence × n\\_prev\\_P',label));
        for ni = 1:3
            ax = subplot(1,3,ni); hold(ax,'on');
            m = string(T.n_prev_P_bin)==npp_labs{ni};
            Tsub = T(m,:); Tsub.y_model = double(Tsub.(y_model));
            hier_scatter(ax, Tsub.conf_z, Tsub.y_plot, NPP_COLS(ni,:), 7, true);
            add_lme_line_and_stats(ax, M.M6, 'conf_z', Tsub.conf_z, Tsub.y_model, NPP_COLS(ni,:));
            decorate(ax,'Confidence (z)',label,...
                sprintf('Prior P blocks: %s',npp_labs{ni}),'M6',label,rev_y);
            text(ax,.97,.96,sprintf('n=%d',height(Tsub)),'Units','normalized',...
                'HorizontalAlignment','right','FontSize',8,'Color',[.5 .5 .5],'VerticalAlignment','top');
        end
        add_interaction_annotation(fig7, M.M6, ...
            {'conf_z:n_prev_P_bin_1','conf_z:n_prev_P_bin_2_'}, ...
            sprintf('M6: conf_z * n\\_prev\\_P\\_bin\n%s',label));
        save_fig(fig7, comp_dir, sprintf('%s_Fig7_M6_nPrevP',tag));
    else
        fprintf('  Fig7 M6 skipped: n_prev_P_bin missing or single-valued.\n');
    end

    % ─── FIG 8: M7 — split by RECENT TRANSITION ─────────────────────────
    if ismember('transition_recent',T.Properties.VariableNames) && ...
            numel(unique(string(T.transition_recent(~isundefined(T.transition_recent))))) > 1
        trans_vals = {'D→D','D→P','P→D','P→P'};
        trans_present = trans_vals(cellfun(@(v) ...
            any(string(T.transition_recent)==v), trans_vals));
        n_tr = numel(trans_present);
        fig8 = figure('Position',[60 60 min(400*n_tr,1600) 480],'Color','w');
        sgtitle_str(fig8, sprintf('%s: Fig 8 — M7: confidence × recent transition',label));
        for ti = 1:n_tr
            ax = subplot(1,n_tr,ti); hold(ax,'on');
            clr_idx = find(strcmp(trans_vals,trans_present{ti}),1);
            if isempty(clr_idx), clr_idx = 1; end
            m = string(T.transition_recent)==trans_present{ti};
            Tsub = T(m,:); Tsub.y_model = double(Tsub.(y_model));
            hier_scatter(ax, Tsub.conf_z, Tsub.y_plot, TRANS_COLS(clr_idx,:), 7, true);
            add_lme_line_and_stats(ax, M.M7, 'conf_z', Tsub.conf_z, Tsub.y_model, TRANS_COLS(clr_idx,:));
            decorate(ax,'Confidence (z)',label,...
                strrep(trans_present{ti},'→','\rightarrow'),'M7',label,rev_y);
            text(ax,.97,.96,sprintf('n=%d',height(Tsub)),'Units','normalized',...
                'HorizontalAlignment','right','FontSize',8,'Color',[.5 .5 .5],'VerticalAlignment','top');
        end
        add_interaction_annotation(fig8, M.M7, ...
            {'conf_z:transition_recent_D.P','conf_z:transition_recent_P.D','conf_z:transition_recent_P.P'}, ...
            sprintf('M7: conf_z * recent transition\n%s',label));
        save_fig(fig8, comp_dir, sprintf('%s_Fig8_M7_transition',tag));
    else
        fprintf('  Fig8 M7 skipped: transition_recent missing or single-valued.\n');
    end

    % Store results
    all_results.(tag) = M;
    fprintf('\n  All figures saved to: %s\n', comp_dir);
end % component loop

% Save model results
save(fullfile(figure_output_folder,'S7_RQ2_RQ3_hierarchical_results.mat'), 'all_results');
fprintf('\n=== S7_RQ2_RQ3 complete. ===\n');
fprintf('Results saved to: %s\n', figure_output_folder);


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

% ── Table/column utilities ─────────────────────────────────────────────────

function col = first_col(T, candidates)
col = '';
for k = 1:numel(candidates)
    if ismember(candidates{k}, T.Properties.VariableNames)
        col = candidates{k}; return;
    end
end
end

function [T, z_col] = ensure_z(T, raw_col)
z_col = [raw_col '_z'];
if ismember(z_col, T.Properties.VariableNames), return; end
T.(z_col) = zscore_ws(T.(raw_col), T.subj_id);
end

function z = zscore_ws(x, subj)
x = double(x); z = nan(size(x));
subs = categories(categorical(subj));
for si = 1:numel(subs)
    m = categorical(subj)==subs{si};
    v = x(m); sd = std(v,'omitnan');
    if sum(~isnan(v))>1 && sd>0
        z(m) = (v - mean(v,'omitnan')) ./ sd;
    end
end
end

% ── Model fitting ──────────────────────────────────────────────────────────

function M = fit_model_hierarchy(T, ~, label)
% Fit M0–M7 with ML (for comparison) and REML (for final estimates).
% Returns struct of fitted models; failed models set to [].

fprintf('\n--- Fitting hierarchy for: %s ---\n', label);

% Stage: need ordinal categorical without 'Ordinal' flag for fitlme interaction
if ismember('stage',T.Properties.VariableNames)
    T.stage = categorical(string(T.stage),{'LN','LE','RN','RE'},'Ordinal',false);
end
if ismember('n_prev_P_bin',T.Properties.VariableNames)
    T.n_prev_P_bin = categorical(string(T.n_prev_P_bin),{'0','1','2+'});
end
if ismember('transition_recent',T.Properties.VariableNames)
    tr = string(T.transition_recent);
    valid = ismember(tr,["D→D","D→P","P→D","P→P"]);
    tr(~valid) = "D→D"; % replace 'first'/other with most common
    T.transition_recent = categorical(tr);
end

specs = {
    'M0', 'y_model ~ conf_z + (1|subj_id)';
    'M1', 'y_model ~ conf_z + outcome + (1|subj_id)';
    'M2', 'y_model ~ conf_z + outcome + block_type + (1|subj_id)';
    'M3', 'y_model ~ conf_z + outcome*block_type + (1|subj_id)';
    'M5', 'y_model ~ conf_z * stage + (1|subj_id)';
    'M6', 'y_model ~ conf_z * n_prev_P_bin + (1|subj_id)';
    'M7', 'y_model ~ conf_z * transition_recent + (1|subj_id)';
};

% M4 uses P-block subset only
M4_formula = 'y_model ~ conf_z + outcome + block_type*false_fb + (1|subj_id)';

M = struct();
for k = 1:size(specs,1)
    nm = specs{k,1}; f = specs{k,2};
    % Check required variables exist
    needed = regexp(f,'[a-zA-Z_][a-zA-Z0-9_]*','match');
    skip = false;
    for ni = 1:numel(needed)
        if ~ismember(needed{ni}, T.Properties.VariableNames) && ...
           ~ismember(needed{ni}, {'y_model','subj_id','conf_z','outcome','block_type','false_fb'})
            skip = true; break;
        end
    end
    if skip
        M.(nm) = []; fprintf('  %s skipped (missing variable)\n', nm);
        continue;
    end
    M.(nm) = fit_safe(T, f, nm, 'ML');
end

% M4: P-block subset
Tp = T(string(T.block_type)=='P' & ~isnan(T.false_fb),:);
has_ffb = ~isempty(Tp) && any(Tp.false_fb==1) && numel(unique(Tp.false_fb))>1;
if has_ffb
    M.M4 = fit_safe(Tp, M4_formula, 'M4', 'ML');
else
    M.M4 = []; fprintf('  M4 skipped (no false FB variation in P blocks)\n');
end
end

function mdl = fit_safe(T, formula, name, method)
try
    mdl = fitlme(T, formula, 'FitMethod', method);
    fprintf('  %s OK  AIC=%.1f  LogL=%.1f\n', name, mdl.ModelCriterion.AIC, mdl.LogLikelihood);
catch ME
    % Try dropping random slope if singular
    f2 = regexprep(formula, '\(1 \+ \w+ \| subj_id\)', '(1|subj_id)');
    try
        mdl = fitlme(T, f2, 'FitMethod', method);
        fprintf('  %s (simplified RE) OK  AIC=%.1f\n', name, mdl.ModelCriterion.AIC);
    catch ME2
        fprintf('  %s FAILED: %s\n', name, ME2.message);
        mdl = [];
    end
end
end

% ── Hierarchical comparison printer ───────────────────────────────────────

function print_hierarchy(M, label)
fprintf('\n--- Model comparisons: %s ---\n', label);
pairs = {'M0','M1','does outcome improve fit?';
         'M1','M2','does block_type improve fit?';
         'M2','M3','does outcome×block_type interaction improve fit?';
         'M2','M5','does stage modulate the conf relationship?';
         'M2','M6','does n_prev_P modulate the conf relationship?';
         'M2','M7','does recent transition modulate the conf relationship?'};
for k = 1:size(pairs,1)
    m0 = M.(pairs{k,1}); m1 = M.(pairs{k,2});
    fprintf('\n%s vs %s — %s\n', pairs{k,1}, pairs{k,2}, pairs{k,3});
    if isempty(m0)||isempty(m1), fprintf('  One model unavailable.\n'); continue; end
    try
        cmp = compare(m0, m1, 'CheckNesting', false);
        fprintf('  χ²(%.0f)=%.3f  p=%.4f  ΔAIC=%.2f\n', ...
            cmp.deltaDF(2), cmp.LRStat(2), cmp.pValue(2), ...
            m1.ModelCriterion.AIC - m0.ModelCriterion.AIC);
    catch ME
        fprintf('  Comparison failed: %s\n', ME.message);
    end
end
if ~isempty(M.M4)
    fprintf('\nM4 (P-block only) fixed effects:\n');
    print_coefs(M.M4);
end
end

function print_coefs(mdl)
if isempty(mdl), return; end
C = mdl.Coefficients;
for k = 1:height(C)
    fprintf('  %-40s  b=%+.4f  SE=%.4f  t=%+.3f  p=%s %s\n', ...
        string(C.Name(k)), C.Estimate(k), C.SE(k), C.tStat(k), ...
        pfmt(C.pValue(k)), pstars(C.pValue(k)));
end
end

% ── CSV saving ────────────────────────────────────────────────────────────

function write_coef_tables(M, comp_dir, tag)
fnames = fieldnames(M);
for k = 1:numel(fnames)
    nm = fnames{k};
    if isempty(M.(nm)), continue; end
    try
        C = M.(nm).Coefficients;
        C.stars = arrayfun(@pstars, C.pValue, 'UniformOutput', false);
        writetable(C, fullfile(comp_dir, sprintf('%s_%s_coefficients.csv',tag,nm)));
    catch
    end
end
end

% ── Plotting helpers ──────────────────────────────────────────────────────

function hier_scatter(ax, x, y, clr, sz, show_scatter)
% Background scatter + quantile-binned means ± 95% CI
x = double(x); y = double(y);
ok = ~isnan(x) & ~isnan(y);
x = x(ok); y = y(ok);
if isempty(x), return; end
if show_scatter
    scatter(ax, x, y, sz, [0.65 0.65 0.65], 'filled', ...
        'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0, 'HandleVisibility','off');
end
% 6 quantile bins
edges = unique(quantile(x, linspace(0,1,7)));
nb = numel(edges)-1;
if nb < 2, return; end
for bi = 1:nb
    if bi < nb, bm = x>=edges(bi) & x<edges(bi+1);
    else,       bm = x>=edges(bi) & x<=edges(bi+1); end
    xv = x(bm); yv = y(bm); yv = yv(~isnan(yv));
    if numel(yv)<2, continue; end
    xm = mean(xv,'omitnan'); ym = mean(yv,'omitnan');
    ci = 1.96*std(yv,'omitnan')/sqrt(numel(yv));
    errorbar(ax, xm, ym, ci, 'o', 'Color',clr, 'MarkerFaceColor',clr, ...
        'MarkerEdgeColor','k', 'MarkerSize',7, 'LineWidth',1.5, ...
        'CapSize',5, 'HandleVisibility','off');
end
end

function add_lme_line_and_stats(ax, mdl, term, x_data, y_data, clr)
% Draw regression line from LME fixed-effect slope + print equation
x_data = double(x_data); y_data = double(y_data);
ok = ~isnan(x_data) & ~isnan(y_data);
if isempty(mdl) || sum(ok) < 5, return; end
C = mdl.Coefficients;
% Find conf_z slope
idx = find(strcmp(string(C.Name), term), 1);
if isempty(idx), idx = find(contains(string(C.Name), 'conf_z') & ...
    ~contains(string(C.Name),':'), 1); end
if isempty(idx), return; end
b  = C.Estimate(idx);
se = C.SE(idx);
p  = C.pValue(idx);
% Intercept
i0 = find(strcmp(string(C.Name),'(Intercept)'),1);
b0 = C.Estimate(i0);
% Regression line over x range
xl = [min(x_data(ok)) max(x_data(ok))];
yl = b0 + b*xl;
lw = 2.0 + (p < 0.05)*0.8; % thicker if significant
ls = ternary_ls(p);
plot(ax, xl, yl, ls, 'Color', clr, 'LineWidth', lw, 'HandleVisibility','off');
% Stats box
eq_str = sprintf('b=%+.3f (SE=%.3f)\nt=%.2f, p=%s %s', b, se, C.tStat(idx), pfmt(p), pstars(p));
text(ax, 0.04, 0.96, eq_str, 'Units','normalized', 'VerticalAlignment','top', ...
    'FontSize', 8, 'BackgroundColor',[1 1 1 0.85], 'EdgeColor',[0.8 0.8 0.8], ...
    'FontName','Arial', 'Interpreter','none');
end

function decorate(ax, xlbl, ylbl, panel_title, model_tag, ~, rev_y)
xlabel(ax, xlbl, 'FontSize',10);
ylabel(ax, ylbl, 'FontSize',10, 'Interpreter','none');
title(ax, panel_title, 'FontSize',10, 'FontWeight','bold');
subtitle(ax, model_tag, 'FontSize',8, 'Color',[0.4 0.4 0.4]);
set(ax, 'TickDir','out','Box','off','LineWidth',1,'TickLength',[0.012 0.012]);
xline(ax, 0, ':', 'Color',[0.6 0.6 0.6], 'LineWidth',0.8, 'HandleVisibility','off');
yline(ax, 0, ':', 'Color',[0.6 0.6 0.6], 'LineWidth',0.8, 'HandleVisibility','off');
if rev_y, set(ax,'YDir','reverse'); end
axis(ax,'square');
end

function add_model_annotation(fig, mdl, term_hint, caption)
% Single text box at bottom of figure with model stats for a key term
if isempty(mdl), return; end
C = mdl.Coefficients;
nm = string(C.Name);
idx = find(contains(nm, term_hint), 1);
if isempty(idx), return; end
b = C.Estimate(idx); p = C.pValue(idx);
txt = sprintf('%s  |  term: %s  b=%+.3f  p=%s %s', ...
    caption, nm(idx), b, pfmt(p), pstars(p));
annotation(fig,'textbox',[0.01 0.01 0.98 0.04], 'String', txt, ...
    'FontSize',8, 'EdgeColor','none', ...
    'BackgroundColor',[0.96 0.96 0.96], 'Interpreter','none');
end

function add_interaction_annotation(fig, mdl, term_hints, caption)
% Annotate interaction terms — try each hint, report the first found
if isempty(mdl), return; end
C = mdl.Coefficients; nm = string(C.Name);
parts = {};
for k = 1:numel(term_hints)
    idx = find(contains(nm, term_hints{k}), 1);
    if isempty(idx), continue; end
    p = C.pValue(idx);
    parts{end+1} = sprintf('%s b=%+.3f p=%s %s', nm(idx), C.Estimate(idx), pfmt(p), pstars(p)); %#ok
end
if isempty(parts), parts = {'(interaction terms not found)'}; end
txt = [caption '   |   ' strjoin(parts,'   ')];
annotation(fig,'textbox',[0.01 0.01 0.98 0.04], 'String', txt, ...
    'FontSize',8, 'EdgeColor','none', ...
    'BackgroundColor',[0.96 0.96 0.96], 'Interpreter','none');
end

function sgtitle_str(fig, str)
sgtitle(fig, str, 'FontSize',12, 'FontWeight','bold', 'Interpreter','none');
end

function save_fig(fig, outdir, fname)
set(fig, 'Color','w');
ax_all = findall(fig,'Type','axes');
for k = 1:numel(ax_all)
    set(ax_all(k),'TickDir','out','Box','off','FontName','Arial','FontSize',10,'LineWidth',1);
end
if exist('apply_fig_style','file'), apply_fig_style(fig); end
exportgraphics(fig, fullfile(outdir, [fname '.pdf']), 'ContentType','vector');
exportgraphics(fig, fullfile(outdir, [fname '.png']), 'Resolution',300);
fprintf('    Saved: %s\n', fname);
% close(fig);
end

% ── Stat formatting ───────────────────────────────────────────────────────

function s = pfmt(p)
if isnan(p),      s = 'n/a';
elseif p < 0.001, s = '<.001';
else,             s = sprintf('=%.3f', p);
end
end

function s = pstars(p)
if isnan(p),      s = '';
elseif p < 0.001, s = '***';
elseif p < 0.01,  s = '**';
elseif p < 0.05,  s = '*';
elseif p < 0.10,  s = '†';
else,             s = 'ns';
end
end

function ls = ternary_ls(p)
if p < 0.05, ls = '-';
else,        ls = '--';
end
end
