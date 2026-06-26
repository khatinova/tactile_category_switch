% ==========================================================================
% plot_RL_params_and_transitions.m
%
% Standalone plotting script. Run AFTER:
%   1. extract_revaligned_alltrialdata_v6.m   → all_trial_data, behav_table
%   2. bayesian_delta_rule_dimshift.m          → model_params_and_sim.mat
%                                                + results struct in workspace
%
% FIGURES PRODUCED
% ────────────────
%  Fig 1  — Nassar H and β: D blocks vs P blocks (paired + distributions)
%  Fig 2  — Noise sensitivity ΔH = H_prob − H_det  (individual + group)
%  Fig 3  — Per-trial model latents (α, ω, |δ|) split by block type
%  Fig 4  — Reversal-aligned ACCURACY by transition type (D→D / D→P / P→D / P→P)
%  Fig 5  — Reversal-aligned CONFIDENCE by transition type
%  Fig 6  — Reversal-aligned CONFIDENCE pooled across ALL blocks (D + P)
%  Fig 7  — Pre vs post reversal accuracy × confidence  (2×2 block-type × rev-state)
%
% EXPECTED WORKSPACE VARIABLES (loaded below if not present)
%  group_T          — long-format behavioural table from extract_revaligned_v6
%  all_trial_data   — per-subject struct with aligned_correct / aligned_confidence
%  results          — per-subject Nassar fit struct from bayesian_delta_rule_dimshift
% ==========================================================================

close all;
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'off');


%% ─── PATHS ──────────────────────────────────────────────────────────────
remote = 0;
switch remote
    case 1
        base_path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch/Salient mod switch KH';
    case 0
        base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH';
    case 2
        base_path = 'Z:\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH';
end
data_path = fullfile(base_path, 'Data');
outpath   = fullfile(base_path, 'Results', 'Simulation results', 'Figures');
if ~exist(outpath,'dir'), mkdir(outpath); end

%% ─── LOAD DATA ───────────────────────────────────────────────────────────
if ~exist('group_T','var') || isempty(group_T)
    fprintf('Loading behav_table.mat...\n');
    load(fullfile(data_path, 'behav_table.mat'), 'group_T');
end
if ~exist('all_trial_data','var') || isempty(all_trial_data)
    fprintf('Loading all_trial_data.mat...\n');
    load(fullfile(data_path, 'all_trial_data.mat'), 'all_trial_data');
end
if ~exist('results','var') || isempty(results)
    fprintf('Loading model_params_and_sim.mat...\n');
    sim_file = fullfile(outpath, 'model_params_and_sim.mat');
    if exist(sim_file,'file')
        load(sim_file);   % loads model_vars, sim_data, params
        fprintf('  Note: ''results'' (Nassar fitted params) not found. Figs 1–3 will be skipped.\n');
        results = [];
    else
        results = [];
    end
end

%% ─── COLOUR SCHEME ───────────────────────────────────────────────────────
CLR_D      = [0.15 0.45 0.70];   % deterministic — blue
CLR_P      = [0.80 0.30 0.10];   % probabilistic — orange
CLR_KH     = [0.15 0.45 0.70];   % Ox cohort
CLR_RR     = [0.80 0.30 0.10];   % Nc cohort
CLR_PRE    = [0.50 0.50 0.50];   % pre-reversal — grey
CLR_POST   = [0.20 0.63 0.17];   % post-reversal — green

% Transition colours
TRANS_TYPES  = {'D→D','D→P','P→D','P→P'};
TRANS_COLORS = {[0.12 0.47 0.71], [0.85 0.33 0.10], [0.47 0.67 0.19], [0.80 0.20 0.60]};

preN   = 30;
postN  = 30;
alignedLen = preN + postN;
rel_ax = -preN : (postN-1);   % −30 … +29

%% ─── SUBJECT LISTS ───────────────────────────────────────────────────────
subj_ids_td = fieldnames(all_trial_data);

if ~isempty(results)
    subj_ids_fit = fieldnames(results);
else
    subj_ids_fit = {};
end

is_kh_td = cellfun(@(s) startsWith(s,'Ox'), subj_ids_td);
is_rr_td = cellfun(@(s) startsWith(s,'Nc'), subj_ids_td);

% =========================================================================
%% FIG 1 — NASSAR H AND β:  D-BLOCKS  vs  P-BLOCKS
% =========================================================================
if ~isempty(results)

    N_fit = numel(subj_ids_fit);
    H_det  = nan(N_fit,1); H_prob  = nan(N_fit,1);
    b_det  = nan(N_fit,1); b_prob  = nan(N_fit,1);
    H_ns   = nan(N_fit,1); H_all   = nan(N_fit,1); b_all = nan(N_fit,1);
    is_kh_fit = false(N_fit,1); is_rr_fit = false(N_fit,1);

    for si = 1:N_fit
        sn = subj_ids_fit{si};
        r  = results.(sn);
        H_all(si)  = r.H_fit;
        b_all(si)  = r.beta_fit;
        is_kh_fit(si) = startsWith(sn,'Ox');
        is_rr_fit(si) = startsWith(sn,'Nc');
        if isfield(r,'H_fit_det'),          H_det(si)  = r.H_fit_det;  end
        if isfield(r,'H_fit_prob'),         H_prob(si) = r.H_fit_prob; end
        if isfield(r,'beta_fit_det'),       b_det(si)  = r.beta_fit_det;  end
        if isfield(r,'beta_fit_prob'),      b_prob(si) = r.beta_fit_prob; end
        if isfield(r,'H_noise_sensitivity'),H_ns(si)   = r.H_noise_sensitivity; end
    end

    fig1 = figure('Position',[50 50 1400 500]);
    sgtitle('Nassar model parameters: Deterministic vs Probabilistic blocks','FontSize',12);

    % ── 1a: H paired ────────────────────────────────────────────────────
    ax1a = subplot(1,4,1); hold(ax1a,'on');
    title(ax1a,'Hazard rate H','FontSize',10);
    ok_h = ~isnan(H_det) & ~isnan(H_prob);
    for si = 1:N_fit
        if ~ok_h(si), continue; end
        lc = ternary(is_kh_fit(si), CLR_KH, CLR_RR);
        plot(ax1a,[1 2],[H_det(si) H_prob(si)],'-o','Color',[lc 0.35],'MarkerSize',5,'MarkerFaceColor',lc);
    end
    % Means ± SEM
    for xi_pair = {1, 2}
        xi = xi_pair{1};
        vals = ternary(xi==1, H_det(ok_h), H_prob(ok_h));
        errorbar(ax1a, xi, mean(vals,'omitnan'), sem(vals), ...
            'ko','MarkerSize',10,'MarkerFaceColor','k','LineWidth',2,'CapSize',8,'HandleVisibility','off');
    end
    [~,p_h] = ttest(H_det(ok_h), H_prob(ok_h));
    text(ax1a,0.5,0.97,sprintf('Paired t-test: p=%.3f',p_h),'Units','normalized', ...
        'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8,'BackgroundColor','w');
    set(ax1a,'XTick',[1 2],'XTickLabel',{'Deterministic','Probabilistic'});
    ylabel(ax1a,'H (hazard rate)'); xlim(ax1a,[0.5 2.5]); ylim(ax1a,[0 max([H_det;H_prob],[],'omitnan')*1.2]);
    add_significance_bracket(ax1a, 1, 2, max([H_det;H_prob],[],'omitnan')*1.1, p_h);

    % ── 1b: β paired ────────────────────────────────────────────────────
    ax1b = subplot(1,4,2); hold(ax1b,'on');
    title(ax1b,'Decision noise β','FontSize',10);
    ok_b = ~isnan(b_det) & ~isnan(b_prob);
    for si = 1:N_fit
        if ~ok_b(si), continue; end
        lc = ternary(is_kh_fit(si), CLR_KH, CLR_RR);
        plot(ax1b,[1 2],[b_det(si) b_prob(si)],'-o','Color',[lc 0.35],'MarkerSize',5,'MarkerFaceColor',lc);
    end
    for xi_pair = {1, 2}
        xi = xi_pair{1};
        vals = ternary(xi==1, b_det(ok_b), b_prob(ok_b));
        errorbar(ax1b, xi, mean(vals,'omitnan'), sem(vals), ...
            'ko','MarkerSize',10,'MarkerFaceColor','k','LineWidth',2,'CapSize',8,'HandleVisibility','off');
    end
    [~,p_b] = ttest(b_det(ok_b), b_prob(ok_b));
    text(ax1b,0.5,0.97,sprintf('p=%.3f',p_b),'Units','normalized', ...
        'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8,'BackgroundColor','w');
    set(ax1b,'XTick',[1 2],'XTickLabel',{'Deterministic','Probabilistic'});
    ylabel(ax1b,'\beta (inverse temperature)'); xlim(ax1b,[0.5 2.5]);
    add_significance_bracket(ax1b, 1, 2, max([b_det;b_prob],[],'omitnan')*1.05, p_b);

    % ── 1c: H scatter D vs P ────────────────────────────────────────────
    ax1c = subplot(1,4,3); hold(ax1c,'on');
    title(ax1c,'H_{det} vs H_{prob} per subject','FontSize',10);
    scatter(ax1c, H_det(ok_h & is_kh_fit), H_prob(ok_h & is_kh_fit), 60, CLR_KH,'filled','MarkerFaceAlpha',0.8,'DisplayName','Ox');
    scatter(ax1c, H_det(ok_h & is_rr_fit), H_prob(ok_h & is_rr_fit), 60, CLR_RR,'filled','MarkerFaceAlpha',0.8,'DisplayName','Nc');
    lims_h = [0, max([H_det(ok_h); H_prob(ok_h)])*1.1];
    plot(ax1c, lims_h, lims_h, 'k--', 'HandleVisibility','off');
    for si = 1:N_fit
        if ok_h(si), text(ax1c,H_det(si)+0.002,H_prob(si),subj_ids_fit{si},'FontSize',6,'Color',[0.5 0.5 0.5]); end
    end
    text(ax1c,0.05,0.97,'Points above diagonal: H_{prob} > H_{det}', ...
        'Units','normalized','VerticalAlignment','top','FontSize',7,'Color',[0.4 0.4 0.4]);
    xlabel(ax1c,'H_{det}'); ylabel(ax1c,'H_{prob}');
    xlim(ax1c,lims_h); ylim(ax1c,lims_h); axis(ax1c,'square');
    legend(ax1c,'Box','off','FontSize',8,'Location','southeast');

    % ── 1d: Noise sensitivity ΔH ────────────────────────────────────────
    ax1d = subplot(1,4,4); hold(ax1d,'on');
    title(ax1d,'\DeltaH = H_{prob} − H_{det}  (noise sensitivity)','FontSize',10);
    ok_ns = ~isnan(H_ns);
    histogram(ax1d, H_ns(ok_ns & is_kh_fit), 6, 'FaceColor',CLR_KH,'FaceAlpha',0.7,'EdgeColor','w','DisplayName','Ox');
    histogram(ax1d, H_ns(ok_ns & is_rr_fit), 6, 'FaceColor',CLR_RR,'FaceAlpha',0.7,'EdgeColor','w','DisplayName','Nc');
    xline(ax1d, 0, 'k--', 'LineWidth',1.5, 'HandleVisibility','off');
    xline(ax1d, mean(H_ns(ok_ns),'omitnan'), 'k-', 'LineWidth',2, 'HandleVisibility','off');
    [~,p_ns] = ttest(H_ns(ok_ns));
    text(ax1d,0.98,0.97,sprintf('Mean=%.3f\nt-test vs 0: p=%.3f',mean(H_ns(ok_ns),'omitnan'),p_ns), ...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',8,'BackgroundColor','w');
    xlabel(ax1d,'\DeltaH'); ylabel(ax1d,'Subjects');
    legend(ax1d,'Box','off','FontSize',8,'Location','best');

    saveas(fig1, fullfile(outpath,'fig_params_D_vs_P.png'));
    fprintf('Fig 1 saved.\n');
end

% =========================================================================
%% FIG 2 — PER-TRIAL MODEL LATENTS BY BLOCK TYPE
%  α_t, ω_t, |δ_t| aligned to reversal separately for D and P blocks
% =========================================================================
if ~isempty(results)

    half_win  = 20;
    x_ali     = -half_win : half_win;
    n_x       = numel(x_ali);
    latent_vars = {'alpha_trial','omega_trial',{'delta_trial','abs_delta'}};
    latent_lbls = {'\alpha_t (learning rate)', '\omega_t (CP probability)', '|\delta_t| (prediction error)'};

    % Matrices: rows = (subject×block), cols = time points, 4 conditions
    mats = struct();
    for lv = 1:3
        mats(lv).D  = NaN(0,n_x);
        mats(lv).P  = NaN(0,n_x);
    end

    for si = 1:numel(subj_ids_fit)
        sn = subj_ids_fit{si};
        r  = results.(sn);
        if ~isfield(all_trial_data, sn), continue; end
        td = all_trial_data.(sn).trial_data;
        if ~isfield(td,'revTrial'), continue; end

        for b = 1:numel(td.revTrial)
            rev = td.revTrial(b);
            if isnan(rev), continue; end

            is_det_b = isfield(r,'block_is_det') && b<=numel(r.block_is_det) && r.block_is_det(b);
            tag = ternary(is_det_b,'D','P');

            bm     = r.block_id == b;
            b_t    = r.trial_id(bm);

            latent_data = {r.alpha_trial(bm), r.omega_trial(bm), abs(r.delta_trial(bm))};

            for lv = 1:3
                row = NaN(1,n_x);
                for xi = 1:n_x
                    t_abs = round(rev) + x_ali(xi);
                    m_idx = find(b_t == t_abs, 1);
                    if ~isempty(m_idx), row(xi) = latent_data{lv}(m_idx); end
                end
                mats(lv).(tag)(end+1,:) = row;
            end
        end
    end

    fig2 = figure('Position',[50 50 1400 400]);
    sgtitle('Model latent variables aligned to reversal: D vs P blocks','FontSize',11);

    for lv = 1:3
        ax = subplot(1,3,lv); hold(ax,'on');
        title(ax, latent_lbls{lv},'FontSize',10);
        plot_ribbon_lc(ax, x_ali, mats(lv).D, CLR_D, '-',  sprintf('Det (n=%d)',size(mats(lv).D,1)));
        plot_ribbon_lc(ax, x_ali, mats(lv).P, CLR_P, '--', sprintf('Prob (n=%d)',size(mats(lv).P,1)));
        xline(ax,0,'k--','LineWidth',1.5,'HandleVisibility','off');
        xlabel(ax,'Trial relative to reversal');
        ylabel(ax, latent_lbls{lv},'Interpreter','tex');
        legend(ax,'Box','off','FontSize',9,'Location','best');
        xlim(ax,[-half_win half_win]);
    end
    saveas(fig2, fullfile(outpath,'fig_latents_D_vs_P_revaligned.png'));
    fprintf('Fig 2 saved.\n');
end

% =========================================================================
%% FIGS 3 & 4 — ACCURACY AND CONFIDENCE BY TRANSITION TYPE
%  Uses group_T (long-format table) with group_T.transition column
% =========================================================================

% Clean up transition column: remove NaN-string, keep only clean transitions
if ismember('transition', group_T.Properties.VariableNames)
    trans_col = string(group_T.transition);
    valid_trans = ismember(trans_col, TRANS_TYPES);
else
    warning('group_T has no transition column — add it: group_T.transition = strcat(string(group_T.prev_block_type),"→",string(group_T.block_type))');
    valid_trans = false(height(group_T),1);
end

% Block type mask
is_det_trial  = strcmp(string(group_T.block_type), 'D');
is_prob_trial = strcmp(string(group_T.block_type), 'P');

% Rev state (pre=0, post=1) already in table as group_T.rev_state
has_revstate = ismember('rev_state', group_T.Properties.VariableNames);

% ── Compute reversal-aligned matrices per transition type ─────────────────
% Strategy: for each subject, each block → determine transition type,
% then extract reversal-aligned accuracy and confidence vectors
% (these are already pre-computed in all_trial_data.aligned_correct etc.)

% We also pull aligned data directly from all_trial_data for clean alignment

trans_acc  = cell(1,numel(TRANS_TYPES));   % one matrix per transition
trans_conf = cell(1,numel(TRANS_TYPES));
trans_acc_pre  = cell(1,numel(TRANS_TYPES));  % pre-rev summary
trans_acc_post = cell(1,numel(TRANS_TYPES));
trans_conf_pre  = cell(1,numel(TRANS_TYPES));
trans_conf_post = cell(1,numel(TRANS_TYPES));
for ti = 1:numel(TRANS_TYPES)
    trans_acc{ti}  = NaN(0, alignedLen);
    trans_conf{ti} = NaN(0, alignedLen);
    trans_acc_pre{ti}  = [];
    trans_acc_post{ti} = [];
    trans_conf_pre{ti}  = [];
    trans_conf_post{ti} = [];
end

for si = 1:numel(subj_ids_td)
    sn = subj_ids_td{si};
    td = all_trial_data.(sn).trial_data;
    if ~isfield(td,'revTrial') || ~isfield(td,'aligned_correct'), continue; end

    [nB,~] = size(td.correct);

    % Determine block types from block_structure if available
    if isfield(td,'block_structure') && ~isempty(td.block_structure)
        bs = td.block_structure;
    else
        % Try to infer from pTrueFB
        bs = repmat('?',1,nB);
        if isfield(td,'pTrueFB')
            for b=1:nB
                pfb = td.pTrueFB(b, ~isnan(td.pTrueFB(b,:)));
                if ~isempty(pfb)
                    bs(b) = ternary(mean(pfb)>=0.99,'D','P');
                end
            end
        end
    end

    for b = 1:nB
        if b > numel(bs), continue; end
        curr_type = upper(char(bs(b)));
        if b == 1
            prev_type = 'none';
            tr_str    = '';
        else
            prev_type = upper(char(bs(b-1)));
            tr_str    = [prev_type '→' curr_type];
        end
        ti = find(strcmp(TRANS_TYPES, tr_str));
        if isempty(ti), continue; end

        % Aligned accuracy
        acc_row  = td.aligned_correct(b,:);
        conf_row = NaN(1,alignedLen);
        if isfield(td,'aligned_confidence')
            conf_row = td.aligned_confidence(b,:);
        end

        trans_acc{ti}(end+1,:)  = acc_row;
        trans_conf{ti}(end+1,:) = conf_row;

        % Summary: pre = positions 1:preN, post = preN+1:end
        trans_acc_pre{ti}(end+1)   = nanmean(acc_row(1:preN));
        trans_acc_post{ti}(end+1)  = nanmean(acc_row(preN+1:end));
        trans_conf_pre{ti}(end+1)  = nanmean(conf_row(1:preN));
        trans_conf_post{ti}(end+1) = nanmean(conf_row(preN+1:end));
    end
end

% ── Fig 3: Reversal-aligned ACCURACY by transition ────────────────────────
fig3 = figure('Position',[50 50 1600 500]);
sgtitle('Reversal-aligned accuracy by block transition type','FontSize',12);

for ti = 1:numel(TRANS_TYPES)
    ax = subplot(1, numel(TRANS_TYPES), ti);
    hold(ax,'on');
    title(ax, TRANS_TYPES{ti}, 'FontSize',11);

    mat = trans_acc{ti};
    if isempty(mat) || all(isnan(mat(:)))
        text(ax,0.5,0.5,'No data','Units','normalized','HorizontalAlignment','center');
        continue;
    end

    % Plot individual block-lines (light)
    for row = 1:size(mat,1)
        plot(ax, rel_ax, mat(row,:), 'Color', [TRANS_COLORS{ti} 0.12], 'LineWidth',0.5);
    end

    % Grand mean ribbon
    plot_ribbon_lc(ax, rel_ax, mat, TRANS_COLORS{ti}, '-', ...
        sprintf('Mean (n=%d)', size(mat,1)));

    xline(ax, 0, 'k--', 'LineWidth',1.5, 'HandleVisibility','off');
    yline(ax, 0.5, 'k:', 'HandleVisibility','off');
    xlabel(ax,'Trial relative to reversal');
    if ti==1, ylabel(ax,'P(correct)'); end
    xlim(ax,[-preN postN-1]); ylim(ax,[0 1]);

    % Add pre/post summary text
    mn_pre  = mean(trans_acc_pre{ti},'omitnan');
    mn_post = mean(trans_acc_post{ti},'omitnan');
    [~,p_t] = ttest(trans_acc_pre{ti}, trans_acc_post{ti});
    text(ax,0.02,0.97,sprintf('Pre=%.2f  Post=%.2f\nt-test: p=%.3f',mn_pre,mn_post,p_t), ...
        'Units','normalized','VerticalAlignment','top','FontSize',8,'BackgroundColor','w');

    % Shade pre and post regions
    ylims = [0 1];
    patch(ax,[-preN 0 0 -preN],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.85 0.85 0.85], ...
        'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
end
saveas(fig3, fullfile(outpath,'fig_accuracy_by_transition.png'));
fprintf('Fig 3 saved.\n');

% ── Fig 4: Reversal-aligned CONFIDENCE by transition ──────────────────────
fig4 = figure('Position',[50 50 1600 500]);
sgtitle('Reversal-aligned confidence by block transition type','FontSize',12);

for ti = 1:numel(TRANS_TYPES)
    ax = subplot(1, numel(TRANS_TYPES), ti);
    hold(ax,'on');
    title(ax, TRANS_TYPES{ti}, 'FontSize',11);

    mat = trans_conf{ti};
    if isempty(mat) || all(isnan(mat(:)))
        text(ax,0.5,0.5,'No data','Units','normalized','HorizontalAlignment','center');
        continue;
    end

    for row = 1:size(mat,1)
        plot(ax, rel_ax, mat(row,:), 'Color', [TRANS_COLORS{ti} 0.12], 'LineWidth',0.5);
    end

    plot_ribbon_lc(ax, rel_ax, mat, TRANS_COLORS{ti}, '-', ...
        sprintf('Mean (n=%d)', size(mat,1)));

    xline(ax, 0, 'k--', 'LineWidth',1.5, 'HandleVisibility','off');
    xlabel(ax,'Trial relative to reversal');
    if ti==1, ylabel(ax,'Confidence (1–10)'); end
    xlim(ax,[-preN postN-1]); ylim(ax,[1 10]);

    mn_pre  = mean(trans_conf_pre{ti},'omitnan');
    mn_post = mean(trans_conf_post{ti},'omitnan');
    [~,p_t] = ttest(trans_conf_pre{ti}, trans_conf_post{ti});
    text(ax,0.02,0.97,sprintf('Pre=%.1f  Post=%.1f\np=%.3f',mn_pre,mn_post,p_t), ...
        'Units','normalized','VerticalAlignment','top','FontSize',8,'BackgroundColor','w');

    patch(ax,[-preN 0 0 -preN],[1 1 10 10],[0.85 0.85 0.85], ...
        'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
end
saveas(fig4, fullfile(outpath,'fig_confidence_by_transition.png'));
fprintf('Fig 4 saved.\n');

% =========================================================================
%% FIG 5 — PRE vs POST ACCURACY & CONFIDENCE: 2×4 SUMMARY BAR PLOT
%  Rows: Accuracy / Confidence
%  Cols: D→D / D→P / P→D / P→P
% =========================================================================
fig5 = figure('Position',[50 50 1200 600]);
sgtitle('Pre vs post reversal: accuracy and confidence by transition','FontSize',12);

for ti = 1:numel(TRANS_TYPES)
    % Row 1: accuracy
    ax_a = subplot(2, numel(TRANS_TYPES), ti);
    hold(ax_a,'on');
    if ti==1, ylabel(ax_a,'P(correct)'); end
    title(ax_a, TRANS_TYPES{ti},'FontSize',11);

    pre_a  = trans_acc_pre{ti}(:);
    post_a = trans_acc_post{ti}(:);
    ok_a   = ~isnan(pre_a) & ~isnan(post_a);

    if sum(ok_a) > 1
        bar(ax_a, [1 2], [mean(pre_a(ok_a)) mean(post_a(ok_a))], ...
            'FaceColor',TRANS_COLORS{ti},'FaceAlpha',0.6,'EdgeColor','none');
        errorbar(ax_a, [1 2], [mean(pre_a(ok_a)) mean(post_a(ok_a))], ...
            [sem(pre_a(ok_a)) sem(post_a(ok_a))], ...
            'k.', 'LineWidth',1.5, 'CapSize',8);
        % Overlay individual points
        jitter_x = 0.15*(rand(sum(ok_a),1)-0.5);
        pre_ok  = pre_a(ok_a);
        post_ok = post_a(ok_a);
        scatter(ax_a, 1+jitter_x, pre_ok,  25, TRANS_COLORS{ti},'filled','MarkerFaceAlpha',0.6);
        scatter(ax_a, 2+jitter_x, post_ok, 25, TRANS_COLORS{ti},'filled','MarkerFaceAlpha',0.6);
        for i=1:numel(pre_ok)
            plot(ax_a,[1+jitter_x(i) 2+jitter_x(i)],[pre_ok(i) post_ok(i)],'-','Color',[0 0 0 0.15],'HandleVisibility','off');
        end
        [~,p_ta] = ttest(pre_a(ok_a), post_a(ok_a));
        add_significance_bracket(ax_a,1,2,0.95,p_ta);
    end
    set(ax_a,'XTick',[1 2],'XTickLabel',{'Pre','Post'});
    ylim(ax_a,[0 1]); yline(ax_a,0.5,'k:','HandleVisibility','off');

    % Row 2: confidence
    ax_c = subplot(2, numel(TRANS_TYPES), ti+numel(TRANS_TYPES));
    hold(ax_c,'on');
    if ti==1, ylabel(ax_c,'Confidence (1–10)'); end

    pre_c  = trans_conf_pre{ti}(:);
    post_c = trans_conf_post{ti}(:);
    ok_c   = ~isnan(pre_c) & ~isnan(post_c);

    if sum(ok_c) > 1
        bar(ax_c, [1 2], [mean(pre_c(ok_c)) mean(post_c(ok_c))], ...
            'FaceColor',TRANS_COLORS{ti},'FaceAlpha',0.6,'EdgeColor','none');
        errorbar(ax_c, [1 2], [mean(pre_c(ok_c)) mean(post_c(ok_c))], ...
            [sem(pre_c(ok_c)) sem(post_c(ok_c))], ...
            'k.', 'LineWidth',1.5, 'CapSize',8);
        jitter_x = 0.15*(rand(sum(ok_c),1)-0.5);
        pre_ok_c  = pre_c(ok_c);
        post_ok_c = post_c(ok_c);
        scatter(ax_c, 1+jitter_x, pre_ok_c,  25, TRANS_COLORS{ti},'filled','MarkerFaceAlpha',0.6);
        scatter(ax_c, 2+jitter_x, post_ok_c, 25, TRANS_COLORS{ti},'filled','MarkerFaceAlpha',0.6);
        for i=1:numel(pre_ok_c)
            plot(ax_c,[1+jitter_x(i) 2+jitter_x(i)],[pre_ok_c(i) post_ok_c(i)],'-','Color',[0 0 0 0.15],'HandleVisibility','off');
        end
        [~,p_tc] = ttest(pre_c(ok_c), post_c(ok_c));
        add_significance_bracket(ax_c,1,2,9.5,p_tc);
    end
    set(ax_c,'XTick',[1 2],'XTickLabel',{'Pre','Post'});
    ylim(ax_c,[1 10]);
end
saveas(fig5, fullfile(outpath,'fig_pre_post_by_transition.png'));
fprintf('Fig 5 saved.\n');

% =========================================================================
%% FIG 6 — REVERSAL-ALIGNED CONFIDENCE POOLED ACROSS ALL BLOCKS (D + P)
% =========================================================================

% Gather all aligned confidence rows regardless of block type
all_conf_rows = NaN(0, alignedLen);
all_acc_rows  = NaN(0, alignedLen);
row_block_type = {};   % 'D' or 'P' per row

for si = 1:numel(subj_ids_td)
    sn = subj_ids_td{si};
    td = all_trial_data.(sn).trial_data;
    if ~isfield(td,'aligned_correct'), continue; end

    [nB,~] = size(td.correct);

    % Determine block types
    if isfield(td,'block_structure') && ~isempty(td.block_structure)
        bs = td.block_structure;
    else
        bs = repmat('?',1,nB);
        if isfield(td,'pTrueFB')
            for b=1:nB
                pfb = td.pTrueFB(b, ~isnan(td.pTrueFB(b,:)));
                if ~isempty(pfb), bs(b) = ternary(mean(pfb)>=0.99,'D','P'); end
            end
        end
    end

    for b = 1:nB
        acc_row = td.aligned_correct(b,:);
        all_acc_rows(end+1,:) = acc_row;

        if isfield(td,'aligned_confidence')
            conf_row = td.aligned_confidence(b,:);
        else
            conf_row = NaN(1,alignedLen);
        end
        all_conf_rows(end+1,:) = conf_row;

        if b <= numel(bs)
            row_block_type{end+1} = upper(char(bs(b)));
        else
            row_block_type{end+1} = '?';
        end
    end
end

is_D_row = strcmp(row_block_type, 'D');
is_P_row = strcmp(row_block_type, 'P');

fig6 = figure('Position',[50 50 1200 700]);
sgtitle('Reversal-aligned performance and confidence — all subjects pooled','FontSize',12);

% ── 6a: Accuracy — pooled ─────────────────────────────────────────────────
ax6a = subplot(2,3,1); hold(ax6a,'on');
title(ax6a,'Accuracy: ALL blocks pooled','FontSize',10);
plot_ribbon_lc(ax6a, rel_ax, all_acc_rows, [0.3 0.3 0.3], '-', ...
    sprintf('All (n=%d)', size(all_acc_rows,1)));
xline(ax6a,0,'k--','LineWidth',1.5,'HandleVisibility','off');
yline(ax6a,0.5,'k:','HandleVisibility','off');
patch(ax6a,[-preN 0 0 -preN],[0 0 1 1],[0.85 0.85 0.85],'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
xlabel(ax6a,'Trial relative to reversal'); ylabel(ax6a,'P(correct)');
xlim(ax6a,[-preN postN-1]); ylim(ax6a,[0 1]);

% ── 6b: Accuracy — D vs P ─────────────────────────────────────────────────
ax6b = subplot(2,3,2); hold(ax6b,'on');
title(ax6b,'Accuracy: D vs P blocks','FontSize',10);
plot_ribbon_lc(ax6b, rel_ax, all_acc_rows(is_D_row,:), CLR_D, '-', ...
    sprintf('Det (n=%d)', sum(is_D_row)));
plot_ribbon_lc(ax6b, rel_ax, all_acc_rows(is_P_row,:), CLR_P, '--', ...
    sprintf('Prob (n=%d)', sum(is_P_row)));
xline(ax6b,0,'k--','LineWidth',1.5,'HandleVisibility','off');
yline(ax6b,0.5,'k:','HandleVisibility','off');
patch(ax6b,[-preN 0 0 -preN],[0 0 1 1],[0.85 0.85 0.85],'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
xlabel(ax6b,'Trial relative to reversal'); ylabel(ax6b,'P(correct)');
xlim(ax6b,[-preN postN-1]); ylim(ax6b,[0 1]);
legend(ax6b,'Box','off','FontSize',9);

% ── 6c: Accuracy — Cohort ─────────────────────────────────────────────────
% Build cohort tag per row
row_cohort = {};
for si = 1:numel(subj_ids_td)
    sn  = subj_ids_td{si};
    td  = all_trial_data.(sn).trial_data;
    if ~isfield(td,'aligned_correct'), continue; end
    nB = size(td.correct,1);
    for b = 1:nB
        row_cohort{end+1} = ternary(startsWith(sn,'Ox'),'Ox','Nc');
    end
end
is_kh_row = strcmp(row_cohort,'Ox');
is_rr_row = strcmp(row_cohort,'Nc');

ax6c = subplot(2,3,3); hold(ax6c,'on');
title(ax6c,'Accuracy: Ox vs Nc cohorts','FontSize',10);
plot_ribbon_lc(ax6c, rel_ax, all_acc_rows(is_kh_row,:), CLR_KH, '-', ...
    sprintf('Ox (n=%d)', sum(is_kh_row)));
plot_ribbon_lc(ax6c, rel_ax, all_acc_rows(is_rr_row,:), CLR_RR, '--', ...
    sprintf('Nc (n=%d)', sum(is_rr_row)));
xline(ax6c,0,'k--','LineWidth',1.5,'HandleVisibility','off');
yline(ax6c,0.5,'k:','HandleVisibility','off');
patch(ax6c,[-preN 0 0 -preN],[0 0 1 1],[0.85 0.85 0.85],'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
xlabel(ax6c,'Trial relative to reversal'); ylabel(ax6c,'P(correct)');
xlim(ax6c,[-preN postN-1]); ylim(ax6c,[0 1]);
legend(ax6c,'Box','off','FontSize',9);

% ── 6d: Confidence — pooled ───────────────────────────────────────────────
ax6d = subplot(2,3,4); hold(ax6d,'on');
title(ax6d,'Confidence: ALL blocks pooled','FontSize',10);
plot_ribbon_lc(ax6d, rel_ax, all_conf_rows, [0.3 0.3 0.3], '-', ...
    sprintf('All (n=%d)', size(all_conf_rows,1)));
xline(ax6d,0,'k--','LineWidth',1.5,'HandleVisibility','off');
patch(ax6d,[-preN 0 0 -preN],[1 1 10 10],[0.85 0.85 0.85],'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
xlabel(ax6d,'Trial relative to reversal'); ylabel(ax6d,'Confidence (1–10)');
xlim(ax6d,[-preN postN-1]); ylim(ax6d,[1 10]);
text(ax6d,0.5,0.02,'Confidence rated after decision, before feedback', ...
    'Units','normalized','HorizontalAlignment','center','FontSize',7,'Color',[0.5 0.5 0.5]);

% ── 6e: Confidence — D vs P ───────────────────────────────────────────────
ax6e = subplot(2,3,5); hold(ax6e,'on');
title(ax6e,'Confidence: D vs P blocks','FontSize',10);
plot_ribbon_lc(ax6e, rel_ax, all_conf_rows(is_D_row,:), CLR_D, '-', ...
    sprintf('Det (n=%d)', sum(is_D_row)));
plot_ribbon_lc(ax6e, rel_ax, all_conf_rows(is_P_row,:), CLR_P, '--', ...
    sprintf('Prob (n=%d)', sum(is_P_row)));
xline(ax6e,0,'k--','LineWidth',1.5,'HandleVisibility','off');
patch(ax6e,[-preN 0 0 -preN],[1 1 10 10],[0.85 0.85 0.85],'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
xlabel(ax6e,'Trial relative to reversal'); ylabel(ax6e,'Confidence (1–10)');
xlim(ax6e,[-preN postN-1]); ylim(ax6e,[1 10]);
legend(ax6e,'Box','off','FontSize',9);

% ── 6f: Confidence — Cohort ───────────────────────────────────────────────
ax6f = subplot(2,3,6); hold(ax6f,'on');
title(ax6f,'Confidence: Ox vs Nc cohorts','FontSize',10);
plot_ribbon_lc(ax6f, rel_ax, all_conf_rows(is_kh_row,:), CLR_KH, '-', ...
    sprintf('Ox (n=%d)', sum(is_kh_row)));
plot_ribbon_lc(ax6f, rel_ax, all_conf_rows(is_rr_row,:), CLR_RR, '--', ...
    sprintf('Nc (n=%d)', sum(is_rr_row)));
xline(ax6f,0,'k--','LineWidth',1.5,'HandleVisibility','off');
patch(ax6f,[-preN 0 0 -preN],[1 1 10 10],[0.85 0.85 0.85],'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
xlabel(ax6f,'Trial relative to reversal'); ylabel(ax6f,'Confidence (1–10)');
xlim(ax6f,[-preN postN-1]); ylim(ax6f,[1 10]);
legend(ax6f,'Box','off','FontSize',9);

saveas(fig6, fullfile(outpath,'fig_revaligned_pooled.png'));
fprintf('Fig 6 saved.\n');

% =========================================================================
%% FIG 7 — TRANSITION SUMMARY: SIDE-BY-SIDE COMPARISON
%  Compact summary: all 4 transitions in one panel per measure,
%  with pre/post error bars
% =========================================================================
fig7 = figure('Position',[50 50 900 700]);
sgtitle('Transition type × reversal state: accuracy and confidence','FontSize',12);

% ── Top: accuracy ribbons overlay ─────────────────────────────────────────
ax7a = subplot(2,1,1); hold(ax7a,'on');
title(ax7a,'Accuracy: all transition types','FontSize',11);
for ti = 1:numel(TRANS_TYPES)
    mat = trans_acc{ti};
    if isempty(mat), continue; end
    plot_ribbon_lc(ax7a, rel_ax, mat, TRANS_COLORS{ti}, ...
        ternary(mod(ti,2)==1,'-','--'), sprintf('%s (n=%d)',TRANS_TYPES{ti},size(mat,1)));
end
xline(ax7a,0,'k--','LineWidth',1.5,'HandleVisibility','off');
yline(ax7a,0.5,'k:','HandleVisibility','off');
patch(ax7a,[-preN 0 0 -preN],[0 0 1 1],[0.85 0.85 0.85],'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
xlabel(ax7a,'Trial relative to reversal'); ylabel(ax7a,'P(correct)');
xlim(ax7a,[-preN postN-1]); ylim(ax7a,[0 1]);
legend(ax7a,'Box','off','FontSize',9,'Location','southeast','NumColumns',2);

% ── Bottom: confidence ribbons overlay ────────────────────────────────────
ax7b = subplot(2,1,2); hold(ax7b,'on');
title(ax7b,'Confidence: all transition types','FontSize',11);
for ti = 1:numel(TRANS_TYPES)
    mat = trans_conf{ti};
    if isempty(mat), continue; end
    plot_ribbon_lc(ax7b, rel_ax, mat, TRANS_COLORS{ti}, ...
        ternary(mod(ti,2)==1,'-','--'), sprintf('%s (n=%d)',TRANS_TYPES{ti},size(mat,1)));
end
xline(ax7b,0,'k--','LineWidth',1.5,'HandleVisibility','off');
patch(ax7b,[-preN 0 0 -preN],[1 1 10 10],[0.85 0.85 0.85],'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
xlabel(ax7b,'Trial relative to reversal'); ylabel(ax7b,'Confidence (1–10)');
xlim(ax7b,[-preN postN-1]); ylim(ax7b,[1 10]);
legend(ax7b,'Box','off','FontSize',9,'Location','southeast','NumColumns',2);

saveas(fig7, fullfile(outpath,'fig_transition_overlay.png'));
fprintf('Fig 7 saved.\n');

% =========================================================================
%% PRINT GROUP STATISTICS
% =========================================================================
fprintf('\n=== GROUP STATISTICS: Accuracy and Confidence by Transition ===\n');
fprintf('%-8s  %-6s  %-6s  %-7s   %-6s  %-6s  %-7s   n\n', ...
    'Trans','AccPre','AccPost','p(acc)','ConfPre','ConfPost','p(conf)');
for ti = 1:numel(TRANS_TYPES)
    pa = trans_acc_pre{ti};  qa = trans_acc_post{ti};
    pc = trans_conf_pre{ti}; qc = trans_conf_post{ti};
    ok_a = ~isnan(pa) & ~isnan(qa);
    ok_c = ~isnan(pc) & ~isnan(qc);
    if sum(ok_a)>1, [~,p_a]=ttest(pa(ok_a),qa(ok_a)); else, p_a=NaN; end
    if sum(ok_c)>1, [~,p_c]=ttest(pc(ok_c),qc(ok_c)); else, p_c=NaN; end
    fprintf('%-8s  %-6.3f  %-6.3f  %-7.3f   %-6.2f  %-6.2f  %-7.3f   n=%d/%d\n', ...
        TRANS_TYPES{ti}, mean(pa,'omitnan'), mean(qa,'omitnan'), p_a, ...
        mean(pc,'omitnan'), mean(qc,'omitnan'), p_c, sum(ok_a), sum(ok_c));
end

fprintf('\n=== REVERSAL COST (pre − post accuracy) BY TRANSITION ===\n');
for ti = 1:numel(TRANS_TYPES)
    d = trans_acc_pre{ti} - trans_acc_post{ti};
    ok = ~isnan(d);
    if sum(ok)<2, continue; end
    [~,p,~,st] = ttest(d(ok));
    fprintf('%s:  ΔAcc=%.3f±%.3f  t(%d)=%.2f  p=%.3f\n', ...
        TRANS_TYPES{ti}, mean(d(ok)), std(d(ok)), st.df, st.tstat, p);
end

fprintf('\nAll figures saved to: %s\n', outpath);


% =========================================================================
%% LOCAL HELPER FUNCTIONS
% =========================================================================

function plot_ribbon_lc(ax, x, mat, clr, ls, lbl)
%PLOT_RIBBON_LC  Mean ± SEM ribbon with light smoothing.
if isempty(mat) || all(isnan(mat(:))), return; end
mn = movmean(mean(mat,1,'omitnan'), 3, 'omitnan');
se = std(mat,0,1,'omitnan') ./ sqrt(max(sum(~isnan(mat),1),1));
fill(ax,[x,fliplr(x)],[mn+se,fliplr(mn-se)],clr, ...
    'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
plot(ax,x,mn,'Color',clr,'LineWidth',2,'LineStyle',ls,'DisplayName',lbl);
end

function s = sem(x)
%SEM  Standard error of the mean (ignoring NaN).
x = x(~isnan(x));
s = std(x) / sqrt(max(numel(x),1));
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function add_significance_bracket(ax, x1, x2, y_top, p_val)
%ADD_SIGNIFICANCE_BRACKET  Draw a bracket with significance star above bars.
if isnan(p_val), return; end
if     p_val < 0.001, sig_str = '***';
elseif p_val < 0.01,  sig_str = '**';
elseif p_val < 0.05,  sig_str = '*';
else,                  sig_str = 'ns';
end
y_bar = y_top * 0.98;
line(ax, [x1 x1 x2 x2], [y_bar*0.97 y_bar y_bar y_bar*0.97], ...
    'Color','k','LineWidth',1,'HandleVisibility','off');
text(ax, mean([x1 x2]), y_bar*1.005, sig_str, ...
    'HorizontalAlignment','center','FontSize',10,'HandleVisibility','off');
end