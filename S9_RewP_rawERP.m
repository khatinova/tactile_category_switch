% =============================================================================
% P9_poster_ERP_behaviour_figures.m
%
% PURPOSE
% -------
% Generates four sets of publication/poster-quality figures:
%
%   FIGURE 1  — Reward Positivity (RewP) difference waves per stage × block type
%               (correct - incorrect ERP, showing the RewP component).
%               Derived from frn_rewp_stage_table (output of S3/S2).
%
%   FIGURE 2  — Grand-average ERP waveforms: correct, incorrect,
%               false-correct, false-incorrect (P blocks have 4 lines,
%               D blocks have 2), split across D and P panels.
%               Below each ERP: the corresponding RewP difference wave.
%
%   FIGURE 3  — First-10-trials accuracy: one axis for D blocks,
%               one for P blocks. Shows how quickly subjects settle
%               into each block type. Individual subjects as faint
%               lines; group mean ± SEM as bold overlay.
%
%   FIGURE 4  — Mock raw EEG traces (5 channels) for poster illustration.
%               Physiologically plausible synthetic data with alpha,
%               theta, beta components and realistic noise structure,
%               plus a vertical event marker.
%
% INPUTS (expected in workspace or loaded below)
% -----------------------------------------------
%   group_feature_table_combined.mat   — per-trial table (group_table_combined)
%   frn_rewp_by_stage_combined.mat     — per-stage FRN/RewP table
%   all_trial_data.mat                 — raw behavioural struct (for Fig 3)
%   behav_table.mat                    — long-format table (group_T, for Fig 3)
%
%   If the per-trial table contains prefrontal_waveform (cell column of
%   per-trial FCz waveforms), these are used to reconstruct grand averages
%   for Fig 2. Otherwise the script uses FRN_mean_amp as a proxy.
%
% OUTPUT FOLDER
% -------------
%   Results/EEG analysis/Figures/Poster_P9/
%
% STYLE
% -----
%   Ticks outside, no top/right spine, Arial, negative-up convention for
%   frontal ERPs, shaded SEM ribbons, stage colour palette consistent with
%   S6/S7.
% =============================================================================
% 
% clear; close all;
% 
% %% ── PATHS ───────────────────────────────────────────────────────────────────
% remote = 0;
% if remote
%     base_path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch';
% else
%     base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch';
% end
% 
% kh_results   = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis');
% feat_dir     = fullfile(kh_results, 'Outcome_feature_tables_v4_merged');
% data_dir     = fullfile(base_path, 'Salient mod switch KH', 'Data');
% outdir       = fullfile(kh_results, 'Figures', 'Poster_P9');
% if ~exist(outdir, 'dir'), mkdir(outdir); end
% 
% %% ── GLOBAL STYLE ────────────────────────────────────────────────────────────
% set(groot, 'defaultAxesTickDir',      'out');
% set(groot, 'defaultAxesBox',          'off');
% set(groot, 'defaultAxesFontSize',     11);
% set(groot, 'defaultAxesFontName',     'Arial');
% set(groot, 'defaultLineLineWidth',    1.8);
% set(groot, 'defaultTextFontName',     'Arial');
% 
% % Colour palette (matches S6/S7)
% CLR_D       = [0.15 0.45 0.70];   % blue   – deterministic
% CLR_P       = [0.80 0.30 0.10];   % orange – probabilistic
% CLR_CORRECT = [0.10 0.60 0.10];   % green
% CLR_INCORR  = [0.80 0.10 0.10];   % red
% CLR_FC      = [0.55 0.25 0.65];   % purple – false-correct (told wrong, were right)
% CLR_FI      = [0.90 0.60 0.10];   % amber  – false-incorrect (told right, were wrong)
% 
% STAGE_NAMES  = {'LN','LE','RN','RE'};
% STAGE_COLORS = [0.12 0.62 0.47;
%                 0.85 0.65 0.00;
%                 0.80 0.27 0.13;
%                 0.40 0.25 0.65];
% 
% ERP_XLIM    = [-200 800];   % ms
% FRN_WIN     = [250 350];    % ms – shaded region
% REWP_WIN    = [250 400];    % ms – positive shaded region
% BL_WIN      = [-200 0];     % ms
% 
% %% ── LOAD DATA ───────────────────────────────────────────────────────────────
% fprintf('Loading data...\n');
% 
% % Combined per-trial feature table
% gt = [];
% cand_gt = {fullfile(feat_dir, 'group_feature_table_combined.mat'), ...
%            fullfile(kh_results, 'Epoched_data', 'group_feature_table_combined.mat'), ...
%            fullfile(kh_results, 'Epoched_data_noisefiltering', 'group_feature_table_combined.mat')};
% for ci = 1:numel(cand_gt)
%     if exist(cand_gt{ci}, 'file')
%         S = load(cand_gt{ci});
%         if isfield(S, 'group_table_combined'), gt = S.group_table_combined;
%         elseif isfield(S, 'group_table'),       gt = S.group_table;
%         elseif isfield(S, 'all_trials_table'),  gt = S.all_trials_table;
%         end
%         if ~isempty(gt), fprintf('  Loaded per-trial table: %s\n', cand_gt{ci}); break; end
%     end
% end
% 
% % FRN/RewP per-stage table
% frn_tbl = [];
% cand_frn = {fullfile(feat_dir, 'frn_rewp_by_stage_combined.mat'), ...
%             fullfile(feat_dir, 'group_feature_table_KH_final.mat')};
% for ci = 1:numel(cand_frn)
%     if exist(cand_frn{ci}, 'file')
%         S = load(cand_frn{ci});
%         if isfield(S, 'frn_rewp_stage_table'),   frn_tbl = S.frn_rewp_stage_table;
%         elseif isfield(S, 'frn_tbl'),            frn_tbl = S.frn_tbl;
%         end
%         if ~isempty(frn_tbl), fprintf('  Loaded FRN/RewP stage table: %s\n', cand_frn{ci}); break; end
%     end
% end
% 
% % Behavioural table (for Fig 3)
% group_T = [];
% cand_beh = {fullfile(data_dir, 'behav_table_June2026_RL.mat'), ...
%             fullfile(data_dir, 'behav_table_June2026.mat'), ...
%             fullfile(data_dir, 'behav_table.mat')};
% for ci = 1:numel(cand_beh)
%     if exist(cand_beh{ci}, 'file')
%         S = load(cand_beh{ci});
%         if isfield(S, 'group_T'), group_T = S.group_T; end
%         if ~isempty(group_T), fprintf('  Loaded behav table: %s\n', cand_beh{ci}); break; end
%     end
% end
% 
% % all_trial_data (for Fig 3 first-10-trials)
% all_trial_data = [];
% cand_atd = {fullfile(data_dir, 'all_trial_data_June2026.mat'), ...
%             fullfile(data_dir, 'all_trial_data.mat')};
% for ci = 1:numel(cand_atd)
%     if exist(cand_atd{ci}, 'file')
%         S = load(cand_atd{ci});
%         if isfield(S, 'all_trial_data'), all_trial_data = S.all_trial_data; end
%         if ~isempty(all_trial_data), fprintf('  Loaded all_trial_data: %s\n', cand_atd{ci}); break; end
%     end
% end

%% ── TIME AXIS ───────────────────────────────────────────────────────────────
% Try to recover from the loaded table or the stage table
t_ax = [];
if ~isempty(frn_tbl) && ismember('t_ax', fieldnames(frn_tbl))
    t_ax = frn_tbl.t_ax;
end
% Search loaded workspaces for a saved t_ax variable
cand_tax = {fullfile(feat_dir, 'group_feature_table_KH_final.mat'), ...
            fullfile(feat_dir, 'group_feature_table_combined.mat'), ...
            fullfile(kh_results, 'Epoched_data', 'group_feature_table_combined.mat')};
if isempty(t_ax)
    for ci = 1:numel(cand_tax)
        if exist(cand_tax{ci}, 'file')
            S = load(cand_tax{ci}, 't_ax');
            if isfield(S, 't_ax') && ~isempty(S.t_ax)
                t_ax = S.t_ax; break;
            end
        end
    end
end
if isempty(t_ax)
    % Construct a plausible default (500 Hz, -200 to 800 ms)
    t_ax = -200 : 2 : 800;
    fprintf('  WARNING: t_ax not found – using default -200:2:800 ms\n');
end

% Time masks
in_erp    = t_ax >= ERP_XLIM(1) & t_ax <= ERP_XLIM(2);
t_plot    = t_ax(in_erp);
frn_mask  = t_ax >= FRN_WIN(1)  & t_ax <= FRN_WIN(2);
rewp_mask = t_ax >= REWP_WIN(1) & t_ax <= REWP_WIN(2);
bl_mask   = t_ax >= BL_WIN(1)   & t_ax <= BL_WIN(2);

n_t = sum(in_erp);

%% ═══════════════════════════════════════════════════════════════════════════
%%  FIGURE 1 — RewP Difference Waves per Stage × Block Type
%% ═══════════════════════════════════════════════════════════════════════════
fprintf('\nBuilding Figure 1: RewP difference waves per stage × block type...\n');

fig1 = figure('Position', [50 50 1400 650], 'Color', 'w');
sgtitle('Reward Positivity (RewP) — correct minus incorrect ERP by stage and feedback type', ...
    'FontSize', 13, 'FontWeight', 'bold');

if ~isempty(frn_tbl) && ismember('diff_wave', frn_tbl.Properties.VariableNames)
    % frn_rewp_stage_table structure:
    %   subj_id, block_type, stage, FRN_amp, RewP_amp, diff_wave (cell: correct-incorrect)
    %   diff_wave is (incorrect - correct); RewP_amp = mean over RewP window of (correct - incorrect)
    %   So RewP wave = -diff_wave

    bt_list = {'D', 'P'};

    for bi = 1:2
        bt  = bt_list{bi};
        clr = ternary_p9(bi==1, CLR_D, CLR_P);

        for si = 1:4
            sg  = STAGE_NAMES{si};
            ax  = subplot(2, 4, (bi-1)*4 + si);
            hold(ax, 'on');

            sel  = string(frn_tbl.block_type) == bt & ...
                   string(frn_tbl.stage)      == sg;
            dw_cell = frn_tbl.diff_wave(sel);
            dw_cell = dw_cell(~cellfun(@isempty, dw_cell));

            if ~isempty(dw_cell)
                % Stack diff waves; each cell is (incorrect - correct), so
                % multiply by -1 to get RewP polarity (correct - incorrect)
                M = cell2mat(cellfun(@(v) -v(:)', dw_cell, 'UniformOutput', false));

                if size(M, 2) ~= numel(t_ax)
                    % Interpolate to match t_ax length if sizes differ
                    t_orig = linspace(t_ax(1), t_ax(end), size(M,2));
                    M_interp = NaN(size(M,1), numel(t_ax));
                    for ri = 1:size(M,1)
                        M_interp(ri,:) = interp1(t_orig, M(ri,:), t_ax, 'linear', NaN);
                    end
                    M = M_interp;
                end

                mn = mean(M(:, in_erp), 1, 'omitnan');
                se = std( M(:, in_erp), 0, 1, 'omitnan') / sqrt(size(M,1));

                % Shade RewP window
                yl = [min(mn-se)-0.3, max(mn+se)+0.3];
                if any(isnan(yl)) || range(yl) < 0.01, yl = [-1 2]; end
                patch(ax, [REWP_WIN(1) REWP_WIN(2) REWP_WIN(2) REWP_WIN(1)], ...
                    [yl(1) yl(1) yl(2) yl(2)], [1 0.92 0.85], ...
                    'EdgeColor','none','FaceAlpha',0.5,'HandleVisibility','off');

                % Shade FRN window
                patch(ax, [FRN_WIN(1) FRN_WIN(2) FRN_WIN(2) FRN_WIN(1)], ...
                    [yl(1) yl(1) yl(2) yl(2)], [0.85 0.90 1.00], ...
                    'EdgeColor','none','FaceAlpha',0.5,'HandleVisibility','off');

                fill(ax, [t_plot, fliplr(t_plot)], [mn+se, fliplr(mn-se)], ...
                    clr, 'FaceAlpha',0.18,'EdgeColor','none','HandleVisibility','off');
                plot(ax, t_plot, mn, 'Color', clr, 'LineWidth', 2.2, ...
                    'DisplayName', sprintf('%s-%s (n=%d)', bt, sg, size(M,1)));

                % Mark RewP amplitude
                rewp_in_erp = rewp_mask(in_erp);
                if any(rewp_in_erp)
                    rewp_val = mean(mn(rewp_in_erp), 'omitnan');
                    text(ax, mean(REWP_WIN), rewp_val, sprintf(' %.2f uV', rewp_val), ...
                        'Color', clr, 'FontSize', 8, 'VerticalAlignment','middle');
                end
            else
                text(ax, 0.5, 0.5, 'No data', 'Units','normalized', ...
                    'HorizontalAlignment','center', 'Color', [0.6 0.6 0.6]);
            end

            xline(ax, 0, 'k:', 'LineWidth', 1, 'HandleVisibility','off');
            yline(ax, 0, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
            xlim(ax, ERP_XLIM);
            xlabel(ax, 'Time (ms)', 'FontSize', 9);
            if si == 1, ylabel(ax, 'RewP amplitude (uV)', 'FontSize', 9); end
            title(ax, sprintf('%s blocks — %s', bt, sg), 'FontSize', 10);
            legend(ax, 'Box','off','FontSize',7,'Location','northeast');
        end
    end

    % Add annotation explaining RewP
    annotation(fig1, 'textbox', [0.01 0.01 0.98 0.05], 'String', ...
        ['RewP (Reward Positivity) = correct minus incorrect ERP. ' ...
         'Blue shading = FRN window (250–350 ms). Orange shading = RewP window (250–400 ms). ' ...
         'Positive deflection in RewP window indicates reward positivity. ' ...
         'Each line is the grand average across subjects; shading = ±1 SEM.'], ...
        'FontSize', 8, 'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);
else
    % Fallback: plot RewP amplitude from the feature table scalar column
    fprintf('  diff_wave column not found – plotting RewP_amp scalar values per stage\n');

    has_rewp = ~isempty(frn_tbl) && ismember('RewP_amp', frn_tbl.Properties.VariableNames);
    has_frn  = ~isempty(frn_tbl) && ismember('FRN_amp',  frn_tbl.Properties.VariableNames);

    measures = {};
    if has_rewp, measures{end+1} = 'RewP_amp'; end
    if has_frn,  measures{end+1} = 'FRN_amp';  end
    if isempty(measures), measures = {'FRN_mean_amp'}; end

    for mi = 1:min(2, numel(measures))
        m = measures{mi};
        if ~ismember(m, frn_tbl.Properties.VariableNames), continue; end
        ax = subplot(2, 2, mi);
        hold(ax, 'on');
        title(ax, strrep(m,'_','\_'), 'FontSize', 11);

        for bt_i = 1:2
            bt = {'D','P'}; clr_bt = {CLR_D, CLR_P};
            for si = 1:4
                sel = string(frn_tbl.block_type)==bt{bt_i} & string(frn_tbl.stage)==STAGE_NAMES{si};
                vals = frn_tbl.(m)(sel);
                vals = vals(~isnan(vals));
                if isempty(vals), continue; end
                xpos = si + (bt_i-1.5)*0.25;
                bar(ax, xpos, mean(vals,'omitnan'), 0.2, 'FaceColor', clr_bt{bt_i}, ...
                    'EdgeColor','none','FaceAlpha',0.75);
                errorbar(ax, xpos, mean(vals,'omitnan'), ...
                    std(vals,'omitnan')/sqrt(numel(vals)), 'k.','CapSize',4,'LineWidth',1.2);
            end
        end
        plot(ax,NaN,NaN,'-','Color',CLR_D,'DisplayName','Deterministic');
        plot(ax,NaN,NaN,'-','Color',CLR_P,'DisplayName','Probabilistic');
        set(ax,'XTick',1:4,'XTickLabel',STAGE_NAMES);
        xlabel(ax,'Stage'); ylabel(ax,strrep(m,'_','\_'));
        legend(ax,'Box','off','Location','best');
        yline(ax, 0,'k--');
    end
end

save_poster_fig_p9(fig1, outdir, 'P9_Fig1_RewP_by_stage_blocktype');
fprintf('  Figure 1 saved.\n');

%% ═══════════════════════════════════════════════════════════════════════════
%%  FIGURE 2 — Grand-Average ERPs × Condition, with RewP below
%% ═══════════════════════════════════════════════════════════════════════════
fprintf('\nBuilding Figure 2: Grand-average ERP waveforms...\n');

% We need per-trial waveform data from the per-trial table.
% Source: prefrontal_waveform column (cell array of per-trial FCz traces)
% If absent, reconstruct from FRN_mean_amp as scalar fallback.

has_waves = ~isempty(gt) && ismember('prefrontal_waveform', gt.Properties.VariableNames);

if ~has_waves
    % Try alternative column names
    wave_cols = {'FCzCz_waveform','prefrontal_signal','FCz_waveform'};
    for wci = 1:numel(wave_cols)
        if ~isempty(gt) && ismember(wave_cols{wci}, gt.Properties.VariableNames)
            gt.prefrontal_waveform = gt.(wave_cols{wci});
            has_waves = true;
            fprintf('  Using waveform column: %s\n', wave_cols{wci});
            break;
        end
    end
end

fig2 = figure('Position', [50 50 1400 900], 'Color', 'w');
sgtitle('Grand-average ERPs: correct, incorrect, false feedback trials', ...
    'FontSize', 13, 'FontWeight', 'bold');

% We'll build 2 columns (D, P) × 2 rows (ERP panel, RewP panel)
% D column: 2 lines (correct, incorrect – true feedback only)
% P column: 4 lines (true-correct, true-incorrect, false-correct, false-incorrect)

bt_specs = {'D', 'P'};

for bi = 1:2
    bt = bt_specs{bi};
    is_P = strcmp(bt, 'P');

    %% ── Collect waveforms ────────────────────────────────────────────────
    if has_waves && ~isempty(gt)
        % Ensure string comparison works
        gt_bt = string(gt.block_type);
        gt_bt(gt_bt == "V") = "P";   % legacy code

        if isnumeric(gt.correct) || islogical(gt.correct)
            gt.correct_num = double(gt.correct);
        else
            cs = lower(string(gt.correct));
            gt.correct_num = double(cs == "1" | cs == "correct" | cs == "true");
        end

        if ~ismember('false_fb', gt.Properties.VariableNames)
            gt.false_fb = false(height(gt), 1);
        end
        gt.false_fb = logical(gt.false_fb);

        extract_erp = @(mask) extract_grand_avg_p9(gt, 'prefrontal_waveform', mask, t_ax, n_t, in_erp);

        mask_tc = gt_bt == bt & gt.correct_num == 1 & ~gt.false_fb;
        mask_ti = gt_bt == bt & gt.correct_num == 0 & ~gt.false_fb;

        [mn_tc, se_tc, n_tc] = extract_erp(mask_tc);
        [mn_ti, se_ti, n_ti] = extract_erp(mask_ti);

        if is_P
            % false_fb = 1 + correct_num = 1 means subject was correct,
            % shown incorrect feedback (false-negative in shown terms)
            mask_fc = gt_bt == bt & gt.correct_num == 1 & gt.false_fb;  % told wrong, were right
            mask_fi = gt_bt == bt & gt.correct_num == 0 & gt.false_fb;  % told right, were wrong
            [mn_fc, se_fc, n_fc] = extract_erp(mask_fc);
            [mn_fi, se_fi, n_fi] = extract_erp(mask_fi);
        end
    else
        % Scalar fallback using prefrontal_mean_amp
        fprintf('  Waveform column not available – using scalar mean amplitude for Fig 2\n');
        mn_tc = []; mn_ti = []; mn_fc = []; mn_fi = [];
        se_tc = []; se_ti = []; se_fc = []; se_fi = [];
        n_tc = 0; n_ti = 0; n_fc = 0; n_fi = 0;
    end

    %% ── ERP panel (top row) ──────────────────────────────────────────────
    ax_erp = subplot(2, 2, bi);
    hold(ax_erp, 'on');

    title(ax_erp, sprintf('%s blocks — grand-average ERP', ...
        ternary_p9(strcmp(bt,'D'),'Deterministic','Probabilistic')), 'FontSize', 11);

    % Shade windows
    plot_window_shade(ax_erp, FRN_WIN,  [0.85 0.90 1.00]);
    plot_window_shade(ax_erp, REWP_WIN, [1.00 0.92 0.85]);

    if ~isempty(mn_tc)
        plot_ribbon_p9(ax_erp, t_plot, mn_tc, se_tc, CLR_CORRECT, '-', ...
            sprintf('True correct (n=%d)', n_tc));
        plot_ribbon_p9(ax_erp, t_plot, mn_ti, se_ti, CLR_INCORR, '--', ...
            sprintf('True incorrect (n=%d)', n_ti));
        if is_P && ~isempty(mn_fc)
            plot_ribbon_p9(ax_erp, t_plot, mn_fc, se_fc, CLR_FC, '-.', ...
                sprintf('False correct* (n=%d)', n_fc));
            plot_ribbon_p9(ax_erp, t_plot, mn_fi, se_fi, CLR_FI, ':', ...
                sprintf('False incorrect† (n=%d)', n_fi));
        end
    else
        % Use frn_rewp_stage_table as fallback to at least show D vs P difference waves
        if ~isempty(frn_tbl) && ismember('diff_wave', frn_tbl.Properties.VariableNames)
            sel = string(frn_tbl.block_type) == bt;
            dw_cell = frn_tbl.diff_wave(sel);
            dw_cell = dw_cell(~cellfun(@isempty, dw_cell));
            if ~isempty(dw_cell)
                M = cell2mat(cellfun(@(v) v(:)', dw_cell, 'UniformOutput', false));
                if size(M,2) == numel(t_ax)
                    mn_dw = mean(M(:,in_erp),1,'omitnan');
                    se_dw = std( M(:,in_erp),0,1,'omitnan')/sqrt(size(M,1));
                    plot_ribbon_p9(ax_erp, t_plot, mn_dw, se_dw, ...
                        ternary_p9(strcmp(bt,'D'),CLR_D,CLR_P), '-', ...
                        sprintf('%s difference wave (n=%d)',bt,size(M,1)));
                end
            end
        else
            text(ax_erp, 0.5, 0.5, 'Waveform data not available', ...
                'Units','normalized','HorizontalAlignment','center','Color',[0.5 0.5 0.5]);
        end
    end

    xline(ax_erp, 0, 'k:', 'LineWidth',1, 'HandleVisibility','off');
    yline(ax_erp, 0, 'k--', 'LineWidth',0.8, 'HandleVisibility','off');
    set(ax_erp, 'YDir', 'reverse');   % EEG negative-up convention
    xlim(ax_erp, ERP_XLIM);
    xlabel(ax_erp, 'Time (ms)');
    ylabel(ax_erp, 'Amplitude (uV)  [negative up]');
    legend(ax_erp, 'Box','off','FontSize',8,'Location','northwest');

    % Annotation: FRN and RewP labels
    yl_erp = ylim(ax_erp);
    text(ax_erp, mean(FRN_WIN),  yl_erp(1)*0.85, 'FRN', ...
        'HorizontalAlignment','center','FontSize',8,'Color',[0.3 0.3 0.8]);
    text(ax_erp, mean(REWP_WIN), yl_erp(2)*0.85, 'RewP', ...
        'HorizontalAlignment','center','FontSize',8,'Color',[0.8 0.4 0.1]);

    if is_P
        annotation(fig2, 'textbox', [0.52 0.03 0.45 0.06], 'String', ...
            ['* False correct: participant was correct (right answer) but shown incorrect feedback (negative). ' ...
             '† False incorrect: participant was incorrect but shown positive feedback.'], ...
            'FontSize', 7, 'EdgeColor','none','BackgroundColor','none');
    end

    %% ── RewP difference wave panel (bottom row) ──────────────────────────
    ax_rewp = subplot(2, 2, 2 + bi);
    hold(ax_rewp, 'on');
    title(ax_rewp, sprintf('%s blocks — RewP wave (correct minus incorrect)', ...
        ternary_p9(strcmp(bt,'D'),'Deterministic','Probabilistic')), 'FontSize', 11);

    plot_window_shade(ax_rewp, FRN_WIN,  [0.85 0.90 1.00]);
    plot_window_shade(ax_rewp, REWP_WIN, [1.00 0.92 0.85]);

    if ~isempty(mn_tc) && ~isempty(mn_ti)
        rewp_true = mn_tc - mn_ti;
        se_comb   = sqrt((se_tc.^2 + se_ti.^2)/2);

        plot_ribbon_p9(ax_rewp, t_plot, rewp_true, se_comb, ...
            ternary_p9(strcmp(bt,'D'),CLR_D,CLR_P), '-', ...
            sprintf('True FB: correct−incorrect (n≈%d)', round((n_tc+n_ti)/2)));

        if is_P && ~isempty(mn_fc) && ~isempty(mn_fi)
            % False feedback RewP: what would be "correct-incorrect" in
            % terms of SHOWN outcome
            rewp_false = mn_fi - mn_fc;   % told-correct minus told-incorrect
            se_false   = sqrt((se_fi.^2 + se_fc.^2)/2);
            plot_ribbon_p9(ax_rewp, t_plot, rewp_false, se_false, CLR_FC, '--', ...
                sprintf('False FB: told-correct−told-incorrect (n≈%d)', round((n_fi+n_fc)/2)));
        end

        % Mark zero-crossing / RewP peak
        rewp_in = rewp_mask(in_erp);
        if any(rewp_in)
            [peak_val, peak_idx] = max(rewp_true(rewp_in));
            rewp_t = t_plot(rewp_in);
            if ~isnan(peak_val)
                plot(ax_rewp, rewp_t(peak_idx), peak_val, 'o', ...
                    'Color', ternary_p9(strcmp(bt,'D'),CLR_D,CLR_P), ...
                    'MarkerFaceColor', ternary_p9(strcmp(bt,'D'),CLR_D,CLR_P), ...
                    'MarkerSize', 7, 'HandleVisibility','off');
            end
        end
    elseif ~isempty(frn_tbl) && ismember('diff_wave', frn_tbl.Properties.VariableNames)
        sel = string(frn_tbl.block_type) == bt;
        dw_cell = frn_tbl.diff_wave(sel);
        dw_cell = dw_cell(~cellfun(@isempty, dw_cell));
        if ~isempty(dw_cell)
            M = cell2mat(cellfun(@(v) -v(:)', dw_cell, 'UniformOutput', false));
            if size(M,2) == numel(t_ax)
                mn_rw = mean(M(:,in_erp),1,'omitnan');
                se_rw = std( M(:,in_erp),0,1,'omitnan')/sqrt(size(M,1));
                plot_ribbon_p9(ax_rewp, t_plot, mn_rw, se_rw, ...
                    ternary_p9(strcmp(bt,'D'),CLR_D,CLR_P), '-', ...
                    sprintf('%s RewP (n=%d)', bt, size(M,1)));
            end
        end
    end

    xline(ax_rewp, 0, 'k:', 'LineWidth',1, 'HandleVisibility','off');
    yline(ax_rewp, 0, 'k--', 'LineWidth',0.8, 'HandleVisibility','off');
    xlim(ax_rewp, ERP_XLIM);
    xlabel(ax_rewp, 'Time (ms)');
    ylabel(ax_rewp, 'RewP amplitude (uV)  [correct − incorrect]');
    legend(ax_rewp, 'Box','off','FontSize',8,'Location','northwest');

    % Label windows on bottom plots
    yl_rw = ylim(ax_rewp);
    % text(ax_rewp, mean(FRN_WIN),  yl_rw(1) + 0.1*range(yl_rw), 'FRN window', ...
    %     'HorizontalAlignment','center','FontSize',7,'Color',[0.3 0.3 0.8]);
    text(ax_rewp, mean(REWP_WIN), yl_rw(2) - 0.1*range(yl_rw), 'RewP window', ...
        'HorizontalAlignment','center','FontSize',7,'Color',[0.8 0.4 0.1]);
end

save_poster_fig_p9(fig2, outdir, 'P9_Fig2_ERP_by_condition');
fprintf('  Figure 2 saved.\n');

%% ═══════════════════════════════════════════════════════════════════════════
%%  FIGURE 3 — First 10 Trials Performance: D vs P blocks
%% ═══════════════════════════════════════════════════════════════════════════
fprintf('\nBuilding Figure 3: First-10-trials accuracy...\n');

fig3 = figure('Position', [50 50 1200 620], 'Color', 'w');
sgtitle('First 10 trials accuracy per block: early-block performance by block type', ...
    'FontSize', 13, 'FontWeight', 'bold');

N_FIRST = 10;

% Build per-subject first-10-trials matrices from all_trial_data or group_T
first10_D = [];   % rows = blocks, length = N_FIRST
first10_P = [];
subj_idx_D = {};  % subject label per row
subj_idx_P = {};

if ~isempty(all_trial_data)
    subj_ids = fieldnames(all_trial_data);
    for si = 1:numel(subj_ids)
        sn  = subj_ids{si};
        atd = all_trial_data.(sn);
        if ~isfield(atd, 'trial_data'), continue; end
        td = atd.trial_data;
        if ~isfield(td, 'correct'), continue; end

        [nB, nT] = size(td.correct);

        % Infer block type
        bs = '';
        if isfield(td,'block_structure') && ~isempty(td.block_structure)
            bs = upper(char(td.block_structure));
            bs(bs == 'V') = 'P';
        elseif isfield(td,'trueFB')
            bs_arr = repmat('D', 1, nB);
            for b = 1:nB
                pfb = td.trueFB(b, ~isnan(td.trueFB(b,:)));
                if ~isempty(pfb) && mean(pfb) < 0.99, bs_arr(b) = 'P'; end
            end
            bs = bs_arr;
        end
        if isempty(bs), bs = repmat('D', 1, nB); end

        for b = 1:nB
            if nT < N_FIRST, continue; end
            row = td.correct(b, 1:N_FIRST);
            if all(isnan(row)), continue; end
            curr_type = char(bs(min(b, numel(bs))));
            if curr_type == 'D'
                first10_D(end+1,:) = double(row);
                subj_idx_D{end+1}  = sn;
            elseif curr_type == 'P'
                first10_P(end+1,:) = double(row);
                subj_idx_P{end+1}  = sn;
            end
        end
    end
elseif ~isempty(group_T)
    % Fallback: use group_T sorted by trial within block
    gt3 = group_T;
    if ~ismember('block_type_clean', gt3.Properties.VariableNames)
        gt3.block_type_clean = string(gt3.block_type);
        gt3.block_type_clean(gt3.block_type_clean == "V") = "P";
    end
    gt3.correct = double(gt3.correct);
    gt3.trial   = double(gt3.trial);
    gt3.block   = double(gt3.block);

    subjs3 = unique(string(gt3.subjID));
    for si = 1:numel(subjs3)
        sn = subjs3(si);
        Ts = gt3(string(gt3.subjID)==sn, :);
        blocks = unique(Ts.block)';
        for b = blocks
            Tb = Ts(Ts.block==b & Ts.trial <= N_FIRST, :);
            if height(Tb) < N_FIRST, continue; end
            [~, ord] = sort(Tb.trial);
            Tb = Tb(ord,:);
            row = Tb.correct(1:N_FIRST)';
            bt = char(Tb.block_type_clean(1));
            if bt == 'D'
                first10_D(end+1,:) = row;
                subj_idx_D{end+1} = char(sn);
            elseif bt == 'P'
                first10_P(end+1,:) = row;
                subj_idx_P{end+1} = char(sn);
            end
        end
    end
end

trial_ax = 1:N_FIRST;

for bt_i = 1:2
    ax3 = subplot(1, 2, bt_i);
    hold(ax3, 'on');

    if bt_i == 1
        mat = first10_D; clr = CLR_D; lbl = 'Deterministic blocks'; subj_lst = subj_idx_D;
    else
        mat = first10_P; clr = CLR_P; lbl = 'Probabilistic blocks'; subj_lst = subj_idx_P;
    end

    if isempty(mat)
        text(ax3, 0.5, 0.5, 'No data', 'Units','normalized', 'HorizontalAlignment','center');
        title(ax3, lbl);
        continue;
    end

    % Individual subject lines (one line per block occurrence)
    for ri = 1:size(mat,1)
        row = mat(ri,:);
        valid = ~isnan(row);
        if sum(valid) < 2, continue; end
        hl = plot(ax3, trial_ax(valid), row(valid), '-', ...
            'Color', clr, 'LineWidth', 0.6, 'HandleVisibility','off');
        try, hl.Color(4) = 0.18; catch, end
    end

    % Group mean ± SEM (per-block, not per-subject – valid as descriptive)
    gm = mean(mat, 1, 'omitnan');
    gs = std( mat, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(mat), 1));

    % Per-subject mean (across all their blocks of this type) for SEM
    subj_ids_3 = unique(string(subj_lst));
    subj_means = NaN(numel(subj_ids_3), N_FIRST);
    for si3 = 1:numel(subj_ids_3)
        sn3 = subj_ids_3(si3);
        row_mask = strcmp(string(subj_lst), sn3);
        if sum(row_mask) > 0
            subj_means(si3,:) = mean(mat(row_mask,:), 1, 'omitnan');
        end
    end
    gm_subj = mean(subj_means, 1, 'omitnan');
    gs_subj = std(subj_means, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(subj_means),1));

    fill(ax3, [trial_ax, fliplr(trial_ax)], ...
        [gm_subj + gs_subj, fliplr(gm_subj - gs_subj)], ...
        clr, 'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off');
    plot(ax3, trial_ax, gm_subj, '-', 'Color', clr, 'LineWidth', 2.8, ...
        'DisplayName', sprintf('Group mean ± SEM (n=%d blocks, n=%d subj)', ...
        size(mat,1), numel(subj_ids_3)));

    yline(ax3, 0.5, 'k--', 'LineWidth',1,'HandleVisibility','off');

    % Per-trial n
    n_per = sum(~isnan(mat),1);
    for ti = 1:N_FIRST
        text(ax3, ti, 0.02, sprintf('%d',n_per(ti)), ...
            'HorizontalAlignment','center','FontSize',7,'Color',[0.5 0.5 0.5]);
    end

    xlim(ax3, [0.5 N_FIRST+0.5]);
    ylim(ax3, [0 1.05]);
    set(ax3, 'XTick', 1:N_FIRST, 'TickDir','out');
    xlabel(ax3, 'Trial number within block');
    ylabel(ax3, 'P(correct)');
    title(ax3, sprintf('%s — first %d trials of each block', lbl, N_FIRST), 'FontSize',11);
    legend(ax3, 'Box','off','FontSize',9,'Location','southeast');

    text(ax3, 0.5, 1.0, 'Faint lines = individual blocks', ...
        'Units','normalized','HorizontalAlignment','center', ...
        'FontSize',8,'Color',[0.5 0.5 0.5],'VerticalAlignment','top');
end

annotation(fig3, 'textbox', [0.01 0.01 0.98 0.04], 'String', ...
    ['Figure 3: First 10 trials of each block. Each faint line is one block (pooled across subjects). ' ...
     'Bold = per-subject-averaged group mean ± 1 SEM. n per trial in grey. ' ...
     'In D blocks, performance should reach ceiling rapidly; in P blocks, noise may slow acquisition.'], ...
    'FontSize', 8, 'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

save_poster_fig_p9(fig3, outdir, 'P9_Fig3_First10trials_DvsP');
fprintf('  Figure 3 saved.\n');

%% ═══════════════════════════════════════════════════════════════════════════
%%  FIGURE 4 — Mock Raw EEG Traces for Poster Illustration
%% ═══════════════════════════════════════════════════════════════════════════
fprintf('\nBuilding Figure 4: Mock raw EEG traces...\n');

rng(42);   % reproducible

% --- Parameters ---
fs       = 500;             % Hz
dur_s    = 4.0;             % seconds of data to show
t_eeg    = (0:1/fs:dur_s - 1/fs);
n_samp   = numel(t_eeg);
event_t  = 2.0;             % outcome event at 2 s
event_s  = round(event_t * fs);

ch_names = {'FCz','Cz','Pz','C3','C4'};
n_ch     = numel(ch_names);

% Colour per channel (muted, readable on poster)
ch_clrs = [0.15 0.35 0.65;   % FCz – frontal
           0.30 0.55 0.75;   % Cz
           0.65 0.30 0.15;   % Pz – parietal
           0.20 0.60 0.30;   % C3 – left central
           0.60 0.25 0.50];  % C4 – right central

% --- Generate physiologically plausible synthetic EEG ---
% Each channel = alpha + theta + beta + 1/f noise + channel-specific noise
alpha_freq = 10 + randn(1, n_ch)*0.5;   % individual alpha peak 8-12 Hz
theta_freq = 6;
beta_freq  = 20;

eeg = zeros(n_ch, n_samp);

for ch = 1:n_ch
    % 1/f noise (pink noise approximation via AR filter)
    pink  = cumsum(randn(1, n_samp+200));
    pink  = pink(201:end);
    pink  = pink / std(pink) * 8;             % ~8 uV std

    % Alpha burst (stronger post-event)
    alpha_env = ones(1, n_samp) * 0.6;
    alpha_env(event_s+1:end) = 0.4 + 0.8 * exp(-(0:n_samp-event_s-1)/(0.5*fs));
    alpha  = alpha_env .* sin(2*pi*alpha_freq(ch)*t_eeg) * 12;

    % Theta component (frontal channels stronger)
    theta_amp = 6 * (ch <= 2) + 3 * (ch > 2);
    theta  = theta_amp .* sin(2*pi*theta_freq*t_eeg + rand*2*pi);

    % Post-event ERP-like deflection (N2/P3 morphology at frontal channels)
    erp = zeros(1, n_samp);
    if ch <= 3
        t_post = t_eeg - event_t;
        % N2: negative peak ~0.25 s post event
        erp = erp - 6 * exp(-(t_post - 0.25).^2 / (2*0.03^2)) .* (t_post >= 0);
        % P300: positive peak ~0.4 s post event (parietal stronger)
        p300_amp = 4 * (ch == 3) + 2 * (ch <= 2);
        erp = erp + p300_amp * exp(-(t_post - 0.40).^2 / (2*0.05^2)) .* (t_post >= 0);
    end

    % Beta spindles (motor cortex)
    beta_amp = 3 * (ch >= 4) + 1 * (ch < 4);
    beta = beta_amp .* sin(2*pi*beta_freq*t_eeg + rand*2*pi);

    % Eye blink artefact (slow, large, only at frontal, ~0.8 s in)
    blink = zeros(1, n_samp);
    if ch <= 2
        blink_t = 0.8;
        blink = 80 * exp(-((t_eeg - blink_t).^2)/(2*0.06^2));
    end

    % Muscle artefact (high-freq burst ~1.5 s, only C3/C4)
    muscle = zeros(1, n_samp);
    if ch >= 4
        mu_t = 1.5;
        mu_dur = 0.15;
        mu_mask = t_eeg >= mu_t & t_eeg <= mu_t + mu_dur;
        muscle(mu_mask) = randn(1, sum(mu_mask)) * 15;
    end

    eeg(ch,:) = pink + alpha + theta + erp + beta + blink + muscle;

    % Light bandpass smoothing for alpha (no edge effects needed for mock)
    eeg(ch,:) = eeg(ch,:) + randn(1, n_samp) * 2;
end

% --- Plot ---
fig4 = figure('Position', [50 50 1300 700], 'Color', 'w');

% Vertical offset for stacked display
offset_step = 120;   % uV spacing
offsets = (0:n_ch-1) * offset_step;

ax4 = axes(fig4, 'Position', [0.08 0.12 0.86 0.78]);
hold(ax4, 'on');

% Shaded epoch window (outcome-locked, e.g. –200 to +800 ms)
epoch_start = event_t - 0.2;
epoch_end   = event_t + 0.8;
yl_mock = [-offset_step/2, n_ch*offset_step + offset_step/2];
patch(ax4, [epoch_start epoch_end epoch_end epoch_start], ...
    [yl_mock(1) yl_mock(1) yl_mock(2) yl_mock(2)], ...
    [0.88 0.94 1.0], 'EdgeColor','none','FaceAlpha',0.6,'HandleVisibility','off');

for ch = 1:n_ch
    trace = eeg(ch,:) + offsets(ch);
    plot(ax4, t_eeg, trace, 'Color', ch_clrs(ch,:), ...
        'LineWidth', 1.0, 'DisplayName', ch_names{ch});
end

% Event marker
xline(ax4, event_t, 'k-', 'LineWidth', 2.5, 'HandleVisibility','off');
text(ax4, event_t + 0.04, yl_mock(2)*0.90, 'Outcome', ...
    'FontSize', 10, 'FontWeight','bold', 'Color', [0 0 0]);

% Annotations for visible artefacts
text(ax4, 0.82, offsets(2) + 55, '← Eye blink', ...
    'FontSize', 8, 'Color', [0.6 0.2 0.2], 'FontStyle','italic');
text(ax4, 1.52, offsets(4) + 55, '← Muscle', ...
    'FontSize', 8, 'Color', [0.5 0.5 0.5], 'FontStyle','italic');
text(ax4, event_t + 0.22, offsets(1) - 40, 'N2↓', ...
    'FontSize', 9, 'Color', ch_clrs(1,:), 'FontWeight','bold');
text(ax4, event_t + 0.40, offsets(3) + 30, 'P300↑', ...
    'FontSize', 9, 'Color', ch_clrs(3,:), 'FontWeight','bold');
text(ax4, event_t + 0.28, offsets(3) - 45, '← Epoch window (–200 to +800 ms)', ...
    'FontSize', 8, 'Color', [0.2 0.4 0.7], 'FontStyle','italic');

% Scale bar (50 uV × 0.5 s)
sb_t = [3.50 4.00]; sb_uv_centre = offsets(3);
plot(ax4, sb_t, [sb_uv_centre sb_uv_centre], 'k-', 'LineWidth', 2.5,'HandleVisibility','off');
plot(ax4, [sb_t(2) sb_t(2)], sb_uv_centre + [-25 25], 'k-', 'LineWidth',2.5,'HandleVisibility','off');
text(ax4, mean(sb_t), sb_uv_centre - 30, '500 ms', ...
    'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
text(ax4, sb_t(2)+0.04, sb_uv_centre, '50 uV', ...
    'HorizontalAlignment','left','FontSize',9,'FontWeight','bold');

% Y ticks = channel names
set(ax4, 'YTick', offsets, 'YTickLabel', ch_names, 'FontSize', 11, ...
    'TickDir','out', 'Box','off', 'YColor','k');
xlim(ax4, [t_eeg(1) t_eeg(end)]);
ylim(ax4, yl_mock);
xlabel(ax4, 'Time (s)', 'FontSize', 12);
title(ax4, 'Illustrative raw EEG traces — 5 channels, tactile category-switch outcome epoch', ...
    'FontSize', 12, 'FontWeight','bold');
legend(ax4, 'Box','off','FontSize',10,'Location','northeast');

% Coloured channel labels on right margin
for ch = 1:n_ch
    text(ax4, t_eeg(end)+0.03, offsets(ch), ch_names{ch}, ...
        'Color', ch_clrs(ch,:), 'FontSize',10,'FontWeight','bold', ...
        'HorizontalAlignment','left','Clipping','off');
end

annotation(fig4, 'textbox', [0.01 0.01 0.98 0.04], 'String', ...
    ['Figure 4: Illustrative raw EEG (synthetic, physiologically plausible). ' ...
     'Channels stacked with 120 uV offset. Blue shading = outcome epoch (−200 to +800 ms). ' ...
     'Alpha (∼10 Hz), theta (∼6 Hz), post-outcome N2 and P300 morphology visible. ' ...
     'Eye blink (FCz/Cz) and muscle artefact (C3/C4) shown for context.'], ...
    'FontSize', 8, 'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

save_poster_fig_p9(fig4, outdir, 'P9_Fig4_Mock_raw_EEG');
fprintf('  Figure 4 saved.\n');

%% ── DONE ────────────────────────────────────────────────────────────────────
fprintf('\n=== P9 COMPLETE ===\n');
fprintf('All figures saved to:\n  %s\n', outdir);
fprintf('Files:\n');
fprintf('  P9_Fig1_RewP_by_stage_blocktype.pdf / .png\n');
fprintf('  P9_Fig2_ERP_by_condition.pdf / .png\n');
fprintf('  P9_Fig3_First10trials_DvsP.pdf / .png\n');
fprintf('  P9_Fig4_Mock_raw_EEG.pdf / .png\n');


%% ═══════════════════════════════════════════════════════════════════════════
%%  LOCAL HELPER FUNCTIONS
%% ═══════════════════════════════════════════════════════════════════════════

function save_poster_fig_p9(fig, outdir, fname)
%SAVE_POSTER_FIG_P9  Export vector PDF + 300-dpi PNG.
    if ~exist(outdir,'dir'), mkdir(outdir); end
    try
        exportgraphics(fig, fullfile(outdir, [fname '.pdf']), 'ContentType','vector');
        exportgraphics(fig, fullfile(outdir, [fname '.png']), 'Resolution',300);
    catch ME_save
        saveas(fig, fullfile(outdir, [fname '.pdf']));
        warning('exportgraphics failed (%s); used saveas instead.', ME_save.message);
    end
end

% ─────────────────────────────────────────────────────────────────────────────
function [mn, se, n] = extract_grand_avg_p9(gt, wave_col, mask, t_ax, n_t, in_erp)
%EXTRACT_GRAND_AVG_P9  Per-subject-averaged grand mean from per-trial waveforms.
%  Returns mn and se (each 1 × n_t), and n = number of subjects.
    mn = NaN(1, n_t);
    se = NaN(1, n_t);
    n  = 0;

    if ~ismember(wave_col, gt.Properties.VariableNames), return; end

    rows = find(mask);
    if isempty(rows), return; end

    subjs = unique(string(gt.subj_id(rows)));
    if isempty(subjs) && ismember('subjID', gt.Properties.VariableNames)
        subjs = unique(string(gt.subjID(rows)));
    end
    if isempty(subjs), return; end

    subj_avg = NaN(numel(subjs), n_t);
    for si = 1:numel(subjs)
        sn = subjs(si);
        if ismember('subj_id', gt.Properties.VariableNames)
            sm = mask & string(gt.subj_id) == sn;
        else
            sm = mask & string(gt.subjID)  == sn;
        end

        waves = gt.(wave_col)(sm);
        waves = waves(~cellfun(@isempty, waves));
        if isempty(waves), continue; end

        % Stack and baseline-correct
        M_raw = cell2mat(cellfun(@(v) v(:)', waves, 'UniformOutput', false));

        % Handle length mismatch
        if size(M_raw,2) ~= numel(t_ax)
            t_orig = linspace(t_ax(1), t_ax(end), size(M_raw,2));
            M_new  = NaN(size(M_raw,1), numel(t_ax));
            for ri = 1:size(M_raw,1)
                M_new(ri,:) = interp1(t_orig, M_raw(ri,:), t_ax, 'linear', NaN);
            end
            M_raw = M_new;
        end

        % Baseline correction (−200 to 0 ms)
        bl = t_ax >= -200 & t_ax <= 0;
        if any(bl)
            bl_mean = mean(M_raw(:, bl), 2, 'omitnan');
            M_raw   = M_raw - bl_mean;
        end

        subj_avg(si,:) = mean(M_raw(:, in_erp), 1, 'omitnan');
    end

    ok = ~all(isnan(subj_avg), 2);
    n  = sum(ok);
    if n == 0, return; end
    mn = mean(subj_avg(ok,:), 1, 'omitnan');
    se = std( subj_avg(ok,:), 0, 1, 'omitnan') ./ sqrt(n);
end

% ─────────────────────────────────────────────────────────────────────────────
function plot_ribbon_p9(ax, t, mn, se, clr, ls, lbl)
%PLOT_RIBBON_P9  Plot mean ± SEM ribbon on ax.
    ok = ~isnan(mn);
    if ~any(ok), return; end
    t_ok  = t(ok); mn_ok = mn(ok); se_ok = se(ok);
    fill(ax, [t_ok, fliplr(t_ok)], [mn_ok+se_ok, fliplr(mn_ok-se_ok)], ...
        clr, 'FaceAlpha',0.18,'EdgeColor','none','HandleVisibility','off');
    plot(ax, t_ok, mn_ok, 'Color',clr,'LineWidth',2.0,'LineStyle',ls, ...
        'DisplayName',lbl);
end

% ─────────────────────────────────────────────────────────────────────────────
function plot_window_shade(ax, win_ms, rgb)
%PLOT_WINDOW_SHADE  Add a shaded background patch for a time window.
    yl = ylim(ax);
    patch(ax, [win_ms(1) win_ms(2) win_ms(2) win_ms(1)], ...
        [yl(1) yl(1) yl(2) yl(2)], rgb, ...
        'EdgeColor','none','FaceAlpha',0.45,'HandleVisibility','off');
end

% ─────────────────────────────────────────────────────────────────────────────
function out = ternary_p9(cond, a, b)
    if cond, out = a; else, out = b; end
end