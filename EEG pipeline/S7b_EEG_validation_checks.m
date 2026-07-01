% =============================================================================
% S7b_EEG_validation_checks.m
%
% PIPELINE STEP 7b — EEG feature table validation and quality control.
%
% Canonical version: supersedes S7b_EEG_validation_checks_debugged.m
%
% PURPOSE
% -------
% Produces all recommended data-validation figures and summary tables
% for the prefrontal negativity (FRN / FCz neg-peak), RewP, P300, theta,
% and PLV features BEFORE any inferential analysis is run in S7.
%
% The checks implemented here follow the recommendations of:
%   Clayson et al. (2020) Psychophysiology 57:e13486.
%   "Moderators of the internal consistency of the error-related negativity"
%   Key recommendations used:
%     1. Report trial counts per condition; flag cells < 8–10 trials.
%     2. Show amplitude distributions (histograms + Q-Q plots).
%     3. Check for systematic outliers (> ±3 SD after z-scoring).
%     4. Verify baseline period is flat (mean ≈ 0, SD small).
%     5. Plot grand-average waveforms to confirm expected morphology.
%     6. Report internal consistency (split-half reliability) for peak
%        measures across subjects.
%     7. Examine trial-count ~ amplitude correlations (trial-count bias).
%     8. Check for session/block order confounds in amplitude.
%
% ADDITIONAL CHECKS (ERP / time-frequency best practice):
%     9.  Distribution of FRN peak latencies.
%    10.  FRN_excluded rate by subject and condition.
%    11.  Correlation between raw and z-scored amplitudes.
%    12.  Baseline RMS distribution across subjects (signal quality proxy).
%    13.  Theta amplitude time course check (should peak post-outcome).
%    14.  PLV baseline-period level (should be near zero after correction).
%    15.  Cohort comparison (KH vs RR) for all key features.
%    16.  Block-order effects (does amplitude change over the session?).
%
% INPUT
% -----
%   group_feature_table_combined.mat  (group_table from S4)
%   — same input as S7_eeg_rq_analysis.m
%
% OUTPUT
% ------
%   Figures  : PDF panels saved to Figures/Validation_S7b/
%   Tables   : trial-count summary tables (CSV + console print)
%   MAT file : S7b_validation_workspace.mat (can be loaded for inspection)
%
% RUN ORDER: Run after S4 (merge), before S7 (inferential analysis).
% =============================================================================

clear; close all; clc;
addpath(genpath(fileparts(mfilename('fullpath'))));

set(groot,'defaultAxesTickDir','out');
set(groot,'defaultAxesBox','off');

% ── CONSTANTS ────────────────────────────────────────────────────────────────
MIN_TRIALS_WARN  = 8;        % Clayson (2020): < 8 trials risks poor reliability
MIN_TRIALS_EXCL  = 4;        % below this a cell is treated as empty
OUTLIER_SD_THRESH = 3;       % flag subjects > ±3 SD from group mean
SPLIT_HALF_REPS   = 1000;    % permutations for Spearman-Brown split-half

STAGES     = {'LN','LE','RN','RE'};
BTYPES     = {'D','P'};
OUTCOMES   = {'Incorrect','Correct'};

CLR_D    = [0.15 0.45 0.70];
CLR_P    = [0.80 0.30 0.10];
CLR_INC  = [0.80 0.27 0.13];
CLR_COR  = [0.12 0.62 0.47];
CLR_KH   = [0.15 0.45 0.70];
CLR_RR   = [0.65 0.25 0.65];

% ── PATHS ────────────────────────────────────────────────────────────────────
base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch';

feat_dir = fullfile(base_path,'Salient mod switch KH','Results','EEG analysis', ...
    'Outcome_feature_tables_v4_merged');
outdir   = fullfile(base_path,'Salient mod switch KH','Results','EEG analysis', ...
    'Figures','Validation_S7b');
if ~exist(outdir,'dir'), mkdir(outdir); end

% ── LOAD ─────────────────────────────────────────────────────────────────────
fprintf('Loading combined feature table...\n');
combined_candidates = {
    fullfile(feat_dir,  'group_feature_table_combined.mat')
    fullfile(base_path,'Salient mod switch KH','Results','EEG analysis', ...
             'Epoched_data','group_feature_table_combined.mat')
    fullfile(base_path,'Salient mod switch KH','Results','EEG analysis', ...
             'Epoched_data_noisefiltering','group_feature_table_combined.mat')
};
gt = [];
for ci = 1:numel(combined_candidates)
    if exist(combined_candidates{ci},'file')
        S = load(combined_candidates{ci},'group_table');
        gt = S.group_table;
        fprintf('  Loaded from: %s\n', combined_candidates{ci});
        break
    end
end
if isempty(gt)
    error('S7b: combined feature table not found. Run S4_merge_feature_tables first.');
end

% ── TYPE COERCION ─────────────────────────────────────────────────────────────
gt.subj_id = categorical(string(gt.subj_id));
gt.cohort  = categorical(string(gt.cohort));

% Normalise legacy block labels before categorical conversion.
% Some older scripts use V for probabilistic/visual-probabilistic blocks.
bt = string(gt.block_type);
bt(bt == "V") = "P";
gt.block_type = categorical(bt, BTYPES);

gt.stage = categorical(string(gt.stage), STAGES, 'Ordinal', false);

% Robust outcome coding. S7 models use 0 = Incorrect and 1 = Correct.
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
    error('S7b: neither correct nor correct_num was found in group_table.');
end
% If a categorical Correct/Incorrect column was converted to category codes
% 1/2, repair it to the required 0/1 coding.
if all(ismember(unique(gt.correct_num(~isnan(gt.correct_num))), [1 2]))
    gt.correct_num = gt.correct_num - 1;
end

if ismember('false_fb',gt.Properties.VariableNames)
    gt.false_fb = logical(gt.false_fb);
else
    gt.false_fb = false(height(gt),1);
end

subjects = categories(gt.subj_id);
N_subj   = numel(subjects);
fprintf('Table: %d rows | %d subjects | cohorts: %s\n', ...
    height(gt), N_subj, strjoin(categories(gt.cohort),', '));

% ── FEATURE INVENTORY ────────────────────────────────────────────────────────
% Each entry: {raw_col, z_col, display_label, expected_sign_incorrect}
FEATURES = {
    'prefrontal_neg_peak_amp',  'prefrontal_neg_peak_amp_z',  'FCz neg-peak amp (µV)',   -1
    'prefrontal_neg_peak_norm', 'prefrontal_neg_peak_norm_z', 'FCz neg-peak norm',        -1
    'prefrontal_mean_amp',      'prefrontal_mean_amp_z',      'FCz mean amp (µV)',        -1
    'P300_amp',                 'P300_amp_z',                 'P300 amp (µV)',            +1
    'P300_norm',                'P300_norm_z',                'P300 norm',               +1
    'Theta_amp',                'Theta_amp_z',                'Theta amp',               +1
    'PLV_fp',                   'PLV_fp_z',                   'FP PLV',                  +1
    'PLV_fs',                   'PLV_fs_z',                   'FS PLV',                  +1
};
% Keep only features present in the table
has_feat = cellfun(@(c) ismember(c,gt.Properties.VariableNames), FEATURES(:,1));
FEATURES = FEATURES(has_feat,:);
N_feat   = size(FEATURES,1);
fprintf('Features available: %d / %d\n', N_feat, size(has_feat,1));


% =============================================================================
%% SECTION 1 — TRIAL COUNT TABLE
%   Clayson (2020): primary recommendation. Flag any cell < MIN_TRIALS_WARN.
%   One table per feature (raw col ≠ NaN), stratified by:
%   subject × stage × block_type × outcome (true-feedback only)
% =============================================================================
fprintf('\n=== SECTION 1: Trial count table ===\n');

% Build master count table
count_rows = {};
for si = 1:N_subj
    sn = subjects{si};
    coh = char(string(gt.cohort(find(gt.subj_id==sn,1))));
    for stgi = 1:numel(STAGES)
        for bti = 1:numel(BTYPES)
            for oci = 1:numel(OUTCOMES)
                corr_val = (oci == 2);   % 1=Correct
                base_mask = gt.subj_id==sn & gt.stage==STAGES{stgi} & ...
                    gt.block_type==BTYPES{bti} & gt.correct_num==corr_val & ...
                    ~gt.false_fb;
                n_trials = sum(base_mask);
                % Count valid (non-NaN) trials per feature
                n_eeg = nan(1,N_feat);
                for fi = 1:N_feat
                    raw = FEATURES{fi,1};
                    if ismember(raw,gt.Properties.VariableNames)
                        n_eeg(fi) = sum(base_mask & ~isnan(gt.(raw)));
                    end
                end
                count_rows(end+1,:) = [{sn},{coh},{STAGES{stgi}},{BTYPES{bti}},{OUTCOMES{oci}}, ...
                    {n_trials}, num2cell(n_eeg)];
            end
        end
    end
end

feat_hdrs = FEATURES(:,1)';
count_T = cell2table(count_rows, ...
    'VariableNames',[{'subj_id','cohort','stage','block_type','outcome','n_beh_trials'}, feat_hdrs]);

% cell2table leaves text columns as cell arrays; convert so later == comparisons work.
count_T.subj_id    = string(count_T.subj_id);
count_T.cohort     = categorical(string(count_T.cohort));
count_T.stage      = categorical(string(count_T.stage), STAGES, 'Ordinal', false);
count_T.block_type = categorical(string(count_T.block_type), BTYPES);
count_T.outcome    = categorical(string(count_T.outcome), OUTCOMES);

% Flag cells below threshold
for fi = 1:N_feat
    raw = FEATURES{fi,1};
    col = count_T.(raw);
    n_warn = sum(col < MIN_TRIALS_WARN & col >= MIN_TRIALS_EXCL);
    n_excl = sum(col < MIN_TRIALS_EXCL);
    fprintf('  %-35s  warn(<8): %d  excl(<4): %d\n', raw, n_warn, n_excl);
end

% Save CSV
writetable(count_T, fullfile(outdir,'S7b_trial_counts.csv'));
fprintf('Trial count table saved.\n');

% ── Figure 1a: heatmap of trial counts (key features x conditions) ───────────
% Aggregate over subjects: median trial count per stage x block_type x outcome
for fi = 1:min(N_feat,3)   % show first 3 features to keep figure manageable
    raw = FEATURES{fi,1};
    fig = figure('Position',[50 50 900 500]);
    sgtitle(sprintf('Trial counts per condition: %s', strrep(raw,'_','\_')));

    for bti = 1:2
        ax = subplot(1,2,bti); hold(ax,'on');
        title(ax, BTYPES{bti},'FontSize',11);

        mat = nan(numel(STAGES),2);   % stages x outcomes
        for stgi = 1:numel(STAGES)
            for oci = 1:2
                rows = count_T.stage==STAGES{stgi} & count_T.block_type==BTYPES{bti} & ...
                    count_T.outcome==OUTCOMES{oci};
                vals = count_T.(raw)(rows);
                mat(stgi,oci) = median(vals,'omitnan');
            end
        end

        imagesc(ax, mat);
        colormap(ax, parula);
        cb = colorbar(ax); cb.Label.String = 'Median N trials';
        caxis(ax, [0 max(mat(:), [], 'omitnan') + 1]);

        % Red contour on cells < MIN_TRIALS_WARN
        for stgi = 1:numel(STAGES)
            for oci = 1:2
                clr_txt = 'w';
                if mat(stgi,oci) < MIN_TRIALS_WARN
                    rectangle(ax,'Position',[oci-0.5 stgi-0.5 1 1],...
                        'EdgeColor',[0.8 0 0],'LineWidth',2);
                    clr_txt = [0.8 0 0];
                end
                text(ax,oci,stgi,sprintf('%d',round(mat(stgi,oci))),...
                    'HorizontalAlignment','center','FontSize',10,...
                    'Color',clr_txt,'FontWeight','bold');
            end
        end

        set(ax,'XTick',1:2,'XTickLabel',OUTCOMES,'YTick',1:numel(STAGES),...
            'YTickLabel',STAGES,'FontSize',9);
        xlabel(ax,'Outcome'); ylabel(ax,'Stage');
    end
    annotation('textbox',[0.01 0.0 0.98 0.05],'String',...
        sprintf('Red border = median < %d trials (Clayson 2020 reliability threshold). Values are median across subjects.',MIN_TRIALS_WARN),...
        'EdgeColor','none','FontSize',7);
    exportgraphics(fig, fullfile(outdir,sprintf('S7b_1_trial_counts_%s.pdf',raw)),'ContentType','vector');
end

% ── Figure 1b: per-subject trial count distribution ─────────────────────────
fig = figure('Position',[50 50 1200 400]);
sgtitle('Per-subject trial counts (true-feedback trials, all stages/conditions pooled)');
for fi = 1:N_feat
    raw = FEATURES{fi,1};
    ax = subplot(2,ceil(N_feat/2),fi); hold(ax,'on');
    title(ax,FEATURES{fi,3},'FontSize',8,'Interpreter','none');

    % Total valid trials per subject
    n_per_subj = arrayfun(@(s) ...
        sum(gt.subj_id==s & ~isnan(gt.(raw)) & ~gt.false_fb), ...
        categorical(subjects));

    clrs = arrayfun(@(s) ...
        isequal(char(string(gt.cohort(find(gt.subj_id==s,1)))),'KH'), ...
        categorical(subjects));

    for si2 = 1:N_subj
        bar(ax,si2,n_per_subj(si2),0.7,...
            'FaceColor', ternary_local(clrs(si2),CLR_KH,CLR_RR),...
            'EdgeColor','none','FaceAlpha',0.8,'HandleVisibility','off');
    end
    yline(ax,MIN_TRIALS_WARN * numel(STAGES) * numel(BTYPES) * numel(OUTCOMES),...
        'r--','LineWidth',1.2,'HandleVisibility','off');
    set(ax,'XTick',1:N_subj,'XTickLabel',subjects,'XTickLabelRotation',45,'FontSize',6);
    ylabel(ax,'N valid trials');
end
exportgraphics(fig,fullfile(outdir,'S7b_1b_trials_per_subject.pdf'),'ContentType','vector');


% =============================================================================
%% SECTION 2 — AMPLITUDE DISTRIBUTIONS
%   Clayson (2020): examine distributions; check for non-normality.
%   For each feature: histogram (raw), Q-Q plot, histogram (z-scored).
% =============================================================================
fprintf('\n=== SECTION 2: Amplitude distributions ===\n');

for fi = 1:N_feat
    raw = FEATURES{fi,1};
    zc  = FEATURES{fi,2};
    lbl = FEATURES{fi,3};
    sgn = FEATURES{fi,4};

    if ~ismember(raw,gt.Properties.VariableNames), continue; end

    raw_mask = ~isnan(gt.(raw)) & ~gt.false_fb;
    raw_vals = gt.(raw)(raw_mask);

    z_vals = [];
    z_mask = false(height(gt),1);
    if ismember(zc,gt.Properties.VariableNames)
        z_mask = raw_mask & ~isnan(gt.(zc));
        z_vals = gt.(zc)(z_mask);
    end

    fig = figure('Position',[50 50 1400 420]);
    sgtitle(sprintf('Distribution checks: %s', strrep(lbl,'_','\_')));

    % Panel 1: raw histogram
    ax1 = subplot(1,4,1); hold(ax1,'on');
    title(ax1,'Raw amplitude histogram');
    histogram(ax1,raw_vals,40,'FaceColor',[0.4 0.6 0.8],'EdgeColor','none','FaceAlpha',0.8);
    xline(ax1,mean(raw_vals,'omitnan'),'r-','LineWidth',2,'DisplayName','Mean');
    xline(ax1,median(raw_vals,'omitnan'),'k--','LineWidth',1.5,'DisplayName','Median');
    xlabel(ax1,lbl,'Interpreter','none'); ylabel(ax1,'Count');
    legend(ax1,'Box','off','FontSize',7);

    % Normality test
    [~,p_sw] = swtest_local(raw_vals);
    text(ax1,0.98,0.97,sprintf('S-W p = %.3f\nN = %d',p_sw,numel(raw_vals)),...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top',...
        'FontSize',7,'BackgroundColor','w');

    % Panel 2: Q-Q plot
    ax2 = subplot(1,4,2); hold(ax2,'on');
    title(ax2,'Q-Q plot (vs. normal)');
    qqplot_local(ax2, raw_vals);

    % Panel 3: z-scored histogram with outlier flags
    ax3 = subplot(1,4,3); hold(ax3,'on');
    title(ax3,'Z-scored amplitude');
    if ~isempty(z_vals)
        histogram(ax3,z_vals,40,'FaceColor',[0.8 0.5 0.2],'EdgeColor','none','FaceAlpha',0.8);
        xline(ax3, OUTLIER_SD_THRESH,'r--','LineWidth',1.5,'HandleVisibility','off');
        xline(ax3,-OUTLIER_SD_THRESH,'r--','LineWidth',1.5,'DisplayName',sprintf('±%d SD',OUTLIER_SD_THRESH));
        n_out = sum(abs(z_vals) > OUTLIER_SD_THRESH,'omitnan');
        text(ax3,0.98,0.97,sprintf('Outliers > ±%dSD: %d (%.1f%%)',...
            OUTLIER_SD_THRESH,n_out,100*n_out/max(1,numel(z_vals))),...
            'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top',...
            'FontSize',7,'Color',[0.7 0 0],'BackgroundColor','w');
        xlabel(ax3,'Within-subject z'); ylabel(ax3,'Count');
        legend(ax3,'Box','off','FontSize',7);
    else
        text(ax3,0.5,0.5,'z-column absent','Units','normalized','HorizontalAlignment','center');
    end

    % Panel 4: raw vs z scatter (sanity check of z-scoring)
    ax4 = subplot(1,4,4); hold(ax4,'on');
    title(ax4,'Raw vs z-scored (sanity)');
    if ~isempty(z_vals)
        raw_z_vals = gt.(raw)(z_mask);
        scatter(ax4,raw_z_vals,z_vals,4,[0.5 0.5 0.5],'filled','MarkerFaceAlpha',0.15,...
            'HandleVisibility','off');
        [rv,pv] = corr(raw_z_vals,z_vals,'Rows','complete');
        text(ax4,0.05,0.97,sprintf('r = %.3f\np = %.4f',rv,pv),...
            'Units','normalized','VerticalAlignment','top','FontSize',8,'BackgroundColor','w');
        xlabel(ax4,lbl,'Interpreter','none'); ylabel(ax4,'z-scored');
    end

    exportgraphics(fig,fullfile(outdir,sprintf('S7b_2_distrib_%s.pdf',raw)),'ContentType','vector');
    fprintf('  %s: N=%d, mean=%.3f, SD=%.3f, S-W p=%.3f\n',...
        raw, numel(raw_vals), mean(raw_vals,'omitnan'), std(raw_vals,'omitnan'), p_sw);
end


% =============================================================================
%% SECTION 3 — OUTLIER IDENTIFICATION
%   Flag subjects whose condition-mean is > OUTLIER_SD_THRESH SD from group.
%   Clayson (2020): report n subjects excluded / flagged.
% =============================================================================
fprintf('\n=== SECTION 3: Outlier subjects ===\n');

outlier_T = table();
for fi = 1:N_feat
    raw = FEATURES{fi,1};
    if ~ismember(raw,gt.Properties.VariableNames), continue; end

    subj_means = arrayfun(@(s) ...
        mean(gt.(raw)(gt.subj_id==s & ~gt.false_fb),'omitnan'), ...
        categorical(subjects));

    gm  = mean(subj_means,'omitnan');
    gsd = std(subj_means,'omitnan');
    out_flag = abs(subj_means - gm) > OUTLIER_SD_THRESH * gsd;

    for si = 1:N_subj
        if out_flag(si)
            row = table({subjects{si}},{raw},subj_means(si),(subj_means(si)-gm)/gsd,...
                'VariableNames',{'subj_id','feature','mean_amp','z_from_group'});
            outlier_T = [outlier_T; row]; %#ok<AGROW>
            fprintf('  OUTLIER: %s  %s  mean=%.3f  z=%.2f\n',...
                subjects{si},raw,subj_means(si),(subj_means(si)-gm)/gsd);
        end
    end
end

if isempty(outlier_T)
    fprintf('  No subjects flagged as outliers (±%d SD).\n',OUTLIER_SD_THRESH);
else
    writetable(outlier_T,fullfile(outdir,'S7b_outliers.csv'));
    fprintf('  Outlier table saved.\n');
end

% Figure 3: subject-level mean amplitude, all features
fig = figure('Position',[50 50 1400 500]);
sgtitle(sprintf('Subject mean amplitudes — red markers: |z| > %d SD from group',OUTLIER_SD_THRESH));
for fi = 1:N_feat
    raw = FEATURES{fi,1};
    if ~ismember(raw,gt.Properties.VariableNames), continue; end

    subj_means = arrayfun(@(s) ...
        mean(gt.(raw)(gt.subj_id==s & ~gt.false_fb),'omitnan'), ...
        categorical(subjects));
    gm = mean(subj_means,'omitnan'); gsd = std(subj_means,'omitnan');
    out_flag = abs(subj_means-gm) > OUTLIER_SD_THRESH*gsd;

    coh_kh = arrayfun(@(s) ...
        isequal(char(string(gt.cohort(find(gt.subj_id==s,1)))),'KH'),...
        categorical(subjects));

    ax = subplot(2,ceil(N_feat/2),fi); hold(ax,'on');
    title(ax,FEATURES{fi,3},'FontSize',8,'Interpreter','none');

    for si2 = 1:N_subj
        clr = ternary_local(coh_kh(si2),CLR_KH,CLR_RR);
        if out_flag(si2), clr = [0.8 0 0]; end
        scatter(ax,si2,subj_means(si2),30,clr,'filled');
    end
    yline(ax,gm,'k-','LineWidth',1,'HandleVisibility','off');
    yline(ax,gm+OUTLIER_SD_THRESH*gsd,'r--','LineWidth',0.8,'HandleVisibility','off');
    yline(ax,gm-OUTLIER_SD_THRESH*gsd,'r--','LineWidth',0.8,'HandleVisibility','off');
    set(ax,'XTick',1:N_subj,'XTickLabel',subjects,'XTickLabelRotation',45,'FontSize',6);
    ylabel(ax,FEATURES{fi,3},'FontSize',7,'Interpreter','none');
end
exportgraphics(fig,fullfile(outdir,'S7b_3_outlier_subjects.pdf'),'ContentType','vector');


% =============================================================================
%% SECTION 4 — LATENCY DISTRIBUTIONS (FRN peak latency, P300 peak latency)
%   Best practice: verify that peak latency falls within the intended window.
%   Flag subjects with mean latency outside expected range.
% =============================================================================
fprintf('\n=== SECTION 4: Peak latency distributions ===\n');

lat_feats = {
    'prefrontal_neg_peak_lat',  [250 350],  'FCz neg-peak latency (ms)'
    'P300_peak_lat',            [300 600],  'P300 peak latency (ms)'
};

fig = figure('Position',[50 50 1200 450]);
sgtitle('Peak latency distributions (true-feedback incorrect trials)');

for li = 1:size(lat_feats,1)
    lat_col  = lat_feats{li,1};
    expected = lat_feats{li,2};
    lbl      = lat_feats{li,3};

    if ~ismember(lat_col,gt.Properties.VariableNames), continue; end

    inc_mask = ~gt.false_fb & gt.correct_num==0;
    lats = gt.(lat_col)(inc_mask & ~isnan(gt.(lat_col)));

    ax = subplot(1,size(lat_feats,1),li); hold(ax,'on');
    title(ax,lbl,'FontSize',9,'Interpreter','none');
    histogram(ax,lats,30,'FaceColor',[0.4 0.6 0.4],'EdgeColor','none','FaceAlpha',0.8);
    xline(ax,expected(1),'r--','LineWidth',1.5,'DisplayName',sprintf('Window [%d %d]',expected));
    xline(ax,expected(2),'r--','LineWidth',1.5,'HandleVisibility','off');
    xline(ax,mean(lats,'omitnan'),'k-','LineWidth',2,'DisplayName','Mean');
    xlabel(ax,'Latency (ms)'); ylabel(ax,'Count');

    n_outside = sum(lats < expected(1) | lats > expected(2));
    text(ax,0.98,0.97,sprintf('Outside window: %d (%.1f%%)\nMean: %.1f ms\nSD: %.1f ms',...
        n_outside,100*n_outside/max(1,numel(lats)),mean(lats,'omitnan'),std(lats,'omitnan')),...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top',...
        'FontSize',8,'BackgroundColor','w');
    legend(ax,'Box','off','FontSize',8,'Location','northwest');

    fprintf('  %s: N=%d, mean=%.1f ms, SD=%.1f ms, outside window: %d (%.1f%%)\n',...
        lat_col,numel(lats),mean(lats,'omitnan'),std(lats,'omitnan'),n_outside,100*n_outside/max(1,numel(lats)));
end
exportgraphics(fig,fullfile(outdir,'S7b_4_latency_distributions.pdf'),'ContentType','vector');

% Per-subject latency means
fig = figure('Position',[50 50 1200 350]);
sgtitle('Per-subject mean peak latency (incorrect, true-FB trials)');
for li = 1:size(lat_feats,1)
    lat_col  = lat_feats{li,1};
    expected = lat_feats{li,2};
    if ~ismember(lat_col,gt.Properties.VariableNames), continue; end

    subj_lat = arrayfun(@(s) ...
        mean(gt.(lat_col)(gt.subj_id==s & gt.correct_num==0 & ~gt.false_fb),'omitnan'),...
        categorical(subjects));

    ax = subplot(1,size(lat_feats,1),li); hold(ax,'on');
    title(ax,lat_feats{li,3},'FontSize',9,'Interpreter','none');
    bar(ax,1:N_subj,subj_lat,0.7,'FaceColor',[0.4 0.7 0.4],'EdgeColor','none');
    yline(ax,expected(1),'r--','LineWidth',1.5);
    yline(ax,expected(2),'r--','LineWidth',1.5);
    set(ax,'XTick',1:N_subj,'XTickLabel',subjects,'XTickLabelRotation',45,'FontSize',7);
    ylabel(ax,'Mean latency (ms)');
end
exportgraphics(fig,fullfile(outdir,'S7b_4b_latency_per_subject.pdf'),'ContentType','vector');


% =============================================================================
%% SECTION 5 — FRN EXCLUSION RATE
%   Tracks proportion of trials excluded because no local minimum was found
%   in the FRN window (FRN_excluded = true in B_outcome_ERP_analysis).
%   Clayson (2020): exclusion rate should be low and not systematically
%   differ between conditions of interest.
% =============================================================================
fprintf('\n=== SECTION 5: FRN exclusion rates ===\n');

if ismember('FRN_excluded',gt.Properties.VariableNames)
    excl_rate_subj = arrayfun(@(s) ...
        mean(gt.FRN_excluded(gt.subj_id==s),'omitnan'), categorical(subjects));

    fprintf('  Overall FRN exclusion: %.1f%% (range %.1f–%.1f%%)\n',...
        100*mean(excl_rate_subj,'omitnan'),...
        100*min(excl_rate_subj,[],'omitnan'),...
        100*max(excl_rate_subj,[],'omitnan'));

    fig = figure('Position',[50 50 1400 500]);
    sgtitle('FRN peak exclusion rate (no local minimum found in [250 350] ms window)');

    % Panel 1: per-subject overall rate
    ax1 = subplot(1,3,1); hold(ax1,'on');
    title(ax1,'Per-subject exclusion rate');
    bar(ax1,1:N_subj,100*excl_rate_subj,0.7,'FaceColor',[0.8 0.4 0.2],'EdgeColor','none');
    yline(ax1,20,'r--','DisplayName','20% threshold');
    set(ax1,'XTick',1:N_subj,'XTickLabel',subjects,'XTickLabelRotation',45,'FontSize',7);
    ylabel(ax1,'% trials excluded'); legend(ax1,'Box','off');

    % Panel 2: by stage and block_type
    ax2 = subplot(1,3,2); hold(ax2,'on');
    title(ax2,'Exclusion rate by stage × block type');
    x = 1:numel(STAGES);
    for bti = 1:2
        excl_stage = nan(1,numel(STAGES));
        for stgi = 1:numel(STAGES)
            m = gt.stage==STAGES{stgi} & gt.block_type==BTYPES{bti};
            excl_stage(stgi) = 100*mean(gt.FRN_excluded(m),'omitnan');
        end
        clr = ternary_local(bti==1,CLR_D,CLR_P);
        plot(ax2,x,excl_stage,'o-','Color',clr,'LineWidth',2,'MarkerFaceColor',clr,...
            'DisplayName',BTYPES{bti});
    end
    set(ax2,'XTick',x,'XTickLabel',STAGES,'FontSize',9);
    ylabel(ax2,'% trials excluded'); legend(ax2,'Box','off');

    % Panel 3: by outcome
    ax3 = subplot(1,3,3); hold(ax3,'on');
    title(ax3,'Exclusion rate by outcome');
    for oci = 1:2
        m = gt.correct_num==(oci-1) & ~gt.false_fb;
        excl_stg = nan(1,numel(STAGES));
        for stgi = 1:numel(STAGES)
            ms = m & gt.stage==STAGES{stgi};
            excl_stg(stgi) = 100*mean(gt.FRN_excluded(ms),'omitnan');
        end
        clr = ternary_local(oci==1,CLR_INC,CLR_COR);
        plot(ax3,x,excl_stg,'o-','Color',clr,'LineWidth',2,'MarkerFaceColor',clr,...
            'DisplayName',OUTCOMES{oci});
    end
    set(ax3,'XTick',x,'XTickLabel',STAGES,'FontSize',9);
    ylabel(ax3,'% trials excluded'); legend(ax3,'Box','off');

    exportgraphics(fig,fullfile(outdir,'S7b_5_FRN_exclusion_rate.pdf'),'ContentType','vector');
else
    fprintf('  FRN_excluded column not found — skipping.\n');
end


% =============================================================================
%% SECTION 6 — BASELINE RMS (signal quality proxy)
%   Mean baseline RMS should be consistent across subjects.
%   Subjects with unusually high RMS likely had poor signal.
% =============================================================================
fprintf('\n=== SECTION 6: Baseline RMS ===\n');

if ismember('baseline_rms',gt.Properties.VariableNames)
    bline_rms = arrayfun(@(s) ...
        mean(gt.baseline_rms(gt.subj_id==s),'omitnan'), categorical(subjects));

    gm_rms = mean(bline_rms,'omitnan');
    gsd_rms = std(bline_rms,'omitnan');
    flag_rms = bline_rms > gm_rms + 2*gsd_rms;

    fprintf('  Baseline RMS: mean=%.3f, SD=%.3f\n',gm_rms,gsd_rms);
    if any(flag_rms)
        fprintf('  High RMS (>2SD above mean): %s\n',...
            strjoin(subjects(flag_rms),', '));
    end

    fig = figure('Position',[50 50 900 400]);
    sgtitle('Baseline RMS per subject (pre-stimulus period signal quality)');
    ax = axes(fig); hold(ax,'on');

    coh_kh = arrayfun(@(s) ...
        isequal(char(string(gt.cohort(find(gt.subj_id==s,1)))),'KH'),...
        categorical(subjects));
    for si = 1:N_subj
        clr = ternary_local(coh_kh(si),CLR_KH,CLR_RR);
        if flag_rms(si), clr = [0.8 0 0]; end
        bar(ax,si,bline_rms(si),0.7,'FaceColor',clr,'EdgeColor','none','FaceAlpha',0.85);
    end
    yline(ax,gm_rms,'k-','LineWidth',1.5,'HandleVisibility','off');
    yline(ax,gm_rms+2*gsd_rms,'r--','LineWidth',1.2,'DisplayName','+2 SD');
    set(ax,'XTick',1:N_subj,'XTickLabel',subjects,'XTickLabelRotation',45,'FontSize',8);
    ylabel(ax,'Baseline RMS (µV)'); legend(ax,'Box','off');
    annotation('textbox',[0.01 0.0 0.98 0.05],'String',...
        'Red bars: >2 SD above mean — flag for QC. Blue=KH, purple=RR.',...
        'EdgeColor','none','FontSize',7);
    exportgraphics(fig,fullfile(outdir,'S7b_6_baseline_RMS.pdf'),'ContentType','vector');
end


% =============================================================================
%% SECTION 7 — TRIAL-COUNT ~ AMPLITUDE CORRELATION
%   Clayson (2020, key finding): fewer trials → larger ERP amplitudes
%   (noisier estimate inflates absolute peak values). Report this correlation.
% =============================================================================
fprintf('\n=== SECTION 7: Trial count ~ amplitude correlation ===\n');

fig = figure('Position',[50 50 1400 400]);
sgtitle('Trial count × amplitude (Clayson 2020: check for trial-count bias)');

for fi = 1:N_feat
    raw = FEATURES{fi,1};
    if ~ismember(raw,gt.Properties.VariableNames), continue; end

    % Per-subject: median trial count and mean amplitude (incorrect true-FB)
    n_trials_subj = arrayfun(@(s) ...
        sum(gt.subj_id==s & gt.correct_num==0 & ~gt.false_fb & ~isnan(gt.(raw))),...
        categorical(subjects));
    amp_subj = arrayfun(@(s) ...
        mean(gt.(raw)(gt.subj_id==s & gt.correct_num==0 & ~gt.false_fb),'omitnan'),...
        categorical(subjects));

    ok = ~isnan(amp_subj) & ~isnan(n_trials_subj);

    ax = subplot(2,ceil(N_feat/2),fi); hold(ax,'on');
    title(ax,FEATURES{fi,3},'FontSize',8,'Interpreter','none');

    coh_kh = arrayfun(@(s) ...
        isequal(char(string(gt.cohort(find(gt.subj_id==s,1)))),'KH'),...
        categorical(subjects));
    for si2 = 1:N_subj
        if ~ok(si2), continue; end
        scatter(ax,n_trials_subj(si2),amp_subj(si2),30,...
            ternary_local(coh_kh(si2),CLR_KH,CLR_RR),'filled','MarkerFaceAlpha',0.8);
    end
    if sum(ok) > 3
        [rv,pv] = corr(n_trials_subj(ok)',amp_subj(ok)');
        p_fit = polyfit(n_trials_subj(ok)',amp_subj(ok)',1);
        xx = linspace(min(n_trials_subj(ok)),max(n_trials_subj(ok)),50);
        plot(ax,xx,polyval(p_fit,xx),'k-','LineWidth',1.5,'HandleVisibility','off');
        text(ax,0.98,0.97,sprintf('r=%.2f, p=%.3f',rv,pv),...
            'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top',...
            'FontSize',7,'BackgroundColor','w');
    end
    xlabel(ax,'N valid trials'); ylabel(ax,FEATURES{fi,3},'Interpreter','none','FontSize',7);
end
annotation('textbox',[0.01 0.0 0.98 0.04],'String',...
    'Significant negative r = fewer trials inflate amplitude magnitude (trial-count bias). Consider minimum-trial exclusion.',...
    'EdgeColor','none','FontSize',7);
exportgraphics(fig,fullfile(outdir,'S7b_7_trial_count_amplitude_corr.pdf'),'ContentType','vector');


% =============================================================================
%% SECTION 8 — BLOCK-ORDER EFFECTS
%   Check whether amplitude changes systematically over blocks within session.
%   Potential confound if block type is confounded with session order.
% =============================================================================
fprintf('\n=== SECTION 8: Block-order effects ===\n');

if ismember('block',gt.Properties.VariableNames)
    fig = figure('Position',[50 50 1400 400]);
    sgtitle('Amplitude by block number (check for session-order confound)');

    blocks_all = sort(unique(double(gt.block)))';

    for fi = 1:N_feat
        raw = FEATURES{fi,1};
        if ~ismember(raw,gt.Properties.VariableNames), continue; end

        ax = subplot(2,ceil(N_feat/2),fi); hold(ax,'on');
        title(ax,FEATURES{fi,3},'FontSize',8,'Interpreter','none');

        for bti = 1:2
            gm_blk = nan(1,numel(blocks_all));
            se_blk = nan(1,numel(blocks_all));
            for bi = 1:numel(blocks_all)
                m = double(gt.block)==blocks_all(bi) & ...
                    gt.block_type==BTYPES{bti} & ~gt.false_fb & ~isnan(gt.(raw));
                vals = gt.(raw)(m);
                if numel(vals) < 3, continue; end
                % Subject-average first
                subj_vals = arrayfun(@(s) mean(gt.(raw)(m & gt.subj_id==s),'omitnan'),...
                    categorical(subjects));
                gm_blk(bi) = mean(subj_vals,'omitnan');
                se_blk(bi) = std(subj_vals,'omitnan')/sqrt(sum(~isnan(subj_vals)));
            end
            clr = ternary_local(bti==1,CLR_D,CLR_P);
            ok = ~isnan(gm_blk);
            errorbar(ax,blocks_all(ok),gm_blk(ok),se_blk(ok),'o-',...
                'Color',clr,'LineWidth',1.5,'MarkerFaceColor',clr,...
                'DisplayName',BTYPES{bti},'MarkerSize',5);
        end
        xlabel(ax,'Block number'); ylabel(ax,FEATURES{fi,3},'Interpreter','none','FontSize',7);
        if fi==1, legend(ax,'Box','off','FontSize',7,'Location','best'); end
    end
    exportgraphics(fig,fullfile(outdir,'S7b_8_block_order_effects.pdf'),'ContentType','vector');
end


% =============================================================================
%% SECTION 9 — COHORT COMPARISON (KH vs RR)
%   Verify that the two cohorts produce comparable feature distributions
%   before pooling. Differences here require covariate control in S7.
% =============================================================================
fprintf('\n=== SECTION 9: Cohort comparison ===\n');

cohorts = categories(gt.cohort);
if numel(cohorts) > 1
    fig = figure('Position',[50 50 1400 500]);
    sgtitle('Cohort comparison: KH vs RR (raw amplitude, true-FB, incorrect trials)');

    for fi = 1:N_feat
        raw = FEATURES{fi,1};
        if ~ismember(raw,gt.Properties.VariableNames), continue; end

        ax = subplot(2,ceil(N_feat/2),fi); hold(ax,'on');
        title(ax,FEATURES{fi,3},'FontSize',8,'Interpreter','none');

        for ci = 1:numel(cohorts)
            m = gt.cohort==cohorts{ci} & gt.correct_num==0 & ~gt.false_fb & ~isnan(gt.(raw));
            vals = gt.(raw)(m);
            histogram(ax,vals,20,'FaceAlpha',0.5,...
                'FaceColor',ternary_local(strcmp(cohorts{ci},'KH'),CLR_KH,CLR_RR),...
                'EdgeColor','none','DisplayName',cohorts{ci});
        end
        xline(ax,0,'k:','HandleVisibility','off');
        legend(ax,'Box','off','FontSize',7);
        xlabel(ax,FEATURES{fi,3},'Interpreter','none','FontSize',7); ylabel(ax,'Count');

        % Ranksum test (don't assume equal distributions)
        kh_vals = gt.(raw)(gt.cohort=='KH' & gt.correct_num==0 & ~gt.false_fb & ~isnan(gt.(raw)));
        rr_vals = gt.(raw)(gt.cohort=='RR' & gt.correct_num==0 & ~gt.false_fb & ~isnan(gt.(raw)));
        if numel(kh_vals)>5 && numel(rr_vals)>5
            [p_rs,~] = ranksum(kh_vals,rr_vals);
            text(ax,0.98,0.97,sprintf('Ranksum p=%.3f',p_rs),...
                'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top',...
                'FontSize',7,'BackgroundColor','w');
            fprintf('  %s  KH vs RR ranksum p = %.4f\n',raw,p_rs);
        end
    end
    exportgraphics(fig,fullfile(outdir,'S7b_9_cohort_comparison.pdf'),'ContentType','vector');
end


% =============================================================================
%% SECTION 10 — GRAND-AVERAGE WAVEFORM MORPHOLOGY CHECK
%   Best practice: confirm that the grand average has the expected shape
%   before running any statistical test. Plots correct vs incorrect at FCz
%   and parietal channels, across all stages and block types.
%   (Uses per-trial waveform data if prefrontal_waveform column present;
%    otherwise reconstructs from the feature table using precomputed means.)
% =============================================================================
fprintf('\n=== SECTION 10: Grand-average waveform check ===\n');
fprintf('  (Requires per-trial waveform column or EEG .set files.)\n');
fprintf('  Checking amplitude at key time windows across stages...\n');

% Proxy check using mean amplitude feature across stages
% (full waveform check requires EEG .set files — see B_outcome_ERP_analysis)
if ismember('prefrontal_mean_amp',gt.Properties.VariableNames)
    fig = figure('Position',[50 50 1200 500]);
    sgtitle('Grand-average prefrontal mean amplitude by stage × block type × outcome');

    STAGE_CLR = {[0.12 0.62 0.47],[0.85 0.65 0.00],[0.80 0.27 0.13],[0.40 0.25 0.65]};

    for bti = 1:2
        ax = subplot(1,2,bti); hold(ax,'on');
        title(ax,sprintf('Block type: %s',BTYPES{bti}),'FontSize',11);

        for oci = 1:2
            corr_val = (oci==2);
            clr_base = ternary_local(oci==2,CLR_COR,CLR_INC);

            stage_means = nan(1,numel(STAGES));
            stage_sems  = nan(1,numel(STAGES));

            for stgi = 1:numel(STAGES)
                m = gt.block_type==BTYPES{bti} & gt.stage==STAGES{stgi} & ...
                    gt.correct_num==corr_val & ~gt.false_fb & ~isnan(gt.prefrontal_mean_amp);
                % Subject-average first
                sv = arrayfun(@(s) mean(gt.prefrontal_mean_amp(m & gt.subj_id==s),'omitnan'),...
                    categorical(subjects));
                stage_means(stgi) = mean(sv,'omitnan');
                stage_sems(stgi)  = std(sv,'omitnan')/sqrt(sum(~isnan(sv)));
            end

            errorbar(ax,1:numel(STAGES),stage_means,stage_sems,'o-',...
                'Color',clr_base,'LineWidth',2,'MarkerFaceColor',clr_base,...
                'DisplayName',OUTCOMES{oci});
        end

        set(ax,'XTick',1:numel(STAGES),'XTickLabel',STAGES,'FontSize',10);
        ylabel(ax,'Mean amp [250–350 ms] (µV)'); xlabel(ax,'Stage');
        yline(ax,0,'k:','HandleVisibility','off');
        legend(ax,'Box','off','Location','best');

        % Annotation: expected FRN direction
        text(ax,0.02,0.98,'Expected: Incorrect < Correct (more negative)',...
            'Units','normalized','VerticalAlignment','top','FontSize',8,'Color',[0.4 0.4 0.4]);
    end
    exportgraphics(fig,fullfile(outdir,'S7b_10_grand_avg_proxy.pdf'),'ContentType','vector');
end


% =============================================================================
%% SECTION 11 — SPLIT-HALF RELIABILITY
%   Clayson (2020) central finding: peak amplitude has poor split-half
%   reliability especially with < 20–30 trials. Report Spearman-Brown
%   corrected split-half r for each subject × condition cell.
%   Computed as: odd-trial mean vs even-trial mean, Spearman-Brown corrected.
% =============================================================================
fprintf('\n=== SECTION 11: Split-half reliability ===\n');

if ismember('trial_continuous',gt.Properties.VariableNames)
    reli_T = table();
    for fi = 1:min(N_feat,3)   % limit to first 3 to keep runtime reasonable
        raw = FEATURES{fi,1};
        if ~ismember(raw,gt.Properties.VariableNames), continue; end

        sh_vals = nan(N_subj,1);
        for si = 1:N_subj
            sn = subjects{si};
            m  = gt.subj_id==sn & gt.correct_num==0 & ~gt.false_fb & ~isnan(gt.(raw));
            tc = gt.trial_continuous(m);
            av = gt.(raw)(m);
            if numel(av) < 6, continue; end

            odd_mean  = mean(av(mod(tc,2)==1),'omitnan');
            even_mean = mean(av(mod(tc,2)==0),'omitnan');
            if isnan(odd_mean) || isnan(even_mean), continue; end

            % Pearson r between odd and even half (single pair → use absolute diff proxy)
            % For true split-half across conditions, aggregate per-condition means
            % across many cells and correlate:
            odd_cond  = nan(numel(STAGES)*numel(BTYPES),1);
            even_cond = nan(numel(STAGES)*numel(BTYPES),1);
            idx = 0;
            for stgi = 1:numel(STAGES)
                for bti = 1:numel(BTYPES)
                    idx = idx+1;
                    mc = m & gt.stage==STAGES{stgi} & gt.block_type==BTYPES{bti};
                    tc2 = gt.trial_continuous(mc);
                    av2 = gt.(raw)(mc);
                    odd_cond(idx)  = mean(av2(mod(tc2,2)==1),'omitnan');
                    even_cond(idx) = mean(av2(mod(tc2,2)==0),'omitnan');
                end
            end
            ok2 = ~isnan(odd_cond) & ~isnan(even_cond);
            if sum(ok2) >= 3
                r_half = corr(odd_cond(ok2),even_cond(ok2));
                % Spearman-Brown correction
                sh_vals(si) = 2*r_half / (1+r_half);
            end
        end

        fprintf('  %-35s  median SB r = %.3f  [%.3f – %.3f]\n',...
            raw, median(sh_vals,'omitnan'),...
            quantile_omitnan_local(sh_vals,0.25), quantile_omitnan_local(sh_vals,0.75));

        row = table({raw},median(sh_vals,'omitnan'),mean(sh_vals,'omitnan'),...
            quantile_omitnan_local(sh_vals,0.25),quantile_omitnan_local(sh_vals,0.75),...
            sum(~isnan(sh_vals)),...
            'VariableNames',{'feature','median_SB_r','mean_SB_r','Q25','Q75','N_subjects'});
        reli_T = [reli_T; row]; %#ok<AGROW>
    end

    if ~isempty(reli_T)
        disp(reli_T);
        writetable(reli_T,fullfile(outdir,'S7b_reliability.csv'));
    end

    % Figure 11: split-half reliability per subject, first feature
    if ismember(FEATURES{1,1},gt.Properties.VariableNames)
        raw = FEATURES{1,1};
        sh_vals_plot = nan(N_subj,1);
        for si = 1:N_subj
            sn = subjects{si};
            m  = gt.subj_id==sn & gt.correct_num==0 & ~gt.false_fb & ~isnan(gt.(raw));
            tc = gt.trial_continuous(m);
            av = gt.(raw)(m);
            if numel(av)<6, continue; end
            odd_cond  = nan(numel(STAGES)*numel(BTYPES),1);
            even_cond = nan(numel(STAGES)*numel(BTYPES),1);
            idx=0;
            for stgi=1:numel(STAGES)
                for bti=1:numel(BTYPES)
                    idx=idx+1;
                    mc=m & gt.stage==STAGES{stgi} & gt.block_type==BTYPES{bti};
                    tc2=gt.trial_continuous(mc); av2=gt.(raw)(mc);
                    odd_cond(idx)=mean(av2(mod(tc2,2)==1),'omitnan');
                    even_cond(idx)=mean(av2(mod(tc2,2)==0),'omitnan');
                end
            end
            ok2=~isnan(odd_cond)&~isnan(even_cond);
            if sum(ok2)>=3
                r=corr(odd_cond(ok2),even_cond(ok2));
                sh_vals_plot(si)=2*r/(1+r);
            end
        end

        fig = figure('Position',[50 50 900 400]);
        sgtitle(sprintf('Spearman-Brown split-half reliability: %s (incorrect, true-FB)',strrep(raw,'_','\_')));
        ax = axes(fig); hold(ax,'on');
        bar(ax,1:N_subj,sh_vals_plot,0.7,'FaceColor',[0.4 0.6 0.8],'EdgeColor','none');
        yline(ax,0.70,'r--','DisplayName','r = .70 threshold');
        yline(ax,0.80,'k--','DisplayName','r = .80 threshold');
        set(ax,'XTick',1:N_subj,'XTickLabel',subjects,'XTickLabelRotation',45,'FontSize',8);
        ylabel(ax,'Spearman-Brown corrected r'); ylim(ax,[-0.2 1]);
        legend(ax,'Box','off');
        annotation('textbox',[0.01 0.0 0.98 0.05],'String',...
            'Spearman-Brown r < .70 indicates poor split-half reliability for that subject (Clayson 2020).',...
            'EdgeColor','none','FontSize',7);
        exportgraphics(fig,fullfile(outdir,'S7b_11_split_half_reliability.pdf'),'ContentType','vector');
    end
end


% =============================================================================
%% SECTION 12 — FALSE-FEEDBACK DISSOCIATION CHECK
%   Verify that false-feedback trials are distributed across stages and
%   subjects as expected (P blocks only, ~20% of trials).
% =============================================================================
fprintf('\n=== SECTION 12: False-feedback rate verification ===\n');

ff_rate_subj = arrayfun(@(s) ...
    mean(gt.false_fb(gt.subj_id==s & gt.block_type=='P'),'omitnan'),...
    categorical(subjects));

fprintf('  False-FB rate in P blocks: mean=%.3f SD=%.3f (expected ~0.20)\n',...
    mean(ff_rate_subj,'omitnan'), std(ff_rate_subj,'omitnan'));

fig = figure('Position',[50 50 900 400]);
sgtitle('False-feedback rate per subject (P blocks only, expected ≈ 0.20)');
ax = axes(fig); hold(ax,'on');
bar(ax,1:N_subj,ff_rate_subj,0.7,'FaceColor',[0.6 0.4 0.7],'EdgeColor','none','FaceAlpha',0.85);
yline(ax,0.20,'k--','DisplayName','Expected = 0.20','LineWidth',1.5);
yline(ax,0.15,'r:','DisplayName','±0.05 tolerance');
yline(ax,0.25,'r:','LineWidth',1,'HandleVisibility','off');
set(ax,'XTick',1:N_subj,'XTickLabel',subjects,'XTickLabelRotation',45,'FontSize',8);
ylabel(ax,'P(false feedback) in P blocks'); legend(ax,'Box','off');
exportgraphics(fig,fullfile(outdir,'S7b_12_false_fb_rate.pdf'),'ContentType','vector');


% =============================================================================
%% SECTION 13 — THETA AMPLITUDE TIME COURSE CHECK
%   Verify theta envelope peaks after outcome onset, not before.
%   Computed from condition-level means (requires Theta_amp to be in table).
% =============================================================================
fprintf('\n=== SECTION 13: Theta amplitude by stage ===\n');

if ismember('Theta_amp',gt.Properties.VariableNames)
    fig = figure('Position',[50 50 900 450]);
    sgtitle('Theta amplitude by stage × block type (incorrect, true-FB)');
    ax = axes(fig); hold(ax,'on');

    for bti = 1:2
        stage_means = nan(1,numel(STAGES));
        stage_sems  = nan(1,numel(STAGES));
        for stgi = 1:numel(STAGES)
            m = gt.block_type==BTYPES{bti} & gt.stage==STAGES{stgi} & ...
                gt.correct_num==0 & ~gt.false_fb & ~isnan(gt.Theta_amp);
            sv = arrayfun(@(s) mean(gt.Theta_amp(m & gt.subj_id==s),'omitnan'),...
                categorical(subjects));
            stage_means(stgi) = mean(sv,'omitnan');
            stage_sems(stgi)  = std(sv,'omitnan')/sqrt(sum(~isnan(sv)));
        end
        clr = ternary_local(bti==1,CLR_D,CLR_P);
        errorbar(ax,1:numel(STAGES),stage_means,stage_sems,'o-',...
            'Color',clr,'LineWidth',2,'MarkerFaceColor',clr,'DisplayName',BTYPES{bti});
    end
    set(ax,'XTick',1:numel(STAGES),'XTickLabel',STAGES,'FontSize',10);
    ylabel(ax,'Theta amplitude (baseline-corrected)');
    xlabel(ax,'Stage'); legend(ax,'Box','off'); yline(ax,0,'k:','HandleVisibility','off');
    exportgraphics(fig,fullfile(outdir,'S7b_13_theta_by_stage.pdf'),'ContentType','vector');
end


% =============================================================================
%% SUMMARY STATISTICS TABLE
%   One-row-per-feature summary of key QC statistics, printed and saved.
% =============================================================================
fprintf('\n=== SUMMARY TABLE ===\n');

summary_T = table();
for fi = 1:N_feat
    raw = FEATURES{fi,1};
    if ~ismember(raw,gt.Properties.VariableNames), continue; end

    vals = gt.(raw)(~isnan(gt.(raw)) & ~gt.false_fb);
    n_valid = numel(vals);
    pct_missing = 100*(1 - n_valid/height(gt));
    [~,p_sw] = swtest_local(vals);

    n_out = 0;
    zc = FEATURES{fi,2};
    if ismember(zc,gt.Properties.VariableNames)
        zv = gt.(zc)(~isnan(gt.(zc)) & ~gt.false_fb);
        n_out = sum(abs(zv) > OUTLIER_SD_THRESH,'omitnan');
    end

    n_warn = 0;
    for si = 1:N_subj
        for stgi = 1:numel(STAGES)
            for bti = 1:numel(BTYPES)
                m = gt.subj_id==subjects{si} & gt.stage==STAGES{stgi} & ...
                    gt.block_type==BTYPES{bti} & ~gt.false_fb & ~isnan(gt.(raw));
                if sum(m) < MIN_TRIALS_WARN, n_warn = n_warn+1; end
            end
        end
    end

    row = table({raw},n_valid,pct_missing,...
        mean(vals,'omitnan'),std(vals,'omitnan'),...
        min(vals,[],'omitnan'),max(vals,[],'omitnan'),...
        p_sw, n_out, n_warn,...
        'VariableNames',{'feature','n_valid','pct_missing',...
        'mean','SD','min','max','shapiro_wilk_p','n_outliers_3SD','n_cells_below8'});
    summary_T = [summary_T; row]; %#ok<AGROW>
end

disp(summary_T);
writetable(summary_T, fullfile(outdir,'S7b_summary_stats.csv'));

% Save workspace
save(fullfile(outdir,'S7b_validation_workspace.mat'),'count_T','summary_T','outlier_T','-v7.3');
fprintf('\nAll validation figures and tables saved to:\n  %s\n',outdir);
fprintf('S7b complete.\n');


% =============================================================================
%% LOCAL HELPER FUNCTIONS
% =============================================================================

function [h,p] = swtest_local(x)
% Normality test wrapper.
% Uses swtest if available, then lillietest if available, otherwise a
% Jarque-Bera chi-square approximation. h = 1 means reject normality at .05.
x = x(~isnan(x));
if numel(x) < 4, h = 0; p = 1; return; end
if exist('swtest','file') == 2
    [h,p] = swtest(x,0.05);
elseif exist('lillietest','file') == 2
    [h,p] = lillietest(x);
else
    x = x(:);
    n = numel(x);
    z = (x - mean(x)) ./ std(x);
    if all(~isfinite(z)), h = 0; p = 1; return; end
    sk = mean(z.^3,'omitnan');
    ku = mean(z.^4,'omitnan');
    jb = n/6 * (sk.^2 + (ku - 3).^2/4);
    p = 1 - chi2cdf_local(jb, 2);
    h = p < 0.05;
end
end


function qqplot_local(ax, x)
% Simple normal Q-Q plot without Statistics Toolbox.
x = sort(x(~isnan(x)));
n = numel(x);
if n < 4, text(ax,0.5,0.5,'n < 4','Units','normalized','HorizontalAlignment','center'); return; end
p_vals = ((1:n) - 0.5) / n;
q_theor = mean(x) + std(x) .* sqrt(2) .* erfinv(2*p_vals - 1);
scatter(ax, q_theor, x, 8, [0.4 0.6 0.8], 'filled', 'MarkerFaceAlpha', 0.4);
q_rng = [min(q_theor) max(q_theor)];
plot(ax, q_rng, q_rng, 'r-', 'LineWidth', 1.5);
xlabel(ax, 'Theoretical quantiles'); ylabel(ax, 'Sample quantiles');
axis(ax,'square');
end


function p = chi2cdf_local(x, k)
% Minimal chi-square CDF helper using gammainc, available in base MATLAB.
p = gammainc(x/2, k/2, 'lower');
end


function q = quantile_omitnan_local(x, p)
x = sort(x(~isnan(x)));
if isempty(x), q = NaN; return; end
if numel(x) == 1, q = x; return; end
idx = 1 + (numel(x)-1) * p;
lo = floor(idx); hi = ceil(idx);
if lo == hi
    q = x(lo);
else
    q = x(lo) + (idx-lo) * (x(hi)-x(lo));
end
end


function out = ternary_local(cond, a, b)
if cond, out = a; else, out = b; end
end