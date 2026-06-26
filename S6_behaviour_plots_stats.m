% ==========================================================================
% S6_behaviour_plots_stats.m  (was: B_plot_task_stage_behaviour.m)
%
% PIPELINE STEP 6 of 7 — behavioural figures + statistics.
%
% CHANGES v2:
%   - Individual subject points are now CONNECTED (spaghetti lines) across
%     stages in all 4-stage panels, and paired between conditions (D vs P,
%     win vs lose, correct vs incorrect) in all 2-condition panels.
%   - All axes are square (axis square applied everywhere).
%   - Key statistics (t-tests / Wilcoxon signed-rank) are printed for all
%     pairwise comparisons, alongside the LME tables.
%   - Normality is formally tested before every t-test using Shapiro-Wilk
%     (swtest, if available on path) or Lilliefors (lillietest, Statistics
%     Toolbox) or standardised KS (fallback).  The test used, its p-value,
%     and the resulting decision (t-test vs Wilcoxon) are printed.
%
% TASK STAGES
%   LN — Learning Naive    : first 20 trials of block
%   LE — Learning Expert   : 20 trials before reversal
%   RN — Reversal Naive    : 20 trials after reversal
%   RE — Reversal Expert   : last 20 trials of block
%
% FIGURES
%   Fig S1  per-stage accuracy + response rates
%   Fig S2  true-FB win/loss rates (P blocks)
%   Fig S3  confidence per stage x block type
%   Fig S4  RT per stage x block type
%   Fig S5  stay behaviour per stage
%   Fig S6  stage profiles by transition type
%   Fig S7  heatmap summary
%   Fig S8  reversal dynamics by experienced uncertainty
% ==========================================================================

%% ─── SETTINGS ─────────────────────────────────────────────────────────────
addpath(genpath(fileparts(mfilename('fullpath'))));
set(groot,'defaultAxesTickDir','out');
set(groot,'defaultAxesBox','off');

STAGE_WIN  = 20;
STAGES     = {'LN','LE','RN','RE'};
STAGE_LBLS = {'LN (block start)','LE (pre-rev)','RN (post-rev)','RE (block end)'};

CLR_D    = [0.15 0.45 0.70];
CLR_P    = [0.80 0.30 0.10];
CLR_STGS = {[0.20 0.63 0.17],[0.12 0.47 0.71],[0.85 0.33 0.10],[0.58 0.40 0.74]};

TRANS_TYPES  = {'D→D','D→P','P→D','P→P'};
TRANS_COLORS = {[0.12 0.47 0.71],[0.85 0.33 0.10],[0.47 0.67 0.19],[0.80 0.20 0.60]};


%% ─── PATHS ─────────────────────────────────────────────────────────────────
remote = 0;
switch remote
    case 1, base_path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch/Salient mod switch KH';
    case 0, base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH';
    case 2, base_path = 'Z:\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH';
end
outpath = fullfile(base_path,'Results','Behav results','Stage Figures');
cd(base_path)

beh_candidates = { ...
    fullfile(base_path,'Data','behav_table_RL.mat'), ...
    fullfile(base_path,'Data','behav_table_June2026_RL.mat'), ...
    fullfile(base_path,'Data','behav_table_June2026.mat'), ...
    fullfile(base_path,'Data','behav_table.mat')};
loaded_beh = false;
for ci = 1:numel(beh_candidates)
    if exist(beh_candidates{ci},'file')
        load(beh_candidates{ci},'group_T');
        fprintf('Loaded behaviour table: %s\n',beh_candidates{ci});
        loaded_beh = true; break;
    end
end
if ~loaded_beh; error('No behaviour table found in %s/Data.',base_path); end

%% ─── PREPROCESS group_T ───────────────────────────────────────────────────
group_T.trial            = double(group_T.trial);
group_T.block            = double(group_T.block);
group_T.revTrial         = double(group_T.revTrial);
group_T.correct          = double(group_T.correct);
group_T.confidence       = double(group_T.confidence);
group_T.RT               = double(group_T.RT);
group_T.RT(group_T.RT > 1) = NaN;
group_T.trueFB           = double(group_T.trueFB);
group_T.goTrial          = double(group_T.goTrial);
group_T.perceivedCorrect = double(group_T.perceivedCorrect);
group_T.block_type       = string(group_T.block_type);
group_T.prev_block_type  = string(group_T.prev_block_type);
group_T.subjID           = string(group_T.subjID);

group_T.block_type_clean = group_T.block_type;
group_T.block_type_clean(group_T.block_type == "V") = "P";

group_T.transition = group_T.prev_block_type + "→" + group_T.block_type_clean;
group_T.transition(group_T.prev_block_type == "NaN" | ...
                   group_T.prev_block_type == "") = "first";

group_T.isHit  = double(group_T.correct==1 & group_T.goTrial==1);
group_T.isFA   = double(group_T.correct==0 & group_T.goTrial==0);
group_T.isMiss = double(group_T.correct==0 & group_T.goTrial==1);
group_T.isCR   = double(group_T.correct==1 & group_T.goTrial==0);

group_T.trueFB_win  = double(group_T.trueFB==1 & group_T.perceivedCorrect==1);
group_T.trueFB_loss = double(group_T.trueFB==1 & group_T.perceivedCorrect==0);
group_T.falseFB     = double(group_T.trueFB==0);


%% ─── ASSIGN STAGE PER TRIAL ─────────────────────────────────────────────────
nRows = height(group_T);
group_T.stage = repmat("",nRows,1);

subj_block_pairs = unique(group_T(:,{'subjID','block'}),'rows');

for pb = 1:height(subj_block_pairs)
    sn  = subj_block_pairs.subjID(pb);
    blk = subj_block_pairs.block(pb);
    mask = group_T.subjID==sn & group_T.block==blk;
    rows = find(mask);
    if isempty(rows), continue; end
    rev_t  = group_T.revTrial(rows(1));
    trials = group_T.trial(rows);
    max_t  = max(trials);

    if isnan(rev_t) || ~isfinite(rev_t)
        for ri = 1:numel(rows)
            t = trials(ri);
            if     t <= STAGE_WIN,             group_T.stage(rows(ri)) = "LN";
            elseif t > max_t - STAGE_WIN,      group_T.stage(rows(ri)) = "RE";
            end
        end
        continue;
    end
    rev_t = round(rev_t);
    LN_start = 1;              LN_end = min(STAGE_WIN,rev_t-1);
    LE_start = max(1,rev_t-STAGE_WIN); LE_end = rev_t-1;
    RN_start = rev_t;          RN_end = min(max_t,rev_t+STAGE_WIN-1);
    RE_start = max(1,max_t-STAGE_WIN+1); RE_end = max_t;

    for ri = 1:numel(rows)
        t = trials(ri);
        if     t>=RE_start && t<=RE_end, group_T.stage(rows(ri)) = "RE";
        elseif t>=RN_start && t<=RN_end, group_T.stage(rows(ri)) = "RN";
        elseif t>=LE_start && t<=LE_end, group_T.stage(rows(ri)) = "LE";
        elseif t>=LN_start && t<=LN_end, group_T.stage(rows(ri)) = "LN";
        end
    end
end

T = group_T(group_T.stage ~= "",:);
fprintf('Valid stage labels: %d / %d rows\n',height(T),nRows);

%% ─── BUILD SUBJECT × STAGE MATRICES ──────────────────────────────────────
all_subjs = unique(T.subjID);
N_subj    = numel(all_subjs);
n_stg     = numel(STAGES);
x_pos     = 1:n_stg;

measures = {'correct','isHit','isFA','isMiss','isCR', ...
            'confidence','RT','trueFB_win','trueFB_loss','falseFB'};

SMAT = struct();
for mi = 1:numel(measures)
    SMAT.(measures{mi}).all = NaN(N_subj,n_stg);
    SMAT.(measures{mi}).D   = NaN(N_subj,n_stg);
    SMAT.(measures{mi}).P   = NaN(N_subj,n_stg);
end

for si = 1:N_subj
    sn = all_subjs(si);
    Ts = T(T.subjID==sn,:);
    for stgi = 1:n_stg
        sg    = STAGES{stgi};
        Tsg   = Ts(Ts.stage==sg,:);
        Tsg_D = Tsg(Tsg.block_type_clean=="D",:);
        Tsg_P = Tsg(Tsg.block_type_clean=="P",:);
        for mi = 1:numel(measures)
            m = measures{mi};
            if height(Tsg)   > 0, SMAT.(m).all(si,stgi) = mean(Tsg.(m),  'omitnan'); end
            if height(Tsg_D) > 0, SMAT.(m).D(si,stgi)   = mean(Tsg_D.(m),'omitnan'); end
            if height(Tsg_P) > 0, SMAT.(m).P(si,stgi)   = mean(Tsg_P.(m),'omitnan'); end
        end
    end
end
fprintf('Stage matrices built: %d subjects x %d stages.\n',N_subj,n_stg);


% =========================================================================
%% FIG S1 — PER-STAGE ACCURACY AND RESPONSE RATES
% =========================================================================
fig_s1 = figure('Position',[50 50 1500 650]);
sgtitle('Task stage performance: accuracy and response rates (subject means +/- SEM)','FontSize',12);

meas_s1 = {'correct','isHit','isFA','isMiss','isCR'};
ylbl_s1 = {'P(correct)','P(Hit)','P(FA)','P(Miss)','P(CR)'};
ylim_s1 = {[0 1],[0 1],[0 1],[0 1],[0 1]};

for mi = 1:numel(meas_s1)
    m = meas_s1{mi};

    % ── Top row: all blocks pooled ────────────────────────────────────
    ax = subplot(2,numel(meas_s1),mi); hold(ax,'on');
    title(ax,ylbl_s1{mi},'FontSize',10);

    mat_all = SMAT.(m).all;
    gp_mean = mean(mat_all,1,'omitnan');
    gp_se   = stage_sem(mat_all);
    n_con   = sum(~isnan(mat_all),1);

    for stgi = 1:n_stg
        bar(ax,stgi,gp_mean(stgi),0.55,'FaceColor',CLR_STGS{stgi}, ...
            'EdgeColor','none','FaceAlpha',0.75);
    end
    errorbar(ax,x_pos,gp_mean,gp_se,'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');

    % CONNECTED spaghetti: one line per subject across the 4 stages
    plot_spaghetti(ax, mat_all, x_pos);

    for stgi = 1:n_stg
        text(ax,stgi,ylim_s1{mi}(2)*0.03,sprintf('n=%d',n_con(stgi)), ...
            'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
    end
    set(ax,'XTick',x_pos,'XTickLabel',STAGES,'FontSize',9);
    ylabel(ax,ylbl_s1{mi}); ylim(ax,ylim_s1{mi});
    axis(ax,'square');
    add_pairwise_brackets(ax,mat_all,x_pos,ylim_s1{mi}(2));

    % ── Bottom row: D vs P ────────────────────────────────────────────
    ax2 = subplot(2,numel(meas_s1),mi+numel(meas_s1)); hold(ax2,'on');
    title(ax2,sprintf('%s: D vs P',ylbl_s1{mi}),'FontSize',9);

    bt_data  = {SMAT.(m).D,  SMAT.(m).P};
    bt_clrs  = {CLR_D,       CLR_P};
    bt_lbls  = {'Det','Prob'};
    bt_xoffs = [-0.18, 0.18];

    for bti = 1:2
        mat_bt = bt_data{bti};
        gm  = mean(mat_bt,1,'omitnan');
        gs  = stage_sem(mat_bt);
        n_b = sum(~isnan(mat_bt),1);
        xp  = x_pos + bt_xoffs(bti);
        bar(ax2,xp,gm,0.3,'FaceColor',bt_clrs{bti},'EdgeColor','none', ...
            'FaceAlpha',0.75,'DisplayName',bt_lbls{bti});
        errorbar(ax2,xp,gm,gs,'k.','LineWidth',1.2,'CapSize',4,'HandleVisibility','off');
        for stgi = 1:n_stg
            text(ax2,xp(stgi),ylim_s1{mi}(2)*0.03,sprintf('n=%d',n_b(stgi)), ...
                'HorizontalAlignment','center','FontSize',6,'Color',bt_clrs{bti}*0.7);
        end
    end
    % CONNECTED paired dots: D↔P per subject at each stage
    plot_paired_dots(ax2, SMAT.(m).D, SMAT.(m).P, x_pos, bt_xoffs(1), bt_xoffs(2), ...
        CLR_D, CLR_P, 1);

    set(ax2,'XTick',x_pos,'XTickLabel',STAGES,'FontSize',9);
    ylabel(ax2,ylbl_s1{mi}); ylim(ax2,ylim_s1{mi});
    axis(ax2,'square');
    legend(ax2,'Box','off','FontSize',8,'Location','best');
end

saveas(fig_s1,fullfile(outpath,'figS1_stage_accuracy.pdf'));
fprintf('Fig S1 saved.\n');


% =========================================================================
%% FIG S2 — TRUE-FB WIN / LOSS / FALSE-FB RATES (P blocks only)
% =========================================================================
fig_s2 = figure('Position',[50 50 1100 500]);
sgtitle({'True-FB win/loss and false-feedback rate per stage', ...
    '(P blocks only)'},'FontSize',12);

meas_s2  = {'trueFB_win','trueFB_loss','falseFB'};
ylbl_s2  = {'P(trueFB win)','P(trueFB loss)','P(false feedback)'};
y_expect = {0.64,0.16,0.20};

for mi = 1:numel(meas_s2)
    m   = meas_s2{mi};
    ax  = subplot(1,3,mi); hold(ax,'on');
    title(ax,ylbl_s2{mi},'FontSize',10);

    mat_P = SMAT.(m).P;
    gm    = mean(mat_P,1,'omitnan');
    gs    = stage_sem(mat_P);
    n_c   = sum(~isnan(mat_P),1);

    for stgi = 1:n_stg
        bar(ax,stgi,gm(stgi),0.55,'FaceColor',CLR_STGS{stgi},'EdgeColor','none','FaceAlpha',0.75);
    end
    errorbar(ax,x_pos,gm,gs,'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');
    yline(ax,y_expect{mi},'k--','LineWidth',1.2,'HandleVisibility','off');
    text(ax,0.02,y_expect{mi}+0.02,sprintf('Expected=%.2f',y_expect{mi}), ...
        'FontSize',7,'Color',[0.3 0.3 0.3]);

    % CONNECTED spaghetti
    plot_spaghetti(ax, mat_P, x_pos);

    for stgi = 1:n_stg
        text(ax,stgi,0.02,sprintf('n=%d',n_c(stgi)), ...
            'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
    end
    set(ax,'XTick',x_pos,'XTickLabel',STAGES,'FontSize',9);
    ylabel(ax,ylbl_s2{mi}); ylim(ax,[0 1]);
    axis(ax,'square');
    add_pairwise_brackets(ax,mat_P,x_pos,0.95);
    subtitle(ax,'P blocks only','FontSize',8,'Color',[0.5 0.5 0.5]);
end

annotation('textbox',[0.01 0.01 0.98 0.07],'String', ...
    ['trueFB win = honest feedback shown as correct. ' ...
    'trueFB loss = honest feedback shown as incorrect. ' ...
    'false feedback = misleading trial (unknown to participant). ' ...
    'Expected at p(trueFB)=0.8: win=0.64, loss=0.16, falseFB=0.20.'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig_s2,fullfile(outpath,'figS2_stage_trueFB.pdf'));
fprintf('Fig S2 saved.\n');


% =========================================================================
%% FIG S3 — CONFIDENCE PER STAGE × BLOCK TYPE
% =========================================================================
fig_s3 = figure('Position',[50 50 1300 550]);
sgtitle('Confidence per stage: overall and D vs P blocks','FontSize',12);

% Panel 1: all blocks pooled
ax_s3a = subplot(1,3,1); hold(ax_s3a,'on');
title(ax_s3a,'Confidence: all blocks','FontSize',10);
mat_all_c = SMAT.confidence.all;
gm = mean(mat_all_c,1,'omitnan');
gs = stage_sem(mat_all_c);
n_c = sum(~isnan(mat_all_c),1);
for stgi = 1:n_stg
    bar(ax_s3a,stgi,gm(stgi),0.55,'FaceColor',CLR_STGS{stgi},'EdgeColor','none','FaceAlpha',0.75);
end
errorbar(ax_s3a,x_pos,gm,gs,'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');

% CONNECTED spaghetti
plot_spaghetti(ax_s3a, mat_all_c, x_pos);

for stgi = 1:n_stg
    text(ax_s3a,stgi,1.2,sprintf('n=%d',n_c(stgi)), ...
        'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
end
set(ax_s3a,'XTick',x_pos,'XTickLabel',STAGES,'FontSize',9);
ylabel(ax_s3a,'Confidence (1-10)'); ylim(ax_s3a,[1 10]);
axis(ax_s3a,'square');
add_pairwise_brackets(ax_s3a,mat_all_c,x_pos,9.5);

% Panel 2: D vs P bars + connected paired dots
ax_s3b = subplot(1,3,2); hold(ax_s3b,'on');
title(ax_s3b,'Confidence: D vs P blocks','FontSize',10);

bt_data  = {SMAT.confidence.D, SMAT.confidence.P};
bt_clrs  = {CLR_D,              CLR_P};
bt_lbls  = {'Det','Prob'};
bt_xoffs = [-0.18, 0.18];

for bti = 1:2
    mat_bt = bt_data{bti};
    gm2 = mean(mat_bt,1,'omitnan');
    gs2 = stage_sem(mat_bt);
    xp  = x_pos + bt_xoffs(bti);
    bar(ax_s3b,xp,gm2,0.3,'FaceColor',bt_clrs{bti},'EdgeColor','none', ...
        'FaceAlpha',0.75,'DisplayName',bt_lbls{bti});
    errorbar(ax_s3b,xp,gm2,gs2,'k.','LineWidth',1.2,'CapSize',4,'HandleVisibility','off');
end
% CONNECTED paired dots: D↔P per subject at each stage
plot_paired_dots(ax_s3b, SMAT.confidence.D, SMAT.confidence.P, x_pos, ...
    bt_xoffs(1), bt_xoffs(2), CLR_D, CLR_P, 1);

set(ax_s3b,'XTick',x_pos,'XTickLabel',STAGES,'FontSize',9);
ylabel(ax_s3b,'Confidence (1-10)'); ylim(ax_s3b,[1 10]);
axis(ax_s3b,'square');
legend(ax_s3b,'Box','off','FontSize',9);

% Panel 3: spaghetti by block type (group means + per-subject lines)
ax_s3c = subplot(1,3,3); hold(ax_s3c,'on');
title(ax_s3c,'Stage means per subject','FontSize',10);
subtitle(ax_s3c,'Lines = within-subject stage means','FontSize',8,'Color',[0.5 0.5 0.5]);
for si = 1:N_subj
    row_D = SMAT.confidence.D(si,:);
    row_P = SMAT.confidence.P(si,:);
    if ~all(isnan(row_D))
        hl = plot(ax_s3c,x_pos,row_D,'o-','Color',CLR_D, ...
            'LineWidth',0.8,'MarkerSize',3,'MarkerFaceColor',CLR_D,'HandleVisibility','off');
        hl.Color(4) = 0.25;
    end
    if ~all(isnan(row_P))
        hl = plot(ax_s3c,x_pos,row_P,'s--','Color',CLR_P, ...
            'LineWidth',0.8,'MarkerSize',3,'MarkerFaceColor',CLR_P,'HandleVisibility','off');
        hl.Color(4) = 0.25;
    end
end
plot(ax_s3c,x_pos,mean(SMAT.confidence.D,1,'omitnan'),'o-', ...
    'Color',CLR_D,'LineWidth',2.5,'MarkerSize',8,'MarkerFaceColor',CLR_D,'DisplayName','Det');
plot(ax_s3c,x_pos,mean(SMAT.confidence.P,1,'omitnan'),'s--', ...
    'Color',CLR_P,'LineWidth',2.5,'MarkerSize',8,'MarkerFaceColor',CLR_P,'DisplayName','Prob');
set(ax_s3c,'XTick',x_pos,'XTickLabel',STAGES,'FontSize',9);
ylabel(ax_s3c,'Confidence (1-10)'); ylim(ax_s3c,[1 10]);
axis(ax_s3c,'square');
legend(ax_s3c,'Box','off','FontSize',9,'Location','best');

saveas(fig_s3,fullfile(outpath,'figS3_stage_confidence.pdf'));
fprintf('Fig S3 saved.\n');


% =========================================================================
%% FIG S4 — REACTION TIME PER STAGE × BLOCK TYPE
% =========================================================================
fig_s4 = figure('Position',[50 50 1300 500]);
sgtitle('Reaction time per stage: overall and D vs P blocks','FontSize',12);

% Panel 1: RT all blocks
ax_s4a = subplot(1,3,1); hold(ax_s4a,'on');
title(ax_s4a,'RT: all blocks','FontSize',10);
mat_rt = SMAT.RT.all;
gm = mean(mat_rt,1,'omitnan') * 1000;
gs = stage_sem(mat_rt) * 1000;
n_c = sum(~isnan(mat_rt),1);
for stgi = 1:n_stg
    bar(ax_s4a,stgi,gm(stgi),0.55,'FaceColor',CLR_STGS{stgi},'EdgeColor','none','FaceAlpha',0.75);
end
errorbar(ax_s4a,x_pos,gm,gs,'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');

% CONNECTED spaghetti (scale to ms)
plot_spaghetti(ax_s4a, mat_rt*1000, x_pos);

for stgi = 1:n_stg
    if ~isnan(gm(stgi))
        text(ax_s4a,stgi,min(gm(~isnan(gm)))*0.85,sprintf('n=%d',n_c(stgi)), ...
            'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
    end
end
set(ax_s4a,'XTick',x_pos,'XTickLabel',STAGES,'FontSize',9);
ylabel(ax_s4a,'RT (ms)');
if any(~isnan(gm))
    ylim(ax_s4a,[0 max(gm+gs,[],'omitnan')*1.3]);
    add_pairwise_brackets(ax_s4a,mat_rt*1000,x_pos,max(gm+gs,[],'omitnan')*1.25);
end
axis(ax_s4a,'square');

% Panel 2: D vs P + connected paired dots
ax_s4b = subplot(1,3,2); hold(ax_s4b,'on');
title(ax_s4b,'RT: D vs P blocks','FontSize',10);

bt_data  = {SMAT.RT.D, SMAT.RT.P};
bt_clrs  = {CLR_D,      CLR_P};
bt_lbls  = {'Det','Prob'};
bt_xoffs = [-0.18, 0.18];

for bti = 1:2
    mat_bt = bt_data{bti};
    gm2 = mean(mat_bt,1,'omitnan') * 1000;
    gs2 = stage_sem(mat_bt) * 1000;
    xp  = x_pos + bt_xoffs(bti);
    bar(ax_s4b,xp,gm2,0.3,'FaceColor',bt_clrs{bti},'EdgeColor','none', ...
        'FaceAlpha',0.75,'DisplayName',bt_lbls{bti});
    errorbar(ax_s4b,xp,gm2,gs2,'k.','LineWidth',1.2,'CapSize',4,'HandleVisibility','off');
end
% CONNECTED paired dots (scale to ms)
plot_paired_dots(ax_s4b, SMAT.RT.D, SMAT.RT.P, x_pos, bt_xoffs(1), bt_xoffs(2), ...
    CLR_D, CLR_P, 1000);

set(ax_s4b,'XTick',x_pos,'XTickLabel',STAGES,'FontSize',9);
ylabel(ax_s4b,'RT (ms)');
axis(ax_s4b,'square');
legend(ax_s4b,'Box','off','FontSize',9);

% Panel 3: correct vs incorrect RT + connected paired dots
ax_s4c = subplot(1,3,3); hold(ax_s4c,'on');
title(ax_s4c,'RT: correct vs incorrect','FontSize',10);

SMAT_RT_corr   = NaN(N_subj,n_stg);
SMAT_RT_incorr = NaN(N_subj,n_stg);
for si = 1:N_subj
    sn = all_subjs(si);
    Ts = T(T.subjID==sn,:);
    for stgi = 1:n_stg
        sg  = STAGES{stgi};
        Tsg = Ts(Ts.stage==sg,:);
        SMAT_RT_corr(si,stgi)   = mean(Tsg.RT(Tsg.correct==1),'omitnan') * 1000;
        SMAT_RT_incorr(si,stgi) = mean(Tsg.RT(Tsg.correct==0),'omitnan') * 1000;
    end
end

bt_data  = {SMAT_RT_corr,   SMAT_RT_incorr};
bt_clrs  = {[0.2 0.6 0.2],  [0.8 0.2 0.2]};
bt_lbls  = {'Correct','Incorrect'};
bt_xoffs = [-0.18, 0.18];

for bti = 1:2
    mat_bt = bt_data{bti};
    gm2 = mean(mat_bt,1,'omitnan');
    gs2 = stage_sem(mat_bt);
    xp  = x_pos + bt_xoffs(bti);
    bar(ax_s4c,xp,gm2,0.3,'FaceColor',bt_clrs{bti},'EdgeColor','none', ...
        'FaceAlpha',0.75,'DisplayName',bt_lbls{bti});
    errorbar(ax_s4c,xp,gm2,gs2,'k.','LineWidth',1.2,'CapSize',4,'HandleVisibility','off');
end
% CONNECTED paired dots: correct↔incorrect per subject at each stage
% (data already in ms)
plot_paired_dots(ax_s4c, SMAT_RT_corr, SMAT_RT_incorr, x_pos, bt_xoffs(1), bt_xoffs(2), ...
    bt_clrs{1}, bt_clrs{2}, 1);

set(ax_s4c,'XTick',x_pos,'XTickLabel',STAGES,'FontSize',9);
ylabel(ax_s4c,'RT (ms)');
axis(ax_s4c,'square');
legend(ax_s4c,'Box','off','FontSize',9,'Location','best');

saveas(fig_s4,fullfile(outpath,'figS4_stage_RT.pdf'));
fprintf('Fig S4 saved.\n');


% =========================================================================
%% FIG S5 — STAY BEHAVIOUR PER STAGE
% =========================================================================
has_stay = ismember('stay_choice',T.Properties.VariableNames) && ...
           ismember('prevCorrect',T.Properties.VariableNames);

if has_stay
    T_stay = T(~isnan(T.stay_choice) & ~isnan(T.prevCorrect),:);
else
    T_stay = T(false(height(T),1),:);
    warning('plot_task_stages: stay_choice / prevCorrect not found. Fig S5 will be empty.');
end

SMAT_winstay  = NaN(N_subj,n_stg);
SMAT_losestay = NaN(N_subj,n_stg);
SMAT_stay_all = NaN(N_subj,n_stg);
SMAT_stay_D   = NaN(N_subj,n_stg);
SMAT_stay_P   = NaN(N_subj,n_stg);

if has_stay
    for si = 1:N_subj
        sn = all_subjs(si);
        Ts_stay = T_stay(T_stay.subjID==sn,:);
        for stgi = 1:n_stg
            sg    = STAGES{stgi};
            Tsg   = Ts_stay(Ts_stay.stage==sg,:);
            if height(Tsg)==0, continue; end
            SMAT_stay_all(si,stgi) = mean(Tsg.stay_choice,'omitnan');
            win_r  = Tsg(Tsg.prevCorrect==1,:);
            lose_r = Tsg(Tsg.prevCorrect==0,:);
            if height(win_r)  > 0, SMAT_winstay(si,stgi)  = mean(win_r.stay_choice,'omitnan');  end
            if height(lose_r) > 0, SMAT_losestay(si,stgi) = mean(lose_r.stay_choice,'omitnan'); end
            Tsg_D = Tsg(Tsg.block_type_clean=="D",:);
            Tsg_P = Tsg(Tsg.block_type_clean=="P",:);
            if height(Tsg_D) > 0, SMAT_stay_D(si,stgi) = mean(Tsg_D.stay_choice,'omitnan'); end
            if height(Tsg_P) > 0, SMAT_stay_P(si,stgi) = mean(Tsg_P.stay_choice,'omitnan'); end
        end
    end
end

fig_s5 = figure('Position',[50 50 1300 500]);
sgtitle('Stay behaviour per stage (stimulus-specific)','FontSize',12);

% Panel 1: P(stay) all — CONNECTED spaghetti
ax_s5a = subplot(1,3,1); hold(ax_s5a,'on');
title(ax_s5a,'P(stay): all trials','FontSize',10);
mat5 = SMAT_stay_all;
gm5 = mean(mat5,1,'omitnan'); gs5 = stage_sem(mat5); n5 = sum(~isnan(mat5),1);
for stgi = 1:n_stg
    bar(ax_s5a,stgi,gm5(stgi),0.55,'FaceColor',CLR_STGS{stgi},'EdgeColor','none','FaceAlpha',0.75);
end
errorbar(ax_s5a,x_pos,gm5,gs5,'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');
plot_spaghetti(ax_s5a, mat5, x_pos);
for stgi = 1:n_stg
    text(ax_s5a,stgi,0.02,sprintf('n=%d',n5(stgi)), ...
        'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
end
if ~has_stay
    text(ax_s5a,0.5,0.5,'stay\_choice not in group\_T', ...
        'Units','normalized','HorizontalAlignment','center','FontSize',9,'Color',[0.6 0.3 0]);
end
yline(ax_s5a,0.5,'k:','HandleVisibility','off');
set(ax_s5a,'XTick',x_pos,'XTickLabel',STAGES);
ylabel(ax_s5a,'P(stay)'); ylim(ax_s5a,[0 1]);
axis(ax_s5a,'square');
add_pairwise_brackets(ax_s5a,mat5,x_pos,0.95);

% Panel 2: win-stay vs lose-stay + connected paired dots
ax_s5b = subplot(1,3,2); hold(ax_s5b,'on');
title(ax_s5b,'Win-stay vs lose-stay','FontSize',10);

bt_data  = {SMAT_winstay,   SMAT_losestay};
bt_clrs  = {[0.2 0.6 0.2],  [0.8 0.2 0.2]};
bt_lbls  = {'Win-stay','Lose-stay'};
bt_xoffs = [-0.18, 0.18];

for bti = 1:2
    mat_bt = bt_data{bti};
    gm2 = mean(mat_bt,1,'omitnan');
    gs2 = stage_sem(mat_bt);
    xp  = x_pos + bt_xoffs(bti);
    bar(ax_s5b,xp,gm2,0.3,'FaceColor',bt_clrs{bti},'EdgeColor','none', ...
        'FaceAlpha',0.75,'DisplayName',bt_lbls{bti});
    errorbar(ax_s5b,xp,gm2,gs2,'k.','LineWidth',1.2,'CapSize',4,'HandleVisibility','off');
end
% CONNECTED paired dots: win↔lose per subject at each stage
plot_paired_dots(ax_s5b, SMAT_winstay, SMAT_losestay, x_pos, bt_xoffs(1), bt_xoffs(2), ...
    bt_clrs{1}, bt_clrs{2}, 1);

yline(ax_s5b,0.5,'k:','HandleVisibility','off');
set(ax_s5b,'XTick',x_pos,'XTickLabel',STAGES);
ylabel(ax_s5b,'P(stay)'); ylim(ax_s5b,[0 1]);
axis(ax_s5b,'square');
legend(ax_s5b,'Box','off','FontSize',9,'Location','best');

% Panel 3: D vs P stay + connected paired dots
ax_s5c = subplot(1,3,3); hold(ax_s5c,'on');
title(ax_s5c,'Stay: D vs P blocks','FontSize',10);

bt_data  = {SMAT_stay_D, SMAT_stay_P};
bt_clrs  = {CLR_D,        CLR_P};
bt_lbls  = {'Det','Prob'};
bt_xoffs = [-0.18, 0.18];

for bti = 1:2
    mat_bt = bt_data{bti};
    gm2 = mean(mat_bt,1,'omitnan');
    gs2 = stage_sem(mat_bt);
    xp  = x_pos + bt_xoffs(bti);
    bar(ax_s5c,xp,gm2,0.3,'FaceColor',bt_clrs{bti},'EdgeColor','none', ...
        'FaceAlpha',0.75,'DisplayName',bt_lbls{bti});
    errorbar(ax_s5c,xp,gm2,gs2,'k.','LineWidth',1.2,'CapSize',4,'HandleVisibility','off');
end
% CONNECTED paired dots: D↔P stay per subject at each stage
plot_paired_dots(ax_s5c, SMAT_stay_D, SMAT_stay_P, x_pos, bt_xoffs(1), bt_xoffs(2), ...
    CLR_D, CLR_P, 1);

yline(ax_s5c,0.5,'k:','HandleVisibility','off');
set(ax_s5c,'XTick',x_pos,'XTickLabel',STAGES);
ylabel(ax_s5c,'P(stay)'); ylim(ax_s5c,[0 1]);
axis(ax_s5c,'square');
legend(ax_s5c,'Box','off','FontSize',9);

saveas(fig_s5,fullfile(outpath,'figS5_stage_stay.pdf'));
fprintf('Fig S5 saved.\n');


% =========================================================================
%% FIG S6 — STAGE PROFILES BY TRANSITION TYPE
% =========================================================================
fig_s6 = figure('Position',[50 50 1600 700]);
sgtitle('Stage profiles by block transition type','FontSize',12);

meas_s6  = {'correct','confidence'};
ylbl_s6  = {'P(correct)','Confidence (1-10)'};
ylims_s6 = {[0 1],[1 10]};

trans_subjs_mat = cell(numel(TRANS_TYPES),2);
for ti = 1:numel(TRANS_TYPES)
    for mi = 1:2
        trans_subjs_mat{ti,mi} = NaN(N_subj,n_stg);
    end
    tr_str = TRANS_TYPES{ti};
    for si = 1:N_subj
        sn = all_subjs(si);
        Ts = T(T.subjID==sn & T.transition==tr_str,:);
        if height(Ts)==0, continue; end
        for stgi = 1:n_stg
            sg  = STAGES{stgi};
            Tsg = Ts(Ts.stage==sg,:);
            if height(Tsg)==0, continue; end
            trans_subjs_mat{ti,1}(si,stgi) = mean(Tsg.correct,   'omitnan');
            trans_subjs_mat{ti,2}(si,stgi) = mean(Tsg.confidence,'omitnan');
        end
    end
end

le_idx = find(strcmp(STAGES,'LE'));
rn_idx = find(strcmp(STAGES,'RN'));

for mi = 1:2
    for ti = 1:numel(TRANS_TYPES)
        ax = subplot(2,numel(TRANS_TYPES),(mi-1)*numel(TRANS_TYPES)+ti);
        hold(ax,'on');
        if mi==1, title(ax,TRANS_TYPES{ti},'FontSize',11); end
        if ti==1, ylabel(ax,ylbl_s6{mi}); end

        mat = trans_subjs_mat{ti,mi};
        gm  = mean(mat,1,'omitnan');
        gs  = stage_sem(mat);
        n_c = sum(~isnan(mat),1);

        % Individual spaghetti (light) — already connected per original S6
        for si = 1:N_subj
            if sum(~isnan(mat(si,:))) >= 2
                hl = plot(ax,x_pos,mat(si,:),'o-', ...
                    'Color',TRANS_COLORS{ti},'LineWidth',0.8,'MarkerSize',3,'HandleVisibility','off');
                hl.Color(4) = 0.15;
            end
        end
        % Group mean ribbon
        fill(ax,[x_pos,fliplr(x_pos)],[gm+gs,fliplr(gm-gs)], ...
            TRANS_COLORS{ti},'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
        plot(ax,x_pos,gm,'o-','Color',TRANS_COLORS{ti},'LineWidth',2.5,'MarkerSize',8, ...
            'MarkerFaceColor',TRANS_COLORS{ti},'DisplayName',sprintf('n=%d',max(n_c)));

        % LE→RN connector
        if ~isnan(gm(le_idx)) && ~isnan(gm(rn_idx))
            plot(ax,[le_idx rn_idx],[gm(le_idx) gm(rn_idx)],'k-','LineWidth',2,'HandleVisibility','off');
        end

        set(ax,'XTick',x_pos,'XTickLabel',STAGE_LBLS,'XTickLabelRotation',20,'FontSize',8);
        ylim(ax,ylims_s6{mi});
        axis(ax,'square');
        for stgi = 1:n_stg
            text(ax,stgi,ylims_s6{mi}(1)+diff(ylims_s6{mi})*0.03, ...
                sprintf('n=%d',n_c(stgi)),'HorizontalAlignment','center', ...
                'FontSize',7,'Color',TRANS_COLORS{ti}*0.7);
        end
        legend(ax,'Box','off','FontSize',8,'Location','best');
    end
end

annotation('textbox',[0.01 0.01 0.98 0.04],'String', ...
    ['Bold line connects LE to RN to highlight reversal cost per transition. ' ...
    'n = subjects contributing to that transition x stage cell.'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig_s6,fullfile(outpath,'figS6_stage_by_transition.pdf'));
fprintf('Fig S6 saved.\n');


% =========================================================================
%% FIG S7 — HEATMAP SUMMARY
% =========================================================================
fig_s7 = figure('Position',[50 50 1400 600]);
sgtitle('Heatmap: group-mean measures x stage x block type','FontSize',12);

hmap_measures = {'correct','isHit','isFA','isMiss','isCR','confidence','RT'};
hmap_lbls     = {'Accuracy','Hit','FA','Miss','CR','Confidence','RT (s)'};
n_hm = numel(hmap_measures);

for bti = 1:2
    bt_tag = ternary_local(bti==1,'D','P');
    bt_lbl = ternary_local(bti==1,'Deterministic','Probabilistic');
    ax_hm  = subplot(1,2,bti);

    hmat = NaN(n_hm,n_stg);
    for mi = 1:n_hm
        m = hmap_measures{mi};
        if isfield(SMAT,m)
            col = mean(SMAT.(m).(bt_tag),1,'omitnan');
            if strcmp(m,'RT')
                r = range(col); if r==0, r=1; end
                col = (col - min(col)) / r;
            end
            hmat(mi,:) = col;
        end
    end

    imagesc(ax_hm,hmat);
    colormap(ax_hm);
    cb = colorbar(ax_hm); cb.Label.String = 'Mean (RT normalised 0-1)';
    set(ax_hm,'XTick',1:n_stg,'XTickLabel',STAGES, ...
              'YTick',1:n_hm,'YTickLabel',hmap_lbls,'FontSize',10);
    title(ax_hm,bt_lbl,'FontSize',11);
    xlabel(ax_hm,'Stage'); ylabel(ax_hm,'Measure');
    axis(ax_hm,'square');

    for mi = 1:n_hm
        for stgi = 1:n_stg
            v = mean(SMAT.(hmap_measures{mi}).(bt_tag)(:,stgi),'omitnan');
            if ~isnan(v)
                text(ax_hm,stgi,mi,sprintf('%.2f',v), ...
                    'HorizontalAlignment','center','VerticalAlignment','middle', ...
                    'FontSize',8,'Color','w','FontWeight','bold');
            end
        end
    end
end

saveas(fig_s7,fullfile(outpath,'figS7_heatmap_stages.pdf'));
fprintf('Fig S7 saved.\n');

% =========================================================================
%% PRINT STAGE STATISTICS — ONE-WAY REPEATED-MEASURES ANOVA
% =========================================================================
fprintf('\n=== ONE-WAY RM-ANOVA: STAGE EFFECT (all blocks pooled) ===\n');
fprintf('Each row = subject, each column = stage.\n');
for mi = 1:numel(measures)
    m = measures{mi};
    mat = SMAT.(m).all;
    ok_subj = ~any(isnan(mat),2);
    if sum(ok_subj) < 3, continue; end
    [p_rm,tbl_rm] = anova1(mat(ok_subj,:),[],'off');
    F_stat = tbl_rm{2,5};
    df_num = tbl_rm{2,3};
    df_den = tbl_rm{3,3};
    if ~isempty(F_stat) && ~isnan(F_stat)
        fprintf('  %-14s  F(%d,%d)=%.2f  p=%.4f\n',m,df_num,df_den,F_stat,p_rm);
    end
end
fprintf('\nAll stage figures saved to:\n  %s\n',outpath);


% =========================================================================
%% STAGE & BLOCK-TYPE STATISTICS: PAIRED T-TESTS WITH NORMALITY CHECKS
%
% Normality of pairwise differences is tested BEFORE every t-test using:
%   1. Shapiro-Wilk (swtest) if available on the MATLAB path
%   2. Lilliefors test (lillietest) from the Statistics Toolbox
%   3. Standardised one-sample KS test (kstest) as a last fallback
%
% Decision rule:
%   If normality NOT rejected (p >= .05) --> paired t-test
%   If normality IS rejected  (p <  .05) --> Wilcoxon signed-rank
%
% All test names, normality p-values, test statistics, df, and final
% p-values are printed so the choice of test is fully transparent.
% =========================================================================
fprintf('\n');
fprintf('=========================================================\n');
fprintf(' PAIRWISE STAGE COMPARISONS — t-test / Wilcoxon\n');
fprintf(' (normality tested on pairwise differences)\n');
fprintf('=========================================================\n');

key_pairs = {[1 3],[2 4],[1 2],[3 4],[2 3]};
key_names = {'LN vs RN','LE vs RE','LN vs LE','RN vs RE','LE vs RN (reversal cost)'};

stage_stat_measures = {'correct','confidence','RT','isHit','isFA','isCR','isMiss'};

for mi = 1:numel(stage_stat_measures)
    m = stage_stat_measures{mi};
    fprintf('\n--- %s ---\n',m);
    for pi = 1:numel(key_pairs)
        a = key_pairs{pi}(1);  b = key_pairs{pi}(2);
        va = SMAT.(m).all(:,a);
        vb = SMAT.(m).all(:,b);
        print_paired_test(key_names{pi}, va, vb, '  ');
    end
end

fprintf('\n');
fprintf('=========================================================\n');
fprintf(' D vs P BLOCK-TYPE COMPARISONS — t-test / Wilcoxon\n');
fprintf(' (within each stage; paired by subject)\n');
fprintf('=========================================================\n');

for mi = 1:numel(stage_stat_measures)
    m = stage_stat_measures{mi};
    fprintf('\n--- %s ---\n',m);
    for stgi = 1:n_stg
        va = SMAT.(m).D(:,stgi);
        vb = SMAT.(m).P(:,stgi);
        lbl = sprintf('D vs P at stage %s',STAGES{stgi});
        print_paired_test(lbl, va, vb, '  ');
    end
end

fprintf('\n');
fprintf('=========================================================\n');
fprintf(' WIN-STAY vs LOSE-STAY — t-test / Wilcoxon\n');
fprintf('=========================================================\n');
if has_stay
    for stgi = 1:n_stg
        va = SMAT_winstay(:,stgi);
        vb = SMAT_losestay(:,stgi);
        lbl = sprintf('Win-stay vs Lose-stay at stage %s',STAGES{stgi});
        print_paired_test(lbl, va, vb, '  ');
    end
else
    fprintf('  stay_choice not available — skipped.\n');
end

fprintf('\n');
fprintf('=========================================================\n');
fprintf(' CORRECT vs INCORRECT RT — t-test / Wilcoxon\n');
fprintf('=========================================================\n');
for stgi = 1:n_stg
    va = SMAT_RT_corr(:,stgi);
    vb = SMAT_RT_incorr(:,stgi);
    lbl = sprintf('Correct vs Incorrect RT at stage %s',STAGES{stgi});
    print_paired_test(lbl, va, vb, '  ');
end


% =========================================================================
%% NEW: CUMULATIVE PRIOR UNCERTAINTY (n_prev_P) + TRANSITION STATISTICS
% =========================================================================
GT = group_T;
GT.subj_id_s        = string(GT.subjID);
GT.block            = double(GT.block);
GT.correct          = double(GT.correct);
GT.block_type_clean = string(GT.block_type);
GT.block_type_clean(GT.block_type_clean=="V") = "P";

GT.prev_block_type = strings(height(GT),1);
GT.transition2     = strings(height(GT),1);
GT.n_prev_P        = nan(height(GT),1);
subs = unique(GT.subj_id_s);
for si = 1:numel(subs)
    sm = GT.subj_id_s == subs(si);
    blks = sort(unique(GT.block(sm)))';
    for bi = 1:numel(blks)
        b = blks(bi); bm = sm & GT.block==b;
        curr = GT.block_type_clean(find(bm,1));
        if bi==1
            prev='first'; trans='first'; nP=0;
        else
            prevb = blks(bi-1);
            prev  = GT.block_type_clean(find(sm & GT.block==prevb,1));
            trans = prev + "->" + curr;
            priorP = 0;
            for pb = blks(1:bi-1)
                if GT.block_type_clean(find(sm & GT.block==pb,1))=="P", priorP=priorP+1; end
            end
            nP = priorP;
        end
        GT.prev_block_type(bm) = prev;
        GT.transition2(bm)     = trans;
        GT.n_prev_P(bm)        = nP;
    end
end

GT.stage = group_T.stage;

lme_tbl = GT(ismember(GT.stage,{'LN','LE','RN','RE'}) & GT.transition2~="first",:);
lme_tbl.subj_id     = categorical(lme_tbl.subj_id_s);
lme_tbl.stage       = categorical(string(lme_tbl.stage),{'LN','LE','RN','RE'});
lme_tbl.transition2 = categorical(lme_tbl.transition2,{'D->D','D->P','P->D','P->P'});
try
    mdl_trans = fitglme(lme_tbl,'correct ~ stage * transition2 + (1|subj_id)', ...
        'Distribution','Binomial','Link','logit');
    fprintf('\n=== LME: accuracy ~ stage * transition ===\n'); disp(mdl_trans.Coefficients);
catch ME
    warning('Transition LME failed: %s',ME.message);
end

lme_tbl.n_prev_P_z = (lme_tbl.n_prev_P - mean(lme_tbl.n_prev_P,'omitnan')) ./ std(lme_tbl.n_prev_P,'omitnan');
try
    mdl_nprevP = fitglme(lme_tbl,'correct ~ stage * n_prev_P_z + (1|subj_id)', ...
        'Distribution','Binomial','Link','logit');
    fprintf('\n=== LME: accuracy ~ stage * n_prev_P ===\n'); disp(mdl_nprevP.Coefficients);
catch ME
    warning('n_prev_P LME failed: %s',ME.message);
end

% ── Fig S8: reversal dynamics by experienced uncertainty ──────────────────
fig_np = figure('Position',[60 60 1100 440]);
sgtitle('Reversal dynamics by experienced uncertainty');

tt = {'D->D','D->P','P->D','P->P'};
trans_rn_vals = cell(4,1);  % collect per-subject RN values per transition
for k = 1:4
    v = [];
    for si = 1:numel(subs)
        m_mask = lme_tbl.subj_id==categorical(subs(si)) & ...
                 lme_tbl.transition2==tt{k} & lme_tbl.stage=='RN';
        if any(m_mask), v(end+1) = mean(lme_tbl.correct(m_mask),'omitnan'); end %#ok<AGROW>
    end
    trans_rn_vals{k} = v;
end

% Panel A: RN accuracy by transition
axA = subplot(1,2,1); hold(axA,'on');
title(axA,'RN accuracy by transition');
for k = 1:4
    v = trans_rn_vals{k};
    if isempty(v), continue; end
    bar(axA,k,mean(v,'omitnan'),0.6,'FaceColor',[0.3 0.5 0.75],'EdgeColor','none','FaceAlpha',0.75);
    errorbar(axA,k,mean(v,'omitnan'),std(v,'omitnan')/sqrt(numel(v)),'k.','LineWidth',1.2,'CapSize',5);
    % individual subject dots
    scatter(axA, k*ones(size(v)), v, 20, [0.3 0.3 0.3], 'filled', ...
        'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
end
set(axA,'XTick',1:4,'XTickLabel',tt);
ylabel(axA,'P(correct) at RN'); ylim(axA,[0 1]);
axis(axA,'square');

% Panel B: RN accuracy by n_prev_P
axB = subplot(1,2,2); hold(axB,'on');
title(axB,'RN accuracy by prior P exposure');
np_levels = unique(lme_tbl.n_prev_P(~isnan(lme_tbl.n_prev_P)))';
for k = 1:numel(np_levels)
    v = [];
    for si = 1:numel(subs)
        m_mask = lme_tbl.subj_id==categorical(subs(si)) & ...
                 lme_tbl.n_prev_P==np_levels(k) & lme_tbl.stage=='RN';
        if any(m_mask), v(end+1) = mean(lme_tbl.correct(m_mask),'omitnan'); end %#ok<AGROW>
    end
    if isempty(v), continue; end
    bar(axB,k,mean(v,'omitnan'),0.6,'FaceColor',[0.80 0.40 0.20],'EdgeColor','none','FaceAlpha',0.75);
    errorbar(axB,k,mean(v,'omitnan'),std(v,'omitnan')/sqrt(numel(v)),'k.','LineWidth',1.2,'CapSize',5);
    scatter(axB, k*ones(size(v)), v, 20, [0.3 0.3 0.3], 'filled', ...
        'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
end
set(axB,'XTick',1:numel(np_levels),'XTickLabel',string(np_levels));
xlabel(axB,'# previous P blocks'); ylabel(axB,'P(correct) at RN'); ylim(axB,[0 1]);
axis(axB,'square');

if exist('save_fig','file')
    save_fig(fig_np,fullfile(outpath,'figS8_reversal_by_uncertainty'));
else
    exportgraphics(fig_np,fullfile(outpath,'figS8_reversal_by_uncertainty.pdf'),'ContentType','vector');
end
fprintf('Fig S8 saved.\n');

% ── RN accuracy: pairwise t-tests across transition types ─────────────────
fprintf('\n');
fprintf('=========================================================\n');
fprintf(' RN ACCURACY BY TRANSITION TYPE — pairwise t-tests\n');
fprintf('=========================================================\n');
trans_pairs = nchoosek(1:4,2);
for pi = 1:size(trans_pairs,1)
    ka = trans_pairs(pi,1); kb = trans_pairs(pi,2);
    va = trans_rn_vals{ka}';
    vb = trans_rn_vals{kb}';
    lbl = sprintf('%s vs %s',tt{ka},tt{kb});
    % independent samples (subjects can appear in only one transition per block)
    if numel(va)<3 || numel(vb)<3
        fprintf('  %s: too few observations (n=%d, n=%d) — skipped\n',lbl,numel(va),numel(vb));
        continue
    end
    [use_t, p_norm_a, test_nm_a] = check_normality_local(va);
    [use_t2, p_norm_b, ~]        = check_normality_local(vb);
    use_t = use_t && use_t2;
    if use_t
        [~,p_t,~,st] = ttest2(va,vb);
        fprintf('  %s: Normality OK (%s, p=%.3f & %.3f) --> independent t-test: t(%d)=%.3f, p=%.4f  (M=%.3f vs %.3f)\n', ...
            lbl, test_nm_a, p_norm_a, p_norm_b, st.df, st.tstat, p_t, mean(va,'omitnan'), mean(vb,'omitnan'));
    else
        [p_w,~,st_w] = ranksum(va,vb,'method','approximate');
        fprintf('  %s: Normality REJECTED (%s, p=%.3f/%.3f) --> Wilcoxon rank-sum: z=%.3f, p=%.4f  (M=%.3f vs %.3f)\n', ...
            lbl, test_nm_a, p_norm_a, p_norm_b, st_w.zval, p_w, mean(va,'omitnan'), mean(vb,'omitnan'));
    end
end


% =========================================================================
%% LOCAL HELPER FUNCTIONS
% =========================================================================

function gs = stage_sem(mat)
%STAGE_SEM  Per-column SEM across subjects (rows), ignoring NaN.
n  = sum(~isnan(mat),1);
gs = std(mat,0,1,'omitnan') ./ sqrt(max(n,1));
end

% ─────────────────────────────────────────────────────────────────────────
function plot_spaghetti(ax, mat, x_pos)
%PLOT_SPAGHETTI  Draw one thin translucent line per subject (row) across
%  the stage positions in x_pos.  Points at NaN stages are skipped.
%  Transparency requires R2019b+; on older MATLAB the lines are just gray.
for si = 1:size(mat,1)
    row   = mat(si,:);
    valid = ~isnan(row);
    if sum(valid) < 2, continue; end
    hl = plot(ax, x_pos(valid), row(valid), 'o-', ...
        'Color',         [0.45 0.45 0.45], ...
        'LineWidth',     0.7, ...
        'MarkerSize',    3, ...
        'MarkerFaceColor',[0.45 0.45 0.45], ...
        'MarkerEdgeColor','none', ...
        'HandleVisibility','off');
    try, hl.Color(4) = 0.25; catch, end   % transparency (R2019b+)
end
end

% ─────────────────────────────────────────────────────────────────────────
function plot_paired_dots(ax, mat_a, mat_b, x_pos, xoff_a, xoff_b, ...
                          clr_a, clr_b, scale)
%PLOT_PAIRED_DOTS  For every subject and every stage, plot a dot for
%  condition A and a dot for condition B, connected by a thin gray line.
%  scale: multiply values before plotting (use 1000 to convert s→ms).
if nargin < 9, scale = 1; end
for si = 1:size(mat_a,1)
    for xi = 1:numel(x_pos)
        xa = x_pos(xi) + xoff_a;
        xb = x_pos(xi) + xoff_b;
        va = mat_a(si,xi) * scale;
        vb = mat_b(si,xi) * scale;

        if ~isnan(va)
            plot(ax,xa,va,'o','Color',clr_a,'MarkerSize',3, ...
                'MarkerFaceColor',clr_a,'MarkerEdgeColor','none', ...
                'HandleVisibility','off');
        end
        if ~isnan(vb)
            plot(ax,xb,vb,'o','Color',clr_b,'MarkerSize',3, ...
                'MarkerFaceColor',clr_b,'MarkerEdgeColor','none', ...
                'HandleVisibility','off');
        end
        if ~isnan(va) && ~isnan(vb)
            hl = plot(ax,[xa xb],[va vb],'-', ...
                'Color',[0.5 0.5 0.5],'LineWidth',0.6, ...
                'HandleVisibility','off');
            try, hl.Color(4) = 0.20; catch, end
        end
    end
end
end

% ─────────────────────────────────────────────────────────────────────────
function [use_ttest, p_norm, test_name] = check_normality_local(d)
%CHECK_NORMALITY_LOCAL  Test whether vector d is plausibly normal.
%  Priority:  Shapiro-Wilk (swtest)  >  Lilliefors (lillietest)  >  KS
%
%  Returns:
%    use_ttest  — true  if normality is NOT rejected  (use t-test)
%                 false if normality IS  rejected      (use Wilcoxon)
%    p_norm     — p-value of the normality test
%    test_name  — string naming the test used
d = d(~isnan(d));
if numel(d) < 4
    use_ttest = true; p_norm = NaN;
    test_name = sprintf('n=%d<4,assumed normal',numel(d));
    return
end
if exist('swtest','file') == 2
    [h_n, p_norm] = swtest(d, 0.05);
    test_name = 'Shapiro-Wilk';
elseif license('test','statistics_toolbox')
    try
        [h_n, p_norm] = lillietest(d);
        test_name = 'Lilliefors';
    catch
        dz = (d - mean(d)) / std(d);
        [h_n, p_norm] = kstest(dz);
        test_name = 'KS(standardised)';
    end
else
    dz = (d - mean(d)) / std(d);
    [h_n, p_norm] = kstest(dz);
    test_name = 'KS(standardised)';
end
use_ttest = (h_n == 0);   % h=0 → fail to reject normality → use t-test
end

% ─────────────────────────────────────────────────────────────────────────
function print_paired_test(label, va, vb, indent)
%PRINT_PAIRED_TEST  Paired t-test or Wilcoxon with normality justification.
%  Prints one line per comparison to the command window.
if nargin < 4, indent = ''; end
ok = ~isnan(va) & ~isnan(vb);
n  = sum(ok);
if n < 3
    fprintf('%s%s: n=%d paired obs — skipped\n', indent, label, n);
    return
end
d_ok  = va(ok) - vb(ok);
va_ok = va(ok);
vb_ok = vb(ok);

[use_t, p_norm, test_name] = check_normality_local(d_ok);

if use_t
    [~, p_t, ~, st] = ttest(va_ok, vb_ok);
    sig = sig_stars(p_t);
    fprintf('%s%s  [n=%d | normality: %s p=%.3f --> paired t-test]\n', ...
        indent, label, n, test_name, p_norm);
    fprintf('%s  t(%d)=%.3f, p=%.4f %s  (M_a=%.4f, M_b=%.4f, mean diff=%.4f)\n', ...
        indent, st.df, st.tstat, p_t, sig, mean(va_ok), mean(vb_ok), mean(d_ok));
else
    [p_w, ~, st_w] = signrank(va_ok, vb_ok, 'method', 'approximate');
    sig = sig_stars(p_w);
    fprintf('%s%s  [n=%d | normality REJECTED: %s p=%.3f --> Wilcoxon signed-rank]\n', ...
        indent, label, n, test_name, p_norm);
    fprintf('%s  z=%.3f, p=%.4f %s  (M_a=%.4f, M_b=%.4f, mean diff=%.4f)\n', ...
        indent, st_w.zval, p_w, sig, mean(va_ok), mean(vb_ok), mean(d_ok));
end
end

% ─────────────────────────────────────────────────────────────────────────
function s = sig_stars(p)
%SIG_STARS  Return *** / ** / * / ns string for a p-value.
if     p < 0.001, s = '***';
elseif p < 0.01,  s = '**';
elseif p < 0.05,  s = '*';
else,             s = 'ns';
end
end

% ─────────────────────────────────────────────────────────────────────────
function add_pairwise_brackets(ax, mat, x_pos, y_top)
%ADD_PAIRWISE_BRACKETS  Draw significance brackets for LN vs RN and LE vs RE.
%  Also prints full t-test / Wilcoxon stats (with normality check) to console.
pairs  = {[1 3],[2 4]};
labels = {'LN vs RN','LE vs RE'};
y_br   = y_top * [0.90 0.96];

for pi = 1:2
    a = pairs{pi}(1); b = pairs{pi}(2);
    va = mat(:,a); vb = mat(:,b);
    ok = ~isnan(va) & ~isnan(vb);
    if sum(ok) < 3, continue; end

    d_ok  = va(ok) - vb(ok);
    va_ok = va(ok); vb_ok = vb(ok);

    [use_t, p_norm, test_name] = check_normality_local(d_ok);

    if use_t
        [~, p, ~, stats] = ttest(va_ok, vb_ok);
        stat_str = sprintf('t(%d)=%.2f',stats.df,stats.tstat);
    else
        [p, ~, stats] = signrank(va_ok, vb_ok, 'method', 'approximate');
        stat_str = sprintf('z=%.2f',stats.zval);
    end

    % Print to console
    test_label = ternary_local(use_t,'t-test','Wilcoxon');
    fprintf('  Bracket [%s]: norm=%s(p=%.3f) --> %s: %s, p=%.4f %s\n', ...
        labels{pi}, test_name, p_norm, test_label, stat_str, p, sig_stars(p));

    % Draw bracket on axis
    if     p < 0.001, star = '***';
    elseif p < 0.01,  star = '**';
    elseif p < 0.05,  star = '*';
    else,             star = 'ns';
    end

    xa = x_pos(a); xb = x_pos(b);
    yb = y_br(pi);
    line(ax,[xa xa xb xb],[yb*0.97 yb yb yb*0.97], ...
        'Color','k','LineWidth',0.8,'HandleVisibility','off');
    text(ax,mean([xa xb]),yb*1.005, ...
        sprintf('%s %s',labels{pi},star), ...
        'HorizontalAlignment','center','FontSize',7,'HandleVisibility','off');
end
end

% ─────────────────────────────────────────────────────────────────────────
function out = ternary_local(cond, a, b)
if cond, out = a; else, out = b; end
end
