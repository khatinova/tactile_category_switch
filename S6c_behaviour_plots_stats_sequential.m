% ==========================================================================
% S6c_behaviour_plots_stats_sequential.m
%
% PIPELINE STEP 6c — behavioural stage + sequential-block figures/statistics.
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
%   Flexible stage windows controlled by STAGE_WIN.
%   LN — Learning Naive    : first STAGE_WIN trials of block
%   LE — Learning Expert   : STAGE_WIN trials before reversal, truncated if needed
%   RN — Reversal Naive    : STAGE_WIN trials after reversal, truncated if needed
%   RE — Reversal Expert   : last STAGE_WIN trials of block
%
% STAGE ASSIGNMENT POLICY
%   This version protects the edge stages LN and RE. LN and RE get their full
%   STAGE_WIN trials whenever the block is long enough. If candidate windows
%   overlap, LE/RN are truncated first rather than stealing trials from LN/RE.
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

% Set this to 20, 15, 10, etc.
% For your requested sensitivity checks, run once with STAGE_WIN = 15 and
% once with STAGE_WIN = 10.
STAGE_WIN  = 15;
STAGE_ASSIGNMENT_POLICY = 'protect_edges';  % protects LN and RE; truncates LE/RN on overlap
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
outpath = fullfile(base_path,'Results','Behav results', ...
    sprintf('S6c_Stage_Sequential_Figures_%02dtrial_%s',STAGE_WIN,STAGE_ASSIGNMENT_POLICY));
if ~exist(outpath,'dir'), mkdir(outpath); end
cd(base_path)
fprintf('S6 stage window: %d trials | assignment policy: %s\n',STAGE_WIN,STAGE_ASSIGNMENT_POLICY);
fprintf('Outputs will be saved to: %s\n',outpath);

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
% Candidate windows:
%   LN = first STAGE_WIN trials in the block
%   LE = STAGE_WIN trials immediately before reversal
%   RN = STAGE_WIN trials starting at reversal
%   RE = last STAGE_WIN trials in the block
%
% Important change from the older S6 scripts:
%   The older scripts assigned overlapping windows with priority
%       RE > RN > LE > LN
%   so LN could be shortened by LE and RN could be shortened by RE.
%
% This version protects LN and RE. If windows overlap, LE and/or RN are
% shortened first. That gives stable edge baselines and edge endpoints.
nRows = height(group_T);
group_T.stage = repmat("",nRows,1);

subj_block_pairs = unique(group_T(:,{'subjID','block'}),'rows');
stage_count_log = table();

for pb = 1:height(subj_block_pairs)
    sn  = subj_block_pairs.subjID(pb);
    blk = subj_block_pairs.block(pb);
    mask = group_T.subjID==sn & group_T.block==blk;
    rows = find(mask);
    if isempty(rows), continue; end

    rev_t  = group_T.revTrial(rows(1));
    trials = group_T.trial(rows);
    max_t  = max(trials);

    % Protected edge windows.
    LN_start = 1;
    LN_end   = min(STAGE_WIN,max_t);
    RE_start = max(1,max_t-STAGE_WIN+1);
    RE_end   = max_t;

    % If a block is shorter than 2*STAGE_WIN, LN and RE cannot both be full.
    % In that rare case, keep both labels non-overlapping by splitting the
    % available block; this avoids double-labelling the same trial.
    if RE_start <= LN_end
        midpoint = floor(max_t/2);
        LN_end   = midpoint;
        RE_start = midpoint + 1;
        warning('%s block %g is too short for full LN and RE windows; split edges as LN=1:%d, RE=%d:%d.', ...
            sn,blk,LN_end,RE_start,RE_end);
    end

    LE_start = NaN; LE_end = NaN; RN_start = NaN; RN_end = NaN;
    if ~isnan(rev_t) && isfinite(rev_t)
        rev_t = round(rev_t);

        % Middle windows are truncated to avoid stealing trials from LN/RE.
        LE_start = max(LN_end+1, rev_t-STAGE_WIN);
        LE_end   = min(rev_t-1, RE_start-1);

        RN_start = max(rev_t, LN_end+1);
        RN_end   = min(rev_t+STAGE_WIN-1, RE_start-1);
    end

    % Assign edges first; they are protected and non-overlapping.
    group_T.stage(rows(trials>=LN_start & trials<=LN_end)) = "LN";
    group_T.stage(rows(trials>=RE_start & trials<=RE_end)) = "RE";

    % Assign middle phases only where valid after truncation.
    if ~isnan(LE_start) && LE_start <= LE_end
        group_T.stage(rows(trials>=LE_start & trials<=LE_end)) = "LE";
    end
    if ~isnan(RN_start) && RN_start <= RN_end
        group_T.stage(rows(trials>=RN_start & trials<=RN_end)) = "RN";
    end

    % Log per-block stage counts so you can verify what was truncated.
    row_log = table(sn,blk,STAGE_WIN, ...
        sum(group_T.stage(rows)=="LN"), sum(group_T.stage(rows)=="LE"), ...
        sum(group_T.stage(rows)=="RN"), sum(group_T.stage(rows)=="RE"), ...
        'VariableNames',{'subjID','block','stage_win','n_LN','n_LE','n_RN','n_RE'});
    stage_count_log = [stage_count_log; row_log]; %#ok<AGROW>
end

T = group_T(group_T.stage ~= "",:);
fprintf('Valid stage labels: %d / %d rows\n',height(T),nRows);
fprintf('Stage count summary across subject-blocks:\n');
disp(varfun(@mean,stage_count_log,'InputVariables',{'n_LN','n_LE','n_RN','n_RE'}));
writetable(stage_count_log, fullfile(outpath,sprintf('stage_count_log_%02dtrial.csv',STAGE_WIN)));
fprintf('Stage-count log saved: %s\n',fullfile(outpath,sprintf('stage_count_log_%02dtrial.csv',STAGE_WIN)));

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
%% S6c ADD-ON — SEQUENTIAL BLOCK / n_prev_P FIGURES FROM E SCRIPT
% =========================================================================
% This section folds the standalone sequential-block behavioural figures into
% the S6 pipeline. It reuses the S6-preprocessed group_T and saves outputs in
% the same S6c folder. Stage assignment above remains flexible via STAGE_WIN
% and protects LN/RE before truncating LE/RN.

fprintf('\n=========================================================\n');
fprintf(' S6c ADD-ON: sequential-block behavioural figures\n');
fprintf('=========================================================\n');

% Keep the S6 style consistent.
set(groot,'defaultAxesTickDir','out');
set(groot,'defaultAxesBox','off');
set(groot,'defaultAxesLineWidth',0.5);
set(groot,'defaultAxesFontSize',10);

% Load all_trial_data when available. These figures need aligned_correct and
% stimulus-level trial_data; if the file is absent, the S6 core figures still run.
all_trial_data = [];
all_trial_candidates = { ...
    fullfile(base_path,'Data','all_trial_data_June2026.mat'), ...
    fullfile(base_path,'Data','all_trial_data.mat')};
loaded_all_trials = false;
for ci = 1:numel(all_trial_candidates)
    if exist(all_trial_candidates{ci},'file')
        S_all = load(all_trial_candidates{ci});
        if isfield(S_all,'all_trial_data')
            all_trial_data = S_all.all_trial_data;
            fprintf('Loaded all_trial_data: %s\n', all_trial_candidates{ci});
            loaded_all_trials = true;
            break
        end
    end
end
if ~loaded_all_trials
    warning('S6c sequential figures skipped because all_trial_data was not found.');
else
    % Optional Nassar/model results. Only used if switch_stims_by_block exists.
    results = [];
    sim_candidates = { ...
        fullfile(base_path,'Results','Simulation results','Figures','nassar_results.mat'), ...
        fullfile(base_path,'Results','Simulation results','Figures','sim_data.mat'), ...
        fullfile(base_path,'Results','EEG analysis','Epoched_data','group_feature_table_combined_v9c_RL.mat')};
    for ci = 1:numel(sim_candidates)
        if exist(sim_candidates{ci},'file')
            try
                tmp = load(sim_candidates{ci});
                if isfield(tmp,'results')
                    results = tmp.results;
                    fprintf('Loaded optional results struct: %s\n', sim_candidates{ci});
                    break
                end
            catch ME_load_results
                warning('Could not load optional results from %s: %s', sim_candidates{ci}, ME_load_results.message);
            end
        end
    end
    if isempty(results)
        fprintf('Optional results struct not found; switched stimuli default to [2 3] unless detectable.\n');
    end

    % Colour scheme: preserve S6 D/P and stage colours, add sequential colours.
    CLR_KH   = CLR_D;
    CLR_RR   = CLR_P;
    CLR_STGS_MAT = vertcat(CLR_STGS{:});
    CLR_SWITCH = [0.84 0.15 0.16];
    CLR_MAINT  = [0.12 0.47 0.71];
    CLR_BOTH   = [0.50 0.50 0.50];
    CLR_NPP = [0.26 0.58 0.78; 0.97 0.58 0.02; 0.70 0.17 0.12];
    TRANS_TYPES_ASCII  = {'D->D','D->P','P->D','P->P'};
    NPP_BIN_LABELS = {'0 prior P','1 prior P','2+ prior P'};

    % Reversal-aligned axis. Keep the ±30-trial alignment used in the pasted E script.
    preN       = 30;
    postN      = 30;
    alignedLen = preN + postN;
    rel_ax     = -preN : (postN-1);
    SWITCH_STIMS_DEFAULT = [2 3];
    MAINT_STIMS_DEFAULT  = [1 4]; %#ok<NASGU>

    % Ensure n_prev_P exists on group_T for consistency with S6/S6c.
    group_T.n_prev_P = zeros(height(group_T), 1);
    all_subjs_s6c = unique(group_T.subjID);
    for si_np = 1:numel(all_subjs_s6c)
        sn_np = all_subjs_s6c(si_np);
        smask_np = group_T.subjID == sn_np;
        blocks_np = sort(unique(group_T.block(smask_np)))';
        block_types_s = strings(numel(blocks_np),1);
        for bi_np = 1:numel(blocks_np)
            b_rows = smask_np & group_T.block == blocks_np(bi_np);
            bt_vals = unique(string(group_T.block_type_clean(b_rows)));
            bt_vals(bt_vals=="V") = "P";
            if isempty(bt_vals), block_types_s(bi_np) = "D"; else, block_types_s(bi_np) = bt_vals(1); end
        end
        for bi_np = 1:numel(blocks_np)
            b_rows = smask_np & group_T.block == blocks_np(bi_np);
            group_T.n_prev_P(b_rows) = sum(block_types_s(1:max(bi_np-1,0)) == "P");
        end
    end
    group_T.n_prev_P_bin = min(group_T.n_prev_P, 2);
    fprintf('S6c n_prev_P range: 0-%d across all participants\n', max(group_T.n_prev_P));

    % Build reversal-aligned matrices and metadata for E1-E11.
    acc_rows     = NaN(0, alignedLen);
    conf_rows    = NaN(0, alignedLen);
    acc_sw_rows  = NaN(0, alignedLen);
    acc_mn_rows  = NaN(0, alignedLen);
    conf_sw_rows = NaN(0, alignedLen);
    conf_mn_rows = NaN(0, alignedLen);
    rt_rows      = NaN(0, alignedLen);

    meta_subj       = {};
    meta_block      = [];
    meta_block_type = {};
    meta_n_prev_P   = [];
    meta_npp_bin    = [];
    meta_transition = {};
    meta_cohort     = {};
    meta_stim_config = {};
    meta_rev_trial  = [];

    subj_ids_td = fieldnames(all_trial_data);
    for si_td = 1:numel(subj_ids_td)
        sn_char = subj_ids_td{si_td};
        if ~isfield(all_trial_data.(sn_char),'trial_data'), continue; end
        td = all_trial_data.(sn_char).trial_data;
        if ~isfield(td,'revTrial') || ~isfield(td,'aligned_correct') || ~isfield(td,'correct'), continue; end
        [nB,~] = size(td.correct);
        is_kh  = startsWith(sn_char, 'Ox');

        bs = '';
        if isfield(td,'block_structure') && ~isempty(td.block_structure)
            bs = upper(char(td.block_structure));
        elseif isfield(td,'trueFB')
            bs_arr = repmat('D',1,nB);
            for b=1:nB
                pfb = td.trueFB(b, ~isnan(td.trueFB(b,:)));
                if ~isempty(pfb) && mean(pfb) < 0.99
                    bs_arr(b) = 'P';
                end
            end
            bs = bs_arr;
        end
        if isempty(bs), bs = repmat('D',1,nB); end
        bs(bs=='V') = 'P';

        n_prev_P_running = 0;
        if ~isempty(results) && isfield(results, sn_char) && isfield(results.(sn_char),'switch_stims_by_block')
            sw_by_block_s = results.(sn_char).switch_stims_by_block;
        else
            sw_by_block_s = {};
        end

        for b = 1:nB
            if b > size(td.aligned_correct,1), continue; end
            rev = td.revTrial(b);
            curr_type = char(bs(min(b,numel(bs))));
            curr_stim_config = get_block_stim_config_lc(td, group_T, sn_char, b);

            if b <= numel(sw_by_block_s) && ~isempty(sw_by_block_s{b})
                sw_b = sw_by_block_s{b};
            else
                sw_b = infer_switch_stims_lc(td,b,SWITCH_STIMS_DEFAULT);
            end
            mn_b = setdiff(1:4, sw_b);

            if b == 1
                trans_str = 'first';
            else
                prev_type = char(bs(min(b-1,numel(bs))));
                trans_str = [prev_type '->' curr_type];
            end

            acc_row = td.aligned_correct(b,:);
            if numel(acc_row) ~= alignedLen
                acc_row = resize_row_lc(acc_row, alignedLen);
            end
            acc_rows(end+1,:) = acc_row; %#ok<AGROW>

            conf_row = NaN(1,alignedLen);
            if isfield(td,'aligned_confidence') && b <= size(td.aligned_confidence,1)
                conf_row = resize_row_lc(td.aligned_confidence(b,:), alignedLen);
            end
            conf_rows(end+1,:) = conf_row; %#ok<AGROW>

            rt_row = NaN(1,alignedLen);
            if isfield(td,'aligned_rt') && b <= size(td.aligned_rt,1)
                rt_row = resize_row_lc(td.aligned_rt(b,:), alignedLen);
                rt_row(rt_row > 1) = NaN;
            elseif isfield(td,'RT') && isfinite(rev)
                rt_row = aligned_from_raw_lc(td.RT(b,:), round(rev), preN, postN);
                rt_row(rt_row > 1) = NaN;
            end
            rt_rows(end+1,:) = rt_row; %#ok<AGROW>

            rel_idx = -preN : (postN-1);
            sw_row = NaN(1,alignedLen); mn_row = NaN(1,alignedLen);
            sw_conf_row = NaN(1,alignedLen); mn_conf_row = NaN(1,alignedLen);
            sf2 = '';
            if isfield(td,'stimType'), sf2='stimType'; elseif isfield(td,'stimID'), sf2='stimID'; end
            if ~isempty(sf2) && isfinite(rev) && rev > 0
                nT_td = size(td.correct, 2);
                for w = 1:alignedLen
                    t_abs = round(rev) + rel_idx(w);
                    if t_abs < 1 || t_abs > nT_td, continue; end
                    s_id = td.(sf2)(b, t_abs);
                    if isnan(s_id), continue; end
                    if ismember(s_id, sw_b)
                        sw_row(w) = td.correct(b, t_abs);
                        if isfield(td,'confidence'), sw_conf_row(w) = td.confidence(b,t_abs); end
                    elseif ismember(s_id, mn_b)
                        mn_row(w) = td.correct(b, t_abs);
                        if isfield(td,'confidence'), mn_conf_row(w) = td.confidence(b,t_abs); end
                    end
                end
            end
            acc_sw_rows(end+1,:) = sw_row; %#ok<AGROW>
            acc_mn_rows(end+1,:) = mn_row; %#ok<AGROW>
            conf_sw_rows(end+1,:) = sw_conf_row; %#ok<AGROW>
            conf_mn_rows(end+1,:) = mn_conf_row; %#ok<AGROW>

            meta_subj{end+1}       = sn_char; %#ok<AGROW>
            meta_block(end+1)      = b; %#ok<AGROW>
            meta_block_type{end+1} = curr_type; %#ok<AGROW>
            meta_n_prev_P(end+1)   = n_prev_P_running; %#ok<AGROW>
            meta_npp_bin(end+1)    = min(n_prev_P_running, 2); %#ok<AGROW>
            meta_transition{end+1} = trans_str; %#ok<AGROW>
            meta_cohort{end+1}     = ternary_local(is_kh, 'KH', 'RR'); %#ok<AGROW>
            meta_stim_config{end+1}= curr_stim_config; %#ok<AGROW>
            meta_rev_trial(end+1)  = rev; %#ok<AGROW>

            if curr_type == 'P', n_prev_P_running = n_prev_P_running + 1; end
        end
    end

    meta_block       = meta_block(:);
    meta_n_prev_P    = meta_n_prev_P(:);
    meta_npp_bin     = meta_npp_bin(:);
    meta_rev_trial   = meta_rev_trial(:);
    meta_subj        = string(meta_subj(:));
    meta_block_type  = string(meta_block_type(:));
    meta_transition  = string(meta_transition(:));
    meta_cohort      = string(meta_cohort(:));
    meta_stim_config = string(meta_stim_config(:));
    meta_stim_config(meta_stim_config == "") = "unknown";

    is_D_row  = meta_block_type == "D";
    is_P_row  = meta_block_type == "P";
    is_KH_row = meta_cohort == "KH";

    assert(size(acc_rows,1) == numel(meta_block), 'acc_rows and metadata row counts differ');
    fprintf('\nS6c sequential rows: %d subject x block rows (D=%d, P=%d)\n', numel(meta_block), sum(is_D_row), sum(is_P_row));
    fprintf('Stim-config metadata: known=%d unknown=%d\n', sum(meta_stim_config ~= "unknown"), sum(meta_stim_config == "unknown"));
%% FIG E1 — SEQUENTIAL BLOCK LEARNING CURVES
%
% RATIONALE: If participants learn "how to learn" across blocks (meta-learning
% or Bayesian prior calibration), we expect reversal cost to shrink and
% recovery rate to grow as a function of block number. This is the
% within-subject learning curve for volatility adaptation. Split by D/P
% type because the Nassar model predicts type-specific H calibration
% (Behrens et al. 2007: learning rate for learning rate updates).
%
% Block sequential index (1–5) preserves the temporal order of experience,
% which n_prev_P alone does not encode (it ignores D-blocks).
% =========================================================================
fprintf('\n--- Building Fig E1: Sequential block learning curves ---\n');

MAX_BLOCK = 5;   % typical maximum number of real task blocks

fig_E1 = figure('Position',[50 50 1400 580]);
sgtitle({'Sequential block learning: reversal-aligned accuracy across blocks', ...
    '(Behrens et al. 2007: learning rate updates as block sequence progresses)'}, ...
    'FontSize',12);

% Pre-allocate subject × block matrices
N_sub = numel(subj_ids_td);
acc_pre_D  = NaN(N_sub, MAX_BLOCK);   % pre-reversal accuracy, D blocks
acc_post_D = NaN(N_sub, MAX_BLOCK);
acc_pre_P  = NaN(N_sub, MAX_BLOCK);
acc_post_P = NaN(N_sub, MAX_BLOCK);

% Also: all reversal-aligned traces split by block number (1..MAX_BLOCK)
acc_by_block = cell(1, MAX_BLOCK);   % each cell: rows = (subj×block)
type_by_block = cell(1, MAX_BLOCK);
for bk = 1:MAX_BLOCK, acc_by_block{bk} = NaN(0,alignedLen); type_by_block{bk} = {}; end

for si = 1:N_sub
    sn_char = subj_ids_td{si};
    td = all_trial_data.(sn_char).trial_data;
    if ~isfield(td,'aligned_correct') || ~isfield(td,'revTrial'), continue; end

    [nB,~] = size(td.correct);
    bs = '';
    if isfield(td,'block_structure') && ~isempty(td.block_structure)
        bs = upper(char(td.block_structure));
    else
        bs = repmat('D',1,nB);
    end
    bs(bs=='V') = 'P';

    for b = 1:min(nB, MAX_BLOCK)
        acc_row = td.aligned_correct(b,:);
        pre_acc  = mean(acc_row(1:preN), 'omitnan');
        post_acc = mean(acc_row(preN+1:end), 'omitnan');
        curr_type = char(bs(min(b,numel(bs))));

        si_idx = find(strcmp(subj_ids_td, sn_char));
        if curr_type == 'D'
            acc_pre_D(si_idx,b)  = pre_acc;
            acc_post_D(si_idx,b) = post_acc;
        else
            acc_pre_P(si_idx,b)  = pre_acc;
            acc_post_P(si_idx,b) = post_acc;
        end

        acc_by_block{b}(end+1,:) = acc_row;
        type_by_block{b}{end+1}  = curr_type;
    end
end

% ── E1a: Accuracy ribbon per block number, all types pooled ───────────────
ax_E1a = subplot(2, 3, [1 2]); hold(ax_E1a,'on');
title(ax_E1a,'Reversal-aligned accuracy per sequential block (all types)', 'FontSize',10);

block_clr = lines(MAX_BLOCK);
for bk = 1:MAX_BLOCK
    mat = acc_by_block{bk};
    if isempty(mat) || all(isnan(mat(:))), continue; end
    plot_ribbon_lc(ax_E1a, rel_ax, mat, block_clr(bk,:), '-', ...
        sprintf('Block %d (n=%d)', bk, size(mat,1)));
end
xline(ax_E1a, 0, 'k--','LineWidth',1.5,'HandleVisibility','off');
yline(ax_E1a, 0.5,'k:','HandleVisibility','off');
patch(ax_E1a,[-preN 0 0 -preN],[0 0 1 1],[0.9 0.9 0.9],'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
xlabel(ax_E1a,'Trial relative to reversal');
ylabel(ax_E1a,'P(correct)');
xlim(ax_E1a,[-preN postN-1]); ylim(ax_E1a,[0.3 1]);
legend(ax_E1a,'Box','off','FontSize',8,'Location','southeast','NumColumns',2);

% ── E1b: Pre vs post accuracy — block number trajectory ───────────────────
ax_E1b = subplot(2,3,3); hold(ax_E1b,'on');
title(ax_E1b,'Reversal cost across block sequence','FontSize',10);

bk_ax = 1:MAX_BLOCK;
for bk = 1:MAX_BLOCK
    mat = acc_by_block{bk};
    if isempty(mat), continue; end
    pre_m  = mean(mat(:,1:preN),2,'omitnan');
    post_m = mean(mat(:,preN+1:end),2,'omitnan');
    cost   = pre_m - post_m;   % positive = performance drop at reversal
    cost   = cost(~isnan(cost));
    if isempty(cost), continue; end
    errorbar(ax_E1b, bk, mean(cost,'omitnan'), sem_lc(cost), ...
        'o-','Color',block_clr(bk,:),'MarkerFaceColor',block_clr(bk,:), ...
        'MarkerSize',8,'LineWidth',1.5,'HandleVisibility','off');
end
yline(ax_E1b,0,'k:','HandleVisibility','off');
xlabel(ax_E1b,'Block number in session');
ylabel(ax_E1b,'Reversal cost: P(correct)_{pre} − P(correct)_{post}');
set(ax_E1b,'XTick',1:MAX_BLOCK,'TickDir','out'); xlim(ax_E1b,[0.5 MAX_BLOCK+0.5]);
subtitle(ax_E1b,'↑ = larger accuracy drop at reversal','FontSize',8,'Color',[0.5 0.5 0.5]);

% ── E1c/d: Pre-reversal accuracy trajectory per block type ─────────────────
for bt_idx = 1:2
    bt_tag = ternary_lc(bt_idx==1, 'D', 'P');
    pre_mat  = ternary_lc(bt_idx==1, acc_pre_D,  acc_pre_P);
    post_mat = ternary_lc(bt_idx==1, acc_post_D, acc_post_P);
    clr_bt   = ternary_lc(bt_idx==1, CLR_D, CLR_P);

    ax_bt = subplot(2,3,3+bt_idx); hold(ax_bt,'on');
    title(ax_bt,sprintf('Pre vs post accuracy — %s blocks',bt_tag),'FontSize',10);

    for bk = 1:MAX_BLOCK
        pre_col  = pre_mat(:,bk);  pre_col  = pre_col(~isnan(pre_col));
        post_col = post_mat(:,bk); post_col = post_col(~isnan(post_col));
        if isempty(pre_col), continue; end

        bk_jitter = (bk-0.12)*[1 1];
        errorbar(ax_bt,bk-0.15,mean(pre_col,'omitnan'),sem_lc(pre_col), ...
            'o','Color',CLR_STGS_MAT(2,:),'MarkerFaceColor',CLR_STGS_MAT(2,:), ...
            'MarkerSize',7,'LineWidth',1.5,'HandleVisibility','off');
        if ~isempty(post_col)
            errorbar(ax_bt,bk+0.15,mean(post_col,'omitnan'),sem_lc(post_col), ...
                's','Color',CLR_STGS_MAT(3,:),'MarkerFaceColor',CLR_STGS_MAT(3,:), ...
                'MarkerSize',7,'LineWidth',1.5,'HandleVisibility','off');
        end
    end
    % Legend entries
    plot(ax_bt,NaN,NaN,'o-','Color',CLR_STGS_MAT(2,:),'MarkerFaceColor',CLR_STGS_MAT(2,:),'DisplayName','Pre-rev (LE)');
    plot(ax_bt,NaN,NaN,'s-','Color',CLR_STGS_MAT(3,:),'MarkerFaceColor',CLR_STGS_MAT(3,:),'DisplayName','Post-rev (RN)');
    xlabel(ax_bt,'Block number in session');
    ylabel(ax_bt,'P(correct)');
    set(ax_bt,'XTick',1:MAX_BLOCK, 'TickDir', 'out'); xlim(ax_bt,[0.5 MAX_BLOCK+0.5]); ylim(ax_bt,[0.4 1]);
    legend(ax_bt,'Box','off','FontSize',8,'Location','best');
end

% ── E1e: Recovery slope across blocks ─────────────────────────────────────
ax_E1e = subplot(2,3,6); hold(ax_E1e,'on');
title(ax_E1e,'Post-reversal recovery slope per block','FontSize',10);
set(ax_E1e,'TickDir','out');
for bk = 1:MAX_BLOCK
    mat = acc_by_block{bk};
    if isempty(mat), continue; end

    slopes = NaN(size(mat,1),1);
    for ri = 1:size(mat,1)
        post_seg = mat(ri, preN+5:end);   % exclude first ~5 post-rev trials (initial confusion)
        post_t   = (1:numel(post_seg));
        ok = ~isnan(post_seg);
        if sum(ok) > 3
            p = polyfit(post_t(ok), post_seg(ok), 1);
            slopes(ri) = p(1);
        end
    end
    slopes = slopes(~isnan(slopes));
    if isempty(slopes), continue; end
    errorbar(ax_E1e, bk, mean(slopes,'omitnan'), sem_lc(slopes), ...
        'o-','Color',block_clr(bk,:),'MarkerFaceColor',block_clr(bk,:), ...
        'MarkerSize',8,'LineWidth',1.5,'HandleVisibility','off');
end
yline(ax_E1e,0,'k:','HandleVisibility','off');
xlabel(ax_E1e,'Block number in session');
ylabel(ax_E1e,'Recovery slope (Δ accuracy / trial)');
set(ax_E1e,'XTick',1:MAX_BLOCK, 'TickDir','out'); xlim(ax_E1e,[0.5 MAX_BLOCK+0.5]);
subtitle(ax_E1e,'↑ = faster post-rev recovery','FontSize',8,'Color',[0.5 0.5 0.5]);

annotation('textbox',[0.01 0.01 0.98 0.04],'String', ...
    ['Fig E1: Sequential learning curve. Nassar (2010) model predicts H calibration improves across blocks; '...
     'Behrens et al. (2007) show learning rate for learning rate updates consolidates with experience. '...
     'Reversal cost = LE accuracy − RN accuracy (pre-reversal expert minus post-reversal naive). '...
     'Recovery slope fitted on trials +5 to +30 post-reversal (avoids initial confusion period).'],...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig_E1, fullfile(outpath,'figE1_sequential_block_learning.pdf'));
saveas(fig_E1, fullfile(outpath,'figE1_sequential_block_learning.png'));
fprintf('Fig E1 saved.\n');

% =========================================================================
%% FIG E2 — REVERSAL COST TRAJECTORY (subject-level)
%
% RATIONALE: This is the within-subject analogue of Behrens et al. (2007)
% Fig 3: the learner's response to volatility should sharpen as they
% accumulate experience. We expect a significant block_number × block_type
% interaction: the D→P transition should show heightened reversal cost
% (more noise confusable with a change-point) while accumulated P experience
% should reduce spurious resets, lowering cost in later blocks.
% =========================================================================
fprintf('--- Building Fig E2: Reversal cost trajectory ---\n');

fig_E2 = figure('Position',[50 50 1300 500]);
sgtitle({'Reversal cost trajectory across sequential blocks', ...
    '(Yu & Dayan 2005: unexpected uncertainty disrupts reliable change-point detection)'}, ...
    'FontSize',12);

% Build subject × block cost matrix per type
cost_D = NaN(N_sub, MAX_BLOCK);
cost_P = NaN(N_sub, MAX_BLOCK);

for si_idx = 1:numel(subj_ids_td)
    sn_char = subj_ids_td{si_idx};
    td = all_trial_data.(sn_char).trial_data;
    if ~isfield(td,'aligned_correct') || ~isfield(td,'revTrial'), continue; end
    [nB,~] = size(td.correct);

    bs = '';
    if isfield(td,'block_structure') && ~isempty(td.block_structure)
        bs = upper(char(td.block_structure)); bs(bs=='V')='P';
    else
        bs = repmat('D',1,nB);
    end

    for b = 1:min(nB,MAX_BLOCK)
        acc_row  = td.aligned_correct(b,:);
        pre_acc  = mean(acc_row(1:preN), 'omitnan');
        post_acc = mean(acc_row(preN+1:min(preN+10,alignedLen)), 'omitnan');
        cost_val = pre_acc - post_acc;
        curr_type = char(bs(min(b,numel(bs))));
        if curr_type == 'D'
            cost_D(si_idx,b) = cost_val;
        else
            cost_P(si_idx,b) = cost_val;
        end
    end
end

% ── E2a: Subject-level trajectories — D blocks ───────────────────────────
ax_E2a = subplot(1,3,1); hold(ax_E2a,'on');
title(ax_E2a,'Reversal cost: Deterministic blocks','FontSize',10);
for si_idx = 1:N_sub
    row = cost_D(si_idx,:);
    ok = ~isnan(row);
    if sum(ok) > 1
        plot(ax_E2a, find(ok), row(ok), 'o-','Color',fade_lc(CLR_D,0.75),'MarkerSize',4,'HandleVisibility','off');
    end
end
% Group mean ± SEM
for bk = 1:MAX_BLOCK
    col = cost_D(:,bk); col=col(~isnan(col));
    if ~isempty(col)
        errorbar(ax_E2a,bk,mean(col),sem_lc(col),'ko','MarkerFaceColor',CLR_D,...
            'MarkerSize',9,'LineWidth',2,'HandleVisibility','off');
    end
end
yline(ax_E2a,0,'k:','HandleVisibility','off');
xlabel(ax_E2a,'Block number in session'); ylabel(ax_E2a,'Cost: P(correct)_{pre} − P(correct)_{post}');
set(ax_E2a,'XTick',1:MAX_BLOCK, 'TickDir','out'); xlim(ax_E2a,[0.5 MAX_BLOCK+0.5]);
subtitle(ax_E2a,'Light = individual; Dark = group mean ± SEM','FontSize',8,'Color',[0.5 0.5 0.5]);

% ── E2b: D blocks ─────────────────────────────────────────────────────────
ax_E2b = subplot(1,3,2); hold(ax_E2b,'on');
title(ax_E2b,'Reversal cost: Probabilistic blocks','FontSize',10);
for si_idx = 1:N_sub
    row = cost_P(si_idx,:);
    ok = ~isnan(row);
    if sum(ok) > 1
        plot(ax_E2b, find(ok), row(ok),'o-','Color',fade_lc(CLR_P,0.75),'MarkerSize',4,'HandleVisibility','off');
    end
end
for bk = 1:MAX_BLOCK
    col = cost_P(:,bk); col=col(~isnan(col));
    if ~isempty(col)
        errorbar(ax_E2b,bk,mean(col),sem_lc(col),'ko','MarkerFaceColor',CLR_P,...
            'MarkerSize',9,'LineWidth',2,'HandleVisibility','off');
    end
end
yline(ax_E2b,0,'k:','HandleVisibility','off');
xlabel(ax_E2b,'Block number in session'); ylabel(ax_E2b,'Reversal cost');
set(ax_E2b,'XTick',1:MAX_BLOCK, 'TickDir','out'); xlim(ax_E2b,[0.5 MAX_BLOCK+0.5]);

% ── E2c: D vs P overlaid ─────────────────────────────────────────────────
ax_E2c = subplot(1,3,3); hold(ax_E2c,'on');
title(ax_E2c,'D vs P reversal cost (group mean ± SEM)','FontSize',10);
for bk = 1:MAX_BLOCK
    col_d = cost_D(:,bk); col_d=col_d(~isnan(col_d));
    col_p = cost_P(:,bk); col_p=col_p(~isnan(col_p));
    if ~isempty(col_d)
        errorbar(ax_E2c,bk-0.15,mean(col_d),sem_lc(col_d),'o-','Color',CLR_D,...
            'MarkerFaceColor',CLR_D,'MarkerSize',8,'LineWidth',1.8,'HandleVisibility','off');
    end
    if ~isempty(col_p)
        errorbar(ax_E2c,bk+0.15,mean(col_p),sem_lc(col_p),'s--','Color',CLR_P,...
            'MarkerFaceColor',CLR_P,'MarkerSize',8,'LineWidth',1.8,'HandleVisibility','off');
    end
end
plot(ax_E2c,NaN,NaN,'o-','Color',CLR_D,'MarkerFaceColor',CLR_D,'DisplayName','Deterministic');
plot(ax_E2c,NaN,NaN,'s--','Color',CLR_P,'MarkerFaceColor',CLR_P,'DisplayName','Probabilistic');
yline(ax_E2c,0,'k:','HandleVisibility','off');
xlabel(ax_E2c,'Block number'); ylabel(ax_E2c,'Reversal cost');
set(ax_E2c,'XTick',1:MAX_BLOCK, 'TickDir','out'); xlim(ax_E2c,[0.5 MAX_BLOCK+0.5]);
legend(ax_E2c,'Box','off','FontSize',9,'Location','best');
subtitle(ax_E2c,'Error = ±1 SEM across subjects','FontSize',8,'Color',[0.5 0.5 0.5]);

saveas(fig_E2, fullfile(outpath,'figE2_reversal_cost_trajectory.pdf'));
saveas(fig_E2, fullfile(outpath,'figE2_reversal_cost_trajectory.png'));
fprintf('Fig E2 saved.\n');

% =========================================================================
%% FIG E3 — SWITCHED vs MAINTAINED: REVERSAL-ALIGNED ACCURACY
%
% RATIONALE: The dimensional shift design means only switched stimuli
% carry a genuine change-point signal at reversal. Maintained stimuli
% have the same correct response rule before and after reversal, so any
% accuracy dip on maintained stimuli reflects response competition /
% attentional capture from the switched-stimulus conflict, not true
% uncertainty about the rule. Dissociating these two contributions is
% critical for interpreting FRN amplitude:
%   FRN on switched trials ↔ genuine reward prediction error (Holroyd & Coles 2002)
%   FRN on maintained trials ↔ conflict monitoring signal (Cavanagh & Frank 2014)
% =========================================================================
fprintf('--- Building Fig E3: Switch vs maintained accuracy ---\n');

fig_E3 = figure('Position',[50 50 1400 580]);
sgtitle({'Switched vs maintained stimulus accuracy: reversal-aligned profiles', ...
    '(Collins & Frank 2013; Cavanagh & Frank 2014)'}, 'FontSize',12);

for bt_idx = 1:3   % 1=D, 2=P, 3=All
    if bt_idx == 1,      bt_str='D'; bt_mask=is_D_row; clr_bt=CLR_D; lbl_bt='Deterministic';
    elseif bt_idx == 2,  bt_str='P'; bt_mask=is_P_row; clr_bt=CLR_P; lbl_bt='Probabilistic';
    else,                bt_str='all'; bt_mask=true(numel(meta_block),1); clr_bt=CLR_BOTH; lbl_bt='All blocks';
    end

    sw_mat = acc_sw_rows(bt_mask,:);
    mn_mat = acc_mn_rows(bt_mask,:);
    all_mat = acc_rows(bt_mask,:);

    ax_top = subplot(2,3,bt_idx); hold(ax_top,'on');
    title(ax_top,lbl_bt,'FontSize',10);
    plot_ribbon_lc(ax_top, rel_ax, sw_mat, CLR_SWITCH, '--', ...
        sprintf('Switched (n=%d)', size(sw_mat,1)));
    plot_ribbon_lc(ax_top, rel_ax, mn_mat, CLR_MAINT, '-', ...
        sprintf('Maintained (n=%d)', size(mn_mat,1)));
    plot_ribbon_lc(ax_top, rel_ax, all_mat, [0.6 0.6 0.6], ':', ...
        'All stimuli');
    xline(ax_top,0,'k--','LineWidth',1.5,'HandleVisibility','off');
    yline(ax_top,0.5,'k:','HandleVisibility','off');
    patch(ax_top,[-preN 0 0 -preN],[0 0 1 1],[0.9 0.9 0.9], ...
        'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
    if bt_idx==1, ylabel(ax_top,'P(correct)'); end
    xlabel(ax_top,'Trial relative to reversal');
    xlim(ax_top,[-preN postN-1]); ylim(ax_top,[0.2 1]);
    legend(ax_top,'Box','off','FontSize',8,'Location','southeast');

    % ── Bottom panel: pre/post comparison per stim type ────────────────
    ax_bot = subplot(2,3,3+bt_idx); hold(ax_bot,'on');
    title(ax_bot,sprintf('%s: pre vs post comparison',lbl_bt),'FontSize',10);

    pre_idx  = 1:preN;
    post_idx = preN+1:preN+15;   % first 15 post-rev trials (RN window)

    sw_pre  = mean(sw_mat(:,pre_idx),  2,'omitnan');
    sw_post = mean(sw_mat(:,post_idx), 2,'omitnan');
    mn_pre  = mean(mn_mat(:,pre_idx),  2,'omitnan');
    mn_post = mean(mn_mat(:,post_idx), 2,'omitnan');

    pairs = {sw_pre, sw_post, CLR_SWITCH, 'Switched';
             mn_pre, mn_post, CLR_MAINT,  'Maintained'};

    x_positions = [1 2;  4 5];
    for pi = 1:2
        pre_v  = pairs{pi,1}; post_v = pairs{pi,2};
        clr_pi = pairs{pi,3}; lbl_pi = pairs{pi,4};
        ok_pi  = ~isnan(pre_v) & ~isnan(post_v);
        xp     = x_positions(pi,:);

        if sum(ok_pi) < 2, continue; end

        % Subject lines
        for sj = 1:numel(pre_v)
            if ok_pi(sj)
                plot(ax_bot, xp, [pre_v(sj) post_v(sj)], '-', ...
                    'Color',fade_lc(clr_pi,0.82),'HandleVisibility','off');
            end
        end
        % Bar + scatter
        bar(ax_bot, xp(1), mean(pre_v(ok_pi)), 0.4,'FaceColor',clr_pi,'FaceAlpha',0.55,'EdgeColor','none','HandleVisibility','off');
        bar(ax_bot, xp(2), mean(post_v(ok_pi)),0.4,'FaceColor',clr_pi,'FaceAlpha',0.85,'EdgeColor','none','HandleVisibility','off');
        errorbar(ax_bot, xp, [mean(pre_v(ok_pi)) mean(post_v(ok_pi))], ...
            [sem_lc(pre_v(ok_pi)) sem_lc(post_v(ok_pi))], ...
            'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');

        [~,p_pi] = ttest(pre_v(ok_pi), post_v(ok_pi));
        add_sig_bracket(ax_bot, xp(1), xp(2), ...
            max([pre_v(ok_pi);post_v(ok_pi)],[],'omitnan')*1.05, p_pi, lbl_pi);
    end

    set(ax_bot,'XTick',[1.5 4.5],'XTickLabel',{'Switched','Maintained'},'FontSize',9, 'TickDir','out');
    ylabel(ax_bot,'P(correct)'); ylim(ax_bot,[0.2 1.1]);
    xline(ax_bot,3,'k:','HandleVisibility','off');
    text(1,0.23,'Pre-rev','HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
    text(2,0.23,'Post-rev','HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
    text(4,0.23,'Pre-rev','HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
    text(5,0.23,'Post-rev','HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
end

annotation('textbox',[0.01 0.01 0.98 0.04],'String', ...
    ['Fig E3: Switched stimuli (Go↔NoGo reversed) vs maintained stimuli (unchanged assignment). '...
     'Post-reversal accuracy dip on switched stimuli = true change-point cost. '...
     'Any dip on maintained = response competition / conflict (Cavanagh & Frank 2014 TICS). '...
     'Post-rev window = first 15 trials post-reversal (Reversal Naive stage).'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig_E3, fullfile(outpath,'figE3_switch_vs_maintained_accuracy.pdf'));
saveas(fig_E3, fullfile(outpath,'figE3_switch_vs_maintained_accuracy.png'));
fprintf('Fig E3 saved.\n');

% =========================================================================
%% FIG E4 — SWITCHED vs MAINTAINED: CONFIDENCE AND RT
%
% RATIONALE: The task order is stimulus → response → CONFIDENCE → feedback.
% Confidence is therefore a pre-outcome belief measure (prospective certainty
% |θ−0.5|) not a post-feedback summary. Collins & Frank (2013) show that
% stimulus-specific value tracking drives confidence on stimulus-type
% transitions. We expect:
%  - Confidence drops sharply on switched stimuli after reversal
%    (the agent's θ for that stimulus has been pushed away from 0.5 by
%     the prior rule, and now faces PE; Boldt & Yeung 2015).
%  - Maintained stimuli: confidence is stable or slightly dips due
%    to generalised uncertainty / increased arousal.
%  - RT on switched stimuli increases post-reversal (response conflict);
%    Cavanagh & Frank 2014 predict frontal theta power tracks this conflict.
% =========================================================================
fprintf('--- Building Fig E4: Switch vs maintained confidence + RT ---\n');

fig_E4 = figure('Position',[50 50 1300 560]);
sgtitle({'Switched vs maintained: confidence (pre-outcome) and RT profiles', ...
    '(Boldt & Yeung 2015; Cavanagh & Frank 2014 — confidence as precision-weighted belief)'}, ...
    'FontSize',12);

for row_i = 1:2
    if row_i==1, data_sw=conf_sw_rows; data_mn=conf_mn_rows; ylbl='Confidence (1–10)'; ylims=[1 10]; lbl_='Confidence';
    else,        data_sw=NaN(size(acc_sw_rows)); data_mn=NaN(size(acc_mn_rows));
        % RT data: reconstruct similarly
        % (rt_rows not split by stim; approximation using overall)
        ylbl='RT (ms)'; ylims=[200 700]; lbl_='RT';
    end

    for bt_idx = 1:3
        if bt_idx==1,     bt_mask=is_D_row; lbl_bt='Det'; clr_bt=CLR_D;
        elseif bt_idx==2, bt_mask=is_P_row; lbl_bt='Prob'; clr_bt=CLR_P;
        else,             bt_mask=true(numel(meta_block),1); lbl_bt='All'; clr_bt=CLR_BOTH;
        end

        ax = subplot(2,3,(row_i-1)*3+bt_idx); hold(ax,'on');
        title(ax,sprintf('%s: %s blocks',lbl_,lbl_bt),'FontSize',10);

        if row_i==1
            sw_dat = conf_sw_rows(bt_mask,:);
            mn_dat = conf_mn_rows(bt_mask,:);
        else
            % For RT we use the all-stimulus RT split post-hoc by block type
            sw_dat = rt_rows(bt_mask,:)*1000;   % convert to ms
            mn_dat = rt_rows(bt_mask,:)*1000;   % same (stim-level RT not pre-split)
        end

        plot_ribbon_lc(ax, rel_ax, sw_dat, CLR_SWITCH, '--', 'Switched');
        if row_i==1
            plot_ribbon_lc(ax, rel_ax, mn_dat, CLR_MAINT,  '-',  'Maintained');
        end
        xline(ax,0,'k--','LineWidth',1.5,'HandleVisibility','off');
        patch(ax,[-preN 0 0 -preN],[ylims(1) ylims(1) ylims(2) ylims(2)], ...
            [0.9 0.9 0.9],'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
        xlabel(ax,'Trial relative to reversal');
        if bt_idx==1, ylabel(ax,ylbl); end
        xlim(ax,[-preN postN-1]); ylim(ax,ylims);
        legend(ax,'Box','off','FontSize',8,'Location','best');
    end
end

annotation('textbox',[0.01 0.01 0.98 0.04],'String', ...
    ['Fig E4: Confidence (top) and RT (bottom) for switched vs maintained stimuli. '...
     'Confidence is rated after the decision but BEFORE feedback — it is a prospective precision '...
     'weight, not a post-hoc evaluation (Boldt & Yeung 2015). RT conflict on switched stimuli '...
     'predicts frontal theta power (Cavanagh & Frank 2014 TICS).'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);
set(findall(gcf,'Type','axes'),'TickDir','out','Box','off');
saveas(fig_E4, fullfile(outpath,'figE4_switch_vs_maintained_conf_RT.pdf'));
saveas(fig_E4, fullfile(outpath,'figE4_switch_vs_maintained_conf_RT.png'));
fprintf('Fig E4 saved.\n');

% =========================================================================
%% FIG E5 — n_prev_P: REVERSAL-ALIGNED ACCURACY AND CONFIDENCE
%
% RATIONALE: n_prev_P is the within-subject index of cumulative noise
% exposure at the time of each reversal. It operationalises the
% "uncertainty history" manipulation embedded in the PPDPD block sequence.
% The Bayesian prediction (Yu & Dayan 2005; Behrens et al. 2007):
%   n_prev_P = 0: participant has no noise history → interprets every
%     prediction error as potentially signal → high reversal cost expected.
%   n_prev_P ≥ 2: well-calibrated noise model → smaller H inflation on P
%     blocks → faster, more selective reversal detection → lower cost.
% Note: this is a within-subject, cross-block comparison, so individual
% differences in overall learning speed are controlled.
% =========================================================================
fprintf('--- Building Fig E5: n_prev_P effect on reversal-aligned performance ---\n');

fig_E5 = figure('Position',[50 50 1400 560]);
sgtitle({'Cumulative probabilistic block exposure (n_{prev P}) and reversal adaptation', ...
    '(Yu & Dayan 2005; Behrens et al. 2007: calibrating uncertainty over history)'}, ...
    'FontSize',12);

npp_bins = 0:2;
npp_labels = NPP_BIN_LABELS;

% ── E5a–c: Accuracy by n_prev_P bin ──────────────────────────────────────
ax_E5a = subplot(2,3,1); hold(ax_E5a,'on');
title(ax_E5a,'All blocks: accuracy by n_{prev P}','FontSize',10);

for ni = 1:3
    bin_mask = meta_npp_bin == npp_bins(ni);
    mat = acc_rows(bin_mask,:);
    if isempty(mat), continue; end
    plot_ribbon_lc(ax_E5a, rel_ax, mat, CLR_NPP(ni,:), '-', ...
        sprintf('%s (n=%d)',npp_labels{ni},size(mat,1)));
end
xline(ax_E5a,0,'k--','LineWidth',1.5,'HandleVisibility','off');
yline(ax_E5a,0.5,'k:','HandleVisibility','off');
patch(ax_E5a,[-preN 0 0 -preN],[0 0 1 1],[0.9 0.9 0.9],'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
xlabel(ax_E5a,'Trial relative to reversal'); ylabel(ax_E5a,'P(correct)');
xlim(ax_E5a,[-preN postN-1]); ylim(ax_E5a,[0.3 1]);
legend(ax_E5a,'Box','off','FontSize',8,'Location','southeast');

% ── E5b: D blocks only ────────────────────────────────────────────────────
ax_E5b = subplot(2,3,2); hold(ax_E5b,'on');
title(ax_E5b,'D blocks: accuracy by n_{prev P}','FontSize',10);
for ni = 1:3
    bin_mask = meta_npp_bin == npp_bins(ni) & is_D_row;
    mat = acc_rows(bin_mask,:);
    if isempty(mat), continue; end
    plot_ribbon_lc(ax_E5b, rel_ax, mat, CLR_NPP(ni,:), '-', ...
        sprintf('%s (n=%d)',npp_labels{ni},size(mat,1)));
end
xline(ax_E5b,0,'k--','LineWidth',1.5,'HandleVisibility','off');
yline(ax_E5b,0.5,'k:','HandleVisibility','off');
patch(ax_E5b,[-preN 0 0 -preN],[0 0 1 1],[0.9 0.9 0.9],'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
xlabel(ax_E5b,'Trial relative to reversal'); ylabel(ax_E5b,'P(correct)');
xlim(ax_E5b,[-preN postN-1]); ylim(ax_E5b,[0.3 1]);
legend(ax_E5b,'Box','off','FontSize',8,'Location','southeast');
subtitle(ax_E5b,sprintf('n_{prev P} = prior P blocks; current block = D',1),'FontSize',8,'Color',[0.5 0.5 0.5]);

% ── E5c: P blocks only ────────────────────────────────────────────────────
ax_E5c = subplot(2,3,3); hold(ax_E5c,'on');
title(ax_E5c,'P blocks: accuracy by n_{prev P}','FontSize',10);
for ni = 1:3
    bin_mask = meta_npp_bin == npp_bins(ni) & is_P_row;
    mat = acc_rows(bin_mask,:);
    if isempty(mat), continue; end
    plot_ribbon_lc(ax_E5c, rel_ax, mat, CLR_NPP(ni,:), '--', ...
        sprintf('%s (n=%d)',npp_labels{ni},size(mat,1)));
end
xline(ax_E5c,0,'k--','LineWidth',1.5,'HandleVisibility','off');
yline(ax_E5c,0.5,'k:','HandleVisibility','off');
patch(ax_E5c,[-preN 0 0 -preN],[0 0 1 1],[0.9 0.9 0.9],'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
xlabel(ax_E5c,'Trial relative to reversal'); ylabel(ax_E5c,'P(correct)');
xlim(ax_E5c,[-preN postN-1]); ylim(ax_E5c,[0.3 1]);
legend(ax_E5c,'Box','off','FontSize',8,'Location','southeast');

% ── E5d–f: Confidence by n_prev_P (same layout) ──────────────────────────
for ni_plot = 1:3
    sp_idx = 3 + ni_plot;
    if ni_plot==1, ax_conf=subplot(2,3,sp_idx); hold(ax_conf,'on'); title(ax_conf,'Confidence: all blocks','FontSize',10); conf_row_use=conf_rows; end
    if ni_plot==2, ax_conf=subplot(2,3,sp_idx); hold(ax_conf,'on'); title(ax_conf,'Confidence: D blocks','FontSize',10); conf_row_use=conf_rows; end
    if ni_plot==3, ax_conf=subplot(2,3,sp_idx); hold(ax_conf,'on'); title(ax_conf,'Confidence: P blocks','FontSize',10); conf_row_use=conf_rows; end
    for ni = 1:3
        if ni_plot==1, bin_mask = meta_npp_bin == npp_bins(ni);
        elseif ni_plot==2, bin_mask = meta_npp_bin == npp_bins(ni) & is_D_row;
        else, bin_mask = meta_npp_bin == npp_bins(ni) & is_P_row; end
        mat = conf_row_use(bin_mask,:);
        if isempty(mat), continue; end
        ls = ternary_lc(ni_plot==3,'--','-');
        plot_ribbon_lc(ax_conf, rel_ax, mat, CLR_NPP(ni,:), ls, ...
            sprintf('%s (n=%d)',npp_labels{ni},size(mat,1)));
    end
    xline(ax_conf,0,'k--','LineWidth',1.5,'HandleVisibility','off');
    patch(ax_conf,[-preN 0 0 -preN],[1 1 10 10],[0.9 0.9 0.9],'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
    xlabel(ax_conf,'Trial relative to reversal');
    if ni_plot==1, ylabel(ax_conf,'Confidence (1–10)'); end
    xlim(ax_conf,[-preN postN-1]); ylim(ax_conf,[1 10]);
    legend(ax_conf,'Box','off','FontSize',7,'Location','best');
end

annotation('textbox',[0.01 0.01 0.98 0.04],'String', ...
    ['Fig E5: n_{prev P} = number of probabilistic blocks experienced before the current block. '...
     'Colours: blue=0 prior P blocks, amber=1, red=2+. '...
     'Prediction (Behrens et al. 2007): more noise history → better calibrated H → lower reversal cost. '...
     'Confidence is pre-outcome (rated after response, before feedback; Boldt & Yeung 2015).'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);
set(findall(gcf,'Type','axes'),'TickDir','out','Box','off');
saveas(fig_E5, fullfile(outpath,'figE5_nPrevP_revaligned.pdf'));
saveas(fig_E5, fullfile(outpath,'figE5_nPrevP_revaligned.png'));
fprintf('Fig E5 saved.\n');

% =========================================================================
%% FIG E6 — n_prev_P: REVERSAL COST AND RECOVERY SCATTER
%
% RATIONALE: Quantifies the n_prev_P effect at the level of summary
% statistics used in the planned LME models. Shows both the group-level
% trend and within-subject structure (each point is one block, paired
% within subjects). A fitted regression line tests the linear hypothesis
% that each additional P-block reduces reversal cost by a fixed amount.
% =========================================================================
fprintf('--- Building Fig E6: n_prev_P scatter + regression ---\n');

fig_E6 = figure('Position',[50 50 1200 500]);
sgtitle({'n_{prev P} effect on reversal cost and recovery rate', ...
    '(Within-subject moderator: cumulative noise exposure)'}, 'FontSize',12);

% Build (block, subject) level data frame
cost_vec  = NaN(numel(meta_block),1);
recov_vec = NaN(numel(meta_block),1);
for ri = 1:numel(meta_block)
    row = acc_rows(ri,:);
    pre_  = mean(row(1:preN),'omitnan');
    post_ = mean(row(preN+1:preN+10),'omitnan');
    cost_vec(ri)  = pre_ - post_;
    post_long = row(preN+5:end);
    t_post    = 1:numel(post_long);
    ok = ~isnan(post_long);
    if sum(ok) > 3
        pf = polyfit(t_post(ok), post_long(ok), 1);
        recov_vec(ri) = pf(1);
    end
end

% ── E6a: n_prev_P vs reversal cost (continuous) ──────────────────────────
ax_E6a = subplot(1,3,1); hold(ax_E6a,'on');
title(ax_E6a,'n_{prev P} vs reversal cost (all blocks)','FontSize',10);
n_pp_cont = meta_n_prev_P(:);
ok_E6 = ~isnan(cost_vec) & ~isnan(n_pp_cont);

% Scatter per cohort
kh_rows_E6 = ok_E6 & is_KH_row;
rr_rows_E6 = ok_E6 & ~is_KH_row;
scatter(ax_E6a, n_pp_cont(kh_rows_E6)+0.05*(rand(sum(kh_rows_E6),1)-0.5), ...
    cost_vec(kh_rows_E6), 30, CLR_KH,'filled','MarkerFaceAlpha',0.4,'DisplayName','Ox (KH)');
scatter(ax_E6a, n_pp_cont(rr_rows_E6)+0.05*(rand(sum(rr_rows_E6),1)-0.5), ...
    cost_vec(rr_rows_E6), 30, CLR_RR,'o','MarkerEdgeAlpha',0.6,'DisplayName','Nc (RR)');

if sum(ok_E6) > 5
    [rv_cost,pv_cost] = corr(n_pp_cont(ok_E6), cost_vec(ok_E6),'Rows','complete','Type','Spearman');
    xi_fit = linspace(min(n_pp_cont(ok_E6)), max(n_pp_cont(ok_E6)),100);
    pf_cost = polyfit(n_pp_cont(ok_E6), cost_vec(ok_E6), 1);
    plot(ax_E6a, xi_fit, polyval(pf_cost,xi_fit),'k-','LineWidth',2,'HandleVisibility','off');
    text(ax_E6a,0.05,0.97,sprintf('ρ=%.2f, p=%.3f',rv_cost,pv_cost), ...
        'Units','normalized','VerticalAlignment','top','FontSize',9,'BackgroundColor','w');
end
yline(ax_E6a,0,'k:','HandleVisibility','off');
xlabel(ax_E6a,'n_{prev P} (cumulative prior P blocks)');
ylabel(ax_E6a,'Reversal cost: acc_{pre} − acc_{post}');
legend(ax_E6a,'Box','off','FontSize',8,'Location','best');
subtitle(ax_E6a,'Spearman ρ (block as obs, nested in subject)','FontSize',8,'Color',[0.5 0.5 0.5]);

% ── E6b: n_prev_P vs recovery rate ───────────────────────────────────────
ax_E6b = subplot(1,3,2); hold(ax_E6b,'on');
title(ax_E6b,'n_{prev P} vs post-rev recovery rate','FontSize',10);
ok_E6b = ~isnan(recov_vec) & ~isnan(n_pp_cont);
scatter(ax_E6b, n_pp_cont(ok_E6b & is_KH_row), recov_vec(ok_E6b & is_KH_row), ...
    30, CLR_KH,'filled','MarkerFaceAlpha',0.4,'DisplayName','Ox');
scatter(ax_E6b, n_pp_cont(ok_E6b & ~is_KH_row), recov_vec(ok_E6b & ~is_KH_row), ...
    30, CLR_RR,'o','MarkerEdgeAlpha',0.6,'DisplayName','Nc');
if sum(ok_E6b) > 5
    [rv_rec,pv_rec] = corr(n_pp_cont(ok_E6b), recov_vec(ok_E6b),'Rows','complete','Type','Spearman');
    xi_fit2 = linspace(min(n_pp_cont(ok_E6b)), max(n_pp_cont(ok_E6b)),100);
    pf_rec = polyfit(n_pp_cont(ok_E6b), recov_vec(ok_E6b), 1);
    plot(ax_E6b, xi_fit2, polyval(pf_rec,xi_fit2),'k-','LineWidth',2,'HandleVisibility','off');
    text(ax_E6b,0.05,0.97,sprintf('ρ=%.2f, p=%.3f',rv_rec,pv_rec), ...
        'Units','normalized','VerticalAlignment','top','FontSize',9,'BackgroundColor','w');
end
yline(ax_E6b,0,'k:','HandleVisibility','off');
xlabel(ax_E6b,'n_{prev P}'); ylabel(ax_E6b,'Recovery slope (Δacc/trial, +5 to +30)');
legend(ax_E6b,'Box','off','FontSize',8,'Location','best');

% ── E6c: Binned means — cleaner group summary ─────────────────────────────
ax_E6c = subplot(1,3,3); hold(ax_E6c,'on');
title(ax_E6c,'Reversal cost by n_{prev P} bin','FontSize',10);

for ni = 1:3
    bin_mask_ni = meta_npp_bin == npp_bins(ni) & ok_E6;
    cost_ni = cost_vec(bin_mask_ni);
    if isempty(cost_ni), continue; end

    bar(ax_E6c, ni, mean(cost_ni,'omitnan'), 0.55, 'FaceColor',CLR_NPP(ni,:), ...
        'EdgeColor','none','FaceAlpha',0.8,'HandleVisibility','off');
    errorbar(ax_E6c, ni, mean(cost_ni,'omitnan'), sem_lc(cost_ni), ...
        'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');

    jx = ni + 0.18*(rand(numel(cost_ni),1)-0.5);
    scatter(ax_E6c, jx, cost_ni, 18, [0.3 0.3 0.3],'filled','MarkerFaceAlpha',0.3,'HandleVisibility','off');
    text(ax_E6c, ni, 0.01, sprintf('n=%d',numel(cost_ni)),...
        'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
end

% Pairwise significance brackets
if sum(meta_npp_bin==0 & ok_E6) > 2 && sum(meta_npp_bin==2 & ok_E6) > 2
    c0 = cost_vec(meta_npp_bin==0 & ok_E6);
    c2 = cost_vec(meta_npp_bin==2 & ok_E6);
    [~,p_npp] = ttest2(c0, c2);
    add_sig_bracket(ax_E6c, 1, 3, max(cost_vec(ok_E6))*1.05, p_npp, '0 vs 2+');
end

yline(ax_E6c,0,'k:','HandleVisibility','off');
set(ax_E6c,'XTick',1:3,'XTickLabel',NPP_BIN_LABELS,'FontSize',9, 'TickDir','out');
ylabel(ax_E6c,'Reversal cost');

annotation('textbox',[0.01 0.01 0.98 0.04],'String', ...
    ['Fig E6: Recovery slope = linear fit to accuracy trials +5 to +30 post-reversal (positive=recovering). '...
     'Spearman ρ used (Nassar latents are non-normally distributed). Note: block is nested within subject; '...
     'formal inference requires LME with random slopes (see C_individual_differences_uncertainty_v1.m).'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig_E6, fullfile(outpath,'figE6_nPrevP_cost_recovery.pdf'));
saveas(fig_E6, fullfile(outpath,'figE6_nPrevP_cost_recovery.png'));
fprintf('Fig E6 saved.\n');

% =========================================================================
%% FIG E7 — n_prev_P × BLOCK_TYPE INTERACTION ON REVERSAL COST
%
% RATIONALE: The key interaction predicted by Yu & Dayan (2005):
%   - In D blocks (no noise), reversal is unambiguous → n_prev_P should
%     not substantially change cost (unless prior P raises H globally).
%   - In P blocks (noise present), n_prev_P calibrates noise model →
%     larger cost reduction with each additional prior P block.
% This 2×3 grid (D/P × n_prev_P 0/1/2+) directly visualises the LME
% interaction term block_type × n_prev_P planned in the EEG analysis.
% ERP prediction: FRN amplitude should show the same interaction
% (larger noise calibration benefit reflected in smaller FRN on P blocks
% with high n_prev_P), while P300 tracks the structural update
% independently of n_prev_P.
% =========================================================================
fprintf('--- Building Fig E7: n_prev_P × block_type interaction ---\n');

fig_E7 = figure('Position',[50 50 1200 560]);
sgtitle({'n_{prev P} × block type interaction on reversal performance', ...
    '(Key interaction for LME models and ERP hypotheses)'}, 'FontSize',12);

bt_tags    = {'D','P'};
bt_labels  = {'Deterministic','Probabilistic'};
bt_masks   = {is_D_row, is_P_row};
bt_clrs    = {CLR_D, CLR_P};

for bt_i = 1:2
    for ni = 1:3
        sp_pos = (bt_i-1)*3 + ni;
        ax = subplot(2,3,sp_pos); hold(ax,'on');

        bin_mask = meta_npp_bin == npp_bins(ni) & bt_masks{bt_i};
        mat = acc_rows(bin_mask,:);

        if isempty(mat) || all(isnan(mat(:)))
            text(ax,0.5,0.5,'No data','Units','normalized','HorizontalAlignment','center');
            title(ax,sprintf('%s | %s',bt_labels{bt_i},npp_labels{ni}),'FontSize',9);
            continue;
        end

        % Individual block-level traces (light)
        for ri2 = 1:size(mat,1)
            plot(ax, rel_ax, mat(ri2,:),'Color',fade_lc(bt_clrs{bt_i},0.90),'LineWidth',0.5,'HandleVisibility','off');
        end

        plot_ribbon_lc(ax, rel_ax, mat, bt_clrs{bt_i}, ...
            ternary_lc(bt_i==2,'--','-'), ...
            sprintf('n=%d blocks', size(mat,1)));

        xline(ax,0,'k--','LineWidth',1.5,'HandleVisibility','off');
        yline(ax,0.5,'k:','HandleVisibility','off');
        patch(ax,[-preN 0 0 -preN],[0 0 1 1],[0.9 0.9 0.9],'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');

        pre_ = mean(mat(:,1:preN),2,'omitnan');
        post_ = mean(mat(:,preN+1:preN+10),2,'omitnan');
        ok_t = ~isnan(pre_) & ~isnan(post_);
        if sum(ok_t) > 1
            [~,p_t] = ttest(pre_(ok_t), post_(ok_t));
            cost_mn = mean(pre_(ok_t)-post_(ok_t),'omitnan');
            sig_str = ternary_lc(p_t<0.001,'***',ternary_lc(p_t<0.01,'**',ternary_lc(p_t<0.05,'*','ns')));
            text(ax,0.02,0.97,sprintf('Cost=%.2f %s',cost_mn,sig_str), ...
                'Units','normalized','VerticalAlignment','top','FontSize',8,'BackgroundColor','w');
        end

        title(ax,sprintf('%s | %s',bt_labels{bt_i},npp_labels{ni}),'FontSize',9);
        if ni==1, ylabel(ax,'P(correct)'); end
        if bt_i==2, xlabel(ax,'Trial relative to reversal'); end
        xlim(ax,[-preN postN-1]); ylim(ax,[0.2 1]);
        legend(ax,'Box','off','FontSize',8,'Location','southeast');
    end
end

% Annotate with colour rectangle indicating block type
for bi2 = 1:2
    annotation('rectangle',[0.01 ternary_lc(bi2==1,0.55,0.08) 0.015 0.40], ...
        'FaceColor',bt_clrs{bi2},'EdgeColor','none','FaceAlpha',0.4);
end
annotation('textbox',[0.01 0.01 0.98 0.04],'String', ...
    ['Fig E7: 2×3 grid shows the block_type × n_{prev P} interaction. '...
     'The reversal cost (text in panel) should decrease with n_{prev P} more strongly in P blocks '...
     '(noise calibration effect) than in D blocks. '...
     'ERP prediction: FRN tracks this interaction; P300 tracks structural update independently of n_{prev P}.'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig_E7, fullfile(outpath,'figE7_nPrevP_by_blocktype.pdf'));
saveas(fig_E7, fullfile(outpath,'figE7_nPrevP_by_blocktype.png'));
fprintf('Fig E7 saved.\n');

% =========================================================================
%% FIG E8 — FIRST ENCOUNTER vs LATER ENCOUNTERS POST-REVERSAL
%          (Switch stimuli only; split by n_prev_P)
%
% RATIONALE: On switched stimuli, the first post-reversal encounter of
% that stimulus is the critical learning event — the agent's θ for that
% stimulus will be maximally inconsistent with the new rule. The Nassar
% model predicts ω peaks at this first encounter (see D_Bayesian_RL_models
% Fig S4), and the FRN should be largest here. Dissociating the FRN into
% first-encounter vs later allows a direct neural test of the ω spike.
% The n_prev_P split tests whether prior noise history blunts this first
% encounter signal (noise makes the first anomalous outcome look less
% diagnostic, reducing the ω-driven update).
% Reference: Nassar et al. (2010) J Neurosci §Fig 4; Cavanagh et al.
% (2012) J Neurosci for FRN indexing learning rate.
% =========================================================================
fprintf('--- Building Fig E8: First vs later post-rev encounter ---\n');

fig_E8 = figure('Position',[50 50 1200 560]);
sgtitle({'First vs later post-reversal encounters on switched stimuli', ...
    '(Tests ω-spike prediction: Nassar et al. 2010, Cavanagh et al. 2012)'}, ...
    'FontSize',12);

% Categorise each reversal-aligned position as:
%   "first_encounter" = the first time the switched stimulus is seen in the post-rev window
%   "later_encounter" = subsequent appearances
% Approximated as: trials 1–4 post-rev = likely 1st encounter of each stim
% (4 stimuli, random order ≈ geometric; E[1st encounter of stim i] ≈ 4 trials)
FIRST_ENC_WIN  = 1:4;   % trials +1 to +4 post-reversal  (position in aligned array: preN+1 to preN+4)
LATER_ENC_WIN  = 5:15;  % trials +5 to +15

first_sw_acc  = NaN(numel(meta_block),1);
later_sw_acc  = NaN(numel(meta_block),1);
first_sw_conf = NaN(numel(meta_block),1);
later_sw_conf = NaN(numel(meta_block),1);

for ri = 1:numel(meta_block)
    sw_row = acc_sw_rows(ri,:);
    sc_row = conf_sw_rows(ri,:);
    first_sw_acc(ri)  = mean(sw_row(preN + FIRST_ENC_WIN), 'omitnan');
    later_sw_acc(ri)  = mean(sw_row(preN + LATER_ENC_WIN), 'omitnan');
    first_sw_conf(ri) = mean(sc_row(preN + FIRST_ENC_WIN), 'omitnan');
    later_sw_conf(ri) = mean(sc_row(preN + LATER_ENC_WIN), 'omitnan');
end

n_pp_vec = meta_npp_bin(:);

for ni = 1:3
    bin_mask_ni = n_pp_vec == npp_bins(ni);
    fa = first_sw_acc(bin_mask_ni);   la = later_sw_acc(bin_mask_ni);
    fc = first_sw_conf(bin_mask_ni);  lc = later_sw_conf(bin_mask_ni);
    ok_a = ~isnan(fa) & ~isnan(la);
    ok_c = ~isnan(fc) & ~isnan(lc);

    % ── Accuracy panel ────────────────────────────────────────────────
    ax_a = subplot(2,3,ni); hold(ax_a,'on');
    title(ax_a,sprintf('Accuracy — %s', npp_labels{ni}),'FontSize',10);

    if sum(ok_a) > 1
        bar(ax_a,[1 2],[mean(fa(ok_a)) mean(la(ok_a))],0.45, ...
            'FaceColor',CLR_NPP(ni,:),'FaceAlpha',0.7,'EdgeColor','none');
        errorbar(ax_a,[1 2],[mean(fa(ok_a)) mean(la(ok_a))], ...
            [sem_lc(fa(ok_a)) sem_lc(la(ok_a))], ...
            'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');
        for sj_i = find(ok_a)'
            plot(ax_a,[1 2]+0.04*(rand-0.5),[fa(sj_i) la(sj_i)], ...
                '-','Color',fade_lc(CLR_NPP(ni,:),0.80),'HandleVisibility','off');
        end
        [~,p_enc] = ttest(fa(ok_a), la(ok_a));
        add_sig_bracket(ax_a, 1, 2, 1.02, p_enc, '');
    end
    set(ax_a,'XTick',[1 2],'XTickLabel',{'First','Later'}, 'TickDir','out');
    ylabel(ax_a,'P(correct on switched)'); ylim(ax_a,[0 1.15]);
    subtitle(ax_a,sprintf('n=%d blocks',sum(ok_a)),'FontSize',8,'Color',[0.5 0.5 0.5]);

    % ── Confidence panel ──────────────────────────────────────────────
    ax_c = subplot(2,3,3+ni); hold(ax_c,'on');
    title(ax_c,sprintf('Confidence — %s', npp_labels{ni}),'FontSize',10);

    if sum(ok_c) > 1
        bar(ax_c,[1 2],[mean(fc(ok_c)) mean(lc(ok_c))],0.45, ...
            'FaceColor',CLR_NPP(ni,:),'FaceAlpha',0.7,'EdgeColor','none');
        errorbar(ax_c,[1 2],[mean(fc(ok_c)) mean(lc(ok_c))], ...
            [sem_lc(fc(ok_c)) sem_lc(lc(ok_c))], ...
            'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');
        [~,p_enc_c] = ttest(fc(ok_c), lc(ok_c));
        add_sig_bracket(ax_c, 1, 2, 9.5, p_enc_c, '');
    end
    set(ax_c,'XTick',[1 2],'XTickLabel',{'First','Later'}, 'TickDir','out');
    ylabel(ax_c,'Confidence (1–10)'); ylim(ax_c,[1 10.5]);
    subtitle(ax_c,sprintf('n=%d blocks',sum(ok_c)),'FontSize',8,'Color',[0.5 0.5 0.5]);
end

annotation('textbox',[0.01 0.01 0.98 0.04],'String', ...
    ['Fig E8: First encounter ≈ trials +1 to +4 post-reversal for switched stimuli (geometric expectation). '...
     'The ω spike (Nassar et al. 2010 Fig S4) predicts: (i) lowest accuracy on first encounter, '...
     '(ii) steepest confidence drop. n_{prev P} modulates this: higher prior noise history '...
     'should blunt the ω spike (FRN prediction: smaller first-encounter FRN for n_{prev P}≥2).'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig_E8, fullfile(outpath,'figE8_first_vs_later_encounter.pdf'));
saveas(fig_E8, fullfile(outpath,'figE8_first_vs_later_encounter.png'));
fprintf('Fig E8 saved.\n');

% =========================================================================
%% FIG E9 — EEG-BEHAVIOUR BRIDGE: PREDICTED ERP PATTERNS
%  (Displays empirical data from group_table_combined if loaded;
%   otherwise shows the behavioural correlates that motivate ERP hypotheses)
%
% RATIONALE: This figure anchors the planned ERP analyses to the observed
% behavioural effects, making the logic of each neural measure explicit.
%   FRN: largest for switched, unexpected outcomes (incorrect on switched
%     post-reversal), modulated by n_prev_P — tests whether noise history
%     reduces the FRN (attenuated PE signal when outcomes are expected to
%     be noisy; Holroyd & Coles 2002 Psych Rev; Gehring & Willoughby 2002).
%   P300: largest for high-ω events (genuine change-point detection);
%     should NOT be substantially modulated by n_prev_P because it tracks
%     structural updates, not noise calibration (Polich 2007 Clin Neurophysiol).
%   Frontal theta: indexes response conflict on switched stimuli post-rev;
%     should correlate with RT increase on those trials (Cavanagh & Frank 2014).
% =========================================================================
fprintf('--- Building Fig E9: EEG-behaviour bridge ---\n');

has_eeg_table = exist('group_table_combined','var') && ...
                ismember('FRN_mean_amp', group_table_combined.Properties.VariableNames);

fig_E9 = figure('Position',[50 50 1400 560]);
sgtitle({'EEG–behaviour bridge: behavioural predictors of ERP components', ...
    '(Holroyd & Coles 2002; Polich 2007; Cavanagh & Frank 2014)'}, 'FontSize',12);

% ── Panel 1: Accuracy × reversal cost → FRN prediction ──────────────────
ax_E9a = subplot(1,4,1); hold(ax_E9a,'on');
title(ax_E9a,'FRN predictor: cost on switched stimuli','FontSize',9);
for ni = 1:3
    bin_mask_ni = meta_npp_bin == npp_bins(ni);
    cost_sw = NaN(sum(bin_mask_ni),1);
    rows_ni = find(bin_mask_ni);
    for ri_ni = 1:numel(rows_ni)
        ri = rows_ni(ri_ni);
        sw_row = acc_sw_rows(ri,:);
        pre_ = mean(sw_row(1:preN),'omitnan');
        post_ = mean(sw_row(preN+1:preN+5),'omitnan');
        cost_sw(ri_ni) = pre_ - post_;
    end
    cost_sw = cost_sw(~isnan(cost_sw));
    if isempty(cost_sw), continue; end
    bar(ax_E9a, ni, mean(cost_sw,'omitnan'), 0.55, 'FaceColor',CLR_NPP(ni,:), ...
        'EdgeColor','none','FaceAlpha',0.8,'HandleVisibility','off');
    errorbar(ax_E9a, ni, mean(cost_sw,'omitnan'), sem_lc(cost_sw), ...
        'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');
    text(ax_E9a,ni,-0.01,sprintf('n=%d',numel(cost_sw)), ...
        'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
end
set(ax_E9a,'XTick',1:3,'XTickLabel',{'0','1','2+'},'FontSize',9, 'TickDir','out');
xlabel(ax_E9a,'n_{prev P}'); ylabel(ax_E9a,'Switch cost (acc_{pre}−acc_{post})');
subtitle(ax_E9a,'↑ predicts larger FRN','FontSize',8,'Color',CLR_SWITCH);

% ── Panel 2: Confidence drop on switched → P300 prediction ──────────────
ax_E9b = subplot(1,4,2); hold(ax_E9b,'on');
title(ax_E9b,'P300 predictor: confidence drop at reversal','FontSize',9);
for ni = 1:3
    bin_mask_ni = meta_npp_bin == npp_bins(ni);
    conf_pre_  = mean(conf_rows(bin_mask_ni, max(1,preN-5):preN), 2,'omitnan');
    conf_post_ = mean(conf_rows(bin_mask_ni, preN+1:preN+5), 2,'omitnan');
    conf_drop  = conf_pre_ - conf_post_;
    conf_drop  = conf_drop(~isnan(conf_drop));
    if isempty(conf_drop), continue; end
    bar(ax_E9b, ni, mean(conf_drop,'omitnan'), 0.55, 'FaceColor',CLR_NPP(ni,:), ...
        'EdgeColor','none','FaceAlpha',0.8,'HandleVisibility','off');
    errorbar(ax_E9b, ni, mean(conf_drop,'omitnan'), sem_lc(conf_drop), ...
        'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');
end
set(ax_E9b,'XTick',1:3,'XTickLabel',{'0','1','2+'},'FontSize',9, 'TickDir','out');
xlabel(ax_E9b,'n_{prev P}'); ylabel(ax_E9b,'Confidence drop at reversal');
subtitle(ax_E9b,'↑ predicts larger P300 (context update)','FontSize',8,'Color',[0.40 0.25 0.65]);

% ── Panel 3: RT increase on switched → theta prediction ──────────────────
ax_E9c = subplot(1,4,3); hold(ax_E9c,'on');
title(ax_E9c,'Frontal θ predictor: RT on switched vs maintained','FontSize',9);

for bt_i = 1:2
    bt_mask_e9 = ternary_lc(bt_i==1, is_D_row, is_P_row);
    clr_e9 = ternary_lc(bt_i==1, CLR_D, CLR_P);
    lbl_e9 = ternary_lc(bt_i==1, 'Det','Prob');

    rt_post_this = mean(rt_rows(bt_mask_e9, preN+1:preN+10),2,'omitnan')*1000;
    rt_post_this = rt_post_this(~isnan(rt_post_this));

    errorbar(ax_E9c, bt_i, mean(rt_post_this,'omitnan'), sem_lc(rt_post_this), ...
        'o','Color',clr_e9,'MarkerFaceColor',clr_e9,'MarkerSize',10,...
        'LineWidth',2,'DisplayName',lbl_e9);
end
set(ax_E9c,'XTick',[1 2],'XTickLabel',{'Det','Prob'}, 'TickDir','out');
xlabel(ax_E9c,'Block type'); ylabel(ax_E9c,'RT post-reversal (ms)');
subtitle(ax_E9c,'↑ RT predicts ↑ frontal theta power','FontSize',8,'Color',[0.85 0.65 0.00]);
legend(ax_E9c,'Box','off','FontSize',8);

% ── Panel 4: Schematic hypothesis matrix ──────────────────────────────────
ax_E9d = subplot(1,4,4); hold(ax_E9d,'off');
text(ax_E9d,0.5,0.95,'Predicted ERP patterns', ...
    'HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold','FontSize',10);
ERP_rows = {'FRN/RewP','P300','Frontal θ','FP-PLV','FS-PLV'};
ERP_cols = {'Switch>Maint','P>D block','↑n_{prev P}'};
ERP_predictions = {'+' '+' '−';    % FRN: switch>maint, P>D, attenuated by noise history
                   '+' '+' '~';    % P300: switch, context update, indep of noise history
                   '+' '+' '~';    % Theta: conflict on switch, P blocks
                   '+' '~' '−';    % PLV_fp: frontal-parietal coupling, switch/reversal
                   '+' '~' '~'};   % PLV_fs: frontal-sensorimotor, response conflict
axis(ax_E9d,'off');
y_start = 0.82;
for ri_t = 1:numel(ERP_rows)
    text(ax_E9d, 0.02, y_start - (ri_t-1)*0.14, ERP_rows{ri_t}, ...
        'FontSize',9,'FontWeight','bold');
    for ci_t = 1:numel(ERP_cols)
        val = ERP_predictions{ri_t,ci_t};
        if strcmp(val,'+'), clr_t=[0 0.5 0];
        elseif strcmp(val,'−'), clr_t=[0.7 0 0];
        else, clr_t=[0.5 0.5 0.5]; end
        text(ax_E9d, 0.3 + (ci_t-1)*0.25, y_start - (ri_t-1)*0.14, val, ...
            'Color',clr_t,'FontSize',9,'FontWeight','bold','HorizontalAlignment','center');
    end
end
for ci_t = 1:numel(ERP_cols)
    text(ax_E9d, 0.3+(ci_t-1)*0.25, y_start+0.08, ERP_cols{ci_t}, ...
        'FontSize',7,'HorizontalAlignment','center','Rotation',15,'Color',[0.3 0.3 0.3]);
end
title(ax_E9d,'ERP hypothesis matrix (+ ↑, − ↓, ~ null)','FontSize',9);

annotation('textbox',[0.01 0.01 0.98 0.05],'String', ...
    ['Fig E9: EEG–behaviour bridge. Left panels show behavioural quantities that drive ERP predictions. '...
     'FRN amplitude increases with genuine PE on switched stimuli (Holroyd & Coles 2002); '...
     'n_{prev P} calibration reduces FRN by dampening expected/unexpected uncertainty integration (Yu & Dayan 2005). '...
     'P300 (context update; Polich 2007) predicted to be independent of noise history. '...
     'Right panel: hypothesis matrix for each ERP component × contrast (+: predicted increase).'],...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig_E9, fullfile(outpath,'figE9_EEG_behaviour_bridge.pdf'));
saveas(fig_E9, fullfile(outpath,'figE9_EEG_behaviour_bridge.png'));
fprintf('Fig E9 saved.\n');


% =========================================================================
%% FIG E10 — BLOCK NUMBER CURVES STRATIFIED BY BLOCK_TYPE
%
% Adds the block_type control requested on top of Fig E1/E2. The original
% Fig E1 pooled D/P when drawing block-number reversal-aligned profiles;
% this figure keeps block number on the columns and block_type on the rows.
% This makes it easier to see whether an apparent sequential-block effect is
% actually a deterministic/probabilistic-block composition effect.
% =========================================================================
fprintf('--- Building Fig E10: Block number curves split by block_type ---\n');

fig_E10 = figure('Position',[40 40 1600 620]);
sgtitle({'Block-number learning curves stratified by current block type', ...
    'Control plot: separates sequential experience from D/P block composition'}, ...
    'FontSize',12);

bt_tags    = {'D','P'};
bt_labels  = {'Deterministic','Probabilistic'};
bt_masks   = {is_D_row, is_P_row};
bt_clrs    = {CLR_D, CLR_P};

for bt_i = 1:2
    for bk = 1:MAX_BLOCK
        ax = subplot(2, MAX_BLOCK, (bt_i-1)*MAX_BLOCK + bk); hold(ax,'on');
        row_mask = (meta_block == bk) & bt_masks{bt_i};
        mat = acc_rows(row_mask,:);
        if isempty(mat) || all(isnan(mat(:)))
            text(ax,0.5,0.5,'No data','Units','normalized','HorizontalAlignment','center');
        else
            plot_ribbon_lc(ax, rel_ax, mat, bt_clrs{bt_i}, ...
                ternary_lc(bt_i==2,'--','-'), sprintf('n=%d',size(mat,1)));
            pre_  = mean(mat(:,1:preN),2,'omitnan');
            post_ = mean(mat(:,preN+1:preN+10),2,'omitnan');
            cost_ = pre_ - post_;
            cost_ = cost_(~isnan(cost_));
            if ~isempty(cost_)
                text(ax,0.02,0.97,sprintf('cost=%.2f',mean(cost_,'omitnan')), ...
                    'Units','normalized','VerticalAlignment','top','FontSize',8,'BackgroundColor','w');
            end
        end
        xline(ax,0,'k--','LineWidth',1.2,'HandleVisibility','off');
        yline(ax,0.5,'k:','HandleVisibility','off');
        patch(ax,[-preN 0 0 -preN],[0 0 1 1],[0.9 0.9 0.9], ...
            'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
        title(ax,sprintf('%s | block %d',bt_labels{bt_i},bk),'FontSize',9);
        if bk==1, ylabel(ax,'P(correct)'); end
        if bt_i==2, xlabel(ax,'Trial relative to reversal'); end
        xlim(ax,[-preN postN-1]); ylim(ax,[0.2 1]);
        legend(ax,'Box','off','FontSize',7,'Location','southeast');
    end
end

annotation('textbox',[0.01 0.01 0.98 0.04],'String', ...
    ['Fig E10: Same reversal-aligned block-number logic as Fig E1, but stratified by current block_type. '...
     'Use this to check whether sequential block effects are confounded by the D/P composition of each block number.'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig_E10, fullfile(outpath,'figE10_block_number_by_blocktype.pdf'));
saveas(fig_E10, fullfile(outpath,'figE10_block_number_by_blocktype.png'));
fprintf('Fig E10 saved.\n');

% =========================================================================
%% FIG E11 — STIM_CONFIG CONTROL PLOTS
%
% Adds the stim_config control requested on top of the existing analyses.
% Stimulus configuration can make some mappings easier/harder to learn, so
% these panels show whether reversal cost / recovery / switched-stimulus cost
% differ by configuration, and whether those differences interact with
% block number or block_type.
% =========================================================================
fprintf('--- Building Fig E11: stim_config control plots ---\n');

cfg_known = meta_stim_config ~= "unknown" & meta_stim_config ~= "";
cfg_preferred = ["a";"b";"c";"d";"e"];
cfg_levels = cfg_preferred(ismember(cfg_preferred, unique(meta_stim_config(cfg_known))));
extra_cfg = setdiff(unique(meta_stim_config(cfg_known)), cfg_preferred, 'stable');
cfg_levels = [cfg_levels; extra_cfg(:)];

fig_E11 = figure('Position',[40 40 1500 820]);
sgtitle({'Stimulus-configuration control analyses', ...
    'Checks whether a/b/c/d/e stimulus mappings confound learning difficulty'}, ...
    'FontSize',12);

if isempty(cfg_levels)
    ax = axes(fig_E11); %#ok<LAXES>
    axis(ax,'off');
    text(ax,0.5,0.55,{'No stim_config metadata found in all_trial_data or group_T.', ...
        'The script will still run; add a stim_config/config column to enable these controls.'}, ...
        'HorizontalAlignment','center','FontSize',12);
else
    cfg_cols = lines(numel(cfg_levels));

    % Panel 1: count matrix by stim_config and block number
    ax1 = subplot(2,3,1); hold(ax1,'on');
    count_mat = zeros(numel(cfg_levels), MAX_BLOCK);
    for ci = 1:numel(cfg_levels)
        for bk = 1:MAX_BLOCK
            count_mat(ci,bk) = sum(meta_stim_config == cfg_levels(ci) & meta_block == bk);
        end
    end
    imagesc(ax1, count_mat);
    colorbar(ax1);
    set(ax1,'XTick',1:MAX_BLOCK,'YTick',1:numel(cfg_levels), ...
        'YTickLabel',cellstr(cfg_levels), 'TickDir','out');
    xlabel(ax1,'Block number'); ylabel(ax1,'stim\_config');
    title(ax1,'Observations per config × block','FontSize',10);

    % Panel 2: reversal cost by stim_config
    ax2 = subplot(2,3,2); hold(ax2,'on');
    title(ax2,'Reversal cost by stim\_config','FontSize',10);
    for ci = 1:numel(cfg_levels)
        m = meta_stim_config == cfg_levels(ci) & ~isnan(cost_vec);
        vals = cost_vec(m);
        if isempty(vals), continue; end
        bar(ax2,ci,mean(vals,'omitnan'),0.55,'FaceColor',cfg_cols(ci,:), ...
            'EdgeColor','none','FaceAlpha',0.75,'HandleVisibility','off');
        errorbar(ax2,ci,mean(vals,'omitnan'),sem_lc(vals),'k.', ...
            'LineWidth',1.5,'CapSize',6,'HandleVisibility','off');
        scatter(ax2,ci + 0.18*(rand(numel(vals),1)-0.5), vals, 16, ...
            [0.25 0.25 0.25],'filled','MarkerFaceAlpha',0.25,'HandleVisibility','off');
        text(ax2,ci,0.01,sprintf('n=%d',numel(vals)), ...
            'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
    end
    yline(ax2,0,'k:','HandleVisibility','off');
    set(ax2,'XTick',1:numel(cfg_levels),'XTickLabel',cellstr(cfg_levels), 'TickDir','out');
    xlabel(ax2,'stim\_config'); ylabel(ax2,'Cost: acc_{pre} - acc_{post}');

    % Panel 3: reversal cost by stim_config and block_type
    ax3 = subplot(2,3,3); hold(ax3,'on');
    title(ax3,'Cost by stim\_config × block\_type','FontSize',10);
    for ci = 1:numel(cfg_levels)
        for bt_i = 1:2
            x = ci + ternary_lc(bt_i==1,-0.16,0.16);
            m = meta_stim_config == cfg_levels(ci) & bt_masks{bt_i} & ~isnan(cost_vec);
            vals = cost_vec(m);
            if isempty(vals), continue; end
            errorbar(ax3,x,mean(vals,'omitnan'),sem_lc(vals), ...
                ternary_lc(bt_i==1,'o','s'), 'Color',bt_clrs{bt_i}, ...
                'MarkerFaceColor',bt_clrs{bt_i},'MarkerSize',7, ...
                'LineWidth',1.5,'HandleVisibility','off');
            text(ax3,x,mean(vals,'omitnan'),sprintf(' %d',numel(vals)), ...
                'FontSize',6,'Color',[0.3 0.3 0.3]);
        end
    end
    plot(ax3,NaN,NaN,'o','Color',CLR_D,'MarkerFaceColor',CLR_D,'DisplayName','D');
    plot(ax3,NaN,NaN,'s','Color',CLR_P,'MarkerFaceColor',CLR_P,'DisplayName','P');
    yline(ax3,0,'k:','HandleVisibility','off');
    set(ax3,'XTick',1:numel(cfg_levels),'XTickLabel',cellstr(cfg_levels), 'TickDir','out');
    xlabel(ax3,'stim\_config'); ylabel(ax3,'Reversal cost');
    legend(ax3,'Box','off','FontSize',8,'Location','best');

    % Panel 4: block-number trajectory within stim_config
    ax4 = subplot(2,3,4); hold(ax4,'on');
    title(ax4,'Cost over block number by stim\_config','FontSize',10);
    for ci = 1:numel(cfg_levels)
        y = NaN(1,MAX_BLOCK); se = NaN(1,MAX_BLOCK);
        for bk = 1:MAX_BLOCK
            m = meta_stim_config == cfg_levels(ci) & meta_block == bk & ~isnan(cost_vec);
            vals = cost_vec(m);
            if isempty(vals), continue; end
            y(bk) = mean(vals,'omitnan');
            se(bk) = sem_lc(vals);
        end
        ok = ~isnan(y);
        if any(ok)
            errorbar(ax4,find(ok),y(ok),se(ok),'o-', ...
                'Color',cfg_cols(ci,:),'MarkerFaceColor',cfg_cols(ci,:), ...
                'LineWidth',1.5,'DisplayName',char(cfg_levels(ci)));
        end
    end
    yline(ax4,0,'k:','HandleVisibility','off');
    set(ax4,'XTick',1:MAX_BLOCK, 'TickDir','out'); xlim(ax4,[0.5 MAX_BLOCK+0.5]);
    xlabel(ax4,'Block number'); ylabel(ax4,'Reversal cost');
    legend(ax4,'Box','off','FontSize',8,'Location','best');

    % Panel 5: switched-stimulus cost by stim_config
    ax5 = subplot(2,3,5); hold(ax5,'on');
    title(ax5,'Switched-stimulus cost by stim\_config','FontSize',10);
    sw_cost_vec = mean(acc_sw_rows(:,1:preN),2,'omitnan') - ...
                  mean(acc_sw_rows(:,preN+1:preN+10),2,'omitnan');
    for ci = 1:numel(cfg_levels)
        m = meta_stim_config == cfg_levels(ci) & ~isnan(sw_cost_vec);
        vals = sw_cost_vec(m);
        if isempty(vals), continue; end
        bar(ax5,ci,mean(vals,'omitnan'),0.55,'FaceColor',cfg_cols(ci,:), ...
            'EdgeColor','none','FaceAlpha',0.75,'HandleVisibility','off');
        errorbar(ax5,ci,mean(vals,'omitnan'),sem_lc(vals),'k.', ...
            'LineWidth',1.5,'CapSize',6,'HandleVisibility','off');
        text(ax5,ci,0.01,sprintf('n=%d',numel(vals)), ...
            'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
    end
    yline(ax5,0,'k:','HandleVisibility','off');
    set(ax5,'XTick',1:numel(cfg_levels),'XTickLabel',cellstr(cfg_levels), 'TickDir','out');
    xlabel(ax5,'stim\_config'); ylabel(ax5,'Switch cost');

    % Panel 6: recovery slope by stim_config
    ax6 = subplot(2,3,6); hold(ax6,'on');
    title(ax6,'Recovery slope by stim\_config','FontSize',10);
    for ci = 1:numel(cfg_levels)
        m = meta_stim_config == cfg_levels(ci) & ~isnan(recov_vec);
        vals = recov_vec(m);
        if isempty(vals), continue; end
        bar(ax6,ci,mean(vals,'omitnan'),0.55,'FaceColor',cfg_cols(ci,:), ...
            'EdgeColor','none','FaceAlpha',0.75,'HandleVisibility','off');
        errorbar(ax6,ci,mean(vals,'omitnan'),sem_lc(vals),'k.', ...
            'LineWidth',1.5,'CapSize',6,'HandleVisibility','off');
    end
    yline(ax6,0,'k:','HandleVisibility','off');
    set(ax6,'XTick',1:numel(cfg_levels),'XTickLabel',cellstr(cfg_levels), 'TickDir','out');
    xlabel(ax6,'stim\_config'); ylabel(ax6,'Recovery slope');
end

annotation('textbox',[0.01 0.01 0.98 0.04],'String', ...
    ['Fig E11: Stimulus-configuration controls. If a/b/c/d/e mappings are unevenly distributed across block number or block_type, '...
     'or if specific configurations show systematically higher reversal cost / lower recovery, include stim_config as a covariate in the LME.'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig_E11, fullfile(outpath,'figE11_stim_config_controls.pdf'));
saveas(fig_E11, fullfile(outpath,'figE11_stim_config_controls.png'));
fprintf('Fig E11 saved.\n');

% Optional compact covariate model for the console log, if Statistics Toolbox is available.
try
    control_T = table(categorical(meta_subj), meta_block, categorical(meta_block_type), ...
        categorical(meta_stim_config), meta_n_prev_P, cost_vec, recov_vec, ...
        'VariableNames', {'subjID','block_number','block_type','stim_config','n_prev_P','reversal_cost','recovery_slope'});
    known_cfg_rows = string(control_T.stim_config) ~= "unknown" & ~isnan(control_T.reversal_cost);
    if exist('fitlme','file') == 2 && sum(known_cfg_rows) > 10 && numel(unique(control_T.stim_config(known_cfg_rows))) > 1
        fprintf('\nControl LME: reversal_cost ~ block_number*block_type + stim_config + (1|subjID)\n');
        lme_cfg = fitlme(control_T(known_cfg_rows,:), ...
            'reversal_cost ~ block_number*block_type + stim_config + (1|subjID)');
        disp(lme_cfg);
    else
        fprintf('\nControl LME skipped: fitlme unavailable or insufficient known stim_config rows.\n');
    end
catch ME_cfg
    fprintf('\nControl LME skipped due to error: %s\n', ME_cfg.message);
end

% =========================================================================
%% PRINT SUMMARY STATISTICS
% =========================================================================
fprintf('\n=== SUMMARY STATISTICS ===\n');
fprintf('n_{prev P} distribution across all (subject × block) observations:\n');
for ni = 1:3
    fprintf('  n_prev_P=%d: %d blocks (%.1f%%)\n', npp_bins(ni), ...
        sum(meta_npp_bin==npp_bins(ni)), 100*mean(meta_npp_bin==npp_bins(ni)));
end

fprintf('\nReversal cost by block type:\n');
for bt_i = 1:2
    bt_mask_stat = is_D_row; if bt_i==2, bt_mask_stat=is_P_row; end
    c_vec = cost_vec(bt_mask_stat);
    c_vec = c_vec(~isnan(c_vec));
    [~,pv,~,st] = ttest(c_vec);
    fprintf('  %s blocks: cost=%.3f±%.3f  t(%d)=%.2f  p=%.4f\n', ...
        bt_tags{bt_i}, mean(c_vec,'omitnan'), std(c_vec,'omitnan'), ...
        st.df, st.tstat, pv);
end

fprintf('\nSwitch vs maintained accuracy (post-reversal first 10 trials):\n');
sw_post_ = mean(acc_sw_rows(:,preN+1:preN+10),2,'omitnan');
mn_post_ = mean(acc_mn_rows(:,preN+1:preN+10),2,'omitnan');
ok_sm = ~isnan(sw_post_) & ~isnan(mn_post_);
if sum(ok_sm) > 1
    [~,p_sm,~,st_sm] = ttest(sw_post_(ok_sm), mn_post_(ok_sm));
    fprintf('  Switched=%.3f  Maintained=%.3f  t(%d)=%.2f  p=%.4f\n', ...
        mean(sw_post_(ok_sm)), mean(mn_post_(ok_sm)), st_sm.df, st_sm.tstat, p_sm);
end

fprintf('\nAll figures saved to: %s\n', outpath);
fprintf('=== E_sequential_block_behaviour_plots.m complete. ===\n');


end  % loaded_all_trials guard for S6c sequential figures
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



function c = fade_lc(clr, amount)
%FADE_LC  Mix an RGB colour with white. amount=0 original, amount=1 white.
if numel(clr) > 3, clr = clr(1:3); end
c = clr + amount*(1-clr);
c = max(0,min(1,c));
end

function row = resize_row_lc(row, n)
%RESIZE_ROW_LC  Pad or trim row vector to requested length.
row = row(:)';
if numel(row) >= n
    row = row(1:n);
else
    row = [row, NaN(1,n-numel(row))];
end
end

function out = aligned_from_raw_lc(x, rev, preN, postN)
%ALIGNED_FROM_RAW_LC  Extract a reversal-aligned row from raw trial vector.
x = x(:)';
rel = -preN:(postN-1);
out = NaN(1,numel(rel));
for ii = 1:numel(rel)
    t = rev + rel(ii);
    if t >= 1 && t <= numel(x), out(ii) = x(t); end
end
end

function sw = infer_switch_stims_lc(td, b, default_sw)
%INFER_SWITCH_STIMS_LC  Detect switched stimuli from pre/post goTrial when possible.
sw = default_sw;
try
    sf = '';
    if isfield(td,'stimType'), sf = 'stimType'; elseif isfield(td,'stimID'), sf = 'stimID'; end
    if isempty(sf) || ~isfield(td,'goTrial') || ~isfield(td,'revTrial'), return; end
    rev = round(td.revTrial(b));
    nT = size(td.correct,2);
    if ~isfinite(rev) || rev <= 1 || rev >= nT, return; end
    stim_vec = td.(sf)(b,:);
    go_vec = td.goTrial(b,:);
    detected = [];
    for ss = 1:4
        pre_go = go_vec(stim_vec==ss & (1:nT)<=rev);
        post_go = go_vec(stim_vec==ss & (1:nT)>rev);
        pre_go = pre_go(~isnan(pre_go)); post_go = post_go(~isnan(post_go));
        if ~isempty(pre_go) && ~isempty(post_go) && round(mean(pre_go)) ~= round(mean(post_go))
            detected(end+1) = ss; %#ok<AGROW>
        end
    end
    if ~isempty(detected), sw = detected; end
catch
    sw = default_sw;
end
end


% =========================================================================

function plot_ribbon_lc(ax, x, mat, clr, ls, lbl)
%PLOT_RIBBON_LC  Mean ± SEM ribbon. Smoothed with 3-trial moving average.
if isempty(mat) || all(isnan(mat(:))), return; end
mn = movmean(mean(mat,1,'omitnan'), 3,'omitnan');
se = std(mat,0,1,'omitnan') ./ sqrt(max(sum(~isnan(mat),1),1));
fill(ax,[x,fliplr(x)],[mn+se,fliplr(mn-se)],clr, ...
    'FaceAlpha',0.18,'EdgeColor','none','HandleVisibility','off');
plot(ax,x,mn,'Color',clr,'LineWidth',2,'LineStyle',ls,'DisplayName',lbl);
end

function s = sem_lc(x)
%SEM_LC  Standard error of the mean, ignoring NaN.
x = x(~isnan(x));
if numel(x) < 2, s = NaN; return; end
s = std(x,'omitnan') / sqrt(numel(x));
end

function add_sig_bracket(ax, x1, x2, y_top, p_val, label_str)
%ADD_SIG_BRACKET  Significance bracket above bars.
if isnan(p_val), return; end
if     p_val < 0.001, sig_str = '***';
elseif p_val < 0.01,  sig_str = '**';
elseif p_val < 0.05,  sig_str = '*';
else,                  sig_str = 'ns';
end
y_bar = y_top * 0.98;
line(ax,[x1 x1 x2 x2],[y_bar*0.97 y_bar y_bar y_bar*0.97], ...
    'Color','k','LineWidth',0.8,'HandleVisibility','off');
if ~isempty(label_str)
    text(ax, mean([x1 x2]), y_bar*1.01, sprintf('%s %s',label_str,sig_str), ...
        'HorizontalAlignment','center','FontSize',7,'HandleVisibility','off');
else
    text(ax, mean([x1 x2]), y_bar*1.01, sig_str, ...
        'HorizontalAlignment','center','FontSize',9,'HandleVisibility','off');
end
end


function cfg = get_block_stim_config_lc(td, group_T, sn_char, b)
%GET_BLOCK_STIM_CONFIG_LC  Return block-level stimulus configuration label.
% Searches common field/column names in all_trial_data first, then group_T.
cfg = 'unknown';
field_candidates = {'stim_config'};

for ii = 1:numel(field_candidates)
    fn = field_candidates{ii};
    if isfield(td, fn)
        val = fetch_block_value_lc(td.(fn), b);
        cfg_try = normalize_stim_config_lc(val);
        if ~strcmp(cfg_try, 'unknown')
            cfg = cfg_try;
            return;
        end
    end
end

try
    if istable(group_T) && ismember('subjID', group_T.Properties.VariableNames) && ...
            ismember('block', group_T.Properties.VariableNames)
        rows = string(group_T.subjID) == string(sn_char) & double(group_T.block) == double(b);
        for ii = 1:numel(field_candidates)
            fn = field_candidates{ii};
            if ismember(fn, group_T.Properties.VariableNames) && any(rows)
                vals = group_T.(fn)(rows);
                if iscell(vals) || isstring(vals) || iscategorical(vals)
                    vals = vals(~ismissing(string(vals)));
                elseif isnumeric(vals)
                    vals = vals(~isnan(vals));
                end
                if ~isempty(vals)
                    cfg_try = normalize_stim_config_lc(vals(1));
                    if ~strcmp(cfg_try, 'unknown')
                        cfg = cfg_try;
                        return;
                    end
                end
            end
        end
    end
catch
    cfg = 'unknown';
end
end

function val = fetch_block_value_lc(x, b)
%FETCH_BLOCK_VALUE_LC  Safely extract one block's value from many MATLAB shapes.
val = [];
try
    if iscell(x)
        if numel(x) >= b, val = x{b}; else, val = x{1}; end
    elseif isstring(x) || iscategorical(x)
        if numel(x) >= b, val = x(b); else, val = x(1); end
    elseif isnumeric(x) || islogical(x)
        if isvector(x)
            if numel(x) >= b, val = x(b); else, val = x(1); end
        elseif size(x,1) >= b
            row = x(b,:); row = row(~isnan(row));
            if ~isempty(row), val = row(1); end
        end
    elseif ischar(x)
        if size(x,1) > 1 && size(x,1) >= b
            val = strtrim(x(b,:));
        elseif isrow(x) && numel(x) >= b && all(ismember(lower(x(~isspace(x))), 'abcde'))
            val = x(b);
        else
            val = x;
        end
    end
catch
    val = [];
end
end

function cfg = normalize_stim_config_lc(val)
%NORMALIZE_STIM_CONFIG_LC  Convert raw config field to 'a'...'e' where possible.
cfg = 'unknown';
if isempty(val), return; end
try
    if iscell(val), val = val{1}; end
    if iscategorical(val), val = string(val); end
    if isnumeric(val) || islogical(val)
        val = double(val);
        if isscalar(val) && isfinite(val) && val >= 1 && val <= 5
            cfg = char('a' + round(val) - 1);
            return;
        end
    end
    str = lower(strtrim(char(string(val))));
    if isempty(str) || any(strcmp(str, {'nan','na','none','missing','unknown'})), return; end
    tok = regexp(str, '([abcde])\s*$', 'tokens', 'once');
    if ~isempty(tok), cfg = tok{1}; end
catch
    cfg = 'unknown';
end
end

function out = ternary_lc(cond, a, b)
%TERNARY_LC  Inline ternary operator.
if cond, out = a; else, out = b; end
end