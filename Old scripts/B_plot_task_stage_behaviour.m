% ==========================================================================
% plot_task_stages.m  —  DEBUGGED VERSION
%
% Stage-based performance, confidence, and RT plots.
% Designed to append to plot_RL_params_and_transitions.m, or run standalone.
%
% BUGS FIXED FROM ORIGINAL
% ─────────────────────────
%  BUG 1 (×5): BT_MATS = {...} pseudo-loop — invalid MATLAB syntax.
%    No 'for', no loop variable assignment, no 'end'. Variables mat_bt,
%    clr_bt, lbl_bt, xoff were never defined. The body code dangled after
%    an unclosed cell array literal.
%    FIX: replaced each occurrence with an explicit 'for bti = 1:2' loop.
%
%  BUG 2: 'for [mat_bt,...] = deal_pairs(...)' in Fig S5 — invalid syntax.
%    deal_pairs() at bottom is an empty placeholder.
%    FIX: replaced with explicit for bti = 1:2 loop, same pattern as Bug 1.
%
%  BUG 3: T.stay_choice and T.prevCorrect accessed without existence check.
%    If these columns are absent from group_T the script errors immediately.
%    FIX: guard with ismember() check; Fig S5 shows an informative message
%    if the columns are missing rather than crashing.
%
%  BUG 4: '[~, rows_ok] = find(~isnan(mat))' in ANOVA section.
%    Called find() with two outputs on a matrix, returning row/col indices.
%    rows_ok (column indices) was never used; the useful line
%    'ok_subj = ~any(isnan(mat),2)' was already present on the next line.
%    FIX: removed the dead find() call entirely.
%
%  BUG 5: anova1 table index wrong.
%    anova1 table columns are {Source,SS,df,MS,F,Prob>F}.
%    Original used tbl{2,5} for F and had no p extraction — now fixed to
%    use the second return value from anova1 directly.
%
% TASK STAGES (defined per block, relative to reversal trial)
% ───────────────────────────────────────────────────────────
%   LN  — Learning Naive:    first 20 trials of block
%   LE  — Learning Expert:   20 trials immediately before reversal
%   RN  — Reversal Naive:    20 trials immediately after reversal
%   RE  — Reversal Expert:   last 20 trials of block
%
% FIGURES PRODUCED
% ────────────────
%  Fig S1 — Per-stage means: accuracy, hit, FA, CR, miss rates
%  Fig S2 — True-FB win/loss rates per stage (P blocks only)
%  Fig S3 — Confidence per stage × block type
%  Fig S4 — Reaction time per stage × block type
%  Fig S5 — Stay behaviour per stage × block type
%  Fig S6 — Stage profiles by transition type (D→D / D→P / P→D / P→P)
%  Fig S7 — Stage × block-type heatmap (all measures in one view)
%
% EXPECTED WORKSPACE
%   group_T   — table from extract_revaligned_alltrialdata_v6.m
%               columns used: subjID, block, trial, revTrial, block_type,
%                             correct, perceivedCorrect, confidence, RT,
%                             goTrial, trueFB, prev_block_type
%               optional:     stay_choice, prevCorrect  (for Fig S5)
% ==========================================================================

%% ─── SETTINGS ─────────────────────────────────────────────────────────────
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'off');

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
% if ~exist(outpath,'dir'), mkdir(outpath); end
cd(base_path)



load(fullfile(base_path, 'Data', 'behav_table_June2026.mat'))


%% ─── PREPROCESS group_T ───────────────────────────────────────────────────
group_T.trial             = double(group_T.trial);
group_T.block             = double(group_T.block);
group_T.revTrial          = double(group_T.revTrial);
group_T.correct           = double(group_T.correct);
group_T.confidence        = double(group_T.confidence);
group_T.RT                = double(group_T.RT);
% Remove implausibly slow RTs (>1 s)
group_T.RT(group_T.RT > 1) = NaN;
group_T.trueFB            = double(group_T.trueFB);
group_T.goTrial           = double(group_T.goTrial);
group_T.perceivedCorrect  = double(group_T.perceivedCorrect);
group_T.block_type        = string(group_T.block_type);
group_T.prev_block_type   = string(group_T.prev_block_type);
group_T.subjID            = string(group_T.subjID);

% Normalise block_type: 'V' (legacy visual probabilistic) → 'P'
group_T.block_type_clean  = group_T.block_type;
group_T.block_type_clean(group_T.block_type == "V") = "P";

% Transition label: prev_block_type → current block_type_clean
group_T.transition = group_T.prev_block_type + "→" + group_T.block_type_clean;
group_T.transition(group_T.prev_block_type == "NaN" | ...
                   group_T.prev_block_type == "")  = "first";

% Hit / FA / Miss / CR
group_T.isHit  = double(group_T.correct==1 & group_T.goTrial==1);
group_T.isFA   = double(group_T.correct==0 & group_T.goTrial==0);
group_T.isMiss = double(group_T.correct==0 & group_T.goTrial==1);
group_T.isCR   = double(group_T.correct==1 & group_T.goTrial==0);

% True-FB decomposition
group_T.trueFB_win  = double(group_T.trueFB==1 & group_T.perceivedCorrect==1);
group_T.trueFB_loss = double(group_T.trueFB==1 & group_T.perceivedCorrect==0);
group_T.falseFB     = double(group_T.trueFB==0);

%% ─── ASSIGN STAGE PER TRIAL ─────────────────────────────────────────────────
% Priority: RE > RN > LE > LN so boundary trials belong to the later stage.
nRows = height(group_T);
group_T.stage = repmat("", nRows, 1);

subj_block_pairs = unique(group_T(:,{'subjID','block'}), 'rows');

for pb = 1:height(subj_block_pairs)
    sn  = subj_block_pairs.subjID(pb);
    blk = subj_block_pairs.block(pb);

    mask = group_T.subjID == sn & group_T.block == blk;
    rows = find(mask);
    if isempty(rows), continue; end

    rev_t  = group_T.revTrial(rows(1));
    trials = group_T.trial(rows);
    max_t  = max(trials);

    if isnan(rev_t) || ~isfinite(rev_t)
        % No reversal info — assign LN and RE only
        for ri = 1:numel(rows)
            t = trials(ri);
            if t <= STAGE_WIN
                group_T.stage(rows(ri)) = "LN";
            elseif t > max_t - STAGE_WIN
                group_T.stage(rows(ri)) = "RE";
            end
        end
        continue;
    end

    rev_t = round(rev_t);

    LN_start = 1;           LN_end = min(STAGE_WIN, rev_t-1);
    LE_start = max(1, rev_t-STAGE_WIN); LE_end = rev_t-1;
    RN_start = rev_t;       RN_end = min(max_t, rev_t+STAGE_WIN-1);
    RE_start = max(1, max_t-STAGE_WIN+1); RE_end = max_t;

    for ri = 1:numel(rows)
        t = trials(ri);
        if     t >= RE_start && t <= RE_end,   group_T.stage(rows(ri)) = "RE";
        elseif t >= RN_start && t <= RN_end,   group_T.stage(rows(ri)) = "RN";
        elseif t >= LE_start && t <= LE_end,   group_T.stage(rows(ri)) = "LE";
        elseif t >= LN_start && t <= LN_end,   group_T.stage(rows(ri)) = "LN";
        end
    end
end

% Only keep labelled rows for stage analyses
T = group_T(group_T.stage ~= "", :);
fprintf('Valid stage labels: %d / %d rows\n', height(T), nRows);

%% ─── BUILD SUBJECT × STAGE MATRICES ──────────────────────────────────────
all_subjs = unique(T.subjID);
N_subj    = numel(all_subjs);
n_stg     = numel(STAGES);
x_pos     = 1:n_stg;

measures = {'correct','isHit','isFA','isMiss','isCR', ...
            'confidence','RT','trueFB_win','trueFB_loss','falseFB'};

SMAT = struct();
for mi = 1:numel(measures)
    SMAT.(measures{mi}).all = NaN(N_subj, n_stg);
    SMAT.(measures{mi}).D   = NaN(N_subj, n_stg);
    SMAT.(measures{mi}).P   = NaN(N_subj, n_stg);
end

for si = 1:N_subj
    sn = all_subjs(si);
    Ts = T(T.subjID == sn, :);
    for stgi = 1:n_stg
        sg    = STAGES{stgi};
        Tsg   = Ts(Ts.stage == sg, :);
        Tsg_D = Tsg(Tsg.block_type_clean == "D", :);
        Tsg_P = Tsg(Tsg.block_type_clean == "P", :);
        for mi = 1:numel(measures)
            m = measures{mi};
            if height(Tsg)   > 0, SMAT.(m).all(si,stgi) = mean(Tsg.(m),   'omitnan'); end
            if height(Tsg_D) > 0, SMAT.(m).D(si,stgi)   = mean(Tsg_D.(m), 'omitnan'); end
            if height(Tsg_P) > 0, SMAT.(m).P(si,stgi)   = mean(Tsg_P.(m), 'omitnan'); end
        end
    end
end

is_kh = cellfun(@(s) startsWith(s,'Ox'), cellstr(all_subjs));
is_rr = cellfun(@(s) startsWith(s,'Nc'), cellstr(all_subjs));
fprintf('Stage matrices built: %d subjects × %d stages.\n', N_subj, n_stg);

% =========================================================================
%% FIG S1 — PER-STAGE ACCURACY AND RESPONSE RATES
% =========================================================================
fig_s1 = figure('Position',[50 50 1500 650]);
sgtitle('Task stage performance: accuracy and response rates (subject means ± SEM)','FontSize',12);

meas_s1 = {'correct','isHit','isFA','isMiss','isCR'};
ylbl_s1 = {'P(correct)','P(Hit)','P(FA)','P(Miss)','P(CR)'};
ylim_s1 = {[0 1],[0 1],[0 1],[0 1],[0 1]};

for mi = 1:numel(meas_s1)
    m = meas_s1{mi};

    % ── Top row: all blocks pooled ─────────────────────────────────────
    ax = subplot(2, numel(meas_s1), mi); hold(ax,'on');
    title(ax, ylbl_s1{mi}, 'FontSize',10);

    mat_all = SMAT.(m).all;
    gp_mean = mean(mat_all, 1, 'omitnan');
    gp_se   = stage_sem(mat_all);
    n_con   = sum(~isnan(mat_all), 1);

    for stgi = 1:n_stg
        bar(ax, stgi, gp_mean(stgi), 0.55, 'FaceColor',CLR_STGS{stgi}, ...
            'EdgeColor','none','FaceAlpha',0.75);
    end
    errorbar(ax, x_pos, gp_mean, gp_se, 'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');
    for stgi = 1:n_stg
        v  = mat_all(:,stgi); v = v(~isnan(v));
        jx = stgi + 0.18*(rand(numel(v),1)-0.5);
        scatter(ax, jx, v, 18, [0.3 0.3 0.3],'filled','MarkerFaceAlpha',0.4,'HandleVisibility','off');
        text(ax, stgi, ylim_s1{mi}(2)*0.03, sprintf('n=%d',n_con(stgi)), ...
            'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
    end
    set(ax,'XTick',x_pos,'XTickLabel',STAGES,'FontSize',9);
    ylabel(ax, ylbl_s1{mi}); ylim(ax, ylim_s1{mi});
    add_pairwise_brackets(ax, mat_all, x_pos, ylim_s1{mi}(2));

    % ── Bottom row: D vs P comparison ─────────────────────────────────
    % BUG FIX: was an invalid BT_MATS = {...} pseudo-loop.
    % Replaced with an explicit for loop over the two block types.
    ax2 = subplot(2, numel(meas_s1), mi + numel(meas_s1)); hold(ax2,'on');
    title(ax2, sprintf('%s: D vs P',ylbl_s1{mi}), 'FontSize',9);

    bt_data  = {SMAT.(m).D,  SMAT.(m).P};
    bt_clrs  = {CLR_D,       CLR_P};
    bt_lbls  = {'Det',       'Prob'};
    bt_xoffs = [-0.18,        0.18];

    for bti = 1:2
        mat_bt = bt_data{bti};
        clr_bt = bt_clrs{bti};
        lbl_bt = bt_lbls{bti};
        xoff   = bt_xoffs(bti);

        gm  = mean(mat_bt, 1, 'omitnan');
        gs  = stage_sem(mat_bt);
        n_b = sum(~isnan(mat_bt), 1);
        xp  = x_pos + xoff;

        bar(ax2, xp, gm, 0.3, 'FaceColor',clr_bt,'EdgeColor','none', ...
            'FaceAlpha',0.75,'DisplayName',lbl_bt);
        errorbar(ax2, xp, gm, gs, 'k.','LineWidth',1.2,'CapSize',4,'HandleVisibility','off');
        for stgi = 1:n_stg
            text(ax2, xp(stgi), ylim_s1{mi}(2)*0.03, sprintf('n=%d',n_b(stgi)), ...
                'HorizontalAlignment','center','FontSize',6,'Color',clr_bt*0.7);
        end
    end

    set(ax2,'XTick',x_pos,'XTickLabel',STAGES,'FontSize',9);
    ylabel(ax2, ylbl_s1{mi}); ylim(ax2, ylim_s1{mi});
    legend(ax2,'Box','off','FontSize',8,'Location','best');
end

saveas(fig_s1, fullfile(outpath,'figS1_stage_accuracy.pdf'));
fprintf('Fig S1 saved.\n');

% =========================================================================
%% FIG S2 — TRUE-FB WIN / LOSS / FALSE-FB RATES (P blocks only)
% =========================================================================
fig_s2 = figure('Position',[50 50 1100 500]);
sgtitle({'True-FB win/loss and false-feedback rate per stage', ...
    '(P blocks only — participants unaware of false FB)'},'FontSize',12);

meas_s2  = {'trueFB_win','trueFB_loss','falseFB'};
ylbl_s2  = {'P(trueFB win)','P(trueFB loss)','P(false feedback)'};
y_expect = {0.64, 0.16, 0.20};

for mi = 1:numel(meas_s2)
    m   = meas_s2{mi};
    ax  = subplot(1,3,mi); hold(ax,'on');
    title(ax, ylbl_s2{mi},'FontSize',10);

    mat_P = SMAT.(m).P;
    gm    = mean(mat_P,1,'omitnan');
    gs    = stage_sem(mat_P);
    n_c   = sum(~isnan(mat_P),1);

    for stgi = 1:n_stg
        bar(ax, stgi, gm(stgi), 0.55, 'FaceColor',CLR_STGS{stgi},'EdgeColor','none','FaceAlpha',0.75);
    end
    errorbar(ax, x_pos, gm, gs, 'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');
    yline(ax, y_expect{mi}, 'k--','LineWidth',1.2,'HandleVisibility','off');
    text(ax, 0.02, y_expect{mi}+0.02, sprintf('Expected=%.2f',y_expect{mi}), ...
        'FontSize',7,'Color',[0.3 0.3 0.3]);

    for stgi = 1:n_stg
        v  = mat_P(:,stgi); v = v(~isnan(v));
        jx = stgi + 0.18*(rand(numel(v),1)-0.5);
        scatter(ax, jx, v, 18, [0.3 0.3 0.3],'filled','MarkerFaceAlpha',0.4,'HandleVisibility','off');
        text(ax, stgi, 0.02, sprintf('n=%d',n_c(stgi)), ...
            'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
    end

    set(ax,'XTick',x_pos,'XTickLabel',STAGES,'FontSize',9);
    ylabel(ax, ylbl_s2{mi}); ylim(ax,[0 1]);
    add_pairwise_brackets(ax, mat_P, x_pos, 0.95);
    subtitle(ax,'P blocks only','FontSize',8,'Color',[0.5 0.5 0.5]);
end

annotation('textbox',[0.01 0.01 0.98 0.07],'String', ...
    ['trueFB win = honest feedback shown as correct. ' ...
    'trueFB loss = honest feedback shown as incorrect. ' ...
    'false feedback = misleading trial (unknown to participant). ' ...
    'Expected at p(trueFB)=0.8: win=0.64, loss=0.16, falseFB=0.20.'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig_s2, fullfile(outpath,'figS2_stage_trueFB.pdf'));
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
    bar(ax_s3a, stgi, gm(stgi), 0.55,'FaceColor',CLR_STGS{stgi},'EdgeColor','none','FaceAlpha',0.75);
end
errorbar(ax_s3a, x_pos, gm, gs,'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');
for stgi = 1:n_stg
    v  = mat_all_c(:,stgi); v = v(~isnan(v));
    jx = stgi + 0.18*(rand(numel(v),1)-0.5);
    scatter(ax_s3a, jx, v, 18, [0.3 0.3 0.3],'filled','MarkerFaceAlpha',0.4,'HandleVisibility','off');
    text(ax_s3a, stgi, 1.2, sprintf('n=%d',n_c(stgi)), ...
        'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
end
set(ax_s3a,'XTick',x_pos,'XTickLabel',STAGES,'FontSize',9);
ylabel(ax_s3a,'Confidence (1–10)'); ylim(ax_s3a,[1 10]);
add_pairwise_brackets(ax_s3a, mat_all_c, x_pos, 9.5);

% Panel 2: D vs P bars
% BUG FIX: was an invalid BT_MATS = {...} pseudo-loop.
ax_s3b = subplot(1,3,2); hold(ax_s3b,'on');
title(ax_s3b,'Confidence: D vs P blocks','FontSize',10);

bt_data  = {SMAT.confidence.D, SMAT.confidence.P};
bt_clrs  = {CLR_D,              CLR_P};
bt_lbls  = {'Det',              'Prob'};
bt_xoffs = [-0.18,               0.18];

for bti = 1:2
    mat_bt = bt_data{bti};
    clr_bt = bt_clrs{bti};
    lbl_bt = bt_lbls{bti};
    xoff   = bt_xoffs(bti);

    gm2 = mean(mat_bt,1,'omitnan');
    gs2 = stage_sem(mat_bt);
    xp  = x_pos + xoff;

    bar(ax_s3b, xp, gm2, 0.3,'FaceColor',clr_bt,'EdgeColor','none', ...
        'FaceAlpha',0.75,'DisplayName',lbl_bt);
    errorbar(ax_s3b, xp, gm2, gs2,'k.','LineWidth',1.2,'CapSize',4,'HandleVisibility','off');
end

set(ax_s3b,'XTick',x_pos,'XTickLabel',STAGES,'FontSize',9);
ylabel(ax_s3b,'Confidence (1–10)'); ylim(ax_s3b,[1 10]);
legend(ax_s3b,'Box','off','FontSize',9);

% Panel 3: spaghetti + group means
ax_s3c = subplot(1,3,3); hold(ax_s3c,'on');
title(ax_s3c,'Stage means per subject','FontSize',10);
subtitle(ax_s3c,'Lines = within-subject stage means','FontSize',8,'Color',[0.5 0.5 0.5]);
for si = 1:N_subj
    row_D = SMAT.confidence.D(si,:);
    row_P = SMAT.confidence.P(si,:);
    if ~all(isnan(row_D))
        plot(ax_s3c, x_pos, row_D,'o-','Color',[CLR_D 0.25],'MarkerSize',4,'HandleVisibility','off');
    end
    if ~all(isnan(row_P))
        plot(ax_s3c, x_pos, row_P,'s--','Color',[CLR_P 0.25],'MarkerSize',4,'HandleVisibility','off');
    end
end
plot(ax_s3c, x_pos, mean(SMAT.confidence.D,1,'omitnan'),'o-', ...
    'Color',CLR_D,'LineWidth',2.5,'MarkerSize',8,'MarkerFaceColor',CLR_D,'DisplayName','Det');
plot(ax_s3c, x_pos, mean(SMAT.confidence.P,1,'omitnan'),'s--', ...
    'Color',CLR_P,'LineWidth',2.5,'MarkerSize',8,'MarkerFaceColor',CLR_P,'DisplayName','Prob');
set(ax_s3c,'XTick',x_pos,'XTickLabel',STAGES,'FontSize',9);
ylabel(ax_s3c,'Confidence (1–10)'); ylim(ax_s3c,[1 10]);
legend(ax_s3c,'Box','off','FontSize',9,'Location','best');

saveas(fig_s3, fullfile(outpath,'figS3_stage_confidence.pdf'));
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
    bar(ax_s4a, stgi, gm(stgi), 0.55,'FaceColor',CLR_STGS{stgi},'EdgeColor','none','FaceAlpha',0.75);
end
errorbar(ax_s4a, x_pos, gm, gs,'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');
for stgi = 1:n_stg
    v  = mat_rt(:,stgi)*1000; v = v(~isnan(v));
    jx = stgi + 0.18*(rand(numel(v),1)-0.5);
    scatter(ax_s4a, jx, v, 18, [0.3 0.3 0.3],'filled','MarkerFaceAlpha',0.4,'HandleVisibility','off');
    if ~isnan(gm(stgi))
        text(ax_s4a, stgi, min(gm(~isnan(gm)))*0.85, sprintf('n=%d',n_c(stgi)), ...
            'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
    end
end
set(ax_s4a,'XTick',x_pos,'XTickLabel',STAGES,'FontSize',9);
ylabel(ax_s4a,'RT (ms)');
if any(~isnan(gm))
    ylim(ax_s4a,[0 max(gm+gs,[],'omitnan')*1.3]);
    add_pairwise_brackets(ax_s4a, mat_rt*1000, x_pos, max(gm+gs,[],'omitnan')*1.25);
end

% Panel 2: D vs P
% BUG FIX: was an invalid BT_MATS = {...} pseudo-loop.
ax_s4b = subplot(1,3,2); hold(ax_s4b,'on');
title(ax_s4b,'RT: D vs P blocks','FontSize',10);

bt_data  = {SMAT.RT.D, SMAT.RT.P};
bt_clrs  = {CLR_D,      CLR_P};
bt_lbls  = {'Det',      'Prob'};
bt_xoffs = [-0.18,       0.18];

for bti = 1:2
    mat_bt = bt_data{bti};
    clr_bt = bt_clrs{bti};
    lbl_bt = bt_lbls{bti};
    xoff   = bt_xoffs(bti);

    gm2 = mean(mat_bt,1,'omitnan') * 1000;
    gs2 = stage_sem(mat_bt) * 1000;
    xp  = x_pos + xoff;

    bar(ax_s4b, xp, gm2, 0.3,'FaceColor',clr_bt,'EdgeColor','none', ...
        'FaceAlpha',0.75,'DisplayName',lbl_bt);
    errorbar(ax_s4b, xp, gm2, gs2,'k.','LineWidth',1.2,'CapSize',4,'HandleVisibility','off');
end

set(ax_s4b,'XTick',x_pos,'XTickLabel',STAGES,'FontSize',9);
ylabel(ax_s4b,'RT (ms)');
legend(ax_s4b,'Box','off','FontSize',9);

% Panel 3: RT by correct vs incorrect
% BUG FIX: was an invalid BT_MATS = {...} pseudo-loop.
ax_s4c = subplot(1,3,3); hold(ax_s4c,'on');
title(ax_s4c,'RT: correct vs incorrect','FontSize',10);

SMAT_RT_corr   = NaN(N_subj, n_stg);
SMAT_RT_incorr = NaN(N_subj, n_stg);
for si = 1:N_subj
    sn = all_subjs(si);
    Ts = T(T.subjID == sn, :);
    for stgi = 1:n_stg
        sg  = STAGES{stgi};
        Tsg = Ts(Ts.stage == sg, :);
        SMAT_RT_corr(si,stgi)   = mean(Tsg.RT(Tsg.correct==1),'omitnan') * 1000;
        SMAT_RT_incorr(si,stgi) = mean(Tsg.RT(Tsg.correct==0),'omitnan') * 1000;
    end
end

bt_data  = {SMAT_RT_corr,        SMAT_RT_incorr};
bt_clrs  = {[0.2 0.6 0.2],       [0.8 0.2 0.2]};
bt_lbls  = {'Correct',           'Incorrect'};
bt_xoffs = [-0.18,                0.18];

for bti = 1:2
    mat_bt = bt_data{bti};
    clr_bt = bt_clrs{bti};
    lbl_bt = bt_lbls{bti};
    xoff   = bt_xoffs(bti);

    gm2 = mean(mat_bt,1,'omitnan');
    gs2 = stage_sem(mat_bt);
    xp  = x_pos + xoff;

    bar(ax_s4c, xp, gm2, 0.3,'FaceColor',clr_bt,'EdgeColor','none', ...
        'FaceAlpha',0.75,'DisplayName',lbl_bt);
    errorbar(ax_s4c, xp, gm2, gs2,'k.','LineWidth',1.2,'CapSize',4,'HandleVisibility','off');
end

set(ax_s4c,'XTick',x_pos,'XTickLabel',STAGES,'FontSize',9);
ylabel(ax_s4c,'RT (ms)');
legend(ax_s4c,'Box','off','FontSize',9,'Location','best');

saveas(fig_s4, fullfile(outpath,'figS4_stage_RT.pdf'));
fprintf('Fig S4 saved.\n');

% =========================================================================
%% FIG S5 — STAY BEHAVIOUR PER STAGE
% =========================================================================
% BUG FIX: T.stay_choice and T.prevCorrect accessed without existence check.
% Added ismember guard; Fig S5 shows informative message if absent.

has_stay = ismember('stay_choice',  T.Properties.VariableNames) && ...
           ismember('prevCorrect',  T.Properties.VariableNames);

if has_stay
    T_stay = T(~isnan(T.stay_choice) & ~isnan(T.prevCorrect), :);
else
    T_stay = T(false(height(T),1), :);   % empty table, same schema
    warning('plot_task_stages: stay_choice / prevCorrect not found in group_T. Fig S5 will be empty. Compute stay_choice as same-action repeat for each stimulus and add to group_T.');
end

SMAT_winstay  = NaN(N_subj, n_stg);
SMAT_losestay = NaN(N_subj, n_stg);
SMAT_stay_all = NaN(N_subj, n_stg);
SMAT_stay_D   = NaN(N_subj, n_stg);
SMAT_stay_P   = NaN(N_subj, n_stg);

if has_stay
    for si = 1:N_subj
        sn      = all_subjs(si);
        Ts_stay = T_stay(T_stay.subjID == sn, :);
        for stgi = 1:n_stg
            sg    = STAGES{stgi};
            Tsg   = Ts_stay(Ts_stay.stage == sg, :);
            if height(Tsg) == 0, continue; end

            SMAT_stay_all(si,stgi) = mean(Tsg.stay_choice,'omitnan');
            win_r  = Tsg(Tsg.prevCorrect==1, :);
            lose_r = Tsg(Tsg.prevCorrect==0, :);
            if height(win_r)  > 0, SMAT_winstay(si,stgi)  = mean(win_r.stay_choice,'omitnan');  end
            if height(lose_r) > 0, SMAT_losestay(si,stgi) = mean(lose_r.stay_choice,'omitnan'); end

            Tsg_D = Tsg(Tsg.block_type_clean == "D", :);
            Tsg_P = Tsg(Tsg.block_type_clean == "P", :);
            if height(Tsg_D) > 0, SMAT_stay_D(si,stgi) = mean(Tsg_D.stay_choice,'omitnan'); end
            if height(Tsg_P) > 0, SMAT_stay_P(si,stgi) = mean(Tsg_P.stay_choice,'omitnan'); end
        end
    end
end

fig_s5 = figure('Position',[50 50 1300 500]);
sgtitle('Stay behaviour per stage (stimulus-specific)','FontSize',12);

% Panel 1: P(stay) all
ax_s5a = subplot(1,3,1); hold(ax_s5a,'on');
title(ax_s5a,'P(stay): all trials','FontSize',10);
mat5 = SMAT_stay_all;
gm5 = mean(mat5,1,'omitnan'); gs5 = stage_sem(mat5); n5 = sum(~isnan(mat5),1);
for stgi = 1:n_stg
    bar(ax_s5a, stgi, gm5(stgi), 0.55,'FaceColor',CLR_STGS{stgi},'EdgeColor','none','FaceAlpha',0.75);
end
errorbar(ax_s5a, x_pos, gm5, gs5,'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');
for stgi = 1:n_stg
    v  = mat5(:,stgi); v = v(~isnan(v));
    jx = stgi + 0.18*(rand(numel(v),1)-0.5);
    scatter(ax_s5a, jx, v, 18, [0.3 0.3 0.3],'filled','MarkerFaceAlpha',0.4,'HandleVisibility','off');
    text(ax_s5a, stgi, 0.02, sprintf('n=%d',n5(stgi)), ...
        'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
end
if ~has_stay
    text(ax_s5a,0.5,0.5,'stay\_choice not in group\_T', ...
        'Units','normalized','HorizontalAlignment','center','FontSize',9,'Color',[0.6 0.3 0]);
end
yline(ax_s5a, 0.5,'k:','HandleVisibility','off');
set(ax_s5a,'XTick',x_pos,'XTickLabel',STAGES); ylabel(ax_s5a,'P(stay)'); ylim(ax_s5a,[0 1]);
add_pairwise_brackets(ax_s5a, mat5, x_pos, 0.95);

% Panel 2: win-stay vs lose-stay
% BUG FIX: was an invalid BT_MATS = {...} pseudo-loop.
ax_s5b = subplot(1,3,2); hold(ax_s5b,'on');
title(ax_s5b,'Win-stay vs lose-stay','FontSize',10);

bt_data  = {SMAT_winstay,      SMAT_losestay};
bt_clrs  = {[0.2 0.6 0.2],     [0.8 0.2 0.2]};
bt_lbls  = {'Win-stay',        'Lose-stay'};
bt_xoffs = [-0.18,              0.18];

for bti = 1:2
    mat_bt = bt_data{bti};
    clr_bt = bt_clrs{bti};
    lbl_bt = bt_lbls{bti};
    xoff   = bt_xoffs(bti);

    gm2 = mean(mat_bt,1,'omitnan');
    gs2 = stage_sem(mat_bt);
    xp  = x_pos + xoff;

    bar(ax_s5b, xp, gm2, 0.3,'FaceColor',clr_bt,'EdgeColor','none', ...
        'FaceAlpha',0.75,'DisplayName',lbl_bt);
    errorbar(ax_s5b, xp, gm2, gs2,'k.','LineWidth',1.2,'CapSize',4,'HandleVisibility','off');
end

yline(ax_s5b, 0.5,'k:','HandleVisibility','off');
set(ax_s5b,'XTick',x_pos,'XTickLabel',STAGES); ylabel(ax_s5b,'P(stay)'); ylim(ax_s5b,[0 1]);
legend(ax_s5b,'Box','off','FontSize',9,'Location','best');

% Panel 3: D vs P
% BUG FIX: was an invalid 'for [mat_bt,...] = deal_pairs(...)' call.
% deal_pairs() at the bottom was an empty placeholder — completely invalid.
ax_s5c = subplot(1,3,3); hold(ax_s5c,'on');
title(ax_s5c,'Stay: D vs P blocks','FontSize',10);

bt_data  = {SMAT_stay_D, SMAT_stay_P};
bt_clrs  = {CLR_D,        CLR_P};
bt_lbls  = {'Det',        'Prob'};
bt_xoffs = [-0.18,         0.18];

for bti = 1:2
    mat_bt = bt_data{bti};
    clr_bt = bt_clrs{bti};
    lbl_bt = bt_lbls{bti};
    xoff   = bt_xoffs(bti);

    gm2 = mean(mat_bt,1,'omitnan');
    gs2 = stage_sem(mat_bt);
    xp  = x_pos + xoff;

    bar(ax_s5c, xp, gm2, 0.3,'FaceColor',clr_bt,'EdgeColor','none', ...
        'FaceAlpha',0.75,'DisplayName',lbl_bt);
    errorbar(ax_s5c, xp, gm2, gs2,'k.','LineWidth',1.2,'CapSize',4,'HandleVisibility','off');
end

yline(ax_s5c, 0.5,'k:','HandleVisibility','off');
set(ax_s5c,'XTick',x_pos,'XTickLabel',STAGES); ylabel(ax_s5c,'P(stay)'); ylim(ax_s5c,[0 1]);
legend(ax_s5c,'Box','off','FontSize',9);

saveas(fig_s5, fullfile(outpath,'figS5_stage_stay.pdf'));
fprintf('Fig S5 saved.\n');

% =========================================================================
%% FIG S6 — STAGE PROFILES BY TRANSITION TYPE
% =========================================================================
fig_s6 = figure('Position',[50 50 1600 700]);
sgtitle('Stage profiles by block transition type','FontSize',12);

meas_s6  = {'correct','confidence'};
ylbl_s6  = {'P(correct)','Confidence (1–10)'};
ylims_s6 = {[0 1],[1 10]};

trans_subjs_mat = cell(numel(TRANS_TYPES), 2);
for ti = 1:numel(TRANS_TYPES)
    for mi = 1:2
        trans_subjs_mat{ti,mi} = NaN(N_subj, n_stg);
    end
    tr_str = TRANS_TYPES{ti};
    for si = 1:N_subj
        sn = all_subjs(si);
        Ts = T(T.subjID==sn & T.transition==tr_str, :);
        if height(Ts)==0, continue; end
        for stgi = 1:n_stg
            sg  = STAGES{stgi};
            Tsg = Ts(Ts.stage==sg, :);
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
        ax = subplot(2, numel(TRANS_TYPES), (mi-1)*numel(TRANS_TYPES)+ti);
        hold(ax,'on');
        if mi==1, title(ax, TRANS_TYPES{ti},'FontSize',11); end
        if ti==1, ylabel(ax, ylbl_s6{mi}); end

        mat = trans_subjs_mat{ti,mi};
        gm  = mean(mat,1,'omitnan');
        gs  = stage_sem(mat);
        n_c = sum(~isnan(mat),1);

        % Individual spaghetti (light)
        for si = 1:N_subj
            if ~all(isnan(mat(si,:)))
                plot(ax, x_pos, mat(si,:),'o-', ...
                    'Color',[TRANS_COLORS{ti} 0.15],'LineWidth',0.8,'MarkerSize',3,'HandleVisibility','off');
            end
        end
        % Group mean ribbon
        fill(ax,[x_pos,fliplr(x_pos)],[gm+gs,fliplr(gm-gs)], ...
            TRANS_COLORS{ti},'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
        plot(ax, x_pos, gm,'o-','Color',TRANS_COLORS{ti},'LineWidth',2.5,'MarkerSize',8, ...
            'MarkerFaceColor',TRANS_COLORS{ti},'DisplayName',sprintf('n=%d',max(n_c)));

        % LE→RN connector to highlight reversal cost
        if ~isnan(gm(le_idx)) && ~isnan(gm(rn_idx))
            plot(ax,[le_idx rn_idx],[gm(le_idx) gm(rn_idx)],'k-','LineWidth',2,'HandleVisibility','off');
        end

        set(ax,'XTick',x_pos,'XTickLabel',STAGE_LBLS,'XTickLabelRotation',20,'FontSize',8);
        ylim(ax, ylims_s6{mi});
        for stgi = 1:n_stg
            text(ax, stgi, ylims_s6{mi}(1)+diff(ylims_s6{mi})*0.03, ...
                sprintf('n=%d',n_c(stgi)), 'HorizontalAlignment','center', ...
                'FontSize',7,'Color',TRANS_COLORS{ti}*0.7);
        end
        legend(ax,'Box','off','FontSize',8,'Location','best');
    end
end

annotation('textbox',[0.01 0.01 0.98 0.04],'String', ...
    ['Bold line connects LE→RN to highlight reversal cost per transition. ' ...
    'n = subjects contributing to that transition×stage cell.'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig_s6, fullfile(outpath,'figS6_stage_by_transition.pdf'));
fprintf('Fig S6 saved.\n');

% =========================================================================
%% FIG S7 — HEATMAP SUMMARY
% =========================================================================
fig_s7 = figure('Position',[50 50 1400 600]);
sgtitle('Heatmap: group-mean measures × stage × block type','FontSize',12);

hmap_measures = {'correct','isHit','isFA','isMiss','isCR','confidence','RT'};
hmap_lbls     = {'Accuracy','Hit','FA','Miss','CR','Confidence','RT (s)'};
n_hm = numel(hmap_measures);

for bti = 1:2
    bt_tag = ternary_local(bti==1,'D','P');
    bt_lbl = ternary_local(bti==1,'Deterministic','Probabilistic');
    ax_hm  = subplot(1,2,bti);

    hmat = NaN(n_hm, n_stg);
    for mi = 1:n_hm
        m = hmap_measures{mi};
        if isfield(SMAT, m)
            col = mean(SMAT.(m).(bt_tag), 1, 'omitnan');
            if strcmp(m,'RT')
                r = range(col); if r==0, r=1; end
                col = (col - min(col)) / r;
            end
            hmat(mi,:) = col;
        end
    end

    imagesc(ax_hm, hmat);
    colormap(ax_hm);
    cb = colorbar(ax_hm); cb.Label.String = 'Mean (RT normalised to 0-1)';
    set(ax_hm,'XTick',1:n_stg,'XTickLabel',STAGES, ...
              'YTick',1:n_hm,'YTickLabel',hmap_lbls,'FontSize',10);
    title(ax_hm, bt_lbl,'FontSize',11);
    xlabel(ax_hm,'Stage'); ylabel(ax_hm,'Measure');

    for mi = 1:n_hm
        for stgi = 1:n_stg
            v = mean(SMAT.(hmap_measures{mi}).(bt_tag)(:,stgi),'omitnan');
            if ~isnan(v)
                text(ax_hm, stgi, mi, sprintf('%.2f',v), ...
                    'HorizontalAlignment','center','VerticalAlignment','middle', ...
                    'FontSize',8,'Color','w','FontWeight','bold');
            end
        end
    end
end

saveas(fig_s7, fullfile(outpath,'figS7_heatmap_stages.pdf'));
fprintf('Fig S7 saved.\n');

% =========================================================================
%% PRINT STAGE STATISTICS
% =========================================================================
fprintf('\n=== STAGE STATISTICS: Group means (subject-averaged) ===\n');
fprintf('%-12s  %-5s  %-8s  %-8s  %-8s  %-8s\n','Measure','Type',STAGES{:});
for mi = 1:numel(measures)
    m = measures{mi};
    for bt = {'all','D','P'}
        gm = mean(SMAT.(m).(bt{1}), 1, 'omitnan');
        fprintf('%-12s  %-5s  %-8.3f  %-8.3f  %-8.3f  %-8.3f\n', m, bt{1}, gm);
    end
end

% BUG FIX: removed dead '[~, rows_ok] = find(~isnan(mat))' call.
% ok_subj is computed correctly via any(isnan,2) on the next line.
fprintf('\n=== ONE-WAY REPEATED-MEASURES ANOVA: STAGE EFFECT ===\n');
fprintf('(Each row = a subject, each column = a stage — all blocks pooled)\n');
for mi = 1:numel(measures)
    m = measures{mi};
    mat = SMAT.(m).all;
    ok_subj = ~any(isnan(mat),2);
    if sum(ok_subj) < 3, continue; end

    % BUG FIX: anova1 returns [p, tbl, stats].
    % Table columns: {Source, SS, df, MS, F, Prob>F}.
    % Original indexed tbl_rm{2,3}/{3,3}/{2,5} incorrectly.
    % Now use the direct p return value and extract F cleanly.
    [p_rm, tbl_rm] = anova1(mat(ok_subj,:), [], 'off');
    % tbl_rm{2,5} = F statistic (row 2 = Groups, col 5 = F)
    % tbl_rm{2,3} = numerator df, tbl_rm{3,3} = denominator df
    F_stat   = tbl_rm{2,5};
    df_num   = tbl_rm{2,3};
    df_den   = tbl_rm{3,3};
    if ~isempty(F_stat) && ~isnan(F_stat)
        fprintf('%-12s  F(%d,%d)=%.2f  p=%.4f\n', m, df_num, df_den, F_stat, p_rm);
    end
end

fprintf('\nAll stage figures saved to:\n  %s\n', outpath);


% =========================================================================
%% LOCAL HELPER FUNCTIONS
% =========================================================================

function gs = stage_sem(mat)
%STAGE_SEM  Per-column SEM across subjects (rows), ignoring NaN.
n  = sum(~isnan(mat), 1);
gs = std(mat, 0, 1, 'omitnan') ./ sqrt(max(n, 1));
end


function add_pairwise_brackets(ax, mat, x_pos, y_top)
%ADD_PAIRWISE_BRACKETS  LN vs RN and LE vs RE significance brackets.
%  Only two theory-driven pairs tested to avoid inflating alpha.
pairs  = {[1 3], [2 4]};
labels = {'LN vs RN', 'LE vs RE'};
y_br   = y_top * [0.92 0.97];

for pi = 1:2
    a = pairs{pi}(1); b = pairs{pi}(2);
    va = mat(:,a); vb = mat(:,b);
    ok = ~isnan(va) & ~isnan(vb);
    if sum(ok) < 3, continue; end
    [~,p] = ttest(va(ok), vb(ok));

    if     p < 0.001, sig = '***';
    elseif p < 0.01,  sig = '**';
    elseif p < 0.05,  sig = '*';
    else,             sig = 'ns';
    end

    xa = x_pos(a); xb = x_pos(b);
    yb = y_br(pi);
    line(ax,[xa xa xb xb],[yb*0.97 yb yb yb*0.97], ...
        'Color','k','LineWidth',0.8,'HandleVisibility','off');
    text(ax, mean([xa xb]), yb*1.005, sprintf('%s %s',labels{pi},sig), ...
        'HorizontalAlignment','center','FontSize',7,'HandleVisibility','off');
end
end


function out = ternary_local(cond, a, b)
if cond, out = a; else, out = b; end
end