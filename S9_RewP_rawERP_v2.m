% P9_poster_ERP_behaviour_figures_DEBUGGED.m
% Debugged from uploaded P9_poster_ERP_behaviour_figures.m
%
% FIGURES PRODUCED
% ----------------
%   Fig 1  RewP difference waves per stage × block type (from frn_rewp_stage_table)
%   Fig 2  Grand-average ERP by condition (correct/incorrect/false-correct/
%          false-incorrect) in D and P panels; RewP difference wave below each.
%   Fig 3  First-10-trials accuracy per block: D vs P axes.
%   Fig 4  Mock raw EEG traces (5 channels) for poster illustration.
%   Fig 5  Parietal P300 ERP: grand-average by correct/incorrect × D/P.
%   Fig 6  RewP modulated by PREVIOUS uncertainty:
%            6a  Peak RewP at LE vs RN by block transition type (D→D/D→P/P→D/P→P)
%            6b  Peak RewP at LE vs RN by n_prev_P bin (0 / 1 / 2+)
%            Both panels: per-subject scatter + group bar with error bars.
%
% WINDOW DISCUSSION (see inline comments)
% ----------------------------------------
%   FRN/RewP windows are set as parameters at the top of the script.
%   After inspecting the grand-average waveforms the DEFAULT is expanded to
%   250–400 ms for both FRN (negativity) and RewP (positivity), matching the
%   broader literature consensus (Sambrook & Goslin 2015 meta-analysis: 200–400 ms).
%   The ORIGINAL 250–350 ms window is kept as a narrow comparison option.
%   See WINDOW_MODE parameter below.
%
% INPUTS
%   group_feature_table_combined.mat   per-trial EEG feature table
%   frn_rewp_by_stage_combined.mat     per-stage FRN/RewP table (+ t_ax)
%   all_trial_data.mat / behav_table.mat  behavioural data for Fig 3
%
% OUTPUT
%   Results/EEG analysis/Figures/Poster_P9/
% =============================================================================

clear; close all;

%% ── WINDOW MODE (key parameter) ─────────────────────────────────────────────
% 'broad'  → FRN/RewP window 200–400 ms  (recommended; covers observed peaks)
% 'narrow' → FRN/RewP window 250–350 ms  (original pipeline window)
WINDOW_MODE = 'broad';

%% ── PATHS ───────────────────────────────────────────────────────────────────
remote = 0;
if remote
    base_path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch';
else
    base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch';
end

kh_results = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis');
feat_dir   = fullfile(kh_results, 'Outcome_feature_tables_v4_merged');
data_dir   = fullfile(base_path, 'Salient mod switch KH', 'Data');
outdir     = fullfile(kh_results, 'Figures', 'Poster_P9');
if ~exist(outdir, 'dir'), mkdir(outdir); end

%% ── GLOBAL STYLE ────────────────────────────────────────────────────────────
set(groot,'defaultAxesTickDir','out');
set(groot,'defaultAxesBox','off');
set(groot,'defaultAxesFontSize',11);
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultLineLineWidth',1.8);
set(groot,'defaultTextFontName','Arial');

CLR_D       = [0.15 0.45 0.70];
CLR_P       = [0.80 0.30 0.10];
CLR_CORRECT = [0.10 0.60 0.10];
CLR_INCORR  = [0.80 0.10 0.10];
CLR_FC      = [0.55 0.25 0.65];
CLR_FI      = [0.90 0.60 0.10];

STAGE_NAMES  = {'LN','LE','RN','RE'};
STAGE_COLORS = [0.12 0.62 0.47;
                0.85 0.65 0.00;
                0.80 0.27 0.13;
                0.40 0.25 0.65];

% Transition colours (D→D, D→P, P→D, P→P)
TRANS_TYPES  = {'D->D','D->P','P->D','P->P'};
TRANS_COLORS = [0.12 0.47 0.71;
                0.85 0.33 0.10;
                0.47 0.67 0.19;
                0.80 0.20 0.60];

% n_prev_P colours (0, 1, 2+)
NPP_COLORS = [0.26 0.58 0.78;
              0.97 0.58 0.02;
              0.70 0.17 0.12];

ERP_XLIM = [-200 800];

% Window parameters
switch WINDOW_MODE
    case 'broad'
        FRN_WIN  = [200 400];   % broader – captures observed negativity peak
        REWP_WIN = [200 400];   % broader – captures observed positivity
        fprintf('Window mode: BROAD  (200–400 ms)\n');
    otherwise
        FRN_WIN  = [250 350];
        REWP_WIN = [250 350];
        fprintf('Window mode: NARROW (250–350 ms)\n');
end
BL_WIN = [-200 0];

%% ── LOAD DATA ───────────────────────────────────────────────────────────────
fprintf('Loading data...\n');

% Per-trial feature table
gt = [];
for cand = {fullfile(feat_dir,'group_feature_table_combined.mat'), ...
            fullfile(kh_results,'Epoched_data','group_feature_table_combined.mat'), ...
            fullfile(kh_results,'Epoched_data_noisefiltering','group_feature_table_combined.mat')}
    if exist(cand{1},'file')
        S = load(cand{1});
        if isfield(S,'group_table_combined'), gt = S.group_table_combined;
        elseif isfield(S,'group_table'),      gt = S.group_table;
        elseif isfield(S,'all_trials_table'), gt = S.all_trials_table; end
        if ~isempty(gt), fprintf('  Per-trial table: %s\n',cand{1}); break; end
    end
end

% FRN/RewP per-stage table
frn_tbl = [];
for cand = {fullfile(feat_dir,'frn_rewp_by_stage_combined.mat'), ...
            fullfile(feat_dir,'group_feature_table_KH_final.mat')}
    if exist(cand{1},'file')
        S = load(cand{1});
        if isfield(S,'frn_rewp_stage_table'), frn_tbl = S.frn_rewp_stage_table;
        elseif isfield(S,'frn_tbl'),          frn_tbl = S.frn_tbl; end
        if ~isempty(frn_tbl), fprintf('  FRN/RewP table: %s\n',cand{1}); break; end
    end
end

% Behavioural table (Fig 3)
group_T = [];
for cand = {fullfile(data_dir,'behav_table_June2026_RL.mat'), ...
            fullfile(data_dir,'behav_table_June2026.mat'), ...
            fullfile(data_dir,'behav_table.mat')}
    if exist(cand{1},'file')
        S = load(cand{1});
        if isfield(S,'group_T'), group_T = S.group_T; end
        if ~isempty(group_T), fprintf('  Behav table: %s\n',cand{1}); break; end
    end
end

% all_trial_data (Fig 3)
all_trial_data = [];
for cand = {fullfile(data_dir,'all_trial_data_June2026.mat'), ...
            fullfile(data_dir,'all_trial_data.mat')}
    if exist(cand{1},'file')
        S = load(cand{1});
        if isfield(S,'all_trial_data'), all_trial_data = S.all_trial_data; end
        if ~isempty(all_trial_data), fprintf('  all_trial_data: %s\n',cand{1}); break; end
    end
end

%% ── TIME AXIS ───────────────────────────────────────────────────────────────
t_ax = [];
for cand = {fullfile(feat_dir,'group_feature_table_KH_final.mat'), ...
            fullfile(feat_dir,'group_feature_table_combined.mat'), ...
            fullfile(kh_results,'Epoched_data','group_feature_table_combined.mat')}
    if exist(cand{1},'file')
        S = load(cand{1},'t_ax');
        if isfield(S,'t_ax') && ~isempty(S.t_ax), t_ax = S.t_ax; break; end
    end
end
if isempty(t_ax)
    t_ax = -200:2:800;
    fprintf('  WARNING: t_ax not found – using default -200:2:800 ms\n');
end

in_erp   = t_ax >= ERP_XLIM(1) & t_ax <= ERP_XLIM(2);
t_plot   = t_ax(in_erp);
n_t      = sum(in_erp);
frn_mask  = t_ax >= FRN_WIN(1)  & t_ax <= FRN_WIN(2);
rewp_mask = t_ax >= REWP_WIN(1) & t_ax <= REWP_WIN(2);
p300_mask = t_ax >= 300 & t_ax <= 600;   % P300 always fixed
bl_mask   = t_ax >= BL_WIN(1)   & t_ax <= BL_WIN(2);

% Normalise per-trial table column names for later use
if ~isempty(gt)
    gt_bt = string(gt.block_type);
    gt_bt(gt_bt=="V") = "P";
    gt.block_type_s = gt_bt;
    if isnumeric(gt.correct)||islogical(gt.correct)
        gt.correct_num = double(gt.correct);
    else
        cs = lower(string(gt.correct));
        gt.correct_num = double(cs=="1"|cs=="correct"|cs=="true");
    end
    if ~ismember('false_fb',gt.Properties.VariableNames)
        gt.false_fb = false(height(gt),1);
    else
        gt.false_fb = logical(gt.false_fb);
    end
    % Prefer waveform columns, but keep the selected names as scalar variables.
    % Do NOT store these in gt as table variables: assigning a char scalar to a
    % table with many rows causes a row-size mismatch.
    wave_col_use = '';
    for wc = {'prefrontal_waveform','FCzCz_waveform','prefrontal_signal'}
        if ismember(wc{1},gt.Properties.VariableNames)
            wave_col_use = wc{1};
            break;
        end
    end

    p300_col_use = '';
    for wc = {'P300_waveform','P300_signal'}
        if ismember(wc{1},gt.Properties.VariableNames)
            p300_col_use = wc{1};
            break;
        end
    end
else
    wave_col_use = '';
    p300_col_use = '';
end

%% ═══════════════════════════════════════════════════════════════════════════
%% FIGURE 1 — RewP Difference Waves per Stage × Block Type
%% ═══════════════════════════════════════════════════════════════════════════
fprintf('\nFig 1: RewP difference waves...\n');

fig1 = figure('Position',[50 50 1400 650],'Color','w');
sgtitle('Reward Positivity (RewP) — correct minus incorrect ERP by stage and block type',...
    'FontSize',13,'FontWeight','bold');

if ~isempty(frn_tbl) && ismember('diff_wave',frn_tbl.Properties.VariableNames)
    for bi = 1:2
        bt  = {'D','P'}; clr = ternary_p9(bi==1,CLR_D,CLR_P);
        for si = 1:4
            ax = subplot(2,4,(bi-1)*4+si); hold(ax,'on');
            sel = string(frn_tbl.block_type)==bt{bi} & string(frn_tbl.stage)==STAGE_NAMES{si};
            dw  = frn_tbl.diff_wave(sel);
            dw  = dw(~cellfun(@isempty,dw));
            if ~isempty(dw)
                M = stack_waves_p9(dw, t_ax, -1);   % invert: FRN table stores inc-cor; RewP = cor-inc
                mn = mean(M(:,in_erp),1,'omitnan');
                se = std( M(:,in_erp),0,1,'omitnan')/sqrt(size(M,1));
                yl = compute_yl(mn,se);
                shade_win(ax,FRN_WIN, yl,[0.85 0.90 1.00]);
                shade_win(ax,REWP_WIN,yl,[1.00 0.92 0.85]);
                plot_ribbon_p9(ax,t_plot,mn,se,clr,'-', ...
                    sprintf('%s-%s (n=%d)',bt{bi},STAGE_NAMES{si},size(M,1)));
                % Mark mean RewP amplitude in window
                rewp_in = rewp_mask(in_erp);
                if any(rewp_in)
                    rv = mean(mn(rewp_in),'omitnan');
                    text(ax,mean(REWP_WIN),rv,sprintf(' %.2f µV',rv),...
                        'Color',clr,'FontSize',8,'VerticalAlignment','middle');
                end
            else
                text(ax,0.5,0.5,'No data','Units','normalized',...
                    'HorizontalAlignment','center','Color',[0.6 0.6 0.6]);
            end
            xline(ax,0,'k:','LineWidth',1,'HandleVisibility','off');
            yline(ax,0,'k--','LineWidth',0.8,'HandleVisibility','off');
            xlim(ax,ERP_XLIM);
            xlabel(ax,'Time (ms)','FontSize',9);
            if si==1, ylabel(ax,'RewP (µV) [correct − incorrect]','FontSize',9); end
            title(ax,sprintf('%s blocks — %s',bt{bi},STAGE_NAMES{si}),'FontSize',10);
            legend(ax,'Box','off','FontSize',7,'Location','northeast');
        end
    end
    annotation(fig1,'textbox',[0.01 0.01 0.98 0.04],'String',...
        ['RewP = correct minus incorrect ERP. Blue shading = FRN window. '...
         'Orange shading = RewP window (' num2str(REWP_WIN(1)) '–' num2str(REWP_WIN(2)) ...
         ' ms; WINDOW_MODE=' WINDOW_MODE '). Grand average ± 1 SEM across subjects.'],...
        'FontSize',8,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);
end
save_fig_p9(fig1,outdir,'P9_Fig1_RewP_by_stage_blocktype');

%% ═══════════════════════════════════════════════════════════════════════════
%% FIGURE 2 — Grand-Average ERPs by Condition + RewP wave below
%% ═══════════════════════════════════════════════════════════════════════════
fprintf('Fig 2: Grand-average ERP by condition...\n');

fig2 = figure('Position',[50 50 1400 900],'Color','w');
sgtitle('Grand-average ERPs: correct, incorrect, false feedback trials',...
    'FontSize',13,'FontWeight','bold');

for bi = 1:2
    bt = {'D','P'}; is_P = bi==2;
    wc = '';
    if ~isempty(gt), wc = wave_col_use; end
    has_w = ~isempty(gt) && ~isempty(wc) && ismember(wc,gt.Properties.VariableNames);

    if has_w
        mk = @(c,fb) gt.block_type_s==bt{bi} & gt.correct_num==c & gt.false_fb==fb;
        [mn_tc,se_tc,n_tc] = extract_ga(gt,wc,mk(1,false),t_ax,n_t,in_erp);
        [mn_ti,se_ti,n_ti] = extract_ga(gt,wc,mk(0,false),t_ax,n_t,in_erp);
        if is_P
            [mn_fc,se_fc,n_fc] = extract_ga(gt,wc,mk(1,true),t_ax,n_t,in_erp);
            [mn_fi,se_fi,n_fi] = extract_ga(gt,wc,mk(0,true),t_ax,n_t,in_erp);
        end
    else
        mn_tc=[]; mn_ti=[]; mn_fc=[]; mn_fi=[]; n_tc=0; n_ti=0; n_fc=0; n_fi=0;
    end

    % ERP panel
    ax_erp = subplot(2,2,bi); hold(ax_erp,'on');
    title(ax_erp,sprintf('%s blocks — grand-average ERP',...
        ternary_p9(bi==1,'Deterministic','Probabilistic')),'FontSize',11);
    shade_win(ax_erp,FRN_WIN, [-3 3],[0.85 0.90 1.00]);
    shade_win(ax_erp,REWP_WIN,[-3 3],[1.00 0.92 0.85]);
    if ~isempty(mn_tc)
        plot_ribbon_p9(ax_erp,t_plot,mn_tc,se_tc,CLR_CORRECT,'-',sprintf('True correct (n=%d)',n_tc));
        plot_ribbon_p9(ax_erp,t_plot,mn_ti,se_ti,CLR_INCORR,'--',sprintf('True incorrect (n=%d)',n_ti));
        if is_P && ~isempty(mn_fc)
            plot_ribbon_p9(ax_erp,t_plot,mn_fc,se_fc,CLR_FC,'-.',sprintf('False correct* (n=%d)',n_fc));
            plot_ribbon_p9(ax_erp,t_plot,mn_fi,se_fi,CLR_FI,':',sprintf('False incorrect† (n=%d)',n_fi));
        end
    elseif ~isempty(frn_tbl) && ismember('diff_wave',frn_tbl.Properties.VariableNames)
        sel = string(frn_tbl.block_type)==bt{bi};
        dw  = frn_tbl.diff_wave(sel); dw=dw(~cellfun(@isempty,dw));
        if ~isempty(dw)
            M  = stack_waves_p9(dw,t_ax,1);
            mn = mean(M(:,in_erp),1,'omitnan');
            se = std( M(:,in_erp),0,1,'omitnan')/sqrt(size(M,1));
            plot_ribbon_p9(ax_erp,t_plot,mn,se,ternary_p9(bi==1,CLR_D,CLR_P),'-',...
                sprintf('%s difference wave (n=%d)',bt{bi},size(M,1)));
        end
    end
    xline(ax_erp,0,'k:','HandleVisibility','off');
    yline(ax_erp,0,'k--','HandleVisibility','off');
    set(ax_erp,'YDir','reverse'); xlim(ax_erp,ERP_XLIM);
    xlabel(ax_erp,'Time (ms)'); ylabel(ax_erp,'Amplitude (µV) [negative up]');
    legend(ax_erp,'Box','off','FontSize',8,'Location','northwest');
    yl_e = ylim(ax_erp);
    text(ax_erp,mean(FRN_WIN),yl_e(1)*0.82,'FRN','HorizontalAlignment','center',...
        'FontSize',8,'Color',[0.3 0.3 0.8]);
    text(ax_erp,mean(REWP_WIN),yl_e(2)*0.82,'RewP','HorizontalAlignment','center',...
        'FontSize',8,'Color',[0.8 0.4 0.1]);

    % RewP panel
    ax_rw = subplot(2,2,2+bi); hold(ax_rw,'on');
    title(ax_rw,sprintf('%s blocks — RewP wave (correct − incorrect)',...
        ternary_p9(bi==1,'Deterministic','Probabilistic')),'FontSize',11);
    shade_win(ax_rw,FRN_WIN, [-2 2],[0.85 0.90 1.00]);
    shade_win(ax_rw,REWP_WIN,[-2 2],[1.00 0.92 0.85]);
    if ~isempty(mn_tc) && ~isempty(mn_ti)
        rw  = mn_tc - mn_ti;
        rw_se = sqrt((se_tc.^2+se_ti.^2)/2);
        plot_ribbon_p9(ax_rw,t_plot,rw,rw_se,ternary_p9(bi==1,CLR_D,CLR_P),'-',...
            sprintf('True FB: correct−incorrect (n≈%d)',round((n_tc+n_ti)/2)));
        if is_P && ~isempty(mn_fc) && ~isempty(mn_fi)
            rw_f    = mn_fi - mn_fc;
            rw_f_se = sqrt((se_fi.^2+se_fc.^2)/2);
            plot_ribbon_p9(ax_rw,t_plot,rw_f,rw_f_se,CLR_FC,'--',...
                sprintf('False FB: told-correct−told-incorrect (n≈%d)',round((n_fi+n_fc)/2)));
        end
        % Peak marker
        ri = rewp_mask(in_erp);
        if any(ri)
            [pv,pi2] = max(rw(ri));
            rt = t_plot(ri); rt = rt(pi2);
            if ~isnan(pv)
                plot(ax_rw,rt,pv,'o','Color',ternary_p9(bi==1,CLR_D,CLR_P),...
                    'MarkerFaceColor',ternary_p9(bi==1,CLR_D,CLR_P),'MarkerSize',7,...
                    'HandleVisibility','off');
                text(ax_rw,rt,pv+0.05,sprintf('%.2f µV',pv),'FontSize',8,...
                    'HorizontalAlignment','center','Color',ternary_p9(bi==1,CLR_D,CLR_P));
            end
        end
    elseif ~isempty(frn_tbl) && ismember('diff_wave',frn_tbl.Properties.VariableNames)
        sel = string(frn_tbl.block_type)==bt{bi};
        dw  = frn_tbl.diff_wave(sel); dw=dw(~cellfun(@isempty,dw));
        if ~isempty(dw)
            M  = stack_waves_p9(dw,t_ax,-1);
            mn = mean(M(:,in_erp),1,'omitnan');
            se = std( M(:,in_erp),0,1,'omitnan')/sqrt(size(M,1));
            plot_ribbon_p9(ax_rw,t_plot,mn,se,ternary_p9(bi==1,CLR_D,CLR_P),'-',...
                sprintf('%s RewP (n=%d)',bt{bi},size(M,1)));
        end
    end
    xline(ax_rw,0,'k:','HandleVisibility','off');
    yline(ax_rw,0,'k--','HandleVisibility','off');
    xlim(ax_rw,ERP_XLIM);
    xlabel(ax_rw,'Time (ms)'); ylabel(ax_rw,'RewP (µV) [correct − incorrect]');
    legend(ax_rw,'Box','off','FontSize',8,'Location','northwest');
    yl_r = ylim(ax_rw);
    text(ax_rw,mean(FRN_WIN),yl_r(1)+0.1*range(yl_r),'FRN window',...
        'HorizontalAlignment','center','FontSize',7,'Color',[0.3 0.3 0.8]);
    text(ax_rw,mean(REWP_WIN),yl_r(2)-0.1*range(yl_r),'RewP window',...
        'HorizontalAlignment','center','FontSize',7,'Color',[0.8 0.4 0.1]);
end

annotation(fig2,'textbox',[0.01 0.01 0.98 0.05],'String',...
    ['* False correct = participant correct, shown incorrect feedback (told wrong). '...
     '† False incorrect = participant incorrect, shown correct feedback (told right). '...
     'Window: ' num2str(REWP_WIN(1)) '–' num2str(REWP_WIN(2)) ' ms (' WINDOW_MODE '). '...
     'Top row negative-up (EEG convention). Peak amplitude labelled.'],...
    'FontSize',8,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);
save_fig_p9(fig2,outdir,'P9_Fig2_ERP_by_condition');

%% ═══════════════════════════════════════════════════════════════════════════
%% FIGURE 3 — First-10-Trials Accuracy: D vs P blocks
%% ═══════════════════════════════════════════════════════════════════════════
fprintf('Fig 3: First-10-trials accuracy...\n');

N_FIRST = 10;
first10_D = []; first10_P = []; sidx_D = {}; sidx_P = {};

if ~isempty(all_trial_data)
    sids = fieldnames(all_trial_data);
    for si = 1:numel(sids)
        sn  = sids{si};
        if ~isfield(all_trial_data.(sn),'trial_data'), continue; end
        td  = all_trial_data.(sn).trial_data;
        if ~isfield(td,'correct'), continue; end
        [nB,nT] = size(td.correct);
        bs = infer_block_structure(td,nB);
        for b = 1:nB
            if nT < N_FIRST, continue; end
            row = double(td.correct(b,1:N_FIRST));
            if all(isnan(row)), continue; end
            ct = char(bs(min(b,numel(bs))));
            if ct=='D',      first10_D(end+1,:)=row; sidx_D{end+1}=sn;
            elseif ct=='P',  first10_P(end+1,:)=row; sidx_P{end+1}=sn; end
        end
    end
elseif ~isempty(group_T)
    gt3 = group_T; gt3.correct=double(gt3.correct);
    gt3.trial=double(gt3.trial); gt3.block=double(gt3.block);
    gt3.bt_clean=string(gt3.block_type); gt3.bt_clean(gt3.bt_clean=="V")="P";
    for sn = unique(string(gt3.subjID))'
        Ts=gt3(string(gt3.subjID)==sn,:);
        for b = unique(Ts.block)'
            Tb=Ts(Ts.block==b & Ts.trial<=N_FIRST,:);
            if height(Tb)<N_FIRST, continue; end
            [~,o]=sort(Tb.trial); Tb=Tb(o,:);
            row=Tb.correct(1:N_FIRST)'; ct=char(Tb.bt_clean(1));
            if ct=='D',     first10_D(end+1,:)=row; sidx_D{end+1}=char(sn);
            elseif ct=='P', first10_P(end+1,:)=row; sidx_P{end+1}=char(sn); end
        end
    end
end

fig3 = figure('Position',[50 50 1200 560],'Color','w');
sgtitle('First 10 trials accuracy per block: D vs P blocks',...
    'FontSize',13,'FontWeight','bold');
trial_ax3 = 1:N_FIRST;

for bi = 1:2
    ax3 = subplot(1,2,bi); hold(ax3,'on');
    if bi==1, mat=first10_D; clr=CLR_D; lbl='Deterministic'; slist=sidx_D;
    else,      mat=first10_P; clr=CLR_P; lbl='Probabilistic'; slist=sidx_P; end

    if isempty(mat)
        text(ax3,0.5,0.5,'No data','Units','normalized','HorizontalAlignment','center');
        title(ax3,lbl); continue
    end
    % Faint individual-block traces
    for ri=1:size(mat,1)
        row=mat(ri,:); ok=~isnan(row);
        if sum(ok)<2, continue; end
        hl=plot(ax3,trial_ax3(ok),row(ok),'-','Color',clr,'LineWidth',0.6,'HandleVisibility','off');
        try,hl.Color(4)=0.18;catch,end
    end
    % Per-subject means for proper SEM
    usubs=unique(string(slist));
    sm=NaN(numel(usubs),N_FIRST);
    for si2=1:numel(usubs)
        rm=strcmp(string(slist),usubs(si2));
        if any(rm), sm(si2,:)=mean(mat(rm,:),1,'omitnan'); end
    end
    gm=mean(sm,1,'omitnan'); gs=std(sm,0,1,'omitnan')./sqrt(sum(~isnan(sm),1));
    fill(ax3,[trial_ax3,fliplr(trial_ax3)],[gm+gs,fliplr(gm-gs)],...
        clr,'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off');
    plot(ax3,trial_ax3,gm,'-','Color',clr,'LineWidth',2.8,...
        'DisplayName',sprintf('Group mean ± SEM (n=%d subj, %d blocks)',numel(usubs),size(mat,1)));
    yline(ax3,0.5,'k--','LineWidth',1,'HandleVisibility','off');
    n_per=sum(~isnan(mat),1);
    for ti=1:N_FIRST
        text(ax3,ti,0.02,sprintf('%d',n_per(ti)),...
            'HorizontalAlignment','center','FontSize',7,'Color',[0.5 0.5 0.5]);
    end
    xlim(ax3,[0.5 N_FIRST+0.5]); ylim(ax3,[0 1.05]);
    set(ax3,'XTick',1:N_FIRST,'TickDir','out');
    xlabel(ax3,'Trial within block'); ylabel(ax3,'P(correct)');
    title(ax3,sprintf('%s blocks — first %d trials',lbl,N_FIRST),'FontSize',11);
    legend(ax3,'Box','off','FontSize',9,'Location','southeast');
end
annotation(fig3,'textbox',[0.01 0.01 0.98 0.04],'String',...
    ['Faint lines = individual blocks pooled across subjects. '...
     'Bold = per-subject-averaged group mean ± 1 SEM. '...
     'n (grey) = number of blocks contributing to each trial position. '...
     'Dashed line = chance (0.5).'],...
    'FontSize',8,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);
save_fig_p9(fig3,outdir,'P9_Fig3_First10trials_DvsP');

%% ═══════════════════════════════════════════════════════════════════════════
%% FIGURE 4 — Mock Raw EEG Traces (poster illustration)
%% ═══════════════════════════════════════════════════════════════════════════
fprintf('Fig 4: Mock raw EEG...\n');

rng(42);
fs=500; dur=4.0; t_eeg=(0:1/fs:dur-1/fs); n_samp=numel(t_eeg); ev_t=2.0; ev_s=round(ev_t*fs);
ch_names={'FCz','Cz','Pz','C3','C4'};
ch_clrs=[0.15 0.35 0.65;0.30 0.55 0.75;0.65 0.30 0.15;0.20 0.60 0.30;0.60 0.25 0.50];
alpha_f=10+randn(1,5)*0.5;
eeg=zeros(5,n_samp);
for ch=1:5
    pink=cumsum(randn(1,n_samp+200)); pink=pink(201:end)/std(pink(201:end))*8;
    ae=ones(1,n_samp)*0.6; ae(ev_s:end)=0.4+0.8*exp(-(0:n_samp-ev_s-1)/(0.5*fs));
    alp=ae.*sin(2*pi*alpha_f(ch)*t_eeg)*12;
    th=(6*(ch<=2)+3*(ch>2)).*sin(2*pi*6*t_eeg+rand*2*pi);
    erp=zeros(1,n_samp);
    if ch<=3
        tp=t_eeg-ev_t;
        erp=erp-6*exp(-(tp-0.25).^2/(2*0.03^2)).*(tp>=0);
        pa=(4*(ch==3)+2*(ch<=2)); erp=erp+pa*exp(-(tp-0.40).^2/(2*0.05^2)).*(tp>=0);
    end
    bt2=(3*(ch>=4)+1*(ch<4)).*sin(2*pi*20*t_eeg+rand*2*pi);
    blink=zeros(1,n_samp); if ch<=2, blink=80*exp(-((t_eeg-0.8).^2)/(2*0.06^2)); end
    musc=zeros(1,n_samp);
    if ch>=4, mm=t_eeg>=1.5&t_eeg<=1.65; musc(mm)=randn(1,sum(mm))*15; end
    eeg(ch,:)=pink+alp+th+erp+bt2+blink+musc+randn(1,n_samp)*2;
end

fig4=figure('Position',[50 50 1300 700],'Color','w');
off=120; offs=(0:4)*off;
ax4=axes(fig4,'Position',[0.08 0.12 0.86 0.78]); hold(ax4,'on');
yl4=[-off/2,5*off+off/2];
patch(ax4,[ev_t-0.2 ev_t+0.8 ev_t+0.8 ev_t-0.2],[yl4(1) yl4(1) yl4(2) yl4(2)],...
    [0.88 0.94 1.0],'EdgeColor','none','FaceAlpha',0.6,'HandleVisibility','off');
for ch=1:5
    plot(ax4,t_eeg,eeg(ch,:)+offs(ch),'Color',ch_clrs(ch,:),'LineWidth',1.0,'DisplayName',ch_names{ch});
end
xline(ax4,ev_t,'k-','LineWidth',2.5,'HandleVisibility','off');
text(ax4,ev_t+0.04,yl4(2)*0.90,'Outcome','FontSize',10,'FontWeight','bold');
text(ax4,0.82,offs(2)+55,'← Eye blink','FontSize',8,'Color',[0.6 0.2 0.2],'FontStyle','italic');
text(ax4,1.52,offs(4)+55,'← Muscle','FontSize',8,'Color',[0.5 0.5 0.5],'FontStyle','italic');
text(ax4,ev_t+0.22,offs(1)-40,'N2↓','FontSize',9,'Color',ch_clrs(1,:),'FontWeight','bold');
text(ax4,ev_t+0.40,offs(3)+30,'P300↑','FontSize',9,'Color',ch_clrs(3,:),'FontWeight','bold');
sb_t=[3.50 4.00]; sb_c=offs(3);
plot(ax4,sb_t,[sb_c sb_c],'k-','LineWidth',2.5,'HandleVisibility','off');
plot(ax4,[sb_t(2) sb_t(2)],sb_c+[-25 25],'k-','LineWidth',2.5,'HandleVisibility','off');
text(ax4,mean(sb_t),sb_c-30,'500 ms','HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
text(ax4,sb_t(2)+0.04,sb_c,'50 µV','HorizontalAlignment','left','FontSize',9,'FontWeight','bold');
set(ax4,'YTick',offs,'YTickLabel',ch_names,'FontSize',11,'TickDir','out','Box','off');
xlim(ax4,[t_eeg(1) t_eeg(end)]); ylim(ax4,yl4);
xlabel(ax4,'Time (s)','FontSize',12);
title(ax4,'Illustrative raw EEG traces — 5 channels, tactile category-switch outcome epoch',...
    'FontSize',12,'FontWeight','bold');
legend(ax4,'Box','off','FontSize',10,'Location','northeast');
annotation(fig4,'textbox',[0.01 0.01 0.98 0.04],'String',...
    'Synthetic, physiologically plausible EEG. Blue shading = outcome epoch (−200 to +800 ms). Alpha (∼10 Hz), theta (∼6 Hz), N2 and P300 morphology visible. Eye blink (FCz/Cz) and muscle (C3/C4) shown for context.',...
    'FontSize',8,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);
save_fig_p9(fig4,outdir,'P9_Fig4_Mock_raw_EEG');

%% ═══════════════════════════════════════════════════════════════════════════
%% FIGURE 5 — Parietal P300 ERP (correct vs incorrect, D and P panels)
%% ═══════════════════════════════════════════════════════════════════════════
fprintf('Fig 5: Parietal P300 ERP...\n');

fig5 = figure('Position',[50 50 1400 900],'Color','w');
sgtitle('Parietal P300 ERP: correct vs incorrect by block type',...
    'FontSize',13,'FontWeight','bold');

for bi = 1:2
    bt = {'D','P'}; is_P2 = bi==2;
    wcp = '';
    if ~isempty(gt), wcp = p300_col_use; end
    has_wp = ~isempty(gt) && ~isempty(wcp) && ismember(wcp,gt.Properties.VariableNames);

    if has_wp
        mk = @(c,fb) gt.block_type_s==bt{bi} & gt.correct_num==c & gt.false_fb==fb;
        [p_tc,ps_tc,np_tc] = extract_ga(gt,wcp,mk(1,false),t_ax,n_t,in_erp);
        [p_ti,ps_ti,np_ti] = extract_ga(gt,wcp,mk(0,false),t_ax,n_t,in_erp);
        if is_P2
            [p_fc,ps_fc,np_fc] = extract_ga(gt,wcp,mk(1,true),t_ax,n_t,in_erp);
            [p_fi,ps_fi,np_fi] = extract_ga(gt,wcp,mk(0,true),t_ax,n_t,in_erp);
        end
    else
        % Scalar fallback using P300_amp from per-trial table
        p_tc=[]; p_ti=[]; p_fc=[]; p_fi=[]; np_tc=0; np_ti=0; np_fc=0; np_fi=0;
    end

    % Grand-average waveform panel (top)
    ax5a = subplot(2,2,bi); hold(ax5a,'on');
    title(ax5a,sprintf('%s blocks — Parietal P300 (Pz/P1/P2)',...
        ternary_p9(bi==1,'Deterministic','Probabilistic')),'FontSize',11);
    shade_win(ax5a,[300 600],[-2 4],[1.00 0.95 0.85]);  % P300 window always fixed
    if ~isempty(p_tc)
        plot_ribbon_p9(ax5a,t_plot,p_tc,ps_tc,CLR_CORRECT,'-',sprintf('True correct (n=%d)',np_tc));
        plot_ribbon_p9(ax5a,t_plot,p_ti,ps_ti,CLR_INCORR,'--',sprintf('True incorrect (n=%d)',np_ti));
        if is_P2 && ~isempty(p_fc)
            plot_ribbon_p9(ax5a,t_plot,p_fc,ps_fc,CLR_FC,'-.',sprintf('False correct* (n=%d)',np_fc));
            plot_ribbon_p9(ax5a,t_plot,p_fi,ps_fi,CLR_FI,':',sprintf('False incorrect† (n=%d)',np_fi));
        end
        % Mark P300 peak in correct trials
        pk_in = p300_mask(in_erp);
        if any(pk_in) && ~isempty(p_tc)
            [pk_v,pk_i] = max(p_tc(pk_in));
            pt_v  = t_plot(pk_in); pt_v = pt_v(pk_i);
            if ~isnan(pk_v)
                plot(ax5a,pt_v,pk_v,'o','Color',CLR_CORRECT,...
                    'MarkerFaceColor',CLR_CORRECT,'MarkerSize',7,'HandleVisibility','off');
                text(ax5a,pt_v,pk_v+0.1,sprintf('%.2f µV',pk_v),...
                    'FontSize',8,'HorizontalAlignment','center','Color',CLR_CORRECT);
            end
        end
    else
        % Fallback: P300_amp scalar per stage
        if ~isempty(gt) && ismember('P300_amp',gt.Properties.VariableNames)
            for ci = 0:1
                msk_p = gt.block_type_s==bt{bi} & gt.correct_num==ci & ~gt.false_fb;
                stage_mu = NaN(1,4); stage_se = NaN(1,4);
                for si4=1:4
                    v=gt.P300_amp(msk_p & string(gt.stage)==STAGE_NAMES{si4});
                    v=v(~isnan(v));
                    if ~isempty(v), stage_mu(si4)=mean(v); stage_se(si4)=std(v)/sqrt(numel(v)); end
                end
                clr5=ternary_p9(ci==1,CLR_CORRECT,CLR_INCORR);
                errorbar(ax5a,1:4,stage_mu,stage_se,'o-','Color',clr5,...
                    'MarkerFaceColor',clr5,'MarkerSize',6,...
                    'DisplayName',ternary_p9(ci==1,'Correct','Incorrect'));
            end
            set(ax5a,'XTick',1:4,'XTickLabel',STAGE_NAMES);
            xlabel(ax5a,'Stage'); ylabel(ax5a,'P300 amplitude (µV)');
            legend(ax5a,'Box','off','FontSize',9,'Location','best');
        else
            text(ax5a,0.5,0.5,'P300 waveform data not available',...
                'Units','normalized','HorizontalAlignment','center','Color',[0.5 0.5 0.5]);
        end
    end
    if ~isempty(p_tc)
        xline(ax5a,0,'k:','HandleVisibility','off');
        yline(ax5a,0,'k--','HandleVisibility','off');
        xlim(ax5a,ERP_XLIM);
        xlabel(ax5a,'Time (ms)'); ylabel(ax5a,'Amplitude (µV)');
        legend(ax5a,'Box','off','FontSize',8,'Location','northwest');
        yl5=ylim(ax5a);
        text(ax5a,450,yl5(2)-0.12*range(yl5),'P300 window',...
            'HorizontalAlignment','center','FontSize',8,'Color',[0.8 0.5 0.1]);
    end

    % P300 difference wave panel (bottom): correct - incorrect
    ax5b = subplot(2,2,2+bi); hold(ax5b,'on');
    title(ax5b,sprintf('%s blocks — P300 difference (correct − incorrect)',...
        ternary_p9(bi==1,'Deterministic','Probabilistic')),'FontSize',11);
    shade_win(ax5b,[300 600],[-2 3],[1.00 0.95 0.85]);
    if ~isempty(p_tc) && ~isempty(p_ti)
        pd   = p_tc - p_ti;
        pd_se = sqrt((ps_tc.^2+ps_ti.^2)/2);
        plot_ribbon_p9(ax5b,t_plot,pd,pd_se,ternary_p9(bi==1,CLR_D,CLR_P),'-',...
            sprintf('Correct−incorrect (n≈%d)',round((np_tc+np_ti)/2)));
        % Peak annotation
        pk2_in = p300_mask(in_erp);
        if any(pk2_in)
            [pk2_v,pk2_i]=max(pd(pk2_in));
            pt2=t_plot(pk2_in); pt2=pt2(pk2_i);
            if ~isnan(pk2_v)
                plot(ax5b,pt2,pk2_v,'o','Color',ternary_p9(bi==1,CLR_D,CLR_P),...
                    'MarkerFaceColor',ternary_p9(bi==1,CLR_D,CLR_P),'MarkerSize',7,'HandleVisibility','off');
                text(ax5b,pt2,pk2_v+0.08,sprintf('%.2f µV',pk2_v),...
                    'FontSize',8,'HorizontalAlignment','center',...
                    'Color',ternary_p9(bi==1,CLR_D,CLR_P));
            end
        end
    end
    xline(ax5b,0,'k:','HandleVisibility','off');
    yline(ax5b,0,'k--','HandleVisibility','off');
    xlim(ax5b,ERP_XLIM);
    xlabel(ax5b,'Time (ms)'); ylabel(ax5b,'P300 difference (µV)');
    legend(ax5b,'Box','off','FontSize',8,'Location','northwest');
end

annotation(fig5,'textbox',[0.01 0.01 0.98 0.05],'String',...
    ['P300 reflects context updating (Polich 2007). Parietal ROI: Pz/P1/P2 (KH) or equivalent (RR). '...
     'P300 window: 300–600 ms (fixed). Bottom row: difference wave (correct minus incorrect). '...
     'P300 should be larger to unexpected, surprising outcomes independent of valence.'],...
    'FontSize',8,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);
save_fig_p9(fig5,outdir,'P9_Fig5_P300_by_condition');
fprintf('  Fig 5 saved.\n');

%% ═══════════════════════════════════════════════════════════════════════════
%% FIGURE 6 — RewP modulated by PREVIOUS uncertainty (LE vs RN)
%%
%% Uses the per-trial table to compute mean RewP amplitude in the RewP window
%% (correct minus incorrect difference) at LE and RN stages, split by:
%%   6a  Block transition type (D→D / D→P / P→D / P→P)
%%   6b  Cumulative prior P exposure (n_prev_P: 0 / 1 / 2+)
%%
%% If per-trial waveforms are available the RewP amplitude is extracted
%% from the waveform; otherwise the scalar RewP_mean_amp column is used.
%% ═══════════════════════════════════════════════════════════════════════════
fprintf('Fig 6: RewP modulated by previous uncertainty (LE vs RN)...\n');

fig6 = figure('Position',[50 50 1500 900],'Color','w');
sgtitle({'RewP modulated by prior uncertainty: LE (pre-reversal) vs RN (post-reversal)',...
    'Previous block type shapes how the brain processes outcomes at the reversal boundary'},...
    'FontSize',13,'FontWeight','bold');

% ── Derive transition and n_prev_P on per-trial table ─────────────────────
if ~isempty(gt)
    % Rebuild previous-block-type columns robustly using block column
    if ~ismember('n_prev_P',gt.Properties.VariableNames) || ...
       ~ismember('transition_cat',gt.Properties.VariableNames)
        fprintf('  Deriving transition and n_prev_P columns...\n');
        gt = derive_prev_uncertainty(gt);
    end

    has_trans = ismember('transition_cat',gt.Properties.VariableNames);
    has_npp   = ismember('n_prev_P',gt.Properties.VariableNames);

    % Determine best source of RewP amplitude
    use_wave = ~isempty(wave_col_use) && ismember(wave_col_use,gt.Properties.VariableNames);
    use_scalar_rewp = ismember('RewP_mean_amp',gt.Properties.VariableNames);
    use_scalar_frn  = ismember('FRN_mean_amp', gt.Properties.VariableNames);

    % Build per-subject × stage × transition RewP (correct - incorrect in window)
    % We extract subject-level RewP contrast at LE and RN for each grouping variable

    rev_stages = {'LE','RN'};   % the two reversal-boundary stages
    if ismember('subj_id',gt.Properties.VariableNames)
        subjs6 = unique(string(gt.subj_id));
    else
        subjs6 = unique(string(gt.subjID));
    end

    % ── Transition analysis (Fig 6a) ──────────────────────────────────────
    if has_trans
        % Build: subjects × stages × transitions matrix of RewP amplitude
        trans_levs = TRANS_TYPES;
        N_subj6    = numel(subjs6);
        N_stg6     = numel(rev_stages);
        N_tr6      = numel(trans_levs);

        rewp_by_trans = NaN(N_subj6, N_stg6, N_tr6);

        for si6 = 1:N_subj6
            sn6 = subjs6(si6);
            sm6 = get_subj_mask(gt,sn6);
            for sti = 1:N_stg6
                sg6 = rev_stages{sti};
                stg_m = sm6 & string(gt.stage)==sg6 & ~gt.false_fb;
                for tri = 1:N_tr6
                    tr_m = stg_m & string(gt.transition_cat)==trans_levs{tri};
                    rewp_by_trans(si6,sti,tri) = ...
                        compute_rewp_contrast(gt,tr_m,use_wave,use_scalar_rewp,...
                        use_scalar_frn,rewp_mask,in_erp,t_ax,wave_col_use);
                end
            end
        end

        % Plot: 2 × N_tr6 panels (LE vs RN per transition)
        for tri = 1:N_tr6
            ax6a = subplot(3,4,tri); hold(ax6a,'on');
            title(ax6a,sprintf('Transition: %s',strrep(trans_levs{tri},'->','→')),...
                'FontSize',10);

            for sti = 1:N_stg6
                vals = rewp_by_trans(:,sti,tri);
                vals = vals(~isnan(vals));
                if isempty(vals), continue; end
                xpos = sti;
                bar(ax6a,xpos,mean(vals,'omitnan'),0.5,...
                    'FaceColor',STAGE_COLORS(sti+1,:),'EdgeColor','none','FaceAlpha',0.8,...
                    'HandleVisibility','off');
                errorbar(ax6a,xpos,mean(vals,'omitnan'),std(vals,'omitnan')/sqrt(numel(vals)),...
                    'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');
                % Subject dots
                jx = xpos + 0.15*(rand(numel(vals),1)-0.5);
                scatter(ax6a,jx,vals,20,STAGE_COLORS(sti+1,:),'filled',...
                    'MarkerFaceAlpha',0.5,'HandleVisibility','off');
                % n
                text(ax6a,xpos,-0.02,sprintf('n=%d',numel(vals)),...
                    'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
            end

            % Connector line (LE to RN) per subject
            for si6 = 1:N_subj6
                v_le = rewp_by_trans(si6,1,tri);
                v_rn = rewp_by_trans(si6,2,tri);
                if ~isnan(v_le) && ~isnan(v_rn)
                    plot(ax6a,[1 2],[v_le v_rn],'-','Color',[0.5 0.5 0.5],...
                        'LineWidth',0.6,'HandleVisibility','off');
                end
            end

            % Paired t-test
            v_le_all = rewp_by_trans(:,1,tri); v_rn_all = rewp_by_trans(:,2,tri);
            ok_pair  = ~isnan(v_le_all) & ~isnan(v_rn_all);
            if sum(ok_pair) > 2
                [~,p_pair] = ttest(v_le_all(ok_pair), v_rn_all(ok_pair));
                add_sig_bracket_p9(ax6a, 1, 2, ...
                    max([v_le_all(ok_pair);v_rn_all(ok_pair)],[],'omitnan')*1.1, ...
                    p_pair);
            end

            set(ax6a,'XTick',[1 2],'XTickLabel',{'LE','RN'},'TickDir','out');
            ylabel(ax6a,'RewP amplitude (µV)','FontSize',9);
            xlabel(ax6a,'Stage','FontSize',9);
            yline(ax6a,0,'k--','LineWidth',0.8,'HandleVisibility','off');
            % Colour the transition label
            idx_tc = find(strcmp(trans_levs,trans_levs{tri}));
            ax6a.Title.Color = TRANS_COLORS(idx_tc,:);
        end
    end

    % ── n_prev_P analysis (Fig 6b) ────────────────────────────────────────
    if has_npp
        npp_bins   = {0, 1, [2 Inf]};
        npp_labels = {'0 prior P','1 prior P','2+ prior P'};
        N_npp      = numel(npp_bins);

        rewp_by_npp = NaN(numel(subjs6), N_stg6, N_npp);

        for si6 = 1:numel(subjs6)
            sn6 = subjs6(si6);
            sm6 = get_subj_mask(gt,sn6);
            for sti = 1:N_stg6
                sg6  = rev_stages{sti};
                stg_m = sm6 & string(gt.stage)==sg6 & ~gt.false_fb;
                for ni = 1:N_npp
                    bin  = npp_bins{ni};
                    if numel(bin)==1
                        npp_m = stg_m & gt.n_prev_P == bin(1);
                    else
                        npp_m = stg_m & gt.n_prev_P >= bin(1) & gt.n_prev_P < bin(2);
                    end
                    rewp_by_npp(si6,sti,ni) = ...
                        compute_rewp_contrast(gt,npp_m,use_wave,use_scalar_rewp,...
                        use_scalar_frn,rewp_mask,in_erp,t_ax,wave_col_use);
                end
            end
        end

        % Summary plot: LE vs RN by n_prev_P (row 2 of fig6)
        for ni = 1:N_npp
            ax6b = subplot(3,4,4+ni); hold(ax6b,'on');
            title(ax6b,npp_labels{ni},'FontSize',10);

            for sti = 1:N_stg6
                vals = rewp_by_npp(:,sti,ni);
                vals = vals(~isnan(vals));
                if isempty(vals), continue; end
                xpos = sti;
                bar(ax6b,xpos,mean(vals,'omitnan'),0.5,...
                    'FaceColor',NPP_COLORS(ni,:),'EdgeColor','none','FaceAlpha',0.8,...
                    'HandleVisibility','off');
                errorbar(ax6b,xpos,mean(vals,'omitnan'),std(vals,'omitnan')/sqrt(numel(vals)),...
                    'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');
                jx=xpos+0.15*(rand(numel(vals),1)-0.5);
                scatter(ax6b,jx,vals,20,NPP_COLORS(ni,:),'filled',...
                    'MarkerFaceAlpha',0.5,'HandleVisibility','off');
                text(ax6b,xpos,-0.02,sprintf('n=%d',numel(vals)),...
                    'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
            end
            for si6=1:numel(subjs6)
                v1=rewp_by_npp(si6,1,ni); v2=rewp_by_npp(si6,2,ni);
                if ~isnan(v1)&&~isnan(v2)
                    plot(ax6b,[1 2],[v1 v2],'-','Color',[0.5 0.5 0.5],...
                        'LineWidth',0.6,'HandleVisibility','off');
                end
            end
            v1a=rewp_by_npp(:,1,ni); v2a=rewp_by_npp(:,2,ni);
            ok_p=~isnan(v1a)&~isnan(v2a);
            if sum(ok_p)>2
                [~,p_p]=ttest(v1a(ok_p),v2a(ok_p));
                add_sig_bracket_p9(ax6b,1,2,max([v1a(ok_p);v2a(ok_p)],[],'omitnan')*1.1,p_p);
            end
            set(ax6b,'XTick',[1 2],'XTickLabel',{'LE','RN'},'TickDir','out');
            ylabel(ax6b,'RewP amplitude (µV)','FontSize',9);
            xlabel(ax6b,'Stage','FontSize',9);
            yline(ax6b,0,'k--','LineWidth',0.8,'HandleVisibility','off');
        end

        % Summary plot: trajectory LE→RN as n_prev_P increases (row 3)
        for sti = 1:N_stg6
            ax6c = subplot(3,4,8+sti); hold(ax6c,'on');
            title(ax6c,sprintf('%s — RewP by n_{prev P}',rev_stages{sti}),'FontSize',10);
            gm_npp = NaN(1,N_npp); gs_npp = NaN(1,N_npp); nn_npp = zeros(1,N_npp);
            for ni=1:N_npp
                vals=rewp_by_npp(:,sti,ni); vals=vals(~isnan(vals));
                if isempty(vals), continue; end
                gm_npp(ni)=mean(vals,'omitnan'); gs_npp(ni)=std(vals,'omitnan')/sqrt(numel(vals));
                nn_npp(ni)=numel(vals);
                scatter(ax6c,ni+0.2*(rand(numel(vals),1)-0.5),vals,...
                    18,NPP_COLORS(ni,:),'filled','MarkerFaceAlpha',0.35,'HandleVisibility','off');
            end
            ok_s=~isnan(gm_npp);
            if any(ok_s)
                fill(ax6c,[find(ok_s),fliplr(find(ok_s))],...
                    [gm_npp(ok_s)+gs_npp(ok_s),fliplr(gm_npp(ok_s)-gs_npp(ok_s))],...
                    STAGE_COLORS(sti+1,:),'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
                plot(ax6c,find(ok_s),gm_npp(ok_s),'o-','Color',STAGE_COLORS(sti+1,:),...
                    'MarkerFaceColor',STAGE_COLORS(sti+1,:),'MarkerSize',8,'LineWidth',2,...
                    'DisplayName',rev_stages{sti});
            end
            for ni=1:N_npp
                if nn_npp(ni)>0 && ~isnan(gm_npp(ni))
                    text(ax6c,ni,gm_npp(ni)-(gs_npp(ni)*1.5+0.05),...
                        sprintf('n=%d',nn_npp(ni)),...
                        'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
                end
            end
            yline(ax6c,0,'k--','LineWidth',0.8,'HandleVisibility','off');
            set(ax6c,'XTick',1:N_npp,'XTickLabel',npp_labels,'XTickLabelRotation',20,'TickDir','out');
            ylabel(ax6c,'RewP amplitude (µV)','FontSize',9);
            xlabel(ax6c,'Prior probabilistic exposure','FontSize',9);
            legend(ax6c,'Box','off','FontSize',8,'Location','best');
        end

        % Scatter: LE RewP vs RN RewP per subject, coloured by n_prev_P
        ax6d = subplot(3,4,11); hold(ax6d,'on');
        title(ax6d,'LE vs RN RewP per subject','FontSize',10);
        all_le = NaN(numel(subjs6),1); all_rn = NaN(numel(subjs6),1); all_ni = NaN(numel(subjs6),1);
        for si6=1:numel(subjs6)
            for ni=1:N_npp
                v1=rewp_by_npp(si6,1,ni); v2=rewp_by_npp(si6,2,ni);
                if ~isnan(v1)&&~isnan(v2)
                    all_le(si6)=v1; all_rn(si6)=v2; all_ni(si6)=ni; break
                end
            end
        end
        ok_sc=~isnan(all_le)&~isnan(all_rn);
        for ni=1:N_npp
            m=ok_sc&all_ni==ni;
            scatter(ax6d,all_le(m),all_rn(m),50,NPP_COLORS(ni,:),'filled','DisplayName',npp_labels{ni});
        end
        if any(ok_sc)
            lims_sc=[min([all_le(ok_sc);all_rn(ok_sc)])-0.1,...
                      max([all_le(ok_sc);all_rn(ok_sc)])+0.1];
            if range(lims_sc)>0
                plot(ax6d,lims_sc,lims_sc,'k--','LineWidth',0.8,'HandleVisibility','off');
                xlim(ax6d,lims_sc); ylim(ax6d,lims_sc);
            end
        else
            text(ax6d,0.5,0.5,'No paired LE/RN data','Units','normalized',...
                'HorizontalAlignment','center','Color',[0.5 0.5 0.5]);
        end
        xlabel(ax6d,'RewP at LE (µV)'); ylabel(ax6d,'RewP at RN (µV)');
        legend(ax6d,'Box','off','FontSize',8,'Location','best');
        axis(ax6d,'square');
        text(ax6d,0.04,0.96,'Points above diagonal: RN > LE','Units','normalized',...
            'VerticalAlignment','top','FontSize',7,'Color',[0.4 0.4 0.4]);

        % Legend for subplot 3,4,12 – colour key
        ax6e=subplot(3,4,12); axis(ax6e,'off');
        text(ax6e,0.05,0.95,'Colour key','FontWeight','bold','FontSize',10,...
            'Units','normalized','VerticalAlignment','top');
        items={'D→D','D→P','P→D','P→P','','0 prior P','1 prior P','2+ prior P'};
        clrs_key=[TRANS_COLORS;0.9 0.9 0.9;NPP_COLORS];
        for ki=1:numel(items)
            if isempty(items{ki}), continue; end
            rectangle(ax6e,'Position',[0.02 0.85-ki*0.10 0.12 0.08],...
                'FaceColor',clrs_key(ki,:),'EdgeColor','none');
            text(ax6e,0.18,0.89-ki*0.10,items{ki},'FontSize',9,'Units','normalized',...
                'VerticalAlignment','middle');
        end
    end
else
    text(0.5,0.5,'Per-trial table not available — cannot compute Fig 6',...
        'Units','normalized','HorizontalAlignment','center','FontSize',12);
end

annotation(fig6,'textbox',[0.01 0.01 0.98 0.04],'String',...
    ['Fig 6: RewP (mean amplitude in ' num2str(REWP_WIN(1)) '–' num2str(REWP_WIN(2)) ' ms window; correct minus incorrect) '...
     'at LE (pre-reversal) and RN (post-reversal). '...
     'Row 1: by block transition (D→D/D→P/P→D/P→P). Row 2: by n_{prev P} bin. '...
     'Row 3: trajectory across n_{prev P}. Lines connect same subject. Brackets = paired t-test. '...
     'Theory: prior uncertainty (n_{prev P}) should attenuate RewP especially at RN (reversal noise).'],...
    'FontSize',8,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);
save_fig_p9(fig6,outdir,'P9_Fig6_RewP_previous_uncertainty');
fprintf('  Fig 6 saved.\n');

%% ── DONE ────────────────────────────────────────────────────────────────────
fprintf('\n=== P9 COMPLETE ===\n');
fprintf('All figures saved to:\n  %s\n',outdir);
for fn={'P9_Fig1_RewP_by_stage_blocktype','P9_Fig2_ERP_by_condition',...
        'P9_Fig3_First10trials_DvsP','P9_Fig4_Mock_raw_EEG',...
        'P9_Fig5_P300_by_condition','P9_Fig6_RewP_previous_uncertainty'}
    fprintf('  %s.pdf / .png\n',fn{1});
end


%% ═══════════════════════════════════════════════════════════════════════════
%% LOCAL HELPER FUNCTIONS
%% ═══════════════════════════════════════════════════════════════════════════

function save_fig_p9(fig,outdir,fname)
    if ~exist(outdir,'dir'), mkdir(outdir); end
    try
        exportgraphics(fig,fullfile(outdir,[fname '.pdf']),'ContentType','vector');
        exportgraphics(fig,fullfile(outdir,[fname '.png']),'Resolution',300);
    catch
        saveas(fig,fullfile(outdir,[fname '.pdf']));
    end
    fprintf('  Saved: %s\n',fname);
end

function [mn,se,n] = extract_ga(gt,wave_col,mask,t_ax,n_t,in_erp)
    mn=NaN(1,n_t); se=NaN(1,n_t); n=0;
    if ~ismember(wave_col,gt.Properties.VariableNames), return; end
    rows = find(mask);
    if isempty(rows), return; end
    if ismember('subj_id',gt.Properties.VariableNames)
        subjs=unique(string(gt.subj_id(rows)));
        sc='subj_id';
    else
        subjs=unique(string(gt.subjID(rows)));
        sc='subjID';
    end
    sm=NaN(numel(subjs),n_t);
    for si=1:numel(subjs)
        sn=subjs(si);
        smask=mask & string(gt.(sc))==sn;
        waves=gt.(wave_col)(smask);
        waves=waves(~cellfun(@isempty,waves));
        if isempty(waves), continue; end
        M=cell2mat(cellfun(@(v)v(:)',waves,'UniformOutput',false));
        if size(M,2)~=numel(t_ax)
            t2=linspace(t_ax(1),t_ax(end),size(M,2));
            M2=NaN(size(M,1),numel(t_ax));
            for ri=1:size(M,1), M2(ri,:)=interp1(t2,M(ri,:),t_ax,'linear',NaN); end
            M=M2;
        end
        bl=t_ax>=-200&t_ax<=0;
        if any(bl), M=M-mean(M(:,bl),2,'omitnan'); end
        sm(si,:)=mean(M(:,in_erp),1,'omitnan');
    end
    ok=~all(isnan(sm),2);
    n=sum(ok);
    if n==0,return;end
    mn=mean(sm(ok,:),1,'omitnan');
    se=std( sm(ok,:),0,1,'omitnan')/sqrt(n);
end

function M = stack_waves_p9(dw_cell,t_ax,sign_mult)
    M=cell2mat(cellfun(@(v)sign_mult*v(:)',dw_cell,'UniformOutput',false));
    if size(M,2)~=numel(t_ax)
        t2=linspace(t_ax(1),t_ax(end),size(M,2));
        M2=NaN(size(M,1),numel(t_ax));
        for ri=1:size(M,1), M2(ri,:)=interp1(t2,M(ri,:),t_ax,'linear',NaN); end
        M=M2;
    end
end



function shade_win(ax,win,yl_hint,rgb)
    yl=ylim(ax);
    if abs(yl(1))<1e-6 && abs(yl(2)-1)<1e-6, yl=yl_hint; end
    patch(ax,[win(1) win(2) win(2) win(1)],[yl(1) yl(1) yl(2) yl(2)],...
        rgb,'EdgeColor','none','FaceAlpha',0.45,'HandleVisibility','off');
end

function yl = compute_yl(mn,se)
    lo = mn - se;
    hi = mn + se;
    vals = [lo(:); hi(:)];
    vals = vals(~isnan(vals));
    if isempty(vals)
        yl = [-1 2];
        return;
    end
    yl = [min(vals)-0.3, max(vals)+0.3];
    if range(yl) < 0.01
        yl = [-1 2];
    end
end

function out = ternary_p9(c,a,b)
    if c, out=a; else, out=b; end
end

function bs = infer_block_structure(td,nB)
    bs=repmat('D',1,nB);
    if isfield(td,'block_structure')&&~isempty(td.block_structure)
        bs=upper(char(td.block_structure));
        bs(bs=='V')='P';
    elseif isfield(td,'trueFB')
        for b=1:nB
            pfb=td.trueFB(b,~isnan(td.trueFB(b,:)));
            if ~isempty(pfb)&&mean(pfb)<0.99, bs(b)='P'; end
        end
    end
end

function gt = derive_prev_uncertainty(gt)
    % Add n_prev_P and transition_cat columns to per-trial table
    gt.n_prev_P      = nan(height(gt),1);
    gt.transition_cat = repmat("first",height(gt),1);
    sc = 'subj_id';
    if ~ismember(sc,gt.Properties.VariableNames), sc='subjID'; end
    bc = 'block';
    if ~ismember(bc,gt.Properties.VariableNames)&&ismember('blocknum',gt.Properties.VariableNames)
        bc='blocknum';
    end
    subjs=unique(string(gt.(sc)));
    for si=1:numel(subjs)
        sn=subjs(si);
        sm=string(gt.(sc))==sn;
        blks=sort(unique(double(gt.(bc)(sm))))';
        bt_arr=strings(1,numel(blks));
        for bi=1:numel(blks)
            b=blks(bi); bm=sm & double(gt.(bc))==b;
            bts=unique(gt.block_type_s(bm));
            bt_arr(bi)=bts(1);
        end
        n_prev_P_run=0;
        for bi=1:numel(blks)
            b=blks(bi); bm=sm & double(gt.(bc))==b;
            gt.n_prev_P(bm)=n_prev_P_run;
            if bi==1
                gt.transition_cat(bm)="first";
            else
                prev=bt_arr(bi-1); curr=bt_arr(bi);
                gt.transition_cat(bm)=prev+"->"+curr;
            end
            if bt_arr(bi)=="P", n_prev_P_run=n_prev_P_run+1; end
        end
    end
    gt.transition_cat=string(gt.transition_cat);
end

function sm = get_subj_mask(gt,sn)
    if ismember('subj_id',gt.Properties.VariableNames)
        sm=string(gt.subj_id)==sn;
    else
        sm=string(gt.subjID)==sn;
    end
end

function rv = compute_rewp_contrast(gt,mask,use_wave,use_scalar_rewp,...
    use_scalar_frn,rewp_mask,in_erp,t_ax,wave_col_use)
% Compute mean RewP amplitude (correct - incorrect) for trials in mask.
% Returns per-subject mean RewP in µV.
    rv = NaN;
    mask_cor = mask & gt.correct_num==1;
    mask_inc = mask & gt.correct_num==0;
    if sum(mask_cor)<2 || sum(mask_inc)<2, return; end

    if use_wave
        wc = wave_col_use;
        wc_cor=gt.(wc)(mask_cor); wc_cor=wc_cor(~cellfun(@isempty,wc_cor));
        wc_inc=gt.(wc)(mask_inc); wc_inc=wc_inc(~cellfun(@isempty,wc_inc));
        if isempty(wc_cor)||isempty(wc_inc), rv=NaN; return; end
        M_cor=cell2mat(cellfun(@(v)v(:)',wc_cor,'UniformOutput',false));
        M_inc=cell2mat(cellfun(@(v)v(:)',wc_inc,'UniformOutput',false));
        if size(M_cor,2)~=numel(t_ax)
            t2=linspace(t_ax(1),t_ax(end),size(M_cor,2));
            M2=NaN(size(M_cor,1),numel(t_ax));
            for ri=1:size(M_cor,1), M2(ri,:)=interp1(t2,M_cor(ri,:),t_ax,'linear',NaN); end
            M_cor=M2;
        end
        if size(M_inc,2)~=numel(t_ax)
            t2=linspace(t_ax(1),t_ax(end),size(M_inc,2));
            M2=NaN(size(M_inc,1),numel(t_ax));
            for ri=1:size(M_inc,1), M2(ri,:)=interp1(t2,M_inc(ri,:),t_ax,'linear',NaN); end
            M_inc=M2;
        end
        bl = t_ax >= -200 & t_ax <= 0;
        if any(bl)
            M_cor = M_cor - mean(M_cor(:,bl),2,'omitnan');
            M_inc = M_inc - mean(M_inc(:,bl),2,'omitnan');
        end
        win_idx = rewp_mask;
        if ~any(win_idx), return; end
        mn_cor = mean(M_cor(:,win_idx), 'omitnan');
        mn_cor = mean(mn_cor, 'omitnan');
        mn_inc = mean(M_inc(:,win_idx), 'omitnan');
        mn_inc = mean(mn_inc, 'omitnan');
        rv = mn_cor - mn_inc;
    elseif use_scalar_rewp
        rv = mean(gt.RewP_mean_amp(mask_cor),'omitnan') - ...
             mean(gt.RewP_mean_amp(mask_inc),'omitnan');
    elseif use_scalar_frn
        % FRN is more negative for incorrect; RewP proxy = negative of FRN difference
        rv = -(mean(gt.FRN_mean_amp(mask_cor),'omitnan') - ...
               mean(gt.FRN_mean_amp(mask_inc),'omitnan'));
    end
end

function add_sig_bracket_p9(ax,x1,x2,y_top,p_val)
    if isnan(p_val)||isnan(y_top), return; end
    if     p_val<0.001, s='***';
    elseif p_val<0.01,  s='**';
    elseif p_val<0.05,  s='*';
    else,               s='ns'; end
    yb=y_top;
    line(ax,[x1 x1 x2 x2],[yb*0.96 yb yb yb*0.96],'Color','k','LineWidth',0.8,'HandleVisibility','off');
    text(ax,mean([x1 x2]),yb*1.01,s,'HorizontalAlignment','center','FontSize',10,'HandleVisibility','off');
end


function plot_ribbon_p9(ax,t,mn,se,clr,ls,lbl)
    % Plot mean +/- SEM ribbon. Uses only finite points so fill() does not fail.
    if isempty(mn) || isempty(se) || isempty(t), return; end
    ok = ~isnan(mn) & ~isnan(se) & ~isnan(t);
    if ~any(ok), return; end
    t_ok = t(ok);
    mn_ok = mn(ok);
    se_ok = se(ok);
    fill(ax,[t_ok fliplr(t_ok)],[mn_ok+se_ok fliplr(mn_ok-se_ok)], ...
        clr,'FaceAlpha',0.18,'EdgeColor','none','HandleVisibility','off');
    plot(ax,t_ok,mn_ok,'Color',clr,'LineStyle',ls,'LineWidth',2.0,'DisplayName',lbl);
end
