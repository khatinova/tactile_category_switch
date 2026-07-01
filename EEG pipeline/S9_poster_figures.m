% S9_poster_figures.m — Generate all poster panels (ENHANCED AESTHETICS + STATS)
% "Uncertainty calibrates behavioural and neural adaptation to category switching"
% Requires: S5 (nassar_results.mat), S5b (model_comparison_RW_HGF.mat),
%           S4 (group_feature_table_combined.mat)
% Focuses on stochasticity (P blocks / false FB) vs volatility (LE→RN reversal)
%
% PLOT DIVERSITY:
%   Fig 1: Behaviour — spaghetti + group CI, dot-CI + subject lines, slope
%   Fig 2: Modelling — BIC scatter + Nassar schematic, D/P surprise ribbons,
%                       paired raincloud for H
%   Fig 3B: FRN — raincloud, paired strip
%   Fig 3C: P300 — grouped bar + jitter, binned ribbon scatter
%   Fig 3D: Theta — violin + box overlay, binned ribbon scatter
%   Fig 3E: PLV — dot-CI pathway comparison, dual ribbon scatter
%   Fig 3F: Neural × n_prev_P — slope + ribbon
%
% All panels annotated with paired t / ANOVA stats + significance stars.
% ═══════════════════════════════════════════════════════════════════════════
close all; clc;
addpath(genpath(fileparts(mfilename('fullpath'))));
putils=fullfile(fileparts(fileparts(mfilename('fullpath'))),'pipeline','utils');
if exist(putils,'dir'),addpath(putils);end
base_path='\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch';
stf=fullfile(base_path,'Salient mod switch KH','Results','EEG analysis','Outcome_feature_tables_v4_merged');
figdir=fullfile(base_path,'Salient mod switch KH','Results','EEG analysis','Figures','Poster_2026');
if ~exist(figdir,'dir'),mkdir(figdir);end

%% LOAD DATA
if exist('gt','var')&&istable(gt), fprintf('Using gt already in workspace.\n');
elseif exist('group_table','var')&&istable(group_table), gt=group_table;
else, load(fullfile(stf,'group_feature_table_combined.mat'),'group_table'); gt=group_table; end
gt.subj_id=categorical(gt.subj_id); gt.block_type=categorical(gt.block_type);
gt.stage=categorical(gt.stage,{'LN','LE','RN','RE'},'Ordinal',true);

% DERIVE false_fb FROM trueFB (matching S8 logic) — this is why the panel was empty
if ~ismember('false_fb',gt.Properties.VariableNames) && ismember('trueFB',gt.Properties.VariableNames)
    gt.false_fb = double(~gt.trueFB);
    fprintf('Derived false_fb from trueFB column.\n');
elseif ismember('false_fb',gt.Properties.VariableNames)
    gt.false_fb = double(gt.false_fb);
end

subj_list=unique(gt.subj_id); N_subj=numel(subj_list);

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

%% MERGE NASSAR LATENTS (same as S8)
if ~ismember('surprise',gt.Properties.VariableNames)||sum(~isnan(gt.surprise))<10
fprintf('Merging Nassar latents...\n');
nrc={fullfile(base_path,'Salient mod switch KH','Results','Simulation results','Figures','nassar_results.mat')};
if exist(nrc{1},'file')
tmp=load(nrc{1},'results');results=tmp.results;sr=fieldnames(results);
ci2={'PE_nassar','PE_unsigned','omega','alpha_nassar','certainty','surprise','theta_nassar'};
for c=ci2,if ~ismember(c{1},gt.Properties.VariableNames),gt.(c{1})=nan(height(gt),1);end,end
sc='subj_id';bc='block';tc='trial';
if ~ismember(bc,gt.Properties.VariableNames),bc='block_number';end
if ~ismember(tc,gt.Properties.VariableNames),tc='trialnum';end
nm=0;
for si=1:numel(sr),sn=sr{si};r=results.(sn);
    if ~isfield(r,'delta_trial'),continue;end
    rows=find(string(gt.(sc))==string(sn));if isempty(rows),continue;end
    bg=double(gt.(bc)(rows));tg=double(gt.(tc)(rows));
    for t=1:numel(r.trial_id)
        m=rows(bg==r.block_id(t)&tg==r.trial_id(t));if isempty(m),continue;end
        gt.PE_nassar(m(1))=r.delta_trial(t);gt.PE_unsigned(m(1))=abs(r.delta_trial(t));
        gt.omega(m(1))=r.omega_trial(t);gt.alpha_nassar(m(1))=r.alpha_trial(t);
        gt.certainty(m(1))=r.certainty_trial(t);gt.surprise(m(1))=r.surprise(t);
        gt.theta_nassar(m(1))=r.theta_trial(t);nm=nm+1;
    end
end
fprintf('  Merged %d trials.\n',nm);
end
end

%% Z-SCORE WITHIN SUBJECT
vz={'PE_nassar','PE_unsigned','omega','alpha_nassar','certainty','surprise','theta_nassar',...
    'prefrontal_mean_norm','P300_norm','Theta_amp','PLV_fp','PLV_fs'};
for f=1:numel(vz),fn=vz{f};fnz=[fn '_z'];
if ~ismember(fn,gt.Properties.VariableNames),continue;end
gt.(fnz)=nan(height(gt),1);
for si=1:numel(subj_list),mask=gt.subj_id==subj_list(si);v=gt.(fn)(mask);sd=std(v,'omitnan');
if sd>0,gt.(fnz)(mask)=(v-mean(v,'omitnan'))/sd;end,end,end

%% ADD n_prev_P IF MISSING
if ~ismember('n_prev_P',gt.Properties.VariableNames)
    gt = add_transition_history_columns_local(gt);
end

% ═══════════════════════════════════════════════════════════════════════════
% COLOUR PALETTE — poster-consistent, colourblind-safe
% ═══════════════════════════════════════════════════════════════════════════
CLR_D=[.15 .45 .70]; CLR_P=[.80 .30 .10]; CLR_T=[.20 .60 .30]; CLR_F=[.75 .20 .55];
CLR_FP=[.15 .45 .70]; CLR_FS=[.60 .20 .55];
CLR_D_LIGHT=[.55 .75 .90]; CLR_P_LIGHT=[.95 .70 .55];
STAGES={'LN','LE','RN','RE'}; STAGE_X=1:4;
fprintf('\n=== POSTER FIGURES (ENHANCED + STATS) ===\n');

%% ═══════════════════════════════════════════════════════════════════════════
%  FIG 1: BEHAVIOUR
%  ═══════════════════════════════════════════════════════════════════════════
fig1=figure('Position',[30 30 1600 480],'Color','w');

% --- 1A: Accuracy × Stage × Block type with subject spaghetti ---
subplot(1,3,1); hold on
for si=1:N_subj
    acc_d=nan(1,4); acc_p=nan(1,4);
    for sti=1:4
        md=gt.subj_id==subj_list(si)&gt.block_type=='D'&gt.stage==STAGES{sti};
        mp=gt.subj_id==subj_list(si)&gt.block_type=='P'&gt.stage==STAGES{sti};
        acc_d(sti)=mean(gt.correct(md),'omitnan');
        acc_p(sti)=mean(gt.correct(mp),'omitnan');
    end
    plot(STAGE_X,acc_d,'-','Color',[CLR_D .12],'LineWidth',.6,'HandleVisibility','off');
    plot(STAGE_X,acc_p,'-','Color',[CLR_P .12],'LineWidth',.6,'HandleVisibility','off');
end
[m_d,se_d]=subj_means_by_stage(gt,'correct','D',subj_list);
[m_p,se_p]=subj_means_by_stage(gt,'correct','P',subj_list);
errorbar(STAGE_X-.04,m_d,se_d,'o-','Color',CLR_D,'LineWidth',2.5,...
    'MarkerFaceColor',CLR_D,'MarkerSize',9,'CapSize',5,'DisplayName','Deterministic');
errorbar(STAGE_X+.04,m_p,se_p,'o-','Color',CLR_P,'LineWidth',2.5,...
    'MarkerFaceColor',CLR_P,'MarkerSize',9,'CapSize',5,'DisplayName','Probabilistic');
set(gca,'XTick',STAGE_X,'XTickLabel',STAGES); xlim([.5 4.5]); ylim([.3 1]);
yline(.5,':','Color',[.5 .5 .5],'LineWidth',1,'HandleVisibility','off');
% Stats: LE vs RN drop for D
[acc_d_subj,~]=subj_vals_by_stage(gt,'correct','D',subj_list);
[~,p_rev_d]=ttest(acc_d_subj(:,2),acc_d_subj(:,3));
text(2.5,.38,sprintf('LE\\rightarrowRN (D): %s',pstar(p_rev_d)),'FontSize',9,'Color',CLR_D);
xlabel('Task stage'); ylabel('P(correct)');
title('A   Accuracy \times stage','FontSize',12);
legend('Box','off','Location','southwest','FontSize',9);

% --- 1B: Confidence × Stage ---
subplot(1,3,2); hold on
for si=1:N_subj
    conf_d=nan(1,4); conf_p=nan(1,4);
    for sti=1:4
        md=gt.subj_id==subj_list(si)&gt.block_type=='D'&gt.stage==STAGES{sti};
        mp=gt.subj_id==subj_list(si)&gt.block_type=='P'&gt.stage==STAGES{sti};
        conf_d(sti)=mean(gt.confidence(md),'omitnan');
        conf_p(sti)=mean(gt.confidence(mp),'omitnan');
    end
    plot(STAGE_X,conf_d,'-','Color',[CLR_D .12],'LineWidth',.6,'HandleVisibility','off');
    plot(STAGE_X,conf_p,'-','Color',[CLR_P .12],'LineWidth',.6,'HandleVisibility','off');
end
[m_d,se_d]=subj_means_by_stage(gt,'confidence','D',subj_list);
[m_p,se_p]=subj_means_by_stage(gt,'confidence','P',subj_list);
errorbar(STAGE_X-.04,m_d,se_d,'s-','Color',CLR_D,'LineWidth',2.5,...
    'MarkerFaceColor',CLR_D,'MarkerSize',9,'CapSize',5,'DisplayName','D');
errorbar(STAGE_X+.04,m_p,se_p,'s-','Color',CLR_P,'LineWidth',2.5,...
    'MarkerFaceColor',CLR_P,'MarkerSize',9,'CapSize',5,'DisplayName','P');
set(gca,'XTick',STAGE_X,'XTickLabel',STAGES); xlim([.5 4.5]);
% Stats: D vs P overall
[conf_d_subj,~]=subj_vals_by_stage(gt,'confidence','D',subj_list);
[conf_p_subj,~]=subj_vals_by_stage(gt,'confidence','P',subj_list);
[~,p_dp]=ttest(mean(conf_d_subj,2,'omitnan'),mean(conf_p_subj,2,'omitnan'));
yl=ylim; text(2.5,yl(1)+.05*range(ylim),sprintf('D vs P: %s',pstar(p_dp)),...
    'FontSize',9,'Color',[.3 .3 .3]);
xlabel('Task stage'); ylabel('Confidence (1–10)');
title('B   Confidence \times stage','FontSize',12);
legend('Box','off','Location','southwest','FontSize',9);

% --- 1C: Reversal cost by n_prev_P ---
subplot(1,3,3); hold on
if ismember('n_prev_P',gt.Properties.VariableNames)
    npps=[0 1 2]; npp_x=1:3; npp_labels={'0','1','2+'};
    rev_cost_d=nan(N_subj,3); rev_cost_p=nan(N_subj,3);
    for si=1:N_subj
        for ni=1:3
            if ni<3, nm=gt.subj_id==subj_list(si)&gt.n_prev_P==npps(ni);
            else, nm=gt.subj_id==subj_list(si)&gt.n_prev_P>=2; end
            for bt_i=1:2
                bt=ternary_pf(bt_i==1,'D','P');
                le_m=nm&gt.stage=='LE'&gt.block_type==bt;
                rn_m=nm&gt.stage=='RN'&gt.block_type==bt;
                le_v=mean(gt.correct(le_m),'omitnan');
                rn_v=mean(gt.correct(rn_m),'omitnan');
                if bt_i==1, rev_cost_d(si,ni)=le_v-rn_v;
                else, rev_cost_p(si,ni)=le_v-rn_v; end
            end
        end
    end
    m_rd=mean(rev_cost_d,'omitnan'); se_rd=std(rev_cost_d,'omitnan')./sqrt(sum(~isnan(rev_cost_d)));
    m_rp=mean(rev_cost_p,'omitnan'); se_rp=std(rev_cost_p,'omitnan')./sqrt(sum(~isnan(rev_cost_p)));
    fill([npp_x fliplr(npp_x)],[m_rd fliplr(m_rp)],[.85 .85 .85],...
        'FaceAlpha',.25,'EdgeColor','none','HandleVisibility','off');
    errorbar(npp_x-.05,m_rd,se_rd,'o-','Color',CLR_D,'LineWidth',2.5,...
        'MarkerFaceColor',CLR_D,'MarkerSize',9,'CapSize',5,'DisplayName','D');
    errorbar(npp_x+.05,m_rp,se_rp,'o-','Color',CLR_P,'LineWidth',2.5,...
        'MarkerFaceColor',CLR_P,'MarkerSize',9,'CapSize',5,'DisplayName','P');
    set(gca,'XTick',npp_x,'XTickLabel',npp_labels); xlim([.5 3.5]);
    % Stats: linear trend
    ok_lin=~isnan(rev_cost_p(:,1))&~isnan(rev_cost_p(:,3));
    if sum(ok_lin)>5
        [~,p_trend]=ttest(rev_cost_p(ok_lin,1),rev_cost_p(ok_lin,3));
        yl=ylim; text(2,yl(2)-.05*range(ylim),sprintf('P 0vs2+: %s',pstar(p_trend)),...
            'FontSize',9,'Color',CLR_P,'HorizontalAlignment','center');
    end
    xlabel('Prior P-block exposure'); ylabel('Reversal cost (LE − RN)');
end
title('C   Prior uncertainty recalibrates cost','FontSize',12);
legend('Box','off','Location','northeast','FontSize',9);
style_save(fig1,fullfile(figdir,'Fig1_behaviour'));

%% ═══════════════════════════════════════════════════════════════════════════
%  FIG 2: COMPUTATIONAL MODELLING
%  2A: Nassar model schematic (text-based introduction)
%  2B: Model comparison BIC scatter (Nassar vs RW — loaded from S5b output)
%  2C: Surprise at reversal (D vs P ribbons)
%  2D: Hazard rate paired raincloud
%  ═══════════════════════════════════════════════════════════════════════════
fig2=figure('Position',[30 30 1800 460],'Color','w');

% Load nassar results
nr_path=fullfile(base_path,'Salient mod switch KH','Results','Simulation results','Figures','nassar_results.mat');
has_nassar=exist(nr_path,'file');
if has_nassar
    tmp=load(nr_path,'results'); results=tmp.results; sr=fieldnames(results);
    N_r=numel(sr);
    H_det=nan(N_r,1); H_prob=nan(N_r,1); bic_n=nan(N_r,1);
    for si=1:N_r, r=results.(sr{si});
        if isfield(r,'H_fit_det'), H_det(si)=r.H_fit_det; end
        if isfield(r,'H_fit_prob'), H_prob(si)=r.H_fit_prob; end
        if isfield(r,'bic'), bic_n(si)=r.bic; end
    end
end

% Load RW BICs from S5b (model_comparison_RW_HGF.mat)
mc_path=fullfile(base_path,'Salient mod switch KH','Results','Simulation results','Figures','model_comparison_RW_HGF.mat');
bic_rw=nan(N_r,1);
if exist(mc_path,'file')
    mc=load(mc_path);
    % S5b saves: res (table with bic_nassar, bic_rw1, bic_rwdual, bic_hgf)
    %            bic_mat (N×4 matrix), model_lbl
    if isfield(mc,'bic_mat')
        % Column 2 = single-alpha RW (the standard comparison)
        bic_rw_raw = mc.bic_mat(:,2);
        % Match subjects by order (S5 and S5b iterate fieldnames(all_trial_data) identically)
        if numel(bic_rw_raw)>=N_r
            bic_rw = bic_rw_raw(1:N_r);
        else
            bic_rw(1:numel(bic_rw_raw)) = bic_rw_raw;
        end
    elseif isfield(mc,'res') && istable(mc.res) && ismember('bic_rw1',mc.res.Properties.VariableNames)
        bic_rw_raw = mc.res.bic_rw1;
        if numel(bic_rw_raw)>=N_r
            bic_rw = bic_rw_raw(1:N_r);
        else
            bic_rw(1:numel(bic_rw_raw)) = bic_rw_raw;
        end
    end
    fprintf('Loaded RW BICs from model_comparison_RW_HGF.mat (%d subjects)\n',sum(~isnan(bic_rw)));
else
    fprintf('WARNING: model_comparison_RW_HGF.mat not found — run S5b first.\n');
    fprintf('  Expected at: %s\n',mc_path);
end

% --- 2A: Nassar model schematic (text description for poster) ---
subplot(1,4,1); axis off; hold on
box_txt = {'\bf Nassar (2010) Change-Point Model',...
    '',...
    '\rm Belief update:  \theta_{t+1} = \theta_t + \alpha_t \cdot \delta_t',...
    '',...
    'Adaptive learning rate:',...
    '  \alpha_t = \omega_t + (1-\omega_t) / n_{eff}',...
    '',...
    'Change-point probability:',...
    '  \omega_t = H \cdot \chi_t / [H\cdot\chi_t + (1-H)]',...
    '',...
    'Surprise = \omega_t \times |\delta_t|',...
    '',...
    '\it Fitted per subject: H (hazard), \beta (softmax)',...
    '\it Latents: \omega, \alpha, surprise, PE'};
text(.05,.95,box_txt,'Units','normalized','VerticalAlignment','top',...
    'FontSize',9,'FontName','Arial','Interpreter','tex');
title('A   Nassar model','FontSize',12);

% --- 2B: Model comparison BIC scatter ---
subplot(1,4,2); hold on; axis square
ok=~isnan(bic_n)&~isnan(bic_rw);
if any(ok)
    winner_clr=repmat(CLR_D,sum(ok),1);
    nassar_wins=bic_n(ok)<bic_rw(ok);
    winner_clr(~nassar_wins,:)=repmat(CLR_P,sum(~nassar_wins),1);
    scatter(bic_rw(ok),bic_n(ok),70,winner_clr,'filled','MarkerFaceAlpha',.75);
    lims=[min([bic_n(ok);bic_rw(ok)])-20 max([bic_n(ok);bic_rw(ok)])+20];
    plot(lims,lims,'k--','LineWidth',1.2,'HandleVisibility','off');
    xlim(lims);ylim(lims);
    n_nw=sum(nassar_wins); n_tot=sum(ok);
    % Wilcoxon signed-rank test on BIC difference
    [p_mc,~]=signrank(bic_n(ok),bic_rw(ok));
    text(.05,.95,{sprintf('Nassar wins: %d/%d',n_nw,n_tot),...
        sprintf('RW wins: %d/%d',n_tot-n_nw,n_tot),...
        sprintf('signrank %s',pstar(p_mc))},...
        'Units','normalized','VerticalAlignment','top','FontSize',9,...
        'BackgroundColor',[1 1 1 .8]);
    xlabel('RW BIC'); ylabel('Nassar BIC');
else
    text(.5,.5,{'RW BIC not available','Run S5b_models_RW_HGF.m first'},...
        'HorizontalAlignment','center','Units','normalized','FontSize',10);
end
title('B   Model comparison','FontSize',12);

% --- 2C: Surprise aligned to reversal — D vs P ribbons (shallow/wide) ---
ax2c=subplot(1,4,3); hold on
% Make this panel wider and shallower
pos=get(ax2c,'Position'); set(ax2c,'Position',[pos(1)-.02 pos(2)+.08 pos(3)+.04 pos(4)-.12]);
half_w=15; x_ali=-half_w:half_w;
all_surp_d=[]; all_surp_p=[];
atd_path=fullfile(base_path,'Salient mod switch KH','Data','all_trial_data.mat');
if exist(atd_path,'file') && has_nassar
load(atd_path,'all_trial_data');
for si=1:N_r, sn=sr{si}; r=results.(sn);
    if ~isfield(all_trial_data,sn),continue;end
    td=all_trial_data.(sn).trial_data;if ~isfield(td,'revTrial'),continue;end
    nB=size(td.correct,1);
    for b=1:nB
        rev=td.revTrial(b);if isnan(rev),continue;end
        bm=r.block_id==b;bt_ids=r.trial_id(bm);bs=r.surprise(bm);
        row=NaN(1,numel(x_ali));
        for xi=1:numel(x_ali),ta=round(rev)+x_ali(xi);mi=find(bt_ids==ta,1);
            if ~isempty(mi),row(xi)=bs(mi);end
        end
        % Determine block type
        if isfield(td,'block_type')
            bt_char=td.block_type(b);
            if iscell(bt_char),bt_char=bt_char{1};end
            if isstring(bt_char),bt_char=char(bt_char);end
        elseif isfield(td,'block_is_det')
            if td.block_is_det(b),bt_char='D';else,bt_char='P';end
        else, bt_char='D'; end
        if upper(bt_char(1))=='D'
            all_surp_d(end+1,:)=row; %#ok
        else
            all_surp_p(end+1,:)=row; %#ok
        end
    end
end
end
% Plot ribbons
if ~isempty(all_surp_d)
    mn_d=mean(all_surp_d,1,'omitnan');se_d=std(all_surp_d,0,1,'omitnan')./sqrt(size(all_surp_d,1));
    fill([x_ali fliplr(x_ali)],[mn_d+se_d fliplr(mn_d-se_d)],CLR_D,...
        'FaceAlpha',.18,'EdgeColor','none','HandleVisibility','off');
    plot(x_ali,mn_d,'Color',CLR_D,'LineWidth',2.5,'DisplayName','Deterministic');
end
if ~isempty(all_surp_p)
    mn_p=mean(all_surp_p,1,'omitnan');se_p=std(all_surp_p,0,1,'omitnan')./sqrt(size(all_surp_p,1));
    fill([x_ali fliplr(x_ali)],[mn_p+se_p fliplr(mn_p-se_p)],CLR_P,...
        'FaceAlpha',.18,'EdgeColor','none','HandleVisibility','off');
    plot(x_ali,mn_p,'Color',CLR_P,'LineWidth',2.5,'DisplayName','Probabilistic');
end
xline(0,'--','Color',[.3 .3 .3],'LineWidth',1.5,'HandleVisibility','off');
% Stats: surprise at reversal trial (t=0) D vs P
if ~isempty(all_surp_d)&&~isempty(all_surp_p)
    idx0=find(x_ali==0);
    [~,p_surp]=ttest2(all_surp_d(:,idx0),all_surp_p(:,idx0));
    yl=ylim; text(1,yl(2)-.05*range(ylim),sprintf('D vs P at rev: %s',pstar(p_surp)),...
        'FontSize',9,'Color',[.3 .3 .3]);
end
xlabel('Trial relative to reversal'); ylabel('Surprise (\omega\times|\delta|)');
title('C   Surprise at reversal','FontSize',12);
legend('Box','off','Location','northeast','FontSize',9);

% --- 2D: H_det vs H_prob — paired raincloud ---
subplot(1,4,4); hold on
ok2=~isnan(H_det)&~isnan(H_prob);
if any(ok2)
    hd=H_det(ok2); hp=H_prob(ok2);
    % Generate jitter FIRST so lines can connect to actual dot positions
    jit1=1+.08*randn(numel(hd),1); jit2=2+.08*randn(numel(hp),1);
    % Paired lines connecting to jittered positions
    for i=1:numel(hd)
        plot([jit1(i) jit2(i)],[hd(i) hp(i)],'-','Color',[.6 .6 .6 .35],'LineWidth',.8,'HandleVisibility','off');
    end
    % Jittered dots
    scatter(jit1,hd,45,CLR_D,'filled','MarkerFaceAlpha',.6,'HandleVisibility','off');
    scatter(jit2,hp,45,CLR_P,'filled','MarkerFaceAlpha',.6,'HandleVisibility','off');
    % Mean + CI
    plot(1,mean(hd),'d','Color',CLR_D,'MarkerFaceColor',CLR_D,'MarkerSize',14,'LineWidth',1.5,'DisplayName','H_{det}');
    plot(2,mean(hp),'d','Color',CLR_P,'MarkerFaceColor',CLR_P,'MarkerSize',14,'LineWidth',1.5,'DisplayName','H_{prob}');
    errorbar(1,mean(hd),std(hd)/sqrt(numel(hd)),'Color',CLR_D,'LineWidth',2,'CapSize',8,'HandleVisibility','off');
    errorbar(2,mean(hp),std(hp)/sqrt(numel(hp)),'Color',CLR_P,'LineWidth',2,'CapSize',8,'HandleVisibility','off');
    % Stats
    [~,p_hp]=ttest(hp-hd);
    dH=mean(hp-hd); d_cohen=dH/std(hp-hd);
    text(.5,.95,sprintf('\\DeltaH=%.3f, d=%.2f\n%s',dH,d_cohen,pstar(p_hp)),...
        'Units','normalized','HorizontalAlignment','center','VerticalAlignment','top',...
        'FontSize',9,'BackgroundColor',[1 1 1 .8]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Deterministic','Probabilistic'});
    ylabel('Fitted hazard rate (H)'); xlim([.4 2.6]);
    legend('Box','off','Location','northwest','FontSize',9);
end
title('D   Noise sensitivity','FontSize',12);
style_save(fig2,fullfile(figdir,'Fig2_modelling'));

%% ═══════════════════════════════════════════════════════════════════════════
%  FIG 3A: RAW ERP — PLACEHOLDER
%  ═══════════════════════════════════════════════════════════════════════════
fig3a=figure('Position',[30 30 600 400],'Color','w');
axis off;
text(.5,.5,{'Fig 3A: Raw ERP waveforms','(insert from exported ERP figure)'},...
    'HorizontalAlignment','center','FontSize',14,'FontWeight','bold');
style_save(fig3a,fullfile(figdir,'Fig3A_ERP_placeholder'));

%% ═══════════════════════════════════════════════════════════════════════════
%  FIG 3B: FRN — RAINCLOUD (volatility) + PAIRED STRIP (stochasticity)
%  ═══════════════════════════════════════════════════════════════════════════
fig3b=figure('Position',[30 30 1100 480],'Color','w');

feat='prefrontal_mean_norm';

% LEFT: Raincloud — FRN at LE vs RN, D and P
subplot(1,2,1); hold on
frn_d_le=nan(N_subj,1);frn_d_rn=nan(N_subj,1);
frn_p_le=nan(N_subj,1);frn_p_rn=nan(N_subj,1);
for si=1:N_subj
    frn_d_le(si)=mean(gt.(feat)(gt.subj_id==subj_list(si)&gt.block_type=='D'&gt.stage=='LE'),'omitnan');
    frn_d_rn(si)=mean(gt.(feat)(gt.subj_id==subj_list(si)&gt.block_type=='D'&gt.stage=='RN'),'omitnan');
    frn_p_le(si)=mean(gt.(feat)(gt.subj_id==subj_list(si)&gt.block_type=='P'&gt.stage=='LE'),'omitnan');
    frn_p_rn(si)=mean(gt.(feat)(gt.subj_id==subj_list(si)&gt.block_type=='P'&gt.stage=='RN'),'omitnan');
end
% Plot as half-violin + dots + mean
grps={frn_d_le,frn_d_rn,frn_p_le,frn_p_rn};
clrs={CLR_D_LIGHT,CLR_D,CLR_P_LIGHT,CLR_P};
xpos=[1 2 3.5 4.5]; xlbls={'LE','RN','LE','RN'};
for gi=1:4
    v=grps{gi}; v=v(~isnan(v));
    if numel(v)<3,continue;end
    % Half-violin (kernel density)
    [kde_f,kde_xi]=ksdensity(v,'NumPoints',50);
    kde_f=kde_f/max(kde_f)*.35;
    fill([xpos(gi)-kde_f fliplr(xpos(gi)-kde_f*0)],[kde_xi fliplr(kde_xi)],...
        clrs{gi},'FaceAlpha',.3,'EdgeColor','none','HandleVisibility','off');
    % Jittered dots
    jx=xpos(gi)+.05+.06*randn(numel(v),1);
    scatter(jx,v,20,clrs{gi},'filled','MarkerFaceAlpha',.5,'HandleVisibility','off');
    % Mean + SEM
    plot(xpos(gi),mean(v),'_','Color',clrs{gi}*.6,'MarkerSize',18,'LineWidth',3,'HandleVisibility','off');
    errorbar(xpos(gi),mean(v),std(v)/sqrt(numel(v)),'Color',clrs{gi}*.6,...
        'LineWidth',2,'CapSize',6,'HandleVisibility','off');
end
set(gca,'XTick',xpos,'XTickLabel',xlbls,'YDir','reverse');
% Stats: paired t LE vs RN within each block type
ok_d=~isnan(frn_d_le)&~isnan(frn_d_rn); [~,p_d_frn]=ttest(frn_d_le(ok_d),frn_d_rn(ok_d));
ok_p=~isnan(frn_p_le)&~isnan(frn_p_rn); [~,p_p_frn]=ttest(frn_p_le(ok_p),frn_p_rn(ok_p));
yl=ylim;
text(1.5,yl(1)+.05*range(ylim),pstar(p_d_frn),'HorizontalAlignment','center','FontSize',10,'Color',CLR_D);
text(4,yl(1)+.05*range(ylim),pstar(p_p_frn),'HorizontalAlignment','center','FontSize',10,'Color',CLR_P);
% Block type labels
text(1.5,yl(2)-.02*range(ylim),'D blocks','HorizontalAlignment','center','FontSize',9,'Color',CLR_D,'FontWeight','bold');
text(4,yl(2)-.02*range(ylim),'P blocks','HorizontalAlignment','center','FontSize',9,'Color',CLR_P,'FontWeight','bold');
ylabel('Prefrontal negativity (FRN)');
title('Volatility: LE \rightarrow RN','FontSize',11);

% RIGHT: Paired strip plot — true vs false FB (stochasticity)
% false_fb derived from ~trueFB (always accurate in this dataset)
subplot(1,2,2); hold on
if ismember('false_fb',gt.Properties.VariableNames)
    gpi=gt(gt.block_type=='P'&gt.correct==0&~isnan(gt.(feat))&~isnan(gt.false_fb),:);
    mt=nan(N_subj,1);mf=nan(N_subj,1);
    for si=1:N_subj
        sm=gpi.subj_id==subj_list(si);
        mt(si)=mean(gpi.(feat)(sm&gpi.false_fb==0),'omitnan');
        mf(si)=mean(gpi.(feat)(sm&gpi.false_fb==1),'omitnan');
    end
    ok_fb=~isnan(mt)&~isnan(mf);
    % Generate jitter FIRST so lines connect to actual dot positions
    n_ok=sum(ok_fb);
    jit_t=ones(n_ok,1)+.06*randn(n_ok,1);
    jit_f=2*ones(n_ok,1)+.06*randn(n_ok,1);
    mt_ok=mt(ok_fb); mf_ok=mf(ok_fb);
    % Paired lines connecting jittered positions
    for i=1:n_ok
        plot([jit_t(i) jit_f(i)],[mt_ok(i) mf_ok(i)],'-','Color',[.6 .6 .6 .4],'LineWidth',.7,'HandleVisibility','off');
    end
    % Scatter at jittered positions
    scatter(jit_t,mt_ok,40,CLR_T,'filled','MarkerFaceAlpha',.6,'HandleVisibility','off');
    scatter(jit_f,mf_ok,40,CLR_F,'filled','MarkerFaceAlpha',.6,'HandleVisibility','off');
    % Grand mean diamonds
    plot(1,mean(mt(ok_fb)),'d','Color',CLR_T,'MarkerFaceColor',CLR_T,'MarkerSize',14,'LineWidth',1.5,'DisplayName','True FB');
    plot(2,mean(mf(ok_fb)),'d','Color',CLR_F,'MarkerFaceColor',CLR_F,'MarkerSize',14,'LineWidth',1.5,'DisplayName','False FB');
    errorbar(1,mean(mt(ok_fb)),std(mt(ok_fb))/sqrt(sum(ok_fb)),'Color',CLR_T,'LineWidth',2,'CapSize',8,'HandleVisibility','off');
    errorbar(2,mean(mf(ok_fb)),std(mf(ok_fb))/sqrt(sum(ok_fb)),'Color',CLR_F,'LineWidth',2,'CapSize',8,'HandleVisibility','off');
    % Stats
    [~,p_fb]=ttest(mt(ok_fb),mf(ok_fb));
    d_fb=(mean(mt(ok_fb))-mean(mf(ok_fb)))/std(mt(ok_fb)-mf(ok_fb));
    text(.5,.05,sprintf('paired t: %s, d=%.2f',pstar(p_fb),d_fb),'Units','normalized',...
        'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
    set(gca,'XTick',[1 2],'XTickLabel',{'True FB','False FB'},'YDir','reverse');
    xlim([.4 2.6]);
    legend('Box','off','Location','northwest','FontSize',9);
else
    text(.5,.5,{'false\_fb column not found','Check that trueFB exists in group\_table'},...
        'HorizontalAlignment','center','Units','normalized','FontSize',10,'Color','r');
end
ylabel('Prefrontal negativity (FRN)');
title('Stochasticity: true vs false FB','FontSize',11);
style_save(fig3b,fullfile(figdir,'Fig3B_FRN'));

%% ═══════════════════════════════════════════════════════════════════════════
%  FIG 3C: P300 — GROUPED BAR + JITTER + BINNED RIBBON
%  ═══════════════════════════════════════════════════════════════════════════
fig3c=figure('Position',[30 30 1100 480],'Color','w');

feat='P300_norm';
% LEFT: Grouped bar with individual dots
subplot(1,2,1); hold on
p3_d_le=nan(N_subj,1);p3_d_rn=nan(N_subj,1);
p3_p_le=nan(N_subj,1);p3_p_rn=nan(N_subj,1);
for si=1:N_subj
    p3_d_le(si)=mean(gt.(feat)(gt.subj_id==subj_list(si)&gt.block_type=='D'&gt.stage=='LE'),'omitnan');
    p3_d_rn(si)=mean(gt.(feat)(gt.subj_id==subj_list(si)&gt.block_type=='D'&gt.stage=='RN'),'omitnan');
    p3_p_le(si)=mean(gt.(feat)(gt.subj_id==subj_list(si)&gt.block_type=='P'&gt.stage=='LE'),'omitnan');
    p3_p_rn(si)=mean(gt.(feat)(gt.subj_id==subj_list(si)&gt.block_type=='P'&gt.stage=='RN'),'omitnan');
end
% Bar means
bar_data=[mean(p3_d_le,'omitnan') mean(p3_d_rn,'omitnan');
          mean(p3_p_le,'omitnan') mean(p3_p_rn,'omitnan')];
hb=bar([1 2],bar_data,.7,'grouped');
hb(1).FaceColor=CLR_D_LIGHT; hb(1).EdgeColor='none'; hb(1).FaceAlpha=.5;
hb(2).FaceColor=CLR_P_LIGHT; hb(2).EdgeColor='none'; hb(2).FaceAlpha=.5;
% Overlay individual dots
xoff_d=hb(1).XEndPoints; xoff_p=hb(2).XEndPoints;
scatter(xoff_d(1)+.04*randn(sum(~isnan(p3_d_le)),1),p3_d_le(~isnan(p3_d_le)),18,CLR_D,'filled','MarkerFaceAlpha',.4,'HandleVisibility','off');
scatter(xoff_p(1)+.04*randn(sum(~isnan(p3_p_le)),1),p3_p_le(~isnan(p3_p_le)),18,CLR_P,'filled','MarkerFaceAlpha',.4,'HandleVisibility','off');
scatter(xoff_d(2)+.04*randn(sum(~isnan(p3_d_rn)),1),p3_d_rn(~isnan(p3_d_rn)),18,CLR_D,'filled','MarkerFaceAlpha',.4,'HandleVisibility','off');
scatter(xoff_p(2)+.04*randn(sum(~isnan(p3_p_rn)),1),p3_p_rn(~isnan(p3_p_rn)),18,CLR_P,'filled','MarkerFaceAlpha',.4,'HandleVisibility','off');
% Error bars
se_d_le=std(p3_d_le,'omitnan')/sqrt(sum(~isnan(p3_d_le)));
se_d_rn=std(p3_d_rn,'omitnan')/sqrt(sum(~isnan(p3_d_rn)));
se_p_le=std(p3_p_le,'omitnan')/sqrt(sum(~isnan(p3_p_le)));
se_p_rn=std(p3_p_rn,'omitnan')/sqrt(sum(~isnan(p3_p_rn)));
errorbar(xoff_d,[bar_data(1,1) bar_data(2,1)],[se_d_le se_p_le],'k.','LineWidth',1.5,'CapSize',5,'HandleVisibility','off');
errorbar(xoff_p,[bar_data(1,2) bar_data(2,2)],[se_d_rn se_p_rn],'k.','LineWidth',1.5,'CapSize',5,'HandleVisibility','off');
set(gca,'XTick',[1 2],'XTickLabel',{'LE (pre-rev)','RN (post-rev)'});
% Stats: LE vs RN paired t within D
ok_d=~isnan(p3_d_le)&~isnan(p3_d_rn); [~,p_d_p3]=ttest(p3_d_le(ok_d),p3_d_rn(ok_d));
ok_p=~isnan(p3_p_le)&~isnan(p3_p_rn); [~,p_p_p3]=ttest(p3_p_le(ok_p),p3_p_rn(ok_p));
yl=ylim;
text(xoff_d(1),.5*(bar_data(1,1)+bar_data(2,1))+.08*range(ylim),pstar(p_d_p3),...
    'HorizontalAlignment','center','FontSize',10,'Color',CLR_D);
text(xoff_p(1),.5*(bar_data(1,2)+bar_data(2,2))+.08*range(ylim),pstar(p_p_p3),...
    'HorizontalAlignment','center','FontSize',10,'Color',CLR_P);
legend(hb,{'D','P'},'Box','off','Location','northwest');
ylabel('P300 amplitude [norm]');
title('Volatility: change-point P300','FontSize',11);

% RIGHT: P300 ~ surprise — binned scatter with CI ribbon
subplot(1,2,2); hold on
if ismember('surprise_z',gt.Properties.VariableNames)&&ismember([feat '_z'],gt.Properties.VariableNames)
    gs=gt(~isnan(gt.([feat '_z']))&~isnan(gt.surprise_z),:);
    [r_d,p_d]=plot_bins_ribbon(gs,'surprise_z',[feat '_z'],'block_type',{'D','P'},{CLR_D,CLR_P},{'D','P'});
    % Print correlation stats
    yl=ylim;
    text(.05,.95,sprintf('D: r=%.2f %s\nP: r=%.2f %s',r_d(1),pstar(p_d(1)),r_d(2),pstar(p_d(2))),...
        'Units','normalized','VerticalAlignment','top','FontSize',9,'BackgroundColor',[1 1 1 .8]);
end
xlabel('Surprise (\omega\times|\delta|) [z]'); ylabel('P300 [z]');
title('P300 ~ model surprise','FontSize',11);
legend('Box','off','Location','southwest','FontSize',9);
style_save(fig3c,fullfile(figdir,'Fig3C_P300'));

%% ═══════════════════════════════════════════════════════════════════════════
%  FIG 3D: THETA — VIOLIN + BOX (volatility) + BINNED SCATTER (PE)
%  ═══════════════════════════════════════════════════════════════════════════
fig3d=figure('Position',[30 30 1100 480],'Color','w');

feat='Theta_amp';
% LEFT: Violin-style with box overlay
subplot(1,2,1); hold on
th_d_le=nan(N_subj,1);th_d_rn=nan(N_subj,1);
th_p_le=nan(N_subj,1);th_p_rn=nan(N_subj,1);
for si=1:N_subj
    th_d_le(si)=mean(gt.(feat)(gt.subj_id==subj_list(si)&gt.block_type=='D'&gt.stage=='LE'),'omitnan');
    th_d_rn(si)=mean(gt.(feat)(gt.subj_id==subj_list(si)&gt.block_type=='D'&gt.stage=='RN'),'omitnan');
    th_p_le(si)=mean(gt.(feat)(gt.subj_id==subj_list(si)&gt.block_type=='P'&gt.stage=='LE'),'omitnan');
    th_p_rn(si)=mean(gt.(feat)(gt.subj_id==subj_list(si)&gt.block_type=='P'&gt.stage=='RN'),'omitnan');
end
grps_th={th_d_le,th_d_rn,th_p_le,th_p_rn};
clrs_th={CLR_D_LIGHT,CLR_D,CLR_P_LIGHT,CLR_P};
xpos_th=[1 2 3.5 4.5]; xlbls_th={'LE','RN','LE','RN'};
for gi=1:4
    v=grps_th{gi}; v=v(~isnan(v));
    if numel(v)<3,continue;end
    % Symmetric violin
    [kde_f,kde_xi]=ksdensity(v,'NumPoints',50);
    kde_f=kde_f/max(kde_f)*.3;
    fill([xpos_th(gi)-kde_f xpos_th(gi)+fliplr(kde_f)],[kde_xi fliplr(kde_xi)],...
        clrs_th{gi},'FaceAlpha',.35,'EdgeColor',clrs_th{gi},'EdgeAlpha',.5,'HandleVisibility','off');
    % Box: IQR
    q=quantile(v,[.25 .5 .75]);
    rectangle('Position',[xpos_th(gi)-.12 q(1) .24 q(3)-q(1)],...
        'EdgeColor',clrs_th{gi}*.6,'LineWidth',1.5,'Curvature',.1);
    plot([xpos_th(gi)-.12 xpos_th(gi)+.12],[q(2) q(2)],...
        'Color',clrs_th{gi}*.4,'LineWidth',2.5,'HandleVisibility','off');
    % Mean diamond
    plot(xpos_th(gi),mean(v),'d','Color',clrs_th{gi}*.5,'MarkerFaceColor',clrs_th{gi},...
        'MarkerSize',10,'HandleVisibility','off');
end
set(gca,'XTick',xpos_th,'XTickLabel',xlbls_th);
% Stats
ok_d=~isnan(th_d_le)&~isnan(th_d_rn); [~,p_d_th]=ttest(th_d_le(ok_d),th_d_rn(ok_d));
ok_p=~isnan(th_p_le)&~isnan(th_p_rn); [~,p_p_th]=ttest(th_p_le(ok_p),th_p_rn(ok_p));
yl=ylim;
text(1.5,yl(2)-.08*range(ylim),sprintf('D: %s',pstar(p_d_th)),'HorizontalAlignment','center','FontSize',9,'Color',CLR_D);
text(4,yl(2)-.08*range(ylim),sprintf('P: %s',pstar(p_p_th)),'HorizontalAlignment','center','FontSize',9,'Color',CLR_P);
text(1.5,yl(1)+.02*range(ylim),'D blocks','HorizontalAlignment','center','FontSize',9,'Color',CLR_D,'FontWeight','bold');
text(4,yl(1)+.02*range(ylim),'P blocks','HorizontalAlignment','center','FontSize',9,'Color',CLR_P,'FontWeight','bold');
ylabel('Frontal theta power');
title('Volatility: LE \rightarrow RN','FontSize',11);

% RIGHT: Theta ~ |PE| — binned scatter with ribbon
subplot(1,2,2); hold on
if ismember('PE_unsigned_z',gt.Properties.VariableNames)&&ismember([feat '_z'],gt.Properties.VariableNames)
    gs=gt(~isnan(gt.([feat '_z']))&~isnan(gt.PE_unsigned_z),:);
    [r_th,p_th]=plot_bins_ribbon(gs,'PE_unsigned_z',[feat '_z'],'block_type',{'D','P'},{CLR_D,CLR_P},{'D','P'});
    text(.05,.95,sprintf('D: r=%.2f %s\nP: r=%.2f %s',r_th(1),pstar(p_th(1)),r_th(2),pstar(p_th(2))),...
        'Units','normalized','VerticalAlignment','top','FontSize',9,'BackgroundColor',[1 1 1 .8]);
end
xlabel('|PE| [z]'); ylabel('Theta [z]');
title('Theta ~ |PE| magnitude','FontSize',11);
legend('Box','off','Location','northwest','FontSize',9);
style_save(fig3d,fullfile(figdir,'Fig3D_theta'));

%% ═══════════════════════════════════════════════════════════════════════════
%  FIG 3E: PLV — PATHWAY COMPARISON
%  ═══════════════════════════════════════════════════════════════════════════
fig3e=figure('Position',[30 30 1200 500],'Color','w');

% LEFT: Dot-CI — pathway × stage × block type
subplot(1,2,1); hold on
plv_fp_d_le=nan(N_subj,1);plv_fp_d_rn=nan(N_subj,1);
plv_fp_p_le=nan(N_subj,1);plv_fp_p_rn=nan(N_subj,1);
plv_fs_d_le=nan(N_subj,1);plv_fs_d_rn=nan(N_subj,1);
plv_fs_p_le=nan(N_subj,1);plv_fs_p_rn=nan(N_subj,1);
for si=1:N_subj
    s=subj_list(si);
    plv_fp_d_le(si)=mean(gt.PLV_fp(gt.subj_id==s&gt.block_type=='D'&gt.stage=='LE'),'omitnan');
    plv_fp_d_rn(si)=mean(gt.PLV_fp(gt.subj_id==s&gt.block_type=='D'&gt.stage=='RN'),'omitnan');
    plv_fp_p_le(si)=mean(gt.PLV_fp(gt.subj_id==s&gt.block_type=='P'&gt.stage=='LE'),'omitnan');
    plv_fp_p_rn(si)=mean(gt.PLV_fp(gt.subj_id==s&gt.block_type=='P'&gt.stage=='RN'),'omitnan');
    plv_fs_d_le(si)=mean(gt.PLV_fs(gt.subj_id==s&gt.block_type=='D'&gt.stage=='LE'),'omitnan');
    plv_fs_d_rn(si)=mean(gt.PLV_fs(gt.subj_id==s&gt.block_type=='D'&gt.stage=='RN'),'omitnan');
    plv_fs_p_le(si)=mean(gt.PLV_fs(gt.subj_id==s&gt.block_type=='P'&gt.stage=='LE'),'omitnan');
    plv_fs_p_rn(si)=mean(gt.PLV_fs(gt.subj_id==s&gt.block_type=='P'&gt.stage=='RN'),'omitnan');
end
xp=[1 2 4 5];
fp_means=[mean(plv_fp_d_le,'omitnan') mean(plv_fp_d_rn,'omitnan') mean(plv_fp_p_le,'omitnan') mean(plv_fp_p_rn,'omitnan')];
fp_se=[std(plv_fp_d_le,'omitnan') std(plv_fp_d_rn,'omitnan') std(plv_fp_p_le,'omitnan') std(plv_fp_p_rn,'omitnan')]./sqrt(N_subj);
fs_means=[mean(plv_fs_d_le,'omitnan') mean(plv_fs_d_rn,'omitnan') mean(plv_fs_p_le,'omitnan') mean(plv_fs_p_rn,'omitnan')];
fs_se=[std(plv_fs_d_le,'omitnan') std(plv_fs_d_rn,'omitnan') std(plv_fs_p_le,'omitnan') std(plv_fs_p_rn,'omitnan')]./sqrt(N_subj);

errorbar(xp-.1,fp_means,fp_se,'o-','Color',CLR_FP,'LineWidth',2.2,...
    'MarkerFaceColor',CLR_FP,'MarkerSize',10,'CapSize',5,'DisplayName','Fronto-parietal');
errorbar(xp+.1,fs_means,fs_se,'s-','Color',CLR_FS,'LineWidth',2.2,...
    'MarkerFaceColor',CLR_FS,'MarkerSize',10,'CapSize',5,'DisplayName','Fronto-somatosensory');
% Shaded reversal columns
yl=ylim;
fill([1.5 2.5 2.5 1.5],[yl(1) yl(1) yl(2) yl(2)],[.9 .9 .9],...
    'FaceAlpha',.15,'EdgeColor','none','HandleVisibility','off');
fill([4.5 5.5 5.5 4.5],[yl(1) yl(1) yl(2) yl(2)],[.9 .9 .9],...
    'FaceAlpha',.15,'EdgeColor','none','HandleVisibility','off');
set(gca,'XTick',xp,'XTickLabel',{'LE','RN','LE','RN'});
% Stats: LE vs RN for each pathway in D
ok_fp=~isnan(plv_fp_d_le)&~isnan(plv_fp_d_rn);
[~,p_fp_d]=ttest(plv_fp_d_le(ok_fp),plv_fp_d_rn(ok_fp));
ok_fs=~isnan(plv_fs_d_le)&~isnan(plv_fs_d_rn);
[~,p_fs_d]=ttest(plv_fs_d_le(ok_fs),plv_fs_d_rn(ok_fs));
text(1.5,yl(2)-.06*range(ylim),sprintf('FP: %s',pstar(p_fp_d)),'HorizontalAlignment','center','FontSize',8,'Color',CLR_FP);
text(1.5,yl(2)-.12*range(ylim),sprintf('FS: %s',pstar(p_fs_d)),'HorizontalAlignment','center','FontSize',8,'Color',CLR_FS);
text(1.5,yl(1)+.02*range(ylim),'D blocks','HorizontalAlignment','center','FontSize',9,'Color',CLR_D,'FontWeight','bold');
text(4.5,yl(1)+.02*range(ylim),'P blocks','HorizontalAlignment','center','FontSize',9,'Color',CLR_P,'FontWeight','bold');
ylabel('Phase-locking value (PLV)');
title('Pathway \times stage \times block type','FontSize',11);
legend('Box','off','Location','northeast','FontSize',9);
xlim([0 6]);

% RIGHT: PLV ~ surprise — dual pathway overlay with ribbons
subplot(1,2,2); hold on
if all(ismember({'PLV_fp_z','PLV_fs_z','surprise_z'},gt.Properties.VariableNames))
    gs_fp=gt(~isnan(gt.PLV_fp_z)&~isnan(gt.surprise_z),:);
    gs_fs=gt(~isnan(gt.PLV_fs_z)&~isnan(gt.surprise_z),:);
    [r_fp,p_fp]=plot_bins_ribbon_single(gs_fp,'surprise_z','PLV_fp_z',CLR_FP,'Fronto-parietal');
    [r_fs,p_fs]=plot_bins_ribbon_single(gs_fs,'surprise_z','PLV_fs_z',CLR_FS,'Fronto-somatosensory');
    text(.05,.95,sprintf('FP: r=%.2f %s\nFS: r=%.2f %s',r_fp,pstar(p_fp),r_fs,pstar(p_fs)),...
        'Units','normalized','VerticalAlignment','top','FontSize',9,'BackgroundColor',[1 1 1 .8]);
end
xlabel('Surprise [z]'); ylabel('PLV [z]');
title('PLV ~ surprise by pathway','FontSize',11);
legend('Box','off','Location','southwest','FontSize',9);
style_save(fig3e,fullfile(figdir,'Fig3E_PLV_pathways'));

%% ═══════════════════════════════════════════════════════════════════════════
%  FIG 3F: NEURAL SIGNATURES × n_prev_P — slope graph with individual traces
%  ═══════════════════════════════════════════════════════════════════════════
fig3f=figure('Position',[30 30 1400 450],'Color','w');

npps_bins=[0 1 2]; npp_x=1:3; npp_labels={'0','1','2+'};
feats_f={'prefrontal_mean_norm','Theta_amp','PLV_fp'};
feat_labels={'FRN at reversal','Theta at reversal','FP-PLV at reversal'};
rev_y=[true false false];

for fi=1:3
    subplot(1,3,fi); hold on
    feat_f=feats_f{fi};
    if ~ismember(feat_f,gt.Properties.VariableNames),continue;end

    for bt_i=1:2
        bt=ternary_pf(bt_i==1,'D','P');
        clr=ternary_pf(bt_i==1,CLR_D,CLR_P);

        subj_traces=nan(N_subj,3);
        for ni=1:3
            for si=1:N_subj
                if ni<3, nm=gt.subj_id==subj_list(si)&gt.n_prev_P==npps_bins(ni);
                else, nm=gt.subj_id==subj_list(si)&gt.n_prev_P>=2; end
                sm=nm&gt.stage=='RN'&gt.block_type==bt&~isnan(gt.(feat_f));
                subj_traces(si,ni)=mean(gt.(feat_f)(sm),'omitnan');
            end
        end
        % Individual traces
        ok_tr=all(~isnan(subj_traces),2);
        for si=find(ok_tr)'
            off=ternary_pf(bt_i==1,-.03,.03);
            plot(npp_x+off,subj_traces(si,:),'-','Color',[clr .10],'LineWidth',.5,'HandleVisibility','off');
        end
        % Group mean + ribbon
        m_vals=mean(subj_traces,1,'omitnan');
        se_vals=std(subj_traces,0,1,'omitnan')./sqrt(sum(~isnan(subj_traces),1));
        off=ternary_pf(bt_i==1,-.06,.06);
        fill([npp_x+off fliplr(npp_x+off)],[m_vals+se_vals fliplr(m_vals-se_vals)],...
            clr,'FaceAlpha',.12,'EdgeColor','none','HandleVisibility','off');
        plot(npp_x+off,m_vals,'-','Color',clr,'LineWidth',2.5,'HandleVisibility','off');
        scatter(npp_x+off,m_vals,80,clr,'filled','MarkerFaceAlpha',.9,'DisplayName',bt);
    end
    set(gca,'XTick',npp_x,'XTickLabel',npp_labels);
    if rev_y(fi), set(gca,'YDir','reverse'); end
    xlabel('Prior P-block exposure'); ylabel(feat_labels{fi});
    title(feat_labels{fi},'FontSize',11);
    xlim([.5 3.5]);
    if fi==1, legend('Box','off','Location','southwest','FontSize',9); end
    % Stat: linear trend test (0 vs 2+ within P blocks)
    % Recompute for P
    subj_p0=nan(N_subj,1);subj_p2=nan(N_subj,1);
    for si=1:N_subj
        nm0=gt.subj_id==subj_list(si)&gt.n_prev_P==0&gt.stage=='RN'&gt.block_type=='P'&~isnan(gt.(feat_f));
        nm2=gt.subj_id==subj_list(si)&gt.n_prev_P>=2&gt.stage=='RN'&gt.block_type=='P'&~isnan(gt.(feat_f));
        subj_p0(si)=mean(gt.(feat_f)(nm0),'omitnan');
        subj_p2(si)=mean(gt.(feat_f)(nm2),'omitnan');
    end
    ok_t=~isnan(subj_p0)&~isnan(subj_p2);
    if sum(ok_t)>5
        [~,p_trend]=ttest(subj_p0(ok_t),subj_p2(ok_t));
        yl=ylim;
        text(.5,.05,sprintf('P 0vs2+: %s',pstar(p_trend)),'Units','normalized',...
            'HorizontalAlignment','center','FontSize',9,'Color',CLR_P);
    end
end
style_save(fig3f,fullfile(figdir,'Fig3F_neural_nPrevP'));

%% ═══════════════════════════════════════════════════════════════════════════
%  FIG 3G: CONFIDENCE × PREFRONTAL NEGATIVITY (from S7d RQ2)
%  Shows: confidence predicts FRN magnitude (metacognitive evaluation signal)
%  Left: Confidence binned scatter → FRN (D vs P, incorrect trials)
%  Right: Confidence × FRN split by stage (LE vs RN) — reversal disrupts
%         the confidence-FRN coupling
%  ═══════════════════════════════════════════════════════════════════════════
fig3g=figure('Position',[30 30 1100 480],'Color','w');

% Prepare confidence z-score if not already done
if ~ismember('conf_z',gt.Properties.VariableNames)
    gt.conf_z=nan(height(gt),1);
    for si=1:N_subj
        sm=gt.subj_id==subj_list(si);
        cv=gt.confidence(sm); sd=std(cv,'omitnan');
        if sd>0, gt.conf_z(sm)=(cv-mean(cv,'omitnan'))/sd; end
    end
end

% Use prefrontal_mean_norm (same as FRN panels) or prefrontal_mean_norm
feat_conf='prefrontal_mean_norm';
if ~ismember(feat_conf,gt.Properties.VariableNames)
    feat_conf='prefrontal_mean_norm';
end
feat_conf_z=[feat_conf '_z'];

% LEFT: Confidence → FRN, incorrect trials, D vs P
subplot(1,2,1); hold on
if ismember('conf_z',gt.Properties.VariableNames) && ismember(feat_conf,gt.Properties.VariableNames)
    % Incorrect trials only (where FRN is most meaningful)
    gs_inc=gt(gt.correct==0 & ~isnan(gt.conf_z) & ~isnan(gt.(feat_conf)),:);
    if ismember(feat_conf_z,gs_inc.Properties.VariableNames)
        y_col=feat_conf_z;
    else
        y_col=feat_conf;
    end
    [r_c,p_c]=plot_bins_ribbon(gs_inc,'conf_z',y_col,'block_type',{'D','P'},{CLR_D,CLR_P},{'D','P'});
    text(.05,.95,sprintf('D: r=%.2f %s\nP: r=%.2f %s',r_c(1),pstar(p_c(1)),r_c(2),pstar(p_c(2))),...
        'Units','normalized','VerticalAlignment','top','FontSize',9,'BackgroundColor',[1 1 1 .8]);
    set(gca,'YDir','reverse');
end
xlabel('Confidence [z]'); ylabel('FRN [z or norm]');
title('Confidence \rightarrow FRN (incorrect trials)','FontSize',11);
legend('Box','off','Location','southwest','FontSize',9);

% RIGHT: Confidence-FRN coupling by stage (LE vs RN) — does reversal disrupt?
subplot(1,2,2); hold on
if ismember('conf_z',gt.Properties.VariableNames) && ismember(feat_conf,gt.Properties.VariableNames)
    % Compute per-subject confidence-FRN correlation at LE vs RN
    r_le=nan(N_subj,1); r_rn=nan(N_subj,1);
    for si=1:N_subj
        s_le=gt.subj_id==subj_list(si)&gt.stage=='LE'&gt.correct==0&~isnan(gt.conf_z)&~isnan(gt.(feat_conf));
        s_rn=gt.subj_id==subj_list(si)&gt.stage=='RN'&gt.correct==0&~isnan(gt.conf_z)&~isnan(gt.(feat_conf));
        if sum(s_le)>5, r_le(si)=corr(gt.conf_z(s_le),gt.(feat_conf)(s_le),'Rows','complete'); end
        if sum(s_rn)>5, r_rn(si)=corr(gt.conf_z(s_rn),gt.(feat_conf)(s_rn),'Rows','complete'); end
    end
    ok_r=~isnan(r_le)&~isnan(r_rn);
    if sum(ok_r)>5
        % Generate jitter first, then connect
        n_ok=sum(ok_r); r_le_ok=r_le(ok_r); r_rn_ok=r_rn(ok_r);
        jit1=1+.06*randn(n_ok,1); jit2=2+.06*randn(n_ok,1);
        for i=1:n_ok
            plot([jit1(i) jit2(i)],[r_le_ok(i) r_rn_ok(i)],'-','Color',[.6 .6 .6 .35],'LineWidth',.7,'HandleVisibility','off');
        end
        scatter(jit1,r_le_ok,40,[.3 .7 .3],'filled','MarkerFaceAlpha',.6,'HandleVisibility','off');
        scatter(jit2,r_rn_ok,40,[.7 .3 .3],'filled','MarkerFaceAlpha',.6,'HandleVisibility','off');
        % Grand means
        plot(1,mean(r_le_ok),'d','Color',[.2 .6 .2],'MarkerFaceColor',[.2 .6 .2],'MarkerSize',14,'LineWidth',1.5,'DisplayName','LE (expert)');
        plot(2,mean(r_rn_ok),'d','Color',[.6 .2 .2],'MarkerFaceColor',[.6 .2 .2],'MarkerSize',14,'LineWidth',1.5,'DisplayName','RN (reversal)');
        errorbar(1,mean(r_le_ok),std(r_le_ok)/sqrt(n_ok),'Color',[.2 .6 .2],'LineWidth',2,'CapSize',8,'HandleVisibility','off');
        errorbar(2,mean(r_rn_ok),std(r_rn_ok)/sqrt(n_ok),'Color',[.6 .2 .2],'LineWidth',2,'CapSize',8,'HandleVisibility','off');
        [~,p_stage]=ttest(r_le_ok,r_rn_ok);
        d_stage=(mean(r_le_ok)-mean(r_rn_ok))/std(r_le_ok-r_rn_ok);
        text(.5,.05,sprintf('%s, d=%.2f',pstar(p_stage),d_stage),'Units','normalized',...
            'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
        yline(0,':k','HandleVisibility','off');
        set(gca,'XTick',[1 2],'XTickLabel',{'LE (expert)','RN (reversal)'}); xlim([.4 2.6]);
        legend('Box','off','Location','northwest','FontSize',9);
    end
end
xlabel('Task stage'); ylabel('Within-subject r(confidence, FRN)');
title('Reversal disrupts confidence-FRN coupling','FontSize',11);
style_save(fig3g,fullfile(figdir,'Fig3G_confidence_FRN'));


fprintf('\n=== ALL POSTER FIGURES COMPLETE ===\n');
fprintf('Saved to: %s\n',figdir);

%% ═══════════════════════════════════════════════════════════════════════════
%  LOCAL HELPER FUNCTIONS
%  ═══════════════════════════════════════════════════════════════════════════

function [m_vals,se_vals]=subj_means_by_stage(gt,feat,bt,subj_list)
stages={'LN','LE','RN','RE'}; N_s=numel(subj_list);
sm=nan(N_s,4);
for si=1:N_s
    for sti=1:4
        mask=gt.subj_id==subj_list(si)&gt.block_type==bt&gt.stage==stages{sti}&~isnan(gt.(feat));
        sm(si,sti)=mean(gt.(feat)(mask),'omitnan');
    end
end
m_vals=mean(sm,1,'omitnan'); se_vals=std(sm,0,1,'omitnan')./sqrt(sum(~isnan(sm),1));
end

function [sm,se_vals]=subj_vals_by_stage(gt,feat,bt,subj_list)
% Returns per-subject matrix (N_subj × 4 stages) + SEM
stages={'LN','LE','RN','RE'}; N_s=numel(subj_list);
sm=nan(N_s,4);
for si=1:N_s
    for sti=1:4
        mask=gt.subj_id==subj_list(si)&gt.block_type==bt&gt.stage==stages{sti}&~isnan(gt.(feat));
        sm(si,sti)=mean(gt.(feat)(mask),'omitnan');
    end
end
se_vals=std(sm,0,1,'omitnan')./sqrt(sum(~isnan(sm),1));
end

function [r_vals,p_vals]=plot_bins_ribbon(T,xcol,ycol,gcol,glevels,colors,labels)
% Binned scatter with SEM ribbon. Returns per-group correlation r and p.
nb=7; r_vals=nan(numel(glevels),1); p_vals=nan(numel(glevels),1);
for gi=1:numel(glevels)
    gm=T.(gcol)==glevels{gi};
    xv=T.(xcol)(gm);yv=T.(ycol)(gm);
    ok=~isnan(xv)&~isnan(yv);xv=xv(ok);yv=yv(ok);
    if numel(xv)<20,continue;end
    [r_vals(gi),p_vals(gi)]=corr(xv,yv,'Rows','complete');
    edges=unique(quantile(xv,linspace(0,1,nb+1)));
    nb_eff=numel(edges)-1;
    bx=nan(1,nb_eff);by=nan(1,nb_eff);bs=nan(1,nb_eff);
    for bi=1:nb_eff
        if bi<nb_eff,bm=xv>=edges(bi)&xv<edges(bi+1);else,bm=xv>=edges(bi)&xv<=edges(bi+1);end
        bx(bi)=mean(xv(bm),'omitnan');by(bi)=mean(yv(bm),'omitnan');
        bs(bi)=std(yv(bm),'omitnan')/sqrt(max(sum(bm),1));
    end
    fill([bx fliplr(bx)],[by+bs fliplr(by-bs)],colors{gi},...
        'FaceAlpha',.15,'EdgeColor','none','HandleVisibility','off');
    plot(bx,by,'-','Color',colors{gi},'LineWidth',2.2,'HandleVisibility','off');
    scatter(bx,by,55,colors{gi},'filled','MarkerFaceAlpha',.85,'DisplayName',labels{gi});
end
end

function [r_val,p_val]=plot_bins_ribbon_single(T,xcol,ycol,clr,lbl)
% Single-group binned scatter with ribbon. Returns correlation r and p.
nb=7;xv=T.(xcol);yv=T.(ycol);ok=~isnan(xv)&~isnan(yv);xv=xv(ok);yv=yv(ok);
r_val=NaN;p_val=NaN;
if numel(xv)<20,return;end
[r_val,p_val]=corr(xv,yv,'Rows','complete');
edges=quantile(xv,linspace(0,1,nb+1));
edges=unique(edges); nb_eff=numel(edges)-1;
bx=nan(1,nb_eff);by=nan(1,nb_eff);bs=nan(1,nb_eff);
for bi=1:nb_eff
    if bi<nb_eff,bm=xv>=edges(bi)&xv<edges(bi+1);else,bm=xv>=edges(bi)&xv<=edges(bi+1);end
    bx(bi)=mean(xv(bm),'omitnan');by(bi)=mean(yv(bm),'omitnan');
    bs(bi)=std(yv(bm),'omitnan')/sqrt(max(sum(bm),1));
end
fill([bx fliplr(bx)],[by+bs fliplr(by-bs)],clr,...
    'FaceAlpha',.15,'EdgeColor','none','HandleVisibility','off');
plot(bx,by,'-','Color',clr,'LineWidth',2.2,'HandleVisibility','off');
scatter(bx,by,55,clr,'filled','MarkerFaceAlpha',.85,'DisplayName',lbl);
end

function style_save(fig,path_no_ext)
% Unified styling + export (PDF vector + PNG 300dpi)
set(fig,'Color','w');
ax_all=findall(fig,'Type','axes');
for k=1:numel(ax_all)
    set(ax_all(k),'TickDir','out','Box','off','FontName','Arial','FontSize',11,...
        'LineWidth',1,'TickLength',[.012 .012]);
    set(ax_all(k).Title,'FontWeight','bold');
end
if exist('apply_fig_style','file'),apply_fig_style(fig);end
[d,~,~]=fileparts(path_no_ext);if ~isempty(d)&&~exist(d,'dir'),mkdir(d);end
exportgraphics(fig,[path_no_ext '.pdf'],'ContentType','vector');
exportgraphics(fig,[path_no_ext '.png'],'Resolution',300);
fprintf('  Saved: %s\n',path_no_ext);
end

function s=pstar(p)
% Return significance string with stars
if p<.001, s=sprintf('p<.001 ***');
elseif p<.01, s=sprintf('p=%.3f **',p);
elseif p<.05, s=sprintf('p=%.3f *',p);
else, s=sprintf('p=%.3f n.s.',p);
end
end

function out=ternary_pf(cond,a,b)
if cond,out=a;else,out=b;end
end

function T = add_transition_history_columns_local(T)
% Adds n_prev_P column: number of P blocks before the current block per subject.
if ~ismember('subj_id',T.Properties.VariableNames),return;end
bc='block';if ~ismember(bc,T.Properties.VariableNames),bc='block_number';end
if ~ismember(bc,T.Properties.VariableNames),T.n_prev_P=nan(height(T),1);return;end
T.n_prev_P=nan(height(T),1);
bt_s=string(T.block_type); bt_s(bt_s=="V")="P";
subs=unique(string(T.subj_id));
for si=1:numel(subs)
    sm=string(T.subj_id)==subs(si);
    blks=unique(T.(bc)(sm)); blks=sort(blks(:))';
    prev_p=0;
    for bi=1:numel(blks)
        bm=sm&T.(bc)==blks(bi);
        T.n_prev_P(bm)=prev_p;
        curr_bt=bt_s(find(bm,1));
        if curr_bt=="P",prev_p=prev_p+1;end
    end
end
end
