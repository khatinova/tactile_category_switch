% S8_model_based_eeg.m - Model-based single-trial EEG analysis
% No stage in models (latents capture it). false_fb as numeric 0/1.
close all; clc;
addpath(genpath(fileparts(mfilename('fullpath'))));
putils=fullfile(fileparts(fileparts(mfilename('fullpath'))),'pipeline','utils');
if exist(putils,'dir'),addpath(putils);end
base_path='\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch';
stf=fullfile(base_path,'Salient mod switch KH','Results','EEG analysis','Outcome_feature_tables_v4_merged');
figdir=fullfile(base_path,'Salient mod switch KH','Results','EEG analysis','Figures','S8_model_based_EEG');
if ~exist(figdir,'dir'),mkdir(figdir);end
if exist('gt','var')&&istable(gt),fprintf('Using gt.\n');
elseif exist('group_table','var')&&istable(group_table),gt=group_table;
else,load(fullfile(stf,'group_feature_table_combined.mat'),'group_table');gt=group_table;end
CLR_D=[.15 .45 .70];CLR_P=[.80 .30 .10];CLR_T=[.20 .60 .30];CLR_F=[.75 .20 .55];


fprintf('\n=== S8 MODEL-BASED EEG ===\n');
gt.subj_id=categorical(gt.subj_id);gt.block_type=categorical(gt.block_type);
if ismember('false_fb',gt.Properties.VariableNames),gt.false_fb=double(~gt.trueFB);end
subj_list=unique(gt.subj_id);

gt = gt(gt.subj_id)

%% MERGE NASSAR LATENTS IF ABSENT
nn={'PE_nassar','omega','surprise'};
has_n=all(ismember(nn,gt.Properties.VariableNames))&&sum(~isnan(gt.PE_nassar))>10;
if ~has_n
fprintf('Merging Nassar latents...\n');
nrc={fullfile(base_path,'Salient mod switch KH','Results','Simulation results','Figures','nassar_results.mat'),...
     fullfile(base_path,'Salient mod switch KH','Data','nassar_results.mat')};
nr='';for ci=1:numel(nrc),if exist(nrc{ci},'file'),nr=nrc{ci};break;end,end
if ~isempty(nr)
tmp=load(nr,'results');results=tmp.results;sr=fieldnames(results);
ci2={'PE_nassar','PE_unsigned','omega','alpha_nassar','certainty','surprise','theta_nassar'};
for c=ci2,if ~ismember(c{1},gt.Properties.VariableNames),gt.(c{1})=nan(height(gt),1);end,end
sc='';bc='';tc='';

for x={'subj_id','subjID'},if ismember(x{1},gt.Properties.VariableNames),sc=x{1};break;end,end

for x={'block','blocknum','block_number'},if ismember(x{1},gt.Properties.VariableNames),bc=x{1};break;end,end

for x={'trial','trialnum'},if ismember(x{1},gt.Properties.VariableNames),tc=x{1};break;end,end
fprintf('  Keys: %s %s %s\n',sc,bc,tc);nm=0;
for si=1:numel(sr)
sn=sr{si};r=results.(sn);if ~isfield(r,'delta_trial'),continue;end
rows=find(string(gt.(sc))==string(sn));if isempty(rows),continue;end
bg=double(gt.(bc)(rows));tg=double(gt.(tc)(rows));
for t=1:numel(r.trial_id)
m=rows(bg==r.block_id(t)&tg==r.trial_id(t));if isempty(m),continue;end
gt.PE_nassar(m(1))=r.delta_trial(t);gt.PE_unsigned(m(1))=abs(r.delta_trial(t));
gt.omega(m(1))=r.omega_trial(t);gt.alpha_nassar(m(1))=r.alpha_trial(t);
gt.certainty(m(1))=r.certainty_trial(t);gt.surprise(m(1))=r.surprise(t);
gt.theta_nassar(m(1))=r.theta_trial(t);nm=nm+1;
end,end
fprintf('  Merged %d trials.\n',nm);
else,warning('nassar_results.mat not found.');end
end

%% Z-SCORE AND NEXT CORRECT
vz={'PE_nassar','PE_unsigned','omega','alpha_nassar','certainty','surprise','theta_nassar','prefrontal_neg_peak_norm','P300_norm','Theta_amp','PLV_fp'};
for f=1:numel(vz),fn=vz{f};fnz=[fn '_z'];
if ~ismember(fn,gt.Properties.VariableNames),continue;end
gt.(fnz)=nan(height(gt),1);
for si=1:numel(subj_list),mask=gt.subj_id==subj_list(si);v=gt.(fn)(mask);sd=std(v,'omitnan');
if sd>0,gt.(fnz)(mask)=(v-mean(v,'omitnan'))/sd;end,end,end
gt.next_correct=nan(height(gt),1);
bc2='block';if ~ismember(bc2,gt.Properties.VariableNames),bc2='block_number';end
for si=1:numel(subj_list),for b=1:10
idx=find(gt.subj_id==subj_list(si)&gt.(bc2)==b);if numel(idx)<2,continue;end
gt.next_correct(idx(1:end-1))=gt.correct(idx(2:end));end,end
fprintf('PE_nassar:%d Theta:%d P300:%d non-NaN\n',sum(~isnan(gt.PE_nassar)),sum(~isnan(gt.Theta_amp)),sum(~isnan(gt.P300_norm)));

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
%% MRQ1: P300 ~ SURPRISE
fprintf('\n=== MRQ1 ===\n');
mdl1a=[];mdl1b=[];mdl1c=[];
if all(ismember({'P300_norm_z','surprise_z','omega_z'},gt.Properties.VariableNames))
gt1=gt(~isnan(gt.P300_norm_z)&~isnan(gt.surprise_z),:);
mdl1a=safe_lme(gt1,'P300_norm_z ~ surprise_z * block_type + (1|subj_id)');disp(mdl1a.Coefficients);
gt1b=gt(~isnan(gt.P300_norm_z)&~isnan(gt.omega_z),:);
mdl1b=safe_lme(gt1b,'P300_norm_z ~ omega_z * block_type + (1|subj_id)');disp(mdl1b.Coefficients);
if ismember('false_fb',gt.Properties.VariableNames)&&ismember('alpha_nassar_z',gt.Properties.VariableNames)
gtp=gt(gt.block_type=='P'&~isnan(gt.P300_norm_z)&~isnan(gt.alpha_nassar_z)&~isnan(gt.false_fb),:);
if height(gtp)>50
    mdl1c=safe_lme(gtp,'alpha_nassar_z ~ P300_norm_z * false_fb + (1|subj_id)');
    disp(mdl1c.Coefficients);end
end
fig1=figure('Position',[50 50 1200 420],'Color','w');
sgtitle('MRQ1: P300 tracks surprise','FontSize',13,'FontWeight','bold');
subplot(1,3,1);hold on;plot_bins(gt1,'surprise_z','P300_norm_z','block_type',{'D','P'},{CLR_D,CLR_P},{'D','P'});
xlabel('Surprise [z]');ylabel('P300 [z]');title('A. P300~Surprise');legend('Box','off','Location','nw');
subplot(1,3,2);hold on;plot_bins(gt1b,'omega_z','P300_norm_z','block_type',{'D','P'},{CLR_D,CLR_P},{'D','P'});
xlabel('\omega [z]');ylabel('P300 [z]');title('B. P300~\omega');
subplot(1,3,3);hold on
if ~isempty(mdl1c),plot_bins(gtp,'P300_norm_z','alpha_nassar_z','false_fb',{0,1},{CLR_T,CLR_F},{'True','False'});end
xlabel('P300 [z]');ylabel('\alpha [z]');title('C. P300\rightarrow\alpha');legend('Box','off','Location','nw');
style_save(fig1,fullfile(figdir,'MRQ1_P300'));
end

%% MRQ2: FRN ~ PE
fprintf('\n=== MRQ2 ===\n');
mdl2a=[];mdl2b=[];mdl2c=[];mdl2d=[];
if all(ismember({'prefrontal_neg_peak_norm_z','PE_nassar_z'},gt.Properties.VariableNames))
gt2=gt(~isnan(gt.prefrontal_neg_peak_norm_z)&~isnan(gt.PE_nassar_z),:);
mdl2a=safe_lme(gt2,'prefrontal_neg_peak_norm_z ~ PE_nassar_z * block_type + correct + (1|subj_id)');disp(mdl2a.Coefficients);
if ismember('PE_unsigned_z',gt.Properties.VariableNames)
gt2b=gt(~isnan(gt.prefrontal_neg_peak_norm_z)&~isnan(gt.PE_unsigned_z),:);
mdl2b=safe_lme(gt2b,'prefrontal_neg_peak_norm_z ~ PE_unsigned_z * block_type + correct + (1|subj_id)');disp(mdl2b.Coefficients);end
if ismember('false_fb',gt.Properties.VariableNames)
gpi=gt(gt.block_type=='P'&gt.correct==0&~isnan(gt.prefrontal_neg_peak_norm_z)&~isnan(gt.false_fb),:);
if height(gpi)>30,mdl2c=safe_lme(gpi,'prefrontal_neg_peak_norm_z ~ false_fb + (1|subj_id)');disp(mdl2c.Coefficients);end
gpa=gt(gt.block_type=='P'&~isnan(gt.prefrontal_neg_peak_norm_z),:);
if ismember('trial',gpa.Properties.VariableNames)&&height(gpa)>50
gpa.trial_z=nan(height(gpa),1);
for si=1:numel(subj_list),m=gpa.subj_id==subj_list(si);tv=double(gpa.trial(m));sd=std(tv,'omitnan');
if sd>0,gpa.trial_z(m)=(tv-mean(tv,'omitnan'))/sd;end,end
mdl2d=safe_lme(gpa,'prefrontal_neg_peak_norm_z ~ trial_z * correct + (1|subj_id)');disp(mdl2d.Coefficients);end,end
fig2=figure('Position',[50 50 1400 420],'Color','w');
sgtitle('MRQ2: FRN ~ PE','FontSize',13,'FontWeight','bold');
subplot(1,4,1);hold on;plot_bins(gt2,'PE_nassar_z','prefrontal_neg_peak_norm_z','block_type',{'D','P'},{CLR_D,CLR_P},{'D','P'});
xlabel('Signed PE [z]');ylabel('FRN [z]');title('A');set(gca,'YDir','reverse');
subplot(1,4,2);hold on
if ~isempty(mdl2b),plot_bins(gt2b,'PE_unsigned_z','prefrontal_neg_peak_norm_z','block_type',{'D','P'},{CLR_D,CLR_P},{'D','P'});end
xlabel('|PE| [z]');ylabel('FRN [z]');title('B');set(gca,'YDir','reverse');
subplot(1,4,3);hold on
if ~isempty(mdl2c)
mt=mean(gpi.prefrontal_neg_peak_norm_z(gpi.false_fb==0),'omitnan');
mf=mean(gpi.prefrontal_neg_peak_norm_z(gpi.false_fb==1),'omitnan');
st=std(gpi.prefrontal_neg_peak_norm_z(gpi.false_fb==0),'omitnan')/sqrt(sum(gpi.false_fb==0));
sf=std(gpi.prefrontal_neg_peak_norm_z(gpi.false_fb==1),'omitnan')/sqrt(sum(gpi.false_fb==1));
bar(1,mt,'FaceColor',CLR_T,'EdgeColor','none','FaceAlpha',.7);bar(2,mf,'FaceColor',CLR_F,'EdgeColor','none','FaceAlpha',.7);
errorbar([1 2],[mt mf],[st sf],'k.','LineWidth',1.5);set(gca,'XTick',[1 2],'XTickLabel',{'True','False'},'YDir','reverse');end
ylabel('FRN [z]');title('C. True vs False');
subplot(1,4,4);hold on
if ~isempty(mdl2d)&&exist('gpa','var'),plot_bins(gpa,'trial_z','prefrontal_neg_peak_norm_z','correct',{0,1},{CLR_F,CLR_T},{'Inc','Cor'});end
xlabel('Trial [z]');ylabel('FRN [z]');title('D. Discount');set(gca,'YDir','reverse');legend('Box','off');
style_save(fig2,fullfile(figdir,'MRQ2_FRN'));
end

%% MRQ3: THETA ~ |PE|
fprintf('\n=== MRQ3 ===\n');
mdl3a=[];mdl3b=[];mdl3c=[];
if all(ismember({'Theta_amp_z','PE_unsigned_z','surprise_z'},gt.Properties.VariableNames))
gt3=gt(~isnan(gt.Theta_amp_z)&~isnan(gt.PE_unsigned_z),:);
mdl3a=safe_lme(gt3,'Theta_amp_z ~ PE_unsigned_z * block_type + correct + (1|subj_id)');disp(mdl3a.Coefficients);
gt3b=gt(~isnan(gt.Theta_amp_z)&~isnan(gt.surprise_z),:);
mdl3b=safe_lme(gt3b,'Theta_amp_z ~ surprise_z * block_type + (1|subj_id)');disp(mdl3b.Coefficients);
gt3c=gt(~isnan(gt.Theta_amp_z)&~isnan(gt.next_correct),:);
if height(gt3c)>50
mdl3c=safe_glme(gt3c,'next_correct ~ Theta_amp_z * block_type + correct + (1|subj_id)');disp(mdl3c.Coefficients);end
fig3=figure('Position',[50 50 1200 420],'Color','w');
sgtitle('MRQ3: Theta ~ |PE|','FontSize',13,'FontWeight','bold');
subplot(1,3,1);hold on;plot_bins(gt3,'PE_unsigned_z','Theta_amp_z','block_type',{'D','P'},{CLR_D,CLR_P},{'D','P'});
xlabel('|PE| [z]');ylabel('\theta [z]');title('A');legend('Box','off','Location','nw');
subplot(1,3,2);hold on;plot_bins(gt3b,'surprise_z','Theta_amp_z','block_type',{'D','P'},{CLR_D,CLR_P},{'D','P'});
xlabel('Surprise [z]');ylabel('\theta [z]');title('B');
subplot(1,3,3);hold on
if ~isempty(mdl3c),plot_bins(gt3c,'Theta_amp_z','next_correct','block_type',{'D','P'},{CLR_D,CLR_P},{'D','P'});end
xlabel('\theta [z]');ylabel('P(correct_{t+1})');title('C');
style_save(fig3,fullfile(figdir,'MRQ3_theta'));
end

%% MRQ4: PLV ~ THETA
fprintf('\n=== MRQ4 ===\n');
mdl4a=[];mdl4b=[];mdl4c=[];mdl4d=[];
if all(ismember({'PLV_fp_z','Theta_amp_z','surprise_z'},gt.Properties.VariableNames))
gt4=gt(~isnan(gt.PLV_fp_z)&~isnan(gt.Theta_amp_z),:);
mdl4a=safe_lme(gt4,'PLV_fp_z ~ Theta_amp_z * block_type + (1|subj_id)');disp(mdl4a.Coefficients);
gt4b=gt(~isnan(gt.PLV_fp_z)&~isnan(gt.surprise_z),:);
mdl4b=safe_lme(gt4b,'PLV_fp_z ~ surprise_z * block_type + (1|subj_id)');disp(mdl4b.Coefficients);
if ismember('PE_unsigned_z',gt.Properties.VariableNames)
gt4c=gt(~isnan(gt.PLV_fp_z)&~isnan(gt.PE_unsigned_z)&~isnan(gt.next_correct),:);
if height(gt4c)>50
mdl4c=safe_lme(gt4c,'PLV_fp_z ~ PE_unsigned_z + block_type + (1|subj_id)');
mdl4d=safe_glme(gt4c,'next_correct ~ PLV_fp_z + PE_unsigned_z + block_type + (1|subj_id)');
end,end
fig4=figure('Position',[50 50 1200 420],'Color','w');
sgtitle('MRQ4: PLV ~ theta','FontSize',13,'FontWeight','bold');
subplot(1,3,1);hold on;plot_bins(gt4,'Theta_amp_z','PLV_fp_z','block_type',{'D','P'},{CLR_D,CLR_P},{'D','P'});
xlabel('\theta [z]');ylabel('PLV [z]');title('A');legend('Box','off','Location','nw');
subplot(1,3,2);hold on;plot_bins(gt4b,'surprise_z','PLV_fp_z','block_type',{'D','P'},{CLR_D,CLR_P},{'D','P'});
xlabel('Surprise [z]');ylabel('PLV [z]');title('B');
ax3=subplot(1,3,3);axis(ax3,'off');title(ax3,'C. Mediation');
if ~isempty(mdl4c)&&~isempty(mdl4d)
ab=mdl4c.Coefficients.Estimate(2);ap=mdl4c.Coefficients.pValue(2);
bb=mdl4d.Coefficients.Estimate(2);bp=mdl4d.Coefficients.pValue(2);
cb=mdl4d.Coefficients.Estimate(3);cpv=mdl4d.Coefficients.pValue(3);
txt={sprintf('|PE|->PLV: b=%.3f %s',ab,pf(ap)),sprintf('PLV->next: OR=%.2f %s',exp(bb),pf(bp)),...
     sprintf('|PE|->next: OR=%.2f %s',exp(cb),pf(cpv)),'',sprintf('Indirect=%.4f',ab*bb)};
text(ax3,.05,.85,txt,'Units','normalized','VerticalAlignment','top','FontSize',10);end
style_save(fig4,fullfile(figdir,'MRQ4_PLV'));
end
fprintf('\n=== S8 COMPLETE ===\n');

%% LOCAL FUNCTIONS
function mdl=safe_lme(T,formula)
try, mdl=fitlme(T,formula,'FitMethod','REML');
catch ME, fprintf('  LME failed: %s\n  Dropping interaction...\n',ME.message);
    f2=regexprep(formula,'\*','+');
    try, mdl=fitlme(T,f2,'FitMethod','REML');
    catch, f3=regexprep(f2,'\+ block_type','');mdl=fitlme(T,f3,'FitMethod','REML');end
end,end

function mdl=safe_glme(T,formula)
try, mdl=fitglme(T,formula,'Distribution','Binomial','Link','logit','FitMethod','Laplace');
catch ME, fprintf('  GLME failed: %s\n  Dropping interaction...\n',ME.message);
    f2=regexprep(formula,'\*','+');
    try, mdl=fitglme(T,f2,'Distribution','Binomial','Link','logit','FitMethod','Laplace');
    catch, f3=regexprep(f2,'\+ block_type','');mdl=fitglme(T,f3,'Distribution','Binomial','Link','logit','FitMethod','Laplace');end
end,end

function plot_bins(T,xcol,ycol,gcol,glevels,colors,labels)
nb=5;
for gi=1:numel(glevels)
    gm=T.(gcol)==glevels{gi};
    xv=T.(xcol)(gm);yv=T.(ycol)(gm);
    ok=~isnan(xv)&~isnan(yv);xv=xv(ok);yv=yv(ok);
    if numel(xv)<20,continue;end
    edges=quantile(xv,linspace(0,1,nb+1));
    bx=nan(1,nb);by=nan(1,nb);bs=nan(1,nb);
    for bi=1:nb
        if bi<nb,bm=xv>=edges(bi)&xv<edges(bi+1);else,bm=xv>=edges(bi)&xv<=edges(bi+1);end
        bx(bi)=mean(xv(bm),'omitnan');by(bi)=mean(yv(bm),'omitnan');
        bs(bi)=std(yv(bm),'omitnan')/sqrt(max(sum(bm),1));
    end
    errorbar(bx,by,bs,'o-','Color',colors{gi},'LineWidth',1.8,...
        'MarkerFaceColor',colors{gi},'MarkerSize',7,'CapSize',4,'DisplayName',labels{gi});
end,end

function style_save(fig,path_no_ext)
set(fig,'Color','w');
ax_all=findall(fig,'Type','axes');
for k=1:numel(ax_all)
    set(ax_all(k),'TickDir','out','Box','off','FontName','Arial','FontSize',10,'LineWidth',1);
end
if exist('apply_fig_style','file'),apply_fig_style(fig);end
[d,~,~]=fileparts(path_no_ext);if ~isempty(d)&&~exist(d,'dir'),mkdir(d);end
exportgraphics(fig,[path_no_ext '.pdf'],'ContentType','vector');
exportgraphics(fig,[path_no_ext '.png'],'Resolution',300);
fprintf('  Saved: %s\n',path_no_ext);
end

function s=pf(p)
if p<.001,s='p<.001';else,s=sprintf('p=.%03d',round(p*1000));end,end
