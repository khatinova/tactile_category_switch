% S8_model_based_eeg.m — Model-based single-trial EEG analysis
% No stage in fixed effects (model latents capture it). false_fb as numeric.
close all; clc;
addpath(genpath(fileparts(mfilename('fullpath'))));
putils = fullfile(fileparts(fileparts(mfilename('fullpath'))),'pipeline','utils');
if exist(putils,'dir'), addpath(putils); end
base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch';
saved_tables_folder = fullfile(base_path,'Salient mod switch KH','Results','EEG analysis','Outcome_feature_tables_v4_merged');
figure_output_folder = fullfile(base_path,'Salient mod switch KH','Results','EEG analysis','Figures','S8_model_based_EEG');
if ~exist(figure_output_folder,'dir'), mkdir(figure_output_folder); end
if exist('gt','var')&&istable(gt), fprintf('Using gt.\n');
elseif exist('group_table','var')&&istable(group_table), gt=group_table;
else, load(fullfile(saved_tables_folder,'group_feature_table_combined.mat'),'group_table'); gt=group_table; end
CLR_D=[.15 .45 .70]; CLR_P=[.80 .30 .10]; CLR_TRUE=[.20 .60 .30]; CLR_FALSE=[.75 .20 .55];
fprintf('\n=== S8: MODEL-BASED EEG ===\nOutput: %s\n\n',figure_output_folder);
gt.subj_id=categorical(gt.subj_id); gt.block_type=categorical(gt.block_type);
if ismember('false_fb',gt.Properties.VariableNames), gt.false_fb=double(gt.false_fb); end
subj_list=unique(gt.subj_id);

%% MERGE NASSAR LATENTS
nassar_needed={'PE_nassar','omega','surprise'};
has_n=all(ismember(nassar_needed,gt.Properties.VariableNames))&&sum(~isnan(gt.PE_nassar))>10;
if ~has_n
    fprintf('Merging Nassar latents...\n');
    nrc={fullfile(base_path,'Salient mod switch KH','Results','Simulation results','Figures','nassar_results.mat'),...
         fullfile(base_path,'Salient mod switch KH','Data','nassar_results.mat')};
    nr_path=''; for ci=1:numel(nrc), if exist(nrc{ci},'file'), nr_path=nrc{ci}; break; end, end
    if ~isempty(nr_path)
        tmp=load(nr_path,'results'); results=tmp.results; subjs_r=fieldnames(results);
        cols_i={'PE_nassar','PE_unsigned','omega','alpha_nassar','certainty','surprise','theta_nassar'};
        for c=cols_i, if ~ismember(c{1},gt.Properties.VariableNames), gt.(c{1})=nan(height(gt),1); end, end
        sc='';bc='';tc='';
        for x={'subj_id','subjID'},if ismember(x{1},gt.Properties.VariableNames),sc=x{1};break;end,end
        for x={'block','blocknum','block_number'},if ismember(x{1},gt.Properties.VariableNames),bc=x{1};break;end,end
        for x={'trial','trialnum'},if ismember(x{1},gt.Properties.VariableNames),tc=x{1};break;end,end
        fprintf('  Keys: %s %s %s\n',sc,bc,tc); nm=0;
        for si=1:numel(subjs_r)
            sn=subjs_r{si}; r=results.(sn); if ~isfield(r,'delta_trial'),continue;end
            rows=find(string(gt.(sc))==string(sn)); if isempty(rows),continue;end
            bg=double(gt.(bc)(rows)); tg=double(gt.(tc)(rows));
            for t=1:numel(r.trial_id)
                m=rows(bg==r.block_id(t)&tg==r.trial_id(t));
                if isempty(m),continue;end
                gt.PE_nassar(m(1))=r.delta_trial(t); gt.PE_unsigned(m(1))=abs(r.delta_trial(t));
                gt.omega(m(1))=r.omega_trial(t); gt.alpha_nassar(m(1))=r.alpha_trial(t);
                gt.certainty(m(1))=r.certainty_trial(t); gt.surprise(m(1))=r.surprise(t);
                gt.theta_nassar(m(1))=r.theta_trial(t); nm=nm+1;
            end
        end
        fprintf('  Merged %d trials.\n',nm);
    else, warning('nassar_results.mat not found.');
    end
end

%% Z-SCORE & NEXT CORRECT
vz={'PE_nassar','PE_unsigned','omega','alpha_nassar','certainty','surprise','theta_nassar','prefrontal_neg_peak_norm','P300_norm','Theta_amp','PLV_fp'};
for f=1:numel(vz), fn=vz{f}; fnz=[fn '_z'];
    if ~ismember(fn,gt.Properties.VariableNames),continue;end
    gt.(fnz)=nan(height(gt),1);
    for si=1:numel(subj_list), mask=gt.subj_id==subj_list(si); v=gt.(fn)(mask); sd=std(v,'omitnan');
        if sd>0, gt.(fnz)(mask)=(v-mean(v,'omitnan'))/sd; end
    end
end
gt.next_correct=nan(height(gt),1);
bc2='block'; if ~ismember(bc2,gt.Properties.VariableNames),bc2='block_number';end
for si=1:numel(subj_list), for b=1:10
    idx=find(gt.subj_id==subj_list(si)&gt.(bc2)==b); if numel(idx)<2,continue;end
    gt.next_correct(idx(1:end-1))=gt.correct(idx(2:end));
end,end
fprintf('PE_nassar non-NaN: %d | Theta_amp non-NaN: %d | P300_norm non-NaN: %d\n',...
    sum(~isnan(gt.PE_nassar)),sum(~isnan(gt.Theta_amp)),sum(~isnan(gt.P300_norm)));

%% MRQ1: P300 ~ SURPRISE
fprintf('\n=== MRQ1 ===\n');
mdl_p3_surprise=[]; mdl_p3_omega=[]; mdl_p3_learn=[];
if all(ismember({'P300_norm_z','surprise_z','omega_z'},gt.Properties.VariableNames))
    gt_m1=gt(~isnan(gt.P300_norm_z)&~isnan(gt.surprise_z),:);
    mdl_p3_surprise=safe_fitlme(gt_m1,'P300_norm_z ~ surprise_z * block_type + (1|subj_id)');
    disp(mdl_p3_surprise.Coefficients);
    gt_m1b=gt(~isnan(gt.P300_norm_z)&~isnan(gt.omega_z),:);
    mdl_p3_omega=safe_fitlme(gt_m1b,'P300_norm_z ~ omega_z * block_type + (1|subj_id)');
    disp(mdl_p3_omega.Coefficients);
    if ismember('false_fb',gt.Properties.VariableNames)&&ismember('alpha_nassar_z',gt.Properties.VariableNames)
        gt_p=gt(gt.block_type=='P'&~isnan(gt.P300_norm_z)&~isnan(gt.alpha_nassar_z)&~isnan(gt.false_fb),:);
        if height(gt_p)>50
            mdl_p3_learn=safe_fitlme(gt_p,'alpha_nassar_z ~ P300_norm_z * false_fb + (1|subj_id)');
            disp(mdl_p3_learn.Coefficients);
        end
    end
    fig1=figure('Position',[50 50 1200 420],'Color','w');
    sgtitle('MRQ1: P300 ~ surprise','FontSize',13,'FontWeight','bold');
    subplot(1,3,1);hold on;plot_bins(gt_m1,'surprise_z','P300_norm_z','block_type',{'D','P'},{CLR_D,CLR_P},{'D','P'});
    xlabel('Surprise [z]');ylabel('P300 [z]');title('A. P300~Surprise');legend('Box','off','Location','nw');
    subplot(1,3,2);hold on;plot_bins(gt_m1b,'omega_z','P300_norm_z','block_type',{'D','P'},{CLR_D,CLR_P},{'D','P'});
    xlabel('\omega [z]');ylabel('P300 [z]');title('B. P300~\omega');
    subplot(1,3,3);hold on
    if ~isempty(mdl_p3_learn),plot_bins(gt_p,'P300_norm_z','alpha_nassar_z','false_fb',{0,1},{CLR_TRUE,CLR_FALSE},{'True','False'});end
    xlabel('P300 [z]');ylabel('\alpha [z]');title('C. P300\rightarrow\alpha');legend('Box','off','Location','nw');
    style_and_save(fig1,fullfile(figure_output_folder,'MRQ1_P300'));
end

%% MRQ2: FRN ~ PE
fprintf('\n=== MRQ2 ===\n');
mdl_frn_pe=[]; mdl_frn_upe=[]; mdl_frn_false=[]; mdl_frn_discount=[];
if all(ismember({'prefrontal_neg_peak_norm_z','PE_nassar_z'},gt.Properties.VariableNames))
    gt_m2=gt(~isnan(gt.prefrontal_neg_peak_norm_z)&~isnan(gt.PE_nassar_z),:);
    mdl_frn_pe=safe_fitlme(gt_m2,'prefrontal_neg_peak_norm_z ~ PE_nassar_z * block_type + correct + (1|subj_id)');
    disp(mdl_frn_pe.Coefficients);
    if ismember('PE_unsigned_z',gt.Properties.VariableNames)
        gt_m2b=gt(~isnan(gt.prefrontal_neg_peak_norm_z)&~isnan(gt.PE_unsigned_z),:);
        mdl_frn_upe=safe_fitlme(gt_m2b,'prefrontal_neg_peak_norm_z ~ PE_unsigned_z * block_type + correct + (1|subj_id)');
        disp(mdl_frn_upe.Coefficients);
    end
    if ismember('false_fb',gt.Properties.VariableNames)
        gt_pi=gt(gt.block_type=='P'&gt.correct==0&~isnan(gt.prefrontal_neg_peak_norm_z)&~isnan(gt.false_fb),:);
        if height(gt_pi)>30
            mdl_frn_false=safe_fitlme(gt_pi,'prefrontal_neg_peak_norm_z ~ false_fb + (1|subj_id)');
            disp(mdl_frn_false.Coefficients);
        end
        gt_pa=gt(gt.block_type=='P'&~isnan(gt.prefrontal_neg_peak_norm_z),:);
        if ismember('trial',gt_pa.Properties.VariableNames)&&height(gt_pa)>50
            gt_pa.trial_z=nan(height(gt_pa),1);
            for si=1:numel(subj_list),m=gt_pa.subj_id==subj_list(si);tv=double(gt_pa.trial(m));sd=std(tv,'omitnan');
                if sd>0,gt_pa.trial_z(m)=(tv-mean(tv,'omitnan'))/sd;end,end
            mdl_frn_discount=safe_fitlme(gt_pa,'prefrontal_neg_peak_norm_z ~ trial_z * correct + (1|subj_id)');
            disp(mdl_frn_discount.Coefficients);
        end
    end
    fig2=figure('Position',[50 50 1400 420],'Color','w');
    sgtitle('MRQ2: FRN ~ PE','FontSize',13,'FontWeight','bold');
    subplot(1,4,1);hold on;plot_bins(gt_m2,'PE_nassar_z','prefrontal_neg_peak_norm_z','block_type',{'D','P'},{CLR_D,CLR_P},{'D','P'});
    xlabel('Signed PE [z]');ylabel('FRN [z]');title('A');set(gca,'YDir','reverse');
    subplot(1,4,2);hold on
    if ~isempty(mdl_frn_upe),plot_bins(gt_m2b,'PE_unsigned_z','prefrontal_neg_peak_norm_z','block_type',{'D','P'},{CLR_D,CLR_P},{'D','P'});end
    xlabel('|PE| [z]');ylabel('FRN [z]');title('B');set(gca,'YDir','reverse');
    subplot(1,4,3);hold on
    if ~isempty(mdl_frn_false)
        mt=mean(gt_pi.prefrontal_neg_peak_norm_z(gt_pi.false_fb==0),'omitnan');
        mf=mean(gt_pi.prefrontal_neg_peak_norm_z(gt_pi.false_fb==1),'omitnan');
        st=std(gt_pi.prefrontal_neg_peak_norm_z(gt_pi.false_fb==0),'omitnan')/sqrt(sum(gt_pi.false_fb==0));
        sf=std(gt_pi.prefrontal_neg_peak_norm_z(gt_pi.false_fb==1),'omitnan')/sqrt(sum(gt_pi.false_fb==1));
        bar(1,mt,'FaceColor',CLR_TRUE,'EdgeColor','none','FaceAlpha',.7);bar(2,mf,'FaceColor',CLR_FALSE,'EdgeColor','none','FaceAlpha',.7);
        errorbar([1 2],[mt mf],[st sf],'k.','LineWidth',1.5);set(gca,'XTick',[1 2],'XTickLabel',{'True','False'},'YDir','reverse');
    end
    ylabel('FRN [z]');title('C. True vs False');
    subplot(1,4,4);hold on
    if ~isempty(mdl_frn_discount)&&exist('gt_pa','var'),plot_bins(gt_pa,'trial_z','prefrontal_neg_peak_norm_z','correct',{0,1},{CLR_FALSE,CLR_TRUE},{'Inc','Cor'});end
    xlabel('Trial [z]');ylabel('FRN [z]');title('D. Discount');set(gca,'YDir','reverse');legend('Box','off');
    style_and_save(fig2,fullfile(figure_output_folder,'MRQ2_FRN'));
end

%% MRQ3: THETA ~ |PE|
fprintf('\n=== MRQ3 ===\n');
mdl_theta_pe=[]; mdl_theta_surp=[]; mdl_theta_next=[];
if all(ismember({'Theta_amp_z','PE_unsigned_z','surprise_z'},gt.Properties.VariableNames))
    gt_m3=gt(~isnan(gt.Theta_amp_z)&~isnan(gt.PE_unsigned_z),:);
    mdl_theta_pe=safe_fitlme(gt_m3,'Theta_amp_z ~ PE_unsigned_z * block_type + correct + (1|subj_id)');
    disp(mdl_theta_pe.Coefficients);
    gt_m3b=gt(~isnan(gt.Theta_amp_z)&~isnan(gt.surprise_z),:);
    mdl_theta_surp=safe_fitlme(gt_m3b,'Theta_amp_z ~ surprise_z * block_type + (1|subj_id)');
    disp(mdl_theta_surp.Coefficients);
    gt_m3c=gt(~isnan(gt.Theta_amp_z)&~isnan(gt.next_correct),:);
    if height(gt_m3c)>50
        mdl_theta_next=safe_fitglme(gt_m3c,'next_correct ~ Theta_amp_z * block_type + correct + (1|subj_id)');
        disp(mdl_theta_next.Coefficients);
    end
    fig3=figure('Position',[50 50 1200 420],'Color','w');
    sgtitle('MRQ3: Theta ~ |PE|','FontSize',13,'FontWeight','bold');
    subplot(1,3,1);hold on;plot_bins(gt_m3,'PE_unsigned_z','Theta_amp_z','block_type',{'D','P'},{CLR_D,CLR_P},{'D','P'});
    xlabel('|PE| [z]');ylabel('\theta [z]');title('A');legend('Box','off','Location','nw');
    subplot(1,3,2);hold on;plot_bins(gt_m3b,'surprise_z','Theta_amp_z','block_type',{'D','P'},{CLR_D,CLR_P},{'D','P'});
    xlabel('Surprise [z]');ylabel('\theta [z]');title('B');
    subplot(1,3,3);hold on
    if ~isempty(mdl_theta_next),plot_bins(gt_m3c,'Theta_amp_z','next_correct','block_type',{'D','P'},{CLR_D,CLR_P},{'D','P'});end
    xlabel('\theta [z]');ylabel('P(correct_{t+1})');title('C');
    style_and_save(fig3,fullfile(figure_output_folder,'MRQ3_theta'));
end

%% MRQ4: PLV ~ THETA / MEDIATION
fprintf('\n=== MRQ4 ===\n');
mdl_plv_theta=[]; mdl_plv_surp=[]; mdl_path_a=[]; mdl_path_bc=[];
if all(ismember({'PLV_fp_z','Theta_amp_z','surprise_z'},gt.Properties.VariableNames))
    gt_m4=gt(~isnan(gt.PLV_fp_z)&~isnan(gt.Theta_amp_z),:);
    mdl_plv_theta=safe_fitlme(gt_m4,'PLV_fp_z ~ Theta_amp_z * block_type + (1|subj_id)');
    disp(mdl_plv_theta.Coefficients);
    gt_m4b=gt(~isnan(gt.PLV_fp_z)&~isnan(gt.surprise_z),:);
    mdl_plv_surp=safe_fitlme(gt_m4b,'PLV_fp_z ~ surprise_z * block_type + (1|subj_id)');
    disp(mdl_plv_surp.Coefficients);
    if ismember('PE_unsigned_z',gt.Properties.VariableNames)
        gt_m4c=gt(~isnan(gt.PLV_fp_z)&~isnan(gt.PE_unsigned_z)&~isnan(gt.next_correct),:);
        if height(gt_m4c)>50
            mdl_path_a=safe_fitlme(gt_m4c,'PLV_fp_z ~ PE_unsigned_z + block_type + (1|subj_id)');
            mdl_path_bc=safe_fitglme(gt_m4c,'next_correct ~ PLV_fp_z + PE_unsigned_z + block_type + (1|subj_id)');
        end
    end
    fig4=figure('Position',[50 50 1200 420],'Color','w');
    sgtitle('MRQ4: PLV ~ theta','FontSize',13,'FontWeight','bold');
    subplot(1,3,1);hold on;plot_bins(gt_m4,'Theta_amp_z','PLV_fp_z','block_type',{'D','P'},{CLR_D,CLR_P},{'D','P'});
    xlabel('\theta [z]');ylabel('PLV [z]');title('A');legend('Box','off','Location','nw');
    subplot(1,3,2);hold on;plot_bins(gt_m4b,'surprise_z','PLV_fp_z','block_type',{'D','P'},{CLR_D,CLR_P},{'D','P'});
    xlabel('Surprise [z]');ylabel('PLV [z]');title('B');
    ax3=subplot(1,3,3);axis(ax3,'off');title(ax3,'C. Mediation');
    if ~isempty(mdl_path_a)&&~isempty(mdl_path_bc)
        ab=mdl_path_a.Coefficients.Estimate(2);ap=mdl_path_a.Coefficients.pValue(2);
        bb=mdl_path_bc.Coefficients.Estimate(2);bp=mdl_path_bc.Coefficients.pValue(2);
        cb=mdl_path_bc.Coefficients.Estimate(3);cp=mdl_path_bc.Coefficients.pValue(3);
        txt={sprintf('|PE|->PLV: b=%.3f %s',ab,pfmt(ap)),sprintf('PLV->next: OR=%.2f %s',exp(bb),pfmt(bp)),...
             sprintf('|PE|->next: OR=%.2f %s',exp(cb),pfmt(cp)),'',sprintf('Indirect=%.4f',ab*bb)};
        text(ax3,.05,.85,txt,'Units','normalized','VerticalAlignment','top','FontSize',10);
    end
    style_and_save(fig4,fullfile(figure_output_folder,'MRQ4_PLV'));
end
fprintf('\n=== S8 COMPLETE ===\n');

%% LOCAL FUNCTIONS
function mdl = safe_fitlme(T, formula)
try
    mdl = fitlme(T, formula, 'FitMethod','REML');
catch ME
    fprintf('  [safe_fitlme] Failed: %s\n  Trying without interaction...\n', ME.message);
    f2 = regexprep(formula, '\*', '+');
    try
        mdl = fitlme(T, f2, 'FitMethod','REML');
    catch ME2
        fprintf('  [safe_fitlme] Still failed: %s\n  Dropping block_type...\n', ME2.message);
        f3 = regexprep(f2, '\+\s*block_type', '');
        mdl = fitlme(T, f3, 'FitMethod','REML');
    end
end
end

function mdl = safe_fitglme(T, formula)
try
    mdl = fitglme(T, formula, 'Distribution','Binomial','Link','logit','FitMethod','Laplace');
catch ME
    fprintf('  [safe_fitglme] Failed: %s\n  Trying without interaction...\n', ME.message);
    f2 = regexprep(formula, '\*', '+');
    try
        mdl = fitglme(T, f2, 'Distribution','Binomial','Link','logit','FitMethod','Laplace');
    catch ME2
        fprintf('  [safe_fitglme] Still failed: %s\n', ME2.message);
        f3 = regexprep(f2, '\+\s*block_type', '');
        mdl = fitglme(T, f3, 'Distribution','Binomial','Link','logit','FitMethod','Laplace');
    end
end
end

function plot_bins(T, xcol, ycol, gcol, glevels, colors, labels)
nb=5;
for gi=1:numel(glevels)
    gm = T.(gcol)==glevels{gi};
    xv=T.(xcol)(gm); yv=T.(ycol)(gm);
    ok=~isnan(xv)&~isnan(yv); xv=xv(ok); yv=yv(ok);
    if numel(xv)<20, continue; end
    edges=quantile(xv,linspace(0,1,nb+1));
    bx=nan(1,nb); by=nan(1,nb); bs=nan(1,nb);
    for bi=1:nb
        if bi<nb, bm=xv>=edges(bi)&xv<edges(bi+1); else, bm=xv>=edges(bi)&xv<=edges(bi+1); end
        bx(bi)=mean(xv(bm),'omitnan'); by(bi)=mean(yv(bm),'omitnan');
        bs(bi)=std(yv(bm),'omitnan')/sqrt(max(sum(bm),1));
    end
    errorbar(bx,by,bs,'o-','Color',colors{gi},'LineWidth',1.8,...
        'MarkerFaceColor',colors{gi},'MarkerSize',7,'CapSize',4,'DisplayName',labels{gi});
end
end

function style_and_save(fig, path_no_ext)
set(fig,'Color','w');
ax_all=findall(fig,'Type','axes');
for k=1:numel(ax_all)
    set(ax_all(k),'TickDir','out','Box','off','FontName','Arial','FontSize',10,'LineWidth',1);
end
if exist('apply_fig_style','file'), apply_fig_style(fig); end
[d,~,~]=fileparts(path_no_ext); if ~isempty(d)&&~exist(d,'dir'),mkdir(d);end
exportgraphics(fig,[path_no_ext '.pdf'],'ContentType','vector');
exportgraphics(fig,[path_no_ext '.png'],'Resolution',300);
fprintf('  Saved: %s\n',path_no_ext);
end

function s = pfmt(p)
if p<.001, s='p<.001'; else, s=sprintf('p=.%03d',round(p*1000)); end
end
