%% EEG preprocessing and epoching - FIXED Feb 2026


%% ================================
% OUTCOME ERP EPOCHING
% ================================

% load behav results
load("\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\Behav results\all_trial_data_v2.mat")

% set file paths for EEG data output, input and figure output
output_folder = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\EEG analysis\Winter 2026 outcomes\Epoch analysis';
input_folder = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\EEG analysis\Winter 2026 outcomes';
figure_output_folder = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\EEG analysis\Pilot2_Jan2026\Figures';


addpath('C:\Users\khatinova\OneDrive - Nexus365\Pre_2026_Folders\Documents\MATLAB\eeglab2025.1.0')

eeglab nogui


epoch_window = [-0.5 1.5];
baseline_win = [-500 0];

valid_subjects = [3:12,14,15,17,18];

all_correct  = [];
all_incorrect = [];
all_rewp     = [];
all_p3_correct = [];
all_p3_incorrect = [];
all_theta = [];

for s = valid_subjects

    subj = sprintf('Ox%02d',s);
    fprintf('\n=== Epoching %s ===\n',subj);

    EEG = pop_loadset([subj '_ICA_pruned.set'], input_folder);

    % -------- outcome codes ----------
    if s < 9
        outcome_codes = {'10','11'};
        correct_codes = [10];
        incorrect_codes = [11];
    else
        outcome_codes = {'31','32','33','34'};
        correct_codes = [31 32];
        incorrect_codes = [33 34];
    end

    % -------- epoch once ----------
    EEG_all = pop_epoch(EEG, outcome_codes, epoch_window);
    EEG_all = pop_rmbase(EEG_all, baseline_win);

    % -------- parse event codes ----------
    nE = length(EEG_all.epoch);
    evtCodes = nan(nE,1);

    for e = 1:nE
        evt = EEG_all.epoch(e).eventtype;
        if iscell(evt)
            match = ismember(evt,outcome_codes);
            if any(match)
                evtCodes(e) = str2double(evt{find(match,1)});
            end
        elseif ischar(evt)
            evtCodes(e) = str2double(evt);
        end
    end

    isCorrect = ismember(evtCodes,correct_codes);
    isIncorrect = ismember(evtCodes,incorrect_codes);

    EEG_correct = pop_select(EEG_all,'trial',find(isCorrect));
    EEG_incorrect = pop_select(EEG_all,'trial',find(isIncorrect));

    % -------- save ----------
    pop_saveset(EEG_correct,[subj '_correct'],output_folder);
    pop_saveset(EEG_incorrect,[subj '_incorrect'],output_folder);
    pop_saveset(EEG_all,[subj '_outcome'],output_folder);


    % =====================================
    % FEEDBACK P3 (Pz 300–600 ms)
    % =====================================
    
    pz_idx = find(strcmpi({EEG_all.chanlocs.labels},'Pz'),1);
    p3_window = [300 600];
    p3_samples = find(EEG_all.times >= p3_window(1) & EEG_all.times <= p3_window(2));
    
    if ~isempty(pz_idx)
    
        p3_cor = mean(mean(EEG_correct.data(pz_idx,p3_samples,:),2,'omitnan'),3,'omitnan');
        p3_inc = mean(mean(EEG_incorrect.data(pz_idx,p3_samples,:),2,'omitnan'),3,'omitnan');
    
        all_p3_correct(end+1)   = p3_cor;
        all_p3_incorrect(end+1) = p3_inc;
    
    end
    
    
    % =====================================
    % MIDFRONTAL THETA (FCz 4–8 Hz)
    % =====================================
    
    fc_idx = find(strcmpi({EEG_all.chanlocs.labels},'FCz'),1);
    
    if ~isempty(fc_idx)
    
        [ersp,~,~,times_tf,freqs_tf] = newtimef( ...
            EEG_all.data(fc_idx,:,:), ...
            EEG_all.pnts, ...
            [EEG_all.xmin*1000 EEG_all.xmax*1000], ...
            EEG_all.srate, ...
            0, ...
            'freqs',[4 8], ...
            'nfreqs',10, ...
            'baseline',[-500 0], ...
            'plotersp','off', ...
            'plotitc','off');
    
        theta_time_idx = find(times_tf >= 200 & times_tf <= 500);
        theta_freq_idx = find(freqs_tf >= 4 & freqs_tf <= 8);
    
        theta_val = mean(mean(ersp(theta_freq_idx,theta_time_idx,:),1),2);
    
        all_theta(end+1) = mean(theta_val,'omitnan');
    
    end






    % -------- diagnostic plot ----------
    chanIdx = find(strcmpi({EEG_all.chanlocs.labels},'FCz'),1);
    cz_idx = find(strcmpi({EEG_all.chanlocs.labels},'Cz'),1);

    erp_cor = squeeze(mean(EEG_correct.data(chanIdx,:,:),3,'omitnan'));
    erp_inc = squeeze(mean(EEG_incorrect.data(chanIdx,:,:),3,'omitnan'));
    cz_erp_cor = squeeze(mean(EEG_correct.data(cz_idx,:,:),3,'omitnan'));
    cz_erp_inc = squeeze(mean(EEG_incorrect.data(cz_idx,:,:),3,'omitnan'));
    rewp = cz_erp_cor - cz_erp_inc;

    % ---- store ----
    all_correct(:,end+1)   = erp_cor;
    all_incorrect(:,end+1) = erp_inc;
    all_rewp(:,end+1)      = rewp;

    time = EEG_all.times; % save once


    figure;
    
    subplot(2,1,1)
    hold on
    plot(EEG_all.times, erp_cor,'LineWidth',1.5, 'Color','g');
    plot(EEG_all.times, erp_inc,'LineWidth',1.5, 'Color', 'r');
    % plot(EEG_all.times, rewp,'LineWidth',1.5, 'Color', 'b');
    xlim([-200 800]);
    title([subj ' RewP (Correct - Incorrect)']);
    xlabel('Time (ms)');
    ylabel('\muV');
    xline(0,'--k');
    yline(0,'--k')
    hold off
    
    
    subplot(2,1,2)
    hold on
    plot(EEG_all.times, rewp,'LineWidth',1.5, 'Color', 'b');
    xlim([-200 800]);
    title([subj ' RewP (Correct - Incorrect)']);
    xlabel('Time (ms)');
    ylabel('\muV');
    xline(0,'--k');
    yline(0,'--k')
    hold off
end

% ================================
% GROUP AVERAGES
% ================================

grand_correct   = mean(all_correct,2,'omitnan');
grand_incorrect = mean(all_incorrect,2,'omitnan');
grand_rewp      = mean(all_rewp,2,'omitnan');

sem_correct   = std(all_correct,[],2,'omitnan')   ./ sqrt(size(all_correct,2));
sem_incorrect = std(all_incorrect,[],2,'omitnan') ./ sqrt(size(all_incorrect,2));
sem_rewp      = std(all_rewp,[],2,'omitnan')      ./ sqrt(size(all_rewp,2));

% ================================
% PLOT
% ================================

figure;

% ---- Correct vs Incorrect ----
subplot(2,1,1)
hold on

% Shaded SEM for correct
fill([time fliplr(time)], ...
     [grand_correct'+sem_correct' fliplr(grand_correct'-sem_correct')], ...
     [0.6 1 0.6], 'EdgeColor','none','FaceAlpha',0.3);

% Shaded SEM for incorrect
fill([time fliplr(time)], ...
     [grand_incorrect'+sem_incorrect' fliplr(grand_incorrect'-sem_incorrect')], ...
     [1 0.6 0.6], 'EdgeColor','none','FaceAlpha',0.3);

plot(time, grand_correct,'g','LineWidth',2);
plot(time, grand_incorrect,'r','LineWidth',2);

xline(0,'--k');
yline(0,'--k');
xlim([-200 800]);
title('Grand Average Cz: Correct vs Incorrect');
xlabel('Time (ms)');
ylabel('\muV');
legend({'Correct SEM','Incorrect SEM','Correct','Incorrect'});
hold off

% ---- RewP ----
subplot(2,1,2)
hold on

fill([time fliplr(time)], ...
     [grand_rewp'+sem_rewp' fliplr(grand_rewp'-sem_rewp')], ...
     [0.7 0.7 1], 'EdgeColor','none','FaceAlpha',0.3);

plot(time, grand_rewp,'b','LineWidth',2);

xline(0,'--k');
yline(0,'--k');
xlim([-200 800]);
title('Grand Average RewP (Correct − Incorrect) at Cz');
xlabel('Time (ms)');
ylabel('\muV');
hold off


%% P3 plot
figure;

bar([mean(all_p3_correct,'omitnan'), mean(all_p3_incorrect,'omitnan')]);
set(gca,'XTickLabel',{'Correct','Incorrect'});
ylabel('\muV');
title('Grand Average P3 (Pz 300–600 ms)');

%% prefrontal theta plot
figure;

bar(mean(all_theta,'omitnan'));
ylabel('Theta Power (dB)');
title('Grand Average Midfrontal Theta (FCz 200–500 ms)');