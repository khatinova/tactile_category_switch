% clear all
% close all
% data_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Data';
% addpath(genpath(data_path));
% cd(data_path);

% % load behav results
% load("\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\Behav results\all_trial_data_v2.mat")
% 
% % set file paths for EEG data output, input and figure output
% output_folder = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\EEG analysis\Winter 2026 outcomes\Epoch analysis\Epoch_analysis_long';
% input_folder = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\EEG analysis\Winter 2026 outcomes';
% figure_output_folder = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\EEG analysis\Pilot2_Jan2026\Figures';
% 
% 
% addpath('C:\Users\khatinova\OneDrive - Nexus365\Pre_2026_Folders\Documents\MATLAB\eeglab2025.1.0')
% 
% eeglab nogui
% 
% % extract specific electrodes
% electrodeList = {'FPZ', 'AFZ', 'FZ', 'FCZ', 'CZ', 'PZ'}; % PZ used for P300, CZ used for RewP calculations
% rev_window = 30; % how many trials around the reversal do you want to look at revaligned signatures?
% valid_participants=[3:8,10:12,14:15];
% rm_baseline = [-100 0];
% epoch_window = [-0.2 0.8];

triggers = struct( ...
    'trial_start',        1, ...
    'block_start',        2, ...
    'stim_1_GO',         11, ...
    'stim_2_GO',         12, ...
    'stim_3_GO',         13, ...
    'stim_4_GO',         14, ...
    'stim_1_NOGO',       21, ...
    'stim_2_NOGO',       22, ...
    'stim_3_NOGO',       23, ...
    'stim_4_NOGO',       24, ...
    'correct_trueFB',    31, ...
    'correct_falseFB',   32, ...
    'incorrect_trueFB',  33, ...
    'incorrect_falseFB', 34, ...
    'confidence_rating', 40, ...
    'response',          50, ...
    'response_go',       51, ...
    'response_nogo',     52, ...
    'response_miss',     53, ...
    'response_FA',       54, ...
    'response_late',     55 ...
);
 
for i=valid_participants
    participant_string       = sprintf('Ox%02d', i);   
    
    data_by_participant{i}   = all_trial_data.(participant_string).trial_data;
end

%%%%%%%%%%%% COMPUTING OPTIONS %%%%%%%%%%%%%
save_outcome_epochs   = 0;
outcome_erp_analysis  = 1;
outcome_ERP_plotting  = 1; % requires you to have outcome_erp_analysis = 1 or have the variables calculated and loaded from directory
valid_participants = [3:12,14:18];
if save_outcome_epochs==1
    
        for i=valid_participants

            
            
            subject_set_file = sprintf('Ox%02d_ICA_pruned.set', i);
            

            EEG1 = pop_loadset(subject_set_file,input_folder);
            
            

            if i < 9

                [EEG2c, outcome_indices] =  pop_epoch(EEG1, {'10', '11'}, epoch_window);
                 

            elseif i > 8 

                [EEG2c, outcome_indices] =  pop_epoch(EEG1, {'31', '32', '33', '34'}, epoch_window);


            end 

            EEG2c = pop_rmbase(EEG2c, rm_baseline);


            save2c = sprintf('Ox%02d_outcome_short', i);

            pop_saveset(EEG2c, save2c, output_folder)
        end


end




if outcome_erp_analysis==1
    
    for participant = valid_participants           

         % load the behavioural data
         beh_data_one_participant=data_by_participant{participant};         
  
         num_blocks = height(beh_data_one_participant.correct);
         num_trials_per_block = width(beh_data_one_participant.correct);
         total_trials = num_blocks*num_trials_per_block; % total number of trials in the behavioural dataset

         beh_correct  = [];
         for b = 1:num_blocks
            r= beh_data_one_participant.correct(b,:);
            beh_correct = [beh_correct, r];
         end
         beh_incorrect = ~beh_correct;
         
       
         % load the eeg data

        fname = sprintf('Ox%02d_outcome_short.set', participant);
        EEG_   = pop_loadset(fname,output_folder);

        nChan    = size(EEG_.data,1);
        nSamples    = size(EEG_.data,2);
        nEEGTrials  = size(EEG_.data, 3);            % total number of EEG recorded trails


        % aligns each of the EEG trials relative to the behaviour trial number
        alignedEEG = robust_align_epochs(EEG_, beh_correct);

        
        trialNum = alignedEEG';

   
                       
        %% ---------- pre-allocate trial-aligned containers ------------------
        EEGstruct.correct   = nan(nChan, nSamples, total_trials);
        EEGstruct.incorrect = nan(nChan, nSamples, total_trials);

        %% ---------- drop every epoch into its trial slot -------------------
     for epoch = 1:length(alignedEEG)

        t = alignedEEG(epoch);   % behavioural trial number
        
        if isnan(t)
            continue   % skip unmatched EEG epochs
        end
        
        % Now classify using BEHAVIOURAL outcome (not EEG trigger)
        if beh_correct(t) == 1
            EEGstruct.correct(:,:,t) = EEG_.data(:,:,epoch);
        else
            EEGstruct.incorrect(:,:,t) = EEG_.data(:,:,epoch);
        end
        
     end

     %% ------------ Find pre and post reversal trial indices -----------------
        rev_cont_idx = beh_data_one_participant.revTrial + ...
                       (0:num_blocks-1) * num_trials_per_block;
        
        globalRevIdx(:, participant) = rev_cont_idx(:);
        
        pre_reversal_trials  = [];
        post_reversal_trials = [];
        
        for r = rev_cont_idx
        
            pre_idx  = max(1, r-rev_window) : r-1;
            post_idx = r : min(total_trials, r+rev_window);
        
            pre_reversal_trials  = [pre_reversal_trials  pre_idx];
            post_reversal_trials = [post_reversal_trials post_idx];
        end



        

        EEGstruct.prerev_correct= EEGstruct.correct(:,:,pre_reversal_trials);
        EEGstruct.postrev_correct = EEGstruct.correct(:,:,post_reversal_trials);
        EEGstruct.prerev_incorrect = EEGstruct.incorrect(:,:,pre_reversal_trials);
        EEGstruct.postrev_incorrect = EEGstruct.incorrect(:,:,post_reversal_trials);
        
        
        
        %% ---------- build map: electrode name → index in this dataset -------
        
        labelList = {EEG_.chanlocs.labels};
        
        % Clean labels (remove whitespace + force uppercase)
        labelList = strtrim(labelList);
        labelList = upper(labelList);
        
        ERPout.Chanlocs{participant} = labelList;
        
        eIdx = nan(1, numel(electrodeList));
        
        for e = 1:numel(electrodeList)
        
            targetLabel = upper(strtrim(electrodeList{e}));
        
            channel_number = find(strcmp(labelList, targetLabel), 1);
        
            if isempty(channel_number)
                warning('Electrode "%s" NOT found for participant %d.', targetLabel, participant);
                eIdx(e) = NaN;
            else
                eIdx(e) = channel_number;
            end
        end

   
        
        %% Loop over all fields to compute mean ERP (chan x time) across trials for each participant
        for f = 1:numel(fieldNames)
            field = fieldNames{f};
            
            % Extract data: chan x time x trials
            data = EEGstruct.(field);
            
            % Average across trials (3rd dimension), ignoring NaNs
            ERPallChan = squeeze(mean(data, 3, 'omitnan')); % chan x time
            
            for e = 1:numel(electrodeList)
                electrodeName = electrodeList{e};
                chanIdx = eIdx(e); % channel index of this electrode
                sig = ERPallChan(chanIdx, :); % 1 x time vector in a single channel
                
                if ~isfield(ERPout, electrodeName)
                    ERPout.(electrodeName) = struct();
                end
                if ~isfield(ERPout.(electrodeName), field)
                    ERPout.(electrodeName).(field) = [];
                end
                
                ERPout.(electrodeName).(field)(:, participant) = sig.';

            end
        end
    end
    % -------------- save everything in one MAT file -------------------------
    ERPout.globalrevIdx = globalRevIdx;
    save('ERPoutcome.mat','ERPout');
   
end


if outcome_ERP_plotting == 1    

        % -------------------------------------------------------
        %  electrodeList = {'FCZ' 'CZ' 'PZ' 'C1' 'C3' 'CP1' 'CP3'};
        %  ERPout.(electrode).(field)  :  time × participants
        % -------------------------------------------------------
        
        SEM = struct();                          % container for standard errors
        
        for e = 1:numel(electrodeList)
            elec = electrodeList{e};
            if ~isfield(ERPout, elec), continue, end
        
            allFields = fieldnames(ERPout.(elec));   % e.g. 'GOstim','pre_trig3',...
            for f = 1:numel(allFields)
                fld  = allFields{f};                 % current ERP field
                A    = ERPout.(elec).(fld);          % time × participants
        
                % SEM:  std over participants  ÷ sqrt(N valid per time point)
                N             = sum(~isnan(A), 2);             % time × 1
                SEM.(elec).(fld) = std(A,[],2,'omitnan') ./ sqrt(N);
            end
        end

        save('SEMoutcome.mat','SEM');

        load('ERPoutcome.mat')
        % load("SEMoutcome.mat")
        timeline = -1000:1:1999;
    % ------------------------------------------------------------
    % Loop over all electrodes and generate the two figures
    % ------------------------------------------------------------
    for e = 1:numel(electrodeList)
        plotted_electrode = electrodeList{e};
    
        %% -------- pull data for this electrode -----------------
        el  = plotted_electrode;         % shorthand
        correct          = ERPout.(el).correct;
        incorrect        = ERPout.(el).incorrect;
    
        pre_correct      = ERPout.(el).pre_correct;
        post_correct     = ERPout.(el).post_correct;
    
        pre_incorrect    = ERPout.(el).pre_incorrect;
        post_incorrect   = ERPout.(el).post_incorrect;
    
    
        %% ---------- FIGURE 1 : stimulus ERPs -------------------
        fig1 = figure('Name',[el ' Outcome ERP'],'Visible','off');
        sgtitle(el)
    
        % (1) grand‑average correct vs incorrect
        subplot(3,2,1); hold on
        plot(timeline, mean(correct,  2,'omitnan'), "Color",'g');
        plot(timeline, mean(incorrect,2,'omitnan'), "Color",'r');
        legend({'correct stim','NoGo stim'}); title([el ' correct vs incorrect']);
        xlabel('Time (ms)'); ylabel('µV'); xlim([-100 800]);
        xline(0); yline(0);
    
        % (1) zoomed in correct vs incorrect
        subplot(3,2,2); hold on
        plot(timeline, mean(correct,  2,'omitnan'), "Color",'g');
        plot(timeline, mean(incorrect,2,'omitnan'), "Color",'r');
        legend({'correct stim','NoGo stim'}); title([el ' correct vs incorrect zoomed in']);
        xlabel('Time (ms)'); ylabel('µV'); xlim([-100 800]);
        xline(0); yline(0);
    
        % (2) individual correct
        subplot(3,2,3); hold on
        plot(timeline, correct);  % one line per participant
        title([el ' correct individuals']); xlim([-100 800]);
        xline(0); yline(0);
    
        % (3) individual incorrect
        subplot(3,2,4); hold on
        plot(timeline, incorrect);
        title([el ' incorrect individuals']); xlim([-100 800]);
        xline(0); yline(0);
        legend()
    
        % % (1) Pre‑correct vs Post‑correct
        % subplot(3,2,5); hold on
        % plot(timeline, mean(pre_correct,2,'omitnan'), 'Color',"b");
        % plot(timeline, mean(post_correct,2,'omitnan'), 'Color',"g");
        % legend({'pre','post'},'location','bestoutside'); title('correct before & after switch');
        % xlabel('Time (ms)'); ylabel('µV'); xlim([-200 800]);
        % xline(0); yline(0);
        % 
        % 
        % % (4) Pre‑incorrect vs Post‑correct
        % subplot(3,2,6); hold on
        % plot(timeline, mean(pre_incorrect,2,'omitnan'), 'Color',"b");
        % plot(timeline, mean(post_incorrect,2,'omitnan'), 'Color',"g");
        % legend({'pre','post'},'location','bestoutside'); title('incorrect - mean');
        % xlim([-200 800]); xline(0); yline(0);
    
        exportgraphics(fig1, fullfile(figure_output_folder, [fig1.Name '.pdf']), ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'none');  % Keeps fonts editable in Illustrator
    
      
    end
    
    fprintf('All electrode figures saved to: %s\n', figure_output_folder);


end


clear; close all; clc;

% -------- Paths --------
input_folder = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\EEG analysis\Winter 2026 outcomes\Epoch analysis';

addpath('C:\Users\khatinova\OneDrive - Nexus365\Pre_2026_Folders\Documents\MATLAB\eeglab2025.1.0')
eeglab nogui

valid_subjects = [3:12,14:18];

electrode = 'FCz';   % change if needed

all_rewp = [];

for s = valid_subjects

    subj = sprintf('Ox%02d',s);
    fprintf('Loading %s\n',subj);

    try
        EEG_correct   = pop_loadset([subj '_correct.set'], input_folder);
        EEG_incorrect = pop_loadset([subj '_incorrect.set'], input_folder);
    catch
        warning('Missing files for %s — skipping',subj);
        continue;
    end

    % Find electrode index
    chanIdx = find(strcmpi({EEG_correct.chanlocs.labels}, electrode),1);
    if isempty(chanIdx)
        warning('%s not found in %s',electrode,subj);
        continue;
    end

    % Compute subject-level ERPs
    erp_cor = squeeze(mean(EEG_correct.data(chanIdx,:,:),3,'omitnan'));
    erp_inc = squeeze(mean(EEG_incorrect.data(chanIdx,:,:),3,'omitnan'));

    % RewP
    rewp = erp_cor - erp_inc;

    all_rewp(:,end+1) = rewp;

end

% -------- Grand Average --------
grand_rewp = mean(all_rewp,2,'omitnan');
sem_rewp   = std(all_rewp,[],2,'omitnan') ./ sqrt(size(all_rewp,2));

time = EEG_correct.times;

% -------- Plot --------
figure;
hold on;

% Shaded SEM
fill([time fliplr(time)], ...
     [grand_rewp'+sem_rewp' fliplr(grand_rewp'-sem_rewp')], ...
     [0.7 0.7 1], 'EdgeColor','none', 'FaceAlpha',0.3);

plot(time, grand_rewp,'b','LineWidth',2);

xline(0,'--k');
yline(0,'k');

xlim([-200 800]);
xlabel('Time (ms)');
ylabel('\muV');
title(['Grand Average RewP (' electrode ')']);
set(gca,'FontSize',12);

fprintf('Grand average computed from %d subjects.\n', size(all_rewp,2));
