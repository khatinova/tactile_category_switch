% clear all
% close all
data_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch';
addpath(genpath(data_path));
cd(data_path);
load("Results\Behav_Results_salient_mod_switch\group_performance_results.mat")

addpath('\\humerus\PHARM_BANERJEE\data\Personal Folders\Klara\EEG_analysis\eeglab2025.0.0')
% eeglab
 
% unique_block_string          = unique(data(2:end,1))
% block_identity_ordered       = [unique_block_string(131:end);unique_block_string(1:130)]
data_titles                    = fieldnames(group_results);
output_folder = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\EEGLAB_STUDY_v4\Epoched_by_trigger';
input_folder = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\EEGLAB_STUDY_v4\full_processed_data';
figure_output_folder = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\EEGLAB_STUDY_v4\ERP Figures\somatosensory ERPs';


valid_participants=2:8;

% organise the behavioural data into a matrix that can be collated with the
% EEG epochs
for i=1:8
    participant_string       = data_titles{i};   
    participant_data         = group_results.(participant_string);
    data_by_participant{i}   = group_results.(participant_string);
end

if ~exist('data_by_participant.mat')
    save('data_by_participant.mat', "data_by_participant")
end

%%%%%%%%%%%% COMPUTING OPTIONS %%%%%%%%%%%%%
save_stimulus_epochs   = 0; % by go and nogo, plus stimuli that switch (triggers 3 and 6)
stimulus_erp_analysis  = 0;
ERP_plotting           = 1; % requires you to have outcome_erp_analysis = 1 or have the variables calculated and loaded from directory


trigSets = { ...
    'Go'           , 1:4 ;  ... % GO   (all pre- & post-reversal)
    'NoGo'         , 5:8 ;  ... % No-GO
    'Postrev_GO'   , 3   ;  ... % GO   after reversal only
    'Postrev_NoGo' , 6   ; ... % No-GO after reversal only
    'ALL_Stim'     , 1:8};     
    baseline = [-200 0];
    epoch_duration = [-1 1];

if save_stimulus_epochs==1


    for s = 1:8
        if s == 1, continue, 
        end          % skip p1
    
        %% Load the subject's dataset ---------------------------------------
        fname  = sprintf('subject%d_trimmed_filtered_CzRef_CAR_ASR_ICA_ICL rejected.set',s);
        EEG_   = pop_loadset(fname, ...
                '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\EEGLAB_STUDY_v4\full_processed_data');
    
        % Numeric vector of event codes
        evt = [EEG_.urevent.type];        % 1 × nEvents
    
        %% Loop over the four trigger groups --------------------------------
        for k = 5 % 1:size(trigSets,1)
            label   = trigSets{k,1};
            trgList = trigSets{k,2};      % e.g. 1:4 or [3]
    
            if numel(trgList) == 1
                trgArg = num2str(trgList);                                    % keep numeric
            else
                trgArg = arrayfun(@num2str, trgList, 'UniformOutput', false);
                % trgArg is now {'1' '2' '3' '4'} – legal for pop_epoch
            end
            
            EEGtmp = pop_epoch(EEG_, trgArg, epoch_duration);                % works!
            EEGtmp = pop_rmbase(EEGtmp, baseline);
        
            % ---------- save ----------
            % saveName = sprintf('p%d_%s_epochs.set',s,label);
            % pop_saveset(EEGtmp,'filename', saveName, 'filepath', output_folder);
        end
    end
end


% extract specific electrodes
electrodeList = {'FPZ', 'AFZ', 'FZ', 'FCZ', 'CZ', 'PZ', 'C1', 'C3', 'CP1', 'CP3'};

if stimulus_erp_analysis==1
    
    for participant = valid_participants                   % skip p1 if needed

       

        % load the behavioural data
         beh_data_one_participant=data_by_participant{participant};  
         % go_indices_check = beh_data_one_participant.goTrial==1;
         % go_indices_check = go_indices_check(:)';
         % nogo_indices_check = beh_data_one_participant.goTrial==0;
         % nogo_indices_check = nogo_indices_check(:)';
         % stim_number = beh_data_one_participant.stimType;
         % stim_number = stim_number(:)';
         num_blocks = height(beh_data_one_participant.correct);
         num_trials_per_block = width(beh_data_one_participant.correct);
         total_trials = num_blocks*num_trials_per_block; % total number of trials in the behavioural dataset

         % load the eeg data

        label = trigSets{5,1}; % load only the data epoched by all stimuli, stored in 5th index of this strcuture
        fname = sprintf('p%d_ALL_stim_epochs.set',participant);
        EEG_   = pop_loadset(fname,output_folder);
   
        %% ---------- convert trigger codes & trial number -------------------
        evtCodes = str2double({EEG_.epoch.eventtype});     % 1×nEpoch
        nChan    = size(EEG_.data,1);
        nSamples    = size(EEG_.data,2);
    
        nEEGTrials  = numel(evtCodes);            % total number of EEG recorded trails
        % ------------------------------------------------------
        % Case-1  (practice at start)   nEEGTrials  > total_trials
        %         → trim the FIRST (nEEGTrials-total_trials) epochs
        %
        % Case-2  (late start)          nEEGTrials  < total_trials
        %         → keep everything; later we'll treat missing early
        %           trials as NaN by shifting indices
        % ------------------------------------------------------
        if nEEGTrials > total_trials                       % practice present
            n_practice = nEEGTrials - total_trials;        % how many to drop
            keepTrials = (n_practice+1) : nEEGTrials;      % practise+1 … end
            EEG_       = pop_select(EEG_ ,'trial', keepTrials);
        
            % After trimming, trials now line up 1…total_trials
            trialOffset = 0;
        
        else                                               % ❷ recording started late
            % keepTrials = 1:nEEGTrials  (no trimming needed)
            trialOffset = total_trials - nEEGTrials;       % missing at the FRONT
        end

        trialNum = trialOffset + (1:nEEGTrials);
        total_pre_trials = 20;
        total_post_trials = 21;
                
        if ~exist(sprintf('p%d_outcome_valid_trials_only.set', participant))
            pop_saveset(EEG_, 'filename', sprintf('p%d_outcome_valid_trials_only.set', participant));
        end  
        %% ---------- pre-allocate trial-aligned containers ------------------
        fieldNames = {'GOstim', 'NoGOstim', 'trig1','trig2', 'trig3', 'trig4', 'trig6', 'trig7'};
        for i = 1:numel(fieldNames)
            fieldNames{end+1} = ['pre_' fieldNames{i}];
            fieldNames{end+1} = ['post_' fieldNames{i}];
        end

        % Initialize structure to hold reversal-aligned data
        nChan = size(EEG_.data,1);
        nSamples = size(EEG_.data,2);

        % Initialize structure with NaNs for all fields
        for f = 1:numel(fieldNames)
            EEGstruct.(fieldNames{f}) = nan(nChan, nSamples, total_trials);
            
            % add fieldnames for reversal aligned trials
            if startsWith(fieldNames{f}, 'pre_')
                EEGstruct.(fieldNames{f}) = nan(nChan, nSamples, total_pre_trials);
            elseif startsWith(fieldNames{f}, 'post_')
                EEGstruct.(fieldNames{f}) = nan(nChan, nSamples, total_post_trials);
            end
        end
            
        %% ---------- drop every epoch into its trial slot -------------------
        for epoch = 1:length(keepTrials)-trialOffset
            t = trialNum(epoch);               % behavioural trial index (accurately offset from the EEG trial)
            code = evtCodes(epoch);            % 1 … 8 trigger
    
            if ismember(code, 1:4)          % GO stimuli
                EEGstruct.GOstim(:,:,t) = EEG_.data(:,:,epoch);
    
            elseif ismember(code, 5:8)      % No-Go stimuli
                EEGstruct.NoGOstim(:,:,t) = EEG_.data(:,:,epoch);
            end

            % postreversal go only
            if code == 3
                EEGstruct.trig3(:,:,t) = EEG_.data(:,:,epoch);
            end
            % pre and postreversal go
            if code == 1
                EEGstruct.trig1(:,:,t) = EEG_.data(:,:,epoch);
            end
            % postreversal nogo only
            if code == 6 
                EEGstruct.trig6(:,:,t) = EEG_.data(:,:,epoch);
            end
            % pre and postreversal nogo
            if code == 4
                EEGstruct.trig4(:,:,t) = EEG_.data(:,:,epoch);
            end
            if code == 7
                EEGstruct.trig7(:,:,t) = EEG_.data(:,:,epoch);
            end
            if code == 2
                EEGstruct.trig2(:,:,t) = EEG_.data(:,:,epoch);
            end


        end

        pre_reversal_signal = [];
        post_reversal_signal = [];

        % get the trials that are before and after a reversal to extract relevant epochs
        for j = 1:total_trials
            for block = 1:num_blocks

                rev_cont_idx = beh_data_one_participant.reversal_trials(block) + ((block-1)*100);
                
                globalRevIdx(block, participant) = rev_cont_idx;

                if j == rev_cont_idx

                    pre_reversal_signal  = [pre_reversal_signal, j-21:j-1];
                    post_reversal_signal = [post_reversal_signal, j:j+20];

                end
            end


        end

        EEGstruct.EEG_prerev_GO = EEGstruct.GOstim(:,:,pre_reversal_signal);
        EEGstruct.EEG_prerev_NoGO = EEGstruct.NoGOstim(:,:,pre_reversal_signal);
        EEGstruct.EEG_postrev_GO = EEGstruct.GOstim(:,:,post_reversal_signal);
        EEGstruct.EEG_postrev_NoGO = EEGstruct.NoGOstim(:,:,post_reversal_signal);
        EEGstruct.pre_trig1  = EEGstruct.trig1(:,:,pre_reversal_signal);
        EEGstruct.post_trig1 = EEGstruct.trig1(:,:,post_reversal_signal); 
        EEGstruct.post_trig3 = EEGstruct.trig3(:,:,post_reversal_signal);  
        EEGstruct.pre_trig4  = EEGstruct.trig4(:,:,pre_reversal_signal);
        EEGstruct.post_trig4 = EEGstruct.trig4(:,:,post_reversal_signal);
        EEGstruct.post_trig6 = EEGstruct.trig6(:,:,post_reversal_signal);
        EEGstruct.pre_trig2 = EEGstruct.trig2(:,:,pre_reversal_signal);
        EEGstruct.pre_trig7 = EEGstruct.trig7(:,:,pre_reversal_signal);
        EEGstruct.post_trig2 = EEGstruct.trig2(:,:,post_reversal_signal);
        EEGstruct.post_trig7 = EEGstruct.trig7(:,:,post_reversal_signal);


        % ---------- build map: electrode name → index in this dataset -------
        labelList  = {EEG_.chanlocs.labels};
        ERPout.Chanlocs{participant} = labelList;

        eIdx = nan(1,numel(electrodeList));
        for e = 1:numel(electrodeList)
            eIdx(e) = find(strcmpi(labelList, electrodeList{e}),1);
            if isempty(eIdx(e))
                error('Electrode "%s" not found in participant %d.', electrodeList{e});
            end
        end
 
        % Loop over all fields to compute mean ERP (chan x time) across trials for each participant
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
    save('ERPstim.mat','ERPout');
    fprintf('Saved %s\n', outputMAT);
end


if ERP_plotting == 1    
   load('ERPstim.mat')
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
        save('SEMstim.mat','SEM');


        % 
        % load('Stimulus_ERPs_ALL.mat')
        % load('Outcome_PreRev_ERPs_ALL.mat')
        % load('Outcome_PostRev_ERPs_ALL.mat')
     
        load("SEMstim.mat")
        timeline = -1000:1:999;
% ------------------------------------------------------------
% Loop over all electrodes and generate the two figures
% ------------------------------------------------------------
for e = 1:numel(electrodeList)
    plotted_electrode = electrodeList{e};

    %% -------- pull data for this electrode -----------------
    el  = plotted_electrode;         % shorthand
    GO          = ERPout.(el).GOstim;
    NoGO        = ERPout.(el).NoGOstim;

    pre_GO      = ERPout.(el).pre_trig2;
    post_NoGO   = ERPout.(el).post_trig6;

    pre_NoGO    = ERPout.(el).pre_trig7;
    post_GO     = ERPout.(el).post_trig3;

    % always‑Go / always‑NoGo (optional, remove if not needed)
    pre_alwG    = ERPout.(el).pre_trig1;
    post_alwG   = ERPout.(el).post_trig1;
    pre_alwN    = ERPout.(el).pre_trig4;
    post_alwN   = ERPout.(el).post_trig4;

    %% ---------- FIGURE 1 : stimulus ERPs -------------------
    fig1 = figure('Name',[el ' ALL ERP plots together'],'Visible','off');
    sgtitle(el)

    % (1) grand‑average GO vs NoGo
    subplot(2,2,1); hold on
    plot(timeline, mean(GO(:,valid_participants),  2,'omitnan'));
    plot(timeline, mean(NoGO(:,valid_participants),2,'omitnan'));
    legend({'Go stim','NoGo stim'}); title([el ' Go vs NoGo']);
    xlabel('Time (ms)'); ylabel('µV'); xlim([-500 1000]);
    xline(0); yline(0);

    % (2) individual GO
    subplot(2,2,3); hold on
    plot(timeline, GO(:,valid_participants));  % one line per participant
    title([el ' GO individuals']); xlim([-100 800]);
    xline(0); yline(0);

    % (3) individual NoGo
    subplot(2,2,4); hold on
    plot(timeline, NoGO(:,valid_participants));
    title([el ' NoGo individuals']); xlim([-100 800]);
    xline(0); yline(0);

    % (4) mean ± SEM ribbon
    sem_GO   = SEM.(el).GOstim;
    sem_NoGO = SEM.(el).NoGOstim;
    subplot(2,2,2); hold on
    shadedErrorBar(timeline, mean(GO(:,valid_participants),2,'omitnan'), sem_GO);
    shadedErrorBar(timeline, mean(NoGO(:,valid_participants),2,'omitnan'), sem_NoGO);
    plot(timeline, mean(GO(:,valid_participants),  2,'omitnan'));
    plot(timeline, mean(NoGO(:,valid_participants),2,'omitnan'));
    title([el ' Go vs NoGo ± SEM']); xlim([-500 1000]);
    xline(0); yline(0);

    exportgraphics(fig1, fullfile(figure_output_folder, [fig1.Name '.pdf']), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'none');  % Keeps fonts editable in Illustrator

    close(fig1);

    %% ---------- FIGURE 2 : reversal‑aligned ----------------
    fig2 = figure('Name',[el ' Reversed stimuli ALL plots'],'Visible','off');
    sgtitle(el)

    % (1) Pre‑GO vs Post‑NoGo
    subplot(2,3,1); hold on
    plot(timeline, mean(pre_GO   (:,valid_participants),2,'omitnan'));
    plot(timeline, mean(post_NoGO(:,valid_participants),2,'omitnan'));
    legend({'pre‑GO','post‑NoGo'}); title('Go→NoGo switch');
    xlabel('Time (ms)'); ylabel('µV'); xlim([-200 800]);
    xline(0); yline(0);

    % (2) Individuals (preGO & postGO)
    subplot(2,3,2); hold on
    plot(timeline, pre_GO  (:,valid_participants),'Color','b');
    plot(timeline, post_GO (:,valid_participants),'Color','r');
    title('Go→NoGo individuals'); xlim([-200 800]); xline(0); yline(0);

    % (3) Mean ± SEM for Go→NoGo
    subplot(2,3,3); hold on
    plot(timeline, pre_alwG(:,valid_participants),'Color',[0 .6 0 .4]);
    plot(timeline, post_alwG(:,valid_participants),'Color', [.8 0 0 .4]);
    title('Go → Go (stim 1) ± SEM'); xlim([-200 800]); xline(0); yline(0);

    % shadedErrorBar(timeline, mean(pre_GO(:,valid_participants),2,'omitnan'),  SEM.(el).pre_trig2);
    % shadedErrorBar(timeline, mean(post_NoGO(:,valid_participants),2,'omitnan'),SEM.(el).post_trig6);
    % plot(timeline, mean(pre_GO   (:,valid_participants),2,'omitnan'));
    % plot(timeline, mean(post_NoGO(:,valid_participants),2,'omitnan'));
    % title('Go→NoGo ± SEM'); xlim([-200 800]); xline(0); yline(0);

    % (4) Pre‑NoGo vs Post‑GO
    subplot(2,3,4); hold on
    plot(timeline, mean(pre_NoGO(:,valid_participants),2,'omitnan'));
    plot(timeline, mean(post_GO (:,valid_participants),2,'omitnan'));
    legend({'pre‑NoGo','post‑GO'}); title('NoGo→Go switch');
    xlim([-200 800]); xline(0); yline(0);

    % (5) Individuals (preNoGo & postNoGo)
    subplot(2,3,5); hold on
    plot(timeline, pre_NoGO(:,valid_participants),'Color',[0 .6 0 .4]);
    plot(timeline, post_GO(:,valid_participants),'Color', [.8 0 0 .4]);
    title('NoGo→Go individuals'); xlim([-200 800]); xline(0); yline(0);

    % (6) Mean ± SEM for NoGo→Go
    subplot(2,3,6); hold on
    plot(timeline, pre_alwN(:,valid_participants),'Color','b');
    plot(timeline, post_alwN(:,valid_participants),'Color', 'r');
    title('Go → Go (stim 1) ± SEM'); xlim([-200 800]); xline(0); yline(0);


    % shadedErrorBar(timeline, mean(post_GO (:,valid_participants),2,'omitnan'), SEM.(el).post_trig3);
    % shadedErrorBar(timeline, mean(pre_NoGO(:,valid_participants),2,'omitnan'), SEM.(el).pre_trig7);
    % plot(timeline, mean(pre_NoGO(:,valid_participants),2,'omitnan'));
    % plot(timeline, mean(post_GO (:,valid_participants),2,'omitnan'));
    % title('NoGo→Go ± SEM'); xlim([-200 800]); xline(0); yline(0);

   exportgraphics(fig2, fullfile(figure_output_folder, [fig2.Name '.pdf']), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'none');  % Keeps fonts editable in Illustrator

    close(fig2);
end

fprintf('All electrode figures saved to: %s\n', figure_output_folder);


end