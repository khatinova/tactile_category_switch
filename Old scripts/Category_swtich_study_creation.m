clear; close all; clc;

% === Setup Paths ===
eeglab_path = 'C:\Users\khatinova\OneDrive - Nexus365\Pre_2026_Folders\Documents\MATLAB\eeglab2025.1.0';
addpath(eeglab_path);
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

data_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Data';
study_filepath = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\EEG analysis\Winter 2026 outcomes';
ica_path = fullfile(study_filepath, 'ICA_decomposed');

subjects = dir(fullfile(data_path,'Ox*'));
subjects = subjects([subjects.isdir]);
valid_subjects = [11,12,14:numel(subjects)];
% 
pop_editoptions('option_storedisk', 0); % one dataset at a time
cd(study_filepath)
preprocess =1;


if preprocess == 1
    for i = valid_subjects
    
        subjID = subjects(i).name;
        fprintf('\n==== Processing %s ====\n', subjID);

        if ~exist([subjID, '_cleaned.set'])
    
            subjPath = fullfile(data_path, subjID);
            curryFile = dir(fullfile(subjPath,'Acquisition*.dat'));
            if isempty(curryFile), continue; end
        
            EEG = loadcurry(fullfile(subjPath,curryFile(1).name), ...
                'KeepTriggerChannel','True','CurryLocations','False');
            EEG.subject = subjID;
        
            %% TRIM 10s AFTER LAST EVENT
            last_latency = max([EEG.event.latency]);
            trim_sample  = min(last_latency + 10*EEG.srate, EEG.pnts);
            EEG = pop_select(EEG,'point',[1 trim_sample]);
        
            %% FILTER
            EEG = pop_eegfiltnew(EEG,'locutoff',0.5);
            EEG = pop_eegfiltnew(EEG,'hicutoff',40);
        
            %% REFERENCE
            cz_index = find(strcmpi({EEG.chanlocs.labels}, 'Cz'));
            EEG = pop_reref(EEG, cz_index, 'keepref','on');  % temporarily restore Cz
            EEG = pop_reref(EEG, []);
        
            %% CLEAN (ASR)
            EEG = pop_clean_rawdata(EEG,'FlatlineCriterion',5,...
                'ChannelCriterion',0.87,'LineNoiseCriterion',4,...
                'Highpass',[0.25 0.75],'BurstCriterion',20,...
                'WindowCriterion',0.25,'BurstRejection','on');
        
            EEG = pop_reref(EEG, []);
        
            EEG.setname = [subjID '_cleaned'];
            % 
            EEG = pop_saveset(EEG, ...
                'filename', [EEG.subject '_cleaned.set'], ...
                'filepath', study_filepath, ...
                'savemode','onefile');

         end

        % if  ~exist(fullfile(study_filepath, 'ICA_decomposed', [subjID, '_ICA.set']))
        % 
            if exist([subjID, '_cleaned.set'])
                EEG = pop_loadset([subjID, '_cleaned.set']);
            end
            EEG = pop_runica(EEG, 'extended', 1);
         % end
       
        % if  ~exist( [subjID, '_ICA_pruned.set'])

            % if exist(fullfile(study_filepath, 'ICA_decomposed', [subjID, '_ICA.set']))
            %     EEG = pop_loadset(fullfile(study_filepath, 'ICA_decomposed', [subjID, '_ICA.set']));
            % end
            % ----- ICLabel -----
            EEG= pop_iclabel(EEG,'default');
            
            % ----- Flag artifactual components -----
            EEG= pop_icflag(EEG, ...
                [0 0.2;      % Brain keep
                 0.5 1;      % Muscle
                 0.5 1;      % Eye
                 0.5 1;      % Heart
                 0.5 1;      % LineNoise
                 0.5 1;      % ChannelNoise
                 0.5 1]);    % Other
            
            % ----- Remove flagged components -----
            EEG= pop_subcomp(EEG, []);
            
            % ----- Re-reference AFTER IC removal -----
            EEG= pop_reref(EEG, [], ...
                'huber', 25, ...
                'interpchan', [], ...
                'refica', 'remove');
    
        
            % [ALLEEG, EEG] = eeg_store(ALLEEG, EEG);
    
            EEG = pop_saveset(EEG, ...
                'filename', [EEG.subject '_ICA_pruned.set'], ...
                'filepath', study_filepath, ...
                'savemode','onefile');
    
        % end

    end
    
    % 
    % [STUDY, ALLEEG] = std_editset([], ALLEEG, ...
    % 'name', 'CategorySwitch.study', ...
    % 'task', 'CategorySwitch', ...
    % 'filename', 'CategorySwitch.study', ...
    % 'filepath', study_filepath);
end



% STUDY = [];
% commands = {};
% 
% block_struct = {}; % for condition data/study metadata
% stim_config = {}; % for condition/study metadata
% task_version = {'v1', 'v2'}; % which taske v1 = 1-8, v2 = 9-onwards
% 
% for k = 1:length(ALLEEG)
%     commands{k} = { ...
%         'index', k, ...
%         'subject', ALLEEG(k).subject, ...
%         'condition', 'continuous'};
% end
% 
% % [STUDY, ALLEEG] = std_editset(STUDY, ALLEEG, ...
% %     'name','SalientModSwitch', ...
% %     'commands', commands);
% [STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG);
% 
% % optional: picard PCA
% plugin_askinstall('picard','picard',1);
% 
% % ----- ICA per subject -----
% for i = 8:length(ALLEEG)
%     fprintf('Running ICA for %s\n', ALLEEG(i).subject);
% 
%     % OPTIONAL but recommended: rank reduction
%     ncomps = ALLEEG(i).nbchan - 1;  % after average reference
%     ALLEEG(i) = pop_runica(ALLEEG(i), 'extended', 1);
% end
% 
% % ---------- SAVE ICA data --------
% ica_path = fullfile(study_filepath, 'ICA_decomposed');
% if ~exist(ica_path,'dir'), mkdir(ica_path); end
% 
% for i = 1:length(ALLEEG)
% 
%     EEG = ALLEEG(i);
% 
%     if isfield(EEG,'icaweights') && ~isempty(EEG.icaweights)
% 
%         EEG.setname = [EEG.subject '_ICA.set'];
% 
%         EEG = pop_saveset(EEG, ...
%             'filename', [EEG.subject '_ICA.set'], ...
%             'filepath', ica_path, ...
%             'savemode','onefile');
% 
%         fprintf('Saved ICA file for %s\n', EEG.subject);
% 
%     else
%         fprintf('No ICA found for %s — skipping save\n', EEG.subject);
%     end
% end
% 
% % ----- Save cleaned datasets -----
% for i = 1:length(ALLEEG)
%     EEG = ALLEEG(i);
%     EEG = pop_saveset(EEG, ...
%         'filename',[EEG.setname '_ICAclean.set'], ...
%         'filepath', study_filepath, ...
%         'savemode','onefile');
%     [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, i);
% end
% 
% % parameters: change windows as needed
% groupA_triggers = {'10','11'};         % subjects 1-8
% groupB_triggers = {'31','32','33','34'};% subjects 9+ (your earlier triggers)
% % epoch windows (in seconds) — adjust per timing difference
% epoch_limits = [-0.2 1.0];      % example for group A (10/11)
% epoch_baseline = [-0.2 0];            % baseline for group A
% 
% 
% %% Epoching
% ALLEEG = pop_epoch(ALLEEG, ...
%     {'31','32','33','34'}, epoch_limits, 'epochinfo','yes');
% 
% ALLEEG = pop_rmbase(ALLEEG, epoch_baseline);
% 
% STUDY = std_maketrialinfo(STUDY, ALLEEG);
% 
% STUDY = std_makedesign(STUDY, ALLEEG, 1, ...
%     'name','OutcomeDesign', ...
%     'variable1','type', ...
%     'values1',{'31','32','33','34'}, ...
%     'vartype1','categorical');
% 
% [STUDY, ALLEEG] = std_precomp(STUDY, ALLEEG, {}, ...
%     'savetrials','on', ...
%     'rmicacomps','on', ...
%     'interp','on', ...
%     'recompute','on', ...
%     'erp','on');
% 
% STUDY = pop_erpparams(STUDY,'topotime',350);
% STUDY = std_erpplot(STUDY, ALLEEG, 'design',1);
