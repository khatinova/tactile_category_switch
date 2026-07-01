clear; close all; clc;
%% ============================================================
% LOOP 1: Preprocess + ICA + Auto-Flag (NO removal yet)
% =============================================================

clear; close all; clc;

% === Setup Paths ===
eeglab_path = 'C:\Users\khatinova\OneDrive - Nexus365\Pre_2026_Folders\Documents\MATLAB\eeglab2025.1.0';
addpath(eeglab_path);
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

data_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Data';
study_filepath = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\EEG analysis\Pilot2_Jan2026';
output_folder_processed = fullfile(study_filepath, 'full_processed_data');

if ~exist(output_folder_processed,'dir')
    mkdir(output_folder_processed);
end

subjects = dir(fullfile(data_path,'Ox*'));
subjects = subjects([subjects.isdir]);

valid_subjects = [4:12,14:18];  % your selection

for i = valid_subjects

    subjID = subjects(i).name;
    fprintf('\n=== Processing %s ===\n', subjID);

    subjPath = fullfile(data_path, subjID);
    curryFile = dir(fullfile(subjPath,'Acquisition*.dat'));
    if isempty(curryFile)
        warning('No Curry file for %s', subjID);
        continue
    end

    %% LOAD
    EEG = loadcurry(fullfile(subjPath,curryFile(1).name), ...
        'KeepTriggerChannel','True','CurryLocations','False');
    EEG.subject = subjID;
    EEG = eeg_checkset(EEG,'loaddata');

    %% TRIM 10s AFTER LAST EVENT
    last_latency = max([EEG.event.latency]);
    trim_sample  = last_latency + (10 * EEG.srate);
    trim_sample  = min(trim_sample, EEG.pnts);
    EEG = pop_select(EEG,'point',[1 trim_sample]);

    %% FILTER 0.5–40 Hz
    EEG = pop_eegfiltnew(EEG,'locutoff',0.5);
    EEG = pop_eegfiltnew(EEG,'hicutoff',40);

    %% Cz → CAR
    cz_index = find(strcmpi({EEG.chanlocs.labels},'Cz'));
    EEG = pop_reref(EEG, cz_index,'keepref','on');
    EEG = pop_reref(EEG, []);

    original_chanlocs = EEG.chanlocs;
    original_nbchan   = EEG.nbchan;

    %% ASR - Mild 1st rejection
    EEG = pop_clean_rawdata( ...
        EEG,'FlatlineCriterion',5, ...
        'ChannelCriterion',0.8,...
        'LineNoiseCriterion',4, ...
        'Highpass',[0.25 0.75], ...
        'BurstCriterion',40,...
        'WindowCriterion',0.25, ...
        'BurstRejection','on', ...
        'Distance','Euclidian',...
        'WindowCriterionTolerances',[-Inf 10] );
    
    %% INTERPOLATE REMOVED CHANNELS (CORRECT WAY)
    if EEG.nbchan < original_nbchan
        fprintf('Interpolating %d removed channels\n', ...
            original_nbchan - EEG.nbchan);
        
        EEG = pop_interp(EEG, original_chanlocs, 'spherical');
    end

    %% Re-reference again after ASR
    EEG = pop_reref(EEG, []);

    %% ICA (Picard recommended, but runica if you prefer)
    EEG = pop_runica(EEG, 'extended', 1);

    %% ICLabel
    EEG = pop_iclabel(EEG,'default');

    %% AUTO FLAG USING YOUR CRITERIA (0.5–1)
    EEG = pop_icflag(EEG, ...
        [0 0.2;     % Brain keep
         0.5 1;       % Muscle reject
         0.5 1;       % Eye reject
         0.5 1;       % Heart reject
         0.5 1;       % LineNoise reject
         0.5 1;       % ChannelNoise reject
         0.5 1]);     % Other reject

    fprintf('Auto-flagged ICs: ');
    EEG.reject.icareject = EEG.reject.gcompreject;
    disp(find(EEG.reject.icareject==1));

    %% SAVE FLAGGED VERSION (NO REMOVAL)
    % EEG.setname = [subjID '_ICAflagged'];
    % pop_saveset(EEG, ...
    %     'filename',[subjID '_ICAflagged.set'], ...
    %     'filepath',output_folder_processed, ...
    %     'savemode','onefile');

    fprintf('\nManual inspection for %s\n', subjID);

    %% Show all ICs in grid (red = auto-flagged)
    % pop_topoplot(EEG,0,1:size(EEG.icaweights,1), ...
    %     'ICs (red = auto-flagged)',0,'electrodes','on');

    %% Detailed properties (grouped)
    % pop_prop(EEG,0,1:size(EEG.icaweights,1),NaN,{'freqrange',[2 40]});

 

    %% Remove flagged ICs (including any you manually changed)
    EEG = pop_subcomp(EEG, [], 0);

    %% ASR2 - Aggressive data rejevtion, based on good quality data portions
    EEG = pop_clean_rawdata( EEG, ...
        'FlatlineCriterion','off', ...
        'ChannelCriterion','off',...
        'LineNoiseCriterion','off', ...
        'Highpass','off', ...
        'BurstCriterion',20,...
        'WindowCriterion',0.25, ...
        'BurstRejection','on', ...
        'Distance','Euclidian',...
        'WindowCriterionTolerances',[-Inf 7] );


    EEG = pop_reref(EEG, [], 'refica', 'remove');
    EEG = pop_runica(EEG, 'extended', 1);

    %% OPTINAL: ICLabel v2
    % EEG = pop_iclabel(EEG,'default');
    % 
    % %% AUTO FLAG USING YOUR CRITERIA (0.5–1)
    % EEG = pop_icflag(EEG, ...
    %     [0 0.2;     % Brain keep
    %      0.5 1;       % Muscle reject
    %      0.5 1;       % Eye reject
    %      0.5 1;       % Heart reject
    %      0.5 1;       % LineNoise reject
    %      0.5 1;       % ChannelNoise reject
    %      0.5 1]);     % Other reject
    % 
    % fprintf('Auto-flagged ICs round 2: ');
    % EEG.reject.icareject = EEG.reject.gcompreject;
    % disp(find(EEG.reject.icareject==1));
    % EEG = pop_subcomp(EEG, [], 0);
    % EEG = pop_reref(EEG, [], 'interpchan', 'on');

    EEG.setname = [subjID '_cleaned'];
    pop_saveset(EEG, ...
        'filename',[subjID '_cleaned.set'], ...
        'filepath',output_folder_processed, ...
        'savemode','onefile');


end

fprintf('\n=== LOOP 2 COMPLETE ===\n');


%% ============================================================
% LOOP 3: Epoching
% =============================================================

files = dir(fullfile(processed_path,'*_cleaned.set'));

for f = 1:length(files)

    EEG = pop_loadset('filename',files(f).name,'filepath',processed_path);
    subjID = EEG.subject;

    fprintf('\nEpoching %s\n', subjID);
    if f < 9
        outcome_trigger = {'10','11'};
    else
        outcome_trigger = {'31', '32', '33', '34'};
    end

    %% Example stimulus-locked epoch
    EEG_ep = pop_epoch(EEG, outcome_trigger, [-0.2 1], 'epochinfo','yes');
    EEG_ep = pop_rmbase(EEG_ep, [-200 0]);

    EEG_ep.setname = [subjID '_outcome'];
    pop_saveset(EEG_ep, ...
        'filename',[subjID '_outcome.set'], ...
        'filepath',processed_path, ...
        'savemode','onefile');

end

fprintf('\n=== LOOP 3 COMPLETE ===\n');
