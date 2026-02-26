%% ===================== SETUP =====================


eeglab_path = 'C:\Users\khatinova\OneDrive - Nexus365\Pre_2026_Folders\Documents\MATLAB\eeglab2025.1.0';
addpath(eeglab_path);
eeglab;

% -------- Paths --------
rr_data_path   = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch RR\Data';
rr_study_path  = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch RR\Results\EEG_analysis';
rr_study_name  = 'CategorySwitchRR.study';

rr_sets_path = fullfile(rr_study_path, 'raw_sets');
if ~exist(rr_sets_path,'dir'), mkdir(rr_sets_path); end
if ~exist(rr_study_path,'dir'), mkdir(rr_study_path); end

%% ===================== FIND SUBJECTS =====================
subjects = dir(fullfile(rr_data_path, 'RR*'));
subjects = subjects([subjects.isdir]);

fprintf('Found %d RR subjects\n', numel(subjects));

ALLEEG = [];
CURRENTSET = 0;

%% ===================== IMPORT MFF FILES =====================
for s = 1:numel(subjects)

    subjID   = subjects(s).name;
    subjPath = fullfile(rr_data_path, subjID);

    fprintf('\n--- Importing %s ---\n', subjID);

    mffFile = dir(fullfile(subjPath, '*.mff'));
    if isempty(mffFile)
        warning('No .mff file found for %s — skipping', subjID);
        continue;
    end

    mffFullPath = fullfile(subjPath, mffFile(1).name);

    % ---------- Import MFF ----------
    EEG = pop_mffimport( ...
        {mffFullPath}, ...
        'code', 'eventtype', ...
        'reference', [] ...
    );

    EEG.subject = subjID;
    EEG.setname = sprintf('%s_raw', subjID);

    % ---------- Save SET ----------
    EEG = pop_saveset(EEG, ...
        'filename', [EEG.setname '.set'], ...
        'filepath', rr_sets_path, ...
        'savemode', 'onefile');

    % ---------- Store ----------
    CURRENTSET = CURRENTSET + 1;
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

end

%% ===================== CREATE RR STUDY =====================
% [STUDY, ALLEEG] = std_editset([], ALLEEG, ...
%     'name', 'categorySwitchRR', ...
%     'task', 'CategorySwitch', ...
%     'filename', rr_study_name, ...
%     'filepath', rr_study_path);

% Assign subject metadata
for i = 1:length(ALLEEG)
    [STUDY, ALLEEG] = std_editset(STUDY, ALLEEG, ...
        'commands', {{ ...
            'index', i, ...
            'subject', ALLEEG(i).subject, ...
            'condition', 'continuous' ...
        }}, ...
        'updatedat', 'on');
end

% Save study
[STUDY, ALLEEG] = pop_savestudy(STUDY, ALLEEG, ...
    'filename', rr_study_name, ...
    'filepath', rr_study_path);

fprintf('\n✅ RR MFF STUDY CREATED:\n%s\n', fullfile(rr_study_path, rr_study_name));

%% ===================== LOAD RR STUDY =====================
[STUDY, ALLEEG] = pop_loadstudy( ...
    'filename', 'RR_MFF.study', ...
    'filepath', rr_study_path);

fprintf('Loaded RR study with %d datasets\n', numel(ALLEEG));

%% ===================== PREPROCESSING =====================
for i = 1:numel(ALLEEG)

    EEG = ALLEEG(i);

    fprintf('\n--- Preprocessing %s (%d/%d) ---\n', EEG.subject, i, numel(ALLEEG));

    % -------- Skip if already done --------
    if isfield(EEG,'etc') && isfield(EEG.etc,'preprocessing_complete') ...
            && EEG.etc.preprocessing_complete
        fprintf('✔ Already preprocessed, skipping\n');
        continue;
    end

    % -------- Load data into memory (critical for .fdt safety) --------
    EEG = pop_loadset('filename', EEG.filename, 'filepath', EEG.filepath);
    EEG = eeg_checkset(EEG,'loaddata');

    %% ---------- STEP 1: TRIM ----------
    latencies = [EEG.event.latency];
    last_latency = max(latencies) + round(0.2 * EEG.srate); % +200 ms
    end_time_sec = last_latency / EEG.srate;

    EEG = pop_select(EEG, 'time', [0 end_time_sec]);
    EEG.setname = [EEG.setname '_trim'];

    %% ---------- STEP 2: FILTER ----------
    EEG = pop_eegfiltnew(EEG, 'locutoff', 1);
    EEG = pop_eegfiltnew(EEG, 'hicutoff', 40);
    EEG.setname = [EEG.setname '_filt'];

    %% ---------- STEP 3: RE-REFERENCE ----------
    % RR datasets: assume Cz reference (change here if needed)
    EEG = pop_reref(EEG, 61, 'keepref','on');
    EEG = pop_reref(EEG, []); % CAR
    EEG.setname = [EEG.setname '_CAR'];

    %% ---------- STEP 4: ASR ----------
    EEG = pop_clean_rawdata(EEG, ...
        'FlatlineCriterion', 5, ...
        'ChannelCriterion', 0.8, ...
        'Highpass', 'off', ...
        'BurstCriterion', 'off', ...
        'WindowCriterion', 'off', ...
        'BurstRejection', 'off');

    EEG.setname = [EEG.setname '_ASR'];

    %% ---------- STEP 5: ICA ----------
    if isempty(EEG.icaweights)
        rankData = rank(double(eeg_getdatact(EEG)));
        EEG = pop_runica(EEG, 'extended', 1, 'pca', rankData);
        EEG.setname = [EEG.setname '_ICA'];
    end

    %% ---------- STEP 6: ICLabel ----------
    EEG = iclabel(EEG);
    EEG.setname = [EEG.setname '_ICL'];

    %% ---------- SAVE WITH NEW FILENAME (SAFE) ----------
    
    new_filename = sprintf('%s_preprocessed.set', EEG.subject);
    new_filepath = output_folder_processed;
    
    EEG.etc.preprocessing_complete = true;
    
    EEG = pop_saveset(EEG, ...
        'filename', new_filename, ...
        'filepath', new_filepath, ...
        'savemode', 'onefile');
    
    % Update EEG struct
    EEG.filename = new_filename;
    EEG.filepath = new_filepath;
    
    % Put back into ALLEEG
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, i);
    
    % 🔴 CRITICAL: update STUDY pointers
    STUDY.datasetinfo(i).filename = new_filename;
    STUDY.datasetinfo(i).filepath = new_filepath;
    
    fprintf('✔ Saved and relinked: %s\n', fullfile(new_filepath,new_filename));
    

end

%% ===================== SAVE STUDY =====================
[STUDY, ALLEEG] = pop_savestudy(STUDY, ALLEEG, ...
    'filename', 'RR_MFF.study', ...
    'filepath', rr_study_path);

fprintf('\n✅ RR STUDY updated with PREPROCESSED datasets\n');
eeglab redraw;

