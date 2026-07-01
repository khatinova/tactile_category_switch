%% ===================== SETUP =====================


eeglab_path = 'C:\Users\khatinova\OneDrive - Nexus365\Pre_2026_Folders\Documents\MATLAB\eeglab2025.1.0';
addpath(eeglab_path);
eeglab;

% -------- Paths --------
rr_data_path   = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch RR\Data\det_or_prob_and_conf';
rr_study_path  = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch RR\Results\EEG_analysis';
rr_study_name  = 'CategorySwitchRR.study';

rr_sets_path = fullfile(rr_study_path, 'raw_sets');
if ~exist(rr_sets_path,'dir'), mkdir(rr_sets_path); end
if ~exist(rr_study_path,'dir'), mkdir(rr_study_path); end

%% ===================== FIND SUBJECTS =====================
subjects = dir(fullfile(rr_data_path, 'Nc*'));
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
    EEG = pop_mffimport(mffFullPath, 'code', 0, 1);

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
[STUDY, ALLEEG] = std_editset([], ALLEEG, ...
    'name', 'categorySwitchRR', ...
    'task', 'CategorySwitch', ...
    'filename', rr_study_name, ...
    'filepath', rr_study_path);

% Assign subject metadata
for i = 1:length(ALLEEG)

    IDsplit = strsplit(ALLEEG(i).subject, '_');
    [STUDY, ALLEEG] = std_editset(STUDY, ALLEEG, ...
        'commands', {{ ...
            'index', i, ...
            'subject',IDsplit{1}, ...
            'condition', IDsplit{3}, ...
            'group', 'RR'...
        }}, ...
        'updatedat', 'on');
end

% Save study
[STUDY, ALLEEG] = pop_savestudy(STUDY, ALLEEG, ...
    'filename', rr_study_name, ...
    'filepath', rr_study_path, ...
    'resavedatasets', 'on');

fprintf('\n✅ RR MFF STUDY CREATED:\n%s\n', fullfile(rr_study_path, rr_study_name));

