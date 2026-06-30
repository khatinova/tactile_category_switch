% % === Setup Paths ===

eeglab_path = 'C:\Users\khatinova\OneDrive - Nexus365\Pre_2026_Folders\Documents\MATLAB\eeglab2025.1.0';
addpath(eeglab_path)
cd(eeglab_path)
% eeglab;

data_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Data';
study_filepath = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\EEG analysis\Pilot2_Jan2026';
study_filename = 'Pilot2_Jan2026.study';
processed_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\EEG analysis\Pilot2_Jan2026\full_processed_data'
raw_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\EEG analysis\Pilot2_Jan2026\raw_sets';
raw_sets_path = fullfile(study_filepath, 'raw_sets');
if ~exist(raw_sets_path, 'dir')
    mkdir(raw_sets_path);
end


% Designated folder where ALL fully processed .set files will be saved
% This folder will contain the final trimmed, filtered, re-referenced, ASR'd, ICA'd, and ICLabeled datasets.
output_folder_processed = fullfile(study_filepath, 'full_processed_data'); 
save_rejected_components_data_path = fullfile(output_folder_processed, "ICA rejection");

% --- Trimming Parameters (from your previous query) ---
% Define trigger event types to identify the end of relevant data
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


% --- Preprocessing Parameters ---
% Filtering (High-pass and Low-pass)
high_pass_filter_hz = 40; % Hz
low_pass_filter_hz = 1; % Hz

% Re-referencing: Cz channel index
cz_channel_index = 61; 


% ASR (Artifact Subspace Reconstruction) parameters
asr_flatline_criterion = 5; % Number of standard deviations for burst rejection
asr_channel_criterion = 0.8; % Max correlation with neighboring channels to be kept
asr_burst_criterion = 'off'; % 'off' means ASR handles bursts
asr_window_criterion = 'off'; % 'off' means ASR handles windows
asr_burst_rejection = 'off'; % 'off' means ASR handles rejection

% ICA (Independent Component Analysis) parameters
icatype     = 'runica'; 
ica_options = 'extended'; 
 
% --- Create output directory if it doesn't exist ---
if ~exist(output_folder_processed, 'dir')
    mkdir(output_folder_processed);
    fprintf('Created output directory for all processed data: %s\n', output_folder_processed);
end
 
% % --- Initialize EEGLAB and Load Study ---
% [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
 
iclabel_thresholds = struct(...
    'Brain', 0.8, ...   % Keep components where Brain prob is high (or reject if below certain threshold and other prob high)
    'Muscle', 0.8, ...  % Reject if Muscle probability is >= this threshold
    'Eye', 0.8, ...     % Reject if Eye probability is >= this threshold
    'Heart', 0.8, ...   % Reject if Heart probability is >= this threshold
    'LineNoise', 0.8, ... % Reject if Line Noise probability is >= this threshold
    'ChannelNoise', 0.8,... % Reject if Channel Noise probability is >= this threshold
    'Other', 0.8 ...    % Reject if Other probability is >= this threshold (often indicates problematic components)
);
% 


%%%%% SET OPTIONS FOR WHICH SCRIPT SEGMENTS TO RUN %%%%%%%%%%
load_data                   = 0; % load initial data into study
num_datasets                = 15;
preprocessing               = 1;
ICArejection                = 1;
epoching                    = 0;


subjects = dir(fullfile(data_path, 'Ox*'));
subjects = subjects([subjects.isdir]);
valid_subjects = [3:12,14:18];
% [STUDY, ALLEEG] = pop_loadstudy( ...
%     'filename', study_filename, ...
%     'filepath', study_filepath);
% 
% CURRENTSET = length(ALLEEG);
% existing_subjects = unique({STUDY.datasetinfo.subject});

if load_data == 1

    % existing_subjects = unique({STUDY.datasetinfo.subject});


    for i = 14:18
    
        subjID = subjects(i).name;
    
        % % ---- SKIP if already in STUDY ----
        % if ismember(subjID, existing_subjects)
        %     fprintf('Skipping %s (already in STUDY)\n', subjID);
        %     continue
        % end
    
        subjPath = fullfile(data_path, subjID);
        curryFile = dir(fullfile(subjPath, 'Acquisition*.dat'));
    
        if isempty(curryFile)
            warning('No Curry file found for %s — skipping', subjID);
            continue
        end
    
        fprintf('Appending NEW subject %s\n', subjID);
    
        EEG = loadcurry( ...
            fullfile(subjPath, curryFile(1).name), ...
            'KeepTriggerChannel', 'True', ...
            'CurryLocations', 'False');
    
        EEG.subject = subjID;
        EEG.setname = sprintf('%s_raw', subjID);
    
        EEG = pop_saveset(EEG, ...
            'filename', sprintf('%s_raw.set', subjID), ...
            'filepath', raw_sets_path);
    
        % ---- APPEND (critical line) ----
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET+1);
    
        % ---- Update STUDY metadata ----
        [STUDY, ALLEEG] = std_editset(STUDY, ALLEEG, ...
            'commands', {{ ...
                'index', CURRENTSET, ...
                'subject', subjID, ...
                'condition', 'continuous' }}, ...
            'updatedat', 'on');
    end
    
    [STUDY, ALLEEG] = pop_savestudy(STUDY, ALLEEG, ...
        'filename', 'Category_Switch_ALL.study', ...
        'filepath', study_filepath);

end
%%%%%%% INITIAL PREPROCESSING %%%%%%%%%%%%
if preprocessing == 1

    % % Subject(s) that used LM/RM reference
    % LMRM_subjects = {'Ox13'};
    % 
    % % Helper function
    % uses_LMRM = @(subj) ismember(subj, LMRM_subjects);


     % % === Load in study created by ERP_analysis_KH STUDY  ===
    % [STUDY, ALLEEG, CURRENTSET] = pop_loadstudy('filename', 'Category_Switch_ALL.study', 'filepath', study_filepath);
    
    for i = valid_subjects
    
        % EEG = ALLEEG(i);
        subj = EEG.subject;
    
        fprintf('\n--- Subject %s (%d/%d) ---\n', subj, i, numel(ALLEEG));
    
        % % ================= SKIP IF ALREADY PREPROCESSED =================
        % if isfield(EEG, 'etc') && isfield(EEG.etc, 'preprocessing_complete') ...
        %         && EEG.etc.preprocessing_complete == true
        %     fprintf('  ✔ Already preprocessed — skipping\n');
        %     continue;
        % end
    
        filename = sprintf('Ox%02d_raw.set', i);
        % Load dataset if not already loaded
        if isempty(EEG.data)
            EEG = pop_loadset('filename', filename, 'filepath',raw_path);
        end

    
        if isempty(EEG)
            fprintf(2, 'WARNING: Failed to load raw dataset %s. Skipping this dataset.\n', full_raw_path);
            continue;
        end
    
        if ischar(EEG.data) || isstring(EEG.data)
            fprintf('  Forcing EEG.data into memory to avoid .fdt errors...\n');
            EEG = eeg_checkset(EEG, 'loaddata');

        end

        % --- Step 1: Re-trim Data ---
        fprintf('  Step 1/6: Re-trimming data based on last trigger event...\n');
        event_latencies = [EEG.event.latency];
        last_trigger_latency = max(event_latencies) + 10000;
  
    
     
        end_time_sec = last_trigger_latency / EEG.srate;
            
        fprintf('  Last trigger event at sample: %d (%.3f seconds). Trimming data...\n', last_trigger_latency, end_time_sec);
        EEG = pop_select(EEG, 'time', [0 end_time_sec]);
        EEG.setname = [EEG.setname '_trimmed']; % Update name for clarity
    
        fprintf('  Trimming complete.\n');
    
        % --- Step 2: Filtering (High-pass and Low-pass) ---
        fprintf('  Step 2/6: Filtering data (High-pass: %.1f Hz, Low-pass: %.1f Hz)...\n', high_pass_filter_hz, low_pass_filter_hz);
        EEG = pop_eegfiltnew(EEG, 'locutoff', 1, 'hicutoff', [], 'plotfreqz', 1);
        EEG = pop_eegfiltnew(EEG, 'locutoff', [], 'hicutoff', 40, 'plotfreqz', 1);
        EEG.setname = [EEG.setname '_filtered']; 
        fprintf('  Filtering complete.\n');
    
        % ================= RE-REFERENCING =================
        % if uses_LMRM(subj)
        %     fprintf('  Using LM/RM reference for %s\n', subj);
        % 
        %     % Find LM and RM
        %     lm_labels = {'LM','M1','LMAST','LMast','MastL'};
        %     rm_labels = {'RM','M2','RMAST','RMast','MastR'};
        % 
        %     lm_idx = find(ismember(lower({EEG.chanlocs.labels}), lower(lm_labels)), 1);
        %     rm_idx = find(ismember(lower({EEG.chanlocs.labels}), lower(rm_labels)), 1);
        % 
        %     if isempty(lm_idx) || isempty(rm_idx)
        %         error('LM/RM channels not found for %s', subj);
        %     end
        % 
        %     EEG = pop_reref(EEG, [rm_idx], 'keepref','on');
        %     EEG.etc.reference = 'RM'; % 'LM/RM';
        % 
        % else
            fprintf('  Using Cz reference for %s\n', subj);
            EEG = pop_reref(EEG, cz_channel_index, 'keepref','on');
            EEG.etc.reference = 'Cz';
        % end

        % Always go to CAR after
        fprintf('  Re-referencing to Common Average Reference (CAR)...\n');
        EEG = pop_reref(EEG, []);
        EEG.etc.reference = [EEG.etc.reference ' -> CAR'];
        fprintf('  Re-referencing complete.\n');

    
        % --- Step 4: ASR ---
        EEG_preASR = EEG;
        
        EEG_clean = pop_clean_rawdata(EEG, ...
            'FlatlineCriterion', asr_flatline_criterion, ...
            'ChannelCriterion', asr_channel_criterion, ...
            'Highpass', 'off', ...
            'BurstCriterion', asr_burst_criterion, ...
            'WindowCriterion', asr_window_criterion, ...
            'BurstRejection', asr_burst_rejection, ...
            'Distance', 'Euclidian', ...
            'fusechannels', 0.5 ...
            );
        
        % Restore and interpolate
        EEG = restore_and_interpolate_after_ASR(EEG_preASR, EEG_clean);

        % --- Step 5: Run ICA (Independent Component Analysis) ---
        fprintf('  Step 5/6: Running ICA using %s algorithm...\n', icatype);
        EEG = pop_reref(EEG, []);
        if ~isfield(EEG, 'icaweights') || isempty(EEG.icaweights)
            %rankData = rank(double(getrank(double(eeg_getdatact(EEG)))));
            EEG = pop_runica(EEG, 'extended', 1);
        else
            fprintf('  ICA already present — skipping\n');
        end
    
        % ================= ICLabel =================
        fprintf('  Step 6/6: Classifying ICA components using ICLabel...\n');
        if ~isfield(EEG.etc,'ic_classification')
            fprintf('  Running ICLabel...\n');
            EEG = iclabel(EEG);
        else
            fprintf('  ICLabel already present — skipping\n');
        end

        fprintf('\n--- Automatic ICLabel component flagging ---\n');

        % Reset rejection flags
        EEG.reject.gcompreject = zeros(1, size(EEG.icaweights,1));
        
        % Extract ICLabel probabilities and class names
        ic_scores  = EEG.etc.ic_classification.ICLabel.classifications;
        ic_classes = EEG.etc.ic_classification.ICLabel.classes;
        
        nComp = size(ic_scores,1);
        
        for comp = 1:nComp
        
            % Initialize probabilities
            pBrain  = 0;
            pMuscle = 0;
            pEye    = 0;
            pHeart  = 0;
            pLine   = 0;
            pChan   = 0;
            pOther  = 0;
        
            % Extract probabilities directly
            for c = 1:length(ic_classes)
                switch ic_classes{c}
                    case 'Brain'
                        pBrain = ic_scores(comp,c);
                    case 'Muscle'
                        pMuscle = ic_scores(comp,c);
                    case 'Eye'
                        pEye = ic_scores(comp,c);
                    case 'Heart'
                        pHeart = ic_scores(comp,c);
                    case 'LineNoise'
                        pLine = ic_scores(comp,c);
                    case 'ChannelNoise'
                        pChan = ic_scores(comp,c);
                    case 'Other'
                        pOther = ic_scores(comp,c);
                end
            end
        
            % === Primary rejection rule ===
            if pMuscle >= 0.8 || ...
               pEye    >= 0.8 || ...
               pHeart  >= 0.8 || ...
               pLine   >= 0.8 || ...
               pChan   >= 0.8 || ...
               pOther  >= 0.8
        
                EEG.reject.gcompreject(comp) = 1;
                continue
            end
        
            % === Secondary safeguard rule ===
            % If Brain probability is very low AND artifact present
            if pBrain < 0.2 && (pMuscle > 0.2 || pEye > 0.2 || pChan > 0.2)
                EEG.reject.gcompreject(comp) = 1;
            end
        
        end
        
        flagged = find(EEG.reject.gcompreject == 1);
        fprintf('Auto-flagged %d components: %s\n', numel(flagged), mat2str(flagged));

        EEG.reject.icareject = EEG.reject.gcompreject;
        % Step 6: Manual IC inspection
        pop_viewprops(EEG, 0); % interactive
        pop_prop(EEG, 0, 1:size(EEG.icaweights,1), NaN, ...
    {'freqrange',[2 40], 'plotmode','condensed'});

        
        % Step 7: Remove rejected comps if any (assumes EEG.reject.gcompreject set)
        if isfield(EEG, 'reject') && isfield(EEG.reject, 'gcompreject')
            rej = find(EEG.reject.gcompreject==1);
            if ~isempty(rej)
                EEG = pop_subcomp(EEG, rej, 0);
                EEG.setname = [EEG.setname '_rejIC'];
            end
        end

            % ---- SAVE ----
        EEG = eeg_checkset(EEG,'loaddata');
        EEG = pop_saveset(EEG, ...
            'filename',[EEG.subject '_preprocessed.set'], ...
            'filepath',output_folder_processed, ...
            'savemode','onefile');
    
        EEG.etc.preprocessing_complete = true;
    
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, i);


        [STUDY, ALLEEG] = pop_savestudy(STUDY, ALLEEG, ...
            'filename', study_filename, ...
            'filepath', study_filepath);
        
        fprintf('\n✔ STUDY updated with preprocessed datasets\n');

       
    end
end

%%%%%%%%%%%%% ICA COMPONENT REJECTION (AUTO TAG + MANUAL) %%%%%%%%%%%%%%

%% ===== Replace your ICArejection == 1 block with the following robust workflow =====

if ICArejection == 1

    % Ensure save paths exist
    if ~exist(output_folder_processed, 'dir'), mkdir(output_folder_processed); end
    if ~exist(save_rejected_components_data_path, 'dir'), mkdir(save_rejected_components_data_path); end
    if ~exist(save_rejected_components_data_path, 'dir'), mkdir(save_rejected_components_data_path); end

    for i = valid_subjects


        EEG = ALLEEG(i); % take dataset as stored in ALLEEG
        current_filename = EEG.filename;
        current_subject = EEG.subject;
        fprintf('\n--- Dataset %d/%d : %s (subject: %s) ---\n', i, num_datasets, current_filename, current_subject);

            % ---------- (D) Auto-tag components via thresholds ----------
            % Prepare or reset reject flags
            EEG.reject.gcompreject = zeros(1, size(EEG.icaweights, 1));
            num_components_tagged = 0;

            if isfield(EEG, 'etc') && isfield(EEG.etc, 'ic_classification') && isfield(EEG.etc.ic_classification, 'ICLabel')
                ic_scores = EEG.etc.ic_classification.ICLabel.classifications;  % comps x classes
                ic_classes = EEG.etc.ic_classification.ICLabel.classes;

                % Iterate components and classes; tag if any class passes threshold (except 'Brain' unless you want to tag low-brain)
                for comp = 1:size(ic_scores,1)
                    % If Brain prob is low and one of artifact probs is high, tag
                    brainProb = ic_scores(comp, strcmp(ic_classes,'Brain'));
                    muscleProb = getprob(ic_scores, ic_classes, comp, 'Muscle');
                    eyeProb    = getprob(ic_scores, ic_classes, comp, 'Eye');
                    heartProb  = getprob(ic_scores, ic_classes, comp, 'Heart');
                    lineProb   = getprob(ic_scores, ic_classes, comp, 'LineNoise');
                    chanProb   = getprob(ic_scores, ic_classes, comp, 'ChannelNoise');
                    otherProb  = getprob(ic_scores, ic_classes, comp, 'Other');

                    % Tagging rules (adjust or add)
                    if muscleProb >= iclabel_thresholds.Muscle || ...
                       eyeProb    >= iclabel_thresholds.Eye    || ...
                       heartProb  >= iclabel_thresholds.Heart  || ...
                       lineProb   >= iclabel_thresholds.LineNoise || ...
                       chanProb   >= iclabel_thresholds.ChannelNoise || ...
                       otherProb  >= iclabel_thresholds.Other
                        EEG.reject.gcompreject(comp) = 1;
                        num_components_tagged = num_components_tagged + 1;
                    end

                    % Optional: tag components with very low brain probability
                    if brainProb < 0.05 && (muscleProb>0.05 || eyeProb>0.05 || chanProb>0.05)
                        EEG.reject.gcompreject(comp) = 1;
                        num_components_tagged = num_components_tagged + 1;
                    end
                end
                fprintf('Auto-tagged %d components for %s based on ICLabel thresholds.\n', num_components_tagged, current_filename);
            else
                fprintf('No ICLabel scores available for %s — skipping auto-tagging.\n', current_filename);
            end

            % Save a quick topoplot grid of auto-tagged components before manual inspection
            try
                if any(EEG.reject.gcompreject)
                    comps_to_plot = find(EEG.reject.gcompreject==1);
                    h = figure('visible','off');
                    pop_topoplot(EEG, 0, comps_to_plot, sprintf('Auto-tagged comps for %s', current_filename), 0, 'electrodes','on');
                    saveas(h, fullfile(save_rejected_components_data_path, [current_filename '_autoTagged_topos.png']));
                    close(h);
                end
            catch topoplot_err
                warning('Could not save auto-tagged topoplot for %s: %s', current_filename, topoplot_err.message);
            end

            % ---------- (E) Manual inspection GUI ----------
            fprintf('Opening manual inspection GUI for %s. Mark components to reject (use the GUI buttons) and close window when done.\n', current_filename);
            try
                % pop_viewprops opens GUI and updates EEG.reject.gcompreject after you interact.
                pop_viewprops(EEG, 0, 1); % show all comps, allow marking
            catch ME_view
                warning('pop_viewprops failed for %s: %s. Attempting pop_prop as fallback.\n', current_filename, ME_view.message);
                try
                    pop_prop(EEG, 0, 1);
                catch
                    warning('Fallback GUI also failed. Continuing without manual GUI for %s.\n', current_filename);
                end
            end

            % After GUI closes, EEG.reject.gcompreject should reflect manual selections.
            if ~isfield(EEG, 'reject') || ~isfield(EEG.reject, 'gcompreject')
                warning('No EEG.reject.gcompreject present after manual inspection for %s — skipping removal.', current_filename);
            else
                rejected_comp_indices = find(EEG.reject.gcompreject == 1);
                fprintf('User flagged %d components for rejection in %s: %s\n', numel(rejected_comp_indices), current_filename, mat2str(rejected_comp_indices));

                if ~isempty(rejected_comp_indices)
                    % Save info & scalp maps
                    subject_rejected_components_data = struct();
                    for k = 1:numel(rejected_comp_indices)
                        c = rejected_comp_indices(k);
                        subject_rejected_components_data(k).comp_idx = c;
                        if isfield(EEG, 'etc') && isfield(EEG.etc, 'ic_classification')
                            subject_rejected_components_data(k).iclabel_scores = EEG.etc.ic_classification.ICLabel.classifications(c,:);
                            subject_rejected_components_data(k).iclabel_classes = EEG.etc.ic_classification.ICLabel.classes;
                        else
                            subject_rejected_components_data(k).iclabel_scores = [];
                            subject_rejected_components_data(k).iclabel_classes = {};
                        end
                        subject_rejected_components_data(k).scalp_map = EEG.icawinv(:, c);
                        % activity saved as a short summary to avoid huge files
                        subject_rejected_components_data(k).act_mean = mean(EEG.icaact(c, :), 2);
                        subject_rejected_components_data(k).act_std  = std(EEG.icaact(c, :), [], 2);
                    end

                    % Save .mat describing rejected comps
                    [~, base_name_no_ext, ~] = fileparts(current_filename);
                    rejected_comps_data_filename = [base_name_no_ext '_rejected_components_data.mat'];
                    save(fullfile(save_rejected_components_data_path, rejected_comps_data_filename), 'subject_rejected_components_data');
                    fprintf('Saved rejected component data to %s\n', fullfile(save_rejected_components_data_path, rejected_comps_data_filename));

                    % Save topoplots for each rejected component
                    for k = 1:numel(rejected_comp_indices)
                        c = rejected_comp_indices(k);
                        fig = figure('visible','off');
                        topoplot(EEG.icawinv(:, c), EEG.chanlocs, 'electrodes','on');
                        title(sprintf('%s - Rejected Comp %d', base_name_no_ext, c));
                        figfile = sprintf('%s_Comp%03d_Topoplot.png', base_name_no_ext, c);
                        saveas(fig, fullfile(save_rejected_components_data_path, figfile));
                        close(fig);
                    end

                    % Actually remove the components from data
                    fprintf('Removing components %s from %s ...\n', mat2str(rejected_comp_indices), current_filename);
                    EEG = pop_subcomp(EEG, rejected_comp_indices, 0); % 0 = don't create new dataset in ALLEEG automatically
                    fprintf('Components removed and dataset updated.\n');

                    % Update setname to reflect rejection
                    EEG.setname = [EEG.setname '_rejIC'];
                else
                    fprintf('No components selected for rejection for %s — leaving dataset intact.\n', current_filename);
                end
            end

            % ---------- (F) Save updated dataset ----------
           
            saved_filename = [EEG.setname '.set'];
            EEG = pop_saveset(EEG, 'filename', saved_filename, 'filepath', processed_data_filepath);
            fprintf('Saved processed dataset to: %s\n', fullfile(processed_data_filepath, saved_filename));
            EEG.etc.ICArej_complete = true;

            % Put back into ALLEEG
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, i);

    end % loop datasets

    fprintf('\n--- Finished ICA Component Tagging and Manual Inspection for all datasets. ---\n');
    fprintf('Processed datasets saved to: %s\n', processed_data_filepath);

end % ICArejection

%% ===== helper function used above =====
function p = getprob(ic_scores, ic_classes, comp_idx, class_name)
    % small helper to safely return 0 if class not found
    loc = find(strcmp(ic_classes, class_name));
    if isempty(loc)
        p = 0;
    else
        p = ic_scores(comp_idx, loc);
    end
end


if epoching == 1
for i = 1:num_datasets
    
    %%%%%%%%%%%%%%%%%%%% --- Define Epoching Parameters --- %% in progress %%%%%%%%%%%%%%%%%%%%
    % Stimulus-locked epochs
    % Get dataset information for the current subject (points to ORIGINAL raw data)
    current_dataset_info = STUDY.datasetinfo(i);
    current_dataset_filename = current_dataset_info.filename;
    current_dataset_filepath = current_dataset_info.filepath;
    current_dataset_subject = current_dataset_info.subject;

    fprintf('\n--- Processing dataset %d/%d: %s (Subject: %s) ---\n', i, num_datasets, current_dataset_filename, current_dataset_subject);
    
    % Load the original raw dataset
    full_raw_path = fullfile(current_dataset_filepath, current_dataset_filename);
    fprintf('Loading raw dataset: %s\n', full_raw_path);
    EEG = pop_loadset('filename', current_dataset_filename, 'filepath', current_dataset_filepath);
    
    % Outcome-locked epochs
    outcome_epoch_correct = 10; 
    outcome_epoch_incorrect = 11;
    outcome_epoch_limits = [-0.5 1.5]; % e.g., -500 ms to 1500 ms relative to outcome
    outcome_baseline_limits = [-200 0]; % e.g., -500 ms to 0 ms relative to outcome onset

    stimulus_epoch_GO = [1,2,3,4];
    stimulus_epoch_NOGO = [5,6,7,8];
    stimulus_epoch_limits = [-0.2 1]; % e.g., -200 ms to 1000 ms relative to stimulus
    stimulus_baseline_limits = [-200 0]; % e.g., -200 ms to 0 ms relative to stimulus onset

    % % --- Epoching: Stimulus-locked ---
    fprintf('  Epoching: Creating stimulus-locked epochs...\n');
    EEG_stim_GO_epochs = pop_epoch( EEG, stimulus_epoch_GO, stimulus_epoch_limits, 'newname', [EEG.setname '_GOStimEpochs'], 'epochinfo', 'yes');
    EEG_stim_GO_epochs = pop_rmbase(EEG_stim_GO_epochs, stimulus_baseline_limits); % Baseline correct
    EEG_stim_NOGO_epochs = pop_epoch( EEG, stimulus_epoch_NOGO, stimulus_epoch_limits, 'newname', [EEG.setname '_NOGOStimEpochs'], 'epochinfo', 'yes');
    EEG_stim_NOGO_epochs = pop_rmbase(EEG_stim_NOGO_epochs, stimulus_baseline_limits); % Baseline correct

    fprintf('  Stimulus-locked epochs created. Number of epochs: %d\n', EEG_stim_epochs.trials);

    % --- Epoching: Outcome-locked ---
    fprintf('  Epoching: Creating outcome-locked epochs...\n');
    EEG_correct_epochs = pop_epoch( EEG, outcome_epoch_correct, outcome_epoch_limits, 'newname', [EEG.setname '_CorrectEpochs'], 'epochinfo', 'yes');
    EEG_correct_epochs = pop_rmbase( EEG_correct_epochs, outcome_baseline_limits); % Baseline correct
    
    EEG_incorrect_epochs = pop_epoch( EEG, outcome_epoch_incorrect, outcome_epoch_limits, 'newname', [EEG.setname '_IncorrectEpochs'], 'epochinfo', 'yes');
    EEG_incorrect_epochs = pop_rmbase( EEG_incorrect_epochs, outcome_baseline_limits); % Baseline correct


    fprintf('  Outcome-locked epochs created. Number of epochs: %d\n', EEG_outcome_epochs.trials);

    % --- Save the Epoched Datasets ---
    % It's good practice to save epoched data into a separate subfolder within your processed data.
    % Create specific subfolders for stimulus-locked and outcome-locked epochs if they don't exist
    stim_epochs_output_path = fullfile(study_filepath, 'Stimulus_Epochs');
    if ~exist(stim_epochs_output_path, 'dir'); mkdir(stim_epochs_output_path); end

    outcome_epochs_output_path = fullfile(study_filepath, 'Outcome_Epochs');
    if ~exist(outcome_epochs_output_path, 'dir'); mkdir(outcome_epochs_output_path); end

    % Save stimulus-locked epochs
    stim_epochs_filename = [EEG_stim_epochs.setname '.set'];
    pop_saveset(EEG_stim_GO_epochs, 'filename', stim_epochs_filename, 'filepath', stim_epochs_output_path);
    pop_saveset(EEG_stim_NOGO_epochs, 'filename', stim_epochs_filename, 'filepath', stim_epochs_output_path);

    fprintf('  Saved stimulus-locked epochs to: %s\n', fullfile(stim_epochs_output_path, stim_epochs_filename));

    % Save outcome-locked epochs
    outcome_epochs_filename = [EEG_outcome_epochs.setname '.set'];
    pop_saveset(EEG_correct_epochs, 'filename', outcome_epochs_filename, 'filepath', outcome_epochs_output_path);
    pop_saveset(EEG_incorrect_epochs, 'filename', outcome_epochs_filename, 'filepath', outcome_epochs_output_path);
    fprintf('  Saved outcome-locked epochs to: %s\n', fullfile(outcome_epochs_output_path, outcome_epochs_filename));

    % --- Save the fully processed CONTINUOUS dataset (as before) ---
    % Note: EEG variable here still holds the continuous data
    final_continuous_filename = [EEG.setname, '.set']; % Use the name of the continuous EEG
    full_output_path_final = fullfile(study_filepath, final_continuous_filename);

    fprintf('  Saving fully processed CONTINUOUS dataset to: %s\n', full_output_path_final);
    pop_saveset(EEG, 'filename', final_continuous_filename, 'filepath', study_filepath);

    % Update ALLEEG for the current dataset with its fully processed CONTINUOUS version
    ALLEEG(i) = EEG;


    % --- Save the fully processed dataset ---
    % The filename will reflect all applied steps 
    % (e.g., OriginalFilename_trimmed_filtered_CzRef_CAR_ASR_ICA_ICL.set)
    final_output_filename = [EEG.setname, '.set'];
    full_output_path = fullfile(output_folder_processed, final_output_filename);

    fprintf('  Saving fully processed dataset to: %s\n', full_output_path);
    pop_saveset(EEG, 'filename', final_output_filename, 'filepath', output_folder_processed);

    % Store the modified EEG back into ALLEEG and clear current EEG for memory management
    ALLEEG(i) = EEG;
    clear EEG;

end % End of loop through datasets
end
fprintf('\n--- All batch preprocessing steps complete! ---\n');

% --- Update and Save the EEGLAB STUDY ---
% This is crucial to make sure your STUDY links to the newly processed datasets.
new_study_filename = 'SalientModSwitch_FullyPreprocessed.study'; % New name for the updated study file
fprintf('Updating and saving STUDY file to: %s\n', fullfile(study_filepath, new_study_filename));
try
    [STUDY, ALLEEG] = pop_savestudy(STUDY, ALLEEG, 'filename', new_study_filename, 'filepath', study_filepath);
    fprintf('STUDY file updated and saved successfully.\n');
catch ME_savestudy
    fprintf(2, 'WARNING: Failed to save the updated STUDY file. Error: %s\n', ME_savestudy.message);
    fprintf('The processed datasets are saved individually, but the STUDY file linkage might not be updated.\n');
end

eeglab redraw; % Refresh EEGLAB GUI

function EEG_out = restore_and_interpolate_after_ASR(EEG_preASR, EEG_postASR)

    fprintf('\n--- Restoring ASR-removed channels (safe method) ---\n');

    orig_labels  = upper(strtrim({EEG_preASR.chanlocs.labels}));
    clean_labels = upper(strtrim({EEG_postASR.chanlocs.labels}));

    removed_labels = setdiff(orig_labels, clean_labels, 'stable');

    fprintf('ASR removed %d channels: %s\n', ...
        numel(removed_labels), strjoin(removed_labels, ', '));

    EEG_out = pop_interp(EEG_postASR, EEG_preASR.chanlocs, 'spherical');

    EEG_out = eeg_checkset(EEG_out);

    % Double-check no NaNs remain
    if any(isnan(EEG_out.data(:)))
        error('NaNs remain after interpolation. Check chanlocs or montage.');
    end

    EEG_out.etc.removed_by_ASR = removed_labels;
    EEG_out.etc.ASR_restored = true;

    fprintf('Channel restoration complete (safe mode).\n\n');
end
