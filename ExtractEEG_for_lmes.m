%% ============================================================
% Build trial-wise EEG + behavioural table and run mixed models
% Uses:
%   - all_trial_data (behaviour)
%   - outcome-locked EEG sets (Ox##_outcome.set)
% ============================================================

clear; clc;

%% ---------------- PATHS ----------------
behav_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\Behav results';
eeg_epoch_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\EEG analysis\Pilot2_Jan2026\Epoch analysis';

addpath(genpath(behav_path));

%% ---------------- LOAD BEHAVIOUR ----------------
load(fullfile(behav_path,'all_trial_data_v2.mat')); % -> all_trial_data

subjects = fieldnames(all_trial_data);
subjects = subjects(startsWith(subjects,'Ox'));

%% ---------------- EEG ROI DEFINITIONS ----------------
roiChan   = 'FCZ';

rewp_win  = [0.200 0.350];   % seconds
theta_win = [0.200 0.500];   % seconds
theta_f   = [4 8];           % Hz

%% ---------------- INITIALISE TABLE ----------------
group_T_EEG = table();

%% ================= MAIN LOOP =====================
for s = 1:numel(subjects)

    subj = subjects{s};
    fprintf('Processing %s\n', subj);

    td = all_trial_data.(subj).trial_data;

    nBlocks = size(td.correct,1);
    nTrials = size(td.correct,2);

    %% ---- behavioural vectors (flattened) ----
    correct     = td.correct(:);
    confidence  = td.confidence(:);
    rt          = td.rt(:);
    stimID      = td.stimID(:);
    goTrial     = td.goTrial(:);

    % block type: deterministic (1) vs probabilistic (0)
    blockType = NaN(nBlocks,1);
    for b = 1:nBlocks
        blockType(b) = mean(td.trueFB(b,:)) == 1;
    end
    blockType = repelem(blockType,nTrials);

    % previous block volatility
    prevBlockType = [NaN; blockType(1:end-nTrials)];

    %% ---- load EEG ----
    eegfile = fullfile(eeg_epoch_path, sprintf('%s_outcome.set',subj));
    if ~exist(eegfile,'file')
        warning('%s missing EEG file, skipping',subj);
        continue;
    end

    EEG = pop_loadset(eegfile);

    % channel index
    chanIdx = find(strcmpi({EEG.chanlocs.labels}, roiChan));
    if isempty(chanIdx)
        error('Channel %s not found in %s', roiChan, subj);
    end

    times = EEG.times/1000; % seconds
    rewp_idx  = times>=rewp_win(1)  & times<=rewp_win(2);
    theta_idx = times>=theta_win(1) & times<=theta_win(2);

    %% ---- single-trial EEG extraction ----
    nEEG = size(EEG.data,3);

    RewP_amp    = NaN(nEEG,1);
    theta_power = NaN(nEEG,1);

    for tr = 1:nEEG
        sig = squeeze(EEG.data(chanIdx,:,tr));

        % RewP mean amplitude
        RewP_amp(tr) = mean(sig(rewp_idx),'omitnan');

        % Theta power (Hilbert)
        theta_sig = bandpass(sig,theta_f,EEG.srate);
        theta_pow = abs(hilbert(theta_sig)).^2;
        theta_power(tr) = mean(theta_pow(theta_idx),'omitnan');
    end

    %% ---- align EEG trials to behavioural trials ----
    total_trials = nBlocks*nTrials;

    if nEEG > total_trials
        keep = (nEEG-total_trials+1):nEEG;
    else
        keep = 1:nEEG;
    end

    RewP_amp    = RewP_amp(keep);
    theta_power = theta_power(keep);

    %% ---- build participant table ----
    Tsub = table();
    Tsub.subject       = categorical(repmat({subj},numel(keep),1));
    Tsub.correct       = correct(1:numel(keep));
    Tsub.confidence    = confidence(1:numel(keep));
    Tsub.rt            = rt(1:numel(keep));
    Tsub.stimID        = stimID(1:numel(keep));
    Tsub.goTrial       = goTrial(1:numel(keep));
    Tsub.blockType     = categorical(blockType(1:numel(keep)),[0 1],{'Prob','Det'});
    Tsub.prevBlockType = categorical(prevBlockType(1:numel(keep)),[0 1],{'Prob','Det'});

    Tsub.RewP_amp      = RewP_amp;
    Tsub.theta_power   = theta_power;

    group_T_EEG = [group_T_EEG; Tsub];
end

%% ---------------- CLEAN ----------------
group_T_EEG = group_T_EEG(~isnan(group_T_EEG.prevBlockType),:);

group_T_EEG.logRT = log(group_T_EEG.rt);
group_T_EEG.confidence_z = NaN(height(group_T_EEG),1);

subs = categories(group_T_EEG.subject);
for s = 1:numel(subs)
    idx = group_T_EEG.subject==subs{s};
    c = group_T_EEG.confidence(idx);
    group_T_EEG.confidence_z(idx) = (c-mean(c,'omitnan'))./std(c,'omitnan');
end

save('group_T_EEG.mat','group_T_EEG','-v7.3');

%% ================= MODELS ======================

fprintf('\nRunning mixed models...\n');

% ---------------- Behaviour ----------------
glme_acc = fitglme(group_T_EEG, ...
    'correct ~ blockType*prevBlockType + (1|subject)', ...
    'Distribution','Binomial');

lme_conf = fitlme(group_T_EEG, ...
    'confidence_z ~ blockType*prevBlockType + (1|subject)');

% ---------------- EEG ----------------
lme_rewp = fitlme(group_T_EEG, ...
    'RewP_amp ~ blockType*prevBlockType + confidence_z + (1|subject)');

lme_theta = fitlme(group_T_EEG, ...
    'theta_power ~ blockType*prevBlockType + confidence_z + (1|subject)');

disp(glme_acc);
disp(lme_conf);
disp(lme_rewp);
disp(lme_theta);

fprintf('\n✔ EEG + behaviour group table built and models complete.\n');
