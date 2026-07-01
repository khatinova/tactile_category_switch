% extract_reversal_aligned.m
% Extract 30 trials pre- and 30 trials post- reversal for each subject & block
% Stores results in subject's trial_data.* and also in all_trial_data

clearvars; close all;

%% --- Configure paths (edit to match your environment) ---
remote =2;
if remote == 1
    base_KH = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch/Salient mod switch KH/Data';
    base_RR = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch/Salient mod switch RR/Data/det_or_prob_and_conf';
    func_path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch/Functions';
elseif remote ==0
    base_KH = '\\Humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Data';
    base_RR = '\\Humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch RR\Data\det_or_prob_and_conf';
    func_path = '\\Humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Functions';
elseif remote == 2
    base_KH = 'Z:\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Data';
    base_RR = 'Z:\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch RR\Data\det_or_prob_and_conf';
    func_path = 'Z:\data\Projects\EEG_projects\Salient_Modality_Switch\Functions';
end

% Output save flag
save_per_subject = true;   % if true, writes trial_data_aligned.mat in each subject folder
outfilename = 'trial_data_aligned.mat';
addpath(func_path);

% Alignment window parameters
preN  = 30;  % trials before reversal (excluded reversal trial)
postN = 30;  % trials including reversal trial (reversal trial is first of post window)
alignedLen = preN +postN; % 60

%% --- Gather subject folders (KH first, then RR) ---
dKH = dir(base_KH);
dKH = dKH([dKH.isdir]);
dKH = dKH(~startsWith({dKH.name}, '.'));

dRR = dir(base_RR);
dRR = dRR([dRR.isdir]);
dRR = dRR(~startsWith({dRR.name}, '.'));

subjIDs = {'Nc01', 'Nc02', 'Nc03', 'Nc04', 'Nc05', 'Nc06', 'Nc07', 'Nc08','Nc09', 'Nc11', 'Nc11'};

block_struct_m = {'DDDD', 'DDDD','DDDD', 'DDDD', 'DDDD', 'DDDD', 'DDDD', 'DDDD', ...
    'VDPDP', 'VPDPD','VDPDP', 'VPDPD','VDPDP', 'VPDPD','VDPDP', 'VDDPP', 'VDDPP', 'PPDDP'};


%% stimuli presented in each block:
    % a: left right 4 vs 2 dots
    % b: top bottom two dots/square
    % c: left/right 3 vs 1 dots
    % d: left/right 3 vs 1 dost v2
    % e: L-shapes
    % f: p-shapes
stimuli = {'abcd', 'abcd','abcd', 'abcd','abcd', 'abcd','abcd', 'abcd', ...
    'aebfd','aebfd','aebcd','aebfd','eabfd','eabfd','eabfd','eabfd','abcde','abcde',...
    'abcd', 'abcd','abcd', 'abcd','abcd', 'abcd','abcd', 'abcd','abcd', 'abcd','abcd'};





% Build combined listing (preserve source)
subject_list = struct('name', {}, 'base', {}, 'source', {});
for i=1:numel(dKH)
    subject_list(end+1).name = dKH(i).name; %#ok<SAGROW>
    subject_list(end).base = base_KH;
    subject_list(end).source = 'KH';
    subject_list(end).block_struct = block_struct_m{i};
    subject_list(end).stimuli = stimuli{i};
end
n = numel(subject_list);
for i=1:numel(dRR)
    % k_idx = i + numel(dKH);
    subject_list(end+1).name = subjIDs{i};
    name = strsplit(dRR(i).name, '_');
    subject_list(end).block_struct = name{end}; %#ok<SAGROW>
    subject_list(end).base = base_RR;
    subject_list(end).source = 'RR';
    subject_list(end).stimuli = 'abcd';
end


fprintf('Found %d KH subjects and %d RR subjects (total %d)\n', numel(dKH), numel(dRR), numel(subject_list));

%% --- Main loop: load trial_data, compute aligned matrices, save into structures ---
all_trial_data = struct();  % will hold every subject's enriched trial_data
all_subj_perflines = []; % aggregate data for all subjetcs across blocks for plotting
group_T = table();

for s = 1:numel(subject_list)
    
    subj_name = subject_list(s).name;
    subj_base = subject_list(s).base;
    subj_src  = subject_list(s).source;
    subj_path = dir(fullfile(subj_base, [subj_name, '*']));
    fprintf('\n[%02d/%02d] Processing %s (%s)\n', s, numel(subject_list), subj_name, subj_src);

    tdpath = dir(fullfile(subj_base, [subj_name, '*']));
    % find trial_data file (try trial_data.mat first, else any P*.mat)
    f = fullfile(tdpath(1).folder, tdpath(1).name, 'trial_data.mat');
    % datafile = fullfile(tdpath(1).name, f(1).name);
    try
        loaded = load(f);
    catch ME
        warning('Could not load %s: %s — skipping', tdpath, ME.message);
        continue;
    end

    % Expect trial_data variable in the file (if not, try to construct)
    if isfield(loaded, 'trial_data')
        trial_data = loaded.trial_data;
    else
        beh_file = dir('Ox*.mat')
        load(beh_file(1).name)
    end

   
    if s < 9 
        % ony restructure the first two subjects which were strored trial_data(b,i).(var) rather than
        % trial_data.(var)(b,i)
        trial_data = restructure_trial_data(trial_data);
    end
     % Validate basic dimensions
    [nBlocks, nTrials] = size(trial_data.correct);
    
    
        
   raw_rows = {};
   % Process each block
   subj_perflines = []; % aggregate aligned performance for all blocks
   block_struct = subject_list(s).block_struct;
   stims = subject_list(s).stimuli;

   if isfield(trial_data, 'stimID') && height(trial_data.stimID) < 2
        trial_data.stimID = reshape(trial_data.stimID, 4,100);
   end

   choice_go = zeros(1,100);


    for b = 1:nBlocks

        % Aggregate block-level variables for the table
        nTrials = size(trial_data.goTrial, 2);

          
        % Preallocate aligned vectors (NaN-padded)
        aligned_correct = NaN(1, alignedLen);
        aligned_Hit     = NaN(1, alignedLen);
        aligned_FA      = NaN(1, alignedLen);
        aligned_Miss    = NaN(1, alignedLen);
        aligned_CR      = NaN(1, alignedLen);
        aligned_rt      = NaN(1, alignedLen);
        aligned_conf    = NaN(1, alignedLen);
        aligned_trials  = NaN(1, alignedLen);   % absolute trial indices (optional)
        aligned_perceivedCorrect  = NaN(1, alignedLen); 

        if b == 1 && nBlocks >= 6 
      
            block_struct = ['P', block_struct];
            stims = ['p', stims];

            % Store aligned results back into trial_data
            trial_data.aligned_correct(b,:)    = aligned_correct;
            trial_data.aligned_Hit(b,:)        = aligned_Hit;
            trial_data.aligned_FA(b,:)         = aligned_FA;
            trial_data.aligned_Miss(b,:)       = aligned_Miss;
            trial_data.aligned_CR(b,:)         = aligned_CR;
            trial_data.aligned_rt(b,:)         = aligned_rt;
            trial_data.aligned_confidence(b,:) = aligned_conf;
            trial_data.aligned_trials(b,:)     = aligned_trials;


            continue
        end

        % re-iterate block-level varaibles after accounting for practice
        % block

        stim_config = stims(b);
        curr_block_struct = block_struct(b);

        if b > 1
            prev_block_type = char(block_struct(b-1));
        else
            prev_block_type = 'NaN';
        end
        
        % get reversal trial for this block (should be scalar)
        % if isnan(trial_data.revTrial(b))
            % no reversal detected → leave NaNs and continue
            % fprintf('block %d: no revTrial found → calculating based on stimuli switch \n', b);
                
            % Get stimulus identity and correctness
        if subject_list(s).source(1) == 'K' % only compute revtrials for KH's data
            % --- PARAMETERS ---
            revMin = 30;
            revMax = 70;
            
            % --- Extract trial-wise vectors ---
            if isfield(trial_data, 'stimType')
                stimTypes = trial_data.stimType(b, :);
            else
                stimTypes = trial_data.stimID(b, :);
            end
            
            goTrial = logical(trial_data.goTrial(b, :));
            
            nTrials = numel(stimTypes);
            
            % Clamp window to available trials
            searchIdx = max(revMin,1) : min(revMax,nTrials);
            
            % --- Define flip stimulus sets ---
            A_flip = [2 3];
            B_flip = [1 4];
            
            revCandidates = [];
            
            % ---------- Mapping A (2/3 flip) ----------
            % Pre-reversal assumption:
            %   stim 2 = Go, stim 3 = NoGo
            idx_A_2_nogo = find( ...
                ismember(stimTypes(searchIdx), 2) & ...
                goTrial(searchIdx) == 0, ...
                1, 'first');
            
            idx_A_3_go = find( ...
                ismember(stimTypes(searchIdx), 3) & ...
                goTrial(searchIdx) == 1, ...
                1, 'first');
            
            if ~isempty(idx_A_2_nogo)
                revCandidates(end+1) = searchIdx(idx_A_2_nogo);
            end
            if ~isempty(idx_A_3_go)
                revCandidates(end+1) = searchIdx(idx_A_3_go);
            end
            
            % ---------- Mapping B (1/4 flip) ----------
            % Pre-reversal assumption:
            %   stim 1 = Go, stim 4 = NoGo
            idx_B_1_nogo = find( ...
                ismember(stimTypes(searchIdx), 1) & ...
                goTrial(searchIdx) == 0, ...
                1, 'first');
            
            idx_B_4_go = find( ...
                ismember(stimTypes(searchIdx), 4) & ...
                goTrial(searchIdx) == 1, ...
                1, 'first');
            
            if ~isempty(idx_B_1_nogo)
                revCandidates(end+1) = searchIdx(idx_B_1_nogo);
            end
            if ~isempty(idx_B_4_go)
                revCandidates(end+1) = searchIdx(idx_B_4_go);
            end
            
            % ---------- Final decision ----------
            if isempty(revCandidates)
                revTrial = NaN;
            else
                revTrial = min(revCandidates);   % earliest valid reversal
            end
            
            trial_data.revTrial(b) = revTrial;

        else 
        %else % if revTrial was stored in trial_data, isolate the current blocks' for calculating alignment
             revTrial = trial_data.revTrial(b);
        end

        correct = trial_data.correct(b,:);
        gonogo = trial_data.goTrial(b,:);

        choice_go(correct & gonogo | ~correct & ~gonogo) = 1;

        trial_data.choice_go(b,:) = choice_go;


        
        for t = 1:nTrials
        
            % Reversal state
            rev_state = double(isfinite(revTrial) && t >= revTrial);
            % Safe access to true feedback
            if isfield(trial_data,'pTrueFB') 
                if height(trial_data.pTrueFB)>1
                    trueFB = logical(trial_data.pTrueFB(b,t));
                    trueFB_block = logical(trial_data.pTrueFB(b,:));
                end
            elseif isfield(trial_data,'trueFB')
                trueFB = trial_data.trueFB(b,t);
                trueFB_block = logical(trial_data.trueFB(b,:));
            else
                trueFB = NaN;
            end
            
            % identify the choice on the previous time the same stimulus was seen
            curr_stim = stimTypes(t);
            curr_choice = choice_go(t);
            % find the previous time the same stimulus was shown BEFORE the
            % current trial
            prev_samestim_trial = find(stimTypes(1:t-1) == curr_stim, 1, 'last');

            if isempty(prev_samestim_trial)
                prevChoice = NaN;
                prevCorrect = NaN;
                stay_response = NaN;
                prevTrueFB = NaN;
            else
                prevChoice = choice_go(prev_samestim_trial);
                prevCorrect = trial_data.correct(b,prev_samestim_trial);
                prevTrueFB = trueFB_block(prev_samestim_trial);

                stay_response = prevChoice == choice_go(t);
            end

        

        
            raw_rows(end+1,:) = { ...
                subj_name, ...
                subj_src, ...
                b, ...
                t, ...
                trial_data.goTrial(b,t), ...
                trial_data.correct(b,t), ...
                trial_data.perceivedCorrect(b,t), ...
                trial_data.confidence(b,t), ...
                trial_data.rt(b,t), ...
                stimTypes(t), ... % the overall definity of stimID, unified across my first and 2nd data collection
                revTrial, ...
                rev_state, ...
                stim_config, ...
                curr_block_struct, ...
                trueFB, ...
                prev_block_type, ...
                prevChoice, ...
                prevCorrect, ...
                choice_go(t), ... % did they press the button?
                stay_response, ...
                prevTrueFB ...
            };
        end


  
        % define desired window indices (1-based)
        % ============================
        % REVERSAL-ALIGNED EXTRACTION
        % ============================
        
        % Relative alignment axis: [-30 ... -1 | 0 ... +29]
        rel_idx = -preN : (postN - 1);      % length = alignedLen (=60)        
        % Loop over reversal-aligned positions
        for w = 1:alignedLen
        
            % Convert relative index to absolute trial index
            abs_trial = revTrial + rel_idx(w);
        
            % Only copy data if trial exists
            if abs_trial >= 1 && abs_trial <= nTrials
        
                go  = trial_data.goTrial(b, abs_trial);
                cor = double(trial_data.correct(b, abs_trial));
                pcorr = double(trial_data.perceivedCorrect(b, abs_trial));
                

                aligned_trials(w)  = abs_trial;
                aligned_correct(w) = cor;
        
                aligned_Hit(w)  = (cor == 1 && go == 1);
                aligned_FA(w)   = (cor == 0 && go == 1);
                aligned_Miss(w) = (cor == 0 && go == 0);
                aligned_CR(w)   = (cor == 1 && go == 0);
        
                aligned_rt(w)   = trial_data.rt(b, abs_trial);
                aligned_perceivedCorrect(w) = pcorr;
        
                if isfield(trial_data, 'confidence')
                    aligned_conf(w) = trial_data.confidence(b, abs_trial);
                end
            end
        end
        
        % Store aligned results back into trial_data
        trial_data.aligned_correct(b,:)    = aligned_correct;
        trial_data.aligned_Hit(b,:)        = aligned_Hit;
        trial_data.aligned_FA(b,:)         = aligned_FA;
        trial_data.aligned_Miss(b,:)       = aligned_Miss;
        trial_data.aligned_CR(b,:)         = aligned_CR;
        trial_data.aligned_rt(b,:)         = aligned_rt;
      
        trial_data.aligned_trials(b,:)     = aligned_trials;
        trial_data.aligned_perceivedCorrect(b,:)    = aligned_perceivedCorrect;
        if isfield(trial_data, 'confidence')
            trial_data.aligned_confidence(b,:) = aligned_conf;
        end
        perfline = [movmean(trial_data.aligned_correct(1:preN), 20, 'omitmissing'), ...
                    movmean(trial_data.aligned_correct(preN+1:60), 20, 'omitmissing')];
        % if you want to plot individual block lines
        % subplot(5,5,s)
        % hold on
        % plot(1:60, perfline)    
        % ylim([0 1])
        % xlabel('Reversal-aligned trial')
        % ylabel('P(correct)')
        % title(subj_name)
        %legend({'block 1', 'block2', 'block3', 'block4'})

        subj_perflines = [subj_perflines; perfline];

    end


    subject_T = cell2table(raw_rows, ...
    'VariableNames', { ...
        'subjID', ...
        'researcher', ...
        'block', ...
        'trial', ...
        'goTrial', ...
        'correct', ...
        'perceivedCorrect', ...
        'confidence', ...
        'RT', ...
        'stimID', ...
        'revTrial', ...
        'rev_state', ...
        'stim_config', ...
        'block_type', ...
        'trueFB', ...
        'prev_block_type', ...
        'prevChoice', ...
        'prevCorrect', ...
        'choice',...
        'stay_choice', ...
        'prevTrueFB' ...
    });




    xline(preN+1,'--r');
    hold off
    
    correct_mean = mean(subj_perflines, 1, 'omitnan');
    correct_std  = std(subj_perflines, 0, 1, 'omitnan');
    correct_sem  = correct_std ./ sqrt(size(subj_perflines,1));


    subplot(3,10,s)
    hold on
    shadedErrorBar(1:60, correct_mean, correct_sem)    
    ylim([0 1])
    xlabel('Reversal-aligned trial')
    ylabel('P(correct)')
    title(subj_name)


    all_subj_perflines = [all_subj_perflines; subj_perflines];



    % Save back into all_trial_data

    safe_field = matlab.lang.makeValidName(subj_name);

    all_trial_data.(safe_field).trial_data = trial_data;
    all_trial_data.(safe_field).source = subj_src;
    all_trial_data.(safe_field).path   = subj_path;

    % optionally save per-subject aligned file
    if save_per_subject
        try
            save(fullfile(subj_path, outfilename), 'trial_data');
            fprintf('  Saved aligned trial_data to %s\n', fullfile(subj_path, outfilename));
        catch ME
            warning('  Could not save aligned file for %s: %s', subj_name, ME.message);
        end
    end

    group_T = [group_T; subject_T];
end

group_T.transition = strcat(string(group_T.prev_block_type), "→", string(group_T.block_type)) ;
group_T.transition(group_T.block == 1) = NaN;
group_T = fit_RW_subjectwise(group_T);



%% --- Post-processing summary & save all_trial_data ---
fprintf('\nDone processing %d subjects. Constructed all_trial_data with %d entries.\n', numel(subject_list), numel(fieldnames(all_trial_data)));

% Save the pooled all_trial_data in current folder (adjust path as needed)
try
    
    save(fullfile(base_KH, 'all_trial_data.mat'),'all_trial_data');
    save(fullfile(base_KH, 'behav_table.mat'),'group_T');
    fprintf('Saved all_trial_data_aligned.mat\n');
catch ME
    warning('Could not save all_trial_data_aligned.mat: %s', ME.message);
end


%% Plot the per subject correct performance

correct_mean = mean(all_subj_perflines, 1, 'omitnan');
correct_std  = std(all_subj_perflines, 0, 1, 'omitnan');
correct_sem  = correct_std ./ sqrt(size(all_subj_perflines,1));

hold on
shadedErrorBar(1:60, correct_mean, correct_sem)    
ylim([0 1])
xlabel('Reversal-aligned trial')
ylabel('P(correct)')
title('All subjects')

