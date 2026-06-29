% =============================================================================
% S1_extract_behaviour.m  (was: A_extract_revaligned_alltrialdata.m)
%
% PIPELINE STEP 1 of 7. Builds the behavioural foundation for everything else.
%
% Extracts 30 trials pre- and 30 trials post-reversal for each subject & block,
% stores results in each subject's trial_data.* and the pooled all_trial_data,
% and assembles the long-format behavioural table group_T.
%
% INPUTS : per-subject trial_data.mat (KH "Ox*" and RR "Nc*" folders)
% OUTPUTS: all_trial_data.mat   (per-subject enriched trial_data struct)
%          behav_table.mat      (group_T: one row per behavioural trial)
%
% VARIABLE CONVENTIONS (see pipeline/PIPELINE.md): the long table uses
%   subj_id (string "Ox03"/"Nc07"), subj (numeric), cohort ("KH"/"RR"),
%   block, trial, block_type (D/P), correct, perceivedCorrect, confidence,
%   RT, revTrial, prev_block_type, transition, trueFB, ...
% Legacy column subjID is kept as an alias for backward compatibility.
% =============================================================================

clearvars; close all;

% Put pipeline utils on the path (figure style, subject-id, stage helpers)
addpath(genpath(fileparts(mfilename('fullpath'))));

%% --- Configure paths (edit to match your environment) ---
remote =0;
if remote == 1
    base_KH = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch/Salient mod switch KH/Data';
    base_RR = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch/Salient mod switch RR/Data';
    func_path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch/Functions';
    base_path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch/Salient mod switch KH'
elseif remote ==0
    base_KH = '\\Humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Data';
    base_RR = '\\Humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch RR\Data';
    func_path = '\\Humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Functions';
    base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH';
    fig_outpath = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\Behav results\Figures';
elseif remote == 2
    base_KH = 'Z:\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Data';
    base_RR = 'Z:\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch RR\Data';
    func_path = 'Z:\data\Projects\EEG_projects\Salient_Modality_Switch\Functions';
    base_path = 'Z:\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH';
end

% Output save flag
save_per_subject = true;   % if true, writes trial_data_aligned.mat in each subject folder
outfilename = 'trial_data_aligned.mat';
addpath(func_path);
cd(base_path)

% Alignment window parameters
preN  = 30;  % trials before reversal (excluded reversal trial)
postN = 30;  % trials including reversal trial (reversal trial is first of post window)
alignedLen = preN +postN; % 60
plot_mode = 'full'; % options: 'aligned', 'full', 'both'

%% --- Gather subject folders (KH first, then RR) ---
dKH = dir(base_KH);
dKH = dKH([dKH.isdir]);
dKH = dKH(~startsWith({dKH.name}, '.'));

dRR = dir(base_RR);
dRR = dRR([dRR.isdir]);
dRR = dRR(~startsWith({dRR.name}, '.'));



block_struct_m = {'DDDD', 'DDDD','DDDD', 'DDDD', 'DDDD', 'DDDD', 'DDDD', 'DDDD', ...
    'VDPDP', 'VPDPD','VDPDP', 'VPDPD','VDPDP', 'VPDPD','VDPDP', 'VDDPP', 'VDDPP',... 
    'PPDDP', 'PPDDP', 'DDPPD'};


%% stimuli presented in each block:
    % a: left right 4 vs 2 dots
    % b: top bottom two dots/square
    % c: left/right 3 vs 1 dots
    % d: left/right 3 vs 1 dost v2
    % e: L-shapes
    % f: p-shapes
stimuli = {'abcd', 'abcd','abcd', 'abcd','abcd', 'abcd','abcd', 'abcd', ...
    'aebfd','aebfd','aebcd','aebfd','eabfd','eabfd','eabfd','eabfd','abcde',...
    'abcde','abcde','abcde'}; % all Ox subjects end



% Build combined listing (preserve source)
subject_list = struct('name', {}, 'base', {}, 'source', {});
for i=1:numel(dKH)
    subject_list(end+1).name = dKH(i).name; %#ok<SAGROW>
    subject_list(end).base = base_KH;
    subject_list(end).source = 'KH';
    if i > numel(block_struct_m)
        subject_list(end).block_struct  = [];
    else
        subject_list(end).block_struct = block_struct_m{i};
    end
    if i > numel(stimuli)
        subject_list(end).stimuli = 'abcde'; % note that this needs to be changed
    else
        subject_list(end).stimuli = stimuli{i};
    end
end
n = numel(subject_list);
% In extract_revaligned_alltrialdata_v4.m, after the KH subject loop,
% add a second pass for RR subjects from ALL subfolders:

rr_subfolders = {'det_or_prob_and_conf', 'deterministic', 'det_to_prob', 'probabilistic'};


for sf = 1:numel(rr_subfolders)
    rr_path = fullfile(base_RR, rr_subfolders{sf});
    if ~exist(rr_path, 'dir'); continue; end

    rr_dirs = dir(fullfile(rr_path, 'Nc*'));
    rr_dirs = rr_dirs([rr_dirs.isdir]);


    for r = 1:numel(rr_dirs)
        nc_label = regexp(rr_dirs(r).name, 'Nc\d+', 'match', 'once');



        
        % Add to subject_list with source='RR'
        subject_list(end+1).name         = nc_label;
        subject_list(end).base           = rr_path;
        subject_list(end).source         = 'RR';
        file_components                  = strsplit(rr_dirs(r).name, '_');
        subject_list(end).block_struct   = file_components{end};  % parsed in loop
        file_components{end};
        subject_list(end).stimuli        = 'abcd';
        if strcmp( 'det_to_prob', rr_subfolders{sf})
            subject_list(end).block_struct = 'DDPP';
        elseif strcmp('probabilistic', rr_subfolders{sf})
            subject_list(end).block_struct = 'PPPP';
        end
    end
end


fprintf('Found %d KH subjects and %d RR subjects (total %d)\n', numel(dKH), numel(dRR), numel(subject_list));

%% --- Main loop: load trial_data, compute aligned matrices, save into structures ---
all_trial_data = struct();  % will hold every subject's enriched trial_data
all_subj_perflines = []; % aggregate data for all subjetcs across blocks for plotting
mean_subject_perflines = [];
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
    
    % ======================================================
    % REMOVE PRACTICE BLOCK (if present)
    % ======================================================
    if nBlocks == 6
        fprintf('Subject %s has 6 blocks → removing first (practice) block\n', subj_name);
    
        fn = fieldnames(trial_data);
    
        for f = 1:numel(fn)
            v = trial_data.(fn{f});
    
            % Only trim block-based variables (block x trial OR block x 1)
            if isnumeric(v) || islogical(v)
                if size(v,1) == 6
                    trial_data.(fn{f}) = v(2:end, :);
                end
            end
        end
    
        % Update block count
        [nBlocks, ~] = size(trial_data.correct);
    end
    
        
   raw_rows = {};
   % Process each block
   subj_perflines = []; % aggregate aligned performance for all blocks
   subj_fullblock_perflines = []; % aggregate non-aligned performance from all blocks
   if ~isempty(subject_list(s).block_struct)
    block_struct = subject_list(s).block_struct;
   else
       block_struct = trial_data.block_order;
   end
   stims = subject_list(s).stimuli;

   if isfield(trial_data, 'stimID') && height(trial_data.stimID) < 2
        trial_data.stimID = reshape(trial_data.stimID, 4,100);
   end

  
for b = 1:nBlocks
        % reset per block
        choice_go = zeros(1, size(trial_data.goTrial,2)); 

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
                
            %% Calculate reversal trial based off of first go/nogo violation
        correct = trial_data.correct(b,:);
        gonogo = trial_data.goTrial(b,:);

        choice_go(correct & gonogo | ~correct & ~gonogo) = 1;

        if isfield(trial_data,'stimType')
                stimTypes = trial_data.stimType(b,:);
        else
                stimTypes = trial_data.stimID(b,:);
        end

        trial_data.choice_go(b,:) = choice_go;
        detected_rev = detect_reversal_KH(stimTypes, gonogo, 30, 75);
        
        % Resolve revTrial based on subject source and ID:
        %   KH subjects Ox01–Ox16  → always use detected_rev
        %   KH subjects Ox17+      → use stored revTrial (fall back to detected if missing/NaN)
        %   RR subjects            → use stored revTrial (fall back to detected if missing/NaN)
        subj_num = str2double(regexp(subj_name, '\d+', 'match', 'once'));

        use_detected = strcmp(subj_src, 'KH') && subj_num <= 16;

        if use_detected
            revTrial = detected_rev;
        else
            % Prefer stored; fall back to detected if unavailable or NaN
            if isfield(trial_data, 'revTrial') && isfinite(trial_data.revTrial(b))
                revTrial = trial_data.revTrial(b);
            else
                warning('Block %d (%s, %s): no valid stored revTrial — falling back to detected (%d)', ...
                    b, subj_name, subj_src, detected_rev);
                revTrial = detected_rev;
            end
        end
        trial_data.revTrial(b) = revTrial;




        % fprintf('  block %d: detected revTrial = %s\n', b, mat2str(revTrial));
        
        % %% Quick visual sanity check for the block (around search window)
        % figure(1000 + s*10 + b); clf;
        % subplot(3,1,1);
        % plot(1:nTrials, goTrial, 'k.-'); hold on;
        % plot(searchIdx, goTrial(searchIdx), 'ro');
        % if isfinite(revTrial)
        %     xline(revTrial,'r--','LineWidth',1.5');
        % end
        % xlabel('trial'); ylabel('goTrial'); title(sprintf('%s block %d goTrial', subj_name, b));
        % 
        % subplot(3,1,2);
        % plot(1:nTrials, stimTypes, '.-'); hold on; ylim([min(stimTypes)-0.5 max(stimTypes)+0.5]);
        % if isfinite(revTrial), xline(revTrial,'r--'); end
        % xlabel('trial'); ylabel('stimID');
        % 
        % subplot(3,1,3);
        % plot(1:nTrials, trial_data.correct(b,:), '.-'); if isfinite(revTrial), xline(revTrial,'r--'); end
        % xlabel('trial'); ylabel('correct');
        % subtitle(sprintf('%s block %d: revTrial=%s', subj_name, b, mat2str(revTrial)));
        % drawnow;

        
        for t = 1:nTrials
        
            % Reversal state
            rev_state = double(isfinite(revTrial) && t >= revTrial);
            
            % Safe access to true feedback
            trueFB = NaN;
            % if isfield(trial_data,'pTrueFB') 
            %     if height(trial_data.pTrueFB)>1
            %         trueFB = logical(trial_data.pTrueFB(b,t))
            %         trueFB_block = logical(trial_data.pTrueFB(b,:));
            %     end
            % else
            if isfield(trial_data,'trueFB')
                trueFB = trial_data.trueFB(b,t);
                trueFB_block = logical(trial_data.trueFB(b,:));
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
            % 0th trial is the first trial following a detected reversal
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

        % --------- ALIGNED PERFORMANCE ---------
        perfline_aligned = [ ...
            movmean(trial_data.aligned_correct(b,1:preN), 20, 'omitmissing'), ...
            movmean(trial_data.aligned_correct(b,preN+1:end), 20, 'omitmissing')];
        
        % --------- FULL BLOCK PERFORMANCE ---------
        if isfinite(revTrial)
            left  = movmean(trial_data.correct(b,1:revTrial), 30, 'omitmissing');
            right = movmean(trial_data.correct(b,revTrial+1:end), 30, 'omitmissing');
            perfline_full = [left, right];
        else
            perfline_full = movmean(trial_data.correct(b,:), 30, 'omitmissing');
        end
        
        % store both
        subj_perflines = [subj_perflines; perfline_aligned];
        subj_fullblock_perflines = [subj_fullblock_perflines; perfline_full];

end
if s <19 % only add block struct for KH's participants, 
    % RR's participants should already have it from convert_RR_to_trial_data_v2
 trial_data.block_structure = block_struct;
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




        subplot(5,10,s)
        hold on
        
        switch plot_mode
        
            case 'aligned'
                % --- ALIGNED PLOT ---
                for r = 1:size(subj_perflines,1)
                    plot(-preN:(postN-1), subj_perflines(r,:), 'Color', [0.6 0.6 0.6])
                end
                xline(0,'--r','LineWidth',1.5)
                xlabel('Trials relative to reversal')
                ylabel('P(correct)')
                title([subj_name ' (aligned)'])
                ylim([0 1])
        
            case 'full'
                % --- FULL BLOCK PLOT ---
            
                cmap = lines(nBlocks);  % distinct colours per block
            
                for r = 1:size(subj_fullblock_perflines,1)
            
                    this_color = cmap(r,:);
            
                    % plot performance line
                    plot(1:length(subj_fullblock_perflines(r,:)), ...
                         subj_fullblock_perflines(r,:), ...
                         'Color', this_color, 'LineWidth', 1.5)
            
                    % plot matching reversal line
                    if isfinite(trial_data.revTrial(r))
                        xline(trial_data.revTrial(r), '--', ...
                            'Color', this_color, 'LineWidth', 1.5)
                    end
                end
            
                xlabel('Trial')
                ylabel('P(correct)')
                title([subj_name ' (full)'])
                ylim([0 1])
        
            case 'both'
                % --- BOTH OVERLAID (use carefully) ---
                subplot(1,2,1)
                hold on
                for r = 1:size(subj_perflines,1)
                    plot(-preN:(postN-1), subj_perflines(r,:), 'Color', [0.6 0.6 0.6])
                end
                xline(0,'--r')
                title('Aligned')
                ylim([0 1])
        
                subplot(1,2,2)
                hold on
                for r = 1:size(subj_fullblock_perflines,1)
                    plot(1:length(subj_fullblock_perflines(r,:)), subj_fullblock_perflines(r,:), 'Color', [0.6 0.6 0.6])
                end
                for b = 1:nBlocks
                    if isfinite(trial_data.revTrial(b))
                        xline(trial_data.revTrial(b),'--r')
                    end
                end
                title('Full block')
                ylim([0 1])
        end


    all_subj_perflines = [all_subj_perflines; subj_perflines];

    mean_subject_perflines = [mean_subject_perflines; mean(subj_perflines)];



    % Save back into all_trial_data

    safe_field = matlab.lang.makeValidName(subj_name);

    all_trial_data.(safe_field).trial_data = trial_data;
    all_trial_data.(safe_field).source = subj_src;
    all_trial_data.(safe_field).path   = subj_path;

    % optionally save per-subject aligned file
    % if save_per_subject
    %     try
    %         save(fullfile(subj_path.folder, sprintf(['%s_' outfilename], subj_name)), 'trial_data');
    %         fprintf('  Saved aligned trial_data to %s\n', fullfile(subj_path, outfilename));
    %     catch ME
    %         warning('  Could not save aligned file for %s: %s', subj_name, ME.message);
    %     end
    % end

    group_T = [group_T; subject_T];
end

exportgraphics(gcf, fullfile(fig_outpath,'Full performance all subjects.pdf'),'ContentType','vector')

group_T.transition = strcat(string(group_T.prev_block_type), "→", string(group_T.block_type)) ;
group_T.transition(group_T.block == 1) = NaN;

% ------------------------------------------------------------------
% UNIFY SUBJECT NAMING (canonical subj_id / subj / cohort).
% group_T currently has subjID (string) and researcher ('KH'/'RR').
% kh_subject_id maps these onto the canonical trio used everywhere
% downstream, keeping subjID as a backward-compatible alias.
% ------------------------------------------------------------------
if ismember('researcher', group_T.Properties.VariableNames) && ...
        ~ismember('cohort', group_T.Properties.VariableNames)
    group_T.cohort = string(group_T.researcher);
end
group_T = kh_subject_id('standardise', group_T);

% Normalise legacy probabilistic-visual blocks: 'V' -> 'P'
group_T.block_type = string(group_T.block_type);
group_T.block_type(group_T.block_type == "V") = "P";

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
figure('Position', [50 50 300 300])
set(gca, 'TickDir', 'out')
correct_mean = mean(mean_subject_perflines, 1, 'omitnan');
correct_std  = std(mean_subject_perflines, 0, 1, 'omitnan');
correct_sem  = correct_std ./ sqrt(size(mean_subject_perflines,1));

hold on
plot(1:60, mean_subject_perflines, 'LineWidth',0.3, 'Color', 'k')
shadedErrorBar(1:60, correct_mean, correct_sem)    
ylim([0 1])
xlabel('Reversal-aligned trial')
ylabel('P(correct)')
title('All subjects')


exportgraphics(gcf, fullfile(fig_outpath,'Revaligned performance all blocks.pdf'), 'ContentType', 'vector')


function revTrial = detect_reversal_KH(stimTypes, goTrial, revMin, revMax)
% DETECT_REVERSAL_KH
% Detect the first trial in [revMin, revMax] that is consistent with
% either of the two valid reversal configurations:
%   Config A: stim 3 → Go AND stim 2 → NoGo  (stim 3/2 swap)
%   Config B: stim 4 → Go AND stim 1 → NoGo  (stim 4/1 swap)
%
% First infers which config applies from the pre-reversal mapping,
% then detects the first post-reversal observation in that config.
%
% Returns the index of the detected first post-reversal trial,
% or NaN if detection fails.

if nargin < 3; revMin = 30; end
if nargin < 4; revMax = 75; end

nTrials = numel(stimTypes);
searchIdx = revMin : min(revMax, nTrials);

% --- Infer pre-reversal config from first 20 trials ---
% Which of stim 1/2 were Go before reversal?
pre_idx = 1:min(revMin-1, nTrials);
pre_go_stims = unique(stimTypes(pre_idx & goTrial(pre_idx)));

% Config A pre-reversal: stim 1 and 2 are Go → stim 3 and 4 are NoGo
% Config B pre-reversal: stim 1 and 4 are Go → stim 2 and 3 are NoGo
% (or similar — adjust based on your actual generateStimGoNoGo logic)

% Detect post-reversal evidence:
% Config A indicator: stim 3 appears as Go OR stim 2 appears as NoGo
cond_A = (ismember(stimTypes(searchIdx), 3) & goTrial(searchIdx)) | ...
         (ismember(stimTypes(searchIdx), 2) & ~goTrial(searchIdx));

% Config B indicator: stim 4 appears as Go OR stim 1 appears as NoGo  
cond_B = (ismember(stimTypes(searchIdx), 4) & goTrial(searchIdx)) | ...
         (ismember(stimTypes(searchIdx), 1) & ~goTrial(searchIdx));

% Either config counts as a valid reversal detection
cond_any = cond_A | cond_B;

revRel = find(cond_any, 1, 'first');

if ~isempty(revRel)
    revTrial = searchIdx(revRel);
    
    % Sanity check: detected trial must be within [revMin, 75]
    if revTrial > 75
        warning('Detected reversal at trial %d which exceeds revMax=75 — check data', revTrial);
        revTrial = NaN;
    end
else
    warning('No reversal detected in [%d, %d] — check stimTypes and goTrial', revMin, revMax);
    revTrial = NaN;
end
end