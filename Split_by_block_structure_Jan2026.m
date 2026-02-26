% extract_reversal_aligned.m
% Extract 30 trials pre- and 30 trials post- reversal for each subject & block
% Stores results in subject's trial_data.* and also in all_trial_data

% clearvars; close all;

%% --- Configure paths (edit to match your environment) ---
remote =0;
if remote == 1
    base_path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch/Salient mod switch KH/Data';
    base_RR = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch/Salient mod switch RR/Data/det_or_prob_and_conf';
    func_path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch';
elseif remote ==0
    load_in_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\Behav results';
    fig_output_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\Behav results\Figures';
    func_path = '\\Humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Functions';
elseif remote == 2
    base_KH = 'Z:\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Data';
    base_RR = 'Z:\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch RR\Data\det_or_prob_and_conf';
    func_path = 'Z:\data\Projects\EEG_projects\Salient_Modality_Switch\Functions';
end

addpath(func_path);
cd(load_in_path)

% Alignment window parameters
preN  = 30;  % trials before reversal (excluded reversal trial)
postN = 30;  % trials including reversal trial (reversal trial is first of post window)
alignedLen = preN +postN; % 60

% for some reason this needs to be done manually...
load("\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\Behav results\all_trial_data_v2.mat")
subject_list = fields(all_trial_data);
block_struct_m = {'DDDD', 'DDDD','DDDD', 'DDDD', 'DDDD', 'DDDD', 'DDDD', 'DDDD', ...
    'VDPDP', 'VPDPD','VDPDP', 'VPDPD','VDPDP', 'VPDPD','VDPDP', 'VDDPP'};

%% --- Main loop: load trial_data, compute aligned matrices, save into structures ---
for s = 1:numel(subject_list)
    
    subj_name = subject_list{s};
    fprintf('\n[%02d/%02d] Processing %s\n', s, numel(subject_list), subj_name);


    trial_data = all_trial_data.(subj_name).trial_data;

     % Validate basic dimensions
    [nBlocks, nTrials] = size(trial_data.correct);

    % Process each block
    subj_perflines = []; % aggregate aligned performance for all blocks
    subj_conflines = [];
    subj_perflines_P = [];
    subj_perflines_D = [];
    subj_perflines_V = [];

    %% ----- Subject-level aggregate containers -----
    
    subj_stats = struct();
    
    % Trial-wise aggregates (collapsed across blocks)
    subj_stats.total_correct            = [];
    subj_stats.total_confidence         = [];
    subj_stats.total_conf_correct       = [];
    subj_stats.total_conf_incorrect     = [];
    subj_stats.total_correct_go         = [];
    subj_stats.total_correct_nogo       = [];
    
    % Aligned aggregates (per block, later averaged)
    subj_stats.aligned_correct          = [];
    subj_stats.aligned_confidence       = [];
    subj_stats.aligned_correct_P        = [];
    subj_stats.aligned_correct_D        = [];
    subj_stats.aligned_correct_V        = [];
    subj_stats.aligned_confidence_P     = [];
    subj_stats.aligned_confidence_D     = [];
    subj_stats.aligned_confidence_V     = [];
    
    % full blocks split by type
    subj_stats.correct_P        = [];
    subj_stats.correct_D        = [];
    subj_stats.correct_V        = [];
    subj_stats.confidence_P     = [];
    subj_stats.confidence_D     = [];
    subj_stats.confidence_V     = [];

    block_struct = [];
    if isfield(trial_data, 'block_structure')
        block_struct = trial_data.block_structure;
    else
        block_struct = block_struct_m{s};
    end
    
    for b = 1:nBlocks

        if b == 1 && nBlocks >= 6
            block_struct = [NaN, block_struct];
                continue
        end

 
        %% ----- Trial-wise metrics (NO SHRINKING) -----
        
        go_vec      = logical(trial_data.goTrial(b,:));
        correct_vec = double(trial_data.correct(b,:));
        
        if isfield(trial_data, 'confidence')
            conf_vec = trial_data.confidence(b,:);
        else
            conf_vec = NaN(size(correct_vec));
        end
        
        % ---- Always preserve trial count ----
        correct_vec = correct_vec(:);
        go_vec      = go_vec(:);
        conf_vec    = conf_vec(:);
        
        % ---- Overall ----
        subj_stats.total_correct    = [subj_stats.total_correct;    correct_vec'];
        subj_stats.total_confidence = [subj_stats.total_confidence; conf_vec'];
        
        % ---- Confidence split by correctness (NaN-masked) ----
        conf_correct = conf_vec;
        conf_correct(correct_vec ~= 1) = NaN;
        
        conf_incorrect = conf_vec;
        conf_incorrect(correct_vec ~= 0) = NaN;
        
        subj_stats.total_conf_correct   = [subj_stats.total_conf_correct;   conf_correct'];
        subj_stats.total_conf_incorrect = [subj_stats.total_conf_incorrect; conf_incorrect'];
        
        % ---- Accuracy split by Go / NoGo (NaN-masked) ----
        correct_go = correct_vec;
        correct_go(go_vec ~= 1) = NaN;
        
        correct_nogo = correct_vec;
        correct_nogo(go_vec ~= 0) = NaN;
        
        subj_stats.total_correct_go   = [subj_stats.total_correct_go;   correct_go'];
        subj_stats.total_correct_nogo = [subj_stats.total_correct_nogo; correct_nogo'];


        
        %% ----- Determine block types (P/D) and modality (Visual) -----
        % % Ensure correct length of block_structure
        % if length(block_struct) ~= nBlocks
        %     warning('%s: block_structure length mismatch (%d vs %d). Treating first block as practice and discarding it.', ...
        %             subj_name, length(block_struct), nBlocks);
        % 
        %     % Discard first (practice) block everywhere
        %     block_struct = block_struct(2:end);
        % 
        %     % Truncate or pad to match remaining blocks
        %     if length(block_struct) > nBlocks-1
        %         block_struct = block_struct(1:nBlocks-1);
        %     elseif length(block_struct) < nBlocks-1
        %         block_struct(end+1:nBlocks-1) = 'D';
        %     end
        % 
        %     % IMPORTANT: adjust nBlocks logic downstream
        %     % You are already skipping b==1, so this aligns perfectly
        % end
        
        % ---- Block type flags ----
        isVisualBlock = (block_struct == 'V');
        
        % Probabilistic if P OR V
        isProbBlock   = (block_struct == 'P');
        
        % Deterministic ONLY if D
        isDetBlock    = (block_struct == 'D');
        
        % ---- Sanity check ----
        if any(isVisualBlock & isDetBlock)
            error('%s: Found visual deterministic block — this should be impossible', subj_name);
        end


             % Number of trials in this block
        nTrials = size(trial_data.goTrial, 2);
        
        perfline = [movmean(trial_data.aligned_correct(1:preN), 5, 'omitmissing'), ...
                    movmean(trial_data.aligned_correct(preN+1:60), 5, 'omitmissing')];

        confline = [movmean(trial_data.aligned_confidence(1:preN), 5, 'omitmissing'), ...
            movmean(trial_data.aligned_confidence(preN+1:60), 5, 'omitmissing')];



        %% ----- Aligned aggregates -----

        
        if isProbBlock(b)
            subj_stats.aligned_correct_P    = [subj_stats.aligned_correct_P; trial_data.aligned_correct(b,:)];
            subj_stats.aligned_confidence_P = [subj_stats.aligned_confidence_P;trial_data.aligned_confidence(b,:)];
            subj_stats.correct_P            = [subj_stats.correct_P;  correct_vec'];
            subj_stats.confidence_P         = [subj_stats.confidence_P;  conf_vec'];

        elseif isDetBlock(b)
            trial_data.aligned_correct(b,:)
            subj_stats.aligned_correct_D    = [subj_stats.aligned_correct_D; trial_data.aligned_correct(b,:)];
            subj_stats.aligned_confidence_D = [subj_stats.aligned_confidence_D; trial_data.aligned_confidence(b,:)];
            subj_stats.correct_D            = [subj_stats.correct_D;  correct_vec'];
            subj_stats.confidence_D         = [subj_stats.confidence_D;  conf_vec'];
        
        elseif isVisualBlock(b)
            subj_stats.aligned_correct_V    = [subj_stats.aligned_correct_V; trial_data.aligned_correct(b,:)];
            subj_stats.aligned_confidence_V = [subj_stats.aligned_confidence_V; trial_data.aligned_confidence(b,:)];
            subj_stats.correct_V            = [subj_stats.correct_V;  correct_vec'];
            subj_stats.confidence_V         = [subj_stats.confidence_V;  conf_vec'];
        end

    end % end block loop

    %% ----- Subject-level summaries -----

    subj_summary = struct();
    
    % Trial-wise
    subj_summary.mean_correct          = mean(subj_stats.total_correct,1,'omitnan');
    subj_summary.mean_confidence       = mean(subj_stats.total_confidence,1,'omitnan');
    subj_summary.mean_conf_correct     = mean(subj_stats.total_conf_correct,1,'omitnan');
    subj_summary.mean_conf_incorrect   = mean(subj_stats.total_conf_incorrect,1,'omitnan');
    subj_summary.mean_correct_go       = mean(subj_stats.total_correct_go,1,'omitnan');
    subj_summary.mean_correct_nogo     = mean(subj_stats.total_correct_nogo,1,'omitnan');
    
    % Aligned
    subj_summary.aligned_correct_mean  = mean(subj_stats.aligned_correct,1,'omitnan');
    subj_summary.aligned_conf_mean     = mean(subj_stats.aligned_confidence,1,'omitnan');
    
    subj_summary.mean_aligned_correct_P     = mean(subj_stats.aligned_correct_P,1,'omitnan');
    subj_summary.mean_aligned_correct_D     = mean(subj_stats.aligned_correct_D,1,'omitnan');
    subj_summary.mean_aligned_correct_V     = mean(subj_stats.aligned_correct_V,1,'omitnan');
    
    subj_summary.mean_aligned_conf_P        = mean(subj_stats.aligned_confidence_P,1,'omitnan');
    subj_summary.mean_aligned_conf_D        = mean(subj_stats.aligned_confidence_D,1,'omitnan');
    subj_summary.mean_aligned_conf_V        = mean(subj_stats.aligned_confidence_V,1,'omitnan');


    all_trial_data.(subj_name).stats = subj_stats;
    all_trial_data.(subj_name).summary = subj_summary;

end
    
save('all_trial_data_PD_split.mat','all_trial_data')

%% --- Post-processing summary & save all_trial_data ---
fprintf('\nExtracted matrix including data from %d subjects \n', numel(subject_list));

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


correct_mean = mean(all_subj_perflines, 1, 'omitnan');
correct_std  = std(all_subj_perflines, 0, 1, 'omitnan');
correct_sem  = correct_std ./ sqrt(size(all_subj_perflines,1));

hold on
shadedErrorBar(1:60, correct_mean, correct_sem)    
ylim([0 1])
xlabel('Reversal-aligned trial')
ylabel('P(correct)')
title('All subjects')