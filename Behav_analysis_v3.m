clear

remote = 1; % am I working on remote desktop (1)?

if remote == 1
    path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch/Salient mod switch KH/Data';
else
    path = '\\Humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Data';
end

addpath(genpath(path))
cd(path)

% Initialize group structure
group_results = struct();

% Get list of subfolders (excluding '.' and '..')
folder_list = dir(path);
folder_list = folder_list([folder_list.isdir]);
folder_list = folder_list(~ismember({folder_list.name}, {'.', '..'}));

for subject = 1:length(folder_list)
    sbj_folder = fullfile(path, folder_list(subject).name);
    cd(sbj_folder)

    % Find .mat file starting with 'P'
    behav_datafile = dir('P*.mat');
    if isempty(behav_datafile)
        warning(['No behav file found in ', sbj_folder]);
        continue
    end

    filename = behav_datafile(1).name;
    load(filename);  % loads trial_data
    disp(['Loaded ', filename, ' for subject ', folder_list(subject).name]);

    % Make safe field name
    
    save_name = strjoin(["p", num2str(subject)], ""); % this will be the field name in the group results

    % Initialize structure
    subj_results = struct( ...
        'subject_ID', filename,...
        'pGo', NaN(4,1), ...
        'pCorrect', NaN(4,1), ...
        'pPerceivedCorrect', NaN(4,1),...
        'pTrueFB', NaN(4,1),...
        'confidence', NaN(4,100),...
        'reversal_trials', NaN(4,1),...
        'trial_data', trial_data);

    % Per-block trial processing
    for block = 1:height(trial_data)
        for trial = 1:100
            pGo{block}(trial)        = trial_data(block, trial).goTrial;
            pCorrect{block}(trial)   = trial_data(block, trial).correct;
            pPerceivedCorrect{block}(trial) = trial_data(block, trial).perceivedcorrect;
            pTrueFB{block}(trial)    = trial_data(block, trial).trueFB;
            confidence{block}(trial) = trial_data(block, trial).confidence;
            RT{block}(trial)         = trial_data(block, trial).rt;
            stimType{block}(trial)   = trial_data(block, trial).stimType;
        end

        % === Summary performance metrics ===
        subj_results.goTrial(block,:) =          pGo{block};
        subj_results.pGo(block)       =          mean(pGo{block}, 'omitnan');
        subj_results.correct(block,:) =          pCorrect{block};
        subj_results.pCorrect(block)  =          mean(pCorrect{block}, 'omitnan');
        subj_results.pPerceivedCorrect(block) =  mean(pPerceivedCorrect{block}, 'omitnan');
        subj_results.pTrueFB(block)   =          mean(pTrueFB{block}, 'omitnan');
        subj_results.confidence(block, :) =      confidence{block};
        subj_results.RT(block,:)      =          RT{block};
        subj_results.stimType(block,:)         = stimType{block};

      % --- Improved Reversal Detection ---
        % Get stimulus identity and correctness
        stimTypes = [trial_data(block,:).stimType];
        goTrial = [trial_data(block,:).goTrial];
        correctness = [trial_data(block,:).correct];
        
        % Find first Go trial for stim 3 (it becomes Go only after reversal)
        s1 = find(stimTypes == 1 & goTrial == 0, 1, 'first');   
        s2 = find(stimTypes == 2 & goTrial == 0, 1, 'first');
        s3 = find(stimTypes == 3 & goTrial == 1, 1, 'first');
        s4 = find(stimTypes == 4 & goTrial == 1, 1, 'first');

        % Collect all possible reversal trials into a cell array
        s = {s1, s2, s3, s4};
        
        % Filter out empty entries and concatenate the rest
        nonEmpty = [s{~cellfun(@isempty, s)}];
        
        % Find minimum
        if ~isempty(nonEmpty)
            inrevwindow = nonEmpty > 29;
            rt = min(nonEmpty(inrevwindow));
            subj_results.reversal_trials(block) = rt;
            fprintf('Reversal trial is: %d\n', rt);
        else
            disp('Reversal trial not detected.');
        end

        % === Trial type classification ===
        idx_correct = pCorrect{block} == 1;
        idx_incorrect = pCorrect{block} == 0;
        idx_go = pGo{block} == 1;
        idx_nogo = pGo{block} == 0;

        Hit  = idx_correct   & idx_go;
        FA   = idx_incorrect & idx_go;
        Miss = idx_incorrect & idx_nogo;
        CR   = idx_correct   & idx_nogo;
        
        % Store binary trial arrays in results struct
        subj_results.Hit(block,:)  = Hit;
        subj_results.FA(block,:)   = FA;
        subj_results.Miss(block,:) = Miss;
        subj_results.CR(block,:)   = CR;
        
        % Store total counts (mean also makes sense for proportions if needed)
        subj_results.Hit_total(block,1)  = sum(Hit, 'omitnan');
        subj_results.FA_total(block,1)   = sum(FA, 'omitnan');
        subj_results.Miss_total(block,1) = sum(Miss, 'omitnan');
        subj_results.CR_total(block,1)   = sum(CR, 'omitnan');
        
        % Optional: also store mean (i.e., proportion of trials)
        subj_results.Hit_mean(block,1)  = mean(Hit, 'omitnan');
        subj_results.FA_mean(block,1)   = mean(FA, 'omitnan');
        subj_results.Miss_mean(block,1) = mean(Miss, 'omitnan');
        subj_results.CR_mean(block,1)   = mean(CR, 'omitnan');
  
        % === Reaction Times ===
        subj_results.RT_correct{block,1}   = RT{block}(1, idx_correct);
        subj_results.RT_incorrect{block,1} = RT{block}(1, idx_incorrect);

        % === Confidence split by performance ===
        subj_results.confidence_correct{block}     = confidence{block}(idx_correct);
        subj_results.confidence_wrong{block}       = confidence{block}(idx_incorrect);
        subj_results.mean_confidence_correct(block,1) = mean(subj_results.confidence_correct{block}, 'omitnan');
        subj_results.mean_confidence_wrong(block,1)   = mean(subj_results.confidence_wrong{block}, 'omitnan');


        % ==== Align performance to reversal ====
        % Mark trials as pre/post reversal
        post_reversal = NaN(1, 100);  % preallocate
        if ~isnan(rt)
            post_reversal(1:rt-1) = 0;
            post_reversal(rt:100) = 1;
        end
        subj_results.post_reversal(block, :) = post_reversal;
        
        if ~isnan(rt)
            align_window = 30;
            mov_window = 10;
            num_aligned_points = 2 * align_window;  % 60 points: 30 before + 30 after
        
            % Initialize aligned matrices if not already done
            if ~isfield(subj_results, 'aligned_correct')
                subj_results.aligned_correct = NaN(4, num_aligned_points);
                subj_results.aligned_Hit = NaN(4, num_aligned_points);
                subj_results.aligned_Miss = NaN(4, num_aligned_points);
                subj_results.aligned_FA = NaN(4, num_aligned_points);
                subj_results.aligned_CR = NaN(4, num_aligned_points);
            end
        
            % Smooth correctness
            smoothed_correct = smoothn(pCorrect{block}, mov_window, 'nan');
        
            % Extract pre and post segments
            pre_start = max(1, rt - align_window);
            pre_data = smoothed_correct(pre_start:rt - 1);
            post_end = min(100, rt + align_window - 1);
            post_data = smoothed_correct(rt:post_end);
        
            % Pad pre and post to ensure length = 30 each
            pre_data = [NaN(1, align_window - length(pre_data)), pre_data];
            post_data = [post_data, NaN(1, align_window - length(post_data))];
        
            subj_results.aligned_correct(block,:) = [pre_data, post_data];
        
            % Align trial types
            trial_types = {'Hit', 'Miss', 'FA', 'CR'};
            for t = 1:length(trial_types)
                trial_name = trial_types{t};
                trial_vector = subj_results.(trial_name)(block, :);
        
                smoothed = smoothn(trial_vector, mov_window, 'nan');
                pre_data = smoothed(max(1, rt - align_window):rt - 1);
                post_data = smoothed(rt:min(rt + align_window - 1, 100));
        
                % Pad both sides
                pre_data = [NaN(1, align_window - length(pre_data)), pre_data];
                post_data = [post_data, NaN(1, align_window - length(post_data))];
                aligned_row = [pre_data, post_data];
        
                % Store
                switch trial_name
                    case 'Hit'
                        subj_results.aligned_Hit(block,:) = aligned_row;
                    case 'Miss'
                        subj_results.aligned_Miss(block,:) = aligned_row;
                    case 'FA'
                        subj_results.aligned_FA(block,:) = aligned_row;
                    case 'CR'
                        subj_results.aligned_CR(block,:) = aligned_row;
                end
            
            end

    
        end
      end
    
        % Save individual subject results
        if ~exist('processed_results', 'dir')
            mkdir('processed_results');
        end
        cd('processed_results')
        save(['performance_results_', filename(1:2), '.mat'], "subj_results")
        cd('..')
    
        % Store in group results
        group_results.(save_name) = subj_results;
end

cd('..');
cd('..');

cd("Behav_Results_salient_mod_switch\")
save('group_performance_results.mat', "group_results")
