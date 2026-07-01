 %% Clear environment

 %% During troubleshooting changed/need to change back:
 % trial number
 % isvisual to always be visual
close all;
clear;

remote = 0; % are you working from the USB stick (1= yes, 0 = no)

if remote == 1
     addpath(genpath('E:\EEG Tasks\Salient_Modality_Switch'))
     addpath('F:\EEG Tasks\Salient_Modality_Switch\Functions')
     addpath('C:\Users\Experimenter\Desktop\Psychtoolbox')
else 
    task_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH';
    addpath('C:\Users\khatinova\OneDrive - Nexus365\Documents\MATLAB\Psychtoolbox\')
    addpath(genpath(task_path))
    cd(task_path)
    addpath('\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Functions');
end
%% Settings
subjID              = input('Enter subject ID: ', 's');
practice_setting    = 1; %input('Is this a practice task?');
EEGrecord           = 0; %input('record EEG?');
Metec_connected     = 0; %input('Is the Braille device connected?');
savee               = 1;  %input('save data?');
conf_settings       = true; % do you want confidence ratings in the REAL task?
date                = datetime('today','Format','ddMMyy');
structure_settings = 'PPDDP';
filename            = sprintf('%s_%s_%s', subjID, datestr(date,'ddmmyy'), structure_settings);

 % is the braille device connected?
KbName('UnifyKeyNames') % to enable the function to read the keyboard name
task_version = matlab.desktop.editor.getActiveFilename;
sca;

% load the Braille device
if Metec_connected == 1
    if remote == 1 
        asm = NET.addAssembly('F:\EEG Tasks\Braille device setup\MVBDLight (1).exe');
    else
        asm = NET.addAssembly('\\humerus\pharm_banerjee\data\Projects\EEG_projects\Braille device setup\MVBDLight (1).exe');
    end
    
    import Metec.MVBDLight.*;
    dev = BrailleDevice();
    
end

order = struct('PDPDP', [0.8, 1, 0.8, 1, 0.8],...
               'DPDPD', [1, 0.8, 1, 0.8, 1],...
               'DDDPP', [1, 1, 1, 0.8, 0.8],...
               'PDDPP', [0.8, 1, 1, 0.8, 1],...
               'PPDDP', [0.8, 0.8, 1, 1, 0.8]); 
stim_config = 'abcde' % choose which stimuli are used in each block, labelled in initialiseStim_v2
realNumBlocks = length(order.(structure_settings));  % e.g. 5
 % optionally going up to 6
hasPractice = practice_setting == 1;

if hasPractice
    numBlocks = realNumBlocks + 1;   % block 1 = practice
else
    numBlocks = realNumBlocks;
end

%% EEG Trigger Code Definitions


% Numeric codes for io64 / AntEEG (0-255)
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
    'confidence_high',   41, ...
    'confidence_low',    42, ...
    'response',          50, ...
    'response_go',       51, ...
    'response_nogo',     52, ...
    'response_miss',     53, ...
    'response_FA',       54 ...
);



nTrials.practice = 8;
nTrials.real    = 100; % 100; % divisible by 4
trialnumber     = nTrials.real;
revlim1         = 30;  
revlim2         = 70;  


% Generate reversal points
reversal_trial = round(rand(1, numBlocks) * (revlim2 - revlim1) + revlim1);

%% Initialize stimulus patterns
stimOFF     = uint8(0); % set the baseline off configuration
stimDuration = 1; % duration in seconds
triggerInterval = 0.001; % repeat trigger every 0.5s
nTriggers = floor(stimDuration / triggerInterval);

%% Build stim_GO_NO_GOmatrix from pre-specified order
stim_GO_NO_GOmatrix = cell(1, numBlocks);

% before block loop
pTrueFB_matrix = false(realNumBlocks, trialnumber);  % pre-allocate boolean matrix
rng('shuffle'); % or rng(subjID) for reproducible per-subject ordering


%% Pre-allocate additional timestamp fields (absolute times in seconds from expStart)

stim_blocks = initializeStimPatterns_v2();

for b = 1:realNumBlocks
    
    cfg = order.(structure_settings);

    shiftType = 'ED';                % 'ID' or 'ED'
    pTrueFB = cfg(b);                % e.g. 0.8 or 1
    isVisual(b) = false;

    if structure_settings(1) == 'V'
        isVisual(1) = true;          % logical
    end

    revTrial = reversal_trial(b);

    % generate stim order and go/nogo vector (1 = Go, 0 = NoGo)
    [stim_order, go_nogo] = generateStimGoNoGo(trialnumber, revTrial);

    % store as numeric 2 x trialnumber matrix (row1:stim, row2:go_nogo)
    stim_GO_NO_GOmatrix{b} = [stim_order; double(go_nogo)];

    
    nTrue = round(pTrueFB * trialnumber);    % number of true-feedback trials
    nFalse = trialnumber - nTrue;

    arr = [true(1, nTrue), false(1, nFalse)];   % vector with exact counts
    arr = arr(randperm(trialnumber));           % shuffle within block

    pTrueFB_matrix(b, :) = arr;                 % store row b

end

% Initialize EEG trigger setup
[portobject, portaddress, holdvalue, triggerLength] = initializeEEG(EEGrecord); % only initialises if EEGrecord is true

% %% Initialize sound
% [pahandle, freq, nrchannels] = initializeAudio();
% volume = 0.5;          % Volume from 0 to 1
% fbreak = 625;              % Frequency in Hz of the start of block
% beepLengthSecs2 = 0.3;       % Length of beep
% volumereward = 1;
% volumepunish = 1;
% rewardFreq = 1000;            % Hz, pleasant tone
% punishFreq = 250;            % Hz, unpleasant tone
% beepLengthSecs = 1;         % length of the outcome sound
% repetitions = 1;
% startCue = 1;
% waitForDeviceStart = 1;
% rewardAmp = 0.5;           % Base amplitude for 1000 Hz
% punishAmp = 5.62;          % ~15 dB louder for 250 Hz
% 
% % ---- Generate Beeps ----
% rewardBeep = rewardAmp * MakeBeep(rewardFreq, beepLengthSecs, freq);
% punishBeep = punishAmp * MakeBeep(punishFreq, beepLengthSecs, freq);
% blockStartBeep = MakeBeep(fbreak, beepLengthSecs2, freq);
% 
% % ---- Normalize to Prevent Clipping ----
% rewardBeep = rewardBeep / max(abs(rewardBeep));
% punishBeep = punishBeep / max(abs(punishBeep));
% blockStartBeep = blockStartBeep / max(abs(blockStartBeep));
% 
% % ---- Wrap for Stereo Playback (2 channels) ----
% rewardBeepStereo = [rewardBeep; rewardBeep];
% punishBeepStereo = [punishBeep; punishBeep];
% blockStartBeepStereo = [blockStartBeep; blockStartBeep];


%% Initialize screen
[window, windowRect, white, black, grey, xCenter, yCenter, customColors] = initializeScreen();

%% Fixation and response keys
fixation = '+';

% set the down key as the go key
keyGo = KbName('DownArrow');
% Get default keyboard device
kbdDevice = [];

KbQueueRelease;
KbQueueCreate(kbdDevice);
KbQueueStart(kbdDevice);
responseWindow = 1.0;   % seconds
outcome_window = 2; % seconds for outcome presentation


% optional metadata
trial_data.subjectID       = subjID;     % scalar
trial_data.date            = date;       % scalar or char
trial_data.block_order     = structure_settings;
trial_data.revTrial        = reversal_trial;


%% Block Loop
expStart = GetSecs;
for block = 1:numBlocks

    block_stim_config = stim_config(block);
    stim_block     = stim_blocks.(block_stim_config);

    if hasPractice && block <= realNumBlocks
        % ---- PRACTICE RUN THROUGH ALL REAL BLOCK TYPES ----
        isPracticeBlock = true;
        realBlockIdx = block;   % map practice block 1–5 to real structure 1–5
    else
        isPracticeBlock = false;
        realBlockIdx = block - realNumBlocks*hasPractice;
    end



    if isPracticeBlock
        % ===== PRACTICE BLOCK (DETERMINISTIC RULE) =====
        
        trialnumber = nTrials.practice;
    
        % Random stimulus order (but rule-based responses)
        stimIDs_block = randi(4, 1, trialnumber);
        
        % ---- DEFINE RULE ----
        if mod(block,2) == 1
            goStimIDs = [3 4];
        else
            goStimIDs = [1 2];
        end
        
        % Deterministic mapping
        goTrials_block = ismember(stimIDs_block, goStimIDs);
        
        trueFB_block = true(1, trialnumber);
        measure_confidence = true;
        
        % ===== SHOW VISUAL RULE SCREEN =====
        % showPracticeRuleScreen(window, stim_block, goStimIDs, black);

    else
        % REAL TASK
        trialnumber = nTrials.real;
        stimIDs_block  = stim_GO_NO_GOmatrix{realBlockIdx}(1,:);
        goTrials_block = logical(stim_GO_NO_GOmatrix{realBlockIdx}(2,:));
        trueFB_block   = pTrueFB_matrix(realBlockIdx, :);

        %isVisualBlock  = isVisual(realBlockIdx);
        measure_confidence = conf_settings;
        
        Screen('TextSize', window, 24);
        DrawFormattedText(window, ['Block ' num2str(block-1) ' - Press any key to begin'], 'center', 'center', black);
        Screen('Flip', window);

    end
    

    % send block_start trigger if requested
    if EEGrecord == 1
        sendTrigger(portobject, portaddress, triggers.block_start, triggerLength, holdvalue);
    end

    WaitSecs(2);

    Screen('TextSize', window, 24);
    DrawFormattedText(window, ['Block ' num2str(block-1) ' - Press any key to begin'], 'center', 'center', black);
    Screen('Flip', window);
    KbStrokeWait;

    Screen('TextSize', window, 100);
    DrawFormattedText(window, fixation, 'center', 'center', black);
    Screen('Flip', window);

    block_baseline = GetSecs;

    for trial = 1:trialnumber
        trial_baseline = GetSecs- expStart;

        % optional trial-level trial_start trigger
        if EEGrecord == 1
            sendTrigger(portobject, portaddress, triggers.trial_start, triggerLength, holdvalue);
        end


        stimID = stimIDs_block(trial);
        goTrial = goTrials_block(trial);
        feedbackIsTrue = trueFB_block(trial);


        % Present fixation
        Screen('TextSize', window, 100);
        DrawFormattedText(window, fixation, 'center', 'center', black);
        Screen('Flip', window);

        WaitSecs(0.5)

        % start stimulus presentation 
        stimonset_abs = GetSecs - expStart;  % absolute onset of stimulus

        % send stimulus trigger (if using EEG)
        if EEGrecord == 1
            if goTrial
                keyString = sprintf('stim_%d_GO', stimID);
            else
                keyString = sprintf('stim_%d_NOGO', stimID);
            end
            code = triggers.(keyString);  % numeric code
        end

        % show stimulus (either visual or tactile)
        % if isVisualBlock
        %     if Metec_connected == 1
        %         braille_code = stim_block(stimID);  
        %         showBrailleStim(window, xCenter, yCenter, braille_code, stimDuration, grey)
        %         Screen('Flip', window);
        % 
        %         if EEGrecord
        %            sendTrigger(portobject, portaddress, code, triggerLength, holdvalue);
        %         end
        %         for k = 1:nTriggers
        %             braille_stim = stim_block(stimID);
        %             dev.Send(braille_stim);
        %             WaitSecs(triggerInterval);
        %         end
        %         dev.Send(stimOFF);
        %         Screen('Flip', window);
        %         stimoffset_abs = GetSecs - expStart;    % absolute offset (end of stim presentation)
        %     else
        %         disp(['Stimulus ' num2str(stimID), ' Go? ' num2str(goTrial)]);
        %         WaitSecs(stimDuration);
        %     end
        % else
            if Metec_connected == 1

                for k = 1:nTriggers
                    if EEGrecord && k == 1
                       sendTrigger(portobject, portaddress, code, triggerLength, holdvalue);
                    end
                    braille_stim = stim_block(stimID);
                    dev.Send(braille_stim);
                    WaitSecs(triggerInterval);
                end
                dev.Send(stimOFF);
                stimoffset_abs = GetSecs - expStart;    % absolute offset (end of stim presentation)
            else
                disp(['Stimulus ' num2str(stimID), ' Go? ' num2str(goTrial)]);
                WaitSecs(stimDuration);
                braille_code = stim_block(stimID);  
                showBrailleStim(window, xCenter, yCenter, braille_code, stimDuration, grey)
            end
        % end
  
        Screen('Flip', window)
        
        % ================= RESPONSE WINDOW =================
            % WaitSecs(0.5)
        % Draw purple fixation (response cue)
        DrawFormattedText(window, fixation, 'center', 'center', customColors.green);
        vbl = Screen('Flip', window);      % exact response onset time
        tStart = GetSecs;

        % ---- defaults ----
        respKey = 'NoResponse';
        rt = NaN;
        respWasGo = false;
        respWasLateGo = false;
        
        lateRT = NaN;
        tEnd  = tStart + responseWindow;
        % tLate = tEnd   + lateWindow;

        % ---- MAIN RESPONSE WINDOW ----
        while GetSecs < tEnd
            [keyIsDown, secs, keyCode] = KbCheck;
            if keyIsDown && keyCode(keyGo)
                rt = secs - tStart;
                respKey = 'Go';
                respWasGo = true;
                if EEGrecord 
                    sendTrigger(portobject, portaddress, triggers.response_go, triggerLength, holdvalue);
                end
                break;
            end
        end
        % ---- EEG triggers ----
        if EEGrecord
            if ~respWasGo
                sendTrigger(portobject, portaddress, triggers.response_nogo, triggerLength, holdvalue);
            end
        
        end
        Screen('Flip', window)
        % % ---- RESPONSE WINDOW END (visual) ----
        % DrawFormattedText(window, fixation, 'center', 'center', black);
        % Screen('Flip', window);       

        % % ---- LATE RESPONSE WINDOW ----
        % if ~respWasGo
        %     while GetSecs < tLate
        %         [keyIsDown, secs, keyCode] = KbCheck;
        %         if keyIsDown && keyCode(keyGo)
        %             lateRT = secs - tStart;
        %             respWasLateGo = true;
        %             if EEGrecord 
        %                 sendTrigger(portobject, portaddress, triggers.response_late, triggerLength, holdvalue);
        %             end
        %             break;
        %         end
        %     end
        % end
        % 

        
        WaitSecs(0.5) 
        % ---- correctness ----
        if goTrial && respWasGo
            correct = 1;
        elseif ~goTrial && ~respWasGo
            correct = 1;
        else
            correct = 0;
        end
        
        % ---- perceived correctness ----
        if feedbackIsTrue
            perceivedCorrect = correct;
        else
            perceivedCorrect = ~correct;
        end
        
        % ---- flags ----
        respFlags.respWasGo     = respWasGo;
        respFlags.respWasMiss   = ~respWasGo && goTrial;
        respFlags.respWasFA     = respWasGo && ~goTrial;
        respFlags.respWasCR     = ~respWasGo && ~goTrial;
        % respFlags.respWasLateGo = respWasLateGo;
        % respFlags.lateRT        = lateRT;


        if EEGrecord
            if respFlags.respWasMiss
                sendTrigger(portobject, portaddress, triggers.response_miss, triggerLength, holdvalue);
            elseif respFlags.respWasFA
                sendTrigger(portobject, portaddress, triggers.response_FA, triggerLength, holdvalue);
            end
            % if respWasLateGo
            %     sendTrigger(portobject, portaddress, triggers.response_late, triggerLength, holdvalue);
            % end
        end
        
        % ---- absolute response time ----
        if ~isnan(rt)
            response_abs = tStart - expStart + rt;
        else
            response_abs = NaN;
        end


        fprintf('RT = %.3f | GoTrial = %d | respWasGo = %d| trueFB = %d\n' , ...
            rt, goTrial, respFlags.respWasGo, double(feedbackIsTrue));

        % after response, small wait and then confidence (if measured)
        % WaitSecs(0.5);

        confidence_onset_abs = NaN;
        confidenceRating = NaN;

        if measure_confidence
            % send eeg trigger for confidence rating if EEGrecord

            confidenceFlip = GetSecs;
            confidence_onset_abs = confidenceFlip - expStart;
            maxConfidenceTime = 4;
            if EEGrecord == 1
                sendTrigger(portobject, portaddress, triggers.confidence_rating, triggerLength, holdvalue);
            end
            [confidenceRating] = confidence_scale_function_v1(window, maxConfidenceTime);            
            
            
            if EEGrecord == 1
                
                if confidenceRating >= 5
                    code = 'confidence_high';
                else; code = 'confidence_low';
                end
                
                sendTrigger(portobject, portaddress, triggers.(code), triggerLength, holdvalue);
            end
        end
        WaitSecs(0.5+rand)



        % send correctness/feedback trigger
        if EEGrecord == 1
            if correct && feedbackIsTrue
                triggerID = 'correct_trueFB';
            elseif correct && ~feedbackIsTrue
                triggerID = 'correct_falseFB';
            elseif ~correct && feedbackIsTrue
                triggerID = 'incorrect_trueFB';
            elseif ~correct && ~feedbackIsTrue
                triggerID = 'incorrect_falseFB';
            end
            
        end
        %show outcome
        if EEGrecord
            sendTrigger(portobject, portaddress, triggers.(triggerID), triggerLength, holdvalue);
        end
        outcome_onset_abs = GetSecs - expStart;
        showOutcome(window, customColors, perceivedCorrect);

        
        WaitSecs(outcome_window)
        Screen('Flip', window);
        WaitSecs(0.5+rand); % small inter-trial interval



        % ---- STORE trial_data (block x trial matrices) ----
        trial_data.isPractice(block, trial)       = isPracticeBlock;
        trial_data.stimID(block, trial)           = stimID;
        trial_data.goTrial(block, trial)          = logical(goTrial);
        trial_data.respKey{block, trial}          = respKey;
        trial_data.rt(block, trial)               = rt;
        trial_data.correct(block, trial)          = double(correct);
        trial_data.perceivedCorrect(block, trial) = logical(perceivedCorrect);
        trial_data.trueFB(block, trial)           = logical(feedbackIsTrue);
        trial_data.stimonset(block, trial)        = stimonset_abs - expStart + expStart; % keep relative if you need; storing absolute already below
        trial_data.stimoffset(block, trial)       = stimoffset_abs - expStart + expStart;
        trial_data.stimonset_abs(block, trial)    = stimonset_abs;
        trial_data.stimoffset_abs(block, trial)   = stimoffset_abs;
        trial_data.response_abs(block, trial)     = response_abs;
        trial_data.outcome_onset_abs(block, trial)= outcome_onset_abs;
        trial_data.confidence_onset_abs(block, trial) = confidence_onset_abs;
        trial_data.confidence(block, trial)       = confidenceRating;
        trial_data.blocknum(block, trial)         = block;
        trial_data.trialnum(block, trial)         = trial;
%        trial_data.isVisual(block, trial)         = isVisualBlock;
        trial_data.respWasGo(block, trial)        = respFlags.respWasGo;
        trial_data.respWasMiss(block, trial)      = respFlags.respWasMiss;
        trial_data.respWasFA(block, trial)        = respFlags.respWasFA;
        trial_data.respWasCR(block, trial)        = respFlags.respWasCR;
        % trial_data.respLate(block, trial)         = respFlags.respWasLateGo;
        % trial_data.respLateRT(block, trial)       = respFlags.lateRT;

    end % trial loop

    if isPracticeBlock
        Screen('TextSize', window, 24);
        DrawFormattedText(window, ...
            'Practice complete.\n\nPress any key to start the main task.', ...
            'center', 'center', black);
        Screen('Flip', window);
        KbStrokeWait;
    else 

        Screen('TextSize', window, 24);
        DrawFormattedText(window, 'End of block!', 'center', 'center', black);
        Screen('Flip', window);
    end
    fprintf('Mean performance block %d: %0.3f\n', block, mean(trial_data.correct(block,:), 'omitnan'));

end % block loop


%% End of experiment
Screen('TextSize', window, 24);
DrawFormattedText(window, 'Thank you for participating!', 'center', 'center', black);
Screen('Flip', window);
WaitSecs(2);
sca;

if savee == 1
    cd('F:\EEG Tasks\Salient_Modality_Switch\Results')
    save(sprintf('%s.mat', filename), 'trial_data', 'task_version');
end


%% ====== AFTER RUN: build flattened events table for BIDS ======
% Convert 4xT fields into a 1xN trials long table (N = numBlocks*trialnumber)
nTrialsTotal = numBlocks * trialnumber;
rows = (1:nTrialsTotal)';

% flatten helper
flatten = @(M) reshape(M', [], 1); % reshape 4x100 -> 400x1 by row-major (block1 then block2 ...)

% Build table columns - use stim onset as the event onset (BIDS style)
onset = flatten(trial_data.stimonset_abs);        % seconds from expStart
duration = flatten(trial_data.stimoffset_abs - trial_data.stimonset_abs); % stim duration in seconds
trial_type = flatten(trial_data.goTrial);        % 1=Go, 0=NoGo (convert below)
stim_id = flatten(trial_data.stimID);
response = flatten(trial_data.respKey);          % cell
response_time = flatten(trial_data.rt);
accuracy = flatten(trial_data.correct);
perceivedCorrect = flatten(trial_data.perceivedCorrect);
trueFB = flatten(trial_data.trueFB);

blocknum_flat = flatten(trial_data.blocknum);
trialnum_flat = flatten(trial_data.trialnum);
%isVisual_flat = flatten(trial_data.isVisual);
confidence_flat = flatten(trial_data.confidence);
response_abs_flat = flatten(trial_data.response_abs);
outcome_onset_flat = flatten(trial_data.outcome_onset_abs);
conf_onset_flat = flatten(trial_data.confidence_onset_abs);
revTrial_per_block = repmat(trial_data.revTrial(:), trialnumber, 1); revTrial_flat = reshape(revTrial_per_block, [], 1);

% Normalize some columns
% trial_type to string labels for readability in events.tsv
trial_type_label = repmat({''}, nTrialsTotal, 1);
trial_type_label(trial_type==1) = {'go'};
trial_type_label(trial_type==0) = {'nogo'};
trial_type_label(isnan(trial_type)) = {'n/a'};

% Create table
events_table = table(onset, duration, trial_type_label, stim_id, response, response_time, accuracy, ...
    perceivedCorrect, trueFB,  blocknum_flat, trialnum_flat,  confidence_flat, ...
    response_abs_flat, outcome_onset_flat, conf_onset_flat, revTrial_flat, ...
    'VariableNames', {'onset','duration','trial_type','stim_id','response','response_time','accuracy',...
    'perceived_correct','true_feedback','block','trial_in_block','confidence',...
    'response_time_abs','outcome_onset_abs','confidence_onset_abs','revTrial'});

% Add BIDS-recommended required columns if missing (subject, session)
events_table.subject = repmat({sprintf('sub-%s', subjID)}, height(events_table), 1);
events_table.session = repmat({'ses-01'}, height(events_table), 1);

% Reorder columns to put onset,duration,trial_type first (BIDS style)
events_table = movevars(events_table, {'subject','session'}, 'Before', 1);
events_table = movevars(events_table, {'onset','duration','trial_type'}, 'After', 2);

% Write events.tsv
events_dir = fullfile(pwd, 'bids_events'); % change as desired
if ~exist(events_dir, 'dir'), mkdir(events_dir); end
events_tsv_name = fullfile(events_dir, sprintf('sub-%s_ses-01_task-CategorySwitch_events.tsv', subjID));
writetable(events_table, events_tsv_name, 'FileType','text','Delimiter','\t','QuoteStrings',false);

% Write a minimal JSON sidecar describing columns (so BIDS validators can read)
json_struct = struct();
json_struct.onset = 'seconds from experiment start (GetSecs)'; 
json_struct.duration = 'stimulus duration (seconds)';
json_struct.trial_type = 'go or nogo trial type';
json_struct.stim_id = 'stimulus identity (integer code)';
json_struct.response = 'response key name (or ''NoResponse'')';
json_struct.response_time = 'reaction time (seconds) relative to response window onset';
json_struct.accuracy = '1 = correct, 0 = incorrect';
json_struct.perceived_correct = 'subject perceived correctness (boolean)';
json_struct.true_feedback = 'whether feedback given was true (boolean)';
json_struct.pTrueFB = 'block-level probability that feedback is true for this trial';
json_struct.block = 'block number';
json_struct.trial_in_block = 'trial number within block';
%json_struct.is_visual = '1 = visual block, 0 = tactile block';
json_struct.confidence = 'confidence rating (scale whatever you used)';
json_struct.response_time_abs = 'absolute response time (seconds from expStart)';
json_struct.outcome_onset_abs = 'absolute outcome onset time (seconds from expStart)';
json_struct.confidence_onset_abs = 'absolute confidence onset (seconds from expStart)';
json_struct.revTrial = 'reversal trial index in that block';

jsonname = strrep(events_tsv_name,'.tsv','.json');
fid = fopen(jsonname,'w');
if fid ~= -1
    fprintf(fid, jsonencode(json_struct, PrettyPrint=true));
    fclose(fid);
else
    warning('Could not write events json sidecar at %s', jsonname);
end

% Save trial_data and events table for later analysis
save(sprintf('%s_trial_data_and_events.mat', filename), 'trial_data', 'events_table', 'expStart', '-v7.3');