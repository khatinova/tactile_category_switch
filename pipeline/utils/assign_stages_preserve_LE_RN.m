function T = assign_stages_preserve_LE_RN(T, block_col, trial_col, rev_trials_vec, stage_names)
% ASSIGN_STAGES_PRESERVE_LE_RN  Label each behavioural trial with a task stage.
%
% Canonical, single source of truth for stage assignment across the pipeline
% (previously duplicated in A_B_Preprocessing_KH.m, B_outcome_ERP_analysis_v9c.m
% and the MOCK script). All scripts should call THIS file.
%
% STAGES (20-trial windows, relative to the within-block reversal trial):
%   LN  Learning Naive   - first 20 trials of the block
%   LE  Learning Expert  - 20 trials immediately BEFORE the reversal
%   RN  Reversal Naive   - 20 trials immediately AFTER (and including) reversal
%   RE  Reversal Expert  - last 20 trials of the block
%
% CRITICAL RULE: LE and RN are protected. Where windows overlap (short blocks),
% the overlap is removed from LN/RE, never from LE/RN.
%
% INPUTS
%   T              : trial table
%   block_col      : name of the block-number column   (e.g. 'block')
%   trial_col      : name of the within-block trial col (e.g. 'trial')
%   rev_trials_vec : vector indexed by SEQUENTIAL block position (1 = first
%                    real block after any practice trimming). rev_trials_vec(b)
%                    is the within-block reversal trial number for block b.
%   stage_names    : {'LN','LE','RN','RE'}
%
% OUTPUT
%   T with added columns: stage (categorical), in_stage_window (logical),
%   stage_overlap_resolved (logical).

STAGE_LEN = 20;

if ~ismember(block_col, T.Properties.VariableNames)
    error('assign_stages_preserve_LE_RN: column "%s" not found.', block_col);
end
if ~ismember(trial_col, T.Properties.VariableNames)
    error('assign_stages_preserve_LE_RN: column "%s" not found.', trial_col);
end

T.stage = categorical(repmat(missing, height(T), 1), stage_names, 'Ordinal', false);
T.in_stage_window        = false(height(T), 1);
T.stage_overlap_resolved = false(height(T), 1);

block_nums = kh_to_numeric(T.(block_col));
trial_nums = kh_to_numeric(T.(trial_col));

if all(isnan(block_nums)); error('block column could not be made numeric.'); end
if all(isnan(trial_nums)); error('trial column could not be made numeric.'); end

unique_blocks = sort(unique(block_nums(~isnan(block_nums)))');

for bi = 1:numel(unique_blocks)
    raw_block  = unique_blocks(bi);
    block_rows = find(block_nums == raw_block);
    if isempty(block_rows); continue; end

    tib = trial_nums(block_rows);
    if ~any(~isnan(tib)); continue; end
    max_trial = max(tib, [], 'omitnan');
    if isnan(max_trial) || max_trial < 1; continue; end

    rev_trial = NaN;
    if ~isempty(rev_trials_vec) && bi <= numel(rev_trials_vec)
        rev_trial = rev_trials_vec(bi);
    end

    raw_ln = 1:min(STAGE_LEN, max_trial);
    raw_re = max(1, max_trial - STAGE_LEN + 1):max_trial;

    if isnan(rev_trial)
        ln_trials = raw_ln; le_trials = []; rn_trials = []; re_trials = raw_re;
    else
        le_start = max(1, rev_trial - STAGE_LEN);
        le_end   = max(1, rev_trial - 1);
        rn_start = rev_trial;
        rn_end   = min(max_trial, rev_trial + STAGE_LEN - 1);

        le_trials = le_start:le_end;
        rn_trials = rn_start:rn_end;
        protected = unique([le_trials, rn_trials]);

        ln_trials = setdiff(raw_ln, protected);
        re_trials = setdiff(raw_re, protected);

        overlap = unique([intersect(raw_ln, protected), intersect(raw_re, protected)]);
        if ~isempty(overlap)
            T.stage_overlap_resolved(block_rows(ismember(tib, overlap))) = true;
        end
    end

    T = set_stage(T, block_rows, tib, ln_trials, 'LN', stage_names);
    T = set_stage(T, block_rows, tib, le_trials, 'LE', stage_names);
    T = set_stage(T, block_rows, tib, rn_trials, 'RN', stage_names);
    T = set_stage(T, block_rows, tib, re_trials, 'RE', stage_names);
end
end

% -------------------------------------------------------------------------
function T = set_stage(T, block_rows, tib, trial_list, label, stage_names)
if isempty(trial_list); return; end
global_rows = block_rows(ismember(tib, trial_list));
if isempty(global_rows); return; end
T.stage(global_rows) = categorical({label}, stage_names);
T.in_stage_window(global_rows) = true;
end
