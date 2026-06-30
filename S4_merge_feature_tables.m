% =============================================================================
% S4_merge_feature_tables.m  (was: merge_group_tables_KH_RR.m)
%
% PIPELINE STEP 4 of 7 — merge the KH and RR per-trial EEG feature tables into
% one combined table for the RQ analysis (S7).
%
% INPUTS
%   KH: group_stage_table_features.mat  (group_table)   <- S3 run with COHORT KH
%   RR: group_stage_table_features_RR.mat (group_table) <- S3 run with COHORT RR
%       (legacy fallback: group_stage_table_RR_v7.mat)
% OUTPUT
%   group_feature_table_combined.mat (group_table)  -> consumed by S7
%
% Subject naming is unified to subj_id ("Ox03"/"Nc07") / subj / cohort via
% pipeline/utils/kh_subject_id.m. Study-metadata columns (feedback_modality,
% stimulus_modality, practice_task) are attached per cohort.
% =============================================================================

clear; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

% -------------------------------------------------------------------------
%% PATHS
% -------------------------------------------------------------------------
remote = 0;
if remote == 1
    base_path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch';
else
    base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch';
end

KH_epoch_folder = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis', 'Outcome_feature_tables_v4_merged');
RR_epoch_folder = fullfile(base_path, 'Salient mod switch RR', 'Results', 'EEG analysis',  'Outcome_feature_tables_v4_merged');
out_folder      = KH_epoch_folder;

% Canonical per-trial feature list to z-score within subject.
features_to_zscore = {'prefrontal_neg_peak_amp','prefrontal_neg_peak_norm','prefrontal_mean_amp', ...
                      'P300_amp','P300_norm','Theta_amp', ...
                      'PLV_fp','PLV_fs','PLV_fp_pairwise','PLV_fs_pairwise'};

% -------------------------------------------------------------------------
%% LOAD + PREP EACH COHORT
% -------------------------------------------------------------------------
group_table_KH = load(fullfile(KH_epoch_folder, 'group_stage_table_features.mat'));
group_table_RR = load(fullfile(KH_epoch_folder, 'group_stage_table_features_RR.mat'));

group_table_KH = group_table_KH.group_table;
group_table_RR = group_table_RR.group_table;


% Force the correct cohort label, then unify subject naming for each cohort.
group_table_KH.cohort = repmat("KH", height(group_table_KH), 1);
group_table_RR.cohort = repmat("RR", height(group_table_RR), 1);
group_table_KH = kh_subject_id('standardise', group_table_KH);
group_table_RR = kh_subject_id('standardise', group_table_RR);

% Within-subject z-scores (per cohort, using subj_id).
group_table_KH = zscore_within_subject(group_table_KH, features_to_zscore);
group_table_RR = zscore_within_subject(group_table_RR, features_to_zscore);

% -------------------------------------------------------------------------
%% STUDY-METADATA COLUMNS (per cohort)
% -------------------------------------------------------------------------
% KH: subjects <=8 had auditory feedback; 9-17 had a visual+tactile first block.
group_table_KH.feedback_modality = repmat("visual", height(group_table_KH), 1);
group_table_KH.feedback_modality(double(group_table_KH.subj) <= 8) = "auditory";

group_table_KH.stimulus_modality = repmat("tactile", height(group_table_KH), 1);
block_idx_KH = get_block_index(group_table_KH);
vt_mask = double(group_table_KH.subj) >= 9 & double(group_table_KH.subj) <= 17 & block_idx_KH == 1;
group_table_KH.stimulus_modality(vt_mask) = "visual+tactile";

group_table_KH.practice_task = strings(height(group_table_KH), 1);
group_table_KH.practice_task(double(group_table_KH.subj) <= 8)                                  = "matching";
group_table_KH.practice_task(double(group_table_KH.subj) >= 9 & double(group_table_KH.subj) <= 16) = "structured";
group_table_KH.practice_task(double(group_table_KH.subj) >= 17)                                 = "criterion";

group_table_RR.feedback_modality = repmat("visual",    height(group_table_RR), 1);
group_table_RR.stimulus_modality = repmat("tactile",   height(group_table_RR), 1);
group_table_RR.practice_task     = repmat("criterion", height(group_table_RR), 1);

% -------------------------------------------------------------------------
%% MERGE (align to KH column set) + SAVE
% -------------------------------------------------------------------------
vars = group_table_KH.Properties.VariableNames;
missing_in_RR = setdiff(vars, group_table_RR.Properties.VariableNames, 'stable');
for m = 1:numel(missing_in_RR)
    group_table_RR.(missing_in_RR{m}) = nan(height(group_table_RR), 1);
end
group_table_RR = group_table_RR(:, vars);

group_table = [group_table_KH; group_table_RR];

% Categorical key columns
group_table.subj_id           = categorical(group_table.subj_id);
group_table.cohort            = categorical(group_table.cohort);
group_table.feedback_modality = categorical(group_table.feedback_modality);
group_table.stimulus_modality = categorical(group_table.stimulus_modality);
group_table.practice_task     = categorical(group_table.practice_task);
group_table = movevars(group_table, 'subj_id', 'Before', 1);

save(fullfile(out_folder, 'group_feature_table_combined.mat'), 'group_table');
fprintf('Saved combined table: %d rows, %d cols (%d KH + %d RR subjects).\n', ...
    height(group_table), width(group_table), ...
    numel(unique(group_table_KH.subj_id)), numel(unique(group_table_RR.subj_id)));

% -------------------------------------------------------------------------
%% (Optional) merge the per-stage FRN/RewP difference-wave tables too
% -------------------------------------------------------------------------
frn_KH = try_load_frn(fullfile(KH_epoch_folder, 'frn_rewp_by_stage_KH.mat'));
frn_RR = try_load_frn(fullfile(RR_epoch_folder, 'frn_rewp_by_stage_RR.mat'));
if ~isempty(frn_KH) || ~isempty(frn_RR)
    frn_rewp_stage_table = [frn_KH; frn_RR];
    save(fullfile(out_folder, 'frn_rewp_by_stage_combined.mat'), 'frn_rewp_stage_table');
    fprintf('Saved combined FRN/RewP per-stage table: %d rows.\n', height(frn_rewp_stage_table));
end

% =========================================================================
%% LOCAL FUNCTIONS
% =========================================================================
function T = load_features(paths)
% Load the first existing file and return its group_table.
if ischar(paths) || isstring(paths); paths = {char(paths)}; end
for i = 1:numel(paths)
    if exist(paths{i}, 'file')
        S = load(paths{i});
        if isfield(S, 'group_table'); T = S.group_table; return; end
        fn = fieldnames(S);
        T = S.(fn{1}); return;   % first variable as fallback
    end
end
error('S4: no feature table found among: %s', strjoin(string(paths), ', '));
end

function T = zscore_within_subject(T, feats)
subj_list = unique(string(T.subj_id));
for f = 1:numel(feats)
    fn = feats{f};
    if ~ismember(fn, T.Properties.VariableNames); continue; end
    fn_z = [fn '_z'];
    T.(fn_z) = nan(height(T), 1);
    for si = 1:numel(subj_list)
        mask = string(T.subj_id) == subj_list(si);
        vals = T.(fn)(mask);
        sd = std(vals, 'omitnan');
        if sd > 0
            T.(fn_z)(mask) = (vals - mean(vals, 'omitnan')) / sd;
        end
    end
end
end

function b = get_block_index(T)
% Return the block index column, tolerating legacy names.
if ismember('block', T.Properties.VariableNames)
    b = kh_to_numeric(T.block);
elseif ismember('block_number', T.Properties.VariableNames)
    b = kh_to_numeric(T.block_number);
elseif ismember('blocknum', T.Properties.VariableNames)
    b = kh_to_numeric(T.blocknum);
else
    b = nan(height(T), 1);
end
end

function frn = try_load_frn(path)
frn = table();
if exist(path, 'file')
    S = load(path);
    if isfield(S, 'frn_rewp_stage_table'); frn = S.frn_rewp_stage_table; end
end
end
