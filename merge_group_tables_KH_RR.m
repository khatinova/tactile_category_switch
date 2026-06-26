remote = 0;
if remote == 1
    base_path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch';
else

    base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch';
end

epoch_file_folder = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis', 'Epoched_data');
RR_epoch_folder   = fullfile(base_path, 'Salient mod switch RR', 'Results', 'EEG analysis', 'Epoched_data');
figure_output_folder = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis', 'Figures', 'RQ_analysis_combined_data');

if ~exist(figure_output_folder,'dir'), mkdir(figure_output_folder); end


%% Join KH's and RR's group_tables. IMPORTANT!
load(fullfile(RR_epoch_folder, "group_stage_table_RR_v7.mat"))


features_to_zscore = {'FCz_neg_peak_amp','FRN_amp','RewP_amp', ...
                      'P300_amp','Theta_amp', ...
                      'PLV_fp','PLV_fs','PLV_fp_pairwise','PLV_fs_pairwise'};


% Add formatted subject ID column (e.g., Ox03)
n_rows = height(group_table);
group_table.subj_id = strings(n_rows,1);

for i = 1:n_rows  
        group_table.subj_id(i) = "Nc" + compose("%02d", double(group_table.subj(i)));
end
subj_list = unique(group_table.subj_id);
for f = 1:numel(features_to_zscore)
    fn   = features_to_zscore{f};
    fn_z = [fn '_z'];
    if ~ismember(fn, group_table.Properties.VariableNames); continue; end
    group_table.(fn_z) = nan(height(group_table), 1);
    for si = 1:numel(subj_list)
        mask = group_table.subj_id == subj_list(si);
        vals = group_table.(fn)(mask);
        mn   = mean(vals, 'omitnan');
        sd   = std(vals,  'omitnan');
        if sd > 0
            group_table.(fn_z)(mask) = (vals - mn) / sd;
        end
    end
end

group_table_RR = group_table;

clear group_table

load(fullfile(epoch_file_folder,'group_stage_table_features.mat'), 'group_table');
n_rows = height(group_table);
group_table.subj_id = strings(n_rows,1);

for i = 1:n_rows

    group_table.subj_id(i) = "Ox" + compose("%02d", double(group_table.subj(i)));

end
group_table_KH = group_table;

% ── CORRECT zero-padded subj_id for both cohorts ─────────────────────────
% compose("%02d", n) → "03", "04" etc. regardless of input type.
% Do NOT use strcat("Ox", string(subj)) — string() drops leading zeros.

group_table_KH.subj_id = "Ox" + compose("%02d", double(group_table_KH.subj));
group_table_RR.subj_id = "Nc" + compose("%02d", double(group_table_RR.subj));

% ── Everything else unchanged ─────────────────────────────────────────────
group_table_KH.feedback_modality = repmat("visual", height(group_table_KH), 1);
is_audio = double(group_table_KH.subj) <= 8;
group_table_KH.feedback_modality(is_audio) = "auditory";
group_table_RR.feedback_modality = repmat("visual", height(group_table_RR), 1);

group_table_KH.stimulus_modality = repmat("tactile", height(group_table_KH), 1);
mask = double(group_table_KH.subj) >= 9 & double(group_table_KH.subj) <= 17 ...
     & group_table_KH.block_number == 1;
group_table_KH.stimulus_modality(mask) = "visual+tactile";
group_table_RR.stimulus_modality = repmat("tactile", height(group_table_RR), 1);

group_table_KH.practice_task = strings(height(group_table_KH), 1);
group_table_KH.practice_task(double(group_table_KH.subj) <= 8)  = "matching";
group_table_KH.practice_task(double(group_table_KH.subj) >= 9  & ...
                              double(group_table_KH.subj) <= 16) = "structured";
group_table_KH.practice_task(double(group_table_KH.subj) >= 17) = "criterion";
group_table_RR.practice_task = repmat("criterion", height(group_table_RR), 1);

% ── Join tables ───────────────────────────────────────────────────────────
vars = group_table_KH.Properties.VariableNames;
group_table_RR  = group_table_RR(:, vars);
group_table = [group_table_KH; group_table_RR];

% ── Convert to categorical ────────────────────────────────────────────────
group_table.subj_id           = categorical(group_table.subj_id);
group_table.feedback_modality = categorical(group_table.feedback_modality);
group_table.stimulus_modality = categorical(group_table.stimulus_modality);
group_table.practice_task     = categorical(group_table.practice_task);

% ── Move subj_id to the first column ─────────────────────────────────────
group_table = movevars(group_table, 'subj_id', 'Before', 1);

% ── Save ──────────────────────────────────────────────────────────────────
save(fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis', ...
     'Epoched_data', 'group_feature_table_combined.mat'), 'group_table');
fprintf('Saved combined table: %d rows, %d columns. First column: %s\n', ...
    height(group_table), width(group_table), group_table.Properties.VariableNames{1});