% =============================================================================
% D_INDIVIDUAL_DIFFERENCES_UNCERTAINTY_v1.m
%
% PURPOSE
% -------
% Quantifies and visualises individual differences in adaptation to
% uncertainty across two timescales:
%
%   SHORT-TERM uncertainty  — probabilistic (P) vs deterministic (D)
%                             feedback within a block; the FRN/RewP signal
%                             is expected to be attenuated in P blocks
%                             because outcome PE is lower in expectation
%                             (Cavanagh et al., 2010, Psychophysiology).
%
%   LONG-TERM uncertainty   — the reversal, which demands model updating
%                             regardless of block type. Tracked by the
%                             P300 (context updating; Polich, 2007) and
%                             frontal theta (conflict; Cavanagh & Frank,
%                             2014).
%
%   HISTORY OF UNCERTAINTY  — how exposure to P blocks before the current
%                             block shapes the neural/behavioural response
%                             to reversals in subsequent blocks
%                             (Behrens et al., 2007, Nat Neurosci).
%
% INPUTS (expected in workspace or loaded below)
% -----------------------------------------------
%   group_table_combined   — output of B_outcome_ERP_analysis_v9b.m
%                            (group_feature_table_combined_v9b.mat)
%   figure_output_folder   — path string for saving figures
%
% COLUMNS USED
% ------------
% From stage_table pipeline (v9b):
%   subj_id, cohort, block, stage (LN/LE/RN/RE), block_type (D/P)
%   correct (0/1), false_fb (logical), fb_shown_correct
%   confidence, trial_continuous, has_eeg_epoch
%
% EEG features (all _z are within-subject z-scores from v9b):
%   FRN_mean_amp_z, FRN_peak_amp_z, FRN_peak_lat
%   RewP_mean_amp_z, RewP_peak_amp_z
%   P300_amp_z, P300_peak_lat
%   Theta_amp_z
%   PLV_fp_z, PLV_fs_z, PLV_fp_pairwise_z, PLV_fs_pairwise_z
%
% From behav_table (KH cohort, joined in v9b):
%   PE, value, alpha_det, alpha_prob, alpha_pre, alpha_post
%   conf_weighted_PE, rev_sensitive_v (truncated name)
%   stay_choice, revTrial, rev_state, transition
%
% ANALYSIS STRUCTURE
% ------------------
% SECTION 0  — Load, type coercion, derived variables
% SECTION 1  — Behavioural individual differences (ID) indices
%              1A: Reversal learning index (post-reversal accuracy trajectory)
%              1B: Context-sensitive learning rate (alpha_prob/alpha_det)
%              1C: Win-stay / lose-switch rates by block type
%              1D: Confidence calibration (metacognitive efficiency proxy)
% SECTION 2  — RQ1: FRN/RewP ~ block_type × correct × stage
%              (short-term uncertainty modulates error signal)
% SECTION 3  — RQ2: FRN ~ confidence × block_type × correct
%              (precision weighting and metacognition)
% SECTION 4  — RQ3: Theta ~ stage × block_type (incorrect trials)
%              (conflict signal and adaptive control)
% SECTION 5  — RQ4/5: PLV pathways ~ stage × block_type
%              (fronto-parietal vs fronto-somatosensory connectivity)
% SECTION 6  — Cross-level individual differences: neural ~ behavioural ID
%              (correlate subject-level neural indices with behavioural)
% SECTION 7  — RL model latents ~ neural markers (PE, alpha, value)
%
% MODELLING NOTES
% ---------------
% All LME models use within-subject z-scored EEG features (from v9b).
% Random effects: (1 + [within-subject predictor] | subj_id) where
% justified. Stage is treated as a nominal nuisance covariate or as a
% focussed LE/RN contrast, not as an ordinal slope, to avoid the
% convergence issues documented in the RQ script (BUG 9).
%
% REFERENCES (cited inline at relevant models)
% ---------------------------------------------
% Behrens et al. (2007) Nat Neurosci 10:1214-1221
% Cavanagh et al. (2010) Psychophysiology 47:395-405
% Cavanagh & Frank (2014) Trends Cogn Sci 18:414-421
% Holroyd & Coles (2002) Psychol Rev 109:679-709
% Polich (2007) Clin Neurophysiol 118:2128-2148
% Yu & Dayan (2005) Neuron 46:681-692
% =============================================================================

%% Results summary:


clear; close all;

% -------------------------------------------------------------------------
%% 0A — PATHS & LOAD
% -------------------------------------------------------------------------
remote = 0;
if remote
    base_path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch';
else
    base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch';
end

KH_epoch_folder = fullfile(base_path, 'Salient mod switch KH', 'Results', ...
    'EEG analysis', 'Epoched_data');
figure_output_folder = fullfile(base_path, 'Salient mod switch KH', 'Results', ...
    'EEG analysis', 'Figures', 'IndivDiff_v1');
if ~exist(figure_output_folder, 'dir'), mkdir(figure_output_folder); end

combined_file = fullfile(KH_epoch_folder, 'group_feature_table_combined_v9c.mat');
if exist(combined_file, 'file')
    S = load(combined_file, 'group_table_combined');
    gt_raw = S.group_table_combined;
    fprintf('Loaded combined table: %d rows, %d cols\n', height(gt_raw), width(gt_raw));
else
    error(['group_feature_table_combined_v9b.mat not found.\n' ...
           'Run B_outcome_ERP_analysis_v9b.m first.']);
end

% -------------------------------------------------------------------------
%% 0B — TYPE COERCION
% -------------------------------------------------------------------------
gt = gt_raw;

% --- identifiers and grouping ---
gt.subj_id    = categorical(string(gt.subj_id));
gt.cohort     = categorical(string(gt.cohort));
gt.block_type = categorical(string(gt.block_type), {'D','P'});
gt.stage      = categorical(string(gt.stage), {'LN','LE','RN','RE'}, 'Ordinal', false);

% --- binary outcome: keep as numeric 0/1 for arithmetic, add categorical ---
if islogical(gt.correct), gt.correct = double(gt.correct); end
gt.correct_cat = categorical(gt.correct, [0 1], {'Incorrect','Correct'});

% --- false_fb: numeric + categorical ---
if islogical(gt.false_fb), gt.false_fb = double(gt.false_fb); end
gt.false_fb_cat = categorical(gt.false_fb, [0 1], {'TrueFB','FalseFB'});

% --- reversal transition stage: LE = last epoch before reversal,
%     RN = first epoch after reversal (both are the critical moments).
%     Only rows with stage ∈ {LE, RN} are assigned. ---
gt.reversal_transition_stage = categorical(repmat("other", height(gt), 1));
gt.reversal_transition_stage(gt.stage == 'LE') = "LE";
gt.reversal_transition_stage(gt.stage == 'RN') = "RN";

% -------------------------------------------------------------------------
%% 0C — DERIVE BLOCK-SEQUENCE VARIABLES
% -------------------------------------------------------------------------
% These mirror the RQ-script derivation but are reproduced here so this
% script is self-contained.
%
% Variables created per subject × block:
%   prev_block_type  — categorical: 'first', 'D', or 'P'
%   block_transition — categorical: 'first', 'DtoD', 'DtoP', 'PtoD', 'PtoP'
%   n_prev_P         — integer: how many previous blocks were type P
%                      (cumulative exposure to probabilistic uncertainty)

gt.prev_block_type  = categorical(repmat("first", height(gt), 1));
gt.block_transition = categorical(repmat("first", height(gt), 1));
gt.n_prev_P         = zeros(height(gt), 1);

subjects = categories(gt.subj_id);
for si = 1:numel(subjects)
    sn  = subjects{si};
    sm  = gt.subj_id == sn;
    bls = sort(unique(gt.block(sm))');

    for bi = 1:numel(bls)
        bm        = sm & gt.block == bls(bi);
        curr_type = string(gt.block_type(find(bm, 1)));

        if bi == 1
            prev_type = "first";
            trans     = "first";
            nprevP    = 0;
        else
            prev_bm   = sm & gt.block == bls(bi-1);
            prev_type = string(gt.block_type(find(prev_bm, 1)));
            trans     = prev_type + "to" + curr_type;
            nprevP    = sum(arrayfun(@(pb) ...
                string(gt.block_type(find(sm & gt.block==pb, 1))) == "P", ...
                bls(1:bi-1)));
        end

        gt.prev_block_type(bm)  = categorical(prev_type);
        gt.block_transition(bm) = categorical(trans);
        gt.n_prev_P(bm)         = nprevP;
    end
end

gt.prev_block_type  = categorical(gt.prev_block_type);
gt.block_transition = categorical(gt.block_transition);

% -------------------------------------------------------------------------
%% 0D — WITHIN-SUBJECT Z-SCORES
% -------------------------------------------------------------------------
% The v9b script z-scores EEG features. We repeat here for behavioural
% variables (confidence, PE, value, alpha indices) that were not z-scored
% in v9b.

behav_to_zscore = {'confidence','PE','value','conf_weighted_PE', ...
                   'alpha_pre','alpha_post','alpha_det','alpha_prob'};

for f = 1:numel(behav_to_zscore)
    fn = behav_to_zscore{f};
    if ~ismember(fn, gt.Properties.VariableNames), continue; end
    fn_z = [fn '_z'];
    gt.(fn_z) = nan(height(gt), 1);
    for si = 1:numel(subjects)
        m    = gt.subj_id == subjects{si};
        vals = gt.(fn)(m);
        mn   = mean(vals, 'omitnan');
        sd   = std(vals,  'omitnan');
        if sd > 0
            gt.(fn_z)(m) = (vals - mn) / sd;
        end
    end
end

% Confidence z-score specifically — used heavily in RQ2
if ismember('confidence_z', gt.Properties.VariableNames)
    gt.conf_z = gt.confidence_z;
elseif ismember('confidence', gt.Properties.VariableNames)
    gt.conf_z = nan(height(gt), 1);
    for si = 1:numel(subjects)
        m    = gt.subj_id == subjects{si};
        vals = gt.confidence(m);
        mn   = mean(vals,'omitnan'); sd = std(vals,'omitnan');
        if sd > 0, gt.conf_z(m) = (vals - mn)/sd; end
    end
end

% Trial position within block (time covariate, z-scored globally)
if ismember('trial_continuous', gt.Properties.VariableNames)
    gt.trial_z = (gt.trial_continuous - mean(gt.trial_continuous,'omitnan')) ./ ...
                  std(gt.trial_continuous,'omitnan');
else
    gt.trial_z = zeros(height(gt), 1);
end

if ismember('block', gt.Properties.VariableNames)
    gt.block_z = (gt.block - mean(gt.block,'omitnan')) ./ ...
                  std(gt.block,'omitnan');
else
    gt.block_z = zeros(height(gt), 1);
end

fprintf('\nFull table after coercion: %d rows, %d subjects\n', ...
    height(gt), numel(subjects));

% =========================================================================
%% SECTION 1 — BEHAVIOURAL INDIVIDUAL DIFFERENCES INDICES
% =========================================================================
% For each subject we compute scalar indices summarising:
%   (a) how quickly they adapt after each reversal
%   (b) how context-sensitively they modulate learning rates
%   (c) win-stay / lose-switch rates per block type
%   (d) confidence calibration (overclaiming on false-negative trials)
%
% These indices are then used in Section 6 to correlate with neural markers.

fprintf('\n=== SECTION 1: Behavioural individual difference indices ===\n');

has_RL = ismember('alpha_pre', gt.Properties.VariableNames) && ...
         ismember('alpha_post', gt.Properties.VariableNames);
has_PE = ismember('PE', gt.Properties.VariableNames);

% Pre-allocate subject-level table
id_table = table();
id_table.subj_id = categorical(subjects);
id_table.cohort  = categorical(repmat("", numel(subjects), 1));
id_table.n_blocks = nan(numel(subjects), 1);

% --- 1A: Post-reversal accuracy — area under the post-reversal accuracy
%     curve relative to chance (0.5), window = first 20 trials after rev.
%     Separately for D and P blocks.
%     Theory: faster recovery = better long-term uncertainty adaptation.
%     Ref: Cools et al. (2002) Cereb Cortex 12:469-476.

id_table.rev_AUC_D       = nan(numel(subjects), 1);   % area D blocks
id_table.rev_AUC_P       = nan(numel(subjects), 1);   % area P blocks
id_table.rev_AUC_ratio   = nan(numel(subjects), 1);   % P/D ratio
id_table.trials_to_crit_D = nan(numel(subjects), 1);   % trials until 3/3 correct post-rev
id_table.trials_to_crit_P = nan(numel(subjects), 1);

POST_REV_WINDOW = 20;   % trials after reversal to examine
CRIT_RUN        = 3;    % consecutive correct to call "criterion"

% --- 1B: Context-sensitive learning rate (alpha_prob / alpha_det) ---
%     Theory: high ratio means the subject appropriately inflates their
%     learning rate under probabilistic feedback, matching the volatility
%     of that environment. Ref: Behrens et al. (2007).
id_table.alpha_ratio    = nan(numel(subjects), 1);   % alpha_prob / alpha_det
id_table.delta_alpha    = nan(numel(subjects), 1);   % alpha_post - alpha_pre
id_table.rev_sensitivity = nan(numel(subjects), 1);  % mean(abs(rev_sensitive_v)) if available

% --- 1C: Win-stay / lose-switch rates ---
%     Computed per block type, true feedback only.
%     Theory: lose-switch in D but not P = appropriate uncertainty sensitivity.
%     Ref: Daw et al. (2006) Nature 441:876-879.
id_table.ws_D = nan(numel(subjects), 1);   % win-stay rate, D blocks
id_table.ls_D = nan(numel(subjects), 1);   % lose-switch rate, D blocks
id_table.ws_P = nan(numel(subjects), 1);   % win-stay rate, P blocks
id_table.ls_P = nan(numel(subjects), 1);   % lose-switch rate, P blocks
id_table.ls_ratio = nan(numel(subjects), 1); % ls_D / ls_P (context sensitivity)

% --- 1D: Pre-outcome confidence as a surprise-scaling index ---
%     TRIAL ORDER: stimulus → response → CONFIDENCE → feedback.
%     Confidence is submitted BEFORE the outcome is seen. It therefore
%     indexes the subject's pre-outcome certainty about their choice,
%     not a post-feedback belief revision.
%
%     The relevant individual-difference index is the mean pre-outcome
%     confidence on trials that subsequently receive NEGATIVE feedback,
%     split by whether that negative feedback was TRUE (deserved: they
%     were wrong) or FALSE (undeserved: they were right but were told
%     wrong). The contrast is:
%
%       conf_surprise_index = mean conf on false-negative trials
%                           - mean conf on true-negative trials
%
%     Interpretation: a POSITIVE score means the subject enters
%     false-negative trials with higher pre-outcome confidence than
%     true-negative trials. This is the expected direction — on
%     false-negative trials the subject was objectively correct, so
%     they should be MORE confident before seeing the (misleading)
%     feedback. A high score therefore indicates that the MISMATCH
%     between pre-outcome certainty and the received negative feedback
%     is large on those trials — indexing the potential magnitude of the
%     prediction error surprise signal the brain has to resolve.
%     It is NOT a measure of post-feedback stubbornness.
%
%     This index is then related to FRN amplitude on false-negative
%     trials in Section 6: subjects with a higher surprise index
%     (more confident going in → larger confidence-feedback mismatch)
%     should show a larger FRN on those trials, because the PE is
%     larger. Ref: Boldt & Yeung (2015) J Neurosci 35:2058 on
%     confidence-weighted PE signals; Meyniel et al. (2015).
%
%     Computed in P blocks only (the only block type with false feedback).
id_table.conf_on_false_neg  = nan(numel(subjects), 1);  % pre-FB conf on false-neg trials
id_table.conf_on_true_neg   = nan(numel(subjects), 1);  % pre-FB conf on true-neg trials
id_table.conf_surprise_index = nan(numel(subjects), 1); % false_neg - true_neg (higher = more surprised)

has_stay = ismember('stay_choice', gt.Properties.VariableNames);
has_conf = ismember('confidence', gt.Properties.VariableNames);
has_rev_state = ismember('rev_state', gt.Properties.VariableNames);
rev_sens_col = 'rev_sensitive_v';   % truncated name from v9b

for si = 1:numel(subjects)
    sn = subjects{si};
    sm = gt.subj_id == sn;
    id_table.cohort(si)   = categorical(string(gt.cohort(find(sm,1))));
    id_table.n_blocks(si) = numel(unique(gt.block(sm)));

    % --- 1A: post-reversal AUC ---
    for bt = {'D','P'}
        bm = sm & string(gt.block_type) == bt{1} & ~gt.false_fb;
        bls_bt = sort(unique(gt.block(bm))');
        auc_vals = nan(1, numel(bls_bt));
        ttc_vals = nan(1, numel(bls_bt));

        for bi = 1:numel(bls_bt)
            bk = bls_bt(bi);
            bkm = bm & gt.block == bk;

            % Find the reversal trial within this block
            if has_rev_state
                rev_rows = find(bkm & gt.rev_state == 1);
            else
                % Fallback: treat middle 30-70 as reversal region
                % and find first post-reversal row heuristically.
                tc_bk = gt.trial_continuous(bkm);
                rev_rows = find(bkm & gt.trial_continuous > min(tc_bk)+30 & ...
                    gt.trial_continuous <= min(tc_bk)+70, 1);
            end
            if isempty(rev_rows), continue; end
            rev_start = min(rev_rows);

            % All rows in window after reversal
            tc_start = gt.trial_continuous(rev_start);
            win_mask = bkm & gt.trial_continuous >= tc_start & ...
                       gt.trial_continuous < tc_start + POST_REV_WINDOW;
            c_win = gt.correct(win_mask);
            if numel(c_win) < 5, continue; end

            % AUC = mean accuracy above chance
            auc_vals(bi) = mean(c_win, 'omitnan') - 0.5;

            % Trials to criterion: first run of CRIT_RUN consecutive correct
            ttc = nan;
            run = 0;
            for ti = 1:numel(c_win)
                if c_win(ti) == 1
                    run = run + 1;
                    if run >= CRIT_RUN
                        ttc = ti - CRIT_RUN + 1;
                        break;
                    end
                else
                    run = 0;
                end
            end
            ttc_vals(bi) = ttc;
        end

        if strcmp(bt{1}, 'D')
            id_table.rev_AUC_D(si) = mean(auc_vals, 'omitnan');
            id_table.trials_to_crit_D(si) = mean(ttc_vals, 'omitnan');
        else
            id_table.rev_AUC_P(si) = mean(auc_vals, 'omitnan');
            id_table.trials_to_crit_P(si) = mean(ttc_vals, 'omitnan');
        end
    end
    id_table.rev_AUC_ratio(si) = id_table.rev_AUC_P(si) / max(id_table.rev_AUC_D(si), 0.01);

    % --- 1B: learning rates ---
    if has_RL
        subj_alpha_det  = mean(gt.alpha_det(sm),  'omitnan');
        subj_alpha_prob = mean(gt.alpha_prob(sm), 'omitnan');
        subj_alpha_pre  = mean(gt.alpha_pre(sm),  'omitnan');
        subj_alpha_post = mean(gt.alpha_post(sm), 'omitnan');
        if subj_alpha_det > 0
            id_table.alpha_ratio(si) = subj_alpha_prob / subj_alpha_det;
        end
        id_table.delta_alpha(si) = subj_alpha_post - subj_alpha_pre;
    end
    if ismember(rev_sens_col, gt.Properties.VariableNames)
        id_table.rev_sensitivity(si) = mean(abs(gt.(rev_sens_col)(sm)), 'omitnan');
    end

    % --- 1C: win-stay / lose-switch ---
    if has_stay
        for bt = {'D','P'}
            bm_true = sm & string(gt.block_type) == bt{1} & ~gt.false_fb;
            % Win-stay: correct on trial t AND stays on trial t+1
            ws_num = sum(bm_true & gt.correct == 1 & gt.stay_choice == 1, 'omitnan');
            ws_den = sum(bm_true & gt.correct == 1, 'omitnan');
            ls_num = sum(bm_true & gt.correct == 0 & gt.stay_choice == 0, 'omitnan');
            ls_den = sum(bm_true & gt.correct == 0, 'omitnan');
            ws = ws_num / max(ws_den, 1);
            ls = ls_num / max(ls_den, 1);
            if strcmp(bt{1}, 'D')
                id_table.ws_D(si) = ws;
                id_table.ls_D(si) = ls;
            else
                id_table.ws_P(si) = ws;
                id_table.ls_P(si) = ls;
            end
        end
        id_table.ls_ratio(si) = id_table.ls_D(si) / max(id_table.ls_P(si), 0.01);
    end

    % --- 1D: Pre-outcome confidence surprise index ---
    % Confidence is recorded BEFORE feedback is shown. So we compare:
    %   conf_on_false_neg: mean pre-outcome confidence on trials where the
    %     subject was CORRECT but received NEGATIVE feedback (false-neg).
    %     These are trials where the subject should tend to be confident
    %     (they got it right), so confidence should be relatively high.
    %   conf_on_true_neg: mean pre-outcome confidence on trials where the
    %     subject was INCORRECT and received NEGATIVE feedback (true-neg).
    %     These are genuinely wrong trials; confidence may be lower.
    %   conf_surprise_index = false_neg - true_neg:
    %     A positive value means the subject enters false-negative trials
    %     more confidently than true-negative ones. This is the expected
    %     direction and indexes the magnitude of the confidence-feedback
    %     mismatch on those trials (i.e. how surprised the FRN should be).
    if has_conf
        fn_mask = sm & string(gt.block_type)=='P' & gt.false_fb==1 & gt.correct==1;
        tn_mask = sm & string(gt.block_type)=='P' & gt.false_fb==0 & gt.correct==0;
        id_table.conf_on_false_neg(si)   = mean(gt.confidence(fn_mask), 'omitnan');
        id_table.conf_on_true_neg(si)    = mean(gt.confidence(tn_mask), 'omitnan');
        id_table.conf_surprise_index(si) = id_table.conf_on_false_neg(si) - ...
                                           id_table.conf_on_true_neg(si);
    end
end

fprintf('Subject-level ID index table created: %d subjects\n', height(id_table));

% =========================================================================
%% SECTION 1E — SUMMARY FIGURE: BEHAVIOURAL ID INDICES
% =========================================================================

fig = figure('Position',[50 50 1400 900]);
sgtitle('Individual differences in uncertainty adaptation — behavioural indices', ...
    'FontSize', 12, 'FontWeight', 'bold');

n_plot = 0;

% 1. Post-reversal AUC: D vs P
ax1 = subplot(3,4,1);
plot_subject_scatter(ax1, id_table.rev_AUC_D, id_table.rev_AUC_P, ...
    id_table.cohort, ...
    'Rev. AUC — D blocks', 'Rev. AUC — P blocks', ...
    'Post-reversal accuracy AUC');
refline_identity(ax1);

% 2. AUC ratio distribution
ax2 = subplot(3,4,2);
plot_id_histogram(ax2, id_table.rev_AUC_ratio, ...
    'AUC ratio (P/D)', 'Context-sensitive reversal adaptation', ...
    id_table.cohort);
xline(ax2, 1, 'k--', 'Equal');

% 3. Trials to criterion: D vs P
ax3 = subplot(3,4,3);
plot_subject_scatter(ax3, id_table.trials_to_crit_D, id_table.trials_to_crit_P, ...
    id_table.cohort, ...
    'Trials to crit (D)', 'Trials to crit (P)', ...
    'Trials to criterion after reversal');
refline_identity(ax3);

% 4. Delta-alpha: post - pre reversal learning rate
if has_RL
    ax4 = subplot(3,4,4);
    plot_id_histogram(ax4, id_table.delta_alpha, ...
        '\Delta\alpha (post - pre reversal)', ...
        'Reversal-driven learning rate shift', id_table.cohort);
    xline(ax4, 0, 'k--', 'No shift');
end

% 5. Alpha ratio (prob/det) — context-sensitive updating
if has_RL
    ax5 = subplot(3,4,5);
    plot_id_histogram(ax5, id_table.alpha_ratio, ...
        '\alpha ratio (P/D)', 'Context-sensitive learning rate', ...
        id_table.cohort);
    xline(ax5, 1, 'k--', 'Equal');
end

% 6. Win-stay: D vs P
if has_stay
    ax6 = subplot(3,4,6);
    plot_subject_scatter(ax6, id_table.ws_D, id_table.ws_P, ...
        id_table.cohort, 'Win-stay D', 'Win-stay P', 'Win-stay rate');
    refline_identity(ax6);

    ax7 = subplot(3,4,7);
    plot_subject_scatter(ax7, id_table.ls_D, id_table.ls_P, ...
        id_table.cohort, 'Lose-switch D', 'Lose-switch P', ...
        'Lose-switch rate (uncertainty sensitivity)');
    refline_identity(ax7);
end

% 7. Pre-outcome confidence surprise index
if has_conf
    ax8 = subplot(3,4,8);
    plot_id_histogram(ax8, id_table.conf_surprise_index, ...
        'Conf surprise index (false-neg - true-neg pre-FB)', ...
        'Pre-outcome confidence: false vs true negative trials', id_table.cohort);
    xline(ax8, 0, 'k--', 'No difference');
end

% 8. Correlation matrix of key ID indices
ax9 = subplot(3,4,[9 10 11 12]);
id_vars = {'rev_AUC_D','rev_AUC_P','trials_to_crit_D','trials_to_crit_P'};
if has_RL
    id_vars = [id_vars, {'alpha_ratio','delta_alpha'}];
end
if has_stay
    id_vars = [id_vars, {'ls_ratio'}];
end
if has_conf
    id_vars = [id_vars, {'conf_surprise_index'}];
end

id_mat = nan(height(id_table), numel(id_vars));
for v = 1:numel(id_vars)
    if ismember(id_vars{v}, id_table.Properties.VariableNames)
        id_mat(:,v) = id_table.(id_vars{v});
    end
end
ok_rows = all(~isnan(id_mat), 2);
if sum(ok_rows) > 3
    R = corr(id_mat(ok_rows,:), 'Rows','complete');
    imagesc(ax9, R, [-1 1]);
    colormap(ax9, redblue_colormap());
    colorbar(ax9);
    set(ax9, 'XTick',1:numel(id_vars), 'XTickLabel', strrep(id_vars,'_','\_'), ...
        'XTickLabelRotation', 30, 'YTick', 1:numel(id_vars), ...
        'YTickLabel', strrep(id_vars,'_','\_'), 'FontSize', 8);
    title(ax9, 'Behavioural ID index correlation matrix');
    % Annotate significant correlations
    for ri = 1:numel(id_vars)
        for ci = 1:numel(id_vars)
            if ri ~= ci && ~isnan(R(ri,ci))
                n_ok = sum(ok_rows);
                t_stat = R(ri,ci) * sqrt((n_ok-2)/(1-R(ri,ci)^2));
                p_val  = 2*(1-tcdf(abs(t_stat), n_ok-2));
                if p_val < 0.05
                    text(ax9, ci, ri, sprintf('%.2f*',R(ri,ci)), ...
                        'HorizontalAlignment','center','FontSize',7,'FontWeight','bold');
                end
            end
        end
    end
end

exportgraphics(fig, fullfile(figure_output_folder,'S1_Behavioural_ID_indices.pdf'), ...
    'ContentType','vector');

% =========================================================================
%% SECTION 2 — RQ1: FRN/RewP × BLOCK TYPE × OUTCOME × STAGE
% =========================================================================
% Short-term uncertainty (block_type D vs P) modulates the outcome
% prediction error signal (FRN/RewP). Key predictions from Holroyd &
% Coles (2002) / Cavanagh et al. (2010):
%   • FRN_mean_amp larger (more negative) in D vs P on incorrect trials
%   • RewP_mean_amp larger in D vs P on correct trials
%   • At reversal stages (LE, RN): both signals may be compressed because
%     the PE itself is ambiguous (reversal = large unexpected error)
%
% RQ1A: Full stage × block_type × correct model.
% RQ1B: LE vs RN focused model (reversal transition contrast).
% RQ1C: Past uncertainty: does prev_block_type modulate FRN/RewP?
% RQ1D: False feedback dissociation (P blocks only).

fprintf('\n=== SECTION 2: RQ1 — FRN/RewP models ===\n');

% Primary EEG DV: FRN_mean_amp_z (always defined, no exclusion issue).
% FRN_peak_amp_z used as secondary/sensitivity check.
FRN_DV  = 'FRN_mean_amp_z';
RewP_DV = 'RewP_mean_amp_z';

% Check that _z columns exist; if not, make them from raw columns.
[gt, FRN_DV]  = ensure_z_column(gt, FRN_DV,  subjects);
[gt, RewP_DV] = ensure_z_column(gt, RewP_DV, subjects);
[gt, ~]       = ensure_z_column(gt, 'P300_amp_z',   subjects);
[gt, ~]       = ensure_z_column(gt, 'Theta_amp_z',  subjects);
[gt, ~]       = ensure_z_column(gt, 'PLV_fp_z',     subjects);
[gt, ~]       = ensure_z_column(gt, 'PLV_fs_z',     subjects);

% Working table: true feedback only, valid FRN
gt_true = gt(~gt.false_fb & ~isnan(gt.(FRN_DV)), :);
gt_true.block_type       = categorical(string(gt_true.block_type),       {'D','P'});
gt_true.stage            = categorical(string(gt_true.stage),            {'LN','LE','RN','RE'}, 'Ordinal', false);
gt_true.correct_cat      = categorical(gt_true.correct, [0 1],           {'Incorrect','Correct'});
gt_true.prev_block_type  = categorical(string(gt_true.prev_block_type));
gt_true.block_transition = categorical(string(gt_true.block_transition));

% Reversal-only subset
gt_rev = gt_true(gt_true.stage == 'LE' | gt_true.stage == 'RN', :);
gt_rev.reversal_transition_stage = categorical(string(gt_rev.stage), {'LE','RN'});

% -------------------------------------------------------------------------
% RQ1A: Full-stage model
% Theory: stage captures position within the block (LN=early no-rev,
%   LE=pre-reversal, RN=post-reversal early, RE=post-reversal late).
%   block_type × correct interaction tests whether the FRN outcome effect
%   is compressed in P vs D blocks (Cavanagh et al., 2010).
% Random effect: random slope for correct given subjects vary in how
%   strongly they differentiate correct/incorrect at the neural level.
% -------------------------------------------------------------------------
fprintf('\n--- RQ1A: Full stage model (FRN) ---\n');
mdl_frn_1a = fitlme(gt_true, ...
    [FRN_DV ' ~ block_type * correct_cat * stage + ' ...
     '(1 + correct_cat | subj_id)'], ...
    'FitMethod','REML','DummyVarCoding','effects');
disp(mdl_frn_1a.Coefficients);
anova_frn_1a = anova(mdl_frn_1a, 'DFMethod','Satterthwaite');
disp(anova_frn_1a);

fprintf('\n--- RQ1A: Full stage model (RewP) ---\n');
mdl_rewp_1a = fitlme(gt_true, ...
    [RewP_DV ' ~ block_type * correct_cat * stage + ' ...
     '(1 + correct_cat | subj_id)'], ...
    'FitMethod','REML','DummyVarCoding','effects');
disp(mdl_rewp_1a.Coefficients);
anova_rewp_1a = anova(mdl_rewp_1a, 'DFMethod','Satterthwaite');
disp(anova_rewp_1a);

% LRT: does adding block_type interaction improve over intercept-only for
% uncertainty? (model comparison per Bates et al. 2015)
mdl_frn_1a_null = fitlme(gt_true, ...
    [FRN_DV ' ~ correct_cat * stage + block_type + ' ...
     '(1 + correct_cat | subj_id)'], ...
    'FitMethod','ML','DummyVarCoding','effects');
mdl_frn_1a_full = fitlme(gt_true, ...
    [FRN_DV ' ~ block_type * correct_cat * stage + ' ...
     '(1 + correct_cat | subj_id)'], ...
    'FitMethod','ML','DummyVarCoding','effects');
lrt_1a_frn = compare(mdl_frn_1a_null, mdl_frn_1a_full);
fprintf('\nLRT: adding block_type interactions to FRN model:\n');
disp(lrt_1a_frn);

% -------------------------------------------------------------------------
% RQ1B: LE vs RN reversal transition (focussed)
% Theory: LE is the last trial before the rule switches (maximum commitment
%   to old rule; PE of incorrect = very large). RN is the first post-switch
%   trial (maximum surprise if still using old rule). This contrast isolates
%   the transition period without the noise of early-stable (LN) or
%   late-stable (RE) epochs. Prediction: FRN largest at LE-incorrect and
%   at RN-incorrect (both are high PE moments), but for different reasons.
% -------------------------------------------------------------------------
fprintf('\n--- RQ1B: LE vs RN reversal transition ---\n');
mdl_frn_1b = fitlme(gt_rev, ...
    [FRN_DV ' ~ block_type * correct_cat * reversal_transition_stage + ' ...
     '(1 + correct_cat | subj_id)'], ...
    'FitMethod','REML','DummyVarCoding','effects');
disp(mdl_frn_1b.Coefficients);
anova_frn_1b = anova(mdl_frn_1b, 'DFMethod','Satterthwaite');
disp(anova_frn_1b);

mdl_rewp_1b = fitlme(gt_rev, ...
    [RewP_DV ' ~ block_type * correct_cat * reversal_transition_stage + ' ...
     '(1 + correct_cat | subj_id)'], ...
    'FitMethod','REML','DummyVarCoding','effects');
disp(mdl_rewp_1b.Coefficients);
anova_rewp_1b = anova(mdl_rewp_1b, 'DFMethod','Satterthwaite');
disp(anova_rewp_1b);

% -------------------------------------------------------------------------
% RQ1C: Past uncertainty (history-of-uncertainty effect)
% Theory: if subjects have been trained in a P block immediately before,
%   their outcome PE signal may be globally dampened even in the subsequent
%   D block, because they have learned to discount feedback. This is the
%   neural signature of second-order uncertainty carry-over.
%   Ref: Bland & Schaefer (2012) Front Neurosci.
% n_prev_P is also entered as a continuous covariate to capture cumulative
%   exposure, not just the immediately preceding block.
% -------------------------------------------------------------------------
fprintf('\n--- RQ1C: Past uncertainty model ---\n');
mdl_frn_1c = fitlme(gt_true, ...
    [FRN_DV ' ~ correct_cat * stage * prev_block_type + ' ...
     'block_type + n_prev_P + ' ...
     '(1 + correct_cat | subj_id)'], ...
    'FitMethod','REML','DummyVarCoding','effects');
disp(mdl_frn_1c.Coefficients);
anova_frn_1c = anova(mdl_frn_1c, 'DFMethod','Satterthwaite');
disp(anova_frn_1c);

% P300 past uncertainty
mdl_p3_1c = fitlme(gt_true, ...
    ['P300_amp_z ~ correct_cat * stage * prev_block_type + ' ...
     'block_type + n_prev_P + ' ...
     '(1 + correct_cat | subj_id)'], ...
    'FitMethod','REML','DummyVarCoding','effects');
disp(mdl_p3_1c.Coefficients);
anova_p3_1c = anova(mdl_p3_1c, 'DFMethod','Satterthwaite');
disp(anova_p3_1c);

% -------------------------------------------------------------------------
% RQ1D: False feedback dissociation (P blocks only)
% Theory: the FRN is sensitive to outcome valence, not ground-truth
%   correctness — the brain generates an error signal based on SHOWN
%   feedback. On false-negative trials (told incorrect when correct),
%   we expect a large FRN despite the subject being objectively correct.
%   This tests the Holroyd & Coles (2002) prediction directly.
%   Cavanagh et al. (2012) J Neurosci showed that FRN amplitude scales
%   with outcome certainty, so on false-feedback trials where subjects
%   are unsure, the FRN may be SMALLER than on true-negative trials.
% -------------------------------------------------------------------------
gt_p = gt(string(gt.block_type)=='P' & ~isnan(gt.(FRN_DV)), :);
gt_p.stage        = categorical(string(gt_p.stage), {'LN','LE','RN','RE'}, 'Ordinal',false);
gt_p.correct_cat  = categorical(gt_p.correct, [0 1], {'Incorrect','Correct'});
gt_p.false_fb_cat = categorical(gt_p.false_fb, [0 1], {'TrueFB','FalseFB'});

fprintf('\n--- RQ1D: False feedback model (P blocks) ---\n');
mdl_frn_1d = fitlme(gt_p, ...
    [FRN_DV ' ~ false_fb_cat * correct_cat * stage + ' ...
     '(1 + correct_cat | subj_id)'], ...
    'FitMethod','REML','DummyVarCoding','effects');
disp(mdl_frn_1d.Coefficients);
anova_frn_1d = anova(mdl_frn_1d, 'DFMethod','Satterthwaite');
disp(anova_frn_1d);

% =========================================================================
%% SECTION 2 PLOTS
% =========================================================================

fig2 = figure('Position',[50 50 1600 900]);
sgtitle('RQ1: FRN/RewP — outcome signal under short-term uncertainty', ...
    'FontSize',11,'FontWeight','bold');

% 2.1 Cell means: FRN by block_type × correct × stage
ax = subplot(2,4,1);
plot_cell_means_lineplot(ax, gt_true, FRN_DV, 'stage', 'block_type', ...
    'correct_cat', {'D','P'}, {'Incorrect','Correct'}, ...
    'FRN — stage × block type × outcome', true);

% 2.2 Cell means: RewP
ax = subplot(2,4,2);
plot_cell_means_lineplot(ax, gt_true, RewP_DV, 'stage', 'block_type', ...
    'correct_cat', {'D','P'}, {'Incorrect','Correct'}, ...
    'RewP — stage × block type × outcome', false);

% 2.3 LE vs RN reversal contrast: FRN
ax = subplot(2,4,3);
plot_rev_contrast(ax, gt_rev, FRN_DV, 'block_type', {'D','P'}, ...
    'FRN: LE vs RN × block type');

% 2.4 LE vs RN: P300
ax = subplot(2,4,4);
plot_rev_contrast(ax, gt_rev, 'P300_amp_z', 'block_type', {'D','P'}, ...
    'P300: LE vs RN × block type');

% 2.5 Past uncertainty: FRN by stage × prev_block_type
ax = subplot(2,4,5);
plot_cell_means_lineplot(ax, gt_true, FRN_DV, 'stage', 'prev_block_type', ...
    [], {'first','D','P'}, {}, ...
    'FRN — stage × prev block type (past uncertainty)', true);

% 2.6 Past uncertainty: P300
ax = subplot(2,4,6);
plot_cell_means_lineplot(ax, gt_true, 'P300_amp_z', 'stage', 'prev_block_type', ...
    [], {'first','D','P'}, {}, ...
    'P300 — stage × prev block type', false);

% 2.7 False feedback: FRN in P blocks
ax = subplot(2,4,7);
plot_cell_means_lineplot(ax, gt_p, FRN_DV, 'stage', 'false_fb_cat', ...
    'correct_cat', {'TrueFB','FalseFB'}, {'Incorrect','Correct'}, ...
    'FRN — false vs true feedback (P blocks)', true);

% 2.8 FRN exclusion rate by block_type (diagnostic)
ax = subplot(2,4,8);
if ismember('FRN_excluded', gt.Properties.VariableNames)
    plot_excl_rate_by_condition(ax, gt, 'FRN_excluded', 'block_type', ...
        'stage', 'FRN_excluded rate by block type & stage');
end

exportgraphics(fig2, fullfile(figure_output_folder,'RQ1_FRN_RewP_uncertainty.pdf'), ...
    'ContentType','vector');

% =========================================================================
%% SECTION 3 — RQ2: FRN ~ CONFIDENCE (PRECISION WEIGHTING)
% =========================================================================
% TRIAL ORDER: stimulus → response → CONFIDENCE → feedback.
% Confidence is a PRE-OUTCOME belief state. It indexes how certain the
% subject is about their response before seeing the outcome. This makes
% it a direct operationalisation of the precision weight on the current
% prediction — exactly the quantity that should scale the PE signal.
%
% Meyniel et al. (2015) PLoS Comp Bio showed that humans weight belief
% updates by their pre-outcome certainty. The FRN should therefore be
% LARGER (more negative) on high-confidence incorrect trials than on
% low-confidence incorrect trials, because a confident prediction that is
% then violated produces a larger PE than an uncertain prediction that is
% violated (Boldt & Yeung, 2015, J Neurosci 35:2058).
%
% This gives three testable predictions, all framed in terms of
% pre-outcome confidence:
%
%  RQ2A (current uncertainty): the confidence × incorrect → FRN slope
%    should be STEEPER in D blocks than P blocks. In D blocks, feedback
%    is reliable, so a high-confidence correct prediction that then gets
%    negative feedback is genuinely informative — the PE should be large.
%    In P blocks, even a confident subject has reason to doubt the
%    feedback, so the brain may down-weight the error signal
%    (Cavanagh et al., 2010).
%
%  RQ2B (past uncertainty): prior exposure to P blocks may globally
%    flatten the confidence-FRN coupling, because the subject has learned
%    that their certainty is not a reliable guide to feedback validity.
%
%  RQ2C (false feedback, P blocks): the key dissociation is between
%    high-confidence-correct trials that receive NEGATIVE feedback
%    (false-neg) versus low-confidence-incorrect trials that receive
%    NEGATIVE feedback (true-neg). On false-neg trials, the subject was
%    confident AND correct, so pre-outcome certainty was high — the
%    confidence-feedback mismatch is maximal. The FRN on false-neg trials
%    should therefore be LARGER for subjects with higher pre-outcome
%    confidence on those trials (conf_surprise_index from Section 1D).
%    This is NOT about what the subject thinks AFTER seeing feedback —
%    it is entirely determined by the confidence submitted beforehand.

fprintf('\n=== SECTION 3: RQ2 — Confidence × FRN ===\n');

gt_conf = gt(~isnan(gt.conf_z) & ~isnan(gt.(FRN_DV)), :);
gt_conf.block_type = categorical(string(gt_conf.block_type), {'D','P'});
gt_conf.stage      = categorical(string(gt_conf.stage), {'LN','LE','RN','RE'}, 'Ordinal',false);
gt_conf.correct_cat = categorical(gt_conf.correct, [0 1], {'Incorrect','Correct'});
gt_conf.prev_block_type = categorical(string(gt_conf.prev_block_type));
gt_conf.false_fb_cat    = categorical(gt_conf.false_fb, [0 1], {'TrueFB','FalseFB'});

gt_conf_true = gt_conf(~gt_conf.false_fb, :);

% RQ2A: current uncertainty — does confidence × correct interaction on FRN
% differ between D and P blocks? Prediction: steeper negative slope
% (conf × incorrect → large FRN) in D blocks than P blocks, because
% feedback is more reliable in D.
fprintf('\n--- RQ2A: FRN ~ conf × block_type × correct ---\n');
mdl_conf_2a = fitlme(gt_conf_true, ...
    [FRN_DV ' ~ conf_z * block_type * correct_cat + stage + ' ...
     '(1 + conf_z | subj_id)'], ...
    'FitMethod','REML','DummyVarCoding','effects');
disp(mdl_conf_2a.Coefficients);
anova_conf_2a = anova(mdl_conf_2a, 'DFMethod','Satterthwaite');
disp(anova_conf_2a);

% RQ2B: past uncertainty — does having experienced P blocks before change
% the confidence-FRN relationship? Prediction: after P exposure,
% confidence may be less tightly coupled to FRN because subjects have
% learned to distrust their own certainty.
fprintf('\n--- RQ2B: FRN ~ conf × prev_block_type × correct ---\n');
mdl_conf_2b = fitlme(gt_conf_true, ...
    [FRN_DV ' ~ conf_z * correct_cat * prev_block_type + ' ...
     'block_type + stage + n_prev_P + ' ...
     '(1 + conf_z | subj_id)'], ...
    'FitMethod','REML','DummyVarCoding','effects');
disp(mdl_conf_2b.Coefficients);
anova_conf_2b = anova(mdl_conf_2b, 'DFMethod','Satterthwaite');
disp(anova_conf_2b);

% RQ2C: false feedback (P blocks only)
% Here conf_z is the pre-outcome confidence on the current trial.
% The key interaction is conf_z × false_fb_cat: on false-negative trials
% (correct but told incorrect), a subject who entered that trial with
% HIGH pre-outcome confidence has a larger confidence-feedback mismatch.
% The FRN should therefore be larger for high-confidence false-neg trials
% than for low-confidence false-neg trials (conf_z × FalseFB slope > 0
% in the negative direction). On true-negative trials, the same subject
% was actually wrong, so their pre-outcome confidence may have been lower
% to begin with — the slope should be shallower.
gt_conf_p = gt_conf(string(gt_conf.block_type)=='P', :);
fprintf('\n--- RQ2C: FRN ~ conf × false_fb × correct (P blocks) ---\n');
mdl_conf_2c = fitlme(gt_conf_p, ...
    [FRN_DV ' ~ conf_z * false_fb_cat * correct_cat + stage + ' ...
     '(1 + conf_z | subj_id)'], ...
    'FitMethod','REML','DummyVarCoding','effects');
disp(mdl_conf_2c.Coefficients);
anova_conf_2c = anova(mdl_conf_2c, 'DFMethod','Satterthwaite');
disp(anova_conf_2c);

% ---- Plots ----
fig3 = figure('Position',[50 50 1400 600]);
sgtitle('RQ2: Confidence × FRN (precision weighting of error signal)','FontSize',11);

ax = subplot(1,3,1);
plot_conf_regression_scatter(ax, gt_conf_true, FRN_DV, 'conf_z', ...
    'correct_cat', {'Incorrect','Correct'}, 'block_type', {'D','P'}, ...
    true, 'FRN ~ confidence: D blocks', 'D');

ax = subplot(1,3,2);
plot_conf_regression_scatter(ax, gt_conf_true, FRN_DV, 'conf_z', ...
    'correct_cat', {'Incorrect','Correct'}, 'block_type', {'D','P'}, ...
    true, 'FRN ~ confidence: P blocks (true FB only)', 'P');

ax = subplot(1,3,3);
plot_conf_regression_scatter(ax, gt_conf_p, FRN_DV, 'conf_z', ...
    'false_fb_cat', {'TrueFB','FalseFB'}, [], {}, ...
    false, 'FRN ~ confidence: P blocks, true vs false FB', []);

exportgraphics(fig3, fullfile(figure_output_folder,'RQ2_Confidence_FRN.pdf'), ...
    'ContentType','vector');

% =========================================================================
%% SECTION 4 — RQ3: FRONTAL THETA ~ STAGE × BLOCK TYPE (INCORRECT TRIALS)
% =========================================================================
% Frontal theta indexes conflict monitoring (Cavanagh & Frank, 2014).
% On incorrect trials, especially post-reversal (RN), the perseverative
% error should drive strong theta if conflict resolution is needed.
% In P blocks, where errors are ambiguous, theta may be attenuated relative
% to D blocks (because the error might just be noise, not true conflict).
% Key test: theta × block_type at RN stage — if P-block reversal errors
% produce LESS theta than D-block reversal errors, subjects are failing to
% register genuine rule conflict under uncertainty.

fprintf('\n=== SECTION 4: RQ3 — Frontal theta ===\n');

gt_th = gt(gt.correct == 0 & ~gt.false_fb & ~isnan(gt.Theta_amp_z), :);
gt_th.block_type      = categorical(string(gt_th.block_type),  {'D','P'});
gt_th.stage           = categorical(string(gt_th.stage),       {'LN','LE','RN','RE'}, 'Ordinal',false);
gt_th.prev_block_type = categorical(string(gt_th.prev_block_type));

% RQ3A: stage × block_type
fprintf('\n--- RQ3A: Theta ~ stage × block_type (incorrect, true FB) ---\n');
mdl_theta_3a = fitlme(gt_th, ...
    ['Theta_amp_z ~ stage * block_type + ' ...
     '(1 + block_type | subj_id)'], ...
    'FitMethod','REML','DummyVarCoding','effects');
disp(mdl_theta_3a.Coefficients);
anova_theta_3a = anova(mdl_theta_3a, 'DFMethod','Satterthwaite');
disp(anova_theta_3a);

% RQ3B: LE vs RN reversal contrast
gt_th_rev = gt_th(gt_th.stage=='LE' | gt_th.stage=='RN', :);
gt_th_rev.reversal_transition_stage = categorical(string(gt_th_rev.stage), {'LE','RN'});
fprintf('\n--- RQ3B: Theta ~ reversal_stage × block_type ---\n');
mdl_theta_3b = fitlme(gt_th_rev, ...
    ['Theta_amp_z ~ reversal_transition_stage * block_type + ' ...
     '(1 + block_type | subj_id)'], ...
    'FitMethod','REML','DummyVarCoding','effects');
disp(mdl_theta_3b.Coefficients);
anova_theta_3b = anova(mdl_theta_3b, 'DFMethod','Satterthwaite');
disp(anova_theta_3b);

% RQ3C: past uncertainty
fprintf('\n--- RQ3C: Theta ~ stage × prev_block_type ---\n');
mdl_theta_3c = fitlme(gt_th, ...
    ['Theta_amp_z ~ stage * prev_block_type + block_type + n_prev_P + ' ...
     '(1 + block_type | subj_id)'], ...
    'FitMethod','REML','DummyVarCoding','effects');
disp(mdl_theta_3c.Coefficients);
anova_theta_3c = anova(mdl_theta_3c, 'DFMethod','Satterthwaite');
disp(anova_theta_3c);

% ---- Plots ----
fig4 = figure('Position',[50 50 1400 600]);
sgtitle('RQ3: Frontal theta (incorrect, true-FB trials) — uncertainty × conflict','FontSize',11);

ax = subplot(1,3,1);
plot_cell_means_lineplot(ax, gt_th, 'Theta_amp_z', 'stage', 'block_type', ...
    [], {'D','P'}, {}, 'Theta ~ stage × block type', false);

ax = subplot(1,3,2);
plot_rev_contrast(ax, gt_th_rev, 'Theta_amp_z', 'block_type', {'D','P'}, ...
    'Theta: LE vs RN × block type');

ax = subplot(1,3,3);
plot_cell_means_lineplot(ax, gt_th, 'Theta_amp_z', 'stage', 'prev_block_type', ...
    [], {'first','D','P'}, {}, 'Theta ~ stage × prev block type (history)', false);

exportgraphics(fig4, fullfile(figure_output_folder,'RQ3_Theta_uncertainty.pdf'), ...
    'ContentType','vector');

% =========================================================================
%% SECTION 5 — RQ4/5: PLV PATHWAYS ~ STAGE × BLOCK TYPE
% =========================================================================
% Fronto-parietal PLV (PLV_fp): indexes coherence between frontal control
% regions and parietal outcome-processing areas. In the context of reversal,
% this should peak at RN when top-down updating is most needed.
% Fronto-somatosensory PLV (PLV_fs): specific to this Braille/tactile task;
% indexes how well the motor/sensory response system is being coordinated
% with frontal control signals. May be differentially engaged by Go/NoGo
% reversals where the motor response itself must be suppressed.
%
% PLV is computed at the condition level (condition-level averages), so
% the unit of analysis here is subject × stage × block_type condition.
% This is the appropriate level given the sliding-window PLV estimation
% (min trials constraint means single-trial PLV is noisy).

fprintf('\n=== SECTION 5: RQ4/5 — PLV pathways ===\n');

% Condition-level aggregation
plv_cond = build_plv_condition_table(gt, subjects);

if ~isempty(plv_cond)
    plv_fp_true = plv_cond(plv_cond.pathway=='fp' & ...
                           plv_cond.false_fb_cat=='TrueFB' & ...
                           ~isnan(plv_cond.PLV_z), :);

    % RQ4A: full stage model
    if height(plv_fp_true) > 20
        fprintf('\n--- RQ4A: FP PLV ~ stage × block_type ---\n');
        mdl_plv_4a = fitlme(plv_fp_true, ...
            ['PLV_z ~ stage * block_type + n_rows_contributing + ' ...
             '(1 + block_type | subj_id)'], ...
            'FitMethod','REML','DummyVarCoding','effects');
        disp(mdl_plv_4a.Coefficients);
        anova_plv_4a = anova(mdl_plv_4a, 'DFMethod','Satterthwaite');
        disp(anova_plv_4a);

        % RQ4B: LE vs RN
        plv_fp_rev = plv_fp_true(plv_fp_true.stage=='LE' | plv_fp_true.stage=='RN', :);
        plv_fp_rev.reversal_transition_stage = categorical(string(plv_fp_rev.stage), {'LE','RN'});
        if height(plv_fp_rev) > 10
            fprintf('\n--- RQ4B: FP PLV ~ LE/RN × block_type ---\n');
            mdl_plv_4b = fitlme(plv_fp_rev, ...
                ['PLV_z ~ reversal_transition_stage * block_type + n_rows_contributing + ' ...
                 '(1 + block_type | subj_id)'], ...
                'FitMethod','REML','DummyVarCoding','effects');
            disp(mdl_plv_4b.Coefficients);
        end

        % RQ4C: past uncertainty
        fprintf('\n--- RQ4C: FP PLV ~ stage × prev_block_type ---\n');
        mdl_plv_4c = fitlme(plv_fp_true, ...
            ['PLV_z ~ stage * prev_block_type + block_type + n_prev_P + n_rows_contributing + ' ...
             '(1 + block_type | subj_id)'], ...
            'FitMethod','REML','DummyVarCoding','effects');
        disp(mdl_plv_4c.Coefficients);

        % RQ5: pathway comparison (fp vs fs)
        plv_long = plv_cond(plv_cond.false_fb_cat=='TrueFB' & ...
                            ~isnan(plv_cond.PLV_z_within_pathway), :);
        if height(plv_long) > 20
            fprintf('\n--- RQ5: Pathway × stage × block_type ---\n');
            mdl_pathway = fitlme(plv_long, ...
                ['PLV_z_within_pathway ~ pathway * stage * block_type + n_rows_contributing + ' ...
                 '(1 + pathway | subj_id)'], ...
                'FitMethod','REML','DummyVarCoding','effects');
            disp(mdl_pathway.Coefficients);
            anova_pathway = anova(mdl_pathway, 'DFMethod','Satterthwaite');
            disp(anova_pathway);
        end
    end

    % ---- PLV Plots ----
    fig5 = figure('Position',[50 50 1400 600]);
    sgtitle('RQ4/5: Fronto-parietal vs fronto-somatosensory PLV','FontSize',11);

    if height(plv_fp_true) > 10
        ax = subplot(1,3,1);
        plot_cell_means_lineplot(ax, plv_fp_true, 'PLV_z', 'stage', ...
            'block_type', [], {'D','P'}, {}, ...
            'FP PLV ~ stage × block type', false);

        if height(plv_fp_rev) > 5
            ax = subplot(1,3,2);
            plot_rev_contrast(ax, plv_fp_rev, 'PLV_z', 'block_type', {'D','P'}, ...
                'FP PLV: LE vs RN × block type');
        end
    end

    plv_fp_past = plv_fp_true;
    if height(plv_fp_past) > 10
        ax = subplot(1,3,3);
        plot_cell_means_lineplot(ax, plv_fp_past, 'PLV_z', 'stage', ...
            'prev_block_type', [], {'first','D','P'}, {}, ...
            'FP PLV ~ stage × prev block type', false);
    end

    exportgraphics(fig5, fullfile(figure_output_folder,'RQ4_5_PLV_pathways.pdf'), ...
        'ContentType','vector');
end

% =========================================================================
%% SECTION 6 — CROSS-LEVEL: NEURAL ~ BEHAVIOURAL INDIVIDUAL DIFFERENCES
% =========================================================================
% For each subject, extract neural index (mean of EEG feature in the
% contrast condition of interest) and correlate with behavioural ID index.
%
% KEY PREDICTIONS (theory-driven, so these are confirmatory):
%
%  a) Larger FRN difference (D-incorrect minus P-incorrect) predicts
%     larger ls_ratio (context-sensitive lose-switch).
%     Rationale: subjects whose error signal discriminates D from P
%     environments are also those who adjust their switching strategy.
%
%  b) Larger theta at RN-incorrect (especially D blocks) predicts
%     smaller trials_to_crit after reversal.
%     Rationale: strong conflict signal at reversal onset drives faster
%     rule updating (Cavanagh & Frank, 2014).
%
%  c) alpha_ratio (prob/det) predicts PLV_fp at RN stage.
%     Rationale: subjects with more volatile internal models show more
%     fronto-parietal coupling at the reversal point, as they need more
%     top-down updating.
%
%  d) conf_surprise_index predicts FRN amplitude on false-negative trials.
%     Rationale: conf_surprise_index = mean pre-outcome confidence on
%     false-neg trials minus mean pre-outcome confidence on true-neg
%     trials (both computed BEFORE feedback is shown). A subject with a
%     high index enters false-neg trials more confidently than true-neg
%     ones, meaning the confidence-feedback mismatch is larger on those
%     trials. The FRN — which indexes the brain's detection of the
%     outcome surprise — should therefore be larger (more negative) for
%     high-index subjects on false-neg trials specifically.
%     Ref: Boldt & Yeung (2015) J Neurosci; Meyniel et al. (2015).

fprintf('\n=== SECTION 6: Neural ~ behavioural individual differences ===\n');

% Build per-subject neural contrast indices
subj_neural = table();
subj_neural.subj_id = categorical(subjects);

for si = 1:numel(subjects)
    sn = subjects{si};

    % FRN contrast: D-incorrect vs P-incorrect (true FB only)
    frn_d_inc = mean(gt.(FRN_DV)(gt.subj_id==sn & string(gt.block_type)=='D' & ...
        gt.correct==0 & ~gt.false_fb), 'omitnan');
    frn_p_inc = mean(gt.(FRN_DV)(gt.subj_id==sn & string(gt.block_type)=='P' & ...
        gt.correct==0 & ~gt.false_fb), 'omitnan');
    subj_neural.FRN_D_inc(si) = frn_d_inc;
    subj_neural.FRN_P_inc(si) = frn_p_inc;
    subj_neural.FRN_context_diff(si) = frn_d_inc - frn_p_inc;  % larger = more discriminating

    % FRN on false-negative vs true-negative
    frn_false_neg = mean(gt.(FRN_DV)(gt.subj_id==sn & string(gt.block_type)=='P' & ...
        gt.correct==1 & gt.false_fb==1), 'omitnan');
    frn_true_neg  = mean(gt.(FRN_DV)(gt.subj_id==sn & string(gt.block_type)=='P' & ...
        gt.correct==0 & gt.false_fb==0), 'omitnan');
    subj_neural.FRN_false_neg(si)  = frn_false_neg;
    subj_neural.FRN_true_neg(si)   = frn_true_neg;
    subj_neural.FRN_fb_diff(si)    = frn_false_neg - frn_true_neg;

    % Theta at RN-incorrect: conflict at reversal
    theta_RN_D = mean(gt.Theta_amp_z(gt.subj_id==sn & gt.stage=='RN' & ...
        string(gt.block_type)=='D' & gt.correct==0 & ~gt.false_fb), 'omitnan');
    theta_RN_P = mean(gt.Theta_amp_z(gt.subj_id==sn & gt.stage=='RN' & ...
        string(gt.block_type)=='P' & gt.correct==0 & ~gt.false_fb), 'omitnan');
    subj_neural.Theta_RN_D(si) = theta_RN_D;
    subj_neural.Theta_RN_P(si) = theta_RN_P;

    % P300 at RN (context updating)
    p300_RN_D = mean(gt.P300_amp_z(gt.subj_id==sn & gt.stage=='RN' & ...
        string(gt.block_type)=='D'), 'omitnan');
    p300_RN_P = mean(gt.P300_amp_z(gt.subj_id==sn & gt.stage=='RN' & ...
        string(gt.block_type)=='P'), 'omitnan');
    subj_neural.P300_RN_D(si) = p300_RN_D;
    subj_neural.P300_RN_P(si) = p300_RN_P;
    subj_neural.P300_RN_context_diff(si) = p300_RN_D - p300_RN_P;

    % PLV_fp at RN
    plv_RN_D = mean(gt.PLV_fp_z(gt.subj_id==sn & gt.stage=='RN' & ...
        string(gt.block_type)=='D' & ~gt.false_fb), 'omitnan');
    plv_RN_P = mean(gt.PLV_fp_z(gt.subj_id==sn & gt.stage=='RN' & ...
        string(gt.block_type)=='P' & ~gt.false_fb), 'omitnan');
    subj_neural.PLV_fp_RN_D(si) = plv_RN_D;
    subj_neural.PLV_fp_RN_P(si) = plv_RN_P;
end

% Join neural and behavioural ID tables
id_combined = innerjoin(subj_neural, id_table, 'Keys','subj_id');

% ---- Cross-level correlation plots ----
fig6 = figure('Position',[50 50 1600 900]);
sgtitle('Neural ~ behavioural individual differences', ...
    'FontSize',11,'FontWeight','bold');

cross_pairs = {
    % { neural_x, beh_y, prediction_direction, panel_label }
    'FRN_context_diff', 'ls_ratio', +1, 'FRN D-P contrast vs lose-switch ratio';
    'Theta_RN_D',       'trials_to_crit_D', -1, 'Theta at RN (D) vs trials to criterion (D)';
    'P300_RN_D',        'rev_AUC_D', +1, 'P300 at RN (D) vs post-rev AUC (D)';
    'P300_RN_context_diff', 'rev_AUC_ratio', +1, 'P300 context diff vs AUC ratio';
};
if has_conf
    cross_pairs(end+1,:) = {'FRN_false_neg', 'conf_surprise_index', -1, ...
        'FRN on false-neg trials vs pre-outcome confidence surprise index'};
end
if has_RL
    cross_pairs(end+1,:) = {'PLV_fp_RN_D', 'alpha_ratio', +1, ...
        'FP PLV at RN (D) vs alpha ratio (prob/det)'};
    cross_pairs(end+1,:) = {'Theta_RN_D', 'delta_alpha', +1, ...
        'Theta at RN (D) vs delta-alpha'};
end

n_panels = min(numel(cross_pairs(:,1)), 8);
for pi = 1:n_panels
    ax = subplot(2, 4, pi);
    x_var  = cross_pairs{pi,1};
    y_var  = cross_pairs{pi,2};
    lbl    = cross_pairs{pi,4};
    if ~ismember(x_var, id_combined.Properties.VariableNames) || ...
       ~ismember(y_var, id_combined.Properties.VariableNames)
        continue;
    end
    x_dat = id_combined.(x_var);
    y_dat = id_combined.(y_var);
    cohort_col = id_combined.cohort;
    plot_corr_scatter(ax, x_dat, y_dat, cohort_col, ...
        strrep(x_var,'_','\_'), strrep(y_var,'_','\_'), lbl);
end

exportgraphics(fig6, fullfile(figure_output_folder,'S6_Neural_Behav_ID.pdf'), ...
    'ContentType','vector');

% =========================================================================
%% SECTION 7 — RL MODEL LATENTS ~ NEURAL MARKERS
% =========================================================================
% Tests the Nassar/Bayesian model predictions directly against EEG signals.
% All single-trial models use LME to handle the hierarchical structure.
%
%  7A: PE (prediction error) ~ FRN_mean_amp
%      Theory: FRN should track the signed or unsigned PE signal generated
%      by the RL model. If the Holroyd-Coles AFC theory is correct, the
%      direction should be: negative PE (loss) → more negative FRN.
%      Ref: Walsh & Anderson (2012) Nat Rev Neurosci 13:590.
%
%  7B: PE ~ P300_amp (context updating)
%      Nassar's omega_trial (learning rate) maps to P300 in our framework.
%      PE × stage interaction tests whether post-reversal PEs drive more
%      P300 than equivalent PEs in stable epochs (as the context is now
%      uncertain, the brain allocates more updating resources).
%
%  7C: conf_weighted_PE ~ Theta_amp (precision-weighted conflict)
%      conf_weighted_PE from the behav_table multiplies PE by confidence,
%      testing whether the confidence-scaled mismatch drives theta more
%      than raw PE alone.

if has_PE
    fprintf('\n=== SECTION 7: RL latents ~ EEG markers ===\n');

    gt_pe = gt(~isnan(gt.PE_z) & ~gt.false_fb, :);
    gt_pe.block_type = categorical(string(gt_pe.block_type), {'D','P'});
    gt_pe.stage      = categorical(string(gt_pe.stage), {'LN','LE','RN','RE'}, 'Ordinal',false);

    % 7A: FRN ~ PE
    fprintf('\n--- 7A: FRN ~ PE_z × block_type × stage ---\n');
    mdl_pe_frn = fitlme(gt_pe, ...
        [FRN_DV ' ~ PE_z * block_type + stage + ' ...
         '(1 + PE_z | subj_id)'], ...
        'FitMethod','REML','DummyVarCoding','effects');
    disp(mdl_pe_frn.Coefficients);
    anova_pe_frn = anova(mdl_pe_frn, 'DFMethod','Satterthwaite');
    disp(anova_pe_frn);

    % 7B: P300 ~ PE × stage
    fprintf('\n--- 7B: P300 ~ PE_z × stage ---\n');
    mdl_pe_p300 = fitlme(gt_pe, ...
        ['P300_amp_z ~ PE_z * stage + block_type + ' ...
         '(1 + PE_z | subj_id)'], ...
        'FitMethod','REML','DummyVarCoding','effects');
    disp(mdl_pe_p300.Coefficients);
    anova_pe_p300 = anova(mdl_pe_p300, 'DFMethod','Satterthwaite');
    disp(anova_pe_p300);

    % 7C: Theta ~ conf_weighted_PE
    if ismember('conf_weighted_PE_z', gt_pe.Properties.VariableNames)
        gt_pe_th = gt_pe(~isnan(gt_pe.conf_weighted_PE_z) & ...
                         ~isnan(gt_pe.Theta_amp_z), :);
        fprintf('\n--- 7C: Theta ~ conf_weighted_PE_z × stage ---\n');
        mdl_cwpe_theta = fitlme(gt_pe_th, ...
            ['Theta_amp_z ~ conf_weighted_PE_z * stage + block_type + ' ...
             '(1 + conf_weighted_PE_z | subj_id)'], ...
            'FitMethod','REML','DummyVarCoding','effects');
        disp(mdl_cwpe_theta.Coefficients);
    end

    % ---- RL-Neural plots ----
    fig7 = figure('Position',[50 50 1400 500]);
    sgtitle('RL latents ~ EEG markers (Nassar model)','FontSize',11);

    ax = subplot(1,3,1);
    plot_pe_erp_scatter(ax, gt_pe, 'PE_z', FRN_DV, 'stage', ...
        'PE ~ FRN by reversal stage');

    ax = subplot(1,3,2);
    plot_pe_erp_scatter(ax, gt_pe, 'PE_z', 'P300_amp_z', 'stage', ...
        'PE ~ P300 by reversal stage');

    if has_RL
        ax = subplot(1,3,3);
        % Alpha_post vs P300 at RN: subjects with higher post-rev alpha
        % should show larger P300 at RN (both index belief updating)
        plot_corr_scatter(ax, id_combined.P300_RN_D, ...
            id_combined.delta_alpha, id_combined.cohort, ...
            'P300 at RN (D)', '\Delta\alpha (post-pre)', ...
            'P300 at reversal vs learning rate shift');
    end

    exportgraphics(fig7, fullfile(figure_output_folder,'S7_RL_neural_markers.pdf'), ...
        'ContentType','vector');
end

% =========================================================================
%% SAVE ALL MODELS AND TABLES
% =========================================================================

save(fullfile(figure_output_folder, 'ID_analysis_workspace_v1.mat'), ...
    'id_table', 'id_combined', 'gt_true', 'gt_rev', 'gt_th', ...
    'mdl_frn_1a', 'mdl_frn_1b', 'mdl_frn_1c', 'mdl_frn_1d', ...
    'mdl_rewp_1a', 'mdl_rewp_1b', ...
    'mdl_conf_2a', 'mdl_conf_2b', 'mdl_conf_2c', ...
    'mdl_theta_3a', 'mdl_theta_3b', 'mdl_theta_3c', ...
    '-v7.3');
fprintf('\nWorkspace saved to %s\n', figure_output_folder);

% Report summary
fprintf('\n=== ANALYSIS COMPLETE ===\n');
fprintf('Subjects: %d | Cohorts: %s\n', numel(subjects), ...
    strjoin(categories(gt.cohort),', '));
fprintf('Figures saved to: %s\n', figure_output_folder);

% =========================================================================
%% LOCAL HELPER FUNCTIONS
% =========================================================================

% -------------------------------------------------------------------------
function [gt, col_z] = ensure_z_column(gt, col_z, subjects)
% If col_z already exists (from v9b), leave it. If the raw column exists
% but _z doesn't, z-score it. Logs what happened.
    if ismember(col_z, gt.Properties.VariableNames)
        return;  % already there from v9b
    end
    raw_col = strrep(col_z, '_z', '');
    if ~ismember(raw_col, gt.Properties.VariableNames)
        warning('Neither %s nor %s found. Column will be NaN.', col_z, raw_col);
        gt.(col_z) = nan(height(gt), 1);
        return;
    end
    fprintf('  z-scoring %s → %s\n', raw_col, col_z);
    gt.(col_z) = nan(height(gt), 1);
    for si = 1:numel(subjects)
        m    = gt.subj_id == subjects{si};
        vals = gt.(raw_col)(m);
        mn   = mean(vals,'omitnan');
        sd   = std(vals,'omitnan');
        if sd > 0
            gt.(col_z)(m) = (vals - mn) / sd;
        end
    end
end

% -------------------------------------------------------------------------
function plv_cond = build_plv_condition_table(gt, subjects)
% Aggregate single-trial PLV columns to condition-level means.
% One row per subject × stage × block_type × false_fb_cat × correct_cat.
% PLV_z_within_pathway is z-scored within subject × pathway for RQ5.

    pathways   = {'fp','fs'};
    feat_names = {'PLV_fp_z','PLV_fs_z'};
    rows = {};

    for pp = 1:2
        feat = feat_names{pp};
        if ~ismember(feat, gt.Properties.VariableNames), continue; end
        gtmp = gt(~isnan(gt.(feat)), :);
        if isempty(gtmp), continue; end

        gtmp.subj_id      = categorical(string(gtmp.subj_id));
        gtmp.stage        = categorical(string(gtmp.stage), {'LN','LE','RN','RE'},'Ordinal',false);
        gtmp.block_type   = categorical(string(gtmp.block_type), {'D','P'});
        gtmp.false_fb_cat = categorical(gtmp.false_fb, [0 1], {'TrueFB','FalseFB'});
        gtmp.correct_cat  = categorical(gtmp.correct, [0 1], {'Incorrect','Correct'});
        gtmp.prev_block_type = categorical(string(gtmp.prev_block_type));

        [G, subj, stage, block_type, false_fb_cat, correct_cat, prev_block_type] = ...
            findgroups(gtmp.subj_id, gtmp.stage, gtmp.block_type, ...
                       gtmp.false_fb_cat, gtmp.correct_cat, gtmp.prev_block_type);

        T = table();
        T.subj_id         = subj;
        T.stage           = stage;
        T.block_type      = block_type;
        T.false_fb_cat    = false_fb_cat;
        T.correct_cat     = correct_cat;
        T.prev_block_type = prev_block_type;
        T.PLV_z           = splitapply(@(x) mean(x,'omitnan'), gtmp.(feat), G);
        T.n_rows_contributing = splitapply(@numel, gtmp.(feat), G);
        T.n_prev_P        = splitapply(@(x) mean(x,'omitnan'), gtmp.n_prev_P, G);
        T.pathway         = categorical(repmat({pathways{pp}}, height(T), 1), {'fp','fs'});
        rows{end+1} = T; %#ok<AGROW>
    end

    if isempty(rows)
        plv_cond = table();
        return;
    end

    plv_cond = vertcat(rows{:});

    % Within-subject×pathway z-score for pathway comparison (RQ5)
    plv_cond.PLV_z_within_pathway = nan(height(plv_cond), 1);
    [Gsp, ~, ~] = findgroups(plv_cond.subj_id, plv_cond.pathway);
    mu_v = splitapply(@(x)mean(x,'omitnan'), plv_cond.PLV_z, Gsp);
    sd_v = splitapply(@(x)std(x,'omitnan'),  plv_cond.PLV_z, Gsp);
    for gi = 1:max(Gsp)
        m = Gsp==gi;
        if sd_v(gi) > 0
            plv_cond.PLV_z_within_pathway(m) = (plv_cond.PLV_z(m) - mu_v(gi)) / sd_v(gi);
        end
    end
end

% -------------------------------------------------------------------------
function plot_cell_means_lineplot(ax, T, feat, xvar, linevar, filtervar, ...
    line_cats, filter_cats, ttl, reverse_y)
% General purpose line plot: xvar on x-axis, linevar as separate lines.
% filtervar/filter_cats: if supplied, only matching filtervar rows are used.
% reverse_y: EEG convention (negativity upward).

    if ~iscategorical(T.(xvar)),    T.(xvar)    = categorical(string(T.(xvar)));    end
    if ~iscategorical(T.(linevar)), T.(linevar) = categorical(string(T.(linevar))); end

    xcats = categories(T.(xvar));
    if ~isempty(line_cats)
        lcats = line_cats;
    else
        lcats = categories(T.(linevar));
    end

    COLORS = [0.15 0.45 0.70; 0.80 0.30 0.10; 0.12 0.62 0.47;
              0.85 0.40 0.00; 0.55 0.25 0.65; 0.10 0.10 0.10];
    LS = {'-o','-s','-^','-d','-v','-p'};

    hold(ax,'on');
    subjs = unique(T.subj_id);

    for li = 1:numel(lcats)
        % Filter to correct linevar level and filtervar subset if given
        if ~isempty(filtervar) && ~isempty(filter_cats)
            lm = T.(linevar)==lcats{li} & ismember(string(T.(filtervar)), filter_cats);
        else
            lm = T.(linevar)==lcats{li};
        end

        ms = nan(1, numel(xcats));
        se = nan(1, numel(xcats));

        for xi = 1:numel(xcats)
            m = lm & T.(xvar)==xcats{xi} & ~isnan(T.(feat));

            % Subject-level means to avoid trial-count weighting
            subj_means = nan(numel(subjs), 1);
            for si2 = 1:numel(subjs)
                sm2 = m & T.subj_id==subjs(si2);
                if any(sm2)
                    subj_means(si2) = mean(T.(feat)(sm2),'omitnan');
                end
            end
            ok = ~isnan(subj_means);
            if sum(ok) < 2, continue; end
            ms(xi) = mean(subj_means(ok));
            se(xi) = std(subj_means(ok)) / sqrt(sum(ok));
        end

        clr = COLORS(mod(li-1, size(COLORS,1))+1, :);
        ls  = LS{mod(li-1, numel(LS))+1};
        errorbar(ax, 1:numel(xcats), ms, se, ls, ...
            'Color', clr, 'LineWidth', 1.8, ...
            'MarkerFaceColor', clr, ...
            'DisplayName', strrep(lcats{li},'_','\_'));
    end

    set(ax, 'XTick',1:numel(xcats), 'XTickLabel', xcats, ...
        'XTickLabelRotation', 20, 'FontSize', 8);
    xlabel(ax, strrep(xvar,'_','\_'), 'FontSize', 8);
    ylabel(ax, strrep(feat,'_','\_'), 'FontSize', 8);
    title(ax, ttl, 'Interpreter','none', 'FontSize', 9);
    legend(ax, 'Box','off', 'FontSize', 7, 'Location','best');
    yline(ax, 0, '--k', 'HandleVisibility','off');
    if reverse_y, set(ax,'YDir','reverse'); end
end

% -------------------------------------------------------------------------
function plot_rev_contrast(ax, T_rev, feat, groupvar, group_cats, ttl)
% Two-point plot: LE (x=1) vs RN (x=2), lines per group_cats.

    if ~iscategorical(T_rev.(groupvar))
        T_rev.(groupvar) = categorical(string(T_rev.(groupvar)));
    end
    if ~iscategorical(T_rev.reversal_transition_stage)
        T_rev.reversal_transition_stage = categorical(string(T_rev.reversal_transition_stage));
    end

    COLORS = [0.15 0.45 0.70; 0.80 0.30 0.10; 0.12 0.62 0.47];
    subjs  = unique(T_rev.subj_id);
    hold(ax,'on');

    for gi = 1:numel(group_cats)
        gm = T_rev.(groupvar) == group_cats{gi};
        ms = nan(1,2); se = nan(1,2);

        for xi = 1:2
            xcat = {'LE','RN'};
            m = gm & T_rev.reversal_transition_stage == xcat{xi} & ~isnan(T_rev.(feat));
            subj_means = nan(numel(subjs),1);
            for si = 1:numel(subjs)
                sm2 = m & T_rev.subj_id==subjs(si);
                if any(sm2)
                    subj_means(si) = mean(T_rev.(feat)(sm2),'omitnan');
                end
            end
            ok = ~isnan(subj_means);
            if sum(ok) < 2, continue; end
            ms(xi) = mean(subj_means(ok));
            se(xi) = std(subj_means(ok)) / sqrt(sum(ok));
        end

        clr = COLORS(gi,:);
        errorbar(ax, [1 2], ms, se, '-o', 'Color',clr, 'LineWidth',1.8, ...
            'MarkerFaceColor',clr, 'DisplayName', group_cats{gi});
    end

    set(ax, 'XTick',[1 2], 'XTickLabel',{'LE','RN'}, 'FontSize',8);
    xlabel(ax,'Reversal stage','FontSize',8);
    ylabel(ax, strrep(feat,'_','\_'),'FontSize',8);
    title(ax, ttl,'Interpreter','none','FontSize',9);
    legend(ax,'Box','off','FontSize',7,'Location','best');
    yline(ax,0,'--k','HandleVisibility','off');
end

% -------------------------------------------------------------------------
function plot_conf_regression_scatter(ax, T, feat, xvar, splitvar, split_cats, ...
    filtervar, filter_vals, negate_y, ttl, filter_val)
% Scatter + regression line for confidence ~ FRN, split by a categorical.

    if ~isempty(filtervar) && ~isempty(filter_vals) && ~isempty(filter_val)
        T = T(string(T.(filtervar))==filter_val, :);
    end

    if ~iscategorical(T.(splitvar)), T.(splitvar) = categorical(string(T.(splitvar))); end
    if ~isempty(split_cats)
        cats = split_cats;
    else
        cats = categories(T.(splitvar));
    end

    COLORS = [0.80 0.20 0.20; 0.20 0.20 0.80; 0.20 0.70 0.20; 0.70 0.40 0.00];
    hold(ax,'on');

    for ci = 1:numel(cats)
        m = T.(splitvar)==cats{ci} & ~isnan(T.(xvar)) & ~isnan(T.(feat));
        x = T.(xvar)(m);
        y = T.(feat)(m);
        if numel(x) < 5, continue; end

        clr = COLORS(mod(ci-1,size(COLORS,1))+1,:);
        scatter(ax, x, y, 8, clr, 'filled', 'MarkerFaceAlpha',0.12, ...
            'HandleVisibility','off');

        p_fit = polyfit(x, y, 1);
        xi_v  = linspace(min(x), max(x), 80);
        [r_v, p_v] = corr(x, y, 'Rows','complete');
        plot(ax, xi_v, polyval(p_fit, xi_v), '-', 'Color', clr, ...
            'LineWidth', 2.0, 'DisplayName', ...
            sprintf('%s (r=%.2f, p=%.3f)', strrep(cats{ci},'_','\_'), r_v, p_v));
    end

    xlabel(ax, strrep(xvar,'_','\_'), 'FontSize',8);
    ylabel(ax, strrep(feat,'_','\_'), 'FontSize',8);
    title(ax, ttl, 'Interpreter','none', 'FontSize',9);
    legend(ax,'Box','off','FontSize',7,'Location','best');
    yline(ax,0,'--k','HandleVisibility','off');
    if negate_y, set(ax,'YDir','reverse'); end
end

% -------------------------------------------------------------------------
function plot_pe_erp_scatter(ax, T, pe_col, erp_col, stage_col, ttl)
% PE vs ERP scatter, one line per stage.

    stages = {'LN','LE','RN','RE'};
    COLORS = [0.12 0.62 0.47; 0.85 0.65 0.00; 0.80 0.27 0.13; 0.40 0.25 0.65];
    hold(ax,'on');

    for si = 1:4
        m = string(T.(stage_col))==stages{si} & ~isnan(T.(pe_col)) & ~isnan(T.(erp_col));
        x = T.(pe_col)(m);
        y = T.(erp_col)(m);
        if numel(x) < 5, continue; end
        clr = COLORS(si,:);
        scatter(ax, x, y, 6, clr, 'filled', 'MarkerFaceAlpha',0.10, ...
            'HandleVisibility','off');
        p_fit = polyfit(x, y, 1);
        xi_v = linspace(min(x), max(x), 60);
        [r_v, p_v] = corr(x, y,'Rows','complete');
        plot(ax, xi_v, polyval(p_fit, xi_v), '-', 'Color',clr, 'LineWidth',1.8, ...
            'DisplayName', sprintf('%s r=%.2f (p=%.3f)', stages{si}, r_v, p_v));
    end

    xlabel(ax, strrep(pe_col,'_','\_'),'FontSize',8);
    ylabel(ax, strrep(erp_col,'_','\_'),'FontSize',8);
    title(ax, ttl, 'Interpreter','none','FontSize',9);
    legend(ax,'Box','off','FontSize',7,'Location','best');
    yline(ax,0,'--k','HandleVisibility','off');
end

% -------------------------------------------------------------------------
function plot_corr_scatter(ax, x, y, cohort_cat, xlbl, ylbl, ttl)
% Subject-level scatter with cohort-coded colour + regression line.

    ok = ~isnan(x) & ~isnan(y);
    cohorts = categories(categorical(cohort_cat));
    COLORS  = [0.15 0.45 0.70; 0.80 0.30 0.10; 0.30 0.70 0.30];
    hold(ax,'on');

    for ci = 1:numel(cohorts)
        cm = ok & cohort_cat==cohorts{ci};
        if sum(cm) < 2, continue; end
        scatter(ax, x(cm), y(cm), 35, COLORS(ci,:), 'filled', ...
            'DisplayName', cohorts{ci});
    end

    if sum(ok) > 3
        p_fit = polyfit(x(ok), y(ok), 1);
        xi_v = linspace(min(x(ok)), max(x(ok)), 80);
        plot(ax, xi_v, polyval(p_fit, xi_v), 'k-', 'LineWidth',1.8, ...
            'HandleVisibility','off');
        [r_v, p_v] = corr(x(ok), y(ok));
        text(ax, 0.05, 0.93, sprintf('r=%.2f, p=%.3f', r_v, p_v), ...
            'Units','normalized','FontSize',8,'FontWeight','bold');
    end

    xlabel(ax, xlbl,'FontSize',8); ylabel(ax, ylbl,'FontSize',8);
    title(ax, ttl,'Interpreter','none','FontSize',9);
    legend(ax,'Box','off','FontSize',7,'Location','best');
end

% -------------------------------------------------------------------------
function plot_subject_scatter(ax, x, y, cohort_cat, xlbl, ylbl, ttl)
% Subject-level scatter, no regression. Used for AUC D vs P comparisons.

    ok = ~isnan(x) & ~isnan(y);
    cohorts = categories(categorical(cohort_cat));
    COLORS  = [0.15 0.45 0.70; 0.80 0.30 0.10; 0.30 0.70 0.30];
    hold(ax,'on');

    for ci = 1:numel(cohorts)
        cm = ok & cohort_cat==cohorts{ci};
        if ~any(cm), continue; end
        scatter(ax, x(cm), y(cm), 40, COLORS(ci,:), 'filled', ...
            'DisplayName', cohorts{ci});
    end

    xlabel(ax, xlbl,'FontSize',8); ylabel(ax, ylbl,'FontSize',8);
    title(ax, ttl,'Interpreter','none','FontSize',9);
    legend(ax,'Box','off','FontSize',7,'Location','best');
end

% -------------------------------------------------------------------------
function plot_id_histogram(ax, x, xlbl, ttl, cohort_cat)
% Histogram of an ID index, overlaid per cohort.

    cohorts = categories(categorical(cohort_cat));
    COLORS  = [0.15 0.45 0.70; 0.80 0.30 0.10; 0.30 0.70 0.30];
    hold(ax,'on');

    for ci = 1:numel(cohorts)
        cm = cohort_cat==cohorts{ci} & ~isnan(x);
        histogram(ax, x(cm), 10, 'FaceColor', COLORS(ci,:), ...
            'FaceAlpha', 0.55, 'DisplayName', cohorts{ci});
    end

    xlabel(ax, xlbl,'FontSize',8);
    ylabel(ax,'Count','FontSize',8);
    title(ax, ttl,'Interpreter','none','FontSize',9);
    legend(ax,'Box','off','FontSize',7);
end

% -------------------------------------------------------------------------
function plot_excl_rate_by_condition(ax, T, excl_col, groupvar, xvar, ttl)
% Bar chart of exclusion rate (e.g. FRN_excluded) by group × x.

    xcats = categories(T.(xvar));
    gcats = categories(T.(groupvar));
    COLORS = [0.15 0.45 0.70; 0.80 0.30 0.10];
    hold(ax,'on');

    x_base = 1:numel(xcats);
    w = 0.35;
    for gi = 1:numel(gcats)
        for xi = 1:numel(xcats)
            m = T.(groupvar)==gcats{gi} & T.(xvar)==xcats{xi};
            if ~any(m), continue; end
            rate = mean(T.(excl_col)(m),'omitnan');
            x_pos = x_base(xi) + (gi-1.5)*w;
            bar(ax, x_pos, rate, w*0.9, 'FaceColor', COLORS(gi,:), 'EdgeColor','none', ...
                'DisplayName', gcats{gi});
        end
    end

    set(ax,'XTick',x_base,'XTickLabel',xcats,'FontSize',8);
    xlabel(ax,'Stage','FontSize',8);
    ylabel(ax,'Exclusion rate','FontSize',8);
    title(ax, ttl,'Interpreter','none','FontSize',9);
    legend(ax, gcats,'Box','off','FontSize',7);
end

% -------------------------------------------------------------------------
function refline_identity(ax)
% Add a diagonal identity line to a scatter plot.
    xl = xlim(ax); yl = ylim(ax);
    mn = min([xl(1) yl(1)]); mx = max([xl(2) yl(2)]);
    plot(ax, [mn mx], [mn mx], 'k--', 'HandleVisibility','off');
end

% -------------------------------------------------------------------------
function cmap = redblue_colormap()
% Red-white-blue diverging colormap for correlation matrices.
    n = 64;
    r = [linspace(0.70, 1.00, n/2), linspace(1.00, 0.20, n/2)]';
    g = [linspace(0.10, 1.00, n/2), linspace(1.00, 0.20, n/2)]';
    b = [linspace(0.20, 1.00, n/2), linspace(1.00, 0.70, n/2)]';
    cmap = [r g b];
end