% ==========================================================================
% E_sequential_block_behaviour_plots.m
%
% STANDALONE BEHAVIOURAL PLOTTING SCRIPT
% Run AFTER:
%   A_extract_revaligned_alltrialdata.m  → all_trial_data, group_T
%   D_Bayesian_RL_models.m               → results struct (Nassar fits)
%
% ──────────────────────────────────────────────────────────────────────────
% SCIENTIFIC RATIONALE
% ──────────────────────────────────────────────────────────────────────────
% This script distinguishes two timescales of uncertainty that your design
% dissociates:
%
%   LONG-TERM (structural) uncertainty — volatility: the reversal embedded
%     in every block forces participants to detect a change-point and update
%     their policy. The hazard rate H in the Nassar (2010) model captures
%     the estimated probability that the current rule has flipped.
%     Reference: Nassar et al. (2010) J Neurosci; Behrens et al. (2007)
%     Nat Neurosci.
%
%   SHORT-TERM (local) uncertainty — noise: in Probabilistic blocks 20% of
%     feedback is false, generating trial-level stochasticity without a
%     change in the underlying rule. This is the "unexpected uncertainty"
%     of Yu & Dayan (2005) that should downweight individual outcomes
%     without triggering full rule revision. Reference: Yu & Dayan (2005)
%     Neuron.
%
%   INTERACTION: A preceding P block (noise history) modulates how the
%     learner interprets the next reversal signal. If noise raises H
%     globally (H_prob > H_det), the learner should be *faster* to detect a
%     genuine reversal but also more prone to spurious resets. This is the
%     key individual-differences question: n_prev_P indexes cumulative noise
%     exposure at the time of each reversal.
%
% STIMULUS TYPES (2×2 dimensional shift):
%   Switched  (stims 2,3): Go↔NoGo assignment reverses at reversal.
%     Only these produce a genuine change-point signal for the Nassar ω.
%   Maintained (stims 1,4): assignment unchanged. Accuracy on maintained
%     stimuli is a pure measure of belief stability / response competition.
%     Reference: Collins & Frank (2013) Psych Rev for stimulus-specific
%     learning; Cavanagh & Frank (2014) TICS for conflict monitoring.
%
% n_prev_P DEFINITION (computed here, not yet in upstream scripts):
%   For each block b of subject s, n_prev_P(b) = number of blocks with
%   block_type == 'P' in blocks 1…(b-1) for that subject. This is a
%   clean within-subject running count of cumulative noise exposure.
%   It is the behavioural counterpart of the neural moderator used in the
%   LME models of C_individual_differences_uncertainty_v1.m.
%
% ──────────────────────────────────────────────────────────────────────────
% FIGURES PRODUCED
% ──────────────────────────────────────────────────────────────────────────
%  Fig E1  — Sequential block learning curves (accuracy × block number,
%             split D/P, with reversal-aligned inset)
%  Fig E2  — Reversal cost trajectory across blocks (LE→RN accuracy drop
%             as a function of block number and D/P type)
%  Fig E3  — Switched vs maintained stimulus accuracy, reversal-aligned
%             (per block type, with pre/post comparison)
%  Fig E4  — Switched vs maintained: confidence and RT profiles
%  Fig E5  — n_prev_P: reversal-aligned accuracy and confidence as a
%             function of cumulative P-block exposure (0 / 1 / 2+ bins)
%  Fig E6  — n_prev_P: reversal cost and recovery (ΔAcc, slope) as a
%             function of n_prev_P (scatter + fitted line per subject)
%  Fig E7  — n_prev_P × block_type interaction on reversal cost
%             (2×3 grid — complements planned LME models)
%  Fig E8  — Switch-stimulus accuracy: first encounter vs later encounters
%             post-reversal, split by n_prev_P (belief updating speed)
%  Fig E9  — EEG placeholder: predicted FRN / P300 pattern by switch-type
%             and n_prev_P (schematic only — fills when EEG table loaded)
%
% ──────────────────────────────────────────────────────────────────────────
% EXPECTED WORKSPACE VARIABLES (loaded automatically if absent):
%   group_T        — long-format trial table (A_extract)
%   all_trial_data — per-subject struct (A_extract)
%   results        — Nassar fit struct (D_Bayesian_RL_models) [optional]
% ==========================================================================

close all;

% =========================================================================
%% §0  PATHS & DATA LOADING
% =========================================================================
remote = 0;
switch remote
    case 1
        base_path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch/Salient mod switch KH';
    case 0
        base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH';
    case 2
        base_path = 'Z:\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH';
end
data_path = fullfile(base_path, 'Data');
eeg_path  = fullfile(base_path, 'Results', 'EEG analysis', 'Epoched_data');
outpath   = fullfile(base_path, 'Results', 'Behav results', 'Sequential_Block_Figures');
if ~exist(outpath, 'dir'), mkdir(outpath); end

% ── Load behavioural table ────────────────────────────────────────────────
if ~exist('group_T','var') || isempty(group_T)
    fprintf('Loading behav_table.mat...\n');
    load(fullfile(data_path, 'behav_table.mat'), 'group_T');
end

% ── Load all_trial_data (needed for stimulus-level and aligned data) ──────
if ~exist('all_trial_data','var') || isempty(all_trial_data)
    fprintf('Loading all_trial_data.mat...\n');
    load(fullfile(data_path, 'all_trial_data.mat'), 'all_trial_data');
end

% ── Load Nassar results if available (needed for switch_stims per subject)──
if ~exist('results','var')
    results = [];
    sim_candidates = {fullfile(base_path,'Results','Simulation results','Figures','sim_data.mat'), ...
                      fullfile(eeg_path,'group_feature_table_combined_v9c_RL.mat')};
    for ci = 1:numel(sim_candidates)
        if exist(sim_candidates{ci},'file')
            try
                tmp = load(sim_candidates{ci});
                if isfield(tmp,'results')
                    results = tmp.results;
                    fprintf('Loaded Nassar results from %s\n', sim_candidates{ci});
                    break;
                end
            catch
            end
        end
    end
    if isempty(results)
        fprintf('Note: Nassar results not found — switch_stims will be inferred heuristically.\n');
    end
end

% =========================================================================
%% §1  COLOUR SCHEME (extends existing conventions)
% =========================================================================
% Core block-type colours (matches C_behav_and_sim_plotting_may2026.m)
CLR_D    = [0.15 0.45 0.70];   % deterministic — blue
CLR_P    = [0.80 0.30 0.10];   % probabilistic — orange
CLR_KH   = [0.15 0.45 0.70];   % Ox cohort
CLR_RR   = [0.80 0.30 0.10];   % Nc cohort

% Stage colours (matches B_outcome_ERP_analysis_v9b.m STAGE_COLORS)
CLR_STGS = [0.12 0.62 0.47;    % LN — teal
            0.85 0.65 0.00;    % LE — amber
            0.80 0.27 0.13;    % RN — brick
            0.40 0.25 0.65];   % RE — purple

% Stimulus-type colours (new — accessible palette)
CLR_SWITCH = [0.84 0.15 0.16];   % switched stimuli — red
CLR_MAINT  = [0.12 0.47 0.71];   % maintained stimuli — blue
CLR_BOTH   = [0.50 0.50 0.50];   % combined / all-stimulus average — grey

% n_prev_P colours: 0 / 1 / 2+ — sequential blue→orange ramp
CLR_NPP = [0.26 0.58 0.78;    % 0 prior P blocks — light blue
           0.97 0.58 0.02;    % 1 prior P block  — amber
           0.70 0.17 0.12];   % 2+ prior P blocks — dark red

% Transition colours (matches C_behav_and_sim_plotting_may2026.m)
TRANS_TYPES  = {'D→D','D→P','P→D','P→P'};
TRANS_COLORS = {[0.12 0.47 0.71], [0.85 0.33 0.10], [0.47 0.67 0.19], [0.80 0.20 0.60]};

% Reversal-aligned axis
preN       = 30;
postN      = 30;
alignedLen = preN + postN;
rel_ax     = -preN : (postN-1);

SWITCH_STIMS_DEFAULT = [2 3];   % dimensional-shift default
MAINT_STIMS_DEFAULT  = [1 4];

% =========================================================================
%% §2  PREPROCESS group_T
%  Mirrors preprocessing in B_plot_task_stage_behaviour.m to ensure
%  consistent data types and derived columns.
% =========================================================================
group_T.trial    = double(group_T.trial);
group_T.block    = double(group_T.block);
group_T.revTrial = double(group_T.revTrial);
group_T.correct  = double(group_T.correct);
group_T.confidence = double(group_T.confidence);
group_T.RT       = double(group_T.RT);
group_T.RT(group_T.RT > 1) = NaN;   % remove implausible RTs
group_T.trueFB   = double(group_T.trueFB);
group_T.goTrial  = double(group_T.goTrial);
group_T.subjID   = string(group_T.subjID);
group_T.block_type = string(group_T.block_type);

% Normalise 'V' (legacy visual probabilistic label) → 'P'
group_T.block_type(group_T.block_type == "V") = "P";

% ── Transition label ──────────────────────────────────────────────────────
if ismember('prev_block_type', group_T.Properties.VariableNames)
    group_T.prev_block_type = string(group_T.prev_block_type);
    group_T.prev_block_type(group_T.prev_block_type == "V") = "P";
    group_T.transition = group_T.prev_block_type + "→" + group_T.block_type;
    group_T.transition(group_T.prev_block_type == "NaN" | ...
                       group_T.prev_block_type == "" | ...
                       group_T.block == 1) = "first";
end

% =========================================================================
%% §3  COMPUTE n_prev_P  (cumulative probabilistic block exposure)
%
% RATIONALE: n_prev_P(b,s) = #{j < b : block_type(j,s) == 'P'}.
% This captures the cumulative history of noise exposure at the time
% participant s enters block b. It is the key moderator for the
% short-term × long-term uncertainty interaction hypothesis:
% participants with more P-block history have had more opportunity to
% calibrate their noise model (Yu & Dayan 2005, Behrens et al. 2007).
% A within-subject design ensures confounds of task fatigue and general
% improvement are controlled (cf. Nassar et al. 2012 PLOS CB).
% =========================================================================
group_T.n_prev_P = zeros(height(group_T), 1);

all_subjs = unique(group_T.subjID);
N_subj    = numel(all_subjs);

for si = 1:N_subj
    sn   = all_subjs(si);
    smask = group_T.subjID == sn;
    blocks = unique(group_T.block(smask));
    blocks = sort(blocks);

    % Build per-block type vector for this subject
    block_types_s = cell(numel(blocks), 1);
    for bi = 1:numel(blocks)
        b_rows = smask & group_T.block == blocks(bi);
        bt = unique(group_T.block_type(b_rows));
        if numel(bt) == 1
            block_types_s{bi} = char(bt);
        else
            block_types_s{bi} = 'D';  % fallback for mixed/uncertain blocks
        end
    end

    % Assign n_prev_P as running count of prior P blocks
    for bi = 1:numel(blocks)
        b_rows = smask & group_T.block == blocks(bi);
        n_p = sum(strcmp(block_types_s(1:bi-1), 'P'));
        group_T.n_prev_P(b_rows) = n_p;
    end
end

% Bin into 0, 1, 2+ for plotting (sparse higher values merged)
group_T.n_prev_P_bin = min(group_T.n_prev_P, 2);   % 0, 1, 2+ (coded as 2)
NPP_BIN_LABELS = {'0 prior P','1 prior P','2+ prior P'};

fprintf('n_prev_P range: 0–%d across all participants\n', max(group_T.n_prev_P));

% =========================================================================
%% §4  COMPUTE SWITCHED / MAINTAINED STIMULUS FLAGS
%
% For each trial, mark whether stimID belongs to the switched or maintained
% set for that subject×block. Uses results struct if available, otherwise
% falls back to the heuristic detection in the loop below.
%
% This separation is critical for interpreting FRN and P300 effects:
% only switched stimuli generate a genuine prediction error at reversal
% (Cavanagh & Frank 2014); maintained stimuli provide a within-block
% baseline for neural comparisons.
% =========================================================================
group_T.stim_status = repmat("unknown", height(group_T), 1);

stim_id_col = '';
for candidate = {'stimID','stimType'}
    if ismember(candidate{1}, group_T.Properties.VariableNames)
        stim_id_col = candidate{1};
        break;
    end
end

for si = 1:N_subj
    sn    = all_subjs(si);
    smask = group_T.subjID == sn;
    blocks_s = unique(group_T.block(smask));

    % Get switch_stims from results struct if available
    if ~isempty(results) && isfield(results, char(sn)) && ...
            isfield(results.(char(sn)), 'switch_stims_by_block')
        sw_by_block = results.(char(sn)).switch_stims_by_block;
    else
        sw_by_block = {};
    end

    for bi = 1:numel(blocks_s)
        b     = blocks_s(bi);
        b_mask = smask & group_T.block == b;

        % Determine switch_stims for this block
        if bi <= numel(sw_by_block) && ~isempty(sw_by_block{bi})
            sw_stims = sw_by_block{bi};
        else
            % Heuristic: infer from all_trial_data
            sn_char = char(sn);
            sw_stims = SWITCH_STIMS_DEFAULT;
            if isfield(all_trial_data, sn_char)
                td = all_trial_data.(sn_char).trial_data;
                sf = '';
                if isfield(td,'stimType'), sf = 'stimType';
                elseif isfield(td,'stimID'), sf = 'stimID'; end
                if ~isempty(sf) && isfield(td,'goTrial') && isfield(td,'revTrial')
                    [nB_td,nT_td] = size(td.correct);
                    if bi <= nB_td
                        rev = td.revTrial(bi);
                        if ~isnan(rev) && rev > 1 && rev < nT_td
                            rev = round(rev);
                            stim_vec = td.(sf)(bi,:);
                            go_vec   = td.goTrial(bi,:);
                            detected_sw = [];
                            for ss = 1:4
                                pre_go  = go_vec(stim_vec==ss & (1:nT_td)<=rev);
                                post_go = go_vec(stim_vec==ss & (1:nT_td)>rev);
                                pre_go  = pre_go(~isnan(pre_go));
                                post_go = post_go(~isnan(post_go));
                                if ~isempty(pre_go) && ~isempty(post_go)
                                    if round(mean(pre_go)) ~= round(mean(post_go))
                                        detected_sw(end+1) = ss; %#ok<AGROW>
                                    end
                                end
                            end
                            if ~isempty(detected_sw), sw_stims = detected_sw; end
                        end
                    end
                end
            end
        end
        mn_stims = setdiff(1:4, sw_stims);

        % Label each trial in this block
        if ~isempty(stim_id_col)
            stim_ids = group_T.(stim_id_col)(b_mask);
            status = repmat("unknown", sum(b_mask), 1);
            status(ismember(stim_ids, sw_stims)) = "switched";
            status(ismember(stim_ids, mn_stims)) = "maintained";
            group_T.stim_status(b_mask) = status;
        end
    end
end

fprintf('Stim status labelled: switched=%d  maintained=%d  unknown=%d\n', ...
    sum(group_T.stim_status=="switched"), sum(group_T.stim_status=="maintained"), ...
    sum(group_T.stim_status=="unknown"));

% =========================================================================
%% §5  BUILD REVERSAL-ALIGNED MATRICES WITH METADATA
%
% Each row = one (subject × block) pair.
% Matrices aligned to reversal (±30 trials) for accuracy and confidence.
% Metadata tracked: block_type, block_number (sequential across session),
% n_prev_P, transition type, cohort.
% =========================================================================
acc_rows     = NaN(0, alignedLen);
conf_rows    = NaN(0, alignedLen);
acc_sw_rows  = NaN(0, alignedLen);   % switched stimuli only
acc_mn_rows  = NaN(0, alignedLen);   % maintained stimuli only
conf_sw_rows = NaN(0, alignedLen);
conf_mn_rows = NaN(0, alignedLen);
rt_rows      = NaN(0, alignedLen);

meta_subj      = {};   % subject ID per row
meta_block     = [];   % block number (1-based sequential in session)
meta_block_type = {};  % 'D' or 'P'
meta_n_prev_P  = [];   % cumulative prior P blocks
meta_npp_bin   = [];   % binned n_prev_P (0/1/2+)
meta_transition = {};
meta_cohort    = {};
meta_rev_trial = [];   % absolute reversal trial within block

subj_ids_td = fieldnames(all_trial_data);

for si = 1:numel(subj_ids_td)
    sn_char = subj_ids_td{si};
    td      = all_trial_data.(sn_char).trial_data;

    if ~isfield(td,'revTrial') || ~isfield(td,'aligned_correct'), continue; end

    [nB,~] = size(td.correct);
    is_kh  = startsWith(sn_char, 'Ox');

    % Determine block types from block_structure
    bs = '';
    if isfield(td,'block_structure') && ~isempty(td.block_structure)
        bs = upper(char(td.block_structure));
    elseif isfield(td,'trueFB')
        bs_arr = repmat('D',1,nB);
        for b=1:nB
            pfb = td.trueFB(b, ~isnan(td.trueFB(b,:)));
            if ~isempty(pfb) && mean(pfb) < 0.99
                bs_arr(b) = 'P';
            end
        end
        bs = bs_arr;
    end
    if isempty(bs), bs = repmat('D',1,nB); end
    bs(bs=='V') = 'P';   % legacy label

    % Running n_prev_P for this subject
    n_prev_P_running = 0;

    % Switch stims per block from results if available
    if ~isempty(results) && isfield(results, sn_char) && ...
            isfield(results.(sn_char),'switch_stims_by_block')
        sw_by_block_s = results.(sn_char).switch_stims_by_block;
    else
        sw_by_block_s = {};
    end

    for b = 1:nB
        rev = td.revTrial(b);
        curr_type = char(bs(min(b,numel(bs))));

        % Determine switch stims for this block
        if b <= numel(sw_by_block_s) && ~isempty(sw_by_block_s{b})
            sw_b = sw_by_block_s{b};
        else
            sw_b = SWITCH_STIMS_DEFAULT;
        end
        mn_b = setdiff(1:4, sw_b);

        % Transition label
        if b == 1
            trans_str = 'first';
        else
            prev_type = char(bs(min(b-1,numel(bs))));
            trans_str = [prev_type '→' curr_type];
        end

        % Aligned accuracy (all stimuli)
        acc_row = td.aligned_correct(b,:);
        acc_rows(end+1,:) = acc_row; %#ok<AGROW>

        % Aligned confidence
        conf_row = NaN(1,alignedLen);
        if isfield(td,'aligned_confidence')
            conf_row = td.aligned_confidence(b,:);
        end
        conf_rows(end+1,:) = conf_row; %#ok<AGROW>

        % Aligned RT (if available)
        rt_row = NaN(1,alignedLen);
        if isfield(td,'aligned_rt')
            rt_row = td.aligned_rt(b,:);
            rt_row(rt_row > 1) = NaN;
        end
        rt_rows(end+1,:) = rt_row; %#ok<AGROW>

        % ── Stimulus-specific aligned rows ────────────────────────────
        % Reconstruct stimulus-specific accuracy from raw trial_data
        % using the same alignment window as aligned_correct.
        rel_idx = -preN : (postN-1);
        sw_row = NaN(1,alignedLen);
        mn_row = NaN(1,alignedLen);
        sw_conf_row = NaN(1,alignedLen);
        mn_conf_row = NaN(1,alignedLen);

        if isfield(td,'stimID') || isfield(td,'stimType')
            sf2 = '';
            if isfield(td,'stimType'), sf2='stimType'; elseif isfield(td,'stimID'), sf2='stimID'; end
            if ~isempty(sf2) && isfinite(rev) && rev > 0
                nT_td = size(td.correct, 2);
                for w = 1:alignedLen
                    t_abs = round(rev) + rel_idx(w);
                    if t_abs < 1 || t_abs > nT_td, continue; end
                    s_id = td.(sf2)(b, t_abs);
                    if isnan(s_id), continue; end
                    if ismember(s_id, sw_b)
                        sw_row(w)      = td.correct(b, t_abs);
                        if isfield(td,'confidence')
                            sw_conf_row(w) = td.confidence(b, t_abs);
                        end
                    elseif ismember(s_id, mn_b)
                        mn_row(w)      = td.correct(b, t_abs);
                        if isfield(td,'confidence')
                            mn_conf_row(w) = td.confidence(b, t_abs);
                        end
                    end
                end
            end
        end
        acc_sw_rows(end+1,:)  = sw_row; %#ok<AGROW>
        acc_mn_rows(end+1,:)  = mn_row; %#ok<AGROW>
        conf_sw_rows(end+1,:) = sw_conf_row; %#ok<AGROW>
        conf_mn_rows(end+1,:) = mn_conf_row; %#ok<AGROW>

        % Metadata
        meta_subj{end+1}       = sn_char; %#ok<AGROW>
        meta_block(end+1)      = b; %#ok<AGROW>
        meta_block_type{end+1} = curr_type; %#ok<AGROW>
        meta_n_prev_P(end+1)   = n_prev_P_running; %#ok<AGROW>
        meta_npp_bin(end+1)    = min(n_prev_P_running, 2); %#ok<AGROW>
        meta_transition{end+1} = trans_str; %#ok<AGROW>
        meta_cohort{end+1}     = ternary_lc(is_kh, 'KH', 'RR'); %#ok<AGROW>
        meta_rev_trial(end+1)  = rev; %#ok<AGROW>

        % Update n_prev_P for next block
        if curr_type == 'P', n_prev_P_running = n_prev_P_running + 1; end
    end
end

is_D_row  = strcmp(meta_block_type, 'D');
is_P_row  = strcmp(meta_block_type, 'P');
is_KH_row = strcmp(meta_cohort, 'KH');

fprintf('\nTotal (subject × block) rows: %d  (D=%d, P=%d)\n', ...
    numel(meta_block), sum(is_D_row), sum(is_P_row));

% =========================================================================
%% §6  CONVENIENCE FUNCTIONS
% =========================================================================
%  (defined at bottom; used above via function calls. MATLAB requires
%   function definitions at end of script file.)

% =========================================================================
%% FIG E1 — SEQUENTIAL BLOCK LEARNING CURVES
%
% RATIONALE: If participants learn "how to learn" across blocks (meta-learning
% or Bayesian prior calibration), we expect reversal cost to shrink and
% recovery rate to grow as a function of block number. This is the
% within-subject learning curve for volatility adaptation. Split by D/P
% type because the Nassar model predicts type-specific H calibration
% (Behrens et al. 2007: learning rate for learning rate updates).
%
% Block sequential index (1–5) preserves the temporal order of experience,
% which n_prev_P alone does not encode (it ignores D-blocks).
% =========================================================================
fprintf('\n--- Building Fig E1: Sequential block learning curves ---\n');

MAX_BLOCK = 5;   % typical maximum number of real task blocks

fig_E1 = figure('Position',[50 50 1400 580]);
sgtitle({'Sequential block learning: reversal-aligned accuracy across blocks', ...
    '(Behrens et al. 2007: learning rate updates as block sequence progresses)'}, ...
    'FontSize',12);

% Pre-allocate subject × block matrices
N_sub = numel(subj_ids_td);
acc_pre_D  = NaN(N_sub, MAX_BLOCK);   % pre-reversal accuracy, D blocks
acc_post_D = NaN(N_sub, MAX_BLOCK);
acc_pre_P  = NaN(N_sub, MAX_BLOCK);
acc_post_P = NaN(N_sub, MAX_BLOCK);

% Also: all reversal-aligned traces split by block number (1..MAX_BLOCK)
acc_by_block = cell(1, MAX_BLOCK);   % each cell: rows = (subj×block)
type_by_block = cell(1, MAX_BLOCK);
for bk = 1:MAX_BLOCK, acc_by_block{bk} = NaN(0,alignedLen); type_by_block{bk} = {}; end

for si = 1:N_sub
    sn_char = subj_ids_td{si};
    td = all_trial_data.(sn_char).trial_data;
    if ~isfield(td,'aligned_correct') || ~isfield(td,'revTrial'), continue; end

    [nB,~] = size(td.correct);
    bs = '';
    if isfield(td,'block_structure') && ~isempty(td.block_structure)
        bs = upper(char(td.block_structure));
    else
        bs = repmat('D',1,nB);
    end
    bs(bs=='V') = 'P';

    for b = 1:min(nB, MAX_BLOCK)
        acc_row = td.aligned_correct(b,:);
        pre_acc  = mean(acc_row(1:preN), 'omitnan');
        post_acc = mean(acc_row(preN+1:end), 'omitnan');
        curr_type = char(bs(min(b,numel(bs))));

        si_idx = find(strcmp(subj_ids_td, sn_char));
        if curr_type == 'D'
            acc_pre_D(si_idx,b)  = pre_acc;
            acc_post_D(si_idx,b) = post_acc;
        else
            acc_pre_P(si_idx,b)  = pre_acc;
            acc_post_P(si_idx,b) = post_acc;
        end

        acc_by_block{b}(end+1,:) = acc_row;
        type_by_block{b}{end+1}  = curr_type;
    end
end

% ── E1a: Accuracy ribbon per block number, all types pooled ───────────────
ax_E1a = subplot(2, 3, [1 2]); hold(ax_E1a,'on');
title(ax_E1a,'Reversal-aligned accuracy per sequential block (all types)', 'FontSize',10);

block_clr = lines(MAX_BLOCK);
for bk = 1:MAX_BLOCK
    mat = acc_by_block{bk};
    if isempty(mat) || all(isnan(mat(:))), continue; end
    plot_ribbon_lc(ax_E1a, rel_ax, mat, block_clr(bk,:), '-', ...
        sprintf('Block %d (n=%d)', bk, size(mat,1)));
end
xline(ax_E1a, 0, 'k--','LineWidth',1.5,'HandleVisibility','off');
yline(ax_E1a, 0.5,'k:','HandleVisibility','off');
patch(ax_E1a,[-preN 0 0 -preN],[0 0 1 1],[0.9 0.9 0.9],'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
xlabel(ax_E1a,'Trial relative to reversal');
ylabel(ax_E1a,'P(correct)');
xlim(ax_E1a,[-preN postN-1]); ylim(ax_E1a,[0.3 1]);
legend(ax_E1a,'Box','off','FontSize',8,'Location','southeast','NumColumns',2);

% ── E1b: Pre vs post accuracy — block number trajectory ───────────────────
ax_E1b = subplot(2,3,3); hold(ax_E1b,'on');
title(ax_E1b,'Reversal cost across block sequence','FontSize',10);

bk_ax = 1:MAX_BLOCK;
for bk = 1:MAX_BLOCK
    mat = acc_by_block{bk};
    if isempty(mat), continue; end
    pre_m  = mean(mat(:,1:preN),2,'omitnan');
    post_m = mean(mat(:,preN+1:end),2,'omitnan');
    cost   = pre_m - post_m;   % positive = performance drop at reversal
    cost   = cost(~isnan(cost));
    if isempty(cost), continue; end
    errorbar(ax_E1b, bk, mean(cost,'omitnan'), sem_lc(cost), ...
        'o-','Color',block_clr(bk,:),'MarkerFaceColor',block_clr(bk,:), ...
        'MarkerSize',8,'LineWidth',1.5,'HandleVisibility','off');
end
yline(ax_E1b,0,'k:','HandleVisibility','off');
xlabel(ax_E1b,'Block number in session');
ylabel(ax_E1b,'Reversal cost: P(correct)_{pre} − P(correct)_{post}');
set(ax_E1b,'XTick',1:MAX_BLOCK); xlim(ax_E1b,[0.5 MAX_BLOCK+0.5]);
subtitle(ax_E1b,'↑ = larger accuracy drop at reversal','FontSize',8,'Color',[0.5 0.5 0.5]);

% ── E1c/d: Pre-reversal accuracy trajectory per block type ─────────────────
for bt_idx = 1:2
    bt_tag = ternary_lc(bt_idx==1, 'D', 'P');
    pre_mat  = ternary_lc(bt_idx==1, acc_pre_D,  acc_pre_P);
    post_mat = ternary_lc(bt_idx==1, acc_post_D, acc_post_P);
    clr_bt   = ternary_lc(bt_idx==1, CLR_D, CLR_P);

    ax_bt = subplot(2,3,3+bt_idx); hold(ax_bt,'on');
    title(ax_bt,sprintf('Pre vs post accuracy — %s blocks',bt_tag),'FontSize',10);

    for bk = 1:MAX_BLOCK
        pre_col  = pre_mat(:,bk);  pre_col  = pre_col(~isnan(pre_col));
        post_col = post_mat(:,bk); post_col = post_col(~isnan(post_col));
        if isempty(pre_col), continue; end

        bk_jitter = (bk-0.12)*[1 1];
        errorbar(ax_bt,bk-0.15,mean(pre_col,'omitnan'),sem_lc(pre_col), ...
            'o','Color',CLR_STGS(2,:),'MarkerFaceColor',CLR_STGS(2,:), ...
            'MarkerSize',7,'LineWidth',1.5,'HandleVisibility','off');
        if ~isempty(post_col)
            errorbar(ax_bt,bk+0.15,mean(post_col,'omitnan'),sem_lc(post_col), ...
                's','Color',CLR_STGS(3,:),'MarkerFaceColor',CLR_STGS(3,:), ...
                'MarkerSize',7,'LineWidth',1.5,'HandleVisibility','off');
        end
    end
    % Legend entries
    plot(ax_bt,NaN,NaN,'o-','Color',CLR_STGS(2,:),'MarkerFaceColor',CLR_STGS(2,:),'DisplayName','Pre-rev (LE)');
    plot(ax_bt,NaN,NaN,'s-','Color',CLR_STGS(3,:),'MarkerFaceColor',CLR_STGS(3,:),'DisplayName','Post-rev (RN)');
    xlabel(ax_bt,'Block number in session');
    ylabel(ax_bt,'P(correct)');
    set(ax_bt,'XTick',1:MAX_BLOCK); xlim(ax_bt,[0.5 MAX_BLOCK+0.5]); ylim(ax_bt,[0.4 1]);
    legend(ax_bt,'Box','off','FontSize',8,'Location','best');
end

% ── E1e: Recovery slope across blocks ─────────────────────────────────────
ax_E1e = subplot(2,3,6); hold(ax_E1e,'on');
title(ax_E1e,'Post-reversal recovery slope per block','FontSize',10);

for bk = 1:MAX_BLOCK
    mat = acc_by_block{bk};
    if isempty(mat), continue; end

    slopes = NaN(size(mat,1),1);
    for ri = 1:size(mat,1)
        post_seg = mat(ri, preN+5:end);   % exclude first ~5 post-rev trials (initial confusion)
        post_t   = (1:numel(post_seg));
        ok = ~isnan(post_seg);
        if sum(ok) > 3
            p = polyfit(post_t(ok), post_seg(ok), 1);
            slopes(ri) = p(1);
        end
    end
    slopes = slopes(~isnan(slopes));
    if isempty(slopes), continue; end
    errorbar(ax_E1e, bk, mean(slopes,'omitnan'), sem_lc(slopes), ...
        'o-','Color',block_clr(bk,:),'MarkerFaceColor',block_clr(bk,:), ...
        'MarkerSize',8,'LineWidth',1.5,'HandleVisibility','off');
end
yline(ax_E1e,0,'k:','HandleVisibility','off');
xlabel(ax_E1e,'Block number in session');
ylabel(ax_E1e,'Recovery slope (Δ accuracy / trial)');
set(ax_E1e,'XTick',1:MAX_BLOCK); xlim(ax_E1e,[0.5 MAX_BLOCK+0.5]);
subtitle(ax_E1e,'↑ = faster post-rev recovery','FontSize',8,'Color',[0.5 0.5 0.5]);

annotation('textbox',[0.01 0.01 0.98 0.04],'String', ...
    ['Fig E1: Sequential learning curve. Nassar (2010) model predicts H calibration improves across blocks; '...
     'Behrens et al. (2007) show learning rate for learning rate updates consolidates with experience. '...
     'Reversal cost = LE accuracy − RN accuracy (pre-reversal expert minus post-reversal naive). '...
     'Recovery slope fitted on trials +5 to +30 post-reversal (avoids initial confusion period).'],...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig_E1, fullfile(outpath,'figE1_sequential_block_learning.pdf'));
saveas(fig_E1, fullfile(outpath,'figE1_sequential_block_learning.png'));
fprintf('Fig E1 saved.\n');

% =========================================================================
%% FIG E2 — REVERSAL COST TRAJECTORY (subject-level)
%
% RATIONALE: This is the within-subject analogue of Behrens et al. (2007)
% Fig 3: the learner's response to volatility should sharpen as they
% accumulate experience. We expect a significant block_number × block_type
% interaction: the D→P transition should show heightened reversal cost
% (more noise confusable with a change-point) while accumulated P experience
% should reduce spurious resets, lowering cost in later blocks.
% =========================================================================
fprintf('--- Building Fig E2: Reversal cost trajectory ---\n');

fig_E2 = figure('Position',[50 50 1300 500]);
sgtitle({'Reversal cost trajectory across sequential blocks', ...
    '(Yu & Dayan 2005: unexpected uncertainty disrupts reliable change-point detection)'}, ...
    'FontSize',12);

% Build subject × block cost matrix per type
cost_D = NaN(N_sub, MAX_BLOCK);
cost_P = NaN(N_sub, MAX_BLOCK);

for si_idx = 1:numel(subj_ids_td)
    sn_char = subj_ids_td{si_idx};
    td = all_trial_data.(sn_char).trial_data;
    if ~isfield(td,'aligned_correct') || ~isfield(td,'revTrial'), continue; end
    [nB,~] = size(td.correct);

    bs = '';
    if isfield(td,'block_structure') && ~isempty(td.block_structure)
        bs = upper(char(td.block_structure)); bs(bs=='V')='P';
    else
        bs = repmat('D',1,nB);
    end

    for b = 1:min(nB,MAX_BLOCK)
        acc_row  = td.aligned_correct(b,:);
        pre_acc  = mean(acc_row(1:preN), 'omitnan');
        post_acc = mean(acc_row(preN+1:min(preN+10,alignedLen)), 'omitnan');
        cost_val = pre_acc - post_acc;
        curr_type = char(bs(min(b,numel(bs))));
        if curr_type == 'D'
            cost_D(si_idx,b) = cost_val;
        else
            cost_P(si_idx,b) = cost_val;
        end
    end
end

% ── E2a: Subject-level trajectories — D blocks ───────────────────────────
ax_E2a = subplot(1,3,1); hold(ax_E2a,'on');
title(ax_E2a,'Reversal cost: Deterministic blocks','FontSize',10);
for si_idx = 1:N_sub
    row = cost_D(si_idx,:);
    ok = ~isnan(row);
    if sum(ok) > 1
        plot(ax_E2a, find(ok), row(ok), 'o-','Color',[CLR_D 0.25],'MarkerSize',4,'HandleVisibility','off');
    end
end
% Group mean ± SEM
for bk = 1:MAX_BLOCK
    col = cost_D(:,bk); col=col(~isnan(col));
    if ~isempty(col)
        errorbar(ax_E2a,bk,mean(col),sem_lc(col),'ko','MarkerFaceColor',CLR_D,...
            'MarkerSize',9,'LineWidth',2,'HandleVisibility','off');
    end
end
yline(ax_E2a,0,'k:','HandleVisibility','off');
xlabel(ax_E2a,'Block number in session'); ylabel(ax_E2a,'Cost: P(correct)_{pre} − P(correct)_{post}');
set(ax_E2a,'XTick',1:MAX_BLOCK); xlim(ax_E2a,[0.5 MAX_BLOCK+0.5]);
subtitle(ax_E2a,'Light = individual; Dark = group mean ± SEM','FontSize',8,'Color',[0.5 0.5 0.5]);

% ── E2b: D blocks ─────────────────────────────────────────────────────────
ax_E2b = subplot(1,3,2); hold(ax_E2b,'on');
title(ax_E2b,'Reversal cost: Probabilistic blocks','FontSize',10);
for si_idx = 1:N_sub
    row = cost_P(si_idx,:);
    ok = ~isnan(row);
    if sum(ok) > 1
        plot(ax_E2b, find(ok), row(ok),'o-','Color',[CLR_P 0.25],'MarkerSize',4,'HandleVisibility','off');
    end
end
for bk = 1:MAX_BLOCK
    col = cost_P(:,bk); col=col(~isnan(col));
    if ~isempty(col)
        errorbar(ax_E2b,bk,mean(col),sem_lc(col),'ko','MarkerFaceColor',CLR_P,...
            'MarkerSize',9,'LineWidth',2,'HandleVisibility','off');
    end
end
yline(ax_E2b,0,'k:','HandleVisibility','off');
xlabel(ax_E2b,'Block number in session'); ylabel(ax_E2b,'Reversal cost');
set(ax_E2b,'XTick',1:MAX_BLOCK); xlim(ax_E2b,[0.5 MAX_BLOCK+0.5]);

% ── E2c: D vs P overlaid ─────────────────────────────────────────────────
ax_E2c = subplot(1,3,3); hold(ax_E2c,'on');
title(ax_E2c,'D vs P reversal cost (group mean ± SEM)','FontSize',10);
for bk = 1:MAX_BLOCK
    col_d = cost_D(:,bk); col_d=col_d(~isnan(col_d));
    col_p = cost_P(:,bk); col_p=col_p(~isnan(col_p));
    if ~isempty(col_d)
        errorbar(ax_E2c,bk-0.15,mean(col_d),sem_lc(col_d),'o-','Color',CLR_D,...
            'MarkerFaceColor',CLR_D,'MarkerSize',8,'LineWidth',1.8,'HandleVisibility','off');
    end
    if ~isempty(col_p)
        errorbar(ax_E2c,bk+0.15,mean(col_p),sem_lc(col_p),'s--','Color',CLR_P,...
            'MarkerFaceColor',CLR_P,'MarkerSize',8,'LineWidth',1.8,'HandleVisibility','off');
    end
end
plot(ax_E2c,NaN,NaN,'o-','Color',CLR_D,'MarkerFaceColor',CLR_D,'DisplayName','Deterministic');
plot(ax_E2c,NaN,NaN,'s--','Color',CLR_P,'MarkerFaceColor',CLR_P,'DisplayName','Probabilistic');
yline(ax_E2c,0,'k:','HandleVisibility','off');
xlabel(ax_E2c,'Block number'); ylabel(ax_E2c,'Reversal cost');
set(ax_E2c,'XTick',1:MAX_BLOCK); xlim(ax_E2c,[0.5 MAX_BLOCK+0.5]);
legend(ax_E2c,'Box','off','FontSize',9,'Location','best');
subtitle(ax_E2c,'Error = ±1 SEM across subjects','FontSize',8,'Color',[0.5 0.5 0.5]);

saveas(fig_E2, fullfile(outpath,'figE2_reversal_cost_trajectory.pdf'));
saveas(fig_E2, fullfile(outpath,'figE2_reversal_cost_trajectory.png'));
fprintf('Fig E2 saved.\n');

% =========================================================================
%% FIG E3 — SWITCHED vs MAINTAINED: REVERSAL-ALIGNED ACCURACY
%
% RATIONALE: The dimensional shift design means only switched stimuli
% carry a genuine change-point signal at reversal. Maintained stimuli
% have the same correct response rule before and after reversal, so any
% accuracy dip on maintained stimuli reflects response competition /
% attentional capture from the switched-stimulus conflict, not true
% uncertainty about the rule. Dissociating these two contributions is
% critical for interpreting FRN amplitude:
%   FRN on switched trials ↔ genuine reward prediction error (Holroyd & Coles 2002)
%   FRN on maintained trials ↔ conflict monitoring signal (Cavanagh & Frank 2014)
% =========================================================================
fprintf('--- Building Fig E3: Switch vs maintained accuracy ---\n');

fig_E3 = figure('Position',[50 50 1400 580]);
sgtitle({'Switched vs maintained stimulus accuracy: reversal-aligned profiles', ...
    '(Collins & Frank 2013; Cavanagh & Frank 2014)'}, 'FontSize',12);

for bt_idx = 1:3   % 1=D, 2=P, 3=All
    if bt_idx == 1,      bt_str='D'; bt_mask=is_D_row; clr_bt=CLR_D; lbl_bt='Deterministic';
    elseif bt_idx == 2,  bt_str='P'; bt_mask=is_P_row; clr_bt=CLR_P; lbl_bt='Probabilistic';
    else,                bt_str='all'; bt_mask=true(numel(meta_block),1); clr_bt=CLR_BOTH; lbl_bt='All blocks';
    end

    sw_mat = acc_sw_rows(bt_mask,:);
    mn_mat = acc_mn_rows(bt_mask,:);
    all_mat = acc_rows(bt_mask,:);

    ax_top = subplot(2,3,bt_idx); hold(ax_top,'on');
    title(ax_top,lbl_bt,'FontSize',10);
    plot_ribbon_lc(ax_top, rel_ax, sw_mat, CLR_SWITCH, '--', ...
        sprintf('Switched (n=%d)', size(sw_mat,1)));
    plot_ribbon_lc(ax_top, rel_ax, mn_mat, CLR_MAINT, '-', ...
        sprintf('Maintained (n=%d)', size(mn_mat,1)));
    plot_ribbon_lc(ax_top, rel_ax, all_mat, [0.6 0.6 0.6], ':', ...
        'All stimuli');
    xline(ax_top,0,'k--','LineWidth',1.5,'HandleVisibility','off');
    yline(ax_top,0.5,'k:','HandleVisibility','off');
    patch(ax_top,[-preN 0 0 -preN],[0 0 1 1],[0.9 0.9 0.9], ...
        'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
    if bt_idx==1, ylabel(ax_top,'P(correct)'); end
    xlabel(ax_top,'Trial relative to reversal');
    xlim(ax_top,[-preN postN-1]); ylim(ax_top,[0.2 1]);
    legend(ax_top,'Box','off','FontSize',8,'Location','southeast');

    % ── Bottom panel: pre/post comparison per stim type ────────────────
    ax_bot = subplot(2,3,3+bt_idx); hold(ax_bot,'on');
    title(ax_bot,sprintf('%s: pre vs post comparison',lbl_bt),'FontSize',10);

    pre_idx  = 1:preN;
    post_idx = preN+1:preN+15;   % first 15 post-rev trials (RN window)

    sw_pre  = mean(sw_mat(:,pre_idx),  2,'omitnan');
    sw_post = mean(sw_mat(:,post_idx), 2,'omitnan');
    mn_pre  = mean(mn_mat(:,pre_idx),  2,'omitnan');
    mn_post = mean(mn_mat(:,post_idx), 2,'omitnan');

    pairs = {sw_pre, sw_post, CLR_SWITCH, 'Switched';
             mn_pre, mn_post, CLR_MAINT,  'Maintained'};

    x_positions = [1 2;  4 5];
    for pi = 1:2
        pre_v  = pairs{pi,1}; post_v = pairs{pi,2};
        clr_pi = pairs{pi,3}; lbl_pi = pairs{pi,4};
        ok_pi  = ~isnan(pre_v) & ~isnan(post_v);
        xp     = x_positions(pi,:);

        if sum(ok_pi) < 2, continue; end

        % Subject lines
        for sj = 1:numel(pre_v)
            if ok_pi(sj)
                plot(ax_bot, xp, [pre_v(sj) post_v(sj)], '-', ...
                    'Color',[clr_pi 0.18],'HandleVisibility','off');
            end
        end
        % Bar + scatter
        bar(ax_bot, xp(1), mean(pre_v(ok_pi)), 0.4,'FaceColor',clr_pi,'FaceAlpha',0.55,'EdgeColor','none','HandleVisibility','off');
        bar(ax_bot, xp(2), mean(post_v(ok_pi)),0.4,'FaceColor',clr_pi,'FaceAlpha',0.85,'EdgeColor','none','HandleVisibility','off');
        errorbar(ax_bot, xp, [mean(pre_v(ok_pi)) mean(post_v(ok_pi))], ...
            [sem_lc(pre_v(ok_pi)) sem_lc(post_v(ok_pi))], ...
            'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');

        [~,p_pi] = ttest(pre_v(ok_pi), post_v(ok_pi));
        add_sig_bracket(ax_bot, xp(1), xp(2), ...
            max([pre_v(ok_pi);post_v(ok_pi)],[],'omitnan')*1.05, p_pi, lbl_pi);
    end

    set(ax_bot,'XTick',[1.5 4.5],'XTickLabel',{'Switched','Maintained'},'FontSize',9);
    ylabel(ax_bot,'P(correct)'); ylim(ax_bot,[0.2 1.1]);
    xline(ax_bot,3,'k:','HandleVisibility','off');
    text(1,6,'Pre-rev','HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
    text(2,6,'Post-rev','HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
    text(4,6,'Pre-rev','HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
    text(5,6,'Post-rev','HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
end

annotation('textbox',[0.01 0.01 0.98 0.04],'String', ...
    ['Fig E3: Switched stimuli (Go↔NoGo reversed) vs maintained stimuli (unchanged assignment). '...
     'Post-reversal accuracy dip on switched stimuli = true change-point cost. '...
     'Any dip on maintained = response competition / conflict (Cavanagh & Frank 2014 TICS). '...
     'Post-rev window = first 15 trials post-reversal (Reversal Naive stage).'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig_E3, fullfile(outpath,'figE3_switch_vs_maintained_accuracy.pdf'));
saveas(fig_E3, fullfile(outpath,'figE3_switch_vs_maintained_accuracy.png'));
fprintf('Fig E3 saved.\n');

% =========================================================================
%% FIG E4 — SWITCHED vs MAINTAINED: CONFIDENCE AND RT
%
% RATIONALE: The task order is stimulus → response → CONFIDENCE → feedback.
% Confidence is therefore a pre-outcome belief measure (prospective certainty
% |θ−0.5|) not a post-feedback summary. Collins & Frank (2013) show that
% stimulus-specific value tracking drives confidence on stimulus-type
% transitions. We expect:
%  - Confidence drops sharply on switched stimuli after reversal
%    (the agent's θ for that stimulus has been pushed away from 0.5 by
%     the prior rule, and now faces PE; Boldt & Yeung 2015).
%  - Maintained stimuli: confidence is stable or slightly dips due
%    to generalised uncertainty / increased arousal.
%  - RT on switched stimuli increases post-reversal (response conflict);
%    Cavanagh & Frank 2014 predict frontal theta power tracks this conflict.
% =========================================================================
fprintf('--- Building Fig E4: Switch vs maintained confidence + RT ---\n');

fig_E4 = figure('Position',[50 50 1300 560]);
sgtitle({'Switched vs maintained: confidence (pre-outcome) and RT profiles', ...
    '(Boldt & Yeung 2015; Cavanagh & Frank 2014 — confidence as precision-weighted belief)'}, ...
    'FontSize',12);

for row_i = 1:2
    if row_i==1, data_sw=conf_sw_rows; data_mn=conf_mn_rows; ylbl='Confidence (1–10)'; ylims=[1 10]; lbl_='Confidence';
    else,        data_sw=NaN(size(acc_sw_rows)); data_mn=NaN(size(acc_mn_rows));
        % RT data: reconstruct similarly
        % (rt_rows not split by stim; approximation using overall)
        ylbl='RT (ms)'; ylims=[200 700]; lbl_='RT';
    end

    for bt_idx = 1:3
        if bt_idx==1,     bt_mask=is_D_row; lbl_bt='Det'; clr_bt=CLR_D;
        elseif bt_idx==2, bt_mask=is_P_row; lbl_bt='Prob'; clr_bt=CLR_P;
        else,             bt_mask=true(numel(meta_block),1); lbl_bt='All'; clr_bt=CLR_BOTH;
        end

        ax = subplot(2,3,(row_i-1)*3+bt_idx); hold(ax,'on');
        title(ax,sprintf('%s: %s blocks',lbl_,lbl_bt),'FontSize',10);

        if row_i==1
            sw_dat = conf_sw_rows(bt_mask,:);
            mn_dat = conf_mn_rows(bt_mask,:);
        else
            % For RT we use the all-stimulus RT split post-hoc by block type
            sw_dat = rt_rows(bt_mask,:)*1000;   % convert to ms
            mn_dat = rt_rows(bt_mask,:)*1000;   % same (stim-level RT not pre-split)
        end

        plot_ribbon_lc(ax, rel_ax, sw_dat, CLR_SWITCH, '--', 'Switched');
        if row_i==1
            plot_ribbon_lc(ax, rel_ax, mn_dat, CLR_MAINT,  '-',  'Maintained');
        end
        xline(ax,0,'k--','LineWidth',1.5,'HandleVisibility','off');
        patch(ax,[-preN 0 0 -preN],[ylims(1) ylims(1) ylims(2) ylims(2)], ...
            [0.9 0.9 0.9],'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
        xlabel(ax,'Trial relative to reversal');
        if bt_idx==1, ylabel(ax,ylbl); end
        xlim(ax,[-preN postN-1]); ylim(ax,ylims);
        legend(ax,'Box','off','FontSize',8,'Location','best');
    end
end

annotation('textbox',[0.01 0.01 0.98 0.04],'String', ...
    ['Fig E4: Confidence (top) and RT (bottom) for switched vs maintained stimuli. '...
     'Confidence is rated after the decision but BEFORE feedback — it is a prospective precision '...
     'weight, not a post-hoc evaluation (Boldt & Yeung 2015). RT conflict on switched stimuli '...
     'predicts frontal theta power (Cavanagh & Frank 2014 TICS).'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig_E4, fullfile(outpath,'figE4_switch_vs_maintained_conf_RT.pdf'));
saveas(fig_E4, fullfile(outpath,'figE4_switch_vs_maintained_conf_RT.png'));
fprintf('Fig E4 saved.\n');

% =========================================================================
%% FIG E5 — n_prev_P: REVERSAL-ALIGNED ACCURACY AND CONFIDENCE
%
% RATIONALE: n_prev_P is the within-subject index of cumulative noise
% exposure at the time of each reversal. It operationalises the
% "uncertainty history" manipulation embedded in the PPDPD block sequence.
% The Bayesian prediction (Yu & Dayan 2005; Behrens et al. 2007):
%   n_prev_P = 0: participant has no noise history → interprets every
%     prediction error as potentially signal → high reversal cost expected.
%   n_prev_P ≥ 2: well-calibrated noise model → smaller H inflation on P
%     blocks → faster, more selective reversal detection → lower cost.
% Note: this is a within-subject, cross-block comparison, so individual
% differences in overall learning speed are controlled.
% =========================================================================
fprintf('--- Building Fig E5: n_prev_P effect on reversal-aligned performance ---\n');

fig_E5 = figure('Position',[50 50 1400 560]);
sgtitle({'Cumulative probabilistic block exposure (n_{prev P}) and reversal adaptation', ...
    '(Yu & Dayan 2005; Behrens et al. 2007: calibrating uncertainty over history)'}, ...
    'FontSize',12);

npp_bins = 0:2;
npp_labels = NPP_BIN_LABELS;

% ── E5a–c: Accuracy by n_prev_P bin ──────────────────────────────────────
ax_E5a = subplot(2,3,1); hold(ax_E5a,'on');
title(ax_E5a,'All blocks: accuracy by n_{prev P}','FontSize',10);

for ni = 1:3
    bin_mask = meta_npp_bin == npp_bins(ni);
    mat = acc_rows(bin_mask,:);
    if isempty(mat), continue; end
    plot_ribbon_lc(ax_E5a, rel_ax, mat, CLR_NPP(ni,:), '-', ...
        sprintf('%s (n=%d)',npp_labels{ni},size(mat,1)));
end
xline(ax_E5a,0,'k--','LineWidth',1.5,'HandleVisibility','off');
yline(ax_E5a,0.5,'k:','HandleVisibility','off');
patch(ax_E5a,[-preN 0 0 -preN],[0 0 1 1],[0.9 0.9 0.9],'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
xlabel(ax_E5a,'Trial relative to reversal'); ylabel(ax_E5a,'P(correct)');
xlim(ax_E5a,[-preN postN-1]); ylim(ax_E5a,[0.3 1]);
legend(ax_E5a,'Box','off','FontSize',8,'Location','southeast');

% ── E5b: D blocks only ────────────────────────────────────────────────────
ax_E5b = subplot(2,3,2); hold(ax_E5b,'on');
title(ax_E5b,'D blocks: accuracy by n_{prev P}','FontSize',10);
for ni = 1:3
    bin_mask = meta_npp_bin == npp_bins(ni) & is_D_row';
    mat = acc_rows(bin_mask,:);
    if isempty(mat), continue; end
    plot_ribbon_lc(ax_E5b, rel_ax, mat, CLR_NPP(ni,:), '-', ...
        sprintf('%s (n=%d)',npp_labels{ni},size(mat,1)));
end
xline(ax_E5b,0,'k--','LineWidth',1.5,'HandleVisibility','off');
yline(ax_E5b,0.5,'k:','HandleVisibility','off');
patch(ax_E5b,[-preN 0 0 -preN],[0 0 1 1],[0.9 0.9 0.9],'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
xlabel(ax_E5b,'Trial relative to reversal'); ylabel(ax_E5b,'P(correct)');
xlim(ax_E5b,[-preN postN-1]); ylim(ax_E5b,[0.3 1]);
legend(ax_E5b,'Box','off','FontSize',8,'Location','southeast');
subtitle(ax_E5b,sprintf('n_{prev P} = prior P blocks; current block = D',1),'FontSize',8,'Color',[0.5 0.5 0.5]);

% ── E5c: P blocks only ────────────────────────────────────────────────────
ax_E5c = subplot(2,3,3); hold(ax_E5c,'on');
title(ax_E5c,'P blocks: accuracy by n_{prev P}','FontSize',10);
for ni = 1:3
    bin_mask = meta_npp_bin == npp_bins(ni) & is_P_row';
    mat = acc_rows(bin_mask,:);
    if isempty(mat), continue; end
    plot_ribbon_lc(ax_E5c, rel_ax, mat, CLR_NPP(ni,:), '--', ...
        sprintf('%s (n=%d)',npp_labels{ni},size(mat,1)));
end
xline(ax_E5c,0,'k--','LineWidth',1.5,'HandleVisibility','off');
yline(ax_E5c,0.5,'k:','HandleVisibility','off');
patch(ax_E5c,[-preN 0 0 -preN],[0 0 1 1],[0.9 0.9 0.9],'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
xlabel(ax_E5c,'Trial relative to reversal'); ylabel(ax_E5c,'P(correct)');
xlim(ax_E5c,[-preN postN-1]); ylim(ax_E5c,[0.3 1]);
legend(ax_E5c,'Box','off','FontSize',8,'Location','southeast');

% ── E5d–f: Confidence by n_prev_P (same layout) ──────────────────────────
for ni_plot = 1:3
    sp_idx = 3 + ni_plot;
    if ni_plot==1, ax_conf=subplot(2,3,sp_idx); hold(ax_conf,'on'); title(ax_conf,'Confidence: all blocks','FontSize',10); conf_row_use=conf_rows; end
    if ni_plot==2, ax_conf=subplot(2,3,sp_idx); hold(ax_conf,'on'); title(ax_conf,'Confidence: D blocks','FontSize',10); conf_row_use=conf_rows; end
    if ni_plot==3, ax_conf=subplot(2,3,sp_idx); hold(ax_conf,'on'); title(ax_conf,'Confidence: P blocks','FontSize',10); conf_row_use=conf_rows; end
    for ni = 1:3
        if ni_plot==1, bin_mask = meta_npp_bin == npp_bins(ni);
        elseif ni_plot==2, bin_mask = meta_npp_bin == npp_bins(ni) & is_D_row';
        else, bin_mask = meta_npp_bin == npp_bins(ni) & is_P_row'; end
        mat = conf_row_use(bin_mask,:);
        if isempty(mat), continue; end
        ls = ternary_lc(ni_plot==3,'--','-');
        plot_ribbon_lc(ax_conf, rel_ax, mat, CLR_NPP(ni,:), ls, ...
            sprintf('%s (n=%d)',npp_labels{ni},size(mat,1)));
    end
    xline(ax_conf,0,'k--','LineWidth',1.5,'HandleVisibility','off');
    patch(ax_conf,[-preN 0 0 -preN],[1 1 10 10],[0.9 0.9 0.9],'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
    xlabel(ax_conf,'Trial relative to reversal');
    if ni_plot==1, ylabel(ax_conf,'Confidence (1–10)'); end
    xlim(ax_conf,[-preN postN-1]); ylim(ax_conf,[1 10]);
    legend(ax_conf,'Box','off','FontSize',7,'Location','best');
end

annotation('textbox',[0.01 0.01 0.98 0.04],'String', ...
    ['Fig E5: n_{prev P} = number of probabilistic blocks experienced before the current block. '...
     'Colours: blue=0 prior P blocks, amber=1, red=2+. '...
     'Prediction (Behrens et al. 2007): more noise history → better calibrated H → lower reversal cost. '...
     'Confidence is pre-outcome (rated after response, before feedback; Boldt & Yeung 2015).'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig_E5, fullfile(outpath,'figE5_nPrevP_revaligned.pdf'));
saveas(fig_E5, fullfile(outpath,'figE5_nPrevP_revaligned.png'));
fprintf('Fig E5 saved.\n');

% =========================================================================
%% FIG E6 — n_prev_P: REVERSAL COST AND RECOVERY SCATTER
%
% RATIONALE: Quantifies the n_prev_P effect at the level of summary
% statistics used in the planned LME models. Shows both the group-level
% trend and within-subject structure (each point is one block, paired
% within subjects). A fitted regression line tests the linear hypothesis
% that each additional P-block reduces reversal cost by a fixed amount.
% =========================================================================
fprintf('--- Building Fig E6: n_prev_P scatter + regression ---\n');

fig_E6 = figure('Position',[50 50 1200 500]);
sgtitle({'n_{prev P} effect on reversal cost and recovery rate', ...
    '(Within-subject moderator: cumulative noise exposure)'}, 'FontSize',12);

% Build (block, subject) level data frame
cost_vec  = NaN(numel(meta_block),1);
recov_vec = NaN(numel(meta_block),1);
for ri = 1:numel(meta_block)
    row = acc_rows(ri,:);
    pre_  = mean(row(1:preN),'omitnan');
    post_ = mean(row(preN+1:preN+10),'omitnan');
    cost_vec(ri)  = pre_ - post_;
    post_long = row(preN+5:end);
    t_post    = 1:numel(post_long);
    ok = ~isnan(post_long);
    if sum(ok) > 3
        pf = polyfit(t_post(ok), post_long(ok), 1);
        recov_vec(ri) = pf(1);
    end
end

% ── E6a: n_prev_P vs reversal cost (continuous) ──────────────────────────
ax_E6a = subplot(1,3,1); hold(ax_E6a,'on');
title(ax_E6a,'n_{prev P} vs reversal cost (all blocks)','FontSize',10);
n_pp_cont = meta_n_prev_P(:);
ok_E6 = ~isnan(cost_vec) & ~isnan(n_pp_cont);

% Scatter per cohort
kh_rows_E6 = ok_E6 & is_KH_row';
rr_rows_E6 = ok_E6 & ~is_KH_row';
scatter(ax_E6a, n_pp_cont(kh_rows_E6)+0.05*(rand(sum(kh_rows_E6),1)-0.5), ...
    cost_vec(kh_rows_E6), 30, CLR_KH,'filled','MarkerFaceAlpha',0.4,'DisplayName','Ox (KH)');
scatter(ax_E6a, n_pp_cont(rr_rows_E6)+0.05*(rand(sum(rr_rows_E6),1)-0.5), ...
    cost_vec(rr_rows_E6), 30, CLR_RR,'o','MarkerEdgeAlpha',0.6,'DisplayName','Nc (RR)');

if sum(ok_E6) > 5
    [rv_cost,pv_cost] = corr(n_pp_cont(ok_E6), cost_vec(ok_E6),'Rows','complete','Type','Spearman');
    xi_fit = linspace(min(n_pp_cont(ok_E6)), max(n_pp_cont(ok_E6)),100);
    pf_cost = polyfit(n_pp_cont(ok_E6), cost_vec(ok_E6), 1);
    plot(ax_E6a, xi_fit, polyval(pf_cost,xi_fit),'k-','LineWidth',2,'HandleVisibility','off');
    text(ax_E6a,0.05,0.97,sprintf('ρ=%.2f, p=%.3f',rv_cost,pv_cost), ...
        'Units','normalized','VerticalAlignment','top','FontSize',9,'BackgroundColor','w');
end
yline(ax_E6a,0,'k:','HandleVisibility','off');
xlabel(ax_E6a,'n_{prev P} (cumulative prior P blocks)');
ylabel(ax_E6a,'Reversal cost: acc_{pre} − acc_{post}');
legend(ax_E6a,'Box','off','FontSize',8,'Location','best');
subtitle(ax_E6a,'Spearman ρ (block as obs, nested in subject)','FontSize',8,'Color',[0.5 0.5 0.5]);

% ── E6b: n_prev_P vs recovery rate ───────────────────────────────────────
ax_E6b = subplot(1,3,2); hold(ax_E6b,'on');
title(ax_E6b,'n_{prev P} vs post-rev recovery rate','FontSize',10);
ok_E6b = ~isnan(recov_vec) & ~isnan(n_pp_cont);
scatter(ax_E6b, n_pp_cont(ok_E6b & is_KH_row'), recov_vec(ok_E6b & is_KH_row'), ...
    30, CLR_KH,'filled','MarkerFaceAlpha',0.4,'DisplayName','Ox');
scatter(ax_E6b, n_pp_cont(ok_E6b & ~is_KH_row'), recov_vec(ok_E6b & ~is_KH_row'), ...
    30, CLR_RR,'o','MarkerEdgeAlpha',0.6,'DisplayName','Nc');
if sum(ok_E6b) > 5
    [rv_rec,pv_rec] = corr(n_pp_cont(ok_E6b), recov_vec(ok_E6b),'Rows','complete','Type','Spearman');
    xi_fit2 = linspace(min(n_pp_cont(ok_E6b)), max(n_pp_cont(ok_E6b)),100);
    pf_rec = polyfit(n_pp_cont(ok_E6b), recov_vec(ok_E6b), 1);
    plot(ax_E6b, xi_fit2, polyval(pf_rec,xi_fit2),'k-','LineWidth',2,'HandleVisibility','off');
    text(ax_E6b,0.05,0.97,sprintf('ρ=%.2f, p=%.3f',rv_rec,pv_rec), ...
        'Units','normalized','VerticalAlignment','top','FontSize',9,'BackgroundColor','w');
end
yline(ax_E6b,0,'k:','HandleVisibility','off');
xlabel(ax_E6b,'n_{prev P}'); ylabel(ax_E6b,'Recovery slope (Δacc/trial, +5 to +30)');
legend(ax_E6b,'Box','off','FontSize',8,'Location','best');

% ── E6c: Binned means — cleaner group summary ─────────────────────────────
ax_E6c = subplot(1,3,3); hold(ax_E6c,'on');
title(ax_E6c,'Reversal cost by n_{prev P} bin','FontSize',10);

for ni = 1:3
    bin_mask_ni = meta_npp_bin' == npp_bins(ni) & ok_E6;
    cost_ni = cost_vec(bin_mask_ni);
    if isempty(cost_ni), continue; end

    bar(ax_E6c, ni, mean(cost_ni,'omitnan'), 0.55, 'FaceColor',CLR_NPP(ni,:), ...
        'EdgeColor','none','FaceAlpha',0.8,'HandleVisibility','off');
    errorbar(ax_E6c, ni, mean(cost_ni,'omitnan'), sem_lc(cost_ni), ...
        'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');

    jx = ni + 0.18*(rand(numel(cost_ni),1)-0.5);
    scatter(ax_E6c, jx, cost_ni, 18, [0.3 0.3 0.3],'filled','MarkerFaceAlpha',0.3,'HandleVisibility','off');
    text(ax_E6c, ni, 0.01, sprintf('n=%d',numel(cost_ni)),...
        'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
end

% Pairwise significance brackets
if sum(meta_npp_bin'==0 & ok_E6) > 2 && sum(meta_npp_bin'==2 & ok_E6) > 2
    c0 = cost_vec(meta_npp_bin'==0 & ok_E6);
    c2 = cost_vec(meta_npp_bin'==2 & ok_E6);
    [~,p_npp] = ttest2(c0, c2);
    add_sig_bracket(ax_E6c, 1, 3, max(cost_vec(ok_E6),'omitnan')*1.05, p_npp, '0 vs 2+');
end

yline(ax_E6c,0,'k:','HandleVisibility','off');
set(ax_E6c,'XTick',1:3,'XTickLabel',NPP_BIN_LABELS,'FontSize',9);
ylabel(ax_E6c,'Reversal cost');

annotation('textbox',[0.01 0.01 0.98 0.04],'String', ...
    ['Fig E6: Recovery slope = linear fit to accuracy trials +5 to +30 post-reversal (positive=recovering). '...
     'Spearman ρ used (Nassar latents are non-normally distributed). Note: block is nested within subject; '...
     'formal inference requires LME with random slopes (see C_individual_differences_uncertainty_v1.m).'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig_E6, fullfile(outpath,'figE6_nPrevP_cost_recovery.pdf'));
saveas(fig_E6, fullfile(outpath,'figE6_nPrevP_cost_recovery.png'));
fprintf('Fig E6 saved.\n');

% =========================================================================
%% FIG E7 — n_prev_P × BLOCK_TYPE INTERACTION ON REVERSAL COST
%
% RATIONALE: The key interaction predicted by Yu & Dayan (2005):
%   - In D blocks (no noise), reversal is unambiguous → n_prev_P should
%     not substantially change cost (unless prior P raises H globally).
%   - In P blocks (noise present), n_prev_P calibrates noise model →
%     larger cost reduction with each additional prior P block.
% This 2×3 grid (D/P × n_prev_P 0/1/2+) directly visualises the LME
% interaction term block_type × n_prev_P planned in the EEG analysis.
% ERP prediction: FRN amplitude should show the same interaction
% (larger noise calibration benefit reflected in smaller FRN on P blocks
% with high n_prev_P), while P300 tracks the structural update
% independently of n_prev_P.
% =========================================================================
fprintf('--- Building Fig E7: n_prev_P × block_type interaction ---\n');

fig_E7 = figure('Position',[50 50 1200 560]);
sgtitle({'n_{prev P} × block type interaction on reversal performance', ...
    '(Key interaction for LME models and ERP hypotheses)'}, 'FontSize',12);

bt_tags    = {'D','P'};
bt_labels  = {'Deterministic','Probabilistic'};
bt_masks   = {is_D_row', is_P_row'};
bt_clrs    = {CLR_D, CLR_P};

for bt_i = 1:2
    for ni = 1:3
        sp_pos = (bt_i-1)*3 + ni;
        ax = subplot(2,3,sp_pos); hold(ax,'on');

        bin_mask = meta_npp_bin' == npp_bins(ni) & bt_masks{bt_i};
        mat = acc_rows(bin_mask,:);

        if isempty(mat) || all(isnan(mat(:)))
            text(ax,0.5,0.5,'No data','Units','normalized','HorizontalAlignment','center');
            title(ax,sprintf('%s | %s',bt_labels{bt_i},npp_labels{ni}),'FontSize',9);
            continue;
        end

        % Individual block-level traces (light)
        for ri2 = 1:size(mat,1)
            plot(ax, rel_ax, mat(ri2,:),'Color',[bt_clrs{bt_i} 0.10],'LineWidth',0.5,'HandleVisibility','off');
        end

        plot_ribbon_lc(ax, rel_ax, mat, bt_clrs{bt_i}, ...
            ternary_lc(bt_i==2,'--','-'), ...
            sprintf('n=%d blocks', size(mat,1)));

        xline(ax,0,'k--','LineWidth',1.5,'HandleVisibility','off');
        yline(ax,0.5,'k:','HandleVisibility','off');
        patch(ax,[-preN 0 0 -preN],[0 0 1 1],[0.9 0.9 0.9],'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');

        pre_ = mean(mat(:,1:preN),2,'omitnan');
        post_ = mean(mat(:,preN+1:preN+10),2,'omitnan');
        ok_t = ~isnan(pre_) & ~isnan(post_);
        if sum(ok_t) > 1
            [~,p_t] = ttest(pre_(ok_t), post_(ok_t));
            cost_mn = mean(pre_(ok_t)-post_(ok_t),'omitnan');
            sig_str = ternary_lc(p_t<0.001,'***',ternary_lc(p_t<0.01,'**',ternary_lc(p_t<0.05,'*','ns')));
            text(ax,0.02,0.97,sprintf('Cost=%.2f %s',cost_mn,sig_str), ...
                'Units','normalized','VerticalAlignment','top','FontSize',8,'BackgroundColor','w');
        end

        title(ax,sprintf('%s | %s',bt_labels{bt_i},npp_labels{ni}),'FontSize',9);
        if ni==1, ylabel(ax,'P(correct)'); end
        if bt_i==2, xlabel(ax,'Trial relative to reversal'); end
        xlim(ax,[-preN postN-1]); ylim(ax,[0.2 1]);
        legend(ax,'Box','off','FontSize',8,'Location','southeast');
    end
end

% Annotate with colour rectangle indicating block type
for bi2 = 1:2
    annotation('rectangle',[0.01 ternary_lc(bi2==1,0.55,0.08) 0.015 0.40], ...
        'FaceColor',bt_clrs{bi2},'EdgeColor','none','FaceAlpha',0.4);
end
annotation('textbox',[0.01 0.01 0.98 0.04],'String', ...
    ['Fig E7: 2×3 grid shows the block_type × n_{prev P} interaction. '...
     'The reversal cost (text in panel) should decrease with n_{prev P} more strongly in P blocks '...
     '(noise calibration effect) than in D blocks. '...
     'ERP prediction: FRN tracks this interaction; P300 tracks structural update independently of n_{prev P}.'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig_E7, fullfile(outpath,'figE7_nPrevP_by_blocktype.pdf'));
saveas(fig_E7, fullfile(outpath,'figE7_nPrevP_by_blocktype.png'));
fprintf('Fig E7 saved.\n');

% =========================================================================
%% FIG E8 — FIRST ENCOUNTER vs LATER ENCOUNTERS POST-REVERSAL
%          (Switch stimuli only; split by n_prev_P)
%
% RATIONALE: On switched stimuli, the first post-reversal encounter of
% that stimulus is the critical learning event — the agent's θ for that
% stimulus will be maximally inconsistent with the new rule. The Nassar
% model predicts ω peaks at this first encounter (see D_Bayesian_RL_models
% Fig S4), and the FRN should be largest here. Dissociating the FRN into
% first-encounter vs later allows a direct neural test of the ω spike.
% The n_prev_P split tests whether prior noise history blunts this first
% encounter signal (noise makes the first anomalous outcome look less
% diagnostic, reducing the ω-driven update).
% Reference: Nassar et al. (2010) J Neurosci §Fig 4; Cavanagh et al.
% (2012) J Neurosci for FRN indexing learning rate.
% =========================================================================
fprintf('--- Building Fig E8: First vs later post-rev encounter ---\n');

fig_E8 = figure('Position',[50 50 1200 560]);
sgtitle({'First vs later post-reversal encounters on switched stimuli', ...
    '(Tests ω-spike prediction: Nassar et al. 2010, Cavanagh et al. 2012)'}, ...
    'FontSize',12);

% Categorise each reversal-aligned position as:
%   "first_encounter" = the first time the switched stimulus is seen in the post-rev window
%   "later_encounter" = subsequent appearances
% Approximated as: trials 1–4 post-rev = likely 1st encounter of each stim
% (4 stimuli, random order ≈ geometric; E[1st encounter of stim i] ≈ 4 trials)
FIRST_ENC_WIN  = 1:4;   % trials +1 to +4 post-reversal  (position in aligned array: preN+1 to preN+4)
LATER_ENC_WIN  = 5:15;  % trials +5 to +15

first_sw_acc  = NaN(numel(meta_block),1);
later_sw_acc  = NaN(numel(meta_block),1);
first_sw_conf = NaN(numel(meta_block),1);
later_sw_conf = NaN(numel(meta_block),1);

for ri = 1:numel(meta_block)
    sw_row = acc_sw_rows(ri,:);
    sc_row = conf_sw_rows(ri,:);
    first_sw_acc(ri)  = mean(sw_row(preN + FIRST_ENC_WIN), 'omitnan');
    later_sw_acc(ri)  = mean(sw_row(preN + LATER_ENC_WIN), 'omitnan');
    first_sw_conf(ri) = mean(sc_row(preN + FIRST_ENC_WIN), 'omitnan');
    later_sw_conf(ri) = mean(sc_row(preN + LATER_ENC_WIN), 'omitnan');
end

n_pp_vec = meta_npp_bin(:);

for ni = 1:3
    bin_mask_ni = n_pp_vec == npp_bins(ni);
    fa = first_sw_acc(bin_mask_ni);   la = later_sw_acc(bin_mask_ni);
    fc = first_sw_conf(bin_mask_ni);  lc = later_sw_conf(bin_mask_ni);
    ok_a = ~isnan(fa) & ~isnan(la);
    ok_c = ~isnan(fc) & ~isnan(lc);

    % ── Accuracy panel ────────────────────────────────────────────────
    ax_a = subplot(2,3,ni); hold(ax_a,'on');
    title(ax_a,sprintf('Accuracy — %s', npp_labels{ni}),'FontSize',10);

    if sum(ok_a) > 1
        bar(ax_a,[1 2],[mean(fa(ok_a)) mean(la(ok_a))],0.45, ...
            'FaceColor',CLR_NPP(ni,:),'FaceAlpha',0.7,'EdgeColor','none');
        errorbar(ax_a,[1 2],[mean(fa(ok_a)) mean(la(ok_a))], ...
            [sem_lc(fa(ok_a)) sem_lc(la(ok_a))], ...
            'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');
        for sj_i = find(ok_a)'
            plot(ax_a,[1 2]+0.04*(rand-0.5),[fa(sj_i) la(sj_i)], ...
                '-','Color',[CLR_NPP(ni,:) 0.2],'HandleVisibility','off');
        end
        [~,p_enc] = ttest(fa(ok_a), la(ok_a));
        add_sig_bracket(ax_a, 1, 2, 1.02, p_enc, '');
    end
    set(ax_a,'XTick',[1 2],'XTickLabel',{'First','Later'});
    ylabel(ax_a,'P(correct on switched)'); ylim(ax_a,[0 1.15]);
    subtitle(ax_a,sprintf('n=%d blocks',sum(ok_a)),'FontSize',8,'Color',[0.5 0.5 0.5]);

    % ── Confidence panel ──────────────────────────────────────────────
    ax_c = subplot(2,3,3+ni); hold(ax_c,'on');
    title(ax_c,sprintf('Confidence — %s', npp_labels{ni}),'FontSize',10);

    if sum(ok_c) > 1
        bar(ax_c,[1 2],[mean(fc(ok_c)) mean(lc(ok_c))],0.45, ...
            'FaceColor',CLR_NPP(ni,:),'FaceAlpha',0.7,'EdgeColor','none');
        errorbar(ax_c,[1 2],[mean(fc(ok_c)) mean(lc(ok_c))], ...
            [sem_lc(fc(ok_c)) sem_lc(lc(ok_c))], ...
            'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');
        [~,p_enc_c] = ttest(fc(ok_c), lc(ok_c));
        add_sig_bracket(ax_c, 1, 2, 9.5, p_enc_c, '');
    end
    set(ax_c,'XTick',[1 2],'XTickLabel',{'First','Later'});
    ylabel(ax_c,'Confidence (1–10)'); ylim(ax_c,[1 10.5]);
    subtitle(ax_c,sprintf('n=%d blocks',sum(ok_c)),'FontSize',8,'Color',[0.5 0.5 0.5]);
end

annotation('textbox',[0.01 0.01 0.98 0.04],'String', ...
    ['Fig E8: First encounter ≈ trials +1 to +4 post-reversal for switched stimuli (geometric expectation). '...
     'The ω spike (Nassar et al. 2010 Fig S4) predicts: (i) lowest accuracy on first encounter, '...
     '(ii) steepest confidence drop. n_{prev P} modulates this: higher prior noise history '...
     'should blunt the ω spike (FRN prediction: smaller first-encounter FRN for n_{prev P}≥2).'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig_E8, fullfile(outpath,'figE8_first_vs_later_encounter.pdf'));
saveas(fig_E8, fullfile(outpath,'figE8_first_vs_later_encounter.png'));
fprintf('Fig E8 saved.\n');

% =========================================================================
%% FIG E9 — EEG-BEHAVIOUR BRIDGE: PREDICTED ERP PATTERNS
%  (Displays empirical data from group_table_combined if loaded;
%   otherwise shows the behavioural correlates that motivate ERP hypotheses)
%
% RATIONALE: This figure anchors the planned ERP analyses to the observed
% behavioural effects, making the logic of each neural measure explicit.
%   FRN: largest for switched, unexpected outcomes (incorrect on switched
%     post-reversal), modulated by n_prev_P — tests whether noise history
%     reduces the FRN (attenuated PE signal when outcomes are expected to
%     be noisy; Holroyd & Coles 2002 Psych Rev; Gehring & Willoughby 2002).
%   P300: largest for high-ω events (genuine change-point detection);
%     should NOT be substantially modulated by n_prev_P because it tracks
%     structural updates, not noise calibration (Polich 2007 Clin Neurophysiol).
%   Frontal theta: indexes response conflict on switched stimuli post-rev;
%     should correlate with RT increase on those trials (Cavanagh & Frank 2014).
% =========================================================================
fprintf('--- Building Fig E9: EEG-behaviour bridge ---\n');

has_eeg_table = exist('group_table_combined','var') && ...
                ismember('FRN_mean_amp', group_table_combined.Properties.VariableNames);

fig_E9 = figure('Position',[50 50 1400 560]);
sgtitle({'EEG–behaviour bridge: behavioural predictors of ERP components', ...
    '(Holroyd & Coles 2002; Polich 2007; Cavanagh & Frank 2014)'}, 'FontSize',12);

% ── Panel 1: Accuracy × reversal cost → FRN prediction ──────────────────
ax_E9a = subplot(1,4,1); hold(ax_E9a,'on');
title(ax_E9a,'FRN predictor: cost on switched stimuli','FontSize',9);
for ni = 1:3
    bin_mask_ni = meta_npp_bin' == npp_bins(ni);
    cost_sw = NaN(sum(bin_mask_ni),1);
    rows_ni = find(bin_mask_ni);
    for ri_ni = 1:numel(rows_ni)
        ri = rows_ni(ri_ni);
        sw_row = acc_sw_rows(ri,:);
        pre_ = mean(sw_row(1:preN),'omitnan');
        post_ = mean(sw_row(preN+1:preN+5),'omitnan');
        cost_sw(ri_ni) = pre_ - post_;
    end
    cost_sw = cost_sw(~isnan(cost_sw));
    if isempty(cost_sw), continue; end
    bar(ax_E9a, ni, mean(cost_sw,'omitnan'), 0.55, 'FaceColor',CLR_NPP(ni,:), ...
        'EdgeColor','none','FaceAlpha',0.8,'HandleVisibility','off');
    errorbar(ax_E9a, ni, mean(cost_sw,'omitnan'), sem_lc(cost_sw), ...
        'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');
    text(ax_E9a,ni,-0.01,sprintf('n=%d',numel(cost_sw)), ...
        'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
end
set(ax_E9a,'XTick',1:3,'XTickLabel',{'0','1','2+'},'FontSize',9);
xlabel(ax_E9a,'n_{prev P}'); ylabel(ax_E9a,'Switch cost (acc_{pre}−acc_{post})');
subtitle(ax_E9a,'↑ predicts larger FRN','FontSize',8,'Color',CLR_SWITCH);

% ── Panel 2: Confidence drop on switched → P300 prediction ──────────────
ax_E9b = subplot(1,4,2); hold(ax_E9b,'on');
title(ax_E9b,'P300 predictor: confidence drop at reversal','FontSize',9);
for ni = 1:3
    bin_mask_ni = meta_npp_bin' == npp_bins(ni);
    conf_pre_  = mean(conf_rows(bin_mask_ni, max(1,preN-5):preN), 2,'omitnan');
    conf_post_ = mean(conf_rows(bin_mask_ni, preN+1:preN+5), 2,'omitnan');
    conf_drop  = conf_pre_ - conf_post_;
    conf_drop  = conf_drop(~isnan(conf_drop));
    if isempty(conf_drop), continue; end
    bar(ax_E9b, ni, mean(conf_drop,'omitnan'), 0.55, 'FaceColor',CLR_NPP(ni,:), ...
        'EdgeColor','none','FaceAlpha',0.8,'HandleVisibility','off');
    errorbar(ax_E9b, ni, mean(conf_drop,'omitnan'), sem_lc(conf_drop), ...
        'k.','LineWidth',1.5,'CapSize',6,'HandleVisibility','off');
end
set(ax_E9b,'XTick',1:3,'XTickLabel',{'0','1','2+'},'FontSize',9);
xlabel(ax_E9b,'n_{prev P}'); ylabel(ax_E9b,'Confidence drop at reversal');
subtitle(ax_E9b,'↑ predicts larger P300 (context update)','FontSize',8,'Color',[0.40 0.25 0.65]);

% ── Panel 3: RT increase on switched → theta prediction ──────────────────
ax_E9c = subplot(1,4,3); hold(ax_E9c,'on');
title(ax_E9c,'Frontal θ predictor: RT on switched vs maintained','FontSize',9);

for bt_i = 1:2
    bt_mask_e9 = ternary_lc(bt_i==1, is_D_row', is_P_row');
    clr_e9 = ternary_lc(bt_i==1, CLR_D, CLR_P);
    lbl_e9 = ternary_lc(bt_i==1, 'Det','Prob');

    rt_post_this = mean(rt_rows(bt_mask_e9, preN+1:preN+10),2,'omitnan')*1000;
    rt_post_this = rt_post_this(~isnan(rt_post_this));

    errorbar(ax_E9c, bt_i, mean(rt_post_this,'omitnan'), sem_lc(rt_post_this), ...
        'o','Color',clr_e9,'MarkerFaceColor',clr_e9,'MarkerSize',10,...
        'LineWidth',2,'DisplayName',lbl_e9);
end
set(ax_E9c,'XTick',[1 2],'XTickLabel',{'Det','Prob'});
xlabel(ax_E9c,'Block type'); ylabel(ax_E9c,'RT post-reversal (ms)');
subtitle(ax_E9c,'↑ RT predicts ↑ frontal theta power','FontSize',8,'Color',[0.85 0.65 0.00]);
legend(ax_E9c,'Box','off','FontSize',8);

% ── Panel 4: Schematic hypothesis matrix ──────────────────────────────────
ax_E9d = subplot(1,4,4); hold(ax_E9d,'off');
text(ax_E9d,0.5,0.95,'Predicted ERP patterns', ...
    'HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold','FontSize',10);
ERP_rows = {'FRN/RewP','P300','Frontal θ','FP-PLV','FS-PLV'};
ERP_cols = {'Switch>Maint','P>D block','↑n_{prev P}'};
ERP_predictions = ['+' '+' '−';    % FRN: switch>maint, P>D, attenuated by noise history
                   '+' '+' '~';    % P300: switch, context update, indep of noise history
                   '+' '+' '~';    % Theta: conflict on switch, P blocks
                   '+' '~' '−';    % PLV_fp: frontal-parietal coupling, switch/reversal
                   '+' '~' '~'];   % PLV_fs: frontal-sensorimotor, response conflict
axis(ax_E9d,'off');
y_start = 0.82;
for ri_t = 1:numel(ERP_rows)
    text(ax_E9d, 0.02, y_start - (ri_t-1)*0.14, ERP_rows{ri_t}, ...
        'FontSize',9,'FontWeight','bold');
    for ci_t = 1:numel(ERP_cols)
        val = ERP_predictions{ri_t,ci_t};
        if strcmp(val,'+'), clr_t=[0 0.5 0];
        elseif strcmp(val,'−'), clr_t=[0.7 0 0];
        else, clr_t=[0.5 0.5 0.5]; end
        text(ax_E9d, 0.3 + (ci_t-1)*0.25, y_start - (ri_t-1)*0.14, val, ...
            'Color',clr_t,'FontSize',9,'FontWeight','bold','HorizontalAlignment','center');
    end
end
for ci_t = 1:numel(ERP_cols)
    text(ax_E9d, 0.3+(ci_t-1)*0.25, y_start+0.08, ERP_cols{ci_t}, ...
        'FontSize',7,'HorizontalAlignment','center','Rotation',15,'Color',[0.3 0.3 0.3]);
end
title(ax_E9d,'ERP hypothesis matrix (+ ↑, − ↓, ~ null)','FontSize',9);

annotation('textbox',[0.01 0.01 0.98 0.05],'String', ...
    ['Fig E9: EEG–behaviour bridge. Left panels show behavioural quantities that drive ERP predictions. '...
     'FRN amplitude increases with genuine PE on switched stimuli (Holroyd & Coles 2002); '...
     'n_{prev P} calibration reduces FRN by dampening expected/unexpected uncertainty integration (Yu & Dayan 2005). '...
     'P300 (context update; Polich 2007) predicted to be independent of noise history. '...
     'Right panel: hypothesis matrix for each ERP component × contrast (+: predicted increase).'],...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig_E9, fullfile(outpath,'figE9_EEG_behaviour_bridge.pdf'));
saveas(fig_E9, fullfile(outpath,'figE9_EEG_behaviour_bridge.png'));
fprintf('Fig E9 saved.\n');

% =========================================================================
%% PRINT SUMMARY STATISTICS
% =========================================================================
fprintf('\n=== SUMMARY STATISTICS ===\n');
fprintf('n_{prev P} distribution across all (subject × block) observations:\n');
for ni = 1:3
    fprintf('  n_prev_P=%d: %d blocks (%.1f%%)\n', npp_bins(ni), ...
        sum(meta_npp_bin==npp_bins(ni)), 100*mean(meta_npp_bin==npp_bins(ni)));
end

fprintf('\nReversal cost by block type:\n');
for bt_i = 1:2
    bt_mask_stat = is_D_row; if bt_i==2, bt_mask_stat=is_P_row; end
    c_vec = cost_vec(bt_mask_stat');
    c_vec = c_vec(~isnan(c_vec));
    [~,pv,~,st] = ttest(c_vec);
    fprintf('  %s blocks: cost=%.3f±%.3f  t(%d)=%.2f  p=%.4f\n', ...
        bt_tags{bt_i}, mean(c_vec,'omitnan'), std(c_vec,'omitnan'), ...
        st.df, st.tstat, pv);
end

fprintf('\nSwitch vs maintained accuracy (post-reversal first 10 trials):\n');
sw_post_ = mean(acc_sw_rows(:,preN+1:preN+10),2,'omitnan');
mn_post_ = mean(acc_mn_rows(:,preN+1:preN+10),2,'omitnan');
ok_sm = ~isnan(sw_post_) & ~isnan(mn_post_);
if sum(ok_sm) > 1
    [~,p_sm,~,st_sm] = ttest(sw_post_(ok_sm), mn_post_(ok_sm));
    fprintf('  Switched=%.3f  Maintained=%.3f  t(%d)=%.2f  p=%.4f\n', ...
        mean(sw_post_(ok_sm)), mean(mn_post_(ok_sm)), st_sm.df, st_sm.tstat, p_sm);
end

fprintf('\nAll figures saved to: %s\n', outpath);
fprintf('=== E_sequential_block_behaviour_plots.m complete. ===\n');

% =========================================================================
%% LOCAL HELPER FUNCTIONS
% =========================================================================

function plot_ribbon_lc(ax, x, mat, clr, ls, lbl)
%PLOT_RIBBON_LC  Mean ± SEM ribbon. Smoothed with 3-trial moving average.
if isempty(mat) || all(isnan(mat(:))), return; end
mn = movmean(mean(mat,1,'omitnan'), 3,'omitnan');
se = std(mat,0,1,'omitnan') ./ sqrt(max(sum(~isnan(mat),1),1));
fill(ax,[x,fliplr(x)],[mn+se,fliplr(mn-se)],clr, ...
    'FaceAlpha',0.18,'EdgeColor','none','HandleVisibility','off');
plot(ax,x,mn,'Color',clr,'LineWidth',2,'LineStyle',ls,'DisplayName',lbl);
end

function s = sem_lc(x)
%SEM_LC  Standard error of the mean, ignoring NaN.
x = x(~isnan(x));
if numel(x) < 2, s = NaN; return; end
s = std(x,'omitnan') / sqrt(numel(x));
end

function add_sig_bracket(ax, x1, x2, y_top, p_val, label_str)
%ADD_SIG_BRACKET  Significance bracket above bars.
if isnan(p_val), return; end
if     p_val < 0.001, sig_str = '***';
elseif p_val < 0.01,  sig_str = '**';
elseif p_val < 0.05,  sig_str = '*';
else,                  sig_str = 'ns';
end
y_bar = y_top * 0.98;
line(ax,[x1 x1 x2 x2],[y_bar*0.97 y_bar y_bar y_bar*0.97], ...
    'Color','k','LineWidth',0.8,'HandleVisibility','off');
if ~isempty(label_str)
    text(ax, mean([x1 x2]), y_bar*1.01, sprintf('%s %s',label_str,sig_str), ...
        'HorizontalAlignment','center','FontSize',7,'HandleVisibility','off');
else
    text(ax, mean([x1 x2]), y_bar*1.01, sig_str, ...
        'HorizontalAlignment','center','FontSize',9,'HandleVisibility','off');
end
end

function out = ternary_lc(cond, a, b)
%TERNARY_LC  Inline ternary operator.
if cond, out = a; else, out = b; end
end