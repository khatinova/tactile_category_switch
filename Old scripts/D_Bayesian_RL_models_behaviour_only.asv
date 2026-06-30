% ==========================================================================
% BAYESIAN DELTA-RULE SIMULATION & FITTING — v2
% Category Switch Braille Go/NoGo Task
%
% MODELS
% ──────
%   Nassar (2010)  — adaptive learning rate via change-point probability ω.
%                    Fitted parameters: H (hazard rate), β (softmax temperature).
%   Rescorla-Wagner — subject-wise RW with separate α_det / α_prob, plus a
%                    confidence-weighted pre/post-reversal variant.
%                    Implemented via fit_RW_subjectwise (see bottom of file).
%
% STRUCTURE
% ─────────
%   §1  Setup & paths
%   §2  Simulation settings + stimulus structure
%   §3  Nassar simulation loop
%   §4  Simulation figures (beliefs, accuracy, learning rate)
%   §5  Fit Nassar to real data (all blocks combined)
%   §6  Fit Nassar separately for D and P blocks
%   §7  Visualise fitted parameters
%   §8  Parameter recovery (fitted → simulate → refit)
%   §9  Model comparison: Nassar BIC vs RW BIC
%   §10 Posterior predictive checks (Nassar only)
%   §11 Add Nassar latents + switch_stims to group_T
%
% KEY DESIGN NOTES
% ─────────────────
%   DIMENSIONAL SHIFT: 2/4 stimuli switch Go/NoGo at reversal; 2 maintained.
%   Only switched stimuli generate a change-point signal post-reversal.
%
%   PER-STIMULUS CP DETECTION: ω is computed independently per stimulus.
%   No cross-stimulus propagation.
%
%   PARTICIPANTS UNAWARE OF PROBABILISTIC FEEDBACK: y_t is derived from the
%   feedback shown on screen (perceivedCorrect), not ground truth.
%   H_prob − H_det therefore indexes noise-induced apparent volatility.
%
%   PROSPECTIVE CERTAINTY |θ−0.5| is the correct confidence predictor
%   (computed before feedback, matching the trial order: stimulus → response
%   → confidence rating → feedback).
%
% Behaviour-only version: does not load EEG feature tables.
% This script fits/simulates the behavioural RL model, saves Nassar results,
% and optionally writes a behaviour-only group_T table enriched with Nassar
% latents via (subjID, block, trial) key.
% ==========================================================================

clear; close all; rng(42);

% =========================================================================
%% §1  SETUP & PATHS
% =========================================================================
remote = 0;
switch remote
    case 1
        base_path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch/Salient mod switch KH';
    case 0
        base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH';
end

outpath   = fullfile(base_path, 'Results', 'Simulation results', 'Figures');
data_path = fullfile(base_path, 'Data');
eeg_path  = fullfile(base_path, 'Results', 'EEG analysis', 'Epoched_data'); %#ok<NASGU>

if ~exist(outpath, 'dir'), mkdir(outpath); end

% -------------------------------------------------------------------------
% Behaviour-only loading
% -------------------------------------------------------------------------
% Do NOT load group_feature_table_combined_v9c.mat here.  That file is EEG-
% feature-specific and is handled by F_RL_EEG_analysis_from_features.m.
%
% This script only needs:
%   all_trial_data : per-subject behavioural trial struct
%   group_T        : long-format behavioural table, when available
% -------------------------------------------------------------------------
load(fullfile(data_path, 'all_trial_data.mat'), 'all_trial_data');

group_T = [];
beh_candidates = { ...
    fullfile(data_path, 'behav_table_June2026.mat'), ...
    fullfile(data_path, 'behav_table.mat'), ...
    fullfile(data_path, 'group_T.mat')};

for ci = 1:numel(beh_candidates)
    if exist(beh_candidates{ci}, 'file')
        tmp_beh = load(beh_candidates{ci});
        if isfield(tmp_beh, 'group_T')
            group_T = tmp_beh.group_T;
            fprintf('Loaded behavioural group_T from %s\n', beh_candidates{ci});
            break;
        elseif isfield(tmp_beh, 'group_table_combined')
            group_T = tmp_beh.group_table_combined;
            fprintf('Loaded behavioural table group_table_combined as group_T from %s\n', beh_candidates{ci});
            break;
        elseif isfield(tmp_beh, 'behav_table')
            group_T = tmp_beh.behav_table;
            fprintf('Loaded behav_table as group_T from %s\n', beh_candidates{ci});
            break;
        end
    end
end

if isempty(group_T)
    warning(['No behavioural group_T table found. Model fitting will still run from all_trial_data, ' ...
             'but §11 will skip writing an enriched behaviour table.']);
else
    if ismember('subj_id', group_T.Properties.VariableNames) && ~ismember('subjID', group_T.Properties.VariableNames)
        group_T.subjID = group_T.subj_id;
    end
    if ismember('blocknum', group_T.Properties.VariableNames) && ~ismember('block', group_T.Properties.VariableNames)
        group_T.block = group_T.blocknum;
    end
    if ismember('trialnum', group_T.Properties.VariableNames) && ~ismember('trial', group_T.Properties.VariableNames)
        group_T.trial = group_T.trialnum;
    end
end

% ── Print all output locations upfront ───────────────────────────────────
fprintf('\n========================================================\n');
fprintf('OUTPUT LOCATIONS — BEHAVIOUR-ONLY RL SCRIPT\n');
fprintf('========================================================\n');
fprintf('Figures (PDFs):\n  %s\n', outpath);
fprintf('\nSaved data:\n');
fprintf('  %s\n', fullfile(outpath,'sim_data.mat'));
fprintf('    → sim_data, model_vars, params\n');
fprintf('  %s\n', fullfile(outpath,'nassar_results.mat'));
fprintf('    → results  (per-subject Nassar fits: H, β, latents, BIC)\n');
fprintf('    → loaded by F_RL_EEG_analysis_from_features.m\n');
fprintf('  %s\n', fullfile(data_path,'behav_table_June2026_RL.mat'));
fprintf('    → behaviour-only group_T + Nassar latents, if group_T is available\n');
fprintf('========================================================\n\n');

% =========================================================================
%% §2  SIMULATION SETTINGS & STIMULUS STRUCTURE
% =========================================================================
N_SUBJECTS = 50;
N_BLOCKS   = 5;
N_TRIALS   = 100;    % N_STIM × N_REPS
N_STIM     = 4;
N_REPS     = 25;
REV_MIN    = 30;
REV_MAX    = 70;

% =========================================================================
%% SUBJECT-SPECIFIC BLOCK ORDERS
% =========================================================================

BLOCK_ORDER_SET = {
    'DDPPD'
    'DPDPD'
    'PPDDP'
    'PDPDD'
    'DPPDD'
    'PDDPP'
};

% Feedback fidelity by block type
P_TRUE_FB_D = 1.0;
P_TRUE_FB_P = 0.8;
BLOCK_IS_DET = P_TRUE_FB >= 0.99;

SOFTMAX_BETA   = 6;
CONF_SCALE_MAX = 10;
CONF_OFFSET    = 1;
CONF_GAIN      = (CONF_SCALE_MAX - CONF_OFFSET) / 0.5;   % = 18
ALPHA_INIT     = 0.5;   % first-encounter learning rate boost

% ── Stimulus feature structure — 2×2 dimensional shift ───────────────────
%   Pre-rev rule:  Go iff Feature1 == 1  → stims 1,2 are Go
%   Post-rev rule: Go iff Feature2 == 1  → stims 1,3 are Go
%
%   Stim | F1 | F2 | Pre-rev | Post-rev | Status
%   ─────┼────┼────┼─────────┼──────────┼─────────────
%     1  |  1 |  1 |   Go    |   Go     | MAINTAINED
%     2  |  1 |  2 |   Go    |  NoGo    | SWITCHED
%     3  |  2 |  1 |  NoGo   |   Go     | SWITCHED
%     4  |  2 |  2 |  NoGo   |  NoGo    | MAINTAINED

STIM_FEATURES = [1 1; 1 2; 2 1; 2 2];
PRE_REV_RULE  = @(s) STIM_FEATURES(s,1) == 1;
POST_REV_RULE = @(s) STIM_FEATURES(s,2) == 1;

% ── Read SWITCH_STIMS from behavioural group_T when available ─────────────
% The behavioural extraction pipeline should preserve switch_stims.  If not,
% the dimensional-shift default [2 3] is used, and fitting itself remains
% stimulus-specific because it reads the actual stimID/stimType sequence.
SWITCH_STIMS = read_switch_stims_from_table_lc(group_T);
if isempty(SWITCH_STIMS)
    warning('switch_stims not found or empty in behavioural group_T. Defaulting to [2 3].');
    SWITCH_STIMS = [2 3];
end
MAINT_STIMS = setdiff(1:size(STIM_FEATURES,1), SWITCH_STIMS);
fprintf('SWITCH_STIMS used for summaries/plots: [%s]\n', num2str(SWITCH_STIMS));
fprintf('MAINT_STIMS: [%s]\n', num2str(MAINT_STIMS));

% Nassar simulation parameters
params.H    = 0.05;
params.beta = SOFTMAX_BETA;

% =========================================================================
%% §3  NASSAR SIMULATION LOOP
% =========================================================================
fprintf('Simulating %d subjects (Nassar model)...\n', N_SUBJECTS);

for subj = 1:N_SUBJECTS
    sl = sprintf('sub%02d', subj);

    subj_structure = BLOCK_ORDER_SET{randi(numel(BLOCK_ORDER_SET))};
    N_BLOCKS = numel(subj_structure);

    reversal_trials = randi([REV_MIN, REV_MAX], 1, N_BLOCKS);

    td = init_trial_data(N_BLOCKS, N_TRIALS);
    mv = init_model_vars(N_BLOCKS, N_TRIALS, N_STIM);

    [N_BLOCKS, N_TRIALS] = size(td.correct);
    for b = 1:N_BLOCKS
        curr_block_type = subj_structure(b);
        td.block_type_label{b} = char(curr_block_type);
        td.block_is_det(b) = curr_block_type == 'D';


        if curr_block_type == 'D'
            trueFB = P_TRUE_FB_D;
        elseif curr_block_type == 'P'
            trueFB = P_TRUE_FB_P;
        else
            error('Unknown block type: %s', curr_block_type);
        end        
        td.trueFB_block(b) = trueFB;

        td.block_structure = subj_structure;
        td.block_type(b) = curr_block_type;
         revTrial = reversal_trials(b);

        stim_order = repmat(1:N_STIM, 1, N_REPS);
        stim_order = stim_order(randperm(N_TRIALS));

        go_nogo = nan(1, N_TRIALS);
        for t = 1:N_TRIALS
            s = stim_order(t);
            go_nogo(t) = double(ternary(t <= revTrial, PRE_REV_RULE(s), POST_REV_RULE(s)));
        end

        nTrue      = round(trueFB * N_TRIALS);
        fb_arr     = [true(1,nTrue), false(1,N_TRIALS-nTrue)];
        trueFB_vec = double(fb_arr(randperm(N_TRIALS)));

        theta     = 0.5 * ones(1, N_STIM);
        n_eff     = ones(1, N_STIM);
        first_obs = true(1, N_STIM);

        for t = 1:N_TRIALS
            s_t     = stim_order(t);
            go_true = go_nogo(t);
            fb_true = trueFB_vec(t);

            p_go    = 1 / (1 + exp(-params.beta * (theta(s_t) - 0.5)));
            resp_go = double(rand < p_go);
            correct = double(resp_go == go_true);
            perc_cor = fb_true * correct + (1-fb_true) * (1-correct);
            y_t = double((resp_go==1 && perc_cor==1) || (resp_go==0 && perc_cor==0));

            % Prospective certainty (before update — matches confidence timing)
            certainty_t = abs(theta(s_t) - 0.5);

            % Nassar learning rate
            pred_prob  = max(theta(s_t)^y_t * (1-theta(s_t))^(1-y_t), 1e-6);
            chi_t      = 0.5 / pred_prob;
            omega_t    = (params.H * chi_t) / (params.H * chi_t + (1-params.H));
            n_eff(s_t) = (1-omega_t) * (n_eff(s_t)+1) + omega_t;
            alpha_t    = omega_t + (1-omega_t) / n_eff(s_t);

            if first_obs(s_t)
                alpha_t       = max(alpha_t, ALPHA_INIT);
                first_obs(s_t) = false;
            end

            delta_t    = y_t - theta(s_t);
            theta(s_t) = max(0.01, min(0.99, theta(s_t) + alpha_t * delta_t));

            conf_raw   = CONF_OFFSET + CONF_GAIN * certainty_t + 0.8*randn;
            confidence = max(1, min(10, round(conf_raw)));

            td.stimID(b,t)           = s_t;
            td.feat1(b,t)            = STIM_FEATURES(s_t,1);
            td.feat2(b,t)            = STIM_FEATURES(s_t,2);
            td.switched(b,t)         = ismember(s_t, SWITCH_STIMS);
            td.goTrial(b,t)          = go_true;
            td.respWasGo(b,t)        = resp_go;
            td.correct(b,t)          = correct;
            td.perceivedCorrect(b,t) = perc_cor;
            td.trueFB(b,t)           = trueFB;
            td.confidence(b,t)       = confidence;
            td.blocknum(b,t)         = b;
            td.trialnum(b,t)         = t;

            mv.theta(b,t,:)   = theta;
            mv.alpha(b,t)     = alpha_t;
            mv.omega(b,t)     = omega_t;
            mv.delta(b,t)     = delta_t;
            mv.certainty(b,t) = certainty_t;
            mv.n_eff(b,t)     = n_eff(s_t);
        end
        td.revTrial(b) = revTrial;
    end

    sim_data.(sl) = td;
    model_vars.(sl) = mv;
end
fprintf('Simulation complete: %d subjects × %d blocks × %d trials\n', ...
    N_SUBJECTS, N_BLOCKS, N_TRIALS);

% =========================================================================
%% §4  SIMULATION FIGURES
% =========================================================================
example_subj = 'sub01';
half_win = 30;
x_ali    = -half_win:half_win;
n_x      = numel(x_ali);

clr_sw = [0.80 0.30 0.10];
clr_mn = [0.15 0.45 0.70];

% ── Fig S1: Belief trajectories for one subject, one block ───────────────
fig1 = figure('Position',[50 50 900 450]);
title('Belief trajectories — sub01, block 1 | solid=maintained, dashed=switched');
hold on;

td  = sim_data.(example_subj);
mv  = model_vars.(example_subj);
rev = td.revTrial(1);
cmap = lines(N_STIM);
for s = 1:N_STIM
    ls  = ternary(ismember(s, SWITCH_STIMS), '--', '-');
    lbl = sprintf('Stim %d (%s)', s, ternary(ismember(s,SWITCH_STIMS),'switch','maintain'));
    plot(1:N_TRIALS, squeeze(mv.theta(1,:,s)), ...
        'Color',cmap(s,:),'LineWidth',1.4,'LineStyle',ls,'DisplayName',lbl);
end
xline(rev,'k--','LineWidth',1.5,'HandleVisibility','off');
yline(0.5,'k:','LineWidth',0.8,'HandleVisibility','off');
xlabel('Trial'); ylabel('P(Go target)');
xlim([1 N_TRIALS]); ylim([0 1]);
legend('Location','east','FontSize',8,'Box','off');
saveas(fig1, fullfile(outpath,'sim_fig1_beliefs.pdf'));

% ── Fig S2: Accuracy aligned to reversal — switch vs maintain ────────────
fig2 = figure('Position',[50 50 700 450]);
hold on;
title('Accuracy aligned to reversal — switched vs maintained stimuli');

sw_acc = NaN(0, n_x);
mn_acc = NaN(0, n_x);

for sl = fieldnames(sim_data)'
    td2 = sim_data.(sl{1});
    for b = 1:N_BLOCKS
        rev2   = td2.revTrial(b);
        sw_row = NaN(1,n_x);
        mn_row = NaN(1,n_x);
        for xi = 1:n_x
            t_abs = round(rev2) + x_ali(xi);
            if t_abs < 1 || t_abs > N_TRIALS, continue; end
            if td2.switched(b,t_abs), sw_row(xi) = td2.correct(b,t_abs);
            else,                     mn_row(xi) = td2.correct(b,t_abs); end
        end
        if any(~isnan(sw_row)), sw_acc(end+1,:) = sw_row; end
        if any(~isnan(mn_row)), mn_acc(end+1,:) = mn_row; end
    end
end

plot_ribbon(gca, x_ali, sw_acc, clr_sw, '--', sprintf('Switched (n=%d)', size(sw_acc,1)));
plot_ribbon(gca, x_ali, mn_acc, clr_mn, '-',  sprintf('Maintained (n=%d)', size(mn_acc,1)));
xline(0,'k--','LineWidth',1.5,'HandleVisibility','off');
yline(0.5,'k:','HandleVisibility','off');
xlabel('Trial relative to reversal'); ylabel('Mean accuracy');
xlim([-half_win half_win]); ylim([0 1]);
legend('Location','south','Box','off','FontSize',9);
saveas(fig2, fullfile(outpath,'sim_fig2_switch_vs_maint.pdf'));

% ── Fig S3: Learning rate aligned to reversal ────────────────────────────
fig3 = figure('Position',[50 50 700 450]);
hold on;
title('Effective learning rate α aligned to reversal');

sw_alpha = NaN(0, n_x);
mn_alpha = NaN(0, n_x);

for sl = fieldnames(sim_data)'
    td2 = sim_data.(sl{1});
    mv2 = model_vars.(sl{1});
    for b = 1:N_BLOCKS
        rev2   = td2.revTrial(b);
        sw_row = NaN(1,n_x);
        mn_row = NaN(1,n_x);
        for xi = 1:n_x
            t_abs = round(rev2) + x_ali(xi);
            if t_abs < 1 || t_abs > N_TRIALS, continue; end
            a_val = mv2.alpha(b, t_abs);
            if isnan(a_val), continue; end
            if td2.switched(b,t_abs), sw_row(xi) = a_val;
            else,                     mn_row(xi) = a_val; end
        end
        if any(~isnan(sw_row)), sw_alpha(end+1,:) = sw_row; end
        if any(~isnan(mn_row)), mn_alpha(end+1,:) = mn_row; end
    end
end

plot_ribbon(gca, x_ali, sw_alpha, clr_sw, '--', 'Switched');
plot_ribbon(gca, x_ali, mn_alpha, clr_mn, '-',  'Maintained');
xline(0,'k--','LineWidth',1.5,'HandleVisibility','off');
xlabel('Trial relative to reversal'); ylabel('\alpha_t');
xlim([-half_win half_win]); ylim([0 1]);
legend('Location','best','Box','off','FontSize',9);

annotation('textbox',[0.12 0.01 0.78 0.08],'String', ...
    ['Note: α appears flat when aligned to reversal trial number because (a) '...
     'stim order is random so the switched stimulus appears ~4 trials after reversal '...
     'on average, (b) pre-reversal n_eff is large so 1/n_eff ≈ H already, and '...
     '(c) averaging across stimuli mixes four independent α trajectories. '...
     'See Fig S4 for stimulus-encounter-aligned α.'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);
saveas(fig3, fullfile(outpath,'sim_fig3_alpha_reversal.pdf'));

% ── Fig S4: α by stimulus encounter (correct alignment) ──────────────────
%
% WHY THIS MATTERS
% Pre-reversal n_eff grows to ~rev/4 ≈ 10-18 per stimulus, so α is already
% near its floor (≈ H). The spike only appears when the switched stimulus
% is first re-encountered after reversal (encounter 0 below). Aligning to
% the reversal trial number smears this across ~4 positions and hides it.

ENC_WIN = 8;
n_enc   = 2*ENC_WIN + 1;
x_enc   = -ENC_WIN:ENC_WIN;

sw_alpha_enc = NaN(0, n_enc);
mn_alpha_enc = NaN(0, n_enc);
sw_omega_enc = NaN(0, n_enc);
mn_omega_enc = NaN(0, n_enc);

for sl = fieldnames(sim_data)'
    td2 = sim_data.(sl{1});
    mv2 = model_vars.(sl{1});
    for b = 1:N_BLOCKS
        rev2 = td2.revTrial(b);
        if isnan(rev2), continue; end
        for s = 1:N_STIM
            s_mask = td2.stimID(b,:) == s;
            [s_trials, sort_i] = sort(find(s_mask));
            s_alpha = mv2.alpha(b, s_trials);
            s_omega = mv2.omega(b, s_trials);
            fp = find(s_trials > rev2, 1);
            if isempty(fp), continue; end
            row_a = NaN(1, n_enc);
            row_o = NaN(1, n_enc);
            for ek = 1:n_enc
                ei = fp + (ek - ENC_WIN - 1);
                if ei >= 1 && ei <= numel(s_trials)
                    row_a(ek) = s_alpha(ei);
                    row_o(ek) = s_omega(ei);
                end
            end
            if ismember(s, SWITCH_STIMS)
                sw_alpha_enc(end+1,:) = row_a;
                sw_omega_enc(end+1,:) = row_o;
            else
                mn_alpha_enc(end+1,:) = row_a;
                mn_omega_enc(end+1,:) = row_o;
            end
        end
    end
end

fig4 = figure('Position',[50 50 1000 450]);
sgtitle({'Learning rate and ω — aligned to first post-reversal encounter per stimulus', ...
    '(x = 0 is the first time each stimulus appears after the reversal)'}, 'FontSize',10);

ax4a = subplot(1,2,1); hold(ax4a,'on');
title(ax4a,'Effective learning rate α');
plot_ribbon(ax4a, x_enc, sw_alpha_enc, clr_sw, '--', sprintf('Switched (n=%d)', size(sw_alpha_enc,1)));
plot_ribbon(ax4a, x_enc, mn_alpha_enc, clr_mn, '-',  sprintf('Maintained (n=%d)', size(mn_alpha_enc,1)));
xline(ax4a,-0.5,'k--','LineWidth',1.5,'HandleVisibility','off');
xlabel(ax4a,'Encounter relative to reversal'); ylabel(ax4a,'\alpha_t');
legend(ax4a,'Box','off','Location','best');

ax4b = subplot(1,2,2); hold(ax4b,'on');
title(ax4b,'Change-point probability ω');
plot_ribbon(ax4b, x_enc, sw_omega_enc, clr_sw, '--', 'Switched');
plot_ribbon(ax4b, x_enc, mn_omega_enc, clr_mn, '-',  'Maintained');
xline(ax4b,-0.5,'k--','LineWidth',1.5,'HandleVisibility','off');
xlabel(ax4b,'Encounter relative to reversal'); ylabel(ax4b,'\omega_t');
legend(ax4b,'Box','off','Location','best');
saveas(fig4, fullfile(outpath,'sim_fig4_alpha_by_encounter.pdf'));

fprintf('Simulation figures saved.\n');
save(fullfile(outpath,'sim_data.mat'), 'sim_data', 'model_vars', 'params');
% results struct saved after §6 (below) once fitting is complete.

% =========================================================================
%% §5  FIT NASSAR TO REAL DATA (all blocks combined)
% =========================================================================
results = fit_all_subjects(all_trial_data, STIM_FEATURES, outpath);

% =========================================================================
%% §6  FIT NASSAR SEPARATELY FOR D AND P BLOCKS
%
% Adds to results.(sn):
%   .H_fit_det / .H_fit_prob     — per-block-type hazard rates
%   .beta_fit_det / .beta_fit_prob
%   .H_noise_sensitivity         — H_prob − H_det (noise-induced volatility)
%   .block_is_det                — logical(1×nBlocks)
%   .block_transition            — cell(1×nBlocks): 'D→D','D→P','P→D','P→P','first'
% =========================================================================
results = fit_subjects_by_blocktype(all_trial_data, STIM_FEATURES, results, outpath);

% Save results struct now that both combined and block-type fits are done.
% F_RL_EEG_analysis.m loads this file.
results_file = fullfile(outpath, 'nassar_results.mat');
save(results_file, 'results');
fprintf('Nassar results saved → %s\n', results_file);

% =========================================================================
%% §7  VISUALISE FITTED PARAMETERS
% =========================================================================
subj_ids  = fieldnames(results);
N_subj    = numel(subj_ids);

H_vals    = nan(N_subj,1);
beta_vals = nan(N_subj,1);
nll_vals  = nan(N_subj,1);
bic_vals  = nan(N_subj,1);
nll_null  = nan(N_subj,1);
H_det_v   = nan(N_subj,1);
H_prob_v  = nan(N_subj,1);
H_ns_v    = nan(N_subj,1);
acc_pre   = nan(N_subj,1);
acc_post  = nan(N_subj,1);
rev_cost  = nan(N_subj,1);
conf_mean = nan(N_subj,1);

for si = 1:N_subj
    sn = subj_ids{si};
    r  = results.(sn);
    H_vals(si)    = r.H_fit;
    beta_vals(si) = r.beta_fit;
    nll_vals(si)  = r.nll;
    bic_vals(si)  = r.bic;
    nll_null(si)  = r.N_obs * log(2);
    if isfield(r,'H_fit_det'),         H_det_v(si)  = r.H_fit_det;           end
    if isfield(r,'H_fit_prob'),        H_prob_v(si) = r.H_fit_prob;          end
    if isfield(r,'H_noise_sensitivity'), H_ns_v(si) = r.H_noise_sensitivity; end
    if isfield(all_trial_data, sn)
        [acc_pre(si), acc_post(si), rev_cost(si), conf_mean(si)] = ...
            extract_beh_summary(all_trial_data.(sn).trial_data);
    end
end

% ── Fig 1: Parameter distributions ───────────────────────────────────────
fig_p1 = figure('Position',[50 50 1100 380]);
sgtitle('Fitted Nassar parameters — empirical data', 'FontSize',11);

ax = subplot(1,3,1); hold(ax,'on');
title(ax,'Combined H (all blocks)');
histogram(ax, H_vals, 10, 'FaceColor',clr_mn,'EdgeColor','w','FaceAlpha',0.8);
xline(ax, mean(H_vals,'omitnan'),   'k-', 'LineWidth',2,'HandleVisibility','off');
xline(ax, median(H_vals,'omitnan'), 'k--','LineWidth',1,'HandleVisibility','off');
xline(ax, 0.02, 'Color',[0.6 0.6 0.6],'LineStyle',':', ...
    'LineWidth',1,'DisplayName','Task true ≈ 0.02');
text(ax,0.98,0.97,sprintf('Mean=%.3f\nMedian=%.3f', ...
    mean(H_vals,'omitnan'),median(H_vals,'omitnan')), ...
    'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
    'FontSize',8,'BackgroundColor','w');
xlabel(ax,'H (hazard rate)'); ylabel(ax,'Subjects');
legend(ax,'Box','off','FontSize',8,'Location','best');

ax = subplot(1,3,2); hold(ax,'on');
title(ax,'H_{det} vs H_{prob}');
ok = ~isnan(H_det_v) & ~isnan(H_prob_v);
scatter(ax, H_det_v(ok), H_prob_v(ok), 60, clr_mn,'filled','MarkerFaceAlpha',0.8);
lims = [0, max([H_det_v(ok); H_prob_v(ok)], [], 'omitnan')+0.02];
plot(ax,lims,lims,'k--','LineWidth',1,'HandleVisibility','off');
[~,p_dp] = ttest(H_prob_v(ok) - H_det_v(ok));
text(ax,0.05,0.97,sprintf('ΔH paired t-test: p=%.3f',p_dp), ...
    'Units','normalized','VerticalAlignment','top','FontSize',8,'BackgroundColor','w');
text(ax,0.05,0.82,'Points above diagonal:', ...
    'Units','normalized','FontSize',7,'Color',[0.5 0.5 0.5],'VerticalAlignment','top');
text(ax,0.05,0.73,'H_{prob} > H_{det}', ...
    'Units','normalized','FontSize',7,'Color',[0.5 0.5 0.5],'VerticalAlignment','top');
xlabel(ax,'H_{det}'); ylabel(ax,'H_{prob}');
xlim(ax,lims); ylim(ax,lims); axis(ax,'square');

ax = subplot(1,3,3); hold(ax,'on');
title(ax,'Noise sensitivity ΔH = H_{prob} − H_{det}');
ok_ns = ~isnan(H_ns_v);
histogram(ax, H_ns_v(ok_ns), 10, 'FaceColor',clr_sw,'EdgeColor','w','FaceAlpha',0.8);
xline(ax,0,'k--','LineWidth',1.5,'HandleVisibility','off');
[~,p_ns] = ttest(H_ns_v(ok_ns));
text(ax,0.98,0.97,sprintf('Mean ΔH=%.3f\nt-test vs 0: p=%.3f', ...
    mean(H_ns_v(ok_ns),'omitnan'),p_ns), ...
    'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
    'FontSize',8,'BackgroundColor','w');
xlabel(ax,'\Delta H (noise sensitivity)'); ylabel(ax,'Subjects');
saveas(fig_p1, fullfile(outpath,'fig1_param_distributions.pdf'));

% ── Fig 2: Parameters vs behaviour ───────────────────────────────────────
fig_p2 = figure('Position',[50 50 1100 420]);
sgtitle('Fitted parameters vs behavioural measures', 'FontSize',11);

pairs = {
    H_vals,    acc_post-acc_pre, 'H',              'Acc change pre→post rev';
    H_vals,    rev_cost,         'H',              'Reversal cost (acc drop)';
    beta_vals, acc_pre,          '\beta',          'Pre-reversal accuracy';
    H_ns_v,    conf_mean,        '\Delta H',       'Mean confidence';
};

for pi = 1:4
    xd  = pairs{pi,1};
    yd  = pairs{pi,2};
    xl  = pairs{pi,3};
    yl  = pairs{pi,4};
    ax  = subplot(1,4,pi); hold(ax,'on');
    title(ax, yl, 'FontSize',8);
    ok2 = ~isnan(xd) & ~isnan(yd);
    scatter(ax, xd(ok2), yd(ok2), 50, clr_mn,'filled','MarkerFaceAlpha',0.7);
    if sum(ok2) > 3
        [rv,pv] = corr(xd(ok2), yd(ok2),'Rows','complete');
        xi2 = linspace(min(xd(ok2)),max(xd(ok2)),100);
        plot(ax,xi2,polyval(polyfit(xd(ok2),yd(ok2),1),xi2),'k-','LineWidth',1.5,'HandleVisibility','off');
        text(ax,0.05,0.97,sprintf('r=%.2f\np=%.3f',rv,pv), ...
            'Units','normalized','VerticalAlignment','top','FontSize',8,'BackgroundColor','w');
    end
    xlabel(ax,xl,'Interpreter','tex'); ylabel(ax,yl,'FontSize',8);
end
saveas(fig_p2, fullfile(outpath,'fig2_param_vs_behaviour.pdf'));

% ── Fig 3: Model fit quality ──────────────────────────────────────────────
fig_p3 = figure('Position',[50 50 850 380]);
sgtitle('Nassar model fit quality', 'FontSize',11);

ax = subplot(1,2,1); hold(ax,'on');
title(ax,'NLL: Nassar vs null model');
scatter(ax, nll_null, nll_vals, 60, clr_mn,'filled','MarkerFaceAlpha',0.8);
lims3 = [min([nll_null;nll_vals])-10, max([nll_null;nll_vals])+10];
plot(ax,lims3,lims3,'k--','LineWidth',1,'HandleVisibility','off');
text(ax,0.05,0.95,'← Model better','Units','normalized', ...
    'VerticalAlignment','top','FontSize',8,'Color',[0 0.5 0]);
pR2 = 1 - nll_vals ./ nll_null;
text(ax,0.98,0.03,sprintf('Mean pseudo-R²=%.3f',mean(pR2,'omitnan')), ...
    'Units','normalized','HorizontalAlignment','right','FontSize',8,'BackgroundColor','w');
xlabel(ax,'Null NLL'); ylabel(ax,'Nassar NLL');
xlim(ax,lims3); ylim(ax,lims3);

ax = subplot(1,2,2); hold(ax,'on');
title(ax,'BIC per subject');
[bic_sorted, sort_idx] = sort(bic_vals,'ascend');
bar(ax,1:N_subj,bic_sorted,'FaceColor',clr_mn,'EdgeColor','none');
set(ax,'XTick',1:N_subj,'XTickLabel',subj_ids(sort_idx),'XTickLabelRotation',45,'FontSize',7);
yline(ax,mean(bic_vals,'omitnan'),'r--','LineWidth',1.5,'HandleVisibility','off');
text(ax,0.02,0.97,sprintf('Mean BIC=%.1f',mean(bic_vals,'omitnan')), ...
    'Units','normalized','VerticalAlignment','top','FontSize',8,'Color','r');
xlabel(ax,'Subject (sorted)'); ylabel(ax,'BIC');
saveas(fig_p3, fullfile(outpath,'fig3_model_fit.pdf'));

% ── Fig 4: Prospective certainty and surprise vs confidence ──────────────
fig_p4 = figure('Position',[50 50 1100 380]);
sgtitle('Model uncertainty signals vs confidence ratings', 'FontSize',11);

measures    = {'certainty_trial','surprise'};
meas_labels = {'|θ−0.5| (prospective, pre-feedback)', 'ω×|δ| (retrospective, post-feedback)'};
r_conf_vals = nan(N_subj,2);

for mi = 1:2
    all_x = []; all_c = [];
    for si = 1:N_subj
        sn = subj_ids{si};
        r  = results.(sn);
        if ~isfield(r,measures{mi}) || ~isfield(all_trial_data,sn), continue; end
        td3 = all_trial_data.(sn).trial_data;
        if ~isfield(td3,'confidence'), continue; end
        [nB3,nT3] = get_n_blocks_trials(td3);
        sx = []; sc = [];
        bm_r = r.block_id;
        for t_idx = 1:numel(r.trial_id)
            b3   = bm_r(t_idx);
            t3   = r.trial_id(t_idx);
            if t3 < 1 || t3 > nT3, continue; end
            cv = td3.confidence(b3,t3);
            xv = r.(measures{mi})(t_idx);
            if isnan(cv) || isnan(xv), continue; end
            sx(end+1) = xv; sc(end+1) = cv;
        end
        if numel(sx) > 5
            [r_conf_vals(si,mi),~] = corr(sx(:),sc(:),'Rows','complete','Type','Spearman');
            all_x = [all_x, sx]; all_c = [all_c, sc];
        end
    end

    ax = subplot(1,4,(mi-1)*2+1); hold(ax,'on');
    scatter(ax,all_x,all_c,3,[0.75 0.75 0.75],'filled','MarkerFaceAlpha',0.12,'HandleVisibility','off');
    edges4 = linspace(min(all_x),max(all_x),11);
    for bi4 = 1:10
        bm4 = all_x >= edges4(bi4) & all_x < edges4(bi4+1);
        if sum(bm4) > 5
            plot(ax,mean(all_x(bm4)),mean(all_c(bm4)),'o','Color',clr_mn, ...
                'MarkerSize',7,'MarkerFaceColor',clr_mn,'HandleVisibility','off');
        end
    end
    [rp,pp] = corr(all_x(:),all_c(:),'Rows','complete','Type','Spearman');
    text(ax,0.98,0.03,sprintf('ρ=%.2f, p=%.3f',rp,pp), ...
        'Units','normalized','HorizontalAlignment','right','FontSize',8,'BackgroundColor','w');
    xlabel(ax,meas_labels{mi},'FontSize',7); ylabel(ax,'Confidence (1-10)');
    title(ax,ternary(mi==1,'Prospective certainty','Retrospective surprise'));

    ax2 = subplot(1,4,(mi-1)*2+2); hold(ax2,'on');
    ok_m = ~isnan(r_conf_vals(:,mi));
    histogram(ax2,r_conf_vals(ok_m,mi),10,'FaceColor',clr_mn,'EdgeColor','w');
    xline(ax2,0,'k--');
    xline(ax2,mean(r_conf_vals(ok_m,mi),'omitnan'),'r-','LineWidth',2,'HandleVisibility','off');
    [~,pt] = ttest(r_conf_vals(ok_m,mi));
    text(ax2,0.98,0.97,sprintf('Mean ρ=%.2f\nt-test: p=%.3f', ...
        mean(r_conf_vals(ok_m,mi),'omitnan'),pt), ...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
        'FontSize',8,'BackgroundColor','w');
    xlabel(ax2,'ρ per subject'); ylabel(ax2,'Subjects');
    title(ax2,ternary(mi==1,'Certainty–conf. r','Surprise–conf. r'));
end

annotation('textbox',[0.01 0.01 0.98 0.07],'String', ...
    ['Trial order: stimulus → response → CONFIDENCE RATING → feedback. '...
     '|θ−0.5| is computed before feedback: the correct causal predictor. '...
     'ω×|δ| requires the outcome (retrospective control). '...
     'If surprise predicts better than certainty, participants may integrate feedback retrospectively.'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);
saveas(fig_p4, fullfile(outpath,'fig4_uncertainty_vs_confidence.pdf'));

% ── Fig 5: Per-subject diagnostic pages ──────────────────────────────────
fprintf('Generating per-subject diagnostic pages...\n');
for si = 1:N_subj
    sn = subj_ids{si};
    r  = results.(sn);
    if ~isfield(all_trial_data,sn), continue; end
    td4 = all_trial_data.(sn).trial_data;
    if ~isfield(td4,'revTrial'), continue; end
    [nB4,~] = get_n_blocks_trials(td4);

    fig_d = figure('Position',[50 50 1400 560],'Visible','off');
    sgtitle(sprintf('%s  H=%.3f  β=%.2f  ΔH=%.3f', sn, r.H_fit, r.beta_fit, ...
        ternary(isfield(r,'H_noise_sensitivity'), r.H_noise_sensitivity, NaN)), ...
        'Interpreter','none','FontSize',10);

    for b = 1:min(nB4,5)
        bm4   = r.block_id == b;
        b_tr  = r.trial_id(bm4);
        b_al  = r.alpha_trial(bm4);
        b_om  = r.omega_trial(bm4);
        b_sp  = r.surprise(bm4);
        b_st  = r.stim_id(bm4);
        rev4  = td4.revTrial(b);

        % SWITCH_STIMS read from group_T in §2 — authoritative source.
        sw_b = SWITCH_STIMS;
        is_det_b = isfield(r,'block_is_det') && b <= numel(r.block_is_det) && r.block_is_det(b);
        blk_lbl  = ternary(is_det_b,'D','P');
        trans_lbl = '';
        if isfield(r,'block_transition') && b <= numel(r.block_transition)
            trans_lbl = r.block_transition{b};
        end

        cmap4 = lines(N_STIM);
        ax = subplot(2,min(nB4,5),b); hold(ax,'on');
        for s = 1:N_STIM
            sm4 = b_st == s;
            scatter(ax,b_tr(sm4),b_al(sm4),15,cmap4(s,:),'filled','MarkerFaceAlpha',0.6,'HandleVisibility','off');
        end
        if ~isnan(rev4), xline(ax,rev4,'k--','HandleVisibility','off'); end
        axis(ax,'square');
        title(ax,sprintf('B%d [%s] %s  α_t',b,blk_lbl,trans_lbl),'FontSize',7);
        xlabel(ax,'Trial'); ylabel(ax,'\alpha_t'); ylim(ax,[0 1]);

        ax2 = subplot(2,min(nB4,5),b+min(nB4,5)); hold(ax2,'on');
        axis(ax2,'square');
        yyaxis(ax2,'left');
        plot(ax2,b_tr,b_om,'Color',[0.5 0.5 0.5],'LineWidth',0.8);
        ylabel(ax2,'\omega_t','Color',[0.5 0.5 0.5]);
        yyaxis(ax2,'right');
        plot(ax2,b_tr,b_sp,'Color',clr_sw,'LineWidth',1.2);
        ylabel(ax2,'Surprise','Color',clr_sw);
        if ~isnan(rev4), xline(ax2,rev4,'k--','HandleVisibility','off'); end
        title(ax2,sprintf('B%d ω_t / surprise',b),'FontSize',7);
        xlabel(ax2,'Trial');
    end
    saveas(fig_d, fullfile(outpath, sprintf('fig5_diag_%s.pdf',sn)));
    close(fig_d);
end
fprintf('Diagnostic pages saved.\n');

% =========================================================================
%% §8  PARAMETER RECOVERY
%
% Correct direction: fitted params → simulate → refit → compare.
% Tests identifiability in the empirical parameter range.
% =========================================================================
fprintf('\n=== §8: Parameter recovery ===\n');

H_fit_all      = nan(N_subj,1);
beta_fit_all   = nan(N_subj,1);
H_recovered    = nan(N_subj,1);
beta_recovered = nan(N_subj,1);

H_grid_r    = logspace(log10(0.005), log10(0.40), 30);
beta_grid_r = linspace(1, 20, 30);

for si = 1:N_subj
    sn = subj_ids{si};
    r  = results.(sn);
    H_fit_all(si)    = r.H_fit;
    beta_fit_all(si) = r.beta_fit;

    if isfield(all_trial_data, sn)
        td_r = all_trial_data.(sn).trial_data;
        [nB_r, nT_r] = get_n_blocks_trials(td_r);
        rev_vec = ternary(isfield(td_r,'revTrial'), td_r.revTrial, ...
            randi([REV_MIN,REV_MAX],1,nB_r));
    else
        nB_r = N_BLOCKS; nT_r = N_TRIALS;
        rev_vec = randi([REV_MIN,REV_MAX],1,nB_r);
    end

    td_syn = simulate_one_subject(r.H_fit, r.beta_fit, nB_r, nT_r, ...
        rev_vec, N_STIM, N_REPS, P_TRUE_FB, PRE_REV_RULE, POST_REV_RULE, ...
        SWITCH_STIMS, ALPHA_INIT);

    obs_syn = pack_observations(td_syn, 'stimID', true);
    nll_g   = nan(numel(H_grid_r), numel(beta_grid_r));
    for hi = 1:numel(H_grid_r)
        for bi2 = 1:numel(beta_grid_r)
            nll_g(hi,bi2) = nassar_nll([H_grid_r(hi), beta_grid_r(bi2)], obs_syn, STIM_FEATURES);
        end
    end
    [~,bl] = min(nll_g(:)); [hi0,bi0] = ind2sub(size(nll_g),bl);
    x0 = [H_grid_r(hi0), beta_grid_r(bi0)];
    opts = optimoptions('fmincon','Display','off','MaxIterations',500);
    try
        xr = fmincon(@(x) nassar_nll(x,obs_syn,STIM_FEATURES), x0,[],[],[],[], ...
            [0.001,0.5],[0.40,25],[],opts);
    catch
        xr = x0;
    end
    H_recovered(si)    = xr(1);
    beta_recovered(si) = xr(2);
    fprintf('  %s: H %.3f→%.3f | β %.2f→%.2f\n', sn, r.H_fit, H_recovered(si), ...
        r.beta_fit, beta_recovered(si));
end

fig_rec = figure('Position',[50 50 850 380]);
sgtitle('Parameter recovery — fitted params as generating params', 'FontSize',11);

ax = subplot(1,2,1); hold(ax,'on');
axis(ax,'square');
title(ax,'Hazard rate H');
ok_r = ~isnan(H_fit_all) & ~isnan(H_recovered);
scatter(ax, H_fit_all(ok_r), H_recovered(ok_r), 60, clr_mn,'filled','MarkerFaceAlpha',0.8);
lims_r = [0, max([H_fit_all(ok_r); H_recovered(ok_r)])+0.02];
plot(ax,lims_r,lims_r,'k--','LineWidth',1.2,'HandleVisibility','off');
[r_h,p_h] = corr(H_fit_all(ok_r), H_recovered(ok_r),'Rows','complete');
text(ax,0.05,0.95,sprintf('r=%.2f, p=%.3f',r_h,p_h), ...
    'Units','normalized','VerticalAlignment','top','FontSize',9,'BackgroundColor','w');
xlabel(ax,'H fitted (real data)'); ylabel(ax,'H recovered (synthetic)');
xlim(ax,lims_r); ylim(ax,lims_r); axis(ax,'square');

ax = subplot(1,2,2); hold(ax,'on');
axis(ax,'square');
title(ax,'Inverse temperature β');
ok_rb = ~isnan(beta_fit_all) & ~isnan(beta_recovered);
scatter(ax, beta_fit_all(ok_rb), beta_recovered(ok_rb), 60, clr_mn,'filled','MarkerFaceAlpha',0.8);
lims_b = [0, max([beta_fit_all(ok_rb); beta_recovered(ok_rb)])+1];
plot(ax,lims_b,lims_b,'k--','LineWidth',1.2,'HandleVisibility','off');
[r_b,p_b] = corr(beta_fit_all(ok_rb), beta_recovered(ok_rb),'Rows','complete');
text(ax,0.05,0.95,sprintf('r=%.2f, p=%.3f',r_b,p_b), ...
    'Units','normalized','VerticalAlignment','top','FontSize',9,'BackgroundColor','w');
xlabel(ax,'\beta fitted (real data)'); ylabel(ax,'\beta recovered (synthetic)');
xlim(ax,lims_b); ylim(ax,lims_b); axis(ax,'square');
saveas(fig_rec, fullfile(outpath,'fig6_param_recovery.pdf'));

% =========================================================================
%% §9  MODEL COMPARISON: NASSAR vs SUBJECT-WISE RW
%
% The RW model (fit_RW_subjectwise) fits separate α_det and α_prob per
% subject plus a confidence-weighted pre/post-reversal variant.
% BIC penalises complexity; Nassar and the basic RW both have k=2 free
% parameters (H/α, β), so ΔBIC = 2×ΔNLL directly.
%
% The RW BIC is derived from the trial-level Bernoulli likelihood of the
% same observed response sequences used for Nassar, ensuring comparability.
%
% NOTE: fit_RW_subjectwise adds columns to group_T. For BIC comparison we
% use the per-trial likelihood computed inside rw_nll (fixed-α version),
% which matches the Nassar likelihood exactly. The richer confidence-
% weighted model is not compared here because it uses a different objective.
% =========================================================================
fprintf('\n=== §9: Model comparison (Nassar vs RW) ===\n');

subjects_comp = fieldnames(all_trial_data);
N_comp        = numel(subjects_comp);
bic_nassar    = nan(N_comp,1);
bic_rw        = nan(N_comp,1);
nll_nassar    = nan(N_comp,1);
nll_rw        = nan(N_comp,1);

alpha_grid = linspace(0.01, 0.60, 30);
beta_grid  = linspace(2, 15, 30);

for si = 1:N_comp
    sn = subjects_comp{si};
    td5 = all_trial_data.(sn).trial_data;
    if isfield(td5,'stimType'), sf = 'stimType';
    elseif isfield(td5,'stimID'), sf = 'stimID';
    else, continue; end
    has_perc = isfield(td5,'perceivedCorrect') && ~all(isnan(td5.perceivedCorrect(:)));
    obs5 = pack_observations(td5, sf, has_perc);
    if isempty(obs5.stim_id), continue; end
    N_obs5 = numel(obs5.stim_id);

    % Nassar BIC (already fitted)
    if isfield(results, sn)
        bic_nassar(si)  = results.(sn).bic;
        nll_nassar(si)  = results.(sn).nll;
    end

    % RW: grid search then refine (fixed α for direct likelihood comparison)
    nll_rw_grid = nan(numel(alpha_grid), numel(beta_grid));
    for ai = 1:numel(alpha_grid)
        for bi2 = 1:numel(beta_grid)
            nll_rw_grid(ai,bi2) = rw_nll([alpha_grid(ai), beta_grid(bi2)], obs5, STIM_FEATURES);
        end
    end
    [~,bl5] = min(nll_rw_grid(:)); [ai0,bi0] = ind2sub(size(nll_rw_grid),bl5);
    x0_rw = [alpha_grid(ai0), beta_grid(bi0)];
    opts5 = optimoptions('fmincon','Display','off','MaxIterations',500);
    try
        xf_rw = fmincon(@(x) rw_nll(x,obs5,STIM_FEATURES), x0_rw,[],[],[],[], ...
            [0.01,0.5],[0.60,30],[],opts5);
    catch
        xf_rw = x0_rw;
    end
    nll_rw(si) = rw_nll(xf_rw, obs5, STIM_FEATURES);
    bic_rw(si) = 2*nll_rw(si) + 2*log(N_obs5);   % k=2 same as Nassar

    fprintf('  %s: Nassar BIC=%.1f  RW BIC=%.1f  ΔBIC=%.1f\n', ...
        sn, bic_nassar(si), bic_rw(si), bic_nassar(si)-bic_rw(si));
end

ok_mc  = ~isnan(bic_nassar) & ~isnan(bic_rw);
d_bic  = bic_nassar(ok_mc) - bic_rw(ok_mc);   % negative = Nassar wins
[~,p_mc,~,st_mc] = ttest(d_bic);

fprintf('\n=== MODEL COMPARISON SUMMARY ===\n');
fprintf('Nassar: mean BIC = %.1f ± %.1f\n', mean(bic_nassar(ok_mc),'omitnan'), std(bic_nassar(ok_mc),'omitnan'));
fprintf('RW:     mean BIC = %.1f ± %.1f\n', mean(bic_rw(ok_mc),'omitnan'),     std(bic_rw(ok_mc),'omitnan'));
fprintf('ΔBIC (Nassar − RW): mean=%.1f  t(%d)=%.2f  p=%.4f\n', ...
    mean(d_bic), st_mc.df, st_mc.tstat, p_mc);
fprintf('Nassar wins in %d/%d subjects\n', sum(d_bic < 0), sum(ok_mc));

fig_mc = figure('Position',[50 50 1100 420]);
sgtitle('Model comparison: Nassar vs Rescorla-Wagner', 'FontSize',11);

ax = subplot(1,3,1); hold(ax,'on');
axis(ax,'square');
title(ax,'BIC per subject');
[~,sort_i5] = sort(bic_nassar(ok_mc),'ascend');
plot(ax,1:sum(ok_mc), bic_nassar(ok_mc(sort_i5)),'o-','Color',clr_mn,'LineWidth',1.5, ...
    'MarkerFaceColor',clr_mn,'DisplayName','Nassar');
plot(ax,1:sum(ok_mc), bic_rw(ok_mc(sort_i5)),    's-','Color',clr_sw,'LineWidth',1.5, ...
    'MarkerFaceColor',clr_sw,'DisplayName','RW');
xlabel(ax,'Subject (sorted by Nassar BIC)'); ylabel(ax,'BIC');
legend(ax,'Box','off','FontSize',8,'Location','northwest');

ax = subplot(1,3,2); hold(ax,'on');
axis(ax,'square');
title(ax,'ΔBIC: Nassar − RW  (< 0 = Nassar wins)');
histogram(ax, d_bic, 10,'FaceColor',clr_mn,'EdgeColor','w','FaceAlpha',0.8);
xline(ax,0,'k--','LineWidth',1.5,'HandleVisibility','off');
xline(ax,-10,'r:','LineWidth',1,'HandleVisibility','off');
xline(ax,10, 'r:','LineWidth',1,'HandleVisibility','off');
xline(ax,mean(d_bic,'omitnan'),'k-','LineWidth',2,'HandleVisibility','off');
text(ax,0.02,0.97,sprintf('Mean ΔBIC=%.1f\nt(%d)=%.2f  p=%.4f\nNassar wins: %d/%d', ...
    mean(d_bic,'omitnan'),st_mc.df,st_mc.tstat,p_mc,sum(d_bic<0),sum(ok_mc)), ...
    'Units','normalized','VerticalAlignment','top','FontSize',8,'BackgroundColor','w');
text(ax,0.5,0.08,'← Nassar better | RW better →','Units','normalized', ...
    'HorizontalAlignment','center','FontSize',8,'Color',[0.5 0.5 0.5]);
xlabel(ax,'ΔBIC'); ylabel(ax,'Subjects');

ax = subplot(1,3,3); hold(ax,'on');
axis(ax,'square');
title(ax,'Nassar BIC vs RW BIC (per subject)');
scatter(ax, bic_rw(ok_mc), bic_nassar(ok_mc), 60, clr_mn,'filled','MarkerFaceAlpha',0.8);
all_b = [bic_nassar(ok_mc); bic_rw(ok_mc)];
lims_mc = [min(all_b,'omitnan')-20, max(all_b,'omitnan')+20];
plot(ax,lims_mc,lims_mc,'k--','HandleVisibility','off');
text(ax,0.05,0.95,'← Nassar better','Units','normalized', ...
    'VerticalAlignment','top','FontSize',8,'Color',[0 0.5 0]);
text(ax,0.95,0.05,'RW better →','Units','normalized', ...
    'HorizontalAlignment','right','FontSize',8,'Color',[0.7 0 0]);
xlabel(ax,'RW BIC'); ylabel(ax,'Nassar BIC');
xlim(ax,lims_mc); ylim(ax,lims_mc); axis(ax,'square');
saveas(fig_mc, fullfile(outpath,'fig7_model_comparison.pdf'));

% =========================================================================
%% §10 POSTERIOR PREDICTIVE CHECKS (Nassar only)
%
% For each subject: simulate N_PPC_SIMS synthetic datasets from fitted H, β.
% Compute 6 summary statistics on real and simulated data.
% Assess whether the model reproduces qualitative features of behaviour.
%
% Statistics:
%   1. Pre-reversal accuracy
%   2. Post-reversal accuracy (first 10 trials)
%   3. Reversal cost (pre − post)
%   4. Recovery rate (slope of accuracy, trials +11 to +30)
%   5. Confidence–accuracy Spearman ρ
%   6. Maintained − switched accuracy gap (post-reversal)
% =========================================================================
fprintf('\n=== §10: Posterior predictive checks ===\n');

N_PPC_SIMS  = 40;
stat_names  = {'Pre-rev accuracy','Post-rev accuracy','Reversal cost', ...
               'Recovery rate','Conf–acc correlation','Maintain–switch gap'};
N_stats     = numel(stat_names);
subjects_ppc = fieldnames(results);
N_ppc        = numel(subjects_ppc);

real_stats      = nan(N_ppc, N_stats);
sim_stats       = nan(N_ppc, N_stats);
sim_stats_ci_lo = nan(N_ppc, N_stats);
sim_stats_ci_hi = nan(N_ppc, N_stats);

for si = 1:N_ppc
    sn = subjects_ppc{si};
    r  = results.(sn);
    if ~isfield(all_trial_data, sn), continue; end
    td6 = all_trial_data.(sn).trial_data;
    if ~isfield(td6,'revTrial'), continue; end
    [nB6, nT6] = get_n_blocks_trials(td6);

    real_stats(si,:) = compute_summary_stats(td6, nB6, nT6, SWITCH_STIMS);

    sim_draws = nan(N_PPC_SIMS, N_stats);
    rev_vec6 = td6.revTrial;
    for draw = 1:N_PPC_SIMS
        td_s = simulate_one_subject(r.H_fit, r.beta_fit, nB6, nT6, rev_vec6, ...
            N_STIM, N_REPS, P_TRUE_FB, PRE_REV_RULE, POST_REV_RULE, SWITCH_STIMS, ALPHA_INIT);
        sim_draws(draw,:) = compute_summary_stats(td_s, nB6, nT6, SWITCH_STIMS);
    end
    sim_stats(si,:)       = mean(sim_draws, 1,'omitnan');
    sim_stats_ci_lo(si,:) = prctile(sim_draws, 5,  1);
    sim_stats_ci_hi(si,:) = prctile(sim_draws, 95, 1);
    fprintf('  %s: pre-acc real=%.2f  sim=%.2f\n', sn, real_stats(si,1), sim_stats(si,1));
end

fprintf('\n%-28s  %-8s  %-8s  %-8s\n','Statistic','r','p','RMSE');
ppc_r = nan(N_stats,1); ppc_p = nan(N_stats,1); ppc_rmse = nan(N_stats,1);
for st = 1:N_stats
    ok_s = ~isnan(real_stats(:,st)) & ~isnan(sim_stats(:,st));
    if sum(ok_s) > 3
        [ppc_r(st),ppc_p(st)] = corr(real_stats(ok_s,st), sim_stats(ok_s,st));
        ppc_rmse(st) = sqrt(mean((real_stats(ok_s,st)-sim_stats(ok_s,st)).^2,'omitnan'));
    end
    fprintf('%-28s  %-8.3f  %-8.4f  %-8.4f\n', stat_names{st}, ppc_r(st), ppc_p(st), ppc_rmse(st));
end

fig_ppc = figure('Position',[50 50 1400 820]);
sgtitle({'Posterior predictive check — Nassar model', ...
    'Each point = one subject | diagonal = perfect prediction'}, 'FontSize',11);

for st = 1:N_stats
    ax = subplot(2,3,st); hold(ax,'on');
    title(ax, stat_names{st},'FontSize',10);
    ok_s = ~isnan(real_stats(:,st)) & ~isnan(sim_stats(:,st));

    % 90% CI bars
    for sj = 1:N_ppc
        if ~ok_s(sj), continue; end
        line(ax,[real_stats(sj,st) real_stats(sj,st)], ...
            [sim_stats_ci_lo(sj,st) sim_stats_ci_hi(sj,st)], ...
            'Color',[0.75 0.75 0.75],'LineWidth',0.8,'HandleVisibility','off');
    end
    scatter(ax, real_stats(ok_s,st), sim_stats(ok_s,st), 60, clr_mn, ...
        'filled','MarkerFaceAlpha',0.8);

    all_v = [real_stats(ok_s,st); sim_stats(ok_s,st)];
    lp = [min(all_v)-0.05*range(all_v), max(all_v)+0.05*range(all_v)];
    if range(lp) < 1e-6, lp = lp + [-0.1 0.1]; end
    plot(ax,lp,lp,'k--','LineWidth',1,'HandleVisibility','off');

    if ~isnan(ppc_r(st))
        border_clr = interp1([0 0.5 1],[0.8 0.2 0.2; 0.9 0.7 0.1; 0.1 0.6 0.1], ...
            min(max(ppc_r(st),0),1));
        set(ax,'XColor',border_clr,'YColor',border_clr,'LineWidth',1.5);
        text(ax,0.04,0.97,sprintf('r=%.2f, p=%.3f\nRMSE=%.3f', ...
            ppc_r(st),ppc_p(st),ppc_rmse(st)), ...
            'Units','normalized','VerticalAlignment','top','FontSize',8, ...
            'BackgroundColor','w','Margin',2);
    end
    xlabel(ax,'Real data','FontSize',9);
    ylabel(ax,sprintf('Simulated (mean of %d)',N_PPC_SIMS),'FontSize',9);
    if ~isempty(lp), xlim(ax,lp); ylim(ax,lp); end
    axis(ax,'square');
end
annotation('textbox',[0.01 0.01 0.98 0.05],'String', ...
    ['PPC: simulate ' num2str(N_PPC_SIMS) ' datasets per subject from fitted H and β. '...
     'Green border = good (r>0.7). Orange = moderate. Red = poor. '...
     'Error bars = 90% simulation CI. A model with low BIC but poor PPC '...
     'is fitting the wrong aspects of behaviour.'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);
saveas(fig_ppc, fullfile(outpath,'fig8_ppc.pdf'));

% =========================================================================
%% §11 ADD NASSAR LATENTS TO BEHAVIOURAL GROUP_T
%
% Behaviour-only output.  The EEG-feature table is deliberately not touched
% here.  The separate script F_RL_EEG_analysis_from_features.m loads the EEG
% combined feature table and joins these same latents there.
% =========================================================================
fprintf('\n=== §11: Joining Nassar latents to behaviour-only group_T ===\n');

if isempty(group_T)
    fprintf('  No group_T was loaded, so behaviour-enriched table was skipped.\n');
    out_file = '';
else
    if ~ismember('switch_stims', group_T.Properties.VariableNames)
        fprintf('  switch_stims column not found in behaviour table; adding default [%s] to all rows.\n', num2str(SWITCH_STIMS));
        group_T.switch_stims = repmat({SWITCH_STIMS}, height(group_T), 1);
    else
        fprintf('  switch_stims column found in behaviour table.\n');
    end

    group_T = add_nassar_to_group_table(group_T, results);

    out_file = fullfile(data_path, 'behav_table_June2026_RL.mat');
    save(out_file, 'group_T');
    save(fullfile(outpath, 'behaviour_RL_table.mat'), 'group_T');

    fprintf('\n========================================================\n');
    fprintf('§11 BEHAVIOUR TABLE SAVE COMPLETE\n');
    fprintf('  File: %s\n', out_file);
    fprintf('  Rows: %d  |  Cols: %d\n', height(group_T), width(group_T));
    fprintf('  Nassar columns added:\n');
    fprintf('    PE_nassar      — δ_t = y_t − θ_t (prediction error)\n');
    fprintf('    omega          — ω_t (change-point probability)\n');
    fprintf('    alpha_nassar   — α_t (effective learning rate)\n');
    fprintf('    certainty      — |θ_t − 0.5| (prospective certainty)\n');
    fprintf('    surprise       — ω_t × |δ_t| (retrospective surprise)\n');
    fprintf('    theta_nassar   — θ_t (current belief state)\n');
    fprintf('========================================================\n');
end

fprintf('\n=== Behaviour-only RL script complete. ===\n');
fprintf('Figures        → %s\n', outpath);
fprintf('Nassar results → %s\n', fullfile(outpath,'nassar_results.mat'));
fprintf('Sim data       → %s\n', fullfile(outpath,'sim_data.mat'));
if ~isempty(out_file), fprintf('Behaviour table → %s\n', out_file); end

% =========================================================================
%% LOCAL FUNCTIONS
% =========================================================================

function td = init_trial_data(nBlocks, nTrials)
fields = {'stimID','feat1','feat2','switched','goTrial','respWasGo','correct', ...
          'perceivedCorrect','trueFB','confidence','blocknum','trialnum'};
for f = fields
    td.(f{1}) = nan(nBlocks, nTrials);
end
td.revTrial = nan(1, nBlocks);
end

function mv = init_model_vars(nBlocks, nTrials, nStim)
mv.theta    = nan(nBlocks, nTrials, nStim);
mv.alpha    = nan(nBlocks, nTrials);
mv.omega    = nan(nBlocks, nTrials);
mv.delta    = nan(nBlocks, nTrials);
mv.certainty = nan(nBlocks, nTrials);
mv.n_eff    = nan(nBlocks, nTrials);
end

function obs = pack_observations(td, stim_field, has_perc)
[nBlocks, nTrials] = get_n_blocks_trials(td);
stim_id  = []; block_id = []; trial_id = []; resp_go = []; y_t_vec = []; is_det = [];

for b = 1:nBlocks
    if isfield(td,'trueFB')
        pfb_b = td.trueFB(b, ~isnan(td.trueFB(b,:)));
        block_det = ~isempty(pfb_b) && mean(pfb_b) >= 0.99;
    else
        block_det = false;
    end

    for t = 1:nTrials
        cor  = td.correct(b,t);
        gotr = td.goTrial(b,t);
        if isnan(cor) || isnan(gotr), continue; end
        stim_val = td.(stim_field)(b,t);
        if isnan(stim_val), continue; end

        rgo = gotr * cor + (1-gotr) * (1-cor);

        if has_perc && isfield(td,'perceivedCorrect') && ~isnan(td.perceivedCorrect(b,t))
            perc_cor = td.perceivedCorrect(b,t);
        else
            perc_cor = cor;
        end

        yt = double((rgo==1 && perc_cor==1) || (rgo==0 && perc_cor==0));

        stim_id(end+1)  = stim_val;
        block_id(end+1) = b;
        trial_id(end+1) = t;
        resp_go(end+1)  = rgo;
        y_t_vec(end+1)  = yt;
        is_det(end+1)   = double(block_det);
    end
end

obs.stim_id  = stim_id;
obs.block_id = block_id;
obs.trial_id = trial_id;
obs.resp_go  = resp_go;
obs.y_t      = y_t_vec;
obs.is_det   = logical(is_det);
end

function [nll, trials_out] = nassar_nll(params, obs, stim_features, return_trials)
if nargin < 4, return_trials = false; end

H    = params(1);
beta = params(2);
N    = numel(obs.stim_id);
N_s  = size(stim_features, 1);

theta      = 0.5 * ones(1, N_s);
n_eff      = ones(1, N_s);
first_obs  = true(1, N_s);
ll         = 0;
prev_block = -1;

if return_trials
    tr_alpha     = nan(1,N); tr_delta    = nan(1,N);
    tr_omega     = nan(1,N); tr_certainty = nan(1,N);
    tr_theta_out = nan(1,N);
end

for t = 1:N
    s_t = obs.stim_id(t);

    if obs.block_id(t) ~= prev_block
        theta      = 0.5 * ones(1, N_s);
        n_eff      = ones(1, N_s);
        first_obs  = true(1, N_s);
        prev_block = obs.block_id(t);
    end

    if isnan(s_t) || s_t < 1 || s_t > N_s, continue; end

    certainty_t = abs(theta(s_t) - 0.5);

    p_go = 1 / (1 + exp(-beta * (theta(s_t) - 0.5)));
    r_go = obs.resp_go(t);
    ll   = ll + log(max(p_go^r_go * (1-p_go)^(1-r_go), 1e-10));

    y_t       = obs.y_t(t);
    pred_prob = max(theta(s_t)^y_t * (1-theta(s_t))^(1-y_t), 1e-6);
    chi_t     = 0.5 / pred_prob;
    omega_t   = (H * chi_t) / (H * chi_t + (1-H));
    n_eff(s_t) = (1-omega_t) * (n_eff(s_t)+1) + omega_t;
    alpha_t   = omega_t + (1-omega_t) / n_eff(s_t);

    if first_obs(s_t)
        alpha_t       = max(alpha_t, 0.5);
        first_obs(s_t) = false;
    end

    delta_t    = y_t - theta(s_t);
    theta(s_t) = max(0.01, min(0.99, theta(s_t) + alpha_t * delta_t));

    if return_trials
        tr_alpha(t)      = alpha_t;
        tr_delta(t)      = delta_t;
        tr_omega(t)      = omega_t;
        tr_certainty(t)  = certainty_t;
        tr_theta_out(t)  = theta(s_t);
    end
end

nll = -ll;
if return_trials
    trials_out.alpha     = tr_alpha;
    trials_out.delta     = tr_delta;
    trials_out.omega     = tr_omega;
    trials_out.certainty = tr_certainty;
    trials_out.theta     = tr_theta_out;
else
    trials_out = [];
end
end

function nll = rw_nll(params, obs, stim_features)
% Fixed learning rate RW — used only for BIC comparison in §9.
% The richer subject-wise RW (fit_RW_subjectwise) is used for group_T columns.
alpha = params(1);
beta  = params(2);
N     = numel(obs.stim_id);
N_s   = size(stim_features, 1);

theta      = 0.5 * ones(1, N_s);
ll         = 0;
prev_block = -1;

for t = 1:N
    s_t = obs.stim_id(t);
    if obs.block_id(t) ~= prev_block
        theta      = 0.5 * ones(1, N_s);
        prev_block = obs.block_id(t);
    end
    if isnan(s_t) || s_t < 1 || s_t > N_s, continue; end

    p_go = 1 / (1 + exp(-beta * (theta(s_t) - 0.5)));
    r_go = obs.resp_go(t);
    ll   = ll + log(max(p_go^r_go * (1-p_go)^(1-r_go), 1e-10));

    y_t        = obs.y_t(t);
    delta_t    = y_t - theta(s_t);
    theta(s_t) = max(0.01, min(0.99, theta(s_t) + alpha * delta_t));
end
nll = -ll;
end

function results = fit_all_subjects(all_trial_data, stim_features, outdir)
% Fit Nassar (H, β) to all subjects across all blocks combined.

subjects  = fieldnames(all_trial_data);
results   = struct();
H_grid    = logspace(log10(0.01), log10(0.30), 25);
beta_grid = linspace(2, 15, 25);

for si = 1:numel(subjects)
    sn = subjects{si};
    fprintf('Fitting %s...\n', sn);
    td = all_trial_data.(sn).trial_data;

    if isfield(td,'stimType'), sf = 'stimType';
    elseif isfield(td,'stimID'), sf = 'stimID';
    else, warning('%s: no stim field, skipping.', sn); continue; end

    has_perc = isfield(td,'perceivedCorrect') && ~all(isnan(td.perceivedCorrect(:)));
    obs = pack_observations(td, sf, has_perc);
    if isempty(obs.stim_id), warning('%s: no valid trials.', sn); continue; end

    % Grid search
    nll_grid = nan(numel(H_grid), numel(beta_grid));
    for hi = 1:numel(H_grid)
        for bi = 1:numel(beta_grid)
            nll_grid(hi,bi) = nassar_nll([H_grid(hi), beta_grid(bi)], obs, stim_features);
        end
    end
    [~,bl] = min(nll_grid(:)); [hi0,bi0] = ind2sub(size(nll_grid),bl);
    x0 = [H_grid(hi0), beta_grid(bi0)];

    opts = optimoptions('fmincon','Display','off','MaxIterations',500);
    try
        x_fit = fmincon(@(x) nassar_nll(x, obs, stim_features), x0,[],[],[],[], ...
            [0.001,0.5],[0.50,30],[],opts);
    catch
        x_fit = x0;
        warning('%s: fmincon failed, using grid best.', sn);
    end

    [nll_fit, trials_out] = nassar_nll(x_fit, obs, stim_features, true);
    N_obs = numel(obs.stim_id);

    % switch_stims comes from group_T (read in §2), not re-derived here.
    results.(sn).subj_id                 = sn;
    results.(sn).H_fit                   = x_fit(1);
    results.(sn).beta_fit                = x_fit(2);
    results.(sn).nll                     = nll_fit;
    results.(sn).bic                     = 2*nll_fit + 2*log(N_obs);
    results.(sn).N_obs                   = N_obs;
    results.(sn).alpha_trial             = trials_out.alpha;
    results.(sn).delta_trial             = trials_out.delta;
    results.(sn).omega_trial             = trials_out.omega;
    results.(sn).certainty_trial         = trials_out.certainty;
    results.(sn).theta_trial             = trials_out.theta;
    results.(sn).surprise                = trials_out.omega .* abs(trials_out.delta);
    results.(sn).block_id                = obs.block_id;
    results.(sn).trial_id                = obs.trial_id;
    results.(sn).stim_id                 = obs.stim_id;
    results.(sn).is_det                  = obs.is_det;

    fprintf('  H=%.3f  β=%.2f  NLL=%.1f  BIC=%.1f\n', x_fit(1), x_fit(2), nll_fit, results.(sn).bic);
end

% Grand-average surprise aligned to reversal (diagnostic figure)
subj_ids2 = fieldnames(results);
half_win2 = 20;  x_ali2 = -half_win2:half_win2;
all_surp = [];
for si2 = 1:numel(subj_ids2)
    sn2 = subj_ids2{si2};
    r2  = results.(sn2);
    if ~isfield(all_trial_data,sn2), continue; end
    td2 = all_trial_data.(sn2).trial_data;
    if ~isfield(td2,'revTrial'), continue; end
    [nB2,~] = get_n_blocks_trials(td2);
    surp_row = NaN(1,numel(x_ali2));
    for b2 = 1:nB2
        rev2 = td2.revTrial(b2);
        if isnan(rev2), continue; end
        bm2 = r2.block_id==b2; b_tr2=r2.trial_id(bm2); b_sp2=r2.surprise(bm2);
        for xi2 = 1:numel(x_ali2)
            t_abs2 = round(rev2)+x_ali2(xi2);
            m2 = find(b_tr2==t_abs2,1);
            if ~isempty(m2), surp_row(xi2)=b_sp2(m2); end
        end
    end
    all_surp(end+1,:) = surp_row;
end

fig_s = figure('Position',[50 50 700 380]);
hold on;
title('Grand-average surprise ω×|δ| aligned to reversal');
if ~isempty(all_surp)
    mn_s = mean(all_surp,1,'omitnan');
    se_s = std(all_surp,0,1,'omitnan') / sqrt(sum(~all(isnan(all_surp),2)));
    fill([x_ali2,fliplr(x_ali2)],[mn_s+se_s,fliplr(mn_s-se_s)], ...
        [0.15 0.45 0.70],'FaceAlpha',0.2,'EdgeColor','none');
    plot(x_ali2, mn_s,'Color',[0.15 0.45 0.70],'LineWidth',2);
end
xline(0,'k--','LineWidth',1.5,'HandleVisibility','off');
xlabel('Trial relative to reversal'); ylabel('\omega_t × |\delta_t|');
saveas(fig_s, fullfile(outdir,'fit_fig_surprise_aligned.pdf'));
end

function results = fit_subjects_by_blocktype(all_trial_data, stim_features, results, outdir)
% Fit H and β separately for D-blocks and P-blocks per subject.

subjects  = fieldnames(results);
H_grid    = logspace(log10(0.005), log10(0.40), 30);
beta_grid = linspace(1, 20, 30);

for si = 1:numel(subjects)
    sn = subjects{si};
    if ~isfield(all_trial_data, sn), continue; end
    td = all_trial_data.(sn).trial_data;
    if isfield(td,'stimType'), sf = 'stimType';
    elseif isfield(td,'stimID'), sf = 'stimID';
    else, continue; end
    has_perc = isfield(td,'perceivedCorrect') && ~all(isnan(td.perceivedCorrect(:)));
    obs_all  = pack_observations(td, sf, has_perc);
    if isempty(obs_all.stim_id), continue; end

    [nB, ~] = get_n_blocks_trials(td);

    % Block type and context variables
    block_is_det     = false(1, nB);
    block_transition = cell(1, nB);
    for b = 1:nB
        if isfield(td,'trueFB')
            pfb_b = td.trueFB(b, ~isnan(td.trueFB(b,:)));
            block_is_det(b) = ~isempty(pfb_b) && mean(pfb_b) >= 0.99;
        else
            block_is_det(b) = true;
        end
    end
    for b = 1:nB
        curr = ternary(block_is_det(b),'D','P');
        if b == 1
            block_transition{b} = 'first';
        else
            prev = ternary(block_is_det(b-1),'D','P');
            block_transition{b} = [prev '→' curr];
        end
    end
    results.(sn).block_is_det    = block_is_det;
    results.(sn).block_transition = block_transition;

    % Fit per block type
    for bt = 1:2
        is_det_fit = (bt == 1);
        type_lbl   = ternary(is_det_fit, 'det', 'prob');
        obs_filt   = filter_obs_by_blocktype(obs_all, is_det_fit);
        if numel(obs_filt.stim_id) < 20
            results.(sn).(['H_fit_'    type_lbl]) = NaN;
            results.(sn).(['beta_fit_' type_lbl]) = NaN;
            results.(sn).(['nll_'      type_lbl]) = NaN;
            continue;
        end

        nll_grid = nan(numel(H_grid), numel(beta_grid));
        for hi = 1:numel(H_grid)
            for bi = 1:numel(beta_grid)
                nll_grid(hi,bi) = nassar_nll([H_grid(hi), beta_grid(bi)], obs_filt, stim_features);
            end
        end
        [~,bl] = min(nll_grid(:)); [hi0,bi0] = ind2sub(size(nll_grid),bl);
        x0 = [H_grid(hi0), beta_grid(bi0)];
        opts = optimoptions('fmincon','Display','off','MaxIterations',500);
        try
            xf = fmincon(@(x) nassar_nll(x, obs_filt, stim_features), x0,[],[],[],[], ...
                [0.005,0.5],[0.40,25],[],opts);
        catch
            xf = x0;
        end

        results.(sn).(['H_fit_'    type_lbl]) = xf(1);
        results.(sn).(['beta_fit_' type_lbl]) = xf(2);
        results.(sn).(['nll_'      type_lbl]) = nassar_nll(xf, obs_filt, stim_features);
        fprintf('  %s [%s]: H=%.3f  β=%.2f\n', sn, upper(type_lbl), xf(1), xf(2));
    end

    results.(sn).H_noise_sensitivity = ...
        results.(sn).H_fit_prob - results.(sn).H_fit_det;
    fprintf('  %s: ΔH=%.3f\n', sn, results.(sn).H_noise_sensitivity);
end

% Block-type comparison figure
subj_ids3 = fieldnames(results);
N_bt = numel(subj_ids3);
H_det_bt = nan(N_bt,1); H_prob_bt = nan(N_bt,1); H_ns_bt = nan(N_bt,1);

for si = 1:N_bt
    sn = subj_ids3{si};
    r  = results.(sn);
    if isfield(r,'H_fit_det'),         H_det_bt(si)  = r.H_fit_det;           end
    if isfield(r,'H_fit_prob'),        H_prob_bt(si) = r.H_fit_prob;          end
    if isfield(r,'H_noise_sensitivity'), H_ns_bt(si) = r.H_noise_sensitivity; end
end

fig_bt = figure('Position',[50 50 1000 380]);
sgtitle('Block-type split: H for D vs P blocks','FontSize',11);
clr_d = [0.15 0.45 0.70]; clr_p = [0.80 0.30 0.10];

ax = subplot(1,3,1); hold(ax,'on');
title(ax,'H_{det} vs H_{prob}');
ok_bt = ~isnan(H_det_bt) & ~isnan(H_prob_bt);
scatter(ax, H_det_bt(ok_bt), H_prob_bt(ok_bt), 60, [0.4 0.4 0.4],'filled','MarkerFaceAlpha',0.8);
lims_bt = [0, max([H_det_bt(ok_bt); H_prob_bt(ok_bt)])+0.02];
plot(ax,lims_bt,lims_bt,'k--','LineWidth',1,'HandleVisibility','off');
[~,p_bt] = ttest(H_prob_bt(ok_bt)-H_det_bt(ok_bt));
text(ax,0.05,0.97,sprintf('ΔH t-test: p=%.3f',p_bt), ...
    'Units','normalized','VerticalAlignment','top','FontSize',8,'BackgroundColor','w');
xlabel(ax,'H_{det}'); ylabel(ax,'H_{prob}');
xlim(ax,lims_bt); ylim(ax,lims_bt); axis(ax,'square');

ax = subplot(1,3,2); hold(ax,'on');
title(ax,'Noise sensitivity ΔH = H_{prob} − H_{det}');
ok_ns2 = ~isnan(H_ns_bt);
histogram(ax,H_ns_bt(ok_ns2),10,'FaceColor',[0.4 0.4 0.4],'EdgeColor','w','FaceAlpha',0.8);
xline(ax,0,'k--','LineWidth',1.5,'HandleVisibility','off');
[~,p_ns2] = ttest(H_ns_bt(ok_ns2));
text(ax,0.98,0.97,sprintf('Mean=%.3f\nt-test p=%.3f', ...
    mean(H_ns_bt(ok_ns2),'omitnan'),p_ns2), ...
    'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
    'FontSize',8,'BackgroundColor','w');
xlabel(ax,'\Delta H'); ylabel(ax,'Subjects');

ax = subplot(1,3,3); hold(ax,'on');
title(ax,'H_{combined} by block transition');
trans_types  = {'first','D→D','D→P','P→D','P→P'};
trans_colors = {[0.5 0.5 0.5], clr_d, clr_p, clr_p*0.7+clr_d*0.3, clr_p};
all_trans_H  = cell(1,5);
for si = 1:N_bt
    sn = subj_ids3{si};
    r  = results.(sn);
    if ~isfield(r,'block_transition') || ~isfield(r,'H_fit'), continue; end
    for b = 1:numel(r.block_transition)
        ti = find(strcmp(trans_types, r.block_transition{b}));
        if ~isempty(ti), all_trans_H{ti}(end+1) = r.H_fit; end
    end
end
for ti = 1:5
    if ~isempty(all_trans_H{ti})
        boxchart(ax, ti*ones(size(all_trans_H{ti})), all_trans_H{ti}, ...
            'BoxFaceColor',trans_colors{ti},'MarkerColor',trans_colors{ti});
    end
end
set(ax,'XTick',1:5,'XTickLabel',trans_types,'XTickLabelRotation',25);
ylabel(ax,'H_{combined}');
saveas(fig_bt, fullfile(outdir,'fig_blocktype_H.pdf'));
end

function obs_filt = filter_obs_by_blocktype(obs, keep_det)
mask = obs.is_det == keep_det;
fn = fieldnames(obs);
for fi = 1:numel(fn)
    obs_filt.(fn{fi}) = obs.(fn{fi})(mask);
end
end

function td_sim = simulate_one_subject(H, beta, nB, nT, rev_trials, ...
    N_STIM, N_REPS, P_TRUE_FB, PRE_REV_RULE, POST_REV_RULE, SWITCH_STIMS, ALPHA_INIT)
td_sim.correct    = nan(nB, nT);
td_sim.goTrial    = nan(nB, nT);
td_sim.confidence = nan(nB, nT);
td_sim.stimID     = nan(nB, nT);
td_sim.switched   = nan(nB, nT);
td_sim.revTrial   = rev_trials(:)';

for b = 1:nB
    rev_b = round(rev_trials(b));
    pfb   = ternary(b <= numel(P_TRUE_FB), P_TRUE_FB(b), 1.0);

    stim_order = repmat(1:N_STIM,1,N_REPS);
    stim_order = stim_order(randperm(nT));
    go_nogo    = nan(1,nT);
    for t = 1:nT
        s = stim_order(t);
        go_nogo(t) = double(ternary(t<=rev_b, PRE_REV_RULE(s), POST_REV_RULE(s)));
    end
    nTrue  = round(pfb * nT);
    fb_vec = double([true(1,nTrue), false(1,nT-nTrue)]);
    fb_vec = fb_vec(randperm(nT));

    theta     = 0.5 * ones(1, N_STIM);
    n_eff     = ones(1, N_STIM);
    first_obs = true(1, N_STIM);

    for t = 1:nT
        s_t      = stim_order(t);
        go_true  = go_nogo(t);
        fb_true  = fb_vec(t);
        certainty_t = abs(theta(s_t) - 0.5);
        p_go     = 1 / (1 + exp(-beta * (theta(s_t) - 0.5)));
        resp_go  = double(rand < p_go);
        correct  = double(resp_go == go_true);
        perc_cor = fb_true * correct + (1-fb_true) * (1-correct);
        y_t = double((resp_go==1 && perc_cor==1)||(resp_go==0 && perc_cor==0));
        pred_prob = max(theta(s_t)^y_t * (1-theta(s_t))^(1-y_t), 1e-6);
        chi_t    = 0.5 / pred_prob;
        omega_t  = (H * chi_t) / (H * chi_t + (1-H));
        n_eff(s_t) = (1-omega_t)*(n_eff(s_t)+1) + omega_t;
        alpha_t  = omega_t + (1-omega_t)/n_eff(s_t);
        if first_obs(s_t), alpha_t = max(alpha_t,ALPHA_INIT); first_obs(s_t)=false; end
        delta_t    = y_t - theta(s_t);
        theta(s_t) = max(0.01, min(0.99, theta(s_t) + alpha_t * delta_t));
        conf_raw   = 1 + 18*certainty_t + 0.8*randn;
        td_sim.correct(b,t)    = correct;
        td_sim.goTrial(b,t)    = go_true;
        td_sim.confidence(b,t) = max(1,min(10,round(conf_raw)));
        td_sim.stimID(b,t)     = s_t;
        td_sim.switched(b,t)   = ismember(s_t, SWITCH_STIMS);
    end
end
end

function stats = compute_summary_stats(td, nB, nT, SWITCH_STIMS)
stats = nan(1,6);
pre_acc = []; post_acc = []; rec_slope = [];
conf_vec = []; acc_vec = [];
sw_acc = []; mn_acc = [];
for b = 1:nB
    rev = td.revTrial(b);
    if isnan(rev) || rev < 5 || rev > nT-5, continue; end
    rev = round(rev);
    acc_row = td.correct(b,:);
    pre_acc(end+1)  = mean(acc_row(1:rev),'omitnan');
    post_acc(end+1) = mean(acc_row(rev+1:min(rev+10,nT)),'omitnan');
    t_rec = (rev+11):min(rev+30,nT);
    if numel(t_rec) > 3
        p_rec = polyfit(t_rec, acc_row(t_rec), 1);
        rec_slope(end+1) = p_rec(1);
    end
    if isfield(td,'confidence')
        c_row = td.confidence(b,:);
        ok_ca = ~isnan(acc_row) & ~isnan(c_row);
        if sum(ok_ca) > 5
            conf_vec = [conf_vec, c_row(ok_ca)];
            acc_vec  = [acc_vec,  acc_row(ok_ca)];
        end
    end
    post_mask = (1:nT) > rev;
    if isfield(td,'switched')
        sw_trials = post_mask & td.switched(b,:)==1;
        mn_trials = post_mask & td.switched(b,:)==0;
    elseif isfield(td,'stimID')
        sw_trials = post_mask & ismember(td.stimID(b,:), SWITCH_STIMS);
        mn_trials = post_mask & ~ismember(td.stimID(b,:), SWITCH_STIMS);
    else
        sw_trials = false(1,nT); mn_trials = false(1,nT);
    end
    if any(sw_trials), sw_acc(end+1) = mean(acc_row(sw_trials),'omitnan'); end
    if any(mn_trials), mn_acc(end+1) = mean(acc_row(mn_trials),'omitnan'); end
end
if ~isempty(pre_acc),  stats(1) = mean(pre_acc,'omitnan');  end
if ~isempty(post_acc), stats(2) = mean(post_acc,'omitnan'); end
if ~isempty(pre_acc) && ~isempty(post_acc)
    stats(3) = mean(pre_acc,'omitnan') - mean(post_acc,'omitnan');
end
if ~isempty(rec_slope), stats(4) = mean(rec_slope,'omitnan'); end
if numel(conf_vec) > 10
    stats(5) = corr(conf_vec(:),acc_vec(:),'Type','Spearman','Rows','complete');
end
if ~isempty(sw_acc) && ~isempty(mn_acc)
    stats(6) = mean(mn_acc,'omitnan') - mean(sw_acc,'omitnan');
end
end

function [acc_pre, acc_post, rev_cost, conf_mean] = extract_beh_summary(td)
acc_pre=NaN; acc_post=NaN; rev_cost=NaN; conf_mean=NaN;
if ~isfield(td,'correct') || ~isfield(td,'revTrial'), return; end
[nB,nT] = get_n_blocks_trials(td);
all_pre = []; all_post = [];
for b = 1:nB
    rev = td.revTrial(b); if isnan(rev), continue; end
    all_pre  = [all_pre,  td.correct(b, 1:round(rev))];
    all_post = [all_post, td.correct(b, round(rev)+1:nT)];
end
acc_pre=mean(all_pre,'omitnan'); acc_post=mean(all_post,'omitnan');
rev_cost=acc_pre-acc_post;
if isfield(td,'confidence'), conf_mean=mean(td.confidence(:),'omitnan'); end
end

function [nBlocks, nTrials] = get_n_blocks_trials(td)
if isfield(td,'correct') && ~isempty(td.correct)
    [nBlocks, nTrials] = size(td.correct);
elseif isfield(td,'goTrial') && ~isempty(td.goTrial)
    [nBlocks, nTrials] = size(td.goTrial);
elseif isfield(td,'stimID') && ~isempty(td.stimID)
    [nBlocks, nTrials] = size(td.stimID);
else
    error('Cannot infer block/trial dimensions from trial_data.');
end
end

function plot_ribbon(ax, x, mat, clr, ls, lbl)
if isempty(mat) || all(isnan(mat(:))), return; end
mn = movmean(mean(mat,1,'omitnan'), 5,'omitnan');
se = std(mat,0,1,'omitnan') ./ sqrt(max(sum(~isnan(mat),1),1));
fill(ax,[x,fliplr(x)],[mn+se,fliplr(mn-se)],clr, ...
    'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
plot(ax,x,mn,'Color',clr,'LineWidth',2,'LineStyle',ls, ...
    'DisplayName',sprintf('%s (n=%d)',lbl,size(mat,1)));
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function group_T = add_nassar_to_group_table(group_T, results)
% Joins Nassar model latents into group_T using (subjID, block, trial) as key.
% Nassar results are indexed per obs-row (packed, no NaN trials).
% group_T has one row per behavioural trial including NaN-epoch rows.

cols = {'PE_nassar','omega','alpha_nassar','certainty','surprise','theta_nassar'};
for c = cols
    group_T.(c{1}) = nan(height(group_T), 1);
end

subj_col = 'subjID';
if ~ismember(subj_col, group_T.Properties.VariableNames) && ismember('subj_id', group_T.Properties.VariableNames)
    subj_col = 'subj_id';
end
subjs = categories(categorical(group_T.(subj_col)));

for si = 1:numel(subjs)
    sn = subjs{si};
    if ~isfield(results, sn), continue; end
    r = results.(sn);

    sm = string(group_T.(subj_col)) == string(sn);
    sub_rows = find(sm);

    % Build lookup: (block, trial) → row index in group_T
    blocks_gt = double(group_T.block(sub_rows));
    trials_gt = double(group_T.trial(sub_rows));

    % Iterate over Nassar obs-level rows
    for t = 1:numel(r.trial_id)
        b_t = r.block_id(t);
        tr_t = r.trial_id(t);
        match = sub_rows(blocks_gt == b_t & trials_gt == tr_t);
        if isempty(match), continue; end
        group_T.PE_nassar(match)      = r.delta_trial(t);
        group_T.omega(match)          = r.omega_trial(t);
        group_T.alpha_nassar(match)   = r.alpha_trial(t);
        group_T.certainty(match)      = r.certainty_trial(t);
        group_T.surprise(match)       = r.surprise(t);
        group_T.theta_nassar(match)   = r.theta_trial(t);
    end
end
end


function switch_stims = read_switch_stims_from_table_lc(group_T)
%READ_SWITCH_STIMS_FROM_TABLE_LC  Robustly read first non-empty switch_stims value.
switch_stims = [];
if isempty(group_T) || ~istable(group_T) || ~ismember('switch_stims', group_T.Properties.VariableNames)
    return;
end
for ri = 1:height(group_T)
    try
        val = group_T.switch_stims(ri);
        if iscell(group_T.switch_stims), val = group_T.switch_stims{ri}; end
        if isempty(val), continue; end
        if isstring(val) || ischar(val)
            nums = regexp(char(val), '\d+', 'match');
            if ~isempty(nums), switch_stims = cellfun(@str2double, nums); end
        elseif isnumeric(val)
            switch_stims = double(val(:)');
            switch_stims = switch_stims(isfinite(switch_stims));
        end
        if ~isempty(switch_stims), return; end
    catch
    end
end
end

function block_types = get_block_types_from_td_lc(td)
    [nB, ~] = size(td.correct);

    if isfield(td,'block_structure') && ~isempty(td.block_structure)
        bs = upper(char(td.block_structure));
        bs(bs == 'V') = 'P';

        if numel(bs) >= nB
            block_types = cellstr(bs(1:nB)');
            return;
        end
    end

    if isfield(td,'trueFB')
        block_types = repmat({'D'}, nB, 1);
        for b = 1:nB
            fb = td.trueFB(b,:);
            fb = fb(~isnan(fb));
            if ~isempty(fb) && mean(fb) < 0.99
                block_types{b} = 'P';
            end
        end
        return;
    end

    block_types = repmat({'unknown'}, nB, 1);
end
