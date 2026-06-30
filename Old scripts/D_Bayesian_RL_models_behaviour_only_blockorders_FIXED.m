% ==========================================================================
% D_BAYESIAN_RL_MODELS_v3.m
% Nassar (2010) adaptive learning rate model — simulation and fitting
% Category Switch Braille Go/NoGo task
%
% KEY CHANGES FROM v2 (behaviour-only FIXED)
% ──────────────────────────────────────────
%  1. EMPIRICALLY GROUNDED SIMULATION
%     The simulation now reproduces the exact structure of the real dataset
%     (all_trial_data_June2026.mat), extracted programmatically:
%
%     Block orders (42 subjects, Ox02 excluded as incomplete):
%       DPDP (n=11), PDPD (n=6), PPDDP (n=5), PDPDP (n=4),
%       PPDPD (n=3), DDPPD (n=3), DDDD (n=2), PDDPP (n=2),
%       DDPP (n=1), PPPP (n=1), DDPD (n=1), DPPPD (n=1),
%       PPPDD (n=1), PDDPD (n=1)
%
%     N blocks per subject: 4 (n=22) or 5 (n=20)
%
%     Switch stims (which 2 of the 4 stims swap Go/NoGo assignment):
%       KH cohort (Ox subjects): [1,4] in 110/128 blocks; [2,3] in 18/128
%         → alternates by block per individual subject schedule
%       RR cohort (Nc subjects): always [2,3]
%
%     Reversal trial (sampled from truncated Normal matching empirical):
%       D blocks: mean=52.5, SD=11.8, range [31,72]
%       P blocks: mean=48.6, SD=12.5, range [30,70]
%       These distributions reflect earlier reversals in P blocks
%       (consistent with Behrens et al. 2007: volatility increases
%       premature updating).
%
%  2. V BLOCKS TREATED AS P BLOCKS
%     'V' (visual probabilistic) blocks have identical feedback
%     fidelity to P blocks (p_trueFB = 0.8) and are treated as P
%     throughout. The distinction is stimulus modality only, not
%     uncertainty structure. This follows your clarification.
%
%  3. switch_stims IS INDEPENDENT OF stim_config
%     The 4 stimuli within a given config can have any two designated as
%     switched — this is determined per block from the empirical data,
%     not inferred from stim_config. Maintained stimuli are the
%     remaining 2.
%
%  4. n_prev_P TRACKED AND STORED PER BLOCK
%     Running count of P (and V) blocks experienced before each block.
%     Stored in sim_data and exported to group_T for LME covariate use
%     (Behrens et al. 2007: accumulated uncertainty experience modulates
%     learning rate priors).
%
% THEORETICAL GROUNDING
% ─────────────────────
%  Dual timescale uncertainty framework (Yu & Dayan 2005, Behrens et al. 2007):
%    Short-term uncertainty: probabilistic feedback (P blocks, p_trueFB=0.8)
%      → inflates apparent PE variance → inflates fitted H
%      → captured by H_noise_sensitivity = H_prob − H_det
%    Long-term uncertainty: reversals (embedded in all blocks)
%      → requires change-point detection via ω_t
%      → modulated by n_prev_P (accumulated P experience)
%
%  Nassar (2010) adaptive delta-rule:
%    θ_{t+1} = θ_t + α_t × δ_t
%    α_t = ω_t + (1−ω_t)/n_eff_t    [learning rate]
%    ω_t = H×χ_t / (H×χ_t + 1−H)   [change-point probability]
%    χ_t = p(y_t | prior) / 0.5      [surprise ratio]
%
%  Confidence (prospective, pre-outcome):
%    certainty_t = |θ_t − 0.5|  [Boldt & Yeung 2015]
%    Rated before feedback, after response. High certainty + neg feedback
%    = maximal FRN (Holroyd & Coles 2002 RPE signal).
%
% ==========================================================================

clear; close all; rng(42);

% =========================================================================
%% §1  PATHS
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
if ~exist(outpath, 'dir'), mkdir(outpath); end

fprintf('\n========================================================\n');
fprintf('D_BAYESIAN_RL_MODELS_v3 — empirically grounded simulation\n');
fprintf('========================================================\n');

% -------------------------------------------------------------------------
% Load real data
% -------------------------------------------------------------------------
load(fullfile(data_path, 'all_trial_data_June2026.mat'), 'all_trial_data');

group_T = [];
beh_candidates = {fullfile(data_path, 'behav_table_June2026.mat')};
for ci = 1:numel(beh_candidates)
    if exist(beh_candidates{ci}, 'file')
        tmp = load(beh_candidates{ci});
        fns = fieldnames(tmp);
        for fi = 1:numel(fns)
            if istable(tmp.(fns{fi}))
                group_T = tmp.(fns{fi});
                fprintf('Loaded group_T from %s (%d rows)\n', beh_candidates{ci}, height(group_T));
                break;
            end
        end
        if ~isempty(group_T), break; end
    end
end

% =========================================================================
%% §2  EMPIRICAL SUBJECT REGISTRY
%
% Built from all_trial_data_June2026.mat (see header).
% Each row = one subject. Block types inferred from mean(trueFB) per block:
%   mean > 0.95 → D (deterministic)
%   mean ≤ 0.95 → P (probabilistic, includes V blocks)
% Switch stims taken directly from td.switch_stims.
% Reversal trials taken directly from td.revTrial.
%
% Ox02 excluded (only 2 blocks, pilot/incomplete).
% =========================================================================

EMPIRICAL_SUBJECTS = build_empirical_registry(all_trial_data);
fprintf('Loaded %d subjects from empirical registry.\n', numel(EMPIRICAL_SUBJECTS));

% Print summary
n4 = sum(cellfun(@(s) s.nblks==4, EMPIRICAL_SUBJECTS));
n5 = sum(cellfun(@(s) s.nblks==5, EMPIRICAL_SUBJECTS));
fprintf('  4-block subjects: %d | 5-block subjects: %d\n', n4, n5);

% =========================================================================
%% §3  MODEL PARAMETERS
% =========================================================================
N_TRIALS       = 100;
N_STIM         = 4;      % 4 stimuli per block (within a given stim_config)
P_TRUE_FB_D    = 1.00;   % deterministic: feedback always valid
P_TRUE_FB_P    = 0.80;   % probabilistic: 80% valid (Behrens et al. 2007)
% V blocks use P_TRUE_FB_P (visual modality, same uncertainty structure)

SOFTMAX_BETA   = 6;
CONF_SCALE_MAX = 10;
CONF_OFFSET    = 1;
CONF_GAIN      = (CONF_SCALE_MAX - CONF_OFFSET) / 0.5;  % = 18 (scales |θ−0.5| to 1–10)
ALPHA_INIT     = 0.50;   % first-encounter learning rate (Nassar 2010)

% Reversal trial sampling — from empirical distributions
% D blocks: mean=52.5, SD=11.8, range [31,72]  (truncated Normal)
% P blocks: mean=48.6, SD=12.5, range [30,70]
REV_D_MU  = 52.5;  REV_D_SIG  = 11.8;  REV_D_LO  = 31;  REV_D_HI  = 72;
REV_P_MU  = 48.6;  REV_P_SIG  = 12.5;  REV_P_LO  = 30;  REV_P_HI  = 70;

% Stimulus feature structure — 2×2 dimensional shift
% Pre-reversal rule:  Go iff Feature1 == 1  → stims 1,2 are Go
% Post-reversal rule: Go iff Feature2 == 1  → stims 1,3 are Go
%
%  Stim | F1 | F2 | Pre-rev  | Post-rev | Category
%  ─────┼────┼────┼──────────┼──────────┼─────────────
%    1  |  1 |  1 |  Go      |  Go      | MAINTAINED
%    2  |  1 |  2 |  Go      |  NoGo    | SWITCHED  (pair [1,4] → stim 2 is one of pair)
%    3  |  2 |  1 |  NoGo    |  Go      | SWITCHED  (pair [2,3] → stim 3 is one of pair)
%    4  |  2 |  2 |  NoGo    |  NoGo    | MAINTAINED
%
% NOTE: which two stims are designated "switched" varies per block/subject.
% Empirically: KH subjects use [1,4] (n=110 blocks) or [2,3] (n=18 blocks);
% RR subjects always use [2,3] (n=60 blocks).
% The switch_stims pair defines the stims that CHANGE their Go/NoGo rule.
% Maintained stims have the same correct response before and after reversal.
STIM_FEATURES = [1 1; 1 2; 2 1; 2 2];  % [F1, F2] for stims 1–4

params.H    = 0.05;
params.beta = SOFTMAX_BETA;

% =========================================================================
%% §4  SIMULATION LOOP — structured to match empirical dataset
%
% Strategy: for each real subject in the registry, simulate ONE synthetic
% dataset with IDENTICAL block structure (order, n_blocks, switch_stims)
% and EMPIRICALLY SAMPLED reversal trials. This means the simulation
% population mirrors the real population in design variability.
% =========================================================================
fprintf('\nRunning simulation (1 synthetic dataset per real subject structure)...\n');

N_SIM_SUBJECTS = numel(EMPIRICAL_SUBJECTS);
sim_data   = struct();
model_vars = struct();

for si = 1:N_SIM_SUBJECTS
    es   = EMPIRICAL_SUBJECTS{si};
    sl   = sprintf('sim_%s', es.sn);  % e.g. sim_Ox09
    nB   = es.nblks;

    td   = init_trial_data(nB, N_TRIALS);
    mv   = init_model_vars(nB, N_TRIALS, N_STIM);

    td.block_structure = es.block_str;  % e.g. 'PDPDP'
    td.block_type_label = es.block_types(:);
    td.cohort = es.cohort;
    td.n_prev_P = zeros(nB,1);

    running_P = 0;
    for b = 1:nB
        % --- Block metadata ---
        curr_type = es.block_types{b};   % 'D' or 'P'
        sw_stims  = es.switch_stims{b};  % e.g. [1 4] or [2 3]
        maint_stims = setdiff(1:N_STIM, sw_stims);

        td.n_prev_P(b) = running_P;

        % Feedback fidelity: D=1.0, P=0.8 (V already recoded as P)
        if strcmp(curr_type, 'D')
            p_trueFB = P_TRUE_FB_D;
        else
            p_trueFB = P_TRUE_FB_P;
        end
        td.trueFB_block(b) = p_trueFB;

        % --- Sample reversal trial from truncated Normal ---
        % Matches the empirical reversal-time distributions observed in
        % real data (D blocks: mu=52.5, P blocks: mu=48.6).
        if strcmp(curr_type, 'D')
            revTrial = sample_truncated_normal(REV_D_MU, REV_D_SIG, REV_D_LO, REV_D_HI);
        else
            revTrial = sample_truncated_normal(REV_P_MU, REV_P_SIG, REV_P_LO, REV_P_HI);
        end

        % --- Stimulus order (balanced, permuted) ---
        stim_order = repmat(1:N_STIM, 1, N_TRIALS/N_STIM);
        stim_order = stim_order(randperm(N_TRIALS));

        % --- Go/NoGo schedule ---
        % Pre-reversal: stims 1,2 are Go (Feature1==1)
        % Post-reversal: stims 1,3 are Go (Feature2==1)
        % BUT: which stims are "switched" is determined by sw_stims.
        % If sw_stims = [1,4]: stims 1 and 4 swap identity at reversal
        % If sw_stims = [2,3]: stims 2 and 3 swap identity at reversal
        % The pre/post rule is constructed so that exactly sw_stims change.
        go_nogo = compute_go_nogo_schedule(stim_order, revTrial, sw_stims, N_TRIALS);

        % --- Feedback validity vector ---
        n_true_fb = round(p_trueFB * N_TRIALS);
        fb_arr = [true(1,n_true_fb), false(1,N_TRIALS-n_true_fb)];
        trueFB_vec = double(fb_arr(randperm(N_TRIALS)));

        % --- Nassar belief update ---
        theta     = 0.5 * ones(1, N_STIM);
        n_eff     = ones(1, N_STIM);
        first_obs = true(1, N_STIM);

        for t = 1:N_TRIALS
            s_t      = stim_order(t);
            go_true  = go_nogo(t);
            fb_true  = trueFB_vec(t);

            % PROSPECTIVE certainty: computed BEFORE response/feedback
            % (trial order: stimulus → response → confidence → feedback)
            % Matches Boldt & Yeung (2015) confidence-as-precision
            certainty_t = abs(theta(s_t) - 0.5);

            % Softmax response
            p_go    = 1 / (1 + exp(-params.beta * (theta(s_t) - 0.5)));
            resp_go = double(rand < p_go);
            correct = double(resp_go == go_true);

            % Perceived correctness (what participants experience)
            perc_cor = fb_true * correct + (1-fb_true) * (1-correct);

            % y_t = outcome that drives belief update (perceived feedback)
            y_t = double((resp_go==1 && perc_cor==1) || (resp_go==0 && perc_cor==0));

            % Nassar (2010) update equations
            pred_prob  = max(theta(s_t)^y_t * (1-theta(s_t))^(1-y_t), 1e-6);
            chi_t      = 0.5 / pred_prob;
            omega_t    = (params.H * chi_t) / (params.H * chi_t + (1-params.H));
            n_eff(s_t) = (1-omega_t) * (n_eff(s_t)+1) + omega_t;
            alpha_t    = omega_t + (1-omega_t) / n_eff(s_t);

            % First-encounter boost (Nassar 2010 §2.4)
            if first_obs(s_t)
                alpha_t        = max(alpha_t, ALPHA_INIT);
                first_obs(s_t) = false;
            end

            delta_t    = y_t - theta(s_t);
            theta(s_t) = max(0.01, min(0.99, theta(s_t) + alpha_t * delta_t));

            % Confidence rating (1–10, pre-outcome certainty + noise)
            % Boldt & Yeung (2015): confidence reflects prospective belief precision
            conf_raw   = CONF_OFFSET + CONF_GAIN * certainty_t + 0.8*randn;
            confidence = max(1, min(10, round(conf_raw)));

            % Store trial-level data
            td.stimID(b,t)           = s_t;
            td.switch_stims{b}       = sw_stims;
            td.is_switch_stim(b,t)   = ismember(s_t, sw_stims);
            td.goTrial(b,t)          = go_true;
            td.respWasGo(b,t)        = resp_go;
            td.correct(b,t)          = correct;
            td.perceivedCorrect(b,t) = perc_cor;
            td.trueFB(b,t)           = fb_true;
            td.p_trueFB(b,t)         = p_trueFB;
            td.confidence(b,t)       = confidence;
            td.blocknum(b,t)         = b;
            td.trialnum(b,t)         = t;

            % Store model variables
            mv.theta(b,t,:)     = theta;
            mv.alpha(b,t)       = alpha_t;
            mv.omega(b,t)       = omega_t;
            mv.delta(b,t)       = delta_t;
            mv.certainty(b,t)   = certainty_t;
            mv.n_eff(b,t)       = n_eff(s_t);
            mv.surprise(b,t)    = omega_t * abs(delta_t);
        end

        td.revTrial(b)    = revTrial;
        td.block_type(b)  = string(curr_type);

        if strcmp(curr_type,'P')
            running_P = running_P + 1;
        end
    end

    sim_data.(sl)   = td;
    model_vars.(sl) = mv;
end

fprintf('Simulation complete: %d subjects, structures matching real data.\n', N_SIM_SUBJECTS);
save(fullfile(outpath,'sim_data_v3.mat'), 'sim_data', 'model_vars', 'params', 'EMPIRICAL_SUBJECTS');

% =========================================================================
%% §5  SIMULATION VALIDATION FIGURES
%
% These figures demonstrate that the simulation reproduces key qualitative
% features of the real data before any fitting. This is prerequisite step 2
% of the parameter recovery framework (simulation sanity checks).
% =========================================================================

clr_D   = [0.15 0.45 0.70];  % blue — deterministic
clr_P   = [0.80 0.30 0.10];  % orange — probabilistic
clr_sw  = [0.75 0.10 0.10];  % red — switched stimuli
clr_mn  = [0.15 0.45 0.70];  % blue — maintained stimuli
CLR_STGS = [0.12 0.62 0.47; 0.85 0.65 0.00; 0.80 0.27 0.13; 0.40 0.25 0.65];

% ── Fig S1: Accuracy reversal-aligned — switch vs maintained, D vs P ──────
fprintf('Generating simulation figures...\n');
half_win = 30; x_ali = -half_win:(half_win-1); n_x = numel(x_ali);

sw_acc_D  = NaN(0,n_x); sw_acc_P  = NaN(0,n_x);
mn_acc_D  = NaN(0,n_x); mn_acc_P  = NaN(0,n_x);
sw_conf_D = NaN(0,n_x); sw_conf_P = NaN(0,n_x);
mn_conf_D = NaN(0,n_x); mn_conf_P = NaN(0,n_x);

snames = fieldnames(sim_data);
for si = 1:numel(snames)
    td  = sim_data.(snames{si});
    mv2 = model_vars.(snames{si});
    nB  = numel(td.revTrial);

    for b = 1:nB
        bt  = char(td.block_type(b));
        rev = td.revTrial(b);
        if isnan(rev), continue; end

        sw_b = td.switch_stims{b};
        row_sw_acc = NaN(1,n_x); row_mn_acc = NaN(1,n_x);
        row_sw_cf  = NaN(1,n_x); row_mn_cf  = NaN(1,n_x);

        for xi = 1:n_x
            t_abs = round(rev) + x_ali(xi);
            if t_abs < 1 || t_abs > 100, continue; end
            s_t = td.stimID(b,t_abs);
            if isnan(s_t), continue; end

            acc_val  = td.correct(b,t_abs);
            conf_val = td.confidence(b,t_abs);

            if ismember(s_t, sw_b)
                row_sw_acc(xi) = acc_val;
                row_sw_cf(xi)  = conf_val;
            else
                row_mn_acc(xi) = acc_val;
                row_mn_cf(xi)  = conf_val;
            end
        end

        if strcmp(bt,'D')
            if any(~isnan(row_sw_acc)), sw_acc_D(end+1,:) = row_sw_acc; end
            if any(~isnan(row_mn_acc)), mn_acc_D(end+1,:) = row_mn_acc; end
            if any(~isnan(row_sw_cf)),  sw_conf_D(end+1,:) = row_sw_cf; end
            if any(~isnan(row_mn_cf)),  mn_conf_D(end+1,:) = row_mn_cf; end
        else
            if any(~isnan(row_sw_acc)), sw_acc_P(end+1,:) = row_sw_acc; end
            if any(~isnan(row_mn_acc)), mn_acc_P(end+1,:) = row_mn_acc; end
            if any(~isnan(row_sw_cf)),  sw_conf_P(end+1,:) = row_sw_cf; end
            if any(~isnan(row_mn_cf)),  mn_conf_P(end+1,:) = row_mn_cf; end
        end
    end
end

fig1 = figure('Position',[50 50 1400 580]);
sgtitle({'Simulated accuracy and confidence — reversal-aligned', ...
    'Switched vs maintained stimuli | D-blocks (blue) vs P-blocks (orange)', ...
    sprintf('N=%d simulated subjects matching empirical block structures', N_SIM_SUBJECTS)}, ...
    'FontSize',11);

% Accuracy rows
ax = subplot(2,4,1); hold(ax,'on');
title(ax,'Accuracy — D blocks, switched','FontSize',9);
plot_ribbon_lc(ax, x_ali, sw_acc_D, clr_sw, '--', sprintf('Switch (n=%d)',size(sw_acc_D,1)));
xline(ax,0,'k--','LineWidth',1.5,'HandleVisibility','off');
yline(ax,0.5,'k:','HandleVisibility','off');
xlabel(ax,'Trial rel. reversal'); ylabel(ax,'P(correct)');
xlim(ax,[-half_win half_win]); ylim(ax,[0 1]);
legend(ax,'Box','off','FontSize',8);

ax = subplot(2,4,2); hold(ax,'on');
title(ax,'Accuracy — D blocks, maintained','FontSize',9);
plot_ribbon_lc(ax, x_ali, mn_acc_D, clr_mn, '-', sprintf('Maint. (n=%d)',size(mn_acc_D,1)));
xline(ax,0,'k--','LineWidth',1.5,'HandleVisibility','off');
yline(ax,0.5,'k:','HandleVisibility','off');
xlabel(ax,'Trial rel. reversal'); ylabel(ax,'P(correct)');
xlim(ax,[-half_win half_win]); ylim(ax,[0 1]);
legend(ax,'Box','off','FontSize',8);

ax = subplot(2,4,3); hold(ax,'on');
title(ax,'Accuracy — P blocks, switched','FontSize',9);
plot_ribbon_lc(ax, x_ali, sw_acc_P, clr_sw, '--', sprintf('Switch (n=%d)',size(sw_acc_P,1)));
xline(ax,0,'k--','LineWidth',1.5,'HandleVisibility','off');
yline(ax,0.5,'k:','HandleVisibility','off');
xlabel(ax,'Trial rel. reversal'); ylabel(ax,'P(correct)');
xlim(ax,[-half_win half_win]); ylim(ax,[0 1]);
legend(ax,'Box','off','FontSize',8);

ax = subplot(2,4,4); hold(ax,'on');
title(ax,'Accuracy — P blocks, maintained','FontSize',9);
plot_ribbon_lc(ax, x_ali, mn_acc_P, clr_mn, '-', sprintf('Maint. (n=%d)',size(mn_acc_P,1)));
xline(ax,0,'k--','LineWidth',1.5,'HandleVisibility','off');
yline(ax,0.5,'k:','HandleVisibility','off');
xlabel(ax,'Trial rel. reversal'); ylabel(ax,'P(correct)');
xlim(ax,[-half_win half_win]); ylim(ax,[0 1]);
legend(ax,'Box','off','FontSize',8);

% Confidence rows
ax = subplot(2,4,5); hold(ax,'on');
title(ax,'Confidence — D blocks, switched','FontSize',9);
plot_ribbon_lc(ax, x_ali, sw_conf_D, clr_sw, '--', 'Switch');
xline(ax,0,'k--','LineWidth',1.5,'HandleVisibility','off');
xlabel(ax,'Trial rel. reversal'); ylabel(ax,'Confidence (1–10)');
xlim(ax,[-half_win half_win]); ylim(ax,[1 10]);
annotation_text(ax, 'Prospective certainty|θ−0.5| (pre-outcome)', 8);

ax = subplot(2,4,6); hold(ax,'on');
title(ax,'Confidence — D blocks, maintained','FontSize',9);
plot_ribbon_lc(ax, x_ali, mn_conf_D, clr_mn, '-', 'Maintained');
xline(ax,0,'k--','LineWidth',1.5,'HandleVisibility','off');
xlabel(ax,'Trial rel. reversal'); ylabel(ax,'Confidence (1–10)');
xlim(ax,[-half_win half_win]); ylim(ax,[1 10]);

ax = subplot(2,4,7); hold(ax,'on');
title(ax,'Confidence — P blocks, switched','FontSize',9);
plot_ribbon_lc(ax, x_ali, sw_conf_P, clr_sw, '--', 'Switch');
xline(ax,0,'k--','LineWidth',1.5,'HandleVisibility','off');
xlabel(ax,'Trial rel. reversal'); ylabel(ax,'Confidence (1–10)');
xlim(ax,[-half_win half_win]); ylim(ax,[1 10]);

ax = subplot(2,4,8); hold(ax,'on');
title(ax,'Confidence — P blocks, maintained','FontSize',9);
plot_ribbon_lc(ax, x_ali, mn_conf_P, clr_mn, '-', 'Maintained');
xline(ax,0,'k--','LineWidth',1.5,'HandleVisibility','off');
xlabel(ax,'Trial rel. reversal'); ylabel(ax,'Confidence (1–10)');
xlim(ax,[-half_win half_win]); ylim(ax,[1 10]);

annotation('textbox',[0.01 0.01 0.98 0.05],'String', ...
    ['Switched stims: 2 of 4 stimuli that change Go/NoGo assignment at reversal (block-specific, from empirical data). ' ...
     'Maintained: the other 2 stims (same rule throughout). ' ...
     'Reversal drop should appear only for switched, not maintained. ' ...
     'P blocks show attenuated recovery due to noisy feedback (p_trueFB=0.80).'], ...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig1, fullfile(outpath,'simv3_fig1_acc_conf_switch_maintain.pdf'));
fprintf('  Fig 1 saved.\n');

% ── Fig S2: Learning rate (α) by block type and stimulus category ─────────
% Aligned to FIRST ENCOUNTER of each stimulus after reversal (not reversal
% trial number), because the α spike only appears when the switched stim is
% first seen post-reversal (Nassar 2010 eq. 6: first-obs boost).
ENC_WIN = 8; n_enc = 2*ENC_WIN+1; x_enc = -ENC_WIN:ENC_WIN;

sw_alpha_D = NaN(0,n_enc); sw_alpha_P = NaN(0,n_enc);
mn_alpha_D = NaN(0,n_enc); mn_alpha_P = NaN(0,n_enc);
sw_omega_D = NaN(0,n_enc); sw_omega_P = NaN(0,n_enc);

for si = 1:numel(snames)
    td  = sim_data.(snames{si});
    mv2 = model_vars.(snames{si});
    nB  = numel(td.revTrial);
    for b = 1:nB
        bt  = char(td.block_type(b));
        rev = td.revTrial(b);
        if isnan(rev), continue; end
        sw_b = td.switch_stims{b};

        for s = 1:N_STIM
            s_mask = td.stimID(b,:) == s;
            s_trials = find(s_mask);
            s_alpha = mv2.alpha(b, s_trials);
            s_omega = mv2.omega(b, s_trials);

            fp = find(s_trials > round(rev), 1);
            if isempty(fp), continue; end

            row_a = NaN(1,n_enc); row_o = NaN(1,n_enc);
            for ek = 1:n_enc
                ei = fp + (ek - ENC_WIN - 1);
                if ei >= 1 && ei <= numel(s_trials)
                    row_a(ek) = s_alpha(ei);
                    row_o(ek) = s_omega(ei);
                end
            end

            is_sw = ismember(s, sw_b);
            if strcmp(bt,'D')
                if is_sw, sw_alpha_D(end+1,:)=row_a; sw_omega_D(end+1,:)=row_o;
                else,     mn_alpha_D(end+1,:)=row_a; end
            else
                if is_sw, sw_alpha_P(end+1,:)=row_a; sw_omega_P(end+1,:)=row_o;
                else,     mn_alpha_P(end+1,:)=row_a; end
            end
        end
    end
end

fig2 = figure('Position',[50 50 1200 500]);
sgtitle({'Learning rate α and ω — stimulus-encounter-aligned', ...
    'x=0: first time the stimulus appears post-reversal', ...
    '(Reversal-trial alignment smears the spike; encounter alignment isolates it)'}, ...
    'FontSize',10);

ax = subplot(1,3,1); hold(ax,'on'); title(ax,'α — D blocks','FontSize',10);
plot_ribbon_lc(ax,x_enc,sw_alpha_D,clr_sw,'--',sprintf('Switch (n=%d)',size(sw_alpha_D,1)));
plot_ribbon_lc(ax,x_enc,mn_alpha_D,clr_mn,'-', sprintf('Maint. (n=%d)',size(mn_alpha_D,1)));
xline(ax,-0.5,'k--','LineWidth',1.5,'HandleVisibility','off');
xlabel(ax,'Encounter rel. reversal'); ylabel(ax,'\alpha_t');
xlim(ax,[-ENC_WIN ENC_WIN]); ylim(ax,[0 1]);
legend(ax,'Box','off','FontSize',9,'Location','best');

ax = subplot(1,3,2); hold(ax,'on'); title(ax,'α — P blocks','FontSize',10);
plot_ribbon_lc(ax,x_enc,sw_alpha_P,clr_sw,'--',sprintf('Switch (n=%d)',size(sw_alpha_P,1)));
plot_ribbon_lc(ax,x_enc,mn_alpha_P,clr_mn,'-', sprintf('Maint. (n=%d)',size(mn_alpha_P,1)));
xline(ax,-0.5,'k--','LineWidth',1.5,'HandleVisibility','off');
xlabel(ax,'Encounter rel. reversal'); ylabel(ax,'\alpha_t');
xlim(ax,[-ENC_WIN ENC_WIN]); ylim(ax,[0 1]);
legend(ax,'Box','off','FontSize',9,'Location','best');

ax = subplot(1,3,3); hold(ax,'on'); title(ax,'\omega — Switch stims: D vs P','FontSize',10);
plot_ribbon_lc(ax,x_enc,sw_omega_D,clr_D,'-', sprintf('D-block (n=%d)',size(sw_omega_D,1)));
plot_ribbon_lc(ax,x_enc,sw_omega_P,clr_P,'--',sprintf('P-block (n=%d)',size(sw_omega_P,1)));
xline(ax,-0.5,'k--','LineWidth',1.5,'HandleVisibility','off');
xlabel(ax,'Encounter rel. reversal'); ylabel(ax,'\omega_t (CP probability)');
xlim(ax,[-ENC_WIN ENC_WIN]); ylim(ax,[0 1]);
legend(ax,'Box','off','FontSize',9,'Location','best');
annotation('textbox',[0.01 0.01 0.98 0.05],'String',...
    ['ω_t is Bayesian estimate of change-point probability (Nassar 2010). ' ...
     'Spike at encounter 0 for switched stims: first post-reversal observation triggers maximal ω. ' ...
     'P blocks: ω spike is broader/slower because noisy feedback reduces χ_t (surprise ratio). ' ...
     'This captures Yu & Dayan (2005) expected vs unexpected uncertainty distinction.'],...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig2, fullfile(outpath,'simv3_fig2_alpha_omega_encounter.pdf'));
fprintf('  Fig 2 saved.\n');

% ── Fig S3: n_prev_P effect on reversal cost ─────────────────────────────
% Key hypothesis (Behrens et al. 2007): accumulated P experience → higher
% prior H estimate → larger initial α spike at reversal → faster recovery.
% Here we show simulated reversal cost (pre−post accuracy) as a function of
% n_prev_P, separately for D and P blocks.
% NOTE: In this simulation H is fixed; n_prev_P effect requires adaptive H
% (see §5 fitting). This figure shows the BASELINE (no n_prev_P effect)
% to confirm the simulation is well-behaved before fitting.

fig3 = figure('Position',[50 50 900 420]);
sgtitle({'n_{prev_P} vs reversal cost — simulated', ...
    '(Fixed H: no systematic effect expected — confirms simulation baseline)'}, ...
    'FontSize',10);

% Compute reversal cost per block per subject
nP_vals_D = []; rc_D = [];
nP_vals_P = []; rc_P = [];

for si = 1:numel(snames)
    td = sim_data.(snames{si});
    nB = numel(td.revTrial);
    for b = 1:nB
        bt  = char(td.block_type(b));
        rev = td.revTrial(b);
        nP  = td.n_prev_P(b);
        if isnan(rev) || rev < 10 || rev > 90, continue; end

        pre_acc  = mean(td.correct(b,max(1,round(rev)-20):round(rev)-1), 'omitnan');
        post_acc = mean(td.correct(b,round(rev):min(100,round(rev)+19)), 'omitnan');
        rc = pre_acc - post_acc;

        if strcmp(bt,'D')
            nP_vals_D(end+1) = nP; rc_D(end+1) = rc;
        else
            nP_vals_P(end+1) = nP; rc_P(end+1) = rc;
        end
    end
end

ax = subplot(1,2,1); hold(ax,'on');
title(ax,'D blocks: reversal cost vs n_{prev_P}','FontSize',10);
scatter(ax,nP_vals_D,rc_D,30,clr_D,'filled','MarkerFaceAlpha',0.4);
if numel(nP_vals_D)>5
    xfit = linspace(min(nP_vals_D),max(nP_vals_D),50);
    pfit = polyfit(nP_vals_D(:),rc_D(:),1);
    plot(ax,xfit,polyval(pfit,xfit),'k-','LineWidth',1.5,'HandleVisibility','off');
    [rv,pv] = corr(nP_vals_D(:),rc_D(:));
    text(ax,0.05,0.97,sprintf('r=%.2f, p=%.3f\n(null: fixed H)',rv,pv), ...
        'Units','normalized','VerticalAlignment','top','FontSize',8,'BackgroundColor','w');
end
xlabel(ax,'n_{prev\_P} (prior P blocks experienced)');
ylabel(ax,'Reversal cost (pre−post accuracy)');
yline(ax,0,'k:','HandleVisibility','off'); ylim(ax,[-0.3 0.8]);

ax = subplot(1,2,2); hold(ax,'on');
title(ax,'P blocks: reversal cost vs n_{prev_P}','FontSize',10);
scatter(ax,nP_vals_P,rc_P,30,clr_P,'filled','MarkerFaceAlpha',0.4);
if numel(nP_vals_P)>5
    xfit = linspace(min(nP_vals_P),max(nP_vals_P),50);
    pfit = polyfit(nP_vals_P(:),rc_P(:),1);
    plot(ax,xfit,polyval(pfit,xfit),'k-','LineWidth',1.5,'HandleVisibility','off');
    [rv,pv] = corr(nP_vals_P(:),rc_P(:));
    text(ax,0.05,0.97,sprintf('r=%.2f, p=%.3f\n(null: fixed H)',rv,pv), ...
        'Units','normalized','VerticalAlignment','top','FontSize',8,'BackgroundColor','w');
end
xlabel(ax,'n_{prev\_P}');
ylabel(ax,'Reversal cost');
yline(ax,0,'k:','HandleVisibility','off'); ylim(ax,[-0.3 0.8]);

annotation('textbox',[0.01 0.01 0.98 0.04],'String',...
    ['n_prev_P = number of P blocks experienced before this block. ' ...
     'Behrens et al. (2007): prior P experience should increase learning rate H prior, ' ...
     'reducing reversal cost over repeated exposure. In the null simulation (fixed H) no effect is expected. ' ...
     'In real data, a significant negative correlation would support the accumulated-uncertainty hypothesis.'],...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig3, fullfile(outpath,'simv3_fig3_nprevP_revCost.pdf'));
fprintf('  Fig 3 saved.\n');

% ── Fig S4: Block-order effects — simulated accuracy profiles ─────────────
% Shows the four most common empirical block orders side-by-side.
% Demonstrates that the simulation captures design variability in
% uncertainty context sequences.

fig4 = figure('Position',[50 50 1400 500]);
sgtitle({'Block-order effects on accuracy — four most common orders', ...
    '(Simulated from empirical block structures)'}, 'FontSize',11);

order_examples = {'DPDP','PDPD','PPDDP','DDPPD'};
order_colors   = {clr_D, clr_P, [0.5 0.1 0.6], [0.1 0.6 0.5]};

for oi = 1:4
    target_order = order_examples{oi};
    ax = subplot(1,4,oi); hold(ax,'on');
    title(ax,sprintf('Order: %s',target_order),'FontSize',10);

    % Find simulated subjects with this block order
    match_snames = {};
    for si2 = 1:N_SIM_SUBJECTS
        es2 = EMPIRICAL_SUBJECTS{si2};
        if strcmp(es2.block_str, target_order)
            match_snames{end+1} = fieldnames_at(sim_data, si2);
        end
    end

    n_match = numel(match_snames);
    if n_match == 0
        text(ax,0.5,0.5,'No subjects','Units','normalized','HorizontalAlignment','center');
        continue;
    end

    % Pool all blocks across subjects, sorted by position in sequence
    nblks_order = numel(target_order);
    for b = 1:nblks_order
        bt_this = target_order(b);
        block_means = NaN(n_match,1);
        for mi = 1:n_match
            sn2 = match_snames{mi};
            td2 = sim_data.(sn2);
            if b > numel(td2.revTrial), continue; end
            block_means(mi) = mean(td2.correct(b,:),'omitnan');
        end
        xval = b;
        clr_b = ternary(bt_this=='D', clr_D, clr_P);
        scatter(ax,xval*ones(n_match,1)+0.1*(rand(n_match,1)-0.5), ...
            block_means,30,clr_b,'filled','MarkerFaceAlpha',0.5,'HandleVisibility','off');
        errorbar(ax,xval,mean(block_means,'omitnan'),std(block_means,'omitnan'), ...
            'Color',clr_b,'LineWidth',2,'CapSize',8,'Marker','o', ...
            'MarkerSize',8,'MarkerFaceColor',clr_b, ...
            'DisplayName',sprintf('Blk%d(%s)',b,bt_this));
    end

    set(ax,'XTick',1:nblks_order,'XTickLabel',cellstr(target_order(:))');
    xlabel(ax,'Block position'); ylabel(ax,'Mean accuracy');
    ylim(ax,[0 1]); yline(ax,0.5,'k:','HandleVisibility','off');
    legend(ax,'Box','off','FontSize',7,'Location','best');
    subtitle(ax,sprintf('n=%d subjects',n_match),'FontSize',8,'Color',[0.5 0.5 0.5]);
end

annotation('textbox',[0.01 0.01 0.98 0.04],'String',...
    ['Blue bars = D blocks, orange = P blocks. ' ...
     'DPDP and PDPD are the most frequent orders (RR cohort). ' ...
     'PPDDP/DDPPD most frequent in KH. ' ...
     'P-blocks should show lower accuracy due to noisy feedback masking reversal signal.'],...
    'FontSize',7,'EdgeColor','none','BackgroundColor',[0.95 0.95 0.95]);

saveas(fig4, fullfile(outpath,'simv3_fig4_block_order_effects.pdf'));
fprintf('  Fig 4 saved.\n');

% =========================================================================
%% §6  FIT NASSAR TO REAL DATA
% =========================================================================
fprintf('\n=== §6: Fitting Nassar model to real data ===\n');

STIM_FEATURES_FIT = [1 1; 1 2; 2 1; 2 2];
results = fit_all_subjects_v3(all_trial_data, STIM_FEATURES_FIT, outpath);

% =========================================================================
%% §7  FIT SEPARATELY FOR D AND P BLOCKS
% =========================================================================
fprintf('\n=== §7: Block-type-specific fitting ===\n');
results = fit_subjects_by_blocktype_v3(all_trial_data, STIM_FEATURES_FIT, results, outpath);

save(fullfile(outpath,'nassar_results_v3.mat'), 'results');
fprintf('Saved nassar_results_v3.mat\n');

% =========================================================================
%% §8  PARAMETER VISUALISATION
% =========================================================================
subj_ids  = fieldnames(results);
N_subj    = numel(subj_ids);

H_all     = nan(N_subj,1); beta_all  = nan(N_subj,1);
H_det     = nan(N_subj,1); H_prob    = nan(N_subj,1);
H_ns      = nan(N_subj,1); bic_all   = nan(N_subj,1);
nP_mean   = nan(N_subj,1);  % mean n_prev_P across blocks for this subject

for si = 1:N_subj
    sn = subj_ids{si};
    r  = results.(sn);
    H_all(si)   = r.H_fit;
    beta_all(si)= r.beta_fit;
    bic_all(si) = r.bic;
    if isfield(r,'H_fit_det'),         H_det(si)  = r.H_fit_det;           end
    if isfield(r,'H_fit_prob'),        H_prob(si) = r.H_fit_prob;          end
    if isfield(r,'H_noise_sensitivity'), H_ns(si) = r.H_noise_sensitivity; end
    if isfield(all_trial_data, sn)
        td_r = all_trial_data.(sn).trial_data;
        nB_r = size(td_r.correct, 1);
        tfb  = td_r.trueFB;
        np_v = 0; nP_arr = zeros(nB_r,1);
        for b = 1:nB_r
            bfb = tfb(b, ~isnan(tfb(b,:)));
            is_P = ~isempty(bfb) && mean(bfb) < 0.95;
            nP_arr(b) = np_v;
            if is_P, np_v = np_v+1; end
        end
        nP_mean(si) = mean(nP_arr);
    end
end

is_kh = cellfun(@(s) startsWith(s,'Ox'), subj_ids);
is_rr = cellfun(@(s) startsWith(s,'Nc'), subj_ids);

fig5 = figure('Position',[50 50 1400 420]);
sgtitle('Fitted Nassar parameters — empirical data','FontSize',12);

ax = subplot(1,4,1); hold(ax,'on');
title(ax,'H_{det} vs H_{prob}','FontSize',10);
ok = ~isnan(H_det) & ~isnan(H_prob);
scatter(ax,H_det(ok&is_kh),H_prob(ok&is_kh),60,clr_D,'filled','MarkerFaceAlpha',0.8,'DisplayName','KH');
scatter(ax,H_det(ok&is_rr),H_prob(ok&is_rr),60,clr_P,'filled','MarkerFaceAlpha',0.8,'DisplayName','RR');
lims = [0, max([H_det(ok);H_prob(ok)])*1.1];
plot(ax,lims,lims,'k--','HandleVisibility','off');
[~,p_h]=ttest(H_prob(ok)-H_det(ok));
text(ax,0.05,0.97,sprintf('Paired t: p=%.3f',p_h),'Units','normalized','VerticalAlignment','top','FontSize',8,'BackgroundColor','w');
xlabel(ax,'H_{det}'); ylabel(ax,'H_{prob}');
xlim(ax,lims); ylim(ax,lims); axis(ax,'square');
legend(ax,'Box','off','FontSize',8);

ax = subplot(1,4,2); hold(ax,'on');
title(ax,'\DeltaH = H_{prob} − H_{det}','FontSize',10);
ok_ns = ~isnan(H_ns);
histogram(ax,H_ns(ok_ns&is_kh),8,'FaceColor',clr_D,'FaceAlpha',0.7,'EdgeColor','w','DisplayName','KH');
histogram(ax,H_ns(ok_ns&is_rr),8,'FaceColor',clr_P,'FaceAlpha',0.7,'EdgeColor','w','DisplayName','RR');
xline(ax,0,'k--','LineWidth',1.5,'HandleVisibility','off');
[~,p_ns]=ttest(H_ns(ok_ns));
text(ax,0.98,0.97,sprintf('Mean ΔH=%.3f\nt-test: p=%.3f',mean(H_ns(ok_ns),'omitnan'),p_ns),...
    'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',8,'BackgroundColor','w');
xlabel(ax,'\DeltaH (noise sensitivity)'); ylabel(ax,'Subjects');
legend(ax,'Box','off','FontSize',8,'Location','best');

ax = subplot(1,4,3); hold(ax,'on');
title(ax,'H_{prob} vs n_{prev\_P}','FontSize',10);
ok_np = ~isnan(H_prob) & ~isnan(nP_mean);
scatter(ax,nP_mean(ok_np&is_kh),H_prob(ok_np&is_kh),60,clr_D,'filled','MarkerFaceAlpha',0.8,'DisplayName','KH');
scatter(ax,nP_mean(ok_np&is_rr),H_prob(ok_np&is_rr),60,clr_P,'filled','MarkerFaceAlpha',0.8,'DisplayName','RR');
if sum(ok_np)>3
    [rv,pv]=corr(nP_mean(ok_np),H_prob(ok_np),'Rows','complete');
    xfit=linspace(min(nP_mean(ok_np)),max(nP_mean(ok_np)),50);
    plot(ax,xfit,polyval(polyfit(nP_mean(ok_np),H_prob(ok_np),1),xfit),'k-','LineWidth',1.5,'HandleVisibility','off');
    text(ax,0.05,0.97,sprintf('r=%.2f, p=%.3f',rv,pv),'Units','normalized','VerticalAlignment','top','FontSize',8,'BackgroundColor','w');
end
xlabel(ax,'Mean n_{prev\_P}'); ylabel(ax,'H_{prob}');
legend(ax,'Box','off','FontSize',8);
subtitle(ax,'Accumulated P experience → higher fitted H?','FontSize',8,'Color',[0.5 0.5 0.5]);

ax = subplot(1,4,4); hold(ax,'on');
title(ax,'BIC per subject','FontSize',10);
[bic_s,idx_s]=sort(bic_all,'ascend');
bar(ax,1:N_subj,bic_s,'FaceColor',0.6*[1 1 1],'EdgeColor','none');
yline(ax,mean(bic_all,'omitnan'),'r--','LineWidth',1.5,'HandleVisibility','off');
text(ax,0.02,0.97,sprintf('Mean BIC=%.1f',mean(bic_all,'omitnan')),'Units','normalized','VerticalAlignment','top','FontSize',8,'Color','r');
set(ax,'XTick',1:5:N_subj); xlabel(ax,'Subject (sorted)'); ylabel(ax,'BIC');

saveas(fig5, fullfile(outpath,'simv3_fig5_params.pdf'));
fprintf('  Fig 5 saved.\n');

% =========================================================================
%% §9  PARAMETER RECOVERY
% =========================================================================
fprintf('\n=== §9: Parameter recovery ===\n');

H_grid_r    = logspace(log10(0.005), log10(0.40), 30);
beta_grid_r = linspace(1, 20, 30);
H_fit_orig  = nan(N_subj,1); beta_fit_orig = nan(N_subj,1);
H_recovered = nan(N_subj,1); beta_recovered= nan(N_subj,1);

for si = 1:N_subj
    sn = subj_ids{si};
    r  = results.(sn);
    H_fit_orig(si)   = r.H_fit;
    beta_fit_orig(si)= r.beta_fit;

    % Find matching empirical structure for this subject
    idx = find(cellfun(@(e) strcmp(e.sn,sn), EMPIRICAL_SUBJECTS),1);
    if isempty(idx)
        warning('%s not in empirical registry, skipping recovery.',sn); continue;
    end
    es_r = EMPIRICAL_SUBJECTS{idx};

    % Simulate synthetic dataset from fitted params, matching THIS subject's
    % exact block structure (block order, n_blocks, switch_stims)
    nB_r = es_r.nblks;
    td_syn = simulate_from_registry_entry(es_r, r.H_fit, r.beta_fit, ...
        N_TRIALS, N_STIM, STIM_FEATURES_FIT, P_TRUE_FB_D, P_TRUE_FB_P, ...
        REV_D_MU, REV_D_SIG, REV_D_LO, REV_D_HI, ...
        REV_P_MU, REV_P_SIG, REV_P_LO, REV_P_HI, ALPHA_INIT);

    obs_syn = pack_observations_v3(td_syn, 'stimID');
    if isempty(obs_syn.stim_id), continue; end

    % Grid search + refine
    nll_g = nan(numel(H_grid_r),numel(beta_grid_r));
    for hi=1:numel(H_grid_r)
        for bi2=1:numel(beta_grid_r)
            nll_g(hi,bi2)=nassar_nll_v3([H_grid_r(hi),beta_grid_r(bi2)],obs_syn,STIM_FEATURES_FIT);
        end
    end
    [~,bl]=min(nll_g(:)); [hi0,bi0]=ind2sub(size(nll_g),bl);
    x0=[H_grid_r(hi0),beta_grid_r(bi0)];
    opts=optimoptions('fmincon','Display','off','MaxIterations',500);
    try
        xr=fmincon(@(x) nassar_nll_v3(x,obs_syn,STIM_FEATURES_FIT),x0,[],[],[],[],...
            [0.001,0.5],[0.40,25],[],opts);
    catch
        xr=x0;
    end
    H_recovered(si)    = xr(1);
    beta_recovered(si) = xr(2);
    fprintf('  %s: H %.3f→%.3f | β %.2f→%.2f\n',sn,r.H_fit,H_recovered(si),r.beta_fit,beta_recovered(si));
end

fig6=figure('Position',[50 50 850 380]);
sgtitle('Parameter recovery — same subject block structure used for synthesis','FontSize',11);

ax=subplot(1,2,1); hold(ax,'on'); axis(ax,'square'); title(ax,'H recovery');
ok_r=~isnan(H_fit_orig)&~isnan(H_recovered);
scatter(ax,H_fit_orig(ok_r&is_kh),H_recovered(ok_r&is_kh),60,clr_D,'filled','MarkerFaceAlpha',0.8,'DisplayName','KH');
scatter(ax,H_fit_orig(ok_r&is_rr),H_recovered(ok_r&is_rr),60,clr_P,'filled','MarkerFaceAlpha',0.8,'DisplayName','RR');
lims_r=[0,max([H_fit_orig(ok_r);H_recovered(ok_r)])+0.02];
plot(ax,lims_r,lims_r,'k--'); [rv,pv]=corr(H_fit_orig(ok_r),H_recovered(ok_r));
text(ax,0.05,0.95,sprintf('r=%.2f, p=%.3f',rv,pv),'Units','normalized','VerticalAlignment','top','FontSize',9,'BackgroundColor','w');
xlabel(ax,'H fitted (real data)'); ylabel(ax,'H recovered (synthetic)');
xlim(ax,lims_r); ylim(ax,lims_r); legend(ax,'Box','off','FontSize',8);

ax=subplot(1,2,2); hold(ax,'on'); axis(ax,'square'); title(ax,'β recovery');
ok_rb=~isnan(beta_fit_orig)&~isnan(beta_recovered);
scatter(ax,beta_fit_orig(ok_rb&is_kh),beta_recovered(ok_rb&is_kh),60,clr_D,'filled','MarkerFaceAlpha',0.8,'DisplayName','KH');
scatter(ax,beta_fit_orig(ok_rb&is_rr),beta_recovered(ok_rb&is_rr),60,clr_P,'filled','MarkerFaceAlpha',0.8,'DisplayName','RR');
lims_b=[0,max([beta_fit_orig(ok_rb);beta_recovered(ok_rb)])+1];
plot(ax,lims_b,lims_b,'k--'); [rv,pv]=corr(beta_fit_orig(ok_rb),beta_recovered(ok_rb));
text(ax,0.05,0.95,sprintf('r=%.2f, p=%.3f',rv,pv),'Units','normalized','VerticalAlignment','top','FontSize',9,'BackgroundColor','w');
xlabel(ax,'\beta fitted'); ylabel(ax,'\beta recovered');
xlim(ax,lims_b); ylim(ax,lims_b); legend(ax,'Box','off','FontSize',8);

saveas(fig6, fullfile(outpath,'simv3_fig6_param_recovery.pdf'));
fprintf('  Fig 6 saved.\n');

% =========================================================================
%% §10 EXPORT: add Nassar latents to group_T
% =========================================================================
fprintf('\n=== §10: Exporting Nassar latents to group_T ===\n');
if ~isempty(group_T)
    group_T = add_nassar_latents_v3(group_T, results);
    save(fullfile(data_path,'behav_table_June2026_RL_v3.mat'), 'group_T');
    fprintf('Saved enriched group_T: %d rows, %d cols\n', height(group_T), width(group_T));
else
    fprintf('  No group_T available — skipping export.\n');
end

fprintf('\n=== v3 complete ===\n');
fprintf('Figures → %s\n', outpath);
fprintf('Results → nassar_results_v3.mat\n');
fprintf('Sim data → sim_data_v3.mat\n');

% =========================================================================
%% LOCAL FUNCTIONS
% =========================================================================

function reg = build_empirical_registry(all_trial_data)
% Build per-subject block registry from real data.
% Returns cell array of structs, one per subject (Ox02 excluded).
% Block type inferred from mean(trueFB) per block: >0.95 = D, else = P.
% V blocks (visual probabilistic) also classified as P.

subjects = fieldnames(all_trial_data);
reg = {};

for si = 1:numel(subjects)
    sn = subjects{si};
    if strcmp(sn,'Ox02'), continue; end  % incomplete (only 2 blocks)

    subj = all_trial_data.(sn);
    td   = subj.trial_data;
    nblks = size(td.correct, 1);

    tfb = td.trueFB;
    rt  = td.revTrial(:);
    sw  = td.switch_stims;  % nblks × 2

    block_types   = cell(nblks,1);
    switch_stims  = cell(nblks,1);
    rev_trials    = nan(nblks,1);

    for b = 1:nblks
        bfb = tfb(b, ~isnan(tfb(b,:)));
        if isempty(bfb),    bt='P';
        elseif mean(bfb)>0.95, bt='D';
        else,               bt='P';  % P and V both classified P
        end
        block_types{b} = bt;

        if size(sw,1) >= b
            sw_b = sw(b,:);
            if any(isnan(sw_b))
                sw_b = [2 3];  % fallback for 3 blocks with NaN
            end
            switch_stims{b} = sort(double(sw_b));
        else
            switch_stims{b} = [2 3];
        end

        if b <= numel(rt) && ~isnan(rt(b))
            rev_trials(b) = rt(b);
        end
    end

    block_str = strjoin(block_types,'');
    cohort = ternary(startsWith(sn,'Ox'),'KH','RR');

    entry.sn          = sn;
    entry.cohort      = cohort;
    entry.nblks       = nblks;
    entry.block_str   = block_str;
    entry.block_types = block_types;
    entry.switch_stims= switch_stims;
    entry.rev_trials  = rev_trials;
    reg{end+1}        = entry;
end
end

function td = simulate_from_registry_entry(es, H, beta, N_TRIALS, N_STIM, ...
    STIM_FEATURES, P_TRUE_FB_D, P_TRUE_FB_P, ...
    REV_D_MU, REV_D_SIG, REV_D_LO, REV_D_HI, ...
    REV_P_MU, REV_P_SIG, REV_P_LO, REV_P_HI, ALPHA_INIT)
% Simulate one subject's dataset matching their empirical block structure.

nB = es.nblks;
td = init_trial_data(nB, N_TRIALS);
td.block_structure  = es.block_str;
td.block_type_label = es.block_types;
td.cohort           = es.cohort;
td.n_prev_P         = zeros(nB,1);

running_P = 0;
for b = 1:nB
    curr_type = es.block_types{b};
    sw_b      = es.switch_stims{b};

    td.n_prev_P(b) = running_P;

    p_trueFB = ternary(strcmp(curr_type,'D'), P_TRUE_FB_D, P_TRUE_FB_P);
    td.trueFB_block(b) = p_trueFB;

    if strcmp(curr_type,'D')
        revTrial = sample_truncated_normal(REV_D_MU, REV_D_SIG, REV_D_LO, REV_D_HI);
    else
        revTrial = sample_truncated_normal(REV_P_MU, REV_P_SIG, REV_P_LO, REV_P_HI);
    end
    td.revTrial(b) = revTrial;
    td.block_type(b) = string(curr_type);

    stim_order = repmat(1:N_STIM, 1, N_TRIALS/N_STIM);
    stim_order = stim_order(randperm(N_TRIALS));

    go_nogo = compute_go_nogo_schedule(stim_order, revTrial, sw_b, N_TRIALS);

    n_true_fb = round(p_trueFB * N_TRIALS);
    fb_arr    = [true(1,n_true_fb), false(1,N_TRIALS-n_true_fb)];
    trueFB_vec= double(fb_arr(randperm(N_TRIALS)));

    theta     = 0.5 * ones(1, N_STIM);
    n_eff     = ones(1, N_STIM);
    first_obs = true(1, N_STIM);

    for t = 1:N_TRIALS
        s_t     = stim_order(t);
        go_true = go_nogo(t);
        fb_true = trueFB_vec(t);
        certainty_t = abs(theta(s_t)-0.5);
        p_go    = 1/(1+exp(-beta*(theta(s_t)-0.5)));
        resp_go = double(rand<p_go);
        correct = double(resp_go==go_true);
        perc_cor= fb_true*correct+(1-fb_true)*(1-correct);
        y_t     = double((resp_go==1&&perc_cor==1)||(resp_go==0&&perc_cor==0));
        pred_prob=max(theta(s_t)^y_t*(1-theta(s_t))^(1-y_t),1e-6);
        chi_t   = 0.5/pred_prob;
        omega_t = (H*chi_t)/(H*chi_t+(1-H));
        n_eff(s_t)=(1-omega_t)*(n_eff(s_t)+1)+omega_t;
        alpha_t = omega_t+(1-omega_t)/n_eff(s_t);
        if first_obs(s_t), alpha_t=max(alpha_t,ALPHA_INIT); first_obs(s_t)=false; end
        delta_t = y_t-theta(s_t);
        theta(s_t)=max(0.01,min(0.99,theta(s_t)+alpha_t*delta_t));
        conf_raw=1+18*certainty_t+0.8*randn;
        td.stimID(b,t)          =s_t;
        td.switch_stims{b}      =sw_b;
        td.is_switch_stim(b,t)  =ismember(s_t,sw_b);
        td.goTrial(b,t)         =go_true;
        td.respWasGo(b,t)       =resp_go;
        td.correct(b,t)         =correct;
        td.perceivedCorrect(b,t)=perc_cor;
        td.trueFB(b,t)          =fb_true;
        td.p_trueFB(b,t)        =p_trueFB;
        td.confidence(b,t)      =max(1,min(10,round(conf_raw)));
        td.blocknum(b,t)        =b;
        td.trialnum(b,t)        =t;
    end

    if strcmp(curr_type,'P'), running_P=running_P+1; end
end
end

function go_nogo = compute_go_nogo_schedule(stim_order, revTrial, sw_stims, N_TRIALS)
% Compute the Go/NoGo schedule given which stims switch at reversal.
% Pre-reversal rule: stims in sw_stims are Go (they will switch to NoGo).
% Post-reversal rule: stims in maint_stims that were NoGo become Go;
%   sw_stims that were Go become NoGo.
%
% This implements the 2×2 dimensional shift:
%   Switched stims: reverse their Go/NoGo assignment at reversal.
%   Maintained stims: keep the same assignment throughout.
%
% The pre-reversal assignment of the switch pair can be either
% Go-pre/NoGo-post or NoGo-pre/Go-post. Here we use the convention
% that the FIRST element of sw_stims is Go pre-reversal and
% the SECOND is NoGo pre-reversal (matching your detect_reversal_KH
% logic: cond_A/cond_B detect either config).

go_nogo = zeros(1, N_TRIALS);
all_stims = 1:4;
maint_stims = setdiff(all_stims, sw_stims);

% Assign pre-reversal Go identity:
% sw_stims(1) = Go, sw_stims(2) = NoGo
% maint_stims(1) = Go, maint_stims(2) = NoGo
% (four stims: two Go, two NoGo, defined by which pair is switch vs maint)
pre_go  = [sw_stims(1), maint_stims(1)];
post_go = [sw_stims(2), maint_stims(1)];  % sw reverses, maint stays

for t = 1:N_TRIALS
    s = stim_order(t);
    if t <= revTrial
        go_nogo(t) = double(ismember(s, pre_go));
    else
        go_nogo(t) = double(ismember(s, post_go));
    end
end
end

function rev = sample_truncated_normal(mu, sigma, lo, hi)
% Sample from a Normal(mu,sigma) truncated to [lo,hi].
% Uses rejection sampling (efficient for wide tails relative to truncation).
for attempt = 1:1000
    v = round(mu + sigma * randn);
    if v >= lo && v <= hi
        rev = v;
        return;
    end
end
rev = round((lo+hi)/2);  % fallback — extremely unlikely
end

function td = init_trial_data(nBlocks, nTrials)
fields = {'stimID','is_switch_stim','goTrial','respWasGo','correct',...
          'perceivedCorrect','trueFB','p_trueFB','confidence','blocknum','trialnum'};
for f = fields
    td.(f{1}) = nan(nBlocks, nTrials);
end
td.revTrial     = nan(1,nBlocks);
td.trueFB_block = nan(1,nBlocks);
td.n_prev_P     = zeros(nBlocks,1);
td.block_type   = strings(nBlocks,1);
td.switch_stims = cell(nBlocks,1);
end

function mv = init_model_vars(nBlocks, nTrials, nStim)
mv.theta    = nan(nBlocks, nTrials, nStim);
mv.alpha    = nan(nBlocks, nTrials);
mv.omega    = nan(nBlocks, nTrials);
mv.delta    = nan(nBlocks, nTrials);
mv.certainty= nan(nBlocks, nTrials);
mv.n_eff    = nan(nBlocks, nTrials);
mv.surprise = nan(nBlocks, nTrials);
end

function obs = pack_observations_v3(td, stim_field)
obs.stim_id=[]; obs.block_id=[]; obs.trial_id=[];
obs.resp_go=[]; obs.y_t=[]; obs.is_det=[];

[nB,nT]=size(td.correct);
for b=1:nB
    bt=char(td.block_type(b)); is_det=strcmp(bt,'D');
    for t=1:nT
        cor=td.correct(b,t); gt=td.goTrial(b,t);
        if isnan(cor)||isnan(gt), continue; end
        sv=td.(stim_field)(b,t); if isnan(sv), continue; end
        rgo=gt*cor+(1-gt)*(1-cor);
        pc=td.perceivedCorrect(b,t); if isnan(pc), pc=cor; end
        yt=double((rgo==1&&pc==1)||(rgo==0&&pc==0));
        obs.stim_id(end+1)=sv;  obs.block_id(end+1)=b;
        obs.trial_id(end+1)=t;  obs.resp_go(end+1)=rgo;
        obs.y_t(end+1)=yt;      obs.is_det(end+1)=double(is_det);
    end
end
end

function [nll,trials_out]=nassar_nll_v3(params, obs, stim_features, return_trials)
if nargin<4, return_trials=false; end
H=params(1); beta=params(2);
N=numel(obs.stim_id); Ns=size(stim_features,1);
theta=0.5*ones(1,Ns); n_eff=ones(1,Ns); first_obs=true(1,Ns);
ll=0; prev_block=-1;
if return_trials
    tr_alpha=nan(1,N); tr_delta=nan(1,N); tr_omega=nan(1,N);
    tr_cert=nan(1,N); tr_theta=nan(1,N);
end
for t=1:N
    s_t=obs.stim_id(t);
    if obs.block_id(t)~=prev_block
        theta=0.5*ones(1,Ns); n_eff=ones(1,Ns); first_obs=true(1,Ns);
        prev_block=obs.block_id(t);
    end
    if isnan(s_t)||s_t<1||s_t>Ns, continue; end
    cert_t=abs(theta(s_t)-0.5);
    p_go=1/(1+exp(-beta*(theta(s_t)-0.5)));
    rgo=obs.resp_go(t);
    ll=ll+log(max(p_go^rgo*(1-p_go)^(1-rgo),1e-10));
    y_t=obs.y_t(t);
    pp=max(theta(s_t)^y_t*(1-theta(s_t))^(1-y_t),1e-6);
    chi_t=0.5/pp; omega_t=(H*chi_t)/(H*chi_t+(1-H));
    n_eff(s_t)=(1-omega_t)*(n_eff(s_t)+1)+omega_t;
    alpha_t=omega_t+(1-omega_t)/n_eff(s_t);
    if first_obs(s_t), alpha_t=max(alpha_t,0.5); first_obs(s_t)=false; end
    delta_t=y_t-theta(s_t);
    theta(s_t)=max(0.01,min(0.99,theta(s_t)+alpha_t*delta_t));
    if return_trials
        tr_alpha(t)=alpha_t; tr_delta(t)=delta_t; tr_omega(t)=omega_t;
        tr_cert(t)=cert_t; tr_theta(t)=theta(s_t);
    end
end
nll=-ll;
if return_trials
    trials_out.alpha=tr_alpha; trials_out.delta=tr_delta;
    trials_out.omega=tr_omega; trials_out.certainty=tr_cert;
    trials_out.theta=tr_theta;
else
    trials_out=[];
end
end

function results=fit_all_subjects_v3(all_trial_data, stim_features, outdir)
subjects=fieldnames(all_trial_data);
results=struct();
H_grid=logspace(log10(0.01),log10(0.30),25);
beta_grid=linspace(2,15,25);

for si=1:numel(subjects)
    sn=subjects{si};
    td=all_trial_data.(sn).trial_data;
    sf=ternary(isfield(td,'stimType'),'stimType','stimID');
    has_perc=isfield(td,'perceivedCorrect')&&~all(isnan(td.perceivedCorrect(:)));
    obs=pack_obs_from_real(td, sf, has_perc);
    if isempty(obs.stim_id), continue; end
    fprintf('  Fitting %s...\n',sn);
    nll_g=nan(numel(H_grid),numel(beta_grid));
    for hi=1:numel(H_grid)
        for bi=1:numel(beta_grid)
            nll_g(hi,bi)=nassar_nll_v3([H_grid(hi),beta_grid(bi)],obs,stim_features);
        end
    end
    [~,bl]=min(nll_g(:)); [hi0,bi0]=ind2sub(size(nll_g),bl);
    x0=[H_grid(hi0),beta_grid(bi0)];
    opts=optimoptions('fmincon','Display','off','MaxIterations',500);
    try
        xf=fmincon(@(x) nassar_nll_v3(x,obs,stim_features),x0,[],[],[],[],...
            [0.001,0.5],[0.50,30],[],opts);
    catch
        xf=x0;
    end
    [nll_fit,tr]=nassar_nll_v3(xf,obs,stim_features,true);
    N_obs=numel(obs.stim_id);
    results.(sn).H_fit=xf(1); results.(sn).beta_fit=xf(2);
    results.(sn).nll=nll_fit; results.(sn).bic=2*nll_fit+2*log(N_obs);
    results.(sn).N_obs=N_obs;
    results.(sn).alpha_trial=tr.alpha; results.(sn).delta_trial=tr.delta;
    results.(sn).omega_trial=tr.omega; results.(sn).certainty_trial=tr.certainty;
    results.(sn).theta_trial=tr.theta;
    results.(sn).surprise=tr.omega.*abs(tr.delta);
    results.(sn).block_id=obs.block_id; results.(sn).trial_id=obs.trial_id;
    results.(sn).stim_id=obs.stim_id; results.(sn).is_det=obs.is_det;
    fprintf('    H=%.3f β=%.2f NLL=%.1f BIC=%.1f\n',xf(1),xf(2),nll_fit,results.(sn).bic);
end
end

function results=fit_subjects_by_blocktype_v3(all_trial_data, stim_features, results, outdir)
subjects=fieldnames(results);
H_grid=logspace(log10(0.005),log10(0.40),30);
beta_grid=linspace(1,20,30);

for si=1:numel(subjects)
    sn=subjects{si};
    td=all_trial_data.(sn).trial_data;
    sf=ternary(isfield(td,'stimType'),'stimType','stimID');
    has_perc=isfield(td,'perceivedCorrect')&&~all(isnan(td.perceivedCorrect(:)));
    obs_all=pack_obs_from_real(td, sf, has_perc);
    if isempty(obs_all.stim_id), continue; end

    tfb=td.trueFB; nB=size(td.correct,1);
    block_is_det=false(1,nB);
    for b=1:nB
        bfb=tfb(b,~isnan(tfb(b,:)));
        block_is_det(b)=~isempty(bfb)&&mean(bfb)>0.95;
    end
    valid_b=obs_all.block_id>=1&obs_all.block_id<=nB;
    obs_all.is_det(valid_b)=block_is_det(obs_all.block_id(valid_b));
    results.(sn).block_is_det=block_is_det;

    for bt=1:2
        is_det_fit=(bt==1); lbl=ternary(is_det_fit,'det','prob');
        obs_f.stim_id  = obs_all.stim_id(obs_all.is_det==double(is_det_fit));
        obs_f.block_id = obs_all.block_id(obs_all.is_det==double(is_det_fit));
        obs_f.trial_id = obs_all.trial_id(obs_all.is_det==double(is_det_fit));
        obs_f.resp_go  = obs_all.resp_go(obs_all.is_det==double(is_det_fit));
        obs_f.y_t      = obs_all.y_t(obs_all.is_det==double(is_det_fit));
        obs_f.is_det   = obs_all.is_det(obs_all.is_det==double(is_det_fit));
        if numel(obs_f.stim_id)<20
            results.(sn).(['H_fit_' lbl])=NaN;
            results.(sn).(['beta_fit_' lbl])=NaN;
            continue;
        end
        nll_g=nan(numel(H_grid),numel(beta_grid));
        for hi=1:numel(H_grid)
            for bi=1:numel(beta_grid)
                nll_g(hi,bi)=nassar_nll_v3([H_grid(hi),beta_grid(bi)],obs_f,stim_features);
            end
        end
        [~,bl]=min(nll_g(:)); [hi0,bi0]=ind2sub(size(nll_g),bl);
        x0=[H_grid(hi0),beta_grid(bi0)];
        opts=optimoptions('fmincon','Display','off','MaxIterations',500);
        try
            xf=fmincon(@(x) nassar_nll_v3(x,obs_f,stim_features),x0,[],[],[],[],...
                [0.005,0.5],[0.40,25],[],opts);
        catch, xf=x0; end
        results.(sn).(['H_fit_'    lbl])=xf(1);
        results.(sn).(['beta_fit_' lbl])=xf(2);
        fprintf('  %s [%s]: H=%.3f β=%.2f\n',sn,upper(lbl),xf(1),xf(2));
    end
    results.(sn).H_noise_sensitivity = ...
        results.(sn).H_fit_prob - results.(sn).H_fit_det;
end
end

function obs=pack_obs_from_real(td, stim_field, has_perc)
obs.stim_id=[]; obs.block_id=[]; obs.trial_id=[];
obs.resp_go=[]; obs.y_t=[]; obs.is_det=[];
[nB,nT]=size(td.correct);
for b=1:nB
    tfb_b=td.trueFB(b,~isnan(td.trueFB(b,:)));
    is_det_b=~isempty(tfb_b)&&mean(tfb_b)>0.95;
    for t=1:nT
        cor=td.correct(b,t); gt=td.goTrial(b,t);
        if isnan(cor)||isnan(gt), continue; end
        sv=td.(stim_field)(b,t); if isnan(sv), continue; end
        rgo=gt*cor+(1-gt)*(1-cor);
        if has_perc&&isfield(td,'perceivedCorrect')&&~isnan(td.perceivedCorrect(b,t))
            pc=td.perceivedCorrect(b,t);
        else
            pc=cor;
        end
        yt=double((rgo==1&&pc==1)||(rgo==0&&pc==0));
        obs.stim_id(end+1)=sv; obs.block_id(end+1)=b;
        obs.trial_id(end+1)=t; obs.resp_go(end+1)=rgo;
        obs.y_t(end+1)=yt; obs.is_det(end+1)=double(is_det_b);
    end
end
end

function group_T=add_nassar_latents_v3(group_T, results)
cols={'PE_nassar','omega','alpha_nassar','certainty','surprise','theta_nassar'};
for c=cols, group_T.(c{1})=nan(height(group_T),1); end
subj_col=ternary(ismember('subjID',group_T.Properties.VariableNames),'subjID','subj_id');
subjs=unique(string(group_T.(subj_col)));
for si=1:numel(subjs)
    sn=char(subjs(si));
    if ~isfield(results,sn), continue; end
    r=results.(sn);
    sm=string(group_T.(subj_col))==string(sn);
    sub_rows=find(sm);
    blocks_gt=double(group_T.block(sub_rows));
    trials_gt=double(group_T.trial(sub_rows));
    for t=1:numel(r.trial_id)
        b_t=r.block_id(t); tr_t=r.trial_id(t);
        match=sub_rows(blocks_gt==b_t&trials_gt==tr_t);
        if isempty(match), continue; end
        group_T.PE_nassar(match)     =r.delta_trial(t);
        group_T.omega(match)         =r.omega_trial(t);
        group_T.alpha_nassar(match)  =r.alpha_trial(t);
        group_T.certainty(match)     =r.certainty_trial(t);
        group_T.surprise(match)      =r.surprise(t);
        group_T.theta_nassar(match)  =r.theta_trial(t);
    end
end
end

function sn=fieldnames_at(sim_data, si)
fn=fieldnames(sim_data); sn=fn{si};
end

function out=ternary(cond,a,b)
if cond, out=a; else, out=b; end
end

function plot_ribbon_lc(ax,x,mat,clr,ls,lbl)
if isempty(mat)||all(isnan(mat(:))), return; end
mn=movmean(mean(mat,1,'omitnan'),3,'omitnan');
se=std(mat,0,1,'omitnan')./sqrt(max(sum(~isnan(mat),1),1));
fill(ax,[x,fliplr(x)],[mn+se,fliplr(mn-se)],clr,'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
plot(ax,x,mn,'Color',clr,'LineWidth',2,'LineStyle',ls,'DisplayName',lbl);
end

function annotation_text(ax,txt,fsize)
text(ax,0.03,0.07,txt,'Units','normalized','FontSize',fsize,'Color',[0.5 0.5 0.5],...
    'VerticalAlignment','bottom','Interpreter','none');
end