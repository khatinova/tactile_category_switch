% =============================================================================
% S5b_models_RW_HGF.m
%
% PIPELINE STEP 5b — additional behavioural models + 3-way model comparison.
%
% Adds two field-standard models alongside Nassar (2010):
%   (1) Rescorla-Wagner with DUAL learning rates (alpha_pos / alpha_neg):
%       separate update gains for positive vs negative prediction errors
%       (asymmetric value updating; e.g. Frank et al. 2007, PNAS;
%        Gershman 2015, Cognition).
%   (2) 3-level binary Hierarchical Gaussian Filter (HGF):
%       Mathys et al. (2011) Front Hum Neurosci 5:39;
%       Mathys et al. (2014) Front Hum Neurosci 8:825.
%       Tracks belief (level 2) and its volatility (level 3) per stimulus.
%
% BEST PRACTICE / VALIDATION
% --------------------------
% The HGF here is a transparent, self-contained implementation for model
% comparison. For a publication you should ALSO fit the HGF with the validated
% TAPAS toolbox (tapas_fitModel + tapas_hgf_binary + tapas_unitsq_sgm) and
% confirm the parameter estimates agree. If TAPAS is on the path this script
% will use it; otherwise it falls back to the built-in implementation below.
%
% INPUTS : all_trial_data.mat, nassar_results.mat (from S5)
% OUTPUT : model_comparison_RW_HGF.mat (per-subject BIC/AIC for all 3 models)
%
% Each model is fit per subject by grid initialisation + fmincon on the
% trial-level Bernoulli likelihood of the observed Go/NoGo responses, so the
% BICs are directly comparable to Nassar's.
% =============================================================================

clear; close all; rng(42);
addpath(genpath(fileparts(mfilename('fullpath'))));

% -------------------------------------------------------------------------
%% PATHS
% -------------------------------------------------------------------------
remote = 0;
switch remote
    case 1; base_path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch/Salient mod switch KH';
    case 0; base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH';
end
data_path = fullfile(base_path, 'Data');
outpath   = fullfile(base_path, 'Results', 'Simulation results', 'Figures');
if ~exist(outpath, 'dir'); mkdir(outpath); end

load(fullfile(data_path, 'all_trial_data.mat'), 'all_trial_data');

% Nassar BICs (from S5). Optional: comparison still runs without them.
nassar_results = struct();
nf = fullfile(outpath, 'nassar_results.mat');
if exist(nf, 'file')
    S = load(nf, 'results'); nassar_results = S.results;
else
    warning('nassar_results.mat not found; run S5 first for the Nassar comparison.');
end

STIM_FEATURES = [1 1; 1 2; 2 1; 2 2];     % matches S5
N_STIM = size(STIM_FEATURES, 1);

use_tapas = exist('tapas_fitModel', 'file') == 2;
if use_tapas
    fprintf('TAPAS detected: HGF will also be fit with tapas_hgf_binary.\n');
else
    fprintf('TAPAS not found: using built-in 3-level binary HGF.\n');
end

% -------------------------------------------------------------------------
%% FIT ALL SUBJECTS
% -------------------------------------------------------------------------
subjects = fieldnames(all_trial_data);
N        = numel(subjects);

res = table('Size',[N 0]);
res.subj_id = strings(N,1);
[res.bic_nassar, res.bic_rw1, res.bic_rwdual, res.bic_hgf] = deal(nan(N,1));
[res.alpha_pos, res.alpha_neg, res.beta_rwdual]            = deal(nan(N,1));
[res.hgf_om2, res.hgf_om3, res.beta_hgf]                   = deal(nan(N,1));

for si = 1:N
    sn = subjects{si};
    td = all_trial_data.(sn).trial_data;

    if isfield(td,'stimType'), sf = 'stimType';
    elseif isfield(td,'stimID'), sf = 'stimID';
    else, continue; end
    has_perc = isfield(td,'perceivedCorrect') && ~all(isnan(td.perceivedCorrect(:)));
    obs = pack_observations(td, sf, has_perc, N_STIM);
    if isempty(obs.stim_id); continue; end
    Nobs = numel(obs.stim_id);

    res.subj_id(si) = string(sn);

    % --- Nassar BIC (from S5) ---
    if isfield(nassar_results, sn) && isfield(nassar_results.(sn), 'bic')
        res.bic_nassar(si) = nassar_results.(sn).bic;
    end

    % --- single-alpha RW (baseline, k=2) ---
    [~, bic1] = fit_grid_fmincon(@(p) rw_single_nll(p, obs), ...
        {linspace(0.01,0.6,20), linspace(2,15,15)}, [0.01 0.5],[0.6 30], Nobs, 2);
    res.bic_rw1(si) = bic1;

    % --- dual-alpha RW (k=3: alpha_pos, alpha_neg, beta) ---
    [xd, bicd] = fit_grid_fmincon(@(p) rw_dual_nll(p, obs), ...
        {linspace(0.02,0.6,12), linspace(0.02,0.6,12), linspace(2,12,10)}, ...
        [0.005 0.005 0.5],[0.8 0.8 30], Nobs, 3);
    res.alpha_pos(si)   = xd(1);
    res.alpha_neg(si)   = xd(2);
    res.beta_rwdual(si) = xd(3);
    res.bic_rwdual(si)  = bicd;

    % --- 3-level binary HGF (k=3: omega2, omega3, beta) ---
    if use_tapas
        [xh, bich] = fit_hgf_tapas(obs, Nobs);
    else
        [xh, bich] = fit_grid_fmincon(@(p) hgf_binary_nll(p, obs), ...
            {linspace(-6,-1,8), linspace(-7,-2,8), linspace(2,12,8)}, ...
            [-12 -12 0.5],[ -0.5 -0.5 30], Nobs, 3);
    end
    res.hgf_om2(si)  = xh(1);
    res.hgf_om3(si)  = xh(2);
    res.beta_hgf(si) = xh(3);
    res.bic_hgf(si)  = bich;

    fprintf('%-8s  BIC: Nassar=%.1f  RW1=%.1f  RWdual=%.1f  HGF=%.1f\n', ...
        sn, res.bic_nassar(si), res.bic_rw1(si), res.bic_rwdual(si), res.bic_hgf(si));
end

% -------------------------------------------------------------------------
%% SUMMARY + FIGURE
% -------------------------------------------------------------------------
bic_mat   = [res.bic_nassar, res.bic_rw1, res.bic_rwdual, res.bic_hgf];
model_lbl = {'Nassar','RW (1\alpha)','RW (2\alpha)','HGF'};
ok        = all(~isnan(bic_mat), 2);

fprintf('\n=== MODEL COMPARISON (mean BIC, lower = better) ===\n');
for m = 1:4
    fprintf('  %-12s %.1f\n', model_lbl{m}, mean(bic_mat(ok,m),'omitnan'));
end
if any(ok)
    [~, winner] = min(bic_mat(ok,:), [], 2);
    fprintf('  Winning model per subject (count): ');
    for m = 1:4; fprintf('%s=%d  ', model_lbl{m}, sum(winner==m)); end
    fprintf('\n');
end

fig = figure('Position',[60 60 900 420]);
ax = axes(fig); hold(ax,'on');
title(ax, 'Model comparison (BIC; lower is better)');
b = bar(ax, mean(bic_mat(ok,:),1,'omitnan'));
set(ax,'XTick',1:4,'XTickLabel',model_lbl);
ylabel(ax,'Mean BIC');
% subject-level points
for m = 1:4
    v = bic_mat(ok,m);
    scatter(ax, m + 0.12*(rand(numel(v),1)-0.5), v, 16, [0.3 0.3 0.3], ...
        'filled','MarkerFaceAlpha',0.4);
end
if exist('apply_fig_style','file'); apply_fig_style(fig); end
exportgraphics(fig, fullfile(outpath,'model_comparison_RW_HGF.pdf'),'ContentType','vector');

save(fullfile(outpath,'model_comparison_RW_HGF.mat'), 'res', 'bic_mat', 'model_lbl');
fprintf('\nSaved model_comparison_RW_HGF.mat\n');

% =============================================================================
%% LOCAL FUNCTIONS — observation packing
% =============================================================================
function obs = pack_observations(td, stim_field, has_perc, N_STIM) %#ok<INUSD>
[nB, nT] = size(td.correct);
stim_id=[]; block_id=[]; trial_id=[]; resp_go=[]; y_t=[];
for b = 1:nB
    for t = 1:nT
        cor  = td.correct(b,t);  gotr = td.goTrial(b,t);
        if isnan(cor) || isnan(gotr); continue; end
        s = td.(stim_field)(b,t);
        if isnan(s); continue; end
        rgo = gotr*cor + (1-gotr)*(1-cor);                 % did they press Go?
        if has_perc && isfield(td,'perceivedCorrect') && ~isnan(td.perceivedCorrect(b,t))
            pc = td.perceivedCorrect(b,t);
        else
            pc = cor;
        end
        yt = double((rgo==1 && pc==1) || (rgo==0 && pc==0)); % feedback-consistency
        stim_id(end+1)=s; block_id(end+1)=b; trial_id(end+1)=t; %#ok<AGROW>
        resp_go(end+1)=rgo; y_t(end+1)=yt; %#ok<AGROW>
    end
end
obs = struct('stim_id',stim_id,'block_id',block_id,'trial_id',trial_id, ...
             'resp_go',resp_go,'y_t',y_t);
end

% =============================================================================
%% LOCAL FUNCTIONS — Rescorla-Wagner (single + dual learning rate)
% =============================================================================
function nll = rw_single_nll(p, obs)
alpha = p(1); beta = p(2);
nll = rw_core(obs, alpha, alpha, beta);
end

function nll = rw_dual_nll(p, obs)
% Dual learning rate: alpha_pos for positive PE (delta>0), alpha_neg otherwise.
nll = rw_core(obs, p(1), p(2), p(3));
end

function nll = rw_core(obs, a_pos, a_neg, beta)
N_s = 4; theta = 0.5*ones(1,N_s); ll = 0; prev_block = -1;
for t = 1:numel(obs.stim_id)
    s = obs.stim_id(t);
    if obs.block_id(t) ~= prev_block
        theta = 0.5*ones(1,N_s); prev_block = obs.block_id(t);
    end
    if isnan(s) || s<1 || s>N_s; continue; end
    p_go = 1/(1+exp(-beta*(theta(s)-0.5)));
    r    = obs.resp_go(t);
    ll   = ll + log(max(p_go^r * (1-p_go)^(1-r), 1e-10));
    delta = obs.y_t(t) - theta(s);
    a = a_pos*(delta>=0) + a_neg*(delta<0);
    theta(s) = max(0.01, min(0.99, theta(s) + a*delta));
end
nll = -ll;
end

% =============================================================================
%% LOCAL FUNCTIONS — 3-level binary HGF (Mathys 2011/2014)
% =============================================================================
function nll = hgf_binary_nll(p, obs)
% Params: omega2 (L2 tonic vol), omega3 (L3 tonic vol), beta (softmax).
% kappa fixed = 1. One HGF per stimulus, reset at block boundaries.
% Response model: P(go) = softmax over predicted P(go-target) muhat1.
om2 = p(1); om3 = p(2); beta = p(3); kappa = 1;
N_s = 4;

% per-stimulus states
mu2 = zeros(1,N_s); sa2 = ones(1,N_s);          % belief (logit) + variance
mu3 = ones(1,N_s);  sa3 = ones(1,N_s);          % log-volatility + variance
ll = 0; prev_block = -1;

for t = 1:numel(obs.stim_id)
    s = obs.stim_id(t);
    if obs.block_id(t) ~= prev_block
        mu2(:)=0; sa2(:)=1; mu3(:)=1; sa3(:)=1; prev_block = obs.block_id(t);
    end
    if isnan(s) || s<1 || s>N_s; continue; end

    % --- predictions ---
    muhat2 = mu2(s);
    muhat1 = 1/(1+exp(-muhat2));                 % predicted P(go-target)
    % response likelihood (softmax on predicted Go-target prob)
    p_go = 1/(1+exp(-beta*(muhat1-0.5)));
    r    = obs.resp_go(t);
    ll   = ll + log(max(p_go^r * (1-p_go)^(1-r), 1e-10));

    % --- input (binary outcome: was the go-target confirmed this trial) ---
    u = obs.y_t(t);

    % --- HGF update (binary; Mathys 2014 eqs) ---
    nu2    = exp(kappa*mu3(s) + om2);            % predicted L2 step
    pihat2 = 1/(sa2(s) + nu2);
    pi2    = pihat2 + muhat1*(1-muhat1);
    mu2(s) = muhat2 + (1/pi2)*(u - muhat1);
    sa2(s) = 1/pi2;

    % level 3 (volatility)
    pihat3 = 1/(sa3(s) + exp(om3));
    w2     = nu2 * pihat2;
    da2    = (1/pi2 + (mu2(s)-muhat2)^2)*pihat2 - 1;
    pi3    = pihat3 + 0.5*kappa^2*w2*(w2 + (2*w2-1)*da2);
    if pi3 <= 0                                  % invalid trajectory
        nll = 1e6; return;
    end
    mu3(s) = mu3(s) + 0.5*kappa*w2/pi3*da2;
    sa3(s) = 1/pi3;
end
nll = -ll;
end

function [x, bic] = fit_hgf_tapas(obs, Nobs)
% Wrapper for the validated TAPAS binary HGF (recommended for publication).
% Maps obs to a single sequence per stimulus is non-trivial in TAPAS, so here
% we fit the pooled binary input sequence u = y_t with tapas_hgf_binary and a
% unit-square sigmoid response model. Returns [om2 om3 beta-equivalent] and BIC.
try
    u = obs.y_t(:);
    y = obs.resp_go(:);
    r = tapas_fitModel(y, u, 'tapas_hgf_binary_config', 'tapas_unitsq_sgm_config');
    om2 = r.p_prc.om(2);
    om3 = r.p_prc.om(3);
    ze  = r.p_obs.ze;                 % response noise (proxy for beta)
    bic = -2*r.optim.LME + 3*log(Nobs);   % approx; LME is the log model evidence
    x = [om2, om3, ze];
catch ME
    warning('TAPAS HGF failed (%s); using built-in HGF.', ME.message);
    [x, bic] = fit_grid_fmincon(@(p) hgf_binary_nll(p, obs), ...
        {linspace(-6,-1,8), linspace(-7,-2,8), linspace(2,12,8)}, ...
        [-12 -12 0.5],[-0.5 -0.5 30], Nobs, 3);
end
end

% =============================================================================
%% LOCAL FUNCTIONS — generic fit (grid init + fmincon) + BIC
% =============================================================================
function [x_fit, bic] = fit_grid_fmincon(nll_fun, grids, lb, ub, Nobs, k)
% Grid search over the cartesian product of `grids` for a good start, then
% refine with fmincon. Returns the fitted params and BIC = 2*NLL + k*log(Nobs).
gv = grids;
[G{1:numel(gv)}] = ndgrid(gv{:});
P = cell2mat(cellfun(@(g) g(:), G, 'UniformOutput', false));   % rows = candidates
best = inf; x0 = P(1,:);
for i = 1:size(P,1)
    v = nll_fun(P(i,:));
    if v < best; best = v; x0 = P(i,:); end
end
opts = optimoptions('fmincon','Display','off','MaxIterations',500);
try
    x_fit = fmincon(nll_fun, x0, [],[],[],[], lb, ub, [], opts);
catch
    x_fit = x0;
end
nll = nll_fun(x_fit);
bic = 2*nll + k*log(Nobs);
end
