function stageT = kh_compute_frn_rewp_by_stage(T, times, opts)
% KH_COMPUTE_FRN_REWP_BY_STAGE  Per-condition FRN/RewP difference-wave features.
%
% WHY THIS EXISTS
% ---------------
% The FRN (feedback-related negativity) and RewP (reward positivity) are NOT
% single-trial quantities. They are *difference waves* between correct and
% incorrect outcomes, so they can only be computed once trials are aggregated
% into a condition (here: subject x block_type x stage). This function builds
% that stage-level table. The per-trial table should therefore carry only the
% single-trial frontocentral measures (FCz negative-peak and FCz/Cz mean),
% NOT FRN/RewP.
%
% SIGN CONVENTIONS (field standard; see references below)
% -------------------------------------------------------
%   FRN_amp  = mean_{FRN window} ( ERP_incorrect - ERP_correct )   -> NEGATIVE
%              (loss - win): the relative negativity to worse-than-expected
%              outcomes. More negative = larger FRN.
%   RewP_amp = mean_{RewP window} ( ERP_correct - ERP_incorrect )  -> POSITIVE
%              (win - loss): the relative positivity to rewards/correct.
%              RewP is the same component as FRN with the opposite sign;
%              reduced/absent for losses.
%
%   References:
%     Holroyd & Coles (2002) Psychological Review 109:679-709.
%     Proudfit (2015) Psychophysiology 52:449-459  (the "reward positivity").
%     Sambrook & Goslin (2015) Psychological Bulletin 141:213-235 (meta-analysis;
%        reward-prediction-error scaling; 250-350 ms window).
%
% INPUTS
% ------
%   T      : per-trial table. Required columns:
%              subj_id  (string/categorical)
%              block_type (categorical/string 'D'/'P')
%              stage    (categorical 'LN'/'LE'/'RN'/'RE')
%              correct  (numeric/logical 0/1, OR categorical Incorrect/Correct)
%              <wave_col> : cell column, each cell = 1 x nTime FCz/Cz waveform
%                           (baseline-corrected), aligned to `times`.
%            Optional:
%              false_fb (logical) - if present, only true-feedback trials used.
%   times  : 1 x nTime vector of epoch time points (ms), matching the waveforms.
%   opts   : struct with optional fields:
%              .wave_col   (default 'FCzCz_signal')
%              .FRN_win    (default [250 300] ms)   FRN difference window
%              .RewP_win   (default [250 350] ms)   RewP difference window
%              .stages     (default {'LN','LE','RN','RE'})
%              .block_types(default {'D','P'})
%              .min_trials (default 3) min trials per cell per outcome
%
% OUTPUT
% ------
%   stageT : one row per subj_id x block_type x stage that has BOTH a correct
%            and an incorrect average available. Columns:
%              subj_id, block_type, stage,
%              n_correct, n_incorrect,
%              FRN_amp, RewP_amp,
%              mean_amp_correct_FRNwin, mean_amp_incorrect_FRNwin,
%              diff_wave (cell: incorrect-correct waveform, for plotting)
%
% The windows default to the user's stated choices (FRN 250-300, RewP 250-350).
% The classic literature window is ~250-350 ms for both; inspect the grand
% averages and adjust opts.FRN_win / opts.RewP_win if your peak differs.

% -------------------------------------------------------------------------
if nargin < 3; opts = struct(); end
def = struct('wave_col','FCzCz_signal', ...
             'FRN_win',[250 300], 'RewP_win',[250 350], ...
             'stages',{{'LN','LE','RN','RE'}}, ...
             'block_types',{{'D','P'}}, 'min_trials',3);
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}); opts.(fn{i}) = def.(fn{i}); end
end

times    = times(:)';
frn_mask = times >= opts.FRN_win(1)  & times <= opts.FRN_win(2);
rewp_mask= times >= opts.RewP_win(1) & times <= opts.RewP_win(2);

% Normalise key columns to comparable types
subj_str = string(T.subj_id);
bt_str   = string(T.block_type);
stg_str  = string(T.stage);
corr_num = local_correct_to_numeric(T.correct);

if ismember('false_fb', T.Properties.VariableNames)
    true_fb = ~logical(T.false_fb);
else
    true_fb = true(height(T), 1);
end

waves = T.(opts.wave_col);

subj_list = unique(subj_str, 'stable');
rows = {};

for si = 1:numel(subj_list)
    sn = subj_list(si);
    for bi = 1:numel(opts.block_types)
        bt = string(opts.block_types{bi});
        for sti = 1:numel(opts.stages)
            stg = string(opts.stages{sti});

            base = subj_str == sn & bt_str == bt & stg_str == stg & true_fb;

            erp_c = local_mean_wave(waves, base & corr_num == 1, opts.min_trials);
            erp_i = local_mean_wave(waves, base & corr_num == 0, opts.min_trials);

            if isempty(erp_c) || isempty(erp_i)
                continue;   % need both outcomes for a difference wave
            end

            diff_wave = erp_i - erp_c;          % incorrect - correct

            FRN_amp  = mean(diff_wave(frn_mask),   'omitnan');   % neg (inc-cor)
            RewP_amp = mean(-diff_wave(rewp_mask), 'omitnan');   % pos (cor-inc)

            row = table();
            row.subj_id    = sn;
            row.block_type = categorical(bt, opts.block_types);
            row.stage      = categorical(stg, opts.stages);
            row.n_correct   = sum(base & corr_num == 1);
            row.n_incorrect = sum(base & corr_num == 0);
            row.FRN_amp     = FRN_amp;
            row.RewP_amp    = RewP_amp;
            row.mean_amp_correct_FRNwin   = mean(erp_c(frn_mask), 'omitnan');
            row.mean_amp_incorrect_FRNwin = mean(erp_i(frn_mask), 'omitnan');
            row.diff_wave   = {diff_wave};
            rows{end+1,1}   = row; %#ok<AGROW>
        end
    end
end

if isempty(rows)
    stageT = table();
else
    stageT = vertcat(rows{:});
end
end

% -------------------------------------------------------------------------
function m = local_mean_wave(waves, mask, min_trials)
% Average non-empty waveforms selected by mask; [] if fewer than min_trials.
sel = waves(mask);
sel = sel(~cellfun(@isempty, sel));
if numel(sel) < min_trials
    m = [];
    return;
end
M = cell2mat(cellfun(@(x) x(:)', sel, 'UniformOutput', false));
m = mean(M, 1, 'omitnan');
end

% -------------------------------------------------------------------------
function y = local_correct_to_numeric(v)
if isnumeric(v) || islogical(v)
    y = double(v(:));
    return;
end
sv = lower(strtrim(string(v)));
y  = nan(numel(sv), 1);
y(sv == "1" | sv == "true"  | sv == "correct")   = 1;
y(sv == "0" | sv == "false" | sv == "incorrect") = 0;
end
