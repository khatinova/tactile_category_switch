function P = kh_poster_palette()
% KH_POSTER_PALETTE  Single source of truth for figure colours across the
% pipeline (behaviour S6 + EEG poster P8). Keeping every variable's colour in
% one place means the poster reads as one coherent story: the same blue always
% means "deterministic / certain", the same orange-red always means
% "probabilistic / uncertain", etc.
%
%   P = kh_poster_palette() returns a struct of RGB triplets (0-1):
%
%   BLOCK TYPE (certainty)
%     P.D      deterministic / certain   (blue)
%     P.P      probabilistic / uncertain (orange-red)
%
%   OUTCOME
%     P.correct, P.incorrect            (green / red)
%
%   TASK STAGE (relative to reversal): LN, LE, RN, RE
%     P.stage.LN/LE/RN/RE  and  P.stage_order = {'LN','LE','RN','RE'}
%     P.stage_list = [4x3] in LN,LE,RN,RE order
%
%   BLOCK TRANSITION (prior -> current uncertainty)
%     P.trans.DD / DP / PD / PP
%     P.trans_order = {'D->D','D->P','P->D','P->P'}
%     P.trans_list  = [4x3]
%
%   CONNECTIVITY PATHWAYS
%     P.fp   fronto-parietal  (teal)
%     P.fs   fronto-sensory   (purple)
%
%   NEUTRALS
%     P.grey, P.lightgrey, P.k
%
% NOTE: the block/stage/transition RGBs match S6_behaviour_plots_stats.m so the
% behaviour and EEG figures are colour-consistent on the poster.

% --- Block type (certainty) -------------------------------------------------
P.D = [0.15 0.45 0.70];     % deterministic / certain  (blue)
P.P = [0.80 0.30 0.10];     % probabilistic / uncertain (orange-red)

% --- Outcome ----------------------------------------------------------------
P.correct   = [0.20 0.60 0.20];   % green
P.incorrect = [0.80 0.20 0.20];   % red

% --- Task stage -------------------------------------------------------------
P.stage.LN = [0.20 0.63 0.17];    % green  (learning naive)
P.stage.LE = [0.12 0.47 0.71];    % blue   (learning expert)
P.stage.RN = [0.85 0.33 0.10];    % orange (reversal naive)
P.stage.RE = [0.58 0.40 0.74];    % purple (reversal expert)
P.stage_order = {'LN','LE','RN','RE'};
P.stage_list  = [P.stage.LN; P.stage.LE; P.stage.RN; P.stage.RE];

% --- Block transition (prior -> current) ------------------------------------
P.trans.DD = [0.12 0.47 0.71];
P.trans.DP = [0.85 0.33 0.10];
P.trans.PD = [0.47 0.67 0.19];
P.trans.PP = [0.80 0.20 0.60];
P.trans_order = {'D->D','D->P','P->D','P->P'};
P.trans_list  = [P.trans.DD; P.trans.DP; P.trans.PD; P.trans.PP];

% --- Connectivity pathways --------------------------------------------------
P.fp = [0.00 0.55 0.55];    % fronto-parietal  (teal)
P.fs = [0.55 0.25 0.60];    % fronto-sensory   (purple)

% --- Neutrals ---------------------------------------------------------------
P.grey      = [0.40 0.40 0.40];
P.lightgrey = [0.80 0.80 0.80];
P.k         = [0.00 0.00 0.00];

end
