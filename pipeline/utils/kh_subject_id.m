function out = kh_subject_id(action, varargin)
% KH_SUBJECT_ID  Unified subject-identifier conventions for the whole pipeline.
%
% THE CANONICAL CONVENTION (use everywhere downstream):
%   subj_id : zero-padded STRING label, cohort-prefixed, e.g. "Ox03", "Nc07"
%             - KH cohort  -> "Ox%02d"
%             - RR cohort  -> "Nc%02d"
%   subj    : numeric participant number (double), e.g. 3, 7
%   cohort  : "KH" or "RR" (categorical/string)
%
% This single function replaces the many ad-hoc conventions found across the
% original scripts (subjID / subj_id / subj / subject; strcat vs compose;
% leading-zero loss from string(double)). Always build labels with compose so
% leading zeros are preserved ("Ox03", NOT "Ox3").
%
% USAGE
% -----
%   id   = kh_subject_id('make', 3, 'KH')        -> "Ox03"
%   id   = kh_subject_id('make', 7, 'RR')        -> "Nc07"
%   n    = kh_subject_id('num',  "Ox03")         -> 3
%   coh  = kh_subject_id('cohort', "Nc07")       -> "RR"
%   T    = kh_subject_id('standardise', T)       -> table with canonical
%                                                   subj_id / subj / cohort
%
% The 'standardise' action repairs a table that may use any of the legacy
% column names (subjID, subj_id, subject, subj) and guarantees the canonical
% trio exists and is consistent.

switch lower(action)

    case 'make'
        % kh_subject_id('make', num, cohort)
        num    = varargin{1};
        cohort = string(varargin{2});
        prefix = local_prefix(cohort);
        out    = prefix + compose("%02d", double(num));

    case 'num'
        % kh_subject_id('num', id)  -> numeric part
        id  = string(varargin{1});
        tok = regexp(id, '\d+', 'match', 'once');
        out = str2double(tok);

    case 'cohort'
        % kh_subject_id('cohort', id)  -> "KH" / "RR"
        id = string(varargin{1});
        if startsWith(id, "Ox")
            out = "KH";
        elseif startsWith(id, "Nc")
            out = "RR";
        else
            out = "unknown";
        end

    case 'standardise'
        out = local_standardise_table(varargin{1});

    otherwise
        error('kh_subject_id: unknown action "%s".', action);
end
end

% -------------------------------------------------------------------------
function prefix = local_prefix(cohort)
switch upper(string(cohort))
    case "KH";  prefix = "Ox";
    case "RR";  prefix = "Nc";
    otherwise;  error('kh_subject_id: unknown cohort "%s" (expected KH or RR).', cohort);
end
end

% -------------------------------------------------------------------------
function T = local_standardise_table(T)
% Guarantee canonical subj_id (string Ox##/Nc##), subj (numeric), cohort.

vn = T.Properties.VariableNames;

% 1) Resolve a source string id from any legacy name.
src_id = strings(height(T), 1);
if ismember('subj_id', vn)
    src_id = string(T.subj_id);
elseif ismember('subjID', vn)
    src_id = string(T.subjID);
elseif ismember('subject', vn)
    src_id = string(T.subject);
end

% 2) If we still have no usable label but a numeric subj + cohort exist, build it.
need_build = all(src_id == "" | ismissing(src_id));
if need_build && ismember('subj', vn)
    if ismember('cohort', vn)
        coh = string(T.cohort);
    elseif ismember('researcher', vn)
        coh = string(T.researcher);
    else
        coh = repmat("KH", height(T), 1);   % default; adjust if needed
    end
    for i = 1:height(T)
        src_id(i) = kh_subject_id('make', double(T.subj(i)), coh(i));
    end
end

% 3) Write canonical columns.
T.subj_id = src_id;
T.subj    = arrayfun(@(s) kh_subject_id('num', s),    src_id);
T.cohort  = arrayfun(@(s) kh_subject_id('cohort', s), src_id);

% 4) Keep a copy under subjID for backward compatibility with old scripts.
T.subjID = src_id;
end
