function x = kh_to_numeric(v)
% KH_TO_NUMERIC  Robustly coerce a table column of any type to a numeric column.
% Handles numeric, logical, categorical, string, char and cell inputs.
% Canonical replacement for the several local copies named local_to_numeric.
if isnumeric(v)
    x = double(v);
elseif islogical(v)
    x = double(v);
elseif iscategorical(v)
    x = str2double(string(v));
elseif isstring(v) || ischar(v)
    x = str2double(string(v));
elseif iscell(v)
    x = str2double(string(v));
else
    error('kh_to_numeric: unsupported column class "%s".', class(v));
end
x = x(:);
end
