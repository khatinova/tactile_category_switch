function apply_fig_style(fig, opts)
% APPLY_FIG_STYLE  Apply consistent, poster/thesis/paper-quality styling to a figure.
%
%   apply_fig_style()            styles the current figure (gcf)
%   apply_fig_style(fig)         styles a specific figure handle
%   apply_fig_style(fig, opts)   styles with overrides (struct)
%
% WHAT IT DOES (applied to every axes in the figure):
%   - Tick direction OUTSIDE the axes (TickDir = 'out')  [requested]
%   - Removes the top/right box so only left + bottom spines remain
%   - Consistent font (Arial) and font sizes
%   - Consistent line width for axes and data lines
%   - Light, even layout suitable for posters, theses and slides
%   - Sets white background and removes the figure menu/tool clutter on export
%
% NOTE ON "EVEN" FIGURES:
%   For panels that must look visually even (e.g. on a poster), call
%   axis(ax,'square') yourself or pass opts.square = true to force every
%   data axes to be square.
%
% COMPANION EXPORT:
%   Use save_fig(fig, fullfile(outdir, name)) below for vector PDF + PNG export
%   with the same styling baked in.
%
% -------------------------------------------------------------------------
% Defaults (override any of these via the opts struct)
% -------------------------------------------------------------------------
if nargin < 1 || isempty(fig); fig = gcf; end
if nargin < 2; opts = struct(); end

def = struct( ...
    'FontName',     'Arial', ...
    'FontSize',     11, ...
    'TitleSize',    12, ...
    'LabelSize',    11, ...
    'AxesLineWidth',1.0, ...
    'LineWidth',    1.8, ...   % data lines (errorbar/plot) bumped to this minimum
    'TickDir',      'out', ...
    'Box',          'off', ...
    'TickLength',   [0.015 0.015], ...
    'square',       false, ...
    'Color',        'w');

fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}); opts.(fn{i}) = def.(fn{i}); end
end

set(fig, 'Color', opts.Color);

% Style every axes in the figure (subplots included), skipping colorbars/legends
ax_all = findall(fig, 'Type', 'axes');
for k = 1:numel(ax_all)
    ax = ax_all(k);

    set(ax, ...
        'TickDir',        opts.TickDir, ...     % ticks on the OUTSIDE
        'Box',            opts.Box, ...         % no top/right spine
        'TickLength',     opts.TickLength, ...
        'LineWidth',      opts.AxesLineWidth, ...
        'FontName',       opts.FontName, ...
        'FontSize',       opts.FontSize, ...
        'Layer',          'top', ...
        'XColor',         [0 0 0], ...
        'YColor',         [0 0 0]);

    % Title / label fonts
    set(ax.Title,  'FontSize', opts.TitleSize, 'FontWeight', 'bold', 'FontName', opts.FontName);
    set(ax.XLabel, 'FontSize', opts.LabelSize, 'FontName', opts.FontName);
    set(ax.YLabel, 'FontSize', opts.LabelSize, 'FontName', opts.FontName);

    % Bump thin data lines up to the minimum line width for visibility
    lines = findobj(ax, 'Type', 'line');
    for li = 1:numel(lines)
        if lines(li).LineWidth < opts.LineWidth && ~strcmp(lines(li).LineStyle, 'none')
            lines(li).LineWidth = opts.LineWidth;
        end
    end

    % Legends inside this axes: no box, matched font
    lg = get(ax, 'Legend');
    if ~isempty(lg)
        set(lg, 'Box', 'off', 'FontName', opts.FontName, 'FontSize', opts.FontSize - 1);
    end

    if opts.square
        axis(ax, 'square');
    end
end
end


function save_fig(fig, path_no_ext, opts)
% SAVE_FIG  Export a styled figure to vector PDF (+ PNG) for publications.
%
%   save_fig(fig, '/path/to/figure_name')        -> figure_name.pdf + .png
%   save_fig(fig, path_no_ext, opts)              with overrides
%
% Uses exportgraphics with vector content for the PDF (crisp at any size) and
% a 300-dpi PNG for quick previews / slides.

if nargin < 3; opts = struct(); end
if ~isfield(opts, 'dpi');  opts.dpi  = 300; end
if ~isfield(opts, 'png');  opts.png  = true; end

apply_fig_style(fig);

[outdir, ~, ~] = fileparts(path_no_ext);
if ~isempty(outdir) && ~exist(outdir, 'dir'); mkdir(outdir); end

exportgraphics(fig, [path_no_ext '.pdf'], 'ContentType', 'vector');
if opts.png
    exportgraphics(fig, [path_no_ext '.png'], 'Resolution', opts.dpi);
end
end
