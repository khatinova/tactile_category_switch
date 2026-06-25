% =============================================================================
% EEG AMPLITUDE NORMALISATION + FEATURE EXTRACTION + RQ ANALYSIS
%
% This version loads each outcome .set file ONCE per subject and reuses the
% same data for baseline RMS, single-trial features, and ERP plotting.
%
% It extracts BOTH:
%   1) FRN  = negative-going frontocentral outcome response
%   2) RewP = positive-going reward positivity / reward-related difference
%
% IMPORTANT INTERPRETATION:
%   - FRN is quantified here as a fixed-window mean amplitude.
%   - RewP is quantified here as a fixed-window mean amplitude.
%   - Peak-picking is intentionally avoided because it is unstable when the
%     waveform is noisy or flat.
%
% -------------------------------------------------------------------------
% clear; close all; clc;
close all
% -------------------------------------------------------------------------
%% PATHS
% -------------------------------------------------------------------------
remote = 0;
if remote == 1
    base_path = '/Volumes/PHARM_BANERJEE/data/Projects/EEG_projects/Salient_Modality_Switch';
else

    base_path = '\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch';
end

epoch_file_folder = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis', 'Epoched_data');
RR_epoch_folder   = fullfile(base_path, 'Salient mod switch RR', 'Results', 'EEG analysis', 'Epoched_data');
figure_output_folder = fullfile(base_path, 'Salient mod switch KH', 'Results', 'EEG analysis', 'Figures', 'RQ_analysis');
% if ~exist(figure_output_folder,'dir'), mkdir(figure_output_folder); end


load(fullfile(epoch_file_folder,'group_stage_table_mICA.mat'), 'group_table');

% -------------------------------------------------------------------------
%% SETTINGS
% -------------------------------------------------------------------------
valid_participants = [3:12, 14:23,27];
ERP_plot_window = [-200 800];
rm_baseline     = [-500 0];

% Outcome windows (ms)
FRN_win   = [120 350];
RewP_win  = [200 400];
P300_win  = [300 600];
Theta_win = [200 500];
PLV_win   = [200 400];

% Frontocentral ROI for outcome evaluation
frontocentral_channels = {'FCz','Cz'};

% Other channel groups
par_channels = {'Pz','P1','P2'};
acc_channels = {'FCz','Fz','AFz'};
som_channels = {'C3','C4','CP3','CP1','C5','CP5'};

MIN_TRIALS_PLV = 5;
MIN_MATCH_RATE = 0.70;

n_rows = height(group_table);

% Feature columns
group_table.FCz_neg_peak_amp   = nan(n_rows,1);   % diagnostic only, not FRN
group_table.FCz_neg_peak_norm  = nan(n_rows,1);
group_table.FRN_amp            = nan(n_rows,1);
group_table.FRN_norm           = nan(n_rows,1);
group_table.RewP_amp           = nan(n_rows,1);
group_table.RewP_norm          = nan(n_rows,1);
group_table.P300_amp           = nan(n_rows,1);
group_table.P300_norm          = nan(n_rows,1);
group_table.Theta_amp          = nan(n_rows,1);
group_table.PLV_fp             = nan(n_rows,1);
group_table.PLV_fs             = nan(n_rows,1);
group_table.PLV_fp_pairwise    = nan(n_rows,1);
group_table.PLV_fs_pairwise    = nan(n_rows,1);
group_table.FCzCz_signal        = cell(n_rows,1);
group_table.P300_signal         = cell(n_rows,1);
group_table.Theta_signal        = cell(n_rows,1);
group_table.FCz_neg_peak_lat  = nan(n_rows,1);
group_table.P300_peak_lat     = nan(n_rows,1);

subj_list = unique(group_table.subj);
subj_baseline_rms = containers.Map('KeyType','double','ValueType','double');

% Storage for the normalised ERP plot, so no subject files need reloading later
raw_correct = [];
raw_incorrect = [];
norm_correct = [];
norm_incorrect = [];
t_global = [];

% -------------------------------------------------------------------------
%% SINGLE MASTER LOOP: load each subject once
% -------------------------------------------------------------------------
fprintf('Processing subjects in one pass...\n');

for si = valid_participants
    
    subj_num = double(si);
    subj     = sprintf('Ox%02d', subj_num);
    subj_rows = group_table.subj == categorical(si);

    fprintf('  %s (%d trials)...\n', subj, sum(subj_rows));

    fname = sprintf('%s_outcome_trimmed.set', subj);
    fpath = fullfile(epoch_file_folder, fname);
    if ~exist(fpath,'file')
        fname = sprintf('%s_outcome.set', subj);
        fpath = fullfile(epoch_file_folder, fname);
    end

    % Load outcome .set ONCE
    EEGp = pop_loadset(fname, epoch_file_folder);
    if isempty(t_global)
        t_global = EEGp.times;
    end

    % Optional theta/phase files
    theta_fname = sprintf('%s_outcome_theta_trimmed.set', subj);
    phase_fname = sprintf('%s_outcome_phase_trimmed.set', subj);
    has_theta = exist(fullfile(epoch_file_folder, theta_fname),'file');
    has_phase = exist(fullfile(epoch_file_folder, phase_fname),'file');

    theta_fname_untrimmed = sprintf('%s_outcome_theta.set', subj);
    phase_fname_untrimmed = sprintf('%s_outcome_phase.set', subj);

    if has_theta
        EEGp_theta = pop_loadset(theta_fname, epoch_file_folder);
    else
        EEGp_theta = pop_loadset(theta_fname_untrimmed, epoch_file_folder);
    end
    if has_phase
        EEGp_phase = pop_loadset(phase_fname, epoch_file_folder);
    else
        EEGp_phase = pop_loadset(theta_fname_untrimmed, epoch_file_folder);
    end

    % Channel indices
    fcz_idx = find(strcmpi({EEGp.chanlocs.labels}, 'FCz'), 1);
    cz_idx  = find(strcmpi({EEGp.chanlocs.labels}, 'Cz'), 1);
    par_idx = find(ismember(lower({EEGp.chanlocs.labels}), lower(par_channels)));
    acc_idx = find(ismember(lower({EEGp.chanlocs.labels}), lower(acc_channels)));
    som_idx = find(ismember(lower({EEGp.chanlocs.labels}), lower(som_channels)));

    bl_mask   = EEGp.times >= rm_baseline(1) & EEGp.times <= rm_baseline(2);
    frn_mask  = EEGp.times >= FRN_win(1)  & EEGp.times <= FRN_win(2);
    rewp_mask = EEGp.times >= RewP_win(1) & EEGp.times <= RewP_win(2);
    p300_mask = EEGp.times >= P300_win(1) & EEGp.times <= P300_win(2);
    th_mask   = EEGp.times >= Theta_win(1) & EEGp.times <= Theta_win(2);
    plv_mask  = EEGp.times >= PLV_win(1)   & EEGp.times <= PLV_win(2);
    plv_bl    = EEGp.times >= -200 & EEGp.times <= -20;

    % Baseline RMS for this subject from the already-loaded EEGp
    if ~isempty(fcz_idx)
        bl_data = squeeze(EEGp.data(fcz_idx, bl_mask, :));
        bline_rms = rms(bl_data(:));
        subj_baseline_rms(subj_num) = bline_rms;
        fprintf('    baseline RMS = %.3f uV\n', bline_rms);
    else
        warning('%s: FCz not found, skipping baseline-dependent measures.', subj);
        bline_rms = NaN;
    end

    % Trial-wise feature extraction
    row_indices = find(subj_rows);

    for ri = 1:numel(row_indices)
        r  = row_indices(ri);
        ep = group_table.epoch(r);
        if isnan(ep) || ep < 1 || ep > EEGp.trials
            continue
        end

        % ----- Diagnostic FCz negative peak, not called FRN -----
        if ~isempty(fcz_idx) && ~isnan(bline_rms)
            sig_fcz = double(EEGp.data(fcz_idx, :, ep));
            sig_fcz = sig_fcz - mean(sig_fcz(bl_mask), 'omitnan');
            group_table.FCz_neg_peak_amp(r)  = min(sig_fcz(frn_mask));
            group_table.FCz_neg_peak_norm(r) = group_table.FCz_neg_peak_amp(r) / bline_rms;
        end

        if ~isempty(fcz_idx) && ~isempty(cz_idx) && ~isnan(bline_rms)
            sig_frn = mean(double(EEGp.data([fcz_idx cz_idx], :, ep)), 1, 'omitnan');
            sig_frn = sig_frn - mean(sig_frn(bl_mask), 'omitnan');

            % FRN_norm and RewP_norm are not microvolts anymore, because
            % they are baseline corrected!!!
        
            group_table.FCzCz_signal{r} = sig_frn;
            group_table.FRN_amp(r)  = mean(sig_frn(frn_mask), 'omitnan');
            group_table.FRN_norm(r) = group_table.FRN_amp(r) / bline_rms;
            group_table.RewP_amp(r)  = mean(sig_frn(rewp_mask), 'omitnan');
            group_table.RewP_norm(r) = group_table.RewP_amp(r) / bline_rms;
        end

        % ----- Diagnostic FCz negative peak, not called FRN -----
        if ~isempty(fcz_idx) && ~isnan(bline_rms)
            sig_fcz = double(EEGp.data(fcz_idx, :, ep));
            sig_fcz = sig_fcz - mean(sig_fcz(bl_mask), 'omitnan');
        
            win_vals = sig_fcz(frn_mask);
            win_t    = EEGp.times(frn_mask);
        
            [pk_amp, pk_ix] = min(win_vals, [], 'omitnan');
            pk_time = win_t(pk_ix);
        
            group_table.FCz_neg_peak_amp(r)  = pk_amp;
            group_table.FCz_neg_peak_norm(r) = pk_amp / bline_rms;
            group_table.FCz_neg_peak_lat(r)  = pk_time;

            % % Debug figure
            % if ismember(ep, [1, 15, 30, 45, 60, 78, 141, 206, 255])
            %     figure; 
            %     plot(EEGp.times, sig_fcz, 'k-', 'LineWidth', 1); hold on;
            %     xline(FRN_win(1), 'b--', 'HandleVisibility','off');
            %     xline(FRN_win(2), 'b--', 'HandleVisibility','off');
            %     xline(0, 'k:', 'HandleVisibility','off');
            %     plot(pk_time, pk_amp, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
            %     yline(0, 'k:');
            %     title(sprintf('%s trial %d | FCz negative peak', subj, ep), 'Interpreter','none');
            %     xlabel('Time (ms)'); ylabel('\muV');
            %     legend({'Baseline-corrected FCz','FRN window','Peak'}, 'Location', 'best');
            % 
            % 
            %     fprintf('%s trial %d: FCz neg peak = %.3f uV at %.1f ms\n', ...
            %         subj, ep, pk_amp, pk_time);
            % end
        end
        
        %% rerun to validate poistioning of P300 window!
       if ~isempty(par_idx) && ~isnan(bline_rms)
            sig_par = mean(double(EEGp.data(par_idx, :, ep)), 1, 'omitnan');
            sig_par = sig_par - mean(sig_par(bl_mask), 'omitnan');
        
            win_vals = sig_par(p300_mask);
            win_t    = EEGp.times(p300_mask);
        
            [pk_amp, pk_ix] = max(win_vals, [], 'omitnan');
            pk_time = win_t(pk_ix);
        
            group_table.P300_signal{r} = sig_par;
            group_table.P300_amp(r)    = pk_amp;
            group_table.P300_norm(r)   = pk_amp / bline_rms;
            group_table.P300_peak_lat(r) = pk_time;
        
            % Debug plot for selected epochs
            % if ismember(ep, [1, 15, 30, 45, 60, 78, 141, 206, 255])
            %     figure; 
            %     plot(EEGp.times, sig_par, 'k-', 'LineWidth', 1); hold on;
            %     xline(P300_win(1), 'b--', 'HandleVisibility','off');
            %     xline(P300_win(2), 'b--', 'HandleVisibility','off');
            %     xline(0, 'k:', 'HandleVisibility','off');
            %     plot(pk_time, pk_amp, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
            %     yline(0, 'k:');
            % 
            %     title(sprintf('%s trial %d | P300 peak', subj, ep), 'Interpreter','none');
            %     xlabel('Time (ms)');
            %     ylabel('\muV');
            %     legend({'Baseline-corrected parietal trace','P300 window','Peak'}, 'Location', 'best');
            %     drawnow;
            % 
            %     fprintf('%s trial %d: P300 peak = %.3f uV at %.1f ms\n', ...
            %         subj, ep, pk_amp, pk_time);
            % end
        end
        
        if ~isempty(EEGp_theta) && ~isempty(acc_idx)
            sig_th = mean(double(EEGp_theta.data(acc_idx, :, ep)), 1, 'omitnan');
            env    = abs(hilbert(sig_th));
            env    = env - mean(env(bl_mask), 'omitnan');
        
            group_table.Theta_signal{r} = env;
            group_table.Theta_amp(r) = mean(env(th_mask), 'omitnan');
        end
    end

    % ----- PLV computed per condition -----
    if ~isempty(EEGp_phase) && ~isempty(acc_idx)
        cond_combos = unique(group_table(subj_rows, {'block_type','stage','correct','false_fb'}), 'rows');
        for c = 1:height(cond_combos)
            bt  = cond_combos.block_type(c);
            st  = cond_combos.stage(c);
            cor = cond_combos.correct(c);
            ff  = cond_combos.false_fb(c);

            cond_mask = subj_rows & ...
                group_table.block_type == bt & group_table.stage == st & ...
                group_table.correct == cor & group_table.false_fb == ff;
            ep_cond = group_table.epoch(cond_mask);
            ep_cond = ep_cond(~isnan(ep_cond) & ep_cond >= 1 & ep_cond <= EEGp_phase.trials);

            if numel(ep_cond) < MIN_TRIALS_PLV; continue; end

            if ~isempty(par_idx)
                plv_fp = cross_trial_plv_scalar(EEGp_phase, acc_idx, par_idx, ep_cond, plv_mask, plv_bl);
                plv_fp_pair = cross_trial_plv_pairwise(EEGp_phase, acc_idx, par_idx, ep_cond, plv_mask, plv_bl);
                group_table.PLV_fp(cond_mask & ismember(group_table.epoch, ep_cond)) = plv_fp;
                group_table.PLV_fp_pairwise(cond_mask & ismember(group_table.epoch, ep_cond)) = plv_fp_pair;
            end
            if ~isempty(som_idx)
                plv_fs = cross_trial_plv_scalar(EEGp_phase, acc_idx, som_idx, ep_cond, plv_mask, plv_bl);
                plv_fs_pair = cross_trial_plv_pairwise(EEGp_phase, acc_idx, som_idx, ep_cond, plv_mask, plv_bl);
                group_table.PLV_fs(cond_mask & ismember(group_table.epoch, ep_cond)) = plv_fs;
                group_table.PLV_fs_pairwise(cond_mask & ismember(group_table.epoch, ep_cond)) = plv_fs_pair;
            end
        end
    end

    % Store ERP traces for plotting without reloading files later
    if ~isempty(fcz_idx) && ~isnan(bline_rms)
        subj_rows_plot = group_table(group_table.subj == categorical(si) & ~group_table.false_fb, :);
        ep_c = clean_ep(subj_rows_plot.epoch(subj_rows_plot.correct==1), EEGp);
        ep_i = clean_ep(subj_rows_plot.epoch(subj_rows_plot.correct==0), EEGp);

        erp_c = bl_correct_erp(EEGp, fcz_idx, ep_c, bl_mask);
        erp_i = bl_correct_erp(EEGp, fcz_idx, ep_i, bl_mask);

        if ~isempty(erp_c)
            raw_correct(end+1,:)  = erp_c;
            norm_correct(end+1,:) = erp_c / bline_rms;
        end
        if ~isempty(erp_i)
            raw_incorrect(end+1,:)  = erp_i;
            norm_incorrect(end+1,:) = erp_i / bline_rms;
        end
    end

    clear EEGp EEGp_theta EEGp_phase
end

% -------------------------------------------------------------------------
%% STEP 2: Z-SCORE FEATURES WITHIN SUBJECT (for LME models)
% -------------------------------------------------------------------------
fprintf('\nZ-scoring features within subject...\n');
features_to_zscore = {'FCz_neg_peak_amp','FRN_amp','RewP_amp','P300_amp','Theta_amp','PLV_fp','PLV_fs','PLV_fp_pairwise','PLV_fs_pairwise'};
for f = 1:numel(features_to_zscore)
    fn = features_to_zscore{f};
    fn_z = [fn '_z'];
    group_table.(fn_z) = nan(n_rows,1);
    for si = 1:numel(subj_list)
        mask = group_table.subj == subj_list(si);
        vals = group_table.(fn)(mask);
        mn = mean(vals,'omitnan');
        sd = std(vals,'omitnan');
        if sd > 0
            group_table.(fn_z)(mask) = (vals - mn) / sd;
        end
    end
end

save(fullfile(epoch_file_folder,'group_stage_table_features.mat'), 'group_table');
fprintf('Saved group_stage_table_features.mat\n');

% -------------------------------------------------------------------------
%% STEP 3: NORMALISED GRAND AVERAGE ERP FIGURE
% -------------------------------------------------------------------------
fprintf('\nPlotting normalised vs raw grand averages...\n');

fig_norm = figure('Position',[50 50 1200 500]);
in_w = t_global >= ERP_plot_window(1) & t_global <= ERP_plot_window(2);
sgtitle('Effect of baseline-RMS normalisation on grand average ERP');

subplot(1,2,1); hold on; title('Raw (uV)');
shaded_erp(t_global(in_w), raw_correct(:,in_w),   [0.1 0.6 0.1], 'Correct');
shaded_erp(t_global(in_w), raw_incorrect(:,in_w), [0.7 0.1 0.1], 'Incorrect');
xline(0,'k--'); yline(0,'k:');
xlabel('Time (ms)'); ylabel('Amplitude (uV)'); legend('Box','off');

subplot(1,2,2); hold on; title('Normalised (signal / baseline RMS)');
for ii = 1:size(norm_correct,1)
    plot(t_global(in_w), norm_correct(ii,in_w), 'Color', [0.75 0.92 0.75], 'LineWidth', 0.5, 'HandleVisibility','off');
end
for ii = 1:size(norm_incorrect,1)
    plot(t_global(in_w), norm_incorrect(ii,in_w), 'Color', [0.96 0.78 0.78], 'LineWidth', 0.5, 'HandleVisibility','off');
end
shaded_erp(t_global(in_w), norm_correct(:,in_w),   [0.1 0.6 0.1], 'Correct');
shaded_erp(t_global(in_w), norm_incorrect(:,in_w), [0.7 0.1 0.1], 'Incorrect');
xline(0,'k--'); yline(0,'k:');
xlabel('Time (ms)'); ylabel('Amplitude (normalised)'); legend('Box','off');

saveas(fig_norm, fullfile(figure_output_folder,'Normalisation_comparison.pdf'));

% =========================================================================
%% LOCAL FUNCTIONS
% =========================================================================

function plv_mean = cross_trial_plv_scalar(EEGp_phase, ref_idx, tgt_idx, ep, plv_mask, bl_mask)
% Legacy ROI-collapsed PLV: average phase within each ROI first.
n_ep = numel(ep);
if n_ep == 0
    plv_mean = NaN; return;
end

if isscalar(ref_idx)
    phi_r = squeeze(double(EEGp_phase.data(ref_idx,:,ep)))';
else
    phi_r = squeeze(angle(mean(exp(1i*double(EEGp_phase.data(ref_idx,:,ep))),1,'omitnan')))';
end
if isscalar(tgt_idx)
    phi_t = squeeze(double(EEGp_phase.data(tgt_idx,:,ep)))';
else
    phi_t = squeeze(angle(mean(exp(1i*double(EEGp_phase.data(tgt_idx,:,ep))),1,'omitnan')))';
end
if isvector(phi_r) && n_ep==1; phi_r=phi_r(:)'; end
if isvector(phi_t) && n_ep==1; phi_t=phi_t(:)'; end

plv_ts  = abs(mean(exp(1i*(phi_r-phi_t)),1,'omitnan'));
plv_bl  = mean(plv_ts(bl_mask),'omitnan');
plv_ts  = plv_ts - plv_bl;
plv_mean = mean(plv_ts(plv_mask),'omitnan');
end

function plv_mean = cross_trial_plv_pairwise(EEGp_phase, ref_idx, tgt_idx, ep, plv_mask, bl_mask)
% Pairwise-averaged PLV: compute every ref x target pair separately.
if isempty(ep)
    plv_mean = NaN; return;
end

ref_idx = ref_idx(:)';
tgt_idx = tgt_idx(:)';
if isempty(ref_idx) || isempty(tgt_idx)
    plv_mean = NaN; return;
end

pair_vals = nan(numel(ref_idx) * numel(tgt_idx), 1);
k = 0;
for r = 1:numel(ref_idx)
    for t = 1:numel(tgt_idx)
        k = k + 1;
        phi_r = squeeze(double(EEGp_phase.data(ref_idx(r),:,ep)))';
        phi_t = squeeze(double(EEGp_phase.data(tgt_idx(t),:,ep)))';
        if isvector(phi_r) && numel(ep)==1; phi_r = phi_r(:)'; end
        if isvector(phi_t) && numel(ep)==1; phi_t = phi_t(:)'; end
        plv_ts = abs(mean(exp(1i*(phi_r - phi_t)),1,'omitnan'));
        plv_bl = mean(plv_ts(bl_mask),'omitnan');
        plv_ts = plv_ts - plv_bl;
        pair_vals(k) = mean(plv_ts(plv_mask),'omitnan');
    end
end

plv_mean = mean(pair_vals,'omitnan');
end

function ep = clean_ep(epoch_col, EEGp)
ep = epoch_col(~isnan(epoch_col) & epoch_col>=1 & epoch_col<=EEGp.trials);
end

function erp = bl_correct_erp(EEGp, ch_idx, ep, bl_mask)
if isempty(ep); erp=[]; return; end
if isscalar(ch_idx)
    raw = squeeze(EEGp.data(ch_idx,:,ep))';
else
    raw = squeeze(mean(EEGp.data(ch_idx,:,ep),1,'omitnan'))';
end
if isvector(raw) && numel(ep)==1; raw=raw(:)'; end
raw = raw - mean(raw(:,bl_mask),2,'omitnan');
erp = mean(raw,1,'omitnan');
end

function shaded_erp(ax_or_t, data_or_ax, clr_or_data, lbl_or_clr, lbl)
% Flexible: shaded_erp(t, data, clr, lbl) called with implicit gca
if nargin == 4
    t = ax_or_t; data = data_or_ax; clr = clr_or_data; lbl = lbl_or_clr;
    ax = gca;
else
    ax = ax_or_t; t = data_or_ax; data = clr_or_data; clr = lbl_or_clr;
end
if isempty(data); return; end
mn = mean(data,1,'omitnan');
se = std(data,0,1,'omitnan')/sqrt(size(data,1));
fill(ax,[t,fliplr(t)],[mn+se,fliplr(mn-se)],clr,'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
plot(ax,t,mn,'Color',clr,'LineWidth',2.5,'DisplayName',sprintf('%s (n=%d)',lbl,size(data,1)));
end
