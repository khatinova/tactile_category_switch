% basic behav plotting

smoothWin = 5;
preN = 30;
postN = 30;
rel_x = -preN:(postN-1);

subjects = fieldnames(all_trial_data);

for s = 9:15 % numel(subjects)

    subj = subjects{s};
    td = all_trial_data.(subj).trial_data;

    [nBlocks, nTrials] = size(td.correct);
        
    % ---- Figure ----
        figure('Position',[100 100 500 1000]); hold on;

    for b = 1:nBlocks

        % ---- Skip practice block if you do elsewhere ----
        if b == 1 && nBlocks >= 6
            continue
        end

        rev = td.revTrial(b);
        if isnan(rev)
            continue
        end

        % ---- Extract aligned window ----
        abs_idx = rev + rel_x;
        valid = abs_idx >= 1 & abs_idx <= nTrials;

        correct = nan(size(rel_x));
        pCorr   = nan(size(rel_x));
        conf    = nan(size(rel_x));

        correct(valid) = td.correct(b, abs_idx(valid));
        pCorr(valid)   = td.perceivedCorrect(b, abs_idx(valid));
        confidence(valid)    = (td.confidence(b, abs_idx(valid))/10);

        % ---- 5-trial moving average (no cross-reversal leakage) ----
        pre  = correct(1:preN);
        post = correct(preN+1:end);

        c_pre  = confidence(1:preN);
        c_post = confidence(preN+1:end);


        smooth_correct = ...
            [ movmean(pre,  smoothWin, 'omitnan'), ...
              movmean(post, smoothWin, 'omitnan') ];

        smooth_confidence = ...
            [ movmean(c_pre,  smoothWin, 'omitnan'), ...
              movmean(c_post, smoothWin, 'omitnan') ];

        subplot(nBlocks,1,b)
        hold on

        % Smoothed line
        plot(rel_x, smooth_confidence, 'k', 'LineWidth', 2);

        % ---- Scatter raw trials ----
        for i = 1:numel(rel_x)
            if isnan(correct(i)); continue; end

            if pCorr(i) == 1
                col = 'g';
            elseif pCorr(i) == 0
                col = 'r';
            else
                continue
            end

            scatter(rel_x(i), correct(i), 50, ...
                'filled', ...
                'MarkerFaceColor', col, ...
                'MarkerEdgeColor', 'k');
        end

        % ---- Axes & labels ----
        xline(0,'--k','Reversal');
        yline(0.5,'--k');

        hold off
        xlim([rel_x(1) rel_x(end)]);
        ylim([-0.05 1.05]);

        %xlabel('Reversal-aligned trial','FontSize',14);
        ylabel('Confidence','FontSize',14);

        title(sprintf('%s — Block %d', subj, b), ...
              'Interpreter','none','FontSize',16);

        grid on

    end

    save_name = sprintf('%s block confidence.pdf', subj)
    exportgraphics(gcf, save_name, 'ContentType','vector')
end
