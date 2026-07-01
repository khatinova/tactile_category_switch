% quick script to analyse correct for go nogo

subjects = fieldnames(group_results);
nSubj    = numel(subjects);

pGoCorrect    = NaN(nSubj,1);
pNoGoCorrect  = NaN(nSubj,1);

for s = 1:nSubj
    subj   = subjects{s};
    data   = group_results.(subj);
    % Concatenate all blocks
    goFlag = data.goTrial(:);          % 400 × 1  (1 = Go, 0 = NoGo)
    corr   = data.correct(:);          % 400 × 1
    
    pGoCorrect  (s) = mean(corr(goFlag==1) ,'omitnan');   % P(correct|Go)
    pNoGoCorrect(s) = mean(corr(goFlag==0) ,'omitnan');   % P(correct|NoGo)
end

% ---------- group statistics -----------------------------------
meanGo    = mean(pGoCorrect   ,'omitnan');
meanNoGo  = mean(pNoGoCorrect ,'omitnan');

semGo     = std(pGoCorrect   ,'omitnan') ./ sqrt(sum(~isnan(pGoCorrect)));
semNoGo   = std(pNoGoCorrect ,'omitnan') ./ sqrt(sum(~isnan(pNoGoCorrect)));

fprintf('\nGroup mean P(correct):  Go = %.3f ± %.3f  |  NoGo = %.3f ± %.3f\n', ...
        meanGo, semGo, meanNoGo, semNoGo);

%% ------------- optional quick bar‑plot ------------------------
figure; hold on
bh = bar([meanGo meanNoGo],'FaceColor','flat');
bh.CData = [0 0.6 0; 0.8 0 0];    % green / red
errorbar(1, meanGo , semGo , 'k','linestyle','none','capsize',6);
errorbar(2, meanNoGo, semNoGo,'k','linestyle','none','capsize',6);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'Go','NoGo'});
ylabel('P(correct)'); title('Proportion correct  ± SEM');
box off; axis square

validIdx = ~isnan(pGoCorrect) & ~isnan(pNoGoCorrect);
[~, pPair, ~, statsPair] = ttest( ...
        pGoCorrect (validIdx), ...
        pNoGoCorrect(validIdx));

fprintf('Paired t‑test  (Go vs NoGo):  t(%d) = %.3f ,  p = %.4f\n', ...
        statsPair.df, statsPair.tstat, pPair);

%% ---------- annotate significance on the bar‑plot ------------
if exist('bh','var')          % bar‑plot already drawn
    hold on
    yStar = max([meanGo+semGo, meanNoGo+semNoGo]) + 0.03;

    plot([1 2],[yStar yStar],'k','linewidth',1.5);       % horizontal line
    if pPair < .001
        starText = '***';
    elseif pPair < .01
        starText = '**';
    elseif pPair < .05
        starText = '*';
    else
        starText = sprintf('p=%.2f',pPair);
    end
    text(1.5, yStar+0.02, starText,'HorizontalAlignment','center', ...
         'FontSize',12,'FontWeight','bold');
    hold off
end