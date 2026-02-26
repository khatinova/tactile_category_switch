load("\\humerus\pharm_banerjee\data\Projects\EEG_projects\Salient_Modality_Switch\Salient mod switch KH\Results\Behav results\all_trial_data_v2.mat")

T = group_T;
T_confidence = T(T.subjID ~= "Ox11",:);
T_tactile = T(T.block_type ~= "V",:);


summary_T = groupsummary(T, 'subjID', 'mean', {'correct', 'confidence', 'perceivedCorrect', 'stay_choice'});


formula1 = 'correct ~ block_type + trial  + (1 + trial|subjID)';
formula2 = 'correct ~ block_type + trial  + block + stim_config + (1 + trial|subjID)';
formula3 = 'correct ~ block_type + trial  + block + stim_config + trueFB + (1 + trial|subjID)';
formula4 = 'correct ~ transition + block_type + trial  + block + stim_config + trueFB + (1 + trial|subjID)';
formula5 = 'confidence ~ block_type + trial + rev_state + trueFB + (1 + trial|subjID)';

glme1 = fitglme(group_T, formula1, 'Distribution','Binomial','Link','logit');
glme2 = fitglme(group_T, formula2, 'Distribution','Binomial','Link','logit');
glme3 = fitglme(group_T, formula3, 'Distribution','Binomial','Link','logit');
glme4 = fitglme(group_T, 'correct ~ transition + block_type + trial  + block + stim_config + trueFB + (1 + trial|subjID)' ...
    , 'Distribution','Binomial','Link','logit')
lme1 = fitlme(T_confidence,'confidence ~ block_type + prev_block_type + trial + rev_state + trueFB + (1 + trial|subjID)')


% probability of repeating the same choice = win-stay? Determined by
% prev-correct positive effect, interaction with prevTrueFB is true/falseFB
glme1 = fitglme(group_T,'stay_choice ~ prevCorrect + prevTrueFB + (1|subjID)', ...
     'Distribution','Binomial','Link','logit')
glme2 = fitglme(group_T,'stay_choice ~ prevCorrect*prevTrueFB + (1|subjID)', ...
     'Distribution','Binomial','Link','logit')
glme3 = fitglme(group_T,'stay_choice ~ prevCorrect*prevTrueFB + (trial|subjID)', ...
     'Distribution','Binomial','Link','logit')
glme4 = fitglme(group_T,'stay_choice ~ prevChoice +prevCorrect*prevTrueFB + block_type + (1 + trial|subjID)', ...
     'Distribution','Binomial','Link','logit')
glme4 = fitglme(group_T,'stay_choice ~ prevCorrect*prevTrueFB + prevChoice + block_type + trial + (1 + trial|subjID)', ...
     'Distribution','Binomial','Link','logit')

glme4 = fitglme(T_tactile,'stay_choice ~ prevCorrect*prevTrueFB*prev_block_type + prevChoice + trial + stimID + (1 + trial|subjID)', ...
     'Distribution','Binomial','Link','logit')

% Results:
% 
% Model information:
%     Number of observations            9800
%     Fixed effects coefficients           5
%     Random effects coefficients         48
%     Covariance parameters                4
% 
% Formula:
%     confidence ~ 1 + trial + rev_state + block_type + trueFB + (1 + trial | subjID)
% 
% Model fit statistics:
%     AIC      BIC      LogLikelihood    Deviance
%     42833    42898    -21408           42815   
% 
% Fixed effects coefficients (95% CIs):
%     Name                    Estimate    SE           tStat     DF      pValue         Lower        Upper   
%     {'(Intercept)' }          6.0152      0.31352    19.186    9795     1.4297e-80       5.4007      6.6298
%     {'trial'       }        0.013507    0.0025647    5.2667    9795     1.4182e-07    0.0084802    0.018535
%     {'rev_state'   }        -0.37363     0.075497    -4.949    9795     7.5851e-07     -0.52162    -0.22564
%     {'block_type_P'}         -1.1334       0.0518    -21.88    9795    1.1971e-103      -1.2349     -1.0318
%     {'trueFB_1'    }         0.15945     0.070512    2.2614    9795       0.023759     0.021235     0.29767
% 
% Random effects covariance parameters (95% CIs):
% Group: subjID (24 Levels)
%     Name1                  Name2                  Type            Estimate    Lower        Upper   
%     {'(Intercept)'}        {'(Intercept)'}        {'std' }           1.484       1.1119      1.9805
%     {'trial'      }        {'(Intercept)'}        {'corr'}         -0.4807     -0.73856    -0.10015
%     {'trial'      }        {'trial'      }        {'std' }        0.010843    0.0079096    0.014864
% 
% Group: Error
%     Name               Estimate    Lower     Upper 
%     {'Res Std'}        2.131       2.1013    2.1611


compare(glme2, glme1)
    % Theoretical Likelihood Ratio Test
    % 
    % Model    DF    AIC      BIC      LogLik    LRStat    deltaDF    pValue
    % glme2    13    56312    56409    -28143                               
    % glme1     7    56266    56318    -28126    34.311    -6         NaN 

compare(glme3, glme1)

    % Theoretical Likelihood Ratio Test
    % 
    % Model    DF    AIC      BIC      LogLik    LRStat    deltaDF    pValue
    % glme3    14    56320    56424    -28146                               
    % glme1     7    56266    56318    -28126    40.01     -7         NaN 



% RW fitted param models:
fitlme(group_T, 'confidence ~ rev_state * alpha_pre + (1|subjID)')

% Model fit statistics:
%     AIC      BIC      LogLikelihood    Deviance
%     48821    48865    -24405           48809   
% 
% Fixed effects coefficients (95% CIs):
%     Name                           Estimate    SE          tStat      DF       pValue        Lower        Upper  
%     {'(Intercept)'        }          5.7608     0.38362     15.017    10896    1.8165e-50       5.0088     6.5127
%     {'rev_state'          }        0.090927    0.062218     1.4614    10896       0.14393    -0.031032    0.21289
%     {'alpha_pre'          }          4.2783      4.6256    0.92492    10896       0.35503      -4.7887     13.345
%     {'rev_state:alpha_pre'}          1.7546     0.77474     2.2648    10896      0.023547      0.23597     3.2732
% 
% Random effects covariance parameters (95% CIs):
% Group: subjID (25 Levels)
%     Name1                  Name2                  Type           Estimate    Lower      Upper 
%     {'(Intercept)'}        {'(Intercept)'}        {'std'}        1.306       0.98796    1.7263
% 
% Group: Error
%     Name               Estimate    Lower     Upper 
%     {'Res Std'}        2.2576      2.2278    2.2878

fitlme(group_T,'alpha_post ~ prev_block_type + (1|subjID)')


