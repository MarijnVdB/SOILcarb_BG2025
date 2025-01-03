% This script provides an application example
% of the PAWN sensitivity analysis approach (Pianosi and Wagener, 2015)
% 
% REFERENCES
%
% Pianosi, F. and Wagener, T. (2015), A simple and efficient method 
% for global sensitivity analysis based on cumulative distribution 
% functions, Env. Mod. & Soft., 67, 1-11.

% This script prepared by Francesca Pianosi and Fanny Sarrazin
% University of Bristol, 2015
% mail to: francesca.pianosi@bristol.ac.uk

% Adapted by Marijn Van de Broek for Van de Broek et al., A microbially-driven 
% and depth-explicit soil organic carbon model constrained by carbon isotopes 
% to reduce equifinality, Biogeosciences

clc
clear all
close all

%% Add paths

% Add the folder '/safe_R1.1' and all subfolders to the current path, as we need functions stored there
addpath(genpath('../'))

createParams = 0;
evaluateModelResults = 1;

%% Define the parameter ranges

DistrFun  = 'unif'  ; % Parameter distribution

distrpar  = {[-29.8 -28.8]; ...  % d13C of aboveground vegetation
            [-28.3 -27.3]; ...   % d13C of belowground vegetation
            [-29.4 -28.4]; ...   % d13C of root exudates
            [0 0.05]; ...        % Fraction of microbial CO2 uptake
            [0.0108 0.0172]; ... % Effect of the atmospheric CO2 concentration on 13C-fractionation by leaves
            [0 1]};              % Dummy parameter

labelparams = {'d13C_AGveg','d13C_BGveg','d13C_rhizodep', ...
    'f_microbialCO2Uptake','CO2_13CfractionationEffect','dummy'};
    
%     save('paramNames_dummy','labelparams')

M  = numel(distrpar); % number of uncertain parameters [ Sm beta alfa Rs Rf ]

%% Inputparameter for SOILcarb are created

if createParams == 1 || evaluateModelResults == 1

    NU = 500 ; % number of samples to estimate unconditional CDF
    NC = 500 ; % number of samples to estimate conditional CDFs
    n  = 50 ;  % number of conditioning points
    
end
    
if createParams == 1

    % Create input/output samples to estimate the unconditional output CDF:
    Xu = AAT_sampling('lhs',M,'unif',distrpar,NU); % matrix (NU,M)

    % Create input/output samples to estimate the conditional output CDFs:
    [ XX, xc ] = pawn_sampling('lhs',M,'unif',distrpar,n,NC);

    % XX is reformated to serve as input values for SOILcarb
    [r_cell c_cell] = size(XX); % The size of the input cell array
    numberofCells = numel(XX); % The total number of cells in the cell array
    [r_input c_input] = size(XX{1}); % The size of the input array in every cell
    XX_out = NaN(r_input*numberofCells, c_input); % The array that contains all input values for SOILcarb

    % The input values in the cell array are ordered in the output matrix
    % (XX_out), column by column
    counter = 0; % To keep track of the rows where the values have to be stored
    for i = 1:c_cell
        for j = 1:r_cell
            XX_out((counter*r_input)+1:(counter+1)*r_input,:) = XX{j,i};
            counter = counter + 1;
        end
    end

    % Xu and XX are combined to be exported to serve as inputs to SOILcarb
    X_out = [Xu; XX_out];

    % The parameters are stored in a table with a header
    T = array2table(X_out);
    T.Properties.VariableNames(1:6) = labelparams;

    % X_out and Xc are exported as .csv, to be imported in R
    writetable(T, 'X_out.csv'); % To be used as input for SOILcarb
    save('xc', 'xc'); % To be imported into this program together with the SOILcarb results

end

%% The SOILcarb output is loaded to evaluate the results

% Col 1: d13C topsoil
% Col 2: d13C subsoil
% Col 3: d14C topsoil
% Col 4: d14C subsoil
% Col 5: diff d13C
% Col 6: diff d14C

if evaluateModelResults == 1
    
    folder = 'PAWN output - isotopes';
    
    % SOILcarb outputs are loaded
    X_in = readmatrix([folder '/output_SAFE.csv']);
    
    % The first row with ID's is removed
    X_in(:,1) = [];
    
    % The xc input
    xc = load([folder '/xc']);
    xc = xc.xc;
    
    % The parameter names
    labelparams = load([folder '/paramNames_dummy']);
    labelparams = labelparams.labelparams;
    
    % The YY and Yu variables are calculated for every assessed variable
    colNum = 1;
    [YY_1,Yu_1] = calculate_YY_Yu(colNum, NU, NC, n, X_in, M);
    
    colNum = 2;
    [YY_2,Yu_2] = calculate_YY_Yu(colNum, NU, NC, n, X_in, M);
    
    colNum = 3;
    [YY_3,Yu_3] = calculate_YY_Yu(colNum, NU, NC, n, X_in, M);
    
    colNum = 4;
    [YY_4,Yu_4] = calculate_YY_Yu(colNum, NU, NC, n, X_in, M);
    
    colNum = 5;
    [YY_5,Yu_5] = calculate_YY_Yu(colNum, NU, NC, n, X_in, M);
    
    colNum = 6;
    [YY_6,Yu_6] = calculate_YY_Yu(colNum, NU, NC, n, X_in, M);

end

%% The outcomes are plotted in a histogram

hfig = figure;
set(hfig, 'units','centimeters', 'position', [3 3 40 20], 'color',[1,1,1])

% d13C - topsoil
subplot(2,3,1)
col = 1;
histogram(X_in(:,col),50);
xlim([min(X_in(:,col)) max(X_in(:,col))])
title('d13C - topsoil')
ylabel('frequency')

% d13C - subsoil
subplot(2,3,2)
col = 2;
histogram(X_in(:,col),50);
xlim([min(X_in(:,col)) max(X_in(:,col))])
title('d13C - subsoil')
ylabel('frequency')

% d13C - diff topsoil/subsoil
subplot(2,3,3)
col = 5;
histogram(X_in(:,col),50);
xlim([min(X_in(:,col)) max(X_in(:,col))])
title('d13C - diff topsoil/subsoil')
ylabel('frequency')

% d14C - topsoil
subplot(2,3,4)
col = 3;
histogram(X_in(:,col),100);
xlim([min(X_in(:,col)) max(X_in(:,col))])
title('d14C - topsoil')
ylabel('frequency')

% d14C - subsoil
subplot(2,3,5)
col = 4;
histogram(X_in(:,col),100);
xlim([min(X_in(:,col)) max(X_in(:,col))])
title('d14C - subsoil')
ylabel('frequency')

% d14C - diff topsoil/subsoil
subplot(2,3,6)
col = 6;
histogram(X_in(:,col),100);
xlim([min(X_in(:,col)) max(X_in(:,col))])
title('d14C - diff topsoil/subsoil')
ylabel('frequency')

%% Calculating and plotting results

plotEverything = 0;

% YY = pawn_model_evaluation(fun_test,XX,rain,evap) ;

% Estimate unconditional and conditional CDFs:
[ YF_1, Fu_1, Fc_1  ] = pawn_cdfs(Yu_1,YY_1) ;
[ YF_2, Fu_2, Fc_2  ] = pawn_cdfs(Yu_2,YY_2) ;
[ YF_3, Fu_3, Fc_3  ] = pawn_cdfs(Yu_3,YY_3) ;
[ YF_4, Fu_4, Fc_4  ] = pawn_cdfs(Yu_4,YY_4) ;
[ YF_5, Fu_5, Fc_5  ] = pawn_cdfs(Yu_5,YY_5) ;
[ YF_6, Fu_6, Fc_6  ] = pawn_cdfs(Yu_6,YY_6) ;

if plotEverything == 1
    % Plot CDFs:
    figure
    for i=1:M
       subplot(1,M,i)
       pawn_plot_cdf(YF, Fu, Fc(i,:),[],'y (max flow)')
    end

    % Further analyze CDF of one input:
    i = 1 ;
    figure;
    pawn_plot_cdf(YF, Fu, Fc(i,:),xc{i},'y (max flow)',labelparams{i}) % same
    % function as before but exploiting more optional input arguments

end

% Compute KS statistics:
KS_1 = pawn_ks(YF_1,Fu_1,Fc_1) ;
KS_2 = pawn_ks(YF_2,Fu_2,Fc_2) ;
KS_3 = pawn_ks(YF_3,Fu_3,Fc_3) ;
KS_4 = pawn_ks(YF_4,Fu_4,Fc_4) ;
KS_5 = pawn_ks(YF_5,Fu_5,Fc_5) ;
KS_6 = pawn_ks(YF_6,Fu_6,Fc_6) ;

if plotEverything == 1

    % Plot KS statistics:
    figure
    for i=1:M
       subplot(1,M,i)
       pawn_plot_kstest(KS(:,i),NC,NU,0.05,xc{i},labelparams{i})
    end

end

% Compute PAWN index by taking a statistic of KSs (e.g. max):
Pi_1 = max(KS_1);
Pi_2 = max(KS_2);
Pi_3 = max(KS_3);
Pi_4 = max(KS_4);
Pi_5 = max(KS_5);
Pi_6 = max(KS_6);

% Plot
hfig = figure;
set(hfig, 'units','centimeters', 'position', [3 3 40 20], 'color',[1,1,1])

subplot(2,3,1)
plot(0:M+1,ones(size(0:M+1)).*Pi_1(end),'--r')
hold on
boxplot1(Pi_1,labelparams)
set(gca,'fontsize',8,'XTickLabelRotation',45)
title('d13C - topsoil')

subplot(2,3,2)
plot(0:M+1,ones(size(0:M+1)).*Pi_2(end),'--r')
hold on
boxplot1(Pi_2,labelparams)
set(gca,'fontsize',8,'XTickLabelRotation',45)
title('d13C - subsoil')

% !!! This data is in the 5th column
subplot(2,3,3)
plot(0:M+1,ones(size(0:M+1)).*Pi_5(end),'--r')
hold on
boxplot1(Pi_5,labelparams)
set(gca,'fontsize',8,'XTickLabelRotation',45)
title('d13C - diff topsoil/subsoil')

subplot(2,3,4)
plot(0:M+1,ones(size(0:M+1)).*Pi_3(end),'--r')
hold on
boxplot1(Pi_3,labelparams)
set(gca,'fontsize',8,'XTickLabelRotation',45)
title('d14C - topsoil')

subplot(2,3,5)
plot(0:M+1,ones(size(0:M+1)).*Pi_4(end),'--r')
hold on
boxplot1(Pi_4,labelparams)
set(gca,'fontsize',8,'XTickLabelRotation',45)
title('d14C - subsoil')

subplot(2,3,6)
plot(0:M+1,ones(size(0:M+1)).*Pi_6(end),'--r')
hold on
boxplot1(Pi_6,labelparams)
set(gca,'fontsize',8,'XTickLabelRotation',45)
title('d14C - diff topsoil/subsoil')

%% Bootstrapping to get robustness

% Use bootstrapping to assess robustness of PAWN indices:
stat = 'max' ; % statistic to be applied to KSs
Nboot = 100  ; % number of boostrap resamples

[T_m_1, T_lb_1, T_ub_1] = pawn_indices(Yu_1,YY_1,stat,[],Nboot);
[T_m_2, T_lb_2, T_ub_2] = pawn_indices(Yu_2,YY_2,stat,[],Nboot);
[T_m_3, T_lb_3, T_ub_3] = pawn_indices(Yu_3,YY_3,stat,[],Nboot);
[T_m_4, T_lb_4, T_ub_4] = pawn_indices(Yu_4,YY_4,stat,[],Nboot);
[T_m_5, T_lb_5, T_ub_5] = pawn_indices(Yu_5,YY_5,stat,[],Nboot);
[T_m_6, T_lb_6, T_ub_6] = pawn_indices(Yu_6,YY_6,stat,[],Nboot);

fontSize = 14;

hfig = figure;
set(hfig, 'units','centimeters', 'position', [3 3 40 20], 'color',[1,1,1])

subplot(2,3,1)
hold on
plot(0:M+1,ones(size(0:M+1)).*T_m_1(end),'--r')
boxplot1(T_m_1,labelparams,[],T_lb_1,T_ub_1)
set(gca,'fontsize',fontSize,'XTickLabelRotation',45)
title('d13C - topsoil')

subplot(2,3,2)
hold on
plot(0:M+1,ones(size(0:M+1)).*T_m_2(end),'--r')
boxplot1(T_m_2,labelparams,[],T_lb_2,T_ub_2)
set(gca,'fontsize',fontSize,'XTickLabelRotation',45)
title('d13C - subsoil')

subplot(2,3,3)
hold on
plot(0:M+1,ones(size(0:M+1)).*T_m_5(end),'--r')
boxplot1(T_m_5,labelparams,[],T_lb_5,T_ub_5)
set(gca,'fontsize',fontSize,'XTickLabelRotation',45)
title('d13C - diff topsoil/subsoil')

subplot(2,3,4)
hold on
plot(0:M+1,ones(size(0:M+1)).*T_m_3(end),'--r')
boxplot1(T_m_3,labelparams,[],T_lb_3,T_ub_3)
set(gca,'fontsize',fontSize,'XTickLabelRotation',45)
title('d14C - topsoil')

subplot(2,3,5)
hold on
plot(0:M+1,ones(size(0:M+1)).*T_m_4(end),'--r')
boxplot1(T_m_4,labelparams,[],T_lb_4,T_ub_4)
set(gca,'fontsize',fontSize,'XTickLabelRotation',45)
title('d14C - subsoil')

subplot(2,3,6)
hold on
plot(0:M+1,ones(size(0:M+1)).*T_m_6(end),'--r')
boxplot1(T_m_6,labelparams,[],T_lb_6,T_ub_6)
set(gca,'fontsize',fontSize,'XTickLabelRotation',45)
title('d14C - diff topsoil/subsoil')


% [T_m_1, T_lb_1, T_ub_1]
% [T_m_2, T_lb_2, T_ub_2]
% [T_m_3, T_lb_3, T_ub_3]
% [T_m_4, T_lb_4, T_ub_4]
% [T_m_5, T_lb_5, T_ub_5]
% [T_m_6, T_lb_6, T_ub_6]

% ------------------------------------------
% The results are exported for plotting in R
% ------------------------------------------

% Results for d13C - topsoil
out_1 = [T_m_1' T_lb_1' T_ub_1'];
T = array2table(out_1);
T.Properties.VariableNames(1:3) = {'mid', 'lower', 'upper'};
writetable(T, 'Sensitivity_d13C_topsoil.csv');

% Results for d13C - subsoil
out_2 = [T_m_2' T_lb_2' T_ub_2'];
T = array2table(out_2);
T.Properties.VariableNames(1:3) = {'mid', 'lower', 'upper'};
writetable(T, 'Sensitivity_d13C_subsoil.csv');

% Results for diff d13C topsoil - subsoil
out_3 = [T_m_5' T_lb_5' T_ub_5'];
T = array2table(out_3);
T.Properties.VariableNames(1:3) = {'mid', 'lower', 'upper'};
writetable(T, 'Sensitivity_diff_d13C_topsoil_subsoil.csv');

%% Convergence analysis:
stat = 'max' ; % statistic to be applied to KSs
% NCb = [ NC/10 NC/2 NC ] ;
% NUb = [ NU/10 NU/2 NU ] ;

NCb = linspace(NC/10, NC, 10);
NUb = linspace(NU/10, NU, 10);

[T_m_n_1, T_lb_n_1, T_ub_n_1] = pawn_convergence( Yu_1, YY_1, stat, NUb, NCb,[],Nboot );
NN = NUb+n*NCb ;

figure
plot_convergence(T_m_n_1,NN, T_lb_n_1, T_ub_n_1, [], 'no of evals', [], labelparams)

%% Step 4: Apply PAWN to sub-region of the output range

% Compute the PAWN index over a sub-range of the output distribution, for
% instance only output values above a given threshold
thres = 50 ;
[ T_m2, T_lb2, T_ub2 ]= pawn_indices( Yu, YY, stat,[], Nboot,[],'above',thres ) ;

% Plot:
figure; boxplot1(T_m2,labelparams,[],T_lb2,T_ub2)
hold on
plot(0:M+1,ones(size(0:M+1)).*Pi(end),'-r')

plot(0:M+1,ones(size(0:M+1)).*0.15,'-r')

