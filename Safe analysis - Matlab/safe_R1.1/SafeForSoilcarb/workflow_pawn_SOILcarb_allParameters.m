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

addpath('export_fig')

%% Add paths

% Add the folder '/safe_R1.1' and all subfolders to the current path, as we need functions stored there
addpath(genpath('../'))

createParams = 0;
evaluateModelResults = 1;
correctoutputs = 0;

%% Define the parameter ranges

if createParams == 1

    DistrFun  = 'unif'  ; % Parameter distribution

    distrpar  = {[3.21e-05 9.97e-05]; ...   % 1. Db0_stage1
                [0.014 0.332]; ...      % 2. Db_eFold_depth_stage1
                [0.851 0.909]; ...     % 3. betaRoots + fineRootBeta
                [0.696 0.905]; ...       % 4. VmaxD_root_R_stage1
                [0.011 0.998]; ...        % 5. VmaxU_BioAv_stage1
                [243.82 999.530]; ...        % 6. Vmax_ads
                [0.0133 0.291]; ...     % 7. Km_ads
                [44.05 943.17]; ...    % 8. VmaxD_M_stage1
                [0.0582 0.999]; ...    % 9. Km_depol_M_stage1
                [0.03279 0.2694]; ...     % 10. kDes_init
                [0.467 0.999]; ...     % 11. advectionRate_polyC_stage1
                [0 1]};             % 12. Dummy parameter
                
    labelparams = {'Db0_stage1',...
        'Db_eFold_depth_stage1',...
        'betaRoots', ...
        'VmaxD_root_R_stage1',...
        'VmaxU_BioAv_stage1',...
        'Vmax_ads',...
        'Km_ads',...
        'VmaxD_M_stage1',...
        'Km_depol_M_stage1',...
        'kDes_init',...
        'advectionRate_polyC_stage1',...
        'dummy'};
    
    save('paramNames_dummy','labelparams')
end

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
    T.Properties.VariableNames = labelparams;

    % X_out and Xc are exported as .csv, to be imported in R
    writetable(T, 'X_out.csv'); % To be used as input for SOILcarb
    save('xc', 'xc'); % To be imported into this program together with the SOILcarb results

end

%% The SOILcarb output is loaded to evaluate the results

if evaluateModelResults == 1
    
    folder = 'PAWN output - all parameters';
    
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
    
    % =======================================================
    % Rows with unrealistic model outcomes are changed to NaN
    
    if correctoutputs == 1
        
        % Unrealistic model outcomes are converted to NaN
        
        % -----------
        % d13C values
        % -----------
        minVal = -35;
        maxVal = -20;
        
        [r c] = find(X_in(:,1) < minVal | X_in(:,2) < minVal);
        X_in(r,:) = NaN;
        
        [r c] = find(X_in(:,1) > maxVal | X_in(:,2) > maxVal);
        X_in(r,:) = NaN;
        
        % -------------------
        % Topsoil d14C values
        % -------------------
        minVal = -150;
        maxVal = 200;
        
        [r c] = find(X_in(:,3) < minVal | X_in(:,3) > maxVal);
        X_in(r,:) = NaN;
        
        % -------------------
        % Subsoil d14C values
        % -------------------
        minVal = -500;
        maxVal = 200;
        
        [r c] = find(X_in(:,4) < minVal | X_in(:,4) > maxVal);
        X_in(r,:) = NaN;
        
        % -------------------
        % Dd13C values
        % -------------------
        minVal = 0;
        maxVal = 5;
        
        [r c] = find(X_in(:,5) < minVal | X_in(:,5) > maxVal);
        X_in(r,:) = NaN;
        
        % -------------------
        % Dd14C values
        % -------------------
        minVal = -600;
        maxVal = 0;
        
        [r c] = find(X_in(:,6) < minVal | X_in(:,6) > maxVal);
        X_in(r,:) = NaN;
        
        % -------------------
        % Carbon stocks
        % -------------------
        minVal = 0;
        maxVal = 30;
        
        [r c] = find(X_in(:,7) < minVal | X_in(:,7) > maxVal);
        X_in(r,:) = NaN;
    
    end
    % =======================================================
    
    % The column with the results at which we will look
%     colNum = 7;
    
    % YY and Yu are reconstructed
    % Yu are the results for the unconditional run, in an array
    Yu_1 = X_in(1:NU,1);
    Yu_2 = X_in(1:NU,2);
    Yu_3 = X_in(1:NU,3);
    Yu_4 = X_in(1:NU,4);
    Yu_5 = X_in(1:NU,5);
    Yu_6 = X_in(1:NU,6);
    Yu_7 = X_in(1:NU,7);
    
    % YY are the results of the conditional runs, in a cell array
    YY_1 = cell(M,n);
    YY_2 = cell(M,n);
    YY_3 = cell(M,n);
    YY_4 = cell(M,n);
    YY_5 = cell(M,n);
    YY_6 = cell(M,n);
    YY_7 = cell(M,n);
    
    % YY is filled
    counter = 1; % To keep track of the rows where the values have to be stored
    for i = 1:n % Columns
        for j = 1:M % Rows
            startRow = NU + NC*(counter-1) + 1;
            endRow = NU + (counter*NC);
            
            YY_1{j,i} = X_in(startRow:endRow, 1);
            YY_2{j,i} = X_in(startRow:endRow, 2);
            YY_3{j,i} = X_in(startRow:endRow, 3);
            YY_4{j,i} = X_in(startRow:endRow, 4);
            YY_5{j,i} = X_in(startRow:endRow, 5);
            YY_6{j,i} = X_in(startRow:endRow, 6);
            YY_7{j,i} = X_in(startRow:endRow, 7);
            
%             YY{j,i} = X_in((counter*NU)+1:(NU*(counter+1)), colNum);
            counter = counter + 1;
        end
    end

end

%% The outcomes are plotted in a histogram

hfig = figure;
set(hfig, 'units','centimeters', 'position', [3 3 15 12])

colNum = 7;

% hold on
p1 = histogram(X_in(:,colNum),50);
% p1 = bar(X_in(:,colNum),1);
ylabel('Frequency')
xlabel(['colNum ' int2str(colNum)])
% xlim([0 ceil(max(X_in(:,colNum)))])
xlim([min(X_in(:,colNum)) max(X_in(:,colNum))])
% xlim([-50 1000])

%% Calculating and plotting results

plotEverything = 0;

% YY = pawn_model_evaluation(fun_test,XX,rain,evap) ;

% Estimate unconditional and conditional CDFs:
[YF_1, Fu_1, Fc_1] = pawn_cdfs(Yu_1, YY_1);
[YF_2, Fu_2, Fc_2] = pawn_cdfs(Yu_2, YY_2);
[YF_3, Fu_3, Fc_3] = pawn_cdfs(Yu_3, YY_3);
[YF_4, Fu_4, Fc_4] = pawn_cdfs(Yu_4, YY_4);
[YF_5, Fu_5, Fc_5] = pawn_cdfs(Yu_5, YY_5);
[YF_6, Fu_6, Fc_6] = pawn_cdfs(Yu_6, YY_6);
[YF_7, Fu_7, Fc_7] = pawn_cdfs(Yu_7, YY_7);

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
KS_1 = pawn_ks(YF_1, Fu_1, Fc_1);
KS_2 = pawn_ks(YF_2, Fu_2, Fc_2);
KS_3 = pawn_ks(YF_3, Fu_3, Fc_3);
KS_4 = pawn_ks(YF_4, Fu_4, Fc_4);
KS_5 = pawn_ks(YF_5, Fu_5, Fc_5);
KS_6 = pawn_ks(YF_6, Fu_6, Fc_6);
KS_7 = pawn_ks(YF_7, Fu_7, Fc_7);

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
Pi_7 = max(KS_7);

% The values of the dummies are subtracted
Pi_1 = Pi_1 - Pi_1(end); Pi_1(end) = [];
Pi_2 = Pi_2 - Pi_2(end); Pi_2(end) = [];
Pi_3 = Pi_3 - Pi_3(end); Pi_3(end) = [];
Pi_4 = Pi_4 - Pi_4(end); Pi_4(end) = [];
Pi_5 = Pi_5 - Pi_5(end); Pi_5(end) = [];
Pi_6 = Pi_6 - Pi_6(end); Pi_6(end) = [];
Pi_7 = Pi_7 - Pi_7(end); Pi_7(end) = [];

% Al Pi values are combined for plotting
allPi = [Pi_1; Pi_2; Pi_3; Pi_4; Pi_5; Pi_6; Pi_7]';

% Values <= are given a value of 0
allPi(allPi < 0) = 0;

% --------
% Plotting
% --------

% Colors
c1 = [27,158,119; ...
217,95,2; ...
117,112,179; ...
231,41,138; ...
102,166,30; ...
230,171,2; ...
166,118,29]./255; ...

hfig = figure;
set(hfig, 'units','centimeters', 'position', [3 3 35 20], 'color', [1 1 1])

b = bar(allPi);

for ii = 1:size(allPi,2)
    b(ii).FaceColor = c1(ii,:);
end

% b(7).LineWidth = 2;

labelparams = {'D_b(0)*',...
        'z_b',...
        '\beta_r*', ...
        'V_{max,POC-r}*',...
        'V_{maxU,mic-r}*',...
        'V_{max,ads}*',...
        'K_{m\_ads}*',...
        'V_{max,DOC-b}*',...
        'K_{m\_DOC-b}*',...
        'k_{deprotect}*',...
        '\nu'};

xLabels = labelparams;

set(gca,'xticklabel',xLabels)
set(gca,'fontsize',16,'XTickLabelRotation',45)
ylabel('PAWN sensitivity index (max(KS))')

set(b, {'DisplayName'}, {'Topsoil \delta^{13}C', 'Subsoil \delta^{13}C', 'Topsoil \Delta^{14}C','Subsoil \Delta^{14}C', '\Delta(\delta^{13}C)', '\Delta(\Delta^{14}C)', 'SOC stock'}')

legend('NumColumns',1, 'Orientation','horizontal', 'Location', 'northwest')

% The plot is saved
% export_fig 'Figure_supplement_sensitivityAllParameters.png' -r300 -nocrop -transparent


%% Bootstrapping to get robustness

% Use bootstrapping to assess robustness of PAWN indices:
stat = 'max' ; % statistic to be applied to KSs
Nboot = 100  ; % number of boostrap resamples
[ T_m, T_lb, T_ub ] = pawn_indices(Yu,YY,stat,[],Nboot);

% Plot:
figure;
hold on 
plot(0:M+1,ones(size(0:M+1)).*T_m(end),'--r')
boxplot1(T_m,labelparams,[],T_lb,T_ub)
set(gca,'fontsize',8)

%% Convergence analysis:
stat = 'max' ; % statistic to be applied to KSs
% NCb = [ NC/10 NC/2 NC ] ;
% NUb = [ NU/10 NU/2 NU ] ;

NCb = linspace(NC/10, NC, 10);
NUb = linspace(NU/10, NU, 10);

[ T_m_n, T_lb_n, T_ub_n ] = pawn_convergence( Yu, YY, stat, NUb, NCb,[],Nboot );
NN = NUb+n*NCb ;

figure
plot_convergence(T_m_n,NN,T_lb_n,T_ub_n,[],'no of evals',[],labelparams)

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

