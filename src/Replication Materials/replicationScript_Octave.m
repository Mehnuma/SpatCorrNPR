% This script contains all the code snippets included in the manuscript
% titled 'SpatCorrNPR: A MATLAB toolbox for Local Polynomial Regression
% with Correlated Errors', describing the use of various modules of the
% SpatCorrNPR toolbox and replication of the practical examples.
%
% Authors: Mehnuma Tabassum, Kris De Brabanter
% Email: mehnuma@iastate.edu
% ********************************************************************** %

%% Codes for the Section "SpatCorrNPR: The MATLAB Toolbox"
% There are 6 modules in this toolbox (5 in Octave), which are represented in the
% following sections.

%% Module-1: Data Import
% The following lines of code can be used to import the data. The first
% dataset (testData.csv) used is without the original function, m values.
data = dlmread('testData.csv',",",1,0);
X_1d = data(:,1);        % for 1-D input
X_2d = data(:,1:2);      % for 2-D input
y = data(:,3);           % The response

% If the user wishes to use a dataset with m values, they can import the
% testDataWithM.csv file. In that case, the data import would be like the
% following:
dataWithM = dlmread('testDataWithM.csv', ",",1,0);
X_1d = dataWithM(:,1);   % for 1-D input
X_2d = dataWithM(:,1:2); % for 2-D input
y = dataWithM(:,3);      % The response
m = dataWithM(:,4);      % The true function

% Note: The first three columns of these two datasets are the same; the
% only difference is that the second dataset contain a column for the true
% funtion values, m.

%% Module-2: Bandwidth Optimization (Table-4)
% Case: 1-D X (Table 4, Column 3)
% The bandwidth optimization for 1-D data is fairly straightforward as it
% is a scalar value. The following lines of codes demonstrate bandwidth
% optimization for the all the objectives and kernels combination. The
% limits for the 1-D case are [0.01, 0.5]. The syntax for calling the
% optimizer is as follows.

% Note: Some results may be different because of using a special
% algorithm to set the seed value for the optimizer search.

% Unadjusted GCV (Epanechnikov Kernel)
optimal_h_gcv_1d_ep = bwOptimization(X_1d, y, [], 'Epanechnikov', ...
    'Unadjusted GCV', 0.01, 0.5, []); % for 1-D X
% Unadjusted GCV (Gaussian Kernel)
optimal_h_gcv_1d_gau = bwOptimization(X_1d, y, [], 'Gaussian', ...
    'Unadjusted GCV', 0.01, 0.5, []); % for 1-D X

% If the optimization criterion is changed to Corrected GCV (GCVc), the R
% matrix needs to be imported, which can be done using the following line:
RMatrix = dlmread('trueR.csv');
% Corrected GCV, GCVc (Epanechnikov Kernel)
optimal_h_gcvc_1d_ep = bwOptimization(X_1d, y, [], 'Epanechnikov', ...
    'Corrected GCV (GCVc)', 0.01, 0.5, RMatrix); % for 1-D X
% Corrected GCV, GCVc (Gaussian Kernel)
optimal_h_gcvc_1d_gau = bwOptimization(X_1d, y, [], 'Gaussian', ...
    'Corrected GCV (GCVc)', 0.01, 0.5, RMatrix); % for 1-D X

% Corrected and Estimated GCV, GCVce (Epanechnikov Kernel)
optimal_h_gcvce_1d_ep = bwOptimization(X_1d, y, [], 'Epanechnikov', ...
    'Corrected and Estimated GCV (GCVce)', 0.01, 0.5, []); % for 1-D X
% Corrected and Estimated GCV, GCVce (Gaussian Kernel)
optimal_h_gcvce_1d_gau = bwOptimization(X_1d, y, [], 'Gaussian', ...
    'Corrected and Estimated GCV (GCVce)', 0.01, 0.5, []); % for 1-D X

% Residual Sum of Squares, RSS (Epanechnikov Kernel)
optimal_h_rss_1d_ep = bwOptimization(X_1d, y, [], 'Epanechnikov', ...
    'Residual Sum of Squares (RSS)', 0.04, 0.5, []); % for 1-D X
% Residual Sum of Squares, RSS (Gaussian Kernel)
optimal_h_rss_1d_gau = bwOptimization(X_1d, y, [], 'Gaussian', ...
    'Residual Sum of Squares (RSS)', 0.04, 0.5, []); % for 1-D X

% The m values are required if the Asymptotic Squared Error (ASE) objective
% is chosen. We can use the following lines in that case:
% Asymptotic Squared Error, ASE (Epanechnikov Kernel)
optimal_h_ase_1d_ep = bwOptimization(X_1d, y, m, 'Epanechnikov', ...
    'Asymptotic Squared Error (ASE)', 0.01, 0.5, []); % for 1-D X
% Asymptotic Squared Error, ASE (Gaussian Kernel)
optimal_h_ase_1d_gau = bwOptimization(X_1d, y, m, 'Gaussian', ...
    'Asymptotic Squared Error (ASE)', 0.01, 0.5, []); % for 1-D X

% Case: 2-D X  (Table 4, Column 4)
% The optimization for the 2-D input is similar, except the bounds are
% required to be 2 × 2 symmetric matrices. Below are the codes for the 2-D
% bandwidth optimization for all objectives and kernel combinations.

% Unadjusted GCV (Epanechnikov Kernel)
optimal_H_gcv_2d_ep = bwOptimization(X_2d, y,[], 'Epanechnikov', 'Unadjusted GCV',...
     [0.08, 0.001; 0.001 0.08], [0.5, 0.1; 0.1 0.5], []); % for 2-D X
% Unadjusted GCV (Gaussian Kernel)
optimal_H_gcv_2d_gau = bwOptimization(X_2d, y,[], 'Gaussian', 'Unadjusted GCV',...
     [0.08, 0.001; 0.001 0.08], [0.5, 0.1; 0.1 0.5], []); % for 2-D X

% The syntax for the Corrected GCV (GCVc) criterion for the 2-D is similar
% to the 1-D case, which is shown below:
% Corrected GCV, GCVc (Epanechnikov Kernel)
optimal_H_gcvc_2d_ep = bwOptimization(X_2d, y,[], 'Epanechnikov', ...
    'Corrected GCV (GCVc)', [0.08, 0.001; 0.001 0.08], ...
    [0.5, 0.1; 0.1 0.5], RMatrix); % for 2-D X
    % The R Matrix can be imported as in Line 57.
% Corrected GCV, GCVc (Gaussian Kernel)
optimal_H_gcvc_2d_gau = bwOptimization(X_2d, y,[], 'Gaussian', ...
    'Corrected GCV (GCVc)', [0.08, 0.001; 0.001 0.08], ...
    [0.5, 0.1; 0.1 0.5], RMatrix); % for 2-D X

% Corrected and Estimated GCV, GCVce (Epanechnikov Kernel)
optimal_H_gcvce_2d_ep = bwOptimization(X_2d, y,[], 'Epanechnikov', ...
    'Corrected and Estimated GCV (GCVce)', [0.08, 0.001; 0.001 0.08], ...
    [0.5, 0.1; 0.1 0.5], []); % for 2-D X
% Corrected and Estimated GCV, GCVce (Gaussian Kernel)
optimal_H_gcvce_2d_gau = bwOptimization(X_2d, y,[], 'Gaussian', ...
    'Corrected and Estimated GCV (GCVce)', [0.08, 0.001; 0.001 0.08], ...
    [0.5, 0.1; 0.1 0.5], []); % for 2-D X

% Residual Sum of Squares, RSS (Epanechnikov Kernel)
optimal_H_rss_2d_ep = bwOptimization(X_2d, y, [], 'Epanechnikov', ...
    'Residual Sum of Squares (RSS)', [0.2, 0.001; 0.001 0.2], ...
    [0.5, 0.1; 0.1 0.5], []); % for 2-D X
% Residual Sum of Squares, RSS (Gaussian Kernel)
optimal_H_rss_2d_gau = bwOptimization(X_2d, y, [], 'Gaussian', ...
    'Residual Sum of Squares (RSS)', [0.2, 0.001; 0.001 0.2], ...
    [0.5, 0.1; 0.1 0.5], []); % for 2-D X

% For the Asymptotic Squared Error (ASE) objective, the user can execute
% the following lines:
% Asymptotic Squared Error, ASE (Epanechnikov Kernel)
optimal_H_ase_2d_ep = bwOptimization(X_2d, y, m, 'Epanechnikov', ...
    'Asymptotic Squared Error (ASE)', [0.08, 0.001; 0.001 0.08], ...
    [0.5, 0.1; 0.1 0.5], []); % for 2-D X
% Asymptotic Squared Error, ASE (Gaussian Kernel)
optimal_H_ase_2d_gau = bwOptimization(X_2d, y, m, 'Gaussian', ...
    'Asymptotic Squared Error (ASE)', [0.08, 0.001; 0.001 0.08], ...
    [0.5, 0.1; 0.1 0.5], []); % for 2-D X

disp("%%%%%%%%%%%%%%%%%%%%%%%% Table 4 Results %%%%%%%%%%%%%%%%%%%%%%%%")
disp('Optimal Bandwidth for 1-D X (Column 3 of Table 4):')

disp('Unadjusted GCV, Epanechnikov kernel: '), disp(optimal_h_gcv_1d_ep);
disp('Unadjusted GCV, Gaussian kernel: '),     disp(optimal_h_gcv_1d_gau);

disp('Corrected GCV, Epanechnikov kernel: '), disp(optimal_h_gcvc_1d_ep);
disp('Corrected GCV, Gaussian kernel: '),    disp(optimal_h_gcvc_1d_gau);

disp('Corrected and Estimated GCV, Epanechnikov kernel: '), disp(optimal_h_gcvce_1d_ep);
disp('Corrected and Estimated GCV, Gaussian kernel: '),     disp(optimal_h_gcvce_1d_gau);

disp('Residual Sum of Squares, Epanechnikov kernel: '), disp(optimal_h_rss_1d_ep);
disp('Residual Sum of Squares, Gaussian kernel: '),    disp(optimal_h_rss_1d_gau);

disp('Asymptotic Squared Error, Epanechnikov kernel: '), disp(optimal_h_ase_1d_ep);
disp('Asymptotic Squared Error, Gaussian kernel: '),    disp(optimal_h_ase_1d_gau);

disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
disp('Optimal Bandwidth for 2-D X (Column 4 of Table 4):')

disp('Unadjusted GCV, Epanechnikov kernel: '), disp(optimal_H_gcv_2d_ep);
disp('Unadjusted GCV, Gaussian kernel: '),     disp(optimal_H_gcv_2d_gau);

disp('Corrected GCV, Epanechnikov kernel: '), disp(optimal_H_gcvc_2d_ep);
disp('Corrected GCV, Gaussian kernel: '),    disp(optimal_H_gcvc_2d_gau);

disp('Corrected and Estimated GCV, Epanechnikov kernel: '), disp(optimal_H_gcvce_2d_ep);
disp('Corrected and Estimated GCV, Gaussian kernel: '),     disp(optimal_H_gcvce_2d_gau);

disp('Residual Sum of Squares, Epanechnikov kernel: '), disp(optimal_H_rss_2d_ep);
disp('Residual Sum of Squares, Gaussian kernel: '),    disp(optimal_H_rss_2d_gau);

disp('Asymptotic Squared Error, Epanechnikov kernel: '), disp(optimal_H_ase_2d_ep);
disp('Asymptotic Squared Error, Gaussian kernel: '),    disp(optimal_H_ase_2d_gau);

disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

%% Module-3: Error Evaluation (Table-5)
% The following four snippets calculate the PRESS values for 1-D and 2-D X
% with Epanechnikov and Gaussian kernels, respectively. The results are
% shown for arbitrary bandwidth values.
% PRESS (1-D, both kernels)
press_1d_ep = evaluate_error(X_1d, y, [], 'PRESS', 0.0824, 'Epanechnikov'); % for 1-D X
press_1d_gau = evaluate_error(X_1d, y, [], 'PRESS', 0.0824, 'Gaussian');     % for 1-D X

% PRESS (2-D, both kernels)
press_2d_ep = evaluate_error(X_2d, y, [], 'PRESS', [0.25 0.1; 0.1 0.25], ...
           'Epanechnikov'); % for 2-D X
press_2d_gau = evaluate_error(X_2d, y, [], 'PRESS', [0.25 0.1; 0.1 0.25], ...
           'Gaussian'); % for 2-D X

% Similarly, we can find the ASE values for the datasets. The additional
% input here is the m values.
% ASE (1-D, both kernels)
ase_1d_ep = evaluate_error(X_1d, y, m, 'ASE', 0.0824, 'Epanechnikov'); % for 1-D X
ase_1d_gau = evaluate_error(X_1d, y, m, 'ASE', 0.0824, 'Gaussian');     % for 1-D X

% ASE (2-D, both kernels)
ase_2d_ep = evaluate_error(X_2d, y, m, 'ASE', [0.25 0.1; 0.1 0.25], ...
        'Epanechnikov'); % for 2-D X
ase_2d_gau = evaluate_error(X_2d, y, m, 'ASE', [0.25 0.1; 0.1 0.25], ...
        'Gaussian'); % for 2-D X

disp("%%%%%%%%%%%%%%%%%%%%%%%% Table 5 Results %%%%%%%%%%%%%%%%%%%%%%%%")
disp('Error Criterion Values for 1-D X (Column 3 of Table 5):')

disp('PRESS Value, Epanechnikov kernel: '), disp(press_1d_ep)
disp('PRESS Value, Gaussian kernel: '), disp(press_1d_gau)

disp('ASE Value, Epanechnikov kernel: '), disp(ase_1d_ep)
disp('ASE Value, Gaussian kernel: '), disp(ase_1d_gau)

disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
disp('Error Criterion Values for 2-D X (Column 4 of Table 5):')

disp('PRESS Value, Epanechnikov kernel: '), disp(press_2d_ep)
disp('PRESS Value, Gaussian kernel: '), disp(press_2d_gau)

disp('ASE Value, Epanechnikov kernel: '), disp(ase_2d_ep)
disp('ASE Value, Gaussian kernel: '), disp(ase_2d_gau)
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

%% Module-4: Plotting
% Plot Raw Data (Figure 1)
% To get a scatterplot of 1-D X, the following code is used:
plotRawData(X_1d,y)   % Figure 1(a)
xlabel('x1')
ylabel('y')
print -dpng Fig1a.png
close

% For 2-D X, the following line produces a 3-D scatterplot.
plotRawData(X_2d,y)   % Figure 1(b)
xlabel('x1')
ylabel('x2')
zlabel('y')
print -dpng Fig1b.png
close

% Plot Objective Function
% The following snippet provide the syntax for plotting the objective
% functions with Epanechnikov kernel for 1-D X and 2-D X,
% respectively.
% Figure 2
plotObjectiveFunction(X_1d, y,[], [0.01, 0.5], 'Epanechnikov', ...
    'Unadjusted GCV', [])   % Figure 2(a)
print -dpng Fig2a.png
close

plotObjectiveFunction(X_1d, y,[], [0.01, 0.5], 'Epanechnikov', ...
    'Corrected and Estimated GCV (GCVce)', [])   % Figure 2(b)
print -dpng Fig2b.png
close

plotObjectiveFunction(X_1d, y,[], [0.01, 5], 'Epanechnikov', ...
    'Residual Sum of Squares (RSS)', [])   % Figure 2(c)
print -dpng Fig2c.png
close

plotObjectiveFunction(X_1d, y, m, [0.01, 5], 'Epanechnikov', ...
    'Asymptotic Squared Error (ASE)', [])    % Figure 2(d)
print -dpng Fig2d.png
close

% Figure 3
plotObjectiveFunction(X_2d, y,[], [0.05 0.05; 2 2], 'Epanechnikov', ...
    'Unadjusted GCV', [])   % Figure 3(a)
colormap jet
print -dpng Fig3a.png
close

plotObjectiveFunction(X_2d, y,[], [0.05 0.05; 2 2], 'Epanechnikov', ...
    'Corrected and Estimated GCV (GCVce)', [])   % Figure 3(b)
colormap jet
print -dpng Fig3b.png
close

plotObjectiveFunction(X_2d, y,[], [0.05 0.05; 2 2], 'Epanechnikov', ...
    'Residual Sum of Squares (RSS)', [])   % Figure 3(c)
print -dpng Fig3c.png
close

plotObjectiveFunction(X_2d, y, m, [0.05 0.05; 2 2], 'Epanechnikov', ...
    'Asymptotic Squared Error (ASE)', [])    % Figure 3(d)
print -dpng Fig3d.png
close

% The syntax for using the Corrected GCV (GCVc) with Epanechnikov kernel
% is as follows: (ONLY AVILABLE FOR THE TOOLBOX, NOT SHOWN IN THE
% MANUSCRIPT)

% RMatrix = dlmread('trueR.csv');
% plotObjectiveFunction(X_1d, y,[], [0.01, 0.5], 'Epanechnikov', ...
%    'Corrected GCV (GCVc)', RMatrix)    % for 1-D X
% plotObjectiveFunction(X_2d, y,[], [0.05 0.05; 2.5 2.5], 'Epanechnikov', ...
%    'Corrected GCV (GCVc)', RMatrix)    % for 2-D X

% Note: These plots are not shown in the manuscript. The plots shown in the
% manuscript are common to the app and the command line toolbox. GCVc
% plotting is only available in the command line toolbox.

% Plot Smoothed Function (Figure 4)
% Below is the syntax for plotting the smoothed function for 1-D X and 2-D
% X, respectively.
plotSmoothedFunction(X_1d, y, 0.0824, 'Epanechnikov', []) % Figure 4(a)
xlabel('x1')
ylabel('y')
print -dpng Fig4a.png
close

plotSmoothedFunction(X_2d, y, [0.233 0.058; 0.058 0.311], 'Epanechnikov', [])    % Figure 4(b)
view([-24, 23])
xlabel('x1')
ylabel('x2')
zlabel('y')
print -dpng Fig4b.png
close

%% Module-5: Covariance Estimation
% Here are some example codes to produce semivariograms and calculate the
% covariance matrix of the local linear regression residuals.
% Warning: These plots might need some time to load compared to MATLAB.
% Figure 5
pkg load statistics
CMat_exp = covarianceEstimation(X_2d, y, [0.233 0.058; 0.058 0.311], ...
    'Epanechnikov', [], [], [], [], [], [], 'Exponential', [], []); % Model: Exponential, Figure 5(a)
print -dpng Fig5a.png
close

CMat_gau = covarianceEstimation(X_2d, y, [0.233 0.058; 0.058 0.311], ...
    'Epanechnikov', [], [], [], [], [], [], 'Gaussian', [], []);   % Model: Gaussian, Figure 5(b)
print -dpng Fig5b.png
close

%% Module-6:  Geoplotting
% Please note that this module uses the MATLAB Mapping Toolbox, which is
% not available to use in Octave. Therefore, the plots in this module can
% only be produced in MATLAB.

% ********************************************************************** %

