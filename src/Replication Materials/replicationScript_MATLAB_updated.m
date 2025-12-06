% This script contains all the code snippets included in the manuscript
% titled 'SpatCorrNPR: A MATLAB toolbox for Local Polynomial Regression 
% with Correlated Errors', describing the use of various modules of the
% SpatCorrNPR toolbox and replication of the practical examples.
%
% Authors: Mehnuma Tabassum, Kris De Brabanter
% Email: mehnuma@iastate.edu
% ********************************************************************** %
%% Codes for the Section "SpatCorrNPR: descriptions of GUI and command line toolbox"
% There are 6 modules in this toolbox, which are represented in the 6
% following sections.

%% Module-1: Data Import
% The following lines of code can be used to import the data. The first 
% dataset (testData.csv) used is without the original function, m values.

data = readmatrix('testData.csv'); 
X_1d = data(:,1);        % for 1-D input 
X_2d = data(:,1:2);      % for 2-D input 
y = data(:,3);           % The response

% If the user wishes to use a dataset with m values, they can import the 
% testDataWithM.csv file. In that case, the data import would be like the 
% following:

dataWithM = readmatrix('testDataWithM.csv'); 
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
% matrix needs to be imported, which can be done using the following lines:
RMatrix = readmatrix('trueR.csv'); 
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

%% Module-4: Plotting
% Plot Raw Data (Figure 1)
% To get a scatterplot of 1-D X, the following code is used:
plotRawData(X_1d,y)   % Figure 1(a)
colormap(brewermap([],"Spectral"))
xlabel('x1')
ylabel('y')
saveas(gcf, 'fig1a.png')
close

% For 2-D X, the following line produces a 3-D scatterplot.
plotRawData(X_2d,y)   % Figure 1(b)
colormap(brewermap([],"Spectral"))
xlabel('x1')
ylabel('x2')
zlabel('y')
saveas(gcf, 'fig1b.png')
close

% Plot Objective Function
% The following snippet provide the syntax for plotting the objective
% functions with Epanechnikov kernel for 1-D X and 2-D X,
% respectively.
% Figure 2
plotObjectiveFunction(X_1d, y,[], [0.01, 0.5], 'Epanechnikov', ...
    'Unadjusted GCV', [])   % Figure 2(a)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig2a.png')
close

plotObjectiveFunction(X_1d, y,[], [0.01, 0.5], 'Epanechnikov', ...
    'Corrected and Estimated GCV (GCVce)', [])   % Figure 2(b)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig2b.png')
close

plotObjectiveFunction(X_1d, y,[], [0.01, 5], 'Epanechnikov', ...
    'Residual Sum of Squares (RSS)', [])   % Figure 2(c)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig2c.png')
close

plotObjectiveFunction(X_1d, y, m, [0.01, 5], 'Epanechnikov', ...
    'Asymptotic Squared Error (ASE)', [])    % Figure 2(d)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig2d.png')
close

% Figure 3
plotObjectiveFunction(X_2d, y,[], [0.05 0.05; 2 2], 'Epanechnikov', ...
    'Unadjusted GCV', [])   % Figure 3(a)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig3a.png')
close

plotObjectiveFunction(X_2d, y,[], [0.05 0.05; 2 2], 'Epanechnikov', ...
    'Corrected and Estimated GCV (GCVce)', [])   % Figure 3(b)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig3b.png')
close

plotObjectiveFunction(X_2d, y,[], [0.05 0.05; 2 2], 'Epanechnikov', ...
    'Residual Sum of Squares (RSS)', [])   % Figure 3(c)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig3c.png')
close

plotObjectiveFunction(X_2d, y, m, [0.05 0.05; 2 2], 'Epanechnikov', ...
    'Asymptotic Squared Error (ASE)', [])    % Figure 3(d)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig3d.png')
close

% % The syntax for using the Corrected GCV (GCVc) with Epanechnikov kernel 
% % is as follows: 
% RMatrix = readmatrix('trueR.csv');
% plotObjectiveFunction(X_1d, y,[], [0.01, 0.5], 'Epanechnikov', ...
%     'Corrected GCV (GCVc)', RMatrix)    % for 1-D X
% plotObjectiveFunction(X_2d, y,[], [0.05 0.05; 2.5 2.5], 'Epanechnikov', ...
%     'Corrected GCV (GCVc)', RMatrix)    % for 2-D X

% Note: These plots are not shown in the manuscript. The plots shown in the
% manuscript are common to the app and the command line toolbox. GCVc
% plotting is only available in the command line toolbox. 

% Plot Smoothed Function (Figure 4)
% Below is the syntax for plotting the smoothed function for 1-D X and 2-D
% X, respectively.

plotSmoothedFunction(X_1d, y, 0.0824, 'Epanechnikov', []) % Figure 4(a)
colormap(brewermap([],"Spectral"))
xlabel('x1')
ylabel('y')
saveas(gcf, 'fig4a.png')
close

plotSmoothedFunction(X_2d, y, [0.233 0.058; 0.058 0.311], 'Epanechnikov', [])
                                                          % Figure 4(b)
colormap(brewermap([],"RdYlGn"))
view([-24, 23])
xlabel('x1')
ylabel('x2')
zlabel('y')                                                        
saveas(gcf, 'fig4b.png')
close

%% Module-5: Covariance Estimation
% Here are some example codes to produce semivariograms and calculate the
% covariance matrix of the local linear regression residuals.

% Warning: These plots might need some time to load compared to MATLAB.

% Figure 5
CMat_exp = covarianceEstimation(X_2d, y, [0.233 0.058; 0.058 0.311], ...
    'Epanechnikov', [], [], [], [], [], [], 'Exponential', [], []);
     % Model: Exponential, Figure 5(a)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig5a.png')
close

CMat_gau = covarianceEstimation(X_2d, y, [0.233 0.058; 0.058 0.311], ...
    'Epanechnikov', [], [], [], [], [], [], 'Gaussian', [], []);
     % Model: Gaussian, Figure 5(b)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig5b.png')
close

%% Module-6: Geoplotting (ONLY available in MATLAB)
% Example-1: The first geoPlot is for the US precipitation data. The 
% following codes are used:
geoData = readmatrix('usPrecipitation.csv');
X = geoData(:,1:2); y = geoData(:,3);
geoPlot(X, y, [3.615 0.0026; 0.0026 2.5873], 'Epanechnikov', [], ...
    'Contiguous US', [], [], 'Positive')  % Figure 6(a)
colormap(brewermap([],"RdBu"))
saveas(gcf, 'fig6a.png')
close

% Example-2: The following example uses a toy dataset from Australia, and 
% defines the perimeter by importing the perimeter data, instead of 
% choosing the geographicRegion. Here is the code:
geoPerimeter = readmatrix('australiaPerimeter.csv');
Lat = geoPerimeter(:,1);
Lon = geoPerimeter(:,2);

geoData = readmatrix('ausData.csv');
X = geoData(:,1:2); y = geoData(:,3);
geoPlot(X, y, [5.2 0.04; 0.04 7.1], 'Epanechnikov', [], 'N/A', Lat, Lon, ...
    'other')   % Figure 6(b)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig6b.png')
close

%% Codes for the Section "Some Practical Examples"
iter = 10;
addpath(append(pwd,'\ooDACE-1.4\ooDACE'))
addpath(append(pwd,'UQLab_Rel2.1.0\UQLab_Rel2.1.0'))
addpath(append(pwd, '\UQLab_Rel2.1.0\UQLab_Rel2.1.0\core'))
uqlab_install

%% Dataset-1 (Texas Wildfire)
% SpatCorrNPR
time_wildfire_spatcorr = nan(1,iter);
n1_wildfire = 178; n2_wildfire = 76; ind_wildfire1 = 1; ind_wildfire2 = 1;
rmse_wildfire_spatcorr = nan(1,iter);
for i=1:iter
    tic
    dataTrain = readmatrix('Train_wildfire.csv');
    dataTest = readmatrix('Test_wildfire.csv');
    X = dataTrain(ind_wildfire1:n1_wildfire*i,1:2); y = dataTrain(ind_wildfire1:n1_wildfire*i,3);
    y = (y-min(y))/(max(y)-min(y));
    xtest = dataTest(ind_wildfire2:n2_wildfire*i,1:2); ytest = dataTest(ind_wildfire2:n2_wildfire*i,3);
    ytest = (ytest-min(ytest))/(max(ytest)-min(ytest));
    H = [1.25 0.01; 0.01 1.55];
    yhat_spatcorr = locpol_spatcorrnpr(X,y,H,'Epanechnikov',xtest);
    rmse_wildfire_spatcorr(i) = sqrt(mean((ytest-yhat_spatcorr).^2));
    time_wildfire_spatcorr(i) = toc;
    ind_wildfire1 = n1_wildfire*i+1;
    ind_wildfire2 = n2_wildfire*i+1;
end

% ooDACE
time_wildfire_oodace = nan(1,iter);
n1_wildfire = 178; n2_wildfire = 76; ind_wildfire1 = 1; ind_wildfire2 = 1;
rmse_wildfire_oodace = nan(1,iter);
for i=1:iter
    tic
    dataTrain = readmatrix('Train_wildfire.csv');
    dataTest = readmatrix('Test_wildfire.csv');
    X = dataTrain(ind_wildfire1:n1_wildfire*i,1:2); y = dataTrain(ind_wildfire1:n1_wildfire*i,3);
    y = (y-min(y))/(max(y)-min(y));
    xtest = dataTest(ind_wildfire2:n2_wildfire*i,1:2); ytest = dataTest(ind_wildfire2:n2_wildfire*i,3);
    ytest = (ytest-min(ytest))/(max(ytest)-min(ytest));
    k = oodacefit(X,y);
    [yhat_oodace, ~] = k.predict(xtest);   
    rmse_wildfire_oodace(i) = sqrt(mean((ytest-yhat_oodace).^2));
    time_wildfire_oodace(i) = toc;
    ind_wildfire1 = n1_wildfire*i+1;
    ind_wildfire2 = n2_wildfire*i+1;
end

% UQLab
time_wildfire_uqlab = nan(1,iter);
n1_wildfire = 178; n2_wildfire = 76; ind_wildfire1 = 1; ind_wildfire2 = 1;
rmse_wildfire_uqlab = nan(1,iter);
for i=1:iter
    tic
    dataTrain = readmatrix('Train_wildfire.csv');
    dataTest = readmatrix('Test_wildfire.csv');
    X = dataTrain(ind_wildfire1:n1_wildfire*i,1:2); y = dataTrain(ind_wildfire1:n1_wildfire*i,3);
    y = (y-min(y))/(max(y)-min(y));
    xtest = dataTest(ind_wildfire2:n2_wildfire*i,1:2); ytest = dataTest(ind_wildfire2:n2_wildfire*i,3);
    ytest = (ytest-min(ytest))/(max(ytest)-min(ytest));
    uqlab;
    MetaOpts.Type = 'Metamodel';
    MetaOpts.MetaType = 'Kriging';
    MetaOpts.ExpDesign.Sampling = 'User';
    MetaOpts.ExpDesign.X = X;
    MetaOpts.ExpDesign.Y = y;
    myKriging = uq_createModel(MetaOpts);
    yhat_uqlab = uq_evalModel(myKriging, xtest);
    rmse_wildfire_uqlab(i) = sqrt(mean((ytest-yhat_uqlab).^2));
    time_wildfire_uqlab(i) = toc;
    ind_wildfire1 = n1_wildfire*i+1;
    ind_wildfire2 = n2_wildfire*i+1;
end

%% Dataset-2 (US Precipitation)
% SpatCorrNPR
time_precip_spatcorr = nan(1,iter);
n1_prec = 738; n2_prec = 315; ind_prec1 = 1; ind_prec2 = 1;
rmse_precip_spatcorr = nan(1,iter);
for i=1:iter
    tic
    dataTrain = readmatrix('Train_precipitation.csv');
    dataTest = readmatrix('Test_precipitation.csv');
    X = dataTrain(ind_prec1:n1_prec*i,1:2); y = dataTrain(ind_prec1:n1_prec*i,3);
    y = (y-min(y))/(max(y)-min(y));
    xtest = dataTest(ind_prec2:n2_prec*i,1:2); ytest = dataTest(ind_prec2:n2_prec*i,3);
    ytest = (ytest-min(ytest))/(max(ytest)-min(ytest));
    H = [3.15 0.001; 0.001 3.63]; 
    yhat_spatcorr = locpol_spatcorrnpr(X,y,H,'Epanechnikov',xtest);
    rmse_precip_spatcorr(i) = sqrt(mean((ytest-yhat_spatcorr).^2));
    time_precip_spatcorr(i) = toc;
    ind_prec1 = n1_prec*i+1;
    ind_prec2 = n2_prec*i+1;
end

% ooDACE
time_precip_oodace = nan(1,iter);
rmse_precip_oodace = nan(1,iter);
n1_prec = 738; n2_prec = 315; ind_prec1 = 1; ind_prec2 = 1;
for i=1:iter
    tic
    dataTrain = readmatrix('Train_precipitation.csv');
    dataTest = readmatrix('Test_precipitation.csv');
    X = dataTrain(ind_prec1:n1_prec*i,1:2); y = dataTrain(ind_prec1:n1_prec*i,3);
    y = (y-min(y))/(max(y)-min(y));
    xtest = dataTest(ind_prec2:n2_prec*i,1:2); ytest = dataTest(ind_prec2:n2_prec*i,3);
    ytest = (ytest-min(ytest))/(max(ytest)-min(ytest));
    k = oodacefit(X,y);
    [yhat_oodace, ~] = k.predict(xtest);
    rmse_precip_oodace(i) = sqrt(mean((ytest-yhat_oodace).^2));
    time_precip_oodace(i) = toc;
    ind_prec1 = n1_prec*i+1;
    ind_prec2 = n2_prec*i+1;
end

% UQLab
time_precip_uqlab = nan(1,iter);
n1_prec = 738; n2_prec = 315; ind_prec1 = 1; ind_prec2 = 1;
rmse_precip_uqlab = nan(1,iter);
rng('default')
for i=1:iter
    tic
    dataTrain = readmatrix('Train_precipitation.csv');
    dataTest = readmatrix('Test_precipitation.csv');
    X = dataTrain(ind_prec1:n1_prec*i,1:2); y = dataTrain(ind_prec1:n1_prec*i,3);
    y = (y-min(y))/(max(y)-min(y));
    xtest = dataTest(ind_prec2:n2_prec*i,1:2); ytest = dataTest(ind_prec2:n2_prec*i,3);
    ytest = (ytest-min(ytest))/(max(ytest)-min(ytest));
    uqlab;
    MetaOpts.Type = 'Metamodel';
    MetaOpts.MetaType = 'Kriging';
    MetaOpts.ExpDesign.Sampling = 'User';
    MetaOpts.ExpDesign.X = X;
    MetaOpts.ExpDesign.Y = y;
    myKriging = uq_createModel(MetaOpts);
    yhat_uqlab = uq_evalModel(myKriging, xtest);
    rmse_precip_uqlab(i) = sqrt(mean((ytest-yhat_uqlab).^2));
    time_precip_uqlab(i) = toc;
    ind_prec1 = n1_prec*i+1;
    ind_prec2 = n2_prec*i+1;
end

%% Dataset-3 (Global Earthquake)
% SpatCorrNPR
time_eq_spatcorr = nan(1,iter);
n1_eq = 2309; n2_eq = 989; ind_eq1 = 1; ind_eq2 = 1;
rmse_eq_spatcorr = nan(1,iter);
for i=1:iter
    tic
    dataTrain = readmatrix('Train_earthquake.csv');
    dataTest = readmatrix('Test_earthquake.csv');
    X = dataTrain(ind_eq1:n1_eq*i,1:2); y = dataTrain(ind_eq1:n1_eq*i,3);
    y = (y-min(y))/(max(y)-min(y));
    xtest = dataTest(ind_eq2:n2_eq*i,1:2); ytest = dataTest(ind_eq2:n2_eq*i,3);
    ytest = (ytest-min(ytest))/(max(ytest)-min(ytest));
    H = [.0255 0.001; 0.001 0.008];
    yhat_spatcorr = locpol_spatcorrnpr(X,y,H,'Epanechnikov',xtest);
    rmse_eq_spatcorr(i) = sqrt(mean((ytest-yhat_spatcorr).^2));
    time_eq_spatcorr(i) = toc;
    ind_eq1 = n1_eq*i+1;
    ind_eq2 = n2_eq*i+1;
end

% ooDACE
time_eq_oodace = nan(1,iter);
n1_eq = 2309; n2_eq = 989; ind_eq1 = 1; ind_eq2 = 1;
rmse_eq_oodace = nan(1,iter);
for i=1:iter
    tic
    dataTrain = readmatrix('Train_earthquake.csv');
    dataTest = readmatrix('Test_earthquake.csv');
    X = dataTrain(ind_eq1:n1_eq*i,1:2); y = dataTrain(ind_eq1:n1_eq*i,3);
    y = (y-min(y))/(max(y)-min(y));
    xtest = dataTest(ind_eq2:n2_eq*i,1:2); ytest = dataTest(ind_eq2:n2_eq*i,3);
    ytest = (ytest-min(ytest))/(max(ytest)-min(ytest));
    k = oodacefit(X,y);
    [yhat_oodace, ~] = k.predict(xtest);   
    rmse_eq_oodace(i) = sqrt(mean((ytest-yhat_oodace).^2));
    time_eq_oodace(i) = toc;
    ind_eq1 = n1_eq*i+1;
    ind_eq2 = n2_eq*i+1;
end

% UQLab
time_eq_uqlab = nan(1,iter);
n1_eq = 2309; n2_eq = 989; ind_eq1 = 1; ind_eq2 = 1;
rmse_eq_uqlab = nan(1,iter);
for i=1:iter
    tic
    dataTrain = readmatrix('Train_earthquake.csv');
    dataTest = readmatrix('Test_earthquake.csv');
    X = dataTrain(ind_eq1:n1_eq*i,1:2); y = dataTrain(ind_eq1:n1_eq*i,3);
    y = (y-min(y))/(max(y)-min(y));
    xtest = dataTest(ind_eq2:n2_eq*i,1:2); ytest = dataTest(ind_eq2:n2_eq*i,3);
    ytest = (ytest-min(ytest))/(max(ytest)-min(ytest));
    uqlab;
    MetaOpts.Type = 'Metamodel';
    MetaOpts.MetaType = 'Kriging';
    MetaOpts.ExpDesign.Sampling = 'User';
    MetaOpts.ExpDesign.X = X;
    MetaOpts.ExpDesign.Y = y;
    myKriging = uq_createModel(MetaOpts);
    yhat_uqlab = uq_evalModel(myKriging, xtest);   
    rmse_eq_uqlab(i) = sqrt(mean((ytest-yhat_uqlab).^2));
    time_eq_uqlab(i) = toc;
    ind_eq1 = n1_eq*i+1;
    ind_eq2 = n2_eq*i+1;
end


%% Display Results
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


disp("%%%%%%%%%%%%%%%%%%%%%%%% Table 7 Results %%%%%%%%%%%%%%%%%%%%%%%%")
disp('RMSE for Dataset-1 (Column 1 of Table 7):')
fprintf('Texas_Wildfire-RMSE-SpatCorrNPR: %0.3f(%0.3f)\n', mean(rmse_wildfire_spatcorr), std(rmse_wildfire_spatcorr));
fprintf('Texas_Wildfire-RMSE-ooDACE: %0.3f(%0.3f)\n', mean(rmse_wildfire_oodace), std(rmse_wildfire_oodace));
fprintf('Texas_Wildfire-RMSE-UQLab: %0.3f(%0.3f)\n', mean(rmse_wildfire_uqlab), std(rmse_wildfire_uqlab));

disp('RMSE for Dataset-2 (Column 2 of Table 7):')
fprintf('usPrecipitation-RMSE-SpatCorrNPR: %0.3f(%0.3f)\n', mean(rmse_precip_spatcorr), std(rmse_precip_spatcorr));
fprintf('usPrecipitation-RMSE-ooDACE: %0.3f(%0.3f)\n', mean(rmse_precip_oodace), std(rmse_precip_oodace));
fprintf('usPrecipitation-RMSE-UQLab: %0.3f(%0.3f)\n', mean(rmse_precip_uqlab), std(rmse_precip_uqlab));

disp('RMSE for Dataset-3 (Column 3 of Table 7):')
fprintf('Global_Earthquake-RMSE-SpatCorrNPR: %0.3f(%0.3f)\n', mean(rmse_eq_spatcorr), std(rmse_eq_spatcorr));
fprintf('Global_Earthquake-RMSE-ooDACE: %0.3f(%0.3f)\n', mean(rmse_eq_oodace), std(rmse_eq_oodace));
fprintf('Global_Earthquake-RMSE-UQLab: %0.3f(%0.3f)\n', mean(rmse_eq_uqlab), std(rmse_eq_uqlab));
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")


disp("%%%%%%%%%%%%%%%%%%%%%%%% Table 8 Results %%%%%%%%%%%%%%%%%%%%%%%%")
disp('Runtimes for Dataset-1 (Column 1 of Table 8):')
fprintf('Texas_Wildfire-Runtime-SpatCorrNPR: %0.3f(%0.3f)\n', mean(time_wildfire_spatcorr), std(time_wildfire_spatcorr));
fprintf('Texas_Wildfire-Runtime-ooDACE: %0.3f(%0.3f)\n', mean(time_wildfire_oodace), std(time_wildfire_oodace));
fprintf('Texas_Wildfire-Runtime-UQLab: %0.3f(%0.3f)\n', mean(time_wildfire_uqlab), std(time_wildfire_uqlab));

disp('Runtimes for Dataset-2 (Column 2 of Table 8):')
fprintf('usPrecipitation-Runtime-SpatCorrNPR: %0.3f(%0.3f)\n', mean(time_precip_spatcorr), std(time_precip_spatcorr));
fprintf('usPrecipitation-Runtime-ooDACE: %0.3f(%0.3f)\n', mean(time_precip_oodace), std(time_precip_oodace));
fprintf('usPrecipitation-Runtime-UQLab: %0.3f(%0.3f)\n', mean(time_precip_uqlab), std(time_precip_uqlab));

disp('Runtimes for Dataset-3 (Column 3 of Table 8):')
fprintf('Global_Earthquake-Runtime-SpatCorrNPR: %0.3f(%0.3f)\n', mean(time_eq_spatcorr), std(time_eq_spatcorr));
fprintf('Global_Earthquake-Runtime-ooDACE: %0.3f(%0.3f)\n', mean(time_eq_oodace), std(time_eq_oodace));
fprintf('Global_Earthquake-Runtime-UQLab: %0.3f(%0.3f)\n', mean(time_eq_uqlab), std(time_eq_uqlab));
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

% ********************************************************************** %
