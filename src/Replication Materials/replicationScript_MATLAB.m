% This script contains all the code snippets included in the manuscript
% titled 'SpatCorrNPR: A MATLAB toolbox for Local Polynomial Regression 
% with Correlated Errors', describing the use of various modules of the
% SpatCorrNPR toolbox and replication of the practical examples.
%
% Authors: Mehnuma Tabassum, Kris De Brabanter
% Email: mehnuma@iastate.edu

% ********************************************************************** %

%% Codes for the Section "SpatCorrNPR: The MATLAB Toolbox"
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
% Plot Raw Data (Figure 7)
% To get a scatterplot of 1-D X, the following code is used:
plotRawData(X_1d,y)   % Figure 7(a)
colormap(brewermap([],"Spectral"))
xlabel('x1')
ylabel('y')
saveas(gcf, 'fig7a.png')
close

% For 2-D X, the following line produces a 3-D scatterplot.
plotRawData(X_2d,y)   % Figure 7(b)
colormap(brewermap([],"Spectral"))
xlabel('x1')
ylabel('x2')
zlabel('y')
saveas(gcf, 'fig7b.png')
close


% Plot Objective Function
% The following snippet provide the syntax for plotting the objective
% functions with Epanechnikov kernel for 1-D X and 2-D X,
% respectively.

% Figure 8
plotObjectiveFunction(X_1d, y,[], [0.01, 0.5], 'Epanechnikov', ...
    'Unadjusted GCV', [])   % Figure 8(a)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig8a.png')
close

plotObjectiveFunction(X_1d, y,[], [0.01, 0.5], 'Epanechnikov', ...
    'Corrected and Estimated GCV (GCVce)', [])   % Figure 8(b)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig8b.png')
close

plotObjectiveFunction(X_1d, y,[], [0.01, 5], 'Epanechnikov', ...
    'Residual Sum of Squares (RSS)', [])   % Figure 8(c)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig8c.png')
close

plotObjectiveFunction(X_1d, y, m, [0.01, 5], 'Epanechnikov', ...
    'Asymptotic Squared Error (ASE)', [])    % Figure 8(d)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig8d.png')
close

% Figure 9
plotObjectiveFunction(X_2d, y,[], [0.05 0.05; 2 2], 'Epanechnikov', ...
    'Unadjusted GCV', [])   % Figure 9(a)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig9a.png')
close

plotObjectiveFunction(X_2d, y,[], [0.05 0.05; 2 2], 'Epanechnikov', ...
    'Corrected and Estimated GCV (GCVce)', [])   % Figure 9(b)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig9b.png')
close

plotObjectiveFunction(X_2d, y,[], [0.05 0.05; 2 2], 'Epanechnikov', ...
    'Residual Sum of Squares (RSS)', [])   % Figure 9(c)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig9c.png')
close

plotObjectiveFunction(X_2d, y, m, [0.05 0.05; 2 2], 'Epanechnikov', ...
    'Asymptotic Squared Error (ASE)', [])    % Figure 9(d)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig9d.png')
close

% % The syntax for using the Corrected GCV (GCVc) with Epanechnikov kernel 
% % is as follows:
% 
% RMatrix = readmatrix('trueR.csv');
% plotObjectiveFunction(X_1d, y,[], [0.01, 0.5], 'Epanechnikov', ...
%     'Corrected GCV (GCVc)', RMatrix)    % for 1-D X
% plotObjectiveFunction(X_2d, y,[], [0.05 0.05; 2.5 2.5], 'Epanechnikov', ...
%     'Corrected GCV (GCVc)', RMatrix)    % for 2-D X

% Note: These plots are not shown in the manuscript. The plots shown in the
% manuscript are common to the app and the command line toolbox. GCVc
% plotting is only available in the command line toolbox. 


% Plot Smoothed Function (Figure 10)
% Below is the syntax for plotting the smoothed function for 1-D X and 2-D
% X, respectively.

plotSmoothedFunction(X_1d, y, 0.0824, 'Epanechnikov', []) % Figure 10(a)
colormap(brewermap([],"Spectral"))
xlabel('x1')
ylabel('y')
saveas(gcf, 'fig10a.png')
close

plotSmoothedFunction(X_2d, y, [0.233 0.058; 0.058 0.311], 'Epanechnikov', [])
                                                          % Figure 10(b)
colormap(brewermap([],"RdYlGn"))
view([-24, 23])
xlabel('x1')
ylabel('x2')
zlabel('y')                                                        
saveas(gcf, 'fig10b.png')
close


%% Module-5: Covariance Estimation
% Here are some example codes to produce semivariograms and calculate the
% covariance matrix of the local linear regression residuals.

% Warning: These plots might need some time to load compared to MATLAB.

% Figure 11
CMat_exp = covarianceEstimation(X_2d, y, [0.233 0.058; 0.058 0.311], ...
    'Epanechnikov', [], [], [], [], [], [], 'Exponential', [], []);
     % Model: Exponential, Figure 11(a)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig11a.png')
close

CMat_gau = covarianceEstimation(X_2d, y, [0.233 0.058; 0.058 0.311], ...
    'Epanechnikov', [], [], [], [], [], [], 'Gaussian', [], []);
     % Model: Gaussian, Figure 11(b)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig11b.png')
close


%% Module-6: Geoplotting (ONLY available in MATLAB)
% Example-1: The first geoPlot is for the US precipitation data. The 
% following codes are used:

geoData = readmatrix('usPrecipitation.csv');
X = geoData(:,1:2); y = geoData(:,3);
geoPlot(X, y, [3.615 0.0026; 0.0026 2.5873], 'Epanechnikov', [], ...
    'Contiguous US', [], [], 'Positive')  % Figure 12(a)
colormap(brewermap([],"RdBu"))
saveas(gcf, 'fig12a.png')
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
    'other')   % Figure 12(b)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig12b.png')
close

%% Codes for the Section "Some Practical Examples"
% In this section, we implement the method by De Brabanter (2018) on some 
% real-life datasets that are not necessarily spatially correlated. We used 
% the birthrate dataset, the Beluga dataset, and the NBA dataset. These
% datasets have been collected from the Stat-LSSVM Package provided by De
% Brabanter et al. (2013). The analysis primarily includes bandwidth 
% optimization and relevant plotting operations.


%% Practical Example-1: The "Birthrate" Dataset
% The "birthrate" dataset contains the monthly birth rate in the US from 
% January 1940 to December 1947. The dataset includes 96 observations with
% two variables (year and birthrate).

% Reading the birthrate dataset
birthData = readmatrix('birthrate.csv');
X_birth = birthData(:,1);
y_birth = birthData(:,2);

% Plot raw birthrate data
plotRawData(X_birth,y_birth)   % Figure 13(a)
colormap(brewermap([],"Spectral"))
box on
xlabel('Year')
ylabel('BirthRate')
saveas(gcf, 'fig13a.png')
close

% Plot objective functions (birthrate data)
plotObjectiveFunction(X_birth, y_birth, [], [0.055, 30], 'Epanechnikov', ...
    'Residual Sum of Squares (RSS)', [])  % Figure 13(b)
colormap(brewermap([],"Spectral"))
legend('Location','southeast')
saveas(gcf, 'fig13b.png')
close


% Bandwidth optimization of birthrate data (Table 6, Column 3)
hoptGCV_birth = bwOptimization(X_birth, y_birth, [], 'Epanechnikov', ...
    'Unadjusted GCV', 0.055, 30, []);

hoptRSS_birth = bwOptimization(X_birth, y_birth, [], 'Epanechnikov', ...
    'Residual Sum of Squares (RSS)', 1.1, 30, []);


disp("%%%%%%%%%%%%%%%%%%%%%%%% Table 6 Results %%%%%%%%%%%%%%%%%%%%%%%%")
disp('Optimal Bandwidth Values for Birthrate Data (Column 3 of Table 6):')

disp('Optimal Bandwidth, Unadjusted GCV: '),          disp(hoptGCV_birth)
disp('Optimal Bandwidth, Residual Sum of Squares: '), disp(hoptRSS_birth)

disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")


% Plot smoothed function (birthrate data)
plotSmoothedFunction(X_birth, y_birth, 0.0844, 'Epanechnikov', [])
                     % Figure 13(c)
colormap(brewermap([],"Spectral"))                     
xlabel('Year')
ylabel('BirthRate')
legend('Location','northwest')
saveas(gcf, 'fig13c.png')
close

plotSmoothedFunction(X_birth, y_birth, 1.1, 'Epanechnikov', [])
                     % Figure 13(d)
colormap(brewermap([],"Spectral"))
xlabel('Year')
ylabel('BirthRate')
saveas(gcf, 'fig13d.png')
close


%% Practical Example-2: The "Beluga" Dataset
% The "Beluga" dataset contains the nursing times (in seconds) of newborn 
% Beluga whale calves. The dataset has two columns: periods (six-hour time 
% periods), and nursing time, and the total number of observations is 228.

% Reading the Beluga dataset
belugaData = readmatrix('beluga.csv');
X_beluga = belugaData(:,1);
y_beluga = belugaData(:,2);

% Plot raw beluga data
plotRawData(X_beluga, y_beluga)   % Figure 14(a)
colormap(brewermap([],"Spectral"))
box on
xlabel('Period')
ylabel('NursingTime')
saveas(gcf, 'fig14a.png')
close

% Plot objective function (beluga data)
plotObjectiveFunction(X_beluga, y_beluga, [], [0.055, 500], 'Epanechnikov', ...
    'Residual Sum of Squares (RSS)', []) % Figure 14(b)
colormap(brewermap([],"Spectral"))
saveas(gcf, 'fig14b.png')
close

% Bandwidth optimization of beluga data (Table 7, Column 3)
hoptGCV_beluga = bwOptimization(X_beluga, y_beluga, [], 'Epanechnikov', ...
                 'Unadjusted GCV', 0.001, 30, []);

hoptRSS_beluga = bwOptimization(X_beluga, y_beluga, [], 'Epanechnikov', ...
                 'Residual Sum of Squares (RSS)', 17.5, 30, []);

disp("%%%%%%%%%%%%%%%%%%%%%%%% Table 7 Results %%%%%%%%%%%%%%%%%%%%%%%%")
disp('Optimal Bandwidth Values for Beluga Data (Column 3 of Table 7):')

disp('Optimal Bandwidth, Unadjusted GCV: '),          disp(hoptGCV_beluga)
disp('Optimal Bandwidth, Residual Sum of Squares: '), disp(hoptRSS_beluga)

disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

% Plot smoothed function (beluga data)
plotSmoothedFunction(X_beluga, y_beluga, 3.21, 'Epanechnikov', ...
    [])  % Figure 14(c)
colormap(brewermap([],"Spectral"))
xlabel('Period')
ylabel('NursingTime')
saveas(gcf, 'fig14c.png')
close

plotSmoothedFunction(X_beluga, y_beluga, 17.5, 'Epanechnikov', ...
    [])  % Figure 14(d)
colormap(brewermap([],"Spectral"))
xlabel('Period')
ylabel('NursingTime')
saveas(gcf, 'fig14d.png')
close


%% Practical Example-3: The "NBA" Dataset
% The NBA dataset contains the mean points scored per minute based on 
% number of minutes played in each game and the player heights (in 
% centimeters). There are 96 observations in this dataset.

% Reading the nba dataset
nbaData = readmatrix('nba.csv');
X_nba = nbaData(:,1:2);
y_nba = nbaData(:,3);

% Plot raw nba data
plotRawData(X_nba,y_nba)   % Figure 15(a)
xlabel('MPG')
ylabel('Height')
zlabel('MPS')
saveas(gcf, 'fig15a.png')
close

% Plot objective function (nba data)
plotObjectiveFunction(X_nba, y_nba, [], [0.055 0.055; 30 30], ...
    'Epanechnikov', 'Residual Sum of Squares (RSS)', [])   % Figure 15(b)
colormap(brewermap([],"Spectral"))
view([65,21])
saveas(gcf, 'fig15b.png')
close

% Bandwidth optimization of nba data (Table 8, Column 3)
optH_GCV = bwOptimization(X_nba, y_nba, [], 'Epanechnikov', ...
           'Unadjusted GCV', [0.055 0.001; 0.001 0.055], [5 0.1; 0.1 5], ...
           []);

optH_RSS = bwOptimization(X_nba, y_nba, [], 'Epanechnikov', ...
           'Residual Sum of Squares (RSS)', [18.5 0.001; 0.001 26], ...
           [20.5 0.1; 0.1 28], []);

% The optH_RSS result is very sensitive and changes everytime we run this.
% This is perhaps due to using the CSA algorithm to set the seed for the
% fminsearch algorithm. We picked one of the results to report.

disp("%%%%%%%%%%%%%%%%%%%%%%%% Table 8 Results %%%%%%%%%%%%%%%%%%%%%%%%")
disp('Optimal Bandwidth Values for NBA Data (Column 3 of Table 8):')

disp('Optimal Bandwidth, Unadjusted GCV: '),          disp(optH_GCV)
disp('Optimal Bandwidth, Residual Sum of Squares: '), disp(optH_RSS)

disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

% Plot smoothed function (nba data)
plotSmoothedFunction(X_nba, y_nba, [4.21 0.001; 0.001 5.0], 'Epanechnikov', [])
                     % Figure 15(c)
colormap(brewermap([],"Spectral"))
xlabel('MPG')
ylabel('Height')
zlabel('MPS')
saveas(gcf, 'fig15c.png')
close

plotSmoothedFunction(X_nba, y_nba, [20.5 0.1; 0.1 26], 'Epanechnikov', [])
                     % Figure 15(d)
colormap(brewermap([],"Spectral"))
xlabel('MPG')
ylabel('Height')
zlabel('MPS')
saveas(gcf, 'fig15d.png')
close


% ********************************************************************** %

%% References
% De Brabanter, K., Cao, F., Gijbels, I., & Opsomer, J. (2018). Local 
% polynomial regression with correlated errors in random design and unknown 
% correlation structure. Biometrika, 105(3), 681-690.
%
% De Brabanter K, Suykens J, De Moor B (2013). “Nonparametric Regression 
% via StatLSSVM." Journal of Statistical Software, 55, 1–21.


