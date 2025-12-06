% This script implements the sensitivity study presented in the Supplementary 
% Materials for the manuscript titled 'SpatCorrNPR: A MATLAB toolbox for Local 
% Polynomial Regression with Correlated Errors'.
%
% Authors: Mehnuma Tabassum, Kris De Brabanter
% Email: mehnuma@iastate.edu
% ********************************************************************** %

%% Code for the Sensitivity Analysis
dataWithM = readmatrix('testDataWithM.csv'); 
X_1d = dataWithM(:,1);   % for 1-D input 
X_2d = dataWithM(:,1:2); % for 2-D input 
y = dataWithM(:,3);      % The response
m = dataWithM(:,4); 

optimal_h1_epa = bwOptimization(X_1d, y, [], 'Epanechnikov', ...
    'Residual Sum of Squares (RSS)', 0.001, 0.03, []);
optimal_h2_epa = bwOptimization(X_1d, y, [], 'Epanechnikov', ...
    'Residual Sum of Squares (RSS)', 0.009, 0.02, []);
optimal_h3_epa = bwOptimization(X_1d, y, [], 'Epanechnikov', ...
    'Residual Sum of Squares (RSS)', 0.1, 1, []);
optimal_h1_gau = bwOptimization(X_1d, y, [], 'Gaussian', ...
    'Residual Sum of Squares (RSS)', 0.001, 0.03, []);
optimal_h2_gau = bwOptimization(X_1d, y, [], 'Gaussian', ...
    'Residual Sum of Squares (RSS)', 0.009, 0.02, []);
optimal_h3_gau = bwOptimization(X_1d, y, [], 'Gaussian', ...
    'Residual Sum of Squares (RSS)', 0.1, 1, []);


optimal_H1_epa = bwOptimization(X_2d, y, [], 'Epanechnikov', ...
    'Residual Sum of Squares (RSS)', [0.2, 0.001; 0.001 0.2], ...
    [0.5, 0.1; 0.1 0.5], []);
optimal_H2_epa = bwOptimization(X_2d, y, [], 'Epanechnikov', ...
    'Residual Sum of Squares (RSS)', [0.25, 0.01; 0.01 0.25], ...
    [0.8, 0.05; 0.05 0.8], []);
optimal_H3_epa = bwOptimization(X_2d, y, [], 'Epanechnikov', ...
    'Residual Sum of Squares (RSS)', [0.32, 0.05; 0.05 0.32], ...
    [0.45, 0.1; 0.1 0.45], []);
optimal_H1_gau = bwOptimization(X_2d, y, [], 'Gaussian', ...
    'Residual Sum of Squares (RSS)', [0.2, 0.001; 0.001 0.2], ...
    [0.5, 0.1; 0.1 0.5], []);
optimal_H2_gau = bwOptimization(X_2d, y, [], 'Gaussian', ...
    'Residual Sum of Squares (RSS)', [0.25, 0.01; 0.01 0.25], ...
    [0.8, 0.05; 0.05 0.8], []);
optimal_H3_gau = bwOptimization(X_2d, y, [], 'Gaussian', ...
    'Residual Sum of Squares (RSS)', [0.32, 0.05; 0.05 0.32], ...
    [0.45, 0.1; 0.1 0.45], []);


disp("%%%%%%%%%%%%%%%%%%%%%%%% Supplemental Table 1 Results %%%%%%%%%%%%%%%%%%%%%%%%")
disp('1-D Optimal Bandwidth (Column 3 of Table 1):')
fprintf('Kernel: Epanechnikov (1D bandwidth-1): %0.4f\n', optimal_h1_epa);
fprintf('Kernel: Epanechnikov (1D bandwidth-2): %0.4f\n', optimal_h2_epa);
fprintf('Kernel: Epanechnikov (1D bandwidth-3): %0.4f\n', optimal_h3_epa);
fprintf('Kernel: Gaussian (1D bandwidth-1): %0.4f\n', optimal_h1_gau);
fprintf('Kernel: Gaussian (1D bandwidth-2): %0.4f\n', optimal_h2_gau);
fprintf('Kernel: Gaussian (1D bandwidth-3): %0.4f\n', optimal_h3_gau);

disp('2-D Optimal Bandwidth (Column 5 of Table 1):')
fprintf('Kernel: Epanechnikov (2D bandwidth-1): %0.4f\n', optimal_H1_epa);
fprintf('Kernel: Epanechnikov (2D bandwidth-2): %0.4f\n', optimal_H2_epa);
fprintf('Kernel: Epanechnikov (2D bandwidth-3): %0.4f\n', optimal_H3_epa);
fprintf('Kernel: Gaussian (2D bandwidth-1): %0.4f\n', optimal_H1_gau);
fprintf('Kernel: Gaussian (2D bandwidth-2): %0.4f\n', optimal_H2_gau);
fprintf('Kernel: Gaussian (2D bandwidth-3): %0.4f\n', optimal_H3_gau);
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")


