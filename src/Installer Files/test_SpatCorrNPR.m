
% This is a test script for the SpatCorrNPR toolbox.
%
% Authors: Mehnuma Tabassum, Kris De Brabanter
% Email: mehnuma@iastate.edu


% Data Import
dataWithM = readmatrix('testDataWithM.csv'); 
X_1d = dataWithM(:,1);   % for 1-D input 
X_2d = dataWithM(:,1:2); % for 2-D input 
y = dataWithM(:,3);      % The response
m = dataWithM(:,4);      % The true function

% Plot Smmothed Function (2-D)
plotSmoothedFunction(X_2d, y, [0.233 0.058; 0.058 0.311], 'Epanechnikov', [])
                                                          % Figure 4(b)
xlabel('x1')
ylabel('x2')
zlabel('y')
view([-13.7619,29.5192])

