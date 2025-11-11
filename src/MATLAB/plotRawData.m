function plotRawData(X,y)
% Creates a 2-D or 3-D scatterplot
%
% Input Arguments:
% X : Input data (1-D or 2-D)
% y : Response
%
% Output:
% 2-D or 3-D scatterplot

% Authors: Mehnuma Tabassum, Kris De Brabanter
% Email: mehnuma@iastate.edu

[~, dimX] = size(X);

if dimX==1
    scatter(X, y, 40, 'r', 'filled')
    xlabel('Input')
    ylabel('Response')
    title('Scatter Plot of Input Data')
    set(gca, 'fontname', 'Arial', 'fontsize', 14)
elseif dimX==2
    scatter3(X(:,1), X(:,2), y, 40, 'r', 'filled')
    xlabel('Input-1')
    ylabel('Input-2')
    zlabel('Response')
    title('Scatter Plot of Input Data')
    set(gca, 'fontname', 'Arial', 'fontsize', 14)
else
    errordlg('Please Provide 1-D or 2-D Input Data')
end