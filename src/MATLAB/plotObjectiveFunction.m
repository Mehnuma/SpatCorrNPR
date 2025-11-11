function plotObjectiveFunction(X,y,m,bwlimit,kernel,plotObjective,R_true)
% Plots the selected objective function based on the objective and the
% bandwidth limits
% X :            Input data (1-D or 2-D)
% y :            Response
% m :            True function value (optional; required for ASE)
% bwlimit:       Bandwidth limit (vector for 1-D X; 2-by-2 matrix for 2-D 
%                X where first row is for lower limits and second row is 
%                for upper limits)
% kernel:        Kernel used for local polynomial regression
% plotObjective: Objective function to plot
% R_true:        True correlation matrix (Only required for GCVc)
%
% Output:
% 2-D or 3-D objective function plot

% Authors: Mehnuma Tabassum, Kris De Brabanter
% Email: mehnuma@iastate.edu

[~, dimX] = size(X);

if isempty(bwlimit)
    errordlg('Please Provide the Bandwidth Limit')
    return
end
if isempty(kernel)
    errordlg('Please Specify a Kernel')
    return
end
if isempty(plotObjective)
    errordlg('Please Specify the Objective')
    return
else
    if strcmp(plotObjective, 'Corrected GCV (GCVc)') && isempty(R_true)
        errordlg('Please Provide the True Correlation Matrix')
        return
    end
end

grain = 30;
if dimX==1
    if (~isvector(bwlimit) || ~(length(bwlimit)==2))
        errordlg('Please Provide a Scalar Bandwidths for 1-D Input');
        return
    else
        h = linspace(bwlimit(1),bwlimit(2),grain)';
        cvhat = zeros(size(h));
        for j=1:size(h,1)
            cvhat(j) = objFunction(X,y,m,h(j),kernel,plotObjective,R_true);
        end
        figure()
        plot(h,cvhat, 'b-.', 'LineWidth', 3)
        xlabel('Bandwidth')
        ylabel('Objective')
        xlim([bwlimit(1),bwlimit(2)])
        title('Objective Function Plot')
        legend(plotObjective)
        set(gca, 'FontName', 'Arial', 'fontsize', 14)
        colormap jet
        shading interp
    end
elseif dimX==2
    [rowH, colH] = size(bwlimit);
    if ~(rowH==2) || ~(colH==2)
        errordlg('Please Provide a 2x2 Bandwidth Matrix for the Limits');
        return
    else
        w = waitbar(0.5,'Fetching your plot ...', 'Name', 'Plots');
        pause(0.5)
        [h1, h2] = meshgrid(bwlimit(1,1):0.5:bwlimit(2,1), bwlimit(1,2):0.5:bwlimit(2,2));
        H = [reshape(h1,numel(h1),1) reshape(h2,numel(h2),1)];
        cvhat = zeros(size(H,1),1);
        for j=1:size(H,1)
            cvhat(j) = objFunction(X,y,m,diag(H(j,:)),kernel,plotObjective,R_true);
        end
        cvhat = reshape(cvhat,size(h1,1),size(h2,2));
        waitbar(0.98, w, 'Almost done ...');
        pause(1)
        close(w)
        figure()
        surfc(h1,h2,cvhat)
        xlabel('Bandwidth-1')
        ylabel('Bandwidth-2')
        zlabel('Objective')
        title('Objective Function Plot')
        set(gca, 'fontname', 'Arial', 'FontSize', 14)
        colormap jet
        shading interp
    end
else
    errordlg('X needs to be either 1-D or 2-D');
    return
end

end

function cvhat = objFunction(X,y,m,H,kernel,objective,R_true)
[N, dimX] = size(X);
switch objective
    case 'Unadjusted GCV'
        [yhat, S] = locpol(X, y, H, kernel);
        num = y-yhat;
        den = 1 -(1/N)*trace(S);
        cvhat = mean((num/den).^2, 'all');

    case 'Corrected GCV (GCVc)'
        [yhat, S] = locpol(X, y, H, kernel);
        num = y-yhat;
        den = 1 -(1/N)*trace(S*R_true);
        cvhat = mean((num/den).^2, 'all');

    case 'Corrected and Estimated GCV (GCVce)'
        if dimX==1
            [yHat, ~]= locpol(X,y,std(X),kernel);
        elseif dimX==2
            [yHat, ~]= locpol(X,y,[std(X(:,1)) 0; 0 std(X(:,2))],kernel);
        end
        res = y-yHat; D = pdist2(X,X);
        dk = (0.001+2*((1:30)-1)*0.005)'; [ndk, ~] = size(dk);
        gamma = zeros(ndk,1);
        for j = 1:ndk
            ind = D>=max(dk(j)-0.005, 0) & D<=dk(j)+0.005;
            [row, column] = find(triu(ind)==1);
            num_el = numel(row);
            rse = zeros(num_el, 1);
            for i = 1:num_el
                rse(i) = (res(column(i)) -res(row(i))).^2;
            end
            gamma(j) = (1/(2*num_el))*sum(rse);
        end
        sigmaSq = mean(res.^2);
        alphaJ = (1./(N*dk)).*(log(sigmaSq) -log(abs(sigmaSq-gamma)));
        alpha = mean(alphaJ); R_exp = exp(-(alpha*N)*D);

        [yhat, S]= locpol(X,y,H,kernel);
        num = y-yhat;
        den = 1 -(1/N)*trace(S*R_exp);
        cvhat = mean((num/den).^2, 'all');

    case 'Residual Sum of Squares (RSS)'
        if strcmp(kernel, 'Epanechnikov')
            if dimX==1
                H = H/2.57132;
            elseif dimX==2
                H = H/2.58;
            end
        elseif strcmp(kernel, 'Gaussian')
            if dimX==1
                H = H/1.16;
            elseif dimX==2
                H = H/1.14;
            end
        end
        yhat = locpol(X, y, H, 'Kernel with K(0)=0');
        cvhat = sum((y-yhat).^2);

    case 'Asymptotic Squared Error (ASE)'
        if isempty(m)
            errordlg('The true function, m must be supplied');
            cvhat=[];
            return
        else
            yhat = locpol(X, y, H, kernel);
            cvhat = sum((m-yhat).^2);
        end
end
end

%% Supporting Functions
function [yhat,L] = locpol(X,y,H,kernel,xt)
[n,d] = size(X);
opts.SYM = true;
if nargin < 5
    L = zeros(n,n);
    if d==1
        for i = 1:n
            Xx = [ones(n,1) X-X(i,:)];
            Wx = kfun1d(X,H,X(i,:),kernel);
            t = Xx'*diag(Wx);
            tt = linsolve(t*Xx + eye(d+1)*1e-6,t,opts);
            L(i,:) = tt(1,:);
        end
    elseif d==2
        for i = 1:n
            Xx = [ones(n,1) X-X(i,:)];
            Wx = kfun2d(X,H,X(i,:),kernel);
            t = Xx'*diag(Wx);
            tt = linsolve(t*Xx + eye(d+1)*1e-6,t,opts);
            L(i,:) = tt(1,:);
        end
    end
    yhat = L*y;

else
    [nt,d] = size(xt);
    L = zeros(nt,n);
    if d==1
        for i = 1:nt
            Xx = [ones(n,1) X-xt(i,:)];
            Wx = kfun1d(X,H,xt(i,:),kernel);
            t = Xx'*diag(Wx);
            tt = linsolve(t*Xx + eye(d+1)*1e-6,t,opts);
            L(i,:) = tt(1,:);
        end
    elseif d==2
        for i = 1:nt
            Xx = [ones(n,1) X-xt(i,:)];
            Wx = kfun2d(X,H,xt(i,:),kernel);
            t = Xx'*diag(Wx);
            tt = linsolve(t*Xx + eye(d+1)*1e-6,t,opts);
            L(i,:) = tt(1,:);
        end
    end
    yhat = L*y;
end
end

function W = kfun1d(X, H, xt, kernel)
opts.SYM = true;
arg = linsolve(H,(X-xt)',opts);
switch kernel
    case 'Epanechnikov'
        W = 0.75*(1-arg.^2).*(abs(arg)<=1);
    case 'Gaussian'
        W = exp(-arg.*arg/2)/sqrt(2*pi);
    case 'Kernel with K(0)=0'
        W = (2/sqrt(pi))*arg.^2.*exp(-arg.^2);
    otherwise, error('No other kernel supported.')
end
end

function W = kfun2d(X, H, xt, kernel)
opts.SYM = true;
arg = linsolve(H,(X-xt)',opts);
switch kernel
    case 'Epanechnikov'
        W = (2/pi)*max(1-diag(arg'*arg),0)/det(H);
    case 'Gaussian'
        W = (1/(2*pi))*(det(H))^(-1)*exp(-0.5*diag(arg'*arg));
    case 'Kernel with K(0)=0'
        W = (1/pi)*(det(H))^(-1)*diag(arg'*arg).*exp(-diag(arg'*arg));
    otherwise, error('No other kernel supported.')
end
end