function plotSmoothedFunction(X,y,bw,kernel,grain)
% Plots the smoothed function based on bandwidth and kernel
%
% Input Arguments
% X :     Input data (Typically, 2-D Spatial dataset)
% y :     Response
% bw:     Bandwidth (scalar for 1-D X, 2-by-2 symmetric matrix for 2-D X)
% kernel: Kernel used for local polynomial regression
% grain:  Plot resoultion
%
% Output
% 2-D or 3-D smoothed regression plot

% Authors: Mehnuma Tabassum, Kris De Brabanter
% Email: mehnuma@iastate.edu

if isempty(bw)
    errordlg('Please Provide a Bandwidth')
    return
end

if isempty(kernel)
    errordlg('Please Specify a Kernel')
    return
end

if isempty(grain)
    grain = 50;
end
[~, dimX] = size(X);
if dimX==1
    if ~isscalar(bw)
        errordlg('Please Provide a Scalar Bandwidth for 1-D Input');
    else
        xt = linspace(min(X), max(X),grain*grain)';
        [yhat,~] = locpol(X,y,bw,kernel,xt);
        plot(xt,yhat, 'b-.', 'LineWidth', 5)
        colormap jet
        shading interp
        hold on
        scatter(X,y, 36, 'Marker', 'o', 'MarkerFaceColor',[0 0 0], 'MarkerEdgeColor',[0 0 0])
        set(gca, 'fontname', 'Arial', 'FontSize', 14)
        legend('Estimated', 'Observed Data')
        xlabel('Input')
        ylabel('Response')
        title('Local Linear Smoothing')
    end
elseif dimX==2
    [rowBW, colBW] = size(bw);
    if ~(rowBW==2) || ~(colBW==2)
        errordlg('Please Provide a 2x2 Bandwidth Matrix for 2-D Input');
    else
        xmin1 = min(X(:,1)); xmax1 = max(X(:,1)); xmin2 = min(X(:,2)); xmax2 = max(X(:,2));
        if xmin1<0
            xmin1=1.05*xmin1;
        else
            xmin1 = 0.98*xmin1;
        end
        if xmax1>0
            xmax1=1.05*xmax1;
        else
            xmax1 = 0.98*xmax1;
        end
        if xmin2<0
            xmin2=1.05*xmin2;
        else
            xmin2 = 0.98*xmin2;
        end
        if xmax2>0
            xmax2=1.05*xmax2;
        else
            xmax2 = 0.98*xmax2;
        end
        w = waitbar(0.5,'Fetching your plot ...', 'Name', 'Plots');
        pause(0.5)
        xrange1 = xmin1:(xmax1-xmin1)/grain:xmax1;
        xrange2 = xmin2:(xmax2-xmin2)/grain:xmax2;
        [XX,YY] = meshgrid(xrange1,xrange2);
        Xt = [reshape(XX,numel(XX),1) reshape(YY,numel(YY),1)];
        [yhat,~] = locpol(X,y,bw,kernel,Xt);
        yhat = reshape(yhat,size(XX,1),size(XX,2));
        waitbar(0.98, w, 'Almost done ...');
        pause(1)
        close(w)
        h1 = surfc(XX,YY,yhat);
        colormap jet
        shading interp
        hold on
        h2 = scatter3(X(:,1),X(:,2),y,36,'Marker','o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0,0,0]);
        set(gca, 'fontname', 'Arial', 'FontSize', 14)
        axP = get(gca,'Position');
        legend([h1(2), h2], {'Estimated', 'Observed Data'},'Location','NorthEast')
        set(gca, 'Position', axP)
        xlabel('Input-1')
        ylabel('Input-2')
        zlabel('Response')
        title('Local Linear Smoothing')
    end
else
    errordlg('Input must be 1-D or 2-D!');
    return
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
