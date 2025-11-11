function err = evaluate_error(X,y,m,criterion,bw,kernel)
% Evaluates the error of the nonparametric regression estimator given the
% criterion
%
% Input Arguments:
% X :        Input data (1-D or 2-D)
% y :        Response
% m :        True function value (optional; required for ASE)
% criterion: Error criterion (PRESS or ASE)
% bw:        Bandwidth (scalar for 1-D X, 2-by-2 symmetric matrix for 2-D X)
% kernel:    Kernel used for local polynomial regression
%
% Output Argument:
% err:       Calcuated error metric

% Authors: Mehnuma Tabassum, Kris De Brabanter
% Email: mehnuma@iastate.edu

[nX, dimX] = size(X);

if isempty(criterion)
    errordlg('Please Provide an Error Criterion')
    return
end

if ~(strcmp(criterion, 'PRESS') || strcmp(criterion, 'ASE'))
    errordlg('Please use the Available Error Criteria')
    return
end

if isempty(bw)
    errordlg('Please Provide a Bandwidth')
    return
end

if dimX==1
    if ~isscalar(bw)
        errordlg('Please Provide a Scalar Bandwidth for 1-D Data');
        return
    end
elseif dimX==2
    [rowBW, colBW]= size(bw);
    if ~(rowBW==2) ||~(colBW==2)
        errordlg('Please Provide a 2x2 Symmetric Bandwidth Matrix for 2-D Input');
        return
    end
    if ~(bw(1,2)==bw(2,1))
        errordlg('Please Provide a 2x2 Symmetric Bandwidth Matrix for 2-D Input');
        return
    end
end

if isempty(kernel)
    errordlg('Please Specify a Kernel')
    return
end

switch criterion
    case 'PRESS'
        yEst = zeros(nX,1);
        ind = (1:nX)';
        for j=1:nX
            Xt = X(ind~=j,:);
            yt = y(ind~=j);
            yEst(j) = locpol(Xt, yt, bw, kernel, X(j, :));
        end
        err = sum((y-yEst).^2);

    case 'ASE'
        if isempty(m)
            errordlg('Please Supply the True Function, m for the ASE Criterion')
        else
            mhat = locpol(X, y, bw, kernel);
            err = mean((m-mhat).^2);
        end
end
end

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
    otherwise, error('No other kernel supported')
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
