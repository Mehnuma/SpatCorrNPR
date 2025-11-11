function CMat = covarianceEstimation(X,y,bw,kernel,nrbins,maxdist,a0,c0,n0,numobs,model,nu,processVariance)
% Calculates the covariance matrix of the residuals based on the fitted
% semivariogram model and plots the empirical semivariogram with fitted
% model
%
% Input Arguments
% X :      Input data (1-D or 2-D)
% y :      Response
% bw:      Bandwidth (scalar for 1-D X, 2-by-2 symmetric matrix for 2-D X)
% kernel:  Kernel used for local polynomial regression
% nrbins:  Number of bins for grouping
% maxdist: Maximum distance in the data
% a0 :     Initial range
% c0 :     Initial sill
% n0 :     Initial nugget
% numobs:  Number of observations per lag distance
% model:   Model to fit
% nu:      Parameter for the Matérn model
%
% Output:
% CMat:    Covariance matrix of the residuals based on the model chosen
% Empirical semivariogram plot with fitted model

% Authors: Mehnuma Tabassum, Kris De Brabanter
% Email: mehnuma@iastate.edu

[nX,~]= size(X);

[yhat,~] = locpol(X,y,bw,kernel);
residuals = y-yhat;

emp_svar = empirical_svar(X, y, bw, kernel, nrbins, maxdist);
if isempty(processVariance)
    processVariance = emp_svar.var;
end

fitted_svar = fit_svar(emp_svar.dist,emp_svar.gammaEst,a0,c0,n0,numobs,model,nu);

if ~isempty(fitted_svar)
    cormat = zeros(nX);
    for ix=1:nX
        for jx=1:nX
            xi = residuals(ix,:); xj = residuals(jx,:);
            if strcmp(model, 'Linear')
                cormat(ix,jx) = exp(-norm(xi-xj)/fitted_svar.range);
            elseif strcmp(model, 'Circular')
                cormat(ix,jx) = 1-(1-(2./pi)*acos(norm(xi-xj)/fitted_svar.range)+2*norm(xi-xj)/(pi*fitted_svar.range)*sqrt(1-(norm(xi-xj).^2)/(fitted_svar.range^2)));
            elseif strcmp(model, 'Spherical')
                cormat(ix,jx) = 1-((3*norm(xi-xj)/(2*fitted_svar.range))-1/2*(norm(xi-xj)./fitted_svar.range).^3);
            elseif strcmp(model, 'Pentaspherical')
                cormat(ix,jx) = 1-(15*norm(xi-xj)./(8*fitted_svar.range)-5/4*(norm(xi-xj)./fitted_svar.range)^3+3/8*(norm(xi-xj)/fitted_svar.range)^5);
            elseif strcmp(model, 'Cubic')
                cormat(ix,jx) = 1-(7*(norm(xi-xj)/fitted_svar.range)^2 -8.75*(norm(xi-xj)/fitted_svar.range)^3 +3.5*(norm(xi-xj)/fitted_svar.range)^5 -0.75*(norm(xi-xj)/fitted_svar.range)^7);
            elseif strcmp(model, 'Exponential')
                cormat(ix,jx) = exp(-norm(xi-xj)/fitted_svar.range);
            elseif strcmp(model, 'Gaussian')
                cormat(ix,jx) = exp(-(norm(xi-xj)).^2/fitted_svar.range^2);
            elseif strcmp(model, 'Matérn')
                cormat(ix,jx) = 1-(1-(1/((2^(nu-1))*gamma(nu)))*(norm(xi-xj)/fitted_svar.range).^nu.*besselk(nu,norm(xi-xj)/fitted_svar.range));
            end
        end
    end
    CMat = processVariance*cormat;
else
    CMat = [];
end
end

function S = empirical_svar(X, y, bw, kernel, nrbins, maxdist)
% Calculates the empirical semivariogram based on the local polynomial
% regression residuals
%
% Input Arguments:
% X :      Input data (1-D or 2-D)
% y :      Response
% bw:      Bandwidth (scalar for 1-D X, 2-by-2 symmetric matrix for 2-D X)
% kernel:  Kernel used for local polynomial regression
% nrbins:  Number of bins for grouping (default=25)
% maxdist: Maximum distance in the data
%
% Output Argument
% S :      A struct containing the lag distances, empirical semivariogram
%          values, and empirical sill

[~,dimX] = size(X);

if isempty(bw)
    errordlg('Please Provide a Bandwidth')
    return
else
    if dimX==1
        if ~isscalar(bw)
            errordlg('Please Provide a Scalar Bandwidth for 1-D Input');
        end
    elseif dimX==2
        [rowBW, colBW] = size(bw);
        if ~(rowBW==2) ||~(colBW==2)
            errordlg('Please provide a 2x2 symmetric bandwidth matrix for 2-D input');
            return
        end
        if ~(bw(1,2)==bw(2,1))
            errordlg('Please provide a 2x2 symmetric bandwidth matrix for 2-D input');
            return
        end
    else
        errordlg('Please Provide Either 1-D or 2-D Input')
    end
end

[yhat,~] = locpol(X,y,bw,kernel);
res = y-yhat;

D = (pdist2(X,X));
if isempty(maxdist)
    minx = min(X,[],1); maxx = max(X,[],1); maxd = sqrt(sum((maxx-minx).^2));
    maxdist = maxd/2;
end
if ~isempty(nrbins)
    dk = linspace(0,maxdist, nrbins+1);
else
    nrbins = 25;
    % dk = (0.001+2*((1:(nrbins+1))-1)*0.005)';
    dk = linspace(0,maxdist, nrbins+1)';
end
[ndk, ~] = size(dk);
gammaEst = zeros(ndk,1);
for j = 1:numel(dk)
    ind = D>=max(dk(j)-0.005, 0) & D<=dk(j)+0.005;
    [row, column] = find(triu(ind)==1);
    num_el = numel(row);
    rse = zeros(num_el, 1);
    for i = 1:num_el
        rse(i) = (res(column(i)) -res(row(i))).^2;
    end
    gammaEst(j) = (1/(2*num_el))*sum(rse);
end
S.dist = dk;
S.gammaEst = gammaEst;
S.var = max(gammaEst);
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

function S = fit_svar(h,gammaexp,a0,c0,n0,numobs,model,nu)
% Fits a parametric model to an empirical semivariogram
%
% Input Arguments:
% h :       Lag distances
% gammaexp: Empirical semivariogram values
% a0 :      Initial range
% c0 :      Initial sill
% n0 :      Initial nugget
% numobs:   Number of observations per lag distance
% model:    Model to fit
% nu:       Parameter for the Matérn model
%
% Output Parameters
% S :      A struct containing the estimated parameters

% check input arguments
if isempty(h)
    errordlg('Please Provide the Lag Distances')
    S = [];
    return
end
if isempty(gammaexp)
    errordlg('Please Provide the Empirical Semivariogram Values')
    S = [];
    return
end
if isempty(a0)
    a0 = max(h)*2/3;
end
if isempty(c0)
    c0 = max(gammaexp);
end
% nugget variance
if isempty(n0)
    nugget = false;
    funnugget = @(b) 0;
    b0 = [a0 c0 0];
else
    nugget = true;
    funnugget = @(b) b(3);
    b0 = [a0 c0 n0];
end
%model
if isempty(model)
    errordlg('Please Specify a Model to fit the Empirical Semivariogram')
    S = [];
    return
end

% variogram function definitions
switch model
    case 'Linear'
        type = 'bounded';
        func = @(b,h) b(2).*(h./b(1));
    case 'Circular'
        type = 'bounded';
        func = @(b,h) b(2).*(1-(2./pi).*acos(h./b(1))+2*h./(pi.*b(1)).*sqrt(1-(h.^2)./(b(1).^2)));
    case 'Spherical'
        type = 'bounded';
        func = @(b,h) b(2).*((3*h./(2.*b(1)))-1/2.*(h./b(1)).^3);
    case 'Pentaspherical'
        type = 'bounded';
        func = @(b,h) b(2).*(15*h./(8.*b(1))-5/4*(h./b(1)).^3+3/8*(h./b(1)).^5);
    case 'Cubic'
        type = 'bounded';
        func = @(b,h) b(2).*(7*(h./b(1)).^2 -8.75*(h./b(1)).^3 +3.5*(h./b(1)).^5 -0.75*(h./b(1)).^7);
    case 'Exponential'
        type = 'unbounded';
        func = @(b,h) b(2).*(1-exp(-h./b(1)));
    case 'Gaussian'
        type = 'unbounded';
        func = @(b,h) b(2).*(1-exp(-(h.^2)./(b(1).^2)));
    case 'Matérn'
        type = 'unbounded';
        if isempty(nu)
            errordlg('Please Provide the Nu Parameter');
            S=[];
            return
        else
            func = @(b,h) b(2)*(1-(1/((2^(nu-1))*gamma(nu))) * (h/b(1)).^nu .* besselk(nu,h/b(1)));
        end
    otherwise
        errordlg('unknown model')
        return
end

% Remove zeros for Matérn
switch model
    case 'Matérn'
        izero = h==0;
        if any(izero)
            flagzerodistances = true;
        else
            flagzerodistances = false;
        end
    otherwise
        flagzerodistances = false;
end

% if model type is unbounded, then the parameter b(1) is r, which is
% approximately range/3.
switch type
    case 'unbounded'
        b0(1) = b0(1)/3;
end

% create weights (by Cressie)
weights = @(b,h) (numobs./varfunc(b,h).^2)./sum(numobs./varfunc(b,h).^2);

% create objective function: weighted least square
if ~isempty(numobs)
    objectfun = @(b)sum(((varfunc(b,h)-gammaexp).^2).*weights);
else
    objectfun = @(b)sum(((varfunc(b,h)-gammaexp).^2));
end

% Upper and Lower Bounds
lb = zeros(size(b0));
if nugget
    ub = [inf max(gammaexp) max(gammaexp)]; %
else
    ub = [inf max(gammaexp) 0];
end

% create options for fminsearch
options = optimset('MaxFunEvals',1000000);

% call solver
[b,~,exitflag] = fminsearchbnd(objectfun,b0,lb,ub,options);

% Results
S.model     = model; % model
S.range     = b(1);
S.sill      = b(2);
S.nugget    = b(2);
S.h         = h; % distance
S.gamma     = gammaexp; % experimental values
S.gammahat  = varfunc(b,h); % estimated values
S.residuals = gammaexp-S.gammahat; % residuals
S.exitflag  = exitflag; % exitflag

switch type
    case 'bounded'
        plot(h,gammaexp,'bs','MarkerSize',10,'MarkerFaceColor','c','LineWidth', 2);
        hold on
        fplot(@(h) funnugget(b) + func(b,h),[0 b(1)], 'Color','r','LineWidth', 4)
        fplot(@(h) funnugget(b) + b(2),[b(1) max(h)], 'Color','r','LineWidth', 4)
    case 'unbounded'
        plot(h,gammaexp,'bs','MarkerSize',10, 'MarkerFaceColor','c','LineWidth', 2);
        hold on
        fplot(@(h) funnugget(b) + func(b,h),[0 max(h)],'Color','r','LineWidth', 4)
end
axis([0 max(h) 0 max(gammaexp)])
xlabel('lag distance, h')
ylabel('Semivariance, \gamma(h)')
title('Semivariogram')
legend('Empirical Semivariogram', strcat('Fitted Model: ', model), 'Location', 'southeast')
set(gca, 'fontname', 'Arial', 'fontsize',14)
hold off

% fitting functions for  fminsearch/bnd
    function gammahat = varfunc(b,h)
        switch type
            case 'bounded'
                I = h<=b(1);
                gammahat     = zeros(size(I));
                gammahat(I)  = funnugget(b) + func(b,h(I));
                gammahat(~I) = funnugget(b) + b(2);
            case 'unbounded'
                gammahat = funnugget(b) + func(b,h);
                if flagzerodistances
                    gammahat(izero) = funnugget(b);
                end
        end
    end
end

function [x,fval,exitflag,output] = fminsearchbnd(fun,x0,LB,UB,options,varargin)
% size checks
xsize = size(x0);
x0 = x0(:);
n=length(x0);
if (nargin<3) || isempty(LB)
    LB = repmat(-inf,n,1);
else
    LB = LB(:);
end
if (nargin<4) || isempty(UB)
    UB = repmat(inf,n,1);
else
    UB = UB(:);
end
if (n~=length(LB)) || (n~=length(UB))
    error 'x0 is incompatible in size with either LB or UB.'
end
% set default options if necessary
if (nargin<5) || isempty(options)
    options = optimset('fminsearch');
end
% stuff into a struct to pass around
params.args = varargin;
params.LB = LB;
params.UB = UB;
params.fun = fun;
params.n = n;
params.xsize = xsize;
params.OutputFcn = [];
params.BoundClass = zeros(n,1);
for i=1:n
    k = isfinite(LB(i)) + 2*isfinite(UB(i));
    params.BoundClass(i) = k;
    if (k==3) && (LB(i)==UB(i))
        params.BoundClass(i) = 4;
    end
end
% transform starting values into their unconstrained
% surrogates. Check for infeasible starting guesses.
x0u = x0;
k=1;
for i = 1:n
    switch params.BoundClass(i)
        case 1
            % lower bound only
            if x0(i)<=LB(i)
                % infeasible starting value. Use bound.
                x0u(k) = 0;
            else
                x0u(k) = sqrt(x0(i) - LB(i));
            end
            % increment k
            k=k+1;
        case 2
            % upper bound only
            if x0(i)>=UB(i)
                % infeasible starting value. use bound.
                x0u(k) = 0;
            else
                x0u(k) = sqrt(UB(i) - x0(i));
            end
            % increment k
            k=k+1;
        case 3
            % lower and upper bounds
            if x0(i)<=LB(i)
                % infeasible starting value
                x0u(k) = -pi/2;
            elseif x0(i)>=UB(i)
                % infeasible starting value
                x0u(k) = pi/2;
            else
                x0u(k) = 2*(x0(i) - LB(i))/(UB(i)-LB(i)) - 1;
                % shift by 2*pi to avoid problems at zero in fminsearch
                % otherwise, the initial simplex is vanishingly small
                x0u(k) = 2*pi+asin(max(-1,min(1,x0u(k))));
            end
            % increment k
            k=k+1;
        case 0
            % unconstrained variable. x0u(i) is set.
            x0u(k) = x0(i);
            % increment k
            k=k+1;
        case 4
            % fixed variable. drop it before fminsearch sees it.
            % k is not incremented for this variable.
    end
end
% if any of the unknowns were fixed, then we need to shorten
% x0u now.
if k<=n
    x0u(k:n) = [];
end
% were all the variables fixed?
if isempty(x0u)
    % undo the variable transformations into the original space
    x = xtransform(x0u,params);
    % final reshape
    x = reshape(x,xsize);
    % stuff fval with the final value
    fval = feval(params.fun,x,params.args{:});
    % fminsearchbnd was not called
    exitflag = 0;
    output.iterations = 0;
    output.funcCount = 1;
    output.algorithm = 'fminsearch';
    output.message = 'All variables were held fixed by the applied bounds';
    % return with no call at all to fminsearch
    return
end
% Check for an outputfcn. If there is any, then substitute my
% own wrapper function.
if ~isempty(options.OutputFcn)
    params.OutputFcn = options.OutputFcn;
    options.OutputFcn = @outfun_wrapper;
end
% now we can call fminsearch, but with our own
% intra-objective function.
[xu,fval,exitflag,output] = fminsearch(@intrafun,x0u,options,params);
% undo the variable transformations into the original space
x = xtransform(xu,params);
% final reshape to make sure the result has the proper shape
x = reshape(x,xsize);
% Use a nested function as the OutputFcn wrapper
    function stop = outfun_wrapper(x,varargin)
        % we need to transform x first
        xtrans = xtransform(x,params);
        % then call the user supplied OutputFcn
        stop = params.OutputFcn(xtrans,varargin{1:(end-1)});
    end
end % mainline end
%=======================================
function fval = intrafun(x,params)
% transform variables, then call original function
% transform
xtrans = xtransform(x,params);
% and call fun
fval = feval(params.fun,reshape(xtrans,params.xsize),params.args{:});
end % sub function intrafun end
% ======================================
function xtrans = xtransform(x,params)
% converts unconstrained variables into their original domains
xtrans = zeros(params.xsize);
% k allows some variables to be fixed, thus dropped from the
% optimization.
k=1;
for i = 1:params.n
    switch params.BoundClass(i)
        case 1
            % lower bound only
            xtrans(i) = params.LB(i) + x(k).^2;
            k=k+1;
        case 2
            % upper bound only
            xtrans(i) = params.UB(i) - x(k).^2;
            k=k+1;
        case 3
            % lower and upper bounds
            xtrans(i) = (sin(x(k))+1)/2;
            xtrans(i) = xtrans(i)*(params.UB(i) - params.LB(i)) + params.LB(i);
            % just in case of any floating point problems
            xtrans(i) = max(params.LB(i),min(params.UB(i),xtrans(i)));
            k=k+1;
        case 4
            % fixed variable, bounds are equal, set it at either bound
            xtrans(i) = params.LB(i);
        case 0
            % unconstrained variable.
            xtrans(i) = x(k);
            k=k+1;
    end
end
end
