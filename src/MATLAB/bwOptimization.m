function H = bwOptimization(X,y,m,kernel,objective,h_lb,h_ub,R)
% Performs bandwidth optimization for 1-D and 2-D spatially correlated
% data.
%
% Input Arguments:
% X :        Input data (1-D or 2-D)
% y :        Response
% m :        True function value (optional; required for ASE)
% kernel:    Kernel used for local polynomial regression
% objective: Objective function to be minimized
% h_lb:      Lower bound for bandwidth (scalar for 1-D X, 2-by-2 symmetric matrix for 2-D X)
% h_ub:      Upper bound for bandwidth (scalar for 1-D X, 2-by-2 symmetric matrix for 2-D X)
% R:         True correlation matrix (Only needed for Corrected CGV, GCVc)
%
% Output Argument:
% H :        Optimal bandwidth (scalar for 1-D X, 2-by-2 symmetric matrix for 2-D X)

% Authors: Mehnuma Tabassum, Kris De Brabanter
% Email: mehnuma@iastate.edu

[n,dimX] = size(X);

if isempty(kernel)
    errordlg('Please Provide a Kernel')
    return
elseif isempty(objective)
    errordlg('Please Provide the Objective')
    return
elseif isempty(h_lb)
    errordlg('Please Provide the Lower Bound')
    return
elseif isempty(h_ub)
    errordlg('Please Provide the Upper Bound')
    return
end

switch objective
    case 'Unadjusted GCV'
        if dimX==1
            if ~isscalar(h_lb) || ~isscalar(h_ub)
                errordlg('Please Provide a Scalar Bandwidth for 1-D Input');
                H = [];
                return
            else
                w = waitbar(0,'Bandwidth optimization starting ...', 'Name', 'Bandwidth Optimization Status');
                pause(2)
                waitbar(0.1, w, 'Bandwidth optimization running ...');
                H_csa = csa(rand(1,3),@(H_in) gcv(X, y, H_in, kernel));
                waitbar(0.65,w,'Still optimizing ...');
                H = fminsearchbnd(@(H_in) gcv(X, y, H_in, kernel), H_csa, h_lb, h_ub);
                waitbar(0.98,w,'Almost done ...');
                pause(2)
                close(w)
            end

        elseif dimX==2
            [row_lb, col_lb] = size(h_lb); [row_ub, col_ub] = size(h_ub);
            if ~(row_lb==2) || ~(col_lb==2) ||~(row_ub==2) || ~(col_ub==2)
                errordlg('Please Provide a 2x2 Bandwidth Matrix for the Limits');
                H = [];
                return
            elseif ~(h_lb(1,2)==h_lb(2,1)) || ~(h_ub(1,2)==h_ub(2,1))
                errordlg('Please Provide a 2x2 Bandwidth Matrix for the Limits');
                H = [];
                return
            else
                w = waitbar(0,'Bandwidth optimization starting ...', 'Name', 'Bandwidth Optimization Status');
                pause(2)
                waitbar(0.1, w, 'Bandwidth optimization running ...');
                H_csa2 = csa(rand(2,2),@(H_in) gcv2(X, y, H_in, kernel));
                waitbar(0.5,w,'Halfway there ...');
                lb2 = [h_lb(1,1); h_lb(2,2)]; ub2 = [h_ub(1,1); h_ub(2,2)];
                x_gcv2 = fminsearchbnd(@(H_in) gcv2(X, y, H_in, kernel), H_csa2, lb2, ub2);                
                H_csa3 = csa(rand(1,3), @(x) gcv3(X, y, x_gcv2(1), x_gcv2(2), x, kernel));
                waitbar(0.8,w,'Finishing soon ...');
                lb3 = h_lb(1,2); ub3 = h_ub(1,2);
                x_gcv3 = fminsearchbnd(@(x) gcv3(X, y, x_gcv2(1), x_gcv2(2), x, kernel), H_csa3, lb3, ub3);
                H = [x_gcv2(1) x_gcv3; x_gcv3 x_gcv2(2)];
                waitbar(0.98,w,'Almost done ...');
                pause(2)
                close(w)
            end
        end

    case 'Corrected GCV (GCVc)'
        if dimX==1
            if ~isscalar(h_lb) || ~isscalar(h_ub)
                errordlg('Please Provide a Scalar Bandwidth for 1-D Input');
                H = [];
                return
            else
                w = waitbar(0,'Bandwidth optimization starting ...', 'Name', 'Bandwidth Optimization Status');
                pause(2)
                waitbar(0.1, w, 'Bandwidth optimization running ...');                
                H_csa = csa(rand(1,3),@(H_in) gcv_c(X, y, H_in, kernel,R));
                waitbar(0.65,w,'Still optimizing ...');
                H = fminsearchbnd(@(H_in) gcv_c(X, y, H_in, kernel,R), H_csa, h_lb, h_ub);
                waitbar(0.98,w,'Almost done ...');
                pause(2)
                close(w)
            end
        elseif dimX==2
            [row_lb, col_lb] = size(h_lb); [row_ub, col_ub] = size(h_ub);
            if ~(row_lb==2) || ~(col_lb==2) ||~(row_ub==2) || ~(col_ub==2)
                errordlg('Please Provide a 2x2 Bandwidth Matrix for the Limits');
                H = [];
                return
            elseif ~(h_lb(1,2)==h_lb(2,1)) || ~(h_ub(1,2)==h_ub(2,1))
                errordlg('Please Provide a 2x2 Bandwidth Matrix for the Limits');
                H = [];
                return
            else
                w = waitbar(0,'Bandwidth optimization starting ...', 'Name', 'Bandwidth Optimization Status');
                pause(2)
                waitbar(0.1, w, 'Bandwidth optimization running ...');                
                H_csa2 = csa(rand(2,2),@(H_in) gcv_c2(X, y, H_in, kernel,R));
                waitbar(0.5,w,'Halfway there ...');
                lb2 = [h_lb(1,1); h_lb(2,2)]; ub2 = [h_ub(1,1); h_ub(2,2)];
                x_gcv_c2 = fminsearchbnd(@(H_in) gcv_c2(X, y, H_in, kernel,R), H_csa2, lb2, ub2);                
                H_csa3 = csa(rand(1,3), @(x) gcv_c3(X, y, x_gcv_c2(1), x_gcv_c2(2), x, kernel,R));
                waitbar(0.8,w,'Finishing soon ...');
                lb3 = h_lb(1,2); ub3 = h_ub(1,2);
                x_gcv_c3 = fminsearchbnd(@(x) gcv_c3(X, y, x_gcv_c2(1), x_gcv_c2(2), x,kernel,R), H_csa3, lb3, ub3);
                H = [x_gcv_c2(1) x_gcv_c3; x_gcv_c3 x_gcv_c2(2)];
                waitbar(0.98,w,'Almost done ...');
                pause(2)
                close(w)
            end
        end

    case 'Corrected and Estimated GCV (GCVce)'
        if dimX==1
            if ~isscalar(h_lb) || ~isscalar(h_ub)
                errordlg('Please Provide a Scalar Bandwidth for 1-D Input');
                H = [];
                return
            else
                [yhat, ~]= locpol(X,y,std(X),kernel);
            end
        elseif dimX==2
            [row_lb, col_lb] = size(h_lb); [row_ub, col_ub] = size(h_ub);
            if ~(row_lb==2) || ~(col_lb==2) ||~(row_ub==2) || ~(col_ub==2)
                errordlg('Please Provide a 2x2 Bandwidth Matrix for the Limits');
                H = [];
                return
            elseif ~(h_lb(1,2)==h_lb(2,1)) || ~(h_ub(1,2)==h_ub(2,1))
                errordlg('Please Provide a 2x2 Bandwidth Matrix for the Limits');
                H = [];
                return
            else
                [yhat, ~]= locpol(X,y,[std(X(:,1)) 0; 0 std(X(:,2))],kernel);
            end
        end
        res = y-yhat; D = (pdist2(X,X));
        dk = (0.001+2*((1:30)-1)*0.005)'; [ndk, ~] = size(dk);
        gamma = zeros(ndk,1);
        for j = 1:numel(dk)
            ind = D>=max(dk(j)-0.005, 0) & D<=dk(j)+0.005;
            [row, column] = find(triu(ind)==1);
            num_el = numel(row);
            rse = zeros(num_el, 1);
            for i = 1:num_el
                rse(i) = (res(column(i)) -res(row(i))).^2;
            end
            gamma = (1/(2*num_el))*sum(rse);
        end
        sigmaSq = mean(res.^2);
        alphaJ = (1./(n*dk)).*(log(sigmaSq) -log(abs(sigmaSq-gamma)));
        alpha = mean(alphaJ);

        if dimX==1
            w = waitbar(0,'Bandwidth optimization starting ...', 'Name', 'Bandwidth Optimization Status');
            pause(2)
            waitbar(0.1, w, 'Bandwidth optimization running ...');            
            H_csa = csa(rand(1,3),@(H_in) gcv_ce(X, y, H_in, kernel,alpha));
            waitbar(0.65,w,'Still optimizing ...');
            H = fminsearchbnd(@(H_in) gcv_ce(X, y, H_in, kernel,alpha), H_csa, h_lb, h_ub);
            waitbar(0.98,w,'Almost done ...');
            pause(2)
            close(w)
        elseif dimX==2
            w = waitbar(0,'Bandwidth optimization starting ...', 'Name', 'Bandwidth Optimization Status');
            pause(2)
            waitbar(0.1, w, 'Bandwidth optimization running ...');            
            H_csa2 = csa(rand(2,2),@(H_in) gcv_ce2(X, y, H_in, kernel, alpha));
            waitbar(0.5,w,'Halfway there ...');
            lb2 = [h_lb(1,1); h_lb(2,2)]; ub2 = [h_ub(1,1); h_ub(2,2)];
            x_gcv_ce2 = fminsearchbnd(@(H_in) gcv_ce2(X, y, H_in, kernel, alpha), H_csa2, lb2, ub2);            
            H_csa3 = csa(rand(1,3), @(x) gcv_ce3(X, y, x_gcv_ce2(1), x_gcv_ce2(2), x, kernel, alpha));
            waitbar(0.8,w,'Finishing soon ...');
            lb3 = h_lb(1,2); ub3 = h_ub(1,2);
            x_gcv_ce3 = fminsearchbnd(@(x) gcv_ce3(X, y, x_gcv_ce2(1), x_gcv_ce2(2), x, kernel, alpha), H_csa3, lb3, ub3);
            H = [x_gcv_ce2(1) x_gcv_ce3; x_gcv_ce3 x_gcv_ce2(2)];
            waitbar(0.98,w,'Almost done ...');
            pause(2)
            close(w)
        end

    case 'Residual Sum of Squares (RSS)'
        if dimX==1
            if ~isscalar(h_lb) || ~isscalar(h_ub)
                errordlg('Please Provide a Scalar Bandwidth for 1-D Input');
                H = [];
                return
            else
                if strcmp(kernel, 'Epanechnikov')
                    w = waitbar(0,'Bandwidth optimization starting ...', 'Name', 'Bandwidth Optimization Status');
                    pause(2)
                    waitbar(0.1, w, 'Bandwidth optimization running ...');                    
                    H_csa = csa(rand(1,3),@(H_in) rss(X, y, H_in, 'Kernel with K(0)=0'));
                    waitbar(0.65,w,'Still optimizing ...');
                    factor = 2.57132;
                    H = fminsearchbnd(@(H_in) rss(X, y, H_in, 'Kernel with K(0)=0'), H_csa, h_lb/factor, h_ub/factor);
                    H=H*factor;
                    waitbar(0.98,w,'Almost done ...');
                    pause(2)
                    close(w)
                elseif strcmp(kernel, 'Gaussian')
                    w = waitbar(0,'Bandwidth optimization starting ...', 'Name', 'Bandwidth Optimization Status');
                    pause(2)
                    waitbar(0.1, w, 'Bandwidth optimization running ...');                   
                    H_csa = csa(rand(1,3),@(H_in) rss(X, y, H_in, 'Kernel with K(0)=0'));
                    waitbar(0.65,w,'Still optimizing ...');
                    factor = 1.16;
                    H = fminsearchbnd(@(H_in) rss(X, y, H_in, 'Kernel with K(0)=0'), H_csa, h_lb/factor, h_ub/factor);
                    H=H*factor;
                    waitbar(0.98,w,'Almost done ...');
                    pause(2)
                    close(w)
                end
            end
        elseif dimX==2
            [row_lb, col_lb] = size(h_lb); [row_ub, col_ub] = size(h_ub);
            if ~(row_lb==2) || ~(col_lb==2) ||~(row_ub==2) || ~(col_ub==2)
                errordlg('Please Provide a 2x2 Bandwidth Matrix for the Limits');
                H = [];
                return
            elseif ~(h_lb(1,2)==h_lb(2,1)) || ~(h_ub(1,2)==h_ub(2,1))
                errordlg('Please Provide a 2x2 Bandwidth Matrix for the Limits');
                H = [];
                return
            else
                if strcmp(kernel, 'Epanechnikov')
                    w = waitbar(0,'Bandwidth optimization starting ...', 'Name', 'Bandwidth Optimization Status');
                    pause(2)
                    waitbar(0.1, w, 'Bandwidth optimization running ...');                    
                    H_csa2 = csa(rand(2,2),@(H_in) rss2(X, y, H_in, 'Kernel with K(0)=0'));
                    waitbar(0.5,w,'Halfway there ...');
                    factor = 2.58;
                    lb2 = [h_lb(1,1); h_lb(2,2)]./factor; ub2 = [h_ub(1,1); h_ub(2,2)]./factor;
                    x_rss2 = fminsearchbnd(@(H_in) rss2(X, y, H_in, 'Kernel with K(0)=0'), H_csa2, lb2, ub2);                    
                    H_csa3 = csa(rand(1,3), @(x) rss3(X, y, x_rss2(1), x_rss2(2), x, 'Kernel with K(0)=0'));
                    waitbar(0.8,w,'Finishing soon ...');
                    lb3 = h_lb(1,2)/factor; ub3 = h_ub(1,2)/factor;
                    x_rss3 = fminsearchbnd(@(x) rss3(X, y, x_rss2(1), x_rss2(2), x, 'Kernel with K(0)=0'), H_csa3, lb3, ub3);
                    H = [x_rss2(1) x_rss3; x_rss3 x_rss2(2)]*factor;
                    waitbar(0.98,w,'Almost done ...');
                    pause(2)
                    close(w)
                elseif strcmp(kernel, 'Gaussian')
                    w = waitbar(0,'Bandwidth optimization starting ...', 'Name', 'Bandwidth Optimization Status');
                    pause(2)
                    waitbar(0.1, w, 'Bandwidth optimization running ...');                    
                    H_csa2 = csa(rand(2,2),@(H_in) rss2(X, y, H_in, 'Kernel with K(0)=0'));
                    waitbar(0.5,w,'Halfway there ...');
                    factor = 1.14;
                    lb2 = [h_lb(1,1); h_lb(2,2)]./factor; ub2 = [h_ub(1,1); h_ub(2,2)]./factor;
                    x_rss2 = fminsearchbnd(@(H_in) rss2(X, y, H_in, 'Kernel with K(0)=0'), H_csa2, lb2, ub2);                    
                    H_csa3 = csa(rand(1,3), @(x) rss3(X, y, x_rss2(1), x_rss2(2), x, 'Kernel with K(0)=0'));
                    waitbar(0.8,w,'Finishing soon ...');
                    lb3 = h_lb(1,2)/factor; ub3 = h_ub(1,2)/factor;
                    x_rss3 = fminsearchbnd(@(x) rss3(X, y, x_rss2(1), x_rss2(2), x, 'Kernel with K(0)=0'), H_csa3, lb3, ub3);
                    H = [x_rss2(1) x_rss3; x_rss3 x_rss2(2)]*factor;
                    waitbar(0.98,w,'Almost done ...');
                    pause(2)
                    close(w)
                end
            end
        end

    case 'Asymptotic Squared Error (ASE)'
        if isempty(m)
            errordlg('The True Function must be Supplied');
            H = [];
            return
        else
            if dimX==1
                if ~isscalar(h_lb) || ~isscalar(h_ub)
                    errordlg('Please Provide a Scalar Bandwidth for 1-D Input');
                    H = [];
                    return
                else
                    w = waitbar(0,'Bandwidth optimization starting ...', 'Name', 'Bandwidth Optimization Status');
                    pause(2)
                    waitbar(0.1, w, 'Bandwidth optimization running ...');                    
                    H_csa = csa(rand(1,3),@(H_in) ase(X,y,m,H_in, kernel));
                    waitbar(0.65,w,'Still optimizing ...');
                    H = fminsearchbnd(@(H_in) ase(X,y,m, H_in, kernel), H_csa, h_lb, h_ub);
                    waitbar(0.98,w,'Almost done ...');
                    pause(2)
                    close(w)
                end
            elseif dimX==2
                [row_lb, col_lb] = size(h_lb); [row_ub, col_ub] = size(h_ub);
                if ~(row_lb==2) || ~(col_lb==2) ||~(row_ub==2) || ~(col_ub==2)
                    errordlg('Please Provide a 2x2 Bandwidth Matrix for the Limits');
                H = [];
                return
                elseif ~(h_lb(1,2)==h_lb(2,1)) || ~(h_ub(1,2)==h_ub(2,1))
                    errordlg('Please Provide a 2x2 Bandwidth Matrix for the Limits');
                H = [];
                return
                else
                    w = waitbar(0,'Bandwidth optimization starting ...', 'Name', 'Bandwidth Optimization Status');
                    pause(2)
                    waitbar(0.1, w, 'Bandwidth optimization running ...');                    
                    H_csa2 = csa(rand(2,2),@(H_in) ase2(X, y, m, H_in, kernel));
                    waitbar(0.5,w,'Halfway there ...');
                    lb2 = [h_lb(1,1); h_lb(2,2)]; ub2 = [h_ub(1,1); h_ub(2,2)];
                    x_ase2 = fminsearchbnd(@(H_in) ase2(X, y, m, H_in, kernel), H_csa2, lb2, ub2);                    
                    H_csa3 = csa(rand(1,3), @(x) ase3(X, y, m, x_ase2(1), x_ase2(2), x, kernel));
                    waitbar(0.8,w,'Finishing soon ...');
                    lb3 = h_lb(1,2); ub3 = h_ub(1,2);
                    x_ase3 = fminsearchbnd(@(x) ase3(X, y, m, x_ase2(1), x_ase2(2), x, kernel), H_csa3, lb3, ub3);
                    H = [x_ase2(1) x_ase3; x_ase3 x_ase2(2)];
                    waitbar(0.98,w,'Almost done ...');
                    pause(2)
                    close(w)
                end
            end
        end
end
end

%% Definitions of the Objective Functions
% Unadjusted GCV
function cvhat = gcv(X,y,H,kernel)
[N, ~] = size(X); [~, dh] = size(H); cvhat = zeros(dh,1);
for i=1:dh
    [yhat, S] = locpol(X, y, H(i), kernel);
    num = y-yhat;
    den = 1 -(1/N)*trace(S);
    cvhat(i) = mean((num/den).^2, 'all');
end
end
function cvhat = gcv2(X,y,H,kernel)
[N, ~] = size(X); [~, dh] = size(H); Hk = cell(dh^2,1); p = 1; cvhat = zeros(dh,1);
for i=1:dh
    for j=1:dh
        Hk{p} = diag([H(1,i) H(2,j)]);
        p = p+1;
    end
end
for i=1:dh
    Hin = Hk{i};
    [pred, S] = locpol(X, y, Hin, kernel);
    num = y-pred;
    den = 1 -(1/N)*trace(S);
    cvhat(i)= mean((num/den).^2, 'all');
end
end
function cvhat = gcv3(X,y,h1, h2,x,kernel)
[N, ~] = size(X); dh = length(x); cvhat = zeros(dh,1); Hk = cell(dh,1); p = 1;
for i=1:dh
    Hk{p} = [h1 x(i); x(i) h2];
    p = p+1;
end
for i=1:dh
    Hin = Hk{i};
    [pred, S] = locpol(X, y, Hin, kernel);
    num = y-pred;
    den = 1 -(1/N)*trace(S);
    cvhat(i)= mean((num/den).^2, 'all');
end
end

% Corrected GCV (GCVc)
function cvhat = gcv_c(X,y,H,kernel,R)
[N, ~] = size(X); [~, dh] = size(H); cvhat = zeros(dh,1);
for i=1:dh
    [yhat, S] = locpol(X, y, H(i), kernel);
    num = y-yhat;
    den = 1 -(1/N)*trace(S*R);
    cvhat(i) = mean((num/den).^2, 'all');
end
end
function cvhat = gcv_c2(X,y,H,kernel,R)
[N, ~] = size(X); [~, dh] = size(H); Hk = cell(dh^2,1); p = 1; cvhat = zeros(dh,1);
for i=1:dh
    for j=1:dh
        Hk{p} = diag([H(1,i) H(2,j)]);
        p = p+1;
    end
end
for i=1:dh
    Hin = Hk{i};
    [pred, S] = locpol(X, y, Hin, kernel);
    num = y-pred;
    den = 1 -(1/N)*trace(S*R);
    cvhat(i)= mean((num/den).^2, 'all');
end
end
function cvhat = gcv_c3(X,y,h1,h2,x,kernel,R)
[N, ~] = size(X); dh = length(x); cvhat = zeros(dh,1); Hk = cell(dh,1); p = 1;
for i=1:dh
    Hk{p} = [h1 x(i); x(i) h2];
    p = p+1;
end
for i=1:dh
    Hin = Hk{i};
    [pred, S] = locpol(X, y, Hin, kernel);
    num = y-pred;
    den = 1 -(1/N)*trace(S*R);
    cvhat(i)= mean((num/den).^2, 'all');
end
end

% Corrected and Estimated GCV (GCVce)
function cvhat = gcv_ce(X,y,H,kernel,alpha)
[N, ~] = size(X); [~, dh] = size(H); cvhat = zeros(dh,1);
distance = pdist2(X,X); R_exp = exp(-(alpha*N)*distance);
for i=1:dh
    [yhat, S] = locpol(X, y, H(i), kernel);
    num = y-yhat;
    den = 1 -(1/N)*trace(S*R_exp);
    cvhat(i) = mean((num/den).^2, 'all');
end
end
function cvhat = gcv_ce2(X,y,H,kernel,alpha)
[N, ~] = size(X); [~, dh] = size(H); Hk = cell(dh^2,1); p = 1; cvhat = zeros(dh,1);
for i=1:dh
    for j=1:dh
        Hk{p} = diag([H(1,i) H(2,j)]);
        p = p+1;
    end
end
distance = pdist2(X,X); R_exp = exp(-(alpha*N)*distance);
for i=1:dh
    Hin = Hk{i};
    [pred, S] = locpol(X, y, Hin, kernel);
    num = y-pred;
    den = 1 -(1/N)*trace(S*R_exp);
    cvhat(i)= mean((num/den).^2, 'all');
end
end
function cvhat = gcv_ce3(X,y,h1,h2,x,kernel,alpha)
[N, ~] = size(X); dh = length(x); Hk = cell(dh,1); p = 1; cvhat = zeros(dh,1);
for i=1:dh
    Hk{p} = [h1 x(i); x(i) h2];
    p = p+1;
end
distance = pdist2(X,X); R_exp = exp(-(alpha*N)*distance);
for i=1:dh
    Hin = Hk{i};
    [pred, S] = locpol(X, y, Hin, kernel);
    num = y-pred;
    den = 1 -(1/N)*trace(S*R_exp);
    cvhat(i)= mean((num/den).^2, 'all');
end
end

% Residual Sum of Squares (RSS) with KDB Kernel
function cvhat = rss(X,y,H,kernel)
[~, dh] = size(H); cvhat = zeros(dh,1);
for i=1:dh
    [yhat, ~] = locpol(X, y, H(i), kernel);
    cvhat(i) = sum((y-yhat).^2);
end
end
function cvhat = rss2(X,y,H,kernel)
[~, dh] = size(H); Hk = cell(dh^2,1); p = 1; cvhat = zeros(dh,1);
for i=1:dh
    for j=1:dh
        Hk{p} = diag([H(1,i) H(2,j)]);
        p = p+1;
    end
end
for i=1:dh
    Hin = Hk{i};
    [yhat, ~] = locpol(X, y, Hin, kernel);
    cvhat(i)= sum((y-yhat).^2);
end
end
function cvhat = rss3(X,y,h1,h2,x,kernel)
dh = length(x); Hk = cell(dh,1); p = 1; cvhat = zeros(dh,1);
for i=1:dh
    Hk{p} = [h1 x(i); x(i) h2];
    p = p+1;
end
for i=1:dh
    Hin = Hk{i};
    [yhat, ~] = locpol(X, y, Hin, kernel);
    cvhat(i)= sum((y-yhat).^2);
end
end

% Asymptotic Sqaured Error (ASE)
function cvhat = ase(X,y,m,H,kernel)
[~, dh] = size(H); cvhat = zeros(dh,1);
for i=1:dh
    [yhat, ~] = locpol(X, y, H(i), kernel);
    cvhat(i) = sum((m-yhat).^2);
end
end
function cvhat = ase2(X,y,m,H,kernel)
[~, dh] = size(H); Hk = cell(dh^2,1); p = 1; cvhat = zeros(dh,1);
for i=1:dh
    for j=1:dh
        Hk{p} = diag([H(1,i) H(2,j)]);
        p = p+1;
    end
end
for i=1:dh
    Hin = Hk{i};
    [mhat, ~] = locpol(X, y, Hin, kernel);
    cvhat(i)= sum((m-mhat).^2);
end
end
function cvhat = ase3(X, y, m, h1, h2, x, kernel)
dh = length(x); Hk = cell(dh,1); p = 1; cvhat = zeros(dh,1);
for i=1:dh
    Hk{p} = [h1 x(i); x(i) h2];
    p = p+1;
end
for i=1:dh
    Hin = Hk{i};
    [mhat, ~] = locpol(X, y, Hin, kernel);
    cvhat(i)= sum((m-mhat).^2);
end
end

%% Supporting Functions
function [pfinal,efinal] = csa(pn,herrfunc,varargin)
%switch length(varargin)
OPT.T0 = 1; OPT.Tac0 = 1; OPT.FEmax = 40; OPT.FTsteps = 50; OPT.etol = 1e-45; OPT.print = 0;
T0 = OPT.T0;  % initial temperature
Tac0 = OPT.Tac0;  % initial temperature
FEmax = OPT.FEmax; % max # of function evaluations
FTsteps = OPT.FTsteps; % # of steps at fix temperature
etol = OPT.etol; % energy tolerance

clear OPT
% initializes M
pdim = size(pn,1);
pnum = size(pn,2);
NT = ceil(FEmax/FTsteps/pnum); % # max. number of cooling cycles
NI = FTsteps; % #steps per temperature
rand('twister',sum(pn(1)*clock))
randn('state',sum(pn(2)*clock))
e0 = feval(herrfunc,pn,varargin{:});
%if any(e0<0), etol = -etol*etol^-2;end
p0 = pn;
[be0,ind] = min(e0);
bp0 = pn(:,ind);
pblty = zeros(1,pnum);
sgnCR = -1;
CR = 0.1;%0.05;
pvar_est = 0.995;
Tac = Tac0;
for k = 1:NT
    % C = progress('init',['Determine initial tuning parameters for simplex...',': # cooling cycle(s) ', num2str(k)]);
    pbltvar = var(pblty,1);
    sgnCR_ant = sgnCR;
    sgnCR = 2*((pbltvar > (pvar_est*(pnum-1)/(pnum^2)))-0.5);
    Tac = Tac + sgnCR*CR*Tac;
    % T schedules
    % T = T0/log(k+1);
    T = T0/k;
    for l = 1:NI
        % choose new coordinates and compute
        % the function to minimize
        r = tan(pi*(rand(pdim,pnum)-0.5));
        % ****************** Wrapping ***************
        %pn = 2*mod((p0 + r * T + 1)/2,1)-1;
        %**************** Non wrapping *************
        pn = abs(p0 + r * T) ;%* (diag(e0)./sum(e0));
        indd = find(abs(pn)>10);
        while(numel(indd))
            r(indd) = tan(pi*(rand(size(indd))-0.5));
            pn = abs(p0 + r * T);%* (diag(e0)./sum(e0));
            indd = find(abs(pn)>15);
        end
        en = feval(herrfunc,pn,varargin{:});
        Esum = sum(exp((e0-max(e0))./Tac));
        for i=1:pnum
            pblty(i) = min(1,exp((e0(i)-max(e0))/Tac)/Esum);
            if (en(i)-e0(i)) < 0
                % accept
                p0(:,i) = pn(:,i); e0(i) = en(i);
                if e0(i) < be0
                    be0 = e0(i);
                    bp0 = p0(:,i);
                end
            else
                r = rand;
                if (pblty(i)) >= r
                    % accept
                    p0(:,i) = pn(:,i); e0(i) = en(i);
                end
            end
        end
        % C = progress(C,l/NI);
        if any(e0<etol), break, end
    end
    if any(e0<etol), break, end
end
% kfinal = ((k-1)*NI + l)*pnum;
efinal = be0;
pfinal = bp0;
end

function [x,fval,exitflag,output] = fminsearchbnd(fun,x0,LB,UB,options,varargin)
% FMINSEARCHBND: FMINSEARCH, but with bound constraints by transformation
% usage: x=FMINSEARCHBND(fun,x0)
% usage: x=FMINSEARCHBND(fun,x0,LB)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB,options)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB,options,p1,p2,...)
% usage: [x,fval,exitflag,output]=FMINSEARCHBND(fun,x0,...)
% 
% arguments:
%  fun, x0, options - see the help for FMINSEARCH
%
%  LB - lower bound vector or array, must be the same size as x0
%
%       If no lower bounds exist for one of the variables, then
%       supply -inf for that variable.
%
%       If no lower bounds at all, then LB may be left empty.
%
%       Variables may be fixed in value by setting the corresponding
%       lower and upper bounds to exactly the same value.
%
%  UB - upper bound vector or array, must be the same size as x0
%
%       If no upper bounds exist for one of the variables, then
%       supply +inf for that variable.
%
%       If no upper bounds at all, then UB may be left empty.
%
%       Variables may be fixed in value by setting the corresponding
%       lower and upper bounds to exactly the same value.
%
% Notes:
%
%  If options is supplied, then TolX will apply to the transformed
%  variables. All other FMINSEARCH parameters should be unaffected.
%
%  Variables which are constrained by both a lower and an upper
%  bound will use a sin transformation. Those constrained by
%  only a lower or an upper bound will use a quadratic
%  transformation, and unconstrained variables will be left alone.
%
%  Variables may be fixed by setting their respective bounds equal.
%  In this case, the problem will be reduced in size for FMINSEARCH.
%
%  The bounds are inclusive inequalities, which admit the
%  boundary values themselves, but will not permit ANY function
%  evaluations outside the bounds. These constraints are strictly
%  followed.
%
%  If your problem has an EXCLUSIVE (strict) constraint which will
%  not admit evaluation at the bound itself, then you must provide
%  a slightly offset bound. An example of this is a function which
%  contains the log of one of its parameters. If you constrain the
%  variable to have a lower bound of zero, then FMINSEARCHBND may
%  try to evaluate the function exactly at zero.
%
%
% Example usage:
% rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2;
%
% fminsearch(rosen,[3 3])     % unconstrained
% ans =
%    1.0000    1.0000
%
% fminsearchbnd(rosen,[3 3],[2 2],[])     % constrained
% ans =
%    2.0000    4.0000
%
% See test_main.m for other examples of use.
%
%
% See also: fminsearch, fminspleas
%
%
% Author: John D'Errico
% E-mail: woodchips@rochester.rr.com
% Release: 4
% Release date: 7/23/06
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
end % sub function xtransform end

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
