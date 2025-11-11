function geoPlot(X,y,bw,kernel,grain,geographicRegion,Latitude,Longitude,responseType)
% Creates plot that overlays a smoothed spatial regression estimates over a
% geographical region (either specified by user or by imported perimeter data)
%
% Input Arguments:
% X :               Input data (Typically, 2-D Spatial dataset)
% y :               Response
% bw:               Bandwidth (scalar for 1-D X, 2-by-2 symmetric matrix for 2-D X)
% kernel:           Kernel used for local polynomial regression
% grain:            Plot resoultion
% geographicRegion: Speicifies the geographic perimeter ('Contiguous US' or 'N/A')
% Latitude:         Latitude of the geographic region (if geographicRegion is 'N/A')
% Longitude:        Longitude of the geographic region (if geographicRegion is 'N/A')
% responseType:     Sign of the response variable ('Positive', 'Negative', or 'other')
%
% Output
% A geoplot containing the smoothed regression function

% Authors: Mehnuma Tabassum, Kris De Brabanter
% Email: mehnuma@iastate.edu

[~, dimX] = size(X);

if isempty(bw)
    errordlg('Please provide the bandwidth');
end

if isempty(grain)
    grain = 50;
end

latlon = [Latitude Longitude];

[rBW, cBW] = size(bw);
if ~(rBW==2) || ~(cBW==2)
    errordlg('Please provide a 2x2 symmetric bandwidth matrix');
    return
elseif ~(bw(1,2)==bw(2,1))
    errordlg('Please provide a 2x2 symmetric bandwidth matrix');
    return
else
    if dimX==1
        errordlg('Please provide 2-D data (latitude and longitude)');
        return
    elseif dimX==2
        if strcmp(geographicRegion, 'Contiguous US')
            geodata = readmatrix('usaPerimeter.csv');
            G.Latitude = geodata(:,1); G.Longitude = geodata(:,2);
        elseif strcmp(geographicRegion, 'N/A') && ~isempty(latlon)
            G.Latitude = latlon(:,1); G.Longitude = latlon(:,2);
        else
            errordlg('Please select from the available options or import the perimeter data');
            return
        end

        ltlim = [min(G.Latitude) max(G.Latitude)]+ [-1 1];
        lnlim = [min(G.Longitude) max(G.Longitude)]+ [-1 1];

        xmin1 = min(X(:,1)); if xmin1<0, xmin1=1.05*xmin1; else xmin1 = 0.98*xmin1; end
        xmax1 = max(X(:,1)); if xmax1>0, xmax1=1.05*xmax1; else xmax1 = 0.98*xmax1; end
        xmin2 = min(X(:,2)); if xmin2<0, xmin2=1.05*xmin2; else xmin2 = 0.98*xmin2; end
        xmax2 = max(X(:,2)); if xmax2>0, xmax2=1.05*xmax2; else xmax2 = 0.98*xmax2; end
        xrange1 = xmin1:(xmax1-xmin1)/grain:xmax1;
        xrange2 = xmin2:(xmax2-xmin2)/grain:xmax2;
        [XX,YY] = meshgrid(xrange1,xrange2);

        Xt = [reshape(XX,numel(XX),1) reshape(YY,numel(YY),1)];
        w = waitbar(0.5,'Fetching your plot ...', 'Name', 'Geoplotting');
        pause(0.5)
        [yhat,~] = locpol(X,y,bw,kernel,Xt);
        if strcmp(responseType, 'Positive')
            yhat = max(yhat, 0);
            yHat = reshape(yhat,size(XX,1),size(XX,2));
        elseif strcmp(responseType, 'Negative')
            yhat = min(yhat, 0);
            yHat = reshape(yhat,size(XX,1),size(XX,2));
        else
            yHat = reshape(yhat,size(XX,1),size(XX,2));
        end

        waitbar(0.98, w, 'Almost done ...');
        pause(1)
        close(w)

        figure('color','w');
        hold on
        worldmap(ltlim, lnlim)
        geoshow(XX, YY, yHat, 'DisplayType','surface', 'zdata', ones(size(YY))*-1, 'cdata', yHat, 'facecolor', 'flat');          % to show coarseness

        ltbox = ltlim([1 2 2 1 1]);
        lnbox = lnlim([1 1 2 2 1]);
        [ltbox, lnbox] = interpm(ltbox, lnbox, 1);
        [xbox, ybox] = mfwdtran(ltbox, lnbox); % Box around data
        [xus, yus] = mfwdtran(G.Latitude, G.Longitude); % USA polygon
        A = polyshape(xbox, ybox,'Simplify',false);
        B = polyshape(xus, yus,'Simplify',false);
        C = subtract(A,B);
        [xmask, ymask] = boundary(C);
        [f,v] = poly2fv(xmask, ymask);
        patch('faces', f, 'vertices', v, 'facecolor', 'w', 'edgecolor', 'none');
        % set(hmask,'EdgeColor','none','linestyle','none');
        colorbar;
        colormap(jet(256));
        set(gca,'CLim',[min(y) max(y)])
        % xlabel('Latitude')
        % ylabel('Longitude')
        hold off
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

function [x, y, z, savepts] = mfwdtran(varargin)
[mstruct, lat, lon, alt, objectType] = parseInputs(varargin{:});
[x, y, z, savepts] = feval( ...
    mstruct.mapprojection, mstruct, lat, lon, alt, objectType, 'forward');
end

function [mstruct, lat, lon, alt, objectType] = parseInputs(varargin)
if (nargin >= 1) && isstruct(varargin{1})
    narginchk(3, 5)
    mstruct = varargin{1};
    varargin(1) = [];
else
    narginchk(2, 4)
    mstruct = gcm;
end
spheroid = mstruct.geoid;
checkellipsoid(spheroid, 'mfwdtran', 'mstruct.geoid');
if isa(spheroid, 'map.geodesy.Spheroid')
    a = spheroid.SemimajorAxis;
else
    a = spheroid(1);
end
if a <= 0
    error(message('map:validate:expectedPositiveSemimajor'))
end
lat = varargin{1};
lon = varargin{2};
% Assign default values as needed.
defaults = { ...
    zeros(size(lat)), ...  % alt
    'none'};               % objectType
varargin(end+1:numel(defaults)+2) = defaults(numel(varargin)-1:end);
alt        = varargin{3};
objectType = varargin{4};
% Ensure non-empty alt even in the case where varargin{3} is []
if isempty(alt)
    alt = zeros(size(lat));
end
end

