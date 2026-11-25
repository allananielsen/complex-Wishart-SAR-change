function ENL = ENLfilter(X,varargin)

% Estimate ENL of polSAR image with one or more polarizations
%   
% ENL = ENLfilter(X,filtersize);

% (c) Copyright 2018-2024
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, http://people.compute.dtu.dk/alan
% 22 Dec 2018 - 8 Nov 2024 /AA: 4-D added

filtersize = 11; % odd
if nargin>1
    filtersize = varargin{1};
end

if ismatrix(X)
    Xf = filter2(fspecial('average',filtersize),X);
    X2f = filter2(fspecial('average',filtersize),X.^2);
elseif isscalar(X)
    error('input should either 2-D or 3-D array')
elseif isvector(X)
    error('input should either 2-D or 3-D array')
elseif ndims(X)==3 % 3-D array
    [~,~,nv] = size(X);
    Xf = nan(size(X));
    X2f = nan(size(X));
    for iv=1:nv
        Xf(:,:,iv) = filter2(fspecial('average',filtersize),X(:,:,iv));
        X2f(:,:,iv) = filter2(fspecial('average',filtersize),X(:,:,iv).^2);
    end
elseif ndims(X)==4 % 4-D array
    [~,~,nv,k] = size(X);
    Xf = nan(size(X));
    X2f = nan(size(X));
    for iv=1:nv
        for ik=1:k
            Xf(:,:,iv,ik) = filter2(fspecial('average',filtersize),X(:,:,iv,ik));
            X2f(:,:,iv,ik) = filter2(fspecial('average',filtersize),X(:,:,iv,ik).^2);
        end
    end
else
    error('input dimension is wrong')
end

Xf2 = Xf.^2; % mean of X squared
ENL = X2f - Xf2; % mean of squared X minus mean of X squared: variance
ENL = Xf2./ENL; % inverse coefficient of variation squared
