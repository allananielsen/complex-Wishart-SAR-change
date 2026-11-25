function simimage = wc_rand(L,nrows,ncols,Sigma)

% Generate observations from a complex Wishart distribution W(p,ENL,Sigma)
%
% Input
%   L     - ENL, equivalent number of looks, must be an integer  
%   nrows - number of rows, defaults to 1
%   ncols - number of cols, defaults to 1
%   Sigma - 3x3 or 2x2 Hermitian variance-covariance matrix,
%           defaults to 3x3 unit matrix
%           (must be specified to generate 2x2 matrices)
%
% Remember, that here we generate complex Wishart distributed matrices.
% We must therefore divide by ENL to get covariance matrices:
%
% L = 10; nrows = 5; ncols = 20;
% Sigma = [1 2+1j*3 4+1j*5;nan 6 7+1j*8;nan nan 9]+6*eye(3)
% X = wc_rand(L,nrows,ncols,Sigma);
% C = X/L;

% Based on work by Knut Conradsen and Morton J. Canty
%
% Reference
%
% Allan A. Nielsen, Henning Skriver and Knut Conradsen (2019).
% "Generation of Sample Complex Wishart Distributed Matrices and Change
% Detection in Polarimetric SAR Data".
% International Conference on Digital Image and Signal Processing,
% DISP 2019, Article no. 031, St Hugh's College, Oxford University,
% England, 29-30 April 2019.
% https://doi.org/10.48550/arXiv.2611.22943.
% https://www2.imm.dtu.dk/pubdb/pubs/7124-full.html

% (c) Copyright 2018
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 23-24 Oct 2018

% s = rng;

if nargin<2
    nrows = 1;
end
if nargin<3
    ncols = nrows;
end
p = 3;
if nargin==4 % Sigma specified
    p = size(Sigma,1);
    if ~(p==1 || p==2 || p==3)
        error('Sigma must be 1x1, 2x2 or 3x3')
    end
    C = chol(Sigma);
end

nobs = nrows*ncols;
simimage = nan(nobs,p*p);
for i=1:nobs
    X = randn(L,p);
    Y = randn(L,p);
    W = 0.5*(X'*X + Y'*Y - 1j*(X'*Y - Y'*X));
    if nargin==4 % Sigma specified
        W = C'*W*C;
    end
    if p==3
        simimage(i,1) = real(W(1,1));
        simimage(i,2) = real(W(1,2));
        simimage(i,3) = imag(W(1,2));
        simimage(i,4) = real(W(1,3));
        simimage(i,5) = imag(W(1,3));
        simimage(i,6) = real(W(2,2));
        simimage(i,7) = real(W(2,3));
        simimage(i,8) = imag(W(2,3));
        simimage(i,9) = real(W(3,3));
    elseif p==2
        simimage(i,1) = real(W(1,1));
        simimage(i,2) = real(W(1,2));
        simimage(i,3) = imag(W(1,2));
        simimage(i,4) = real(W(2,2));
    else % p==2
        simimage(i,1) = real(W(1,1));
    end
end

simimage = reshape(simimage,nrows,ncols,p*p);

% rng(s);
