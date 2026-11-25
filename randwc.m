function [W,C,X,Y] = randwc(n,nobs,Sigma,varargin)

% Generate observations from a complex Wishart distribution W(p,n,Sigma)
%
% [W,C,X,Y] = randwc(n,nobs,Sigma,varargin)
%
% Input
%   n     - number of observations (ENL)
%   nobs  - number of realisations of W, defaults to one
%   Sigma - p by p symmetric variance-covariance matrix,
%           defaults to 3x3 unit matrix,
%           p is given by dimension of Sigma
%   1     - 1 or any input: if nargin>3 and Sigma==1, plot histogram
%           generated real gamma variable
%
% Output
%   W     - complex Wishart distributed matrix/matrices
%   C     - Cholesky decomposition of Sigma
%   X     - normal distributed aux variable(s) for last observation
%   Y     - normal distributed aux variable(s) for last observation
%
% If Sigma==1, randwc generates (scalar) gamma distributed W with
% shape parameter alpha=n and scale parameter theta=1
% (actually, for any dimensionality of Sigma, all diagonal elements in W
% divided by the corresponding element in Sigma will be gamma
% distributed with shape parameter alpha=n and scale parameter theta=1)
%
% Remember, that here we generate complex Wishart distributed matrices.
% We must therefore divide by n to get covariance matrices:
%
% n = 10; nobs = 20;
% % Sigma = [1 2 4;nan 6 7;nan nan 9]+6*eye(3)
% Sigma = [1 2 4 5;nan 6 7 8 ;nan nan 9 10;nan nan nan 11]+12*eye(4)
% W = randwc(n,nobs,Sigma);
% Wcov = W/n;
%
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

% (c) Copyright 2024
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 31 Dec 2024

% s = rng;
% rng(19520427) % seed for rand (and randwc)
fprintf("\nYou might want to use rng('default')\n\n");

if nargin<1
    error('randwc: randwc(n,nobs,Sigma,1)')
end
if nargin<2
    nobs = 1;
end
p = 3;
C = nan;
if nargin>2 % Sigma specified
    p = size(Sigma,1);
    C = chol(Sigma);
end

W = nan(p,p,nobs);
for ix=1:nobs
    X = randn(n,p);
    Y = randn(n,p);
    W(:,:,ix) = 0.5*(X'*X + Y'*Y - 1j*(X'*Y - Y'*X));
    if nargin>2 % Sigma specified
        W(:,:,ix) = C'*W(:,:,ix)*C;
    end
end

if nargin>3
    if isscalar(Sigma)==1 && Sigma==1
        fprintf('gamma distribution with shape parameter alpha=%d (and scale parameter theta=1)\n',n)
        fprintf('mode is %d\n',n-1)
        figure('Name','randwc (gamma distr.)')
        if size(W(:),1)<100*1000-1
            h = histogram(W(:),'Normalization','pdf');
        else
            h = histogram(W(:),'Normalization','pdf','EdgeColor','none');
        end
        grid on, grid minor
        hold on
        xx = 0:0.1:h.BinLimits(2);
        f = gampdf(xx,n);%    f = gampdf(xx,n,1);
        plot(xx,f,'LineWidth',2)
        hold off
    else
        fprintf('*** randwc: request for histogram plot ignored ***\n')
    end
end
