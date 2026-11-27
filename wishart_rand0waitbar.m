function simimage = wishart_rand0waitbar(L,nrows,ncols,Sigma,varargin)

% Generate observations from a complex Wishart distribution W(p,ENL,Sigma)
%
% same as wishart_rand but no waitbar
%
% simimage = wishart_rand(L,nrows,ncols,Sigma,varargin)
%
% Input
%   L     - ENL, equivalent number of looks, must be an integer  
%   nrows - number of rows, defaults to 1
%   ncols - number of cols, defaults to 1
%   Sigma - 3x3 or 2x2 Hermitian variance-covariance matrix,
%           defaults to 3x3 unit matrix
%           (must be specified to generate 2x2 matrices)
%           or
%           row vector containing upper right elements of Hermitian matrix
%   k     - number of "time points" generated (optional, defaults to 1)
%   if nargin>5 plot 1-D histograms of generated data
%
% Remember, that here we generate covariance matrices.
% To obtain complex Wishart distributed matrices we must therefore multiply
% by ENL:
%
% L = 10; nrows = 5; ncols = 20;
% Sigma = [1 2+1j*3 4+1j*5;nan 6 7+1j*8;nan nan 9]+6*eye(3)
% C = wishart_rand(L,nrows,ncols,Sigma,4,true);
% W = L*C;
% 
% or (where [5 2 3 4] represents the Hermitian matrix [5 2+1i*3;2-1i*3 4])
%
% C = wishart_rand(L,nrows,ncols,[5 2 3 4]);
% wc2mat(C)
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
% https://www2.imm.dtu.dk/pubdb/pubs/7124-full.html

% (c) Copyright 2018-2024
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 23-24 Oct 2018
% 2 Aug 2024 /AA: histograms added
% 8 Dec 2024 /AA: Sigma as row vector (upper right matrix, p=1) added
% 16 Jan 2025 /AA: k added

% s = rng;
% warning("You might want to use rng('default')");

if nargin<1
    error('wishart_rand: wishart_rand(L,nrows,ncols,Sigma,varargin)')
end
if nargin<2
    nrows = 1;
end
if nargin<3
    ncols = nrows;
end
p = 3;
if nargin>3 % Sigma specified
    % if isvector(squeeze(Sigma))
    %     Sigma = wc2mat(squeeze(Sigma));
    % end
    p = size(Sigma,1);
    if ~(p==2 || p==3 || p==1)
        error('wishart_rand: Sigma must be 2x2 or 3x3 or row vector with UR matrix')
    end
    if p==1
        if numel(Sigma)==9
            p = 3;
        elseif numel(Sigma)==4
            p = 2;
        else
            error('wishart_rand: numel of Sigma must be 9 or 4')
        end
        Sigma = wc2mat(Sigma);
    end
    C = chol(Sigma);
end
k = 1;
if nargin>4 % k (number of "time points") specified
    k = varargin{1};
end

nobs = nrows*ncols;
simimage = nan(nobs,p*p,k);
% hw = waitbar(0,'Processing...','Name','wishart_rand');
for ik = 1:k
    for i=1:nobs
        X = randn(L,p);
        Y = randn(L,p);
        W = 0.5*(X'*X + Y'*Y - 1j*(X'*Y - Y'*X));
        if nargin>3 % Sigma specified
            W = C'*W*C;
        end
        if p==3
            simimage(i,1,ik) = real(W(1,1));
            simimage(i,2,ik) = real(W(1,2));
            simimage(i,3,ik) = imag(W(1,2));
            simimage(i,4,ik) = real(W(1,3));
            simimage(i,5,ik) = imag(W(1,3));
            simimage(i,6,ik) = real(W(2,2));
            simimage(i,7,ik) = real(W(2,3));
            simimage(i,8,ik) = imag(W(2,3));
            simimage(i,9,ik) = real(W(3,3));
        else
            simimage(i,1,ik) = real(W(1,1));
            simimage(i,2,ik) = real(W(1,2));
            simimage(i,3,ik) = imag(W(1,2));
            simimage(i,4,ik) = real(W(2,2));
        end
    end
    % waitbar(ik/k,hw);
end
% close(hw)

simimage = reshape(simimage,nrows,ncols,p*p,k);
simimage = simimage/L;

% rng(s);

if nargin>5
    ax = nan(1,p*p);
    figure('Name','Generated matrix data')
    for ix=1:p*p
        ax(ix) = subplot_tight(p,p,ix,[0.06,0.06]);
        h = histogram(simimage(:,:,ix,:));
        if nrows*ncols>100*100+1 % or some other 'high' number
            h.EdgeColor = 'none';
        end
        grid on, grid minor
    end
    % linkaxes(ax,'xy')
    if p==3
        linkaxes([ax(1) ax(6) ax(9)],'xy') % real elements
        % xlim([900 1350])
        linkaxes([ax(2) ax(3) ax(4) ax(5) ax(7) ax(8)],'xy') % complex
        % xlim([900 1350])
    else
        linkaxes([ax(1) ax(4)],'xy') % real elements
        % xlim([900 1350])
        linkaxes([ax(2) ax(3)],'xy') % complex
        % xlim([900 1350])
    end
end
