function simimage = wishart_rand_scene(L,SAR,varargin)

% Generate observations from polSAR imagery W(ENL,SAR)
%
% simimage = wishart_rand_scene(L,SAR,varargin)
%
% Based on work by Knut Conradsen and Morton J. Canty
%
% Input
%   L     - ENL, equivalent number of looks, must be an integer  
%   SAR   - 3-D, nrows x ncols polSAR scene with 2x2 or 3x3 cov. mat.
%   k     - number of realizations (optional, defaults to 1)
%   if nargin>3 plot 1-D histograms of generated data
%
% Remember, that here we generate covariance matrices.
% To obtain complex Wishart distributed matrices we must therefore multiply
% by ENL:
%
% C = wishart_rand_scene(L,SAR); % polSAR covariance matrix representation
% or
% C = wishart_rand_scene(L,SAR,4); % polSAR covariance matrix representation
% W = L*C; % % complex Wishart distributed
%
% see wishart_rand
%
% References
%
% Allan A. Nielsen, Henning Skriver and Knut Conradsen (2019).
% "Generation of Sample Complex Wishart Distributed Matrices and Change
% Detection in Polarimetric SAR Data".
% International Conference on Digital Image and Signal Processing,
% DISP 2019, Article no. 031, St Hugh's College, Oxford University,
% England, 29-30 April 2019.
% https://doi.org/10.48550/arXiv.2611.22943.
% https://www2.imm.dtu.dk/pubdb/pubs/7124-full.html
%
% K. Conradsen, H. Skriver and A. A. Nielsen (2024).
% “Omnibus change detection in block diagonal covariance matrix polSAR data
% illustrated with simulated and Sentinel-1 data,”
% IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing
% pp. 16359–16376,
% https://doi.org/10.1109/JSTARS.2024.3453442.
%
% A. A. Nielsen, H.Skriver and K. Conradsen (2024).
% "Multi-Frequency Omnibus Change Detection in Covariance Matrix PolSAR
% Data,"
% in prep.

% (c) Copyright 2018-2024
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 23-24 Oct 2018
% 17 Jul 2024
% 8 Aug 2024 /AA: k added
% 26 Jan 2025 /AA: simimage divided by L

% s = rng;
fprintf("\nYou might want to use rng('default')\n");
fprintf("\n***********************************************************\n");
fprintf("*** wishart_rand_scene generates covariance matrix data ***\n")
fprintf("*** multiply by L to get Wishart distributed data       ***\n");
fprintf("***********************************************************\n");

if nargin<2
    error('wishart_rand_scene(L,SAR), SAR must be 3-D nrows x ncols with 2x2 or 3x3 cov. mat.')
end
if nargin>4
    error('wishart_rand_scene(L,SAR,k,1), SAR must be 3-D nrows x ncols with 2x2 or 3x3 cov. mat.')
end
k = 1;
if nargin>2
    k = varargin{1};
end

if ismatrix(SAR) | ndims(SAR)>4 % SAR must be 3-D or 4-D
% if ndims(SAR)<3 | ndims(SAR)>4 % SAR must be 3-D or 4-D
    error('wishart_rand_scene: SAR must be 3-D nrows x ncols with 2x2 or 3x3 cov. mat.')
end
if ndims(SAR)==4 % SAR is 4-D
    % [~,~,~,nt] = size(SAR);
    fprintf('wishart_rand_scene: SAR is 4-D, we use first ''time point'' only\n')
    SAR = squeeze(SAR(:,:,:,1));
end
[nrows,ncols,pp] = size(SAR); % SAR is 3-D
if ~(pp==4 || pp==9)
    error('wishart_rand_scene: cov. mat. must be 2x2 or 3x3')
end
p = sqrt(pp);

nobs = nrows*ncols;
SAR = reshape(SAR,nobs,pp);

simimage = nan(nobs,pp,k);
hw = waitbar(0,'Processing...','Name','wishart_rand_scene');
for ik = 1:k
    for ix=1:nobs
        X = randn(L,p);
        Y = randn(L,p);
        W = 0.5*(X'*X + Y'*Y - 1j*(X'*Y - Y'*X));
        C = chol(wc2mat(reshape(SAR(ix,:),1,1,pp)));
        W = C'*W*C;
        if p==3
            simimage(ix,1,ik) = real(W(1,1));
            simimage(ix,2,ik) = real(W(1,2));
            simimage(ix,3,ik) = imag(W(1,2));
            simimage(ix,4,ik) = real(W(1,3));
            simimage(ix,5,ik) = imag(W(1,3));
            simimage(ix,6,ik) = real(W(2,2));
            simimage(ix,7,ik) = real(W(2,3));
            simimage(ix,8,ik) = imag(W(2,3));
            simimage(ix,9,ik) = real(W(3,3));
        else
            simimage(ix,1,ik) = real(W(1,1));
            simimage(ix,2,ik) = real(W(1,2));
            simimage(ix,3,ik) = imag(W(1,2));
            simimage(ix,4,ik) = real(W(2,2));
        end
    end
    waitbar(ik/k,hw);
end
close(hw)

simimage = reshape(simimage,nrows,ncols,pp,k);
simimage = simimage/L;
% if k==1
%     simimage = squeeze(simimage);
% end

% rng(s);

if nargin>3
    if k>1
        fprintf('wishart_rand_scene: ignore request for histograms for k>1\n')
        return
    end
    ax = nan(1,pp);
    figure('Name','Generated matrix data')
    for ix=1:pp
        ax(ix) = subplot_tight(p,p,ix,[0.06,0.06]);
        h = histogram(simimage(:,:,ix));
        % if nrows*ncols>100*100+1 % or some other 'high' number
        %     h.EdgeColor = 'none';
        % end
        grid on
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
