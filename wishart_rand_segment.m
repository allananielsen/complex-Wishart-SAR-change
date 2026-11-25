function [SARg,idx] = wishart_rand_segment(SAR,L,segments)

% generates new covariance matrices based on mean values in segments in polSAR data
%
% Input
%   SAR      - original SAR scene
%   L        - ENL, equivalent number of looks
%   segments - image with segmentation of SAR, background is zero (0)
%
% Output
%   SARg     - generated SAR scene
%   idx      - values for segments
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

% (c) Copyright 2025
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 25-26 Jan 2025

fprintf("\nYou might want to use rng('default')\n");
fprintf("\n********************************************************\n");
fprintf("*** wishart_segment generates covariance matrix data ***\n")
fprintf("*** multiply by L to get Wishart distributed data    ***\n");
fprintf("********************************************************\n\n");

if nargin<3
    error('wishart_rand_segment(SAR,L,segments);')
end

if not(ndims(SAR)==4 || ndims(SAR)==3)
    % error('for now SAR must be 4-D')
    error('SAR must be 4-D or 3-D')
end
if ndims(SAR)==3
    [nr,nc,nv] = size(SAR);
    nt = 1;
    SAR = reshape(SAR,nr,nc,nv,nt);
end

[nr,nc,nv,nt] = size(SAR);
if not(nv==4 ||nv==9)
    error('SAR must have 2x2 or 3x3 covmats')
end
nobs = nr*nc;
[nrs,ncs] = size(segments);
if nrs~=nr || ncs~=nc
    error('SAR and segments must have same #rows and #cols');
end

idx = unique(segments); % contains segment-ids incl possible background 0
% if any(idx)
%     fprintf('segments file contains 0, background\n')
% end

SARv = reshape(SAR,nobs,nv,nt);
SARg = SARv;

for ix=1:numel(idx)
    if idx(ix)==0
        fprintf('segments file contains 0, background, which we skip\n')
        continue
    end
    roi = segments==idx(ix);
    nmobs = sum(roi(:));
    roiv = reshape(roi,nobs,1);
    aux = SARv(roiv,:,:);
    maux = mean(aux,1);
    SARmg = nan(nmobs,1,nv,nt);
    for it=1:nt
        SARmg(:,:,:,it) = wishart_rand0waitbar(L,nmobs,1,wc2mat(maux(:,:,it)));
        SARg(roiv,:,it) = SARmg(:,:,:,it);
    end
end
SARg = reshape(SARg,nr,nc,nv,nt);
