function [lo,c,ilast,nsets] =  wc_when_loewner(SAR,wc,pval)

%  wc_when_loewner takes struct wc from wcRjl as input and outputs occurrence of change
%
% c =  wc_when_loewner(SAR,wc,0.9999);
%
% input
%   SAR is polSAR data (4-D)
%   wc is struct from wcRjl (based on SAR)
%   pval is probability threshold for change prob, default 0.99
% output
%   lo is Loewner order with mean values since last change
%   c(:,:,i) is array showing if lnQ and lnR shows change after time point i
%   ilast is time point for last change
%   nsets is number og time points in SAR
%
% see also omnibus_loewner
%
% References
%
% Allan A. Nielsen, Henning Skriver and Knut Conradsen (2019).
% "The Loewner Order and Direction of Detected Change in Sentinel-1 and
% Radarsat-2 Data".
% IEEE Geoscience and Remote Sensing Letters 17(2), 242-246 (2020).
% DOI:10.1109/LGRS.2019.2918636.
%
% Allan A. Nielsen (2020).
% "Fast matrix based computation of eigenvalues and the Loewner order in
% PolSAR data".
% IEEE Geoscience and Remote Sensing Letters 17(10), 1727-1731.
% DOI:10.1109/LGRS.2019.2952202.
%
% Allan A. Nielsen, Henning Skriver and Knut Conradsen (2024).
% "Multi-Frequency Omnibus Change Detection in Covariance Matrix PolSAR
% Data".
% IEEE Journal of Selected Topics in Applied Earth Observations and Remote
% Sensing

% wcWhenRj detects change based on Rj only
% wc_when and wcWhenQRj detect change based on both Q and Rj

% (c) Copyright 2016-2024
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, http://people.compute.dtu.dk/alan
% 2 May 2016 - 30 Jan 2019
% 17 Aug 2019 (Loewner order added)
% 30 Jul 2024 (variable aux introduced)

if ~isstruct(wc)
    error('wc must be struct as output from wcRjl')
end
if isfield(wc,'sumMat')
    nsets = size(wc,2);
else
    nsets = size(wc,2) + 1;
end
k = nsets;

if ismatrix(wc(nsets-1).P) % two time points only: copy lnQ to lnR
    wc(nsets-1).P(:,:,2) = wc(nsets-1).P;
    wc(nsets-1).P(:,:,1) = wc(nsets-1).P(:,:,2);
end

if nargin<3
    pval = 0.99;
end

sizec = size(wc(1).P);
sizec(3) = sizec(3) - 1;
c = false(sizec);
[nr,nc,~] = size(c);
lotf = false(2,nr,nc,k-1);
ilast = zeros(nr,nc);
h = waitbar(0,'Processing...','Name',' wc_when_loewner');
for cc=1:sizec(2)
    waitbar(cc/sizec(2),h)
    for rr=1:sizec(1)
        time = 0;
        ell = 1;
        while and(ell<k,time<k)
            if wc(ell).P(rr,cc,1)>pval % check Q
                for j=2:(k-ell+1)
                    time = time + 1;
                    if wc(ell).P(rr,cc,j)>pval % check Rjl
                        c(rr,cc,time) = true;
%                         fprintf('ell=%i, j=%i, time=%i, P=%f\n',ell,j,time,1-wc(ell).P(rr,cc,j))
                        ell = ell + j - 1;
                        % aux = loewner_order(mean(SAR(rr,cc,:, ...
                        lotf(:,rr,cc,time) = loewner_order_tf(mean(SAR(rr,cc,:, ...
                            (ilast(rr,cc)+1):time),4),SAR(rr,cc,:,time+1));
                        % lotf(1,rr,cc,time) = aux.posdef;
                        % lotf(2,rr,cc,time) = aux.negdef;
                        ilast(rr,cc) = time;
                        break;
                    end
                end
            else
                break;
            end
        end
    end
end
close(h)
lo = struct;
for ix=1:(k-1)
%     lo(ix) = struct('posdef',{false(nr,nc)},'negdef',{false(nr,nc)});
    lo(ix).posdef = squeeze(lotf(1,:,:,ix));
    lo(ix).negdef = squeeze(lotf(2,:,:,ix));
end
