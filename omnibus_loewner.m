function [lo,chng,ilast] = omnibus_loewner(SAR,wc,pval)

% Omnibus Loewner with mean values since last change
% (see also wcWhen_loewner)
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
% IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing

% wcWhenRj detects change based on Rj only
% wc_when and wcWhenQRj detect change based on both Q and Rj

% (c) Copyright 2019-2024
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 15-16 Aug 2019
% 30 Jul 2024 (variable aux introduced)

if nargin<3
    pval = 0.99;
end

[chng,k] = wcWhenQRj(wc,pval);
[nr,nc,~] = size(chng);
lotf = false(2,nr,nc,k-1);
ilast = zeros(nr,nc);
h = waitbar(0,'Processing...','Name','omnibus_loewner');
for ix=1:(k-1)
    for ic=1:nc
        for ir=1:nr
            if chng(ir,ic,ix)
                % lotf(:,ir,ic,ix) = loewner_order_tf(mean(SAR(ir,ic,:, ...
                aux = loewner_order(mean(SAR(ir,ic,:, ...
                    (ilast(ir,ic)+1):ix),4),SAR(ir,ic,:,ix+1));
                lotf(1,ir,ic,ix) = aux.posdef;
                lotf(2,ir,ic,ix) = aux.negdef;
                ilast(ir,ic) = ix;
            end
        end
    end
    waitbar(ix/(k-1),h)
end
close(h)
for ix=1:(k-1)
    lo(ix).posdef = squeeze(lotf(1,:,:,ix));
    lo(ix).negdef = squeeze(lotf(2,:,:,ix));
end
