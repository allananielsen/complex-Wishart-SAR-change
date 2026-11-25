function [pQRj,roi] = wc_lnRtable(wc,roi,varargin)

% wc_lnRtable takes struct from wcRjl as input and outputs table of
% mean -2*rho*lnQ/lnRj-values for region of interest, roi
%
% it also plots histograms of 
% mean -2*rho*lnQ/lnRj-values for region of interest, roi,
% and theoretical chi-squared distributions
%
% wc_lnRtable(wc,roi);
% wc_lnRtable(wc,roi,txt);
%
% input wc is struct from wcRjl, roi is region of interest, txt is for plot
% if number of input parameters is >3, set rho = 1 and omega2 = 0

% (c) Copyright 2016-2025
% Allan Aasbjerg Nielsen, PhD
% alan@dtu.dk, http://people.compute.dtu.dk/alan
% 14-17 Jan 2016
% 24 Nov 2023 /AA: grid on and linkaxes added
%                  if large dataset, no edge colour on bars
% 16 Jan 2025 /AA: isnan(roi) added

if ~isstruct(wc)
    error('wc must be struct as output from wcRjl')
end
if isfield(wc,'sumMat')
    nsets = size(wc,2);
else
    nsets = size(wc,2) + 1;
end

if nargin==1
    roi = NaN;
end

if ~ismatrix(roi)
    figure, roi = roipolysarrgb(roi);
elseif isnan(roi)
    [nr,nc,~] = size(wc(1).P(:,:,1));
    roi = true(nr,nc);
end

if size(wc(nsets-1).lnR) ~= size(roi)
    error('wc_lnRtable: image and roi must be same size')
end

txt = "";
if nargin>2
    txt = varargin{1};
end

if isscalar(wc(nsets-1).f) % two time points only: copy lnQ to lnR etc.
    wc(nsets-1).lnR(:,:,2) = wc(nsets-1).lnR;
    wc(nsets-1).lnR(:,:,1) = wc(nsets-1).lnR(:,:,2);
    wc(nsets-1).omega2(2) = wc(nsets-1).omega2;
    wc(nsets-1).omega2(1) = wc(nsets-1).omega2(2);
    wc(nsets-1).rho(2) = wc(nsets-1).rho;
    wc(nsets-1).rho(1) = wc(nsets-1).rho(2);
    wc(nsets-1).f(2) = wc(nsets-1).f;
    wc(nsets-1).f(1) = wc(nsets-1).f(2);
end

k = nsets;
pQRj = nan(k,k-1);

figure('Name',txt)
for ell=1:(k-1) % Q
    if nargin>3
        wc(ell).omega2(1) = 0;
        wc(ell).rho(1) = 1;
    end
    lnRaux = wc(ell).lnR(:,:,1);
    lnRaux = -2*wc(ell).rho(1)*lnRaux(roi);
    % if k==2, subplot(1,1,ell), else subplot(k,k-1,ell), end
    hQ(ell) = subplot(k,k-1,ell);
    if size(lnRaux(:),1)<100*1000-1
        h = histogram(lnRaux,'Normalization','pdf');
    else
        h = histogram(lnRaux,'Normalization','pdf','EdgeColor','none');
    end
    grid on, grid minor
    hold on
%     xx = h.BinLimits(1):0.1:h.BinLimits(2);
    xx = 0:0.1:h.BinLimits(2);
    f = (1-wc(ell).omega2(1))*chi2pdf(xx,wc(ell).f(1)) ...
        +  wc(ell).omega2(1) *chi2pdf(xx,wc(ell).f(1)+4);
    plot(xx,f,'LineWidth',1.5)
    % hold off
    pQRj(1,ell) = mean(lnRaux);
end
% linkaxes(hQ,'xy')
% if k==2, k = 1; end
ix = 1;
for ell=1:(k-1)
    m = 2;
    for r=(1+ell):k
        if nargin>3
            wc(ell).omega2(m) = 0;
            wc(ell).rho(m) = 1;
        end
        n = (k-1)*(r-1)+1+ell-1;
%         fprintf('r=%i, ell=%i, n=%i\n',r,ell,n)
        lnRaux = wc(ell).lnR(:,:,m);
        lnRaux = -2*wc(ell).rho(m)*lnRaux(roi);
        hR(ix) = subplot(k,k-1,n);
        if size(lnRaux(:),1)<100*1000-1
            h = histogram(lnRaux,'Normalization','pdf');
        else
            h = histogram(lnRaux,'Normalization','pdf','EdgeColor','none');
        end
        grid on, grid minor
        hold on
%         xx = h.BinLimits(1):0.1:h.BinLimits(2);
        xx = 0:0.1:h.BinLimits(2);
        f = (1-wc(ell).omega2(m))*chi2pdf(xx,wc(ell).f(m)) ...
            +  wc(ell).omega2(m) *chi2pdf(xx,wc(ell).f(m)+4);
        plot(xx,f,'LineWidth',1.5)
        % hold off
        pQRj(r,ell) = mean(lnRaux);
        ix = ix + 1;
        m = m + 1;
    end
end
% linkaxes([hQ(end),hR],'xy')
linkaxes(hR,'xy')
% xlim([0 45])
% ylim([0 .085])
fprintf('you might want to use xlim or ylim, fx xlim([0 50])\n')
