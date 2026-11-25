function [plnQ,roi] = wc_lnRhisto1(wc,roi,varargin)

% wc_lnRhisto takes f, lnQ etc from wc_Q or similar as input plots histograms
% of -2*rho*lnQ values for region of interest, roi,
% and theoretical chi-squared distributions
%
% plnQ contains the number of degrees of freedom,
% mean and median of -2*rho*lnQ values for region of interest, roi
%
% wc_lnRhisto(wcQ,roi); % wcQ is output from wc_Q
% wc_lnRhisto(wcQ,roi,txt);
% wc_lnRhisto(wcQR,roi,txt); % wcQR is output from wcRjl
% wc_lnRhisto(wcQR(2),roi,txt,2);
%
% wc is a struct from wc_Q, wcRjl or similar, roi is region of interest
% third input parameter is text for histogram plot
% fourth input parameter is index for lnR_j if input wc is from wcRjl
%   (ignored if input wc is from wc_Q)

% (c) Copyright 2023-2025
% Allan Aasbjerg Nielsen, PhD
% alan@dtu.dk, http://people.compute.dtu.dk/alan
% 1-24 Nov 2023
% 16 Jan 2025 /AA: isnan(roi) added

if nargin<2
    error('wc_lnRhisto: wrong input')
end
if nargin>5
    error('wc_lnRhisto: wrong input')
end

if ~ismatrix(roi)
    figure, roi = roipolysarrgb(roi);
elseif isnan(roi)
    [nr,nc,~] = size(wc(1).P(:,:,1));
    roi = true(nr,nc);
end

idx = 1;
if nargin==4
    idx = varargin{2};
end

wcnew = wc(1);
if isfield(wc(1),'lnR')
    wcnew.lnQ = wc(1).lnR(:,:,idx);
    wcnew.f = wc(1).f(idx);
    wcnew.rho = wc(1).rho(idx);
    wcnew.omega2 = wc(1).omega2(idx);
end
clear wc

if size(wcnew.lnQ) ~= size(roi)
    error('wc_lnRhisto: image and roi must be same size')
end

txt = "";
if nargin>2
    txt = varargin{1};
end
if nargin>4
    wcnew.rho = 1;
    wcnew.omega2 = 0;
end

lnRaux = -2*wcnew.rho*wcnew.lnQ(roi);
% figure('Name',txt)
if size(lnRaux(:),1)<100*1000-1
    h = histogram(lnRaux,'Normalization','pdf');
else
    h = histogram(lnRaux,'Normalization','pdf','EdgeColor','none');
end
% title("$"+latex(str2sym(txt))+"$",Interpreter="latex") % doesn't work?
hold on
% xx = h.BinLimits(1):0.1:h.BinLimits(2);
xx = 0:0.1:h.BinLimits(2);
f = (1-wcnew.omega2)*chi2pdf(xx,wcnew.f) ...
    +  wcnew.omega2*chi2pdf(xx,wcnew.f+4);
plot(xx,f,'LineWidth',2)
grid on, grid minor
hold off
plnQ = [wcnew.f, mean(lnRaux), median(lnRaux)];
fprintf('%s: df = %f, mean = %f, median = %f\n', ...
    txt, wcnew.f, mean(lnRaux), median(lnRaux))
