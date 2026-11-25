function [plnQ,roi] = wc_lnQhisto(wc,roi,varargin)

% wc_lnQhisto takes f, lnQ etc from wc_q or similar as input plots histograms
% of mean -2*rho*lnQ values for region of interest, roi,
% and theoretical chi-squared distributions
%
% plnQ contains the number of degrees of freedom,
% mean and median of -2*rho*lnQ values for region of interest, roi
%
% wclnQhisto(wcQ,roi); % wcQ is output from wc_Q
% wclnQhisto(wcQ,roi,txt);
% wc_lnQhisto(wcQR,roi,txt); % wcQR is output from wcRjl
% wc_lnQhisto(wcQR(2),roi,txt,3);
%
% wc is a struct from w_Q or similar, roi is region of interest
% third input parameter is text for histogram plot
% fourth input parameter is index for lnR_j if input wc is from wcRjl
%   (ignored if input wc is from wc_Q)

% (c) Copyright 2023
% Allan Aasbjerg Nielsen, PhD
% alan@dtu.dk, http://people.compute.dtu.dk/alan
% 1-8 Nov 2023

if nargin<2
    error('wc_lnQhisto: wrong input')
end
if nargin>5
    error('wc_lnQhisto: wrong input')
end

if ~ismatrix(roi)
    figure, roi = roipolysarrgb(roi);
end

txt = "";
if nargin>2
    txt = varargin{1};
end

idx = 1;
if nargin==4
    idx = varargin{2};
end
if isfield(wc(1),'lnR')
    wcnew = wc(1);
    wcnew.lnQ = wc(1).lnR(:,:,idx);
    wcnew.f = wc(1).f(idx);
    wcnew.rho = wc(1).rho(idx);
    wcnew.omega2 = wc(1).omega2(idx);
else
    wcnew = wc;
end
clear wc

if size(wcnew.lnQ) ~= size(roi)
    error('wc_lnQhisto: image and roi must be same size')
end

if nargin>4
    wcnew.rho = 1;
    wcnew.omega2 = 0;
end

lnRaux = -2*wcnew.rho*wcnew.lnQ(roi);
figure('Name',txt)
h = histogram(lnRaux,'Normalization','pdf');
% title("$"+latex(str2sym(txt))+"$",Interpreter="latex") % doesn't work?
hold on
% xx = h.BinLimits(1):0.1:h.BinLimits(2);
xx = 0:0.1:h.BinLimits(2);
f = (1-wcnew.omega2)*chi2pdf(xx,wcnew.f) ...
    +  wcnew.omega2*chi2pdf(xx,wcnew.f+4);
plot(xx,f,'LineWidth',2)
grid on
hold off
plnQ = [wcnew.f, mean(lnRaux), median(lnRaux)];
fprintf('%s: df = %f, mean = %f, median = %f\n', ...
    txt, wcnew.f, mean(lnRaux), median(lnRaux))
