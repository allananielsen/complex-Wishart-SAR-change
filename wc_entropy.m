function pP = wc_entropy(wc,roi,varargin)

% wc_Phisto takes P from struct output from wc_Q or similar as input and outputs histogram
% of P-values, i.e., change probabilities for region of interest, roi
%
% wc_Phisto(wcQ.P,roi); % wcQ is output from wc_Q
% wc_Phisto(wcQ,roi,txt);
% wc_Phisto(wcQR(2),roi,txt); % wcQR is output from wcRjl
% wc_Phisto(P,roi);
% wc_Phisto(wc,nan,txt,1,0:0.01:1);
%
% input wc is (P) from wc_Q, wcRjl or similar, roi (logical) is region of interest
%   if roi is colour image, draw ROI with colour image as backdrop
% third input parameter is text for histogram plot
% fourth input parameter is index for lnR_j if input wc is from wcRjl
%   (ignored if input wc is from wc_Q)
%
% pP contains mean, median and entropy measures Shannon, Gini, Rényi
% entropy measures flatness of histogram (flat is 1, one peak is 0)

% mean (1-P)-values, i.e., no-change probabilities for region of interest, roi

% (c) Copyright 2023-2025
% Allan Aasbjerg Nielsen, PhD
% alan@dtu.dk, http://people.compute.dtu.dk/alan
% 1-8 Nov 2023
% 16 Jan 2025 /AA: isnan (roi) added
% 5 Feb 2025 /AA: bins added

if nargin<1
    error('wc_Phisto: wrong input')
end
if nargin>5
    error('wc_Phisto: wrong input')
end

if nargin==1
    roi = NaN;
end

idx = 1;
if nargin>3
    idx = varargin{2};
end

if isstruct(wc)
    % if isfield(wc(1),'P')
        wcnew.P = wc(1).P(:,:,idx);
    % end
else
    wcnew.P = wc;
end
clear wc

if ~ismatrix(roi)
    figure, roi = roipolysarrgb(roi);
elseif isnan(roi)
    [nr,nc,~] = size(wcnew.P);
    roi = true(nr,nc);
end

txt = "";
if nargin>2
    txt = varargin{1};
end

bins = 0:0.05:1;
if nargin>4
    bins = varargin{3};
end

if size(wcnew.P) ~= size(roi)
    error('wc_Phisto: image and roi must be same size')
end
        
% change probability
Paux = wcnew.P(roi);
% % no-change probability
% Paux = 1 - P(roi);
% % figure,imshow(Paux.*(1-double(roi)),[0 1])
% figure('Name',txt)
figure(19520816)
if isnan(bins)
    hh = histogram(Paux,'Visible','off');
else
    hh = histogram(Paux,bins,'Visible','off');
end
% hh = histogram(Paux,bins);
% hh = histogram(Paux); hh.EdgeColor = 'none';
% % hh = histogram(Paux,0:0.05:1,'Normalization','pdf');
% grid on, grid minor
% if size(Paux(:),1)<100*1000-1
%     hold on
%     yy = size(hh.Data,1)/hh.NumBins;
%     plot([0 1],[yy yy],'LineWidth',2)
%     % hold off
% end
aux = hh.Values/hh.NumBins;
probs = aux/nansum(aux);
shannon = -nansum(probs.*log(probs))/log(hh.NumBins);
gini = nansum(probs.*(1-probs));
gini = gini*hh.NumBins/(hh.NumBins-1); % norm to max 1
renyi = -log(nansum(probs.*probs))/log(hh.NumBins);
pP = [mean(Paux), median(Paux), shannon, gini, renyi];
% fprintf('%s: mean = %f, median = %f, Shannon = %f, Gini = %f Rényi = %f\n', ...
%     txt, mean(Paux), median(Paux), shannon, gini, renyi)

close(19520816)
