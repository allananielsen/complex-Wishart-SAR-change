function [pE,pP,roi] = wc_Ptable(wc,roi,varargin)

% wc_Ptable takes struct from wcRjl as input and plots histograms
% and outputs tables of entropies (pE) of change P-histograms, and
% mean change P-values (pP) for a region of interest, roi
%
% wc_Ptable(wc,roi);
% wc_Ptable(wc,roi,txt);
%
% input wc is struct from wcRjl, roi (logical) is region of interest
%   if roi is colour image, draw ROI with colour image as backdrop
%   txt is for plot
% if number of input parameters is >3, set rho = 1 and omega2 = 0

% (c) Copyright 2016-2025
% Allan Aasbjerg Nielsen, PhD
% alan@dtu.dk, http://people.compute.dtu.dk/alan
% 14-17 Jan 2016
% Nov 2023 /AA: now outputs and draws P (and not 1-P)
% 24 Nov 2023 /AA: grid on and linkaxes added
%                  if large dataset, no edge colour on bars
% 30 Jul 2024 /AA: rho=1, omega2=0 fixed
% 20 Dec 2024 /AA: Shannon entropy pE added
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
    figure, roi = roipoly(roi);%roipolysarrgb(roi);
elseif isnan(roi)
    [nr,nc,~] = size(wc(1).P(:,:,1));
    roi = true(nr,nc);
end

if size(wc(nsets-1).P) ~= size(roi)
    error('wc_Ptable: image and roi must be same size')
end

txt = "";
if nargin>2
    txt = varargin{1};
end

k = nsets;
pP = nan(k,k-1);
pE = pP; % Shannon entropy

if nargin<4
    
    if ismatrix(wc(nsets-1).P) % two time points only: copy lnQ to lnR
        wc(nsets-1).P(:,:,2) = wc(nsets-1).P;
        wc(nsets-1).P(:,:,1) = wc(nsets-1).P(:,:,2);
    end
    
    figure('Name',txt)
    for ell=1:(k-1) % Q
        % % no-change probability
        % Paux = 1 - wc(ell).P(:,:,1);
        % change probability
        Paux = wc(ell).P(:,:,1);
        Paux = Paux(roi);
        % figure,imshow(Paux.*(1-double(roi)),[0 1])
        hQ(ell) = subplot(k,k-1,ell);
        hh = histogram(Paux,0:0.05:1);
        % hh = histogram(Paux); hh.EdgeColor = 'none';
        grid on, grid minor
        if size(Paux(:),1)<100*1000-1
            hold on
            yy = size(hh.Data,1)/hh.NumBins;
            plot([0 1],[yy yy],'LineWidth',1.5)
            % hold off
        end
        pP(1,ell) = mean(Paux);
        aux = hh.Values/hh.NumBins;
        probs = aux/nansum(aux);
        pE(1,ell) = -nansum(probs.*log(probs))/log(hh.NumBins);
    end
    ix = 1;
    for ell=1:(k-1)
        m = 2;
        for r=(1+ell):k
            n = (k-1)*(r-1)+1+ell-1;
            %         fprintf('r=%i, ell=%i, n=%i\n',r,ell,n)
            % Paux = 1 - wc(ell).P(:,:,m);
            Paux = wc(ell).P(:,:,m);
            Paux = Paux(roi);
            % figure,imshow(Paux.*(1-double(roi)),[0 1])
            hR(ix) = subplot(k,k-1,n);
            hh = histogram(Paux,0:0.05:1);
            % hh = histogram(Paux); hh.EdgeColor = 'none';
            grid on, grid minor
            if size(Paux(:),1)<100*1000-1
                hold on
                yy = size(hh.Data,1)/hh.NumBins;
                plot([0 1],[yy yy],'LineWidth',1.5)
                % hold off
            end
            pP(r,ell) = mean(Paux);
            aux = hh.Values/hh.NumBins;
            probs = aux/nansum(aux);
            pE(r,ell) = -nansum(probs.*log(probs))/log(hh.NumBins);
            ix = ix + 1;
            m = m + 1;
        end
    end
    linkaxes(hQ,'xy')
    linkaxes(hR,'xy')
    % linkaxes([hQ,hR],'xy')
  
else % rho = 1, omega2 = 0

    % if ismatrix(wc(nsets-1).P) % two time points only: copy lnQ to lnR
    %     wc(nsets-1).P(:,:,2) = wc(nsets-1).P;
    %     wc(nsets-1).P(:,:,1) = wc(nsets-1).P(:,:,2);
    % end
    if isscalar(wc(nsets-1).f) % two time points only: copy lnQ to lnR
        wc(nsets-1).lnR(:,:,2) = wc(nsets-1).lnR;
        wc(nsets-1).lnR(:,:,1) = wc(nsets-1).lnR(:,:,2);
        wc(nsets-1).f(2) = wc(nsets-1).f;
        wc(nsets-1).f(1) = wc(nsets-1).f(2);
    end
    
    figure('Name',txt)
    % change probability
    for ell=1:(k-1) % Q
        Paux = wc(ell).lnR(:,:,1);
        Paux = Paux(roi);
        Paux = gammainc(-Paux,0.5*wc(ell).f(1));
        % figure,imshow(Paux.*(1-double(roi)),[0 1])
        hQ(ell) = subplot(k,k-1,ell);
        hh = histogram(Paux,0:0.05:1);
        % hh = histogram(Paux); hh.EdgeColor = 'none';
        grid on, grid minor
        if size(Paux(:),1)<100*1000-1
            hold on
            yy = size(hh.Data,1)/hh.NumBins;
            plot([0 1],[yy yy],'LineWidth',1.5)
            % hold off
        end
        pP(1,ell) = mean(Paux);
        aux = hh.Values/hh.NumBins;
        probs = aux/nansum(aux);
        pE(1,ell) = -nansum(probs.*log(probs))/log(hh.NumBins);
    end
    ix = 1;
    for ell=1:(k-1)
        m = 2;
        for r=(1+ell):k
            n = (k-1)*(r-1)+1+ell-1;
            %         fprintf('r=%i, ell=%i, n=%i\n',r,ell,n)
            Paux = wc(ell).lnR(:,:,m);
            Paux = Paux(roi);
            Paux = gammainc(-Paux,0.5*wc(ell).f(m));
            % figure,imshow(Paux.*(1-double(roi)),[0 1])
            hR(ix) = subplot(k,k-1,n);
            hh = histogram(Paux,0:0.05:1);
            % hh = histogram(Paux); hh.EdgeColor = 'none';
            grid on, grid minor
            if size(Paux(:),1)<100*1000-1
                hold on
                yy = size(hh.Data,1)/hh.NumBins;
                plot([0 1],[yy yy],'LineWidth',1.5)
                % hold off
            end
            pP(r,ell) = mean(Paux);
            aux = hh.Values/hh.NumBins;
            probs = aux/nansum(aux);
            pE(r,ell) = -nansum(probs.*log(probs))/log(hh.NumBins);
            ix = ix + 1;
            m = m + 1;
        end
    end
    linkaxes(hQ,'xy')
    linkaxes(hR,'xy')
    % linkaxes([hQ,hR],'xy')

end
fprintf('you might want to use xlim or ylim, fx ylim([0 100])\n')

% aux = hh.Values/hh.NumBins;
% probs = aux/nansum(aux);
% shannon = -nansum(probs.*log(probs))/log(hh.NumBins);
% % gini = nansum(probs.*(1-probs));
% % gini = gini*hh.NumBins/(hh.NumBins-1); % norm to max 1
% % renyi = -log(nansum(probs.*probs))/log(hh.NumBins);
% % pP = [mean(Paux), median(Paux), shannon, gini, renyi];

return

matrix2latex(pP,3) % 3 digits after the decimal point

figure%('Name','Grass')
plot((3.5:1:7.5).',pP(1,:),'o-k','LineWidth',2);
grid on, grid minor
hold on
plot((3.5:1:7.5).',pP(2:end,:),'o-','LineWidth',2);
axis0 = axis;
% axis([2 9 axis0(3:4)])
axis([2 9 0 axis0(4)])
plot([2 9],[0.01 0.01],'--')
plot([2 9],[0.05 0.05],'--')
xlabel('Month')
ylabel('No-change prob.')
legend('Q^{(l)}','R_j^{(1)}','R_j^{(2)}','R_j^{(3)}','R_j^{(4)}','R_j^{(5)}',...
    'Location','NorthEast')
%     'Location','SouthEast')
%     'Location','NorthWest')
%     'Location','SouthWest')
