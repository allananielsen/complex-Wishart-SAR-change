function [cmax,nsets] = wcMaxQ(wc,pval)

% wcMaxQ takes struct wc from wcRjl as input and outputs time point before max change
%
% c = wcMaxQ(wc,0.9999);
%
% input wc is struct from wcRjl, pval is probability threshold for change prob, default 0.99
% output cmax is a matrix with elements i if max change occurs after time point i
%
% Because the mapping from -2*rho*lnR to P saturates for high values of
% -2*rho*lnR, i.e., very significant change, use wcMaxQ and not wcMaxP
% if you can (i.e., if lnR and rho are stored in wc)

% (c) Copyright 2019
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, http://people.compute.dtu.dk/alan
% 5 Feb 2019

if ~isstruct(wc)
    error('wc must be struct as output from wcRjl')
end
if ~and(isfield(wc,'lnR'),isfield(wc,'rho'))
    error('fields lnR and rho not in struct wc, use wcMaxP instead')
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
    wc(nsets-1).lnR(:,:,2) = wc(nsets-1).lnR;
    wc(nsets-1).lnR(:,:,1) = wc(nsets-1).lnR(:,:,2);
    wc(nsets-1).rho(2) = wc(nsets-1).rho;
    wc(nsets-1).rho(1) = wc(nsets-1).rho(2);
end

if nargin<2
    pval = 0.99;
end

sizec = size(wc(1).P);
sizec(3) = sizec(3) - 1;
cmax = zeros(sizec(1),sizec(2));

h = waitbar(0,'Processing...','Name','wcMaxQ');
for cc=1:sizec(2)
    waitbar(cc/sizec(2),h)
    for rr=1:sizec(1)
        time = 0;
        ell = 1;
        m2rlnRmax = -Inf;
        while and(ell<k,time<k)
            if wc(ell).P(rr,cc,1)>pval % check Q
                for j=2:(k-ell+1)
                    time = time + 1;
                    if wc(ell).P(rr,cc,j)>pval && -2*wc(ell).rho(j)*wc(ell).lnR(rr,cc,j)>m2rlnRmax % check Rjl
                        m2rlnRmax = -2*wc(ell).rho(j)*wc(ell).lnR(rr,cc,j);
                        cmax(rr,cc) = time;
%                         fprintf('ell=%i, j=%i, time=%i, P=%f\n',ell,j,time,1-wc(ell).P(rr,cc,j))
                        ell = ell + j - 1;
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

figure
imshow(cmax,[0 nsets-1])
% title('cMaxP')
colormap(gca,jet0)
ha = colorbar('southoutside');
set(ha,'XTick',0:(nsets-1),'Position',[0.10 0.06 0.80 0.025]) % [left bottom width height]
% set(h,'XTick',0:(nsets-1),'Position',[0.22 0.91 0.7 0.025],'Color','White')

% set(ha,'XTick',0:(nsets-1))

return
