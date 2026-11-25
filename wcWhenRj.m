function [c,nsets] = wcWhenRj(wc,pval)

% wcWhenRj takes struct wc from wcRjl as input and outputs occurrence of change
%
% c = wcWhenRj(wc,0.9999);
%
% input wc is struct from wcRjl, pval is probability threshold for change prob, default 0.99
% output c(:,:,i) is array showing if lnR shows change after time point i
%
% wcWhenRj detects change based on Rj only (this code used to be named simply wcWhen)
% wcWhenQRj detects change based on both Q and Rj

% (c) Copyright 2016
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, http://people.compute.dtu.dk/alan
% 2 May 2016

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

if nargin<2
    pval = 0.99;
end

sizec = size(wc(1).P);
sizec(3) = sizec(3) - 1;
c = false(sizec);

h = waitbar(0,'Processing...','Name','wcWhenRj');
for cc=1:sizec(2)
    waitbar(cc/sizec(2),h)
    for rr=1:sizec(1)
        time = 0;
        ell = 1;
        while and(ell<k,time<k)
            for j=2:(k-ell+1) % no check of lnQ
                time = time + 1;
                if wc(ell).P(rr,cc,j)>pval
                    c(rr,cc,time) = true;
%                     fprintf('ell=%i, j=%i, time=%i, P=%f\n',ell,j,time,1-wc(ell).P(rr,cc,j))
                    ell = ell + j - 1;
                    break;
                end
            end
        end
    end
end
close(h)

return

% Change visualizations

% RGB examples
caux = single(c(:,:,1:3));
% caux = single(c(:,:,end-2:end));
% caux = single(c(:,:,[1 3 5]));
% caux = modefilter(single(caux));
figure, imshow(caux), axis off
% clear caux

% change at time point 2 and on with no previous changes
figure, imshow(c(:,:,1))
for k = 2:(nsets-1)%(nsets-1):-1:2
    figure, imshow(~any(c(:,:,1:k-1),3) & c(:,:,k))
end
% % as in:
% figure, imshow(~c(:,:,1) &  c(:,:,2))
% figure, imshow(~c(:,:,1) & ~c(:,:,2) &  c(:,:,3))
% figure, imshow(~c(:,:,1) & ~c(:,:,2) & ~c(:,:,3) &  c(:,:,4))
% figure, imshow(~c(:,:,1) & ~c(:,:,2) & ~c(:,:,3) & ~c(:,:,4) & c(:,:,5))

% first change
[maxch,change1] = max(c,[],3);
change1(maxch==0) = 0;
figure, imshow(change1,[]), axis off
colormap(gca,jet0)
h = colorbar('southoutside');
set(h,'XTick',0:(nsets-1))
% set(h,'XTick',h.Limits(1):h.Limits(2))

% most recent/last change
[maxch,changel] = max(c(:,:,(nsets-1):-1:1),[],3);
changel = nsets - changel; 
changel(maxch==0) = 0;
figure, imshow(changel,[]), axis off
colormap(gca,jet0)
h = colorbar('southoutside');
set(h,'XTick',0:(nsets-1))
% set(h,'XTick',h.Limits(1):h.Limits(2))

% # of changes
sumc = sum(c,3);
figure, imshow(sumc,[]), axis off
colormap(gca,jet0)
h = colorbar('southoutside');
set(h,'XTick',0:(nsets-1))
% set(h,'XTick',h.Limits(1):h.Limits(2))
% sumcmode = modefilter(sumc);
% figure
% imshow(sumcmode,[])
% colormap(gca,jet0)
% h = colorbar;
% set(h,'XTick',0:(nsets-1))
% % set(h,'XTick',h.Limits(1):h.Limits(2))

% % cm = logical(modefilter(c));

% % second change: wrong!!!
% ix = (change1(:)==1);
% [nrows,ncols,aux] = size(c);
% c0 = reshape(c,nrows*ncols,aux);
% c0(ix,:) = false;
% [maxch,change2] = max(c0,[],2);
% change2(maxch==0) = 0;
% change2 = reshape(change2,nrows,ncols);
% figure
% imshow(change2,[])
% colormap(gca,jet0)
% h = colorbar;
% set(h,'XTick',0:(nsets-1))
% % set(h,'XTick',h.Limits(1):h.Limits(2))

% change at time point 2 and on
for k = 1:(nsets-1)%(nsets-1):-1:1
    figure, imshow(c(:,:,k))
end

% grassidx2
%    0.0000    0.0000    0.0000    0.0000    0.0000
%     0.2808       NaN       NaN       NaN       NaN
%     0.0112    0.0638       NaN       NaN       NaN
%     0.0056    0.0319    0.1244       NaN       NaN
%     0.1094    0.2386    0.3879    0.4311       NaN
%     0.0000    0.0000    0.0000    0.0000    0.0000
wc(1).P(1,1,1) = 1-0.0000;
wc(1).P(1,1,2) = 1-0.2808;
wc(1).P(1,1,3) = 1-0.0112;
wc(1).P(1,1,4) = 1-0.0056;
wc(1).P(1,1,5) = 1-0.1096;
wc(1).P(1,1,6) = 1-0.0000;

wc(2).P(1,1,1) = 1-0.0000;
wc(2).P(1,1,2) = 1-0.0638;
wc(2).P(1,1,3) = 1-0.0319;
wc(2).P(1,1,4) = 1-0.2386;
wc(2).P(1,1,5) = 1-0.0000;

wc(3).P(1,1,1) = 1-0.0000;
wc(3).P(1,1,2) = 1-0.1244;
wc(3).P(1,1,3) = 1-0.3879;
wc(3).P(1,1,4) = 1-0.0000;

wc(4).P(1,1,1) = 1-0.0000;
wc(4).P(1,1,2) = 1-0.4311;
wc(4).P(1,1,3) = 1-0.0000;

wc(5).P(1,1,1) = 1-0.0000;

% ryeidx
%  0.0000    0.0001    0.0112    0.0080    0.0208
%     0.4249       NaN       NaN       NaN       NaN
%     0.0026    0.0131       NaN       NaN       NaN
%     0.1400    0.2599    0.3344       NaN       NaN
%     0.0002    0.0023    0.0237    0.0242       NaN
%     0.1758    0.1292    0.0669    0.0468    0.0208
wc(1).P(1,2,1) = 1-0.0000;
wc(1).P(1,2,2) = 1-0.4249;
wc(1).P(1,2,3) = 1-0.0026;
wc(1).P(1,2,4) = 1-0.1400;
wc(1).P(1,2,5) = 1-0.0002;
wc(1).P(1,2,6) = 1-0.1758;

wc(2).P(1,2,1) = 1-0.0001;
wc(2).P(1,2,2) = 1-0.0131;
wc(2).P(1,2,3) = 1-0.2599;
wc(2).P(1,2,4) = 1-0.0023;
wc(2).P(1,2,5) = 1-0.1292;

wc(3).P(1,2,1) = 1-0.0112;
wc(3).P(1,2,2) = 1-0.3344;
wc(3).P(1,2,3) = 1-0.0237;
wc(3).P(1,2,4) = 1-0.0669;

wc(4).P(1,2,1) = 1-0.0080;
wc(4).P(1,2,2) = 1-0.0242;
wc(4).P(1,2,3) = 1-0.0468;

wc(5).P(1,2,1) = 1-0.0208;

% maskwood
% 0.3494    0.4218    0.4230    0.4108    0.4523
%     0.3925       NaN       NaN       NaN       NaN
%     0.5469    0.6235       NaN       NaN       NaN
%     0.4121    0.4085    0.4913       NaN       NaN
%     0.4229    0.4289    0.4629    0.4567       NaN
%     0.4116    0.4099    0.4095    0.4154    0.4523
wc(1).P(1,3,1) = 1-0.3494;
wc(1).P(1,3,2) = 1-0.3925;
wc(1).P(1,3,3) = 1-0.5469;
wc(1).P(1,3,4) = 1-0.4121;
wc(1).P(1,3,5) = 1-0.4229;
wc(1).P(1,3,6) = 1-0.4116;

wc(2).P(1,3,1) = 1-0.4218;
wc(2).P(1,3,2) = 1-0.6235;
wc(2).P(1,3,3) = 1-0.4085;
wc(2).P(1,3,4) = 1-0.4289;
wc(2).P(1,3,5) = 1-0.4099;

wc(3).P(1,3,1) = 1-0.4230;
wc(3).P(1,3,2) = 1-0.4913;
wc(3).P(1,3,3) = 1-0.4629;
wc(3).P(1,3,4) = 1-0.4095;

wc(4).P(1,3,1) = 1-0.4108;
wc(4).P(1,3,2) = 1-0.4567;
wc(4).P(1,3,3) = 1-0.4154;

wc(5).P(1,3,1) = 1-0.4523;
