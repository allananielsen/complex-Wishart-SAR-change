function pe = polsar_cardano(X)

% Calculates PolSAR eigenvalues, entropy and anisotropy
%
% SAR input images X are images of the (row-wise) upper-right elements of
% the covariance or coherency matrix, e.g. for full polarimetry:
% a vector with nine real (i.e., non-complex) elements, namely
% [ShhShh* Re(ShhShv*) Im(ShhShv*) Re(ShhSvv*) Im(ShhSvv*) ...
%  ShvShv* Re(ShvSvv*) Im(ShvSvv*) SvvSvv*]
%
% Input
%   X       - 3-D array of (1x1,) 2x2 or 3x3 complex covariance matrix
%             (polarimetric SAR image in covariance matrix form,
%             reading order is alphabetical (as read by freadenvisar):
%               [ShhShh* Re(ShhShv*) Im(ShhShv*) Re(ShhSvv*) Im(ShhSvv*) ...
%                        ShvShv*                 Re(ShvSvv*) Im(ShvSvv*)
%                                                SvvSvv*]
%
% Output
%   pe - struct with fields
%           eigvals
%           entropy (Shannon entropy)
%           anisotropy
%           anisotropy12
%
% See also polsar_entropy, H_alpha and van_zyl
%
% References
%
% Gerolamo Cardano, "Ars Magna" (1545).
%
% Allan A. Nielsen, Henning Skriver and Knut Conradsen (2019).
% "The Loewner Order and Direction of Detected Change in Sentinel-1 and
% Radarsat-2 Data".
% IEEE Geoscience and Remote Sensing Letters 17(2), 242-246 (2020).
% DOI:10.1109/LGRS.2019.2918636.
%
% Allan A. Nielsen (2019).
% "Fast matrix based computation of eigenvalues and the Loewner order
% in PolSAR data".
% IEEE Geoscience and Remote Sensing Letters 17(10), 1727-1731 (2020).
% DOI:10.1109/LGRS.2019.2952202.
%
% Allan A. Nielsen (2022).
% "Fast Matrix Based Computation of Eigenvalues in PolSAR Data".
% Accepted for IEEE International Geoscience and Remote Sensing Symposium,
% IGARSS 2020, Hawaii, USA, 19-24 July 2020.

% (c) Copyright 2022
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 16 Feb 2022

warning('polsar_eigvals is nearly twice (~1.8) as fast as polsar_cardano for 3x3')

pX = size(X);
nrows = pX(1);
ncols = pX(2);
if  ndims(X)==3
    nvars = pX(3);
elseif ismatrix(X)
    nvars = 1;
else
    error('polsar_cardano: wrong sizes of input images');
end

if nvars==9 % full pol
    p = 3;
elseif nvars==4 % dual pol
    p = 2;
elseif nvars==2 % dual pol, diagonal only
    p = 2;
elseif nvars==3 % full pol, diagonal only
    p = 3;
elseif nvars==5 % full pol, azimuthal symmetry
    p = 3;
elseif nvars==1 % single channel
    p = 1;
else
    error('polsar_cardano: wrong dimensionality of input images')
end

if isa(X,'double')
    pe.eigvals = nan(nrows,ncols,p,'double');
else
    pe.eigvals = nan(nrows,ncols,p,'single');
end

if nvars==4 % dual pol
    disc = sqrt((X(:,:,1)-X(:,:,4)).^2+4*(X(:,:,2).^2+X(:,:,3).^2));
    pe.eigvals(:,:,1) = 0.5*(X(:,:,1)+X(:,:,4) + disc);
%     pe.eigvals(:,:,2) = 0.5*(X(:,:,1)+X(:,:,4) - disc);
    B = X(:,:,1)+X(:,:,4); % trace
    clear X
    pe.eigvals(:,:,2) = pe.eigvals(:,:,1) - disc;
    clear disc
	p1 = pe.eigvals(:,:,1)./B;
	p2 = pe.eigvals(:,:,2)./B;
    clear B
	pe.entropy = (-p1.*log(p1)-p2.*log(p2))/log(2); % log base 2
    clear p1 p2
    pe.anisotropy12 = (pe.eigvals(:,:,1)-pe.eigvals(:,:,2))./ ...
        (pe.eigvals(:,:,1)+pe.eigvals(:,:,2));
%     pe.gini = (p1.*(1-p1)+p2.*(1-p2))*2;
% 	pe.energy = p1.*p1+p2.*p2;
%     clear p1 p2
%     pe.renyi = -log(pe.energy)/log(2); % Rényi or collision entropy
elseif nvars==9% full pol
	a2 = X(:,:,2).^2+X(:,:,3).^2;
	rho2 = X(:,:,4).^2+X(:,:,5).^2;
	b2 = X(:,:,7).^2+X(:,:,8).^2;
	% minus determinant
    D = -(X(:,:,6).*(X(:,:,1).*X(:,:,9)-rho2)-X(:,:,1).*b2-X(:,:,9).*a2 + ...
		2.0*(X(:,:,4).*(X(:,:,2).*X(:,:,7)-X(:,:,3).*X(:,:,8)) + ...
             X(:,:,5).*(X(:,:,3).*X(:,:,7)+X(:,:,2).*X(:,:,8))));
	C = X(:,:,6).*(X(:,:,1)+X(:,:,9))+X(:,:,1).*X(:,:,9)-rho2-a2-b2;
    clear a2 rho2 b2
	% minus trace
    B = -(X(:,:,1)+X(:,:,6)+X(:,:,9));
    clear X
    % A = 1;
% 	Q = (2.0/27.0*(B./A).^3-1.0/3.0*B.*C./A.^2+D./A)/2.0;
	Q = (2.0/27.0*B.^3-1.0/3.0*B.*C+D)/2.0;
    clear D
%   P = (-1.0/3.0*(B./A).^2+C./A)/3.0;
    P = (-1.0/3.0*B.^2+C)/3.0;
    clear C
%     Q2P3 = Q.^2 + P.^3;
%     if (Q2P3>0) % negation of casus irreducibilis
%         imtool(Q2P3)
%     end
%     clear Q2P3
%   % casus irreducibilis: Q2P3 must be non-positive
    sqrtQ2P3 = 1i*sqrt(-Q.^2 - P.^3);
%     sqrtQ2P3 = sqrt(Q.^2 + P.^3);
    c_alpha = -Q + sqrtQ2P3;
%     c_beta  = -Q - sqrtQ2P3;
    clear Q sqrtQ2P3
    c_alpha = cuberoot(c_alpha,1); % do not calculate third cube root
%     c_beta  = cuberoot(c_beta);
    pe.eigvals(:,:,1) = real(c_alpha(:,:,1) - P./c_alpha(:,:,1));
    pe.eigvals(:,:,2) = real(c_alpha(:,:,2) - P./c_alpha(:,:,2));
%     pe.eigvals(:,:,3) = real(c_alpha(:,:,3) - P./c_alpha(:,:,3));
	pe.eigvals(:,:,3) = -sum(pe.eigvals(:,:,1:2),3); % "eigvals" sum to zero
%     pe.eigvals(:,:,1) = real(c_alpha(:,:,1) - P.*c_alpha(:,:,1)./(c_alpha(:,:,1).*conj(c_alpha(:,:,1))));
%     pe.eigvals(:,:,2) = real(c_alpha(:,:,2) - P.*c_alpha(:,:,2)./(c_alpha(:,:,2).*conj(c_alpha(:,:,2))));
%     pe.eigvals(:,:,3) = real(c_alpha(:,:,3) - P.*c_alpha(:,:,3)./(c_alpha(:,:,3).*conj(c_alpha(:,:,3))));
    clear P c_alpha
    pe.eigvals = pe.eigvals - B/3.0;
%     pe.eigvals = sort(pe.eigvals,3,'descend');
	p1 = -pe.eigvals(:,:,1)./B;
	p2 = -pe.eigvals(:,:,2)./B;
	p3 = -pe.eigvals(:,:,3)./B;
    clear B
	pe.entropy = (-p1.*log(p1)-p2.*log(p2)-p3.*log(p3))/log(3); % log base 3
    pe.anisotropy = (pe.eigvals(:,:,2)-pe.eigvals(:,:,3))./ ...
        (pe.eigvals(:,:,2)+pe.eigvals(:,:,3));
    pe.anisotropy12 = (pe.eigvals(:,:,1)-pe.eigvals(:,:,2))./ ...
        (pe.eigvals(:,:,1)+pe.eigvals(:,:,2));
%     pe.gini = (p1.*(1-p1)+p2.*(1-p2)+p3.*(1-p3))*3/2;
% 	pe.energy = p1.*p1+p2.*p2+p3.*p3;
%     clear p1 p2 p3
%     pe.renyi = -log(pe.energy)/log(3); % Rényi or collision entropy
%     pe.cos3t = theta;
elseif nvars==2 % dual pol diag only
    pe.eigvals = sort(X,3,'descend');
    B = X(:,:,1)+X(:,:,2); % trace
    clear X
	p1 = pe.eigvals(:,:,1)./B;
	p2 = pe.eigvals(:,:,2)./B;
    clear B
	pe.entropy = (-p1.*log(p1)-p2.*log(p2))/log(2); % log base 2
    clear p1 p2
    pe.anisotropy12 = (pe.eigvals(:,:,1)-pe.eigvals(:,:,2))./ ...
        (pe.eigvals(:,:,1)+pe.eigvals(:,:,2));
%     pe.gini = (p1.*(1-p1)+p2.*(1-p2))*2;
% 	pe.energy = p1.*p1+p2.*p2;
%     clear p1 p2
%     pe.renyi = -log(pe.energy)/log(2); % Rényi or collision entropy
elseif nvars==3 % full pol diag only
    pe.eigvals = sort(X,3,'descend');
    B = X(:,:,1)+X(:,:,2)+X(:,:,3); % trace
    clear X
%     pe.anisotropy = (pe.eigvals(:,:,2)-pe.eigvals(:,:,3))./ ...
%         (pe.eigvals(:,:,2)+pe.eigvals(:,:,3));
	p1 = pe.eigvals(:,:,1)./B;
	p2 = pe.eigvals(:,:,2)./B;
	p3 = pe.eigvals(:,:,3)./B;
    clear B
	pe.entropy = (-p1.*log(p1)-p2.*log(p2)-p3.*log(p3))/log(3); % log base 3
    clear p1 p2 p3
    pe.anisotropy = (pe.eigvals(:,:,2)-pe.eigvals(:,:,3))./ ...
        (pe.eigvals(:,:,2)+pe.eigvals(:,:,3));
    pe.anisotropy12 = (pe.eigvals(:,:,1)-pe.eigvals(:,:,2))./ ...
        (pe.eigvals(:,:,1)+pe.eigvals(:,:,2));
%     pe.gini = (p1.*(1-p1)+p2.*(1-p2)+p3.*(1-p3))*3/2;
% 	pe.energy = p1.*p1+p2.*p2+p3.*p3;
%     clear p1 p2 p3
%     pe.renyi = -log(pe.energy)/log(3); % Rényi or collision entropy
elseif nvars==5 % azimuthal symmetry
    disc = sqrt((X(:,:,1)-X(:,:,5)).^2+4*(X(:,:,2).^2+X(:,:,3).^2));
    pe.eigvals(:,:,1) = X(:,:,4);
    pe.eigvals(:,:,2) = 0.5*(X(:,:,1)+X(:,:,5) + disc);
%     pe.eigvals(:,:,3) = 0.5*(X(:,:,1)+X(:,:,5) - disc);
    B = X(:,:,1)+X(:,:,4)+X(:,:,5); % trace
    clear X
    pe.eigvals(:,:,3) = pe.eigvals(:,:,2) - disc;
    clear disc
    pe.eigvals = sort(pe.eigvals,3,'descend'); % to be on the safe side
	p1 = pe.eigvals(:,:,1)./B;
	p2 = pe.eigvals(:,:,2)./B;
	p3 = pe.eigvals(:,:,3)./B;
    clear B
	pe.entropy = (-p1.*log(p1)-p2.*log(p2)-p3.*log(p3))/log(3); % log base 3
    clear p1 p2 p3
    pe.anisotropy = (pe.eigvals(:,:,2)-pe.eigvals(:,:,3))./ ...
        (pe.eigvals(:,:,2)+pe.eigvals(:,:,3));
    pe.anisotropy12 = (pe.eigvals(:,:,1)-pe.eigvals(:,:,2))./ ...
        (pe.eigvals(:,:,1)+pe.eigvals(:,:,2));
%     pe.gini = (p1.*(1-p1)+p2.*(1-p2)+p3.*(1-p3))*3/2;
% 	pe.energy = p1.*p1+p2.*p2+p3.*p3;
%     clear p1 p2 p3
%     pe.renyi = -log(pe.energy)/log(3); % Rényi or collision entropy
else % scalar, single channel power data
    warning('polsar_cardano: eigvals and entropy for scalar are meaningless')
    pe.eigvals = X;
    pe.entropy = nan;
    pe.anisotropy = nan;
    clear X
end

% pp = pe.eigvals./repmat(sum(pe.eigvals,3),1,1,p);
% pe.entropy = -sum(pp.*log(pp),3)/log(p);
% % if p==3
% %     pe.entropy = -sum(pp.*log(pp),3)/log(3);
% % elseif p==2
% %     pe.entropy = -sum(pp.*log(pp),3)/log(2);
% % end

while 0

cd Foulum
L62 = flipud(freadenvisar('fl062_l/fl062_l','b'));
L63 = flipud(freadenvisar('fl063_l/fl063_l','b'));
L64 = flipud(freadenvisar('fl064_l/fl064_l','b'));
L65 = flipud(freadenvisar('fl065_l/fl065_l','b'));
L68 = flipud(freadenvisar('fl068_l/fl068_l','b'));
L74 = flipud(freadenvisar('fl074_l/fl074_l','b'));
cd ..

SAR = cat(4,L62,L63,L64,L65,L68,L74);
clear L62 L63 L64 L65 L68 L74
[nr,nc,~,p] = size(SAR);
L = 13;

% SARcoh = cov2coh(SAR);

for ix=1:p
    img(ix).rgb = log(SAR(:,:,[6 1 9],ix));
%     img(ix).rgb = log(SARcoh(:,:,[6 9 1],ix));
end
% same stretch for R, for G and for B
for ie = 1:3
    ev = [];
    for ix = 1:p
        ev = cat(3,ev,img(ix).rgb(:,:,ie));
    end
    mm = mean(reshape(ev,nr*nc,p));
    ss = std(reshape(ev,nr*nc,p));
    for ix = 1:p
        img(ix).rgb(:,:,ie) = (img(ix).rgb(:,:,ie)-mm(ix))/ss(ix);
    end
end
clear ev
%
figure('name','RGB, same stretch')
nstd = 2;%1.5;
for ix = 1:p
    h(ix) = subplot_tight(2,3,ix,[0.005,0.005]);
    imshow(uint8(255*((img(ix).rgb+nstd)/(2*nstd))))
end
linkaxes(h,'xy')

tic
for ix = 1:p
    pe(ix) = polsar_cardano(SAR(:,:,:,ix));
end
toc
%
figure('name','Eigvals')
for ix = 1:6
    ax(ix) = subplot_tight(2,3,ix,[0.005,0.005]);
    imshowrgb(log(pe(ix).eigvals),[1 2 3],1.5)
end
linkaxes(ax,'xy')

% same stretch
for ix = 1:p
%     pe(ix) = polsar_cardano(SAR(:,:,:,ix));
    pe(ix).eigvals = log(pe(ix).eigvals);
end
for ie = 1:3
    ev = [];
    for ix = 1:p
        ev = cat(3,ev,pe(ix).eigvals(:,:,ie));
    end
    mm = mean(reshape(ev,nr*nc,p))
    ss = std(reshape(ev,nr*nc,p))
    for ix = 1:p
        pe(ix).eigvals(:,:,ie) = (pe(ix).eigvals(:,:,ie)-mm(ix))/ss(ix);
    end
end
clear ev
%
figure('name','Eigvals, same stretch')
nstd = 1.5;
for ix = 1:6
    h(ix) = subplot_tight(2,3,ix,[0.005,0.005]);
    imshow(uint8(255*((pe(ix).eigvals+nstd)/(2*nstd))))
end
linkaxes(h,'xy')

% pairwise Loewner
for ix=1:p-1
    lop(ix) = loewner_order(SAR(:,:,:,ix),SAR(:,:,:,ix+1));
end
figure
for ix=1:p-1
    ax(ix) = subplot_tight(2,3,ix,[0.005,0.005]);
    imshow(255*cat(3,~lop(ix).negdef,~lop(ix).posdef,zeros(nr,nc)));
end
linkaxes(ax,'xy')

tic
wc = wcRjl(SAR,L);
toc

pval = 0.99;

% pairwise comparisons
figure('Name','Pairwise Loewner')
for ix=1:p-1
    ixpair = 2; if ix==p-1; ixpair = 1; end
    axp(ix) = subplot_tight(2,3,ix,[0.005,0.005]);
    imshow(255*cat(3,~lop(ix).negdef & wc(ix).P(:,:,ixpair)>pval, ...
        ~lop(ix).posdef & wc(ix).P(:,:,ixpair)>pval,zeros(nr,nc)));
end
linkaxes(axp,'xy')

% omnibus comparisons
[change,nsets] = wcWhenQRj(wc,pval);

for ix=1:(p-1)
%     lo(ix) = loewner_order(mean(SAR(:,:,:,1:(ix)),4),SAR(:,:,:,ix+1));
    lo(ix) = loewner_order(SAR(:,:,:,ix),SAR(:,:,:,ix+1));
end
%
figure('Name','Omnibus Loewner')
for ix=1:(p-1)
    axo(ix) = subplot_tight(2,3,ix,[0.005,0.005]);
    imshow(255*cat(3,~lo(ix).negdef & change(:,:,ix), ...
        ~lo(ix).posdef & change(:,:,ix),zeros(nr,nc)))
end
linkaxes(axo,'xy')

% Correct omnibus comparisons
tic,[loomni,change] = omnibus_loewner(SAR,wc,pval);toc
% or
tic,[loomni,change,ilast] = wcWhen_loewner(SAR,wc,pval);toc
figure('Name','Omnibus Loewner (with mean values)')
for ix=1:(p-1)
    axl(ix) = subplot_tight(2,3,ix,[0.005,0.005]);
    imshow(255*cat(3,~loomni(ix).negdef & change(:,:,ix), ...
        ~loomni(ix).posdef & change(:,:,ix),zeros(nr,nc)))
end
linkaxes(axl,'xy')
linkaxes([axp,axo,axl],'xy')

end
