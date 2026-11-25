function pe = polsar_eigvals(X)

% Calculates PolSAR eigenvalues, entropy and anisotropy
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
% IEEE Geoscience and Remote Sensing Letters (Early Access).
% DOI:10.1109/LGRS.2019.2952202.
%
% Allan A. Nielsen (2020).
% "Fast Matrix Based Computation of Eigenvalues in PolSAR Data".
% Accepted for IEEE International Geoscience and Remote Sensing Symposium,
% IGARSS 2020, Hawaii, USA, 19-24 July 2020.

% (c) Copyright 2019
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 6 Aug 2019

pX = size(X);
nrows = pX(1);
ncols = pX(2);
if  ndims(X)==3
    nvars = pX(3);
elseif ismatrix(X)
    nvars = 1;
else
    error('polsar_eigvals: wrong sizes of input images');
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
    error('polsar_eigvals: wrong dimensionality of input images')
end

pe.eigvals = nan(nrows,ncols,p);

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
    pe.anisotropy12 = (pe.eigvals(:,:,1)-pe.eigvals(:,:,2))./ ...
        (pe.eigvals(:,:,1)+pe.eigvals(:,:,2));
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
%     casus irreducibilis: Q2P3 must be non-positive
    smp = sqrt(-P);
    theta = acos(Q./P./smp);
    clear Q P
	pe.eigvals(:,:,1) =  2*smp.*cos( theta    /3.0);
	pe.eigvals(:,:,2) = -2*smp.*cos((theta+pi)/3.0);
% 	pe.eigvals(:,:,3) = -2*smp.*cos((theta-pi)/3.0);
    clear smp
	pe.eigvals(:,:,3) = -sum(pe.eigvals(:,:,1:2),3);
    pe.eigvals = pe.eigvals - B/3.0;
	p1 = -pe.eigvals(:,:,1)./B;
	p2 = -pe.eigvals(:,:,2)./B;
	p3 = -pe.eigvals(:,:,3)./B;
    clear B
	pe.entropy = (-p1.*log(p1)-p2.*log(p2)-p3.*log(p3))/log(3); % log base 3
    pe.anisotropy = (pe.eigvals(:,:,2)-pe.eigvals(:,:,3))./ ...
        (pe.eigvals(:,:,2)+pe.eigvals(:,:,3));
    pe.anisotropy12 = (pe.eigvals(:,:,1)-pe.eigvals(:,:,2))./ ...
        (pe.eigvals(:,:,1)+pe.eigvals(:,:,2));
elseif nvars==2 % dual pol diag only
    pe.eigvals = sort(X,3,'descend');
    B = X(:,:,1)+X(:,:,2); % trace
    clear X
	p1 = pe.eigvals(:,:,1)./B;
	p2 = pe.eigvals(:,:,2)./B;
    clear B
	pe.entropy = (-p1.*log(p1)-p2.*log(p2))/log(2); % log base 2
    pe.anisotropy12 = (pe.eigvals(:,:,1)-pe.eigvals(:,:,2))./ ...
        (pe.eigvals(:,:,1)+pe.eigvals(:,:,2));
elseif nvars==3 % full pol diag only
    pe.eigvals = sort(X,3,'descend');
    B = X(:,:,1)+X(:,:,2)+X(:,:,3); % trace
    clear X
	p1 = pe.eigvals(:,:,1)./B;
	p2 = pe.eigvals(:,:,2)./B;
	p3 = pe.eigvals(:,:,3)./B;
    clear B
	pe.entropy = (-p1.*log(p1)-p2.*log(p2)-p3.*log(p3))/log(3); % log base 3
    pe.anisotropy = (pe.eigvals(:,:,2)-pe.eigvals(:,:,3))./ ...
        (pe.eigvals(:,:,2)+pe.eigvals(:,:,3));
    pe.anisotropy12 = (pe.eigvals(:,:,1)-pe.eigvals(:,:,2))./ ...
        (pe.eigvals(:,:,1)+pe.eigvals(:,:,2));
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
    pe.anisotropy = (pe.eigvals(:,:,2)-pe.eigvals(:,:,3))./ ...
        (pe.eigvals(:,:,2)+pe.eigvals(:,:,3));
    pe.anisotropy12 = (pe.eigvals(:,:,1)-pe.eigvals(:,:,2))./ ...
        (pe.eigvals(:,:,1)+pe.eigvals(:,:,2));
else % scalar, single channel power data
    warning('polsar_eigvals: eigvals and entropy for scalar are meaningless')
    pe.eigvals = X;
    pe.entropy = nan;
    pe.anisotropy = nan;
    clear X
end
