function pe = polsar_eigvals_only(X)

% Calculates PolSAR eigenvalues
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
%
% References
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
% Allan A. Nielsen (2020).
% "Fast Matrix Based Computation of Eigenvalues in PolSAR Data".
% IEEE International Geoscience and Remote Sensing Symposium, IGARSS 2021,
% Brussels, Belgium.

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
else
    error('polsar_eigvals: wrong dimensionality of input images')
end

if isa(X,'double')
    pe.eigvals = nan(nrows,ncols,p,'double');
else
    pe.eigvals = nan(nrows,ncols,p,'single');
end

if nvars==4 % dual pol
    disc = sqrt((X(:,:,1)-X(:,:,4)).^2+4*(X(:,:,2).^2+X(:,:,3).^2));
    pe.eigvals(:,:,1) = 0.5*(X(:,:,1)+X(:,:,4) + disc);
    clear X
    pe.eigvals(:,:,2) = pe.eigvals(:,:,1) - disc;
    clear disc
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
    smp = sqrt(-P);
    theta = acos(Q./P./smp);
%     clear Q P
	pe.eigvals(:,:,1) =  2*smp.*cos( theta    /3.0);
	pe.eigvals(:,:,2) = -2*smp.*cos((theta+pi)/3.0);
% 	pe.eigvals(:,:,3) = -2*smp.*cos((theta-pi)/3.0);
    clear smp
	pe.eigvals(:,:,3) = -sum(pe.eigvals(:,:,1:2),3);
    pe.eigvals = pe.eigvals - B/3.0;
    clear B
end
