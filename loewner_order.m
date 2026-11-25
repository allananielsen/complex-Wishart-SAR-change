function lo = loewner_order(X,Y)

% Loewner order of (1x1,) 2x2 or 3x3 Hermitian matrices X and Y: is X > Y?
%
% Loewner order: X > Y if (X-Y) is positive definite
%
% lo = loewner_order(X,Y)
%
% SAR input images are images of the (row-wise) upper-right elements of
% the covariance or coherency matrix, e.g. for full polarimetry:
% a vector with nine real (i.e., non-complex) elements, namely
% [ShhShh* Re(ShhShv*) Im(ShhShv*) Re(ShhSvv*) Im(ShhSvv*) ...
%  ShvShv* Re(ShvSvv*) Im(ShvSvv*) SvvSvv*]
%
% Input
%   X       - 3-D array of (1x1,) 2x2 or 3x3 complex covariance matrix
%             (polarimetric SAR image in covariance matrix form,
%             reading order is alphabetical (as read by freadenvisar):
%               hhhh (real), hhhv (complex), hhvv (complex),
%                            hvhv (real),    hvvv (complex),
%                                            vvvv (real)
%   Y       - 3-D array of (1x1,) 2x2 or 3x3 complex covariance matrix
%
% Output
%   lo - struct with fields
%           posdef	- X-Y is positive definite (logical)
%           negdef	- X-Y is negative definite (logical)
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
% IEEE Geoscience and Remote Sensing Letters (Early Access).
% DOI:10.1109/LGRS.2019.2952202.
%
% Allan A. Nielsen (2020).
% "Fast Matrix Based Computation of Eigenvalues in PolSAR Data".
% Accepted for IEEE International Geoscience and Remote Sensing Symposium,
% IGARSS 2020, Hawaii, USA, 19-24 July 2020.

% (c) Copyright 2018-2019
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 3-21 Sep 2018
% 13 Apr 2019 (fast, matrix based 2x2 and 3x3 azimuthal symmetry)
% 6-9 Jul 2019 (fast, matrix based full 3x3)
% 2 Aug 2019 (finally, pivot based)

if nargin<2
    help loewner_order
    error('loewner_order: wrong input');
end
if nargout>1
    help loewner_order
    error('loewner_order: wrong output');
end

pX = size(X);
pY = size(Y);
if ~all(pX==pY)
    error('loewner_order: X and Y must have same dimensions');
end
if  ndims(X)==3
    nvars = pX(3);
elseif ismatrix(X)
    warning('scalar input, no reason to use loewner_order')
    nvars = 1;
else
    error('loewner_order: wrong sizes of input images');
end

X = X-Y;
clear Y

if nvars==9 % full pol
    d1 = X(:,:,1);
    d2 = wishart_det(X(:,:,[1 2 3 6]));
    d3 = wishart_det(X);
    lo.posdef = d1>0 & d2>0 & d3>0;
    lo.negdef = d1<0 & d2>0 & d3<0;
elseif nvars==4 % dual pol
    d1 = X(:,:,1);
    d2 = wishart_det(X);
    lo.posdef = d1>0 & d2>0;
    lo.negdef = d1<0 & d2>0;
elseif nvars==2 % dual pol, diagonal only
    lo.posdef = X(:,:,1)>0 & X(:,:,2)>0;
    lo.negdef = X(:,:,1)<0 & X(:,:,2)<0;
elseif nvars==1 % single channel
    lo.posdef = X>0;
    lo.negdef = X<0;
elseif nvars==5 % full pol, azimuthal symmetry
    d1 = X(:,:,1);
    d2 = X(:,:,1).*X(:,:,4);
    d3 = wishart_det(X);
    lo.posdef = d1>0 & d2>0 & d3>0;
    lo.negdef = d1<0 & d2>0 & d3<0;
elseif nvars==3 % full pol, diagonal only
    lo.posdef = X(:,:,1)>0 & X(:,:,2)>0 & X(:,:,3)>0;
    lo.negdef = X(:,:,1)<0 & X(:,:,2)<0 & X(:,:,3)<0;
else
    error('loewner_order: wrong dimensionality of input images')
end
