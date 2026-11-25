function [det,negdet] = wishart_det(X)

% Calculate determinant of Hermitian 3x3 or 2x2 matrix
%
% [Xdet,Xnegdet] = wishart_det(X);
%
% Calculate determinants of complex covariance polarimetric SAR data
%
%   full polarimetry
%   full polarimetry, azimutal symmetry
%   full polarimetry, diagonal case
%   dual polarimetry
%   dual polarimetry, diagonal case
%   single polarization
%
% (to read the data from disk, see freadenvisar)
%
% ShhShh* ShhShv* ShhSvv*     k    a   rho
% ShvShh* ShvShv* ShvSvv*  =  a*   ksi b
% SvvShh* SvvShv* SvvSvv*     rho* b*  zeta
%
% Input
%   X       - covariance matrix input at first time point;
%             must be 3-D with either
%               9 (full polarimetry),
%               5 (full, azimuthal symmetry),
%               3 (full, diagonal only),
%               4 (dual polarimetry),
%               2 (dual, diagonal only), or
%               1 band(s);
%             reading order is alphabetical in indexes with
%             real parts before imaginary parts (as read by
%             freadenvisar), fx. for full polarimetry:
%               hhhh (real), hhhv (complex), hhvv (complex),
%               hvhv (real), hvvv (complex), vvvv (real)
%
% Output
%   Xdet    - determinant of X
%   Xnegdet - identifies row and col number of pixels with negative
%             determinants

% (c) Copyright 2012-2024
% Allan Aasbjerg Nielsen
% alan@dtu.dk, people.compute.dtu.dk/alan
% 23 Jan 2012
% 2 Aug 2024 /AA: two 3x3 blocks in 6x6 covmat and two 2x2 blocks in 4x4

if ischar(X), X = freadenvisar(X); end

pX = size(X);
% nrows = pX(1);
% ncols = pX(2);
if  ndims(X)==3
    nvars = pX(3);
elseif ismatrix(X)
    nvars = 1;
else
    error('wishart_det: wrong dimensionality of input');
end

if nvars==9 % fully polarimetric data
%     X(:,:,6) = 2*X(:,:,6);
%     X(:,:,[2 3 7 8]) = sqrt(2)*X(:,:,[2 3 7 8]);
%     Y(:,:,6) = 2*Y(:,:,6);
%     Y(:,:,[2 3 7 8]) = sqrt(2)*Y(:,:,[2 3 7 8]);
%     fprintf('wishart_det: full polarimetry\n');
    % ShhShh* ShhShv* ShhSvv*    k     a    rho   | 1   2,3     4,5
    % ShvShh* ShvShv* ShvSvv* =  a*    ksi  b     |     6       7,8
    % SvvShh* SvvShv* SvvSvv*    rho*  b*   zeta  |             9
    % 
    % det = k*ksi*zeta+a*b*rho*+rho*a**b*-rho^2*ksi-b^2*k-zeta*a^2
    %
    % a*b*rho* = a**b**rho = 2*Re(a*b*rho*)
    det = X(:,:,1).*X(:,:,6).*X(:,:,9) ...
        + 2*(X(:,:,4).*(X(:,:,2).*X(:,:,7)-X(:,:,3).*X(:,:,8)) ...
           + X(:,:,5).*(X(:,:,3).*X(:,:,7)+X(:,:,2).*X(:,:,8))) ...
        - X(:,:,6).*(X(:,:,4).^2+X(:,:,5).^2) ...
        - X(:,:,1).*(X(:,:,7).^2+X(:,:,8).^2) ...
        - X(:,:,9).*(X(:,:,2).^2+X(:,:,3).^2);
elseif nvars==5 % full polarimetry, azimuthal symmetry
%     X(:,:,4) = 2*X(:,:,4);
%     Y(:,:,4) = 2*Y(:,:,4);
%     fprintf('wishart_det: full polarimetry, azimuthal symmetry\n');
    % ShhShh* ShhShv* ShhSvv*    k     0    rho  |  1           2,3
    % ShvShh* ShvShv* ShvSvv* =  0     ksi  0    |     4
    % SvvShh* SvvShv* SvvSvv*    rho*  0    zeta |              5
    % 
    % det = ksi*(k*zeta-rho^2)
    det = X(:,:,4).*(X(:,:,1).*X(:,:,5) - X(:,:,2).^2 - X(:,:,3).^2);
elseif nvars==3 % full polarimetry, diagonal only
%     X(:,:,2) = 2*X(:,:,2);
%     Y(:,:,2) = 2*Y(:,:,2);
%     fprintf('wishart_det: full polarimetry, diagonal only\n');
    % ShhShh* ShhShv* ShhSvv*    k     0    0    |  1
    % ShvShh* ShvShv* ShvSvv* =  0     ksi  0    |     2
    % SvvShh* SvvShv* SvvSvv*    0     0    zeta |        3
    % 
    % det = k*ksi*zeta
    det = X(:,:,1).*X(:,:,2).*X(:,:,3);
elseif nvars==4 % dual polarimetry
%     fprintf('wishart_det: dual polarimetry\n');
    % ShhShh* ShhSvv*    k     rho   | 1   2,3
    % SvvShh* SvvSvv*    rho*  ksi   |     4
    %
    % or
    %
    % ShhShh* ShhShv*    k     a     | 1   2,3
    % ShvShh* ShvShv*    a*    zeta  |     4
    %
    % or
    %
    % ShvShv* ShvSvv*    ksi   b     | 1   2,3
    % SvvShv* SvvSvv*    b*    zeta  |     4
    %
    % det  = k*ksi-rho2;
    det = X(:,:,1).*X(:,:,4) - X(:,:,2).^2 - X(:,:,3).^2;
elseif nvars==2 % dual polarimetry, diagonal only
%     fprintf('wishart_det: dual polarimetry, diagonal only\n');
    % ShhShh* ShhSvv*    k     0     |  1
    % SvvShh* SvvSvv*    0     zeta  |     2
    %
    % det  = k*zeta
    det = X(:,:,1).*X(:,:,2);
elseif nvars==1 % one channel data only
%     fprintf('wishart_det: one channel data\n');
    det = X(:,:,1);
elseif nvars==18 % two 3x3 blocks in 6x6 covmat
    % ShhShh* ShhShv* ShhSvv*    k     a    rho   | 1   2,3     4,5
    % ShvShh* ShvShv* ShvSvv* =  a*    ksi  b     |     6       7,8
    % SvvShh* SvvShv* SvvSvv*    rho*  b*   zeta  |             9
    % 
    % det = k*ksi*zeta+a*b*rho*+rho*a**b*-rho^2*ksi-b^2*k-zeta*a^2
    %
    % a*b*rho* = a**b**rho = 2*Re(a*b*rho*)
    det = X(:,:,1).*X(:,:,6).*X(:,:,9) ...
        + 2*(X(:,:,4).*(X(:,:,2).*X(:,:,7)-X(:,:,3).*X(:,:,8)) ...
           + X(:,:,5).*(X(:,:,3).*X(:,:,7)+X(:,:,2).*X(:,:,8))) ...
        - X(:,:,6).*(X(:,:,4).^2+X(:,:,5).^2) ...
        - X(:,:,1).*(X(:,:,7).^2+X(:,:,8).^2) ...
        - X(:,:,9).*(X(:,:,2).^2+X(:,:,3).^2);
    det1 = X(:,:,9+1).*X(:,:,9+6).*X(:,:,9+9) ...
        + 2*(X(:,:,9+4).*(X(:,:,9+2).*X(:,:,9+7)-X(:,:,9+3).*X(:,:,9+8)) ...
           + X(:,:,9+5).*(X(:,:,9+3).*X(:,:,9+7)+X(:,:,9+2).*X(:,:,9+8))) ...
        - X(:,:,9+6).*(X(:,:,9+4).^2+X(:,:,9+5).^2) ...
        - X(:,:,9+1).*(X(:,:,9+7).^2+X(:,:,9+8).^2) ...
        - X(:,:,9+9).*(X(:,:,9+2).^2+X(:,:,9+3).^2);
    det = det.*det1;
    clear det1
elseif nvars==8 % two 2x2 blocks in 4x4 covmat
    det = X(:,:,1).*X(:,:,4) - X(:,:,2).^2 - X(:,:,3).^2;
    det1 = X(:,:,4+1).*X(:,:,4+4) - X(:,:,4+2).^2 - X(:,:,4+3).^2;
    det = det.*det1;
    clear det1
else
    error('wishart_det: input 18, 8, (9, 5, 3), (4, 2) or 1 band(s), real (and complex)');
%     fprintf('%s\n', '*** order of input images is ***very*** important ***')
%     fprintf('%s\n', '*** hhhh (real), hhhv (complex), hhvv (complex),  ***')
%     fprintf('%s\n', '*** hvhv (real), hvvv (complex), vvvv (real)      ***')
end

if nargout>1
    negdet.r = [];
    negdet.c = [];
    if min(det(:))<=0
%         warning('*** det <= 0, bad *** - set to eps = 2.2204e-16')
        warning('*** det <= 0, bad ***')
        fprintf('min(det(:)) = %g\n', min(det(:)));
        imtool(det<=0,[])
        [negdet.r,negdet.c] = find(det<=0);
%         det(det<=0) = eps;
    end
end
