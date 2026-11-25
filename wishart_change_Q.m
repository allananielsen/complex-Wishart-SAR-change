function [lnQ,f,rho,omega2,P,negdet] = wishart_change_Q(X,L,varargin)

% difference/change detection in k (1x1,) 2x2 or 3x3 complex Wishart distributed matrices
% used for change detection in polarimetric SAR images
% 
% [lnQ,f,rho,omega2,P,negdet] = wishart_change_Q(X,L,pol);
%
% Examples:
% [lnQ,f,rho,omega2,P,negdet] = wishart_change_Q(cat(4,L62,L63,L64,L65,L68,L74),13);
%
% [lnQ6,f6,rho6,omega26,P6] = wishart_change_Q(cat(4,L62,L63,L64,L65,L68,L74),13,'f');
% [lnQ5,f5,rho5,omega25,P5] = wishart_change_Q(cat(4,L63,L64,L65,L68,L74),13,'f');
% [lnQ4,f4,rho4,omega24,P4] = wishart_change_Q(cat(4,L64,L65,L68,L74),13,'f');
% [lnQ3,f3,rho3,omega23,P3] = wishart_change_Q(cat(4,L65,L68,L74),13,'f');
% [lnQ2,f2,rho2,omega22,P2] = wishart_change_Q(cat(4,L68,L74),13,'f');
%
% Wishart distribution based test for equality of several complex 1x1, 2x2 or
% 3x3 variance-covariance matrices (e.g., for change or edge detection).
%
% |ShhShh* ShhShv* ShhSvv*|
% |ShvShh* ShvShv* ShvSvv*|
% |SvvShh* SvvShv* SvvSvv*|_i, i = 1,...,k
%
% |ShhShh* ShhShv*|
% |ShvShh* ShvShv*|_i, i = 1,...,k
%
% |ShvShv* ShvSvv*|
% |SvvShv* SvvSvv*|_i, i = 1,...,k
%
% |ShhShh* ShhSvv*|
% |SvvShh* SvvSvv*|_i, i = 1,...,k
%
% SAR input images are images of the (row-wise) upper-right elements of
% the covariance or coherency matrix, e.g. for full polarimetry:
% a vector with nine real (i.e., non-complex) elements, namely
% [ShhShh* Re(ShhShv*) Im(ShhShv*) Re(ShhSvv*) Im(ShhSvv*) ...
%   ShvShv* Re(ShvSvv*) Im(ShvSvv*) SvvSvv*]
%
% Input
%   X       - 4-D array series of k 1x1, 2x2 or 3x3 comples covariance matrices
%             (polarimetric SAR images in covariance matrix form,
%             reading order is alphabetical (as read by freadenvisar):
%               hhhh (real), hhhv (complex), hhvv (complex),
%                            hvhv (real),    hvvv (complex),
%                                            vvvv (real)
%   L       - k (or one) number of looks for X
%   pol     - 'full', 'azim', 'dual', 'fdiag', 'ddiag' or 'single'
%             (needed only when confusion is possible in which case it
%             defaults to a value corresponding to the lowest k)
%
% Output
%   lnQ     - logarithm of test quantity Q
%   f       - number of degrees of freedom
%   rho
%   omega2
%   P       - probabilities of finding smaller value of -2*rho*lnQ
%   negdet  - identifies row and col number of pixels with non-positive
%             determinants and positive lnQ
%
% Reference
%
% Knut Conradsen, Allan Aasbjerg Nielsen and Henning Skriver (2016):
% "Determining the Points of Change in Time Series of Polarimetric SAR Data".
% IEEE Transactions on Geoscience and Remote Sensing 54(5), 3007-3024.
% http://www.imm.dtu.dk/pubdb/p.php?6825
% http://orbit.dtu.dk/en/publications/determining-the-points-of-change-in-time-series-of-polarimetric-sar-data(4255ddb4-1fd9-4d66-9aad-e8049a976d75).html

% (c) Copyright 2012-2019
% Allan Aasbjerg Nielsen, PhD
% alan@dtu.dk, http://people.compute.dtu.dk/alan
% 3 Jun 2019

% warning('F-distribution for single pol?')

if nargin>3
    help wishart_change_Q
    error('wishart_change_Q: wrong input');
end
if nargout>6
    help wishart_change_Q
    error('wishart_change_Q: wrong output');
end

pX = size(X);
nrows = pX(1);
ncols = pX(2);
nv = pX(3);
nsets = pX(4);
if nsets<2
    error('one set only? possibly wrong polarization')
end
k = nsets;
% fprintf('%s\n', '*** order of input images is ***very*** important ***')

if isscalar(L)
    L = repmat(L,nsets,1);
    n = L(1);
else
    error('number of looks must be equal for all matrices')
end
if length(L)~=nsets
    error('wishart_change_Rj: X and L do not match')
end
for j=1:nsets
    X(:,:,:,j) = L(j)*X(:,:,:,j);
end

if nv==1
    p = 1;
    f = p^2;
    sump = p*(2*p^2-1);
    sump2 = p^2*(p^2-1);
elseif nv==2
    p1 = 1;
    p2 = 1;
    p = p1+p2;
    f = p1^2+p2^2;
    sump = p1*(2*p1^2-1) + p2*(2*p2^2-1);
    sump2 = p1^2*(p1^2-1) + p2^2*(p2^2-1);
elseif nv==3
    p1 = 1;
    p2 = 1;
    p3 = 1;
    p = p1+p2+p3;
    f = p1^2+p2^2+p3^2;
    sump = p1*(2*p1^2-1) + p2*(2*p2^2-1) + p3*(2*p3^2-1);
    sump2 = p1^2*(p1^2-1) + p2^2*(p2^2-1) + p3^2*(p3^2-1);
elseif nv==4
    p = 2;
    f = p^2;
    sump = p*(2*p^2-1);
    sump2 = p^2*(p^2-1);
elseif nv==5
    p1 = 1;
    p2 = 2;
    p = p1+p2;
    f = p1^2+p2^2;
    sump = p1*(2*p1^2-1) + p2*(2*p2^2-1);
    sump2 = p1^2*(p1^2-1) + p2^2*(p2^2-1);
else
    p = 3;
    f = p^2;
    sump = p*(2*p^2-1);
    sump2 = p^2*(p^2-1);
end
f = (k-1)*f;

cst = p*k*log(k);
df = f/2; % this is 0.5*f; for fewer multiplications in call of gammainc
dfp4 = df + 2; % this is 0.5*(f+4)

rho = 1 - (k/n-1/(k*n))*sump/(6*f);
omega2 = -f/4*(1-1/rho)^2 + sump2*(k/n^2-1/n^2/k^2)/(24*rho^2);

logDet1 = zeros(nrows,ncols);
for i=1:nsets
    logDet1 = logDet1 + log(wishart_det(X(:,:,:,i)));
%     logDet1 = logDet1 + log(wishart_det(squeeze(X(:,:,:,i))));
end
lnQ = n*(cst + logDet1 - k*log(wishart_det(sum(X,4))));
if nargout>4
    if ~isreal(lnQ)
        warning('*** some lnQ-values are complex, only real part is used ***')
        fprintf('max(imag(lnQ(:))) = %g\n', max(imag(lnQ(:))));
        lnQ = real(lnQ);
    end
    if max(lnQ(:))>0
        warning('*** lnQ > 0, bad ***')
        fprintf('max(lnQ(:)) = %g\n', max(lnQ(:)));
        if nargout>5
            [negdet.rQ,negdet.cQ] = find(lnQ>0);
            imtool(lnQ>0)
        end
    else
        negdet = nan;
    end
    lnQ(lnQ>0) = 0;
    P = (1 - omega2)*gammainc(-rho*lnQ,df  ) + ...
             omega2 *gammainc(-rho*lnQ,dfp4);
    if ~isreal(P)
        warning('*** some P-values are complex, only real part is output, max(abs(imag(P(:)))): ***')
        max(abs(imag(P(:))))
        P = real(P);
    end
end
