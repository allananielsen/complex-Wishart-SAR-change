function [lnR,f,rho,omega2,P,sumMat,logDet] = wishart_change_Rj_updt(X,Xnew,L,varargin)

% Update ln(Rj) for change detection in polarimetric SAR data
%
% detect difference/change in k 2x2 or 3x3 complex Wishart distributed matrices
% used for change detection in polarimetric SAR images
%
% [lnR,f,rho,omega2,P] = wishart_change_Rj_updt(X,Xnew,L,pol,logDetLXnew);
%
% Examples:
% [lnR,f,rho,omega2,P] = wishart_change_Rj_updt(cat(4,L62,L63,L64,L65,L68),L74,13,'full');
%
% Wishart distribution based test for equality of several complex 2x2 or
% 3x3 variance-covariance matrices (e.g., for change or edge detection)
% with test for when first inequality occurs. If an inequality occurs,
% the prerequsites for the following tests are broken.
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
% Input
%   X       - 4-D array series of k 2x2 or 3x3 comples covariance matrices
%             (polarimetric SAR images in covariance matrix form,
%             reading order is alphabetical (as read by freadenvisar):
%               hhhh (real), hhhv (complex), hhvv (complex),
%                            hvhv (real),    hvvv (complex),
%                                            vvvv (real)
%           or
%             struct wc (see wcRjl)
%   Xnew    - calculate additional Rj with this 3-D data set 
%   L       - k (or one) number of looks for X
%   pol     - 'full', 'azim', 'dual', 'fdiag' or 'ddiag'
%             (needed only when confusion is possible in which case it
%             defaults to a value corresponding to the lowest k)
%   logDet  - logarithm of the determinant of the new data in Xnew
%
% Output
%   lnR     - logarithm of test quantities Q and Rj
%             ln(Q) is lnR(:,:,1), ln(R_2) in the article is lnR(:,:,2) etc.
%   f       - number of degrees of freedom for Q and Rj test statistics
%   rho
%   omega2
%   P       - probabilities of finding smaller value of
%             -2*rho(i)*lnR(:,:,i), i=1,...,k
%
% Reference
%
% Knut Conradsen, Allan Aasbjerg Nielsen and Henning Skriver (2016):
% "Determining the Points of Change in Time Series of Polarimetric SAR Data".
% IEEE Transactions on Geoscience and Remote Sensing 54(5), 3007-3024.
% http://www.imm.dtu.dk/pubdb/p.php?6825
% http://orbit.dtu.dk/en/publications/determining-the-points-of-change-in-time-series-of-polarimetric-sar-data(4255ddb4-1fd9-4d66-9aad-e8049a976d75).html

% (c) Copyright 2016-2019
% Allan Aasbjerg Nielsen, PhD
% alan@dtu.dk, http://people.compute.dtu.dk/alan
% 8 Apr-7 Dec 2016
% 28 Aug 2019 /diag only cases: better P 

if nargin<3 || nargin>5
    help wishart_change_Rj_updt
    error('wishart_change_Rj_updt: wrong input');
end
if nargout>7
    help wishart_change_Rj_updt
    error('wishart_change_Rj_updt: wrong output');
end

calcDet = true;

if isstruct(X) %-----------------------------------------------------------
    
if ~isfield(X,'sumMat')
    error('struct must contain field ''sumMat''')
end

pX = size(X(1).sumMat);
% nrows = pX(1);
% ncols = pX(2);
if ismatrix(X(1).sumMat)
    pX(3) = 1;
end
nv = pX(3);
pY = size(Xnew);
if ismatrix(Xnew)
    pY(3) = 1;
end
if pX(3)~=pY(3)
    error('wishart_change_Rj_updt: update with one image Xnew at a time only');
end
if pY(1)~=pX(1) || pY(2)~=pX(2)
    error('wishart_change_Rj_updt: existing data in struct wc and update image Xnew must be same size');
end

% size(X.P)
if ismatrix(X.P)
    if isempty(X.P)
        nsets = 1;
    else
        nsets = 2;
    end
else
    nsets = size(X.P,3);
end
n = L;
Xnew = L*Xnew;

if nv==9
    pol = 'full';
elseif nv==5
    pol = 'azim';
elseif nv==4
    pol = 'dual';
elseif nv==3
    pol = 'fdiag';
elseif nv==2
    pol = 'ddiag';
else
    pol = 'single';
%     nv = 1;
% %     error('input must contain a multiplicative of 2, 3, 4, 5 or 9 variables');
end
if nargin>3
    pol = varargin{1};
    if strncmp(pol,'s',1)
        nv = 1;
    elseif strncmp(pol,'dd',2)
        nv = 2;
    elseif strncmp(pol,'fd',2)
        nv = 3;
    elseif strncmp(pol,'d',1)
        nv = 4;
    elseif strncmp(pol,'a',1)
        nv = 5;
    elseif strncmp(pol,'f',1)
        nv = 9;
    end
%     if rem(nvars,nv)~=0
%         error('wishart_change_Rj_updt: number of input frames and polarization do not match')
%     end
end
if nargin>4
    calcDet = false;
    logDetLXnew = varargin{2};
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
else % if nv==9
    p = 3;
    f = p^2;
    sump = p*(2*p^2-1);
    sump2 = p^2*(p^2-1); 
end

j = nsets+1;
cst = p*(j*log(j)-(j-1)*log(j-1));
sumMat = X.sumMat + Xnew;
logDet = log(wishart_det(sumMat));
if calcDet
    logDetLXnew = log(wishart_det(Xnew));
end
clear Xnew
lnR = n*(cst + (j-1)*(X.logDet) + logDetLXnew - j*logDet);

else % isstruct(X) --------------------------------------------------------
    
pX = size(X);
nrows = pX(1);
ncols = pX(2);
if ismatrix(X)
    nv = 1;
else
    nv = pX(3);
end
if ndims(X)==3
    pX(4) = 1;
end
nsets = pX(4);
pY = size(Xnew);
nrowY = pY(1);
ncolY = pY(2);
if ismatrix(Xnew)
    nvY = 1;
else
    nvY = pY(3);
end
if nv~=nvY
    error('wishart_change_Rj_updt: X and Xnew are not same mode or more than one Xnew image')
end
if nrows~=nrowY || ncols~=ncolY
    error('wishart_change_Rj_updt: X and Xnew don''t match')
end

if nv==9
    pol = 'full';
elseif nv==5
    pol = 'azim';
elseif nv==4
    pol = 'dual';
elseif nv==3
    pol = 'fdiag';
elseif nv==2
    pol = 'ddiag';
else
    pol = 'single';
%     nv = 1;
% %     error('input must contain a multiplicative of 2, 3, 4, 5 or 9 variables');
end
if nargin>3
    pol = varargin{1};
    if strncmp(pol,'s',1)
        nv = 1;
    elseif strncmp(pol,'dd',2)
        nv = 2;
    elseif strncmp(pol,'fd',2)
        nv = 3;
    elseif strncmp(pol,'d',1)
        nv = 4;
    elseif strncmp(pol,'a',1)
        nv = 5;
    elseif strncmp(pol,'f',1)
        nv = 9;
    end
%     if rem(nvars,nv)~=0
%         error('wishart_change_Rj_updt: number of input frames and polarization do not match')
%     end
end
% calcDet = true;
if nargin>4
    calcDet = false;
    logDetLXnew = varargin{2};
end

% % fprintf('%s\n', '*** order of input images is ***very*** important ***')
% 
% if isscalar(L)
%     L = repmat(L,nsets,1);
% elseif ~isscalar(unique(L))
%     error('wishart_change_Rj_updt: number of looks must be equal for all matrices')
% end
% if length(L)~=nsets
%     error('wishart_change_Rj_updt: wishart_change_Rj: X and L do not match')
% end
for i=1:nsets
    X(:,:,:,i) = L(1)*X(:,:,:,i);
end
% N = sum(L);
Xnew = L(end)*Xnew;

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
else % if nv==9
    p = 3;
    f = p^2;
    sump = p*(2*p^2-1);
    sump2 = p^2*(p^2-1); 
end

j = nsets+1;
cst = p*(j*log(j)-(j-1)*log(j-1));
sumMat = sum(X,4);
clear X
Det1 = wishart_det(sumMat);
sumMat = sumMat + Xnew;
logDet = log(wishart_det(sumMat));
if calcDet
    logDetLXnew = log(wishart_det(Xnew));
end
clear Xnew
n = L(1);
lnR = n*(cst + (j-1)*log(Det1) + logDetLXnew - j*logDet);
clear Det1

end % isstruct(X) ---------------------------------------------------------

% j = nsets;
rho = 1 - (1+1/(j*(j-1)))/(6*f*n)*sump;
omega2 = -f/4*(1-1/rho)^2 + (1+(2*j-1)/(j^2*(j-1)^2))*sump2/(24*n^2*rho^2);
df = f/2; % this is 0.5*f; for fewer multiplications in call of gammainc
dfp4 = df + 2; % this is 0.5*(f+4)

if nargout>4
    P = (1 - omega2)*gammainc(-rho*lnR,df  ) + ...
             omega2 *gammainc(-rho*lnR,dfp4);
    if ~isreal(P)
        warning('*** some P-values are complex, only real part is output, max(abs(imag(P(:)))): ***')
        max(abs(imag(P(:))))
        P = real(P);
    end
end

