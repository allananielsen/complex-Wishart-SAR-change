function [lnR,ff,rho,omega2,P,sumMat,logDet] = wishart_change_Rj(X,L)

% Calculate ln(Q) and ln(Rj) for change detection in polarimetric SAR data
%
% detect difference/change in k 2x2 or 3x3 complex Wishart distributed matrices
% used for change detection in polarimetric SAR images
%
% [lnR,f,rho,omega2,P,sumMat,logDet] = wishart_change_Rj(X,L,pol);
%
% Examples:
% [lnR,ff,rho,omega2,P] = wishart_change_Rj(cat(4,L62,L63,L64,L65,L68,L74),13);
% lnQ = sum(lnR(:,:,2:end),3); % this equals lnR(:,:,1)
%
% [lnR6,f6,rho6,omega26,P6] = wishart_change_Rj(cat(4,L62,L63,L64,L65,L68,L74),13,'f');
% [lnR5,f5,rho5,omega25,P5] = wishart_change_Rj(cat(4,L63,L64,L65,L68,L74),13,'f');
% [lnR4,f4,rho4,omega24,P4] = wishart_change_Rj(cat(4,L64,L65,L68,L74),13,'f');
% [lnR3,f3,rho3,omega23,P3] = wishart_change_Rj(cat(4,L65,L68,L74),13,'f');
% [lnR2,f2,rho2,omega22,P2] = wishart_change_Rj(cat(4,L68,L74),13,'f');
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
% SAR input images are images of the (row-wise) upper-right elements of
% the covariance or coherency matrix, e.g. for full polarimetry:
% a vector with nine real (i.e., non-complex) elements, namely
% [ShhShh* Re(ShhShv*) Im(ShhShv*) Re(ShhSvv*) Im(ShhSvv*) ...
%   ShvShv* Re(ShvSvv*) Im(ShvSvv*) SvvSvv*]
%
% Input
%   X       - 4-D array series of k 2x2 or 3x3 complex covariance matrices
%             (polarimetric SAR images in covariance matrix form,
%             reading order is alphabetical (as read by freadenvisar):
%               hhhh (real), hhhv (complex), hhvv (complex),
%                            hvhv (real),    hvvv (complex),
%                                            vvvv (real)
%   L       - k (or one) number of looks for X
%
% Output
%   lnR     - logarithm of test quantities Q and Rj
%             ln(Q) is lnR(:,:,1), ln(R_2) in the article is lnR(:,:,2) etc.
%   ff      - number of degrees of freedom for Q and Rj test statistics
%   rho
%   omega2
%   P       - probabilities of finding smaller value of
%             -2*rho(i)*lnR(:,:,i), i=1,...,k
%   sumMat  - sum of the L*X (for efficient updating)
%   logDet  - logarithm of the determinant of sumMat
%
% lnQ (see wishart_change_fullk.m or wishart_change_dualk.m) is the sum of
% the lnR's, i.e., lnQ = sum(lnR(:,:,2:end),3)
%
% In the call
%
%   [lnR,f,rho,omega2,P] = wishart_change_Rj(cat(4,L63,L64,L65,L68,L74),13);
%
% lnR(:,:,2) equals lnQ from the call
%
%   [lnQ,f,rho,omega2,P] = wishart_change(L63,L64,13);
%
% References
%
% Knut Conradsen, Allan Aasbjerg Nielsen and Henning Skriver (2016).
% "Determining the Points of Change in Time Series of Polarimetric SAR Data".
% IEEE Transactions on Geoscience and Remote Sensing 54(5), 3007-3024.
% http://www2.imm.dtu.dk/pubdb/p.php?6825
% http://orbit.dtu.dk/en/publications/determining-the-points-of-change-in-time-series-of-polarimetric-sar-data(4255ddb4-1fd9-4d66-9aad-e8049a976d75).html
% DOI:10.1109/TGRS.2015.2510160
%
% Allan A. Nielsen, Knut Conradsen, Henning Skriver and Morton J. Canty (2017).
% "Visualization of and software for omnibus test based change detected in
% a time series of polarimetric SAR data".
% Canadian Journal of Remote Sensing 43(6), 582-592.
% http://www2.imm.dtu.dk/pubdb/p.php?id=7027
% DOI:10.1080/07038992.2017.1394182
%
% Allan A. Nielsen, Morton J. Canty, Henning Skriver and Knut Conradsen (2019).
% "Cloud based spatio-temporal analysis of change in sequences of Sentinel
% images." ESA Conference on Big Data from Space (BiDS), 105-108.
% Munich, Germany, 19-21 February 2019.
% http://www2.imm.dtu.dk/pubdb/p.php?id=7129
% DOI:10.2760/848593

% (c) Copyright 2014-2019
% Allan Aasbjerg Nielsen, PhD
% alan@dtu.dk, http://people.compute.dtu.dk/alan
% 6 Dec 2016 - 11 Apr 2019

if nargin>2
    help wishart_change_Rj
    error('wishart_change_Rj: wrong input');
end
if nargout>7
    help wishart_change_Rj
    error('wishart_change_Rj: wrong output');
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
% N = sum(L);

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

fQ = (k-1)*f;
dfQ = fQ/2; % this is 0.5*f; for fewer multiplications in call of gammainc
dfQp4 = dfQ + 2; % this is 0.5*(f+4)

df = f/2; % this is 0.5*f; for fewer multiplications in call of gammainc
dfp4 = df + 2; % this is 0.5*(f+4)

ff = nan(nsets,1);
ff(1) = fQ;
for j=2:nsets
    ff(j) = f;
end

rho = nan(nsets,1);
omega2 = nan(nsets,1);
rho(1) = 1 - (k/n-1/(n*k))/(6*ff(1))*sump;
omega2(1) = -ff(1)/4*(1-1/rho(1))^2 + ...
    (k/n^2-1/n^2/k^2)*sump2/(24*rho(1)^2);
for j=2:nsets
    rho(j) = 1 - (1+1/(j*(j-1)))/(6*ff(j)*n)*sump;
    omega2(j) = -f/4*(1-1/rho(j))^2 + ...
        (1+(2*j-1)/(j^2*(j-1)^2))*sump2/(24*n^2*rho(j)^2);
end

sumMat = nan(nrows,ncols,nv);
lnR = nan(nrows,ncols,nsets);

for j=2:nsets
    cst = p*(j*log(j)-(j-1)*log(j-1));
    sumMat = sum(X(:,:,:,1:(j-1)),4);
    logDet = log(wishart_det(sumMat));
    lnR(:,:,j) = cst + (j-1)*logDet;
    logDet = log(wishart_det(X(:,:,:,j)));
    lnR(:,:,j) = lnR(:,:,j) + logDet;
    sumMat = sumMat + X(:,:,:,j);
    logDet = log(wishart_det(sumMat));
    lnR(:,:,j) = n*(lnR(:,:,j) - j*logDet);
%     aux = n*(lnR(:,:,j) - j*logDet);
%     aux(aux>0) = 0;
%     lnR(:,:,j) = aux;
end

clear X

if nsets==2
    lnR = lnR(:,:,2); % index "2" is not a mistake!
    omega2 = omega2(1);
    rho = rho(1);
    ff = ff(1);
else
    lnR(:,:,1) = sum(lnR(:,:,2:end),3); % this is ln(Q)
%     aux = sum(lnR(:,:,2:end),3); % this is ln(Q)
%     aux(aux>0) = 0;
%     lnR(:,:,1) = aux;
end

if nargout>4
    if nsets==2
        if ~isreal(lnR)
            warning('*** some lnR-values are complex, only real part is used ***')
            fprintf('max(imag(lnR)) = %g\n', max(imag(lnR(:))));
            lnR = real(lnR);
        end
        P = (1 - omega2)*gammainc(-rho*lnR,dfQ  ) + ...
                 omega2 *gammainc(-rho*lnR,dfQp4);
    else
        P = nan(nrows,ncols,nsets);
        for j=[2:nsets 1]
            if j==1, df = dfQ; dfp4 = dfQp4; end
            if ~isreal(lnR(:,:,j))
                aux = imag(lnR(:,:,j));
                warning('*** some lnR-values are complex, only real part is used ***')
                fprintf('max(imag(lnR(:,:,%i))) = %g\n', j, max(aux(:)));
                lnR(:,:,j) = real(lnR(:,:,j));
            end
            P(:,:,j) = (1 - omega2(j))*gammainc(-rho(j)*lnR(:,:,j),df  ) + ...
                            omega2(j) *gammainc(-rho(j)*lnR(:,:,j),dfp4);
        end
    end
    if ~isreal(P) % some lnR/lnQ values may be > 0
        warning('*** some P-values are complex, only real part is output, max(imag(P(:))): ***')
        max(imag(P(:)))
        P = real(P);
    end
end
