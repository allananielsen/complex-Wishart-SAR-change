function wcQ = wc_Q(X,L,varargin)

% difference/change detection in k (1x1,) 2x2 or 3x3 complex Wishart distributed matrices
% used for change detection in polarimetric SAR images
% 
% wcQ = wc_Q(X,L)
% wcQ = wc_Q(X,L)
% wcQ = wc_Q(X,L,1) % rho=1, omega2=0
%
% Examples:
% wcQ = wc_Q(cat(4,L62,L63,L64,L65,L68,L74),13);
%
% wcQ6 = wc_Q(cat(4,L62,L63,L64,L65,L68,L74),13);
% wcQ5 = wc_Q(cat(4,L63,L64,L65,L68,L74),13);
% wcQ4 = wc_Q(cat(4,L64,L65,L68,L74),13);
% wcQ3 = wc_Q(cat(4,L65,L68,L74),13);
% wcQ2 = wc_Q(cat(4,L68,L74),13);
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
%
% Output: struct with
%   lnQ     - logarithm of test quantity Q
%   f       - number of degrees of freedom
%   rho
%   omega2
%   negdet  - identifies row and col number of pixels with non-positive
%             determinants and positive lnQ
%   P       - probabilities of finding smaller value of -2*rho*lnQ (change prob)
%   L       - k (or one) number of looks for X
%
% see also wishart_change_Q
%
% Reference
%
% Knut Conradsen, Allan Aasbjerg Nielsen and Henning Skriver (2016):
% "Determining the Points of Change in Time Series of Polarimetric SAR Data".
% IEEE Transactions on Geoscience and Remote Sensing 54(5), 3007-3024.
% http://www.imm.dtu.dk/pubdb/p.php?6825
% http://orbit.dtu.dk/en/publications/determining-the-points-of-change-in-time-series-of-polarimetric-sar-data(4255ddb4-1fd9-4d66-9aad-e8049a976d75).html

% (c) Copyright 2012-2025
% Allan Aasbjerg Nielsen, PhD
% alan@dtu.dk, http://people.compute.dtu.dk/alan
%  3 Jun 2019
% 27 Oct 2023 /AA: output is now struct wcQ
%  2 Aug 2024 /AA: two 3x3 in 6x6 covmt and two 2x2 in 4x4 covmat
% 14 Feb 2025 /AA: if k==1 && nv==1: F-test

if nargin>3
    help wc_Q
    error('wc_Q: wrong input');
end
if nargout>1
    help wc_Q
    error('wc_Q: wrong output');
end

pX = size(X);
nrows = pX(1);
ncols = pX(2);
nv = pX(3);
nsets = pX(4);
if nsets<2
    error('wc_Q: one set only? possibly wrong polarization')
end
k = nsets;
% fprintf('%s\n', '*** order of input images is ***very*** important ***')

if isscalar(L)
    L = repmat(L,nsets,1);
    n = L(1);
else
    error('wc_Q: number of looks must be equal for all matrices')
end
if length(L)~=nsets
    error('wc_Q: X and L do not match')
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
elseif nv==9
    p = 3;
    f = p^2;
    sump = p*(2*p^2-1);
    sump2 = p^2*(p^2-1);
elseif nv==18 % two 3x3 blocks in 6x6 covmat
    p1 = 3;
    p2 = 3;
    p = p1+p2;
    f = p1^2+p2^2;
    sump = p1*(2*p1^2-1) + p2*(2*p2^2-1);
    sump2 = p1^2*(p1^2-1) + p2^2*(p2^2-1);
elseif nv==8 % two 2x2 blocks in 4x4 covmat
    p1 = 2;
    p2 = 2;
    p = p1+p2;
    f = p1^2+p2^2;
    sump = p1*(2*p1^2-1) + p2*(2*p2^2-1);
    sump2 = p1^2*(p1^2-1) + p2^2*(p2^2-1);
else
    error('wrong number of input variables')
end
f = (k-1)*f;

cst = p*k*log(k);
df = f/2; % this is 0.5*f; for fewer multiplications in call of gammainc
dfp4 = df + 2; % this is 0.5*(f+4)

rho = 1 - (k/n-1/(k*n))*sump/(6*f);
omega2 = -f/4*(1-1/rho)^2 + sump2*(k/n^2-1/n^2/k^2)/(24*rho^2);
if nargin==3
    rho = 1;
    omega2 = 0;
end

logDet1 = zeros(nrows,ncols);
for i=1:nsets
    logDet1 = logDet1 + log(wishart_det(X(:,:,:,i)));
%     logDet1 = logDet1 + log(wishart_det(squeeze(X(:,:,:,i))));
end
wcQ.lnQ = n*(cst + logDet1 - k*log(wishart_det(sum(X,4))));
clear logDet1
wcQ.f = f; %clear f
wcQ.rho = rho; %clear rho
wcQ.omega2 = omega2; %clear omega2
% if nargout>4
    if ~isreal(wcQ.lnQ)
        warning('*** some lnQ-values are complex, only real part is used ***')
        fprintf('max(imag(lnQ(:))) = %g\n', max(imag(wcQ.lnQ(:))));
        wcQ.lnQ = real(wcQ.lnQ);
    end
    if max(wcQ.lnQ(:))>0
        warning('*** lnQ > 0, bad ***')
        fprintf('max(lnQ(:)) = %g\n', max(wcQ.lnQ(:)));
        if nargout>5
            [wcQ.negdet.rQ,wcQ.negdet.cQ] = find(wcQ.lnQ>0);
            imtool(wcQ.lnQ>0)
        end
    else
        wcQ.negdet = nan;
    end
    wcQ.lnQ(wcQ.lnQ>0) = 0;
    wcQ.P = (1 - omega2)*gammainc(-rho*wcQ.lnQ,df  ) + ...
             omega2 *gammainc(-rho*wcQ.lnQ,dfp4);
    if ~isreal(wcQ.P)
        warning('*** some P-values are complex, only real part is output, max(abs(imag(P(:)))): ***')
        max(abs(imag(wcQ.P(:))))
        wcQ.P = real(wcQ.P);
    end
% end
wcQ.L = L(1);

%% F-test
if k==2 && nv==1
    wcQ.Fstat = X(:,:,:,2) ./ X(:,:,:,1);
    % wcQ.PFstat1 = fcdf(wcQ.Fstat,2*wcQ.L,2*wcQ.L);
    % wcQ.PFstat2 = fcdf(wcQ.Fstat,2*L(1),2*L(2));
    wcQ.PFstat = betainc(L(2)./(L(2)+L(1)*wcQ.Fstat),L(2),L(1)); %,'upper');
    wcQ.comment1 = 'Threshold example: wcQ.PFstat > 0.9995 | wcQ.PFstat < 1-0.9995';
    wcQ.comment2 = '    equivalent to: wcQ.P > 0.999';
end
if false
    figure
    if size(wcQ.Fstat(:))<100*1000-1 % or some other big number
        h = histogram(wcQ.Fstat,'Normalization','pdf');
    else
        h = histogram(wcQ.Fstat,'Normalization','pdf','EdgeColor','none');
    end
    % title("$"+latex(str2sym(txt))+"$",Interpreter="latex") % doesn't work?
    grid on, grid minor
    hold on
    % xx = h.BinLimits(1):0.1:h.BinLimits(2);
    xx = 0:0.1:h.BinLimits(2);
    f = fpdf(xx,2*wcQ.L,2*wcQ.L);
    plot(xx,f,'LineWidth',2)
    hold off
    fprintf('\n*** you may want to use xlim([0 10]) or similar ***\n')
end
