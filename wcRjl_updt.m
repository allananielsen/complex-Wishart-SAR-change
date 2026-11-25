function wc = wcRjl_updt(wc,Xnew,L,varargin)

% Wishart based change detection in time series of k polarimetric SAR images
%
% wc = wcRjl_updt(wc,L74,13,'full');
% wcupdt = wcRjl_updt(wc,L74,13,'full');
%
% wcRjl_updt updates the struct wc (which must contain wc.sumMat)
%
% Input
%   wc      - existing struct (from previous runs)
%   Xnew    - new 3-D data set similar to the ones in X (same size, modality
%             i.e. pol etc.)
%             2x2 or 3x3 comples covariance matrices
%             (polarimetric SAR images in covariance matrix form,
%             reading order is alphabetical (as read by freadenvisar):
%               hhhh (real), hhhv (complex), hhvv (complex),
%                            hvhv (real),    hvvv (complex),
%                                            vvvv (real)
%   L       - k (or one) number of looks for X
%   pol     - 'full', 'azim', 'dual', 'fdiag' or 'ddiag'
%             (needed only when confusion is possible in which case it
%             defaults to a value corresponding to the lowest k)
%
% Output
%   wc      - updated structure with
%             lnR     - logarithm of test quantities Q and Rj
%                       ln(Q) is lnR(:,:,1), ln(R_2) is lnR(:,:,2) etc.
%             f       - number of degrees of freedom for Q and Rj test
%                       statistics
%             rho     - aux variable
%             omega2  - aux variable
%             P       - probabilities of finding smaller value of
%                       -2*rho(i)*lnR(:,:,i), i=1,...,k
%             sumMat  - updated sum of matrices
%             logDet  - logarithm of determinants of sumMat
%             L       - equivalent number of looks, ENL (from input)
%
% Reference
%
% Knut Conradsen, Allan Aasbjerg Nielsen and Henning Skriver (2016):
% "Determining the Points of Change in Time Series of Polarimetric SAR Data".
% IEEE Transactions on Geoscience and Remote Sensing 54(5), 3007-3024.
% http://www.imm.dtu.dk/pubdb/p.php?6825
% http://orbit.dtu.dk/en/publications/determining-the-points-of-change-in-time-series-of-polarimetric-sar-data(4255ddb4-1fd9-4d66-9aad-e8049a976d75).html

% (c) Copyright 2016-2021
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, http://people.compute.dtu.dk/alan
% 11 Apr-5 Dec 2016 - 12 Mar 2019
% 28 Aug 2019 /diag only cases: better P 
% 10 Mar 2021: wc.L added

if nargin<3 || nargin>5
    help wcRjl_updt
    error('wcRjl_updt: wrong input');
end
if nargout>1
    help wcRjl_updt
    error('wcRjl_updt: wrong output');
end

if ~isstruct(wc)
    error('wcRjl_updt: first input argument (wc) must be a structure');
end
if ~isfield(wc,'sumMat')
    error('wcRjl_updt: struct wc must contain sumMat')
end
if ~isfield(wc,'logDet')
    error('wcRjl_updt: struct wc must contain logDet')
%     calcDet = true;
end

pX = size(wc(1).sumMat);
% nrows = pX(1);
% ncols = pX(2);
if ismatrix(wc(1).sumMat)
    pX(3) = 1;
end
nv = pX(3);
pY = size(Xnew);
if ismatrix(Xnew)
    pY(3) = 1;
end
if pX(3)~=pY(3)
    error('wcRj_updt: wrong dimensionality of update data');
end
if pY(1)~=pX(1) || pY(2)~=pX(2)
    error('wcRj_updt: existing data in struct wc and update image Xnew must be same size');
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
%     error('wcRjl_updt: input must contain a multplicative of 2, 3, 4, 5 or 9 variables');
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
%         error('wcRjl_updt: number of input frames and polarization do not match')
%     end
end
if nargin>4 % time points for each acquisition
    error('wcRjl_updt: not implemented yet')
end

if ismatrix(Xnew)
    pY(3) = 1;
end
if pY(3)>nv
    error('wcRjl_updt: update by one data set at a time only');
end

% if isscalar(L)
%     L = repmat(L,nsets,1);
%     n = L(1);
% else
%     error('wcRjl_updt: number of looks must be equal for all matrices')
% end
% if length(L)~=nsets
%     error('wcRjl_updt: X and L do not match')
% end

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

nsets = size(wc,2);
k = nsets;
if isscalar(wc(nsets-1).f) % two time points only: copy lnQ to lnR etc.
    wc(nsets-1).lnR(:,:,2) = wc(nsets-1).lnR;
    wc(nsets-1).lnR(:,:,1) = wc(nsets-1).lnR(:,:,2);
    wc(nsets-1).f(2) = wc(nsets-1).f;
    wc(nsets-1).f(1) = wc(nsets-1).f(2);
    wc(nsets-1).rho(2) = wc(nsets-1).rho;
    wc(nsets-1).rho(1) = wc(nsets-1).rho(2);
    wc(nsets-1).omega2(2) = wc(nsets-1).omega2;
    wc(nsets-1).omega2(1) = wc(nsets-1).omega2(2);
    wc(nsets-1).P(:,:,2) = wc(nsets-1).P;
    wc(nsets-1).P(:,:,1) = wc(nsets-1).P(:,:,2);
end
wc(nsets-1).f = wc(nsets-1).f';
wc(nsets-1).rho = wc(nsets-1).rho';
wc(nsets-1).omega2 = wc(nsets-1).omega2';

wcnew = struct('lnR',[],'f',nan,'rho',nan,'omega2',nan,'P',[],...
    'sumMat',[],'logDet',[],'L',nan);
wc = [wc, wcnew];
clear wcnew
i = nsets+1;
LXnew = L(1)*Xnew;
logDetLXnew = log(wishart_det(LXnew));
wc(i).sumMat = LXnew;
clear LXnew
wc(i).logDet = logDetLXnew;
for i=1:nsets
%     L = L(1);
    j = k-i+2;
    if i==nsets, j = k-i+1; end % two time points only: do not copy lnR to lnQ
    [wc(i).lnR(:,:,j),wc(i).f(j),wc(i).rho(j),wc(i).omega2(j),wc(i).P(:,:,j),wc(i).sumMat,wc(i).logDet] = ...
        wishart_change_Rj_updt(wc(i),Xnew,L(1),pol,logDetLXnew);
end
clear logDetLXnew
for i=1:(nsets-1)
%     wc(i).lnR(:,:,1) = sum(wc(i).lnR(:,:,2:end),3); % this is the updated ln(Q)
    wc(i).lnR(:,:,1) = wc(i).lnR(:,:,1) + wc(i).lnR(:,:,end);
end

% omega2, rho f and P for lnQ
n = L(1);
ff = nan(nsets,1);
rhoQ = nan(nsets,1);
omega2Q = nan(nsets,1);

% f = (k-1)*f;
% 
% cst = p*k*log(k);
% df = f/2; % this is 0.5*f; for fewer multiplications in call of gammainc
% dfp4 = df + 2; % this is 0.5*(f+4)

for i=1:nsets
    k = nsets-i+2;
    ff(i) = (k-1)*f;
    rhoQ(i) = 1 - (k/n-1/(k*n))*sump/(6*ff(i));
    omega2Q(i) = -ff(i)/4*(1-1/rhoQ(i))^2 + sump2*(k/n^2-1/n^2/k^2)/(24*rhoQ(i)^2);
%     wc(i).omega2(1) = omega2Q(i);
%     wc(i).rho(1) = rhoQ(i);
end

for i=1:nsets
    wc(i).omega2(1) = omega2Q(i);
    wc(i).rho(1) = rhoQ(i);
end

for i=1:(nsets-1)
    wc(i).f(1) = wc(i).f(1) + f;
end
for i=1:(nsets-1)
    wc(i).P(:,:,1) = (1 - omega2Q(i))*gammainc(-rhoQ(i)*wc(i).lnR(:,:,1), wc(i).f(1)   /2) + ...
                                  omega2Q(i) *gammainc(-rhoQ(i)*wc(i).lnR(:,:,1),(wc(i).f(1)+4)/2);
    if ~isreal(wc(i).P(:,:,1))
        warning('*** some P-values are complex, only real part is output, max(abs(imag(P(:)))): ***')
        aux = wc(i).P(:,:,1);
        max(abs(imag(aux(:))))
        wc(i).P(:,:,1) = real(aux);
    end                              
end

wc(nsets).L = L(1);

